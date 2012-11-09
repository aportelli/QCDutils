#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <qcd_arg_parse.h>
#include <latan/latan_io.h>
#include <latan/latan_mat.h>
#include <latan/latan_math.h>
#include <latan/latan_mass.h>
#include <latan/latan_minimizer.h>
#include <latan/latan_models.h>
#include <latan/latan_plot.h>
#include <latan/latan_rand.h>
#include <latan/latan_statistics.h>

#define MAX_REL_ERR 0.25
#define MAX_TRY 0

typedef enum
{
    PP = 0,\
    AP = 1,\
    AA = 2
} corr_no;

static const strbuf c_name[3] = {"PP","AP","AA"};

static void correct_range(size_t range[2], const mat *prop, const mat *err,\
                          const qcd_options *opt);

static double fm_AA_func(const mat *x, const mat *p, void *vnt)
{
    double res,m,f,t;
    size_t nt;
    
    t = mat_get(x,0,0);
    m = mat_get(p,0,0);
    f = mat_get(p,1,0);
    
    if (vnt)
    {
        nt = *((size_t *)(vnt));
    }
    else
    {
        nt = 0;
    }
    res = 0.5*SQ(f)*m*(exp(-m*t)+exp(-m*(nt-t)));
    
    return res;
}

static double fm_AP_func(const mat *x, const mat *p, void *vnt)
{
    double res,m,f,t,z;
    size_t nt;
    
    t = mat_get(x,0,0);
    m = mat_get(p,0,0);
    f = mat_get(p,1,0);
    z = mat_get(p,2,0);
    
    if (vnt)
    {
        nt = *((size_t *)(vnt));
    }
    else
    {
        nt = 0;
    }
    res = fabs(0.5*f*z*(exp(-m*t)-exp(-m*(nt-t))));
    
    return res;
}

static double fm_PP_func(const mat *x, const mat *p, void *vnt)
{
    double res,m,t,z;
    size_t nt;
    
    t = mat_get(x,0,0);
    m = mat_get(p,0,0);
    z = mat_get(p,2,0);
    
    if (vnt)
    {
        nt = *((size_t *)(vnt));
    }
    else
    {
        nt = 0;
    }
    res = 0.5*z*z/m*(exp(-m*t)+exp(-m*(nt-t)));
    
    return res;
}

static fit_model fm_pseudosc2 =
{
    "combined pseudoscalar PP/AP",
    {&fm_PP_func,&fm_AP_func},
    &npar_3,
    1,
    2
};

static fit_model fm_pseudosc3 =
{
    "combined pseudoscalar PP/AP/AA",
    {&fm_PP_func,&fm_AP_func,&fm_AA_func},
    &npar_3,
    1,
    3
};

static void correct_range(size_t range[2], const mat *prop, const mat *err,\
                          const qcd_options *opt)
{
    size_t inrange[2],nt;
    strbuf buf;
    double rerr;
    
    inrange[0] = range[0];
    inrange[1] = range[1];
    nt         = nrow(prop);
    
    if (inrange[0] <= nt/2)
    {
        strbufcpy(buf,"");
        for (range[1]=range[0];range[1]<=inrange[1];range[1]++)
        {
            rerr = mat_get(err,range[1],0)/fabs(mat_get(prop,range[1],0));
            if (rerr > MAX_REL_ERR)
            {
                strbufcpy(buf," (rcut)");
                break;
            }
        }
        range[1]--;
    }
    else if (inrange[0] >= nt/2)
    {
        range[1] = inrange[1];
        strbufcpy(buf,"");
        for (range[0]=range[1];range[0]>=inrange[0];range[0]--)
        {
            rerr = mat_get(err,range[0],0)/fabs(mat_get(prop,range[0],0));
            if (rerr > MAX_REL_ERR)
            {
                strbufcpy(buf," (lcut)");
                break;
            }
        }
        range[0]++;
    }
    qcd_printf(opt,"[%d,%d]%s ",(int)range[0],(int)range[1],buf);
}

int main(int argc, char* argv[])
{
    /*              parsing arguments           */
    /********************************************/
    double latspac_nu;
    qcd_options *opt;
    strbuf name[3],manf_name,unit,sink[3],source[3];
    size_t binsize,nboot,propdim[2],nt;
    char *cpt1,*cpt2;
    corr_no lastc;
    fit_model *fm_pt;
    
    opt = qcd_arg_parse(argc,argv,A_PLOT|A_SAVE_RS|A_LOAD_RG|A_PROP_NAME\
                        |A_PROP_LOAD|A_LATSPAC|A_QCOMP|A_FIT,2);
    cpt1 = strtok(opt->ss,"/");
    printf("%s\n",cpt1);
    if (cpt1)
    {
        cpt2 = strchr(cpt1,':');
        if (cpt2 == NULL)
        {
            fprintf(stderr,"error: sink/source option %s is invalid\n",\
                    cpt1);
            exit(EXIT_FAILURE);
        }
        strbufcpy(source[PP],cpt2+1);
        *cpt2 = '\0';
        strbufcpy(sink[PP],cpt1);
    }
    else
    {
        fprintf(stderr,"error: sink/source option %s is invalid\n",\
                opt->ss);
        exit(EXIT_FAILURE);
    }
    cpt1 = strtok(NULL,"/");
    if (cpt1)
    {
        cpt2 = strchr(cpt1,':');
        if (cpt2 == NULL)
        {
            fprintf(stderr,"error: sink/source option %s is invalid\n",\
                    cpt1);
            exit(EXIT_FAILURE);
        }
        strbufcpy(source[AP],cpt2+1);
        *cpt2 = '\0';
        strbufcpy(sink[AP],cpt1);
    }
    else
    {
        strbufcpy(sink[AP],sink[PP]);
        strbufcpy(source[AP],source[PP]);
    }
    cpt1 = strtok(NULL,"/");
    if (cpt1)
    {
        cpt2 = strchr(cpt1,':');
        if (cpt2 == NULL)
        {
            fprintf(stderr,"error: sink/source option %s is invalid\n",\
                    cpt1);
            exit(EXIT_FAILURE);
        }
        strbufcpy(source[AA],cpt2+1);
        *cpt2 = '\0';
        strbufcpy(sink[AA],cpt1);
    }
    else
    {
        strbufcpy(sink[AA],sink[PP]);
        strbufcpy(source[AA],source[PP]);
    }
    cpt1 = strtok(opt->channel[0],"/");
    if (cpt1)
    {
        sprintf(name[PP],"%s_%s_%s_%s",cpt1,opt->quark[0],sink[PP],source[PP]);
    }
    else
    {
        fprintf(stderr,"error: channel option %s is invalid\n",\
                opt->channel[0]);
        exit(EXIT_FAILURE);
    }
    cpt1 = strtok(NULL,"/");
    if (cpt1)
    {
        sprintf(name[AP],"%s_%s_%s_%s",cpt1,opt->quark[0],sink[AP],source[AP]);
    }
    else
    {
        fprintf(stderr,"error: channel option %s is invalid\n",\
                opt->channel[0]);
        exit(EXIT_FAILURE);
    }
    cpt1 = strtok(NULL,"/");
    if (cpt1)
    {
        sprintf(name[AA],"%s_%s_%s_%s",cpt1,opt->quark[0],sink[AA],source[AA]);
        lastc = AA;
        fm_pt = &fm_pseudosc3;
    }
    else
    {
        lastc = AP;
        fm_pt = &fm_pseudosc2;
    }
    strbufcpy(manf_name,opt->manf_name);
    binsize    = opt->binsize;
    nboot      = opt->nboot;
    latspac_nu = opt->latspac_nu;
    if (opt->have_latspac)
    {
        strbufcpy(unit," (MeV)");
    }
    else
    {
        strbufcpy(unit,"");
    }
    latan_set_verb(opt->latan_verb);
    minimizer_set_alg(opt->minimizer);
    mat_ar_loadbin(NULL,propdim,manf_name,name[PP],1);
    nt = propdim[0];
    io_set_fmt(opt->latan_fmt);
    io_init();
    
    /*              loading datas               */
    /********************************************/
    size_t ndat,nbdat;
    mat **prop[3];
    corr_no c;

    ndat    = (size_t)get_nfile(manf_name);
    nbdat   = ndat/binsize + ((ndat%binsize == 0) ? 0 : 1);
    
    for (c=PP;c<=lastc;c++)
    {
        prop[c] = mat_ar_create(nbdat,propdim[0],propdim[1]);
        qcd_printf(opt,"-- loading %s datas from %s...\n",name[c],manf_name);
        mat_ar_loadbin(prop[c],NULL,manf_name,name[c],binsize);
    }
    
    /*                propagator                */
    /********************************************/
    rs_sample *s_mprop[3];
    mat *mprop[3],*sigmprop[3];
    
    for (c=PP;c<=lastc;c++)
    {
        s_mprop[c]  = rs_sample_create(propdim[0],propdim[1],nboot);
        sigmprop[c] = mat_create(propdim[0],propdim[1]);
        qcd_printf(opt,"-- resampling %s mean propagator...\n",name[c]);
        randgen_set_state(opt->state);
        resample(s_mprop[c],prop[c],nbdat,&rs_mean,BOOT,NULL);
        mprop[c] = rs_sample_pt_cent_val(s_mprop[c]);
        rs_sample_varp(sigmprop[c],s_mprop[c]);
        mat_eqsqrt(sigmprop[c]);
    }
    
    /*           effective mass                 */
    /********************************************/
    rs_sample *s_effmass[3];
    mat *tem,*em[3],*sigem[3];
    size_t emdim[2];
    
    get_effmass_size(emdim,mprop[PP],1,EM_ACOSH);
    
    tem = mat_create(emdim[0],1);
    
    for (c=PP;c<=lastc;c++)
    {
        s_effmass[c] = rs_sample_create(emdim[0],emdim[1],nboot);
        sigem[c]     = mat_create(emdim[0],emdim[1]);
        qcd_printf(opt,"-- resampling %s effective mass...\n",name[c]);
        rs_sample_effmass(s_effmass[c],tem,s_mprop[c],1,EM_ACOSH);
        em[c] = rs_sample_pt_cent_val(s_effmass[c]);
        rs_sample_varp(sigem[c],s_effmass[c]);
        mat_eqsqrt(sigem[c]);
    }
    
    /*                  fit mass                */
    /********************************************/
    fit_data *d;
    rs_sample *s_fit;
    mat *fit,*limit,*sigfit,*scanres_t,*scanres_chi2,*scanres_mass,\
    *scanres_masserr;
    size_t npar,nti,tibeg,range[2],ta;
    size_t i,j;
    strbuf buf,range_info,latan_path;
    double pref_i,mass_i;
    
    d        = fit_data_create(nt,1,(size_t)(lastc+1));
    tibeg    = (size_t)(opt->range[0][0]);
    range[0] = 0;
    range[1] = 0;
    npar     = fit_model_get_npar(fm_pt,&nt);
    
    s_fit    = rs_sample_create(npar,1,nboot);
    fit      = rs_sample_pt_cent_val(s_fit);
    limit    = mat_create(npar,2);
    sigfit   = mat_create(npar,1);
    
    /** print operation **/
    if (!opt->do_range_scan)
    {
        qcd_printf(opt,"-- fitting and resampling...\n");
    }
    else
    {
        qcd_printf(opt,"-- scanning ranges [ti,%u] from ti= %u\n",\
                   opt->range[0][1],opt->range[0][0]);
        opt->nmanrange   = 1;
    }
    
    /** check ranges **/
    strbufcpy(range_info,"");
    qcd_printf(opt,"%-20s: ","corrected range(s)");
    fit_data_fit_all_points(d,false);
    for (i=0;i<opt->nmanrange;i++)
    {
        sprintf(buf,"_%u_%u",opt->range[i][0],opt->range[i][1]);
        strbufcat(range_info,buf);
        range[0] = (size_t)(opt->range[i][0]);
        range[1] = (size_t)(opt->range[i][1]);
        correct_range(range,mprop[PP],sigmprop[PP],opt);
        fit_data_fit_range(d,range[0],range[1],true);
    }
    qcd_printf(opt,"\n");
    
    nti             = MAX(range[1] - 1 - tibeg,1);
    scanres_t       = mat_create(nti,1);
    scanres_chi2    = mat_create(nti,1);
    scanres_mass    = mat_create(nti,1);
    scanres_masserr = mat_create(nti,1);
    
    /** set model **/
    fit_data_set_model(d,fm_pt,&nt);
    
    /** set correlation filter **/
    if (opt->corr == NO_COR)
    {
        for (i=0;i<nt;i++)
            for (j=0;j<nt;j++)
            {
                if (i != j)
                {
                    fit_data_set_data_cor(d,i,j,false);
                }
            }
    }
    
    /** set initial parameter values **/
    ta     = nt/8-(size_t)(mat_get(tem,0,0));
    mass_i = mat_get(em[PP],ta,0);
    if (latan_isnan(mass_i))
    {
        mass_i = 0.3;
    }
    mat_set(fit,0,0,mass_i);
    pref_i = mat_get(mprop[PP],ta,0)/(exp(-mass_i*ta)+exp(-mass_i*(nt-ta)));
    if (latan_isnan(pref_i))
    {
        pref_i = 1.0;
    }
    mat_set(fit,1,0,sqrt(pref_i));
    mat_set(fit,2,0,sqrt(pref_i));
    qcd_printf(opt,"%-22smass= %e -- prefactor_0= %e\n","initial parameters: ",\
               mat_get(fit,0,0),pref_i);
    
    /** set parameter limits **/
    mat_cst(limit,latan_nan());
    mat_set(limit,0,0,0.0);
    
    /** positive AP correlator **/
    rs_sample_eqabs(s_mprop[AP]);
    
    /** set x **/
    for (i=0;i<nt;i++)
    {
        fit_data_set_x(d,i,0,(double)(i));
    }
    
    /** regular correlator fit... **/
    if (!opt->do_range_scan)
    {
        latan_set_warn(false);
        rs_data_fit(s_fit,limit,NULL,s_mprop,d,NO_COR,NULL);
        latan_set_warn(true);
        rs_data_fit(s_fit,limit,NULL,s_mprop,d,opt->corr,NULL);
        rs_sample_varp(sigfit,s_fit);
        mat_eqsqrt(sigfit);
        if (fit_data_get_chi2pdof(d) > 2.0)
        {
            fprintf(stderr,"warning: bad final fit (chi^2/dof= %.2e)\n",\
                    fit_data_get_chi2pdof(d));
        }
        qcd_printf(opt,"-- results:\n");
        qcd_printf(opt,"%-10s= %.8f +/- %.8e %s\n","mass",\
                   mat_get(fit,0,0)/latspac_nu,       \
                   mat_get(sigfit,0,0)/latspac_nu,unit);
        qcd_printf(opt,"%-10s= %.8f +/- %.8e %s\n","decay",\
                   mat_get(fit,1,0)/latspac_nu,       \
                   mat_get(sigfit,1,0)/latspac_nu,unit);
        qcd_printf(opt,"%-10s= %d\n","dof",fit_data_get_dof(d));
        qcd_printf(opt,"%-10s= %e\n","chi^2/dof",fit_data_get_chi2pdof(d));
        if (opt->do_save_rs_sample)
        {
            sprintf(latan_path,"%s_pseudosc_fit%s_%s.boot:%s_pseudosc_fit%s_%s",\
                    opt->quark[0],range_info,manf_name,opt->quark[0],\
                    range_info,manf_name);
            rs_sample_save_subsamp(latan_path,'w',s_fit,0,0,1,0);
        }
    }
    /** ...or fit range scanning **/
    else
    {
        qcd_printf(opt,"\n%-5s %-12s a*M_%-8s %-12s","ti/a","chi^2/dof",\
                   opt->quark[0],"error");
        for (i=tibeg;i<range[1]-1;i++)
        {
            latan_set_warn(false);
            rs_data_fit(s_fit,limit,NULL,s_mprop,d,NO_COR,NULL);
            latan_set_warn(true);
            rs_data_fit(s_fit,limit,NULL,s_mprop,d,opt->corr,NULL);
            rs_sample_varp(sigfit,s_fit);
            mat_eqsqrt(sigfit);
            mat_set(scanres_t,i-tibeg,0,(double)(i));
            mat_set(scanres_chi2,i-tibeg,0,fit_data_get_chi2pdof(d));
            mat_set(scanres_mass,i-tibeg,0,mat_get(fit,0,0));
            mat_set(scanres_masserr,i-tibeg,0,mat_get(sigfit,0,0));
            qcd_printf(opt,"\n% -4d % -.5e % -.5e % -.5e",(int)(i),\
                       fit_data_get_chi2pdof(d),mat_get(fit,0,0), \
                       mat_get(sigfit,0,0));
            fit_data_fit_point(d,i,false);
        }
        qcd_printf(opt,"\n\n");
    }
    
    /*                  plot                    */
    /********************************************/
    if (opt->do_plot)
    {
        mat *pr_t,*mbuf,*em_i,*sigem_i,*par,*ft[3],*comp[3];
        plot *p;
        strbuf key,dirname,color;
        size_t maxt,t,npoint;
        double dmaxt,nmass;
        
        pr_t    = mat_create(nt,1);
        mbuf    = mat_create(1,1);
        em_i    = mat_create(nrow(em[PP]),1);
        sigem_i = mat_create(nrow(em[PP]),1);
        par     = mat_create(2,1);
        
        if (!opt->do_range_scan)
        {
            maxt   = nt;
            dmaxt  = (double)maxt;
            npoint = fit_data_fit_point_num(d);
            
            /** chi^2 plot **/
            p = plot_create();
            i = 0;
            for (c=PP;c<=lastc;c++)
            {
                ft[c]   = mat_create(npoint,1);
                comp[c] = mat_create(npoint,1);
                for (t=0;t<nt;t++)
                {
                    if (fit_data_is_fit_point(d,t))
                    {
                        mat_set(ft[c],i%npoint,0,(double)(t)+0.33*(double)(c));
                        mat_set(comp[c],i%npoint,0,mat_get(d->chi2_comp,i,0));
                        i++;
                    }
                }
            }
            plot_set_scale_manual(p,-1.0,dmaxt,-5.0,5.0);
            plot_add_plot(p,"0.0 lt -1 lc rgb 'black' notitle","");
            plot_add_plot(p,"1.0 lt -1 lc rgb 'black' notitle","");
            plot_add_plot(p,"-1.0 lt -1 lc rgb 'black' notitle","");
            plot_add_plot(p,"2.0 lt -1 lc rgb 'dark-gray' notitle","");
            plot_add_plot(p,"-2.0 lt -1 lc rgb 'dark-gray' notitle","");
            plot_add_plot(p,"3.0 lt -1 lc rgb 'gray' notitle","");
            plot_add_plot(p,"-3.0 lt -1 lc rgb 'gray' notitle","");
            plot_add_plot(p,"4.0 lt -1 lc rgb 'light-gray' notitle","");
            plot_add_plot(p,"-4.0 lt -1 lc rgb 'light-gray' notitle","");
            plot_set_ylabel(p,"standard deviations");
            for (c=PP;c<=lastc;c++)
            {
                sprintf(color,"%d",(int)(c)+1);
                plot_add_points(p,ft[c],comp[c],c_name[c],color,"impulses");
            }
            plot_disp(p);
            if (opt->do_save_plot)
            {
                sprintf(dirname,"%s_dev",opt->save_plot_dir);
                plot_save(dirname,p);
            }
            plot_destroy(p);
        
            for (c=PP;c<=lastc;c++)
            {
                /** propagator plot **/
                p = plot_create();
                fit_data_fit_all_points(d,true);
                plot_set_scale_ylog(p);
                plot_set_scale_xmanual(p,0,dmaxt);
                sprintf(key,"%s %s propagator",opt->quark[0],c_name[c]);
                mat_eqabs(mprop[c]);
                mat_set_step(pr_t,0.0,1.0);
                plot_add_fit(p,d,c,mbuf,0,fit,0,dmaxt,1000,false,\
                             PF_FIT|PF_DATA,key,"","rgb 'red'","rgb 'red'");
                plot_disp(p);
                if (opt->do_save_plot)
                {
                    sprintf(dirname,"%s_prop_%s",opt->save_plot_dir,c_name[c]);
                    plot_save(dirname,p);
                }
                plot_destroy(p);
                
                /** effective mass plot **/
                p = plot_create();
                mat_eqmuls(em[c],1.0/latspac_nu);
                mat_eqmuls(sigem[c],1.0/latspac_nu);
                nmass = mat_get(fit,0,0)/latspac_nu;
                plot_add_hlineerr(p,nmass, mat_get(sigfit,0,0)/latspac_nu,\
                                  "rgb 'red'");
                sprintf(key,"%s %s effective mass",opt->quark[0],c_name[c]);
                plot_add_dat(p,tem,em[c],NULL,sigem[c],key,"rgb 'blue'");
                plot_disp(p);
                if (opt->do_save_plot)
                {
                    sprintf(dirname,"%s_em_%s",opt->save_plot_dir,c_name[c]);
                    plot_save(dirname,p);
                }
                plot_destroy(p);
            }
        }
        else
        {
            /* chi^2 plot */
            p = plot_create();
            plot_set_scale_manual(p,0,(double)(nt/2),0,5.0);
            plot_add_hline(p,1.0,"rgb 'black'");
            plot_add_dat(p,scanres_t,scanres_chi2,NULL,NULL,"chi^2/dof",\
                         "rgb 'blue'");
            plot_disp(p);
            if (opt->do_save_plot)
            {
                sprintf(dirname,"%s_chi2",opt->save_plot_dir);
                plot_save(dirname,p);
            }
            plot_destroy(p);
            
            /* mass plot */
            p = plot_create();
            plot_set_scale_xmanual(p,0,(double)(nt/2));
            sprintf(key,"a*M_%s",opt->quark[0]);
            plot_add_dat(p,scanres_t,scanres_mass,NULL,scanres_masserr,key,\
                         "rgb 'red'");
            plot_disp(p);
            if (opt->do_save_plot)
            {
                sprintf(dirname,"%s_mass",opt->save_plot_dir);
                plot_save(dirname,p);
            }
            plot_destroy(p);
        }
        
        mat_destroy(em_i);
        mat_destroy(sigem_i);
        mat_destroy(pr_t);
        mat_destroy(mbuf);
        mat_destroy(par);
        for (c=PP;c<=lastc;c++)
        {
            mat_destroy(ft[c]);
            mat_destroy(comp[c]);
        }
    }
    
    /*              desallocation               */
    /********************************************/
    free(opt);
    io_finish();
    mat_ar_destroy(prop[0],nbdat);
    mat_ar_destroy(prop[1],nbdat);
    for (c=PP;c<=lastc;c++)
    {
        rs_sample_destroy(s_mprop[c]);
        mat_destroy(sigmprop[c]);
        rs_sample_destroy(s_effmass[c]);
        mat_destroy(sigem[c]);
    }
    mat_destroy(tem);
    fit_data_destroy(d);
    rs_sample_destroy(s_fit);
    mat_destroy(limit);
    mat_destroy(sigfit);
    mat_destroy(scanres_t);
    mat_destroy(scanres_chi2);
    mat_destroy(scanres_mass);
    mat_destroy(scanres_masserr);
    
    return EXIT_SUCCESS;
}
