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

static void correct_range(size_t range[2], const mat *prop, const mat *err,\
                          const qcd_options *opt);

static void correct_range(size_t range[2], const mat *prop, const mat *err,\
                          const qcd_options *opt)
{
    size_t inrange[2],nt;
    strbuf buf;
    double rerr;
    
    inrange[0] = range[0];
    inrange[1] = range[1];
    nt         = (opt->nt) ? opt->nt : nrow(prop);
    
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
    bool is_split;
    int emtype;
    strbuf prop_name[2],full_name,manf_name,model,unit;
    size_t binsize,nboot,propdim[2],nt,nstate;
    fit_model *fm_pt;
    void *fmpar_pt;
    
    opt = qcd_arg_parse(argc,argv,A_PLOT|A_SAVE_RS|A_LOAD_RS|A_LOAD_RG  \
                        |A_PROP_NAME|A_PROP_LOAD|A_LATSPAC|A_QCOMP|A_FIT\
                        |A_MODEL,2,2);
    is_split = (strlen(opt->channel[1]) != 0);
    
    if (is_split)
    {
        sprintf(full_name,"%s_%s-%s_%s",opt->channel[0],opt->quark[0],\
                opt->channel[1],opt->quark[1]);
    }
    else
    {
        sprintf(full_name,"%s_%s",opt->channel[0],opt->quark[0]);
    }
    sprintf(prop_name[0],"%s_%s_%s_%s",opt->channel[0],opt->quark[0],opt->sink,\
            opt->source);
    if (is_split)
    {
        sprintf(prop_name[1],"%s_%s_%s_%s",opt->channel[1],opt->quark[1],\
                opt->sink,opt->source);
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
    strbufcpy(model,opt->model);
    minimizer_set_alg(opt->minimizer);
    if (opt->do_load_rs_sample)
    {
        rs_sample_load(NULL,NULL,propdim,opt->load_rs_fname[0]);
    }
    else
    {
        mat_ar_loadbin(NULL,propdim,manf_name,prop_name[0],1);
    }
    nt = (opt->nt) ? opt->nt : propdim[0];
    if (strbufcmp(model,"cosh") == 0)
    {
        fm_pt    = &fm_cosh;
        fmpar_pt = &nt;
        emtype   = EM_ACOSH;
        nstate   = 1;
    }
    else if (strbufcmp(model,"cosh_ex") == 0)
    {
        fm_pt    = &fm_cosh_ex;
        fmpar_pt = &nt;
        emtype   = EM_ACOSH;
        nstate   = 2;
    }
    else if (strbufcmp(model,"expdec_ex") == 0)
    {
        fm_pt    = &fm_expdec_ex;
        fmpar_pt = NULL;
        emtype   = EM_LOG;
        nstate   = 2;
    }
    else
    {
        strbufcpy(opt->model,"expdec");
        fm_pt    = &fm_expdec;
        fmpar_pt = NULL;
        emtype   = EM_LOG;
        nstate   = 1;
    }
    io_set_fmt(opt->latan_fmt);
    io_init();
    
    /*              loading datas               */
    /********************************************/
    mat **prop[2];
    size_t ndat,nbdat;
    
    if (opt->do_load_rs_sample)
    {
        ndat    = 0;
        nbdat   = 0;
        prop[0] = NULL;
        prop[1] = NULL;
    }
    else
    {
        ndat    = (size_t)get_nfile(manf_name);
        nbdat   = ndat/binsize + ((ndat%binsize == 0) ? 0 : 1);
        
        prop[0] = mat_ar_create(nbdat,propdim[0],propdim[1]);
        prop[1] = mat_ar_create(nbdat,propdim[0],propdim[1]);
        
        qcd_printf(opt,"-- loading %s datas from %s...\n",prop_name[0],manf_name);
        mat_ar_loadbin(prop[0],NULL,manf_name,prop_name[0],binsize);
        if (is_split)
        {
            qcd_printf(opt,"-- loading %s datas from %s...\n",prop_name[1],\
                       manf_name);
            mat_ar_loadbin(prop[1],NULL,manf_name,prop_name[1],binsize);
        }
    }
    
    /*                propagator                */
    /********************************************/
    rs_sample *s_tmp,*s_mprop;
    mat *mprop,*sigmprop;

    s_mprop  = rs_sample_create(propdim[0],propdim[1],nboot);
    s_tmp    = rs_sample_create(propdim[0],propdim[1],nboot);
    sigmprop = mat_create(propdim[0],propdim[1]);
    
    if (opt->do_load_rs_sample)
    {
        qcd_printf(opt,"-- loading %s propagator from %s...\n",prop_name[0],\
                   opt->load_rs_fname[0]);
        rs_sample_load(s_mprop,NULL,NULL,opt->load_rs_fname[0]);
        if (is_split)
        {
            qcd_printf(opt,"-- loading %s propagator from %s...\n",\
                       prop_name[1],opt->load_rs_fname[1]);
            rs_sample_load(s_tmp,NULL,NULL,opt->load_rs_fname[1]);
            rs_sample_eqdivp(s_mprop,s_tmp);
        }
    }
    else
    {
        qcd_printf(opt,"-- resampling %s mean propagator...\n",prop_name[0]);
        randgen_set_state(opt->state);
        resample(s_mprop,prop[0],nbdat,&rs_mean,BOOT,NULL);
        if (is_split)
        {
            qcd_printf(opt,"-- resampling %s mean propagator...\n",prop_name[1]);
            randgen_set_state(opt->state);
            resample(s_tmp,prop[1],nbdat,&rs_mean,BOOT,NULL);
            rs_sample_eqdivp(s_mprop,s_tmp);
        }
    }
    mprop = rs_sample_pt_cent_val(s_mprop);
    rs_sample_varp(sigmprop,s_mprop);
    mat_eqsqrt(sigmprop);
        
    /*           effective mass                 */
    /********************************************/
    rs_sample *s_effmass;
    mat *tem,*em,*sigem;
    size_t emdim[2];
    
    get_effmass_size(emdim,mprop,1,emtype);
    
    s_effmass = rs_sample_create(emdim[0],emdim[1],nboot);
    tem       = mat_create(emdim[0],1);
    sigem     = mat_create(emdim[0],emdim[1]);
    
    qcd_printf(opt,"-- resampling %s effective mass...\n",full_name);
    rs_sample_effmass(s_effmass,tem,s_mprop,1,emtype);
    em = rs_sample_pt_cent_val(s_effmass);
    rs_sample_varp(sigem,s_effmass);
    mat_eqsqrt(sigem);
    
    /*                  fit mass                */
    /********************************************/
    fit_data *d;
    rs_sample *s_mass;
    mat *mass,*limit,*sigmass,*scanres_t,*scanres_chi2,*scanres_mass,\
        *scanres_masserr;
    size_t npar,nti,ta,tibeg,range[2];
    size_t i;
    strbuf buf,range_info,latan_path;
    double pref_i,mass_i;

    d        = fit_data_create(propdim[0],1,1);
    tibeg    = (size_t)(opt->range[0][0]); 
    range[0] = 0;
    range[1] = 0;
    npar     = fit_model_get_npar(fm_pt,fmpar_pt);
    
    s_mass   = rs_sample_create(npar,1,nboot);
    mass     = rs_sample_pt_cent_val(s_mass);
    limit    = mat_create(npar,2);
    sigmass  = mat_create(npar,1);

    /** print operation **/
    if (!opt->do_range_scan)
    {
        qcd_printf(opt,"-- fitting and resampling %s mass...\n",full_name);
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
        correct_range(range,mprop,sigmprop,opt);
        fit_data_fit_range(d,range[0],range[1],true);
    }
    qcd_printf(opt,"\n");
    
    nti             = MAX(range[1] - 1 - tibeg,1);
    scanres_t       = mat_create(nti,1);
    scanres_chi2    = mat_create(nti,1);
    scanres_mass    = mat_create(nti,1);
    scanres_masserr = mat_create(nti,1);

    /** set model **/
    fit_data_set_model(d,fm_pt,fmpar_pt);
    
    /** set initial parameter values **/
    ta     = nt/8-(size_t)(mat_get(tem,0,0));
    mass_i = mat_get(em,ta,0);
    if (latan_isnan(mass_i))
    {
        mass_i = 0.1;
    }
    for (i=0;i<nstate;i++)
    {
        mat_set(mass,i,0,((double)(i+1))*mass_i);
    }
    pref_i = log(mat_get(mprop,ta,0)/(exp(-mass_i*ta)+exp(-mass_i*(nt-ta))));
    if (latan_isnan(pref_i))
    {
        pref_i = 1.0;
    }
    pref_i -= log((double)(nstate));
    for (i=nstate;i<npar;i++)
    {
        mat_set(mass,i,0,pref_i);
    }
    qcd_printf(opt,"%-22smass= %e -- prefactor_0= %e\n","initial parameters: ",\
               mat_get(mass,0,0),exp(pref_i));
    for (i=1;i<nstate;i++)
    {
        qcd_printf(opt,"%-22sE_%d= %e -- prefactor_%d= %e\n","",(int)(i),\
                   mat_get(mass,i,0),(int)(i),exp(pref_i));
    }
    
    /** set parameter limits **/
    mat_cst(limit,latan_nan());
    if (!is_split)
    {
        mat_set(limit,0,0,0.0);
        for (i=1;i<nstate;i++)
        {
            mat_set(limit,i,0,mat_get(mass,i-1,0));
        }
    }
    
    /** set x **/
    for (i=0;i<propdim[0];i++)
    {
        fit_data_set_x(d,i,0,(double)(i)+opt->tshift);
    }
    
    /** regular correlator fit... **/
    if (!opt->do_range_scan)
    {
        latan_set_warn(false);
        rs_data_fit(s_mass,limit,NULL,&s_mprop,d,NO_COR,NULL);
        latan_set_warn(true);
        rs_data_fit(s_mass,limit,NULL,&s_mprop,d,opt->corr,NULL);
        if (nstate == 2)
        {
            if (fabs(mat_get(mass,0,0)) > fabs(mat_get(mass,1,0)))
            {
                fprintf(stderr,"warning: mass and excited state may have flip during fit, check result.\n");
                rs_sample *s_buf[2];
                
                s_buf[0] = rs_sample_create(1,1,nboot);
                s_buf[1] = rs_sample_create(1,1,nboot);
                rs_sample_get_subsamp(s_buf[0],s_mass,0,0,0,0);
                rs_sample_get_subsamp(s_buf[1],s_mass,1,0,1,0);
                rs_sample_set_subsamp(s_mass,s_buf[0],1,0,1,0);
                rs_sample_set_subsamp(s_mass,s_buf[1],0,0,0,0);
                rs_sample_get_subsamp(s_buf[0],s_mass,2,0,2,0);
                rs_sample_get_subsamp(s_buf[1],s_mass,3,0,3,0);
                rs_sample_set_subsamp(s_mass,s_buf[0],3,0,3,0);
                rs_sample_set_subsamp(s_mass,s_buf[1],2,0,2,0);
                rs_sample_destroy(s_buf[0]);
                rs_sample_destroy(s_buf[1]);
            }
        }
        rs_sample_varp(sigmass,s_mass);
        mat_eqsqrt(sigmass);
        if (fit_data_get_chi2pdof(d) > 2.0)
        {
            fprintf(stderr,"warning: bad final fit (chi^2/dof= %.2e)\n",\
                    fit_data_get_chi2pdof(d));
        }
        qcd_printf(opt,"-- results:\n");
        qcd_printf(opt,"%-10s= %.8f +/- %.8e %s\n","mass",\
                   mat_get(mass,0,0)/latspac_nu,       \
                   mat_get(sigmass,0,0)/latspac_nu,unit);
        for (i=1;i<nstate;i++)
        {
            sprintf(buf,"E_%d",(int)(i));
            qcd_printf(opt,"%-10s= %.8f +/- %.8e %s\n",buf, \
                       mat_get(mass,i,0)/latspac_nu,     \
                       mat_get(sigmass,i,0)/latspac_nu,unit);
        }
        qcd_printf(opt,"%-10s= %d\n","dof",fit_data_get_dof(d));
        qcd_printf(opt,"%-10s= %e\n","chi^2/dof",fit_data_get_chi2pdof(d));
        if (opt->do_save_rs_sample)
        {
            sprintf(latan_path,"%s_mass_fit%s_%s.boot:%s_mass_fit%s_%s",\
                    full_name,range_info,manf_name,full_name,range_info, \
                    manf_name);
            rs_sample_save_subsamp(latan_path,'w',s_mass,0,0,0,0);
        }
    }
    /** ...or fit range scanning **/
    else
    {
        qcd_printf(opt,"\n%-5s %-12s a*M_%-8s %-12s","ti/a","chi^2/dof",\
                   full_name,"error");
        for (i=tibeg;i<range[1]-1;i++)
        {
            latan_set_warn(false);
            rs_data_fit(s_mass,limit,NULL,&s_mprop,d,NO_COR,NULL);
            latan_set_warn(true);
            rs_data_fit(s_mass,limit,NULL,&s_mprop,d,opt->corr,NULL);
            rs_sample_varp(sigmass,s_mass);
            mat_eqsqrt(sigmass);
            mat_set(scanres_t,i-tibeg,0,(double)(i));
            mat_set(scanres_chi2,i-tibeg,0,fit_data_get_chi2pdof(d));
            mat_set(scanres_mass,i-tibeg,0,mat_get(mass,0,0));
            mat_set(scanres_masserr,i-tibeg,0,mat_get(sigmass,0,0));
            qcd_printf(opt,"\n% -4d % -.5e % -.5e % -.5e",(int)(i),\
                       fit_data_get_chi2pdof(d),mat_get(mass,0,0), \
                       mat_get(sigmass,0,0));
            fit_data_fit_point(d,i,false);
        }
        qcd_printf(opt,"\n\n");
    }
    
    /*                  plot                    */
    /********************************************/
    if (opt->do_plot|opt->do_save_plot)
    {
        mat *mbuf,*em_i,*sigem_i,*prop_pt,*mass_pt,*par,*ft,*comp;
        plot *p;
        strbuf key,xlabel,dirname;
        size_t maxt,j,t;
        double dmaxt,nmass,corr_prop;
        
        mbuf    = mat_create(1,1);
        em_i    = mat_create(nrow(em),1);
        sigem_i = mat_create(nrow(em),1);
        par     = mat_create(2,1);
        ft      = mat_create(fit_data_fit_point_num(d),1);
        comp    = mat_create(fit_data_fit_point_num(d),1);
        
        sprintf(xlabel,"a*t%+.2f",opt->tshift);
        if (emtype == EM_ACOSH)
        {
            maxt = MIN(propdim[0],nt);
        }
        else if (emtype == EM_LOG)
        {
            maxt = MIN(propdim[0],nt/2);
        }
        else
        {
            maxt = propdim[0];
        }
        dmaxt = (double)maxt;
        if (!opt->do_range_scan)
        {
            /** chi^2 plot **/
            p = plot_create();
            i = 0;
            for (t=0;t<propdim[0];t++)
            {
                if (fit_data_is_fit_point(d,t))
                {
                    mat_set(ft,i,0,(double)(t)+opt->tshift);
                    mat_set(comp,i,0,mat_get(d->chi2_comp,i,0));
                    i++;
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
            plot_set_xlabel(p,xlabel);
            plot_set_ylabel(p,"standard deviations");
            plot_add_points(p,ft,comp,"","rgb 'red'","impulses");
            if (opt->do_plot)
            {
                plot_disp(p);
            }
            if (opt->do_save_plot)
            {
                sprintf(dirname,"%s_dev",opt->save_plot_dir);
                plot_save(dirname,p);
            }
            plot_destroy(p);
            
            /** propagator plot **/
            p = plot_create();
            fit_data_fit_all_points(d,true);
            if (is_split)
            {
                plot_set_scale_ymanual(p,0.92,1.08);
            }
            else
            {
                plot_set_scale_ylog(p);
            }
            plot_set_scale_xmanual(p,0,dmaxt);
            plot_set_xlabel(p,xlabel);
            sprintf(key,"%s propagator",full_name);
            mat_eqabs(mprop);
            plot_add_fit(p,d,0,mbuf,0,mass,0,dmaxt,1000,false,\
                         PF_FIT|PF_DATA,key,"","rgb 'red'","rgb 'red'");
            if (opt->do_plot)
            {
                plot_disp(p);
            }
            if (opt->do_save_plot)
            {
                sprintf(dirname,"%s_prop",opt->save_plot_dir);
                plot_save(dirname,p);
            }
            plot_destroy(p);
            
            /** effective mass plot **/
            p = plot_create();
            mat_eqmuls(em,1.0/latspac_nu);
            mat_eqmuls(sigem,1.0/latspac_nu);
            nmass = mat_get(mass,0,0)/latspac_nu;
            plot_add_hlineerr(p,nmass, mat_get(sigmass,0,0)/latspac_nu,\
                              "rgb 'red'");
            if (is_split)
            {
                plot_add_hline(p,0.0,"rgb 'black'");
            }
            for (i=1;i<nstate;i++)
            {
                nmass = mat_get(mass,i,0)/latspac_nu;
                plot_add_hlineerr(p,nmass,mat_get(sigmass,i,0)/latspac_nu,\
                                  "rgb 'red'");
            }
            plot_set_scale_manual(p,0.0,dmaxt,0.0,2.0*nmass);
            plot_set_xlabel(p,xlabel);
            sprintf(key,"%s effective mass",full_name);
            plot_add_dat(p,tem,em,NULL,sigem,key,"rgb 'blue'");
            plot_add_dat(p,tem,em,NULL,NULL,"","rgb 'blue'");
            for (i=1;i<nstate;i++)
            {
                if (emtype == EM_ACOSH)
                {
                    fm_pt = &fm_cosh;
                }
                else if (emtype == EM_LOG)
                {
                    fm_pt = &fm_expdec;
                }
                for (j=0;j<=rs_sample_get_nsample(s_mprop);j++)
                {
                    prop_pt = (j == 0) ? rs_sample_pt_cent_val(s_mprop)\
                                       : rs_sample_pt_sample(s_mprop,j-1);
                    mass_pt = (j == 0) ? rs_sample_pt_cent_val(s_mass)\
                                       : rs_sample_pt_sample(s_mass,j-1);
                    mat_set(par,0,0,mat_get(mass_pt,i-1,0));
                    mat_set(par,1,0,mat_get(mass_pt,nstate+i-1,0));
                    for (t=0;t<propdim[0];t++)
                    {
                        mat_set(mbuf,0,0,(double)(t)-opt->tshift);
                        corr_prop = mat_get(prop_pt,t,0)-\
                                    fit_model_eval(fm_pt,0,mbuf,par,fmpar_pt);
                        mat_set(prop_pt,t,0,corr_prop);
                    }
                }
                rs_sample_effmass(s_effmass,tem,s_mprop,1,emtype);
                rs_sample_varp(sigem,s_effmass);
                mat_eqsqrt(sigem);
                plot_add_dat(p,tem,em,NULL,sigem,"","rgb 'blue'");
                plot_add_dat(p,tem,em,NULL,NULL,"","rgb 'blue'");
            }
            if (opt->do_plot)
            {
                plot_disp(p);
            }
            if (opt->do_save_plot)
            {
                sprintf(dirname,"%s_em",opt->save_plot_dir);
                plot_save(dirname,p);
            }
            plot_destroy(p);
        }
        else
        {
            /* chi^2 plot */
            p = plot_create();
            plot_set_scale_manual(p,0,dmaxt,0,5.0);
            plot_add_hline(p,1.0,"rgb 'black'");
            plot_add_dat(p,scanres_t,scanres_chi2,NULL,NULL,"chi^2/dof",\
                         "rgb 'blue'");
            plot_set_xlabel(p,xlabel);
            if (opt->do_plot)
            {
                plot_disp(p);
            }
            if (opt->do_save_plot)
            {
                sprintf(dirname,"%s_chi2",opt->save_plot_dir);
                plot_save(dirname,p);
            }
            plot_destroy(p);
            
            /* mass plot */
            p = plot_create();
            plot_set_scale_xmanual(p,0,dmaxt);
            sprintf(key,"a*M_%s",full_name);
            plot_add_dat(p,scanres_t,scanres_mass,NULL,scanres_masserr,key,\
                         "rgb 'red'");
            plot_set_xlabel(p,xlabel);
            if (opt->do_plot)
            {
                plot_disp(p);
            }
            if (opt->do_save_plot)
            {
                sprintf(dirname,"%s_mass",opt->save_plot_dir);
                plot_save(dirname,p);
            }
            plot_destroy(p);
        }

        mat_destroy(em_i);
        mat_destroy(sigem_i);
        mat_destroy(mbuf);
        mat_destroy(par);
        mat_destroy(ft);
        mat_destroy(comp);
    }
    
    /*              desallocation               */
    /********************************************/
    free(opt);
    io_finish();
    mat_ar_destroy(prop[0],nbdat);
    mat_ar_destroy(prop[1],nbdat);
    rs_sample_destroy(s_mprop);
    rs_sample_destroy(s_tmp);
    mat_destroy(sigmprop);
    rs_sample_destroy(s_effmass);
    mat_destroy(tem);
    mat_destroy(sigem);
    fit_data_destroy(d);
    rs_sample_destroy(s_mass);
    mat_destroy(limit);
    mat_destroy(sigmass);
    mat_destroy(scanres_t);
    mat_destroy(scanres_chi2);
    mat_destroy(scanres_mass);
    mat_destroy(scanres_masserr);
    
return EXIT_SUCCESS;
}
