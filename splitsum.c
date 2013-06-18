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
    int emtype;
    strbuf prop_name[2],full_name[2],manf_name,model,unit;
    size_t binsize,nboot,propdim[2],nt,nstate;
    fit_model *fm_pt;
    void *fmpar_pt;
    
    opt = qcd_arg_parse(argc,argv,A_PLOT|A_SAVE_RS|A_LOAD_RG|A_PROP_NAME\
                        |A_PROP_LOAD|A_LATSPAC|A_QCOMP|A_FIT|A_MODEL,2,0);
    

    sprintf(full_name[0],"%s_%s",opt->channel[0],opt->quark[0]);
    sprintf(full_name[1],"%s_%s",opt->channel[1],opt->quark[1]);
    sprintf(prop_name[0],"%s_%s_%s_%s",opt->channel[0],opt->quark[0],opt->sink,\
            opt->source);
    sprintf(prop_name[1],"%s_%s_%s_%s",opt->channel[1],opt->quark[1],\
            opt->sink,opt->source);
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
    mat_ar_loadbin(NULL,propdim,manf_name,prop_name[0],1);
    nt = (opt->nt) ? opt->nt : propdim[0];
    if (strbufcmp(model,"cosh") == 0)
    {
        fm_pt    = &fm_cosh_splitsum;
        fmpar_pt = &nt;
        emtype   = EM_ACOSH;
        nstate   = 1;
    }
    else if (strbufcmp(model,"cosh_ex") == 0)
    {
        fm_pt    = &fm_cosh_ex_splitsum;
        fmpar_pt = &nt;
        emtype   = EM_ACOSH;
        nstate   = 2;
    }
    else if (strbufcmp(model,"expdec_ex") == 0)
    {
        fm_pt    = &fm_expdec_ex_splitsum;
        fmpar_pt = NULL;
        emtype   = EM_LOG;
        nstate   = 2;
    }
    else
    {
        strbufcpy(opt->model,"expdec");
        fm_pt    = &fm_expdec_splitsum;
        fmpar_pt = NULL;
        emtype   = EM_LOG;
        nstate   = 1;
    }
    io_set_fmt(opt->latan_fmt);
    io_init();
    
    /*              loading datas               */
    /********************************************/
    size_t ndat,nbdat;
    mat **prop[2];
    
    ndat    = (size_t)get_nfile(manf_name);
    nbdat   = ndat/binsize + ((ndat%binsize == 0) ? 0 : 1);
    
    prop[0] = mat_ar_create(nbdat,propdim[0],propdim[1]);
    prop[1] = mat_ar_create(nbdat,propdim[0],propdim[1]);
    
    qcd_printf(opt,"-- loading %s datas from %s...\n",prop_name[0],manf_name);
    mat_ar_loadbin(prop[0],NULL,manf_name,prop_name[0],binsize);
    qcd_printf(opt,"-- loading %s datas from %s...\n",prop_name[1],manf_name);
    mat_ar_loadbin(prop[1],NULL,manf_name,prop_name[1],binsize);
    
    /*                propagator                */
    /********************************************/
    rs_sample *s_mprop[2];
    mat *mprop[2],*sigmprop[2];
    
    s_mprop[0]  = rs_sample_create(propdim[0],propdim[1],nboot);
    s_mprop[1]  = rs_sample_create(propdim[0],propdim[1],nboot);
    sigmprop[0] = mat_create(propdim[0],propdim[1]);
    sigmprop[1] = mat_create(propdim[0],propdim[1]);
    
    qcd_printf(opt,"-- resampling %s mean propagator...\n",prop_name[0]);
    randgen_set_state(opt->state);
    resample(s_mprop[0],prop[0],nbdat,&rs_mean,BOOT,NULL);
    qcd_printf(opt,"-- resampling %s mean propagator...\n",prop_name[1]);
    randgen_set_state(opt->state);
    resample(s_mprop[1],prop[1],nbdat,&rs_mean,BOOT,NULL);  
    mprop[0] = rs_sample_pt_cent_val(s_mprop[0]);
    mprop[1] = rs_sample_pt_cent_val(s_mprop[1]);
    rs_sample_varp(sigmprop[0],s_mprop[0]);
    mat_eqsqrt(sigmprop[0]);
    rs_sample_varp(sigmprop[1],s_mprop[1]);
    mat_eqsqrt(sigmprop[1]);
    
    /*           effective mass                 */
    /********************************************/
    rs_sample *s_effmass[2];
    mat *tem,*em[2],*sigem[2];
    size_t emdim[2];
    
    get_effmass_size(emdim,mprop[0],1,emtype);
    
    s_effmass[0] = rs_sample_create(emdim[0],emdim[1],nboot);
    s_effmass[1] = rs_sample_create(emdim[0],emdim[1],nboot);
    tem          = mat_create(emdim[0],1);
    sigem[0]     = mat_create(emdim[0],emdim[1]);
    sigem[1]     = mat_create(emdim[0],emdim[1]);
    
    qcd_printf(opt,"-- resampling %s effective mass...\n",prop_name[0]);
    rs_sample_effmass(s_effmass[0],tem,s_mprop[0],1,emtype);
    qcd_printf(opt,"-- resampling %s effective mass...\n",prop_name[1]);
    rs_sample_effmass(s_effmass[1],tem,s_mprop[1],1,emtype);
    em[0] = rs_sample_pt_cent_val(s_effmass[0]);
    em[1] = rs_sample_pt_cent_val(s_effmass[1]);
    rs_sample_varp(sigem[0],s_effmass[0]);
    mat_eqsqrt(sigem[0]);
    rs_sample_varp(sigem[1],s_effmass[1]);
    mat_eqsqrt(sigem[1]);
    
    /*                  fit mass                */
    /********************************************/
    fit_data *d;
    rs_sample *s_par,*s_mass[2],*s_av_mass,*s_d_mass,*s_d_sqmass;
    mat *par,*limit,*sigpar,*scanres_t,*scanres_chi2,*scanres_av_mass,\
        *scanres_av_mass_err,*scanres_d_mass,*scanres_d_mass_err;
    size_t npar,nti,tibeg,ta,range[2];
    size_t i,j;
    strbuf buf,range_info,latan_path;
    double pref_i,m_i,dm_i;
    
    d        = fit_data_create(propdim[0],1,2);
    tibeg    = (size_t)(opt->range[0][0]); 
    range[0] = 0;
    range[1] = 0;
    npar     = fit_model_get_npar(fm_pt,fmpar_pt);
    
    s_par      = rs_sample_create(npar,1,nboot);
    s_mass[0]  = rs_sample_create(1,1,nboot);
    s_mass[1]  = rs_sample_create(1,1,nboot);
    s_av_mass  = rs_sample_create(1,1,nboot);
    s_d_mass   = rs_sample_create(1,1,nboot);
    s_d_sqmass = rs_sample_create(1,1,nboot);
    par        = rs_sample_pt_cent_val(s_par);
    limit      = mat_create(npar,2);
    sigpar     = mat_create(npar,1);
    
    /** print operation **/
    if (!opt->do_range_scan)
    {
        qcd_printf(opt,"-- fitting and resampling %s,%s mass sum/difference...\n",\
                   full_name[0],full_name[1]);
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
        correct_range(range,mprop[0],sigmprop[0],opt);
        fit_data_fit_range(d,range[0],range[1],true);
    }
    qcd_printf(opt,"\n");
    
    nti                 = MAX(range[1] - 1 - tibeg,1);
    scanres_t           = mat_create(nti,1);
    scanres_chi2        = mat_create(nti,1);
    scanres_av_mass     = mat_create(nti,1);
    scanres_av_mass_err = mat_create(nti,1);
    scanres_d_mass      = mat_create(nti,1);
    scanres_d_mass_err  = mat_create(nti,1);
    
    /** set model **/
    fit_data_set_model(d,fm_pt,fmpar_pt);
    
    /** set correlation filter **/
    if (opt->corr == NO_COR)
    {
        for (i=0;i<propdim[0];i++) 
        for (j=0;j<propdim[0];j++)
        {
            if (i != j)
            {
                fit_data_set_data_cor(d,i,j,false);
            }
        }
    }
    
    /** set initial parameter values **/
    ta  = nt/8-(size_t)(mat_get(tem,0,0));
    m_i = mat_get(em[0],ta,0);
    if (latan_isnan(m_i))
    {
        m_i = 0.1;
    }
    dm_i = 0.01*m_i;
    mat_set(par,0,0,m_i);
    mat_set(par,1,0,dm_i);
    for (i=2;i<2*nstate;i+=2)
    {
        mat_set(par,i,0,((double)(i/2+2))*m_i);
        mat_set(par,i+1,0,((double)(i/2+2))*m_i);
    }
    pref_i = log(mat_get(mprop[0],ta,0)/(exp(-m_i*ta)+exp(-m_i*(nt-ta))));
    if (latan_isnan(pref_i))
    {
        pref_i = 1.0;
    }
    pref_i -= log((double)(nstate));
    for (i=2*nstate;i<npar;i+=2)
    {
        mat_set(par,i,0,pref_i);
        mat_set(par,i+1,0,pref_i);
    }
    qcd_printf(opt,"%-22sav_mass= %e -- d_mass= %e -- prefactor0,1_0= %e\n",\
               "initial parameters: ",mat_get(par,0,0),mat_get(par,1,0),\
               exp(pref_i));
    for (i=2;i<nstate+1;i++)
    {
        qcd_printf(opt,"%-22sE0,1_%d= %e -- prefactor0,1_%d= %e\n","",\
                   (int)(i/2),mat_get(par,i,0),(int)(i/2),exp(pref_i));
    }
    
    /** set parameter limits **/
    mat_cst(limit,latan_nan());
    mat_set(limit,0,0,0.5*m_i);
    mat_set(limit,1,0,-0.5*m_i);
    mat_set(limit,1,1,0.5*m_i);
    for (i=2;i<2*nstate;i+=2)
    {
        mat_set(limit,i,0,mat_get(par,i-2,0));
        mat_set(limit,i+1,0,mat_get(par,i-2,0));
    }
    
    /** set x **/
    for (i=0;i<propdim[0];i++)
    {
        fit_data_set_x(d,i,0,(double)(i));
    }
    
    /** regular correlator fit... **/
    if (!opt->do_range_scan)
    {
        latan_set_warn(false);
        rs_data_fit(s_par,limit,NULL,s_mprop,d,NO_COR,NULL);
        latan_set_warn(true);
        rs_data_fit(s_par,limit,NULL,s_mprop,d,opt->corr,NULL);
        rs_sample_varp(sigpar,s_par);
        mat_eqsqrt(sigpar);
        if (fit_data_get_chi2pdof(d) > 2.0)
        {
            fprintf(stderr,"warning: bad final fit (chi^2/dof= %.2e)\n",\
                    fit_data_get_chi2pdof(d));
        }
        qcd_printf(opt,"-- results:\n");
        qcd_printf(opt,"%-10s= %.8f +/- %.8e %s\n","av_mass",\
                   mat_get(par,0,0)/latspac_nu,             \
                   mat_get(sigpar,0,0)/latspac_nu,unit);
        qcd_printf(opt,"%-10s= %.8f +/- %.8e %s\n","d_mass", \
                   mat_get(par,1,0)/latspac_nu,             \
                   mat_get(sigpar,1,0)/latspac_nu,unit);
        for (i=2;i<2*nstate;i+=2)
        {
            sprintf(buf,"E%d_0",(int)(i/2));
            qcd_printf(opt,"%-10s= %.8f +/- %.8e %s\n",buf, \
                       mat_get(par,i,0)/latspac_nu,        \
                       mat_get(sigpar,i,0)/latspac_nu,unit);
            sprintf(buf,"E%d_1",(int)(i/2));
            qcd_printf(opt,"%-10s= %.8f +/- %.8e %s\n",buf, \
                       mat_get(par,i+1,0)/latspac_nu,        \
                       mat_get(sigpar,i+1,0)/latspac_nu,unit);
        }
        qcd_printf(opt,"%-10s= %d\n","dof",fit_data_get_dof(d));
        qcd_printf(opt,"%-10s= %e\n","chi^2/dof",fit_data_get_chi2pdof(d));
        rs_sample_get_subsamp(s_av_mass,s_par,0,0,0,0);
        rs_sample_get_subsamp(s_d_mass,s_par,1,0,1,0);
        rs_sample_muls(s_mass[0],s_d_mass,0.5);
        rs_sample_eqadd(s_mass[0],s_av_mass);
        rs_sample_muls(s_mass[1],s_d_mass,-0.5);
        rs_sample_eqadd(s_mass[1],s_av_mass);
        if (opt->do_save_rs_sample)
        {
            sprintf(latan_path,"%s_%s_av_mass_fit%s_%s.boot:%s_%s_av_mass_fit%s_%s",\
                    full_name[0],full_name[1],range_info,manf_name,\
                    full_name[0],full_name[1],range_info,manf_name);
            rs_sample_save(latan_path,'w',s_av_mass);
            sprintf(latan_path,"%s_mass_fit%s_%s.boot:%s_mass_fit%s_%s",      \
                    full_name[0],range_info,manf_name,full_name[0],range_info,\
                    manf_name);
            rs_sample_save(latan_path,'w',s_mass[0]);
            sprintf(latan_path,"%s_mass_fit%s_%s.boot:%s_mass_fit%s_%s",      \
                    full_name[1],range_info,manf_name,full_name[1],range_info,\
                    manf_name);
            rs_sample_save(latan_path,'w',s_mass[1]);
            sprintf(latan_path,"%s_%s_split_fit%s_%s.boot:%s_%s_split_fit%s_%s",\
                    full_name[0],full_name[1],range_info,manf_name,\
                    full_name[0],full_name[1],range_info,manf_name);
            rs_sample_save(latan_path,'w',s_d_mass);
        }
        rs_sample_mulp(s_d_sqmass,s_d_mass,s_av_mass);
        rs_sample_eqmuls(s_d_sqmass,2.0);
        if (opt->do_save_rs_sample)
        {
            sprintf(latan_path,"%s_%s_sqsplit_fit%s_%s.boot:%s_%s_sqsplit_fit%s_%s",\
                    full_name[0],full_name[1],range_info,manf_name,\
                    full_name[0],full_name[1],range_info,manf_name);
            rs_sample_save(latan_path,'w',s_d_sqmass);
        }
    }
    /** ...or fit range scanning **/
    else
    {
        qcd_printf(opt,"\n%-5s %-12s a*av_M_%-5s %-12s a*d_M_%-6s %-12s",\
                   "ti/a","chi^2/dof",full_name[0],"error",full_name[0],"error");
        for (i=tibeg;i<range[1]-1;i++)
        {
            rs_data_fit(s_par,limit,NULL,s_mprop,d,opt->corr,NULL);
            rs_sample_varp(sigpar,s_par);
            mat_eqsqrt(sigpar);
            mat_set(scanres_t,i-tibeg,0,(double)(i));
            mat_set(scanres_chi2,i-tibeg,0,fit_data_get_chi2pdof(d));
            mat_set(scanres_av_mass,i-tibeg,0,mat_get(par,0,0));
            mat_set(scanres_av_mass_err,i-tibeg,0,mat_get(sigpar,0,0));
            mat_set(scanres_d_mass,i-tibeg,0,mat_get(par,1,0));
            mat_set(scanres_d_mass_err,i-tibeg,0,mat_get(sigpar,1,0));
            qcd_printf(opt,"\n% -4d % -.5e % -.5e % -.5e % -.5e % -.5e",  \
                       (int)(i),fit_data_get_chi2pdof(d),mat_get(par,0,0),\
                       mat_get(sigpar,0,0),mat_get(par,1,0),              \
                       mat_get(sigpar,1,0));
            fit_data_fit_point(d,i,false);
        }
        qcd_printf(opt,"\n\n");
    }
    
    /*                  plot                    */
    /********************************************/
    if (opt->do_plot|opt->do_save_plot)
    {
        mat *pr_t,*mbuf,*em_i,*sigem_i,*prop_pt,*mass_pt,*mpar,*ft,*comp;
        plot *p;
        strbuf key,color[2],dirname;
        size_t maxt,k,t,npt;
        double dmaxt,nmass,corr_prop;
        
        npt     = fit_data_fit_point_num(d);
        
        pr_t    = mat_create(nt,1);
        mbuf    = mat_create(1,1);
        em_i    = mat_create(nrow(em[0]),1);
        sigem_i = mat_create(nrow(em[0]),1);
        mpar    = mat_create(2,1);
        ft      = mat_create(npt,1);
        comp    = mat_create(npt,1);
        
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
            strbufcpy(color[0],"rgb 'red'");
            strbufcpy(color[1],"rgb 'blue'");
            
            /** chi^2 plot **/
            for (k=0;k<2;k++)
            {
                p = plot_create();
                i = 0;
                for (t=0;t<propdim[0];t++)
                {
                    if (fit_data_is_fit_point(d,t))
                    {
                        mat_set(ft,i,0,(double)(t));
                        mat_set(comp,i,0,mat_get(d->chi2_comp,i+k*npt,0));
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
                plot_set_ylabel(p,"standard deviations");
                plot_add_points(p,ft,comp,"",color[k],"impulses");
                if (opt->do_plot)
                {
                    plot_disp(p);
                }
                if (opt->do_save_plot)
                {
                    sprintf(dirname,"%s_%d_dev",opt->save_plot_dir,(int)k);
                    plot_save(dirname,p);
                }
                plot_destroy(p);
            }
            
            /** propagator plot **/
            for (k=0;k<2;k++)
            {
                p = plot_create();
                fit_data_fit_all_points(d,true);
                plot_set_scale_ylog(p);
                plot_set_scale_xmanual(p,0,dmaxt);
                sprintf(key,"%s propagator",full_name[k]);
                mat_eqabs(mprop[k]);
                mat_set_step(pr_t,0.0,1.0);
                plot_add_fit(p,d,k,mbuf,0,par,0,dmaxt,1000,false,\
                             PF_FIT|PF_DATA,key,"",color[k],color[k]);
                if (opt->do_plot)
                {
                    plot_disp(p);
                }
                if (opt->do_save_plot)
                {
                    sprintf(dirname,"%s_%d_prop",opt->save_plot_dir,(int)k);
                    plot_save(dirname,p);
                }
                plot_destroy(p);
            }
            
            /** effective mass plot **/
            rs_sample_set_subsamp(s_par,s_mass[0],0,0,0,0);
            rs_sample_set_subsamp(s_par,s_mass[1],1,0,1,0);
            rs_sample_varp(sigpar,s_par);
            mat_eqsqrt(sigpar);
            for (k=0;k<2;k++)
            {
                p = plot_create();
                mat_eqmuls(em[k],1.0/latspac_nu);
                mat_eqmuls(sigem[k],1.0/latspac_nu);
                nmass = mat_get(par,k,0)/latspac_nu;
                plot_add_hlineerr(p,nmass, mat_get(sigpar,k,0)/latspac_nu,\
                                  color[k]);
                for (i=2;i<2*nstate;i+=2)
                {
                    nmass = mat_get(par,i+k,0)/latspac_nu;
                    plot_add_hlineerr(p,nmass,mat_get(sigpar,i+k,0)/latspac_nu,\
                                      color[k]);
                }
                plot_set_scale_manual(p,0.0,dmaxt,0.0,1.5*nmass);
                sprintf(key,"%s effective energies",full_name[k]);
                plot_add_dat(p,tem,em[k],NULL,sigem[k],key,"rgb 'blue'");
                plot_add_dat(p,tem,em[k],NULL,NULL,"","rgb 'blue'");
                for (i=2;i<2*nstate;i+=2)
                {
                    if (emtype == EM_ACOSH)
                    {
                        fm_pt = &fm_cosh;
                    }
                    else if (emtype == EM_LOG)
                    {
                        fm_pt = &fm_expdec;
                    }
                    for (j=0;j<=rs_sample_get_nsample(s_mprop[k]);j++)
                    {
                        prop_pt = (j==0) ? rs_sample_pt_cent_val(s_mprop[k]) \
                                         : rs_sample_pt_sample(s_mprop[k],j-1);
                        mass_pt = (j==0) ? rs_sample_pt_cent_val(s_par) \
                                         : rs_sample_pt_sample(s_par,j-1);
                        mat_set(mpar,0,0,mat_get(mass_pt,i-2+k,0));
                        mat_set(mpar,1,0,mat_get(mass_pt,2*nstate+i-2+k,0));
                        for (t=0;t<propdim[0];t++)
                        {
                            mat_set(mbuf,0,0,(double)(t));
                            corr_prop = mat_get(prop_pt,t,0)               \
                                        - fit_model_eval(fm_pt,0,mbuf,mpar,\
                                                         fmpar_pt);
                            mat_set(prop_pt,t,0,corr_prop);
                        }
                    }
                    rs_sample_effmass(s_effmass[k],tem,s_mprop[k],1,emtype);
                    rs_sample_varp(sigem[k],s_effmass[k]);
                    mat_eqsqrt(sigem[k]);
                    plot_add_dat(p,tem,em[k],NULL,sigem[k],"","rgb 'blue'");
                    plot_add_dat(p,tem,em[k],NULL,NULL,"","rgb 'blue'");
                }
                if (opt->do_plot)
                {
                    plot_disp(p);
                }
                if (opt->do_save_plot)
                {
                    sprintf(dirname,"%s_%d_em",opt->save_plot_dir,(int)k);
                    plot_save(dirname,p);
                }
                plot_destroy(p);
            }
        }
        else
        {
            /* chi^2 plot */
            p = plot_create();
            plot_set_scale_manual(p,0,dmaxt,0,5.0);
            plot_add_hline(p,1.0,"rgb 'black'");
            plot_add_dat(p,scanres_t,scanres_chi2,NULL,NULL,"chi^2/dof",\
                         "rgb 'blue'");
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
            
            /* average mass plot */
            p = plot_create();
            plot_set_scale_xmanual(p,0,dmaxt);
            sprintf(key,"a*av_M_%s",full_name[0]);
            plot_add_dat(p,scanres_t,scanres_av_mass,NULL,scanres_av_mass_err,\
                         key,"rgb 'red'");
            if (opt->do_plot)
            {
                plot_disp(p);
            }
            if (opt->do_save_plot)
            {
                sprintf(dirname,"%s_av_mass",opt->save_plot_dir);
                plot_save(dirname,p);
            }
            plot_destroy(p);
            
            /* difference of masses plot */
            p = plot_create();
            plot_set_scale_xmanual(p,0,dmaxt);
            sprintf(key,"a*d_M_%s",full_name[0]);
            plot_add_dat(p,scanres_t,scanres_d_mass,NULL,scanres_d_mass_err,\
                         key,"rgb 'red'");
            if (opt->do_plot)
            {
                plot_disp(p);
            }
            if (opt->do_save_plot)
            {
                sprintf(dirname,"%s_dmass",opt->save_plot_dir);
                plot_save(dirname,p);
            }
            plot_destroy(p);
        }
        
        mat_destroy(em_i);
        mat_destroy(sigem_i);
        mat_destroy(pr_t);
        mat_destroy(mbuf);
        mat_destroy(mpar);
        mat_destroy(ft);
        mat_destroy(comp);
    }
    
    /*              desallocation               */
    /********************************************/
    free(opt);
    io_finish();
    mat_ar_destroy(prop[0],nbdat);
    mat_ar_destroy(prop[1],nbdat);
    rs_sample_destroy(s_mprop[0]);
    rs_sample_destroy(s_mprop[1]);
    mat_destroy(sigmprop[0]);
    mat_destroy(sigmprop[1]);
    rs_sample_destroy(s_effmass[0]);
    rs_sample_destroy(s_effmass[1]);
    mat_destroy(tem);
    mat_destroy(sigem[0]);
    mat_destroy(sigem[1]);
    fit_data_destroy(d);
    rs_sample_destroy(s_par);
    rs_sample_destroy(s_mass[0]);
    rs_sample_destroy(s_mass[1]);
    rs_sample_destroy(s_av_mass);
    rs_sample_destroy(s_d_mass);
    rs_sample_destroy(s_d_sqmass);
    mat_destroy(sigpar);
    mat_destroy(scanres_t);
    mat_destroy(scanres_chi2);
    mat_destroy(scanres_av_mass);
    mat_destroy(scanres_av_mass_err);
    mat_destroy(scanres_d_mass);
    mat_destroy(scanres_d_mass_err);
    
    return EXIT_SUCCESS;
}
