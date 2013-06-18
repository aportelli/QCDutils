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

#define MAX_REL_ERR 0.4

int main(int argc, char* argv[])
{
    /*              parsing arguments           */
    /********************************************/
    double latspac_nu;
    qcd_options *opt;
    strbuf manf_name,unit,AP_name,PP_name,quark;
    size_t binsize,nboot;
    char *cpt;
    
    opt = qcd_arg_parse(argc,argv,A_PLOT|A_SAVE_RS|A_LOAD_RG|A_PROP_NAME\
                        |A_PROP_LOAD|A_LATSPAC|A_QCOMP|A_FIT,1,0);
    cpt = strchr(opt->channel[0],'/');
    if (cpt == NULL)
    {
        fprintf(stderr,"error: channel option %s is invalid\n",\
                opt->channel[0]);
        exit(EXIT_FAILURE);
    }
    sprintf(PP_name,"%s_%s_%s_%s",cpt+1,opt->quark[0],opt->sink,opt->source);
    *cpt = '\0';
    sprintf(AP_name,"%s_%s_%s_%s",opt->channel[0],opt->quark[0],opt->sink,\
            opt->source);
    strbufcpy(quark,opt->quark[0]);
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
    io_set_fmt(opt->latan_fmt);
    io_init();
    
    /*              loading datas               */
    /********************************************/
    size_t ndat,nbdat,dim[2];
    mat **prop_AP, **prop_PP;
    
    ndat  = (size_t)get_nfile(manf_name);
    nbdat = ndat/binsize + ((ndat%binsize == 0) ? 0 : 1);
    
    mat_ar_loadbin(NULL,dim,manf_name,PP_name,1);
    
    prop_AP = mat_ar_create(nbdat,dim[0],dim[1]);
    prop_PP = mat_ar_create(nbdat,dim[0],dim[1]);

    qcd_printf(opt,"-- loading %s datas from %s...\n",AP_name,manf_name);
    mat_ar_loadbin(prop_AP,NULL,manf_name,AP_name,binsize);
    qcd_printf(opt,"-- loading %s datas from %s...\n",PP_name,manf_name);
    mat_ar_loadbin(prop_PP,NULL,manf_name,PP_name,binsize);
    
    /*                propagator                */
    /********************************************/
    rs_sample *s_mprop_AP, *s_mprop_PP;
    
    s_mprop_AP = rs_sample_create(dim[0],dim[1],nboot);
    s_mprop_PP = rs_sample_create(dim[0],dim[1],nboot);
    
    qcd_printf(opt,"-- resampling %s mean propagator...\n",AP_name);
    randgen_set_state(opt->state);
    resample(s_mprop_AP,prop_AP,nbdat,&rs_mean,BOOT,NULL);
    qcd_printf(opt,"-- resampling %s mean propagator...\n",PP_name);
    randgen_set_state(opt->state);
    resample(s_mprop_PP,prop_PP,nbdat,&rs_mean,BOOT,NULL);
    
    /*            PCAC effective mass           */
    /********************************************/
    rs_sample *s_effmass_pcac;
    mat *effmass_pcac,*sigem;
    size_t nt;
    
    nt = dim[0];
    
    s_effmass_pcac = rs_sample_create(nt-2,1,nboot);
    sigem          = mat_create(nt-2,1);
    
    qcd_printf(opt,"-- resampling %s PCAC effective mass...\n",quark);
    rs_sample_effmass_PCAC(s_effmass_pcac,s_mprop_AP,s_mprop_PP);
    effmass_pcac = rs_sample_pt_cent_val(s_effmass_pcac);
    rs_sample_varp(sigem,s_effmass_pcac);
    
    /*                fit mass                  */
    /********************************************/
    fit_data *d;
    size_t t,i;
    rs_sample *s_mass_pcac;
    mat *mass_pcac,*sigmass;
    strbuf latan_path;
    
    d           = fit_data_create(nt-2,1,1);
    s_mass_pcac = rs_sample_create(1,1,nboot);
    sigmass     = mat_create(1,1);
    
    qcd_printf(opt,"-- fitting and resampling %s PCAC mass...\n",quark);
    fit_data_fit_all_points(d,false);
    for (t=0;t<nt-2;t++)
    {
        fit_data_set_x(d,t,0,(double)(t+1));
        if ((t>=nt/4-1)&&(t<=3*nt/4-1))
        {
            fit_data_fit_point(d,t,true);
        }
    }
    if (opt->nmanrange > 0)
    {
        fit_data_fit_all_points(d,false);
        fit_data_fit_range(d,opt->range[0][0]-1,opt->range[0][1]-1,true);
    }
    fit_data_set_model(d,&fm_const,NULL);
    mat_set(rs_sample_pt_cent_val(s_mass_pcac),0,0,\
            mat_get(effmass_pcac,nt/4-1,0));
    for (i=0;i<rs_sample_get_nsample(s_mass_pcac);i++)
    {
        mat_cp(rs_sample_pt_sample(s_mass_pcac,i),\
               rs_sample_pt_cent_val(s_mass_pcac));
    }
    rs_data_fit(s_mass_pcac,NULL,NULL,&s_effmass_pcac,d,opt->corr,NULL);
    mass_pcac = rs_sample_pt_cent_val(s_mass_pcac);
    if (opt->do_save_rs_sample)
    {
        sprintf(latan_path,"%s_PCAC_mass_fit_%s.boot:%s_PCAC_mass_fit_%s",\
                quark,manf_name,quark,manf_name);
        rs_sample_save(latan_path,'w',s_mass_pcac);
    }
    rs_sample_varp(sigmass,s_mass_pcac);
    mat_eqsqrt(sigmass);
    mat_eqsqrt(sigem);

    /*      switching to right units            */
    /********************************************/
    mat_eqmuls(effmass_pcac,1.0/latspac_nu);
    mat_eqmuls(sigem,1.0/latspac_nu);
    mat_eqmuls(mass_pcac,1.0/latspac_nu);
    mat_eqmuls(sigmass,1.0/latspac_nu);
    
    /*              result output               */
    /********************************************/
    strbuf buf;
    
    sprintf(buf,"m_%s_PCAC",quark);
    qcd_printf(opt,"%-10s= %.8f +/- %.8e %s\n",buf,mat_get(mass_pcac,0,0),\
               mat_get(sigmass,0,0),unit);
    qcd_printf(opt,"%-10s= %d\n","dof",fit_data_get_dof(d));
    qcd_printf(opt,"%-10s= %e\n","chi^2/dof",fit_data_get_chi2pdof(d));
    
    /*                  plot                    */
    /********************************************/
    if (opt->do_plot)
    {
        plot *p;
        mat *pr_t;
        strbuf key,ylabel;
        
        p    = plot_create();
        pr_t = mat_create(nt-2,1);
        
        sprintf(key,"%s PCAC effective mass",quark);
        sprintf(ylabel,"mass%s",unit);
        plot_set_xlabel(p,"time");
        plot_set_ylabel(p,ylabel);
        plot_set_scale_xmanual(p,2.0,nt-1);
        plot_add_hlineerr(p,mat_get(mass_pcac,0,0),mat_get(sigmass,0,0),\
                          "rgb 'red'");
        mat_set_step(pr_t,1.0,1.0);
        plot_add_dat(p,pr_t,effmass_pcac,NULL,sigem,key,"rgb 'blue'");
        plot_disp(p);   
        if (opt->do_save_plot)
        {
            plot_save(opt->save_plot_dir,p);
        }
        plot_destroy(p);
        mat_destroy(pr_t);
    }
    
    /*              desallocation               */
    /********************************************/
    FREE(opt);
    io_finish();
    mat_ar_destroy(prop_AP,nbdat);
    mat_ar_destroy(prop_PP,nbdat);
    rs_sample_destroy(s_mprop_AP);
    rs_sample_destroy(s_mprop_PP);
    rs_sample_destroy(s_effmass_pcac);
    mat_destroy(sigem);
    fit_data_destroy(d);
    rs_sample_destroy(s_mass_pcac);
    mat_destroy(sigmass);
    
    return EXIT_SUCCESS;
}
