#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <qcd_arg_parse.h>
#include <latan/latan_hadron.h>
#include <latan/latan_io.h>
#include <latan/latan_mat.h>
#include <latan/latan_math.h>
#include <latan/latan_mass.h>
#include <latan/latan_minimizer.h>
#include <latan/latan_models.h>
#include <latan/latan_plot.h>
#include <latan/latan_rand.h>
#include <latan/latan_statistics.h>

#define NBOOT 2000

int main(int argc, char* argv[])
{
    /*              parsing arguments           */
    /********************************************/
    ss_no source,sink;
    channel_no ch;
    double latspac_nu;
    qcd_options *opt;
    strbuf spec_name,part_name,manf_name,unit;
    size_t binsize;
    
    opt = qcd_arg_parse(argc,argv,A_PARTICLE|A_PROP_LOAD|A_LATSPAC|A_FIT\
                        |A_SAVE_RS|A_PLOT|A_CHANNEL|A_LOAD_RG);
    source = opt->source;
    sink = opt->sink;
    binsize = opt->binsize;
    strbufcpy(spec_name,opt->spec_name);
    strbufcpy(part_name,opt->part_name);
    strbufcpy(manf_name,opt->manf_name);
    latspac_nu = opt->latspac_nu;
    if (opt->have_latspac)
    {
        strbufcpy(unit," (MeV)");
    }
    else
    {
        strbufcpy(unit,"");
    }
    for (ch=0;ch<NCHANNEL;ch++)
    {
        channel_id_set(ch,opt->channel_id[ch]);
    }
    latan_set_verb(opt->latan_verb);
    minimizer_set_alg(opt->minimizer);
    
    /*          identifying particle            */
    /********************************************/
    spectrum *s;
    hadron *h;
    
    if (strcmp(spec_name,"qcd") == 0)
    {
        s = spectrum_create_qcd();
    }
    else if (strcmp(spec_name,"qcdqed") == 0)
    {
        s = spectrum_create_qcdqed();
    }
    else
    {
        s = NULL;
        fprintf(stderr,"error: spectrum %s unknown\n",spec_name);
        return EXIT_FAILURE;
    }
    h = spectrum_get(s,part_name);
    
    /*              loading datas               */
    /********************************************/
    size_t ndat,nbdat,nt;
    mat **prop;
    
    ndat    = (size_t)get_nfile(manf_name);
    nbdat   = ndat/binsize;
    hadron_prop_load_nt(&nt,h,source,sink,manf_name);
    
    prop = mat_ar_create(nbdat,nt,1);

    io_init();
    qcd_printf(opt,"-- loading %s datas from %s...\n",h->name,manf_name);
    hadron_prop_load_bin(prop,h,source,sink,manf_name,binsize);

    /*      resampling mean propagator          */
    /********************************************/
    rs_sample *s_mprop;
    strbuf sample_name;
    mat *mprop;

    s_mprop = rs_sample_create(nt,NBOOT);
    
    sprintf(sample_name,"%s_prop_%s",h->name,manf_name);
    rs_sample_set_name(s_mprop,sample_name);
    qcd_printf(opt,"-- resampling %s mean propagator...\n",h->name);
    randgen_set_state(opt->state);
    resample(s_mprop,prop,nbdat,1,&rs_mean,BOOT,NULL);
    mprop = rs_sample_pt_cent_val(s_mprop);
    
    /*  computing error on mean propagator      */
    /********************************************/
    mat *sigmprop;
    
    sigmprop = mat_create(nt,1);
    
    qcd_printf(opt,"-- estimating %s mean propagator variance...\n",h->name);
    rs_sample_varp(sigmprop,s_mprop);
    
    /*      resampling effective mass           */
    /********************************************/
    rs_sample *s_effmass;
    mat *em;

    s_effmass = rs_sample_create(nt-2,NBOOT);
    
    sprintf(sample_name,"%s_effmass_%s",h->name,manf_name);
    rs_sample_set_name(s_effmass,sample_name);
    qcd_printf(opt,"-- resampling %s effective mass...\n",h->name);
    randgen_set_state(opt->state);
    resample(s_effmass,prop,nbdat,1,&rs_effmass,BOOT,&(h->parity));
    em = rs_sample_pt_cent_val(s_effmass);
    
    /*  computing error on effective mass       */
    /********************************************/
    mat *sigem;
    
    sigem = mat_create(nt-2,1);
    
    qcd_printf(opt,"-- computing %s effective mass variance...\n",h->name);
    rs_sample_varp(sigem,s_effmass);
    mat_eqsqrt(sigem);
    
    /*              fitting mass                */
    /********************************************/
    fit_data *d;
    rs_sample *s_mass;
    mat *mass;
    
    d = fit_data_create(nt,1);
    s_mass = rs_sample_create(2,NBOOT);
    
    sprintf(sample_name,"%s_mass_fit_%s",h->name,manf_name);
    rs_sample_set_name(s_mass,sample_name);
    qcd_printf(opt,"-- fitting and resampling %s mass...\n",h->name);
    fit_data_mass_fit_tune(d,rs_sample_pt_cent_val(s_mass),mprop,em,sigem,\
                           h->parity);
    rs_data_fit(s_mass,s_mprop,d);
    mass = rs_sample_pt_cent_val(s_mass);
    if (opt->do_save_rs_sample)
    {
        rs_sample_save(s_mass->name,'w',s_mass);
    }
    
    /*      computing error on mass             */
    /********************************************/
    mat *sigmass;
    
    sigmass = mat_create(2,1);
    
    qcd_printf(opt,"-- estimating %s mass variance...\n",h->name);
    rs_sample_varp(sigmass,s_mass);
    mat_eqsqrt(sigmass);
    mat_eqsqrt(sigmprop);
    
    /*      switching to right units            */
    /********************************************/
    mat_eqmuls(em,1.0/latspac_nu);
    mat_eqmuls(sigem,1.0/latspac_nu);
    
    /*              result output               */
    /********************************************/
    qcd_printf(opt,"M_%s\t\t= %.8f +/- %.8e %s\n",h->name,                  \
               mat_get(mass,0,0)/latspac_nu,mat_get(sigmass,0,0)/latspac_nu,\
               unit);
    qcd_printf(opt,"dof\t\t= %d\n",fit_data_get_dof(d));
    qcd_printf(opt,"chi^2/dof\t= %e\n",fit_data_get_chi2pdof(d));
    
    /*                  plot                    */
    /********************************************/
    if (opt->do_plot)
    {
        mat *em_t;
        plot *p;
        strbuf key,plotcmd,plotfmt;
        size_t maxt;
        double dmaxt,shift;
        
        maxt  = (h->parity == EVEN) ? (nt/2) : (nt-1);
        dmaxt = (double)maxt;
        shift = (h->parity == EVEN) ? (0.0) : (-DRATIO(nt,2));
        
        em_t = mat_create(nt-2,1);
        
        /* propagator plot */
        p = plot_create();
        
        plot_set_scale_ylog(p);
        plot_set_scale_xmanual(p,shift,dmaxt+shift);
        sprintf(key,"%s propagator",h->name);
        plot_add_dat_yerr(p,fit_data_pt_x(d),mprop,sigmprop,key,"rgb 'red'");
        fit_model_get_plot_fmt(plotfmt,fit_data_pt_model(d));
        sprintf(plotcmd,plotfmt,mat_get(mass,0,0),mat_get(mass,1,0));
        strcat(plotcmd," t 'fit' lc rgb 'red'");
        plot_add_plot(p,plotcmd);
        plot_disp(p);   
        
        plot_destroy(p);
        
        /* effective mass plot */
        p = plot_create();
        
        plot_set_scale_xmanual(p,0.0,dmaxt);
        plot_add_hlineerr(p,mat_get(mass,0,0)/latspac_nu,\
                          mat_get(sigmass,0,0)/latspac_nu,"0",\
                          "rgb 'red'","rgb 'light-red'");
        sprintf(key,"%s effective mass",h->name);
        mat_set_step(em_t,1.0,1.0);
        plot_add_dat_yerr(p,em_t,em,sigem,key,"rgb 'blue'");
        plot_disp(p);
        
        plot_destroy(p);
    }
    
    /*              desallocation               */
    /********************************************/
    FREE(opt);
    spectrum_destroy(s); 
    /* io_finish(); */ /* unknown memory leak here */
    mat_ar_destroy(prop,nbdat);
    rs_sample_destroy(s_mprop);
    mat_destroy(sigmprop);
    rs_sample_destroy(s_effmass);
    mat_destroy(sigem);
    fit_data_destroy(d);
    rs_sample_destroy(s_mass);
    mat_destroy(sigmass);
    
    return EXIT_SUCCESS;
}
