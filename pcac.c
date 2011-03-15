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
    double latspac_nu;
    quark_no q1,q2;
    channel_no ch;
    qcd_options *opt;
    strbuf manf_name,unit;
    size_t binsize;
    
    opt = qcd_arg_parse(argc,argv,A_QCOMP|A_LATSPAC|A_CHANNEL|A_PROP_LOAD\
                        |A_SAVE_RS|A_PLOT|A_LOAD_RG|A_FIT);
    strbufcpy(manf_name,opt->manf_name);
    source = opt->source;
    sink = opt->sink;
    binsize = opt->binsize;
    q1 = opt->qcomp[0];
    q2 = opt->qcomp[1];
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
    for (ch=0;ch<NCHANNEL;ch++)
    {
        channel_id_set(ch,opt->channel_id[ch]);
    }
    minimizer_set_alg(opt->minimizer);
    io_set_fmt(opt->latan_fmt);
    io_init();
    
    /*          creating PCAC "particle"        */
    /********************************************/
    hadron *h_AP,*h_PP;
    
    h_AP = hadron_create();
    h_PP = hadron_create();
    
    hadron_set_2q_nomix(h_AP,"AP",EVEN,ch_AP,q1,q2);
    hadron_set_2q_nomix(h_PP,"PP",ODD,ch_PP,q1,q2);
    
    /*              loading datas               */
    /********************************************/
    size_t ndat,nbdat,nt;
    mat **prop_AP;
    mat **prop_PP;
    
    ndat  = (size_t)get_nfile(manf_name);
    nbdat = ndat/binsize;
    hadron_prop_load_nt(&nt,h_AP,source,sink,manf_name);
    
    prop_AP = mat_ar_create(nbdat,nt,1);
    prop_PP = mat_ar_create(nbdat,nt,1);

    qcd_printf(opt,"-- loading %s datas from %s...\n",h_AP->name,manf_name);
    hadron_prop_load_bin(prop_AP,h_AP,source,sink,manf_name,binsize);
    qcd_printf(opt,"-- loading %s datas from %s...\n",h_PP->name,manf_name);
    hadron_prop_load_bin(prop_PP,h_PP,source,sink,manf_name,binsize);
    
    /*      resampling mean propagator          */
    /********************************************/
    rs_sample *s_mprop_AP, *s_mprop_PP;
    strbuf sample_name;
    
    s_mprop_AP = rs_sample_create(nt,NBOOT);
    s_mprop_PP = rs_sample_create(nt,NBOOT);
    
    sprintf(sample_name,"%s_prop_AP_%s",opt->qcomp_str,manf_name);
    rs_sample_set_name(s_mprop_AP,sample_name);
    qcd_printf(opt,"-- resampling AP %s mean propagator...\n",opt->qcomp_str);
    randgen_set_state(opt->state);
    resample(s_mprop_AP,prop_AP,nbdat,&rs_mean,BOOT,NULL);
    sprintf(sample_name,"%s_prop_PP_%s",opt->qcomp_str,manf_name);
    rs_sample_set_name(s_mprop_AP,sample_name);
    qcd_printf(opt,"-- resampling PP %s mean propagator...\n",opt->qcomp_str);
    randgen_set_state(opt->state);
    resample(s_mprop_PP,prop_PP,nbdat,&rs_mean,BOOT,NULL);
    
    /*      resampling PCAC effective mass      */
    /********************************************/
    rs_sample *s_effmass_pcac;
    mat *effmass_pcac;
    
    s_effmass_pcac = rs_sample_create(nt-2,NBOOT);
    
    sprintf(sample_name,"%s_effmass_PCAC_%s",opt->qcomp_str,manf_name);
    rs_sample_set_name(s_effmass_pcac,sample_name);
    qcd_printf(opt,"-- resampling %s PCAC effective mass...\n",opt->qcomp_str);
    rs_sample_effmass_PCAC(s_effmass_pcac,s_mprop_AP,s_mprop_PP);
    effmass_pcac = rs_sample_pt_cent_val(s_effmass_pcac);
    
    /* computing variance on PCAC effective mass*/
    /********************************************/
    mat *sigem;
    
    sigem = mat_create(nt-2,1);
    
    qcd_printf(opt,"-- estimating %s PCAC effective mass variance...\n",\
               opt->qcomp);
    rs_sample_varp(sigem,s_effmass_pcac);
    
    /*              fit mass                    */
    /********************************************/
    
    fit_data *d;
    size_t t,i;
    rs_sample *s_mass_pcac;
    mat *mass_pcac;
    
    d = fit_data_create(nt-2,1);
    s_mass_pcac = rs_sample_create(1,NBOOT);
    
    sprintf(sample_name,"%s_massfit_PCAC_%s.boot",opt->qcomp_str,manf_name);
    rs_sample_set_name(s_mass_pcac,sample_name);
    qcd_printf(opt,"-- fitting and resampling %s PCAC mass...\n",opt->qcomp_str);
    fit_data_fit_all_points(d,false);
    for (t=0;t<nt-2;t++)
    {
        fit_data_set_x(d,t,0,(double)(t+1));
        if ((t>=nt/4-1)&&(t<=3*nt/4-1))
        {
            fit_data_fit_point(d,t,true);
        }
    }
    fit_data_set_model(d,&fm_const,NULL);
    mat_set(rs_sample_pt_cent_val(s_mass_pcac),0,0,\
            mat_get(effmass_pcac,nt/4-1,0));
    for (i=0;i<rs_sample_get_nsample(s_mass_pcac);i++)
    {
        mat_cp(rs_sample_pt_sample(s_mass_pcac,i),\
               rs_sample_pt_cent_val(s_mass_pcac));
    }
    rs_data_fit(s_mass_pcac,s_effmass_pcac,d,NO_COR);
    mass_pcac = rs_sample_pt_cent_val(s_mass_pcac);
    if (opt->do_save_rs_sample)
    {
        rs_sample_save(s_mass_pcac->name,'w',s_mass_pcac);
    }
    
    /*      computing variance on PCAC mass     */
    /********************************************/
    mat *sigmass;
    
    sigmass = mat_create(1,1);
    
    qcd_printf(opt,"-- estimating %s PCAC mass variance...\n",opt->qcomp_str);
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
    if (opt->qcomp_str[0] == opt->qcomp_str[1])
    {
        qcd_printf(opt,"m_%c\t\t",opt->qcomp_str[0]);
    }
    else
    {
        qcd_printf(opt,"m_%s\t\t",opt->qcomp_str);
    }
    qcd_printf(opt,"= %.8f +/- %.8e %s\n",mat_get(mass_pcac,0,0),\
               mat_get(sigmass,0,0),unit);
    qcd_printf(opt,"dof\t\t= %d\n",fit_data_get_dof(d));
    qcd_printf(opt,"chi^2/dof\t= %e\n",fit_data_get_chi2pdof(d));
    
    /*                  plot                    */
    /********************************************/
    if (opt->do_plot)
    {
        plot *p;
        mat *pr_t;
        strbuf key,ylabel;
        
        p    = plot_create();
        pr_t = mat_create(nt-2,1);
        
        sprintf(key,"%s PCAC effective mass",opt->qcomp_str);
        sprintf(ylabel,"mass%s",unit);
        plot_set_xlabel(p,"time");
        plot_set_ylabel(p,ylabel);
        plot_set_scale_xmanual(p,2.0,nt-1);
        plot_add_hlineerr(p,mat_get(mass_pcac,0,0),mat_get(sigmass,0,0),"0",\
                          "rgb 'red'","rgb 'light-red'");
        mat_set_step(pr_t,1.0,1.0);
        plot_add_dat_yerr(p,pr_t,effmass_pcac,sigem,key,"rgb 'blue'");
        plot_disp(p);   
        
        plot_destroy(p);
        mat_destroy(pr_t);
    }
    
    /*              desallocation               */
    /********************************************/
    FREE(opt);
    hadron_destroy(h_AP);
    hadron_destroy(h_PP);
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
