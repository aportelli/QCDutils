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
#define MAX_REL_ERR 0.4


int main(int argc, char* argv[])
{
    /*              parsing arguments           */
    /********************************************/
    ss_no source,sink;
    channel_no ch;
    double latspac_nu;
    qcd_options *opt;
    bool is_split;
    strbuf spec_name,part_name[2],full_name,manf_name,unit;
    size_t binsize;
    
    opt = qcd_arg_parse(argc,argv,A_PARTICLE|A_PROP_LOAD|A_LATSPAC|A_FIT\
                        |A_SAVE_RS|A_PLOT|A_CHANNEL|A_LOAD_RG);
    is_split   = (strlen(opt->part_name[1]) != 0);
    strbufcpy(spec_name,opt->spec_name);
    strbufcpy(part_name[0],opt->part_name[0]);
    if (is_split)
    {
        strbufcpy(part_name[1],opt->part_name[1]);
        sprintf(full_name,"%s-%s",part_name[0],part_name[1]);
    }
    else
    {
        strbufcpy(full_name,part_name[0]);
    }
    strbufcpy(manf_name,opt->manf_name);
    source     = opt->source;
    sink       = opt->sink;
    binsize    = opt->binsize;
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
    io_set_fmt(opt->latan_fmt);
    io_init();
    
    /*          identifying particle            */
    /********************************************/
    spectrum *s;
    hadron *h[2];

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
    h[0] = spectrum_get(s,part_name[0]);
    if (is_split)
    {
        h[1] = spectrum_get(s,part_name[1]);
    }
    
    /*              loading datas               */
    /********************************************/
    size_t ndat,nbdat,nt;
    mat **prop[2];
    
    ndat    = (size_t)get_nfile(manf_name);
    nbdat   = ndat/binsize + ((ndat%binsize == 0) ? 0 : 1);
    hadron_prop_load_nt(&nt,h[0],source,sink,manf_name);
    
    prop[0] = mat_ar_create(nbdat,nt,1);
    prop[1] = mat_ar_create(nbdat,nt,1);
    
    qcd_printf(opt,"-- loading %s datas from %s...\n",h[0]->name,manf_name);
    hadron_prop_load_bin(prop[0],h[0],source,sink,manf_name,binsize);
    if (is_split)
    {
        qcd_printf(opt,"-- loading %s datas from %s...\n",h[1]->name,manf_name);
        hadron_prop_load_bin(prop[1],h[1],source,sink,manf_name,binsize);
    }

    /*                propagator                */
    /********************************************/
    rs_sample *s_tmp,*s_mprop;
    strbuf sample_name;
    mat *mprop,*sigmprop;

    s_mprop  = rs_sample_create(nt,NBOOT);
    s_tmp    = rs_sample_create(nt,NBOOT);
    sigmprop = mat_create(nt,1);
    
    sprintf(sample_name,"%s_prop_%s",full_name,manf_name);
    rs_sample_set_name(s_mprop,sample_name);
    qcd_printf(opt,"-- resampling %s mean propagator...\n",h[0]->name);
    randgen_set_state(opt->state);
    resample(s_mprop,prop[0],nbdat,&rs_mean,BOOT,NULL);
    if (is_split)
    {
        qcd_printf(opt,"-- resampling %s mean propagator...\n",h[1]->name);
        randgen_set_state(opt->state);
        resample(s_tmp,prop[1],nbdat,&rs_mean,BOOT,NULL);
        rs_sample_eqdivp(s_mprop,s_tmp);
        h[0]->parity = EVEN;
    }
    mprop = rs_sample_pt_cent_val(s_mprop);
    rs_sample_varp(sigmprop,s_mprop);
    mat_eqsqrt(sigmprop);
    
    /*           effective mass                 */
    /********************************************/
    rs_sample *s_effmass;
    mat *em,*sigem;
    
    em = NULL;
    
    s_effmass = rs_sample_create(nt-2,NBOOT);
    sigem     = mat_create(nt-2,1);
    
    sprintf(sample_name,"%s_effmass_%s",full_name,manf_name);
    rs_sample_set_name(s_effmass,sample_name);
    qcd_printf(opt,"-- resampling %s effective mass...\n",full_name);
    rs_sample_effmass(s_effmass, s_mprop,h[0]->parity);
    em = rs_sample_pt_cent_val(s_effmass);
    qcd_printf(opt,"-- computing %s effective mass variance...\n",full_name);
    rs_sample_varp(sigem,s_effmass);
    mat_eqsqrt(sigem);
    
    /*              fitting mass                */
    /********************************************/
    fit_data *d;
    rs_sample *s_mass;
    mat *mass,*sigmass,*scanres_t,*scanres_chi2,*scanres_mass,*scanres_masserr;
    size_t npar,nti,tibeg,inrmin,rmin,inrmax,rmax,best_t;
    size_t i;
    bool first_elim;
    strbuf buf,range_info;
    double rerr,uc_mass,uc_mass_err,dev;

    npar       = 2;
    d          = fit_data_create(nt,1,1);
    rmax       = 0;
    rmin       = 0;
    tibeg      = (size_t)(opt->rscan_begin); 
    
    first_elim = true;

    s_mass          = rs_sample_create(npar,NBOOT);
    sigmass         = mat_create(npar,1);
   
    mass = rs_sample_pt_cent_val(s_mass);
    fit_data_mass_fit_tune(d,mass,mprop,em,sigem,h[0]->parity);
    if (opt->rscan_begin < 0)
    {
        qcd_printf(opt,"-- fitting and resampling %s mass...\n",full_name);
    }
    else
    {
        qcd_printf(opt,"-- scanning ranges [ti,%d] from ti= %d\n",(int)(nt/2),\
                   opt->rscan_begin);
        opt->nmanrange   = 1;
        opt->range[0][0] = (unsigned int)tibeg;
        opt->range[0][1] = (unsigned int)nt/2;
    }
    strbufcpy(range_info,"");
    if (opt->nmanrange >= 1)
    {
        qcd_printf(opt,"automatic fit range(s) overridden with manual range(s) : ");
        fit_data_fit_all_points(d,false);
        for (i=0;i<opt->nmanrange;i++)
        {
            
            inrmin = (size_t)opt->range[i][0];
            inrmax = (size_t)opt->range[i][1];
            if (inrmin <= nt/2)
            {
                rmin = inrmin;
                strbufcpy(buf,"");
                for (rmax=rmin;rmax<=inrmax;rmax++)
                {
                    rerr = mat_get(sigmprop,rmax,0)/fabs(mat_get(mprop,\
                                                                 rmax,0));
                    if (rerr > MAX_REL_ERR)
                    {
                        strbufcpy(buf," (rcut)");
                        break;
                    }
                }
                rmax--;
            }
            else if (inrmin >= nt/2)
            {
                rmax = inrmax;
                strbufcpy(buf,"");
                for (rmin=rmax;rmin>=inrmin;rmin--)
                {
                    rerr = mat_get(sigmprop,rmin,0)/fabs(mat_get(mprop,\
                                                                 rmin,0));
                    if (rerr > MAX_REL_ERR)
                    {
                        strbufcpy(buf," (lcut)");
                        break;
                    }
                }
                rmin++;
            }
            qcd_printf(opt,"[%d,%d]%s ",(int)rmin,(int)rmax,buf);
            fit_data_fit_range(d,rmin,rmax,true);
            sprintf(buf,"_%d_%d",(int)inrmin,(int)inrmax);
            strcat(range_info,buf);
        }
        printf("\n");
    }
    nti             = MAX(rmax - 1 - tibeg,1);
    scanres_t       = mat_create(nti,1);
    scanres_chi2    = mat_create(nti,1);
    scanres_mass    = mat_create(nti,1);
    scanres_masserr = mat_create(nti,1);
    if (opt->rscan_begin < 0)
    {
        sprintf(sample_name,"%s_mass_fit%s_%s.boot",full_name,range_info,\
                manf_name);
        rs_sample_set_name(s_mass,sample_name);
        rs_data_fit(s_mass,NULL,&s_mprop,d,NO_COR,NULL);
        rs_sample_varp(sigmass,s_mass);
        mat_eqsqrt(sigmass);
        uc_mass     = mat_get(mass,0,0);
        uc_mass_err = mat_get(sigmass,0,0);
        rs_data_fit(s_mass,NULL,&s_mprop,d,opt->corr,NULL);
        if (opt->do_save_rs_sample)
        {
            rs_sample_save(s_mass->name,'w',s_mass);
        }
        rs_sample_varp(sigmass,s_mass);
        mat_eqsqrt(sigmass);
        mat_eqmuls(em,1.0/latspac_nu);
        mat_eqmuls(sigem,1.0/latspac_nu);
        dev = fabs(mat_get(mass,0,0)-uc_mass)/uc_mass_err;
        if (dev > 1.0)
        {
            fprintf(stderr,"warning: correlated mass moved %.1f sigmas away from the uncorrelated mass\n",dev);
        }
        qcd_printf(opt,"M_%-10s= %.8f +/- %.8e %s\n",full_name,\
                   mat_get(mass,0,0)/latspac_nu,             \
                   mat_get(sigmass,0,0)/latspac_nu,unit);
        qcd_printf(opt,"%-12s= %d\n","dof",fit_data_get_dof(d));
        qcd_printf(opt,"%-12s= %e\n","chi^2/dof",fit_data_get_chi2pdof(d));
    }
    else
    {
        printf("%-5s %-12s a*M_%-8s %-12s","ti/a","chi^2/dof",full_name,\
               "error");
        best_t = 0;
        for (i=tibeg;i<rmax-1;i++)
        {
            rs_data_fit(s_mass,NULL,&s_mprop,d,opt->corr,NULL);
            rs_sample_varp(sigmass,s_mass);
            mat_eqsqrt(sigmass);
            mat_set(scanres_t,i-tibeg,0,(double)(i));
            mat_set(scanres_chi2,i-tibeg,0,fit_data_get_chi2pdof(d));
            mat_set(scanres_mass,i-tibeg,0,mat_get(mass,0,0));
            mat_set(scanres_masserr,i-tibeg,0,mat_get(sigmass,0,0));
            if ((i>tibeg)&&(best_t==0)\
                &&(fit_data_get_chi2pdof(d)>mat_get(scanres_chi2,i-tibeg-1,0)))
            {
                best_t = i-1;
                printf(" BEST TIME");
            }
            printf("\n%5d % .5e % .5e % .5e",(int)(i),fit_data_get_chi2pdof(d),\
                   mat_get(mass,0,0),mat_get(sigmass,0,0));
            fit_data_fit_point(d,i,false);
        }
        printf("\n\n");
    }
    
    /*                  plot                    */
    /********************************************/
    if (opt->do_plot)
    {
        mat *em_t,*pr_t;
        plot *p;
        strbuf key,plotcmd;
        size_t maxt;
        double dmaxt,shift,abs_mass;
        
        maxt  = nt-1;
        dmaxt = (double)maxt;
        shift = (h[0]->parity == EVEN) ? 0.0 : -DRATIO(nt,2.0);
        
        em_t = mat_create(nt-2,1);
        pr_t = mat_create(nt,1);
        
        if (opt->rscan_begin < 0)
        {
            /* propagator plot */
            p = plot_create();

            plot_set_scale_ylog(p);
            plot_set_scale_xmanual(p,0,dmaxt);
            sprintf(key,"%s propagator",full_name);
            mat_eqabs(mprop);
            mat_set_step(pr_t,0.0,1.0);
            plot_add_dat(p,pr_t,mprop,NULL,sigmprop,key,"rgb 'red'");
            switch (h[0]->parity)
            {
                case EVEN:
                    sprintf(plotcmd,"exp(-%e*x+%e)",\
                            mat_get(mass,0,0),mat_get(mass,1,0));
                    break;
                case ODD:
                    sprintf(plotcmd,"cosh(%e*(x-%e))*exp(%e)",                   \
                            mat_get(mass,0,0),DRATIO(nt,2),mat_get(mass,1,0));
                    break;
            }

            strcat(plotcmd," t 'fit' lc rgb 'red'");
            plot_add_plot(p,plotcmd);
            plot_disp(p);   
        
            plot_destroy(p);
            
            /* effective mass plot */
            p = plot_create();
            
            abs_mass = fabs(mat_get(mass,0,0))/latspac_nu;
            plot_set_scale_manual(p,0.0,nt-1,0.0,2.0*abs_mass);
            plot_add_hlineerr(p,abs_mass, mat_get(sigmass,0,0)/latspac_nu,\
                              "rgb 'red'");
            sprintf(key,"%s effective mass",full_name);
            mat_set_step(em_t,1.0,1.0);
            plot_add_dat(p,em_t,em,NULL,sigem,key,"rgb 'blue'");
            plot_disp(p);
            
            plot_destroy(p);
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
            
            plot_destroy(p);
            
            /* mass plot */
            p = plot_create();
            
            plot_set_scale_xmanual(p,0,(double)(nt/2));
            sprintf(key,"a*M_%s",full_name);
            plot_add_dat(p,scanres_t,scanres_mass,NULL,scanres_masserr,key,\
                         "rgb 'red'");
            plot_disp(p);
            
            plot_destroy(p);
        }

        mat_destroy(em_t);
        mat_destroy(pr_t);
    }
    
    /*              desallocation               */
    /********************************************/
    free(opt);
    spectrum_destroy(s); 
    io_finish();
    mat_ar_destroy(prop[0],nbdat);
    mat_ar_destroy(prop[1],nbdat);
    rs_sample_destroy(s_mprop);
    rs_sample_destroy(s_tmp);
    mat_destroy(sigmprop);
    rs_sample_destroy(s_effmass);
    mat_destroy(sigem);
    fit_data_destroy(d);
    rs_sample_destroy(s_mass);
    mat_destroy(sigmass);
    mat_destroy(scanres_t);
    mat_destroy(scanres_chi2);
    mat_destroy(scanres_mass);
    mat_destroy(scanres_masserr);
    
return EXIT_SUCCESS;
}
