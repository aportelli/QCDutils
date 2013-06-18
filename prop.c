#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <qcd_arg_parse.h>
#include <latan/latan_io.h>
#include <latan/latan_mat.h>
#include <latan/latan_plot.h>
#include <latan/latan_rand.h>
#include <latan/latan_statistics.h>

int main(int argc, char* argv[])
{
    /*              parsing arguments           */
    /********************************************/
    qcd_options *opt;
    strbuf prop_name,full_name,manf_name;
    size_t binsize,nboot;
     
    opt = qcd_arg_parse(argc,argv,A_SAVE_RS|A_LOAD_RG|A_PROP_NAME|A_PROP_LOAD\
                        |A_QCOMP|A_PLOT,1,0);
    sprintf(prop_name,"%s_%s_%s_%s",opt->channel[0],opt->quark[0],opt->sink,\
            opt->source);
    sprintf(full_name,"%s_%s",opt->channel[0],opt->quark[0]);
    strbufcpy(manf_name,opt->manf_name);
    binsize = opt->binsize;
    nboot   = opt->nboot;
    latan_set_verb(opt->latan_verb);
    io_set_fmt(opt->latan_fmt);
    io_init();
    
    /*              loading datas               */
    /********************************************/
    size_t ndat,nbdat,dim[2],nt;
    mat **prop;
    
    ndat    = (size_t)get_nfile(manf_name);
    nbdat   = ndat/binsize + ((ndat%binsize == 0) ? 0 : 1);
    
    mat_ar_loadbin(NULL,dim,manf_name,prop_name,1);
    prop = mat_ar_create(nbdat,dim[0],dim[1]);
    nt   = dim[0];
    qcd_printf(opt,"-- loading %s datas from %s...\n",prop_name,manf_name);
    mat_ar_loadbin(prop,NULL,manf_name,prop_name,binsize);
    
    /*                 propagator               */
    /********************************************/
    rs_sample *s_mprop;
    mat *mprop,*sig;
    
    s_mprop = rs_sample_create(dim[0],dim[1],nboot);
    sig     = mat_create(dim[0],dim[1]);
    
    qcd_printf(opt,"-- resampling %s mean propagator...\n",full_name);
    randgen_set_state(opt->state);
    resample(s_mprop,prop,nbdat,&rs_mean,BOOT,NULL);
    mprop = rs_sample_pt_cent_val(s_mprop);
    rs_sample_varp(sig,s_mprop);
    mat_eqsqrt(sig);
    
    /*              result output               */
    /********************************************/
    size_t t;
    strbuf latan_path;
    plot *p;
    mat *tvec;
    
    qcd_printf(opt,"\n%-4s %-12s %-12s\n","t","prop","error");
    for (t=0;t<nt;t++)
    {
        qcd_printf(opt,"% -4d % -.5e % -.5e\n",(int)t,mat_get(mprop,t,0),\
                   mat_get(sig,t,0));
    }
    qcd_printf(opt,"\n");
    if (opt->do_save_rs_sample)
    {
        sprintf(latan_path,"%s_prop_%s.boot:%s_prop_%s",full_name,manf_name,\
                full_name,manf_name);
        rs_sample_save(latan_path,'w',s_mprop);
    }
    if (opt->do_plot|opt->do_save_plot)
    {
        p    = plot_create();
        tvec = mat_create(nt,1);
        mat_set_step(tvec,0.0,1.0);
        plot_set_scale_ylog(p);
        plot_add_dat_yerr(p,tvec,mprop,sig,full_name,"rgb 'red'");
        if (opt->do_plot)
        {
            plot_disp(p);
        }
        if (opt->do_save_plot)
        {
            plot_save(opt->save_plot_dir,p);
        }
        mat_destroy(tvec);
        plot_destroy(p);
    }
    
    /*              desallocation               */
    /********************************************/
    FREE(opt);
    io_finish();
    mat_ar_destroy(prop,nbdat);
    rs_sample_destroy(s_mprop);
    mat_destroy(sig);
    
    return EXIT_SUCCESS;
}
