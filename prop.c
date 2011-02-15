#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <qcd_arg_parse.h>
#include <latan/latan_hadron.h>
#include <latan/latan_io.h>
#include <latan/latan_mat.h>
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
    qcd_options *opt;
    strbuf spec_name,part_name,manf_name;
    size_t binsize;
     
    opt = qcd_arg_parse(argc,argv,A_PARTICLE|A_PROP_LOAD|A_SAVE_RS\
                        |A_CHANNEL|A_LOAD_RG);
    strbufcpy(spec_name,opt->spec_name);
    strbufcpy(part_name,opt->part_name);
    strbufcpy(manf_name,opt->manf_name);
    source = opt->source;
    sink = opt->sink;
    binsize = opt->binsize;
    for (ch=0;ch<NCHANNEL;ch++)
    {
        channel_id_set(ch,opt->channel_id[ch]);
    }
    latan_set_verb(opt->latan_verb);
    
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
    size_t ndat,nt;
    mat **prop;
    
    ndat    = (size_t)get_nfile(manf_name);
    hadron_prop_load_nt(&nt,h,source,sink,manf_name);
    
    prop = mat_ar_create(ndat,nt,1);

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
    resample(s_mprop,prop,ndat,1,&rs_mean,BOOT,NULL);
    mprop = rs_sample_pt_cent_val(s_mprop);
    
    /*  computing error on mean propagator      */
    /********************************************/
    mat *sig;
    
    sig = mat_create(nt,1);
    
    qcd_printf(opt,"-- estimating %s mean propagator error...\n",h->name);
    rs_sample_varp(sig,s_mprop);
    mat_eqsqrt(sig);
    
    /*              result output               */
    /********************************************/
    size_t t,maxt;
    
    maxt = (h->parity == EVEN) ? (nt/2) : (nt-1);
    qcd_printf(opt,"\nt\tprop\t\terror\n");
    for (t=0;t<=maxt;t++)
    {
        qcd_printf(opt,"%d\t%e\t%e\n",(int)t,mat_get(mprop,t,0),mat_get(sig,t,0));
    }
    qcd_printf(opt,"\n");
    if (opt->do_save_rs_sample)
    {
        rs_sample_save(s_mprop->name,'w',s_mprop);
    }
    
    /*              desallocation               */
    /********************************************/
    FREE(opt);
    spectrum_destroy(s);
    /* io_finish(); */ /* unknown memory leak here */
    mat_ar_destroy(prop,ndat);
    rs_sample_destroy(s_mprop);
    mat_destroy(sig);
    
    return EXIT_SUCCESS;
}
