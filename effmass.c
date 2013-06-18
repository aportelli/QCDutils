#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <qcd_arg_parse.h>
#include <latan/latan_io.h>
#include <latan/latan_mat.h>
#include <latan/latan_mass.h>
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
    int emtype;
     
    opt = qcd_arg_parse(argc,argv,A_SAVE_RS|A_LOAD_RS|A_LOAD_RG|A_PROP_NAME|\
                        A_PROP_LOAD|A_QCOMP|A_PLOT|A_MODEL,1,1);
    sprintf(prop_name,"%s_%s_%s_%s",opt->channel[0],opt->quark[0],opt->sink,\
            opt->source);
    sprintf(full_name,"%s_%s",opt->channel[0],opt->quark[0]);
    if (opt->do_load_rs_sample)
    {
        strbufcpy(manf_name,opt->load_rs_fname[0]);
    }
    else
    {
        strbufcpy(manf_name,opt->manf_name);
    }
    binsize = opt->binsize;
    nboot   = opt->nboot;
    latan_set_verb(opt->latan_verb);
    if (strbufcmp(opt->model,"cosh") == 0)
    {
        emtype   = EM_ACOSH;
    }
    else
    {
        strbufcpy(opt->model,"expdec");
        emtype   = EM_LOG;
    }
    io_set_fmt(opt->latan_fmt);
    io_init();
    
    /*              loading datas               */
    /********************************************/
    size_t ndat,nbdat,dim[2];
    mat **prop;
    
    if (opt->do_load_rs_sample)
    {
        ndat  = 0;
        nbdat = 0;
        rs_sample_load(NULL,NULL,dim,opt->load_rs_fname[0]);
        prop  = NULL;
    }
    else
    {
        ndat    = (size_t)get_nfile(manf_name);
        nbdat   = ndat/binsize + ((ndat%binsize == 0) ? 0 : 1);
        mat_ar_loadbin(NULL,dim,manf_name,prop_name,1);
        prop = mat_ar_create(nbdat,dim[0],dim[1]);
        qcd_printf(opt,"-- loading %s datas from %s...\n",prop_name,manf_name);
        mat_ar_loadbin(prop,NULL,manf_name,prop_name,binsize);
    }
    
    
    /*                 propagator               */
    /********************************************/
    rs_sample *s_mprop;
    
    s_mprop = rs_sample_create(dim[0],dim[1],nboot);
    
    if (opt->do_load_rs_sample)
    {
        qcd_printf(opt,"-- loading %s propagator from %s...\n",prop_name,\
                   opt->load_rs_fname[0]);
        rs_sample_load(s_mprop,NULL,NULL,opt->load_rs_fname[0]);
    }
    else
    {
        qcd_printf(opt,"-- resampling %s mean propagator...\n",full_name);
        randgen_set_state(opt->state);
        resample(s_mprop,prop,nbdat,&rs_mean,BOOT,NULL);
    }
    
    /*           effective mass                 */
    /********************************************/
    rs_sample *s_effmass;
    mat *tem,*em,*sigem;
    size_t emdim[2];
    
    get_effmass_size(emdim,rs_sample_pt_cent_val(s_mprop),1,emtype);
    
    s_effmass = rs_sample_create(emdim[0],emdim[1],nboot);
    tem       = mat_create(emdim[0],1);
    sigem     = mat_create(emdim[0],emdim[1]);
    
    qcd_printf(opt,"-- resampling %s effective mass...\n",full_name);
    rs_sample_effmass(s_effmass,tem,s_mprop,1,emtype);
    em = rs_sample_pt_cent_val(s_effmass);
    rs_sample_eqabs(s_effmass);
    rs_sample_varp(sigem,s_effmass);
    mat_eqsqrt(sigem);
    
    /*              result output               */
    /********************************************/
    size_t t_i;
    strbuf latan_path;
    plot *p;
    
    qcd_printf(opt,"\n%-4s %-12s %-12s\n","t","effmass","error");
    for (t_i=0;t_i<emdim[0];t_i++)
    {
        qcd_printf(opt,"% -4d % -.5e % -.5e\n",(int)mat_get(tem,t_i,0),\
                   mat_get(em,t_i,0),mat_get(sigem,t_i,0));
    }
    qcd_printf(opt,"\n");
    if (opt->do_save_rs_sample)
    {
        sprintf(latan_path,"%s_effmass_%s.boot:%s_effmass_%s",full_name,manf_name,\
                full_name,manf_name);
        rs_sample_save(latan_path,'w',s_effmass);
    }
    if (opt->do_plot|opt->do_save_plot)
    {
        p    = plot_create();
        plot_add_dat_yerr(p,tem,em,sigem,full_name,"rgb 'red'");
        plot_add_dat(p,tem,em,NULL,NULL,"","rgb 'red'");
        if (opt->do_plot)
        {
            plot_disp(p);
        }
        if (opt->do_save_plot)
        {
            plot_save(opt->save_plot_dir,p);
        }
        plot_destroy(p);
    }
    
    /*              desallocation               */
    /********************************************/
    FREE(opt);
    io_finish();
    mat_ar_destroy(prop,nbdat);
    rs_sample_destroy(s_mprop);
    rs_sample_destroy(s_effmass);
    mat_destroy(tem);
    mat_destroy(sigem);
    
    return EXIT_SUCCESS;
}
