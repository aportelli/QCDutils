#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include "config.h"
#include "qcd_arg_parse.h"
#include <argtable2.h>
#include <latan/latan_math.h>
#include <latan/latan_nunits.h>

#define A_MAX_NERROR 20
#define DEF_NBOOT 100

qcd_options * qcd_arg_parse(int argc, char* argv[], unsigned int argset_flag,
                            const int max_npart)
{
    /*          argument definitions            */
    /********************************************/
    void** a_table;
    size_t a_table_size,i;
    struct arg_lit*  a_help          = NULL;
    struct arg_lit*  a_ver           = NULL;
    struct arg_int*  a_verb          = NULL;
    struct arg_str*  a_fmt           = NULL;
    struct arg_int*  a_nboot         = NULL;
    struct arg_lit*  a_save_rs       = NULL;
    struct arg_str*  a_load_rg       = NULL;
    struct arg_lit*  a_plot          = NULL;
    struct arg_dbl*  a_latspac_fm    = NULL;
    struct arg_str*  a_qcomp         = NULL;
    struct arg_str*  a_channel       = NULL;
    struct arg_str*  a_ss            = NULL;
    struct arg_int*  a_binsize       = NULL;
    struct arg_str*  a_minimizer     = NULL;
    struct arg_str*  a_range         = NULL;
    struct arg_lit*  a_uncorr        = NULL;
    struct arg_lit*  a_rscan         = NULL;
    struct arg_str*  a_model         = NULL;
    struct arg_file* a_manf          = NULL;
    struct arg_end*  a_end           = NULL;
    strbuf help_msg,ver_msg,verb_msg,fmt_msg,nboot_msg,save_rs_msg,load_rg_msg,\
           plot_msg,latspac_fm_msg,qcomp_msg,channel_msg,ss_msg,               \
           binsize_msg,minimizer_msg,range_msg,uncorr_msg,rscan_msg,           \
           manf_msg,model_msg;
    strbuf defmin,deffmt;
    int j;
    
    minalg_id_get(defmin,minimizer_get_alg());
    switch (io_get_fmt())
    {
        case IO_ASCII:
            strbufcpy(deffmt,"ASCII");
            break;
        case IO_XML:
            strbufcpy(deffmt,"XML");
            break;
        default:
            fprintf(stderr,"error: default LatAnalyze format unknown\n");
            exit(EXIT_FAILURE);
            break;
    }
    
    sprintf(help_msg        ,"display this help"                            );
    sprintf(ver_msg         ,"display QCDutils version"                     );
    sprintf(verb_msg        ,"verbosity level (default: 1)"                 );
    sprintf(fmt_msg         ,"I/O format (default: %s)",deffmt              );
    sprintf(nboot_msg       ,"number of bootstrap samples (default: %d)",\
            DEF_NBOOT);
    sprintf(save_rs_msg     ,"save resampled samples"                       );
    sprintf(load_rg_msg     ,"use saved random generator state"             );
    sprintf(plot_msg        ,"show plot"                                    );
    sprintf(latspac_fm_msg  ,"lattice spacing (in fm)"                      );
    sprintf(qcomp_msg       ,"quark composition"                            );
    sprintf(channel_msg     ,"channel"                                      );
    sprintf(ss_msg          ,"particle source and sink (e.g 2:2)"           );
    sprintf(binsize_msg     ,"data binning size (default: 1)"               );
    sprintf(minimizer_msg   ,"minimizer (default: %s)",defmin               );
    sprintf(range_msg       ,"manual fit range"                             );
    sprintf(uncorr_msg      ,"use time-uncorrelated chi^2"                  );
    sprintf(rscan_msg       ,"perform fit range scan"                       );
    sprintf(manf_msg        ,"LatAnalyze data manifest"                     );
    sprintf(model_msg       ,"model name"                                   );

    a_help  = arg_lit0(NULL,"help",help_msg);
    a_ver   = arg_lit0(NULL,"version",ver_msg);
    a_verb  = arg_int0("v","verb",NULL,verb_msg);
    a_fmt   = arg_str0(NULL,"format","{ascii|xml}",fmt_msg);
    a_nboot = arg_int0("n","nboot",NULL,nboot_msg);
    if (argset_flag & A_SAVE_RS)
    {
        a_save_rs       = arg_lit0(NULL,"save",save_rs_msg);
    }
    if (argset_flag & A_LOAD_RG)
    {
        a_load_rg       = arg_str0("r","load_rg","FNAME",load_rg_msg);
    }
    if (argset_flag & A_PLOT)
    {   
        a_plot        = arg_lit0(NULL,"plot",plot_msg);
    }
    if (argset_flag & A_LATSPAC)
    {
        a_latspac_fm    = arg_dbl0("L","lat_spac",NULL,latspac_fm_msg);
    }
    if (argset_flag & A_QCOMP)
    {
        a_qcomp         = arg_strn("q","quarks","QUARKS",1,max_npart,qcomp_msg);
    }
    if (argset_flag & A_PROP_LOAD)
    {
        a_channel       = arg_strn("c","channel","CHNAME",1,max_npart,\
                                   channel_msg);
        a_ss            = arg_str1("S","sink_source","SI:SO",\
                                   ss_msg);
        a_binsize       = arg_int0("b","bin",NULL,binsize_msg);
        a_manf          = arg_file1(NULL,NULL,"<manifest file>",manf_msg);
    }
    if (argset_flag & A_FIT)
    {
        a_minimizer     = arg_str0("M","minimizer","ID",minimizer_msg);
        a_range         = arg_strn("R","range","[min,max]",1,MAX_RANGES,\
                                   range_msg);
        a_uncorr        = arg_lit0(NULL,"uncorr",uncorr_msg);
        a_rscan         = arg_lit0("f","range_scan",rscan_msg);
    }
    if (argset_flag & A_MODEL)
    {
        a_model         = arg_str0("m","model","NAME",model_msg);
    }
    a_end   = arg_end(A_MAX_NERROR);
    
    a_table_size = 5;
    QCD_MALLOC(a_table,void**,a_table_size);
    a_table[0] = a_help;
    a_table[1] = a_ver;
    a_table[2] = a_verb;
    a_table[3] = a_fmt;
    a_table[4] = a_nboot;
    i = 5;
    if (argset_flag & A_SAVE_RS)
    {
        a_table_size += 1;
        QCD_REALLOC(a_table,a_table,void**,a_table_size);
        a_table[i] = a_save_rs;
        i += 1;
    }
    if (argset_flag & A_LOAD_RG)
    {
        a_table_size += 1;
        QCD_REALLOC(a_table,a_table,void**,a_table_size);
        a_table[i] = a_load_rg;
        i += 1;
    }
    if (argset_flag & A_PLOT)
    {
        a_table_size += 1;
        QCD_REALLOC(a_table,a_table,void**,a_table_size);
        a_table[i] = a_plot;
        i += 1;
    }
    if (argset_flag & A_LATSPAC)
    {
        a_table_size += 1;
        QCD_REALLOC(a_table,a_table,void**,a_table_size);
        a_table[i] = a_latspac_fm;
        i += 1;
    }
    if (argset_flag & A_QCOMP)
    {
        a_table_size += 1;
        QCD_REALLOC(a_table,a_table,void**,a_table_size);
        a_table[i] = a_qcomp;
        i += 1;
    }
    if (argset_flag & A_PROP_NAME)
    {
        a_table_size += 2;
        QCD_REALLOC(a_table,a_table,void**,a_table_size);
        a_table[i]   = a_channel;
        a_table[i+1] = a_ss;
        i += 2;
    }
    if (argset_flag & A_PROP_LOAD)
    {
        a_table_size += 2;
        QCD_REALLOC(a_table,a_table,void**,a_table_size);
        a_table[i]   = a_binsize;
        a_table[i+1] = a_manf;
        i += 2;
    }
    if (argset_flag & A_FIT)
    {
        a_table_size += 4;
        QCD_REALLOC(a_table,a_table,void**,a_table_size);
        a_table[i]   = a_minimizer;
        a_table[i+1] = a_range;
        a_table[i+2] = a_uncorr;
        a_table[i+3] = a_rscan;
        i += 4;
    }
    if (argset_flag & A_MODEL)
    {
        a_table_size += 1;
        QCD_REALLOC(a_table,a_table,void**,a_table_size);
        a_table[i] = a_model;
        i += 1;
    }
    a_table_size += 1;
    QCD_REALLOC(a_table,a_table,void**,a_table_size);
    a_table[i] = a_end;

    /*          check arg table allocation      */
    /********************************************/
    if (arg_nullcheck(a_table) != 0)
    {
        fprintf(stderr,"error: argument table allocation failed\n");
        exit(EXIT_FAILURE);
    }
    
    /*      option structure allocation         */
    /********************************************/
    qcd_options *opt;
    
    QCD_MALLOC(opt,qcd_options *,1);
    
    /*          set default options             */
    /********************************************/
    if (argset_flag & A_LOAD_RG)
    {
        opt->have_randgen_state = false;
        randgen_init_from_time();
        randgen_get_state(opt->state);
    }
    if (argset_flag & A_LATSPAC)
    {
        opt->latspac_fm = 1.0;
        opt->latspac_nu = 1.0;
        opt->have_latspac = false;
    }
    if (argset_flag & A_QCOMP)
    {
        strbufcpy(opt->quark[0],"");
        strbufcpy(opt->quark[1],"");
    }
    if (argset_flag & A_PROP_NAME)
    {
        strbufcpy(opt->channel[0],"");
        strbufcpy(opt->channel[1],"");
    }
    if (argset_flag & A_PROP_LOAD)
    {
        opt->binsize = 1;
    }
    if (argset_flag & A_FIT)
    {
        opt->minimizer     = minalg_no_get(defmin);
        opt->corr          = DATA_COR;
        opt->do_range_scan = false;
    }
    if (argset_flag & A_MODEL)
    {
        strbufcpy(opt->model,"");
    }
    opt->nboot      = DEF_NBOOT;
    opt->latan_verb = QUIET;
    opt->qcd_verb   = true;
    opt->latan_fmt  = io_get_fmt();
    
    /*          argument parsing                */
    /********************************************/
    int nerror;
    FILE *fpt;
    char *cpt;
    
    cpt    = NULL;
    fpt    = NULL;
    nerror = arg_parse(argc,argv,a_table);
    
    if ((a_help->count > 0)||(argc == 1))
    {
        if (a_help->count > 0)
        {
            fpt = stdout;
        }
        else
        {
            fpt = stderr;
        }
        fprintf(fpt,"usage: %s",argv[0]);
        arg_print_syntax(fpt,a_table,"\n");
        arg_print_glossary_gnu(fpt,a_table);
        exit(EXIT_SUCCESS);
    }
    if (a_ver->count > 0)
    {
        printf("%s is part of %s v%s\n",argv[0],PACKAGE_NAME,PACKAGE_VERSION);
        printf("Report any bugs or issues to %s\n",PACKAGE_BUGREPORT);
        exit(EXIT_SUCCESS);
    }
    if (nerror > 0)
    {
        arg_print_errors(stdout,a_end,argv[0]);
        printf("Try '%s --help' for more information.\n",argv[0]);
        exit(EXIT_FAILURE);
    }
    if (a_verb->count > 0)
    {
        opt->qcd_verb   = (a_verb->ival[0] > 0);
        opt->latan_verb = MAX(0,a_verb->ival[0]-1);
    }
    if (a_fmt->count > 0)
    {
        if (strcmp(a_fmt->sval[0],"ascii") == 0)
        {
            opt->latan_fmt = IO_ASCII;
        }
        else if (strcmp(a_fmt->sval[0],"xml") == 0)
        {
            opt->latan_fmt = IO_XML;
        }
        else
        {
            fprintf(stderr,"error: I/O format %s unknown\n",a_fmt->sval[0]);
            exit(EXIT_FAILURE);
        }
    }
    if (a_nboot->count > 0)
    {
        opt->nboot = a_nboot->ival[0];
    }
    if (argset_flag & A_SAVE_RS)
    {
        opt->do_save_rs_sample = (a_save_rs->count > 0);
    }
    if (argset_flag & A_LOAD_RG)
    {
        if (a_load_rg->count > 0)
        {
            opt->have_randgen_state = true;
            randgen_load_state(opt->state,a_load_rg->sval[0]);
        }
    }
    if (argset_flag & A_PLOT)
    {
        opt->do_plot = (a_plot->count > 0);
    }
    if (argset_flag & A_LATSPAC)
    {
        if (a_latspac_fm->count > 0)
        {
            opt->latspac_fm = a_latspac_fm->dval[0];
            opt->latspac_nu = opt->latspac_fm*NU_FM;
            opt->have_latspac = true;
        }
    }
    if (argset_flag & A_QCOMP)
    {
        for (j=0;j<a_qcomp->count;j++)
        {
            strbufcpy(opt->quark[j],a_qcomp->sval[j]);
        }
    }
    if (argset_flag & A_PROP_NAME)
    {
        for (j=0;j<a_channel->count;j++)
        {
            strbufcpy(opt->channel[j],a_channel->sval[j]);
        }
        if (argset_flag & A_QCOMP)
        {
            if ((a_qcomp->count == 2)&&(a_channel->count == 1))
            {
                strbufcpy(opt->channel[1],opt->channel[0]);
            }
            else if ((a_qcomp->count == 1)&&(a_channel->count == 2))
            {
                strbufcpy(opt->quark[1],opt->quark[0]);
            }
        }
        cpt = strchr(a_ss->sval[0],':');
        if (cpt == NULL)
        {
            fprintf(stderr,"error: sink/source option %s is invalid\n",\
                    a_ss->sval[0]);
            exit(EXIT_FAILURE);
        }
        strbufcpy(opt->source,cpt+1);
        *cpt = '\0';
        strbufcpy(opt->sink,a_ss->sval[0]);
    }
    if (argset_flag & A_PROP_LOAD)
    {
        if (a_binsize->count)
        {
            opt->binsize = (size_t)(a_binsize->ival[0]);
        }
        strbufcpy(opt->manf_name,a_manf->filename[0]);
    }
    if (argset_flag & A_FIT)
    {
        if (a_minimizer->count > 0)
        {
            opt->minimizer = minalg_no_get(a_minimizer->sval[0]);
        }
        opt->nmanrange = a_range->count;
        for (i=0;i<opt->nmanrange;i++)
        {
            if(sscanf(a_range->sval[i],"[%u,%u]",opt->range[i],\
                      opt->range[i]+1) <= 0)
            {
                fprintf(stderr,"error: range option %s invalid\n",\
                        a_range->sval[i]);
                exit(EXIT_FAILURE);
            }
        }
        if (a_uncorr->count > 0)
        {
            opt->corr = NO_COR;
        }
        opt->do_range_scan = (a_rscan->count > 0);
    }
    if (argset_flag & A_MODEL)
    {
        if (a_model->count > 0)
        {
            strbufcpy(opt->model,a_model->sval[0]);
        }
    }
    
    /*              desallocation               */
    /********************************************/
    arg_freetable(a_table,a_table_size);
    free(a_table);

    return opt;
}

void qcd_printf(const qcd_options *opt, const strbuf fmt, ...)
{
    va_list args;
    
    if (opt->qcd_verb)
    {
        va_start(args,fmt);
        vprintf(fmt,args);
        va_end(args);
    }
}
