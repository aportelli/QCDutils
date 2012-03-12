#ifndef QCD_ARG_PARSE_H_
#define QCD_ARG_PARSE_H_

#include <sys/types.h>
#include <latan/latan_fit.h>
#include <latan/latan_hadron.h>
#include <latan/latan_io.h>
#include <latan/latan_minimizer.h>
#include <latan/latan_rand.h>

#define MAX_RANGES 4

#define QCD_MALLOC(pt,typ,size)\
{\
    pt = (typ)(malloc((size_t)(size)*sizeof(*pt)));\
    if (pt == NULL)\
    {\
        fprintf(stderr,"error: memory allocation failed (%s:%d)\n",__FILE__,\
                __LINE__);\
        abort();\
    }\
}

#define QCD_REALLOC(pt,pt_old,typ,size)\
{\
    pt = (typ)(realloc(pt_old,(size_t)(size)*sizeof(*pt)));\
    if (pt == NULL)\
    {\
        fprintf(stderr,"error: memory reallocation failed (%s:%d)\n",__FILE__,\
                __LINE__);\
        abort();\
    }\
}

enum
{
    A_DEFAULT   = 0,
    A_PLOT      = 1 << 0,
    A_SAVE_RS   = 1 << 1,
    A_LOAD_RG   = 1 << 2,
    A_PROP_LOAD = 1 << 3,
    A_LATSPAC   = 1 << 4,
    A_PARTICLE  = 1 << 5,
    A_CHANNEL   = 1 << 6,
    A_QCOMP     = 1 << 7,
    A_FIT       = 1 << 8
};

typedef struct
{
    bool qcd_verb;
    int latan_verb;
    io_fmt_no latan_fmt;
    double latspac_fm;
    double latspac_nu;
    bool have_latspac;
    quark_no qcomp[2];
    strbuf qcomp_str;
    strbuf spec_name;
    strbuf part_name[2];
    strbuf manf_name;
    rg_state state;
    ss_no source;
    ss_no sink;
    size_t binsize;
    strbuf channel_id[NCHANNEL];
    strbuf quark_id[NQUARK];
    bool do_plot;
    bool do_save_rs_sample;
    bool have_randgen_state;
    cor_flag corr;
    minalg_no minimizer;
    unsigned int range[4][2];
    size_t nmanrange;
    int rscan_begin;
} qcd_options;

qcd_options * qcd_arg_parse(int argc, char* argv[], int argset_flag);
void qcd_printf(const qcd_options *opt, const strbuf fmt, ...);

#endif
