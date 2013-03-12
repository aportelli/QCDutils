#ifndef QCD_ARG_PARSE_H_
#define QCD_ARG_PARSE_H_

#include <sys/types.h>
#include <latan/latan_fit.h>
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
    A_PROP_NAME = 1 << 3,
    A_PROP_LOAD = 1 << 4,
    A_LATSPAC   = 1 << 5,
    A_QCOMP     = 1 << 6,
    A_FIT       = 1 << 7,
    A_MODEL     = 1 << 8
};

typedef struct qcd_options_s
{
    bool qcd_verb;
    bool have_latspac;
    bool do_plot;
    bool do_save_plot;
    bool do_save_rs_sample;
    bool do_range_scan;
    bool have_randgen_state;
    int latan_verb;
    io_fmt_no latan_fmt;
    double latspac_fm;
    double latspac_nu;
    strbuf manf_name;
    size_t binsize;
    rg_state state;
    strbuf source;
    strbuf sink;
    strbuf ss;
    strbuf channel[2];
    strbuf quark[2];
    strbuf model;
    strbuf save_plot_dir;
    cor_flag corr;
    minalg_no minimizer;
    unsigned int range[MAX_RANGES][2];
    double tshift;
    size_t nmanrange;
    size_t nboot;
} qcd_options;

qcd_options * qcd_arg_parse(int argc, char* argv[], unsigned int argset_flag,
                            const int max_npart);
void qcd_printf(const qcd_options *opt, const strbuf fmt, ...);

#endif
