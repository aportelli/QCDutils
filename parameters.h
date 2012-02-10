#ifndef PARAMETERS_H_
#define PARAMETERS_H_

#include <latan/latan_fit.h>

typedef struct fit_init_s
{
    size_t ind;
    double value;
} fit_init;

typedef struct ens_s
{
    size_t T;
    size_t L;
    strbuf beta;
    strbuf ud;
    strbuf s;
    strbuf dataset;
    strbuf dir;
} ens;

typedef struct fit_param_s
{
    int M_ud_deg;
    int s_M_ud_deg;
    int M_s_deg;
    int s_M_s_deg;
    int a_deg;
    int s_with_aM_s;
    int with_umd;
    int s_with_umd;
    int have_umd;
    int with_qed_fvol;
    int s_with_qed_fvol;
    int with_ext_a;
    int verb;
    int plotting;
    int q_dim;
    fit_init *init_param;
    ens *point;
    double M_ud;
    double M_s;
    double M_scale;
    double q_target[2];
    strbuf analyze;
    strbuf q_name;
    strbuf scale_part;
    strbuf ud_name;
    strbuf s_name;
    strbuf *dataset;
    strbuf dataset_cat;
    strbuf manifest;
    strbuf *beta;
    size_t ndataset;
    size_t nbeta;
    size_t nens;
    size_t nsample;
    size_t ninit_param;
    const fit_model *model;
} fit_param;

#define IS_ANALYZE(param,prog) (strcmp((param)->analyze,prog)==0)

int ind_beta(const strbuf beta, const fit_param *param);
fit_param * fit_param_parse(const strbuf fname);
void fit_param_destroy(fit_param *param);

#endif
