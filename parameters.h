#ifndef PARAMETERS_H_
#define PARAMETERS_H_

#include <latan/latan_fit.h>

typedef struct
{
    size_t ind;
    double value;
} fit_init;

typedef struct
{
    int M_ud_deg;
    int s_M_ud_deg;
    int M_s_deg;
    int s_M_s_deg;
    int a_deg;
    int with_umd;
    int s_with_umd;
    int with_qed_fvol;
    int s_with_qed_fvol;
    int with_ext_a;
    int verb;
    int plotting;
    int q_dim;
    fit_init *init_param;
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
void parse_fit_param(fit_param *param, const strbuf fname);

#endif
