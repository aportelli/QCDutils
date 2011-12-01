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
    int M_s_deg;
    int a_deg;
    int with_umd;
    int with_qed_fvol;
    int with_ext_a;
    int verb;
    int q_dim;
    fit_init *init_param;
    double M_ud;
    double M_s;
    double M_scale;
    strbuf analyze;
    strbuf q_name;
    strbuf scale_part;
    strbuf ud_name;
    strbuf s_name;
    strbuf dataset;
    strbuf manifest;
    strbuf *beta;
    size_t nbeta;
    size_t nens;
    size_t ninit_param;
    const fit_model *model;
} fit_param;

#define IS_ANALYZE(param,prog) (strcmp((param)->analyze,prog)==0)

void parse_fit_param(fit_param *param, const strbuf fname);
void add_beta(fit_param *param, const strbuf new_beta);
int ind_beta(const strbuf beta, const fit_param *param);

#endif