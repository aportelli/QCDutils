#ifndef PARAMETERS_H_
#define PARAMETERS_H_

#include <latan/latan_fit.h>

typedef struct
{
    size_t ind;
    double value;
} fit_param;

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
    fit_param *init_param;
    double M_ud;
    double M_ud_ex;
    double M_s;
    double M_s_ex;
    double M_scale;
    strbuf analyze;
    strbuf q_name;
    strbuf scale_part;
    strbuf ud_name;
    strbuf s_name;
    strbuf suffix;
    strbuf manifest;
    strbuf *beta;
    size_t nbeta;
    size_t nens;
    size_t ninit_param;
    const fit_model *model;
} ex_param;

#define IS_ANALYZE(param,prog) (strcmp((param)->analyze,prog)==0)

void parse_ex_param(ex_param *param, const strbuf fname);
void add_beta(ex_param *param, const strbuf new_beta);
int ind_beta(const strbuf beta, const ex_param *param);

#endif