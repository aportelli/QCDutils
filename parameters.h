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
    int Msq_pi_deg;
    int Msq_K_deg;
    int a_deg;
    int with_umd;
    int with_fvol;
    int verb;
    int ex_dim;
    int q_dim;
    fit_param *init_param;
    double M_scale[2];
    double Mpi_cut;
    double MpiL_cut;
    strbuf analyze;
    strbuf q_name;
    strbuf scale_part;
    strbuf pi_name;
    strbuf K_name;
    strbuf suffix;
    strbuf manifest;
    strbuf *beta;
    size_t nbeta;
    size_t nens;
    size_t ninit_param;
    const fit_model *model;
} ex_param;

void parse_ex_param(ex_param *param, const strbuf fname);
void add_beta(ex_param *param, const strbuf new_beta);
int ind_beta(const strbuf beta, const ex_param *param);

#endif