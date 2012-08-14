#ifndef PARAMETERS_H_
#define PARAMETERS_H_

#include <latan/latan_fit.h>
#include <latan/latan_statistics.h>

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

enum
{
    AN_NOTHING = 0,
    AN_PHYPT   = 1 << 0,
    AN_SCALE   = 1 << 1
};

typedef struct fit_param_s
{
    int with_const;
    int M_ud_deg;
    int s_M_ud_deg;
    int M_s_deg;
    int s_M_s_deg;
    int a_deg;
    int with_aalpha;
    int with_a2M_ud;
    int s_with_a2M_ud;
    int with_a2M_s;
    int s_with_a2M_s;
    int umd_deg;
    int s_umd_deg;
    int have_umd;
    int with_udumd;
    int with_sumd;
    int alpha_deg;
    int s_alpha_deg;
    int have_alpha;
    int with_udalpha;
    int with_salpha;
    int with_qed_fvol;
    int s_with_qed_fvol;
    int with_ext_a;
    int with_ext_M_umd;
    int verb;
    int correlated;
    int save_result;
    int plot;
    int warn_missing_data;
    int scale_model;
    int q_dim;
    unsigned int analyze_flag;
    fit_init *init_param;
    ens *point;
    double M_ud;
    double M_s;
    double M_umd_val;
    double M_scale;
    double alpha;
    double q_target[2];
    strbuf analyze;
    strbuf model;
    strbuf s_model;
    strbuf q_name;
    strbuf s_manifest;
    strbuf scale_part;
    strbuf ud_name;
    strbuf s_name;
    strbuf umd_name;
    strbuf *dataset;
    strbuf dataset_cat;
    strbuf manifest;
    strbuf *beta;
    strbuf save_plot;
    strbuf M_umd;
    size_t ndataset;
    size_t nbeta;
    size_t nens;
    size_t nsample;
    size_t ninit_param;
    size_t *save_param;
    size_t nsave_param;
    fit_model fm;
    rs_sample *s_a;
    mat *a_err;
    rs_sample *s_M_umd;
    rs_sample *s_ex;
    mat *ex_err;
} fit_param;

#define IS_AN(param,prog_flag) ((param)->analyze_flag & (prog_flag))

int ind_beta(const strbuf beta, const fit_param *param);
fit_param * fit_param_parse(const strbuf fname);
void fit_param_destroy(fit_param *param);

#endif
