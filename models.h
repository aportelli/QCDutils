#ifndef MODELS_H_
#define MODELS_H_

#include <latan/latan_fit.h>
#include <latan/latan_nunits.h>

/* prototypes */
double a_error_chi2_ext(mat *p, void *vd);
double fm_phypt_a_taylor_func(const mat *X, const mat *p, void *vparam);
size_t fm_phypt_a_taylor_npar(void* vparam);
void   fm_phypt_a_taylor_pstr(strbuf str, const size_t i,   \
                              const mat *x_ex, const mat *p,\
                              void *vparam);
double fm_scaleset_taylor_func(const mat *X, const mat *p, void *vparam);
size_t fm_scaleset_taylor_npar(void* vparam);
void   fm_scaleset_taylor_pstr(strbuf str, const size_t i,   \
                               const mat *x_ex, const mat *p,\
                               void *vparam);
size_t fm_comb_phyt_scale_npar(void* vparam);

extern const fit_model fm_phypt_a_taylor;
extern const fit_model fm_scaleset_taylor;
extern const fit_model fm_comb_phyt_scale_taylor;

#endif
