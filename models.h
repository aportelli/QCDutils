#ifndef MODELS_H_
#define MODELS_H_

#include <latan/latan_fit.h>
#include <latan/latan_nunits.h>
#include <latan/latan_statistics.h>

/* chi^2 extension for lattice spacing errors */
double a_error_chi2_ext(const mat *p, void *vd);

/* model functions */
/** physical point Taylor expansion **/
double fm_phypt_taylor_func(const mat *X, const mat *p, void *vparam);
double fm_scaleset_taylor_func(const mat *X, const mat *p, void *vparam);

/* model number of parameter */
/** physical point Taylor expansion **/
size_t fm_phypt_taylor_npar(void* vparam);
size_t fm_scaleset_taylor_npar(void* vparam);
/** combined models **/
size_t fm_comb_phypt_taylor_scaleset_taylor_npar(void* vparam);

#endif
