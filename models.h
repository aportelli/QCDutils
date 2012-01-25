#ifndef MODELS_H_
#define MODELS_H_

#include <latan/latan_fit.h>
#include <latan/latan_nunits.h>

#define DMSQ_K (SQ(NU_M_K_p)-SQ(NU_M_K_0))

/* prototypes */
extern const fit_model fm_phypt_a_taylor;
extern const fit_model fm_scaleset_taylor;
extern const fit_model fm_comb_phyt_scale_taylor;

#endif