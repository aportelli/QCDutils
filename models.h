#ifndef MODELS_H_
#define MODELS_H_

#include <latan/latan_fit.h>

#define SQ_NU_M_Kchi (0.5*(SQ(NU_M_K_p)+SQ(NU_M_K_0)-SQ(NU_M_pi_0_miso)))
#define NU_DMSQ_pi   (SQ(NU_M_pi_p_miso)-SQ(NU_M_pi_0_miso))

/* stages */
enum
{
    st_cst  = 1 << 0,
    st_pi   = 1 << 1,
    st_K    = 1 << 2,
    st_a    = 1 << 3,
    st_umd  = 1 << 4,
    st_V    = 1 << 5
};

/* prototypes */
extern const fit_model fm_phyptfit_ex_taylor;
extern const fit_model fm_scaleset_ex_taylor;
extern const fit_model fm_Zset_ex_taylor;

#endif