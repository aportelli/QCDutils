#ifndef DATA_LOADER_H_
#define	DATA_LOADER_H_

#include "parameters.h"
#include <latan/latan_statistics.h>

/* indices */
#define N_EX_VAR 10
enum
{
    i_ud    = 0,
    i_s     = 1,
    i_bind  = 2,
    i_vind  = 3,
    i_a     = 4,
    i_umd   = 5,
    i_alpha = 6,
    i_Linv  = 7,
    i_ToL   = 8,
    i_fvM   = 9
};

/* functions */
void data_load(rs_sample *s_x[N_EX_VAR], rs_sample *s_q[2], fit_param *param);

#endif
