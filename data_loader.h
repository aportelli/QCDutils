#ifndef DATA_LOADER_H_
#define	DATA_LOADER_H_

#include "parameters.h"
#include <latan/latan_statistics.h>

/* indices */
#define N_EX_VAR 8
enum
{
    i_ud    = 0,
    i_s     = 1,
    i_bind  = 2,
    i_a     = 3,
    i_umd   = 4,
    i_alpha = 5,
    i_Linv  = 6,
    i_fvM   = 7
};

/* functions */
void data_load(rs_sample *s_x[N_EX_VAR], rs_sample *s_q[2], fit_param *param);

#endif
