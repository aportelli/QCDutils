#ifndef DATA_LOADER_H_
#define	DATA_LOADER_H_

#include "parameters.h"
#include <latan/latan_statistics.h>

/* indices */
#define N_EX_VAR 7
enum
{
    i_pi   = 0,
    i_K    = 1,
    i_bind = 2,
    i_ainv = 3,
    i_umd  = 4,
    i_Linv = 5,
    i_MpiL = 6
};

/* functions */
double get_cst_ainv(double beta);
double get_e0(double beta);
void data_load(rs_sample *s_x[N_EX_VAR], rs_sample *s_q, strbuf beta,\
               const ex_param *param);

#endif

