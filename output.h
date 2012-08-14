#ifndef OUTPUT_H_
#define OUTPUT_H_

#include <stdio.h>
#include <latan/latan_statistics.h>
#include "parameters.h"
#include "data_loader.h"

typedef enum
{
    Q     = 0,
    SCALE = 1
} plot_flag;

void plot_fit(const mat *fit, fit_data *d, fit_param *param, const plot_flag f);
void plot_chi2_comp(const fit_data *d, const fit_param *param, const size_t k,\
                    const strbuf title);
void print_result(const rs_sample *s_fit, fit_param *param);
void fprint_table(FILE* stream, rs_sample *s_x[N_EX_VAR], rs_sample *s_q[2],\
                  const mat *res, const fit_param *param, const plot_flag f);

#endif
