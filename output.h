#ifndef OUTPUT_H_
#define OUTPUT_H_

#include <latan/latan_statistics.h>
#include "parameters.h"

typedef enum
{
    Q     = 0,
    SCALE = 1
} plot_flag;

void plot_fit(const mat *fit, const mat *fit_var, fit_data *d,\
              fit_param *param, const plot_flag f);
void plot_chi2_comp(const fit_data *d, const fit_param *param, const size_t k,\
                    const strbuf title);
void print_result(const rs_sample *s_fit, fit_param *param);

#endif
