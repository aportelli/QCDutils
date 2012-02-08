#ifndef OUTPUT_H_
#define OUTPUT_H_

#include <latan/latan_statistics.h>
#include "parameters.h"

typedef enum
{
    Q     = 0,
    SCALE = 1
} plot_flag;

void plot_fit(const mat *fit, const mat *fit_var, fit_data *d, \
              fit_param *param, plot_flag f);
void print_result(const rs_sample *s_fit, fit_param *param);

#endif
