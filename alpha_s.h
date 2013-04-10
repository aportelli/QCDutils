#ifndef QED_FVOL_TABFUNC_H_
#define QED_FVOL_TABFUNC_H_

#define LAMBDA_MSBAR_2500MEV 335.869

double alpha_s_msbar(const double mu, const double Lambda,          \
                     const unsigned int nf, const unsigned int nloop);

#endif
