#include "alpha_s.h"
#include <math.h>
#include <latan/latan_math.h>

#define ZETA_3 1.2020569031595942854
#define MAX_MSBAR_LOOP 4

double alpha_s_msbar(const double mu, const double Lambda,          \
                     const unsigned int nf, const unsigned int nloop)
{
    double res,logml,loglogml,beta[MAX_MSBAR_LOOP];
    
    if (nloop > MAX_MSBAR_LOOP)
    {
        fprintf(stderr,"error: MSbar alpha_s is not known at %d loops\n",nloop);
        abort();
    }
    
    if (latan_isinf(mu))
    {
        res = 0.0;
    }
    else
    {
        logml    = log(SQ(mu/Lambda));
        loglogml = log(logml);
        beta[0]  = (nloop > 0) ? 11.0 - 2.0/3.0*nf                    : 0.0;
        beta[1]  = (nloop > 1) ? 102.0 - 38.0/3.0*nf                  : 0.0;
        beta[2]  = (nloop > 2) ? 2857.0/2.0-5033.0/18.0*nf                 \
                                 +325.0/54.0*SQ(nf)                   : 0.0;
        beta[3]  = (nloop > 3) ? 149753.0/6.0+3564.0*ZETA_3                \
                                 -(1078361.0/162.0+6508.0/27.0*ZETA_3)*nf  \
                                 +(50065.0/162.0+6472.0/81.0*ZETA_3)*SQ(nf)\
                                 +1093.0/729.0*nf*nf*nf               : 0.0;
        res  = 1.0/(beta[0]*logml);
        res -= beta[1]/(pow(beta[0],3.0)*pow(logml,2.0))*loglogml;
        res += pow(beta[1],2.0)/(pow(beta[0],5.0)*pow(logml,3.0))                 \
               *(pow(loglogml,2.0)-loglogml-1.0+beta[0]*beta[2]/pow(beta[1],2.0));
        res += pow(beta[1],3.0)/(pow(beta[0],7.0)*pow(logml,4.0))                 \
               *(-pow(loglogml,3.0)+5.0/2.0*pow(loglogml,2.0)+2.0*loglogml-0.5    \
                 -3.0*beta[0]*beta[2]/pow(beta[1],2.0)*loglogml                   \
                 +pow(beta[0],2.0)*beta[3]/(2.0*pow(beta[1],3.0)));
        res *= 4.0*C_PI;
    }
    
    return res;
}
