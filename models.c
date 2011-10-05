#include "models.h"
#include "data_loader.h"
#include "parameters.h"
#include <math.h>
#include <string.h>
#include <latan/latan_nunits.h>
#include <latan/latan_mass.h>
#include <latan/latan_math.h>

#define ud_s_taylor(res,p,M_ud,M_s,param)\
{\
    int d_;\
    size_t i_;\
    if (param->M_ud_deg > 0)\
    {\
        res += mat_get(p,ST_pi_I(0,param),0)*M_ud;\
        for (d_=2;d_<=param->M_ud_deg;d_++)\
        {\
            i_    = (size_t)(d_) - 1;\
            res  += mat_get(p,ST_pi_I(i_,param),0)*pow(M_ud,(double)(d_));\
        }\
    }\
    if (param->M_s_deg > 0)\
    {\
        res += mat_get(p,ST_K_I(0,param),0)*M_s;\
        for (d_=2;d_<=param->M_s_deg;d_++)\
        {\
            i_   = (size_t)(d_) - 1;\
            res += mat_get(p,ST_K_I(i_,param),0)*pow(M_s,(double)(d_));\
        }\
    }\
}

#define ud_s_taylor_pstr(str,p,x_str,M_ud_ex,M_s_ex,param)\
{\
    int d_;\
    size_t i_;\
    strbuf buf_;\
    if (param->M_ud_deg > 0)\
    {\
        for (d_=1;d_<=param->M_ud_deg;d_++)\
        {\
            i_ = (size_t)(d_) - 1;\
            sprintf(buf_,"+%e*(%s-%e)**%d",mat_get(p,ST_pi_I(i_,param),0),\
                    x_str[i_ud],M_ud_ex,d_);\
            strcat(str,buf_);\
        }\
    }\
    if (param->M_s_deg > 0)\
    {\
        for (d_=1;d_<=param->M_s_deg;d_++)\
        {\
            i_ = (size_t)(d_) - 1;\
            sprintf(buf_,"+%e*(%s-%e)**%d",mat_get(p,ST_K_I(i_,param),0),\
                    x_str[i_s],M_s_ex,d_);\
            strcat(str,buf_);\
        }\
    }\
}

#define ST_pi_I(i,param)  (i+1)
#define ST_K_I(i,param)   (ST_pi_I(i,param)+(param)->M_ud_deg)
#define ST_a_I(i,param)   (ST_K_I(i,param)+(param)->M_s_deg)
#define ST_umd_I(i,param) (ST_a_I(i,param)+(((param)->a_deg > 0) ? 1 : 0))
#define ST_qedfv_I(i,param)  (ST_umd_I(i,param)+((param)->with_umd ? 1 : 0))

static double fm_phyptfit_taylor_func(const mat *X, const mat *p,\
                                         void *vparam)
{
    double res,M_ud,M_s,a,umd,Linv;
    ex_param *param;
    
    param = (ex_param *)vparam;
    res   = 0.0;
    a     = mat_get(X,i_ainv,0);
    M_ud  = mat_get(X,i_ud,0) - SQ(param->M_ud[0]);
    M_s   = mat_get(X,i_s,0) - SQ(param->M_s[0]);
    umd   = mat_get(X,i_umd,0);
    Linv  = mat_get(X,i_Linv,0); 
    
    res += mat_get(p,0,0);
    ud_s_taylor(res,p,M_ud,M_s,param);
    if (param->a_deg > 0)
    {
        res += mat_get(p,ST_a_I(0,param),0)*pow(a,(double)(param->a_deg));
    }
    if (param->with_umd)
    {
        res += mat_get(p,ST_umd_I(0,param),0)*umd;
    }
    if (param->with_qed_fvol > 0)
    {
        res += mat_get(p,ST_qedfv_I(0,param),0)*SQ(Linv);
    }
    
    return res;
}

static size_t fm_phyptfit_taylor_npar(void* vparam)
{
    ex_param *param;
    size_t npar;
    
    param = (ex_param *)vparam;
    
    npar  = 1;
    npar += param->M_ud_deg;
    npar += param->M_s_deg;
    if (param->a_deg > 0)
    {
        npar++;
    }
    if (param->with_umd)
    {
        npar++;
    }
    if (param->with_qed_fvol)
    {
        npar ++;
    }
    
    return npar;
}

static void fm_phyptfit_taylor_pstr(strbuf str, const size_t i,   \
                                       const mat *x_ex, const mat *p,\
                                       void *vparam)
{
    ex_param *param;
    strbuf buf,x_str[N_EX_VAR];
    double M_ud_phi,M_s_phi;
    size_t j;
    
    param    = (ex_param *)vparam;
    M_ud_phi = SQ(param->M_ud[0]);
    M_s_phi  = SQ(param->M_s[0]);
    
    for (j=0;j<N_EX_VAR;j++)
    {
        if (i == j)
        {
            strbufcpy(x_str[j],"x");
        }
        else
        {
            sprintf(x_str[j],"%e",mat_get(x_ex,j,0));
        }
    }
    sprintf(str,"%e+",mat_get(p,0,0));
    ud_s_taylor_pstr(str,p,x_str,M_ud_phi,M_s_phi,param);
    if (param->a_deg > 0)
    {
        sprintf(buf,"+%e*%s**%d",mat_get(p,ST_a_I(0,param),0),x_str[i_ainv],\
                param->a_deg);
        strcat(str,buf);
    }
    if (param->with_umd)
    {
        sprintf(buf,"+%e*%s",mat_get(p,ST_umd_I(0,param),0),x_str[i_umd]);
        strcat(str,buf);
    }
    if (param->with_qed_fvol > 0)
    {
        sprintf(buf,"+%e*%s**2",mat_get(p,ST_qedfv_I(0,param),0),x_str[i_Linv]);
        strcat(str,buf);
    }
}

const fit_model fm_phyptfit_taylor =
{
    "physical point Taylor expansion",
    &fm_phyptfit_taylor_func,
    &fm_phyptfit_taylor_npar,
    &fm_phyptfit_taylor_pstr,
    N_EX_VAR
};

#undef ST_pi_I
#define ST_pi_I(i,param) (i+(param)->nbeta)

static double fm_scaleset_taylor_func(const mat *X, const mat *p,\
                                         void *vparam)
{
    double res,M_ud,M_s,umd,M_scale,Linv;
    ex_param *param;
    size_t bind;
    
    param   = (ex_param *)vparam;
    res     = 0.0;
    bind    = (size_t)(mat_get(X,i_bind,0));
    M_scale = param->M_scale[0];
    Linv    = mat_get(X,i_Linv,0);
    M_ud    = mat_get(X,i_ud,0)/mat_get(p,bind,0)\
                - SQ(param->M_ud[0])/SQ(M_scale);
    M_s     = mat_get(X,i_s,0)/mat_get(p,bind,0) \
                - SQ(param->M_s[0])/SQ(M_scale);
    umd     = mat_get(X,i_umd,0);
    
    res += 1.0;
    ud_s_taylor(res,p,M_ud,M_s,param);
    if (param->with_umd)
    {
        res += mat_get(p,ST_umd_I(0,param),0)*umd;
    }
    if (param->with_qed_fvol > 0)
    {
        res += mat_get(p,ST_qedfv_I(0,param),0)*SQ(Linv)/mat_get(p,bind,0);
        
    }
    res *= mat_get(p,bind,0);
    
    return res;
}

static size_t fm_scaleset_taylor_npar(void* vparam)
{
    ex_param *param;
    size_t npar;
    
    param = (ex_param *)vparam;
    
    npar  = param->nbeta;
    npar += param->M_ud_deg;
    npar += param->M_s_deg;
    if (param->with_umd)
    {
        npar++;
    }
    if (param->with_qed_fvol)
    {
        npar += 1;
    }
    
    return npar;
}

static void fm_scaleset_taylor_pstr(strbuf str, const size_t i,   \
                                       const mat *x_ex, const mat *p,\
                                       void *vparam)
{
    ex_param *param;
    strbuf buf,x_str[N_EX_VAR],fac_str;
    double M_ud_phi,M_s_phi,M_scale;
    size_t bind;
    size_t j;
    
    param    = (ex_param *)vparam;
    bind     = (size_t)(mat_get(x_ex,i_bind,0));
    M_scale  = param->M_scale[0];
    M_ud_phi = SQ(param->M_ud[0])/SQ(M_scale);
    M_s_phi  = SQ(param->M_s[0])/SQ(M_scale);
    
    for (j=0;j<N_EX_VAR;j++)
    {
        if (i == j)
        {
            sprintf(x_str[j],"(x/%e)",mat_get(p,bind,0));
        }
        else
        {
            sprintf(x_str[j],"%e",mat_get(x_ex,j,0)/mat_get(p,bind,0));
        }
    }
    sprintf(fac_str,"%e",mat_get(p,bind,0));
    sprintf(str,"%s*(1.0",fac_str);
    ud_s_taylor_pstr(str,p,x_str,M_ud_phi,M_s_phi,param);
    if (param->with_umd)
    {
        sprintf(buf,"+%e*%s",mat_get(p,ST_umd_I(0,param),0),x_str[i_umd]);
        strcat(str,buf);
    }
    if (param->with_qed_fvol > 0)
    {
        sprintf(buf,"+%e*%s**2/%e",mat_get(p,ST_qedfv_I(0,param),0),\
                x_str[i_Linv],mat_get(p,bind,0));
        
    }
    strcat(str,")");
}

const fit_model fm_scaleset_taylor =
{
    "physical point Taylor expansion",
    &fm_scaleset_taylor_func,
    &fm_scaleset_taylor_npar,
    &fm_scaleset_taylor_pstr,
    N_EX_VAR
};
