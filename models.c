#include "models.h"
#include "data_loader.h"
#include "parameters.h"
#include <math.h>
#include <string.h>
#include <latan/latan_nunits.h>
#include <latan/latan_mass.h>
#include <latan/latan_math.h>

/* physical point Taylor expansion models */
static double fm_phypt_a_taylor_func(const mat *X, const mat *p, void *vparam);
static size_t fm_phypt_a_taylor_npar(void* vparam);
static void   fm_phypt_a_taylor_pstr(strbuf str, const size_t i,   \
                                     const mat *x_ex, const mat *p,\
                                     void *vparam);
static double fm_scaleset_taylor_func(const mat *X, const mat *p, void *vparam);
static size_t fm_scaleset_taylor_npar(void* vparam);
static void   fm_scaleset_taylor_pstr(strbuf str, const size_t i,   \
                                      const mat *x_ex, const mat *p,\
                                      void *vparam);
static size_t fm_comb_phyt_scale_npar(void* vparam);

#define ud_s_taylor(res,p,s,M_ud,M_ud_deg,M_s,M_s_deg,param)\
{\
    int d_;\
    size_t i_;\
    if ((M_ud_deg) > 0)\
    {\
        res += mat_get(p,ST_pi_I(0,param)+(s),0)*M_ud;\
        for (d_=2;d_<=(M_ud_deg);d_++)\
        {\
            i_    = (size_t)(d_) - 1;\
            res  += mat_get(p,ST_pi_I(i_,param)+(s),0)*pow(M_ud,(double)(d_));\
        }\
    }\
    if ((M_s_deg) > 0)\
    {\
        res += mat_get(p,ST_K_I(0,param)+(s),0)*M_s;\
        for (d_=2;d_<=(M_s_deg);d_++)\
        {\
            i_   = (size_t)(d_) - 1;\
            res += mat_get(p,ST_K_I(i_,param)+(s),0)*pow(M_s,(double)(d_));\
        }\
    }\
}

#define ud_s_taylor_pstr(str,p,s,x_str,M_ud_ex,M_ud_deg,M_s_ex,M_s_deg,param)\
{\
    int d_;\
    size_t i_;\
    strbuf buf_;\
    if ((M_ud_deg) > 0)\
    {\
        for (d_=1;d_<=(M_ud_deg);d_++)\
        {\
            i_ = (size_t)(d_) - 1;\
            sprintf(buf_,"+%e*(%s-%e)**%d",mat_get(p,ST_pi_I(i_,param)+(s),0),\
                    x_str[i_ud],M_ud_ex,d_);\
            strcat(str,buf_);\
        }\
    }\
    if ((M_s_deg) > 0)\
    {\
        for (d_=1;d_<=(M_s_deg);d_++)\
        {\
            i_ = (size_t)(d_) - 1;\
            sprintf(buf_,"+%e*(%s-%e)**%d",mat_get(p,ST_K_I(i_,param)+(s),0),\
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

static double fm_phypt_a_taylor_func(const mat *X, const mat *p, void *vparam)
{
    double res,M_ud,M_s,a,dimfac,umd,Linv;
    size_t s,bind;
    fit_param *param;
    
    param = (fit_param *)vparam;
    res   = 0.0;
    bind  = (size_t)(mat_get(X,i_bind,0));
    if (IS_ANALYZE(param,"phypt"))
    {
        a = mat_get(X,i_a,0);
        s = 0;
    }
    else if (IS_ANALYZE(param,"comb_phypt_scale"))
    {
        a = mat_get(p,bind,0);
        s = fm_scaleset_taylor_npar(param);
    }
    else
    {
        fprintf(stderr,"error: this model should not ne used in program %s\n",\
                param->analyze);
        exit(EXIT_FAILURE);
    }
    dimfac = (!param->plotting) ? a : 1.0;
    M_ud   = mat_get(X,i_ud,0)/SQ(dimfac) - SQ(param->M_ud);
    M_s    = mat_get(X,i_s,0)/SQ(dimfac) - SQ(param->M_s);
    umd    = mat_get(X,i_umd,0)/SQ(dimfac) - DMSQ_K;
    Linv   = mat_get(X,i_Linv,0)/dimfac; 
    
    res += mat_get(p,s,0);
    ud_s_taylor(res,p,s,M_ud,param->M_ud_deg,M_s,param->M_s_deg,param);
    if (param->a_deg > 0)
    {
        res += mat_get(p,ST_a_I(0,param)+s,0)*pow(a,(double)(param->a_deg));
    }
    if (param->with_umd)
    {
        res += mat_get(p,ST_umd_I(0,param)+s,0)*umd;
    }
    if (param->with_qed_fvol)
    {
        res += mat_get(p,ST_qedfv_I(0,param)+s,0)*SQ(Linv);
    }
    res *= pow(dimfac,param->q_dim);
    
    return res;
}

static size_t fm_phypt_a_taylor_npar(void* vparam)
{
    fit_param *param;
    size_t npar;
    
    param = (fit_param *)vparam;
    
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
        npar++;
    }
    
    return npar;
}

static void fm_phypt_a_taylor_pstr(strbuf str, const size_t i,   \
                                   const mat *x_ex, const mat *p,\
                                   void *vparam)
{
    fit_param *param;
    strbuf buf,x_str[N_EX_VAR];
    double M_ud_phi,M_s_phi;
    size_t s;
    size_t j;
    
    param    = (fit_param *)vparam;
    if (IS_ANALYZE(param,"phypt"))
    {
        s = 0;
    }
    else if (IS_ANALYZE(param,"comb_phypt_scale"))
    {
        s = fm_scaleset_taylor_npar(param);
    }
    M_ud_phi = SQ(param->M_ud);
    M_s_phi  = SQ(param->M_s);
    
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
    sprintf(str,"%e+",mat_get(p,s,0));
    ud_s_taylor_pstr(str,p,s,x_str,M_ud_phi,param->M_ud_deg,M_s_phi,\
                     param->M_s_deg,param);
    if (param->a_deg > 0)
    {
        sprintf(buf,"+%e*%s**%d",mat_get(p,ST_a_I(0,param)+s,0),x_str[i_a],\
                param->a_deg);
        strcat(str,buf);
    }
    if (param->with_umd)
    {
        sprintf(buf,"+%e*(%s-%e)",mat_get(p,ST_umd_I(0,param)+s,0),\
                x_str[i_umd],DMSQ_K);
        strcat(str,buf);
    }
    if (param->with_qed_fvol > 0)
    {
        sprintf(buf,"+%e*%s**2",mat_get(p,ST_qedfv_I(0,param)+s,0),\
                x_str[i_Linv]);
        strcat(str,buf);
    }
}

const fit_model fm_phypt_a_taylor =
{
    "physical point Taylor expansion",
    {&fm_phypt_a_taylor_func},
    &fm_phypt_a_taylor_npar,
    {&fm_phypt_a_taylor_pstr},
    N_EX_VAR,
    1
};

#undef ST_pi_I
#undef ST_K_I
#undef ST_a_I
#undef ST_umd_I
#undef ST_qedfv_I
#define ST_pi_I(i,param)  (i+(param)->nbeta)
#define ST_K_I(i,param)   (ST_pi_I(i,param)+(param)->s_M_ud_deg)
#define ST_umd_I(i,param)   (ST_K_I(i,param)+(param)->s_M_s_deg)
#define ST_qedfv_I(i,param)  (ST_umd_I(i,param)+((param)->s_with_umd ? 1 : 0))

static double fm_scaleset_taylor_func(const mat *X, const mat *p, void *vparam)
{
    double res,a,M_ud,M_s,umd,M_scale,Linv;
    fit_param *param;
    size_t bind;
    
    param   = (fit_param *)vparam;
    res     = 0.0;
    bind    = (size_t)(mat_get(X,i_bind,0));
    a       = mat_get(p,bind,0);
    M_scale = param->M_scale;
    Linv    = mat_get(X,i_Linv,0);
    M_ud    = mat_get(X,i_ud,0)/SQ(a*M_scale)-SQ(param->M_ud)/SQ(M_scale);
    M_s     = mat_get(X,i_s,0)/SQ(a*M_scale)-SQ(param->M_s)/SQ(M_scale);
    umd     = mat_get(X,i_umd,0)/SQ(a*M_scale)-DMSQ_K/SQ(M_scale);
    
    res += 1.0;
    ud_s_taylor(res,p,0,M_ud,param->s_M_ud_deg,M_s,param->s_M_ud_deg,param);
    if (param->s_with_umd)
    {
        res += mat_get(p,ST_umd_I(0,param),0)*umd;
    }
    if (param->s_with_qed_fvol > 0)
    {
        res += mat_get(p,ST_qedfv_I(0,param),0)*SQ(Linv)/mat_get(p,bind,0);
        
    }
    res *= SQ(a*M_scale);
    
    return res;
}

static size_t fm_scaleset_taylor_npar(void* vparam)
{
    fit_param *param;
    size_t npar;
    
    param = (fit_param *)vparam;
    
    npar  = param->nbeta;
    npar += param->s_M_ud_deg;
    npar += param->s_M_s_deg;
    if (param->s_with_umd)
    {
        npar++;
    }
    if (param->s_with_qed_fvol)
    {
        npar++;
    }
    
    return npar;
}

static void fm_scaleset_taylor_pstr(strbuf str, const size_t i,   \
                                    const mat *x_ex, const mat *p,\
                                    void *vparam)
{
    fit_param *param;
    strbuf buf,x_str[N_EX_VAR],fac_str;
    double M_ud_phi,M_s_phi,M_scale,a;
    size_t bind;
    size_t j;
    
    param    = (fit_param *)vparam;
    bind     = (size_t)(mat_get(x_ex,i_bind,0));
    a        = mat_get(p,bind,0);
    M_scale  = param->M_scale;
    M_ud_phi = SQ(param->M_ud)/SQ(M_scale);
    M_s_phi  = SQ(param->M_s)/SQ(M_scale);
    
    for (j=0;j<N_EX_VAR;j++)
    {
        if (i == j)
        {
            sprintf(x_str[j],"(x/%e)",SQ(a*M_scale));
        }
        else
        {
            sprintf(x_str[j],"%e",mat_get(x_ex,j,0)/SQ(a*M_scale));
        }
    }
    sprintf(fac_str,"%e",mat_get(p,bind,0));
    sprintf(str,"%s*(1.0",fac_str);
    ud_s_taylor_pstr(str,p,0,x_str,M_ud_phi,param->s_M_ud_deg,M_s_phi,\
                     param->s_M_s_deg,param);
    if (param->s_with_umd)
    {
        sprintf(buf,"+%e*(%s-%e)",mat_get(p,ST_umd_I(0,param),0),x_str[i_umd],\
                DMSQ_K);
        strcat(str,buf);
    }
    if (param->s_with_qed_fvol > 0)
    {
        sprintf(buf,"+%e*%s**2/%e",mat_get(p,ST_qedfv_I(0,param),0),\
                x_str[i_Linv],mat_get(p,bind,0));
        
    }
    strcat(str,")");
}

const fit_model fm_scaleset_taylor =
{
    "physical point Taylor expansion",
    {&fm_scaleset_taylor_func},
    &fm_scaleset_taylor_npar,
    {&fm_scaleset_taylor_pstr},
    N_EX_VAR,
    1
};

static size_t fm_comb_phyt_scale_npar(void* vparam)
{
    return fm_phypt_a_taylor_npar(vparam)+fm_scaleset_taylor_npar(vparam);
}

const fit_model fm_comb_phyt_scale_taylor =
{
    "physical point Taylor expansion",
    {&fm_scaleset_taylor_func,&fm_phypt_a_taylor_func},
    &fm_comb_phyt_scale_npar,
    {&fm_scaleset_taylor_pstr,&fm_phypt_a_taylor_pstr},
    N_EX_VAR,
    2
};
