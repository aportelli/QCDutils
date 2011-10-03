#include "models.h"
#include "data_loader.h"
#include "parameters.h"
#include <math.h>
#include <string.h>
#include <latan/latan_nunits.h>
#include <latan/latan_mass.h>
#include <latan/latan_math.h>

#define ST_pi_I(i,param)  (i+1)
#define ST_K_I(i,param)   (ST_pi_I(i,param)+(param)->Msq_pi_deg)
#define ST_a_I(i,param)   (ST_K_I(i,param)+(param)->Msq_K_deg)
#define ST_umd_I(i,param) (ST_a_I(i,param)+(param)->a_deg)
#define ST_fv_I(i,param)  (ST_umd_I(i,param)+((param)->with_umd ? 1 : 0))

static double fm_phyptfit_ex_taylor_func(const mat *X, const mat *p,\
                                         void *vparam)
{
    double res,Msq_pi,Msq_K,a,umd,Linv;
    ex_param *param;
    int d;
    size_t i;
    
    param   = (ex_param *)vparam;
    res     = 0.0;
    a       = mat_get(X,i_ainv,0);
    Msq_pi  = mat_get(X,i_pi,0) - SQ(NU_M_pi_p_miso);
    Msq_K   = mat_get(X,i_K,0) - SQ_NU_M_Kchi;
    umd     = mat_get(X,i_umd,0);
    Linv    = mat_get(X,i_Linv,0); 
    
    res += mat_get(p,0,0);
    if (param->Msq_pi_deg > 0)
    {
        res += mat_get(p,ST_pi_I(0,param),0)*Msq_pi;
        for (d=2;d<=param->Msq_pi_deg;d++)
        {
            i     = (size_t)(d) - 1;
            res  += mat_get(p,ST_pi_I(i,param),0)*pow(Msq_pi,(double)(d));
        }
    }
    if (param->Msq_K_deg > 0)
    {
        res += mat_get(p,ST_K_I(0,param),0)*Msq_K;
        for (d=2;d<=param->Msq_K_deg;d++)
        {
            i    = (size_t)(d) - 1;
            res += mat_get(p,ST_K_I(i,param),0)*pow(Msq_K,(double)(d));
        }
    }
    if (param->a_deg > 0)
    {
        res += mat_get(p,ST_a_I(0,param),0)*a;
        for (d=2;d<=param->a_deg;d++)
        {
            i    = (size_t)(d) - 1;
            res += mat_get(p,ST_a_I(i,param),0)*pow(a,(double)(d));
        }
    }
    if (param->with_umd)
    {
        res += mat_get(p,ST_umd_I(0,param),0)*umd;
    }
    if (param->with_fvol > 0)
    {
        res += mat_get(p,ST_fv_I(0,param),0)*SQ(Linv);
    }
    
    return res;
}

static size_t fm_phyptfit_ex_taylor_npar(void* vparam)
{
    ex_param *param;
    size_t npar;
    
    param = (ex_param *)vparam;
    
    npar = 1;
    npar += param->Msq_pi_deg;
    npar += param->Msq_K_deg;
    npar += param->a_deg;
    if (param->with_umd)
    {
        npar++;
    }
    if (param->with_fvol)
    {
        npar += 1;
    }
    
    return npar;
}

static void fm_phyptfit_ex_taylor_pstr(strbuf str, const size_t i,   \
                                       const mat *x_ex, const mat *p,\
                                       void *vparam)
{
    ex_param *param;
    strbuf buf,x_str[N_EX_VAR];
    double Msq_pi_phi,Msq_K_phi;
    size_t j;
    int d;
    
    param      = (ex_param *)vparam;
    Msq_pi_phi = SQ(NU_M_pi_p);
    Msq_K_phi  = SQ_NU_M_Kchi;
    
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
    sprintf(str,"%e",mat_get(p,0,0));
    if (param->Msq_pi_deg > 0)
    {
        for (d=1;d<=param->Msq_pi_deg;d++)
        {
            j = (size_t)(d) - 1;
            sprintf(buf,"+%e*(%s-%e)**%d",mat_get(p,ST_pi_I(j,param),0),\
                    x_str[i_pi],Msq_pi_phi,d);
            strcat(str,buf);
        }
    }
    if (param->Msq_K_deg > 0)
    {
        for (d=1;d<=param->Msq_K_deg;d++)
        {
            j = (size_t)(d) - 1;
            sprintf(buf,"+%e*(%s-%e)**%d",mat_get(p,ST_K_I(j,param),0),\
                    x_str[i_K],Msq_K_phi,d);
            strcat(str,buf);
        }
    }
    if (param->a_deg > 0)
    {
        for (d=1;d<=param->a_deg;d++)
        {
            j = (size_t)(d) - 1;
            sprintf(buf,"+%e*%s**%d",mat_get(p,ST_a_I(j,param),0),\
                    x_str[i_ainv],d);
            strcat(str,buf);
        }
    }
    if (param->with_umd)
    {
        sprintf(buf,"+%e*%s",mat_get(p,ST_umd_I(0,param),0),x_str[i_umd]);
        strcat(str,buf);
    }
    if (param->with_fvol > 0)
    {
        sprintf(buf,"+%e*%s**2",mat_get(p,ST_fv_I(0,param),0),x_str[i_Linv]);
        strcat(str,buf);
    }
}

const fit_model fm_phyptfit_ex_taylor =
{
    "physical point Taylor expansion",
    &fm_phyptfit_ex_taylor_func,
    &fm_phyptfit_ex_taylor_npar,
    &fm_phyptfit_ex_taylor_pstr,
    N_EX_VAR
};

#undef ST_pi_I
#define ST_pi_I(i,param) (i+(param)->nbeta)

static double fm_scaleset_ex_taylor_func(const mat *X, const mat *p,\
                                         void *vparam)
{
    double res,Msq_pi,Msq_K,umd,M_scale,Linv;
    ex_param *param;
    int d;
    size_t i,bind;
    
    param   = (ex_param *)vparam;
    res     = 0.0;
    bind    = (size_t)(mat_get(X,i_bind,0));
    M_scale = param->M_scale[0];
    Linv    = mat_get(X,i_Linv,0);
    Msq_pi  = mat_get(X,i_pi,0)/mat_get(p,bind,0) - SQ(NU_M_pi_p)/SQ(M_scale);
    Msq_K   = mat_get(X,i_K,0)/mat_get(p,bind,0) - SQ_NU_M_Kchi/SQ(M_scale);
    umd     = mat_get(X,i_umd,0);
    
    res += 1.0;
    if (param->Msq_pi_deg > 0)
    {
        res += mat_get(p,ST_pi_I(0,param),0)*Msq_pi;
        for (d=2;d<=param->Msq_pi_deg;d++)
        {
            i     = (size_t)(d) - 1;
            res  += mat_get(p,ST_pi_I(i,param),0)*pow(Msq_pi,(double)(d));
        }
    }
    if (param->Msq_K_deg > 0)
    {
        res += mat_get(p,ST_K_I(0,param),0)*Msq_K;
        for (d=2;d<=param->Msq_K_deg;d++)
        {
            i    = (size_t)(d) - 1;
            res += mat_get(p,ST_K_I(i,param),0)*pow(Msq_K,(double)(d));
        }
    }
    if (param->with_umd)
    {
        res += mat_get(p,ST_umd_I(0,param),0)*umd;
    }
    if (param->with_fvol > 0)
    {
        res += mat_get(p,ST_fv_I(0,param),0)*SQ(Linv)/mat_get(p,bind,0);
        
    }
    res *= mat_get(p,bind,0);
    
    return res;
}

static size_t fm_scaleset_ex_taylor_npar(void* vparam)
{
    ex_param *param;
    size_t npar;
    
    param = (ex_param *)vparam;
    
    npar  = param->nbeta;
    npar += param->Msq_pi_deg;
    npar += param->Msq_K_deg;
    if (param->with_umd)
    {
        npar++;
    }
    if (param->with_fvol)
    {
        npar += 1;
    }
    
    return npar;
}

static void fm_scaleset_ex_taylor_pstr(strbuf str, const size_t i,   \
                                       const mat *x_ex, const mat *p,\
                                       void *vparam)
{
    ex_param *param;
    strbuf buf,x_str[N_EX_VAR],fac_str;
    double Msq_pi_phi,Msq_K_phi,M_scale;
    size_t bind;
    size_t j;
    int d;
    
    param      = (ex_param *)vparam;
    bind       = (size_t)(mat_get(x_ex,i_bind,0));
    M_scale    = param->M_scale[0];
    Msq_pi_phi = SQ(NU_M_pi_p)/SQ(M_scale);
    Msq_K_phi  = SQ_NU_M_Kchi/SQ(M_scale);
    
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
    if (param->Msq_pi_deg > 0)
    {
        for (d=1;d<=param->Msq_pi_deg;d++)
        {
            j = (size_t)(d) - 1;
            sprintf(buf,"+%e*(%s-%e)**%d",mat_get(p,ST_pi_I(j,param),0),\
                    x_str[i_pi],Msq_pi_phi,d);
            strcat(str,buf);
        }
    }
    if (param->Msq_K_deg > 0)
    {
        for (d=1;d<=param->Msq_K_deg;d++)
        {
            j = (size_t)(d) - 1;
            sprintf(buf,"+%e*(%s-%e)**%d",mat_get(p,ST_K_I(j,param),0),\
                    x_str[i_K],Msq_K_phi,d);
            strcat(str,buf);
        }
    }
    if (param->a_deg > 0)
    {
        for (d=1;d<=param->a_deg;d++)
        {
            j = (size_t)(d) - 1;
            sprintf(buf,"+%e*%s**%d",mat_get(p,ST_a_I(j,param),0),\
                    x_str[i_ainv],d);
            strcat(str,buf);
        }
    }
    if (param->with_umd)
    {
        sprintf(buf,"+%e*%s",mat_get(p,ST_umd_I(0,param),0),x_str[i_umd]);
        strcat(str,buf);
    }
    if (param->with_fvol > 0)
    {
        sprintf(buf,"+%e*%s**2/%e",mat_get(p,ST_fv_I(0,param),0),x_str[i_Linv],\
                mat_get(p,bind,0));
        
    }
    strcat(str,")");
}

const fit_model fm_scaleset_ex_taylor =
{
    "physical point Taylor expansion",
    &fm_scaleset_ex_taylor_func,
    &fm_scaleset_ex_taylor_npar,
    &fm_scaleset_ex_taylor_pstr,
    N_EX_VAR
};
