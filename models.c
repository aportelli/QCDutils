#include "models.h"
#include "data_loader.h"
#include "parameters.h"
#include "qed_fvol_tabfunc.h"
#include <math.h>
#include <string.h>
#include <latan/latan_nunits.h>
#include <latan/latan_mass.h>
#include <latan/latan_math.h>
#include <latan/latan_tabfunc.h>

/* polynom macro */
#define polynom(res,p,i0,x,deg)\
if ((deg) > 0)\
{\
    int d_;\
    size_t i_;\
    res += mat_get(p,i0,0)*(x);\
    for (d_=2;d_<=(deg);d_++)\
    {\
        i_    = (size_t)(d_) - 1;\
        res  += mat_get(p,(i0)+i_,0)*pow(x,(double)(d_));\
    }\
}\

#define I_ud(i)         (i+1)
#define I_s(i)          (I_ud(i)+(param)->M_ud_deg)
#define I_discr(i)      (I_s(i)+(param)->M_s_deg)
#define I_dis_a2M_ud(i) (I_discr(i)+(param)->a_deg)
#define I_dis_a2M_s(i)  (I_dis_a2M_ud(i)+((param)->with_a2M_ud ? 1 : 0))
#define I_umd(i)        (I_dis_a2M_s(i)+((param)->with_a2M_s ? 1 : 0))
#define I_udumd(i)      (I_umd(i)+(param)->umd_deg) 
#define I_sumd(i)       (I_udumd(i)+((param)->with_udumd ? 1 : 0))
#define I_alpha(i)      (I_sumd(i)+((param)->with_sumd ? 1 : 0))
#define I_udalpha(i)    (I_alpha(i)+(param)->alpha_deg)
#define I_salpha(i)     (I_udalpha(i)+((param)->with_udalpha ? 1 : 0))
#define I_qedfv(i)      (I_salpha(i)+((param)->with_salpha ? 1 : 0))
#define I_a(i)          (I_qedfv(i)+param->with_qed_fvol)

double a_error_chi2_ext(const mat *p, void *vd)
{
    size_t i,s;
    double a,a_err,res;
    fit_param *param;
    fit_data *d;
    
    res   = 0.0;
    d     = (fit_data *)(vd);
    param = (fit_param *)(d->model_param); 
    s     = fit_data_get_sample_counter(d);
    
    for (i=0;i<param->nbeta;i++)
    {
        if (s == 0)
        {
            a = mat_get(rs_sample_pt_cent_val(param->a),i,0);
        }
        else
        {
            a = mat_get(rs_sample_pt_sample(param->a,s-1),i,0);
        }
        a_err = mat_get(param->a_err,i,0);
        res  += SQ(mat_get(p,I_a(i),0) - a)/SQ(a_err);
        d->matperf += 5.0;
    }
    
    return res;
}

double fm_phypt_taylor_func(const mat *X, const mat *p, void *vparam)
{
    double res,M_ud,M_s,a,dimfac,umd,Linv,a2mud,a2ms,alpha,dalpha;
    size_t s,bind;
    fit_param *param;
    
    param = (fit_param *)vparam;
    res   = 0.0;
    bind  = (size_t)(mat_get(X,i_bind,0));
    
    /* what is the lattice spacing ? */
    if (IS_AN(param,AN_PHYPT))
    {
        if (IS_AN(param,AN_SCALE))
        {
            /** a fit parameter from the scale setting model **/
            a = (!param->plotting) ? mat_get(p,bind,0) : mat_get(X,i_a,0);
            s = fm_scaleset_taylor_npar(param);
        }
        else
        {
            /** a global parameter with error (with external scale samples) **/
            if (param->with_ext_a)
            {
                a = (!param->plotting) ? mat_get(p,I_a(bind),0) :\
                                         mat_get(X,i_a,0);
            }
            /** a x coordinate (with the ratio method) **/
            else
            {
                a = mat_get(X,i_a,0);
            }
            s = 0;
        }
    }
    else
    {
        fprintf(stderr,"error: this model should not ne used in program %s\n",\
                param->analyze);
        exit(EXIT_FAILURE);
    }
    /* x values */
    dimfac = (!param->plotting) ? a : 1.0;
    M_ud   = mat_get(X,i_ud,0)/SQ(dimfac) - SQ(param->M_ud);
    M_s    = mat_get(X,i_s,0)/SQ(dimfac) - SQ(param->M_s);
    umd    = mat_get(X,i_umd,0)/SQ(dimfac) - param->M_umd;
    Linv   = mat_get(X,i_Linv,0)/dimfac;
    a2mud  = SQ(a)*mat_get(X,i_ud,0)/SQ(dimfac);
    a2ms   = SQ(a)*mat_get(X,i_s,0)/SQ(dimfac);
    alpha  = mat_get(X,i_alpha,0);
    dalpha = alpha - param->alpha;
    
    /* constant term (extrapolated quantity) */
    res += mat_get(p,s,0);
    /* Taylor expansion around the physical isospin masses */
    polynom(res,p,I_ud(0)+s,M_ud,param->M_ud_deg);
    polynom(res,p,I_s(0)+s,M_s,param->M_s_deg);
    /* m_u-m_d isospin breaking terms */
    polynom(res,p,I_umd(0)+s,umd,param->umd_deg);
    /* m_q*(m_u-m_d) isospin breaking terms */
    if (param->with_udumd)
    {
        res += mat_get(p,I_udumd(0)+s,0)*M_ud*umd;
    }
    if (param->with_sumd)
    {
        res += mat_get(p,I_sumd(0)+s,0)*M_s*umd;
    }
    /* alpha isospin breaking terms */
    polynom(res,p,I_alpha(0)+s,dalpha,param->alpha_deg);
    /* m_q*alpha isospin breaking terms */
    if (param->with_udalpha)
    {
        res += mat_get(p,I_udalpha(0)+s,0)*M_ud*dalpha;
    }
    if (param->with_salpha)
    {
        res += mat_get(p,I_salpha(0)+s,0)*M_s*dalpha;
    }
    /* discretization effects */
    polynom(res,p,I_discr(0)+s,a,param->a_deg);
    if (param->with_a2M_ud)
    {
        res += mat_get(p,I_dis_a2M_ud(0)+s,0)*a2mud;
    }
    if (param->with_a2M_s)
    {
        res += mat_get(p,I_dis_a2M_s(0)+s,0)*a2ms;
    }
    /* QED finite volume effect */
    polynom(res,p,I_qedfv(0)+s,alpha*Linv,param->with_qed_fvol);
    /* dimensional factor */
    res *= pow(dimfac,param->q_dim);
    
    return res;
}

size_t fm_phypt_taylor_npar(void* vparam)
{
    fit_param *param;
    size_t npar;
    
    param = (fit_param *)vparam;
    
    npar  = 1;
    npar += param->M_ud_deg;
    npar += param->M_s_deg;
    npar += param->umd_deg;
    if (param->with_udumd)
    {
        npar++;
    }
    if (param->with_sumd)
    {
        npar++;
    }
    npar += param->alpha_deg;
    if (param->with_udalpha)
    {
        npar++;
    }
    if (param->with_salpha)
    {
        npar++;
    }
    npar += param->a_deg;
    if (param->with_a2M_ud) 
    {
        npar++;
    }
    if (param->with_a2M_s) 
    {
        npar++;
    }
    npar += param->with_qed_fvol;
    if (IS_AN(param,AN_PHYPT)&&!IS_AN(param,AN_SCALE)&&(param->with_ext_a))
    {
        npar += param->nbeta;
    }
    
    return npar;
}

#undef I_ud
#undef I_s
#undef I_discr
#undef I_dis_a2M_ud
#undef I_dis_a2M_s
#undef I_umd
#undef I_udumd
#undef I_sumd
#undef I_alpha
#undef I_udalpha
#undef I_salpha
#undef I_qedfv

#define I_ud(i)         (i+(param)->nbeta)
#define I_s(i)          (I_ud(i)+(param)->s_M_ud_deg)
#define I_dis_a2M_ud(i) (I_s(i)+(param)->s_M_s_deg)
#define I_dis_a2M_s(i)  (I_dis_a2M_ud(i)+(((param)->s_with_a2M_ud) ? 1 : 0))
#define I_umd(i)        (I_dis_a2M_s(i)+(((param)->s_with_a2M_s) ? 1 : 0))
#define I_qedfv(i)      (I_umd(i)+(param)->s_umd_deg)

double fm_scaleset_taylor_func(const mat *X, const mat *p, void *vparam)
{
    double res,a,M_ud,M_s,umd,M_scale,Linv,a2mud,a2ms;
    fit_param *param;
    size_t bind;
    
    param   = (fit_param *)vparam;
    res     = 0.0;
    bind    = (size_t)(mat_get(X,i_bind,0));
    M_scale = param->M_scale;
    a       = mat_get(p,bind,0);
    a2ms    = mat_get(X,i_s,0)/M_scale;
    a2mud   = mat_get(X,i_ud,0)/M_scale;
    Linv    = mat_get(X,i_Linv,0)/(a*M_scale);
    M_ud    = mat_get(X,i_ud,0)/SQ(a*M_scale)-SQ(param->M_ud)/SQ(M_scale);
    M_s     = mat_get(X,i_s,0)/SQ(a*M_scale)-SQ(param->M_s)/SQ(M_scale);
    umd     = mat_get(X,i_umd,0)/SQ(a*M_scale)-param->M_umd/SQ(M_scale);
    
    res += 1.0;
    polynom(res,p,I_ud(0),M_ud,param->s_M_ud_deg);
    polynom(res,p,I_s(0),M_s,param->s_M_s_deg);
    polynom(res,p,I_umd(0),umd,param->s_umd_deg);
    if (param->s_with_qed_fvol > 0)
    {
        res += mat_get(p,I_qedfv(0),0)*SQ(Linv);
    }
    if (param->s_with_a2M_ud)
    {
        res += mat_get(p,I_dis_a2M_ud(0),0)*a2mud;
    }
    if (param->s_with_a2M_s)
    {
        res += mat_get(p,I_dis_a2M_s(0),0)*a2ms;
    }
    res *= a*M_scale;
    
    return res;
}

size_t fm_scaleset_taylor_npar(void* vparam)
{
    fit_param *param;
    size_t npar;
    
    param = (fit_param *)vparam;
    
    npar  = param->nbeta;
    npar += param->s_M_ud_deg;
    npar += param->s_M_s_deg;
    npar += param->s_umd_deg;
    if (param->s_with_qed_fvol)
    {
        npar++;
    }
    if (param->s_with_a2M_ud)
    {
        npar++;
    }
    if (param->s_with_a2M_s)
    {
        npar++;
    }
    
    return npar;
}


size_t fm_comb_phypt_taylor_scaleset_taylor_npar(void* vparam)
{
    return fm_phypt_taylor_npar(vparam) + fm_scaleset_taylor_npar(vparam);
}
