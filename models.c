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

static void polynom(double *res, const mat *par, const size_t i0,\
                    const double x, const int deg)
{
    int d;
    size_t i;
    
    if ((deg) > 0)
    {
        *res += mat_get(par,i0,0)*(x);
        for (d=2;d<=(deg);d++)
        {
            i    = (size_t)(d) - 1;
            *res += mat_get(par,(i0)+i,0)*pow(x,(double)(d));
        }
    }
}

#define I_ud(i)      (i+((param)->with_const ? 1 : 0))
#define I_s(i)       (I_ud(i)+(param)->M_ud_deg)
#define I_a2(i)      (I_s(i)+(param)->M_s_deg)
#define I_a2ud(i)    (I_a2(i)+((param)->with_a2 ? 1 : 0))
#define I_a2s(i)     (I_a2ud(i)+((param)->with_a2ud ? 1 : 0))
#define I_umd(i)     (I_a2s(i)+((param)->with_a2s ? 1 : 0))
#define I_udumd(i)   (I_umd(i)+(param)->umd_deg) 
#define I_sumd(i)    (I_udumd(i)+(param)->with_udumd)
#define I_a2umd(i)   (I_sumd(i)+(param)->with_sumd)
#define I_alpha(i)   (I_a2umd(i)+((param)->with_a2umd ? 1 : 0))
#define I_udalpha(i) (I_alpha(i)+(param)->alpha_deg)
#define I_salpha(i)  (I_udalpha(i)+(param)->with_udalpha)
#define I_aalpha(i)  (I_salpha(i)+(param)->with_salpha)
#define I_qedfv(i)   (I_aalpha(i)+((param)->with_aalpha ? 1 : 0))
#define LAST_IND     (I_qedfv(0)+param->with_qed_fvol-1-param->with_qed_fvol_monopmod)
#define I_a_val(i)   (LAST_IND+1+i)

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
            a = mat_get(rs_sample_pt_cent_val(param->s_a),i,0);
        }
        else
        {
            a = mat_get(rs_sample_pt_sample(param->s_a,s-1),i,0);
        }
        a_err = mat_get(param->a_err,i,0);
        res  += SQ(mat_get(p,I_a_val(i),0) - a)/SQ(a_err);
        d->matperf += 5.0;
    }
    
    return res;
}

double fm_phypt_taylor_func(const mat *X, const mat *p, void *vparam)
{
    double res,lo,ex,buf,M_ud,M_s,M_fvol,a,dimfac,umd,Linv,a2mud,a2ms,alpha,d;
    size_t s,bind;
    fit_param *param;
    
    param = (fit_param *)vparam;
    bind  = (size_t)(mat_get(X,i_bind,0));
    
    /* what is the lattice spacing ? */
    if (IS_AN(param,AN_PHYPT))
    {
        if (IS_AN(param,AN_SCALE))
        {
            /** a fit parameter from the scale setting model **/
            a = (!param->scale_model) ? mat_get(p,bind,0) : mat_get(X,i_a,0);
            s = fm_scaleset_taylor_npar(param);
        }
        else
        {
            /** a global parameter with error (with external scale samples) **/
            if (param->with_ext_a)
            {
                a = (!param->scale_model) ? mat_get(p,I_a_val(bind),0) :\
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
    dimfac = (!param->scale_model) ? a : 1.0;
    M_ud   = mat_get(X,i_ud,0)/SQ(dimfac);
    M_s    = mat_get(X,i_s,0)/SQ(dimfac);
    umd    = mat_get(X,i_umd,0)/SQ(dimfac);
    Linv   = mat_get(X,i_Linv,0)/dimfac;
    a2mud  = SQ(a)*mat_get(X,i_ud,0)/SQ(dimfac);
    a2ms   = SQ(a)*mat_get(X,i_s,0)/SQ(dimfac);
    alpha  = mat_get(X,i_alpha,0);
    M_fvol = mat_get(X,i_fvM,0)/dimfac;
    d      = (double)param->q_dim;
    
    res = 0.0;
    /* Taylor/Pade mass expansion */
    lo = (param->with_const) ? mat_get(p,s,0) : 0.0;
    ex = 0.0;
    polynom(&ex,p,I_ud(0)+s,M_ud-SQ(param->M_ud),param->M_ud_deg);
    polynom(&ex,p,I_s(0)+s,M_s-SQ(param->M_s),param->M_s_deg);
    res += (param->with_pade) ? lo/(1.0-ex/lo) : lo+ex;
    /* discretization effects */
    res += (param->with_a2)   ? mat_get(p,I_a2(0)+s,0)*SQ(a)   : 0.0;
    res += (param->with_a2ud) ? mat_get(p,I_a2ud(0)+s,0)*a2mud : 0.0;
    res += (param->with_a2s)  ? mat_get(p,I_a2s(0)+s,0)*a2ms   : 0.0;
    /* Taylor/Pade m_u-m_d expansion */
    lo = (param->umd_deg) ? mat_get(p,I_umd(0)+s,0) : 0.0;
    ex = 0.0;
    polynom(&ex,p,I_udumd(0)+s,M_ud-SQ(param->M_ud),param->with_udumd);
    polynom(&ex,p,I_sumd(0)+s,M_s-SQ(param->M_s),param->with_sumd);
    ex  += (param->with_a2umd) ? mat_get(p,I_a2umd(0)+s,0)*SQ(a) : 0.0;
    res += (param->with_umd_pade) ? umd*lo/(1.0-ex/lo) : umd*lo+umd*ex;
    /* Taylor/Pade alpha expansion */
    lo = (param->alpha_deg) ? mat_get(p,I_alpha(0)+s,0) : 0.0;
    ex = 0.0;
    polynom(&ex,p,I_udalpha(0)+s,M_ud-SQ(param->M_ud),param->with_udalpha);
    polynom(&ex,p,I_salpha(0)+s,M_s-SQ(param->M_s),param->with_salpha);
    ex  += (param->with_aalpha) ? mat_get(p,I_aalpha(0)+s,0)*a : 0.0;
    /** QED finite volume effects **/
    if (param->have_alpha)
    {
        buf  = 0.0;
        if (param->with_qed_fvol_monopmod)
        {
            
            if (param->with_qed_fvol >= 1)
            {
                buf += -0.5*QED_FVOL_KAPPA*Linv*d*pow(M_fvol,d-1.0);
            }
            if (param->with_qed_fvol >= 2)
            {
                buf += mat_get(p,I_qedfv(0)+s,0)*SQ(Linv)*pow(M_fvol,d-2.0);
            }
            buf *= (double)param->qed_fvol_monopmod_sign;
            res += alpha*buf;
        }
        else
        {
            polynom(&buf,p,I_qedfv(0)+s,Linv/M_fvol,param->with_qed_fvol);
            buf *= pow(M_fvol,d);
            ex  += buf;
        }
    }
    res += (param->with_alpha_pade) ? alpha*lo/(1.0-ex/lo) : alpha*lo+alpha*ex;
    
    /* dimensional factor */
    res *= pow(dimfac,param->q_dim);
    
    return res;
}

size_t fm_phypt_taylor_npar(void* vparam)
{
    fit_param *param;
    size_t npar;
    
    param = (fit_param *)vparam;
    
    npar  = LAST_IND + 1;
    if (IS_AN(param,AN_PHYPT)&&!IS_AN(param,AN_SCALE)&&(param->with_ext_a))
    {
        npar += param->nbeta;
    }
    
    return npar;
}

#undef I_ud
#undef I_s
#undef I_a2
#undef I_aalpha
#undef I_a2ud
#undef I_a2s
#undef I_umd
#undef I_udumd
#undef I_sumd
#undef I_alpha
#undef I_udalpha
#undef I_salpha
#undef I_qedfv
#undef LAST_IND
#undef I_a

#define I_ud(i)   (i+(param)->nbeta)
#define I_s(i)    (I_ud(i)+(param)->s_M_ud_deg)
#define I_a2ud(i) (I_s(i)+(param)->s_M_s_deg)
#define I_a2s(i)  (I_a2ud(i)+(((param)->s_with_a2ud) ? 1 : 0))
#define I_umd(i)  (I_a2s(i)+(((param)->s_with_a2s) ? 1 : 0))
#define LAST_IND  (I_umd(0)+(param)->s_umd_deg-1)

double fm_scaleset_taylor_func(const mat *X, const mat *p, void *vparam)
{
    double res,a,M_ud,M_s,umd,M_scale,a2mud,a2ms;
    fit_param *param;
    size_t bind;
    
    param   = (fit_param *)vparam;
    res     = 0.0;
    bind    = (size_t)(mat_get(X,i_bind,0));
    M_scale = param->M_scale;
    a       = mat_get(p,bind,0);
    a2ms    = mat_get(X,i_s,0)/M_scale;
    a2mud   = mat_get(X,i_ud,0)/M_scale;
    M_ud    = mat_get(X,i_ud,0)/SQ(a*M_scale)-SQ(param->M_ud)/SQ(M_scale);
    M_s     = mat_get(X,i_s,0)/SQ(a*M_scale)-SQ(param->M_s)/SQ(M_scale);
    umd     = mat_get(X,i_umd,0)/SQ(a*M_scale)-param->M_umd_val/SQ(M_scale);
    
    res += 1.0;
    polynom(&res,p,I_ud(0),M_ud,param->s_M_ud_deg);
    polynom(&res,p,I_s(0),M_s,param->s_M_s_deg);
    polynom(&res,p,I_umd(0),umd,param->s_umd_deg);
    if (param->s_with_a2ud)
    {
        res += mat_get(p,I_a2ud(0),0)*a2mud;
    }
    if (param->s_with_a2s)
    {
        res += mat_get(p,I_a2s(0),0)*a2ms;
    }
    res *= a*M_scale;
    
    return res;
}

size_t fm_scaleset_taylor_npar(void* vparam)
{
    fit_param *param;
    
    param = (fit_param *)vparam;
    
    return LAST_IND + 1;
}

#undef I_ud
#undef I_s
#undef I_umd
#undef LAST_IND

size_t fm_comb_phypt_taylor_scaleset_taylor_npar(void* vparam)
{
    return fm_phypt_taylor_npar(vparam) + fm_scaleset_taylor_npar(vparam);
}
