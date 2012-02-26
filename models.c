#include "models.h"
#include "data_loader.h"
#include "parameters.h"
#include <math.h>
#include <string.h>
#include <latan/latan_nunits.h>
#include <latan/latan_mass.h>
#include <latan/latan_math.h>

/* physical point Taylor expansion models */
#define ud_s_taylor(res,p,s,M_ud,M_ud_deg,M_s,M_s_deg,param)\
{\
    int d_;\
    size_t i_;\
    if ((M_ud_deg) > 0)\
    {\
        res += mat_get(p,I_ud(0)+(s),0)*M_ud;\
        for (d_=2;d_<=(M_ud_deg);d_++)\
        {\
            i_    = (size_t)(d_) - 1;\
            res  += mat_get(p,I_ud(i_)+(s),0)*pow(M_ud,(double)(d_));\
        }\
    }\
    if ((M_s_deg) > 0)\
    {\
        res += mat_get(p,I_s(0)+(s),0)*M_s;\
        for (d_=2;d_<=(M_s_deg);d_++)\
        {\
            i_   = (size_t)(d_) - 1;\
            res += mat_get(p,I_s(i_)+(s),0)*pow(M_s,(double)(d_));\
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
            sprintf(buf_,"+%e*(%s-%e)**%d",mat_get(p,I_ud(i_)+(s),0),\
                    x_str[i_ud],M_ud_ex,d_);\
            strcat(str,buf_);\
        }\
    }\
    if ((M_s_deg) > 0)\
    {\
        for (d_=1;d_<=(M_s_deg);d_++)\
        {\
            i_ = (size_t)(d_) - 1;\
            sprintf(buf_,"+%e*(%s-%e)**%d",mat_get(p,I_s(i_)+(s),0),\
                    x_str[i_s],M_s_ex,d_);\
            strcat(str,buf_);\
        }\
    }\
}

#define I_ud(i)         (i+1)
#define I_s(i)          (I_ud(i)+(param)->M_ud_deg)
#define I_discr(i)      (I_s(i)+(param)->M_s_deg)
#define I_dis_a2M_ud(i) (I_discr(i)+(((param)->a_deg > 0) ? 1 : 0))
#define I_dis_a2M_s(i)  (I_dis_a2M_ud(i)+((param)->with_a2M_ud ? 1 : 0))
#define I_umd(i)        (I_dis_a2M_s(i)+((param)->with_a2M_s ? 1 : 0))
#define I_qedfv(i)      (I_umd(i)+((param)->with_umd ? 1 : 0))
#define I_a(i)          (I_qedfv(i)+((param)->with_qed_fvol ? 1 : 0))

double a_error_chi2_ext(mat *p, void *vd)
{
    size_t i,s;
    double a,a_err,res;
    fit_param *param;
    
    res   = 0.0;
    param = (fit_param *)(((fit_data *)vd)->model_param); 
    s     = fit_data_get_sample_counter((fit_data *)vd);
    
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
    }
    
    return res;
}

double fm_phypt_a_taylor_func(const mat *X, const mat *p, void *vparam)
{
    double res,M_ud,M_ud2,M_s,a,dimfac,umd,Linv,a2mud,a2ms;
    size_t s,bind;
    fit_param *param;
    
    param = (fit_param *)vparam;
    res   = 0.0;
    bind  = (size_t)(mat_get(X,i_bind,0));
    
    /* what is the lattice spacing ? */
    if (IS_ANALYZE(param,"phypt"))
    {
        /** a global parameter with error (with external scale samples) **/
        if (param->with_ext_a)
        {
            a = (!param->plotting) ? mat_get(p,I_a(bind),0) : mat_get(X,i_a,0);
        }
        /** a x coordinate (with the ratio method) **/
        else
        {
            a = mat_get(X,i_a,0);
        }
        s = 0;
    }
    else if (IS_ANALYZE(param,"comb_phypt_scale"))
    {
        /** a fit parameter from the scale setting model **/
        a = (!param->plotting) ? mat_get(p,bind,0) : mat_get(X,i_a,0);
        s = fm_scaleset_taylor_npar(param);
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
    M_ud2  = mat_get(X,i_ud,0)/SQ(dimfac);
    M_s    = mat_get(X,i_s,0)/SQ(dimfac) - SQ(param->M_s);
    umd    = mat_get(X,i_umd,0)/SQ(dimfac) - param->M_umd;
    Linv   = mat_get(X,i_Linv,0)/dimfac;
    a2mud  = SQ(a)*mat_get(X,i_ud,0)/SQ(dimfac);
    a2ms   = SQ(a)*mat_get(X,i_s,0)/SQ(dimfac);
    
    /* constant term (extrapolated quantity) */
    res += mat_get(p,s,0);
    /* Taylor expansion around the physical masses in m_ud and m_s */
    ud_s_taylor(res,p,s,M_ud,param->M_ud_deg,M_s,param->M_s_deg,param);
    /* discretization effects */
    if (param->a_deg > 0)
    {
        res += mat_get(p,I_discr(0)+s,0)*pow(a,(double)(param->a_deg));
    }
    if (param->with_a2M_ud)
    {
        res += mat_get(p,I_dis_a2M_ud(0)+s,0)*a2mud;
    }
    if (param->with_a2M_s)
    {
        res += mat_get(p,I_dis_a2M_s(0)+s,0)*a2ms;
    }
    /* mass isospin breaking effects */
    if (param->with_umd)
    {
        res += mat_get(p,I_umd(0)+s,0)*umd;
    }
    /* QED finite volume effect */
    if (param->with_qed_fvol)
    {
        res += mat_get(p,I_qedfv(0)+s,0)*SQ(Linv);
    }
    /* dimensional factor */
    res *= pow(dimfac,param->q_dim);
    
    return res;
}

size_t fm_phypt_a_taylor_npar(void* vparam)
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
    if (param->with_a2M_ud) 
    {
        npar++;
    }
    if (param->with_a2M_s) 
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
    if (IS_ANALYZE(param,"phypt")&&(param->with_ext_a))
    {
        npar += param->nbeta;
    }
    
    return npar;
}

void fm_phypt_a_taylor_pstr(strbuf str, const size_t i,   \
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
        sprintf(buf,"+%e*%s**%d",mat_get(p,I_discr(0)+s,0),x_str[i_a],\
                param->a_deg);
        strcat(str,buf);
    }
    if (param->with_a2M_ud)
    {
        sprintf(buf,"+%e*%s**2*%s",mat_get(p,I_dis_a2M_ud(0)+s,0),x_str[i_a],\
                x_str[i_ud]);
        strcat(str,buf);
    }
    if (param->with_a2M_s)
    {
        sprintf(buf,"+%e*%s**2*%s",mat_get(p,I_dis_a2M_s(0)+s,0),x_str[i_a],\
                x_str[i_s]);
        strcat(str,buf);
    }
    if (param->with_umd)
    {
        sprintf(buf,"+%e*(%s-%e)",mat_get(p,I_umd(0)+s,0),\
                x_str[i_umd],param->M_umd);
        strcat(str,buf);
    }
    if (param->with_qed_fvol)
    {
        sprintf(buf,"+%e*%s**2",mat_get(p,I_qedfv(0)+s,0),\
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

#undef I_ud
#undef I_s
#undef I_discr
#undef I_dis_a2M_ud
#undef I_dis_a2M_s
#undef I_umd
#undef I_qedfv
#define I_ud(i)         (i+(param)->nbeta)
#define I_s(i)          (I_ud(i)+(param)->s_M_ud_deg)
#define I_dis_a2M_ud(i) (I_s(i)+(param)->s_M_s_deg)
#define I_dis_a2M_s(i)  (I_dis_a2M_ud(i)+(((param)->s_with_a2M_ud) ? 1 : 0))
#define I_umd(i)        (I_dis_a2M_s(i)+(((param)->s_with_a2M_s) ? 1 : 0))
#define I_qedfv(i)      (I_umd(i)+((param)->s_with_umd ? 1 : 0))

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
    ud_s_taylor(res,p,0,M_ud,param->s_M_ud_deg,M_s,param->s_M_s_deg,param);
    if (param->s_with_umd)
    {
        res += mat_get(p,I_umd(0),0)*umd;
    }
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
    if (param->s_with_umd)
    {
        npar++;
    }
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

void fm_scaleset_taylor_pstr(strbuf str, const size_t i,   \
                             const mat *x_ex, const mat *p,\
                             void *vparam)
{
    fit_param *param;
    strbuf buf,x_str[N_EX_VAR];
    double M_ud_phi,M_s_phi,M_umd_phi,M_scale,a,fac;
    size_t bind;
    size_t j;
    
    param     = (fit_param *)vparam;
    bind      = (size_t)(mat_get(x_ex,i_bind,0));
    a         = mat_get(p,bind,0);
    M_scale   = param->M_scale;
    M_ud_phi  = SQ(param->M_ud)/SQ(M_scale);
    M_s_phi   = SQ(param->M_s)/SQ(M_scale);
    M_umd_phi = param->M_umd/SQ(M_scale);
    
    for (j=0;j<N_EX_VAR;j++)
    {
        if (j == i_Linv)
        {
            fac = a*M_scale;
        }
        else
        {
            fac = SQ(a*M_scale);
        }
        if (i == j)
        {
            sprintf(x_str[j],"(x/%e)",fac);
        }
        else
        {
            sprintf(x_str[j],"%e",mat_get(x_ex,j,0)/fac);
        }
    }
    sprintf(str,"%e*(1.0",a*M_scale);
    ud_s_taylor_pstr(str,p,0,x_str,M_ud_phi,param->s_M_ud_deg,M_s_phi,\
                     param->s_M_s_deg,param);
    if (param->s_with_umd)
    {
        sprintf(buf,"+%e*(%s-%e)",mat_get(p,I_umd(0),0),x_str[i_umd],\
                M_umd_phi);
        strcat(str,buf);
    }
    if (param->s_with_qed_fvol > 0)
    {
        sprintf(buf,"+%e*%s**2",mat_get(p,I_qedfv(0),0),x_str[i_Linv]);
        strcat(str,buf);
    }
    if (param->s_with_a2M_ud)
    {
        sprintf(buf,"+%e*%e*%s",mat_get(p,I_dis_a2M_ud(0),0),SQ(a)*M_scale,\
                x_str[i_ud]);
        strcat(str,buf);
    }
    if (param->s_with_a2M_s)
    {
        sprintf(buf,"+%e*%e*%s",mat_get(p,I_dis_a2M_s(0),0),SQ(a)*M_scale,\
                x_str[i_s]);
        strcat(str,buf);
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

size_t fm_comb_phyt_scale_npar(void* vparam)
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
