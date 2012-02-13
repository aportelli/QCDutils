#include "output.h"
#include <math.h>
#include <string.h>
#include <latan/latan_math.h>
#include <latan/latan_plot.h>
#include "data_loader.h"
#include "models.h"

#define PLOT_ADD_FIT(obj,kx,ky,title,color)\
plot_add_fit(p[kx],d,ky,phy_pt,kx,fit,true,obj,title,"",color,color)

#define PLOT_ADD_EX(kx,s)\
plot_add_point(p[kx],mat_get(phy_pt,kx,0),param->q_target[0],      \
               -1.0,param->q_target[1],"target","rgb 'dark-blue'");\
plot_add_point(p[kx],mat_get(phy_pt,kx,0),mat_get(fit,s,0),                   \
               -1.0,sqrt(mat_get(fit_var,s,0)),"physical point","rgb 'black'");

#define PLOT_DISP(kx)\
plot_set_title(p[kx],gtitle);\
plot_set_xlabel(p[kx],xlabel);\
plot_set_ylabel(p[kx],ylabel);\
plot_disp(p[kx]);


void plot_fit(const mat *fit, const mat *fit_var, fit_data *d, \
              fit_param *param, plot_flag f)
{
    plot *p[N_EX_VAR];
    double *xb[N_EX_VAR] = {NULL,NULL,NULL,NULL,NULL,NULL};
    double b_int[2],dbind,M_scale,a;
    size_t bind,k,phy_ind,s;
    strbuf color,gtitle,title,xlabel,ylabel;
    mat *phy_pt;
    
    phy_pt = mat_create(N_EX_VAR,1);
    for (k=0;k<N_EX_VAR;k++)
    {
        p[k] = plot_create();
    }
    
    if (IS_ANALYZE(param,"comb_phypt_scale"))
    {
        phy_ind = 1;
        s       = fm_scaleset_taylor_npar(param);
    }
    else
    {
        phy_ind = 0;
        s       = 0;
    }
    
    param->plotting = 1;
    if (f == Q)
    {
        sprintf(gtitle,"quantity: %s -- scale: %s -- datasets: %s -- ensembles: %s",
                param->q_name,param->scale_part,param->dataset_cat,param->manifest);
        mat_set(phy_pt,i_ud,0,SQ(param->M_ud));
        mat_set(phy_pt,i_s,0,SQ(param->M_s));
        mat_set(phy_pt,i_umd,0,param->M_umd);
        mat_set(phy_pt,i_bind,0,0.0);
        mat_set(phy_pt,i_a,0,0.0);
        mat_set(phy_pt,i_Linv,0,0.0);
        for (bind=0;bind<param->nbeta;bind++)
        {
            dbind      = (double)(bind);
            b_int[0]   = dbind - 0.1;
            b_int[1]   = dbind + 0.1;
            xb[i_bind] = b_int;
            fit_data_fit_region(d,xb);
            sprintf(color,"%d",1+(int)bind);
            sprintf(title,"beta = %s",param->beta[bind]);
            PLOT_ADD_FIT(PF_DATA,i_ud,phy_ind,title,color);
            PLOT_ADD_FIT(PF_DATA,i_s,phy_ind,title,color);
            PLOT_ADD_FIT(PF_DATA,i_umd,phy_ind,title,color);
            PLOT_ADD_FIT(PF_DATA,i_a,phy_ind,title,color);
            PLOT_ADD_FIT(PF_DATA,i_Linv,phy_ind,title,color);
            fit_data_fit_all_points(d,true);
        }
        PLOT_ADD_FIT(PF_FIT,i_ud,phy_ind,"","rgb 'black'");
        PLOT_ADD_FIT(PF_FIT,i_s,phy_ind,"","rgb 'black'");
        PLOT_ADD_FIT(PF_FIT,i_umd,phy_ind,"","rgb 'black'");
        PLOT_ADD_FIT(PF_FIT,i_a,phy_ind,"","rgb 'black'");
        PLOT_ADD_FIT(PF_FIT,i_Linv,phy_ind,"","rgb 'black'");
        switch (param->q_dim) 
        {
            case 0:
                strbufcpy(ylabel,param->q_name);
                break;
            case 1:
                sprintf(ylabel,"%s (MeV)",param->q_name);
                break;
            default:
                sprintf(ylabel,"%s^%d (MeV^%d)",param->q_name,param->q_dim,\
                        param->q_dim);
                break;
        }
        
        sprintf(xlabel,"M_%s^2 (MeV^2)",param->ud_name);
        PLOT_ADD_EX(i_ud,s);
        PLOT_DISP(i_ud);
        sprintf(xlabel,"M_%s^2 (MeV^2)",param->s_name);
        PLOT_ADD_EX(i_s,s);
        PLOT_DISP(i_s);
        strbufcpy(xlabel,"a (MeV^-1)");
        PLOT_ADD_EX(i_a,s);
        PLOT_DISP(i_a);
        if (param->have_umd)
        {
            strbufcpy(xlabel,"M_K_p^2 - M_K_0^2 (MeV^2)");
            PLOT_ADD_EX(i_umd,s);
            PLOT_DISP(i_umd);
        }
        strbufcpy(xlabel,"1/L (MeV)");
        PLOT_ADD_EX(i_Linv,s);
        PLOT_DISP(i_Linv);
    }
    else if (f == SCALE)
    {
        sprintf(gtitle,"scale setting: %s -- datasets: %s -- ensembles: %s",
                param->scale_part,param->dataset_cat,param->manifest);
        for (bind=0;bind<param->nbeta;bind++)
        {
            dbind      = (double)(bind);
            b_int[0]   = dbind - 0.1;
            b_int[1]   = dbind + 0.1;
            xb[i_bind] = b_int;
            M_scale    = param->M_scale;
            a          = mat_get(fit,bind,0);
            fit_data_fit_region(d,xb);
            mat_set(phy_pt,i_ud,0,SQ(a*param->M_ud));
            mat_set(phy_pt,i_s,0,SQ(a*param->M_s));
            mat_set(phy_pt,i_umd,0,SQ(a)*param->M_umd);
            mat_set(phy_pt,i_bind,0,bind);
            mat_set(phy_pt,i_a,0,a);
            mat_set(phy_pt,i_Linv,0,0.0);
            sprintf(color,"%d",1+(int)bind);
            sprintf(title,"beta = %s",param->beta[bind]);
            PLOT_ADD_FIT(PF_DATA|PF_FIT,i_ud,0,title,color);
            PLOT_ADD_FIT(PF_DATA|PF_FIT,i_s,0,title,color);
            PLOT_ADD_FIT(PF_DATA|PF_FIT,i_umd,0,title,color);
            PLOT_ADD_FIT(PF_DATA|PF_FIT,i_Linv,0,title,color);
            fit_data_fit_all_points(d,true);
        }
        sprintf(ylabel,"(a*M_%s)^2",param->scale_part);
        sprintf(xlabel,"(a*M_%s)^2",param->ud_name);
        PLOT_DISP(i_ud);
        sprintf(xlabel,"(a*M_%s)^2",param->s_name);
        PLOT_DISP(i_s);
        if (param->have_umd)
        {
            strbufcpy(xlabel,"(a*M_K_p)^2 - (a*M_K_0)^2");
            PLOT_DISP(i_umd);
        }
        strbufcpy(xlabel,"a/L");
        PLOT_DISP(i_Linv);
    }
    param->plotting = 0;
    
    mat_destroy(phy_pt);
    for (k=0;k<N_EX_VAR;k++)
    {
        plot_destroy(p[k]);
    }
}
#undef ADD_PLOT


#define PRINT_PAR(name)\
{\
    printf("%10s = %e (%.4f%%)\n",name,mat_get(fit,i,0),            \
           sqrt(mat_get(fit_var,i,0))/fabs(mat_get(fit,i,0))*100.0);\
    i++;\
}

void print_result(const rs_sample *s_fit, fit_param *param)
{
    size_t i;
    int j;
    mat *fit,*fit_var;
    strbuf buf;
    
    i   = 0;
    fit = rs_sample_pt_cent_val(s_fit);
    
    fit_var = mat_create(fit_model_get_npar(param->model,param),1);
    
    rs_sample_varp(fit_var,s_fit);
    printf("\nfit parameters (m_ud:%d m_s:%d m_u-m_d:%d a:%d) :\n",     \
           param->M_ud_deg,param->M_s_deg,param->with_umd,param->a_deg);
    if (IS_ANALYZE(param,"scaleset")||IS_ANALYZE(param,"comb_phypt_scale"))
    {
        for(j=0;j<(int)param->nbeta;j++)
        {
            sprintf(buf,"a_%s",param->beta[i]);
            PRINT_PAR(buf);
        }
        if (param->s_M_ud_deg > 0)
        {
            for (j=0;j<param->s_M_ud_deg;j++)
            {
                sprintf(buf,"s_p_ud_%d",j+1);
                PRINT_PAR(buf);
            }
        }
        if (param->s_M_s_deg > 0)
        {
            for (j=0;j<param->s_M_s_deg;j++)
            {
                sprintf(buf,"s_p_s_%d",j+1);
                PRINT_PAR(buf);
            }
        }
        if (param->s_with_aM_s)
        {
            PRINT_PAR("s_aM_s");
        }
        if (param->s_with_umd)
        {
            PRINT_PAR("s_p_miso");
        }
        if (param->s_with_qed_fvol)
        {
            PRINT_PAR("s_p_fvol_L");
        }
    }
    if (IS_ANALYZE(param,"phypt")||IS_ANALYZE(param,"comb_phypt_scale"))
    {
        PRINT_PAR(param->q_name);
        if (param->M_ud_deg > 0)
        {
            for (j=0;j<param->M_ud_deg;j++)
            {
                sprintf(buf,"p_ud_%d",j+1);
                PRINT_PAR(buf);
            }
        }
        if (param->M_s_deg > 0)
        {
            for (j=0;j<param->M_s_deg;j++)
            {
                sprintf(buf,"p_s_%d",j+1);
                PRINT_PAR(buf);
            }
        }
        if (param->a_deg > 0)
        {
            PRINT_PAR("p_a");
        }
        if (param->with_umd)
        {
            PRINT_PAR("p_miso");
        }
        if (param->with_qed_fvol)
        {
            PRINT_PAR("p_fvol_L");
        }
    }
    printf("\n");
    
    mat_destroy(fit_var);
}
#undef PRINT_PAR
