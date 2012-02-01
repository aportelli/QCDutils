#include "data_loader.h"
#include "models.h"
#include "parameters.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include <latan/latan_fit.h>
#include <latan/latan_io.h>
#include <latan/latan_math.h>
#include <latan/latan_minimizer.h>
#include <latan/latan_plot.h>

typedef enum
{
    Q     = 0,
    SCALE = 1
} plot_flag;

static void plot_fit(const mat *fit, const mat *fit_var, fit_data *d, \
                     fit_param *param, plot_flag f);
static void print_result(const rs_sample *s_fit, fit_param *param);

#define PLOT_ADD_FIT(obj,kx,ky,title,color)\
plot_add_fit(p[kx],d,ky,phy_pt,kx,fit,true,obj,title,"",color,color)

#define PLOT_ADD_EX(kx,s)\
plot_add_point(p[kx],mat_get(phy_pt,kx,0),param->q_target[0],\
               -1.0,param->q_target[1],"target","rgb 'dark-blue'");\
plot_add_point(p[kx],mat_get(phy_pt,kx,0),mat_get(fit,s,0),\
               -1.0,sqrt(mat_get(fit_var,s,0)),"physical point","rgb 'black'");

#define PLOT_DISP(kx)\
plot_set_title(p[kx],gtitle);\
plot_set_xlabel(p[kx],xlabel);\
plot_set_ylabel(p[kx],ylabel);\
plot_disp(p[kx]);


static void plot_fit(const mat *fit, const mat *fit_var, fit_data *d, \
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
        mat_set(phy_pt,i_umd,0,DMSQ_K);
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
        plot_print(p[i_a]);
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
        if (param->M_ud_deg > 0)
        {
            sprintf(xlabel,"M_%s^2 (MeV^2)",param->ud_name);
            PLOT_ADD_EX(i_ud,s);
            PLOT_DISP(i_ud);
        }
        if (param->M_s_deg > 0)
        {
            sprintf(xlabel,"M_%s^2 (MeV^2)",param->s_name);
            PLOT_ADD_EX(i_s,s);
            PLOT_DISP(i_s);
        }
        if (param->a_deg > 0)
        {
            strbufcpy(xlabel,"a (MeV^-1)");
            PLOT_ADD_EX(i_a,s);
            PLOT_DISP(i_a);
        }
        if (param->with_umd)
        {
            strbufcpy(xlabel,"M_K_p^2 - M_K_0^2 (MeV^2)");
            PLOT_ADD_EX(i_umd,s);
            PLOT_DISP(i_umd);
        }
        if (param->with_qed_fvol)
        {
            strbufcpy(xlabel,"1/L (MeV)");
            PLOT_ADD_EX(i_Linv,s);
            PLOT_DISP(i_Linv);
        }
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
            mat_set(phy_pt,i_umd,0,SQ(a)*DMSQ_K);
            mat_set(phy_pt,i_bind,0,bind);
            mat_set(phy_pt,i_a,0,0.0);
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
        if (param->s_M_ud_deg > 0)
        {
            sprintf(xlabel,"(a*M_%s)^2",param->ud_name);
            PLOT_DISP(i_ud);
        }
        if (param->s_M_s_deg > 0)
        {
            sprintf(xlabel,"(a*M_%s)^2",param->s_name);
            PLOT_DISP(i_s);
        }
        if (param->s_with_umd)
        {
            strbufcpy(xlabel,"(a*M_K_p)^2 - (a*M_K_0)^2");
            PLOT_DISP(i_umd);
        }
        if (param->s_with_qed_fvol)
        {
            strbufcpy(xlabel,"a/L");
            PLOT_DISP(i_Linv);
        }
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
static void print_result(const rs_sample *s_fit, fit_param *param)
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

enum 
{
    AINV = 0,
    A    = 1
};

#define SCALE_DATA(unit)\
{\
    int d_;\
    if (unit == AINV)\
    {\
        rs_sample_eqdivp(s_x[i_ud],s_x[i_a]);\
        rs_sample_eqdivp(s_x[i_ud],s_x[i_a]);\
        rs_sample_eqdivp(s_x[i_s],s_x[i_a]);\
        rs_sample_eqdivp(s_x[i_s],s_x[i_a]);\
        rs_sample_eqdivp(s_x[i_umd],s_x[i_a]);\
        rs_sample_eqdivp(s_x[i_umd],s_x[i_a]);\
        rs_sample_eqdivp(s_x[i_Linv],s_x[i_a]);\
        for (d_=0;d_<param->q_dim;d_++)\
        {\
            rs_sample_eqdivp(s_q[1],s_x[i_a]);\
        }\
    }\
    else if (unit == A)\
    {\
        rs_sample_eqmulp(s_x[i_ud],s_x[i_a]);\
        rs_sample_eqmulp(s_x[i_ud],s_x[i_a]);\
        rs_sample_eqmulp(s_x[i_s],s_x[i_a]);\
        rs_sample_eqmulp(s_x[i_s],s_x[i_a]);\
        rs_sample_eqmulp(s_x[i_umd],s_x[i_a]);\
        rs_sample_eqmulp(s_x[i_umd],s_x[i_a]);\
        rs_sample_eqmulp(s_x[i_Linv],s_x[i_a]);\
        for (d_=0;d_<param->q_dim;d_++)\
        {\
            rs_sample_eqmulp(s_q[1],s_x[i_a]);\
        }\
    }\
}

/* main program */
int main(int argc, char *argv[])
{
    /*              argument parsing            */
    /********************************************/
    fit_param *param;
    
    param = (fit_param *)malloc(sizeof(fit_param));
    
    if (argc != 2)
    {
        fprintf(stderr,"usage: %s <parameter file>\n",argv[0]);
        return EXIT_FAILURE;
    }
    parse_fit_param(param,argv[1]);
    
    /*              global settings             */
    /********************************************/
    latan_set_verb(param->verb);
    
    /*              data loading                */
    /********************************************/
    rs_sample *s_x[N_EX_VAR],*s_q[2];
    strbuf beta;
    size_t i;

    for (i=0;i<N_EX_VAR;i++)
    {
        s_x[i] = rs_sample_create(param->nens,param->nsample);
    }
    s_q[0] = rs_sample_create(param->nens,param->nsample);
    s_q[1] = rs_sample_create(param->nens,param->nsample);
    printf("-- loading data...\n");
    data_load(s_x,s_q,beta,param);
    printf("%-11s : %d\n","datasets",(int)param->ndataset);
    printf("%-11s : %d\n","points",(int)param->nens);
    printf("%-11s : %d\n","betas",(int)param->nbeta);
    printf("%-11s : %d\n","samples",(int)param->nsample);
    printf("\n");

    /*              data fitting                */
    /********************************************/
    fit_data *d;
    size_t npar,nydim,bind;
    rs_sample *s_fit,*s_tmp,**s_pt;
    mat *fit,*fit_var;
    strbuf resf_name, chi2f_name;
    bool use_x_var[N_EX_VAR] = {false,false,false,false,false,false};
    FILE *chi2f;
    
    nydim   = IS_ANALYZE(param,"comb_phypt_scale") ? 2 : 1;
    d       = fit_data_create(param->nens,N_EX_VAR,nydim);
    npar    = fit_model_get_npar(param->model,param);
    s_fit   = rs_sample_create(npar,param->nsample);
    s_pt    = IS_ANALYZE(param,"phypt") ? s_q+1 : s_q;
    s_tmp   = rs_sample_create(1,param->nsample);
    fit_var = mat_create(npar,1);
    
    if (IS_ANALYZE(param,"phypt")||IS_ANALYZE(param,"comb_phypt_scale"))
    {
        sprintf(chi2f_name,"%s_%s_%s.chi2",param->q_name,param->scale_part,\
                param->dataset_cat);
    }
    else if (IS_ANALYZE(param,"scaleset"))
    {
        sprintf(chi2f_name,"scale_%s_%s.chi2",param->scale_part,\
                param->dataset_cat);
    }
    chi2f = fopen(chi2f_name,"w");
    fit_data_fit_all_points(d,true);
    fit_data_set_model(d,param->model,param);
    mat_cst(rs_sample_pt_cent_val(s_fit),0.0001);
    for (i=0;i<param->ninit_param;i++)
    {
        mat_set(rs_sample_pt_cent_val(s_fit),param->init_param[i].ind,0,\
                param->init_param[i].value);
    }
    minimizer_set_alg(MIN_MIGRAD);
    fit = rs_sample_pt_cent_val(s_fit);
    printf("-- pre-fit...\n");
    rs_data_fit(s_fit,s_x,s_pt,d,NO_COR,use_x_var);
    printf("chi^2/dof = %e\n",fit_data_get_chi2pdof(d));
    fprintf(chi2f,"uncorrelated : %e\n",fit_data_get_chi2pdof(d));
    rs_sample_varp(fit_var,s_fit);
    if (IS_ANALYZE(param,"comb_phypt_scale"))
    {
        for (i=0;i<param->nens;i++)
        {
            bind = (size_t)(mat_get(rs_sample_pt_cent_val(s_x[i_bind]),i,0));
            rs_sample_get_subsamp(s_tmp,s_fit,bind,bind);
            rs_sample_set_subsamp(s_x[i_a],s_tmp,i,i);
        }
    }
    print_result(s_fit,param);
    use_x_var[i_ud]  = (param->M_ud_deg != 0)\
                       ||((param->s_M_ud_deg != 0)&&!IS_ANALYZE(param,"phypt"));
    use_x_var[i_s]   = (param->M_s_deg  != 0)\
                       ||((param->s_M_s_deg != 0)&&!IS_ANALYZE(param,"phypt"));
    use_x_var[i_umd] = (param->with_umd != 0)\
                       ||((param->s_with_umd != 0)&&!IS_ANALYZE(param,"phypt"));
    if (IS_ANALYZE(param,"phypt")||IS_ANALYZE(param,"comb_phypt_scale"))
    {
        SCALE_DATA(AINV);
        fit_data_set_covar_from_sample(d,s_x,s_pt,NO_COR,use_x_var);
        plot_fit(fit,fit_var,d,param,Q);
        SCALE_DATA(A);
        fit_data_set_covar_from_sample(d,s_x,s_pt,NO_COR,use_x_var);
    }
    if (IS_ANALYZE(param,"scaleset")||IS_ANALYZE(param,"comb_phypt_scale"))
    {
        plot_fit(fit,fit_var,d,param,SCALE);
    }
    printf("-- fitting and resampling %s...\n",param->q_name);
    rs_data_fit(s_fit,s_x,s_pt,d,X_COR|XDATA_COR|DATA_COR,use_x_var);
    rs_sample_varp(fit_var,s_fit);
    /*plot_fit(fit,d,param);*/
    
    /*              result output               */
    /********************************************/
    /** terminal/file output **/
    /*** chi^2 ***/
    printf("chi^2/dof = %e\n",fit_data_get_chi2pdof(d));
    fprintf(chi2f,"correlated   : %e\n",fit_data_get_chi2pdof(d));
    
    /*** parameters ***/
    print_result(s_fit,param);

    /*** extrapolation ***/
    if (IS_ANALYZE(param,"phypt"))
    {
        printf("extrapolation :\n");
        printf("%10s = %f +/- %e MeV^%d\n",param->q_name,mat_get(fit,0,0),\
               sqrt(mat_get(fit_var,0,0)),param->q_dim);
        sprintf(resf_name,"%s_%s%s.boot",param->q_name,param->scale_part,\
                param->dataset_cat);
        rs_sample_save_subsamp(resf_name,'w',s_fit,0,0);
    }
    else if (IS_ANALYZE(param,"scaleset"))
    {
        size_t j;
        
        sprintf(resf_name,"scale_%s_%s",param->scale_part,param->dataset_cat);
        rs_sample_set_name(s_fit,resf_name);
        rs_sample_eqmuls(s_fit,1.0/SQ(param->M_scale));
        rs_sample_eqsqrt(s_fit);
        printf("scales :\n\n");
        for (j=0;j<param->nbeta;j++)
        {
            rs_sample_varp(fit_var,s_fit);
            printf("beta = %s\n",param->beta[j]);
            printf("a    = %f +/- %e fm\n",mat_get(fit,j,0)/NU_FM,\
                   sqrt(mat_get(fit_var,j,0))/NU_FM);
            rs_sample_eqinvp(s_fit);
            rs_sample_varp(fit_var,s_fit);
            printf("a^-1 = %f +/- %e MeV\n",mat_get(fit,j,0),\
                   sqrt(mat_get(fit_var,j,0)));
            printf("\n");
            sprintf(resf_name,"scale_%s_%s_%s.boot",param->beta[j],\
                    param->scale_part,param->dataset_cat);
            rs_sample_save_subsamp(resf_name,'w',s_fit,j,j);
            rs_sample_eqinvp(s_fit);
        }
    }
    
    /*              desallocation               */
    /********************************************/
    for (i=0;i<N_EX_VAR;i++)
    {
        rs_sample_destroy(s_x[i]);
    }
    rs_sample_destroy(s_q[0]);
    rs_sample_destroy(s_q[1]);
    fit_data_destroy(d);
    rs_sample_destroy(s_fit);
    rs_sample_destroy(s_tmp);
    mat_destroy(fit_var);
    fclose(chi2f);
    free(param->beta);
    free(param->init_param);
    free(param);
    
    return EXIT_SUCCESS;
}
