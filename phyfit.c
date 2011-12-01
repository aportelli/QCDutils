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

void plot_fit(const mat *fit, fit_data *d, fit_param *param);
void print_result(const rs_sample *s_fit, fit_param *param);

#define ADD_PLOT(obj,title,color)\
if (param->M_ud_deg > 0)\
{\
    plot_add_fit(p[i_ud],d,i_ud,phy_pt,fit,true,obj,title,"",color,color);\
}\
if (param->M_s_deg > 0)\
{\
    plot_add_fit(p[i_s],d,i_s,phy_pt,fit,true,obj,title,"",color,color);\
}\
if (param->a_deg > 0)\
{\
    plot_add_fit(p[i_ainv],d,i_ainv,phy_pt,fit,true,obj,title,"",color,color);\
}\
if (param->with_umd)\
{\
    plot_add_fit(p[i_umd],d,i_umd,phy_pt,fit,true,obj,title,"",color,color);\
}\
if (param->with_qed_fvol)\
{\
    plot_add_fit(p[i_Linv],d,i_Linv,phy_pt,fit,true,obj,title,"",color,color);\
}

void plot_fit(const mat *fit, fit_data *d, fit_param *param)
{
    plot *p[N_EX_VAR];
    double *xb[N_EX_VAR] = {NULL,NULL,NULL,NULL,NULL,NULL};
    double b_int[2],dbind,M_scale;
    size_t bind,k;
    strbuf color,title,xlabel,ylabel;
    mat *phy_pt;
    
    phy_pt = mat_create(N_EX_VAR,1);
    for (k=0;k<N_EX_VAR;k++)
    {
        p[k] = plot_create();
    }
    
    mat_set(phy_pt,i_bind,0,0.0);
    mat_set(phy_pt,i_ainv,0,0.0);
    mat_set(phy_pt,i_umd,0,DMSQ_K);
    mat_set(phy_pt,i_Linv,0,0.0);
    if (IS_ANALYZE(param,"phypt"))
    {
        mat_set(phy_pt,i_ud,0,SQ(param->M_ud));
        mat_set(phy_pt,i_s,0,SQ(param->M_s));
        for (bind=0;bind<param->nbeta;bind++)
        {
            dbind      = (double)(bind);
            b_int[0]   = dbind - 0.1;
            b_int[1]   = dbind + 0.1;
            xb[i_bind] = b_int;
            fit_data_fit_region(d,xb);
            sprintf(color,"%d",1+(int)bind);
            sprintf(title,"beta = %s",param->beta[bind]);
            ADD_PLOT(PF_DATA,title,color);
            fit_data_fit_all_points(d,true);
        }
        ADD_PLOT(PF_FIT,"","rgb 'black'");
        switch (param->q_dim) 
        {
            case 0:
                strbufcpy(ylabel,param->q_name);
                break;
            case 1:
                sprintf(ylabel,"%s (MeV)",param->q_name);
                break;
            default:
                sprintf(xlabel,"%s^%d (MeV^%d)",param->q_name,param->q_dim,\
                        param->q_dim);
                break;
        }
    }
    else if (IS_ANALYZE(param,"scaleset"))
    {
        for (bind=0;bind<param->nbeta;bind++)
        {
            dbind      = (double)(bind);
            b_int[0]   = dbind - 0.1;
            b_int[1]   = dbind + 0.1;
            xb[i_bind] = b_int;
            M_scale    = param->M_scale;
            fit_data_fit_region(d,xb);
            mat_set(phy_pt,i_ud,0,\
                    SQ(param->M_ud)*mat_get(fit,bind,0)/SQ(M_scale));
            mat_set(phy_pt,i_s,0,\
                    SQ(param->M_s)*mat_get(fit,bind,0)/SQ(M_scale));
            mat_set(phy_pt,i_bind,0,bind);
            sprintf(color,"%d",1+(int)bind);
            sprintf(title,"beta = %s",param->beta[bind]);
            ADD_PLOT(PF_DATA|PF_FIT,title,color);
            fit_data_fit_all_points(d,true);
        }
        sprintf(ylabel,"a*%s",param->q_name);
    }
    if (param->M_ud_deg > 0)
    {
        if (IS_ANALYZE(param,"phypt"))
        {
            sprintf(xlabel,"M_%s^2 (MeV^2)",param->ud_name);
        }
        else if (IS_ANALYZE(param,"scaleset"))
        {
            sprintf(xlabel,"(a*M_%s)^2",param->ud_name);
        }
        plot_set_xlabel(p[i_ud],xlabel);
        plot_set_ylabel(p[i_ud],ylabel);
        plot_disp(p[i_ud]);

    }
    if (param->M_s_deg > 0)
    {
        if (IS_ANALYZE(param,"phypt"))
        {
            sprintf(xlabel,"M_%s^2 (MeV^2)",param->s_name);
        }
        else if (IS_ANALYZE(param,"scaleset"))
        {
            sprintf(xlabel,"(a*M_%s)^2",param->s_name);
        }
        plot_set_xlabel(p[i_s],xlabel);
        plot_set_ylabel(p[i_s],ylabel);
        plot_disp(p[i_s]);
    }
    if (param->a_deg > 0)
    {
        sprintf(xlabel,"a^-%d (MeV)",param->a_deg);
        plot_set_xlabel(p[i_ainv],xlabel);
        plot_set_ylabel(p[i_ainv],ylabel);
        plot_disp(p[i_ainv]);
    }
    if (param->with_umd)
    {
        strbufcpy(xlabel,"m_u-m_d(PCAC) (MeV)");
        plot_set_xlabel(p[i_umd],xlabel);
        plot_set_ylabel(p[i_umd],ylabel);
        plot_disp(p[i_umd]);
    }
    if (param->with_qed_fvol)
    {
        strbufcpy(xlabel,"L^-1 (MeV)");
        plot_set_xlabel(p[i_Linv],xlabel);
        plot_set_ylabel(p[i_Linv],ylabel);
        plot_disp(p[i_Linv]);
    }
    mat_destroy(phy_pt);
    for (k=0;k<N_EX_VAR;k++)
    {
        plot_destroy(p[k]);
    }
}
#undef ADD_PLOT
#undef PLOT_DISP

#define PRINT_PAR(name)\
{\
    printf("%20s = %e (%.4f%%)\n",name,mat_get(fit,i,0),            \
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
    printf("\nfit parameters (pi:%d K:%d a:%d m_u-m_d:%d) :\n",     \
           param->M_ud_deg,param->M_s_deg,param->a_deg,param->with_umd);
    if (IS_ANALYZE(param,"phypt"))
    {
        PRINT_PAR(param->q_name);
    }
    else if (IS_ANALYZE(param,"scaleset"))
    {
        for(j=0;j<(int)param->nbeta;j++)
        {
            sprintf(buf,"aM_%s_%s",param->scale_part,param->beta[i]);
            PRINT_PAR(buf);
        }
    }
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
    printf("\n");
    
    mat_destroy(fit_var);
}
#undef PRINT_PAR
                      
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
    rs_sample *s_x[N_EX_VAR],*s_q;
    strbuf fset,sfname,beta;
    size_t nset,nsample;
    size_t i;
    
    nset = (size_t)get_nfile(param->manifest);

    get_firstfname(fset,param->manifest);
    sprintf(sfname,"%s/%s_%s.boot%s",fset,param->q_name,param->dataset,\
            (io_get_fmt() == IO_XML) ? ".xml" : "");
    rs_sample_load_nsample(&nsample,sfname,"");
    for (i=0;i<N_EX_VAR;i++)
    {
        s_x[i] = rs_sample_create(nset,nsample);
    }
    s_q = rs_sample_create(nset,nsample);
    printf("-- loading data...\n");
    data_load(s_x,s_q,beta,param);
    

    /*              data fitting                */
    /********************************************/
    fit_data *d;
    size_t npar;
    rs_sample *s_fit;
    mat *fit,*fit_var;
    strbuf resf_name, chi2f_name;
    bool use_x_var[N_EX_VAR] = {false,false,false,false,false,false};
    FILE *chi2f;
    
    d       = fit_data_create(nset,N_EX_VAR);
    npar    = fit_model_get_npar(param->model,param);
    s_fit   = rs_sample_create(npar,(size_t)nsample);
    fit_var = mat_create(npar,1);
    
    if (IS_ANALYZE(param,"phypt"))
    {
        sprintf(chi2f_name,"%s_%s_%s.chi2",param->q_name,param->scale_part,\
                param->dataset);
    }
    else if (IS_ANALYZE(param,"scaleset"))
    {
        sprintf(chi2f_name,"scale_%s_%s.chi2",param->scale_part,param->dataset);
    }
    chi2f = fopen(chi2f_name,"w");
    fit_data_fit_all_points(d,true);
    fit_data_set_model(d,param->model,param);
    mat_cst(rs_sample_pt_cent_val(s_fit),1.0);
    
    for (i=0;i<param->ninit_param;i++)
    {
        mat_set(rs_sample_pt_cent_val(s_fit),param->init_param[i].ind,0,\
                param->init_param[i].value);
    }
    minimizer_set_alg(MIN_MIGRAD);
    fit = rs_sample_pt_cent_val(s_fit);
    printf("-- pre-fit...\n");
    rs_x_data_fit(s_fit,s_x,s_q,d,NO_COR,use_x_var);
    printf("chi^2/dof = %e\n",fit_data_get_chi2pdof(d));
    fprintf(chi2f,"uncorrelated : %e\n",fit_data_get_chi2pdof(d));
    rs_sample_varp(fit_var,s_fit);
    print_result(s_fit,param);
    plot_fit(fit,d,param);
    printf("-- fitting and resampling %s...\n",param->q_name);
    use_x_var[i_ud]   = (param->M_ud_deg != 0);
    use_x_var[i_s]    = (param->M_s_deg  != 0);
    use_x_var[i_umd]  = (param->with_umd != 0);
    use_x_var[i_ainv] = (IS_ANALYZE(param,"phypt"));
    rs_x_data_fit(s_fit,s_x,s_q,d,X_COR|XDATA_COR,use_x_var);
    rs_sample_varp(fit_var,s_fit);
    plot_fit(fit,d,param);
    
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
                param->dataset);
        rs_sample_save_subsamp(resf_name,'w',s_fit,0,0);
    }
    else if (IS_ANALYZE(param,"scaleset"))
    {
        size_t j;
        
        sprintf(resf_name,"scale_%s_%s",param->scale_part,param->dataset);
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
                    param->scale_part,param->dataset);
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
    rs_sample_destroy(s_q);
    fit_data_destroy(d);
    rs_sample_destroy(s_fit);
    mat_destroy(fit_var);
    fclose(chi2f);
    free(param->beta);
    free(param->init_param);
    free(param);
    
    return EXIT_SUCCESS;
}
