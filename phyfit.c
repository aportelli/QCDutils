#include "data_loader.h"
#include "models.h"
#include "output.h"
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
            mat_eqdivp(res[1],rs_sample_pt_cent_val(s_x[i_a]));\
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
            mat_eqmulp(res[1],rs_sample_pt_cent_val(s_x[i_a]));\
        }\
    }\
}

/* main program */
int main(int argc, char *argv[])
{
    /*              argument parsing            */
    /********************************************/
    fit_param *param;
    
    if (argc != 2)
    {
        fprintf(stderr,"usage: %s <parameter file>\n",argv[0]);
        return EXIT_FAILURE;
    }
    param = fit_param_parse(argv[1]);
    
    /*              global settings             */
    /********************************************/
    latan_set_verb(param->verb);
    
    /*              data loading                */
    /********************************************/
    rs_sample *s_x[N_EX_VAR],*s_q[2];
    size_t i;

    for (i=0;i<N_EX_VAR;i++)
    {
        s_x[i] = rs_sample_create(param->nens,param->nsample);
    }
    s_q[0] = rs_sample_create(param->nens,param->nsample);
    s_q[1] = rs_sample_create(param->nens,param->nsample);
    printf("-- loading data...\n");
    data_load(s_x,s_q,param);
    printf("%-11s : %d\n","datasets",(int)param->ndataset);
    printf("%-11s : %d\n","points",(int)param->nens);
    printf("%-11s : %d\n","betas",(int)param->nbeta);
    printf("%-11s : %d\n","samples",(int)param->nsample);
    printf("\n");

    /*              data fitting                */
    /********************************************/
    fit_data *d;
    size_t npar,nydim,bind,s;
    size_t j;
    rs_sample *s_fit,*s_tmp,**s_pt;
    mat *fit,*fit_var,**res;
    strbuf resf_name, chi2f_name;
    bool use_x_var[N_EX_VAR] = {false,false,false,false,false,false};
    FILE *chi2f,*tablef;
    
    nydim   = IS_ANALYZE(param,"comb_phypt_scale") ? 2 : 1;
    d       = fit_data_create(param->nens,N_EX_VAR,nydim);
    npar    = fit_model_get_npar(param->model,param);
    s_fit   = rs_sample_create(npar,param->nsample);
    s_pt    = IS_ANALYZE(param,"phypt") ? s_q + 1 : s_q;
    s_tmp   = rs_sample_create(1,param->nsample);
    fit_var = mat_create(npar,1);
    res     = mat_ar_create(2,param->nens,1);
    
    if (IS_ANALYZE(param,"phypt")||IS_ANALYZE(param,"comb_phypt_scale"))
    {
        sprintf(chi2f_name,"%s_%s_%s.chi2",param->q_name,param->scale_part,\
                param->dataset_cat);
    }
    else if (IS_ANALYZE(param,"scaleset"))
    {
        sprintf(chi2f_name,"a_%s_%s.chi2",param->scale_part,\
                param->dataset_cat);
    }
    chi2f = fopen(chi2f_name,"w");
    fit_data_fit_all_points(d,true);
    fit_data_set_model(d,param->model,param);
    mat_cst(rs_sample_pt_cent_val(s_fit),0.0001);
    if (IS_ANALYZE(param,"comb_phypt_scale"))
    {
        s = fm_scaleset_taylor_npar(param);
    }
    else
    {
        s = 0;
    }
    if (IS_ANALYZE(param,"phypt")||IS_ANALYZE(param,"comb_phypt_scale"))
    {
        mat_set(rs_sample_pt_cent_val(s_fit),s,0,param->q_target[0]);
    }
    if (IS_ANALYZE(param,"phypt")&&(param->with_ext_a))
    {
        fit_data_set_chi2_ext(d,&a_error_chi2_ext);
        mat_set_subm(rs_sample_pt_cent_val(s_fit),   \
                     rs_sample_pt_cent_val(param->a),\
                     npar-param->nbeta,0,npar-1,0);
        fit_data_set_ndumbpar(d,param->nbeta);
    }
    for (i=0;i<param->ninit_param;i++)
    {
        mat_set(rs_sample_pt_cent_val(s_fit),param->init_param[i].ind,0,\
                param->init_param[i].value);
    }
    minimizer_set_alg(MIN_MIGRAD);
    fit = rs_sample_pt_cent_val(s_fit);
    for (i=0;i<param->nens;i++)
    for (j=i+1;j<param->nens;j++)
    {
        if (strcmp(param->point[i].dir,param->point[j].dir) != 0)
        {
            fit_data_set_data_cor(d,i,j,false);
        }
    }
    
    /* uncorrelated fit without x errors to find a good initial value */
    printf("-- pre-fit...\n");
    rs_data_fit(s_fit,s_x,s_pt,d,NO_COR,use_x_var);
    if (IS_ANALYZE(param,"phypt"))
    {
        fit_residual(res[1],d,0,fit);
    }
    else if (IS_ANALYZE(param,"scaleset"))
    {
        fit_residual(res[0],d,0,fit);
    }
    else if (IS_ANALYZE(param,"comb_phypt_scale"))
    {
        fit_residual(res[0],d,0,fit);
        fit_residual(res[1],d,1,fit);
    }
    /** save chi^2/dof **/
    printf("chi^2/dof = %e\n",fit_data_get_chi2pdof(d));
    fprintf(chi2f,"uncorrelated : %e\n",fit_data_get_chi2pdof(d));
    /** compute errors **/
    rs_sample_varp(fit_var,s_fit);
    /** save lattice spacings **/
    if (IS_ANALYZE(param,"comb_phypt_scale"))
    {
        for (i=0;i<param->nens;i++)
        {
            bind = (size_t)(mat_get(rs_sample_pt_cent_val(s_x[i_bind]),i,0));
            rs_sample_get_subsamp(s_tmp,s_fit,bind,bind);
            rs_sample_set_subsamp(s_x[i_a],s_tmp,i,i);
        }
    }
    else if (IS_ANALYZE(param,"phypt")&&(param->with_ext_a))
    {
        for (i=0;i<param->nens;i++)
        {
            bind = (size_t)(mat_get(rs_sample_pt_cent_val(s_x[i_bind]),i,0));
            rs_sample_get_subsamp(s_tmp,s_fit,npar-param->nbeta+bind,\
                                  npar-param->nbeta+bind);
            rs_sample_set_subsamp(s_x[i_a],s_tmp,i,i);
        }
    }
    /** save tables **/
    if (IS_ANALYZE(param,"phypt")||IS_ANALYZE(param,"comb_phypt_scale"))
    {
        tablef = fopen("qcd_phyfit_q.dat","w");
        SCALE_DATA(AINV);
        fprint_table(tablef,s_x,s_q,res[1],param,Q);
        SCALE_DATA(A);
        fclose(tablef);
    }
    if (IS_ANALYZE(param,"scaleset")||IS_ANALYZE(param,"comb_phypt_scale"))
    {
        tablef = fopen("qcd_phyfit_scale.dat","w");
        fprint_table(tablef,s_x,s_q,res[0],param,SCALE);
        fclose(tablef);
    }
    /** print results **/   
    print_result(s_fit,param);
    if (IS_ANALYZE(param,"comb_phypt_scale"))
    {
        plot_chi2_comp(d,param,0,"scale setting");
        plot_chi2_comp(d,param,1,"physical point fit");
    }
    /** display plots **/
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
    
    /* real fit */
    printf("-- fitting and resampling %s...\n",param->q_name);
    rs_data_fit(s_fit,s_x,s_pt,d,X_COR|XDATA_COR|DATA_COR,use_x_var);
    /** save chi^2/dof **/
    printf("chi^2/dof = %e\n",fit_data_get_chi2pdof(d));
    fprintf(chi2f,"correlated   : %e\n",fit_data_get_chi2pdof(d));
    /** compute errors **/
    rs_sample_varp(fit_var,s_fit);
    /** save lattice spacings in 'comb_phypt_scale' **/
    if (IS_ANALYZE(param,"comb_phypt_scale"))
    {
        for (i=0;i<param->nens;i++)
        {
            bind = (size_t)(mat_get(rs_sample_pt_cent_val(s_x[i_bind]),i,0));
            rs_sample_get_subsamp(s_tmp,s_fit,bind,bind);
            rs_sample_set_subsamp(s_x[i_a],s_tmp,i,i);
        }
    }
    /** display plots **/
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
    
    /*              result output               */
    /********************************************/
        
    /* parameters */
    print_result(s_fit,param);

    /* extrapolation */
    if (IS_ANALYZE(param,"phypt")||IS_ANALYZE(param,"comb_phypt_scale"))
    {
        printf("extrapolation :\n");
        printf("%10s = %f +/- %e MeV^%d\n",param->q_name,mat_get(fit,s,0),\
               sqrt(mat_get(fit_var,0,0)),param->q_dim);
        sprintf(resf_name,"%s_%s_%s.boot",param->q_name,param->scale_part,\
                param->dataset_cat);
        rs_sample_save_subsamp(resf_name,'w',s_fit,s,s);
    }
    else if (IS_ANALYZE(param,"scaleset")||IS_ANALYZE(param,"comb_phypt_scale"))
    {
        size_t k;
        
        printf("scales :\n\n");
        for (k=0;k<param->nbeta;k++)
        {
            sprintf(resf_name,"a_%s_%s_%s.boot",param->beta[k],\
                    param->scale_part,param->dataset_cat);
            rs_sample_save_subsamp(resf_name,'w',s_fit,k,k);
            rs_sample_varp(fit_var,s_fit);
            printf("beta = %s\n",param->beta[k]);
            printf("a    = %f +/- %e fm\n",mat_get(fit,k,0)/NU_FM,\
                   sqrt(mat_get(fit_var,k,0))/NU_FM);
            rs_sample_eqinvp(s_fit);
            rs_sample_varp(fit_var,s_fit);
            printf("a^-1 = %f +/- %e MeV\n",mat_get(fit,k,0),\
                   sqrt(mat_get(fit_var,k,0)));
            printf("\n");
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
    fit_param_destroy(param);
    
    return EXIT_SUCCESS;
}
