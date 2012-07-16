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
            mat_eqdivp(res[s],rs_sample_pt_cent_val(s_x[i_a]));\
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
            mat_eqmulp(res[s],rs_sample_pt_cent_val(s_x[i_a]));\
        }\
    }\
}

/* main program */
int main(int argc, char *argv[])
{
    /*               I/O init                   */
    /********************************************/
    io_init();
    
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
        s_x[i] = rs_sample_create(param->nens,1,param->nsample);
    }
    s_q[0] = rs_sample_create(param->nens,1,param->nsample);
    s_q[1] = rs_sample_create(param->nens,1,param->nsample);
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
    size_t npar,nydim,bind,s,s_par;
    size_t j;
    rs_sample *s_fit,*s_tmp,**s_pt;
    mat *fit,*fit_var,**res;
    strbuf chi2f_name;
    bool use_x_var[N_EX_VAR] = {false,false,false,false,false,false,false};
    FILE *chi2f,*tablef;
    
    
    nydim   = (IS_AN(param,AN_PHYPT)&&IS_AN(param,AN_SCALE)) ? 2 : 1;
    d       = fit_data_create(param->nens,N_EX_VAR,nydim);
    npar    = fit_model_get_npar(&param->fm,param);
    s_fit   = rs_sample_create(npar,1,param->nsample);
    s       = (IS_AN(param,AN_PHYPT)&&IS_AN(param,AN_SCALE)) ? 1 : 0;
    s_pt    = (IS_AN(param,AN_PHYPT)&&!IS_AN(param,AN_SCALE)) ? s_q + 1 : s_q;
    s_tmp   = rs_sample_create(1,1,param->nsample);
    fit_var = mat_create(npar,1);
    res     = mat_ar_create(2,param->nens,1);
    
    
    if (IS_AN(param,AN_PHYPT))
    {
        sprintf(chi2f_name,"%s_%s_%s_%s_%s.chi2",param->q_name,        \
                param->scale_part,param->s_manifest,param->dataset_cat,\
                param->manifest);
    }
    else if (IS_AN(param,AN_SCALE))
    {
        sprintf(chi2f_name,"a_%s_%s.chi2",param->scale_part,param->s_manifest);
    }
    fit_data_fit_all_points(d,true);
    fit_data_set_model(d,&param->fm,param);
    mat_cst(rs_sample_pt_cent_val(s_fit),0.00001);
    if (IS_AN(param,AN_PHYPT)&&IS_AN(param,AN_SCALE))
    {
        s_par = fm_scaleset_taylor_npar(param);
    }
    else
    {
        s_par = 0;
    }
    if (IS_AN(param,AN_PHYPT)&&!IS_AN(param,AN_SCALE)&&(param->with_ext_a))
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
        if (strbufcmp(param->point[i].dir,param->point[j].dir) != 0)
        {
            fit_data_set_data_cor(d,i,j,false);
        }
    }
    
    /* uncorrelated fit without x errors */
    if (param->correlated)
    {
        printf("-- pre-fit...\n");
    }
    else
    {
        printf("-- fitting and resampling %s...\n",param->q_name);
    }
    rs_data_fit(s_fit,NULL,s_x,s_pt,d,NO_COR,use_x_var);
    if (IS_AN(param,AN_PHYPT))
    {
        fit_residual(res[s],d,0,fit);
    }
    if (IS_AN(param,AN_SCALE))
    {
        fit_residual(res[0],d,0,fit);
    }
    /** save chi^2/dof **/
    printf("chi^2/dof = %e\n",fit_data_get_chi2pdof(d));
    if (param->save_result)
    {
        chi2f = fopen(chi2f_name,"w");
        fprintf(chi2f,"uncorrelated: %e %e %d\n",fit_data_get_chi2pdof(d),\
                fit_data_get_chi2(d),(int)fit_data_get_dof(d));
        fclose(chi2f);
    }
    /** compute errors **/
    rs_sample_varp(fit_var,s_fit);
    /** save lattice spacings **/
    if (IS_AN(param,AN_PHYPT))
    {
        if (IS_AN(param,AN_SCALE))
        {
            for (i=0;i<param->nens;i++)
            {
                bind = (size_t)(mat_get(rs_sample_pt_cent_val(s_x[i_bind]),\
                                        i,0));
                rs_sample_get_subsamp(s_tmp,s_fit,bind,0,bind,0);
                rs_sample_set_subsamp(s_x[i_a],s_tmp,i,0,i,0);
            }
        }
        else if (param->with_ext_a)
        {
            for (i=0;i<param->nens;i++)
            {
                bind = (size_t)(mat_get(rs_sample_pt_cent_val(s_x[i_bind]),\
                                        i,0));
                rs_sample_get_subsamp(s_tmp,s_fit,npar-param->nbeta+bind,0,\
                                      npar-param->nbeta+bind,0);
                rs_sample_set_subsamp(s_x[i_a],s_tmp,i,0,i,0);
            }
        }
    }
    /** save tables **/
    if (IS_AN(param,AN_PHYPT))
    {
        tablef = fopen("qcd_phyfit_q.dat","w");
        SCALE_DATA(AINV);
        fprint_table(tablef,s_x,s_q,res[s],param,Q);
        SCALE_DATA(A);
        fclose(tablef);
    }
    if (IS_AN(param,AN_SCALE))
    {
        tablef = fopen("qcd_phyfit_scale.dat","w");
        fprint_table(tablef,s_x,s_q,res[0],param,SCALE);
        fclose(tablef);
    }
    /** print results **/   
    print_result(s_fit,param);
    /** display plots **/
    if (param->plot)
    {
        if (IS_AN(param,AN_PHYPT))
        {
            plot_chi2_comp(d,param,s,"physical point fit");
            SCALE_DATA(AINV);
            fit_data_set_covar_from_sample(d,s_x,s_pt,NO_COR,use_x_var);
            plot_fit(fit,fit_var,d,param,Q);
            SCALE_DATA(A);
        }
        if (IS_AN(param,AN_SCALE))
        {
            plot_chi2_comp(d,param,0,"scale setting");
            plot_fit(fit,fit_var,d,param,SCALE);
        }
    }
    
    /* real fit */
    if (param->correlated)
    {
        use_x_var[i_ud]  = (((param->M_ud_deg != 0)||(param->with_udumd)       \
                             ||(param->with_udalpha))&&IS_AN(param,AN_PHYPT))  \
                           ||((param->s_M_ud_deg != 0)&&IS_AN(param,AN_SCALE));
        use_x_var[i_s]   = (((param->M_s_deg != 0)||(param->with_sumd)         \
                             ||(param->with_salpha))&&IS_AN(param,AN_PHYPT))   \
                           ||((param->s_M_ud_deg != 0)&&IS_AN(param,AN_SCALE));
        use_x_var[i_umd] = (((param->umd_deg != 0)&&IS_AN(param,AN_PHYPT))     \
                           ||((param->s_umd_deg != 0)&&IS_AN(param,AN_SCALE))) \
                           &&param->have_umd;
        printf("-- fitting and resampling %s...\n",param->q_name);
        rs_data_fit(s_fit,NULL,s_x,s_pt,d,X_COR|XDATA_COR|DATA_COR,use_x_var);
        if (IS_AN(param,AN_PHYPT))
        {
            fit_residual(res[s],d,0,fit);
        }
        if (IS_AN(param,AN_SCALE))
        {
            fit_residual(res[0],d,0,fit);
        }
        printf("chi^2/dof = %e\n",fit_data_get_chi2pdof(d));
        /** save chi^2/dof **/
        if (param->save_result)
        {
            chi2f = fopen(chi2f_name,"a");
            fprintf(chi2f,"correlated: %e %e %d\n",fit_data_get_chi2pdof(d),\
                    fit_data_get_chi2(d),(int)fit_data_get_dof(d));
            fclose(chi2f);
        }
        /** compute errors **/
        rs_sample_varp(fit_var,s_fit);
        /** save lattice spacings **/
        if (IS_AN(param,AN_PHYPT))
        {
            if (IS_AN(param,AN_SCALE))
            {
                for (i=0;i<param->nens;i++)
                {
                    bind = (size_t)(mat_get(rs_sample_pt_cent_val(s_x[i_bind]),\
                                            i,0));
                    rs_sample_get_subsamp(s_tmp,s_fit,bind,0,bind,0);
                    rs_sample_set_subsamp(s_x[i_a],s_tmp,i,0,i,0);
                }
            }
            else if (param->with_ext_a)
            {
                for (i=0;i<param->nens;i++)
                {
                    bind = (size_t)(mat_get(rs_sample_pt_cent_val(s_x[i_bind]),\
                                            i,0));
                    rs_sample_get_subsamp(s_tmp,s_fit,npar-param->nbeta+bind,0,\
                                          npar-param->nbeta+bind,0);
                    rs_sample_set_subsamp(s_x[i_a],s_tmp,i,0,i,0);
                }
            }
        }
        /** save tables **/
        if (IS_AN(param,AN_PHYPT))
        {
            tablef = fopen("qcd_phyfit_q.dat","w");
            SCALE_DATA(AINV);
            fprint_table(tablef,s_x,s_q,res[s],param,Q);
            SCALE_DATA(A);
            fclose(tablef);
        }
        if (IS_AN(param,AN_SCALE))
        {
            tablef = fopen("qcd_phyfit_scale.dat","w");
            fprint_table(tablef,s_x,s_q,res[0],param,SCALE);
            fclose(tablef);
        }
        /** print results **/
        print_result(s_fit,param);
        /** display plots **/
        if (param->plot)
        {
            if (IS_AN(param,AN_PHYPT))
            {
                SCALE_DATA(AINV);
                fit_data_set_covar_from_sample(d,s_x,s_pt,NO_COR,use_x_var);
                plot_fit(fit,fit_var,d,param,Q);
                SCALE_DATA(A);
                fit_data_set_covar_from_sample(d,s_x,s_pt,NO_COR,use_x_var);
            }
            if (IS_AN(param,AN_SCALE))
            {
                plot_fit(fit,fit_var,d,param,SCALE);
            }
        }
    }
    
    /*              result output               */
    /********************************************/
    strbuf res_path;
    char mode;
    
    if (IS_AN(param,AN_PHYPT))
    {
        printf("extrapolation :\n");
        printf("%10s = %f +/- %e MeV^%d\n",param->q_name,mat_get(fit,s_par,0),\
               sqrt(mat_get(fit_var,s_par,0)),param->q_dim);
        sprintf(res_path,"%s_%s_%s_%s_%s.boot%c%s_%s_%s_%s_%s",param->q_name,\
                param->scale_part,param->s_manifest,param->dataset_cat,      \
                param->manifest,LATAN_PATH_SEP,param->q_name,                \
                param->scale_part,param->s_manifest,param->dataset_cat,      \
                param->manifest);
        if (param->save_result)
        {
            rs_sample_save_subsamp(res_path,'w',s_fit,s_par,0, \
                                   rs_sample_get_nrow(s_fit)-1,0);
        }
    }
    if (IS_AN(param,AN_SCALE))
    {
        printf("scales :\n\n");
        for (i=0;i<param->nbeta;i++)
        {
            if (param->save_result)
            {
                sprintf(res_path,"a_%s_%s.boot%ca_%s_%s_%s",              \
                        param->scale_part,param->s_manifest,LATAN_PATH_SEP,\
                        param->beta[i],param->scale_part,param->s_manifest);
                mode = (i == 0) ? 'w' : 'a';
                rs_sample_save_subsamp(res_path,mode,s_fit,i,0,i,0);
            }
            rs_sample_varp(fit_var,s_fit);
            printf("beta = %s\n",param->beta[i]);
            printf("a    = %f +/- %e fm\n",mat_get(fit,i,0)/NU_FM,\
                   sqrt(mat_get(fit_var,i,0))/NU_FM);
            rs_sample_eqinvp(s_fit);
            rs_sample_varp(fit_var,s_fit);
            printf("a^-1 = %f +/- %e MeV\n",mat_get(fit,i,0),\
                   sqrt(mat_get(fit_var,i,0)));
            printf("\n");
            rs_sample_eqinvp(s_fit);
        }
    }
    
    /*               I/O finish                 */
    /********************************************/
    io_finish();
    
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
    mat_ar_destroy(res,2);
    mat_destroy(fit_var);
    fit_param_destroy(param);
    
    return EXIT_SUCCESS;
}
