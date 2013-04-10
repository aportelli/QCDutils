#include "config.h"
#include "alpha_s.h"
#include "data_loader.h"
#include "models.h"
#include "output.h"
#include "parameters.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#ifdef HAVE_MPI
#include <mpi.h>
#endif
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
        rs_sample_eqdivp(s_x[i_fvM],s_x[i_a]);\
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
        rs_sample_eqmulp(s_x[i_fvM],s_x[i_a]);\
        for (d_=0;d_<param->q_dim;d_++)\
        {\
            rs_sample_eqmulp(s_q[1],s_x[i_a]);\
            mat_eqmulp(res[s],rs_sample_pt_cent_val(s_x[i_a]));\
        }\
    }\
}

/* analysis routine */
static void analysis(fit_param *param)
{
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
    mpi_printf("-- loading data...\n");
    data_load(s_x,s_q,param);
    mpi_printf("%-11s : %d\n","datasets",(int)param->ndataset);
    mpi_printf("%-11s : %d\n","points",(int)param->nens);
    mpi_printf("%-11s : %d\n","betas",(int)param->nbeta);
    mpi_printf("%-11s : %d\n","samples",(int)param->nsample);
    mpi_printf("\n");
    
    /*              data fitting                */
    /********************************************/
    fit_data *d;
    size_t npar,nydim,bind,s;
    size_t j;
    rs_sample *s_fit,*s_tmp,**s_pt;
    mat *fit,*fit_var,*limit,**res,*phy_pt,*ex,*ex_err;
    strbuf chi2f_name;
    bool use_x_var[N_EX_VAR] = {false,false,false,false,false,false,false};
    FILE *chi2f,*tablef;
    double buf;
    
    nydim   = (IS_AN(param,AN_PHYPT)&&IS_AN(param,AN_SCALE)) ? 2 : 1;
    d       = fit_data_create(param->nens,N_EX_VAR,nydim);
    npar    = fit_model_get_npar(&param->fm,param);
    s_fit   = rs_sample_create(npar,1,param->nsample);
    fit     = rs_sample_pt_cent_val(s_fit);
    s       = (IS_AN(param,AN_PHYPT)&&IS_AN(param,AN_SCALE)) ? 1 : 0;
    s_pt    = (IS_AN(param,AN_PHYPT)&&!IS_AN(param,AN_SCALE)) ? s_q + 1 : s_q;
    s_tmp   = rs_sample_create(1,1,param->nsample);
    ex      = rs_sample_pt_cent_val(param->s_ex);
    ex_err  = param->ex_err;
    fit_var = mat_create(npar,1);
    limit   = mat_create(npar,2);
    phy_pt  = mat_create(N_EX_VAR,1);
    res     = mat_ar_create(2,param->nens,1);
    

    sprintf(chi2f_name,"%s.chi2",param->result_file);
    fit_data_fit_all_points(d,true);
    fit_data_set_model(d,&param->fm,param);
    mat_cst(rs_sample_pt_cent_val(s_fit),1.0e-15);
    mat_cst(limit,latan_nan());
    if (IS_AN(param,AN_PHYPT)&&!IS_AN(param,AN_SCALE)&&\
        (strbufcmp(param->with_ext_a,"") != 0))
    {
        fit_data_set_chi2_ext(d,&a_error_chi2_ext);
        mat_set_subm(rs_sample_pt_cent_val(s_fit),   \
                     rs_sample_pt_cent_val(param->s_a),\
                     npar-param->nbeta,0,npar-1,0);
        fit_data_set_ndumbpar(d,param->nbeta);
    }
    if (IS_AN(param,AN_PHYPT))
    {
        for (i=npar-param->nbeta;i<npar;i++)
        {
            mat_set(limit,i,0,0.0);
        }
    }
    for (i=0;i<param->ninit_param;i++)
    {
        mat_set(rs_sample_pt_cent_val(s_fit),param->init_param[i].ind,0,\
                param->init_param[i].value);
    }
    for (i=0;i<param->nlimit_param;i++)
    {
        mat_set(limit,param->limit_param[i].ind,0,\
                param->limit_param[i].value[0]);
        mat_set(limit,param->limit_param[i].ind,1,\
                param->limit_param[i].value[1]);
    }
    minimizer_set_alg(MIN_MIGRAD);
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
        mpi_printf("-- pre-fit...\n");
    }
    else
    {
        mpi_printf("-- fitting and resampling %s...\n",param->q_name);
    }
    rs_data_fit(s_fit,limit,s_x,s_pt,d,NO_COR,use_x_var);
    if (IS_AN(param,AN_PHYPT))
    {
        fit_residual(res[s],d,0,fit);
    }
    if (IS_AN(param,AN_SCALE))
    {
        fit_residual(res[0],d,0,fit);
    }
    /** save chi^2/dof **/
    mpi_printf("chi^2/dof = %e\n",fit_data_get_chi2pdof(d));
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
        else if (strbufcmp(param->with_ext_a,"") != 0)
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
    if (param->nproc == 1)
    {
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
    }
    /** print results **/
    print_result(s_fit,param);
    if (IS_AN(param,AN_PHYPT))
    {
        mat_set(phy_pt,i_ud,0,SQ(param->M_ud));
        mat_set(phy_pt,i_s,0,SQ(param->M_s));
        mat_set(phy_pt,i_umd,0,param->M_umd_val);
        mat_set(phy_pt,i_alpha,0,param->alpha);
        mat_set(phy_pt,i_bind,0,0.0);
        mat_set(phy_pt,i_a,0,0.0);
        mat_set(phy_pt,i_Linv,0,0.0);
        mat_set(phy_pt,i_fvM,0,param->qed_fvol_mass);
        param->scale_model = 1;
        mat_set(phy_pt,i_umd,0,param->M_umd_val);
        buf = fit_data_model_xeval(d,s,phy_pt,rs_sample_pt_cent_val(s_fit));
        mat_set(rs_sample_pt_cent_val(param->s_ex),0,0,buf);
        for (i=0;i<param->nsample;i++)
        {
            if (strbufcmp(param->with_ext_M_umd,"") != 0)
            {
                mat_set(phy_pt,i_umd,0,                                    \
                        mat_get(rs_sample_pt_sample(param->s_M_umd,i),0,0));
            }
            buf = fit_data_model_xeval(d,s,phy_pt,                  \
                                       rs_sample_pt_sample(s_fit,i));
            mat_set(rs_sample_pt_sample(param->s_ex,i),0,0,buf);
        }
        param->scale_model = 0;
        rs_sample_varp(ex_err,param->s_ex);
        mat_eqsqrt(ex_err);
        mpi_printf("extrapolation :\n");
        mpi_printf("%10s = %f +/- %e ( %4.1f%% ) MeV^%d\n",param->q_name,\
                   mat_get(ex,0,0),mat_get(ex_err,0,0),                \
                   mat_get(ex_err,0,0)/fabs(mat_get(ex,0,0))*100.0,    \
                   param->q_dim);
        if ((!latan_isnan(param->q_target[0]))  \
            &&(!latan_isnan(param->q_target[1])))
        {
            mpi_printf("\n");
            mpi_printf("compatible with target within %4.2f sigmas\n\n",    \
                       fabs(mat_get(ex,0,0)-param->q_target[0])             \
                       /sqrt(SQ(mat_get(ex_err,0,0))+SQ(param->q_target[1])));
        }
    }
    /** display plots **/
    if (param->plot)
    {
        mpi_printf("\n-- generating plots...\n");
        if (IS_AN(param,AN_PHYPT))
        {
            plot_chi2_comp(d,param,s,"physical point fit");
            SCALE_DATA(AINV);
            fit_data_set_covar_from_sample(d,s_x,s_pt,NO_COR,use_x_var);
            plot_fit(s_fit,d,param,Q);
            SCALE_DATA(A);
        }
        if (IS_AN(param,AN_SCALE))
        {
            plot_chi2_comp(d,param,0,"scale setting");
            plot_fit(s_fit,d,param,SCALE);
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
                            ||((param->s_umd_deg != 0)&&IS_AN(param,AN_SCALE)))\
        &&param->have_umd;
        mpi_printf("-- fitting and resampling %s...\n",param->q_name);
        rs_data_fit(s_fit,limit,s_x,s_pt,d,X_COR|XDATA_COR|DATA_COR,use_x_var);
        if (IS_AN(param,AN_PHYPT))
        {
            fit_residual(res[s],d,0,fit);
        }
        if (IS_AN(param,AN_SCALE))
        {
            fit_residual(res[0],d,0,fit);
        }
        mpi_printf("chi^2/dof = %e\n",fit_data_get_chi2pdof(d));
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
            else if (strbufcmp(param->with_ext_a,"") != 0)
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
        if (param->nproc == 1)
        {
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
        }
        /** print results **/
        print_result(s_fit,param);
        if (IS_AN(param,AN_PHYPT))
        {
            param->scale_model = 1;
            mat_set(phy_pt,i_umd,0,param->M_umd_val);
            buf = fit_data_model_xeval(d,s,phy_pt,rs_sample_pt_cent_val(s_fit));
            mat_set(rs_sample_pt_cent_val(param->s_ex),0,0,buf);
            for (i=0;i<param->nsample;i++)
            {
                if (strbufcmp(param->with_ext_M_umd,"") != 0)
                {
                    mat_set(phy_pt,i_umd,0,                                    \
                            mat_get(rs_sample_pt_sample(param->s_M_umd,i),0,0));
                }
                buf = fit_data_model_xeval(d,s,phy_pt,                  \
                                           rs_sample_pt_sample(s_fit,i));
                mat_set(rs_sample_pt_sample(param->s_ex,i),0,0,buf);
            }
            param->scale_model = 0;
            rs_sample_varp(ex_err,param->s_ex);
            mat_eqsqrt(ex_err);
            mpi_printf("extrapolation :\n");
            mpi_printf("%10s = %f +/- %e ( %4.1f%% ) MeV^%d\n",param->q_name,\
                       mat_get(ex,0,0),mat_get(ex_err,0,0),                \
                       mat_get(ex_err,0,0)/fabs(mat_get(ex,0,0))*100.0,    \
                       param->q_dim);
            if ((!latan_isnan(param->q_target[0]))  \
                &&(!latan_isnan(param->q_target[1])))
            {
                mpi_printf("\n");
                mpi_printf("compatible with target within %4.2f sigmas\n\n",  \
                           fabs(mat_get(ex,0,0)-param->q_target[0])             \
                           /sqrt(SQ(mat_get(ex_err,0,0))+SQ(param->q_target[1])));
            }
        }
        /** display plots **/
        if (param->plot)
        {
            mpi_printf("\n-- generating plots...\n");
            if (IS_AN(param,AN_PHYPT))
            {
                SCALE_DATA(AINV);
                fit_data_set_covar_from_sample(d,s_x,s_pt,NO_COR,use_x_var);
                plot_fit(s_fit,d,param,Q);
                SCALE_DATA(A);
                fit_data_set_covar_from_sample(d,s_x,s_pt,NO_COR,use_x_var);
            }
            if (IS_AN(param,AN_SCALE))
            {
                plot_fit(s_fit,d,param,SCALE);
            }
        }
    }
    
    /*              result output               */
    /********************************************/
    strbuf res_path;
    char mode;
    
    if (IS_AN(param,AN_PHYPT))
    {
        if (param->save_result)
        {
            sprintf(res_path,"%s.boot%c%s",param->result_file,LATAN_PATH_SEP,\
                    param->result_file);
            rs_sample_save(res_path,'w',param->s_ex);
        }
    }
    if (IS_AN(param,AN_SCALE))
    {
        mpi_printf("scales :\n\n");
        for (i=0;i<param->nbeta;i++)
        {
            if (param->save_result)
            {
                sprintf(res_path,"%s.boot%c%s_%s",param->result_file,   \
                        LATAN_PATH_SEP,param->result_file,param->beta[i]);
                mode = (i == 0) ? 'w' : 'a';
                rs_sample_save_subsamp(res_path,mode,s_fit,i,0,i,0);
            }
            rs_sample_varp(fit_var,s_fit);
            mpi_printf("beta = %s\n",param->beta[i]);
            mpi_printf("a    = %f +/- %e fm\n",mat_get(fit,i,0)/NU_FM,\
                       sqrt(mat_get(fit_var,i,0))/NU_FM);
            rs_sample_eqinvp(s_fit);
            rs_sample_varp(fit_var,s_fit);
            mpi_printf("a^-1 = %f +/- %e MeV\n",mat_get(fit,i,0),\
                       sqrt(mat_get(fit_var,i,0)));
            mpi_printf("\n");
            rs_sample_eqinvp(s_fit);
        }
    }
    if (param->save_all_param)
    {
        for (i=0;i<npar;i++)
        {
            mode = ((i == 0)&&!param->save_result) ? 'w' : 'a';
            sprintf(res_path,"%s.boot%c%s_p_%d",param->result_file,\
                    LATAN_PATH_SEP,param->result_file,(int)i);
            rs_sample_save_subsamp(res_path,mode,s_fit,i,0,i,0);
        }
    }
    else
    {
        for (i=0;i<param->nsave_param;i++)
        {
            mode = ((i == 0)&&!param->save_result) ? 'w' : 'a';
            sprintf(res_path,"%s.boot%c%s_p_%d",param->result_file,\
                    LATAN_PATH_SEP,param->result_file,             \
                    (int)param->save_param[i]);
            rs_sample_save_subsamp(res_path,mode,s_fit,param->save_param[i],0,\
                                   param->save_param[i],0);
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
    mat_ar_destroy(res,2);
    mat_destroy(fit_var);
    mat_destroy(phy_pt);
    mat_destroy(limit);
}
                     
/* main program */
int main(int argc, char *argv[])
{
    /* initialization */
    io_init();

#ifdef HAVE_MPI
    int status;
    
    status = MPI_Init(&argc,&argv);
    if (status != MPI_SUCCESS)
    {
        MPI_Abort(MPI_COMM_WORLD,status);
    }
#endif
    
    /* argument parsing */
    fit_param *param;
    int proc,nproc;
    bool active;
    strbuf prefix;
    
    if (argc < 2)
    {
        fprintf(stderr,"usage: %s <par_file_1> [<par_file_2> ...]\n",argv[0]);
        return EXIT_FAILURE;
    }
#ifdef HAVE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD,&proc);
    MPI_Comm_size(MPI_COMM_WORLD,&nproc);
    if (nproc < argc - 1)
    {
        fprintf(stderr,"error: %d processes for %d parameter files\n",nproc,\
                argc-1);
        return EXIT_FAILURE;
    }
    if (nproc > 1)
    {
        if (proc+1 < argc)
        {
            active = true;
            mpi_printf("*** parameter file %s\n",argv[proc+1]);
            param       = fit_param_parse(argv[proc+1]);
            param->plot = 0;
            sprintf(prefix,"[%d] ",proc);
            latan_set_use_car_ret(false);
            latan_set_msg_prefix(prefix);
            strbufcpy(param->save_plot,"");
        }
        else
        {
            active = false;
            mpi_printf("*** inactive process\n");
            param = NULL;
        }
    }
    else
    {
        active = true;
        param = fit_param_parse(argv[1]);
    }
#else
    proc   = 0;
    nproc  = 1;
    param  = fit_param_parse(argv[1]);
    active = true;
#endif
    
    /* analysis */
    param->nproc = nproc;
    if (active)
    {
        analysis(param);
    }
    
    /* finalization */
    io_finish();
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    fit_param_destroy(param);
    
    return EXIT_SUCCESS;
}
