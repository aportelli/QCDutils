#include "data_loader.h"
#include "parameters.h"
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <sys/types.h>
#include <latan/latan_io.h>
#include <latan/latan_nunits.h>

/* beta -> a function */
double get_cst_ainv(double beta)
{
    double ainv;

    if (beta == 3.31)
    {
        ainv = 1697.0;
    }
    else if (beta == 3.5)
    {
        ainv = 2131.0;
    }
    else if (beta == 3.61)
    {
        ainv = 2561.0;
    }
    else
    {
        fprintf(stderr,"error: beta %f unknown\n",beta);
        exit(EXIT_FAILURE);
    }

    return ainv;
}

void data_load(rs_sample *s_x[N_EX_VAR], rs_sample *s_q, strbuf beta,\
               const ex_param *param)
{
    strbuf *field,*ens,*ens_ar,M_str,ext,sf_name,pi_lbl,K_lbl;
    int lc,nf,i;
    size_t nsample,nens;
    size_t ens_ind;
    int T_ens,L_ens;
    rs_sample *s_tmp;
    mat **x_err,*q_err;
    FILE *table_dump,*outf[2];

    nsample = rs_sample_get_nsample(s_q);
    nens    = rs_sample_get_nrow(s_q);
    ens_ind = 0;
    field   = NULL;
    
    s_tmp   = rs_sample_create(1,nsample);
    x_err   = mat_ar_create(N_EX_VAR,nens,1);
    q_err   = mat_create(nens,1);
    ens_ar  = (strbuf *)malloc(nens*sizeof(strbuf));
    
    switch (param->ex_dim)
    {
        case 0:
            sprintf(pi_lbl,"(a*M_%s)^2",param->pi_name);
            sprintf(K_lbl,"(a*M_%s)^2",param->K_name);
            strbufcpy(M_str,"Msq");
            break;
        case 1:
            sprintf(pi_lbl,"M_%s",param->pi_name);
            sprintf(K_lbl,"M_%s",param->K_name);
            strbufcpy(M_str,"M");
            break;
        case 2:
            sprintf(pi_lbl,"M_%s^2",param->pi_name);
            sprintf(K_lbl,"M_%s^2",param->K_name);
            strbufcpy(M_str,"Msq");
            break;
        default:
            fprintf(stderr,"error: pi/K mass dimension must be 0,1 or 2\n");
            abort();
            break;
    }
    switch (io_get_fmt())
    {
        case IO_XML:
            strbufcpy(ext,".xml");
            break;
        case IO_ASCII:
            strbufcpy(ext,"");
            break;
        default:
            fprintf(stderr,"error: I/O format unknown\n");
            abort();
            break;
    }
    
    /* data loading */
    BEGIN_FOR_LINE_TOK(field,param->manifest,"_",nf,lc)
    if (field[0][0] != '#')
    {
        ens = ens_ar + ens_ind;
        sprintf(*ens,"%s_%s_%s_%s_%s",field[0],field[1],field[2],field[3],\
                field[4]);
        T_ens = atoi(field[0]);
        L_ens = atoi(field[1]);
        strbufcpy(beta,field[2]);
        sprintf(sf_name,"%s/%s_%s%s.boot%s",*ens,M_str,param->pi_name,\
                param->suffix,ext);
        rs_sample_load_subsamp(s_tmp,sf_name,"",0,0);
        rs_sample_set_subsamp(s_x[i_pi],s_tmp,ens_ind,ens_ind);
        if (param->ex_dim != 1)
        {
            rs_sample_eqsqrt(s_tmp);
        }
        rs_sample_eqmuls(s_tmp,-(double)(L_ens));
        rs_sample_eqexpp(s_tmp);
        rs_sample_set_subsamp(s_x[i_MpiL],s_tmp,ens_ind,ens_ind);
        sprintf(sf_name,"%s/%s_%s%s.boot%s",*ens,M_str,param->K_name,\
                param->suffix,ext);
        rs_sample_load_subsamp(s_tmp,sf_name,"",0,0);
        rs_sample_set_subsamp(s_x[i_K],s_tmp,ens_ind,ens_ind);
        rs_sample_cst(s_tmp,ind_beta(beta,param));
        rs_sample_set_subsamp(s_x[i_bind],s_tmp,ens_ind,ens_ind);
        if (param->ex_dim != 0)
        {
            sprintf(sf_name,"./scale_%s_%s%s.boot%s",beta,param->scale_part,\
                    param->suffix,ext);
            if (access(sf_name,R_OK) == 0)
            {
                rs_sample_load_subsamp(s_tmp,sf_name,"",0,0);
            }
            else
            {
                rs_sample_cst(s_tmp,get_cst_ainv(atof(beta)));
            }
            rs_sample_set_subsamp(s_x[i_ainv],s_tmp,ens_ind,ens_ind);
        }
        else
        {
            rs_sample_cst(s_x[i_ainv],1.0);
        }
        if (param->with_umd)
        {
            sprintf(sf_name,"%s/dm_ud%s.boot%s",*ens,param->suffix,ext);
            rs_sample_load_subsamp(s_tmp,sf_name,"",0,0);
            rs_sample_set_subsamp(s_x[i_umd],s_tmp,ens_ind,ens_ind);
        }
        if (param->with_fvol)
        {
            rs_sample_cst(s_tmp,1.0/((double)L_ens));
            rs_sample_set_subsamp(s_x[i_Linv],s_tmp,ens_ind,ens_ind);
        }
        sprintf(sf_name,"%s/%s%s.boot%s",*ens,param->q_name,param->suffix,ext);
        rs_sample_load_subsamp(s_tmp,sf_name,"",0,0);
        rs_sample_set_subsamp(s_q,s_tmp,ens_ind,ens_ind);
        ens_ind++;
    }
    END_FOR_LINE_TOK(field);

    /* scaling */
    for (i=1;i<=param->ex_dim;i++)
    {
        rs_sample_eqmulp(s_x[i_pi],s_x[i_ainv]);
        rs_sample_eqmulp(s_x[i_K],s_x[i_ainv]);
        /*rs_sample_eqmulp(s_x[i_umd],s_x[i_ainv]);*/
    }
    if (param->ex_dim > 0)
    {
        rs_sample_eqmulp(s_x[i_Linv],s_x[i_ainv]);
    }
    for (i=1;i<=param->q_dim;i++)
    {
        rs_sample_eqmulp(s_q,s_x[i_ainv]);
    }
    if (param->ex_dim != 0)
    {
        rs_sample_eqinvp(s_x[i_ainv]);
    }
    
    /* computing errors */
    rs_sample_varp(x_err[i_pi],s_x[i_pi]);
    mat_eqsqrt(x_err[i_pi]);
    rs_sample_varp(x_err[i_K],s_x[i_K]);
    mat_eqsqrt(x_err[i_K]);
    rs_sample_varp(x_err[i_ainv],s_x[i_ainv]);
    mat_eqsqrt(x_err[i_ainv]);
    rs_sample_varp(x_err[i_umd],s_x[i_umd]);
    mat_eqsqrt(x_err[i_umd]);
    rs_sample_varp(x_err[i_Linv],s_x[i_Linv]);
    mat_eqsqrt(x_err[i_Linv]);
    rs_sample_varp(x_err[i_MpiL],s_x[i_MpiL]);
    mat_eqsqrt(x_err[i_MpiL]);
    rs_sample_varp(q_err,s_q);
    mat_eqsqrt(q_err);
    
    /* display */
#define PRINT_DLABEL(name) fprintf(outf[i],"%-12s  ",name)
#define PRINT_DLABEL_WERR(name) fprintf(outf[i],"%-12s %-12s  ",name,"error")
#define PRINT_D(value) fprintf(outf[i],"% .5e  ",value)
#define PRINT_D_WERR(value,err) fprintf(outf[i],"% .5e % .5e  ",value,err)
#define PRINT_X(ind) PRINT_D(mat_get(rs_sample_pt_cent_val(s_x[ind]),ens_ind,0))
#define PRINT_X_WERR(ind) \
PRINT_D_WERR(mat_get(rs_sample_pt_cent_val(s_x[ind]),ens_ind,0),\
             mat_get(x_err[ind],ens_ind,0))
#define PRINT_Q_WERR \
PRINT_D_WERR(mat_get(rs_sample_pt_cent_val(s_q),ens_ind,0),\
             mat_get(q_err,ens_ind,0))
    
    
    table_dump = fopen("phyfit_table.dat","w");
    outf[0]    = stdout;
    outf[1]    = table_dump;
    printf("\n");
    for (i=0;i<2;i++)
    {
        fprintf(outf[i],"#%29s  ","ensemble");
        PRINT_DLABEL_WERR(pi_lbl);
        PRINT_DLABEL_WERR(K_lbl);
        if (param->ex_dim > 0)
        {
            PRINT_DLABEL_WERR("a");
        }
        if (param->with_umd)
        {
            PRINT_DLABEL_WERR("a*(m_u-m_d)");
        }
        if (param->with_fvol)
        {
            PRINT_DLABEL_WERR("1/L");
            PRINT_DLABEL_WERR("M_pi*L");
        }
        PRINT_DLABEL_WERR(param->q_name);
        fprintf(outf[i],"\n");
        for(ens_ind=0;ens_ind<nens;ens_ind++)
        {
            fprintf(outf[i],"%30s  ",ens_ar[ens_ind]);
            PRINT_X_WERR(i_pi);
            PRINT_X_WERR(i_K);
            if (param->ex_dim > 0)
            {
                PRINT_X_WERR(i_ainv);
            }
            if (param->with_umd)
            {
                PRINT_X_WERR(i_umd);
            }
            if (param->with_fvol)
            {
                PRINT_X_WERR(i_Linv);
                PRINT_X_WERR(i_MpiL);
            }
            PRINT_Q_WERR;
            fprintf(outf[i],"\n");
        }
        fprintf(outf[i],"\n");
    }
    fclose(table_dump);

#undef PRINT_DLABEL
#undef PRINT_DLABEL_WERR
#undef PRINT_D
#undef PRINT_D_WERR
#undef PRINT_X
#undef PRINT_X_WERR
#undef PRINT_Q_WERR
    
    rs_sample_destroy(s_tmp);
    mat_ar_destroy(x_err,N_EX_VAR);
    mat_destroy(q_err);
    free(ens_ar);
}
