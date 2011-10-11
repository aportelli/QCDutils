#include "data_loader.h"
#include "parameters.h"
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <sys/types.h>
#include <latan/latan_io.h>
#include <latan/latan_math.h>
#include <latan/latan_nunits.h>

void data_load(rs_sample *s_x[N_EX_VAR], rs_sample *s_q, strbuf beta,\
               const ex_param *param)
{
    strbuf *field,*ens,*ens_ar,M_str,ext,sf_name,ud_lbl,s_lbl;
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
    
    if (IS_ANALYZE(param,"phypt"))
    {
        sprintf(ud_lbl,"M_%s^2",param->ud_name);
        sprintf(s_lbl,"M_%s^2",param->s_name);
    }
    else if (IS_ANALYZE(param,"scaleset"))
    {
        sprintf(ud_lbl,"(a*M_%s)^2",param->ud_name);
        sprintf(s_lbl,"(a*M_%s)^2",param->s_name);
    }
    strbufcpy(M_str,"Msq");
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
        sprintf(sf_name,"%s/%s_%s%s.boot%s",*ens,M_str,param->ud_name,\
                param->suffix,ext);
        rs_sample_load_subsamp(s_tmp,sf_name,"",0,0);
        rs_sample_set_subsamp(s_x[i_ud],s_tmp,ens_ind,ens_ind);
        sprintf(sf_name,"%s/%s_%s%s.boot%s",*ens,M_str,param->s_name,\
                param->suffix,ext);
        rs_sample_load_subsamp(s_tmp,sf_name,"",0,0);
        rs_sample_set_subsamp(s_x[i_s],s_tmp,ens_ind,ens_ind);
        rs_sample_cst(s_tmp,ind_beta(beta,param));
        rs_sample_set_subsamp(s_x[i_bind],s_tmp,ens_ind,ens_ind);
        if (param->with_umd)
        {
            sprintf(sf_name,"%s/dm_ud%s.boot%s",*ens,param->suffix,ext);
            rs_sample_load_subsamp(s_tmp,sf_name,"",0,0);
            rs_sample_set_subsamp(s_x[i_umd],s_tmp,ens_ind,ens_ind);
        }
        if (param->with_qed_fvol)
        {
            rs_sample_cst(s_tmp,1.0/((double)L_ens));
            rs_sample_set_subsamp(s_x[i_Linv],s_tmp,ens_ind,ens_ind);
        }
        if (IS_ANALYZE(param,"phypt"))
        {
            if (param->with_ext_a)
            {
                sprintf(sf_name,"./scale_%s_%s%s.boot%s",beta,\
                        param->scale_part,param->suffix,ext);
                rs_sample_load_subsamp(s_tmp,sf_name,"",0,0);
                rs_sample_set_subsamp(s_x[i_ainv],s_tmp,ens_ind,ens_ind);
            }
            else
            {
                sprintf(sf_name,"%s/%s_%s%s.boot%s",*ens,M_str,\
                        param->scale_part,param->suffix,ext);
                rs_sample_load_subsamp(s_tmp,sf_name,"",0,0);
                rs_sample_set_subsamp(s_x[i_ainv],s_tmp,ens_ind,ens_ind);
            }
        }
        sprintf(sf_name,"%s/%s%s.boot%s",*ens,param->q_name,param->suffix,ext);
        rs_sample_load_subsamp(s_tmp,sf_name,"",0,0);
        rs_sample_set_subsamp(s_q,s_tmp,ens_ind,ens_ind);
        ens_ind++;
    }
    END_FOR_LINE_TOK(field);

    /* scaling */
    if (IS_ANALYZE(param,"phypt"))
    {
        if (!param->with_ext_a)
        {
            rs_sample_eqmuls(s_x[i_ainv],1.0/SQ(param->M_scale));
            rs_sample_eqsqrt(s_x[i_ainv]);
            rs_sample_eqinvp(s_x[i_ainv]);
        }
        rs_sample_eqmulp(s_x[i_ud],s_x[i_ainv]);
        rs_sample_eqmulp(s_x[i_ud],s_x[i_ainv]);
        rs_sample_eqmulp(s_x[i_s],s_x[i_ainv]);
        rs_sample_eqmulp(s_x[i_s],s_x[i_ainv]);
        rs_sample_eqmulp(s_x[i_umd],s_x[i_ainv]);
        rs_sample_eqmulp(s_x[i_Linv],s_x[i_ainv]);
        for (i=1;i<=param->q_dim;i++)
        {
            rs_sample_eqmulp(s_q,s_x[i_ainv]);
        }
        rs_sample_eqinvp(s_x[i_ainv]);
    }
    
    /* computing errors */
    rs_sample_varp(x_err[i_ud],s_x[i_ud]);
    mat_eqsqrt(x_err[i_ud]);
    rs_sample_varp(x_err[i_s],s_x[i_s]);
    mat_eqsqrt(x_err[i_s]);
    rs_sample_varp(x_err[i_ainv],s_x[i_ainv]);
    mat_eqsqrt(x_err[i_ainv]);
    rs_sample_varp(x_err[i_umd],s_x[i_umd]);
    mat_eqsqrt(x_err[i_umd]);
    rs_sample_varp(x_err[i_Linv],s_x[i_Linv]);
    mat_eqsqrt(x_err[i_Linv]);
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
#define PRINT_CV_WERR(s,s_err) \
PRINT_D_WERR(mat_get(rs_sample_pt_cent_val(s),ens_ind,0),\
             mat_get(s_err,ens_ind,0))
    
    
    table_dump = fopen("phyfit_table.dat","w");
    outf[0]    = stdout;
    outf[1]    = table_dump;
    printf("\n");
    for (i=0;i<2;i++)
    {
        fprintf(outf[i],"#%29s  ","ensemble");
        PRINT_DLABEL_WERR(ud_lbl);
        PRINT_DLABEL_WERR(s_lbl);
        if (IS_ANALYZE(param,"phypt"))
        {
            PRINT_DLABEL_WERR("a");
        }
        if (param->with_umd)
        {
            PRINT_DLABEL_WERR("m_u-m_d");
        }
        if (param->with_qed_fvol)
        {
            PRINT_DLABEL_WERR("1/L");
        }
        PRINT_DLABEL_WERR(param->q_name);
        fprintf(outf[i],"\n");
        for(ens_ind=0;ens_ind<nens;ens_ind++)
        {
            fprintf(outf[i],"%30s  ",ens_ar[ens_ind]);
            PRINT_X_WERR(i_ud);
            PRINT_X_WERR(i_s);
            if (IS_ANALYZE(param,"phypt"))
            {
                PRINT_X_WERR(i_ainv);
            }
            if (param->with_umd)
            {
                PRINT_X_WERR(i_umd);
            }
            if (param->with_qed_fvol)
            {
                PRINT_X_WERR(i_Linv);
            }
            PRINT_CV_WERR(s_q,q_err);
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
