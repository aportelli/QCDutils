#define _POSIX_C_SOURCE 199506L /* strtok_r is used here */

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

#define ATOI(str) ((int)strtol(str,(char **)NULL,10))

void data_load(rs_sample *s_x[N_EX_VAR], rs_sample *s_q[2], fit_param *param)
{
    strbuf M_str,Msq_str,ext,sf_name,ud_lbl,s_lbl;
    size_t nsample,nens;
    size_t ens_ind,bind,d;
    int i;
    rs_sample *s_tmp;
    ens *ens_pt;
    mat **x_err,*q_err[2];
    FILE *table_dump,*outf[2];

    nsample  = param->nsample;
    nens     = param->nens;
    
    s_tmp    = rs_sample_create(1,nsample);
    x_err    = mat_ar_create(N_EX_VAR,nens,1);
    q_err[0] = mat_create(nens,1);
    q_err[1] = mat_create(nens,1);
    
    if (IS_ANALYZE(param,"phypt"))
    {
        sprintf(ud_lbl,"M_%s^2",param->ud_name);
        sprintf(s_lbl,"M_%s^2",param->s_name);
    }
    else if (IS_ANALYZE(param,"scaleset")||IS_ANALYZE(param,"comb_phypt_scale"))
    {
        sprintf(ud_lbl,"(a*M_%s)^2",param->ud_name);
        sprintf(s_lbl,"(a*M_%s)^2",param->s_name);
    }
    strbufcpy(M_str,"M");
    strbufcpy(Msq_str,"Msq");
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
    for(ens_ind=0;ens_ind<nens;ens_ind++) 
    {
        ens_pt = param->point + ens_ind;
        bind   = ind_beta(ens_pt->beta,param);
        for (d=0;d<param->ndataset;d++)
        {
            if (strcmp(param->dataset[d],ens_pt->dataset) == 0)
            {
                /* scale setting quantity */
                sprintf(sf_name,"%s/M_%s_%s.boot%s",ens_pt->dir,\
                        param->scale_part,ens_pt->dataset,ext);
                rs_sample_load_subsamp(s_tmp,sf_name,"",0,0);
                rs_sample_set_subsamp(s_q[0],s_tmp,ens_ind,ens_ind);
                /* main quantity (same as the previous one in 'scaleset') */
                sprintf(sf_name,"%s/%s_%s.boot%s",ens_pt->dir,param->q_name,\
                        param->dataset[d],ext);
                rs_sample_load_subsamp(s_tmp,sf_name,"",0,0);
                rs_sample_set_subsamp(s_q[1],s_tmp,ens_ind,ens_ind);
                /* m_ud fixing quantity */
                sprintf(sf_name,"%s/%s_%s_%s.boot%s",ens_pt->dir,Msq_str,\
                        param->ud_name,param->dataset[d],ext);
                rs_sample_load_subsamp(s_tmp,sf_name,"",0,0);
                rs_sample_set_subsamp(s_x[i_ud],s_tmp,ens_ind,ens_ind);
                /* m_s fixing quantity */
                sprintf(sf_name,"%s/%s_%s_%s.boot%s",ens_pt->dir,Msq_str,\
                        param->s_name,param->dataset[d],ext);
                rs_sample_load_subsamp(s_tmp,sf_name,"",0,0);
                rs_sample_set_subsamp(s_x[i_s],s_tmp,ens_ind,ens_ind);
                /* beta index */
                rs_sample_cst(s_tmp,bind);
                rs_sample_set_subsamp(s_x[i_bind],s_tmp,ens_ind,ens_ind);
                /* m_u - m_d fixing quantity if possible */
                sprintf(sf_name,"%s/%s_%s.boot%s",ens_pt->dir,\
                        param->umd_name,param->dataset[d],ext);
                if (access(sf_name,R_OK) == 0)
                {
                    rs_sample_load_subsamp(s_tmp,sf_name,"",0,0);
                    rs_sample_set_subsamp(s_x[i_umd],s_tmp,ens_ind,ens_ind);
                    param->have_umd = 1;
                }
                /* spatial extent */
                rs_sample_cst(s_tmp,1.0/((double)ens_pt->L));
                rs_sample_set_subsamp(s_x[i_Linv],s_tmp,ens_ind,ens_ind);
                /* lattice spacing */
                if (IS_ANALYZE(param,"phypt"))
                {
                    if (param->with_ext_a)
                    {
                        sprintf(sf_name,"./a_%s_%s_%s.boot%s",ens_pt->beta,\
                                param->scale_part,param->dataset_cat,ext);
                        rs_sample_load_subsamp(s_tmp,sf_name,"",0,0);
                        rs_sample_set_subsamp(s_x[i_a],s_tmp,ens_ind,ens_ind);
                        rs_sample_set_subsamp(param->a,s_tmp,bind,bind);
                    }
                    else
                    {
                        sprintf(sf_name,"%s/%s_%s_%s.boot%s",ens_pt->dir,M_str,\
                                param->scale_part,param->dataset[d],ext);
                        rs_sample_load_subsamp(s_tmp,sf_name,"",0,0);
                        rs_sample_eqmuls(s_tmp,1.0/param->M_scale);
                        rs_sample_set_subsamp(s_x[i_a],s_tmp,ens_ind,ens_ind);
                    }
                }
                else if (IS_ANALYZE(param,"scaleset")          \
                         ||IS_ANALYZE(param,"comb_phypt_scale"))
                {
                    rs_sample_cst(s_tmp,1.0);
                    rs_sample_set_subsamp(s_x[i_a],s_tmp,ens_ind,ens_ind);
                }
            }
        }
    }
    
    /* computing errors */
    rs_sample_varp(x_err[i_ud],s_x[i_ud]);
    mat_eqsqrt(x_err[i_ud]);
    rs_sample_varp(x_err[i_s],s_x[i_s]);
    mat_eqsqrt(x_err[i_s]);
    rs_sample_varp(x_err[i_a],s_x[i_a]);
    mat_eqsqrt(x_err[i_a]);
    rs_sample_varp(x_err[i_umd],s_x[i_umd]);
    mat_eqsqrt(x_err[i_umd]);
    rs_sample_varp(x_err[i_Linv],s_x[i_Linv]);
    mat_eqsqrt(x_err[i_Linv]);
    rs_sample_varp(q_err[0],s_q[0]);
    mat_eqsqrt(q_err[0]);
    rs_sample_varp(q_err[1],s_q[1]);
    mat_eqsqrt(q_err[1]);
    rs_sample_varp(param->a_err,param->a);
    mat_eqsqrt(param->a_err);
    
    /* display */
#define PRINT_DLABEL(name) fprintf(outf[i],"%-12s  ",name)
#define PRINT_DLABEL_WERR(name) fprintf(outf[i],"%-12s %-12s  ",name,"error")
#define PRINT_D(value) fprintf(outf[i],"% .5e  ",value)
#define PRINT_D_WERR(value,err) fprintf(outf[i],"% .5e % .5e  ",value,err)
#define PRINT_X(ind,dim)\
{\
    double fac;\
    fac = pow(1.0/mat_get(rs_sample_pt_cent_val(s_x[i_a]),ens_ind,0),dim);\
    PRINT_D(mat_get(rs_sample_pt_cent_val(s_x[ind]),ens_ind,0)*fac);\
}
#define PRINT_X_WERR(ind,dim)\
{\
    double fac;\
    fac = pow(1.0/mat_get(rs_sample_pt_cent_val(s_x[i_a]),ens_ind,0),dim);\
    PRINT_D_WERR(mat_get(rs_sample_pt_cent_val(s_x[ind]),ens_ind,0)*fac,\
                 mat_get(x_err[ind],ens_ind,0)*fac);\
}
#define PRINT_CV_WERR(s,s_err,dim)\
{\
    double fac;\
    fac = pow(1.0/mat_get(rs_sample_pt_cent_val(s_x[i_a]),ens_ind,0),dim);\
    PRINT_D_WERR(mat_get(rs_sample_pt_cent_val(s),ens_ind,0)*fac,\
                 mat_get(s_err,ens_ind,0)*fac);\
}
    
    
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
            PRINT_DLABEL_WERR(param->umd_name);
        }
        if (param->with_qed_fvol)
        {
            PRINT_DLABEL_WERR("1/L");
        }
        PRINT_DLABEL_WERR(param->q_name);
        fprintf(outf[i],"\n");
        for(ens_ind=0;ens_ind<nens;ens_ind++)
        {
            fprintf(outf[i],"%30s  ",param->point[ens_ind].dir);
            PRINT_X_WERR(i_ud,2);
            PRINT_X_WERR(i_s,2);
            if (IS_ANALYZE(param,"phypt"))
            {
                PRINT_X_WERR(i_a,0);
            }
            if (param->with_umd)
            {
                PRINT_X_WERR(i_umd,2);
            }
            if (param->with_qed_fvol)
            {
                PRINT_X_WERR(i_Linv,1);
            }
            PRINT_CV_WERR(s_q[1],q_err[1],param->q_dim);
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
    mat_destroy(q_err[0]);
    mat_destroy(q_err[1]);
}
