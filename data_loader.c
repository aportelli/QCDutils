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
    strbuf M_str,Msq_str,ext,sf_name;
    size_t nsample,nens;
    size_t ens_ind,bind,d;
    rs_sample *s_tmp;
    ens *ens_pt;

    nsample  = param->nsample;
    nens     = param->nens;
    
    s_tmp    = rs_sample_create(1,nsample);
    
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
                rs_sample_varp(param->a_err,param->a);
                mat_eqsqrt(param->a_err);
            }
        }
    }
}
