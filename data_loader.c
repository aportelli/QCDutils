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

#define ATOF(str) (strtod(str,(char **)NULL))
#define ATOI(str) ((int)strtol(str,(char **)NULL,10))

void data_load(rs_sample *s_x[N_EX_VAR], rs_sample *s_q[2], fit_param *param)
{
    strbuf M_str,Msq_str,ext,sf_name;
    size_t nsample,nens;
    size_t ens_ind,bind,vind,d,i,j;
    rs_sample *s_tmp;
    FILE *ef;
    double e;
    ens *ens_pt;

    nsample  = param->nsample;
    nens     = param->nens;
    
    s_tmp    = rs_sample_create(1,1,nsample);
    
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
    for (i=0;i<param->nbeta;i++)
    {
        rs_sample_cst(param->s_vol_av[i],0.0);
    }
    for(ens_ind=0;ens_ind<nens;ens_ind++)
    {
        ens_pt = param->point + ens_ind;
        bind   = (size_t)ind_beta(ens_pt->beta,param);
        vind   = (size_t)ind_volume((unsigned int)ens_pt->L,(int)bind,param);
        for (d=0;d<param->ndataset;d++)
        {
            if (strbufcmp(param->dataset[d],ens_pt->dataset) == 0)
            {
                /* scale setting quantity */
                if (IS_AN(param,AN_SCALE))
                {
                    sprintf(sf_name,"%s/%s_%s_%s.boot%s",ens_pt->dir,M_str,\
                            param->scale_part,ens_pt->dataset,ext);
                    rs_sample_load_subsamp(s_tmp,sf_name,0,0,0,0);
                    rs_sample_set_subsamp(s_q[0],s_tmp,ens_ind,0,ens_ind,0);
                }
                /* main quantity */
                if (IS_AN(param,AN_PHYPT))
                {
                    sprintf(sf_name,"%s/%s_%s.boot%s",ens_pt->dir,\
                            param->q_name,param->dataset[d],ext);
                    rs_sample_load_subsamp(s_tmp,sf_name,0,0,0,0);
                    rs_sample_set_subsamp(s_q[1],s_tmp,ens_ind,0,ens_ind,0);
                    rs_sample_subsamp_eqadd(param->s_vol_av[bind],s_tmp,\
                                            vind,0,vind,0);
                }
                /* m_ud fixing quantity */
                sprintf(sf_name,"%s/%s_%s_%s.boot%s",ens_pt->dir,Msq_str,\
                        param->ud_name,param->dataset[d],ext);
                rs_sample_load_subsamp(s_tmp,sf_name,0,0,0,0);
                rs_sample_set_subsamp(s_x[i_ud],s_tmp,ens_ind,0,ens_ind,0);
                /* m_s fixing quantity */
                sprintf(sf_name,"%s/%s_%s_%s.boot%s",ens_pt->dir,Msq_str,\
                        param->s_name,param->dataset[d],ext);
                rs_sample_load_subsamp(s_tmp,sf_name,0,0,0,0);
                rs_sample_set_subsamp(s_x[i_s],s_tmp,ens_ind,0,ens_ind,0);
                /* beta index */
                rs_sample_cst(s_tmp,bind);
                rs_sample_set_subsamp(s_x[i_bind],s_tmp,ens_ind,0,ens_ind,0);
                /* volume index */
                rs_sample_cst(s_tmp,vind);
                rs_sample_set_subsamp(s_x[i_vind],s_tmp,ens_ind,0,ens_ind,0);
                /* m_u - m_d fixing quantity if possible */
                sprintf(sf_name,"%s/%s_%s.boot%s",ens_pt->dir,\
                        param->umd_name,param->dataset[d],ext);
                if ((param->umd_deg > 0)||(param->with_udumd)        \
                    ||(param->with_sumd)||(access(sf_name,R_OK) == 0)\
                    ||(param->have_umd))
                {
                    rs_sample_load_subsamp(s_tmp,sf_name,0,0,0,0);
                    rs_sample_set_subsamp(s_x[i_umd],s_tmp,ens_ind,0,ens_ind,0);
                    param->have_umd = 1;
                    if (latan_isnan(param->M_umd_val))
                    {
                        fprintf(stderr,"error: %s data found but no physical value specified (parameter M_umd)\n",param->umd_name);
                        exit(EXIT_FAILURE);
                    }
                }
                /* electric charge if possible */
                sprintf(sf_name,"%s/emcoupl_%s",ens_pt->dir,param->dataset[d]);
                if ((param->alpha_deg > 0)||(param->with_udalpha)      \
                    ||(param->with_salpha)||(access(sf_name,R_OK) == 0)\
                    ||(param->have_alpha))
                {
                    ef = fopen(sf_name,"r");
                    if (!ef)
                    {
                        fprintf(stderr,"error: impossible to load charge file %s\n",
                                sf_name);
                        exit(EXIT_FAILURE);
                    }
                    fscanf(ef,"%lf",&e);
                    e = SQ(e)/(4.0*C_PI);
                    fclose(ef);
                    param->have_alpha = 1;
                    if (latan_isnan(param->alpha))
                    {
                        param->alpha = NU_ALPHA_EM;
                    }
                }
                else
                {
                    e = 0.0;
                }
                rs_sample_cst(s_tmp,e);
                rs_sample_set_subsamp(s_x[i_alpha],s_tmp,ens_ind,0,ens_ind,0);
                /* spatial extent */
                rs_sample_cst(s_tmp,1.0/((double)ens_pt->L));
                rs_sample_set_subsamp(s_x[i_Linv],s_tmp,ens_ind,0,ens_ind,0);
                rs_sample_set_subsamp(param->s_vol_Linv[bind],s_tmp,\
                                      vind,0,vind,0);
                rs_sample_cst(s_tmp,((double)ens_pt->T)/((double)ens_pt->L));
                rs_sample_set_subsamp(s_x[i_ToL],s_tmp,ens_ind,0,ens_ind,0);
                /* lattice spacing */
                if (IS_AN(param,AN_PHYPT)&&!IS_AN(param,AN_SCALE))
                {
                    if (strbufcmp(param->with_ext_a,"") != 0)
                    {
                        sprintf(sf_name,"./%s.boot%s%c%s_%s",param->with_ext_a,\
                                ext,LATAN_PATH_SEP,param->with_ext_a,          \
                                ens_pt->beta);
                        rs_sample_load_subsamp(s_tmp,sf_name,0,0,0,0);
                        rs_sample_set_subsamp(s_x[i_a],s_tmp,ens_ind,0,\
                                              ens_ind,0);
                        rs_sample_set_subsamp(param->s_a,s_tmp,bind,0,bind,0);
                    }
                    else
                    {
                        sprintf(sf_name,"%s/%s_%s_%s.boot%s",ens_pt->dir,M_str,\
                                param->scale_part,param->dataset[d],ext);
                        rs_sample_load_subsamp(s_tmp,sf_name,0,0,0,0);
                        rs_sample_eqmuls(s_tmp,1.0/param->M_scale);
                        rs_sample_set_subsamp(s_x[i_a],s_tmp,ens_ind,0,\
                                              ens_ind,0);
                    }
                }
                else if (IS_AN(param,AN_SCALE))
                {
                    rs_sample_cst(s_tmp,1.0);
                    rs_sample_set_subsamp(s_x[i_a],s_tmp,ens_ind,0,ens_ind,0);
                }
                /* mass for QED FV model */
                if (param->with_qed_fvol)
                {
                    sprintf(sf_name,"%s/%s_%s_%s.boot%s",ens_pt->dir,M_str,\
                            param->qed_fvol_mass_name,ens_pt->dataset,ext);
                    rs_sample_load_subsamp(s_tmp,sf_name,0,0,0,0);
                    rs_sample_set_subsamp(s_x[i_fvM],s_tmp,ens_ind,0,\
                                          ens_ind,0);
                }
            }
        }
        if ((param->verb > 0)&&(param->nproc == 1))
        {
            printf("[");
            for (i=0;i<60*(ens_ind+1)/nens;i++)
            {
                printf("=");
            }
            for (i=60*(ens_ind+1)/nens;i<60;i++)
            {
                printf(" ");
            }
            printf("]  %d/%d\r",(int)ens_ind+1,(int)nens);
            fflush(stdout);
        }
    }
    for (i=0;i<param->nbeta;i++)
    for (j=0;j<param->nvol[i];j++)
    {
        rs_sample_subsamp_eqmuls(param->s_vol_av[i],                   \
                                 1.0/((double)(param->nenspvol[i][j])),\
                                 j,0,j,0);
    }
    
    if (IS_AN(param,AN_PHYPT)&&!IS_AN(param,AN_SCALE)&&\
        (strbufcmp(param->with_ext_a,"") != 0))
    {
        rs_sample_varp(param->a_err,param->s_a);
        mat_eqsqrt(param->a_err);
    }
    printf("\n");
    
    rs_sample_destroy(s_tmp);
}
