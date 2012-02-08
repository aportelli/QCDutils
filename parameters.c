#define _POSIX_C_SOURCE 199506L /* strtok_r is used here */

#include "models.h"
#include "parameters.h"
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <latan/latan_io.h>
#include <latan/latan_mass.h>

#define ATOF(str) (strtod(str,(char **)NULL))
#define ATOI(str) ((int)strtol(str,(char **)NULL,10))

#define GET_PARAM_I(param,name)\
if (strcmp(field[0],#name) == 0)\
{\
    (param)->name = ATOI(field[1]);\
}
#define GET_PARAM_D(param,name)\
if (strcmp(field[0],#name) == 0)\
{\
    (param)->name = ATOF(field[1]);\
}
#define GET_PARAM_S(param,name)\
if (strcmp(field[0],#name) == 0)\
{\
    strbufcpy((param)->name,field[1]);\
}

static int ind_dataset(const strbuf dataset, const fit_param *param);
static void add_dataset(fit_param *param, const strbuf new_dataset);
static void add_beta(fit_param *param, const strbuf new_beta);

static int ind_dataset(const strbuf dataset, const fit_param *param)
{
    size_t i;
    
    for (i=0;i<param->ndataset;i++)
    {
        if (strcmp(dataset,param->dataset[i]) == 0)
        {
            return (int)i;
        }
    }
    
    return -1;
}

static void add_dataset(fit_param *param, const strbuf new_dataset)
{
    size_t max_len;
    
    if (ind_dataset(new_dataset,param) < 0)
    {
        param->ndataset++;
        param->dataset = (strbuf *)realloc(param->dataset,               \
                                           param->ndataset*sizeof(strbuf));
        strbufcpy(param->dataset[param->ndataset-1],new_dataset);
        if (strlen(param->dataset_cat) > 0)
        {
            strncat(param->dataset_cat,"_",1);
        }
        max_len = STRING_LENGTH - strlen(param->dataset_cat) - 1;
        strncat(param->dataset_cat,new_dataset,max_len);
    }
}

static void add_beta(fit_param *param, const strbuf new_beta)
{
    if (ind_beta(new_beta,param) < 0)
    {
        param->nbeta++;
        param->beta = (strbuf *)realloc(param->beta,               \
                                        param->nbeta*sizeof(strbuf));
        strbufcpy(param->beta[param->nbeta-1],new_beta);
    }
}

int ind_beta(const strbuf beta, const fit_param *param)
{
    size_t i;
    
    for (i=0;i<param->nbeta;i++)
    {
        if (strcmp(beta,param->beta[i]) == 0)
        {
            return (int)i;
        }
    }
    
    return -1;
}

fit_param * fit_param_parse(const strbuf fname)
{
    fit_param *param;
    strbuf *field,ens_dir,test_fname;
    int nf,lc;
    int i;
    size_t j;
    double dbuf[2];
    
    field = NULL;
    param = (fit_param *)malloc(sizeof(fit_param));
    
    param->M_ud            = -1.0;
    param->M_ud_deg        = 0;
    param->s_M_ud_deg      = 0;
    param->M_s             = -1.0;
    param->M_s_deg         = 0;
    param->s_M_s_deg       = 0;
    param->M_scale         = -1.0;
    param->a_deg           = 0;
    param->with_umd        = 0;
    param->s_with_umd      = 0;
    param->have_umd        = 0;
    param->with_qed_fvol   = 0;
    param->s_with_qed_fvol = 0;
    param->with_ext_a      = 0;
    param->q_dim           = 0;
    param->q_target[0]     = -1.0;
    param->q_target[1]     = 0.0;
    param->verb            = 0;
    param->plotting        = 0;
    param->dataset         = NULL;
    param->ndataset        = 0;
    param->beta            = NULL;
    param->nbeta           = 0;
    param->init_param      = NULL;
    param->ninit_param     = 0;
    param->point           = NULL;
    param->nens            = 0;
    param->nsample         = 0;
    strbufcpy(param->analyze,"");
    strbufcpy(param->q_name,"");
    strbufcpy(param->scale_part,"");
    strbufcpy(param->ud_name,"");
    strbufcpy(param->s_name,"");
    strbufcpy(param->manifest,"");
    BEGIN_FOR_LINE_TOK(field,fname," \t",nf,lc)
    {
        if (field[0][0] != '#')
        {
            GET_PARAM_I(param,M_ud_deg);
            GET_PARAM_I(param,s_M_ud_deg);
            GET_PARAM_D(param,M_ud);
            GET_PARAM_I(param,M_s_deg);
            GET_PARAM_I(param,s_M_s_deg);
            GET_PARAM_D(param,M_s);
            GET_PARAM_D(param,M_scale);
            GET_PARAM_I(param,a_deg);
            GET_PARAM_I(param,with_umd);
            GET_PARAM_I(param,s_with_umd);
            GET_PARAM_I(param,with_qed_fvol);
            GET_PARAM_I(param,s_with_qed_fvol);
            GET_PARAM_I(param,with_ext_a);
            GET_PARAM_I(param,q_dim);
            GET_PARAM_I(param,verb);
            GET_PARAM_S(param,analyze);
            GET_PARAM_S(param,q_name);
            GET_PARAM_S(param,scale_part);
            GET_PARAM_S(param,ud_name);
            GET_PARAM_S(param,s_name);
            GET_PARAM_S(param,manifest);
            if ((strcmp(field[0],"q_target") == 0)&&(nf >= 2))
            {
                param->q_target[0] = ATOF(field[1]);
                param->q_target[1] = ATOF(field[2]);
            }
            if ((strcmp(field[0],"dataset") == 0)&&(nf >= 2))
            {
                for (i=1;i<nf;i++)
                {
                    add_dataset(param,field[i]);
                }
            }
            if ((strcmp(field[0],"init_param") == 0)&&(nf >= 3))
            {
                param->ninit_param++;
                param->init_param = (fit_init *)realloc(param->init_param,    \
                                         param->ninit_param*sizeof(fit_init));
                param->init_param[param->ninit_param-1].ind   =\
                    (size_t)ATOI(field[1]);
                param->init_param[param->ninit_param-1].value = ATOF(field[2]);
            }
        }
    }
    END_FOR_LINE_TOK(field);
    if (IS_ANALYZE(param,"phypt"))
    {
        param->model   = &fm_phypt_a_taylor;
    }
    else if (IS_ANALYZE(param,"scaleset"))
    {
        param->q_dim           = 0;
        param->s_M_ud_deg      = param->M_ud_deg;
        param->s_M_s_deg       = param->M_s_deg;
        param->s_with_umd      = param->with_umd;
        param->s_with_qed_fvol = param->with_qed_fvol;
        param->a_deg           = 0;
        param->with_ext_a      = 0;
        param->model           = &fm_scaleset_taylor;
        sprintf(param->q_name,"Msq_%s",param->scale_part);
    }
    else if (IS_ANALYZE(param,"comb_phypt_scale"))
    {
        param->model = &fm_comb_phyt_scale_taylor;
    }
    else
    {
        fprintf(stderr,"error: analysis program %s unknown\n",param->analyze);
        abort();
    }
    if (param->M_ud < 0)
    {
        get_mass(dbuf,param->ud_name);
        param->M_ud = dbuf[0];
    }
    if (param->M_s < 0)
    {
        get_mass(dbuf,param->s_name);
        param->M_s = dbuf[0];
    }
    if (param->M_scale < 0)
    {
        get_mass(dbuf,param->scale_part);
        param->M_scale = dbuf[0];
    }
    BEGIN_FOR_LINE_TOK(field,param->manifest,"_ \t",nf,lc)
    {
        if ((nf>0)&&(field[0][0] != '#'))
        {
            sprintf(ens_dir,"%s_%s_%s_%s_%s",field[0],field[1],field[2],field[3],\
                    field[4]);
            for (j=0;j<param->ndataset;j++)
            {
                sprintf(test_fname,"%s/%s_%s.boot%s",ens_dir,param->q_name,\
                        param->dataset[j],(io_get_fmt()==IO_XML)?".xml":"");
                if (access(test_fname,R_OK) == 0)
                {
                    param->nens++;
                    param->point  = (ens *)realloc(param->point,           \
                                                   param->nens*sizeof(ens));
                    param->point[param->nens-1].T = ATOI(field[0]);
                    param->point[param->nens-1].L = ATOI(field[1]);
                    strbufcpy(param->point[param->nens-1].beta,field[2]);
                    strbufcpy(param->point[param->nens-1].ud,field[3]);
                    strbufcpy(param->point[param->nens-1].s,field[4]);
                    strbufcpy(param->point[param->nens-1].dataset,\
                              param->dataset[j]);
                    strbufcpy(param->point[param->nens-1].dir,ens_dir);
                    add_beta(param,field[2]);
                    if (param->nsample == 0)
                    {
                        rs_sample_load_nsample(&(param->nsample),test_fname,"");
                    }
                }
                else
                {
                    fprintf(stderr,"warning: no data found for dataset %s in ensemble %s (unable to read file %s)\n",\
                            param->dataset[j],ens_dir,test_fname);
                }
            }
        }
    }
    END_FOR_LINE_TOK(field);
    
    return param;
}

void fit_param_destroy(fit_param *param)
{
    free(param->init_param);
    free(param->point);
    free(param->dataset);
    free(param->beta);
    free(param);
}
