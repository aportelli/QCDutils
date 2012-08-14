#define _POSIX_C_SOURCE 199506L /* strtok_r is used here */

#include "parameters.h"
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <latan/latan_io.h>
#include <latan/latan_mass.h>
#include "models.h"
#include "data_loader.h"

#define ATOF(str) (strtod(str,(char **)NULL))
#define ATOI(str) ((int)strtol(str,(char **)NULL,10))

#define GET_PARAM_I(param,name)\
if (strbufcmp(field[0],#name) == 0)\
{\
    (param)->name = ATOI(field[1]);\
    continue;\
}
#define GET_PARAM_D(param,name)\
if (strbufcmp(field[0],#name) == 0)\
{\
    (param)->name = ATOF(field[1]);\
    continue;\
}
#define GET_PARAM_S(param,name)\
if (strbufcmp(field[0],#name) == 0)\
{\
    strbufcpy((param)->name,field[1]);\
    continue;\
}
#define CHECK_MODEL(param,m,s_m)\
if (IS_AN(param,AN_PHYPT))\
{\
    if (strbufcmp((param)->model,#m) == 0)\
    {\
        (param)->fm.func[s] = &fm_phypt_##m##_func;\
        strbufcat((param)->fm.name,#m);\
        if (IS_AN(param,AN_SCALE))\
        {\
            (param)->fm.npar = &fm_comb_phypt_##m##_scaleset_##s_m##_npar;\
            strbufcat((param)->fm.name,"/");\
        }\
        else\
        {\
            (param)->fm.npar = &fm_phypt_##m##_npar;\
        }\
    }\
}\
if (IS_AN(param,AN_SCALE))\
{\
    if (strbufcmp((param)->s_model,#s_m) == 0)\
    {\
        (param)->fm.func[0] = &fm_scaleset_##s_m##_func;\
        strbufcat((param)->fm.name,#s_m);\
        if (!IS_AN(param,AN_PHYPT))\
        {\
            (param)->fm.npar = &fm_scaleset_##s_m##_npar;\
        }\
    }\
}

static int ind_dataset(const strbuf dataset, const fit_param *param);
static void add_dataset(fit_param *param, const strbuf new_dataset);
static void add_beta(fit_param *param, const strbuf new_beta);

static int ind_dataset(const strbuf dataset, const fit_param *param)
{
    size_t i;
    
    for (i=0;i<param->ndataset;i++)
    {
        if (strbufcmp(dataset,param->dataset[i]) == 0)
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
        if (strbufcmp(beta,param->beta[i]) == 0)
        {
            return (int)i;
        }
    }
    
    return -1;
}

fit_param * fit_param_parse(const strbuf fname)
{
    fit_param *param;
    strbuf *field,ens_dir,test_fname,name;
    int nf,lc;
    int i;
    size_t s;
    size_t j;
    double dbuf[2];
    
    field = NULL;
    param = (fit_param *)malloc(sizeof(fit_param));
    
    /* initialization */
    param->analyze_flag      = AN_NOTHING;
    param->with_const        = 0;
    param->M_ud              = latan_nan();
    param->M_ud_deg          = 0;
    param->s_M_ud_deg        = 0;
    param->M_s               = latan_nan();
    param->M_s_deg           = 0;
    param->s_M_s_deg         = 0;
    param->M_umd_val         = latan_nan();
    param->alpha             = latan_nan();
    param->M_scale           = latan_nan();
    param->a_deg             = 0;
    param->with_aalpha       = 0;
    param->with_a2M_ud       = 0;
    param->s_with_a2M_ud     = 0;
    param->with_a2M_s        = 0;
    param->s_with_a2M_s      = 0;
    param->umd_deg           = 0;
    param->s_umd_deg         = 0;
    param->have_umd          = 0;
    param->with_udumd        = 0;
    param->with_sumd         = 0;
    param->alpha_deg         = 0;
    param->s_alpha_deg       = 0;
    param->have_alpha        = 0;
    param->with_udalpha      = 0;
    param->with_salpha       = 0;
    param->with_qed_fvol     = 0;
    param->s_with_qed_fvol   = 0;
    param->with_ext_a        = 0;
    param->with_ext_M_umd    = 0;
    param->q_dim             = 0;
    param->q_target[0]       = latan_nan();
    param->q_target[1]       = latan_nan();
    param->verb              = 0;
    param->correlated        = 0;
    param->save_result       = 0;
    param->plot              = 0;
    param->warn_missing_data = 0;
    param->scale_model       = 0;
    param->dataset           = NULL;
    param->ndataset          = 0;
    param->beta              = NULL;
    param->nbeta             = 0;
    param->init_param        = NULL;
    param->ninit_param       = 0;
    param->save_param        = NULL;
    param->nsave_param       = 0;
    param->point             = NULL;
    param->nens              = 0;
    param->nsample           = 0;
    strbufcpy(param->analyze,"");
    strbufcpy(param->model,"");
    strbufcpy(param->s_model,"");
    strbufcpy(param->q_name,"");
    strbufcpy(param->scale_part,"");
    strbufcpy(param->ud_name,"");
    strbufcpy(param->s_name,"");
    strbufcpy(param->umd_name,"");
    strbufcpy(param->manifest,"");
    strbufcpy(param->s_manifest,"");
    strbufcpy(param->dataset_cat,"");
    strbufcpy(param->save_plot,"");
    
    /* parse parameter file */
    BEGIN_FOR_LINE_TOK(field,fname," \t",nf,lc)
    {
        if (field[0][0] != '#')
        {
            GET_PARAM_I(param,with_const);
            GET_PARAM_I(param,M_ud_deg);
            GET_PARAM_I(param,s_M_ud_deg);
            GET_PARAM_D(param,M_ud);
            GET_PARAM_I(param,M_s_deg);
            GET_PARAM_I(param,s_M_s_deg);
            GET_PARAM_D(param,M_s);
            GET_PARAM_S(param,M_umd);
            GET_PARAM_D(param,M_scale);
            GET_PARAM_D(param,alpha);
            GET_PARAM_D(param,with_aalpha);
            GET_PARAM_I(param,a_deg);
            GET_PARAM_I(param,with_a2M_ud);
            GET_PARAM_I(param,s_with_a2M_ud);
            GET_PARAM_I(param,with_a2M_s);
            GET_PARAM_I(param,s_with_a2M_s);
            GET_PARAM_I(param,umd_deg);
            GET_PARAM_I(param,s_umd_deg);
            GET_PARAM_I(param,with_udumd);
            GET_PARAM_I(param,with_sumd);
            GET_PARAM_I(param,alpha_deg);
            GET_PARAM_I(param,s_alpha_deg);
            GET_PARAM_I(param,with_udalpha);
            GET_PARAM_I(param,with_salpha);
            GET_PARAM_I(param,with_qed_fvol);
            GET_PARAM_I(param,s_with_qed_fvol);
            GET_PARAM_I(param,with_ext_a);
            GET_PARAM_I(param,with_ext_M_umd);
            GET_PARAM_I(param,q_dim);
            GET_PARAM_I(param,verb);
            GET_PARAM_I(param,correlated);
            GET_PARAM_I(param,save_result);
            GET_PARAM_I(param,plot);
            GET_PARAM_I(param,warn_missing_data);
            GET_PARAM_S(param,analyze);
            GET_PARAM_S(param,model);
            GET_PARAM_S(param,s_model);
            GET_PARAM_S(param,q_name);
            GET_PARAM_S(param,s_manifest);
            GET_PARAM_S(param,scale_part);
            GET_PARAM_S(param,ud_name);
            GET_PARAM_S(param,s_name);
            GET_PARAM_S(param,umd_name);
            GET_PARAM_S(param,manifest);
            GET_PARAM_S(param,save_plot);
            if ((strbufcmp(field[0],"q_target") == 0)&&(nf >= 2))
            {
                param->q_target[0] = ATOF(field[1]);
                param->q_target[1] = ATOF(field[2]);
                continue;
            }
            if ((strbufcmp(field[0],"dataset") == 0)&&(nf >= 2))
            {
                for (i=1;i<nf;i++)
                {
                    if (field[i][0] == '#')
                    {
                        break;
                    }
                    add_dataset(param,field[i]);
                }
                continue;
            }
            if ((strbufcmp(field[0],"init_param") == 0)&&(nf >= 3))
            {
                param->ninit_param++;
                param->init_param = (fit_init *)realloc(param->init_param,    \
                                         param->ninit_param*sizeof(fit_init));
                param->init_param[param->ninit_param-1].ind   =\
                    (size_t)ATOI(field[1]);
                param->init_param[param->ninit_param-1].value = ATOF(field[2]);
                continue;
            }
            if ((strbufcmp(field[0],"save_param") == 0)&&(nf >= 2))
            {
                param->nsave_param++;
                param->save_param = (size_t *)realloc(param->save_param,    \
                                            param->nsave_param*sizeof(size_t));
                param->save_param[param->nsave_param-1]=(size_t)ATOI(field[1]);
                continue;
            }
            fprintf(stderr,"warning: line %d of %s ignored (syntax error)\n",\
                    lc,fname);
        }
    }
    END_FOR_LINE_TOK(field);
    
    /* set analyze program flag */
    if (strbufcmp(param->analyze,"phypt") == 0) 
    {
        param->analyze_flag = AN_PHYPT;
    }
    else if (strbufcmp(param->analyze,"scaleset") == 0)
    {
        param->analyze_flag = AN_SCALE;
    }
    else if (strbufcmp(param->analyze,"comb_phypt_scaleset") == 0)
    {
        param->analyze_flag = AN_PHYPT|AN_SCALE;
    }
    else
    {
        fprintf(stderr,"error: analyze program '%s' unknown\n",param->analyze);
        exit(EXIT_FAILURE);
    }
    
    /* choose the model */
    s               = (IS_AN(param,AN_PHYPT)&&IS_AN(param,AN_SCALE)) ? 1 : 0;
    param->fm.nxdim = N_EX_VAR;
    param->fm.nydim = (IS_AN(param,AN_PHYPT)&&IS_AN(param,AN_SCALE)) ? 2 : 1;
    strbufcpy(param->fm.name,"");
    CHECK_MODEL(param,taylor,taylor);
    if ((param->fm.func[0] == NULL)||((s==1)&&(param->fm.func[s] == NULL))
        ||(param->fm.npar == NULL))
    {
        fprintf(stderr,"error: model initialization failed\n");
        exit(EXIT_FAILURE);
    }

    /* set extrapolation masses */
    if (latan_isnan(param->M_ud))
    {
        get_mass(dbuf,param->ud_name);
        param->M_ud = dbuf[0];
    }
    if (latan_isnan(param->M_s))
    {
        get_mass(dbuf,param->s_name);
        param->M_s = dbuf[0];
    }
    if (latan_isnan(param->M_scale))
    {
        get_mass(dbuf,param->scale_part);
        param->M_scale = dbuf[0];
    }
    
    /* parse ensemble informations */
    if (IS_AN(param,AN_SCALE))
    {
        if (!IS_AN(param,AN_PHYPT)&&(strlen(param->manifest) == 0))
        {
            strbufcpy(param->manifest,param->s_manifest);
        }
        else
        {
            strbufcpy(param->s_manifest,param->manifest);
        }
        if (!IS_AN(param,AN_PHYPT))
        {
            sprintf(param->q_name,"M_%s",param->scale_part);
        }
    }
    BEGIN_FOR_LINE_TOK(field,param->manifest,"_ \t",nf,lc)
    {
        if ((nf>0)&&(field[0][0] != '#'))
        {
            strbufcpy(ens_dir,_line);
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
                        rs_sample_load(NULL,&(param->nsample),NULL,test_fname);
                    }
                }
                else if (param->warn_missing_data)
                {
                    fprintf(stderr,"warning: no data found for dataset %s in ensemble %s (unable to read file %s)\n",\
                            param->dataset[j],ens_dir,test_fname);
                }
            }
        }
    }
    END_FOR_LINE_TOK(field);
    
    /* create samples */
    if (param->with_ext_a)
    {
        param->s_a    = rs_sample_create(param->nbeta,1,param->nsample);
        param->a_err  = mat_create(param->nbeta,1);
    }
    if (param->with_ext_M_umd)
    {
        param->s_M_umd = rs_sample_create(1,1,param->nsample);
    }
    param->s_ex   = rs_sample_create(1,1,param->nsample);
    param->ex_err = mat_create(1,1);
    
    /* take care of M_umd input */
    if (param->with_ext_M_umd)
    {
        sprintf(name,"%s_%s_%s_%s_%s.boot",param->M_umd,param->scale_part,\
                param->s_manifest,param->dataset_cat,param->manifest);
        rs_sample_load(param->s_M_umd,NULL,NULL,name);
        param->M_umd_val = mat_get(rs_sample_pt_cent_val(param->s_M_umd),0,0);
    }
    else
    {
        param->M_umd_val = ATOF(param->M_umd);
    }
    
    return param;
}

void fit_param_destroy(fit_param *param)
{
    free(param->init_param);
    free(param->save_param);
    free(param->point);
    free(param->dataset);
    free(param->beta);
    rs_sample_destroy(param->s_a);
    mat_destroy(param->a_err);
    rs_sample_destroy(param->s_M_umd);
    rs_sample_destroy(param->s_ex);
    mat_destroy(param->ex_err);
    free(param);
}
