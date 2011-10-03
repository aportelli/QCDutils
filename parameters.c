#include "models.h"
#include "parameters.h"
#include <stdio.h>
#include <string.h>
#include <latan/latan_io.h>
#include <latan/latan_mass.h>


#define GET_PARAM_I(param,name)\
if (strcmp(field[0],#name) == 0)\
{\
    (param)->name = atoi(field[1]);\
}
#define GET_PARAM_D(param,name)\
if (strcmp(field[0],#name) == 0)\
{\
    (param)->name = atof(field[1]);\
}
#define GET_PARAM_S(param,name)\
if (strcmp(field[0],#name) == 0)\
{\
    strbufcpy((param)->name,field[1]);\
}

void parse_ex_param(ex_param *param, const strbuf fname)
{
    strbuf *field;
    int nf,lc;
    
    field = NULL;
    
    param->Msq_pi_deg  = 0;
    param->Msq_K_deg   = 0;
    param->a_deg       = 0;
    param->with_umd    = 0;
    param->with_fvol   = 0;
    param->q_dim       = 0;
    param->verb        = 0;
    param->beta        = NULL;
    param->nbeta       = 0;
    param->init_param  = NULL;
    param->ninit_param = 0;
    param->nens        = 0;
    param->Mpi_cut     = -1.0;
    param->MpiL_cut    = -1.0;
    strbufcpy(param->analyze,"");
    strbufcpy(param->q_name,"");
    strbufcpy(param->scale_part,"");
    strbufcpy(param->pi_name,"");
    strbufcpy(param->K_name,"");
    strbufcpy(param->suffix,"");
    strbufcpy(param->manifest,"");
    BEGIN_FOR_LINE_TOK(field,fname," \t",nf,lc)
    {
        if (field[0][0] != '#')
        {
            GET_PARAM_I(param,Msq_pi_deg);
            GET_PARAM_I(param,Msq_K_deg);
            GET_PARAM_I(param,a_deg);
            GET_PARAM_I(param,with_umd);
            GET_PARAM_I(param,with_fvol);
            GET_PARAM_I(param,q_dim);
            GET_PARAM_I(param,verb);
            GET_PARAM_D(param,Mpi_cut);
            GET_PARAM_D(param,MpiL_cut);
            GET_PARAM_S(param,analyze);
            GET_PARAM_S(param,q_name);
            GET_PARAM_S(param,scale_part);
            GET_PARAM_S(param,pi_name);
            GET_PARAM_S(param,K_name);
            GET_PARAM_S(param,suffix);
            GET_PARAM_S(param,manifest);
            if ((strcmp(field[0],"init_param") == 0)&&(nf >= 3))
            {
                param->ninit_param++;
                param->init_param = (fit_param *)realloc(param->init_param,    \
                                         param->ninit_param*sizeof(fit_param));
                param->init_param[param->ninit_param-1].ind   =\
                    (size_t)atoi(field[1]);
                param->init_param[param->ninit_param-1].value = atof(field[2]);
            }
        }
    }
    END_FOR_LINE_TOK(field);
    BEGIN_FOR_LINE_TOK(field,param->manifest,"_",nf,lc)
    {
        add_beta(param,field[2]);
        param->nens++;
    }
    END_FOR_LINE_TOK(field);
    if (strcmp(param->analyze,"phypt") == 0)
    {
        param->ex_dim  = 2;
        param->model   = &fm_phyptfit_ex_taylor;
    }
    else if (strcmp(param->analyze,"scaleset") == 0)
    {
        param->ex_dim  = 0;
        param->q_dim   = 0;
        param->a_deg   = 0;
        param->model   = &fm_scaleset_ex_taylor;
        sprintf(param->q_name,"Msq_%s",param->scale_part);
        get_mass(param->M_scale,param->scale_part);
    }
    else
    {
        fprintf(stderr,"error: analysis program %s unknown\n",param->analyze);
        abort();
    }
}

void add_beta(ex_param *param, const strbuf new_beta)
{
    if (ind_beta(new_beta,param) < 0)
    {
        param->nbeta++;
        param->beta = (strbuf *)realloc(param->beta,               \
                                        param->nbeta*sizeof(strbuf));
        strbufcpy(param->beta[param->nbeta-1],new_beta);
    }
}

int ind_beta(const strbuf beta, const ex_param *param)
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


