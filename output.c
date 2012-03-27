#include "output.h"
#include <math.h>
#include <string.h>
#include <latan/latan_math.h>
#include <latan/latan_plot.h>
#include "data_loader.h"
#include "models.h"

#define PLOT_ADD_FIT(obj,kx,ky,title,color)\
plot_add_fit(p[kx],d,ky,phy_pt,kx,fit,x_range[kx][0],x_range[kx][1],1000,true,\
             obj,title,"",color,color)

#define PLOT_ADD_EX(kx,s)\
if ((!latan_isnan(param->q_target[0]))&&(!latan_isnan(param->q_target[1])))\
{\
    plot_add_datpoint(p[kx],mat_get(phy_pt,kx,0),param->q_target[0],      \
                      -1.0,param->q_target[1],"target","rgb 'dark-blue'");\
}\
plot_add_datpoint(p[kx],mat_get(phy_pt,kx,0),mat_get(fit,s,0),     \
                  -1.0,sqrt(mat_get(fit_var,s,0)),"physical point",\
                  "rgb 'black'");

#define PLOT_DISP(kx)\
plot_set_title(p[kx],gtitle);\
plot_set_xlabel(p[kx],xlabel);\
plot_set_ylabel(p[kx],ylabel);\
plot_disp(p[kx]);


void plot_fit(const mat *fit, const mat *fit_var, fit_data *d,\
              fit_param *param, const plot_flag f)
{
    plot *p[N_EX_VAR];
    double *xb[N_EX_VAR] = {NULL,NULL,NULL,NULL,NULL,NULL};
    double x_range[N_EX_VAR][2],b_int[2],dbind,a;
    size_t bind,k,phy_ind,s;
    strbuf color,gtitle,title,xlabel,ylabel;
    mat *phy_pt,*x_k;
    
    phy_pt = mat_create(N_EX_VAR,1);
    x_k    = mat_create(param->nens,1);
    for (k=0;k<N_EX_VAR;k++)
    {
        p[k] = plot_create();
    }
    
    param->plotting = 1;
    if (IS_AN(param,AN_PHYPT)&&IS_AN(param,AN_SCALE))
    {
        phy_ind = 1;
        s       = fm_scaleset_taylor_npar(param);
    }
    else
    {
        phy_ind = 0;
        s       = 0;
    }
    for (k=0;k<N_EX_VAR;k++)
    {
        fit_data_get_x_k(x_k,d,k);
        if ((k == i_a)||(k == i_ud)||(k == i_Linv))
        {
            x_range[k][0] = 0.0;
        }
        else
        {
            x_range[k][0] = mat_get_min(x_k)-0.40*fabs(mat_get_min(x_k));
        }
        x_range[k][1] = mat_get_max(x_k)+0.40*fabs(mat_get_min(x_k));
    }
    if (f == Q)
    {
        sprintf(gtitle,"quantity: %s -- scale: %s -- datasets: %s -- ensembles: %s",
                param->q_name,param->scale_part,param->dataset_cat,param->manifest);
        mat_set(phy_pt,i_ud,0,SQ(param->M_ud));
        mat_set(phy_pt,i_s,0,SQ(param->M_s));
        mat_set(phy_pt,i_umd,0,param->M_umd);
        mat_set(phy_pt,i_bind,0,0.0);
        mat_set(phy_pt,i_a,0,0.0);
        mat_set(phy_pt,i_Linv,0,0.0);
        for (bind=0;bind<param->nbeta;bind++)
        {
            dbind      = (double)(bind);
            b_int[0]   = dbind - 0.1;
            b_int[1]   = dbind + 0.1;
            xb[i_bind] = b_int;
            fit_data_fit_region(d,xb);
            sprintf(color,"%d",1+(int)bind);
            sprintf(title,"beta = %s",param->beta[bind]);
            PLOT_ADD_FIT(PF_DATA,i_ud,phy_ind,title,color);
            PLOT_ADD_FIT(PF_DATA,i_s,phy_ind,title,color);
            PLOT_ADD_FIT(PF_DATA,i_umd,phy_ind,title,color);
            PLOT_ADD_FIT(PF_DATA,i_a,phy_ind,title,color);
            PLOT_ADD_FIT(PF_DATA,i_Linv,phy_ind,title,color);
            fit_data_fit_all_points(d,true);
        }
        PLOT_ADD_FIT(PF_FIT,i_ud,phy_ind,"","rgb 'black'");
        PLOT_ADD_FIT(PF_FIT,i_s,phy_ind,"","rgb 'black'");
        PLOT_ADD_FIT(PF_FIT,i_umd,phy_ind,"","rgb 'black'");
        PLOT_ADD_FIT(PF_FIT,i_a,phy_ind,"","rgb 'black'");
        PLOT_ADD_FIT(PF_FIT,i_Linv,phy_ind,"","rgb 'black'");
        switch (param->q_dim) 
        {
            case 0:
                strbufcpy(ylabel,param->q_name);
                break;
            case 1:
                sprintf(ylabel,"%s (MeV)",param->q_name);
                break;
            default:
                sprintf(ylabel,"%s^%d (MeV^%d)",param->q_name,param->q_dim,\
                        param->q_dim);
                break;
        }
        
        sprintf(xlabel,"M_%s^2 (MeV^2)",param->ud_name);
        PLOT_ADD_EX(i_ud,s);
        PLOT_DISP(i_ud);
        sprintf(xlabel,"M_%s^2 (MeV^2)",param->s_name);
        PLOT_ADD_EX(i_s,s);
        PLOT_DISP(i_s);
        strbufcpy(xlabel,"a (MeV^-1)");
        PLOT_ADD_EX(i_a,s);
        PLOT_DISP(i_a);
        if (param->have_umd)
        {
            sprintf(xlabel,"%s (MeV^2)",param->umd_name);
            PLOT_ADD_EX(i_umd,s);
            PLOT_DISP(i_umd);
        }
        strbufcpy(xlabel,"1/L (MeV)");
        PLOT_ADD_EX(i_Linv,s);
        PLOT_DISP(i_Linv);
    }
    else if (f == SCALE)
    {
        sprintf(gtitle,"scale setting: %s -- datasets: %s -- ensembles: %s",
                param->scale_part,param->dataset_cat,param->manifest);
        for (bind=0;bind<param->nbeta;bind++)
        {
            dbind      = (double)(bind);
            b_int[0]   = dbind - 0.1;
            b_int[1]   = dbind + 0.1;
            xb[i_bind] = b_int;
            a          = mat_get(fit,bind,0);
            fit_data_fit_region(d,xb);
            mat_set(phy_pt,i_ud,0,SQ(a*param->M_ud));
            mat_set(phy_pt,i_s,0,SQ(a*param->M_s));
            mat_set(phy_pt,i_umd,0,SQ(a)*param->M_umd);
            mat_set(phy_pt,i_bind,0,bind);
            mat_set(phy_pt,i_a,0,a);
            mat_set(phy_pt,i_Linv,0,0.0);
            sprintf(color,"%d",1+(int)bind);
            sprintf(title,"beta = %s",param->beta[bind]);
            PLOT_ADD_FIT(PF_DATA|PF_FIT,i_ud,0,title,color);
            PLOT_ADD_FIT(PF_DATA|PF_FIT,i_s,0,title,color);
            PLOT_ADD_FIT(PF_DATA|PF_FIT,i_umd,0,title,color);
            PLOT_ADD_FIT(PF_DATA|PF_FIT,i_Linv,0,title,color);
            fit_data_fit_all_points(d,true);
        }
        sprintf(ylabel,"(a*M_%s)^2",param->scale_part);
        sprintf(xlabel,"(a*M_%s)^2",param->ud_name);
        PLOT_DISP(i_ud);
        sprintf(xlabel,"(a*M_%s)^2",param->s_name);
        PLOT_DISP(i_s);
        if (param->have_umd)
        {
            sprintf(xlabel,"a^2*%s",param->umd_name);
            PLOT_DISP(i_umd);
        }
        strbufcpy(xlabel,"a/L");
        PLOT_DISP(i_Linv);
    }
    param->plotting = 0;
    
    mat_destroy(phy_pt);
    mat_destroy(x_k);
    for (k=0;k<N_EX_VAR;k++)
    {
        plot_destroy(p[k]);
    }
}
#undef ADD_PLOT

void plot_chi2_comp(const fit_data *d, const fit_param *param, const size_t k,\
                    const strbuf title)
{
    strbuf *tmpfname,plotcmd;
    size_t nbeta,nens,nydim,bind;
    size_t i;
    FILE **tmpf;
    plot *p;
    
    nbeta    = param->nbeta;
    nens     = param->nens;
    nydim    = fit_data_get_nydim(d);
    
    tmpf     = (FILE **)malloc(nbeta*nydim*sizeof(FILE *));
    tmpfname = (strbuf *)malloc(nbeta*nydim*sizeof(strbuf));
    p        = plot_create();
    
    for (i=0;i<nbeta;i++)
    {
        sprintf(tmpfname[i],".qcd_phyfit_tmp_%d",(int)i);
        tmpf[i] = fopen(tmpfname[i],"w");
    }
    for (i=0;i<nens;i++)
    {
        bind = ind_beta(param->point[i].beta,param);
        fprintf(tmpf[bind],"%s %e %e\n",param->point[i].dir,(double)(i),\
                mat_get(d->chi2_comp,i+k*nens,0));
    }
    for (i=0;i<nbeta;i++)
    {
        fclose(tmpf[i]);
        sprintf(plotcmd,"'%s' u 2:3:xtic(1) t '%s' w impulse",tmpfname[i],\
                param->beta[i]);
        plot_add_plot(p,plotcmd);
    }
    plot_add_head(p,"set xtics rotate by -90 font 'courier, 10'");
    plot_set_scale_manual(p,-1.0,(double)(param->nens+1),-5.0,5.0);
    plot_add_plot(p,"0.0 lt -1 lc rgb 'black' notitle");
    plot_add_plot(p,"1.0 lt -1 lc rgb 'black' notitle");
    plot_add_plot(p,"-1.0 lt -1 lc rgb 'black' notitle");
    plot_add_plot(p,"2.0 lt -1 lc rgb 'dark-gray' notitle");
    plot_add_plot(p,"-2.0 lt -1 lc rgb 'dark-gray' notitle");
    plot_add_plot(p,"3.0 lt -1 lc rgb 'gray' notitle");
    plot_add_plot(p,"-3.0 lt -1 lc rgb 'gray' notitle");
    plot_add_plot(p,"4.0 lt -1 lc rgb 'light-gray' notitle");
    plot_add_plot(p,"-4.0 lt -1 lc rgb 'light-gray' notitle");
    plot_set_ylabel(p,"standard deviations");
    plot_set_title(p,title);
    plot_disp(p);
    for (i=0;i<nbeta;i++)
    {
        remove(tmpfname[i]);
    }
    
    free(tmpf);
    free(tmpfname);
    plot_destroy(p);
}

#define PRINT_PAR(name)\
{\
    printf("%12s = % e (%4.0f%%)\n",name,mat_get(fit,i,0),            \
           sqrt(mat_get(fit_var,i,0))/fabs(mat_get(fit,i,0))*100.0);\
    i++;\
}
#define PRINT_SCALE(name)\
{\
    double sig_,p_;\
    sig_ = sqrt(mat_get(fit_var,i,0));\
    p_   = mat_get(fit,i,0);\
    printf("%12s = % e (%4.1f%%) [%12s^-1 = %5.0f(%4.0f) MeV]\n",name,\
           p_,sig_/fabs(p_)*100.0,name,1.0/p_,sig_/SQ(p_));\
    i++;\
}

void print_result(const rs_sample *s_fit, fit_param *param)
{
    size_t i;
    int j;
    mat *fit,*fit_var;
    strbuf buf;
    
    i   = 0;
    fit = rs_sample_pt_cent_val(s_fit);
    
    fit_var = mat_create(fit_model_get_npar(&(param->fm),param),1);
    
    rs_sample_varp(fit_var,s_fit);
    printf("\nfit parameters :\n");
    if (IS_AN(param,AN_SCALE))
    {
        for(j=0;j<(int)param->nbeta;j++)
        {
            sprintf(buf,"a_%s",param->beta[j]);
            PRINT_SCALE(buf);
        }
        if (param->s_M_ud_deg > 0)
        {
            for (j=0;j<param->s_M_ud_deg;j++)
            {
                sprintf(buf,"s_p_ud_%d",j+1);
                PRINT_PAR(buf);
            }
        }
        if (param->s_M_s_deg > 0)
        {
            for (j=0;j<param->s_M_s_deg;j++)
            {
                sprintf(buf,"s_p_s_%d",j+1);
                PRINT_PAR(buf);
            }
        }
        if (param->s_with_a2M_ud)
        {
            PRINT_PAR("s_a2M_ud");
        }
        if (param->s_with_a2M_s)
        {
            PRINT_PAR("s_a2M_s");
        }
        if (param->s_umd_deg > 0)
        {
            for (j=0;j<param->s_umd_deg;j++)
            {
                sprintf(buf,"s_p_umd_%d",j+1);
                PRINT_PAR(buf);
            }
        }
        if (param->s_with_qed_fvol)
        {
            PRINT_PAR("s_p_fvol_L");
        }
    }
    if (IS_AN(param,AN_PHYPT))
    {
        PRINT_PAR(param->q_name);
        if (param->M_ud_deg > 0)
        {
            for (j=0;j<param->M_ud_deg;j++)
            {
                sprintf(buf,"p_ud_%d",j+1);
                PRINT_PAR(buf);
            }
        }
        if (param->M_s_deg > 0)
        {
            for (j=0;j<param->M_s_deg;j++)
            {
                sprintf(buf,"p_s_%d",j+1);
                PRINT_PAR(buf);
            }
        }
        if (param->a_deg > 0)
        {
            PRINT_PAR("p_a");
        }
        if (param->with_a2M_ud)
        {
            PRINT_PAR("p_a2M_ud");
        }
        if (param->with_a2M_s)
        {
            PRINT_PAR("p_a2M_s");
        }
        if (param->umd_deg > 0)
        {
            for (j=0;j<param->umd_deg;j++)
            {
                sprintf(buf,"p_umd_%d",j+1);
                PRINT_PAR(buf);
            }
        }
        if (param->with_qed_fvol)
        {
            PRINT_PAR("p_fvol_L");
        }
        if (IS_AN(param,AN_PHYPT)&&!IS_AN(param,AN_SCALE)&&(param->with_ext_a))
        {
            for(j=0;j<(int)param->nbeta;j++)
            {
                sprintf(buf,"corr_a_%s",param->beta[j]);
                PRINT_SCALE(buf);
            }
        }
    }
    printf("\n");
    
    mat_destroy(fit_var);
}
#undef PRINT_PAR

#define PRINT_DLABEL(name) fprintf(stream,"%-12s  ",name)
#define PRINT_DLABEL_WERR(name) fprintf(stream,"%-12s %-12s  ",name,"error")
#define PRINT_D(value) fprintf(stream,"% .5e  ",value)
#define PRINT_D_WERR(value,err) fprintf(stream,"% .5e % .5e  ",value,err)
#define PRINT_X(ind,dim)\
PRINT_D(mat_get(rs_sample_pt_cent_val(s_x[ind]),ens_ind,0))
#define PRINT_CV_WERR(s,s_err)\
PRINT_D_WERR(mat_get(rs_sample_pt_cent_val(s),ens_ind,0),\
             mat_get(s_err,ens_ind,0))
#define PRINT_X_WERR(ind) PRINT_CV_WERR(s_x[ind],x_err[ind])
void fprint_table(FILE* stream, rs_sample *s_x[N_EX_VAR], rs_sample *s_q[2],\
                  const mat *res, const fit_param *param, const plot_flag f)
{
    mat **x_err,*q_err;
    const rs_sample *s_q_pt;
    size_t nens,ydim;
    size_t ens_ind;
    
    nens     = param->nens;
    ydim     = (f == SCALE) ? 0   : 1;
    s_q_pt   = s_q[ydim];
    
    x_err = mat_ar_create(N_EX_VAR,nens,1);
    q_err = mat_create(nens,1);
    
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
    rs_sample_varp(q_err,s_q_pt);
    mat_eqsqrt(q_err);
    
    /* display */
    fprintf(stream,"#%29s  ","ensemble");
    PRINT_DLABEL_WERR(param->ud_name);
    PRINT_DLABEL_WERR(param->s_name);
    if (f == Q)
    {
        PRINT_DLABEL_WERR("a");
    }
    if (param->have_umd)
    {
        PRINT_DLABEL_WERR(param->umd_name);
    }
    PRINT_DLABEL_WERR("1/L");
    if (f == Q)
    {
        PRINT_DLABEL_WERR(param->q_name);
    }
    else if (f == SCALE)
    {
        PRINT_DLABEL_WERR(param->scale_part);
    }
    PRINT_DLABEL("residuals");
    fprintf(stream,"\n");
    for(ens_ind=0;ens_ind<nens;ens_ind++)
    {
        fprintf(stream,"%30s  ",param->point[ens_ind].dir);
        PRINT_X_WERR(i_ud);
        PRINT_X_WERR(i_s);
        if (f == Q)
        {
            PRINT_X_WERR(i_a);
        }
        if (param->have_umd)
        {
            PRINT_X_WERR(i_umd);
        }
        PRINT_X_WERR(i_Linv);
        PRINT_CV_WERR(s_q_pt,q_err);
        PRINT_D(mat_get(res,ens_ind,0));
        fprintf(stream,"\n");
    }
    fprintf(stream,"\n");

    mat_ar_destroy(x_err,N_EX_VAR);
    mat_destroy(q_err);
}

#undef PRINT_DLABEL
#undef PRINT_DLABEL_WERR
#undef PRINT_D
#undef PRINT_D_WERR
#undef PRINT_X
#undef PRINT_X_WERR
#undef PRINT_Q_WERR
