#include "config.h"
#include "output.h"
#include <math.h>
#include <stdarg.h>
#include <string.h>
#ifdef HAVE_MPI
#include <mpi.h>
#endif
#include <latan/latan_math.h>
#include <latan/latan_plot.h>
#include "data_loader.h"
#include "models.h"

#define MOD_PLOT_NPT 1000

#define PLOT_ADD_FIT(obj,kx,ky,title,color)\
plot_add_fit(p[kx],d,ky,phy_pt,kx,fit,x_range[kx][0],x_range[kx][1],\
             MOD_PLOT_NPT,true,obj,title,"",color,color);\

#define PLOT_ADD_PB(kx,ky,color)\
plot_add_fit_predband(p[kx],d,ky,phy_pt,kx,s_fit,x_range[kx][0],x_range[kx][1],\
                      MOD_PLOT_NPT/4,color);

#define PLOT_ADD_EX(kx,s)\
if ((!latan_isnan(param->q_target[0]))&&(!latan_isnan(param->q_target[1])))\
{\
    plot_add_datpoint(p[kx],mat_get(phy_pt,kx,0),param->q_target[0],      \
                      -1.0,param->q_target[1],"target","rgb 'dark-blue'");\
}\
plot_add_datpoint(p[kx],mat_get(phy_pt,kx,0),\
                  mat_get(rs_sample_pt_cent_val(param->s_ex),s,0), \
                  -1.0,mat_get(param->ex_err,s,0),"physical point",\
                  "rgb 'black'");

#define PLOT_DISP(kx,name)\
plot_set_title(p[kx],gtitle);\
plot_set_xlabel(p[kx],xlabel);\
plot_set_ylabel(p[kx],ylabel);\
plot_disp(p[kx]);\
if (strlen(param->save_plot))\
{\
    strbuf dirname;\
    sprintf(dirname,"%s_%s",param->save_plot,name);\
    plot_save(dirname,p[kx]);\
}

/* MPI friendly printf function */
void mpi_printf(const strbuf fmt, ...)
{
    va_list args;
    int proc,nproc;
    strbuf buf;
    
#ifdef HAVE_MPI
    MPI_Comm_size(MPI_COMM_WORLD,&nproc);
    if (nproc > 1)
    {
        MPI_Comm_rank(MPI_COMM_WORLD,&proc);
        sprintf(buf,"[%d] %s",proc,fmt);
    }
    else
    {
        strbufcpy(buf,fmt);
    }
#else
    strbufcpy(buf,fmt);
#endif
    va_start(args,fmt);
    vprintf(buf,args);
    va_end(args);
}

void plot_fit(const rs_sample *s_fit, fit_data *d, fit_param *param, const plot_flag f)
{
    plot *p[N_EX_VAR];
    double *xb[N_EX_VAR] = {NULL,NULL,NULL,NULL,NULL,NULL,NULL};
    double x_range[N_EX_VAR][2],b_int[2],dbind,a;
    size_t bind,vind,eind,k,phy_ind,s;
    strbuf color,gtitle,title,xlabel,ylabel;
    mat *phy_pt,*x_k,*fit,*cordat,**vol_av_corr,*yerrtmp;
    ens *ept;
    
    phy_pt     = mat_create(N_EX_VAR,1);
    x_k        = mat_create(param->nens,1);
    cordat     = mat_create(param->nens,1);
    MALLOC(vol_av_corr,mat **,param->nbeta);
    for (bind=0;bind<param->nbeta;bind++)
    {
        vol_av_corr[bind] = mat_create(param->nvol[bind],1);
    }

    for (k=0;k<N_EX_VAR;k++)
    {
        p[k] = plot_create();
    }
    
    param->scale_model = 1;
    fit = rs_sample_pt_cent_val(s_fit);
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
        if (k == i_vind)
        {
            fit_data_get_x_k(x_k,d,i_Linv);
        }
        else
        {
            fit_data_get_x_k(x_k,d,k);
        }
        if ((k == i_a)||(k == i_ud)||(k == i_Linv)||(k == i_alpha)\
            ||(k == i_vind))
        {
            x_range[k][0] = 0.0;
        }
        else
        {
            x_range[k][0] = mat_get_min(x_k)-0.15*fabs(mat_get_min(x_k));
        }
        x_range[k][1] = mat_get_max(x_k)+0.15*fabs(mat_get_min(x_k));
        plot_set_scale_xmanual(p[k],x_range[k][0],x_range[k][1]);
    }
    if (f == Q)
    {
        sprintf(gtitle,"quantity: %s -- scale: %s -- datasets: %s -- ensembles: %s",\
                param->q_name,param->scale_part,param->dataset_cat,param->manifest);
        mat_set(phy_pt,i_ud,0,SQ(param->M_ud));
        mat_set(phy_pt,i_s,0,SQ(param->M_s));
        mat_set(phy_pt,i_umd,0,param->M_umd_val);
        mat_set(phy_pt,i_alpha,0,param->alpha);
        mat_set(phy_pt,i_bind,0,0.0);
        mat_set(phy_pt,i_vind,0,0.0);
        mat_set(phy_pt,i_a,0,0.0);
        mat_set(phy_pt,i_Linv,0,0.0);
        mat_set(phy_pt,i_fvM,0,param->qed_fvol_mass);
        /* regular plots */
        PLOT_ADD_FIT(PF_FIT,i_ud,phy_ind,"","rgb 'black'");
        PLOT_ADD_PB(i_ud,phy_ind,"rgb 'black'");
        PLOT_ADD_FIT(PF_FIT,i_s,phy_ind,"","rgb 'black'");
        PLOT_ADD_PB(i_s,phy_ind,"rgb 'black'");
        PLOT_ADD_FIT(PF_FIT,i_umd,phy_ind,"","rgb 'black'");
        PLOT_ADD_PB(i_umd,phy_ind,"rgb 'black'");
        PLOT_ADD_FIT(PF_FIT,i_alpha,phy_ind,"","rgb 'black'");
        PLOT_ADD_PB(i_alpha,phy_ind,"rgb 'black'");
        PLOT_ADD_FIT(PF_FIT,i_a,phy_ind,"","rgb 'black'");
        PLOT_ADD_PB(i_a,phy_ind,"rgb 'black'");
        PLOT_ADD_FIT(PF_FIT,i_Linv,phy_ind,"","rgb 'black'");
        PLOT_ADD_PB(i_Linv,phy_ind,"rgb 'black'");
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
            PLOT_ADD_FIT(PF_DATA,i_alpha,phy_ind,title,color);
            PLOT_ADD_FIT(PF_DATA,i_a,phy_ind,title,color);
            PLOT_ADD_FIT(PF_DATA,i_Linv,phy_ind,title,color);
            fit_data_fit_all_points(d,true);
        }
        /* volume averages plot */
        plot_add_fit(p[i_vind],d,phy_ind,phy_pt,i_Linv,fit,x_range[i_Linv][0],\
                     x_range[i_Linv][1],MOD_PLOT_NPT,true,PF_FIT,"","",       \
                     "rgb 'black'","rgb 'black'");
        plot_add_fit_predband(p[i_vind],d,phy_ind,phy_pt,i_Linv,s_fit,\
                              x_range[i_Linv][0],x_range[i_Linv][1],  \
                              MOD_PLOT_NPT/4,"rgb 'black'");
        fit_partresidual(cordat,d,phy_ind,phy_pt,i_Linv,fit);
        for(bind=0;bind<param->nbeta;bind++)
        {
            mat_zero(vol_av_corr[bind]);
        }
        for(eind=0;eind<param->nens;eind++)
        {
            ept  = param->point + eind;
            bind = (size_t)ind_beta(ept->beta,param);
            vind = (size_t)ind_volume((unsigned int)ept->L,(int)bind,param);
            mat_inc(vol_av_corr[bind],vind,0,mat_get(cordat,eind,0));
        }
        for (bind=0;bind<param->nbeta;bind++)
        for (vind=0;vind<param->nvol[bind];vind++)
        {
            mat_set(vol_av_corr[bind],vind,0,               \
                    mat_get(vol_av_corr[bind],vind,0)       \
                    /((double)(param->nenspvol[bind][vind])));
        }
        for(bind=0;bind<param->nbeta;bind++)
        {
            yerrtmp = mat_create(param->nvol[bind],1);
            
            rs_sample_varp(yerrtmp,param->s_vol_av[bind]);
            mat_eqsqrt(yerrtmp);
            sprintf(color,"%d",1+(int)bind);
            sprintf(title,"beta = %s",param->beta[bind]);
            plot_add_dat_yerr(p[i_vind],                                     \
                              rs_sample_pt_cent_val(param->s_vol_Linv[bind]),\
                              vol_av_corr[bind],yerrtmp,title,color);
            
            mat_destroy(yerrtmp);
        }

        /* display plots */
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
        PLOT_DISP(i_ud,"ud");
        sprintf(xlabel,"M_%s^2 (MeV^2)",param->s_name);
        PLOT_ADD_EX(i_s,s);
        PLOT_DISP(i_s,"s");
        strbufcpy(xlabel,"a (MeV^-1)");
        PLOT_ADD_EX(i_a,s);
        PLOT_DISP(i_a,"a");
        if (param->have_umd)
        {
            sprintf(xlabel,"%s (MeV^2)",param->umd_name);
            PLOT_ADD_EX(i_umd,s);
            PLOT_DISP(i_umd,"umd");
        }
        if (param->have_alpha)
        {
            strbufcpy(xlabel,"alpha");
            PLOT_ADD_EX(i_alpha,s);
            PLOT_DISP(i_alpha,"alpha");
        }
        strbufcpy(xlabel,"1/L (MeV)");
        PLOT_ADD_EX(i_Linv,s);
        PLOT_DISP(i_Linv,"Linv");
        PLOT_ADD_EX(i_vind,s);
        PLOT_DISP(i_vind,"Linv_av");
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
            mat_set(phy_pt,i_umd,0,SQ(a)*param->M_umd_val);
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
        PLOT_DISP(i_ud,"ud");
        sprintf(xlabel,"(a*M_%s)^2",param->s_name);
        PLOT_DISP(i_s,"s");
        if (param->have_umd)
        {
            sprintf(xlabel,"a^2*%s",param->umd_name);
            PLOT_DISP(i_umd,"umd");
        }
        strbufcpy(xlabel,"a/L");
        PLOT_DISP(i_Linv,"Linv");
    }
    param->scale_model = 0;
    
    mat_destroy(phy_pt);
    mat_destroy(x_k);
    mat_destroy(cordat);
    for (bind=0;bind<param->nbeta;bind++)
    {
        mat_destroy(vol_av_corr[bind]);
    }
    free(vol_av_corr);
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
        bind = (size_t)ind_beta(param->point[i].beta,param);
        fprintf(tmpf[bind],"%s %e %e\n",param->point[i].dir,(double)(i),\
                mat_get(d->chi2_comp,i+k*nens,0));
    }
    for (i=0;i<nbeta;i++)
    {
        fclose(tmpf[i]);
        sprintf(plotcmd,"u 2:3:xtic(1) t '%s' w impulse",param->beta[i]);
        plot_add_plot(p,plotcmd,tmpfname[i]);
    }
    plot_add_head(p,"set xtics rotate by -90 font 'courier, 10'");
    plot_set_scale_manual(p,-1.0,(double)(param->nens+1),-5.0,5.0);
    plot_add_plot(p,"0.0 lt -1 lc rgb 'black' notitle","");
    plot_add_plot(p,"1.0 lt -1 lc rgb 'black' notitle","");
    plot_add_plot(p,"-1.0 lt -1 lc rgb 'black' notitle","");
    plot_add_plot(p,"2.0 lt -1 lc rgb 'dark-gray' notitle","");
    plot_add_plot(p,"-2.0 lt -1 lc rgb 'dark-gray' notitle","");
    plot_add_plot(p,"3.0 lt -1 lc rgb 'gray' notitle","");
    plot_add_plot(p,"-3.0 lt -1 lc rgb 'gray' notitle","");
    plot_add_plot(p,"4.0 lt -1 lc rgb 'light-gray' notitle","");
    plot_add_plot(p,"-4.0 lt -1 lc rgb 'light-gray' notitle","");
    plot_set_ylabel(p,"standard deviations");
    plot_set_title(p,title);
    plot_disp(p);
    
    free(tmpf);
    free(tmpfname);
    plot_destroy(p);
}

#define PRINT_PAR(name)\
{\
    mpi_printf("%16s = % e ( %4.0f%% )\n",name,mat_get(fit,i,0),          \
               sqrt(mat_get(fit_var,i,0))/fabs(mat_get(fit,i,0))*100.0);\
    i++;\
}
#define PRINT_SCALE(name)\
{\
    double sig_,p_;\
    sig_ = sqrt(mat_get(fit_var,i,0));\
    p_   = mat_get(fit,i,0);\
    mpi_printf("%16s = % e ( %4.1f%% ) [%12s^-1 = %5.0f(%4.0f) MeV]\n",name,\
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
    mpi_printf("\n");
    mpi_printf("fit parameters :\n");
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
        if (param->s_with_a2ud)
        {
            PRINT_PAR("s_a2ud");
        }
        if (param->s_with_a2s)
        {
            PRINT_PAR("s_a2s");
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
        if (param->with_const)
        {
            PRINT_PAR("const");
        }
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
        if (param->with_a2 > 0)
        {
            PRINT_PAR("p_a2");
        }
        if (param->with_alpha_sa > 0)
        {
            PRINT_PAR("p_alpha_sa");
        }
        if (param->with_a2ud)
        {
            PRINT_PAR("p_a2ud");
        }
        if (param->with_a2s)
        {
            PRINT_PAR("p_a2s");
        }
        if (param->umd_deg)
        {
            PRINT_PAR("p_umd");
        }
        if (param->with_udumd)
        {
            for (j=0;j<param->with_udumd;j++)
            {
                sprintf(buf,"p_udumd_%d",j+1);
                PRINT_PAR(buf);
            }
        }
        if (param->with_sumd)
        {
            for (j=0;j<param->with_sumd;j++)
            {
                sprintf(buf,"p_sumd_%d",j+1);
                PRINT_PAR(buf);
            }
        }
        if (param->with_a2umd)
        {
            PRINT_PAR("p_a2umd");
        }
        if (param->with_alpha_saumd)
        {
            PRINT_PAR("p_alpha_saumd");
        }
        if (param->alpha_deg)
        {
            PRINT_PAR("p_alpha");
        }
        if (param->with_udalpha)
        {
            for (j=0;j<param->with_udalpha;j++)
            {
                sprintf(buf,"p_udalpha_%d",j+1);
                PRINT_PAR(buf);
            }
        }
        if (param->with_salpha)
        {
            for (j=0;j<param->with_salpha;j++)
            {
                sprintf(buf,"p_salpha_%d",j+1);
                PRINT_PAR(buf);
            }
        }
        if (param->with_aalpha)
        {
            PRINT_PAR("p_aalpha");
        }
        if (param->with_qed_fvol)
        {
            if (param->with_qed_fvol_monopmod)
            {
                if (param->with_qed_fvol == 2)
                {
                    PRINT_PAR("p_Linv_2");
                }
            }
            else
            {
                for (j=0;j<param->with_qed_fvol;j++)
                {
                    sprintf(buf,"p_Linv_%d",j+1);
                    PRINT_PAR(buf);
                }
            }
        }
        if (IS_AN(param,AN_PHYPT)&&!IS_AN(param,AN_SCALE)&&\
            (strbufcmp(param->with_ext_a,"") != 0))
        {
            for(j=0;j<(int)param->nbeta;j++)
            {
                sprintf(buf,"corr_a_%s",param->beta[j]);
                PRINT_SCALE(buf);
            }
        }
    }
    mpi_printf("\n");
    
    mat_destroy(fit_var);
}
#undef PRINT_PAR

#define PRINT_DLABEL(name) fprintf(stream,"%-12s  ",name)
#define PRINT_DLABEL_WERR(name) fprintf(stream,"%-12s %-12s  ",name,"error")
#define PRINT_D(value) fprintf(stream,"% .5e  ",value)
#define PRINT_D_WERR(value,err) fprintf(stream,"% .5e % .5e  ",value,err)
#define PRINT_X(ind)\
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
    fprintf(stream,"#%49s  ","ensemble");
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
    if (param->have_alpha)
    {
        PRINT_DLABEL("alpha");
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
        fprintf(stream,"%50s  ",param->point[ens_ind].dir);
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
        if (param->have_alpha)
        {
            PRINT_X(i_alpha);
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
