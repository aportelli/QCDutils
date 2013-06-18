#define _POSIX_C_SOURCE 199506L /* strtok_r is used here */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <latan/latan_io.h>
#include <latan/latan_math.h>
#include <latan/latan_plot.h>
#include <latan/latan_statistics.h>

#define ATOF(str) (strtod(str,(char **)NULL))
#define ATOI(str) ((int)strtol(str,(char **)NULL,10))
#define STAT_ERR_COL "rgb 'red'"
#define SYS_ERR_COL  "rgb 'blue'"
#define TOT_ERR_COL  "rgb '#444444'"
#define TARG_ERR_COL "rgb 'green'"
#define NPBIN 40

static void load_res(rs_sample *s_res, double chi2_val[2],\
                     const strbuf latan_path)
{
    char *pt;
    strbuf buf,*field;
    size_t nf,lc;
    
    field = NULL;

    strbufcpy(buf,latan_path);
    rs_sample_load_subsamp(s_res,latan_path,0,0,0,0);
    pt = strstr(buf,".boot");
    strncpy(pt,".chi2\0",6);
    BEGIN_FOR_LINE_TOK(field,buf," \t",nf,lc)
    {
        if ((nf >= 4)&&(strbufcmp(field[0],"correlated:") == 0)) 
        {
            chi2_val[0] = ATOF(field[2]);
            chi2_val[1] = ATOF(field[3]);
        }
    }
    END_FOR_LINE_TOK(field)
}

static void plot_res_hist(const mat *res, const mat *w, const double final,\
                          const double stat_err, const double sys_err[2],  \
                          const double target, const double target_err,    \
                          const size_t nbin)
{
    double xmin,xmax,ymax,l,tot_err[2];
    mat *hist,*hist_flat;
    plot *p;
    size_t nres;
    
    hist      = mat_create(nbin,1);
    hist_flat = mat_create(nbin,1);
    
    nres       = nrow(res);
    tot_err[0] = sqrt(SQ(stat_err)+SQ(sys_err[0]));
    tot_err[1] = sqrt(SQ(stat_err)+SQ(sys_err[1]));
    
    /* compute histograms */
    xmin  = MIN(target-target_err,final-tot_err[0]);
    xmin  = MIN(xmin,mat_get_min(res));
    xmax  = MAX(target+target_err,final+tot_err[1]);
    xmax  = MAX(xmax,mat_get_max(res));
    l     = xmax-xmin;
    xmin -= 0.25*l;
    xmax += 0.25*l;
    histogram(hist,res,w,xmin,xmax,nbin);
    histogram(hist_flat,res,NULL,xmin,xmax,nbin);
    
    /* make plot */
    p    = plot_create();
    ymax = MAX(mat_get_max(hist)*DRATIO(nbin,mat_elsum(w,NULL)*(xmax-xmin)),\
               mat_get_max(hist_flat)*DRATIO(nbin,nres*(xmax-xmin)));
    ymax = 1.2*ymax;
    plot_set_scale_manual(p,xmin-0.2*l,xmax+0.2*l,0.0,ymax);
    plot_add_vlineaerr(p,final,tot_err,TOT_ERR_COL);
    if ((stat_err > sys_err[0])&&(stat_err > sys_err[1]))
    {
        plot_add_vlineerr(p,final,stat_err,STAT_ERR_COL);
        plot_add_vlineaerr(p,final,sys_err,SYS_ERR_COL);
    }
    else
    {
        plot_add_vlineaerr(p,final,sys_err,SYS_ERR_COL);
        plot_add_vlineerr(p,final,stat_err,STAT_ERR_COL);
    }
    plot_add_vlineerr(p,target,target_err,TARG_ERR_COL);
    plot_add_vline(p,final-tot_err[0],TOT_ERR_COL);
    plot_add_vline(p,final+tot_err[1],TOT_ERR_COL);
    plot_add_vline(p,final-stat_err,STAT_ERR_COL);
    plot_add_vline(p,final+stat_err,STAT_ERR_COL);
    plot_add_vline(p,final-sys_err[0],SYS_ERR_COL);
    plot_add_vline(p,final+sys_err[1],SYS_ERR_COL);
    plot_add_vline(p,final,"rgb 'black'");
    plot_add_histogram(p,hist_flat,xmin,xmax,(double)(nres),true,"",\
                       "rgb 'dark-gray'");
    plot_add_histogram(p,hist,xmin,xmax,mat_elsum(w,NULL),true,"","");
    plot_disp(p);
    plot_save("qcd_error",p);
    plot_destroy(p);
    
    mat_destroy(hist);
    mat_destroy(hist_flat);
}

static void plot_p_hist(const mat *w)
{
    double ymax,cl[2];
    mat *phist;
    plot *p;
    size_t nres;
    
    phist = mat_create(NPBIN,1);
    
    nres = nrow(w);
    
    histogram(phist,w,NULL,0.0,1.0,NPBIN);
    p     = plot_create();
    ymax  = mat_get_max(phist);
    ymax  = 1.2*ymax;
    cl[1] = 0.0;
    plot_set_scale_manual(p,0.0,1.0,0.0,ymax);
    cl[0] = 0.954499736;
    plot_add_vlineaerr(p,1.0,cl,"rgb 'dark-gray'");
    cl[0] = 0.682689492;
    plot_add_vlineaerr(p,1.0,cl,"rgb 'black'");
    plot_add_histogram(p,phist,0.0,1.0,(double)(nres),false,"","rgb 'red'");
    plot_disp(p);
    plot_destroy(p);
    
    mat_destroy(phist);
}

int main(int argc, char *argv[])
{
    rs_sample *s_res,*s_med,*s_mean;
    mat *chi2_val,*w,*med_var,*mean_var,*res,*med,*mean;
    size_t nsample,nres,nbin;
    size_t j,s,count;
    double target,target_err,stat_err,sys_err[2],tot_err[2],med_s,mean_s,final,\
           res_var;
    int i;
    
    if (argc <= 5)
    {
        fprintf(stderr,"usage: %s <nbin> <target> <target_err> <result_1> <result_2> ...\n",\
                argv[0]);
        return EXIT_FAILURE;
    }
    
    rs_sample_load(NULL,&nsample,NULL,argv[4]);
    nres       = (size_t)(argc) - 4;
    nbin       = ATOF(argv[1]);
    count      = 0;
    target     = ATOF(argv[2]);
    target_err = ATOF(argv[3]);
    
    
    s_res     = rs_sample_create(nres,1,nsample);
    s_med     = rs_sample_create(1,1,nsample);
    s_mean    = rs_sample_create(1,1,nsample);
    chi2_val  = mat_create(nres,2);
    w         = mat_create(nres,1);
    med_var   = mat_create(1,1);
    mean_var  = mat_create(1,1);
    
    res  = rs_sample_pt_cent_val(s_res);
    med  = rs_sample_pt_cent_val(s_med);
    mean = rs_sample_pt_cent_val(s_mean);
    
    io_init();
    
    /* load results */
    printf("-- loading results...\n");
#ifdef _OPENMP
    #pragma omp parallel
#endif
    {
        rs_sample *s_res_i;
        
        s_res_i = rs_sample_create(1,1,nsample);
#ifdef _OPENMP
        #pragma omp for
#endif
        for (i=0;i<(int)(nres);i++) 
        {
            double chi2_val_i[2];
            
            load_res(s_res_i,chi2_val_i,argv[i+4]);
            rs_sample_set_subsamp(s_res,s_res_i,i,0,i,0);
            mat_set(chi2_val,i,0,chi2_val_i[0]);
            mat_set(chi2_val,i,1,chi2_val_i[1]);
            mat_set(w,i,0,chi2_pvalue(chi2_val_i[0],(size_t)chi2_val_i[1]));
#ifdef _OPENMP
            #pragma omp critical
#endif
            {
                count++;
                printf("[");
                for (j=0;j<60*count/nres;j++)
                {
                    printf("=");
                }
                for (j=60*count/nres;j<60;j++)
                {
                    printf(" ");
                }
                printf("]  %d/%d\r",(int)count,(int)nres);
                fflush(stdout);
            }
        }
        rs_sample_destroy(s_res_i);
    }
    printf("\n%d sample files successfully loaded\n",(int)nres);
    
    /* compute median and mean */
    printf("-- resampling median and mean...\n");
    med_s  = mat_elpercentile(res,w,50.0);
    mean_s = mat_elmean(res,w);
    mat_set(med,0,0,med_s);
    mat_set(mean,0,0,mean_s);
    for (s=0;s<nsample;s++)
    {
        med_s  = mat_elpercentile(rs_sample_pt_sample(s_res,s),w,50.0);
        mean_s = mat_elmean(rs_sample_pt_sample(s_res,s),w);
        mat_set(rs_sample_pt_sample(s_med,s),0,0,med_s);
        mat_set(rs_sample_pt_sample(s_mean,s),0,0,mean_s);
    }
    
    /* display mean result */
    rs_sample_varp(mean_var,s_mean);
    res_var    = mat_elvar(res,w);
    final      = mat_get(mean,0,0);
    stat_err   = sqrt(mat_get(mean_var,0,0));
    sys_err[0] = sqrt(res_var);
    sys_err[1] = sqrt(res_var);
    tot_err[0] = sqrt(SQ(stat_err)+SQ(sys_err[0]));
    tot_err[1] = sqrt(SQ(stat_err)+SQ(sys_err[1]));
    printf("result (mean): %f (%f)[stat.](-%f,+%f)[sys.](-%f,+%f)[tot.]\n",final,\
           stat_err,sys_err[0],sys_err[1],tot_err[0],tot_err[1]);
    plot_res_hist(res,w,final,stat_err,sys_err,target,target_err,nbin);
    
    /* display median result */
    rs_sample_varp(med_var,s_med);
    conf_int(sys_err,res,w,1.0);
    final      = mat_get(med,0,0);
    stat_err   = sqrt(mat_get(med_var,0,0));
    sys_err[0] = final - sys_err[0];
    sys_err[1] = sys_err[1] - final;
    tot_err[0] = sqrt(SQ(stat_err)+SQ(sys_err[0]));
    tot_err[1] = sqrt(SQ(stat_err)+SQ(sys_err[1]));
    printf("result (med) : %f (%f)[stat.](-%f,+%f)[sys.](-%f,+%f)[tot.]\n",final,\
           stat_err,sys_err[0],sys_err[1],tot_err[0],tot_err[1]);
    plot_res_hist(res,w,final,stat_err,sys_err,target,target_err,nbin);
    
    io_finish();
    
    rs_sample_destroy(s_res);
    rs_sample_destroy(s_med);
    mat_destroy(chi2_val);
    mat_destroy(w);
    mat_destroy(med_var);
    mat_destroy(mean_var);
}
