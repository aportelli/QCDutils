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
#define STAT_ERR_COL "rgb 'red'"
#define SYS_ERR_COL  "rgb 'blue'"
#define TOT_ERR_COL  "rgb '#444444'"
#define TARG_ERR_COL "rgb 'green'"

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

int main(int argc, char *argv[])
{
    rs_sample *s_res,*s_med;
    mat *chi2_val,*w,*hist,*hist_flat,*phist,*med_var,*res,*med;
    size_t nsample,nres,nbin;
    size_t j,s,count;
    double target,target_err,xmin,xmax,l,ymax,stat_err,sys_err[2],tot_err[2],\
           med_s,final,cl[2];
    plot *p;
    int i;
    
    if (argc <= 4)
    {
        fprintf(stderr,"usage: %s <target> <target_err> <result_1> <result_2> ...\n",\
                argv[0]);
        return EXIT_FAILURE;
    }
    
    rs_sample_load(NULL,&nsample,NULL,argv[3]);
    nres       = (size_t)(argc) - 3;
    nbin       = 10;
    count      = 0;
    target     = ATOF(argv[1]);
    target_err = ATOF(argv[2]);
    
    
    s_res     = rs_sample_create(nres,1,nsample);
    s_med     = rs_sample_create(1,1,nsample);
    chi2_val  = mat_create(nres,2);
    w         = mat_create(nres,1);
    hist      = mat_create(nbin,1);
    hist_flat = mat_create(nbin,1);
    phist     = mat_create(40,1);
    med_var   = mat_create(1,1);
    
    res = rs_sample_pt_cent_val(s_res);
    med = rs_sample_pt_cent_val(s_med);
    
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
            
            load_res(s_res_i,chi2_val_i,argv[i+3]);
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
    
    /* compute median */
    printf("-- resampling median...\n");
    med_s = mat_elpercentile(res,w,50.0);
    mat_set(med,0,0,med_s);
    for (s=0;s<nsample;s++)
    {
        med_s = mat_elpercentile(rs_sample_pt_sample(s_res,s),w,50.0);
        mat_set(rs_sample_pt_sample(s_med,s),0,0,med_s);
    }
    
    /* display result */
    rs_sample_varp(med_var,s_med);
    conf_int(sys_err,res,w,1.0);
    final      = mat_get(med,0,0);
    stat_err   = sqrt(mat_get(med_var,0,0));
    sys_err[0] = final - sys_err[0];
    sys_err[1] = sys_err[1] - final;
    tot_err[0] = sqrt(SQ(stat_err)+SQ(sys_err[0]));
    tot_err[1] = sqrt(SQ(stat_err)+SQ(sys_err[1]));
    printf("result: %f (%f)[stat.](-%f,+%f)[sys.](-%f,+%f)[tot.]\n",final,\
           stat_err,sys_err[0],sys_err[1],tot_err[0],tot_err[1]);
    
    /* compute histograms */
    xmin = MIN(target-target_err,final-tot_err[0]);
    xmax = MAX(target+target_err,final+tot_err[1]);
    l    = xmax-xmin;
    xmin = xmin - 0.2*l;
    xmax = xmax + 0.2*l;
    histogram(hist,rs_sample_pt_cent_val(s_res),w,xmin,xmax,nbin);
    histogram(hist_flat,rs_sample_pt_cent_val(s_res),NULL,xmin,xmax,nbin);
    histogram(phist,w,NULL,0.0,1.0,40);
    
    /* make plot */
    /** result histogram **/
    p    = plot_create();
    ymax = MAX(mat_get_max(hist)*DRATIO(nbin,mat_elsum(w)*(xmax-xmin)),\
               mat_get_max(hist_flat)*DRATIO(nbin,nres*(xmax-xmin)));
    ymax = 1.2*ymax;
    plot_set_scale_manual(p,xmin-0.3*l,xmax+0.3*l,0.0,ymax);
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
    plot_add_histogram(p,hist,xmin,xmax,mat_elsum(w),true,"","");
    plot_disp(p);
    plot_save("qcd_error",p);
    plot_destroy(p);
    /** p-value histogram **/
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
    
    io_finish();
    
    rs_sample_destroy(s_res);
    rs_sample_destroy(s_med);
    mat_destroy(chi2_val);
    mat_destroy(w);
    mat_destroy(hist);
    mat_destroy(phist);
}
