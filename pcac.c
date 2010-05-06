#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <qcd_arg_parse.h>
#include <latan/latan_hadron.h>
#include <latan/latan_io.h>
#include <latan/latan_mat.h>
#include <latan/latan_nunits.h>
#include <latan/latan_plot.h>
#include <latan/latan_rand.h>
#include <latan/latan_statistics.h>

#define NBOOT 2000

int main(int argc, char* argv[])
{
	/*				parsing arguments			*/
	/********************************************/
	int source,sink;
	double latspac_nu;
	quark_no q1,q2;
	qcd_options opt;
	stringbuf manf_name,unit;
	
	opt = qcd_arg_parse(argc,argv,A_QCOMP|A_LATSPAC);
	strcpy(manf_name,opt->manf_name);
	source = opt->source;
	sink = opt->sink;
	q1 = quark_no_get(opt->qcomp[0]);
	q2 = quark_no_get(opt->qcomp[1]);
	latspac_nu = opt->latspac_nu;
	if (opt->have_latspac)
	{
		strcpy(unit," (MeV)");
	}
	else
	{
		strcpy(unit,"");
	}
	
	/*			creating PCAC "particle"		*/
	/********************************************/
	hadron h_AP,h_PP;
	
	h_AP = hadron_create();
	h_PP = hadron_create();
	
	hadron_set_2q_nomix(h_AP,"AP",ODD,ch_AP,q1,q2);
	hadron_set_2q_nomix(h_PP,"PP",ODD,ch_PP,q1,q2);
	
	/*				loading datas				*/
	/********************************************/
	size_t ndat,nt;
	mat* prop;
	mat* prop_AP;
	mat* prop_PP;
	
	ndat	= (size_t)get_nfile(manf_name);
	nt		= (size_t)hadron_getnt(h_AP,source,sink,manf_name);
	
	prop = mat_ar_create(2*ndat,nt,1);
	prop_AP = prop;
	prop_PP = prop + ndat;
	
	printf("-- loading %s datas from %s...\n",h_AP->name,manf_name);
	hadron_prop(prop_AP,h_AP,source,sink,manf_name);
	printf("-- loading %s datas from %s...\n",h_PP->name,manf_name);
	hadron_prop(prop_PP,h_PP,source,sink,manf_name);
	
	
	/*		resampling mean propagator			*/
	/********************************************/
	rs_sample s_effmass_pcac;
	stringbuf sample_name;
	mat effmass_pcac;
	
	sprintf(sample_name,"effmass_PCAC_%s",manf_name);
	s_effmass_pcac = rs_sample_create_boot(nt-2,NBOOT,sample_name);
	
	printf("-- resampling PCAC effective mass...\n");
	randgen_init_from_time();
	resample(s_effmass_pcac,prop,ndat,2,&rs_effmass_PCAC,NULL);
	effmass_pcac = rs_sample_get_cent_val(s_effmass_pcac);
	
	/*	computing error on mean propagator		*/
	/********************************************/
	mat sig;
	
	sig = mat_create(nt-2,1);
	
	printf("-- computing PCAC effective mass error...\n");
	rs_sample_varp(sig,s_effmass_pcac);
	mat_eqsqrt(sig);
	
	/*				result output				*/
	/********************************************/
	size_t t;
	
	mat_eqmuls(effmass_pcac,1.0/latspac_nu);
	mat_eqmuls(sig,1.0/latspac_nu);
	printf("\nt\tm_PCAC%s\t\terror%s\n",unit,unit);
	for (t=0;t<nt-2;t++)
	{
		printf("%d\t%e\t\t%e\n",(int)t+1,mat_get(effmass_pcac,t,0),\
			   mat_get(sig,t,0));
	}
	printf("\n");
	if (opt->do_save_rs_sample)
	{
		rs_sample_save(s_effmass_pcac,s_effmass_pcac->name);
	}
	
	/*					plot					*/
	/********************************************/
	if (opt->do_plot)
	{
		plot p;
		stringbuf key;
		
		p = plot_create();
		
		sprintf(key,"PCAC %s effective mass",opt->qcomp);
		plot_set_scale_xmanual(p,0.0,nt/2);
		plot_add_daterr(p,effmass_pcac,sig,1.0,1.0,key);
		plot_disp(p);	
		
		plot_destroy(p);
	}
	
	/*				desallocation				*/
	/********************************************/
	FREE(opt);
	hadron_destroy(h_AP);
	mat_ar_destroy(prop,2*ndat);
	rs_sample_destroy(s_effmass_pcac);
	
	return EXIT_SUCCESS;
}