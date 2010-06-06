#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <qcd_arg_parse.h>
#include <latan/latan_hadron.h>
#include <latan/latan_io.h>
#include <latan/latan_mat.h>
#include <latan/latan_plot.h>
#include <latan/latan_rand.h>
#include <latan/latan_statistics.h>

#define NBOOT 2000

int main(int argc, char* argv[])
{
	/*				parsing arguments			*/
	/********************************************/
	ss_no source,sink;
	channel_no ch;
	qcd_options opt;
	stringbuf spec_name,part_name,manf_name;
	size_t binsize;
	 
	opt = qcd_arg_parse(argc,argv,A_PARTICLE|A_PROP_LOAD|A_SAVE_RS\
						|A_CHANNEL|A_LOAD_RG);
	strcpy(spec_name,opt->spec_name);
	strcpy(part_name,opt->part_name);
	strcpy(manf_name,opt->manf_name);
	source = opt->source;
	sink = opt->sink;
	binsize = opt->binsize;
	for (ch=0;ch<NCHANNEL;ch++)
	{
		channel_id_set(ch,opt->channel_id[ch]);
	}
	latan_set_verb(opt->latan_verb);
	
	/*			identifying particle			*/
	/********************************************/
	spectrum s;
	hadron h;
	
	if (strcmp(spec_name,"qcd") == 0)
	{
		s = spectrum_create_qcd();
	}
	else if (strcmp(spec_name,"qcdqed") == 0)
	{
		s = spectrum_create_qcdqed();
	}
	else
	{
		s = NULL;
		fprintf(stderr,"error: spectrum %s unknown\n",spec_name);
		return EXIT_FAILURE;
	}
	h = spectrum_get(s,part_name);
	
	/*				loading datas				*/
	/********************************************/
	size_t ndat,nt;
	mat* prop;
	
	ndat	= (size_t)get_nfile(manf_name);
	nt		= (size_t)hadron_getnt(h,source,sink,manf_name);
	
	prop = mat_ar_create(ndat,nt,1);
	
	printf("-- loading %s datas from %s...\n",h->name,manf_name);
	hadron_propbin(prop,h,source,sink,manf_name,binsize);
	
	/*		resampling mean propagator			*/
	/********************************************/
	rs_sample s_mprop;
	stringbuf sample_name;
	mat mprop;
	
	sprintf(sample_name,"%s_prop_%s",h->name,manf_name);
	s_mprop = rs_sample_create_boot(nt,NBOOT,sample_name);
	
	printf("-- resampling %s mean propagator...\n",h->name);
	randgen_set_state(opt->state);
	resample(s_mprop,prop,ndat,1,&rs_mean,NULL);
	mprop = rs_sample_pt_cent_val(s_mprop);
	
	/*	computing error on mean propagator		*/
	/********************************************/
	mat sig;
	
	sig = mat_create(nt,1);
	
	printf("-- estimating %s mean propagator error...\n",h->name);
	rs_sample_varp(sig,s_mprop);
	mat_eqsqrt(sig);
	
	/*				result output				*/
	/********************************************/
	size_t t,maxt;
	
	maxt = (h->parity == EVEN) ? (nt/2) : (nt-1);
	printf("\nt\tprop\t\terror\n");
	for (t=0;t<=maxt;t++)
	{
		printf("%d\t%e\t%e\n",(int)t,mat_get(mprop,t,0),mat_get(sig,t,0));
	}
	printf("\n");
	if (opt->do_save_rs_sample)
	{
		rs_sample_save(s_mprop,s_mprop->name);
	}
	
	
	/*					plot					*/
	/********************************************/
	/*
	if (opt->do_plot)
	{
		plot p;
		stringbuf key;
		const double dmaxt = (double)maxt;
		
		p = plot_create();
		
		sprintf(key,"%s propagator",h->name);
		plot_set_scale_ylog(p);
		plot_set_scale_xmanual(p,0.0,dmaxt);
		plot_add_daterr(p,mprop,sig,0.0,1.0,key);
		plot_disp(p);	
		
		plot_destroy(p);
	}
	*/
	
	/*				desallocation				*/
	/********************************************/
	FREE(opt);
	spectrum_destroy(s);
	mat_ar_destroy(prop,ndat);
	rs_sample_destroy(s_mprop);
	mat_destroy(sig);
	
	return EXIT_SUCCESS;
}