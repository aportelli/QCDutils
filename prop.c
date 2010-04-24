#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <latan/hadron.h>
#include <latan/io.h>
#include <latan/mat.h>
#include <latan/plot.h>
#include <latan/rand.h>
#include <latan/spectrum.h>
#include <latan/statistics.h>

#define A_SPEC_NAME argv[1]
#define A_PART_NAME argv[2]
#define A_SOURCE argv[3]
#define A_SINK argv[4]
#define A_MANF_NAME argv[5]

#define NBOOT 2000

int main(int argc, char* argv[])
{
	/*				parsing arguments			*/
	/********************************************/
	int source,sink;
	stringbuf spec_name,part_name,manf_name;
	
	if (argc != 6)
	{
		fprintf(stderr,														\
				"usage: %s (qcd|qcdqed) particle source sink manifest\n",	\
				argv[0]);
		return EXIT_FAILURE;
	}
	
	strcpy(spec_name,A_SPEC_NAME);
	strcpy(part_name,A_PART_NAME);
	source	= atoi(A_SOURCE);
	sink	= atoi(A_SINK);
	strcpy(manf_name,A_MANF_NAME);
	
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
		abort();
	}
	h = spectrum_get(s,part_name);
	
	/*				loading datas				*/
	/********************************************/
	size_t ndat,nt;
	mat* prop;
	
	ndat	= (size_t)get_nfile(manf_name);
	nt		= (size_t)hadron_getnt(h,source,sink,manf_name);
	
	prop = mat_create_ar(ndat,nt,1);
	
	printf("-- loading %s datas from %s...\n",h->name,manf_name);
	hadron_prop(prop,h,source,sink,manf_name);
	
	/*		resampling mean propagator			*/
	/********************************************/
	rs_sample s_mprop;
	stringbuf sample_name;
	mat mprop;
	
	sprintf(sample_name,"%s_prop_%s",h->name,manf_name);
	s_mprop = rs_sample_create_boot(nt,NBOOT,sample_name);
	
	printf("-- resampling %s mean propagator...\n",h->name);
	randgen_init_from_time();
	resample(s_mprop,prop,ndat,&rs_mean,NULL);
	mprop = rs_sample_get_cent_val(s_mprop);
	
	/*	computing error on mean propagator		*/
	/********************************************/
	mat sig;
	
	sig = mat_create(nt,1);
	
	printf("-- computing %s mean propagator error...\n",h->name);
	rs_sample_varp(sig,s_mprop);
	mat_eqsqrt(sig);
	
	/*				result output				*/
	/********************************************/
	size_t t,maxt;
	
	maxt = (h->parity == EVEN) ? (nt/2) : (nt-1);
	printf("\nt\tprop\t\terror\n");
	for (t=0;t<maxt;t++)
	{
		printf("%d\t%e\t%e\n",(int)t,mat_get(mprop,t,0),mat_get(sig,t,0));
	}
	printf("\n");
	
	/*					plot					*/
	/********************************************/
	plot p;
	stringbuf key;
	const double dmaxt = (double)maxt;
	
	p = plot_create();
	
	sprintf(key,"%s propagator",h->name);
	plot_set_scale_ylog(p);
	plot_set_scale_xmanual(p,0.0,dmaxt);
	plot_add_daterr(p,mprop,sig,0.0,1.0,key);
	plot_disp(p);
	
	/*				desallocation				*/
	/********************************************/
	spectrum_destroy(s);
	mat_destroy_ar(prop,ndat);
	rs_sample_destroy(s_mprop);
	plot_destroy(p);
	
	return EXIT_SUCCESS;
}