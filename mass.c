#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <qcd_arg_parse.h>
#include <latan/latan_hadron.h>
#include <latan/latan_io.h>
#include <latan/latan_mat.h>
#include <latan/latan_math.h>
#include <latan/latan_mass.h>
#include <latan/latan_minimizer.h>
#include <latan/latan_models.h>
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
	double latspac_nu;
	qcd_options opt;
	stringbuf spec_name,part_name,manf_name,unit;
	size_t binsize;
	
	opt = qcd_arg_parse(argc,argv,A_PARTICLE|A_PROP_LOAD|A_LATSPAC|A_PROP_FIT\
						|A_SAVE_RS|A_PLOT|A_CHANNEL|A_LOAD_RG);
	strcpy(spec_name,opt->spec_name);
	strcpy(part_name,opt->part_name);
	strcpy(manf_name,opt->manf_name);
	source = opt->source;
	sink = opt->sink;
	binsize = opt->binsize;
	latspac_nu = opt->latspac_nu;
	if (opt->have_latspac)
	{
		strcpy(unit," (MeV)");
	}
	else
	{
		strcpy(unit,"");
	}
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
	
	qcd_printf(opt,"-- loading %s datas from %s...\n",h->name,manf_name);
	hadron_propbin(prop,h,source,sink,manf_name,binsize);
	
	/*		resampling mean propagator			*/
	/********************************************/
	rs_sample s_mprop;
	stringbuf sample_name;
	mat mprop;
	
	sprintf(sample_name,"%s_prop_%s",h->name,manf_name);
	s_mprop = rs_sample_create_boot(nt,NBOOT,sample_name);
	
	qcd_printf(opt,"-- resampling %s mean propagator...\n",h->name);
	randgen_set_state(opt->state);
	resample(s_mprop,prop,ndat,1,&rs_mean,NULL);
	mprop = rs_sample_pt_cent_val(s_mprop);
	
	/*	computing error on mean propagator		*/
	/********************************************/
	mat sigmprop;
	
	sigmprop = mat_create(nt,1);
	
	qcd_printf(opt,"-- estimating %s mean propagator variance...\n",h->name);
	rs_sample_varp(sigmprop,s_mprop);
	
	/*		resampling effective mass			*/
	/********************************************/
	rs_sample s_effmass;
	mat em;
	
	sprintf(sample_name,"%s_effmass_%s",h->name,manf_name);
	s_effmass = rs_sample_create_boot(nt-2,NBOOT,sample_name);
	
	qcd_printf(opt,"-- resampling %s effective mass...\n",h->name);
	randgen_set_state(opt->state);
	resample(s_effmass,prop,ndat,1,&rs_effmass,&(h->parity));
	em = rs_sample_pt_cent_val(s_effmass);
	
	/*	computing error on effective mass		*/
	/********************************************/
	mat sigem;
	
	sigem = mat_create(nt-2,1);
	
	qcd_printf(opt,"-- computing %s effective mass variance...\n",h->name);
	rs_sample_varp(sigem,s_effmass);
	mat_eqsqrt(sigem);
	
	/*				fitting mass				*/
	/********************************************/
	fit_data d;
	rs_sample s_mass;
	mat mass;
	size_t i;
	
	d = fit_data_create(nt,1);
	sprintf(sample_name,"%s_massfit_%s",h->name,manf_name);
	s_mass = rs_sample_create_boot(2,NBOOT,sample_name);
	
	qcd_printf(opt,"-- fitting and resampling %s mass...\n",h->name);
	fit_data_set_var(d,sigmprop);
	fit_data_mass_fit_tune(d,rs_sample_pt_cent_val(s_mass),mprop,em,sigem,\
						   h->parity);
	for (i=0;i<rs_sample_get_nsample(s_mass);i++)
	{
		mat_cp(rs_sample_pt_sample(s_mass,i),rs_sample_pt_cent_val(s_mass));
	}
	randgen_set_state(opt->state);
	resample(s_mass,prop,ndat,1,&rs_data_fit,d);
	mass = rs_sample_pt_cent_val(s_mass);
	if (opt->do_save_rs_sample)
	{
		rs_sample_save(s_mass,s_mass->name);
	}
	
	/*		computing error on mass				*/
	/********************************************/
	mat sigmass;
	
	sigmass = mat_create(2,1);
	
	qcd_printf(opt,"-- estimating %s mass variance...\n",h->name);
	rs_sample_varp(sigmass,s_mass);
	mat_eqsqrt(sigmass);
	mat_eqsqrt(sigmprop);
	
	/*		switching to right units			*/
	/********************************************/
	mat_eqmuls(em,1.0/latspac_nu);
	mat_eqmuls(sigem,1.0/latspac_nu);
	
	/*				result output				*/
	/********************************************/
	qcd_printf(opt,"M_%s\t\t= %.8f +/- %.8e %s\n",h->name,					\
			   mat_get(mass,0,0)/latspac_nu,mat_get(sigmass,0,0)/latspac_nu,\
			   unit);
	qcd_printf(opt,"dof\t\t= %d\n",fit_data_get_dof(d));
	qcd_printf(opt,"chi^2/dof\t= %e\n",fit_data_get_chi2pdof(d));
	
	/*					plot					*/
	/********************************************/
	if (opt->do_plot)
	{
		mat em_t;
		plot p;
		stringbuf key,plotcmd,plotfmt;
		size_t maxt;
		double dmaxt,shift;
		
		maxt  = (h->parity == EVEN) ? (nt/2) : (nt-1);
		dmaxt = (double)maxt;
		shift = (h->parity == EVEN) ? (0.0) : (-DRATIO(nt,2));
		
		em_t = mat_create(nt-2,1);
		
		/* propagator plot */
		p = plot_create();
		
		plot_set_scale_ylog(p);
		plot_set_scale_xmanual(p,shift,dmaxt+shift);
		sprintf(key,"%s propagator",h->name);
		plot_add_daterr(p,fit_data_pt_x(d),mprop,sigmprop,key);
		fit_model_get_plot_fmt(plotfmt,fit_data_pt_model(d));
		sprintf(plotcmd,plotfmt,mat_get(mass,0,0),mat_get(mass,1,0));
		strcat(plotcmd," t \"fit\"");
		plot_add_plot(p,plotcmd);
		plot_disp(p);	
		
		plot_destroy(p);
		
		/* effective mass plot */
		p = plot_create();
		
		plot_set_scale_xmanual(p,0.0,dmaxt);
		plot_add_hlineerr(p,mat_get(mass,0,0)/latspac_nu,\
						  mat_get(sigmass,0,0)/latspac_nu,"0",\
						  "rgb \"red\"","rgb \"light-red\"");
		sprintf(key,"%s effective mass",h->name);
		mat_set_step(em_t,1.0,1.0);
		plot_add_daterr(p,em_t,em,sigem,key);
		plot_disp(p);
		
		plot_destroy(p);
	}
	
	/*				desallocation				*/
	/********************************************/
	FREE(opt);
	spectrum_destroy(s);
	mat_ar_destroy(prop,ndat);
	rs_sample_destroy(s_mprop);
	mat_destroy(sigmprop);
	rs_sample_destroy(s_effmass);
	mat_destroy(sigem);
	fit_data_destroy(d);
	rs_sample_destroy(s_mass);
	mat_destroy(sigmass);
	
	return EXIT_SUCCESS;
}