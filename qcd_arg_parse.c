#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <sys/types.h>
#include "config.h"
#include "qcd_arg_parse.h"
#include <argtable2.h>
#include <latan/latan_io.h>
#include <latan/latan_math.h>
#include <latan/latan_nunits.h>

#define A_MAX_NERROR 20

qcd_options qcd_arg_parse(int argc, char* argv[], int argset_flag)
{
	/*			argument definitions			*/
	/********************************************/
	void** a_table;
	size_t a_table_size,i;
	struct arg_lit*	a_help			= NULL;
	struct arg_lit* a_ver			= NULL;
	struct arg_int* a_verb			= NULL;
	struct arg_lit* a_save_rs		= NULL;
	struct arg_str* a_load_rg		= NULL;
	struct arg_lit* a_noplot		= NULL;
	struct arg_dbl* a_latspac_fm	= NULL;
	struct arg_str* a_qcomp			= NULL;
	struct arg_str* a_spec			= NULL;
	struct arg_str* a_part			= NULL;
	struct arg_str* a_channel		= NULL;
	struct arg_str* a_ss			= NULL;
	struct arg_int* a_binsize		= NULL;
	struct arg_file* a_manf			= NULL;
	struct arg_end* a_end			= NULL;
	stringbuf help_msg,ver_msg,verb_msg,save_rs_msg,load_rg_msg,noplot_msg,latspac_fm_msg,\
	qcomp_msg,spec_msg,part_msg,channel_msg,ss_msg,binsize_msg,  \
	manf_msg;
	
	sprintf(help_msg		,"display this help"							);
	sprintf(ver_msg			,"display QCDutils version"						);
	sprintf(verb_msg		,"verbosity level (default: 1)"					);
	sprintf(save_rs_msg		,"save resampled samples"						);
	sprintf(load_rg_msg		,"use saved random generator state"				);
	sprintf(noplot_msg		,"disable plot"									);
	sprintf(latspac_fm_msg	,"lattice spacing (in fm)"						);
	sprintf(qcomp_msg		,"quark composition"							);
	sprintf(spec_msg		,"particle spectrum (default: qcd)"				);
	sprintf(part_msg		,"particle name"								);
	sprintf(channel_msg		,"custom channel label"							);
	sprintf(ss_msg			,"particle source and sink"						);
	sprintf(binsize_msg		,"data binning size (default: 1)"				);
	sprintf(manf_msg		,"particle LatAnalyze data manifest"			);
	
	a_help	= arg_lit0(NULL,"help",help_msg);
	a_ver	= arg_lit0(NULL,"version",ver_msg);
	a_verb	= arg_int0("v","verb",NULL,verb_msg);
	if (argset_flag & A_SAVE_RS)
	{
		a_save_rs		= arg_lit0(NULL,"save-rs",save_rs_msg);
	}
	if (argset_flag & A_LOAD_RG)
	{
		a_load_rg		= arg_str0("r","load-rg","FNAME",load_rg_msg);
	}
	if (argset_flag & A_PLOT)
	{	
		a_noplot		= arg_lit0(NULL,"no-plot",noplot_msg);
	}
	if (argset_flag & A_LATSPAC)
	{
		a_latspac_fm	= arg_dbl0("L","lat-spac",NULL,latspac_fm_msg);
	}
	if (argset_flag & A_QCOMP)
	{
		a_qcomp			= arg_str1("m","qu-comp","Q1Q2",qcomp_msg);
	}
	if (argset_flag & A_PARTICLE)
	{
		a_spec			= arg_str0("s","spectrum","{qcd|qcdqed}",spec_msg);
		a_part			= arg_str1("p","particle","NAME",part_msg);
	}
	if (argset_flag & A_CHANNEL)
	{
		a_channel		= arg_strn("c","channel","CH_TYPE:CH_LABEL",0,\
								   NCHANNEL,channel_msg);
	}
	if (argset_flag & A_PROP_LOAD)
	{
		a_ss			= arg_str1("S","source-sink","SOURCE_TYPESINK_TYPE",\
								   ss_msg);
		a_binsize		= arg_int0("b","bin",NULL,binsize_msg);
		a_manf			= arg_file1(NULL,NULL,"<manifest file>",manf_msg);
	}
	a_end	= arg_end(A_MAX_NERROR);
	
	a_table_size = 3;
	MALLOC_ERRVAL(a_table,void**,a_table_size,NULL);
	a_table[0] = a_help;
	a_table[1] = a_ver;
	a_table[2] = a_verb;
	i = 3;
	if (argset_flag & A_SAVE_RS)
	{
		a_table_size++;
		REALLOC_ERRVAL(a_table,a_table,void**,a_table_size,NULL);
		a_table[i] = a_save_rs;
		i++;
	}
	if (argset_flag & A_LOAD_RG)
	{
		a_table_size++;
		REALLOC_ERRVAL(a_table,a_table,void**,a_table_size,NULL);
		a_table[i] = a_load_rg;
		i++;
	}
	if (argset_flag & A_PLOT)
	{
		a_table_size++;
		REALLOC_ERRVAL(a_table,a_table,void**,a_table_size,NULL);
		a_table[i] = a_noplot;
		i++;
	}
	if (argset_flag & A_LATSPAC)
	{
		a_table_size++;
		REALLOC_ERRVAL(a_table,a_table,void**,a_table_size,NULL);
		a_table[i] = a_latspac_fm;
		i++;
	}
	if (argset_flag & A_QCOMP)
	{
		a_table_size++;
		REALLOC_ERRVAL(a_table,a_table,void**,a_table_size,NULL);
		a_table[i] = a_qcomp;
		i++;
	}
	if (argset_flag & A_PARTICLE)
	{
		a_table_size += 2;
		REALLOC_ERRVAL(a_table,a_table,void**,a_table_size,NULL);
		a_table[i] = a_spec;
		a_table[i+1] = a_part;
		i += 2;
	}
	if (argset_flag & A_CHANNEL)
	{
		a_table_size++;
		REALLOC_ERRVAL(a_table,a_table,void**,a_table_size,NULL);
		a_table[i] = a_channel;
		i++;
	}
	if (argset_flag & A_PROP_LOAD)
	{
		a_table_size += 3;
		REALLOC_ERRVAL(a_table,a_table,void**,a_table_size,NULL);
		a_table[i] = a_ss;
		a_table[i+1] = a_binsize;
		a_table[i+2] = a_manf;
		i += 3;
	}
	a_table_size++;
	REALLOC_ERRVAL(a_table,a_table,void**,a_table_size,NULL);
	a_table[i] = a_end;
	
	/*			check arg table allocation		*/
	/********************************************/
	if (arg_nullcheck(a_table) != 0)
	{
		fprintf(stderr,"error: argument table allocation failed\n");
		exit(EXIT_FAILURE);
	}
	
	/*		option structure allocation			*/
	/********************************************/
	qcd_options opt;
	
	MALLOC_ERRVAL(opt,qcd_options,1,NULL);
	
	/*			set default options				*/
	/********************************************/
	if (argset_flag & A_PARTICLE)
	{
		a_spec->sval[0] = "qcd";
	}
	if (argset_flag & A_CHANNEL)
	{
		channel_no ch;
		
		for (ch=0;ch<NCHANNEL;ch++)
		{
			channel_id_get(opt->channel_id[ch],ch);
		}
	}
	if (argset_flag & A_PROP_LOAD)
	{
		a_binsize->ival[0] = 1;
	}
	a_verb->ival[0] = 1;
	
	/*			argument parsing				*/
	/********************************************/
	
	int nerror;
	
	nerror = arg_parse(argc,argv,a_table);
	if ((a_help->count > 0)||(argc == 1))
	{
		printf("usage: %s",argv[0]);
        arg_print_syntax(stdout,a_table,"\n");
        arg_print_glossary_gnu(stdout,a_table);
		exit(EXIT_SUCCESS);
	}
	if (a_ver->count > 0)
	{
		printf("%s is part of %s v%s\n",argv[0],PACKAGE_NAME,PACKAGE_VERSION);
		printf("Report any bugs or issues to %s\n",PACKAGE_BUGREPORT);
		exit(EXIT_SUCCESS);
	}
	if (nerror > 0)
	{
        arg_print_errors(stdout,a_end,argv[0]);
        printf("Try '%s --help' for more information.\n",argv[0]);
        exit(EXIT_FAILURE);
	}
	opt->qcd_verb = (a_verb->ival[0] > 0);
	opt->latan_verb = MAX(0,a_verb->ival[0]-1);
	if (argset_flag & A_SAVE_RS)
	{
		opt->do_save_rs_sample = (a_save_rs->count > 0);
	}
	if (argset_flag & A_LOAD_RG)
	{
		if (a_load_rg->count > 0)
		{
			opt->have_randgen_state = true;
			randgen_load_state(opt->state,a_load_rg->sval[0]);
		}
		else
		{
			opt->have_randgen_state = false;
			randgen_init_from_time();
			randgen_get_state(opt->state);
		}
	}
	if (argset_flag & A_PLOT)
	{
		opt->do_plot = !(a_noplot->count > 0);
	}
	if (argset_flag & A_LATSPAC)
	{
		if (a_latspac_fm->count == 1)
		{
			opt->latspac_fm = a_latspac_fm->dval[0];
			opt->latspac_nu = opt->latspac_fm*NU_FM;
			opt->have_latspac = true;
		}
		else
		{
			opt->latspac_fm = 1.0;
			opt->latspac_nu = 1.0;
			opt->have_latspac = false;
		}
	}
	if (argset_flag & A_QCOMP)
	{
		strcpy(opt->qcomp,a_qcomp->sval[0]);
	}
	if (argset_flag & A_PARTICLE)
	{
		strcpy(opt->spec_name,a_spec->sval[0]);
		strcpy(opt->part_name,a_part->sval[0]);
	}
	if (argset_flag & A_CHANNEL)
	{
		char* pch = NULL;
		stringbuf label;
		stringbuf new_id;
		size_t col_ind;
		int i;
		
		for(i=0;i<a_channel->count;i++)
		{
			pch = strchr(a_channel->sval[i],':');
			if (pch == NULL)
			{
				fprintf(stderr,"error: channel option %s invalid\n",\
						a_channel->sval[i]);
				exit(EXIT_FAILURE);
			}
			col_ind = (size_t)(pch-a_channel->sval[i]);
			strncpy(label,a_channel->sval[i],col_ind);
			strcpy(new_id,a_channel->sval[i]+col_ind+1);
			strcpy(opt->channel_id[channel_no_get(label)],new_id);
		}
	}
	if (argset_flag & A_PROP_LOAD)
	{
		opt->source = ss_no_get(a_ss->sval[0][0]);
		opt->sink = ss_no_get(a_ss->sval[0][1]);
		opt->binsize = (size_t)(a_binsize->ival[0]);
		strcpy(opt->manf_name,a_manf->filename[0]);
	}
	
	/*				desallocation				*/
	/********************************************/
	arg_freetable(a_table,a_table_size);
	FREE(a_table);

	return opt;
}

void qcd_printf(const qcd_options opt, const stringbuf fmt, ...)
{
	va_list args;
	
	if (opt->qcd_verb)
	{
		va_start(args,fmt);
		vprintf(fmt,args);
		va_end(args);
	}
}