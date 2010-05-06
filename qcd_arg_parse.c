#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include "config.h"
#include "qcd_arg_parse.h"
#include <argtable2.h>
#include <latan/latan_nunits.h>

#define A_MAX_NERROR 20
#define A_TABLE_DEF_SIZE 8

qcd_options qcd_arg_parse(int argc, char* argv[], int argset_flag)
{
	/*			argument definitions			*/
	/********************************************/
	void** a_table;
	size_t a_table_size,i;
	struct arg_lit*	a_help;
	struct arg_lit* a_ver;
	struct arg_lit* a_save;
	struct arg_lit* a_noplot;
	struct arg_dbl* a_latspac_fm;
	struct arg_str* a_qcomp;
	struct arg_str* a_spec;
	struct arg_str* a_part;
	struct arg_int* a_source;
	struct arg_int* a_sink;
	struct arg_file* a_manf;
	struct arg_end* a_end;
	stringbuf help_msg,ver_msg,save_msg,noplot_msg,latspac_fm_msg,qcomp_msg,\
	spec_msg,part_msg,source_msg,sink_msg,manf_msg;
	
	sprintf(help_msg		,"display this help"             		);
	sprintf(ver_msg			,"display QCDutils version"         	);
	sprintf(save_msg		,"save resampled samples"				);
	sprintf(noplot_msg		,"disable plot"							);
	sprintf(latspac_fm_msg	,"lattice spacing (in fm)"				);
	sprintf(qcomp_msg		,"quark composition"					);
	sprintf(spec_msg		,"particle spectrum (default: qcd)"		);
	sprintf(part_msg		,"particle name"	             		);
	sprintf(source_msg		,"particle source"	             		);
	sprintf(sink_msg		,"particle sink"	            		);
	sprintf(manf_msg		,"particle LatAnalyze data manifest"	);
	
	a_help			= arg_lit0(NULL,"help",help_msg);
	a_ver			= arg_lit0(NULL,"version",ver_msg);
	a_save			= arg_lit0(NULL,"save",save_msg);
	a_noplot		= arg_lit0(NULL,"noplot",noplot_msg);
	a_latspac_fm	= arg_dbl0("L","latspac",NULL,latspac_fm_msg);
	a_qcomp			= arg_str1(NULL,"qcomp","q1q2",qcomp_msg);
	a_spec			= arg_str0("s","spectrum","{qcd|qcdqed}",spec_msg);
	a_part			= arg_str1("p","particle","NAME",part_msg);
	a_source		= arg_int1("A","source",NULL,source_msg);
	a_sink			= arg_int1("B","sink",NULL,sink_msg);
	a_manf			= arg_file1(NULL,NULL,"<manifest file>",manf_msg);
	a_end			= arg_end(A_MAX_NERROR);
	
	a_table_size = A_TABLE_DEF_SIZE;
	MALLOC_ERRVAL(a_table,void**,a_table_size,NULL);
	a_table[0] = a_help;
	a_table[1] = a_ver;
	a_table[2] = a_save;
	a_table[3] = a_noplot;
	i = 4;
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
	a_table[i] = a_source;
	a_table[i+1] = a_sink;
	a_table[i+2] = a_manf;
	a_table[i+3] = a_end;
	
	/*			check arg table allocation		*/
	/********************************************/
	if (arg_nullcheck(a_table) != 0)
	{
		fprintf(stderr,"error: argument table allocation failed\n");
		exit(EXIT_FAILURE);
	}
	
	/*			set default options				*/
	/********************************************/
	a_spec->sval[0] = "qcd";
	
	/*			argument parsing				*/
	/********************************************/
	qcd_options opt;
	int nerror;
	
	MALLOC_ERRVAL(opt,qcd_options,1,NULL);
	
	nerror = arg_parse(argc,argv,a_table);
	if ((a_help->count > 0)||(argc == 1))
	{
		printf("usage: %s",argv[0]);
        arg_print_syntax(stdout,a_table,"\n");
        arg_print_glossary(stdout,a_table,"  %-30s %s\n");
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
	opt->do_save_rs_sample = (a_save->count > 0);
	opt->do_plot = !(a_noplot->count > 0);
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
			fprintf(stderr,"warning: no lattice spacing specified, results will be in lattice unit\n");
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
	opt->source = a_source->ival[0];
	opt->sink = a_sink->ival[0];
	strcpy(opt->manf_name,a_manf->filename[0]);
	
	/*				desallocation				*/
	/********************************************/
	arg_freetable(a_table,a_table_size);
	FREE(a_table);

	return opt;
}