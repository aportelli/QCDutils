#ifndef QCD_ARG_PARSE_H_
#define QCD_ARG_PARSE_H_

#include <latan/latan_hadron.h>

enum
{
	A_DEFAULT	= 0x0,
	A_LATSPAC	= 0x1,
	A_PARTICLE	= 0x2,
	A_QCOMP		= 0x4,
	A_PROP_FIT	= 0x8,
	A_CHI_FIT	= 0x10
};

typedef struct
{
	double latspac_fm;
	double latspac_nu;
	bool have_latspac;
	stringbuf qcomp;
	stringbuf spec_name;
	stringbuf part_name;
	stringbuf manf_name;
	int source;
	int sink;
	stringbuf channel_id[NCHANNEL];
	stringbuf quark_id[NQUARK];
	bool do_plot;
	bool do_save_rs_sample;
}* qcd_options;

qcd_options qcd_arg_parse(int argc, char* argv[], int argset_flag);

#endif