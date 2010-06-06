#ifndef QCD_ARG_PARSE_H_
#define QCD_ARG_PARSE_H_

#include <latan/latan_hadron.h>
#include <latan/latan_rand.h>

enum
{
	A_DEFAULT	= 0x0,
	A_PLOT		= 0x1,
	A_SAVE_RS	= 0x2,
	A_LOAD_RG	= 0x4,
	A_PROP_LOAD	= 0x8,
	A_LATSPAC	= 0x10,
	A_PARTICLE	= 0x20,
	A_CHANNEL	= 0x40,
	A_QCOMP		= 0x80,
	A_PROP_FIT	= 0x100,
	A_CHI_FIT	= 0x200
};

typedef struct
{
	int qcd_verb;
	int latan_verb;
	double latspac_fm;
	double latspac_nu;
	bool have_latspac;
	stringbuf qcomp;
	stringbuf spec_name;
	stringbuf part_name;
	stringbuf manf_name;
	randgen_state state;
	ss_no source;
	ss_no sink;
	size_t binsize;
	stringbuf channel_id[NCHANNEL];
	stringbuf quark_id[NQUARK];
	bool do_plot;
	bool do_save_rs_sample;
	bool have_randgen_state;
}* qcd_options;

qcd_options qcd_arg_parse(int argc, char* argv[], int argset_flag);
void qcd_printf(const qcd_options opt, const stringbuf fmt, ...);

#endif