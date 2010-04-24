#ifndef QCD_ARG_PARSE_H_
#define QCD_ARG_PARSE_H_

#include <latan/spectrum.h>

typedef struct
{
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

qcd_options qcd_arg_parse(int argc, char* argv[]);

#endif