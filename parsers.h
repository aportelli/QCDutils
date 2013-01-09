#ifndef QCD_PARSERS_H_
#define QCD_PARSERS_H_

#include <latan/latan_io.h>

#define ATOI(str) ((int)strtol(str,(char **)NULL,10))
#define ATOF(str) (strtod(str,(char **)NULL))

typedef struct ukhadron_par_s
{
    int quark[2];
    int ss[2];
    strbuf type;
    endian_no in_endian;
} ukhadron_par;

typedef void parse_func(strbuf, io_fmt_no, void *);

parse_func parse_bmw_dynqcd;
parse_func parse_ukhadron_mes;
parse_func parse_ukhadron_bar;

#endif
