#define _POSIX_C_SOURCE 199506L /* strtok_r is used here */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <libgen.h>
#include <sys/types.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <latan/latan_io.h>
#include <latan/latan_math.h>
#include <parsers.h>

#define MAX_PROP_NT 1024

void parse_bmw_dynqcd(strbuf meas_fname, io_fmt_no fmt, void *nothing __dumb)
{
    strbuf out_fname,prop_name,rsource,rsink;
    bool is_in_prop;
    strbuf channel;
    int q1,q2;
    size_t nt;
    int nf,lc;
    mat *prop;
    double prop_buf[MAX_PROP_NT];
    strbuf *field;
    
    field      = NULL;
    is_in_prop = false;
    nt         = 0;
    q1         = 0;
    q2         = 0;
    strbufcpy(rsource,"");
    strbufcpy(rsink,"");
    
    switch (fmt)
    {
        case IO_XML:
            sprintf(out_fname,"latan.%s.xml",basename(meas_fname));
            break;
        case IO_ASCII:
            sprintf(out_fname,"latan.%s.dat",basename(meas_fname));
            break;
        default:
            fprintf(stderr,"format unknown\n");
            exit(EXIT_FAILURE);
    }
    if (access(out_fname,F_OK) != 0)
    {
        BEGIN_FOR_LINE_TOK(field,meas_fname," ",nf,lc)
        {
            if (strbufcmp(field[0],"START_PROP") == 0)
            {
                strbufcpy(channel,field[1]);
                q1         = ATOI(field[4]);
                q2         = ATOI(field[6]);
                strbufcpy(rsource,field[8]);
                strbufcpy(rsink,field[10]);
                is_in_prop = true;
                nt         = 0;
            }
            else if (is_in_prop)
            {
                if (strbufcmp(field[0],"END_PROP") == 0)
                {
                    prop = mat_create_from_ar(prop_buf,nt,1);
                    
                    sprintf(prop_name,"%s%c%s_%d%d_%s_%s",         \
                            out_fname,LATAN_PATH_SEP,channel,q1,q2,\
                            rsink,rsource);
                    mat_save(prop_name,'a',prop);
                    mat_destroy(prop);
                    is_in_prop = false;
                }
                else if (sscanf(field[1],"%lf",prop_buf+nt) > 0)
                {
                    nt++;
                }
                else
                {
                    fprintf(stderr,"error (%s:%d): error while parsing propagator\n",\
                            meas_fname,lc);
                    exit(EXIT_FAILURE);
                }
            }
        }
        END_FOR_LINE_TOK(field)
    }
}

void parse_ukhadron_mes(strbuf meas_fname, io_fmt_no fmt, void *par_v)
{
    ukhadron_par *par;
    int ibuf[10],mom[3],mom_i,gsrc,gsink,t;
    int nmom,ngsrc,ngsink,nt;
    double *dbuf;
    mat *prop;
    FILE *in_f;
    size_t rcount,ind,x[3],dim[3];
    strbuf out_fname,prop_name;
    
    par  = (ukhadron_par *)par_v;
    
    switch (fmt)
    {
        case IO_XML:
            sprintf(out_fname,"latan.%s.xml",basename(meas_fname));
            break;
        case IO_ASCII:
            sprintf(out_fname,"latan.%s.dat",basename(meas_fname));
            break;
        default:
            fprintf(stderr,"format unknown\n");
            exit(EXIT_FAILURE);
    }
    in_f = fopen(meas_fname,"r");
    if (in_f == NULL)
    {
        fprintf(stderr,"error: cannot open file %s\n",meas_fname);
        abort();
    }
    rcount = fread(ibuf,sizeof(int),1,in_f);
    if (rcount != 1)
    {
        fprintf(stderr,"error (%s:%lx): reading momentum number failed\n",\
                meas_fname,ftell(in_f));
        abort();
    }
    nmom   = latan_conv_endianness_i(ibuf[0],par->in_endian);
    for (mom_i=0;mom_i<nmom;mom_i++)
    {
        rcount = fread(ibuf,sizeof(int),10,in_f);
        if (rcount != 10)
        {
            fprintf(stderr,"error (%s:%lx): reading correlator header failed\n",\
                    meas_fname,ftell(in_f));
            abort();
        }
        mom[0] = latan_conv_endianness_i(ibuf[4],par->in_endian);
        mom[1] = latan_conv_endianness_i(ibuf[5],par->in_endian);
        mom[2] = latan_conv_endianness_i(ibuf[6],par->in_endian);
        ngsrc  = latan_conv_endianness_i(ibuf[7],par->in_endian);
        ngsink = latan_conv_endianness_i(ibuf[8],par->in_endian);
        nt     = latan_conv_endianness_i(ibuf[9],par->in_endian);
        dim[0] = (size_t)nt;
        dim[1] = (size_t)ngsink;
        dim[2] = (size_t)ngsrc;
        dbuf   = (double *)malloc(2*nt*ngsrc*ngsink*sizeof(double));
        prop   = mat_create(nt,1);
        rcount = fread(dbuf,sizeof(double),2*nt*ngsrc*ngsink,in_f);
        if (rcount != (size_t)(2*nt*ngsrc*ngsink))
        {
            fprintf(stderr,"error (%s:%lx): reading correlator failed\n",\
                    meas_fname,ftell(in_f));
            abort();
        }
        for (gsrc=0;gsrc<ngsrc;gsrc++)
        for (gsink=0;gsink<ngsink;gsink++)
        {
            x[1] = (size_t)gsink;
            x[2] = (size_t)gsrc;
            for (t=0;t<nt;t++)
            {
                x[0] = (size_t)t;
                ind  = coord_to_rowmaj(x,dim,3);
                mat_set(prop,(size_t)t,0,                                     \
                        latan_conv_endianness_d(dbuf[2*ind],par->in_endian));
            }
            sprintf(prop_name,"%s%cG%dG%d_%d%d%d_%d%d_%d_%d",out_fname,\
                    LATAN_PATH_SEP,gsink,gsrc,mom[0],mom[1],mom[2],    \
                    par->quark[0],par->quark[1],par->ss[0],par->ss[1]);
            mat_save(prop_name,'a',prop);
        }
        free(dbuf);
        mat_destroy(prop);
    }
    fclose(in_f);
}

void parse_ukhadron_bar(strbuf meas_fname, io_fmt_no fmt, void *par_v)
{
    ukhadron_par *par;
    int ibuf[7],mom[3],mom_i,op,t;
    int nmom,nop,nt;
    double *dbuf;
    mat *prop;
    FILE *in_f;
    size_t rcount,ind,x[2],dim[2];
    strbuf out_fname,prop_name;
    
    par  = (ukhadron_par *)par_v;
    
    switch (fmt)
    {
        case IO_XML:
            sprintf(out_fname,"latan.%s.xml",basename(meas_fname));
            break;
        case IO_ASCII:
            sprintf(out_fname,"latan.%s.dat",basename(meas_fname));
            break;
        default:
            fprintf(stderr,"format unknown\n");
            exit(EXIT_FAILURE);
    }
    in_f = fopen(meas_fname,"r");
    if (in_f == NULL)
    {
        fprintf(stderr,"error: cannot open file %s\n",meas_fname);
        abort();
    }
    rcount = fread(ibuf,sizeof(int),1,in_f);
    if (rcount != 1)
    {
        fprintf(stderr,"error (%s:%lx): reading momentum number failed\n",\
                meas_fname,ftell(in_f));
        abort();
    }
    nmom   = latan_conv_endianness_i(ibuf[0],par->in_endian);
    for (mom_i=0;mom_i<nmom;mom_i++)
    {
        rcount = fread(ibuf,sizeof(int),7,in_f);
        if (rcount != 7)
        {
            fprintf(stderr,"error (%s:%lx): reading correlator header failed\n",\
                    meas_fname,ftell(in_f));
            abort();
        }
        mom[0] = latan_conv_endianness_i(ibuf[2],par->in_endian);
        mom[1] = latan_conv_endianness_i(ibuf[3],par->in_endian);
        mom[2] = latan_conv_endianness_i(ibuf[4],par->in_endian);
        nop    = latan_conv_endianness_i(ibuf[5],par->in_endian);
        nt     = latan_conv_endianness_i(ibuf[6],par->in_endian);
        dim[0] = (size_t)nt;
        dim[1] = (size_t)nop;
        dbuf   = (double *)malloc(2*nt*nop*sizeof(double));
        prop   = mat_create(nt,1);
        rcount = fread(dbuf,sizeof(double),2*nt*nop,in_f);
        if (rcount != (size_t)(2*nt*nop))
        {
            fprintf(stderr,"error (%s:%lx): reading correlator failed\n",\
                    meas_fname,ftell(in_f));
            abort();
        }
        for (op=0;op<nop;op++)
        {
            x[1] = (size_t)op;
            for (t=0;t<nt;t++)
            {
                x[0] = (size_t)t;
                ind  = coord_to_rowmaj(x,dim,2);
                mat_set(prop,(size_t)t,0,                                   \
                        latan_conv_endianness_d(dbuf[2*ind],par->in_endian));
            }
            sprintf(prop_name,"%s%cBAR%d_%d%d%d_%d%d_%d_%d",out_fname,\
                    LATAN_PATH_SEP,op,mom[0],mom[1],mom[2],           \
                    par->quark[0],par->quark[1],par->ss[0],par->ss[1]);
            mat_save(prop_name,'a',prop);
        }
        free(dbuf);
        mat_destroy(prop);
    }
    fclose(in_f);
}
