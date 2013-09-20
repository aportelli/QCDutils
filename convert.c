#include <parsers.h>

#define MALLOC(pt,typ,size)\
{\
    pt = (typ)(malloc((size_t)(size)*sizeof(*pt)));\
    if (pt == NULL)\
    {\
        fprintf(stderr,"error: memory allocation failed\n");\
        exit(EXIT_FAILURE);\
    }\
}

#define FREE(pt)\
{\
    if (pt != NULL)\
    {\
        free(pt);\
        pt = NULL;\
    }\
}

int main(int argc, char *argv[])
{
    strbuf man_fname, *meas_fname, fname;
    io_fmt_no io_fmt;
    int i,m,count,dumb;
    int nmeas;
    parse_func *parser;
    void *par;
    ukhadron_par uk_par;
    
    if (argc >= 4)
    {
        strbufcpy(man_fname,argv[1]);
        io_fmt = io_get_fmt();
        if (strbufcmp(argv[2],"xml") == 0)
        {
            io_fmt = IO_XML;
        }
        else if (strbufcmp(argv[2],"ascii") == 0)
        {
            io_fmt = IO_ASCII;
        }
        else
        {
            fprintf(stderr,"format unknown\n");
            return EXIT_FAILURE;
        }
        io_set_fmt(io_fmt);
        if (strbufcmp(argv[3],"bmw_dynqcd") == 0)
        {
            parser = &parse_bmw_dynqcd;
        }
        else if ((strbufcmp(argv[3],"ukqcd_mes") == 0)  \
                 ||(strbufcmp(argv[3],"ukqcd_bar") == 0)\
                 ||(strbufcmp(argv[3],"ukqcd_rarekaon") == 0)\
                 ||(strbufcmp(argv[3],"ukqcd_hvp") == 0))
        {
            if (strbufcmp(argv[3],"ukqcd_mes") == 0)
            {
                parser = &parse_ukhadron_mes;
            }
            else if (strbufcmp(argv[3],"ukqcd_bar") == 0)
            {
                parser = &parse_ukhadron_bar;
            }
            else if (strbufcmp(argv[3],"ukqcd_rarekaon") == 0)
            {
                parser = &parse_ukhadron_rarekaon;
            }
            else if (strbufcmp(argv[3],"ukqcd_hvp") == 0)
            {
                parser = &parse_ukhadron_hvp;
            }
            par    = &uk_par;
            if (argc >= 8)
            {
                uk_par.quark[0] = ATOI(argv[4]);
                uk_par.quark[1] = ATOI(argv[5]);
                uk_par.ss[0]    = ATOI(argv[6]);
                uk_par.ss[1]    = ATOI(argv[7]);
                if (argc == 9)
                {
                    if (strbufcmp(argv[8],"le") == 0)
                    {
                        uk_par.in_endian = LE;
                    }
                    else if (strbufcmp(argv[8],"be") == 0)
                    {
                        uk_par.in_endian = BE;
                    }
                    else
                    {
                        fprintf(stderr,"error: endianness option \"%s\" unknown\n",
                                argv[8]);
                        return EXIT_FAILURE;
                    }
                }
                else
                {
                    uk_par.in_endian = latan_get_endianness();
                }
            }
            else
            {
                fprintf(stderr,"error: wrong options for UKhadron parser\n");
                fprintf(stderr,"expected: <q1_no> <q2_no> <sink_no> <source_no> [{le|be}]\n");
                return EXIT_FAILURE;
            }
        }
        else
        {
            fprintf(stderr,"error: parser unknown\n");
            return EXIT_FAILURE;
        }
    }
    else
    {
        fprintf(stderr,"usage: %s <manifest> {xml|ascii} <parser> [<parser_arg_1> ...]\n",\
                argv[0]);
        return EXIT_FAILURE;
    }
    
    m          = 0;
    count      = 0;
    nmeas      = get_nfile(man_fname);
    
    MALLOC(meas_fname,strbuf *,nmeas);
    
    io_init();
    BEGIN_FOR_LINE(fname,man_fname,dumb)
    {
        strbufcpy(meas_fname[m],fname);
        m++;
    }
    END_FOR_LINE
#ifdef _OPENMP
    switch (io_fmt)
    {
        case IO_XML:
            omp_set_num_threads(3*omp_get_num_procs());
            break;
        case IO_ASCII:
            omp_set_num_threads(omp_get_num_procs());
            break;
        default:
            fprintf(stderr,"format unknown\n");
            exit(EXIT_FAILURE);
    }
#pragma omp parallel for
#endif
    for (m=0;m<nmeas;m++)
    {
        parser(meas_fname[m],io_fmt,par);
#ifdef _OPENMP
#pragma omp barrier
#pragma omp critical
#endif
        {
            count++;
            printf("[");
            for (i=0;i<60*count/nmeas;i++)
            {
                printf("=");
            }
            for (i=60*count/nmeas;i<60;i++)
            {
                printf(" ");
            }
            printf("]  %d/%d\r",count,nmeas);
            fflush(stdout);
        }
    }
    printf("\n");
    io_finish();
    
    FREE(meas_fname);
    
    return EXIT_SUCCESS;
}
