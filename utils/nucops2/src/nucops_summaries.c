#include <stdlib.h>
#include <string.h>
#include "argparser.h"
#include "sequence_base.h"
#include "textparsing.h"
#include "vecops.h"

args_t* nucops_fastasummary_init_args(int argc, char** argv) {
    args_t* result;
    result = args_alloc();
    args_add(result, "input", 'i', "str");
    args_add(result, "output", 'o', "str");
    args_add(result, "minlen", 'm', "");
    args_add(result, "maxlen", 'M', "");
    args_add(result, "nseq", 'n', "");
    args_add(result, "totlen", 't', "");
    args_add(result, "GC", 'g', "");
    args_add(result, "Nxx", 'N', "int");
    args_add(result, "Lxx", 'L', "int");
    args_add_help(result,NULL,"Positional argument","A single positional value is accepted, and should be the file name.", "treated as -i");
    args_parse(result, argc, argv);
    return result;
}

int nucops_fastasummary_a(args_t* args) {
    char* inputfn;
    char* outputfn;
    PF_t* f_in;
    PF_t* f_out;
    long long resume_at;

    int defaultprint;
    int _m, _M, _t, _n, _g, _N, _L;
    int Nxx, Lxx;
    size_t A,C,G,T,N,tmp1,tmp2,tmp3,tmp4,tmp5;
    int64_t* all_lengths;
    size_t nalloclengths, curseqid;

    contig_s** inseqs;
    size_t minlen;
    size_t maxlen;
    size_t totlen;
    size_t nseq, nseqtot;
    size_t i;
    size_t curlen;
    int failure_code;
    int stdin_input;

    stdin_input = 0;
    inputfn = args_getstr(args, "input", 0, args_getstr(args, NULL, 0, "stdin"));
    if (strcmp(inputfn,"stdin")==0){
        stdin_input = 1;
        args_report_warning(NULL,"When reading stdin, all contigs are copied to memory before being parsed - please ensure that your file fits within half of available RAM before doing this\n");
    }
    outputfn = args_getstr(args, "output", 0, "stdout");
    defaultprint = 1;
    _m = args_ispresent(args, "minlen"); if (_m)defaultprint = 0;
    _M = args_ispresent(args, "maxlen"); if (_M)defaultprint = 0;
    _t = args_ispresent(args, "totlen"); if (_t)defaultprint = 0;
    _n = args_ispresent(args, "nseq"); if (_n)defaultprint = 0;
    _g = args_ispresent(args, "GC"); if (_g)defaultprint = 0;
    _N = args_ispresent(args, "Nxx"); if (_N)defaultprint = 0;
    _L = args_ispresent(args, "Lxx"); if (_L)defaultprint = 0;
    Nxx = args_getint(args,"Nxx",0,50);
    if(Nxx>99){
       args_report_warning(NULL, "N%d is not valid and has been transformed into to N99",Nxx);
       Nxx=99;
    }
    if(Nxx<1){
       args_report_warning(NULL, "N%d is not valid and has been transformed into to N1",Nxx);
       Nxx=1;
    }
    Lxx = args_getint(args,"Lxx",0,50);
    if(Lxx>99){
       args_report_warning(NULL, "L%d is not valid and has been transformed into to L99",Lxx);
       Lxx=99;
    }
    if(Lxx<1){
       args_report_warning(NULL, "L%d is not valid and has been transformed into to L1",Lxx);
       Lxx=1;
    }
    args_report_info(NULL, "Parsing of arguments finished.\n");

    nalloclengths = 32;
    all_lengths = malloc(sizeof(int64_t)*nalloclengths);
    failure_code = 0;
    f_in = NULL;
    f_out = NULL;
    if (!failure_code && PFopen(&f_in, inputfn, "rb")) {
        args_report_error(NULL, "The input file could not be opened.\n");
        failure_code = 1;
    }

    if (!failure_code && PFopen(&f_out, outputfn, "wb")) {
        args_report_error(NULL, "The output file could not be opened.\n");
        failure_code = 2;
    }
    if (!failure_code){
        args_report_info(NULL, "Input and output files successfully opened.\n");
        resume_at = 0;
        if(stdin_input){
            inseqs = contigs_from_fastxPF_limited(f_in, &nseq, 0, &resume_at);
        }
        else {
            inseqs = contigs_from_fastxPF_limited(f_in, &nseq,0xC0000000, &resume_at);
        }
        if (nseq < 1) {
            args_report_error(NULL, "The input file does not contain any sequences.\n");
            failure_code = 1;
        }
    }
    if (!failure_code){
        nseqtot = 0;
        curseqid = 0;
        do {
            nseqtot += nseq;
            args_report_progress(NULL, _LLD_ " sequences read\r", nseqtot);
            curlen = contig_length(inseqs[0]);
            minlen = curlen;
            maxlen = minlen;
            A = C = G = T = N = 0;
            if(_g || defaultprint){
                for (i = 0;i < nseq;i++) {
                    contig_countACGTN(inseqs[i],&tmp1,&tmp2,&tmp3,&tmp4,&tmp5);
                    A += tmp1;
                    C += tmp2;
                    G += tmp3;
                    T += tmp4;
                    N += tmp5;
                }
            }
            totlen = 0;
            for (i = 0;i < nseq;i++) {
                curlen = contig_length(inseqs[i]);
                totlen += curlen;
                if (curlen > maxlen) maxlen = curlen;
                if (curlen < minlen) minlen = curlen;
                contig_free(inseqs[i]);
                all_lengths[curseqid] = (int64_t)curlen;
                curseqid++;
                if(curseqid == nalloclengths){
                    if(nalloclengths < 0x10000000) nalloclengths*=2;
                    else nalloclengths += 0x10000000;
                    all_lengths = realloc(all_lengths, sizeof(int64_t)*nalloclengths);
                }
            }
            if(inseqs)free(inseqs);
            inseqs = contigs_from_fastxPF_limited(f_in, &nseq,0xC0000000, &resume_at);
        } while(resume_at>0);
        if(nseqtot != curseqid){
            args_report_warning(NULL,"Sequences were not properly counted, somehow...");
        }
        args_report_progress(NULL, _LLD_ "\nFile read completely\n", nseqtot);
        if (defaultprint) {
            PFprintf(f_out, "Shortest sequence:\t" _LLD_ "\n", (int64_t)minlen);
            PFprintf(f_out, "Longest sequence:\t" _LLD_ "\n", (int64_t)maxlen);
            PFprintf(f_out, "Total length:\t" _LLD_ "\n", (int64_t)totlen);
            PFprintf(f_out, "Number of sequences:\t" _LLD_ "\n", (int64_t)nseqtot);
            PFprintf(f_out, "GC content:\t%.4f\n", (double)(G+C)/(double)(A+G+C+T+N) );
            vec_sorti64(all_lengths,nseqtot);
            curseqid = nseqtot;
            curlen = 0;
            while(curlen < (totlen * Nxx)/100 && curseqid > 0){
              curseqid--;
              curlen += all_lengths[curseqid];
            }
            PFprintf(f_out, "N" _LLD_ ":\t" _LLD_ "\n", Nxx, all_lengths[curseqid] );
            curseqid = nseqtot;
            curlen = 0;
            while(curlen < (totlen * Lxx)/100 && curseqid > 0){
              curseqid--;
              curlen += all_lengths[curseqid];
            }
            PFprintf(f_out, "L" _LLD_ ":\t" _LLD_ "\n", Lxx, (int64_t)(nseqtot-curseqid) );
        }
        else {
            i = 0;
            if (_m) { PFprintf(f_out, _LLD_, minlen);i = 1; }
            if (_M) { if (i)PFputc('\t', f_out); PFprintf(f_out, _LLD_, (int64_t)maxlen);i = 1; }
            if (_t) { if (i)PFputc('\t', f_out); PFprintf(f_out, _LLD_, (int64_t)totlen);i = 1; }
            if (_n) { if (i)PFputc('\t', f_out); PFprintf(f_out, _LLD_, (int64_t)nseqtot);i = 1; }
            if (_g) { if (i)PFputc('\t', f_out); PFprintf(f_out, "%.4f", (double)(G+C)/(double)(A+G+C+T+N) );i = 1; }
            if (_N || _L) {vec_sorti64(all_lengths,nseqtot);}
            if (_N) {
                curseqid = nseqtot;
                curlen = 0;
                while(curlen < (totlen * Nxx)/100 && curseqid > 0){
                   curseqid--;
                   curlen += all_lengths[curseqid];
                }
                if (i)PFputc('\t', f_out);
                PFprintf(f_out, _LLD_ , all_lengths[curseqid] );
                i = 1;
            }
            if (_L) {
                curseqid = nseqtot;
                curlen = 0;
                while(curlen < (totlen * Lxx)/100 && curseqid > 0){
                   curseqid--;
                   curlen += all_lengths[curseqid];
                }
                if (i)PFputc('\t', f_out);
                PFprintf(f_out, _LLD_ , (int64_t)(nseqtot-curseqid) );
                i = 1;
            }

            if (i) PFputc('\n', f_out);
        }
    }
    free(all_lengths);
    PFclose(f_in);
    PFclose(f_out);
    if(inseqs)free(inseqs);
    return failure_code;
}


int nucops_fastasummary(int argc, char** argv) {
    int result;
    args_t* args;

    args = nucops_fastasummary_init_args(argc, argv);
    result = 0;
    if (!args_ispresent(args, "help")) {
        result = nucops_fastasummary_a(args);
        if (result != 0) {
            args_report_error(args, "Fastasummary failed with code <%d>\n", result);
        }
    }
    args_free(args);
    return result;
}
