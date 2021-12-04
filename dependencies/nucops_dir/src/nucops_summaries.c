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
    args_add(result, "avglen", 'a', "");
    args_add(result, "medlen", 'e', "");
    args_add(result, "lenhistogram", 'H', "int");
    args_add(result, "lenquantiles", 'Q', "int");
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
    int _m, _M, _t, _n, _g, _N, _L, _a, _e, _h, _Q, nbins;
    int Nxx, Lxx;
    size_t A,C,G,T,N,tmp1,tmp2,tmp3,tmp4,tmp5;
    int64_t* all_lengths;
    size_t nalloclengths, curseqid;

    contig_s** inseqs;
    size_t minlen;
    size_t maxlen;
    size_t totlen, totlensum;
    size_t nseq, nseqtot;
    size_t i, j;
    size_t curlen, binsize;
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
    _a = args_ispresent(args, "avglen"); if (_a)defaultprint = 0;
    _e = args_ispresent(args, "medlen"); if (_e)defaultprint = 0;
    _t = args_ispresent(args, "totlen"); if (_t)defaultprint = 0;
    _n = args_ispresent(args, "nseq"); if (_n)defaultprint = 0;
    _g = args_ispresent(args, "GC"); if (_g)defaultprint = 0;
    _N = args_ispresent(args, "Nxx"); if (_N)defaultprint = 0;
    _L = args_ispresent(args, "Lxx"); if (_L)defaultprint = 0;
    _h = args_ispresent(args, "lenhistogram"); if (_h)defaultprint = 0;
    nbins = args_getint(args, "lenhistogram", 0, 20);
    _Q = args_ispresent(args, "lenquantiles"); if (_Q)defaultprint = 0;
    Nxx = args_getint(args,"Nxx",0,50);
    if(Nxx>99){
       args_report_warning(NULL, "N%d is not valid and has been transformed into to N99\n",Nxx);
       Nxx=99;
    }
    if(Nxx<1){
       args_report_warning(NULL, "N%d is not valid and has been transformed into to N1\n",Nxx);
       Nxx=1;
    }
    Lxx = args_getint(args,"Lxx",0,50);
    if(Lxx>99){
       args_report_warning(NULL, "L%d is not valid and has been transformed into to L99\n",Lxx);
       Lxx=99;
    }
    if(Lxx<1){
       args_report_warning(NULL, "L%d is not valid and has been transformed into to L1\n",Lxx);
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
    inseqs = NULL;
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
        minlen = 0;
        maxlen = 0;
        totlensum = 0;
        do {
            nseqtot += nseq;
            args_report_progress(NULL, _LLD_ " sequences read\r", nseqtot);
            curlen = contig_length(inseqs[0]);
            if(minlen==0){
                minlen = curlen;
                maxlen = minlen;
            }
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
            totlensum += totlen;
            if(inseqs)free(inseqs);
            inseqs = NULL;
            if(resume_at>0)
                inseqs = contigs_from_fastxPF_limited(f_in, &nseq,0xC0000000, &resume_at);
        } while(resume_at>0);
        if (inseqs) {
            nseqtot += nseq;
            args_report_progress(NULL, _LLD_ " sequences read\r", nseqtot);
            curlen = contig_length(inseqs[0]);
            if (minlen == 0) {
                minlen = curlen;
                maxlen = minlen;
            }
            A = C = G = T = N = 0;
            if (_g || defaultprint) {
                for (i = 0;i < nseq;i++) {
                    contig_countACGTN(inseqs[i], &tmp1, &tmp2, &tmp3, &tmp4, &tmp5);
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
                if (curseqid == nalloclengths) {
                    if (nalloclengths < 0x10000000) nalloclengths *= 2;
                    else nalloclengths += 0x10000000;
                    all_lengths = realloc(all_lengths, sizeof(int64_t)*nalloclengths);
                }
            }
            totlensum += totlen;
            if (inseqs)free(inseqs);
        }
        inseqs = NULL;
        if(nseqtot != curseqid){
            args_report_warning(NULL,"Sequences were not properly counted, somehow...");
        }
        args_report_progress(NULL, _LLD_ "\nFile read completely\n", nseqtot);
        if (defaultprint) {
            PFprintf(f_out, "Shortest sequence:\t" _LLD_ "\n", (int64_t)minlen);
            PFprintf(f_out, "Longest sequence:\t" _LLD_ "\n", (int64_t)maxlen);
            PFprintf(f_out, "Total length:\t" _LLD_ "\n", (int64_t)totlensum);
            PFprintf(f_out, "Number of sequences:\t" _LLD_ "\n", (int64_t)nseqtot);
            PFprintf(f_out, "GC content:\t%.4f\n", (double)(G+C)/(double)(A+G+C+T+N) );
            vec_sorti64(all_lengths,nseqtot);
            curseqid = nseqtot;
            curlen = 0;
            while(curlen < (totlensum * Nxx)/100 && curseqid > 0){
              curseqid--;
              curlen += all_lengths[curseqid];
            }
            PFprintf(f_out, "N" _LLD_ ":\t" _LLD_ "\n", Nxx, all_lengths[curseqid] );
            curseqid = nseqtot;
            curlen = 0;
            while(curlen < (totlensum * Lxx)/100 && curseqid > 0){
              curseqid--;
              curlen += all_lengths[curseqid];
            }
            PFprintf(f_out, "L" _LLD_ ":\t" _LLD_ "\n", Lxx, (int64_t)(nseqtot-curseqid) );
        }
        else {
            i = 0;
            if (_m) { PFprintf(f_out, _LLD_, minlen);i = 1; }
            if (_M) { if (i)PFputc('\t', f_out); PFprintf(f_out, _LLD_, (int64_t)maxlen);i = 1; }
            if (_t) { if (i)PFputc('\t', f_out); PFprintf(f_out, _LLD_, (int64_t)totlensum);i = 1; }
            if (_n) { if (i)PFputc('\t', f_out); PFprintf(f_out, _LLD_, (int64_t)nseqtot);i = 1; }
            if (_g) { if (i)PFputc('\t', f_out); PFprintf(f_out, "%.4f", (double)(G+C)/(double)(A+G+C+T+N) );i = 1; }
            if (_a) { if (i)PFputc('\t', f_out); PFprintf(f_out, _LLD_, veci64_avg(all_lengths, nseqtot)); i = 1; }
            if (_N || _L || _e || _h || _Q) {vec_sorti64(all_lengths,nseqtot);}
            if (_e) { if (i)PFputc('\t', f_out); PFprintf(f_out, _LLD_, (all_lengths[nseqtot / 2] + all_lengths[(nseqtot + 1) / 2]) / 2); i = 1; }
            if (_N) {
                curseqid = nseqtot;
                curlen = 0;
                while(curlen < (totlensum * Nxx)/100 && curseqid > 0){
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
                while(curlen < (totlensum * Lxx)/100 && curseqid > 0){
                   curseqid--;
                   curlen += all_lengths[curseqid];
                }
                if (i)PFputc('\t', f_out);
                PFprintf(f_out, _LLD_ , (int64_t)(nseqtot-curseqid) );
                i = 1;
            }
            if (_h) {
                if (i) PFputc('\n', f_out);
                if (nseqtot > nbins) {
                    minlen = all_lengths[(nseqtot) / nbins];
                    maxlen = all_lengths[(nseqtot * (nbins-1)) / nbins];
                    if (maxlen - minlen < nbins) {
                        nbins = maxlen - minlen + 1;
                        minlen = all_lengths[(nseqtot) / nbins];
                        maxlen = all_lengths[(nseqtot * (nbins - 1)) / nbins];
                    }
                    binsize = (maxlen - minlen + 1) / nbins;
                    curlen = minlen;
                    j = 0;
                    for (i = 0; i < nbins-1; i++) {
                        curlen += binsize;
                        nseq = 0;
                        while (all_lengths[j] <= curlen && j<nseqtot) {
                            nseq++;
                            j++;
                        }
                        PFprintf(f_out, "<=" _LLD_ "\t" _LLD_ "\n", curlen, nseq);
                    }
                    PFprintf(f_out, " >" _LLD_ "\t" _LLD_ "\n", curlen, nseqtot-j);
                }
                else {
                    binsize = 0;
                    PFprintf(f_out, _LLD_ , all_lengths[0]);
                    for (j = 1;j < nseqtot;j++) {
                        PFprintf(f_out, "\t" _LLD_, all_lengths[j]);
                    }
                    PFprintf(f_out, "\n");
                }
                i = 0;
            }
            if (_Q) {
                if (i) PFputc('\n', f_out);
                PFprintf(f_out, "P1\t" _LLD_, all_lengths[nseqtot / 100]);
                PFprintf(f_out, "\nP5\t" _LLD_, all_lengths[nseqtot / 20]);
                PFprintf(f_out, "\nQ1\t" _LLD_, all_lengths[nseqtot / 4]);
                PFprintf(f_out, "\nQ3\t" _LLD_, all_lengths[(nseqtot * 3) / 4]);
                PFprintf(f_out, "\nP95\t" _LLD_, all_lengths[(nseqtot*19) / 20]);
                PFprintf(f_out, "\nP99\t" _LLD_, all_lengths[(nseqtot*99) / 100]);
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


args_t* nucops_seqsummary_init_args(int argc, char** argv) {
    args_t* result;
    result = args_alloc();
    args_add(result, "input", 'i', "str");
    args_add(result, "output", 'o', "str");
    args_add(result, "length", 'l', "");
    args_add(result, "order", 'x', "");
    args_add(result, "name", 'n', "");
    args_add(result, "ACGT", 'c', "");
    args_add(result, "quality", 'q', "");
    args_add(result, "minquality", 'm', "");
    args_add(result, "maxquality", 'M', "");
    args_add_help(result, NULL, "Positional argument", "A single positional value is accepted, and should be the file name.", "treated as -i");
    args_parse(result, argc, argv);
    return result;
}

int nucops_seqsummary_a(args_t* args) {
    char* inputfn;
    char* outputfn;
    PF_t* f_in;
    PF_t* f_out;
    
    int defaultprint;
    int _m, _M, _q, _n, _c, _l, _x;
    size_t A, C, G, T, N;
    size_t seqid;
    contig_s* curcontig;
    int failure_code;
    int seemslike_fastq;
    int firstcol;
    
    inputfn = args_getstr(args, "input", 0, args_getstr(args, NULL, 0, "stdin"));
    if (strcmp(inputfn, "stdin") == 0) {
        args_report_error(NULL, "Using stdin as input is not supported\n");
        return 3;
    }
    outputfn = args_getstr(args, "output", 0, "stdout");
    defaultprint = 1;
    _m = args_ispresent(args, "minquality"); if (_m)defaultprint = 0;
    _M = args_ispresent(args, "maxquality"); if (_M)defaultprint = 0;
    _q = args_ispresent(args, "quality"); if (_q)defaultprint = 0;
    _x = args_ispresent(args, "order"); if (_x)defaultprint = 0;
    _n = args_ispresent(args, "name"); if (_n)defaultprint = 0;
    _c = args_ispresent(args, "ACGT"); if (_c)defaultprint = 0;
    _l = args_ispresent(args, "length"); if (_l)defaultprint = 0;

    if (defaultprint) {
        _l = 1;
        _c = 1;
        _n = 1;
        _x = 1;
    }

    args_report_info(NULL, "Parsing of arguments finished.\n");

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
    if (!failure_code) {
        seemslike_fastq = 0;
        args_report_info(NULL, "Input and output files successfully opened.\n");
        curcontig = contig_from_fastxPF(f_in);
        if (contig_qscores(curcontig))seemslike_fastq = 1;
        if (defaultprint && seemslike_fastq) {
            _m = 1;
            _M = 1;
            _q = 1;
        }
        firstcol = 1;
        seqid = 0;
        if (_x) {
            if (firstcol) {
                firstcol = 0;
                PFprintf(f_out, "order");
            }
            else {
                PFprintf(f_out, "\torder");
            }
        }
        if (_n) {
            if (firstcol) {
                firstcol = 0;
                PFprintf(f_out, "name");
            }
            else {
                PFprintf(f_out, "\tname");
            }
        }
        if (_l) {
            if (firstcol) {
                firstcol = 0;
                PFprintf(f_out, "length");
            }
            else {
                PFprintf(f_out, "\tlength");
            }
        }
        if (_c) {
            if (firstcol) {
                firstcol = 0;
                PFprintf(f_out, "A\tC\tG\tT\tN");
            }
            else {
                PFprintf(f_out, "\tA\tC\tG\tT\tN");
            }
        }
        if (_m) {
            if (firstcol) {
                firstcol = 0;
                PFprintf(f_out, "min_Quality");
            }
            else {
                PFprintf(f_out, "\tmin_Quality");
            }
        }
        if (_q) {
            if (firstcol) {
                firstcol = 0;
                PFprintf(f_out, "avg_Quality");
            }
            else {
                PFprintf(f_out, "\tavg_Quality");
            }
        }
        if (_M) {
            if (firstcol) {
                firstcol = 0;
                PFprintf(f_out, "max_Quality");
            }
            else {
                PFprintf(f_out, "\tmax_Quality");
            }
        }
        PFprintf(f_out, "\n");
        while (curcontig) {
            firstcol = 1;
            if (_x) {
                if (firstcol) {
                    firstcol = 0;
                    PFprintf(f_out, _LLD_, (long long)seqid);
                }
                else {
                    PFprintf(f_out, "\t" _LLD_ , (long long)seqid);
                }
            }
            if (_n) {
                if (firstcol) {
                    firstcol = 0;
                    PFprintf(f_out, "%s", contig_nameptr(curcontig));
                }
                else {
                    PFprintf(f_out, "\t%s", contig_nameptr(curcontig));
                }
            }
            if (_l) {
                if (firstcol) {
                    firstcol = 0;
                    PFprintf(f_out, _LLD_ , (long long)contig_length(curcontig));
                }
                else {
                    PFprintf(f_out, "\t" _LLD_, (long long)contig_length(curcontig));
                }
            }
            if (_c) {
                contig_countACGTN(curcontig, &A, &C, &G, &T, &N);
                if (firstcol) {
                    firstcol = 0;
                    PFprintf(f_out, _LLD_ "\t" _LLD_ "," _LLD_ "\t" _LLD_ "\t" _LLD_ , (long long)A, (long long)C, (long long)G, (long long)T, (long long)N);
                }
                else {
                    PFprintf(f_out, "\t" _LLD_ "\t" _LLD_ "\t" _LLD_ "\t" _LLD_ "\t" _LLD_, (long long)A, (long long)C, (long long)G, (long long)T, (long long)N);
                }
            }
            if (_m) {
                if (firstcol) {
                    firstcol = 0;
                    PFprintf(f_out, "%d", contig_min_qscore(curcontig));
                }
                else {
                    PFprintf(f_out, "\t%d", contig_min_qscore(curcontig));
                }
            }
            if (_q) {
                if (firstcol) {
                    firstcol = 0;
                    PFprintf(f_out, "%f", contig_avg_qscore(curcontig));
                }
                else {
                    PFprintf(f_out, "\t%f", contig_avg_qscore(curcontig));
                }
            }
            if (_M) {
                if (firstcol) {
                    firstcol = 0;
                    PFprintf(f_out, "%d", contig_max_qscore(curcontig));
                }
                else {
                    PFprintf(f_out, "\t%d", contig_max_qscore(curcontig));
                }
            }
            PFprintf(f_out, "\n");
            seqid++;
            contig_free(curcontig);
            curcontig = contig_from_fastxPF(f_in);
            if (seqid>10000 && seqid % 10000 == 1) {
                args_report_progress(NULL, "\r" _LLD_ "k sequences parsed", (long long)(seqid/1000));
            }
        }
        if (seqid > 10000) {
            args_report_progress(NULL, "\n" _LLD_ " sequences parsed\n", (long long)seqid);
        }
        else {
            args_report_progress(NULL, _LLD_ " sequences parsed\n", (long long)seqid);
        }
    }
    PFclose(f_in);
    PFclose(f_out);
    return failure_code;
}

int nucops_seqsummary(int argc, char** argv) {
    int result;
    args_t* args;

    args = nucops_seqsummary_init_args(argc, argv);
    result = 0;
    if (!args_ispresent(args, "help")) {
        result = nucops_seqsummary_a(args);
        if (result != 0) {
            args_report_error(args, "Seqsummary failed with code <%d>\n", result);
        }
    }
    args_free(args);
    return result;
}
