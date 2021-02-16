#include "osportstd.h"
#include "sequence_base.h"
#include "argparser.h"
#include "textparsing.h"
#include <math.h>
#include <string.h>

args_t* nucops_winsplit_init_args(int argc, char** argv) {
    args_t* result;
    result = args_alloc();
    args_add(result, "input", 'i', "str");
    args_add(result, "output", 'o', "str");
    args_add(result, "window", 'w', "int");
    args_add(result, "step", 's', "int");
    args_add(result, "progressive", 'p', "");
    args_parse(result, argc, argv);
    return result;
}

int nucops_winsplit_a(args_t* args) {
    char* inputfn;
    char* outputfn;
    char* newname;
    contig_s** original;
    contig_s** result;
    size_t i, j, ncontigs, nnew, nnew_extra, nalloc, namelen, idlen;
    size_t winsize;
    size_t step;
    size_t contiglen;
    /* initialize relevant parameters */
    inputfn = args_getstr(args, "input", 0, "stdin");
    outputfn = args_getstr(args, "output", 0, "stdout");
    winsize = (size_t)args_getint(args, "window", 0, 1000);
    step = (size_t)args_getint(args, "step", 0, 1000);

    original = contigs_from_fasta(inputfn, &ncontigs);
    if (!original) return 1;

    nnew = 0;
    nalloc = 16;
    result = (contig_s**)malloc(sizeof(contig_s*)*nalloc);
    newname = NULL;
    for (i = 0;i < ncontigs;i++) {
        contiglen = contig_length(original[i]);
        if (contiglen >= winsize) {
            namelen = strlen(contig_nameptr(original[i]));
            idlen = (size_t) log10((double)contiglen)+1;
            newname = (char*)realloc(newname, namelen+idlen+3);
            memcpy(newname, contig_nameptr(original[i]), namelen);
            newname[namelen] = '+';
            nnew_extra = (contiglen - winsize) / step + 1;
            if (nnew + nnew_extra > nalloc) {
                nalloc *= 2;
                if (nnew + nnew_extra > nalloc) nalloc = (nnew + nnew_extra)*2;
                result = (contig_s**)realloc(result, sizeof(contig_s*)*nalloc);
            }
            for (j = 0;j < nnew_extra;j++) {
#ifndef _WIN32
                sprintf(newname + namelen + 1, _LLD_, step*j);
#else
                sprintf_s(newname + namelen + 1, idlen+1, _LLD_, (int64_t)step*j);
#endif
                newname[namelen + idlen + 1] = 0;
                result[nnew + j] = contig_subsection(original[i], step*j, winsize);
                contig_rename(result[nnew + j], newname);
            }
            nnew += nnew_extra;
        }
        contig_free(original[i]);
    }
    if(newname) free(newname);
    free(original);
    contigs_to_fasta(outputfn, result, nnew);
    for(i=0;i<nnew;i++)
        contig_free(result[i]);
    free(result);
    return 0;
}

int nucops_winsplit_p(args_t* args) {
    /* not yet implemented */
    return -1;
}

int nucops_winsplit(int argc, char** argv) {
    args_t* args;
    int result;
    args = nucops_winsplit_init_args(argc, argv);

    result = 0;
    if (!args_ispresent(args, "help")) {
        if (args_ispresent(args, "progressive")) {
            result = nucops_winsplit_p(args);
        }
        else {
            result = nucops_winsplit_a(args);
        }
        if (result != 0) args_report_error(args, "winsplit failed with code <%d>\n", result);
    }
    args_free(args);
    return result;
}