#include "argparser.h"
#include "sequence_base.h"

args_t* nucops_concat_init_args(int argc, char** argv) {
    args_t* result;
    result = args_alloc();

    args_add(result, "input", 'i', "str");
    args_add(result, "output", 'o', "str");
    args_parse(result, argc, argv);
    return result;
}

int nucops_concat_a(args_t* args) {
    contig_s** contigs;
    contig_s* result;
    char* inputfn;
    char* outputfn;
    size_t i, ncontigs;

    inputfn = args_getstr(args, "input", 0, "stdin");
    outputfn = args_getstr(args, "output", 0, "stdout");
    
    contigs = contigs_from_fasta(inputfn, &ncontigs);
    result = contigs_concatenate(contigs, ncontigs);
    contig_rename(result, inputfn);
    contigs_to_fasta(outputfn, &result, 1);
    contig_free(result);
    for (i = 0;i < ncontigs;i++) {
        contig_free(contigs[i]);
    }
    free(contigs);
    return 0;
}

int nucops_concat(int argc, char** argv) {
    int result;
    args_t* args;

    args = nucops_concat_init_args(argc, argv);
    result = 0;
    if (!args_ispresent(args, "help")) {
        result = nucops_concat_a(args);
        if (result != 0) {
            args_report_error(args, "Concat failed with code <%d>\n", result);
        }
    }
    args_free(args);
    return result;
}