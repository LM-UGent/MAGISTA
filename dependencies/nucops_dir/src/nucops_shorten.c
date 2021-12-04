#include "argparser.h"
#include "sequence_base.h"

args_t* nucops_shorten_init_args(int argc, char** argv) {
    args_t* result;
    result = args_alloc();

    args_add(result, "input", 'i', "str");
    args_add(result, "output", 'o', "str");
    args_add(result, "kmer", 'k', "int");
    args_add(result, "window", 'w', "int");
    args_add(result, "minimizers", 'm', "");
    args_add(result, "maximizers", 'M', "");
    args_parse(result, argc, argv);
    return result;
}

int nucops_shorten_a(args_t* args) {
    contig_s** contigs;
    contig_s* result;
    char* inputfn;
    char* outputfn;
    size_t i, ncontigs;
    size_t window, kmerlen;
    size_t* positions;
    size_t npositions;
    size_t curcontig;
    int mode;

    inputfn = args_getstr(args, "input", 0, "stdin");
    outputfn = args_getstr(args, "output", 0, "stdout");
    window = args_getint(args, "window", 0, 1000);
    kmerlen = args_getint(args, "kmer", 0, 15);
    mode = 0;
    if (args_ispresent(args, "minimizers"))mode = 1;
    if (args_ispresent(args, "maximizers"))mode = 2;
    if (mode == 0) {
        args_report_info(NULL, "Mode not set - using minimizers.");
        args_report_info(NULL, "You can set mode with the '-m' or '-M' option.");
        mode = 1;
    }

    contigs = contigs_from_fasta(inputfn, &ncontigs);
    curcontig = 0;
    for (i = 0;i < ncontigs;i++) {
        switch(mode) {
            case 0:
            case 1:
                positions = contig_minimizer_positions(contigs[i], window, kmerlen, &npositions);
                break;
            case 2:
                positions = contig_maximizer_positions(contigs[i], window, kmerlen, &npositions);
                break;
            default:
                positions = NULL;
                break;
        }
        result = NULL;
        if (positions) {
            result = contig_kmer_sequence_from_positions(contigs[i], positions, npositions, kmerlen);
            free(positions);
        }
        contig_free(contigs[i]);
        if (result) {
            contigs[curcontig] = result;
            curcontig++;
        }
    }
    contigs_to_fasta(outputfn, contigs, curcontig);
    free(contigs);
    return 0;
}

int nucops_shorten(int argc, char** argv) {
    int result;
    args_t* args;

    args = nucops_shorten_init_args(argc, argv);
    result = 0;
    if (!args_ispresent(args, "help")) {
        result = nucops_shorten_a(args);
        if (result != 0) {
            args_report_error(args, "Shorten failed with code <%d>\n", result);
        }
    }
    args_free(args);
    return result;
}