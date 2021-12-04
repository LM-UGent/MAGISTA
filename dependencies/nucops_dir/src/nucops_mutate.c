#include "argparser.h"
#include "sequence_base.h"
#include "randgen.h"
#include "vecops.h"

args_t* nucops_mutate_init_args(int argc, char** argv) {
    args_t* result;
    result = args_alloc();

    args_add(result, "input", 'i', "str");
    args_add(result, "output", 'o', "str");
    args_add(result, "snps", 's', "int");
    args_add(result, "random_seed", 'r', "int");
    args_parse(result, argc, argv);
    return result;
}

void nucops_mutate_add_snps(contig_s** contigs,size_t ncontigs, size_t count, randgen_s* rng) {
    size_t i;
    size_t totlen;
    uint64_t* snppos;
    size_t contigid;
    size_t contigoffset;
    char possible_nucs[4] = { 'A','C','G','T' };
    char nuc, randnuc;
    totlen = 0;
    for (i = 0;i < ncontigs;i++) {
        totlen += contig_length(contigs[i]);
    }
    snppos = malloc(sizeof(size_t)*count);
    for (i = 0;i < count;i++) {
        snppos[i] = (uint64_t)randgen_uniform_i64(rng, 0, totlen - 1);
    }
    vec_sorti64((int64_t*)snppos, count);
    contigid = 0;
    contigoffset = 0;
    for (i = 0;i < count;i++) {
        while ((size_t)(snppos[i]) - contigoffset > contig_length(contigs[contigid])) {
            contigoffset += contig_length(contigs[contigid]);
            contigid++;
        }
        nuc = contig_getnuc(contigs[contigid], snppos[i]);
        randnuc = nuc;
        while(randnuc == nuc)
            randnuc = possible_nucs[randgen_uniform_i64(rng, 0, 3)];
        contig_setnuc(contigs[contigid], snppos[i], randnuc);
    }
    free(snppos);
}

int nucops_mutate_prog(args_t* args) {
    contig_s** contigs;
    char* inputfn;
    char* outputfn;
    size_t i, ncontigs;
    size_t snp_count;
    randgen_s* rng;

    rng = randgen_alloc();
    if (args_ispresent(args, "random_seed")) {
        randgen_seed(rng, args_getint(args, "random_seed", 0, 2056));
    }
    inputfn = args_getstr(args, "input", 0, "stdin");
    outputfn = args_getstr(args, "output", 0, "stdout");
    snp_count = (size_t)args_getint(args, "snps", 0, 1);
    contigs = contigs_from_fasta(inputfn, &ncontigs);
    if (args_ispresent(args, "snps")) {
        nucops_mutate_add_snps(contigs, ncontigs, snp_count, rng);
    }
    contigs_to_fasta(outputfn, contigs, ncontigs);
    for (i = 0;i < ncontigs;i++) {
        contig_free(contigs[i]);
    }
    free(contigs);
    randgen_free(rng);
    return 0;
}

int nucops_mutate(int argc, char** argv) {
    int result;
    args_t* args;

    args = nucops_mutate_init_args(argc, argv);
    result = 0;
    if (!args_ispresent(args, "help")) {
        result = nucops_mutate_prog(args);
        if (result != 0) {
            args_report_error(args, "Mutate failed with code <%d>\n", result);
        }
        args_free(args);
    }
    return result;
}