#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "argparser.h"
#include "genosig.h"
#include "vecops.h"
#include "osportstd.h"
#include "datamap.h"


args_t* cmd_generate_args(int argc, char** argv) {
    args_t* result;
    result = args_alloc();
    args_add(result, "input", 'i', "str");
    args_add_help(result, "input", "INPUT FILE", "Specify the file based on which the signature should be generated", "input file (default: stdin)");
    args_add(result, "type", 't', "str,int");
    args_add_help(result, "type", "SIGNATURE TYPE",
        "Specify the signature type to be generated (fasta files only)."
        "Currently, the following signature types are available:\n"
        "  [freq] - frquency profile (2nd arg: kmer size)\n"
        "  [karl] - karlin signature (2nd arg: kmer size)\n"
        "  [logkarl] - base-10 logarithm of the kalrin signature\n"
        "  [mmz] - markov model z-scores (2nd arg: kmer size)\n"
        "  [GCskew] - GC skew\n",
        "signature type");
    args_add(result, "canonical", 'c', "");
    args_add_help(result, "canonical", "CANONICAL K-MERS", "When this argument is present, k-mers which  arethe reverse complement of "
        "\"smaller\" k-mers are removed", "keep only canonical k-mer (default: no)");
    args_add(result, "output", 'o', "str");
    args_add_help(result, "output", "OUTPUT FILE", "Specify the output file name", "ouput file (default: stdout)");
    args_add(result, "readable", 'r', "");
    args_add_help(result, "readable", "HUMAN-READABLE OUTPUT", "Generate human-readbale output instead of binary output", "generate human-readable output");
    args_add(result, "minlen", 'm', "int");
    args_add_help(result, "minlen", "MINIMUM CONTIG LENGTH",
        "When using fasta-formatted files as input, ignores contigs shorter than the specified length",
        "minimum contig length (default: 2000)");
    args_parse(result, argc, argv);
    return result;
}

void _exporttxtsig(genosig_t** sig, PF_t* output) {
    genosig_t** sigs;
    size_t i, count;
    if (genosig_info_ismuxed(*sig)) {
        sigs = genosig_demux(*sig, &count);
        for (i = 0;i<count;i++) {
            _exporttxtsig(sigs + i, output);
        }
        *sig = genosig_multiplex(sigs, count);
    }
    else {
        genosig_savetxt(*sig, output);
    }
}

uint64_t _revcompindex(uint64_t kmerindex, int kmerlen) {
    int i;
    uint64_t result;
    result = 0;
    for (i = 0;i < kmerlen;i++) {
        result <<= 2;
        result += 3 - ((kmerindex >> (2 * i)) % 4);
    }
    return result;
}

void _write_kmerheader(args_t* args, PF_t* output) {
    uint64_t i;
    uint64_t maxi;
    int kmerlen;
    char* cursequence;
    int iscanonical;
    int j;

    kmerlen = args_getint(args, "type", 0, 4);
    maxi = ((uint64_t)1) << (2 * kmerlen);

    iscanonical = args_ispresent(args, "canonical");
    cursequence = (char*)malloc(kmerlen + 1);
    cursequence[kmerlen] = 0;
    for (i = 0;i < maxi;i++) {
        if (!iscanonical || i <= _revcompindex(i, kmerlen)) {
            for (j = 0;j < kmerlen;j++) {
                switch ((i >> (2 * (kmerlen - j - 1))) % 4) {
                case 0:
                    cursequence[j] = 'A'; break;
                case 1:
                    cursequence[j] = 'C'; break;
                case 2:
                    cursequence[j] = 'G'; break;
                case 3:
                    cursequence[j] = 'T'; break;
                default:
                    cursequence[j] = 'N'; break;
                }
            }
            PFprintf(output, ",%s", cursequence);
        }
    }
}

void _exportsig(args_t* args, genosig_t** psig) {
    PF_t* output;
    if (args_ispresent(args, "readable")) {
        if (PFopen(&output, args_getstr(args, "output", 0, "stdout"), "wb") != 0) {
            args_report_warning(NULL, "Could not open file : <%s>\n", args_getstr(args, "output", 0, "stdout"));
            args_report_warning(NULL, "Output will be printed to stdout\n");
            PFclose(output);
            PFopen(&output, "stdout", "wb");
        }
        if (strcmp(args_getstr(args,"type",0,""),"GC")==0 || strcmp(args_getstr(args,"type",0,""),"len")==0)
            PFprintf(output,",value");
        else if (strcmp(args_getstr(args,"type",0,""),"GCskew")==0)
            PFprintf(output,",value,sd");
        else
            _write_kmerheader(args, output);
        PFprintf(output, "\n");
        _exporttxtsig(psig, output);
        PFclose(output);
    }
    else genosig_export(*psig, args_getstr(args, "output", 0, "stdout"));
}

int cmd_generate(int argc, char** argv) {
    args_t* args;
    PF_t* input;
    genosig_t* sig;
    genosig_t* nsig;
    genosig_t** allsigs;
    nucseq** fastafile;
    char* tmpstr;
    char* input_file_name;
    size_t ifn_len;
    size_t nseqs, i, j;
    size_t minseqlen;
    char** names;
    size_t nnames;
    size_t tmpsz;
    size_t badcount;
    args = cmd_generate_args(argc, argv);
    if (args_is_helpmode(args)) {
        args_free(args);
        return 0;
    }
    input = NULL;
    input_file_name = args_getstr(args, "input", 0, "stdin");
    nnames = 0;
    names = NULL;
    sig = genosig_importextra(input_file_name, GENOSIG_IMPORT_NOFASTA);
    args_report_info(NULL, "Bork\n");
    if (PFopen(&input, input_file_name, "rb") != 0) {
        args_report_error(NULL, "Could not open file : <%s>\n", input_file_name);
    }
    else if (!sig) {
        /* if that fail, try to read it as a fasta file */
        minseqlen = args_getint(args, "minlen", 0, 2000);
        fastafile = nucseq_array_from_fasta(input, &nseqs, 1, minseqlen, &badcount);
        if (!fastafile) {
            args_report_error(NULL, "Could not parse file : <%s>\n", PFgetbasenameptr(input));
            sig = NULL;
        }
        else {
            args_report_info(NULL, "File has %d sequences of sufficient (> %d) length\n", (int)nseqs, (int)minseqlen);
            tmpstr = args_getstr(args, "type", 0, "freq");
            args_report_info(NULL, "Selected signature type: \"%s\"\n", tmpstr );
            if (args_checkforstr(args, NULL, "global")) {
                nnames = 1;
                names = malloc(sizeof(char*) * nnames);
                ifn_len = strlen(input_file_name);
                names[0] = (char*)malloc(ifn_len + 1);
                memcpy(names[0], input_file_name, ifn_len);
                names[0][ifn_len] = 0;
                for (j = 0;j < ifn_len;j++) {
                    if (names[0][j] == '.') names[0][j] = '_';
                    if (names[0][j] == ' ') names[0][j] = '_';
                    if (names[0][j] == ',') names[0][j] = '_';
                    if (names[0][j] == '\t') names[0][j] = '_';
                }

                sig = genosig_fullgenome(fastafile, nseqs, 0, COPYLVL_INTEGRATE);
                
                if (strcmp(tmpstr, "freq") == 0) {
                    args_report_info(NULL, "Computing frequency signature\n");
                    genosig_kmerfreq(sig, args_getint(args, "type", 0, 4));
                }
                else if (strcmp(tmpstr, "karl") == 0) {
                    args_report_info(NULL, "Computing karlin signature of length :%d \n");
                    genosig_karlinsig(sig, args_getint(args, "type", 0, 4));
                }
                else if (strcmp(tmpstr, "logkarl") == 0) {
                    args_report_info(NULL, "Computing Base-10 logarithmic karlin signature of length :%d \n");
                    genosig_logkarlinsig(sig, args_getint(args, "type", 0, 4));
                }
                else if (strcmp(tmpstr, "mmz") == 0) {
                    args_report_info(NULL, "Computing markov model z-score signature\n");
                    genosig_mmz(sig, args_getint(args, "type", 0, 4));
                }
                else if (strcmp(tmpstr, "GC") == 0) {
                    args_report_info(NULL, "Computing GC content\n");
                    genosig_GC(sig,0);
                }
                else if (strcmp(tmpstr, "GCskew") == 0) {
                    args_report_info(NULL, "Computing GC skew\n");
                    genosig_GCskew(sig,0);
                }

                sig = genosig_linkname(sig, names[0]);
            }
            else {
                args_report_info(NULL, "Reading File as multiple (%d) signatures\n", (int)nseqs);
                allsigs = malloc(sizeof(genosig_t*)*nseqs);
                nnames = nseqs;
                names = (char**)malloc(sizeof(char*) * nnames);
                ifn_len = strlen(input_file_name);
                for (i = 0;i < nseqs;i++) {
                    allsigs[i] = genosig_fullgenome(fastafile + i, 1, 0, COPYLVL_FULL);

                    tmpsz = ifn_len + 1 + strlen(fastafile[i]->name);
                    names[i] = (char*)malloc(tmpsz + 1);
                    memcpy(names[i], input_file_name, ifn_len);
                    names[i][ifn_len] = 0;
                    memcpy(names[i] + ifn_len + 1, fastafile[i]->name, tmpsz - ifn_len - 1);
                    names[i][tmpsz] = 0;
                    for (j = 0;j < tmpsz;j++) {
                        if (names[i][j] == '.') names[i][j] = '_';
                        if (names[i][j] == ' ') names[i][j] = '_';
                        if (names[i][j] == ',') names[i][j] = '_';
                        if (names[i][j] == '\t') names[i][j] = '_';
                    }
                    names[i][ifn_len] = '.';

                    free(fastafile[i]);
                    if (strcmp(tmpstr, "freq") == 0) {
                        args_report_info(NULL, "Computing frequency signature %d\n", i);
                        genosig_kmerfreq(allsigs[i], args_getint(args, "type", 0, 4));
                    }
                    else if (strcmp(tmpstr, "karl") == 0) {
                        args_report_info(NULL, "Computing karlin signature %d\n", i);
                        genosig_karlinsig(allsigs[i], args_getint(args, "type", 0, 4));
                    }
                    else if (strcmp(tmpstr, "logkarl") == 0) {
                        args_report_info(NULL, "Computing Base-10 logarithmic karlin signature %d\n", i);
                        genosig_logkarlinsig(allsigs[i], args_getint(args, "type", 0, 4));
                    }
                    else if (strcmp(tmpstr, "mmz") == 0) {
                        args_report_info(NULL, "Computing markov model z-scoresignature %d\n", i);
                        genosig_mmz(allsigs[i], args_getint(args, "type", 0, 4));
                    }
                    else if (strcmp(tmpstr, "GC") == 0) {
                        args_report_info(NULL, "Computing GC content\n");
                        genosig_GC(allsigs[i],0);
                    }
                    else if (strcmp(tmpstr, "GCskew") == 0) {
                        args_report_info(NULL, "Computing GC skew\n");
                        genosig_GCskew(allsigs[i],0);
                    }
                    allsigs[i] = genosig_linkname(allsigs[i], names[i]);
                    args_report_info(NULL, "Signature name set to: %s\n", names[i]);
                }
                free(fastafile);
                sig = genosig_multiplex(allsigs, nseqs);
            }
            args_report_info(NULL, "Signatures successfully generated\n", (int)nseqs);
        }
    }
    if (sig) {
        if (args_ispresent(args, "canonical")) {
            sig = genosig_keepnonredundant(sig, args_getint(args, "type", 0, 4));
        }

        if (args_checkforstr(args, NULL, "average")) {
            args_report_info(NULL, "Computing average signature...\n");
            nsig = genosig_avgsig(sig, 0);
            if (nsig != sig) genosig_free(sig);
            if (nsig) {
                args_report_info(NULL, "Writing output...\n");
                _exportsig(args, &nsig);
                genosig_free(nsig);
            }
        }
        if (args_checkforstr(args, NULL, "global")) {
            args_report_info(NULL, "Returning computed signature without modifications...\n");
            _exportsig(args, &sig);
            genosig_free(sig);
        }
        else if (args_checkforstr(args, NULL, "contigs")) {
            args_report_info(NULL, "Returning computed signature without modifications...\n");
            _exportsig(args, &sig);
            genosig_free(sig);
        }
        if (args_checkforstr(args, NULL, "normdist")) {
            args_report_info(NULL, "Computing mean and standard deviation...\n");
            nsig = genosig_normparamsig(sig, 0);
            if (nsig != sig) genosig_free(sig);
            if (nsig) {
                args_report_info(NULL, "Writing output...\n");
                _exportsig(args, &nsig);
                genosig_free(nsig);
            }
        }
    }

    if (names) {
        for (i = 0;i < nnames;i++) {
            free(names[i]);
        }
        free(names);
    }

    PFclose(input);
    args_free(args);
    return 0;
}

int cmd_rename(int argc, char** argv) {
    genosig_t* sig;
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <filename> <new signature name>\n", argv[0]);
    }
    else {
        sig = genosig_importextra(argv[1], GENOSIG_IMPORT_NOFASTA);
        if (sig) {
            genosig_linkname(sig, argv[2]);
            genosig_export(sig, argv[1]);
        }
    }
    return 0;
}

int cmd_merge(int argc, char** argv) {
    genosig_t* sig;
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <filename> <new signature name>\n", argv[0]);
    }
    else {
        sig = genosig_importextra(argv[1], GENOSIG_IMPORT_NOFASTA);
        if (sig) {
            genosig_linkname(sig, argv[2]);
            genosig_export(sig, argv[1]);
        }
    }
    return 0;
}

int main(int argc, char** argv) {
    int nf;
    int32_t cmdid;
    DM32_t* cmds;

    if (argc < 2) {
        fprintf(stderr, "Usage: %s <command> [options]\n", argv[0]);
        return 1;
    }
    cmds = new_DM32(0, 0);
    DM32_append(cmds, "generate", (int)strlen("generate"), 1);
    DM32_append(cmds, "rename", (int)strlen("rename"), 2);
    DM32_sort(cmds);
    cmdid = DM32_get(cmds, argv[1], (int)strlen(argv[1]), &nf);
    if (!nf) {
        if (cmdid == 1) cmd_generate(argc - 1, argv + 1);
        else if (cmdid == 2) cmd_rename(argc - 1, argv + 1);
    }
    else {
        fprintf(stderr, "Unrecognized command: <%s>\n", argv[1]);
    }
    free_DM32(cmds);
    return 0;
}
