#include <stdlib.h>
#include <string.h>
#include "argparser.h"
#include "sequence_base.h"

args_t* nucops_select_init_args(int argc, char** argv) {
    args_t* result;
    result = args_alloc();
    
    args_add(result, "input", 'i', "str");
    args_add(result, "output", 'o', "str");
    args_add(result, "shorter_than", 's', "int");
    args_add(result, "longer_than", 'l', "int");
    args_add(result, "name_contains", 'n', "str");
    args_add(result, "ids_only", '@', "str");
    args_add(result, "avgq_above" , 'a', "int");
    args_add(result, "avgq_below", 'b', "int");
    args_add(result, "minq_above" , 'A', "int");
    args_add(result, "maxq_below", 'B', "int");
    args_add(result, "maxN", 'N', "int");
    args_add(result, "invert_criteria", '!', "");
    
    args_parse(result, argc, argv);
    return result;
}

int nucops_select_a(args_t* args) {
    contig_s** contigs;
    char* inputfn;
    char* outputfn;
    size_t shorter_than;
    size_t longer_than;
    char* name_contains;
    PF_t* f;
    size_t i, ncontigs, j;
    int avgq_above;
    int avgq_below;
    int minq_above;
    int maxq_below;
    int maxN;
    int inverse_select;
    long long resume_at;
    int passid;
    size_t A,C,G,T,N;

    inputfn = args_getstr(args, "input", 0, "stdin");
    outputfn = args_getstr(args, "output", 0, "stdout");
    longer_than = args_getint(args, "longer_than", 0, 0);
    shorter_than = args_getint(args, "shorter_than", 0, 0);
    name_contains = args_getstr(args, "name_contains", 0, NULL);
    inverse_select = (args_ispresent(args, "invert_criteria")?1:0);
    avgq_above = args_getint(args, "avgq_above", 0, -2);
    avgq_below = args_getint(args, "avgq_below", 0, -2);
    minq_above = args_getint(args, "minq_above", 0, -2);
    maxq_below = args_getint(args, "maxq_below", 0, -2);
    maxN = args_getint(args, "max?", 0, -1);
    args_report_info(NULL,"Arguments were parsed sucessfully\n");
    
    resume_at = 0;
    passid = 0;
    do {
        contigs = contigs_from_fastx_limited(inputfn, &ncontigs, 0xC0000000, &resume_at);
        args_report_progress(NULL,"Processing up to %dGb\n",(int)(resume_at/1000000000));
        if (longer_than > 0) {
            for (i = 0; i < ncontigs;i++) {
                if (contigs[i]){
                    if ( ((contig_length(contigs[i]) <= longer_than)?1:0) != inverse_select ) {
                        contig_free(contigs[i]);
                        contigs[i] = NULL;
                    }
                }
            }
            args_report_info(NULL,"Short sequences removed\n");
        }
        if (shorter_than > 0) {
            for (i = 0; i < ncontigs;i++) {
                if (contigs[i]){
                    if ( ((contig_length(contigs[i]) >= shorter_than)?1:0) != inverse_select) {
                        contig_free(contigs[i]);
                        contigs[i] = NULL;
                    }
                }
            }
            args_report_info(NULL,"Long sequences removed\n");
        }
        if ( avgq_above > -2) {
            for (i = 0; i < ncontigs;i++) {
                if (contigs[i]) {
                    if ( ((contig_avg_qscore(contigs[i]) > (float)avgq_above )?1:0) != inverse_select) {
                        contig_free(contigs[i]);
                        contigs[i] = NULL;
                    }
                }
            }
            args_report_info(NULL,"Sequences with low average quality removed\n");
        }
        if (maxN >= 0) {
            for (i = 0; i < ncontigs;i++) {
                if (contigs[i]) {
                    contig_countACGTN(contigs[i], &A,&C,&G,&T,&N);
                    if ( ((N>(size_t)maxN)?1:0) != inverse_select) {
                        contig_free(contigs[i]);
                        contigs[i] = NULL;
                    }
                }
            }
            args_report_info(NULL,"Sequences with at least %d N's removed\n", maxN+1);
        }

        if (avgq_below > -2) {
            for (i = 0; i < ncontigs;i++) {
                if (contigs[i]) {
                    if ( ((contig_avg_qscore(contigs[i]) < (float)avgq_below )?1:0) != inverse_select) {
                        contig_free(contigs[i]);
                        contigs[i] = NULL;
                    }
                }
            }
            args_report_info(NULL,"Sequences with high average quality removed\n");
        }
        if (minq_above > -2) {
            for (i = 0; i < ncontigs;i++) {
                if (contigs[i]) {
                    if ( ((contig_min_qscore(contigs[i]) > (float)minq_above )?1:0) != inverse_select) {
                        contig_free(contigs[i]);
                        contigs[i] = NULL;
                    }
                }
            }
            args_report_info(NULL,"Sequences with minimal quality below threshold removed\n");
        }
        if ( maxq_below > -2) {
            for (i = 0; i < ncontigs;i++) {
                if (contigs[i]) {
                    if ( ((contig_max_qscore(contigs[i]) < (float)maxq_below )?1:0) != inverse_select) {
                        contig_free(contigs[i]);
                        contigs[i] = NULL;
                    }
                }
            }
            args_report_info(NULL,"Sequences with maximal quality above threshold removed\n");
        }

        if (name_contains) {
            for (i = 0; i < ncontigs;i++) {
                if (contigs[i]) {
                    if ( ((strstr(contig_nameptr(contigs[i]),name_contains)==NULL)?1:0) != inverse_select) {
                        contig_free(contigs[i]);
                        contigs[i] = NULL;
                    }
                }
            }
            args_report_info(NULL,"Sequences with names not containing target string removed\n");
        }
        j = 0;
        for (i = 0;i < ncontigs;i++) {
            if (i != j) {
                contigs[j] = contigs[i];
            }
            if (contigs[i])j++;
        }
        args_report_info(NULL,"Reordering remaining contigs done.\n");
        if(passid = 0){
            if (PFopen(&f, outputfn, "wb") == 0) {
                if (args_ispresent(args, "ids_only")) {
                    for (i = 0;i < j;i++) {
                        PFprintf(f,"%s\n",contig_nameptr(contigs[i]));
                    }
                }
                else {
                    contigs_to_fastxPF(f, contigs, j);
                }
            }
        }
        else {
            if (PFopen(&f, outputfn, "ab") == 0) {
                if (args_ispresent(args, "ids_only")) {
                    for (i = 0;i < j;i++) {
                        PFprintf(f,"%s\n",contig_nameptr(contigs[i]));
                    }
                }
                else {
                    contigs_to_fastxPF(f, contigs, j);
                }
            }
        }
        passid++;
        PFclose(f);
        
        args_report_info(NULL,"Output written sucessfuly.\n");
        for (i = 0;i < j;i++) {
            contig_free(contigs[i]);
        }
        free(contigs);
    } while(resume_at != 0);
    return 0;
}

int nucops_select(int argc, char** argv) {
    int result;
    args_t* args;

    args = nucops_select_init_args(argc, argv);
    result = 0;
    if (!args_ispresent(args, "help")) {
        result = nucops_select_a(args);
        if (result != 0) {
            args_report_error(args, "Select failed with code <%d>\n", result);
        }
    }
    args_free(args);
    return result;
}
