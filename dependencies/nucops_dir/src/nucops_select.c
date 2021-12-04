#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "randgen.h"
#include "argparser.h"
#include "sequence_base.h"
#include "vecops.h"
#include "datamap.h"
#include "textparsing.h"

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
    args_add(result, "keep_frac", 'k', "float,str");
    args_add_help(result, "keep_frac", "RANDOM FRACTION",
        "Used to specify the fractions of the output selected sequences that should be kept. Values between 0 and 1 indicate the fraction of sequences to be kept,"
        " whereas values above 1 indicate the total number of base pairs to be kept.\n The second argument should be 'bps'(default) or 'seqs' depending on whether "
        " the the desired result should be the fraction of the total length of the sequences or the total count of the sequences, respectively"
        "\nNote : '%','k','M','G', and 'T' are accepted suffixes. ", "fraction of the selected sequences to be kept");
    args_add(result, "ignore_suffix", 'p', "");
    args_add(result, "random_seed", 'r', "int");
    args_add(result, "sequence_list", 'e', "str");
    
    args_parse(result, argc, argv);
    return result;
}

int nucops_select_a(args_t* args) {
    contig_s** contigs;
    contig_s** hashed_contigs;
    randgen_s* rng;
    DM64_t* nameset;
    char* inputfn;
    char* outputfn;
    char* seqlistfn;
    char* line;
    char* tmpstr;
    size_t shorter_than;
    size_t longer_than;
    char* name_contains;
    PF_t* f;
    PF_t* names;
    size_t i, ncontigs, totalcontigs, prevtotalcontigs, j, k;
    int avgq_above;
    int avgq_below;
    int minq_above;
    int maxq_below;
    int maxN;
    int inverse_select;
    long long resume_at;
    long long totalseqs;
    double genomefrac;
    int fractype;
    int passid;
    int nf;
    int errflag;
    size_t A,C,G,T,N;
    uint64_t* shuffled_table;
    uint64_t tmpval, tmpmaxval, tmpcount;
    int64_t lastUID;
    int hasseqlist;
    int ignore_suffix;
    long long randomseed;
    size_t count;
    
    args_report_info(NULL, "nucops select is now runnning\n");

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
    maxN = args_getint(args, "maxN", 0, -1);
    ignore_suffix = args_ispresent(args, "ignore_suffix");

    randomseed = -1;
    rng = randgen_alloc();
    if (args_ispresent(args, "random_seed")) {
        randomseed = (long long)((unsigned int)args_getint(args, "random_seed", 0, 1));
        randgen_seed(rng, (uint64_t)randomseed);
    }
    else {
        randgen_seed(rng, (uint64_t)time(NULL));
    }
    genomefrac = args_getdouble(args, "keep_frac", 0, -1.);
    fractype = 0;
    if (strcmp(args_getstr(args,"keep_frac", 0, "bps"),"seqs")==0) {
        fractype = 1;
    }
    
    if(strcmp(inputfn,"stdin") == 0) {
        if(genomefrac > 0 && genomefrac != 1){
            args_report_error(NULL,"Random subsampling is not supported when input is stdin!\n");
            randgen_free(rng);
            return 1;
        }
        args_report_warning(NULL,"Please note that when stdin is used, the entire input should be able to fit within the computter's RAM\n");
    }
    hasseqlist = args_ispresent(args, "sequence_list");
    args_report_info(NULL, "Using sequence list: %s\n",hasseqlist?"yes":"no");
    seqlistfn = "";
    if (hasseqlist) {
        args_report_info(NULL, "Reads to keep are in:\n");
        seqlistfn = args_getstr(args, "sequence_list", 0, "stdin");
        args_report_info(NULL, "  %s\n", seqlistfn);
    }
    if (strcmp(inputfn, seqlistfn) == 0 || strcmp(outputfn, seqlistfn) == 0 || strcmp(inputfn, outputfn) == 0) {
        args_report_error(NULL, "Cannot use the same file (%s) in multiple arguments\n");
        randgen_free(rng);
        return 1;
    }
    
    args_report_info(NULL,"Arguments were parsed sucessfully\n");
    
    resume_at = 0;
    passid = 0;
    totalcontigs = 0;
    hashed_contigs = NULL;
    lastUID = 0;

    nameset = NULL;
    if (hasseqlist) {
        nameset = new_DM64(DM_ALGORITHM_BASICSORTEDLIST, 0);
        count = 0;
        if (!PFopen(&names, seqlistfn, "rb")) {
            while (line = PFreadline(names)) {
                DM64_append(nameset, line, (int)strlen(line), count);
                count++;
                free(line);
            }
        }
        PFclose(names);
    }

    totalseqs = 0;
    do {
        args_report_progress(NULL,"Reading input... (pass %d)\n",passid);
        if(strcmp(inputfn,"stdin") == 0){
            contigs = contigs_from_fastx(inputfn, &ncontigs);
            resume_at = 0;
        }
        else{
            contigs = contigs_from_fastx_limited(inputfn, &ncontigs, 0xC0000000, &resume_at);
        }
        totalseqs += ncontigs;
        args_report_progress(NULL,"Processing up to %dGb\n",(int)(resume_at/1000000000));
        for (i = 0;i < ncontigs;i++) {
            contig_setUID((contigs[i]),lastUID);
            contig_toggle_suffix(contigs[i], !ignore_suffix);
            lastUID++;
        }
        if (nameset) {
            for (i = 0; i < ncontigs;i++) {
                if (contigs[i]) {
                    (void)DM64_get(nameset, contig_nameptr(contigs[i]), (int)strlen(contig_nameptr(contigs[i])), &nf);
                    if (nf) {
                        k = 0;
                        tmpstr = contig_nameptr(contigs[i]);
                        while (tmpstr[k] && tmpstr[k] != ' ') {
                            k++;
                        }
                        (void)DM64_get(nameset, contig_nameptr(contigs[i]), (int)k, &nf);
                    }
                    if ((nf ? 1 : 0) != inverse_select) {
                        contig_free(contigs[i]);
                        contigs[i] = NULL;
                    }
                }
            }
            args_report_info(NULL, "Sequences with names not in the list file removed\n");
        }

        if (longer_than > 0) {
            count = 0;
            for (i = 0; i < ncontigs;i++) {
                if (contigs[i]){
                    if ( ((contig_length(contigs[i]) <= longer_than)?1:0) != inverse_select ) {
                        contig_free(contigs[i]);
                        contigs[i] = NULL;
                        count++;
                    }
                }
            }
            args_report_info(NULL, _LLD_ " short sequences removed\n", (long long)count );
        }
        if (shorter_than > 0) {
            count = 0;
            for (i = 0; i < ncontigs;i++) {
                if (contigs[i]){
                    if ( ((contig_length(contigs[i]) >= shorter_than)?1:0) != inverse_select) {
                        contig_free(contigs[i]);
                        contigs[i] = NULL;
                        count++;
                    }
                }
            }
            args_report_info(NULL, _LLD_ " long sequences removed\n", (long long)count);
        }
        if ( avgq_above > -2) {
            count = 0;
            for (i = 0; i < ncontigs;i++) {
                if (contigs[i]) {
                    if ( ((contig_avg_qscore(contigs[i]) < (float)avgq_above )?1:0) != inverse_select) {
                        contig_free(contigs[i]);
                        contigs[i] = NULL;
                        count++;
                    }
                }
            }
            args_report_info(NULL, _LLD_ " sequences with low average quality removed\n", (long long)count);
        }
        
        if (maxN >= 0) {
            count = 0;
            for (i = 0; i < ncontigs;i++) {
                if (contigs[i]) {
                    contig_countACGTN(contigs[i], &A,&C,&G,&T,&N);
                    if ( ((N>(size_t)maxN)?1:0) != inverse_select) {
                        contig_free(contigs[i]);
                        contigs[i] = NULL;
                        count++;
                    }
                }
            }
            args_report_info(NULL, _LLD_ " sequences with at least %d N's removed\n", (long long)count, maxN+1);
        }

        if (avgq_below > -2) {
            count = 0;
            for (i = 0; i < ncontigs;i++) {
                if (contigs[i]) {
                    if ( ((contig_avg_qscore(contigs[i]) > (float)avgq_below )?1:0) != inverse_select) {
                        contig_free(contigs[i]);
                        contigs[i] = NULL;
                        count++;
                    }
                }
            }
            args_report_info(NULL, _LLD_ " sequences with high average quality removed\n", (long long)count);
        }
        if (minq_above > -2) {
            count = 0;
            for (i = 0; i < ncontigs;i++) {
                if (contigs[i]) {
                    if ( ((contig_min_qscore(contigs[i]) < (float)minq_above )?1:0) != inverse_select) {
                        contig_free(contigs[i]);
                        contigs[i] = NULL;
                        count++;
                    }
                }
            }
            args_report_info(NULL, _LLD_ " sequences with minimal quality below threshold removed\n");
        }
        if ( maxq_below > -2) {
            count = 0;
            for (i = 0; i < ncontigs;i++) {
                if (contigs[i]) {
                    if ( ((contig_max_qscore(contigs[i]) > (float)maxq_below )?1:0) != inverse_select) {
                        contig_free(contigs[i]);
                        contigs[i] = NULL;
                        count++;
                    }
                }
            }
            args_report_info(NULL, _LLD_ " sequences with maximal quality above threshold removed\n", (long long)count);
        }

        if (name_contains) {
            count = 0;
            for (i = 0; i < ncontigs;i++) {
                if (contigs[i]) {
                    if ( ((strstr(contig_nameptr(contigs[i]),name_contains)==NULL)?1:0) != inverse_select) {
                        contig_free(contigs[i]);
                        contigs[i] = NULL;
                        count++;
                    }
                }
            }
            args_report_info(NULL, _LLD_ " sequences with names not containing target string removed\n", (long long)count);
        }
        j = 0;
        for (i = 0;i < ncontigs;i++) {
            if (i != j) {
                contigs[j] = contigs[i];
            }
            if (contigs[i])j++;
        }
        ncontigs = j;
        prevtotalcontigs = totalcontigs;
        totalcontigs += ncontigs;
        args_report_info(NULL,"Reordering remaining (" _LLD_ ") contigs done.\n", (long long)j);

        if (genomefrac > 0 && (genomefrac != 1 || fractype == 1)) {
            /* discard the nucleotide & qscore data to reduce memory usage (note: a hash is kept)*/
            hashed_contigs = realloc(hashed_contigs,sizeof(contig_s*)*(totalcontigs+ncontigs));
            errflag = 0;
            for (i = 0;i < ncontigs;i++) {
                contig_discard_data(contigs[i]);
                hashed_contigs[prevtotalcontigs + i] = contigs[i];
                if (hashed_contigs[prevtotalcontigs + i] == NULL) {
                    errflag = 1;
                }
            }
            args_report_info(NULL, "Contigs were hashed and unneeded sequence information was discarded.\n");
        }
        else {
            if(passid == 0){
                if (PFopen(&f, outputfn, "wb") == 0) {
                    if (args_ispresent(args, "ids_only")) {
                        args_report_info(args, "Only reporting ids\n");
                        for (i = 0;i < ncontigs;i++) {
                            PFprintf(f,"%s\n",contig_nameptr(contigs[i]));
                        }
                    }
                    else {
                        contigs_to_fastxPF(f, contigs, ncontigs);
                    }
                }
            }
            else {
                if (PFopen(&f, outputfn, "ab") == 0) {
                    if (args_ispresent(args, "ids_only")) {
                        args_report_info(args, "Only reporting ids\n");
                        for (i = 0;i < ncontigs;i++) {
                            PFprintf(f,"%s\n",contig_nameptr(contigs[i]));
                        }
                    }
                    else {
                        contigs_to_fastxPF(f, contigs, ncontigs);
                    }
                }
            }
            PFclose(f);
            args_report_info(NULL,"Output written sucessfuly.\n");
            for (i = 0;i < ncontigs;i++) {
                contig_free(contigs[i]);
            }
        }
        passid++;
        free(contigs);
    } while(resume_at != 0);
    if (totalcontigs == 0) {
        args_report_warning(NULL, "No contigs match selection criteria!");
    }
    if (genomefrac > 0 && (genomefrac != 1 || fractype == 1)) {
        args_report_info(NULL,"Final number of contigs: " _LLD_ "\n", totalcontigs);
        shuffled_table = malloc(sizeof(uint64_t)*totalcontigs);
        tmpmaxval = 0;
        errflag = 0;
        for (i = 0;i < totalcontigs;i++) {
            shuffled_table[i] = (uint64_t)i;
            if (hashed_contigs[i])
                tmpmaxval += contig_length(hashed_contigs[i]);
            else
                errflag = 1;
        }
        if (errflag) {
            args_report_error(NULL, "Some sequences were NULL. This can only be the results of a bug.\n");
            exit(1);
        }
        args_report_info(NULL,"Total length of qualifying contigs: " _LLD_ "\n", (uint64_t)tmpmaxval);
        if (randomseed != 0) {
            for (i = 0;i < totalcontigs;i++) {
                k = totalcontigs - i;
                j = (size_t)randgen_uniform_i64(rng,0,(int64_t)k);
                tmpval = shuffled_table[k - 1];
                shuffled_table[k - 1] = shuffled_table[j];
                shuffled_table[j] = tmpval;
            }
            args_report_info(NULL,"Shuffling contigs done\n");
        }
        else {
            args_report_info(NULL,"Contigs were not shuffled, since random seed was set to zero\n");
        }
        if (fractype == 1) {
            tmpmaxval = totalcontigs;
        }
        if (genomefrac < 1) {
            tmpmaxval = (uint64_t)(((double)tmpmaxval)*genomefrac);
        }
        else {
            tmpmaxval = (uint64_t)genomefrac;
        }
        tmpval = 0;
        i = 0;
        tmpcount = 0;
        if (fractype == 0) {
            args_report_info(NULL, "Target number of bases: " _LLD_ "\n", tmpmaxval);
        }
        else if (fractype == 1) {
            args_report_info(NULL, "Target number of sequences: " _LLD_ "\n", tmpmaxval);
        }
        while (tmpval < tmpmaxval && i < totalcontigs) {
            k = shuffled_table[i];
            if (fractype == 1) tmpval++;
            else tmpval += contig_length(hashed_contigs[k]);
            tmpcount += contig_length(hashed_contigs[k]);
            i++;
        }
        ncontigs = i;
        args_report_info(NULL, "Kept " _LLD_ " sequences totalling " _LLD_ " bases\n", ncontigs, tmpcount);
        while (i < totalcontigs) {
            k = shuffled_table[i];
            contig_free(hashed_contigs[k]);
            hashed_contigs[k] = NULL;
            i++;
        }
        shuffled_table = realloc(shuffled_table,sizeof(uint64_t)*ncontigs);
        j = 0;
        for (i = 0;i < totalcontigs;i++) {
            if (i != j) {
                hashed_contigs[j] = hashed_contigs[i];
            }
            if (hashed_contigs[i])j++;
        }
        totalcontigs = j;
        for (i = 0;i < totalcontigs;i++) {
            shuffled_table[i] = contig_UID(hashed_contigs[i]) ;
            contig_free(hashed_contigs[i]);
        }
        free(hashed_contigs);
        vec_sorti64((int64_t*)shuffled_table,totalcontigs);
        
        resume_at = 0;
        k = 0;
        j = 0;
        lastUID = 0;
        passid = 0;
        do {
            args_report_progress(NULL,"Reading input... (pass %d)\n",passid);
            contigs = contigs_from_fastx_limited(inputfn, &ncontigs, 0xC0000000, &resume_at);
            args_report_progress(NULL,"Re-processing up to %dGb\n",(int)(resume_at/1000000000));
            for (i = 0;i < ncontigs;i++) {
                contig_setUID((contigs[i]),lastUID);
                lastUID++;
            }
            for (i = 0;i < ncontigs;i++) {
                if(k == totalcontigs || contig_UID(contigs[i]) != shuffled_table[k]) {
                    contig_free(contigs[i]);
                    contigs[i] = NULL;
                }
                else{
                    k++;
                }
            }
            j = 0;
            for (i = 0;i < ncontigs;i++) {
                if (i != j) {
                    contigs[j] = contigs[i];
                }
                if (contigs[i])j++;
            }
            ncontigs = j;
            if(passid == 0){
                if (PFopen(&f, outputfn, "wb") == 0) {
                    if (args_ispresent(args, "ids_only")) {
                        args_report_info(args, "Only reporting ids\n");
                        for (i = 0;i < ncontigs;i++) {
                            PFprintf(f,"%s\n",contig_nameptr(contigs[i]));
                        }
                    }
                    else {
                        contigs_to_fastxPF(f, contigs, ncontigs);
                    }
                }
            }
            else {
                if (PFopen(&f, outputfn, "ab") == 0) {
                    if (args_ispresent(args, "ids_only")) {
                        args_report_info(args, "Only reporting ids\n");
                        for (i = 0;i < ncontigs;i++) {
                            PFprintf(f,"%s\n",contig_nameptr(contigs[i]));
                        }
                    }
                    else {
                        contigs_to_fastxPF(f, contigs, ncontigs);
                    }
                }
            }
            PFclose(f);
            passid++;
            for (i = 0;i < ncontigs;i++) {
                contig_free(contigs[i]);
            }
            free(contigs);
        } while(resume_at != 0);
    }
    randgen_free(rng);
    if (nameset)free_DM64(nameset);
    args_report_info(NULL,"Output written sucessfuly.\n");
    
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
            args_report_error(args, "splitfile failed with code <%d>\n", result);
        }
    }
    args_free(args);
    return result;
}


args_t* nucops_splitfile_init_args(int argc, char** argv) {
    args_t* result;
    result = args_alloc();

    args_add(result, "input", 'i', "str");
    args_add(result, "output", 'o', "str");
    args_add(result, "max_bases_per_file", 'b', "int");
    args_add(result, "max_sequences_per_file", 's', "int");
    args_add(result, "by_prefix", 'p', "str");
    args_add(result, "keep_names", 'n', "");
    args_add(result, "ids_only", '@', "str");
    args_add(result, "by_table", 'e', "str");
    
    args_add_help(result, "by_prefix", "SPLIT FILE BY SEQUENCE PREFIX",
        "If the input file contains sequences with multiple parts (eg. >experiment:run:cell), separated by the specified string (in this case, ':'), each distinct"
        " prefix will produce a separate file. If all sequences have the exact same first-level prefix, the second-level prefix will be used instead (and so on)."
        " The sequence names will not contain the prefix unless the '--keep_names' flag is present",
        "Split file by sequence prefix");
    args_add_help(result, "keep_names", "DO NOT MODIFY SEQUENCE NAMES",
        "Use this flag to ensure the sequence names stay the same as in the original file",
        "Do not modify sequence names");
    args_add_help(result, "output", "OUTPUT PREFIX",
        "Specify what should be prepended to output fiel names",
        "prefix for outputted files");
    args_add_help(result, "by_table", "SPLIT BY TABLE",
        "Specify the file which contains a two column-table (assumed to be tab-separated). The first columns should contain the read name, and the second column should contain the taerget file.",
        "Assign each read to a file according to a table");

    args_parse(result, argc, argv);
    return result;
}

#define MAX_SIMULTOPEN 50 - 3
int nucops_splitfile_a(args_t* args) {
    contig_s** contigs;
    char* inputfn;
    char* outputfn;
    char* separator;
    char sepchar;
    char* curgroupname;
    char* tmpstr;
    PF_t* outfile;
    PF_t* tablefile;
    PF_t** openfiles;
    char* tablefilename;
    char* basename;
    char** outfiles;
    char* filestatus;
    long long resume_at;
    size_t passid;
    size_t totalseqs, ncontigs, i, j, totalcontigs, nlines;
    int64_t maxseqs;
    int64_t maxbases;
    int64_t lastUID;
    int64_t curbases;
    int64_t curseqs;
    int64_t fileid;
    int fileUID;
    int _n, _at;
    int doabort;
    int opennewfile;
    size_t ignored_lines;
    DM64_t* read2fileid;
    DM64_t* filename2fileid;
    int nf;
    size_t noutfiles, nallocfiles;
    size_t nexttoclose;

    args_report_info(NULL, "nucops select is now runnning\n");

    inputfn = args_getstr(args, "input", 0, "stdin");
    outputfn = args_getstr(args, "output", 0, NULL);
    maxbases = args_getint(args, "max_bases_per_file", 0, -1);
    maxseqs = args_getint(args, "max_sequences_per_file", 0, -1);
    _n = args_ispresent(args, "keep_names");
    separator = args_getstr(args, "by_prefix", 0, NULL);
    _at = args_ispresent(args, "ids_only");

    if (strcmp(inputfn, "stdin") == 0) {
        args_report_warning(NULL, "Please note that when stdin is used, the entire input should be able to fit within the computer's RAM\n");
    }

    if (separator != NULL) {
        args_report_error(NULL, "by_prefix splitting is currently not implemented - aborting\n");
        return 1;
    }
    tablefile = NULL;
    tablefilename = NULL;
    if (args_ispresent(args, "by_table")) {
        tablefilename = args_getstr(args, "by_table", 0, NULL);
        if (!tablefilename) {
            args_report_error(NULL, "Please specify a table file - aborting\n");
            return 1;
        }
    }

    args_report_info(NULL, "Arguments were parsed sucessfully\n");

    read2fileid = NULL;
    filename2fileid = NULL;
    openfiles = NULL;
    outfiles = NULL;
    filestatus = NULL;
    nexttoclose = 0;
    if (tablefilename) {
        if (PFopen(&tablefile, tablefilename, "rb") != 0) {
            args_report_error(NULL, "Could not open table file\n");
        }
        else {
            read2fileid = new_DM64(0, 0);
            filename2fileid = new_DM64(0, 0);
            sepchar = 0;
            ignored_lines = 0;
            nallocfiles = 16;
            noutfiles = 0;
            outfiles = malloc(sizeof(char*) * nallocfiles);
            nlines = 0;
            while (tmpstr = PFreadline(tablefile)) {
                if (tmpstr[0] != '#') {
                    i = 0;
                    if (sepchar == 0) {
                        sepchar = '\t';
                        while (tmpstr[i] && tmpstr[i] != sepchar)i++;
                        if (tmpstr[i] == 0) {
                            i = 0;
                            sepchar = '\t';
                            while (tmpstr[i] && tmpstr[i] != sepchar)i++;

                        }
                    }
                    while (tmpstr[i] && tmpstr[i] != sepchar)i++;
                    if (tmpstr[i] != 0) {
                        tmpstr[i] = 0;
                        i++;
                        fileid = DM64_get(filename2fileid, tmpstr + i, (int)strlen(tmpstr + i), &nf);
                        if (nf) {
                            fileid = (int64_t)noutfiles;
                            DM64_assign(filename2fileid, tmpstr + i, (int)strlen(tmpstr + i), fileid);
                            if (noutfiles == nallocfiles) {
                                nallocfiles *= 2;
                                outfiles = realloc(outfiles, sizeof(char*) * nallocfiles);
                            }
                            noutfiles++;
                            outfiles[fileid] = memcpyalloc(tmpstr + i, strlen(tmpstr + i) + 1);
                        }
                        DM64_append(read2fileid, tmpstr, (int)strlen(tmpstr), fileid);
                    }
                    else {
                        ignored_lines++;
                    }
                    nlines++;
                    if (nlines % 1000000 == 0) {
                        args_report_info(NULL, _LLD_ "M entries parsed\n", (nlines-ignored_lines)/1000000);
                    }
                }
                free(tmpstr);
            }
            i = 0;
            if (ignored_lines > 0) {
                args_report_warning(NULL, "Some lines (" _LLD_ " out of " _LLD_ ") from the input table were improperly formatted\n", ignored_lines, nlines);
                args_report_warning(NULL, "  Note: comment lines (beginning with #) are not counted towards the total\n");
            }
            openfiles = malloc(sizeof(PF_t*)*noutfiles);
            filestatus = malloc(noutfiles);
            while (i < MAX_SIMULTOPEN && i < noutfiles) {
                PFopen(&(openfiles[i]), outfiles[i], "wb");
                filestatus[i] = 3; /*created and open*/
                i++;
            }
            while (i < noutfiles) {
                openfiles[i] = NULL;
                filestatus[i] = 0; 
                i++;
            }
            nexttoclose = 0;
            args_report_info(NULL, "Finished indexing read assignments\n", ignored_lines, nlines);
            args_report_info(NULL, "Potentially " _LLD_ " files will be created\n", noutfiles);
        }
        PFclose(tablefile);
        if (!read2fileid) {
            return 2;
        }
    }

    resume_at = 0;
    passid = 0;
    
    totalseqs = 0;
    totalcontigs = 0;
    lastUID = 0;
    /* parse the file */
    curbases = 0;
    curseqs = 0;
    if (outputfn != NULL) {
        basename = text_catstr(outputfn, "_");
    }
    else if (inputfn != NULL) {
        basename = memcpyalloc(inputfn, strlen(inputfn)+1);
        text_rmextention(basename);
        tmpstr = text_catstr(basename, "_");
        free(basename);
        basename = tmpstr;
        tmpstr = NULL;
    }
    else {
        basename = memcpyalloc("output_", strlen("output_")+1);
    }
    curgroupname = NULL;
    doabort = 0;
    outfile = NULL;
    if (separator == NULL) {
        tmpstr = text_catint(basename, 0);
        curgroupname = text_catstr(tmpstr, ".fastx");
        free(tmpstr);
        if (PFopen(&outfile, curgroupname, "wb")) {
            args_report_error(NULL, "Could not open target file - aborting\n");
            doabort = 1;
        }
        fileUID = 0;
    }
    if (!doabort) {
        do {
            args_report_progress(NULL, "Reading input... (pass %d)\n", passid);
            if (strcmp(inputfn, "stdin") == 0) {
                contigs = contigs_from_fastx(inputfn, &ncontigs);
                resume_at = 0;
            }
            else {
                contigs = contigs_from_fastx_limited(inputfn, &ncontigs, 0xC0000000, &resume_at);
            }
            totalseqs += ncontigs;
            args_report_progress(NULL, "Processing up to %dGb\n", (int)(resume_at / 1000000000));
            for (i = 0;i < ncontigs;i++) {
                contig_setUID((contigs[i]), lastUID);
                lastUID++;
            }
            totalcontigs += ncontigs;
            
            if (separator != NULL) {
                args_report_warning(NULL, "Impossible branch - if you get this error message, a mysterious bug has occured.\n");
            }
            else {
                for (i = 0;i < ncontigs;i++) {
                    nf = 1;
                    if (read2fileid) {
                        tmpstr = contig_nameptr(contigs[i]);
                        nf = 0;
                        if (!tmpstr)nf = 1;
                        else {
                            fileid = DM64_get(read2fileid, tmpstr, (int)strlen(tmpstr), &nf);
                            if (nf) {
                                nf = 0;
                                j = 0;
                                tmpstr = memcpyalloc(tmpstr, strlen(tmpstr) + 1);
                                while (tmpstr[j] && tmpstr[j] != ' ')j++;
                                tmpstr[j] = 0;
                                fileid = DM64_get(read2fileid, tmpstr, (int)strlen(tmpstr), &nf);
                                free(tmpstr);
                            }
                        }
                        if(!nf){
                            if (filestatus[fileid] < 2) {
                                /* file is not created, nor is it open */
                                if (filestatus[nexttoclose] / 2 == 1) {
                                    filestatus[nexttoclose] -= 2;
                                    PFclose(openfiles[nexttoclose]);
                                    openfiles[nexttoclose] = NULL;
                                    if (filestatus[fileid] == 0) {
                                        PFopen(&(openfiles[fileid]), outfiles[fileid], "wb");
                                        filestatus[fileid] |= 1;
                                    }
                                    else if (filestatus[fileid] == 1) {
                                        /* file hasb been created, but it is closed */
                                        PFopen(&(openfiles[fileid]), outfiles[fileid], "ab");
                                    }
                                    filestatus[fileid] |= 2;
                                    j = 0;
                                    while (j < MAX_SIMULTOPEN && filestatus[nexttoclose] / 2 == 0) {
                                        nexttoclose++;
                                        if (nexttoclose == noutfiles) nexttoclose = 0;
                                        j++;
                                    }
                                }
                            }
                            contigs_to_fastxPF(openfiles[fileid], &(contigs[i]), 1);
                        }
                    }
                    if (nf) {
                        opennewfile = 0;
                        if (maxseqs > 0) {
                            if (curseqs >= maxseqs) {
                                args_report_info(NULL, "Reached sequence limit.\n");
                                opennewfile = 1;
                                curseqs = 0;
                                curbases = 0;
                            }
                        }
                        if (maxbases > 0) {
                            if (curbases >= maxbases && curseqs > 0) {
                                args_report_info(NULL, "Reached base limit.\n");
                                opennewfile = 1;
                                curseqs = 0;
                                curbases = 0;
                            }
                        }
                        if (opennewfile) {
                            free(curgroupname);
                            PFclose(outfile);
                            fileUID++;
                            tmpstr = text_catint(basename, fileUID);
                            curgroupname = text_catstr(tmpstr, ".fastx");
                            args_report_info(NULL, "New file: %s\n", curgroupname);

                            free(tmpstr);
                            if (PFopen(&outfile, curgroupname, "wb")) {
                                args_report_error(NULL, "Could not open target file - aborting\n");
                                doabort = 1;
                                break;
                            }
                        }
                        contigs_to_fastxPF(outfile, &(contigs[i]), 1);
                        curseqs++;
                        curbases += contig_length(contigs[i]);
                    }
                    contig_free(contigs[i]);
                    contigs[i] = NULL;
                }
            }
            passid++;
            for (i = 0;i < ncontigs;i++) {
                if (contigs[i])contig_free(contigs[i]);
                contigs[i] = NULL;
            }
            free(contigs);
        } while (resume_at != 0 && !doabort);

        if(!doabort) args_report_info(NULL, "Output written sucessfuly.\n");
    }
    if (openfiles) {
        for (i = 0;i < noutfiles;i++) {
            if (openfiles[i])PFclose(openfiles[i]);
            openfiles[i] = NULL;
        }
        free(openfiles);
    }
    if (outfiles) {
        for (i = 0;i < noutfiles;i++) {
            if (outfiles[i])free(outfiles[i]);
            outfiles[i] = NULL;
        }
        free(outfiles);
    }
    if (filestatus) free(filestatus);
    if (read2fileid)free_DM64(read2fileid);
    if (filename2fileid)free_DM64(filename2fileid);
    if (curgroupname)free(curgroupname);
    if (basename)free(basename);
    return doabort;
}


int nucops_splitfile(int argc, char** argv) {
    int result;
    args_t* args;

    args = nucops_splitfile_init_args(argc, argv);
    result = 0;
    if (!args_ispresent(args, "help")) {
        result = nucops_splitfile_a(args);
        if (result != 0) {
            args_report_error(args, "splitfile failed with code <%d>\n", result);
        }
    }
    args_free(args);
    return result;
}

