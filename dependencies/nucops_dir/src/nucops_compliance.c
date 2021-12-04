#include "argparser.h"
#include "sequence_base.h"
#include "vecops.h"
#include "osportstd.h"
#include "nucops_compliance.h"
#include "textparsing.h"
#include <string.h>
#include <stdio.h>

args_t* nucops_compliance_init_args(int argc, char** argv) {
    args_t* result;
    result = args_alloc();

    args_add_help(result, NULL, "POSITIONAL ARGUMENTS", "Only one positional argument is examined. Positional argument should be the name of the file to examine.\n", "File to examine");
    args_add(result, "input", 'i', "str");
    args_add_help(result, NULL, "INPUT FILE", "should be a list of files to examine\n", "File to examine");
    args_add(result, "format", 'f', "str,str");
    args_add_help(result, "format", "EXPECTED INPUT FORMAT", "Specify the expected input format\nThe following formats are available:\n fasta\n fastq\n","expected input format");
    args_add(result, "mode", 'm', "str,str");
    args_add_help(result, "mode", "EXPECTED INPUT MODE", "Specify the expected input mode\nThe following formats are available:\n nucleotide\n", "expected input mode");
    args_add(result, "corrected", 'c', "str");
    args_add_help(result, "corrected", "CORRECTED OUTPUT", "Used to specify the name of the corrected ouput (can be stdout).\nCorrected output is not produced if this flag is absent.\n","Output file name.");
    args_add(result, "output", 'o', "str");
    args_add_help(result, "output", "CORRECTED OUTPUT", "Used to specify the name of the corrected ouput (can be stdout).\nCorrected output is not produced if this flag is absent.\n", "''");
    args_add(result, "lax", 'l', "");
    args_add_help(result, "lax", "LAX NUCLEOTIDE MODE", "Enable lax mode - in this mode, sequences may also contain B,D,H,M,R,K,S,Y,U,v and W without raising errors\n", "View B,D,H,M,R,Y,K,S,U,V,W as valid nucleotides");

    args_parse(result, argc, argv);
    return result;
}

#define MODE_NUCLEOTIDE 0
#define MODE_LAXNUCLEOTIDE 2

int correct_invalid_chars(char* line, size_t lineid) {
    size_t i;
    int invalids_found;
    i = 0;
    invalids_found = 0;
    while (line[i]) {
        if (line[i] < ' ' || line[i] > '~') {
            args_report_info(NULL, "found \\%d @ loc:%d in line \"%s\" \n",(int)line[i],(int)i,line);
            invalids_found = 1;
            line[i] = '_';
        }
        i++;
    }
    if (invalids_found) {
        fprintf(stderr, "Non-printable character(s) found");
        fprintf(stderr, " (line " _LLD_ ")\n", lineid);
    }
    return invalids_found;
}

int correct_ACGTN(char* line, size_t lineid) {
    size_t i;
    int invalids_found;
    size_t shorten_by;
    i = 0;
    invalids_found = 0;
    shorten_by = 0;
    while (line[i]) {
        if (shorten_by > 0)line[i - shorten_by] = line[i];
        if (line[i] != 'A' && line[i] != 'a' && line[i] != 'C' && line[i] != 'c' && line[i] != 'G' && line[i] != 'g' && line[i] != 'T' && line[i] != 't' && line[i] != 'N' && line[i] != 'n') {
            if (line[i] != '_')line[i] = 'N';
            else shorten_by++;
            invalids_found = 1;
        }
        i++;
    }
    if (invalids_found) {
        fprintf(stderr, "Sequence contains characters otehr than A,C,G,T or N");
        fprintf(stderr, " (line " _LLD_ ")\n", lineid);
    }
    return invalids_found;
}
int check_ACGTNjunk(char* line, size_t lineid) {
    size_t i;
    size_t invalids_found;
    size_t shorten_by;
    int invalids_are_alt;
    i = 0;
    invalids_found = 0;
    shorten_by = 0;
    invalids_are_alt=1;
    while (line[i]) {
        if (shorten_by > 0)line[i - shorten_by] = line[i];
        if (line[i] != 'A' && line[i] != 'a' && line[i] != 'C' && line[i] != 'c' && line[i] != 'G' && line[i] != 'g' && line[i] != 'T' && line[i] != 't' && line[i] != 'N' && line[i] != 'n') {
            invalids_found++;
            if(line[i] != 'K' && line[i] != 'k' && line[i] != 'R' && line[i] != 'r' && line[i] != 'S' && line[i] != 's' && line[i] != 'W' && line[i] != 'w' && line[i] != 'Y' && line[i] != 'y' && line[i] != 'M' && line[i] != 'm' && line[i] != 'H' && line[i] != 'h' && line[i] != 'B' && line[i] != 'b' && line[i] != 'D' && line[i] != 'd' && line[i] != 'U' && line[i] != 'u' && line[i] != 'V' && line[i] != 'v')invalids_are_alt=0;
        }
        i++;
    }
    if (invalids_found) {
        if (invalids_found > i / 2 && !invalids_are_alt) {
            fprintf(stderr, "Sequence line contains junk");
            fprintf(stderr, " (line " _LLD_ ")\n", lineid);
            return 1;
        }
    }
    return 0;
}
int correct_ACGTNKRYSW(char* line, size_t lineid) {
    size_t i;
    int invalids_found;
    size_t shorten_by;
    i = 0;
    invalids_found = 0;
    shorten_by = 0;
    while (line[i]) {
        if (shorten_by > 0)line[i - shorten_by] = line[i];
        if (line[i] != 'A' && line[i] != 'a' && line[i] != 'C' && line[i] != 'c' && line[i] != 'G' && line[i] != 'g' && line[i] != 'T' && line[i] != 't' && line[i] != 'N' && line[i] != 'n' &&
            line[i] != 'K' && line[i] != 'k' && line[i] != 'R' && line[i] != 'r' && line[i] != 'S' && line[i] != 's' && line[i] != 'W' && line[i] != 'w' && line[i] != 'Y' && line[i] != 'y' && line[i] != 'M' && line[i] != 'm' && line[i] != 'H' && line[i] != 'h' && line[i] != 'B' && line[i] != 'b' && line[i] != 'D' && line[i] != 'd' && line[i] != 'U' && line[i] != 'u' && line[i] != 'V' && line[i] != 'v') {
            if (line[i] != '_')line[i] = 'N';
            else shorten_by++;
            invalids_found = 1;
        }
        i++;
    }
    if (invalids_found) {
        fprintf(stderr, "Sequence contains characters other than A,B,C,D,G,H,R,S,T,U,V,M,Y,W or N");
        fprintf(stderr, " (line " _LLD_ ")\n", lineid);
    }
    return invalids_found;
}
int check_ACGTNKRYSWjunk(char* line, size_t lineid) {
    size_t i;
    size_t invalids_found;
    size_t shorten_by;
    i = 0;
    invalids_found = 0;
    shorten_by = 0;
    while (line[i]) {
        if (shorten_by > 0)line[i - shorten_by] = line[i];
        if (line[i] != 'A' && line[i] != 'a' && line[i] != 'C' && line[i] != 'c' && line[i] != 'G' && line[i] != 'g' && line[i] != 'T' && line[i] != 't' && line[i] != 'N' && line[i] != 'n' &&
            line[i] != 'K' && line[i] != 'k' && line[i] != 'R' && line[i] != 'r' && line[i] != 'S' && line[i] != 's' && line[i] != 'W' && line[i] != 'w' && line[i] != 'Y' && line[i] != 'y' && line[i] != 'M' && line[i] != 'm' && line[i] != 'H' && line[i] != 'h' && line[i] != 'B' && line[i] != 'b' && line[i] != 'D' && line[i] != 'd' && line[i] != 'U' && line[i] != 'u' && line[i] != 'V' && line[i] != 'v') {
            invalids_found++;
        }
        i++;
    }
    if (invalids_found) {
        if (invalids_found > i / 2) {
            fprintf(stderr, "Sequence line contains junk");
            fprintf(stderr, " (line " _LLD_ ")\n", lineid);
            return 1;
        }
    }
    return 0;
}

int local__is_fasta_compliant(PF_s* infile, PF_s* outfile, int mode) {
    int res;
    char* line;
    char* previous_line;
    char* curseq;
    int last_line_type;
    int line_type;
    size_t lasti;
    size_t lastlinelen;
    size_t lineid;
    size_t linelen;
    size_t curseqalloc;
    size_t curseqlen;
    size_t ncols;
    size_t mlseq_found;
    size_t nseq;
    
    res = 1;
    last_line_type = 0;
    line_type = 0;
    lineid = 0;
    curseqalloc = 150;
    curseqlen = 0;
    curseq = malloc(curseqalloc);
    ncols = 0;
    lastlinelen = 0;
    mlseq_found = 0;
    previous_line = NULL;
    nseq=0;
    while (line = PFreadline(infile)) {
        if(lineid%100000==0){
           args_report_progress(NULL,_LLD_ " lines (" _LLD_ " sequences) parsed\n",lineid,nseq);
        }
        lineid++;
        last_line_type = line_type;
        if (line[0] == 0) {
            fprintf(stderr, "Empty line");
            fprintf(stderr, " (line " _LLD_ ")\n", lineid);
            free(line);
            res = 0;
        }
        else{
            if (correct_invalid_chars(line, lineid)) res = 0;
            if (last_line_type == 0) {
                if (line[0] != '>') {
                    fprintf(stderr, "Missing expected fasta header");
                    fprintf(stderr, " (line " _LLD_ ")\n", lineid);
                    res = 0;
                }
                else {
                    line_type = 1;
                    nseq++;
                }
            }
            else if (last_line_type == 1) {
                if (line[0] == '>') {
                    fprintf(stderr, "Sequence data is missing");
                    fprintf(stderr, " (line " _LLD_ ")\n", lineid);
                    res = 0;
                }
                else {
                    line_type = 2;
                }
            }
            else if (last_line_type == 2) {
                if (line[0] == '>') {
                    if (outfile) {
                        lasti = 0;
                        if (ncols > 0 && curseqlen > ncols) {
                            while (curseqlen > lasti + ncols) {
                                PFwrite(curseq + lasti, 1, ncols, outfile);
                                PFputc('\n',outfile);
                                lasti += ncols;
                            }
                        }
                        PFprintf(outfile, "%s\n", curseq+lasti);
                    }
                    line_type = 1;
                    nseq++;
                }
                else {
                    if(mode == MODE_NUCLEOTIDE){
                      if (check_ACGTNjunk(line,lineid)) {
                        line_type = 1;
                        res = 0;
                        linelen = strlen(line);
                        line = realloc(line, linelen + 2);
                        linelen++;
                        while (linelen > 0) {
                            line[linelen] = line[linelen - 1];
                        }
                        line[0] = '>';
                      }
                    }
                    else if (mode == MODE_LAXNUCLEOTIDE) {
                      if (check_ACGTNKRYSWjunk(line,lineid)) {
                        line_type = 1;
                        res = 0;
                        linelen = strlen(line);
                        line = realloc(line, linelen + 2);
                        linelen++;
                        while (linelen > 0) {
                            line[linelen] = line[linelen - 1];
                        }
                        line[0] = '>';
                      }

                    }
                    line_type = 2;
                }
            }

            if (line_type == 1) {
                curseqlen = 0;
            }
            else if (line_type == 2) {
                if (last_line_type == 1) {
                    if (outfile && previous_line) {
                        PFprintf(outfile, "%s\n", previous_line);
                    }
                }
                if (mode == MODE_NUCLEOTIDE) {
                    if (correct_ACGTN(line, lineid))res = 0;
                }
                else if (mode == MODE_LAXNUCLEOTIDE) {
                    if (correct_ACGTNKRYSW(line, lineid))res = 0;
                }
                linelen = strlen(line);

                if (mlseq_found && linelen > ncols) {
                    fprintf(stderr, "Line is longer than expected");
                    fprintf(stderr, " (line " _LLD_ ")\n", lineid - 1);
                    res = 0;
                }
                if (last_line_type == 2) {
                    if (lastlinelen < linelen || (mlseq_found && lastlinelen < ncols)) {
                        fprintf(stderr, "Line is shorter than expected");
                        fprintf(stderr, " (line " _LLD_ ")\n", lineid - 1);
                        res = 0;
                    }
                    if (!mlseq_found) {
                        ncols = lastlinelen;
                        mlseq_found = 1;
                    }
                }
                if (curseqlen + linelen + 1 >= curseqalloc) {
                    curseqalloc = curseqlen + linelen + 2;
                    curseq = realloc(curseq, curseqalloc);
                }
                memcpy(curseq + curseqlen, line, linelen+1);
                curseqlen += linelen;
                lastlinelen = linelen;
            }
            if (previous_line)free(previous_line);
            previous_line = line;
        }
    }
    if (previous_line && line_type==2) {
        if (outfile) {
            lasti = 0;
            if (ncols > 0 && curseqlen > ncols) {
                while (curseqlen > lasti + ncols) {
                    PFwrite(curseq + lasti, 1, ncols, outfile);
                    PFputc('\n', outfile);
                    lasti += ncols;
                }
            }
            PFprintf(outfile, "%s\n", curseq + lasti);
        }
    }
    if (previous_line) free(previous_line);
    free(curseq);

    args_report_info(NULL, "Finished parsing: %s\n",res==1?"success":"failure");

    return res;
}

#define MODE_NUCREADS_LONG  0
#define MODE_NUCREADS_FIXED 1

int local__is_fastq_compliant(PF_s* infile, PF_s* outfile, int mode) {
    int res;
    char* line;
    char* previous_line;
    char* curseq;
    char* qualscoreline;
    char* readname;
    int last_line_type;
    int line_type;
    size_t lastlinelen;
    size_t lineid;
    size_t linelen;
    size_t curseqalloc;
    size_t curseqlen;
    size_t ncols;
    size_t mlseq_found;
    size_t read_length;
    size_t qscorlen;
    int unfixable;

    res = 1;
    last_line_type = 0;
    line_type = 0;
    lineid = 0;
    curseqalloc = 150;
    curseqlen = 0;
    curseq = calloc(curseqalloc,1);
    qualscoreline = calloc(curseqalloc,1);
    ncols = 0;
    lastlinelen = 0;
    mlseq_found = 0;
    previous_line = NULL;
    read_length = 0;
    unfixable = 0;
    readname = NULL;
    qscorlen = 0;
    while (line = PFreadline(infile)) {
        lineid++;
        last_line_type = line_type;
        if (line[0] == 0) {
            fprintf(stderr, "Empty line");
            fprintf(stderr, " (line " _LLD_ ")\n", lineid);
            free(line);
            res = 0;
        }
        else {
            if (correct_invalid_chars(line, lineid)) res = 0;
            if (last_line_type == 0) {
                if (line[0] != '@') {
                    fprintf(stderr, "Missing expected fasta header");
                    fprintf(stderr, " (line " _LLD_ ")\n", lineid);
                    res = 0;
                }
                else {
                    line_type = 1;
                }
            }
            else if (last_line_type == 1) {
                if (line[0] == '@') {
                    fprintf(stderr, "Sequence data is missing");
                    fprintf(stderr, " (line " _LLD_ ")\n", lineid);
                    res = 0;
                }
                else {
                    line_type = 2;
                }
            }
            else if (last_line_type == 2) {
                if (line[0] == '@') {
                    fprintf(stderr, "Missing quality Scores");
                    fprintf(stderr, " (line " _LLD_ ")\n", lineid);
                    res = 0;
                    if (outfile) {
                        memset(qualscoreline, '!', curseqlen);
                        PFprintf(outfile, "%s\n+\n%s\n", curseq, qualscoreline);
                    }
                    line_type = 1;
                }
                else if(line[0]=='+'){
                    line_type = 3;
                }
                if (line_type != 2) {
                    if (mode == MODE_NUCREADS_FIXED) {
                        if (read_length == 0)read_length = curseqlen;
                        if (curseqlen != read_length && !unfixable) {
                            unfixable = 1;
                            fprintf(stderr, "Unfixable : length of read differs from previous values");
                            fprintf(stderr, " (line " _LLD_ ")\n", lineid);
                            res = 0;
                        }
                    }
                }
            }
            else if (last_line_type == 3) {
                line_type = 4;
            }
            else if (last_line_type == 4) {
                if (qscorlen != curseqlen) {
                    res = 0;
                    linelen = strlen(line);
                    if (line[0] == '@') {
                        if (linelen + qscorlen == curseqlen) {
                            memcpy(qualscoreline + qscorlen, line, linelen + 1);
                            fprintf(stderr, "Unfixable : multiline quality score with subsequent line beginning with '@'");
                            fprintf(stderr, " (line " _LLD_ ")\n", lineid);
                            if(!unfixable) unfixable = 2;
                        }
                        else {
                            fprintf(stderr, "Quality score line is of insufficient length");
                            fprintf(stderr, " (line " _LLD_ ")\n", lineid);
                            memset(qualscoreline + qscorlen, '!', curseqlen - qscorlen);
                            line_type = 1;
                        }
                    }
                }
                else {
                    line_type = 1;
                }
            }

            if (line_type == 1) {
                if (outfile && curseqlen > 0) {
                    PFprintf(outfile, "%s\n+\n%s\n", curseq, qualscoreline);
                }
                linelen = strlen(line);
                if (readname!=NULL)free(readname);
                readname = malloc(linelen+1);
                memcpy(readname, line, linelen);
                curseqlen = 0;
            }
            else if (line_type == 2) {
                if (last_line_type == 1) {
                    if (outfile && previous_line) {
                        PFprintf(outfile, "%s\n", previous_line);
                    }
                }
                if (mode == MODE_NUCREADS_FIXED || mode == MODE_NUCREADS_LONG) {
                    if (correct_ACGTN(line, lineid))res = 0;
                }
                linelen = strlen(line);

                if (mlseq_found && linelen > ncols) {
                    fprintf(stderr, "Line is longer than expected");
                    fprintf(stderr, " (line " _LLD_ ")\n", lineid - 1);
                    res = 0;
                }
                if (last_line_type == 2) {
                    if (lastlinelen < linelen || (mlseq_found && lastlinelen < ncols)) {
                        fprintf(stderr, "Line is shorter than expected");
                        fprintf(stderr, " (line " _LLD_ ")\n", lineid - 1);
                        res = 0;
                    }
                    if (!mlseq_found) {
                        fprintf(stderr, "Multiline sequence found");
                        fprintf(stderr, " (line " _LLD_ ")\n", lineid - 1);
                        res = 0;
                        ncols = lastlinelen;
                        mlseq_found = 1;
                    }
                }
                if (curseqlen + linelen + 1 >= curseqalloc) {
                    curseqalloc = curseqlen + linelen + 2;
                    curseq = realloc(curseq, curseqalloc);
                    qualscoreline = realloc(qualscoreline, curseqalloc);
                }
                memcpy(curseq + curseqlen, line, linelen + 1);
                curseqlen += linelen;
                lastlinelen = linelen;
            }
            else if (line_type == 3) {
                qscorlen = 0;
                memset(qualscoreline, '!', curseqlen);
                qualscoreline[curseqlen] = 0;
            }
            else if (line_type == 4) {
                linelen = strlen(line);
                if (linelen + qscorlen == curseqlen) {
                    memcpy(qualscoreline + qscorlen, line, linelen + 1);
                    qscorlen += linelen;
                }
                else {
                    if (linelen + qscorlen < curseqlen) {
                        fprintf(stderr, "Quality score line is of incorrect length (too short)");
                        fprintf(stderr, " (line " _LLD_ ")\n", lineid);
                        memcpy(qualscoreline + qscorlen, line, linelen);
                        qscorlen += linelen;
                    }
                    else {
                        fprintf(stderr, "Quality score line is of incorrect length (too long)");
                        fprintf(stderr, " (line " _LLD_ ")\n", lineid);
                        memcpy(qualscoreline + qscorlen, line, curseqlen - qscorlen);
                        qscorlen += linelen;
                        line_type = 0;
                    }
                }
            }
            if (previous_line)free(previous_line);
            previous_line = line;
        }
    }
    if (previous_line && line_type == 4) {
        if (outfile) {
            PFprintf(outfile, "%s\n+\n%s\n", curseq, qualscoreline);
        }
    }
    if (readname)free(readname);
    if (previous_line) free(previous_line);
    free(curseq);
    free(qualscoreline);

    args_report_info(NULL, "Finished parsing: %s\n", res == 1 ? "success" : "failure");
    if (unfixable == 1) {
        fprintf(stderr, "Detected at least one unfixable error of type 1: length of read differs from previous values");
    }
    if (unfixable == 2) {
        fprintf(stderr, "Detected at least one unfixable error of type 2: multiline quality score with subsequent line beginning with '@'");
    }

    return res;
}


int nucops_compliance_a(args_t* args) {
    PF_s* infile;
    PF_s* outfile;
    char* inputfile;
    char* outputfile;
    char* filetypestr;
    char* modestr;
    int input_is_stdin;
    int mode;
    int res;

    inputfile = args_getstr(args, "input", 0, NULL);
    if(!inputfile) inputfile = args_getstr(args, NULL, 0, "stdin");

    outputfile = args_getstr(args, "corrected", 0, NULL);
    if(!outputfile) outputfile = args_getstr(args, "output", 0, NULL);

    filetypestr = args_getstr(args, "format", 0, NULL);
    modestr = args_getstr(args, "mode", 0, NULL);

    if (filetypestr == NULL) {
        args_report_info(NULL, "Attempting identify file format based on file name : %s\n", inputfile);
        if (endswith(".fasta", inputfile))filetypestr = "fasta";
        if (endswith(".fna", inputfile))filetypestr = "fna";
        if (endswith(".fa", inputfile))filetypestr = "fasta";
        if (endswith(".fastq", inputfile))filetypestr = "fastq";
        if (endswith(".fq", inputfile))filetypestr = "fastq";
        if (!filetypestr) {
            args_report_error(NULL, "Could not determine input format from file name. Please specify with -f\n");
            return 1;
        }
    }
    if (strcmp(filetypestr, "fna") == 0) {
        if (!modestr)modestr = "nucleotide";
        filetypestr = "fasta";
    }
    if (strcmp(filetypestr, "fasta")!=0 && strcmp(filetypestr, "fastq")!=0) {
        args_report_error(NULL, "unknown input format\n");
        return 1;
    }
    if (!modestr)modestr = "nucleotide";
    if (strcmp(modestr, "nucleotide") == 0) {
        if(args_ispresent(args,"lax")){
            args_report_info(NULL, "Lax nucleotide mode is enabled");
            mode = MODE_LAXNUCLEOTIDE;
        }
        else { 
            mode = MODE_NUCLEOTIDE;
        }
    }
    else {
        args_report_error(NULL, "unknown input mode\n");
        return 1;
    }


    outfile = NULL;
    if (PFopen(&infile, inputfile, "rb")) {
        args_report_error(NULL, "Input file could not be read - aborting\n");
        return 1;
    }
    if (outputfile) {
        res = PFopen(&outfile, outputfile, "wb");
        if (res) {
            args_report_warning(NULL, "Output file could not be created (error code: %d) - corrected file will not be produced\n",res);
            res = 0;
            PFclose(outfile);
            outfile = NULL;
            outputfile = NULL;
        }
    }
    if (strcmp(inputfile, "stdin")) {
        input_is_stdin = 1;
    }
    else {
        input_is_stdin = 0;
    }

    args_report_info(NULL, "Arguments parsed successfully\n", inputfile);
    res = 0;
    if (strcmp(filetypestr, "fasta") == 0) {
        if (local__is_fasta_compliant(infile, outfile, mode)) {
            if(outputfile && strcmp(outputfile,"stdout")==0) fprintf(stderr, "Input file was compliant\n");
            else printf("Input file was compliant\n");
        }
        else {
            if (outputfile && strcmp(outputfile, "stdout") == 0) fprintf(stderr, "Input file was not fully compliant\n");
            else fprintf(stdout, "Input file was not fully compliant\n");
        }
        if (outputfile) {
            if(strcmp(outputfile, "stdout") == 0) fprintf(stderr, "Output is stored in %s\n",outputfile);
            else fprintf(stdout, "Output is stored in %s\n", outputfile);
        }
    }
    else if (strcmp(filetypestr, "fastq") == 0) {
        if (local__is_fastq_compliant(infile, outfile, mode)) {
            if (outputfile && strcmp(outputfile, "stdout") == 0) fprintf(stderr, "Input file was compliant\n");
            else printf("Input file was compliant\n");
        }
        else {
            if (outputfile && strcmp(outputfile, "stdout") == 0) fprintf(stderr, "Input file was not fully compliant\n");
            else fprintf(stdout, "Input file was not fully compliant\n");
        }
        if (outputfile) {
            if (strcmp(outputfile, "stdout") == 0) fprintf(stderr, "Output is stored in %s\n", outputfile);
            else fprintf(stdout, "Output is stored in %s\n", outputfile);
        }
    }
    else {
        args_report_error(NULL, "Unknown file format\n");
        res = 1;
    }

    if (outputfile) {
        PFclose(outfile);
    }
    if (!input_is_stdin) {
        PFclose(infile);
    }
    args_report_info(NULL, "Done.\n", inputfile);
    return res;
}

int nucops_compliance(int argc, char** argv) {
    int result;
    args_t* args;

    args = nucops_compliance_init_args(argc, argv);
    result = 0;
    if (!args_ispresent(args, "help")) {
        result = nucops_compliance_a(args);
        if (result != 0) {
            args_report_error(args, "Compliance failed with code <%d>\n", result);
        }
    }
    if(args)args_free(args);
    return result;
}
