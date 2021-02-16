#include <stdlib.h>
#include <string.h>
#include "argparser.h"
#include "sequence_base.h"
#include "textparsing.h"

args_t* nucops_fq2fa_init_args(int argc, char** argv) {
    args_t* result;
    result = args_alloc();

    args_add(result, "input", 'i', "str");
    args_add(result, "output", 'o', "str");
    args_add(result, "quality_threshold", 'q', "int");
    args_add(result, "base", 'b', "int");
    
    args_parse(result, argc, argv);
    return result;
}

/* redo this properly at some point */
int nucops_fq2fa_a(args_t* args) {
    char* inputfn;
    char* outputfn;
    size_t q_min;
    size_t base;
    PF_t* f_in;
    PF_t* f_out;

    size_t seqlen;
    size_t linelen;
    size_t i, pos;
    char* sequence;
    char* line;
    int failure_code;
    int curmode;
    
    inputfn = args_getstr(args, "input", 0, "stdin");
    outputfn = args_getstr(args, "output", 0, "stdout");
    q_min = args_getint(args, "quality_threshold", 0, 20);
    base = args_getint(args, "base", 0, 33);
    
    failure_code = 0;
    f_in = NULL;
    f_out = NULL;
    if (PFopen(&f_in, inputfn, "rb")) {
        failure_code = 1;
    }
    if (!failure_code && PFopen(&f_out, outputfn, "wb")) {
        failure_code = 2;
    }
    if (!failure_code) {
        sequence = NULL;
        curmode = -1;
        seqlen = 0;
        while (line = PFreadline(f_in)) {
            if(line[0]){
             curmode++;
             if(curmode==4) curmode=0;
             if(curmode==0) line[0]='>';
             if(curmode<2) PFprintf(f_out,"%s\n",line);
            }
            /*
            if (line[0] == '@' && (curmode==-1 || curmode == 3)) {
                if (curmode >= 0) {
                    PFwrite(sequence, 1, seqlen + 1, f_out);
                    PFputc('\n', f_out);
                }
                if (sequence)free(sequence);
                sequence = NULL;
                line[0] = '>';
                PFprintf(f_out, "%s\n", line);
                curmode = 1;
                linelen = 0;
                seqlen = 0;
            }
            else if (line[0] == '+' && curmode==1) {
                curmode = 2;
                pos = 0;
            }
            else if (curmode == 1) {
                linelen = strlen(line);
                sequence = (char*)realloc(sequence, linelen + seqlen + 1);
                memcpy(sequence + seqlen, line, linelen);
                seqlen += linelen;
                sequence[seqlen] = 0;
            }
            else if (curmode == 2) {
                linelen = strlen(line);
                for (i = 0;i < linelen;i++) {
                    if (line[i] < base + q_min) {
                        sequence[pos] = 'N';
                    }
                    pos++;
                }
               if(pos >= seqlen) curmode=3;
            }
*/
            free(line);
        }
        if (0 && curmode >= 0) {
            PFwrite(sequence, 1, seqlen + 1, f_out);
            PFputc('\n', f_out);
        }
        if (sequence)free(sequence);
    }
    PFclose(f_in);
    PFclose(f_out);
    return failure_code;
}

int nucops_fq2fa(int argc, char** argv) {
    int result;
    args_t* args;

    args = nucops_fq2fa_init_args(argc, argv);
    result = 0;
    if (!args_ispresent(args, "help")) {
        result = nucops_fq2fa_a(args);
        if (result != 0) {
            args_report_error(args, "fq2fa failed with code <%d>\n", result);
        }
    }
    args_free(args);
    return result;
}
