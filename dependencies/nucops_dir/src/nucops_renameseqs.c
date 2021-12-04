#include "argparser.h"
#include "osportstd.h"
#include "sequence_base.h"
#include <string.h>
#include "textparsing.h"

args_t* nucops_renameseq_init_args(int argc, char** argv) {
    args_t* result;
    result = args_alloc();

    args_add(result, "input", 'i', "str");
    args_add(result, "output", 'o', "str");

    args_add(result, "format", 'f', "str");
    args_add(result, "keywords", 'k', "");

    args_parse(result, argc, argv);
    return result;
}

int nucops_renameseq_a(args_t* args) {
    contig_s* contig;
    char* inputfn;
    char* outputfn;
    char* outputstr;
    char* fmtstr;
    char* fmtstr_parsed;
    char* tmpstr1;
    char* tmpstr2;
    PF_s* fin;
    PF_s* fout;
    int res, nokey, openkey;
    int ninvalid;
    uint16_t kwid;
    size_t nalloc, nalloc2, nalloc3, bracketend, i, j, k, jmax, uid;

    inputfn = args_getstr(args, "input", 0, "stdin");
    outputfn = args_getstr(args, "output", 0, "stdout");
    res = 0;

    if (PFopen(&fin, inputfn, "rb") != 0) {
        args_report_error(NULL, "Could not open %s for reading\n", inputfn);
        res = 1;
    }
    if (PFopen(&fout, outputfn, "wb") != 0) {
        args_report_error(NULL, "Could not open %s for writing\n", outputfn);
        if(!res) res = 2;
    }
    
    fmtstr = args_getstr(args, "format", 0, "{file}.{id}_len.{length}");

    nalloc = 256;
    fmtstr_parsed = calloc(nalloc, 1);
    nalloc2 = 256;
    tmpstr2 = calloc(nalloc2, 1);
    
    tmpstr1 = fmtstr;
    i = 0;
    nokey = 1;
    ninvalid = 0;
    /* parse the format string */
    while (*tmpstr1) {
        if (*tmpstr1 < ' ' || *tmpstr1 > '~' ) {
            args_report_warning(NULL, "Invalid character \\%d encountered in format string, was ignored\n", (int)(unsigned char)(*tmpstr1));
        }
        else if (*tmpstr1 == '{') {
            if (i >= nalloc-4) {
                fmtstr_parsed = realloc(fmtstr_parsed, 2 * nalloc);
                memset(fmtstr_parsed + i, 0, (2 * nalloc) - i);
                nalloc *= 2;
            }
            kwid = 0;
            bracketend = 1;
            while (tmpstr1[bracketend] && tmpstr1[bracketend] != '}')bracketend++;
            if (tmpstr1[bracketend])openkey = 0;
            else openkey = 1;
            if (bracketend >= nalloc2) {
                nalloc2 = bracketend+1;
                tmpstr2 = realloc(tmpstr2, nalloc2);
            }
            tmpstr1 = tmpstr1 + 1;
            bracketend--;
            memcpy(tmpstr2, tmpstr1, bracketend);
            tmpstr2[bracketend] = 0;
            kwid = 0;
            if (strcmp(tmpstr2, "sp") == 0 || strcmp(tmpstr2, "space") == 0) {
                fmtstr_parsed[i] = ' ';
                i++;
            }
            else if (strcmp(tmpstr2, "o") == 0 || strcmp(tmpstr2, "open_curly") == 0) {
                fmtstr_parsed[i] = '{';
                i++;
            }
            else if (strcmp(tmpstr2, "c") == 0 || strcmp(tmpstr2, "close_curly") == 0) {
                fmtstr_parsed[i] = '}';
                i++;
            }
            else if (strcmp(tmpstr2, "q") == 0 || strcmp(tmpstr2, "quote") == 0) {
                fmtstr_parsed[i] = '\'';
                i++;
            }
            else if (strcmp(tmpstr2, "qq") == 0 || strcmp(tmpstr2, "double_quote") == 0) {
                fmtstr_parsed[i] = '"';
                i++;
            }
            else if (strcmp(tmpstr2, "bq") == 0 || strcmp(tmpstr2, "back_quote") == 0) {
                fmtstr_parsed[i] = '`';
                i++;
            }
            else if (strcmp(tmpstr2, "s") == 0 || strcmp(tmpstr2, "slash") == 0 || strcmp(tmpstr2,"stroke") == 0) {
                fmtstr_parsed[i] = '/';
                i++;
            }
            else if (strcmp(tmpstr2, "bs") == 0 || strcmp(tmpstr2, "backslash") == 0 || strcmp(tmpstr2, "backstroke") == 0) {
                fmtstr_parsed[i] = '\\';
                i++;
            }
            else if (strcmp(tmpstr2, "bang") == 0 || strcmp(tmpstr2, "b") == 0) {
                fmtstr_parsed[i] = '!';
                i++;
            }
            else if (strcmp(tmpstr2, "hash") == 0 || strcmp(tmpstr2, "pound") == 0 || strcmp(tmpstr2, "p") == 0 || strcmp(tmpstr2, "h") == 0) {
                fmtstr_parsed[i] = '#';
                i++;
            }
            else if (strcmp(tmpstr2, "dollar") == 0 || strcmp(tmpstr2, "d") == 0) {
                fmtstr_parsed[i] = '$';
                i++;
            }
            else if (strcmp(tmpstr2, "colon") == 0 || strcmp(tmpstr2, "n") == 0) {
                fmtstr_parsed[i] = ':';
                i++;
            }
            else if (strcmp(tmpstr2, "semicolon") == 0 || strcmp(tmpstr2, "sn") == 0) {
                fmtstr_parsed[i] = ';';
                i++;
            }
            else if (strcmp(tmpstr2, "less_than") == 0 || strcmp(tmpstr2, "lt") == 0) {
                fmtstr_parsed[i] = '<';
                i++;
            }
            else if (strcmp(tmpstr2, "greater_than") == 0 || strcmp(tmpstr2, "gt") == 0) {
                fmtstr_parsed[i] = '>';
                i++;
            }
            else if (strcmp(tmpstr2, "percent") == 0 || strcmp(tmpstr2, "pc") == 0) {
                fmtstr_parsed[i] = '%';
                i++;
            }
            else if (strcmp(tmpstr2, "file") == 0) {
                j = strlen(PFgetbasenameptr(fin));
                if (i + j + 1 >= nalloc) {
                    fmtstr_parsed = realloc(fmtstr_parsed, i + j + nalloc + 2);
                    memset(fmtstr_parsed + i, 0, j + nalloc + 1);
                    nalloc = i + j + nalloc + 2;
                }
                memcpy(fmtstr_parsed + i, PFgetbasenameptr(fin), j);
                jmax = j;
                j = j + i;
                while (j > i) {
                    j--;
                    if (fmtstr_parsed[j] == '.')break;
                }
                if (j == 0)j = jmax;
                j -= i;
                i += j;
            }
            else if (strcmp(tmpstr2, "name") == 0) {
                kwid = 1;
            }
            else if (strcmp(tmpstr2, "_name_") == 0) {
                kwid = 2;
            }
            else if (strcmp(tmpstr2, ".name.") == 0) {
                kwid = 3;
            }
            else if (strcmp(tmpstr2, "prefix") == 0) {
                kwid = 4;
            }
            else if (strcmp(tmpstr2, "suffix") == 0) {
                kwid = 5;
            }
            else if (strcmp(tmpstr2, "id") == 0) {
                kwid = 6;
            }
            else if (strcmp(tmpstr2, "length") == 0) {
                kwid = 7;
            }
            else if (strcmp(tmpstr2, "quality") == 0) {
                kwid = 8;
            }
            else if (strcmp(tmpstr2, "gc") == 0) {
                kwid = 9;
            }
            else if (strcmp(tmpstr2, "pacbio1") == 0) {
                kwid = 10;
            }
            else if (strcmp(tmpstr2, "pacbio2") == 0) {
                kwid = 11;
            }
            else if (strcmp(tmpstr2, "pacbio3") == 0) {
                kwid = 12;
            }
            else {
                args_report_info(NULL, "Invalid keyword %s\n", tmpstr2);
                ninvalid++;
            }
            if (kwid > 0) {
                args_report_info(NULL, "Detected keyword: %s (%d)\n", tmpstr2, (int)kwid);
                fmtstr_parsed[i] = 0x10; /* DLE character */
                i++;
                memcpy(fmtstr_parsed+i, &kwid, sizeof(kwid));
                i += sizeof(kwid);
                nokey = 0;
            }
            if (openkey) {
                ninvalid++;
            }
            else {
                tmpstr1 += bracketend;
            }
        }
        else {
            if (i == nalloc - 1) {
                fmtstr_parsed = realloc(fmtstr_parsed, 2 * nalloc);
                memset(fmtstr_parsed + i, 0, nalloc+1);
                nalloc *= 2;
            }
            fmtstr_parsed[i] = *tmpstr1;
            i++;
        }
        tmpstr1 = tmpstr1+1;
    }
    fmtstr_parsed[i] = 0;
    
    if (ninvalid > 0) {
        args_report_error(NULL, "Invalid keywords detected (use '-q none' too see which)\n");
        if(!res)res = 3;
    }

    if (!res && nokey) {
        args_report_warning(NULL, "Format string does not contain a variable part - all sequence names will be the same!\n");
    }
    if (!res) {
        nalloc3 = nalloc;
        outputstr = malloc(nalloc3);
        uid = 0;
        while (contig = contig_from_fastxPF(fin)) {
            uid++;
            if (uid % 10000 == 0) {
                args_report_progress(NULL, "%dk sequences renamed\n", (int)(uid / 1000));
            }
            tmpstr1 = fmtstr_parsed;
            i = 0;
            while (*tmpstr1) {
                if (*tmpstr1 == 0x10) {
                    memcpy(&kwid, tmpstr1 + 1, sizeof(kwid));
                    switch (kwid) {
                        /* name */
                        case 1:
                        case 2:
                        case 3:
                        case 4:
                        case 5:
                            j = strlen(contig_nameptr(contig));
                            if (j + 1 >= nalloc2) {
                                nalloc2 = j + 2;
                                tmpstr2 = realloc(tmpstr2, nalloc2);
                            }
                            memcpy(tmpstr2, contig_nameptr(contig), j);
                            tmpstr2[j] = 0;
                            if (kwid == 2) {
                                jmax = j;
                                for (j = 0;j < jmax;j++) {
                                    if (tmpstr2[j] <= ' ' || tmpstr2[j] > '~') tmpstr2[j] = '_';
                                }
                            }
                            else if (kwid == 3) {
                                jmax = j;
                                for (j = 0;j < jmax;j++) {
                                    if (tmpstr2[j] <= ' ' || tmpstr2[j] > '~') tmpstr2[j] = '.';
                                }
                            }
                            else if (kwid == 4) {
                                jmax = j;
                                j = 0;
                                while (j < jmax && tmpstr2[j] != ' ') {
                                    j++;
                                }
                                tmpstr2[j] = 0;
                            }
                            else if (kwid == 5) {
                                jmax = j;
                                j = 0;
                                while (j < jmax && tmpstr2[j] != ' ') {
                                    j++;
                                }
                                if (tmpstr2[j] == ' ')j++;
                                k = j;
                                while (j < jmax) {
                                    tmpstr2[j - k] = tmpstr2[j];
                                    j++;
                                }
                                tmpstr2[j - k] = 0;
                                j -= k;
                            }
                            break;
#ifdef _WIN32
                        case 6:
                            sprintf_s(tmpstr2, nalloc2, _LLD_, (long long)uid);
                            j = strlen(tmpstr2);
                            break;
                        case 7:
                            sprintf_s(tmpstr2, nalloc2, _LLD_, (long long)contig_length(contig));
                            j = strlen(tmpstr2);
                            break;
                        case 8:
                            sprintf_s(tmpstr2, nalloc2, "%f", (double)contig_avg_qscore(contig));
                            j = strlen(tmpstr2);
                            break;
                        case 9:
                            sprintf_s(tmpstr2, nalloc2, "%f", (double)contig_GC(contig));
                            j = strlen(tmpstr2);
                            break;
#else
                        case 6:
                            sprintf(tmpstr2, _LLD_, (long long)uid);
                            j = strlen(tmpstr2);
                            break;
                        case 7:
                            sprintf(tmpstr2, _LLD_, (long long)contig_length(contig));
                            j = strlen(tmpstr2);
                            break;
                        case 8:
                            sprintf(tmpstr2, "%f", (double)contig_avg_qscore(contig));
                            j = strlen(tmpstr2);
                            break;
                        case 9:
                            sprintf(tmpstr2, "%f", (double)contig_GC(contig));
                            j = strlen(tmpstr2);
                            break;
#endif
                        case 10:
                        case 11:
                        case 12:
                            j = strlen(contig_nameptr(contig));
                            if (j + 1 >= nalloc2) {
                                nalloc2 = j + 2;
                                tmpstr2 = realloc(tmpstr2, nalloc2);
                            }
                            memcpy(tmpstr2, contig_nameptr(contig), j);
                            tmpstr2[j] = 0;
                            jmax = j;
                            j = 0;
                            if (kwid > 10) {
                                while (tmpstr2[j] && tmpstr2[j] != '/')j++;
                                if (tmpstr2[j] == '/')j++;
                            }
                            if (kwid > 11) {
                                while (tmpstr2[j] && tmpstr2[j] != '/')j++;
                                if (tmpstr2[j] == '/')j++;
                            }
                            k = 0;
                            while (tmpstr2[j] && tmpstr2[j] != '/') {
                                tmpstr2[k] = tmpstr2[j];
                                j++;
                                k++;
                            }
                            tmpstr2[k] = 0;
                            j = k;
                            break;
                        default:
                            tmpstr2[0] = 0;
                            j = 0;
                            break;
                    }
                    tmpstr1 = tmpstr1 + 1 + sizeof(kwid);
                    if (nalloc3 <= i + j + 1) {
                        outputstr = realloc(outputstr, nalloc3 + j + i + 1);
                        memset(outputstr + i, 0, j + nalloc3 + 1);
                        nalloc3 = nalloc3 + j + i + 1;
                    }
                    memcpy(outputstr + i, tmpstr2, j);
                    i += j;
                    outputstr[i] = 0;
                }
                else {
                    if (nalloc3 <= i+1) {
                        outputstr = realloc(outputstr, nalloc3 + i + 1);
                        memset(outputstr + i, 0, nalloc3+1);
                        nalloc3 = nalloc3 + i + 1;
                    }
                    outputstr[i] = *tmpstr1;
                    i++;
                    tmpstr1 = tmpstr1 + 1;
                }
            }
            outputstr[i] = 0;
            contig_rename(contig, outputstr);
            contigs_to_fastxPF(fout, &contig, 1);
            contig_free(contig);
        }
        free(outputstr);
    }
    free(tmpstr2);
    free(fmtstr_parsed);

    PFclose(fin);
    PFclose(fout);

    return res;
}

int nucops_renameseq_b(args_t* args) {
    char* outputfn;
    PF_s* f;
    int res;
    
    outputfn = args_getstr(args, "output", 0, "stdout");
    res = 1;
    if (PFopen(&f, outputfn, "wb") == 0) {
        PFprintf(f, "Keywords should be surrounded by curly brackets (eg. {keyword})\n");
        PFprintf(f, "The following keywords are accepted:\n");
        PFprintf(f, " Sequence-specific\n");
        PFprintf(f, "   file       the name of the the input file\n");
        PFprintf(f, "   name       the original name of the sequence\n");
        PFprintf(f, "   _name_     the original name of the sequence with ' ' replaced by '_'\n");
        PFprintf(f, "   .name.     the original name of the sequence with ' ' replaced by '.'\n");
        PFprintf(f, "   prefix     the first word (until the first ' ') of original sequence name\n");
        PFprintf(f, "   suffix     the  original sequence name except the first word\n");
        PFprintf(f, "   id         the id of the sequence - ids are assiged to each sequence according\n");
        PFprintf(f, "              to the order in which they appear.\n");
        PFprintf(f, "   pacbio1    For pacbio reads - run id\n");
        PFprintf(f, "   pacbio2    For pacbio reads - read id\n");
        PFprintf(f, "   pacbio3    For pacbio reads - subread id\n");
        PFprintf(f, "   length     the length of the sequence\n");
        PFprintf(f, "   quality    the average quality of the sequence\n");
        PFprintf(f, "   gc         the gc content of the sequence\n\n");
        PFprintf(f, " Character keywords (with no alternative) :\n");
        PFprintf(f, "   o / open_curly               parsed as '{'\n");
        PFprintf(f, "   c / close_curly              parsed as '}'\n\n");
        PFprintf(f, " Character keywords (added as convience for use in terminals) :\n");
        PFprintf(f, "   sp / space                   parsed as ' '\n");
        PFprintf(f, "   q / quote                    parsed as a single quotation mark\n");
        PFprintf(f, "   qq / double_quote            parsed as '\"'\n");
        PFprintf(f, "   bq / back_quote              parsed as '`'\n");
        PFprintf(f, "   b / bang                     parsed as '!'\n");
        PFprintf(f, "   n / colon                    parsed as ':'\n");
        PFprintf(f, "   sn / semicolon               parsed as ';'\n");
        PFprintf(f, "   p / h / pound / hash         parsed as '#'\n");
        PFprintf(f, "   lt / less_than               parsed as '<'\n");
        PFprintf(f, "   gt / greater_than            parsed as '>'\n");
        PFprintf(f, "   s / slash / stroke           parsed as '/'\n");
        PFprintf(f, "   bs / backslash / backstroke  parsed as '\\'\n");
        PFprintf(f, "   d / dollar                   parsed as '$'\n");
        PFprintf(f, "   pc / percent                 parsed as '%'\n");
        res = 0;
    }
    PFclose(f);
    return res;
}

int nucops_renameseqs(int argc, char** argv) {
    int result;
    args_t* args;

    args = nucops_renameseq_init_args(argc, argv);
    result = 0;
    if (!args_ispresent(args, "help")) {
        if (args_ispresent(args, "keywords")) {
            result = nucops_renameseq_b(args);
        }
        else {
            result = nucops_renameseq_a(args);
        }
        if (result != 0) {
            args_report_error(args, "renameseq failed with code <%d>\n", result);
        }
    }
    args_free(args);
    return result;
}