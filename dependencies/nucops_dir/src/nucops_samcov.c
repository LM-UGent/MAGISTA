#include "osportstd.h"
#include "bnlib_fileops.h"
#include "sequence_base.h"
#include "argparser.h"
#include "datamap.h"
#include <string.h>

int SAMFLAG_multi(int x) { return ((x & 1) != 0); }
int SAMFLAG_aligned(int x) { return ((x & 2) != 0); }
int SAMFLAG_unmapped(int x) { return ((x & 4) != 0); }
int SAMFLAG_next_unmapped(int x) { return ((x & 8) != 0); }
int SAMFLAG_is_revcomp(int x) { return ((x & 16) != 0); }
int SAMFLAG_next_is_revcomp(int x) { return ((x & 32) != 0); }
int SAMFLAG_is_first(int x) { return ((x & 64) != 0); }
int SAMFLAG_is_last(int x) { return ((x & 128) != 0); }
int SAMFLAG_is_secondary(int x) { return ((x & 256) != 0); }
int SAMFLAG_is_bad(int x) { return ((x & 512) != 0); }
int SAMFLAG_is_PCRdup(int x) { return ((x & 1024) != 0); }
int SAMFLAG_is_supplementary(int x) { return ((x & 2048) != 0); }

size_t CIGAR_length_bowtie(char* cigar) {
    int sep_count;
    size_t lasti;
    size_t i;
    size_t cigarstrlen;
    size_t result;
    char* septypes;
    int* seplens;
    cigarstrlen = strlen(cigar);
    sep_count = 0;
    for (i = 0;i < cigarstrlen;i++) {
        if (cigar[i] == 'M' || cigar[i] == 'I' || cigar[i] == 'D' || cigar[i] == 'N' || cigar[i] == 'S' || cigar[i] == 'H' || cigar[i] == 'P')
            sep_count++;
    }
    septypes = malloc(sep_count);
    seplens = malloc(sizeof(int)*sep_count);
    sep_count = 0;
    lasti = 0;
    for (i = 0;i < cigarstrlen;i++) {
        if (cigar[i] == 'M' || cigar[i] == 'I' || cigar[i] == 'D' || cigar[i] == 'N' || cigar[i] == 'S' || cigar[i] == 'H' || cigar[i] == 'P'){
            septypes[sep_count] = cigar[i];
            cigar[i] = 0;
            seplens[sep_count] = atoi(cigar+lasti);
            cigar[i] = septypes[sep_count];
            sep_count++;
        }
        if (cigar[i]<'0' || cigar[i]>'9')lasti = i + 1;
    }
    result = 0;
    for (i = 0;i < sep_count;i++) {
        /*under normal circumstances, I shouldn't increase the reference position, but it does with bowtie alignments */
        if (septypes[i] == 'M' || septypes[i] == 'I' || septypes[i] == 'D' || septypes[i] == 'N') {
            result += seplens[i];
        }
    }
    free(septypes);
    free(seplens);
    return result;
}

size_t CIGAR_length_normal(char* cigar) {
    int sep_count;
    size_t lasti;
    size_t i;
    size_t cigarstrlen;
    size_t result;
    char* septypes;
    int* seplens;
    cigarstrlen = strlen(cigar);
    sep_count = 0;
    for (i = 0;i < cigarstrlen;i++) {
        if (cigar[i] == 'M' || cigar[i] == 'I' || cigar[i] == 'D' || cigar[i] == 'N' || cigar[i] == 'S' || cigar[i] == 'H' || cigar[i] == 'P')
            sep_count++;
    }
    septypes = malloc(sep_count);
    seplens = malloc(sizeof(int)*sep_count);
    sep_count = 0;
    lasti = 0;
    for (i = 0;i < cigarstrlen;i++) {
        if (cigar[i] == 'M' || cigar[i] == 'I' || cigar[i] == 'D' || cigar[i] == 'N' || cigar[i] == 'S' || cigar[i] == 'H' || cigar[i] == 'P') {
            septypes[sep_count] = cigar[i];
            cigar[i] = 0;
            seplens[sep_count] = atoi(cigar + lasti);
            cigar[i] = septypes[sep_count];
            sep_count++;
        }
        if (cigar[i]<'0' || cigar[i]>'9')lasti = i + 1;
    }
    result = 0;
    for (i = 0;i < sep_count;i++) {
        if (septypes[i] == 'M' || septypes[i] == 'D' || septypes[i] == 'N') {
            result += seplens[i];
        }
    }
    free(septypes);
    free(seplens);
    return result;
}

args_t* nucops_samcov_init_args(int argc, char** argv) {
    args_t* result;
    result = args_alloc();

    args_add(result, "input", 'i', "str,str");
    args_add(result, "sam", 's', "str");
    args_add(result, "ref", 'r', "str");
    args_add(result, "output", 'o', "str");
    args_add(result, "bowtie", 'b', "");
    args_add(result, "detailed", 'd', "");
    args_add(result, "separator", ' ', "str");

    args_add_help(result, "input", "INPUT FILE(S)", "To be used instead of or in conjunction with --ref and --sam. If the --sam flag is present, then The first value is the name of the sam file and the second is the name of the reference. Otherwise, it becomes an alias for --ref.", "Specify input file(s)");
    args_add_help(result, "sam", "SAM FILE", "Specify the name of the SAM file to be analysed", "Specify the SAM file");
    args_add_help(result, "ref", "REFERENCE FILE", "Specify the name of the reference file to be analysed", "Specify the reference file");
    args_add_help(result, "output", "OUTPUT FILE", "Specify the name of the output file to be analysed", "Specify the output file");
    args_add_help(result, "bowtie", "BOWTIE-COMPATIBLE MODE", "Bowtie doesn't produce files which conform to the SAM specification. In this mode, it is not possible to compute error rates, but what is computed should correspond to reality.", "bowtie-compatible-mode");
    args_add_help(result, "detailed", "DETAILED ERROR PROFILE", "Output the detailed error profile. Currently limited to outputting all SNPs.", "Output detailed error profile");
    args_add_help(result, "separator", "CHOICE OF SEPARATOR", "specify which format should be used. Currently,  only 'tsv' and 'csv' are valid inputs", "Specify output format ('tsv' or 'csv')");
    args_parse(result, argc, argv);
    return result;
}

int nucops_samcov_a(args_t* args) {
    char* sfile;
    char* rfile;
    char* ofile;
    char* line;
    PF_s* samfile;
    PF_s* output;
    PF_s* reffile;
    contig_s** contigs;
    contig_s* query_seq;
    contig_s* reference_seq;
    size_t i,j,j0,ncontigs;
    int unparsedi,unparsedu;
    int err;
    DM64_t* dm;
    int btmode;
    int nf;
    
    char* refseqname;
    char* cigar;
    char* seqptr;
    char* tmpstr;
    char separator;
    int flags;
    int64_t seqstart, seqlen;
    int64_t TLEN, ntemplates, PNEXT;
    int64_t cigarlen;
    
    size_t Msum,Tsum;
    int64_t readrange[2];
    size_t totlen, totaln;
    size_t M, I, D, R, Mtot, Itot, Dtot, Rtot, MIDRtot;
    size_t missing_count, missing_length;
    monoerr_s error_profile;
    int isdetailedprofile;

    unparsedu = unparsedi = 0;
    memset(&error_profile, 0, sizeof(monoerr_s));

    separator = '\t';
    if (args_ispresent(args, "separator")) {
        tmpstr = args_getstr(args, "separator", 0, NULL);
        if (!tmpstr) {
            args_report_error(args, "please specify an output format after --separator\n");
            return 1;
        }
        if (strcmp(tmpstr, "tsv")==0)separator = '\t';
        else if (strcmp(tmpstr, "csv")==0)separator = ',';
        else {
            args_report_error(args, "%s is not a valid file format\n", tmpstr);
            return 1;
        }
    }
    /* Figure out which input file is the reference and which is the sam file. */
    sfile = args_getstr(args, "sam", 0, NULL);
    if (!sfile) {
        args_getstr(args, "input", unparsedi, NULL);
        if (sfile) unparsedi++;
    }
    rfile = args_getstr(args, "ref", 0, NULL);
    if (!rfile) {
        args_getstr(args, "input", unparsedi, NULL);
        if (rfile) unparsedi++;
    }
    if (!sfile) {
        args_getstr(args, NULL, unparsedu, NULL);
        if (sfile) unparsedu++;
    }
    if (!rfile) {
        args_getstr(args, NULL, unparsedu, NULL);
        if (rfile) unparsedu++;
    }
    if (unparsedu == 1 && unparsedi == 0 && !rfile && sfile) {
        /* special case: When niether --input nor --sam nor --ref are specified, and there is only one positional argument, that argument is the reference */
        rfile = sfile;
        sfile = "stdin";
    }
    if (!sfile) {
        if (!rfile) {
            args_report_error(args, "Please specify a SAM and REF file (one may be stdin)\n");
            return 1;
        }
        sfile = "stdin";
    }
    if (!rfile) {
        rfile = "stdin";
    }

    if (strcmp(sfile, rfile) == 0) {
        args_report_error(args, "SAM and REF files have the same name, aborting\n");
        return 1;
    }
    err = 0;

    btmode = args_ispresent(args, "bowtie");
    isdetailedprofile = args_ispresent(args, "detailed");

    ofile = args_getstr(args, "output", 0, args_getstr(args, NULL, unparsedu, "stdout"));

    args_report_info(args, "Finished parsing arguments\n");
    args_report_info(args, " SAM:%s\n", sfile);
    args_report_info(args, " REF:%s\n", rfile);
    args_report_info(args, " OUT:%s\n", ofile);
    if (btmode) {
        args_report_info(args, " Bowtie-suited mode is enabled\n");
    }

    if (PFopen(&samfile, sfile, "rb") != 0) {
        args_report_error(args, "Failed to open SAM file (%s)\n", sfile);
        if (!err)err = 2;
    }
    if (PFopen(&output, ofile, "wb") != 0) {
        args_report_error(args, "Failed to open OUT file (%s)\n", ofile);
        if (!err)err = 3;
    }
    if (PFopen(&reffile, rfile, "rb") != 0) {
        args_report_error(args, "Failed to open SAM file (%s)\n", rfile);
        if (!err)err = 4;
    }
    contigs = NULL;
    ncontigs = 0;
    if (!err) {
        args_report_progress(args, "Reading reference\n");
        contigs = contigs_from_fastxPF_limited(reffile, &ncontigs, 0, NULL);
    }
    if (!contigs || ncontigs == 0) {
        args_report_error(args, "Failed read sequences from REF file (%s)\n", rfile);
        if (!err)err = 5;
    }
    dm = new_DM64(0, 0);
    if (!err) {
        missing_count = 0;
        missing_length = 0;
        args_report_progress(args, "Registering sequences\n");
        for (i = 0;i < ncontigs;i++) {
            contig_toggle_suffix(contigs[i], 0);
            DM64_append(dm, contig_nameptr(contigs[i]), (int)strlen(contig_nameptr(contigs[i])), (int64_t)i);
        }
        Msum = 0;
        Tsum = 0;
        ntemplates = 0;
        i = 0;
        args_report_progress(args, "Parsing SAM\n", (int64_t)i);
        M = I = D = R = Mtot = Itot = Dtot = Rtot = 0;
        while (line = PFreadline(samfile)) {
            if (line[0] != '@') {
                j0 = 0;
                j = 0;
                while (line[j] && line[j] != '\t')j++;
                if (line[j]) { line[j] = 0; j++; }
                /* ignore read name for now*/

                j0 = j; while (line[j] && line[j] != '\t')j++;
                if (line[j]) { line[j] = 0; j++; }
                flags = atoi(line + j0);

                j0 = j; while (line[j] && line[j] != '\t')j++;
                if (line[j]) { line[j] = 0; j++; }
                refseqname = line + j0;

                j0 = j; while (line[j] && line[j] != '\t')j++;
                if (line[j]) { line[j] = 0; j++; }
                seqstart = atoll(line + j0) - 1;

                j0 = j; while (line[j] && line[j] != '\t')j++;
                if (line[j]) { line[j] = 0; j++; }
                /* ignore match quality for now */

                j0 = j; while (line[j] && line[j] != '\t')j++;
                if (line[j]) { line[j] = 0; j++; }
                cigar = line + j0;

                j0 = j; while (line[j] && line[j] != '\t')j++;
                if (line[j]) { line[j] = 0; j++; }
                /* ignore the rnext field */

                j0 = j; while (line[j] && line[j] != '\t')j++;
                if (line[j]) { line[j] = 0; j++; }
                PNEXT = atoll(line + j0) - 1;

                j0 = j; while (line[j] && line[j] != '\t')j++;
                if (line[j]) { line[j] = 0; j++; }
                TLEN = atoll(line + j0);

                j0 = j; while (line[j] && line[j] != '\t')j++;
                if (line[j]) { line[j] = 0; j++; }
                seqptr = line + j0;
                seqlen = strlen(line + j0);

                if (cigar[0] != '*') {
                    reference_seq = contigs[DM64_get(dm, refseqname, (int)strlen(refseqname), &nf)];
                    query_seq = NULL;
                    if (btmode) {
                        cigarlen = CIGAR_length_bowtie(cigar);
                    }
                    else {
                        cigarlen = CIGAR_length_normal(cigar);
                        if (seqptr[0] == '*') {
                            missing_count++;
                            missing_length += cigarlen;
                        }
                        else {
                            query_seq = contig_from_str(seqptr, "query", CONTIGFMT_BIN);
                            contig_check_ANI_CIGAR(query_seq, reference_seq, cigar, seqstart, SAMFLAG_is_revcomp(flags), &M, &I, &D, &R);
                            Mtot += M;
                            Itot += I;
                            Dtot += D;
                            Rtot += R;
                            if (isdetailedprofile) {
                                contig_check_monerrors_CIGAR(query_seq, reference_seq, cigar, seqstart, SAMFLAG_is_revcomp(flags), &error_profile, 0);
                            }
                        }
                    }
                    if (TLEN == 0)TLEN = cigarlen;

                    Msum += cigarlen;
                    if (SAMFLAG_is_first(flags) || !SAMFLAG_multi(flags)) {
                        Tsum += (TLEN >= 0) ? TLEN : -TLEN;
                        if (!SAMFLAG_is_secondary(flags) && !SAMFLAG_is_supplementary(flags))
                            ntemplates++;
                    }

                    if (SAMFLAG_is_revcomp(flags) && btmode ) {
                        readrange[0] = seqstart - cigarlen + 1;
                        readrange[1] = seqstart + 1;
                    }
                    else {
                        readrange[0] = seqstart;
                        readrange[1] = seqstart + cigarlen;
                    }
                    contig_add_read(reference_seq, readrange[0], readrange[1]);
                    if (query_seq)contig_free(query_seq);
                }
                i++;
                if (i % 1000 == 0) {
                    args_report_progress(args, "SAM contained " _LLD_ " sequences\r", (int64_t)i);
                }
            }
            free(line);
        }
        args_report_progress(args, "SAM contained " _LLD_ " sequences\n", (int64_t)i);
        if (missing_count > 0) {
            args_report_warning(NULL, _LLD_ " sequences, accounting for " _LLD_ " alignment length were missing and ignored for computing error rates\n", missing_count, missing_length);
        }
        PFprintf(output, "Average insert length%c%f\n", separator, ((double)Tsum)/((double)ntemplates));
        totlen = 0;
        totaln = 0;
        args_report_progress(args, "Evaluating coverage\n", (int64_t)i);
        for (i = 0;i < ncontigs;i++) {
            totlen += contig_length(contigs[i]);
            totaln += contig_covered_positions(contigs[i]);
        }
        PFprintf(output, "Genome fraction%c%f\n", separator, ((double)totaln) / ((double)totlen));
        PFprintf(output, "Genome coverage%c%f\n", separator, ((double)Tsum) / ((double)totlen));
        MIDRtot = Rtot + Mtot + Itot + Dtot;
        if (MIDRtot > 0) {
            PFprintf(output, "Identity rate%c%f\n", separator, ((double)Mtot) / ((double)MIDRtot));
            PFprintf(output, "Mismatch rate%c%f\n", separator, ((double)Rtot) / ((double)MIDRtot));
            PFprintf(output, "Insertion rate%c%f\n", separator, ((double)Itot) / ((double)MIDRtot));
            PFprintf(output, "Deletion rate%c%f\n", separator, ((double)Dtot) / ((double)MIDRtot));
            PFprintf(output, "Matches%c" _LLD_ "\n", separator, Mtot);
            PFprintf(output, "Mismatches%c" _LLD_ "\n", separator, Rtot);
            PFprintf(output, "Insertions%c" _LLD_ "\n", separator, Itot);
            PFprintf(output, "Deletions%c" _LLD_ "\n", separator, Dtot);
        }
        if (isdetailedprofile) {
            PFprintf(output, "A->C%c" _LLD_ "\n", separator, error_profile.AtoC);
            PFprintf(output, "A->G%c" _LLD_ "\n", separator, error_profile.AtoG);
            PFprintf(output, "A->T%c" _LLD_ "\n", separator, error_profile.AtoT);
            PFprintf(output, "A->N%c" _LLD_ "\n", separator, error_profile.badA);
            PFprintf(output, "C->A%c" _LLD_ "\n", separator, error_profile.CtoA);
            PFprintf(output, "C->G%c" _LLD_ "\n", separator, error_profile.CtoG);
            PFprintf(output, "C->T%c" _LLD_ "\n", separator, error_profile.CtoT);
            PFprintf(output, "C->N%c" _LLD_ "\n", separator, error_profile.badC);
            PFprintf(output, "G->A%c" _LLD_ "\n", separator, error_profile.GtoA);
            PFprintf(output, "G->C%c" _LLD_ "\n", separator, error_profile.GtoC);
            PFprintf(output, "G->T%c" _LLD_ "\n", separator, error_profile.GtoT);
            PFprintf(output, "G->N%c" _LLD_ "\n", separator, error_profile.badG);
            PFprintf(output, "T->A%c" _LLD_ "\n", separator, error_profile.TtoA);
            PFprintf(output, "T->C%c" _LLD_ "\n", separator, error_profile.TtoC);
            PFprintf(output, "T->G%c" _LLD_ "\n", separator, error_profile.TtoG);
            PFprintf(output, "T->N%c" _LLD_ "\n", separator, error_profile.badT);

            PFprintf(output, "del.A%c" _LLD_ "\n", separator, error_profile.delA);
            PFprintf(output, "del.C%c" _LLD_ "\n", separator, error_profile.delC);
            PFprintf(output, "del.G%c" _LLD_ "\n", separator, error_profile.delG);
            PFprintf(output, "del.T%c" _LLD_ "\n", separator, error_profile.delT);
            PFprintf(output, "ins.A%c" _LLD_ "\n", separator, error_profile.insA);
            PFprintf(output, "ins.C%c" _LLD_ "\n", separator, error_profile.insC);
            PFprintf(output, "ins.G%c" _LLD_ "\n", separator, error_profile.insG);
            PFprintf(output, "ins.T%c" _LLD_ "\n", separator, error_profile.insT);
        }
    }
    for (i = 0;i < ncontigs;i++) {
        contig_free(contigs[i]);
    }
    if (ncontigs > 0)free(contigs);
    free_DM64(dm);
    PFclose(samfile);
    PFclose(output);
    return (int)err;
}
            
 
int nucops_samcov(int argc, char** argv) {
    int result;
    args_t* args;

    args = nucops_samcov_init_args(argc, argv);
    result = 0;
    if (!args_ispresent(args, "help")) {
        result = nucops_samcov_a(args);
        if (result != 0) {
            args_report_error(args, "samcov failed with code <%d>\n", result);
        }
    }
    args_free(args);
    return result;
}