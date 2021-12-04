
#ifndef DEF_SEQBASE
#define DEF_SEQBASE

#include <stdint.h>
#include "osportstd.h"

#define CONTIGFMT_BIN   0
#define CONTIGFMT_CHR   1
#define CONTIGFMT_HASHED    2
#define CONTIGFMT_NODATA    3

typedef struct contig_s contig_s;
typedef struct read_s read_s;

contig_s* contig_alloc();
void contig_free(contig_s* target);

contig_s* contig_from_str(char* src, char* name, int fmt);
contig_s* contig_from_fastaPF(PF_t* src);
contig_s* contig_from_fastxPF(PF_t* src);
char* contig_as_fastastr(contig_s* target, size_t linelen);
char* contig_as_fastqstr(contig_s* target, size_t linelen);
char* contig_as_fastxstr(contig_s* target, size_t linelen);
void contig_discard_data(contig_s* target);
int contig_same_name(contig_s* A, contig_s* B);
int contig_same_content(contig_s* A, contig_s* B);
int64_t contig_UID(contig_s* target);
void contig_setUID(contig_s* target, int64_t newuid);
size_t contig_length(contig_s* target);
void contig_rename(contig_s* target, char* new_name);
void contig_add_read(contig_s* target, int64_t first_position, int64_t post_position);
size_t contig_covered_positions(contig_s* target);
float contig_coverage(contig_s* target);
float contig_covered_coverage(contig_s* target);
float contig_position_coverage(contig_s* target, size_t pos);
void contig_toggle_suffix(contig_s* target, int truefalse);
void contig_remove_qscores(contig_s* target);
char* contig_nameptr(contig_s* target);
void contig_append_contig(contig_s* target, contig_s* extra);
contig_s* contig_subsection(contig_s* target, size_t start, size_t len);
void contig_overwrite(contig_s* target, int64_t offset, contig_s* src, size_t start, size_t len);
void contig_setnuc(contig_s* target, size_t pos, char new_nuc);
char contig_getnuc(contig_s* target, size_t pos);
size_t* contig_minimizer_positions(contig_s* target, size_t winsize, size_t kmerlen, size_t* p_nminimizers);
size_t* contig_maximizer_positions(contig_s* target, size_t winsize, size_t kmerlen, size_t* p_nmaximizers);
contig_s* contig_kmer_sequence_from_positions(contig_s* target, size_t* positions, size_t npositions, size_t kmerlen);
void contig_countACGTN(contig_s* target, size_t* p_A, size_t* p_C, size_t* p_G, size_t* p_T, size_t* p_N);
float contig_GC(contig_s* target);
char* contig_qscores(contig_s* target);
int contig_min_qscore(contig_s* target);
float contig_avg_qscore(contig_s* target);
int contig_max_qscore(contig_s* target);


contig_s** contigs_from_fasta(char* filename, size_t* p_amount);
contig_s** contigs_from_fastx(char* filename, size_t* p_amount);
contig_s** contigs_from_fasta_limited(char* filename, size_t* p_amount, size_t limit, long long* p_resume);
contig_s** contigs_from_fastx_limited(char* filename, size_t* p_amount, size_t limit, long long* p_resume);
contig_s** contigs_from_fastaPF_limited(PF_t* f, size_t* p_amount, size_t limit, long long* p_resume);
contig_s** contigs_from_fastxPF_limited(PF_t* f, size_t* p_amount, size_t limit, long long* p_resume);

int contig_positions_are_identical(contig_s* A, contig_s* B, size_t pos1, size_t pos2);
int contig_positions_are_complementary(contig_s* A, contig_s* B, size_t pos1, size_t pos2);
double contig_check_ANI_CIGAR(contig_s* query, contig_s* reference, char* cigar, size_t refstart, int query_is_rc, size_t* p_M, size_t* p_I, size_t* p_D, size_t* p_R);
typedef struct monoerr_s {
    size_t AtoC, AtoG, AtoT, delA;
    size_t CtoA, CtoG, CtoT, delC;
    size_t GtoA, GtoC, GtoT, delG;
    size_t TtoA, TtoC, TtoG, delT;
    size_t insA, insC, insG, insT;
    size_t badA, badC, badG, badT;
} monoerr_s;
void contig_check_monerrors_CIGAR(contig_s* query, contig_s* reference, char* cigar, size_t refstart, int query_is_rc, monoerr_s* p_error_profile, int reset_profile);
void contig_check_qualest_CIGAR(contig_s* query, contig_s* reference, char* cigar, size_t refstart, int query_is_rc, size_t* counts, int minQ, int maxQ);


contig_s* contigs_concatenate(contig_s** contigs, size_t amount);
void contigs_to_fasta(char* filename, contig_s** contigs, size_t amount);
void contigs_to_fastaPF(PF_t* file, contig_s** contigs, size_t amount);
void contigs_to_fastq(char* filename, contig_s** contigs, size_t amount);
void contigs_to_fastqPF(PF_t* file, contig_s** contigs, size_t amount);
void contigs_to_fastx(char* filename, contig_s** contigs, size_t amount);
void contigs_to_fastxPF(PF_t* file, contig_s** contigs, size_t amount);



#endif
