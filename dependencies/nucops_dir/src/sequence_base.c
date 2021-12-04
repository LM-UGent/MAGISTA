#include <stdlib.h>
#include <string.h>
#include "osportstd.h"
#include "sequence_base.h"
#include "textparsing.h"

struct contig_s {
    char* nucleotides;
    float* coverage;
    size_t length;
    size_t nalloc;
    char* name;
    char* qscore;
    size_t namelen;
    int64_t uid;
    int fmt;
    char discared;
    char suffix_state;
};

struct read_s {
    char* nucleotides;
    char* quality;
    size_t length;
    size_t id;
};

contig_s* contig_alloc() {
    contig_s* result;
    result = (contig_s*)calloc(1, sizeof(contig_s));
    result->nalloc = 16;
    result->length = 0;
    result->nucleotides = (char*)malloc(result->nalloc);
    result->coverage = NULL;
    result->qscore = NULL;
    result->discared = 0;
    result->uid = -1;
    result->suffix_state = 0;
    return result;
}
contig_s* contig_cpyalloc(contig_s* src) {
    contig_s* result;
    result = (contig_s*)calloc(1, sizeof(contig_s));
    result->nalloc = src->nalloc;
    result->length = src->length;
    result->nucleotides = memcpyalloc(src->nucleotides,result->nalloc);
    if (src->qscore)
        result->qscore = memcpyalloc(src->qscore, result->nalloc);
    else
        result->qscore = NULL;
    if (src->coverage)
        result->coverage = memcpyalloc(src->coverage, result->nalloc*sizeof(float));
    else
        result->coverage = NULL;
    result->discared = src->discared;
    result->uid = src->uid;
    result->suffix_state = src->suffix_state;
    result->fmt = src->fmt;
    result->name = memcpyalloc(src->name, src->namelen + 1);
    result->namelen = src->namelen;
    return result;
}
void contig_free(contig_s* target) {
    if (target->name)free(target->name);
    if (target->nucleotides)free(target->nucleotides);
    if (target->coverage)free(target->coverage);
    if (target->qscore)free(target->qscore);
    free(target);
}

void local__contig_chr2bin(contig_s* target) {
    size_t i;
    for (i = 0;i < target->length;i++) {
        switch (target->nucleotides[i]) {
            case 'A':
            case 'a':
                target->nucleotides[i] = 0;
                break;
            case 'C':
            case 'c':
                target->nucleotides[i] = 1;
                break;
            case 'G':
            case 'g':
                target->nucleotides[i] = 2;
                break;
            case 'T':
            case 't':
            case 'U':
            case 'u':
                target->nucleotides[i] = 3;
                break;
            default:
                target->nucleotides[i] = -1;
                break;
        }
    }
    target->fmt = CONTIGFMT_BIN;
}
void local__contig_bin2chr(contig_s* target) {
    size_t i;
    char matches[5] = { 'N','A','C','G','T' };
    for (i = 0;i < target->length;i++) {
        target->nucleotides[i] = matches[target->nucleotides[i] + 1];
    }
    target->fmt = CONTIGFMT_CHR;
}

void contig_discard_data(contig_s* target) {
    uint64_t hash_value;
    if (target->fmt == CONTIGFMT_CHR) local__contig_chr2bin(target);
    if (target->fmt == CONTIGFMT_BIN) {
        if (target->nucleotides) {
            hash_value = data_hash(target->nucleotides, target->length);
            target->nucleotides = realloc(target->nucleotides, sizeof(uint64_t));
            memcpy(target->nucleotides, &hash_value, sizeof(uint64_t));
        }
        if (target->coverage) {
            hash_value = data_hash((char*)(target->coverage), target->length * sizeof(float));
            target->coverage = realloc(target->coverage, sizeof(uint64_t));
            memcpy(target->coverage, &hash_value, sizeof(uint64_t));
        }
        if (target->qscore) {
            hash_value = data_hash(target->qscore, target->length);
            target->qscore = realloc(target->qscore, sizeof(uint64_t));
            memcpy(target->qscore, &hash_value, sizeof(uint64_t));
        }
        target->discared = 1;
        target->fmt = CONTIGFMT_HASHED;
    }
}

contig_s* contig_reverse_complement(contig_s* source) {
    contig_s* result;
    size_t i, j;
    result = (contig_s*)calloc(1, sizeof(contig_s));
    result->nalloc = source->length;
    result->length = source->length;
    result->nucleotides = (char*)malloc(result->nalloc);
    result->fmt = source->fmt;
    if (source->fmt == CONTIGFMT_BIN) {
        j = source->length;
        for (i = 0;i < source->length;i++) {
            j--;
            result->nucleotides[j] = 3 - source->nucleotides[i];
        }
    }
    else if( source->fmt == CONTIGFMT_CHR){
        j = source->length;
        for (i = 0;i < source->length;i++) {
            j--;
            if(source->nucleotides[i]=='A') source->nucleotides[j]='T';
            else if(source->nucleotides[i] == 'a') source->nucleotides[j] = 't';
            else if (source->nucleotides[i] == 'C') source->nucleotides[j] = 'G';
            else if (source->nucleotides[i] == 'c') source->nucleotides[j] = 'g';
            else if (source->nucleotides[i] == 'G') source->nucleotides[j] = 'C';
            else if (source->nucleotides[i] == 'g') source->nucleotides[j] = 'c';
            else if (source->nucleotides[i] == 'T') source->nucleotides[j] = 'A';
            else if (source->nucleotides[i] == 't') source->nucleotides[j] = 'a';
        }
    }
    if (source->coverage) {
        result->coverage = (float*)malloc(result->nalloc * sizeof(float));
        j = source->length;
        for (i = 0;i < source->length;i++) {
            j--;
            result->coverage[j] = source->coverage[i];
        }
    }
    if (source->qscore) {
        result->qscore = (char*)malloc(result->nalloc);
        j = source->length;
        for (i = 0;i < source->length;i++) {
            j--;
            result->qscore[j] = source->qscore[i];
        }
    }
    return result;
}
contig_s* contig_from_str(char* src, char* name, int fmt) {
    contig_s* result;
    contig_s* tmp;
    size_t i, j;
    int is_rc;
    if (!src) return NULL;
    is_rc = 0;
    if (src[0] == '3' && src[1] == '\'') {
        is_rc = 1;
        src += 2;
    }
    else if (src[0] == '^') {
        is_rc = 1;
        src += 2;
    }
    if (src[0] == '5' && src[1] == '\'') {
        src += 2;
    }
    result = (contig_s*)calloc(1, sizeof(contig_s));
    result->nalloc = strlen(src);
    if (fmt != CONTIGFMT_NODATA) {
        result->nucleotides = (char*)malloc(result->nalloc);
    }
    else {
        result->nucleotides = NULL;
    }
    result->length = result->nalloc;

    j = 0;
    if (fmt == CONTIGFMT_BIN || fmt == CONTIGFMT_HASHED) {
        for (i = 0;i < result->length;i++) {
            if (src[i] > ' ') {
                switch (src[i]) {
                case 'A':
                case 'a':
                    result->nucleotides[j] = 0;
                    break;
                case 'C':
                case 'c':
                    result->nucleotides[j] = 1;
                    break;
                case 'G':
                case 'g':
                    result->nucleotides[j] = 2;
                    break;
                case 'T':
                case 't':
                case 'U':
                case 'u':
                    result->nucleotides[j] = 3;
                    break;
                default:
                    result->nucleotides[j] = -1;
                    break;
                }
                j++;
            }
        }
    }
    else if(fmt == CONTIGFMT_CHR){
        for (i = 0;i < result->length;i++) {
            if (src[i] > ' ' && src[i] <='~') {
                result->nucleotides[j] = src[i];
                j++;
            }
        }
    }
    else if (fmt == CONTIGFMT_NODATA) {
        for (i = 0;i < result->length;i++) {
            if (src[i] > ' ' && src[i] <= '~') {
                j++;
            }
        }
    }
    result->length = j;

    if (fmt == CONTIGFMT_HASHED) {
        contig_discard_data(result);
    }
    if (is_rc) {
        tmp = contig_reverse_complement(result);
        contig_free(result);
        result = tmp;
    }
    if (name) {
        result->namelen = strlen(name);
        result->name = (char*)malloc(result->namelen + 1);
        memcpy(result->name, name, result->namelen + 1);
    }
    result->fmt = fmt;
    return result;
}

int local__contig_get_toggle_suffix(contig_s* target) {
    size_t i;
    if (target->suffix_state == 0) {
        for (i = 0;i < target->namelen;i++) {
            if (target->name[i] == ' ')target->suffix_state = 2;
        }
        if (target->suffix_state == 0)target->suffix_state = 1;
    }
    if (target->suffix_state == 2) {
        return 0;
    }
    else if (target->suffix_state == 3) {
        return 1;
    }
    return 1;
}


char* contig_as_fastastr(contig_s* target, size_t linelen) {
    size_t curpos;
    char *result;
    size_t reslen;
    size_t i;
    int suffixtog;
    reslen = 2; /* a fasta sequence begins with '>' and covers at least two lines */
    if (target->name) reslen += target->namelen;
    else reslen += strlen("unnamed");
    reslen += target->length;
    if (linelen > 0) {
        reslen += target->length / linelen;
    }
    result = (char*)malloc(reslen + 1);
    result[0] = '>';
    if (target->name) {
        memcpy(result + 1, target->name, target->namelen);
        curpos = target->namelen + 1;
    }
    else {
        memcpy(result + 1, "unnamed", strlen("unnamed"));
        curpos = 1 + strlen("unnamed");
    }
    suffixtog = local__contig_get_toggle_suffix(target);
    contig_toggle_suffix(target, 1);
    result[curpos] = '\n';
    curpos++;
    if (target->fmt == CONTIGFMT_BIN) {
        for (i = 0;i < target->length;i++) {
            if (target->nucleotides[i] == 0)result[curpos] = 'A';
            else if (target->nucleotides[i] == 1)result[curpos] = 'C';
            else if (target->nucleotides[i] == 2)result[curpos] = 'G';
            else if (target->nucleotides[i] == 3)result[curpos] = 'T';
            else result[curpos] = 'N';
            curpos++;
            if (linelen > 0 && i % linelen == linelen - 1) {
                result[curpos] = '\n';
                curpos++;
            }
        }
    }
    else if (target->fmt == CONTIGFMT_CHR) {
        for (i = 0;i < target->length;i++) {
            result[curpos] = target->nucleotides[i];
            curpos++;
            if (linelen > 0 && i % linelen == linelen - 1) {
                result[curpos] = '\n';
                curpos++;
            }
        }
    }
    else {
        for (i = 0;i < target->length;i++) {
            result[curpos] = 'N';
            curpos++;
            if (linelen > 0 && i % linelen == linelen - 1) {
                result[curpos] = '\n';
                curpos++;
            }
        }
    }
    contig_toggle_suffix(target, suffixtog);
    result[reslen] = 0;
    return result;
}
char* contig_as_fastqstr(contig_s* target, size_t linelen) {
    size_t curpos;
    char *result;
    size_t reslen;
    size_t i;
    int suffixtog;
    reslen = 4; /* a fasta sequence begins with '>' and covers at least four lines */
    if (target->name) reslen += target->namelen;
    else reslen += strlen("unnamed");
    /* seq */
    reslen += target->length;
    if (linelen > 0) {
        reslen += target->length / linelen;
    }
    /* "+\n" and qscores */
    reslen += 2;
    reslen += target->length;
    if (linelen > 0) {
        reslen += target->length / linelen;
    }
    result = (char*)malloc(reslen + 1);
    result[0] = '>';
    if (target->name) {
        memcpy(result + 1, target->name, target->namelen);
        curpos = target->namelen + 1;
    }
    else {
        memcpy(result + 1, "unnamed", strlen("unnamed"));
        curpos = 1 + strlen("unnamed");
    }
    result[curpos] = '\n';
    curpos++;
    suffixtog = local__contig_get_toggle_suffix(target);
    contig_toggle_suffix(target, 1);
    if (target->fmt == CONTIGFMT_BIN) {
        for (i = 0;i < target->length;i++) {
            if (target->nucleotides[i] == 0)result[curpos] = 'A';
            else if (target->nucleotides[i] == 1)result[curpos] = 'C';
            else if (target->nucleotides[i] == 2)result[curpos] = 'G';
            else if (target->nucleotides[i] == 3)result[curpos] = 'T';
            else result[curpos] = 'N';
            curpos++;
            if (linelen > 0 && i % linelen == linelen - 1) {
                result[curpos] = '\n';
                curpos++;
            }
        }
    }
    else if (target->fmt == CONTIGFMT_CHR) {
        for (i = 0;i < target->length;i++) {
            result[curpos] = target->nucleotides[i];
            curpos++;
            if (linelen > 0 && i % linelen == linelen - 1) {
                result[curpos] = '\n';
                curpos++;
            }
        }
    }
    else {
        for (i = 0;i < target->length;i++) {
            result[curpos] = 'N';
            curpos++;
            if (linelen > 0 && i % linelen == linelen - 1) {
                result[curpos] = '\n';
                curpos++;
            }
        }
    }
    result[curpos] = '+';
    result[curpos+1] = '\n';
    curpos += 2;
    if (target->qscore){
        if (target->fmt == CONTIGFMT_BIN) {
            for (i = 0;i < target->length;i++) {
                result[curpos] = 33 + target->qscore[i];
                if (result[curpos] < 0) result[curpos] = 126;
                if (result[curpos] < 33) result[curpos] = 33;
                if (result[curpos] > 126) result[curpos] = 126;
                curpos++;
                if (linelen > 0 && i % linelen == linelen - 1) {
                    result[curpos] = '\n';
                    curpos++;
                }
            }
        }
        else if (target->fmt == CONTIGFMT_CHR) {
            for (i = 0;i < target->length;i++) {
                result[curpos] = target->qscore[i];
                curpos++;
                if (linelen > 0 && i % linelen == linelen - 1) {
                    result[curpos] = '\n';
                    curpos++;
                }
            }
        }
        else {
            for (i = 0;i < target->length;i++) {
                result[curpos] = 33;
                curpos++;
                if (linelen > 0 && i % linelen == linelen - 1) {
                    result[curpos] = '\n';
                    curpos++;
                }
            }
        }
    }
    else {
        for (i = 0;i < target->length;i++) {
            result[curpos] = 33;
            curpos++;
            if (linelen > 0 && i % linelen == linelen - 1) {
                result[curpos] = '\n';
                curpos++;
            }
        }
    }
    contig_toggle_suffix(target, suffixtog);
    result[reslen] = 0;
    return result;
}
char* contig_as_fastxstr(contig_s* target, size_t linelen) {
    size_t curpos;
    char *result;
    size_t reslen;
    size_t i;
    int suffixtog;
    reslen = 4; /* a fastq sequence begins with '>' and covers at least four lines */
    if (target->name) reslen += target->namelen;
    else reslen += strlen("unnamed");
    /* seq */
    reslen += target->length;
    if (linelen > 0) {
        reslen += target->length / linelen;
    }
    /* "\n+\n" and qscores */
    if(target->qscore){
        reslen += 3;
        reslen += target->length;
        if (linelen > 0) {
            reslen += target->length / linelen;
        }
    }
    suffixtog = local__contig_get_toggle_suffix(target);
    contig_toggle_suffix(target,1);
    result = (char*)malloc(reslen + 1);
    if (target->qscore) {
        result[0] = '@';
    }
    else {
        result[0] = '>';
    }
    if (target->name) {
        memcpy(result + 1, target->name, target->namelen);
        curpos = target->namelen + 1;
    }
    else {
        memcpy(result + 1, "unnamed", strlen("unnamed"));
        curpos = 1 + strlen("unnamed");
    }
    result[curpos] = '\n';
    curpos++;
    if (target->fmt == CONTIGFMT_BIN) {
        for (i = 0;i < target->length;i++) {
            if (target->nucleotides[i] == 0)result[curpos] = 'A';
            else if (target->nucleotides[i] == 1)result[curpos] = 'C';
            else if (target->nucleotides[i] == 2)result[curpos] = 'G';
            else if (target->nucleotides[i] == 3)result[curpos] = 'T';
            else result[curpos] = 'N';
            curpos++;
            if (linelen > 0 && i % linelen == linelen - 1) {
                result[curpos] = '\n';
                curpos++;
            }
        }
    }
    else if (target->fmt == CONTIGFMT_CHR) {
        for (i = 0;i < target->length;i++) {
            result[curpos] = target->nucleotides[i];
            curpos++;
            if (linelen > 0 && i % linelen == linelen - 1) {
                result[curpos] = '\n';
                curpos++;
            }
        }
    }
    else {
        for (i = 0;i < target->length;i++) {
            result[curpos] = 'N';
            curpos++;
            if (linelen > 0 && i % linelen == linelen - 1) {
                result[curpos] = '\n';
                curpos++;
            }
        }
    }
    if (target->qscore){
        result[curpos] = '\n';
        result[curpos+1] = '+';
        result[curpos+2] = '\n';
        curpos += 3;
        if (target->fmt == CONTIGFMT_BIN) {
            for (i = 0;i < target->length;i++) {
                result[curpos] = 33 + target->qscore[i];
                if (result[curpos] < 0) result[curpos] = 126;
                if (result[curpos] < 33) result[curpos] = 33;
                if (result[curpos] > 126) result[curpos] = 126;
                curpos++;
                if (linelen > 0 && i % linelen == linelen - 1) {
                    result[curpos] = '\n';
                    curpos++;
                }
            }
        }
        else if (target->fmt == CONTIGFMT_CHR) {
            for (i = 0;i < target->length;i++) {
                result[curpos] = target->qscore[i];
                curpos++;
                if (linelen > 0 && i % linelen == linelen - 1) {
                    result[curpos] = '\n';
                    curpos++;
                }
            }
        }
        else {
            for (i = 0;i < target->length;i++) {
                result[curpos] = 33;
                curpos++;
                if (linelen > 0 && i % linelen == linelen - 1) {
                    result[curpos] = '\n';
                    curpos++;
                }
            }
        }
    }
    contig_toggle_suffix(target, suffixtog);
    result[curpos] = 0;
    return result;
}


int contig_same_name(contig_s* A, contig_s* B){
    if(A->name && B->name){
        if (A->name[0] != 0 && B->name[0] != 0) {
            if (strcmp(A->name,B->name)!=0) return 0;
        }
        else {
            if (A->name[0] != 0 || B->name[0] != 0) return 0;
            if (*((uint64_t*)(A->name+1)) != *((uint64_t*)(B->name+1))) return 0;
        }
    }
    else {
        if (A->name && !(B->name)) return 0;
        if (B->name && !(A->name)) return 0;
    }
    return 1;
}

int contig_same_content(contig_s* A, contig_s* B){
    uint64_t hash_value_A;
    uint64_t hash_value_B;
    if (A->length != B->length) return 0;
    if (A->nucleotides && B->nucleotides){
        if(A->discared) hash_value_A = *((uint64_t*)(A->nucleotides));
        else hash_value_A = data_hash(A->nucleotides, A->length);
        if(B->discared) hash_value_B = *((uint64_t*)(B->nucleotides));
        else hash_value_B = data_hash(B->nucleotides, B->length);
        if (hash_value_A != hash_value_B) return 0;
    }
    else {
        if (A->nucleotides && !(B->nucleotides)) return 0;
        if (B->nucleotides && !(A->nucleotides)) return 0;
    }
    if (A->coverage && B->coverage){
        if(A->discared) hash_value_A = *((uint64_t*)(A->coverage));
        else hash_value_A = data_hash((char*)(A->coverage), A->length*sizeof(float));
        if(B->discared) hash_value_B = *((uint64_t*)(B->coverage));
        else hash_value_B = data_hash((char*)(B->coverage), B->length*sizeof(float));
        if (hash_value_A != hash_value_B) return 0;
    }
    else {
        if (A->coverage && !(B->coverage)) return 0;
        if (B->coverage && !(A->coverage)) return 0;
    }
    if (A->qscore && B->qscore){
        if(A->discared) hash_value_A = *((uint64_t*)(A->qscore));
        else hash_value_A = data_hash(A->qscore, A->length);
        if(B->discared) hash_value_B = *((uint64_t*)(B->qscore));
        else hash_value_B = data_hash(B->qscore, B->length);
        if (hash_value_A != hash_value_B) return 0;
    }
    else {
        if (A->qscore && !(B->qscore)) return 0;
        if (B->qscore && !(A->qscore)) return 0;
    }
    return 1;
}

int64_t contig_UID(contig_s* target){
    return target->uid;
}
void contig_setUID(contig_s* target, int64_t newuid){
    target->uid = newuid;
}

size_t contig_length(contig_s* target) {
    return target->length;
}
void contig_rename(contig_s* target, char* new_name) {
    if (target->name)free(target->name);
    if (new_name) {
        target->namelen = strlen(new_name);
        target->name = (char*)malloc(target->namelen + 1);
        memcpy(target->name, new_name, target->namelen + 1);
    }
    else {
        target->namelen = 0;
        target->name = NULL;
    }
}
void contig_toggle_suffix(contig_s* target, int truefalse) {
    size_t i;
    if (target->suffix_state == 0) {
        for (i = 0;i < target->namelen;i++) {
            if (target->name[i] == ' ')target->suffix_state = 2;
        }
        if (target->suffix_state == 0)target->suffix_state = 1;
    }
    if (target->suffix_state == 2 && truefalse == 0) {
        for (i = 0;i < target->namelen;i++) {
            if (target->name[i] == ' ')target->name[i] = '\0';
        }
        target->suffix_state = 3;
    }
    else if (target->suffix_state == 3 && truefalse == 1) {
        for (i = 0;i < target->namelen;i++) {
            if (target->name[i] == '\0')target->name[i] = ' ';
        }
        target->suffix_state = 2;
    }
}
void contig_remove_qscores(contig_s* target){
    if(target->qscore)free(target->qscore);
    target->qscore=NULL;
}
void contig_add_qscores(contig_s* target) {
    if (!target->qscore && (target->fmt == CONTIGFMT_BIN || target->fmt == CONTIGFMT_CHR)) {
        target->qscore = malloc(sizeof(char)*target->length);
    }
}
void contig_add_read(contig_s* target, int64_t first_position, int64_t post_position) {
    size_t i;
    size_t mini, maxi;
    int64_t tmp;
    if (!target->coverage) {
        target->coverage = calloc(target->length,sizeof(float));
    }
    if (first_position > post_position) {
        tmp = first_position;
        first_position = post_position;
        post_position = tmp;
    }
    if (first_position < 0) {
        first_position = 0;
    }
    if (post_position > (int64_t)(target->length)) {
        post_position = (int64_t)(target->length - 1);
    }
    mini = (size_t)first_position;
    maxi = (size_t)post_position;
    for (i = mini;i < maxi;i++) {
        target->coverage[i] += 1.0;
    }
}
size_t contig_covered_positions(contig_s* target) {
    size_t i;
    size_t result;
    if (!target->coverage)return 0;
    result = 0;
    for (i = 0;i < target->length;i++) {
        result += ((target->coverage[i] > 0) ? 1 : 0);
    }
    return result;
}
float contig_coverage(contig_s* target) {
    size_t i;
    double result;
    if (!target->coverage)return 0;
    result = 0;
    for (i = 0;i < target->length;i++) {
        result += target->coverage[i];
    }
    return ((float)result)/((float)(target->length));
}
float contig_covered_coverage(contig_s* target) {
    size_t i;
    size_t valid;
    double result;
    if (!target->coverage)return 0;
    result = 0;
    valid = 0;
    for (i = 0;i < target->length;i++) {
        if (target->coverage[i] > 0) {
            valid++;
            result += target->coverage[i];
        }
    }
    return ((float)result) / ((float)(valid));
}
float contig_position_coverage(contig_s* target, size_t pos) {
    if (pos >= target->length)return 0.0;
    if (target->coverage == NULL)return 0.0;
    return target->coverage[pos];
}

char* contig_nameptr(contig_s* target) {
    return target->name;
}
void contig_append_contig(contig_s* target, contig_s* extra) {
    size_t newlen;
    size_t i;
    newlen = target->length + extra->length;
    if (target->length == 0) {
        target->fmt = extra->fmt;
    }
    if (target->fmt != extra->fmt) {
        target->fmt = CONTIGFMT_NODATA;
        if (target->nucleotides) free(target->nucleotides);
        if (target->coverage) free(target->coverage);
        if (target->qscore) free(target->qscore);
        target->nucleotides = NULL;
        target->coverage = NULL;
        target->qscore = NULL;
        target->nalloc = 0;
    }
    else {
        if (target->fmt == CONTIGFMT_BIN || target->fmt == CONTIGFMT_CHR) {
            if (newlen > target->nalloc) {
                target->nalloc = (target->length + extra->length) * 2;
                target->nucleotides = (char*)realloc(target->nucleotides, target->nalloc);
                if (target->coverage) {
                    target->coverage = (float*)realloc(target->coverage, target->nalloc * sizeof(float));
                }
                if (target->qscore) {
                    target->qscore = (char*)realloc(target->qscore, target->nalloc);
                }
            }
            memcpy(target->nucleotides + target->length, extra->nucleotides, extra->length);
            if (target->coverage) {
                if (extra->coverage) {
                    memcpy(target->coverage + target->length * sizeof(float), extra->coverage, extra->length * sizeof(float));
                }
                else {
                    for (i = 0;i < extra->length;i++) {
                        target->coverage[target->length + i] = 1.0f;
                    }
                }
            }
            if (target->qscore) {
                if (extra->qscore) {
                    memcpy(target->coverage + target->length, extra->qscore, extra->length);
                }
                else {
                    for (i = 0;i < extra->length;i++) {
                        target->qscore[target->length + i] = -1;
                    }
                }
            }
        }
    }
    target->length = newlen;
}
contig_s* contig_subsection(contig_s* target, size_t start, size_t len) {
    contig_s* result;
    if (start + len > target->length) len = target->length - start;
    result = (contig_s*)calloc(1, sizeof(contig_s));
    result->nalloc = len;
    result->length = len;
    result->fmt = target->fmt;
    if (result->fmt == CONTIGFMT_BIN || result->fmt == CONTIGFMT_CHR) {
        result->nucleotides = memcpyalloc(target->nucleotides + start, len);
        if (target->coverage) result->coverage = memcpyalloc(target->coverage + start * sizeof(float), len * sizeof(float));
        if (target->qscore) result->qscore = memcpyalloc(target->qscore + start, len);
    }
    else {
        result->fmt = CONTIGFMT_NODATA;
        result->nucleotides = NULL;
        result->coverage = NULL;
        result->qscore = NULL;
        result->nalloc = 0;
    }
    return result;
}
void contig_overwrite(contig_s* target, int64_t offset, contig_s* src, size_t start, size_t len) {
    size_t i, maxi;
    if (target->fmt == src->fmt && (target->fmt == CONTIGFMT_BIN || target->fmt == CONTIGFMT_CHR)) {
        if (start + len > target->length) len = target->length - start;
        if (offset + len > target->nalloc) {
            target->nalloc = (size_t)(offset + (int64_t)len) * 2;
            target->nucleotides = (char*)realloc(target->nucleotides, target->nalloc);
            if (target->coverage) {
                target->coverage = (float*)realloc(target->coverage, target->nalloc * sizeof(float));
            }
            if (target->qscore) {
                target->qscore = (char*)realloc(target->qscore, target->nalloc);
            }
        }
        memcpy(target->nucleotides + offset, src->nucleotides + start, len);
        if (target->coverage) {
            if (src->coverage)
                memcpy(target->coverage + offset * sizeof(float), src->coverage + start * sizeof(float), len * sizeof(float));
            else {
                maxi = (size_t)(offset + (int64_t)len);
                for (i = target->length;i < maxi;i++) target->coverage[i] = 1.0f;
            }
        }
        if (target->qscore) {
            if (src->qscore)
                memcpy(target->qscore + offset, src->coverage + start, len);
            else {
                maxi = (size_t)(offset + (int64_t)len);
                for (i = target->length;i < maxi;i++) target->qscore[i] = -1;
            }
        }
    }
}
void contig_setnuc(contig_s* target, size_t pos, char new_nuc) {
    if (target->length > pos) {
        if (target->fmt = CONTIGFMT_BIN) {
            switch (new_nuc) {
            case 'A':
            case 'a':
                target->nucleotides[pos] = 0;
                break;
            case 'C':
            case 'c':
                target->nucleotides[pos] = 1;
                break;
            case 'G':
            case 'g':
                target->nucleotides[pos] = 2;
                break;
            case 'T':
            case 't':
                target->nucleotides[pos] = 3;
                break;
            default:
                target->nucleotides[pos] = -1;
                break;
            }
        }
        else {
            target->nucleotides[pos] = new_nuc;
        }

    }
}
char contig_getnuc(contig_s* target, size_t pos) {
    const char nuc[5] = { 'N','A','C','G','T' };
    if (!target || pos >= target->length)
        return 'N';
    if (target->fmt == CONTIGFMT_BIN) {
        if (target->nucleotides[pos] >= -1 && target->nucleotides[pos] < 4)
            return nuc[target->nucleotides[pos] + 1];
        else
            return 'N';
    }
    else
        return target->nucleotides[pos];
}
static inline int _contig_compare_positions(contig_s* target, size_t position1, size_t position2, size_t kmerlen) {
    size_t k;
    int delta;
    int result = 0;
    for (k = 0;k < kmerlen;k++) {
        if (target->nucleotides[position1 + k] == -1) return -2;
        if (target->nucleotides[position2 + k] == -1) return 2;
        else {
            delta = ((int)(target->nucleotides[position2 + k])) - ((int)(target->nucleotides[position1 + k]));
            if (delta > 0) return 1;
            if (delta < 0) return -1;
        }
    }
    return result;
}
size_t* contig_minimizer_positions(contig_s* target, size_t winsize, size_t kmerlen, size_t* p_nminimizers) {
    size_t curmin;
    size_t i, maxi, j, j0;
    size_t* result;
    size_t* tmp;
    size_t nchanges;
    int cmpres;

    if (target->fmt == CONTIGFMT_CHR)local__contig_chr2bin(target);
    if (target->fmt != CONTIGFMT_BIN) return NULL;

    maxi = target->length - kmerlen + 1;
    tmp = (size_t*)malloc(sizeof(size_t)*maxi);
    curmin = 0;
    for (i = 1;i < maxi;i++) {
        if (i - curmin >= winsize) {
            /* too far from the current minimizer - search for the next one */
            curmin++;
            j0 = curmin;
            for (j = 0; j < winsize;j++) {
                cmpres = _contig_compare_positions(target, curmin, j0 + j, kmerlen);
                if ( cmpres == -1 || cmpres == -2) {
                    curmin = j0 + j;
                }
            }
        }
        else {
            cmpres = _contig_compare_positions(target, curmin, i, kmerlen);
            if ( cmpres == -1 || cmpres == -2) {
                curmin = i;
            }
        }
        tmp[i] = curmin;
    }
    nchanges = 0;
    for (i = winsize;i < maxi;i++) {
        if (tmp[i] != tmp[i - 1])nchanges++;
    }
    result = (size_t*)malloc(sizeof(size_t)*(nchanges + 1));
    result[0] = tmp[winsize-1];
    j = 0;
    for (i = winsize;i < maxi;i++) {
        if (tmp[i] != tmp[i - 1]) {
            j++;
            result[j] = tmp[i];
        }
    }
    free(tmp);
    *p_nminimizers = nchanges + 1;
    return result;
}
size_t* contig_maximizer_positions(contig_s* target, size_t winsize, size_t kmerlen, size_t* p_nmaximizers) {
    size_t curmin;
    size_t i, maxi, j, j0;
    size_t* result;
    size_t* tmp;
    size_t nchanges;
    int cmpres;

    if (target->fmt == CONTIGFMT_CHR)local__contig_chr2bin(target);
    if (target->fmt != CONTIGFMT_BIN) return NULL;

    maxi = target->length - kmerlen + 1;
    tmp = (size_t*)malloc(sizeof(size_t)*maxi);
    curmin = 0;
    for (i = 1;i < maxi;i++) {
        if (i - curmin >= winsize) {
            /* too far from the current minimizer - search for the next one */
            curmin++;
            j0 = curmin;
            for (j = 0; j < winsize;j++) {
                cmpres = _contig_compare_positions(target, curmin, j0 + j, kmerlen);
                if (cmpres == 1 || cmpres == -2) {
                    curmin = j0 + j;
                }
            }
        }
        else {
            cmpres = _contig_compare_positions(target, curmin, i, kmerlen);
            if (cmpres == 1 || cmpres == -2) {
                curmin = i;
            }
        }
        tmp[i] = curmin;
    }
    nchanges = 0;
    for (i = winsize;i < maxi;i++) {
        if (tmp[i] != tmp[i - 1])nchanges++;
    }
    result = (size_t*)malloc(sizeof(size_t)*(nchanges + 1));
    result[0] = tmp[winsize - 1];
    j = 0;
    for (i = winsize;i < maxi;i++) {
        if (tmp[i] != tmp[i - 1]) {
            j++;
            result[j] = tmp[i];
        }
    }
    free(tmp);
    *p_nmaximizers = nchanges + 1;
    return result;
}
contig_s* contig_kmer_sequence_from_positions(contig_s* target, size_t* positions, size_t npositions, size_t kmerlen) {
    contig_s* result;
    size_t i,j;
    
    if (target->fmt == CONTIGFMT_CHR)local__contig_chr2bin(target);
    if (target->fmt != CONTIGFMT_BIN) return NULL;

    if (!target || !positions)return NULL;

    result = contig_alloc();
    result->nucleotides = (char*)malloc(kmerlen*npositions);
    result->length = kmerlen*npositions;
    result->nalloc = kmerlen*npositions;
    if (target->name) {
        result->name = text_join(&(target->name), "", "", "_kmers", 1);
        result->namelen = strlen(result->name);
    }
    else {
        result->name = NULL;
        result->namelen = 0;
    }
    result->fmt = target->fmt;
    for (i = 0;i < npositions;i++) {
        if (positions[i] + kmerlen <= target->length) {
            memcpy(result->nucleotides + i*kmerlen, target->nucleotides + positions[i], kmerlen);
        }
        else {
            memset(result->nucleotides + i*kmerlen, -1, kmerlen);
        }
    }
    if (target->coverage) {
        result->coverage = (float*)malloc(sizeof(float)*kmerlen*npositions);
        memcpy(result->nucleotides + i*kmerlen, target->nucleotides + positions[i], kmerlen);
        for (i = 0;i < npositions;i++) {
            if (positions[i] + kmerlen <= target->length) {
                memcpy(result->coverage + (i*kmerlen), target->coverage + positions[i], kmerlen*sizeof(float));
            }
            else {
                for (j = 0;j < kmerlen;j++) {
                    result->coverage[i*kmerlen + j] = 0.0;
                }
            }
        }
    }
    return result;
}

void contig_countACGTN(contig_s* target, size_t* p_A, size_t* p_C, size_t* p_G, size_t* p_T, size_t* p_N) {
    size_t pos;
    size_t counts[5];
    size_t v;
    counts[0] = counts[1] = counts[2] = counts[3] = counts[4] = 0;
    if (target->fmt == CONTIGFMT_BIN) {
        for (pos = 0; pos < target->length; pos++) {
            v = target->nucleotides[pos] + 1;
            if (v < 5)counts[v]++;
            else counts[0]++;
        }
    }
    else if (target->fmt == CONTIGFMT_CHR) {
        for (pos = 0;pos < target->length;pos++) {
            switch (target->nucleotides[pos]) {
            case 'A':
            case 'a':
                counts[1]++;
                break;
            case 'C':
            case 'c':
                counts[2]++;
                break;
            case 'G':
            case 'g':
                counts[3]++;
                break;
            case 'T':
            case 't':
            case 'U':
            case 'u':
                counts[4]++;
                break;
            default:
                counts[0]++;
                break;
            }
        }
    }
    else {
        counts[0] = target->length;
    }
    *p_N = counts[0];
    *p_A = counts[1];
    *p_C = counts[2];
    *p_G = counts[3];
    *p_T = counts[4];
}
float contig_GC(contig_s* target) {
    size_t GC;
    size_t pos;
    GC = 0;
    if (target->fmt == CONTIGFMT_BIN) {
        for (pos = 0; pos < target->length; pos++) {
            if (target->nucleotides[pos] == 2 || target->nucleotides[pos] == 3)GC++;;
        }
    }
    else if (target->fmt == CONTIGFMT_CHR) {
        for (pos = 0;pos < target->length;pos++) {
            switch (target->nucleotides[pos]) {
            case 'C':
            case 'c':
            case 'G':
            case 'g':
                GC++;
                break;
            default:
                break;
            }
        }
    }
    return ((float)GC) / ((float)target->length);
}

char* contig_qscores(contig_s* target){
    return target->qscore;
}
int contig_min_qscore(contig_s* target){
    size_t i;
    int res;
    res = 256;
    if(target->qscore && (target->fmt==CONTIGFMT_BIN || target->fmt==CONTIGFMT_CHR)) {
        for (i = 0;i < target->length;i++) {
            if ((int)(target->qscore[i]) < res) res = (int)(target->qscore[i]);
        }
        if (target->fmt == CONTIGFMT_CHR) {
            res -= 33;
        }
    }
    else return -1;
    return res;
}
float contig_avg_qscore(contig_s* target){
    size_t i;
    double res;
    res = 0;
    if (target->qscore && (target->fmt == CONTIGFMT_BIN || target->fmt == CONTIGFMT_CHR)) {
        if (target->fmt == CONTIGFMT_BIN) {
            for (i = 0;i < target->length;i++) {
                res += (double)((int)(target->qscore[i]));
            }
        }
        else if (target->fmt == CONTIGFMT_CHR) {
            for (i = 0;i < target->length;i++) {
                res += (double)((int)(target->qscore[i])-33);
            }
        }
        return (float)(res / ((double)(target->length)));
    }
    else return -1;
    return (float)res;
}
int contig_max_qscore(contig_s* target){
    size_t i;
    int res;
    res = -1;
    if (target->qscore && (target->fmt == CONTIGFMT_BIN || target->fmt == CONTIGFMT_CHR)) {
        for (i = 0;i < target->length;i++) {
            if ((int)(target->qscore[i]) > res) res = (int)(target->qscore[i]);
        }
        if (target->fmt == CONTIGFMT_CHR) res -= 33;
    }
    else return -1;
    return res;
}

contig_s* contig_from_fastaPF_fmt(PF_t* src, int fmt) {
    char* line;
    contig_s* curcontig;
    contig_s* tmpcontig;
    int64_t loc;
    curcontig = NULL;
    loc = PFtell(src);
    line = PFreadline(src);
    if (line && line[0] == '>') {
        curcontig = contig_alloc();
        contig_rename(curcontig, line + 1);
        free(line);
        loc = PFtell(src);
        line = PFreadline(src);
        while (line && line[0] != '>') {
            tmpcontig = contig_from_str(line, NULL, fmt);
            contig_append_contig(curcontig, tmpcontig);
            contig_free(tmpcontig);
            free(line);

            loc = PFtell(src);
            line = PFreadline(src);
        }
        free(line);
        PFseek(src, loc, SEEK_SET);
    }
    return curcontig;
}
contig_s* contig_from_fastaPF(PF_t* src) {
    return contig_from_fastaPF_fmt(src, CONTIGFMT_CHR);
}

contig_s** contigs_from_fastaPF_limited(PF_t* f, size_t* p_amount, size_t limit, long long* p_resume) {
    char* line;
    contig_s** result;
    contig_s* curcontig;
    contig_s* tmpcontig;
    size_t nres, nallocres, curalloc;
    long long filepos;
    result = NULL;
    nres = 0;
    nallocres = 0;
    filepos = 0;
    if (p_resume) {
        PFseek(f,*p_resume,SEEK_SET);
        filepos = *p_resume;
        *p_resume=0;
    }
    nallocres = 16;
    result = (contig_s**)malloc(sizeof(contig_s*)*nallocres);
    curcontig = NULL;
    curalloc = 0;
    while (line = PFreadline(f)) {
        curalloc += strlen(line);
        if(limit > 0 && curalloc > limit) {
            free(line);
            break;
        }
        if (line[0] == '>') {
            if(p_resume) *p_resume = filepos;
            curcontig = contig_alloc();
            if (nres == nallocres) {
                nallocres *= 2;
                result = (contig_s**)realloc(result, sizeof(contig_s*)*nallocres);
            }
            contig_rename(curcontig, line + 1);
            result[nres] = curcontig;
            nres++;
        }
        else if (curcontig) {
            tmpcontig = contig_from_str(line, NULL, CONTIGFMT_CHR);
            contig_append_contig(curcontig, tmpcontig);
            contig_free(tmpcontig);
        }
        free(line);
        filepos = PFtell(f);
    }
    if(!line && p_resume) *p_resume = 0;
    *p_amount = nres;
    result = (contig_s**)realloc(result, sizeof(contig_s*)*nres);
    return result;
}

contig_s** contigs_from_fasta_limited(char* filename, size_t* p_amount, size_t limit, long long* p_resume) {
    PF_t* f;
    contig_s** result;
    result = NULL;
    f = NULL;
    if (PFopen(&f, filename, "rb") == 0) {
        result = contigs_from_fastaPF_limited(f, p_amount, limit, p_resume);
    }
    PFclose(f);
    return result;
}
contig_s** contigs_from_fasta(char* filename, size_t* p_amount) {
    return contigs_from_fasta_limited(filename, p_amount, 0, NULL);
}

contig_s* contig_from_fastxPF_fmt(PF_t* src, int fmt) {
    char* line;
    contig_s* curcontig;
    contig_s* tmpcontig;
    int64_t loc;
    size_t maxi, ql;
    curcontig = NULL;
    loc = PFtell(src);
    line = PFreadline(src);
    while (line && line[0] != '>' && line[0] != '@') {
        free(line);
        line = PFreadline(src);
    }
    if (line && line[0] == '>') {
        curcontig = contig_alloc();
        contig_rename(curcontig, line + 1);
        free(line);
        loc = PFtell(src);
        line = PFreadline(src);
        while (line && line[0] != '>') {
            tmpcontig = contig_from_str(line, NULL, fmt);
            contig_append_contig(curcontig, tmpcontig);
            contig_free(tmpcontig);
            free(line);

            loc = PFtell(src);
            line = PFreadline(src);
        }
        free(line);
        PFseek(src, loc, SEEK_SET);
    }
    else if (line && line[0] == '@') {
        curcontig = contig_alloc();
        contig_rename(curcontig, line + 1);
        free(line);
        line = PFreadline(src);
        while (line && line[0] != '+') {
            tmpcontig = contig_from_str(line, NULL, fmt);
            contig_append_contig(curcontig, tmpcontig);
            contig_free(tmpcontig);
            free(line);

            loc = PFtell(src);
            line = PFreadline(src);
        }
        free(line);
        ql = 0;
        contig_add_qscores(curcontig);
        line = PFreadline(src);
        while (line && ql < curcontig->length) {
            maxi = strlen(line);
            if (line[maxi - 1] == '\n')maxi--;
            if (line[maxi - 1] == '\r')maxi--;
            if (maxi + ql > curcontig->length) {
                maxi = curcontig->length - ql;
            }
            memcpy(curcontig->qscore + ql, line, maxi);
            ql += maxi;
            free(line);
            if (ql == curcontig->length) {
                line = NULL;
            }
            else {
                line = PFreadline(src);
            }
        }
        if (line)free(line);
    }
    return curcontig;
}
contig_s* contig_from_fastxPF(PF_t* src) {
    return contig_from_fastxPF_fmt(src, CONTIGFMT_CHR);
}

contig_s** contigs_from_fastxPF_limited(PF_t* f, size_t* p_amount, size_t limit, long long* p_resume) {
    char* line;
    contig_s** result;
    contig_s* curcontig;
    contig_s* tmpcontig;
    size_t nres, nallocres, nqcs, linelen, i, curalloc;
    long long filepos;
    int state;
    result = NULL;
    nres = 0;
    nallocres = 0;
    state = 0;
    nqcs = 0;
    filepos = 0;
    if (p_resume) {
        PFseek(f,*p_resume,SEEK_SET);
        filepos = *p_resume;
        *p_resume = 0;
    }
    nallocres = 16;
    result = (contig_s**)malloc(sizeof(contig_s*)*nallocres);
    curcontig = NULL;
    curalloc = 0;
    while (line = PFreadline(f)) {
        curalloc += strlen(line);
        if(limit > 0 && curalloc > limit) {
            free(line);
            break;
        }
        if ( (line[0] == '@' || line[0] == '>')  && state == 0) {
            if(p_resume) *p_resume = filepos;
            curcontig = contig_alloc();
            if (nres == nallocres) {
                nallocres *= 2;
                result = (contig_s**)realloc(result, sizeof(contig_s*)*nallocres);
            }
            contig_rename(curcontig, line + 1);
            result[nres] = curcontig;
            nres++;
        }
        else if (curcontig) {
            if(state==0){
                if(line[0]=='+') {
                    state=1;
                    nqcs = 0;
                    curcontig->qscore = malloc(curcontig->nalloc);
                }
                else {
                    tmpcontig = contig_from_str(line, NULL, CONTIGFMT_CHR);
                    contig_append_contig(curcontig, tmpcontig);
                    contig_free(tmpcontig);
                }
            }
            else{
               linelen = strlen(line);
               if (nqcs + linelen >= curcontig->length) {
                   linelen = curcontig->length - nqcs;
                   state = 0;
               }
               memcpy(curcontig->qscore+nqcs,line,linelen);
               i = nqcs;
               nqcs += linelen;
               if(state==0) nqcs = curcontig->length;
               /*
               while(i < nqcs){
                   curcontig->qscore[i] -= 33;
                   i++;
               }
               */
            }
        }
        free(line);
        filepos = PFtell(f);
    }
    if(!line && p_resume) *p_resume = 0;
    if (p_resume && *p_resume != 0) {
        nres--;
        contig_free(result[nres]);
    }
    *p_amount = nres;
    result = (contig_s**)realloc(result, sizeof(contig_s*)*nres);
    return result;
}
contig_s** contigs_from_fastx_limited(char* filename, size_t* p_amount, size_t limit, long long* p_resume) {
    PF_t* f;
    contig_s** result;
    result = NULL;
    f = NULL;
    if (PFopen(&f, filename, "rb") == 0) {
        result = contigs_from_fastxPF_limited(f, p_amount, limit, p_resume);
    }
    else {
        if (p_amount)*p_amount = 0;
    }
    PFclose(f);
    return result;
}
contig_s** contigs_from_fastx(char* filename, size_t* p_amount) {
    return contigs_from_fastx_limited(filename, p_amount, 0, NULL);
}

int contig_positions_are_identical(contig_s* A, contig_s* B, size_t pos1, size_t pos2) {
    if (A->fmt == CONTIGFMT_HASHED || A->fmt == CONTIGFMT_NODATA)return 0;
    if (B->fmt == CONTIGFMT_HASHED || B->fmt == CONTIGFMT_NODATA)return 0;
    if (A->fmt == B->fmt) return (A->nucleotides[pos1] == B->nucleotides[pos2]);
    if (A->fmt == CONTIGFMT_BIN && B->fmt == CONTIGFMT_CHR) {
        if (A->nucleotides[pos1] == 0) {
            if (B->nucleotides[pos2] == 'A' || B->nucleotides[pos2] == 'a') return 1;
            else return 0;
        }
        if (A->nucleotides[pos1] == 1) {
            if (B->nucleotides[pos2] == 'C' || B->nucleotides[pos2] == 'c') return 1;
            else return 0;
        }
        if (A->nucleotides[pos1] == 2) {
            if (B->nucleotides[pos2] == 'G' || B->nucleotides[pos2] == 'g') return 1;
            else return 0;
        }
        if (A->nucleotides[pos1] == 3) {
            if (B->nucleotides[pos2] == 'T' || B->nucleotides[pos2] == 't' || B->nucleotides[pos2] == 'U' || B->nucleotides[pos2] == 'u') return 1;
            else return 0;
        }
    }
    if (B->fmt == CONTIGFMT_BIN && A->fmt == CONTIGFMT_CHR) {
        if (B->nucleotides[pos2] == 0) {
            if (A->nucleotides[pos1] == 'A' || A->nucleotides[pos1] == 'a') return 1;
            else return 0;
        }
        if (B->nucleotides[pos2] == 1) {
            if (A->nucleotides[pos1] == 'C' || A->nucleotides[pos1] == 'c') return 1;
            else return 0;
        }
        if (B->nucleotides[pos2] == 2) {
            if (A->nucleotides[pos1] == 'G' || A->nucleotides[pos1] == 'g') return 1;
            else return 0;
        }
        if (B->nucleotides[pos2] == 3) {
            if (A->nucleotides[pos1] == 'T' || A->nucleotides[pos1] == 't' || A->nucleotides[pos1] == 'U' || A->nucleotides[pos1] == 'u') return 1;
            else return 0;
        }
    }
    return 0;
}
int contig_positions_are_complementary(contig_s* A, contig_s* B, size_t pos1, size_t pos2) {
    if (A->fmt == CONTIGFMT_HASHED || A->fmt == CONTIGFMT_NODATA)return 0;
    if (B->fmt == CONTIGFMT_HASHED || B->fmt == CONTIGFMT_NODATA)return 0;
    if (A->fmt == CONTIGFMT_BIN && B->fmt == CONTIGFMT_BIN) return ((A->nucleotides[pos1]) == 3-(B->nucleotides[pos2]));
    if (A->fmt == CONTIGFMT_CHR && B->fmt == CONTIGFMT_CHR) {
        if (A->nucleotides[pos1] == 'T' || A->nucleotides[pos1] == 't' || A->nucleotides[pos1] == 'U' || A->nucleotides[pos1] == 'u') {
            if (B->nucleotides[pos2] == 'A' || B->nucleotides[pos2] == 'a') return 1;
            else return 0;
        }
        if (A->nucleotides[pos1] == 'G' || A->nucleotides[pos1] == 'g') {
            if (B->nucleotides[pos2] == 'C' || B->nucleotides[pos2] == 'c') return 1;
            else return 0;
        }
        if (A->nucleotides[pos1] == 'C' || A->nucleotides[pos1] == 'c') {
            if (B->nucleotides[pos2] == 'G' || B->nucleotides[pos2] == 'g') return 1;
            else return 0;
        }
        if (A->nucleotides[pos1] == 'A' || A->nucleotides[pos1] == 'a') {
            if (B->nucleotides[pos2] == 'T' || B->nucleotides[pos2] == 't' || B->nucleotides[pos2] == 'U' || B->nucleotides[pos2] == 'u') return 1;
            else return 0;
        }
    }
    if (A->fmt == CONTIGFMT_BIN && B->fmt == CONTIGFMT_CHR) {
        if (A->nucleotides[pos1] == 3) {
            if (B->nucleotides[pos2] == 'A' || B->nucleotides[pos2] == 'a') return 1;
            else return 0;
        }
        if (A->nucleotides[pos1] == 2) {
            if (B->nucleotides[pos2] == 'C' || B->nucleotides[pos2] == 'c') return 1;
            else return 0;
        }
        if (A->nucleotides[pos1] == 1) {
            if (B->nucleotides[pos2] == 'G' || B->nucleotides[pos2] == 'g') return 1;
            else return 0;
        }
        if (A->nucleotides[pos1] == 0) {
            if (B->nucleotides[pos2] == 'T' || B->nucleotides[pos2] == 't' || B->nucleotides[pos2] == 'U' || B->nucleotides[pos2] == 'u') return 1;
            else return 0;
        }
    }
    if (B->fmt == CONTIGFMT_BIN && A->fmt == CONTIGFMT_CHR) {
        if (B->nucleotides[pos2] == 3) {
            if (A->nucleotides[pos1] == 'A' || A->nucleotides[pos1] == 'a') return 1;
            else return 0;
        }
        if (B->nucleotides[pos2] == 2) {
            if (A->nucleotides[pos1] == 'C' || A->nucleotides[pos1] == 'c') return 1;
            else return 0;
        }
        if (B->nucleotides[pos2] == 1) {
            if (A->nucleotides[pos1] == 'G' || A->nucleotides[pos1] == 'g') return 1;
            else return 0;
        }
        if (B->nucleotides[pos2] == 0) {
            if (A->nucleotides[pos1] == 'T' || A->nucleotides[pos1] == 't' || A->nucleotides[pos1] == 'U' || A->nucleotides[pos1] == 'u') return 1;
            else return 0;
        }
    }
    return 0;
}

double contig_check_ANI_CIGAR(contig_s* query, contig_s* reference, char* cigar, size_t refstart, int query_is_rc, size_t* p_M, size_t* p_I, size_t* p_D, size_t* p_R) {
    size_t lasti;
    size_t i, i0, k;
    size_t cigarstrlen;
    char septype;
    size_t seplen;
    size_t M, I, D, R;
    size_t qstart;
    
    cigarstrlen = strlen(cigar);
    i0 = 0;
    
    M = I = D = R = 0;
    qstart = 0;
    i = 0;
    while (i < cigarstrlen) {
        lasti = i;
        while (cigar[i] && cigar[i] >= '0' && cigar[i] <= '9')i++;
        if (cigar[i] == 0)break;
        septype = cigar[i];
        cigar[i] = 0;
        seplen = (size_t)atoi(cigar + lasti);
        cigar[i] = septype;
        if (septype == '=') {
            seplen = contig_length(query);
            septype = 'M';
        }
        if (septype == '*') {
            break;
        }
        for (k = 0;k < seplen;k++) {
            switch (septype) {
                case 'M':
                    if (qstart >= query->length || refstart >= reference->length)R++;
                    else {
                        if (contig_positions_are_identical(query, reference, qstart, refstart))M++;
                        else R++;
                    }
                    refstart++;
                    qstart++;
                    break;
                case 'I':
                    I++;
                case 'S':
                    qstart++;
                    break;
                case 'D':
                    D++;
                case 'N':
                    refstart++;
                    break;
                default:
                    break;
            }
        }
        i++;
    }
    if (p_M)*p_M = M;
    if (p_I)*p_I = I;
    if (p_D)*p_D = D;
    if (p_R)*p_R = R;
    return ((double)M) / ((double)(M + I + R + D));
}
void contig_check_monerrors_CIGAR(contig_s* query, contig_s* reference, char* cigar, size_t refstart, int query_is_rc, monoerr_s* p_error_profile, int reset_profile) {
    size_t lasti;
    size_t i, i0, k;
    size_t cigarstrlen;
    char septype;
    char refnuc, qnuc;
    size_t seplen;
    size_t qstart;
    int f;
    static size_t M = 0, I = 0, D = 0, R = 0;

    if (!p_error_profile) return;
    cigarstrlen = strlen(cigar);
    i0 = 0;

    qstart = 0;
    i = 0;
    f = 0;
    if(reset_profile)
        memset(p_error_profile, 0, sizeof(monoerr_s));
    while (i < cigarstrlen) {
        lasti = i;
        while (cigar[i] && cigar[i] >= '0' && cigar[i] <= '9')i++;
        if (cigar[i] == 0)break;
        septype = cigar[i];
        cigar[i] = 0;
        seplen = (size_t)atoi(cigar + lasti);
        cigar[i] = septype;
        if (septype == '=') {
            seplen = contig_length(query);
            septype = 'M';
        }
        if (septype == '*') {
            break;
        }
        for (k = 0;k < seplen;k++) {
            switch (septype) {
            case 'M':
                if (qstart >= query->length || refstart >= reference->length)R++;
                else {
                    if (contig_positions_are_identical(query, reference, qstart, refstart)) {
                        M++;
                    }
                    else {
                        R++;
                        refnuc = contig_getnuc(reference, refstart);
                        if (refnuc >= 'a') refnuc -= 'a' - 'A';
                        qnuc = contig_getnuc(query, qstart);
                        if (qnuc >= 'a') qnuc -= 'a' - 'A';
                        if (refnuc == 'A') {
                            if (qnuc == 'C')  p_error_profile->AtoC++;
                            else if (qnuc == 'G') p_error_profile->AtoG++;
                            else if (qnuc == 'T') p_error_profile->AtoT++;
                            else if (qnuc == 'N') p_error_profile->badA++;
                            else {
                                f++;
                            }
                        }
                        else if (refnuc == 'C') {
                            if (qnuc == 'A')  p_error_profile->CtoA++;
                            else if (qnuc == 'G') p_error_profile->CtoG++;
                            else if (qnuc == 'T') p_error_profile->CtoT++;
                            else if (qnuc == 'N') p_error_profile->badC++;
                            else {
                                f++;
                            }
                        }
                        else if (refnuc == 'G') {
                            if (qnuc == 'A')  p_error_profile->GtoA++;
                            else if (qnuc == 'C') p_error_profile->GtoC++;
                            else if (qnuc == 'T') p_error_profile->GtoT++;
                            else if (qnuc == 'N') p_error_profile->badG++;
                            else {
                                f++;
                            }
                        }
                        else if (refnuc == 'T') {
                            if (qnuc == 'A')  p_error_profile->TtoA++;
                            else if (qnuc == 'C') p_error_profile->TtoC++;
                            else if (qnuc == 'G') p_error_profile->TtoG++;
                            else if (qnuc == 'N') p_error_profile->badT++;
                            else {
                                f++;
                            }
                        }
                        else {
                            f++;
                        }
                    }
                }
                refstart++;
                qstart++;
                break;
            case 'I':
                I++;
                qnuc = contig_getnuc(query, qstart);
                if (qnuc >= 'a') qnuc -= 'a' - 'A';
                if (qnuc == 'A') p_error_profile->insA++;
                else if (qnuc == 'C') p_error_profile->insC++;
                else if (qnuc == 'G') p_error_profile->insG++;
                else if (qnuc == 'T') p_error_profile->insT++;
                else {
                    f++;
                }
            case 'S':
                qstart++;
                break;
            case 'D':
                D++;
                refnuc = contig_getnuc(reference, refstart);
                if (refnuc >= 'a') refnuc -= 'a' - 'A';
                if (refnuc == 'A') p_error_profile->delA++;
                else if (refnuc == 'C') p_error_profile->delC++;
                else if (refnuc == 'G') p_error_profile->delG++;
                else if (refnuc == 'T') p_error_profile->delT++;
            case 'N':
                refstart++;
                break;
            default:
                break;
            }
        }
        i++;
    }
}

void contig_check_qualest_CIGAR(contig_s* query, contig_s* reference, char* cigar, size_t refstart, int query_is_rc, size_t* counts, int minQ, int maxQ) {
    size_t lasti;
    size_t i, i0, k;
    size_t cigarstrlen;
    char septype;
    size_t seplen;
    size_t M, I, D, R;
    size_t qstart;

    cigarstrlen = strlen(cigar);
    i0 = 0;

    M = I = D = R = 0;
    qstart = 0;
    i = 0;
    if (maxQ < minQ || maxQ>128 || minQ<0)return;
    for (i = 0;i < (size_t)(maxQ - minQ);i++) {
        counts[i - (size_t)minQ] = 0;
    }
    if (query->qscore == NULL)return;
    while (i < cigarstrlen) {
        lasti = i;
        while (cigar[i] && cigar[i] >= '0' && cigar[i] <= '9')i++;
        if (cigar[i] == 0)break;
        septype = cigar[i];
        cigar[i] = 0;
        seplen = (size_t)atoi(cigar + lasti);
        cigar[i] = septype;
        if (septype == '=') {
            seplen = contig_length(query);
            septype = 'M';
        }
        if (septype == '*') {
            break;
        }
        for (k = 0;k < seplen;k++) {
            switch (septype) {
            case 'M':
                if (qstart >= query->length || refstart >= reference->length)R++;
                else {
                    if (contig_positions_are_identical(query, reference, qstart, refstart))M++;
                    else {
                        counts[query->qscore[qstart] - minQ]++;
                        R++;
                    }
                }
                refstart++;
                qstart++;
                break;
            case 'I':
                counts[query->qscore[refstart] - minQ]++;
                I++;
            case 'S':
                qstart++;
                break;
            case 'D':
                if (qstart + 1 < query->length) {
                    counts[query->qscore[qstart+1] - minQ]++;
                }
                D++;
            case 'N':
                refstart++;
                break;
            default:
                break;
            }
        }
        i++;
    }
}


contig_s* contigs_concatenate(contig_s** contigs, size_t amount) {
    contig_s* result;
    size_t i;
    int has_coverage;
    result = (contig_s*)calloc(1, sizeof(contig_s));
    result->length = 0;
    result->nalloc = 0;

    has_coverage = 0;
    for (i = 0; i < amount; i++) {
        result->nalloc += contigs[i]->length;
        if (contigs[i]->coverage) has_coverage = 1;
    }
    result->length = 0;
    result->nucleotides = (char*)malloc(result->nalloc);
    if (has_coverage) {
        result->coverage = (float*)calloc(result->nalloc, sizeof(float));
    }
    for (i = 0; i < amount; i++) {
        contig_append_contig(result, contigs[i]);
    }

    return result;
}
void contigs_to_fastaPF(PF_t* file, contig_s** contigs, size_t amount) {
    size_t i;
    char* tmpstr;
    for (i = 0;i < amount;i++) {
        tmpstr = contig_as_fastastr(contigs[i], 0);
        PFprintf(file, "%s\n", tmpstr);
        free(tmpstr);
    }
}
void contigs_to_fastqPF(PF_t* file, contig_s** contigs, size_t amount) {
    size_t i;
    char* tmpstr;
    for (i = 0;i < amount;i++) {
        tmpstr = contig_as_fastqstr(contigs[i], 0);
        PFprintf(file, "%s\n", tmpstr);
        free(tmpstr);
    }
}
void contigs_to_fastxPF(PF_t* file, contig_s** contigs, size_t amount) {
    size_t i;
    char* tmpstr;
    for (i = 0;i < amount;i++) {
        tmpstr = contig_as_fastxstr(contigs[i], 0);
        PFprintf(file, "%s\n", tmpstr);
        free(tmpstr);
    }
}

void contigs_to_fasta(char* filename, contig_s** contigs, size_t amount) {
    PF_t* f;
    f = NULL;
    if (PFopen(&f, filename, "wb") == 0) {
        contigs_to_fastaPF(f, contigs, amount);
    }
    PFclose(f);
}
void contigs_to_fastq(char* filename, contig_s** contigs, size_t amount) {
    PF_t* f;
    f = NULL;
    if (PFopen(&f, filename, "wb") == 0) {
        contigs_to_fastqPF(f, contigs, amount);
    }
    PFclose(f);
}
void contigs_to_fastx(char* filename, contig_s** contigs, size_t amount) {
    PF_t* f;
    f = NULL;
    if (PFopen(&f, filename, "wb") == 0) {
        contigs_to_fastxPF(f, contigs, amount);
    }
    PFclose(f);
}

read_s* read_alloc() {
    read_s* result;
    result = (read_s*)calloc(1, sizeof(read_s));
    return result;
}

/* TODO
struct contig_alignment_s {
    contig_s* reference;
    contig_s* consensus;
    contig_s** unaligned_queue;
    size_t* queue_id;
    char* queue_free_flag;
    size_t nunaligned;
    size_t nalloc_unaligned;
    char** aligned;
    size_t* alnid;
    size_t naligned;
    size_t nalloc_aligned;
    size_t lastid;
};
typedef struct contig_alignment_s conaln_s;

conaln_s* conaln_alloc() { return NULL; }
void conaln_free(conaln_s* target){
    if (target->reference)free(target->reference);
}
size_t conaln_add(conaln_s* target, contig_s* contig){}
void conaln_addref(conaln_s* target, contig_s* contig){
    target->reference = contig_cpyalloc(contig);
}
void conaln_consensus(conaln_s* target, contig_s* contig){}
void conaln_consensus_as_ref(conaln_s* target) {}
void local__contig_global_alignment(contig_s* ref, contig_s* query) {
    char* block;
    size_t p;
    size_t pos_extra;
    pos_extra = 
    while()
}
void conaln_aln2ref(conaln_s* target) {

}
void conaln_stat(conaln_s* target, char* stat){}
contig_s* get_consensus(conaln_s* target){}
*/


