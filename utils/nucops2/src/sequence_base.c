#include <stdlib.h>
#include <string.h>
#include "osportstd.h"
#include "sequence_base.h"
#include "textparsing.h"

#define CONTIGFMT_BIN   0
#define CONTIGFMT_CHR   1

struct contig_s {
    char* nucleotides;
    float* coverage;
    size_t length;
    size_t nalloc;
    char* name;
    char* qscore;
    size_t namelen;
    int fmt;
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
    result->qscore = NULL;
    return result;
}
void contig_free(contig_s* target) {
    if (target->nucleotides)free(target->nucleotides);
    if (target->coverage)free(target->coverage);
    if (target->qscore)free(target->qscore);
    free(target);
}
contig_s* contig_reverse_complement(contig_s* source) {
    contig_s* result;
    size_t i, j;
    result = (contig_s*)calloc(1, sizeof(contig_s));
    result->nalloc = source->length;
    result->length = source->length;
    result->nucleotides = (char*)malloc(result->nalloc);
    if (source->fmt == CONTIGFMT_BIN) {
        j = source->length;
        for (i = 0;i < source->length;i++) {
            j--;
            result->nucleotides[j] = 3 - source->nucleotides[i];
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
contig_s* contig_from_str(char* src, char* name) {
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
    result->nucleotides = (char*)malloc(result->nalloc);
    result->length = result->nalloc;

    j = 0;
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
    return result;
}

char* contig_as_fastastr(contig_s* target, size_t linelen) {
    size_t curpos;
    char *result;
    size_t reslen;
    size_t i;
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
    result[curpos] = '\n';
    curpos++;
    if (target->fmt == 0) {
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
    result[reslen] = 0;
    return result;
}
char* contig_as_fastqstr(contig_s* target, size_t linelen) {
    size_t curpos;
    char *result;
    size_t reslen;
    size_t i;
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
    if (target->fmt == 0) {
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
    result[curpos] = '+';
    result[curpos+1] = '\n';
    curpos += 2;
    if (target->qscore){
        for (i = 0;i < target->length;i++) {
            result[curpos] = 33 + target->qscore[i];
            if (result[curpos]<0) result[curpos]=126;
            if (result[curpos]<33) result[curpos]=33;
            if (result[curpos]>126) result[curpos]=126;
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
    result[reslen] = 0;
    return result;
}
char* contig_as_fastxstr(contig_s* target, size_t linelen) {
    size_t curpos;
    char *result;
    size_t reslen;
    size_t i;
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
    if (target->fmt == 0) {
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
    if (target->qscore){
        result[curpos] = '\n';
        result[curpos+1] = '+';
        result[curpos+2] = '\n';
        curpos += 3;
        for (i = 0;i < target->length;i++) {
            result[curpos] = 33 + target->qscore[i];
            if (result[curpos]<0) result[curpos]=126;
            if (result[curpos]<33) result[curpos]=33;
            if (result[curpos]>126) result[curpos]=126;
            curpos++;
            if (linelen > 0 && i % linelen == linelen - 1) {
                result[curpos] = '\n';
                curpos++;
            }
        }
    }
    else {
        result[curpos]='\n';
        curpos++;
    }
    result[reslen] = 0;
    return result;
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
    else
        target->name = NULL;
}
void contig_remove_qscores(contig_s* target){
    if(target->qscore)free(target->qscore);
    target->qscore=NULL;
}
char* contig_nameptr(contig_s* target) {
    return target->name;
}
void contig_append_contig(contig_s* target, contig_s* extra) {
    size_t newlen;
    size_t i;
    newlen = target->length + extra->length;
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
    target->length = newlen;
}
contig_s* contig_subsection(contig_s* target, size_t start, size_t len) {
    contig_s* result;
    if (start + len > target->length) len = target->length - start;
    result = (contig_s*)calloc(1, sizeof(contig_s));
    result->nalloc = len;
    result->length = len;
    result->nucleotides = memcpyalloc(target->nucleotides + start, len);
    if (target->coverage) result->coverage = memcpyalloc(target->coverage + start * sizeof(float), len * sizeof(float));
    if (target->qscore) result->qscore = memcpyalloc(target->qscore + start, len);
    return result;
}
void contig_overwrite(contig_s* target, int64_t offset, contig_s* src, size_t start, size_t len) {
    size_t i, maxi;
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
    *p_A = 0;
    *p_C = 0;
    *p_G = 0;
    *p_T = 0;
    *p_N = 0;
    for (pos = 0; pos<target->length; pos++) {
        switch (target->nucleotides[pos] + 1) {
        case 0:
            *p_N += 1;
            break;
        case 1:
            *p_A += 1;
            break;
        case 2:
            *p_C += 1;
            break;
        case 3:
            *p_G += 1;
            break;
        case 4:
            *p_T += 1;
            break;
        default:
            *p_N += 1;
            break;
        }
    }
}

char* contig_qscores(contig_s* target){
    return target->qscore;
}
int contig_min_qscore(contig_s* target){
    size_t i;
    int res;
    res = 256;
    if(target->qscore) {
        for(i=0;i<target->length;i++){
            if( (int)(target->qscore[i]) < res) res = (int)(target->qscore[i]);
        }
    }
    else return -1;
    return res;
}
float contig_avg_qscore(contig_s* target){
    size_t i;
    int64_t res;
    res = 0;
    if(target->qscore) {
        for(i=0;i<target->length;i++){
            res += (int)(target->qscore[i]);
        }
        return (float)(((double)res)/((double)(target->length)));
    }
    else return -1;
    return res;
}
int contig_max_qscore(contig_s* target){
    size_t i;
    int res;
    res = -1;
    if(target->qscore) {
        for(i=0;i<target->length;i++){
            if( (int)(target->qscore[i]) > res) res = (int)(target->qscore[i]);
        }
    }
    else return -1;
    return res;
}

contig_s* contig_from_fastaPF(PF_t* src) {
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
            tmpcontig = contig_from_str(line, NULL);
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
            tmpcontig = contig_from_str(line, NULL);
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
                    tmpcontig = contig_from_str(line, NULL);
                    contig_append_contig(curcontig, tmpcontig);
                    contig_free(tmpcontig);
                }
            }
            else{
               linelen = strlen(line);
               if (nqcs + linelen >= curcontig->length) {
                   state = 0;
               }
               memcpy(curcontig->qscore+nqcs,line,linelen);
               i = nqcs;
               nqcs += linelen;
               if(state==0) nqcs = curcontig->length;
               while(i < nqcs){
                   curcontig->qscore[i] -= 33;
                   i++;
               }
            }
        }
        free(line);
        filepos = PFtell(f);
    }
    if(!line && p_resume) *p_resume = 0;
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
    PFclose(f);
    return result;
}
contig_s** contigs_from_fastx(char* filename, size_t* p_amount) {
    return contigs_from_fastx_limited(filename, p_amount, 0, NULL);
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

