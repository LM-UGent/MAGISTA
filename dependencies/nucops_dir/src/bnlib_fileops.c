#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>

/*#include "bnlib_portability.h"*/
/*#include "bnlib_diagnostics.h"*/
#include "bnlib_fileops.h"

/* nucops currently does not actually include bnlib - so some things have to be defined rather than included */
#ifndef MSVC_COMPILER
#define MSVC_COMPILER 1
#define OTHER_COMPILER 0
#endif
#ifdef _WIN32
#define PORT__COMPILER MSVC_COMPILER
#ifdef _WIN64
#define PORT__BITS 64
#else 
#define PORT__BITS 32
#endif
#else
#define PORT__COMPILER OTHER_COMPILER
#define PORT__BITS 64
#endif
#define dgn_malloc malloc
#define dgn_calloc calloc
#define dgn_free free
#define dgn_realloc realloc

#include <stdarg.h>


#if PORT__COMPILER == MSVC_COMPILER
#include <Windows.h>
#endif


struct portablefile_s {
    FILE* filepointer;
    char* filename;
    size_t fnlen, basenamestart;
    char* prevbuffer;
    char* curbuffer;
    char* nextbuffer;
    size_t pbs, cbs, nbs, tbs;
    long long pbp, cbp, nbp;
    size_t cursor;
    int closeable, isstdstream;
};

static inline size_t local__compatiblefread(void* buffer, size_t eltsize, size_t eltcount, FILE* f){
    size_t result;
    #if PORT__COMPILER == MSVC_COMPILER
    result = _fread_nolock(buffer, eltsize, eltcount, f);
    (void)ferror(f);
    #else
    result = fread(buffer, eltsize, eltcount, f);
    #endif
    return result;
}
static inline long long local__compatibleftell(FILE* f){
    long long result;
    #if PORT__COMPILER == MSVC_COMPILER
    #ifdef _WIN64
    return _ftelli64(f);
    #else
    return (long long)ftell(f);
    #endif
    #else
    return (long long)ftell(f);
    #endif
    return result;
}
static errcode_t local__PFnext(PF_s* f){
    char* tmpptr;
    if(!f || !f->filepointer) return 2;
    f->cursor = 0;
    tmpptr = f->prevbuffer;
    f->prevbuffer = f->curbuffer;
    f->curbuffer = f->nextbuffer;
    f->nextbuffer = tmpptr;
    f->pbs = f->cbs;
    f->pbp = f->cbp;
    if(!f->curbuffer) f->curbuffer = dgn_malloc(f->tbs+1);
    if(f->cbs==0){
        f->curbuffer[f->tbs]=0;
        f->cbp = local__compatibleftell(f->filepointer);
        f->cbs = local__compatiblefread(f->curbuffer,1,f->tbs,f->filepointer);
        f->curbuffer[f->cbs]=0;
    }
    else {
        f->cbp = f->nbp;
        f->cbs = f->nbs;
    }
    if(f->cbs < f->tbs){
        f->curbuffer[f->cbs]=0;
        f->nbs = 0;
        f->nbp = 0;
        return -1;
    }
    if(!f->nextbuffer) f->nextbuffer = dgn_malloc(f->tbs+1);
    f->nbp = local__compatibleftell(f->filepointer);
    f->nbs = local__compatiblefread(f->nextbuffer,1,f->tbs,f->filepointer);
    f->nextbuffer[f->nbs]=0;
    return 0;
}
/* adapted std calls*/
errcode_t PFopen(PF_s** f, const char* path, const char* mode){
    return PFopen_spec(f, path, mode, PFBUFFER);
}
errcode_t PFclose(PF_s* f) {
    if (f) {
        if (f->filepointer && f->closeable)fclose(f->filepointer);
        f->closeable = 0;
        if (f->filename)dgn_free(f->filename);
        memset(f, 0, sizeof(PF_s));
        if (f->prevbuffer) dgn_free(f->prevbuffer);
        if (f->curbuffer) dgn_free(f->curbuffer);
        if (f->nextbuffer) dgn_free(f->nextbuffer);
        f->prevbuffer = f->curbuffer = f->nextbuffer = NULL;
        dgn_free(f);
        return 0;
    }
    return 1;
}
size_t PFread(void* buffer, size_t element_size, size_t element_count, PF_s* f) {
    size_t copysize;
    size_t resultsize;
    size_t curstart;
    size_t remaining;
    int eof_reached;
    int nexterr;
    char* cbuf;
    if (!f || !(f->filepointer)) {
        return 0;
    }
    cbuf = (char*)buffer;
    resultsize = element_size*element_count;
    curstart = 0;
    remaining = resultsize;
    if (f->cursor + remaining < f->cbs) {
        copysize = remaining;
    }
    else {
        copysize = f->cbs - f->cursor;
    }
    if(copysize>0){
        memcpy(cbuf, f->curbuffer + f->cursor, copysize);
        f->cursor += copysize;
        remaining -= copysize;
        curstart = copysize;
    }
    eof_reached = 0;
    while (remaining>0 && !eof_reached) {
        if(f->cbs < f->tbs) nexterr=1;
        else nexterr = local__PFnext(f);
        if(nexterr != 0) eof_reached = 1;
        else {
            if (f->cursor + remaining < f->cbs) {
                copysize = remaining;
            }
            else {
                copysize = f->cbs - f->cursor;
            }
            if(copysize>0){
                memcpy(cbuf + curstart, f->curbuffer, copysize);
                f->cursor += copysize;
                remaining -= copysize;
                curstart += copysize;
            }
        }
    }
    return curstart;
}
size_t PFwrite(const void* buffer, size_t element_size, size_t element_count, PF_s* f) {
    if (f && f->filepointer) {
        return fwrite(buffer, element_size, element_count, f->filepointer);
    }
    return 0;
}
long long PFtell(PF_s* f) {
    if(!f || !(f->filepointer) || f->isstdstream) return -1;
    return (long long)(f->cursor + f->cbp);
}
errcode_t PFseek(PF_s* f, long long offset, int origin) {
    errcode_t res;
    if (!f || !(f->filepointer))return 1;
    if (offset >= f->cbp && offset < f->cbp + (long long)(f->cbs)){
        f->cursor = ((size_t)offset) - f->cbp;
        return 0;
    }
    f->pbs = 0;
    f->cbs = 0;
    f->nbs = 0;
#if PORT__COMPILER == MSVC_COMPILER
#if PORT__BITS == 64
    res = _fseeki64(f->filepointer, offset, origin);
#else
    res = fseek(f->filepointer, (long)offset, origin);
#endif
#else
    res = fseek(f->filepointer, offset, origin);
#endif
    if(!res){
        res = local__PFnext(f);
    }
    return res;
}
errcode_t PFrewind(PF_s* f) {
    if (!f)return 1;
    f->pbs = 0;
    f->cbs = 0;
    f->nbs = 0;
    f->cursor = 0;
    rewind(f->filepointer);
    return local__PFnext(f);
}
int PFgetc(PF_s* f) {
    if (!f)return -1;
    if (f->cursor < f->cbs) {
        f->cursor++;
        return f->curbuffer[f->cursor - 1];
    }
    if(local__PFnext(f))return -1;
    f->cursor=1;
    return f->curbuffer[0];
}
errcode_t PFputc(int c, PF_s* f) {
    if (!f)return -1;
    return fputc(c, f->filepointer);
}
errcode_t PFprintf(PF_s* f, const char* fmt, ...) {
    errcode_t res;
    va_list argp;
    va_start(argp, fmt);
    if (!f || !(f->filepointer))return -1;
    res = (errcode_t)vfprintf(f->filepointer, fmt, argp);
    va_end(argp);
    return res;
}
errcode_t PFvprintf(PF_s* f, const char* fmt, va_list argp) {
    return vfprintf(f->filepointer, fmt, argp);
}

/* extra functions */
errcode_t PFopen_spec(PF_s** f, const char* path, const char* mode, size_t buffersize){
    errcode_t outval;
    PF_s* result;
    size_t i, j;
    int readable;

    result = dgn_calloc(1,sizeof(PF_s));
    result->tbs = buffersize;
    result->fnlen = strlen(path);
    result->filename = dgn_malloc(result->fnlen + 1);
    memcpy(result->filename, path, result->fnlen + 1);
    j = result->fnlen;
    readable = 0;
    for (i = result->fnlen;i > 0;i--) {
        j--;
        if (result->filename[j] == '/' || result->filename[j] == '\\') {
            result->basenamestart = j+1;
            break;
        }
    }
    if (i == 0)result->basenamestart = 0;
    outval = 0;
    #if PFSTDIOENABLE == 1
    if (strcmp(path, "stdout") == 0) {
        result->filepointer = stdout;
        result->isstdstream = 1;
    }
    else if (strcmp(path, "stdin") == 0) {
        result->filepointer = stdin;
        result->isstdstream = 1;
        readable = 1;
    }
    else if (strcmp(path, "stderr") == 0) {
        result->filepointer = stderr;
        result->isstdstream = 1;
    }
    else
    #endif
    if(1){
        #if PORT__COMPILER == MSVC_COMPILER
        outval = (errcode_t)fopen_s(&(result->filepointer), path, mode);
        #else
        result->filepointer = fopen(path, mode);
        if (!result->filepointer) outval = 2;
        #endif
        if(result->filepointer){
            result->closeable=1;
        }
    }
    i = 0;
    while(mode[i]) {
        if (mode[i] == 'r')readable = 1;
        i++;
    }
    if(result->filepointer && readable){
        local__PFnext(result);
    }
    *f = result;
    return outval;
}

/* Works like read but does not move cursor - cannot return more characters than is available in two buffers */
size_t PFpoll(void* buffer, size_t element_size, size_t element_count, PF_s* f) {
    size_t copysize;
    size_t resultsize;
    size_t curstart;
    size_t remaining;
    char* cbuf;
    if (!f || !(f->filepointer)) {
        return 0;
    }
    cbuf = (char*)buffer;
    resultsize = element_size*element_count;
    curstart = 0;
    remaining = resultsize;
    if (f->cursor + remaining < f->cbs) {
        copysize = remaining;
    }
    else {
        copysize = f->cbs - f->cursor;
    }
    if(copysize>0){
        memcpy(cbuf, f->curbuffer + f->cursor, copysize);
        remaining -= copysize;
        curstart = copysize;
    }
    if (remaining>0) {
        if (remaining < f->cbs) {
            copysize = remaining;
        }
        else {
            copysize = f->cbs;
        }
        if(copysize>0){
            memcpy(cbuf, f->curbuffer + f->cursor, copysize);
            curstart = copysize;
        }
    }
    return curstart;
}
errcode_t PFprintx(PF_s* f, char c, size_t ntimes) {
    errcode_t res;
    char* tmp;
    if (!f)return -1;
    else res = 0;
    tmp = dgn_malloc(ntimes + 1);
    memset(tmp, c, ntimes);
    tmp[ntimes] = 0;
    res = (errcode_t)fprintf(f->filepointer, tmp);
    dgn_free(tmp);
    return res;
}
int PFseemslikeASCIItext(PF_s* f, size_t bytes_to_test) {
    char* buffer;
    size_t i;
    int result;
    if (!f || !(f->filepointer)) return 2;
    buffer = (char*)dgn_malloc(bytes_to_test);
    result = 1;
    if(f->curbuffer==NULL || f->cbs==0)local__PFnext(f);
    buffer = f->curbuffer;
    if(bytes_to_test == 0 || bytes_to_test > f->cbs) bytes_to_test = f->cbs;
    for (i = 0;i < bytes_to_test;i++) {
        if (buffer[i] == '\t' || buffer[i] == '\n' || buffer[i] == '\r') continue;
        if (buffer[i] < ' ') {
            result = 0;
            break;
        }
        if (result > '~') {
            result = 0;
            break;
        }
    }
    return result;
}
char* PFgetfullpathptr(PF_s* f) {
    return f->filename;
}
char* PFgetbasenameptr(PF_s* f) {
    return f->filename + f->basenamestart;
}
FILE* PFgetFILE(PF_s* f) {
    if (!f)return NULL;
    return f->filepointer;
}
size_t PFnumlines(PF_s* f, int countempty) {
    long long start;
    size_t nread;
    size_t i;
    size_t nlines;
    int stop;
    int nonempty;
    char buffer[0x10000] = { 0 };
    if (!f || !f->filepointer) return 0;
    start = local__compatibleftell(f->filepointer);
    rewind(f->filepointer);
    stop = 0;
    nonempty = 0;
    nlines = 0;
    while (!stop) {
        nread = fread(buffer, 1, 0x10000 - 1, f->filepointer);
        buffer[nread] = 0;
        for (i = 0;i < nread;i++) {
            if (buffer[i] == '\n') {
                if (countempty || nonempty)nlines++;
                nonempty = 0;
            }
            else if (nonempty == 0 && buffer[i] != '\r') {
                nonempty = 1;
            }
        }
        if (nread < 0x10000 - 1)stop = 1;
    }
    if (nonempty)nlines++;
#if PORT__COMPILER == MSVC_COMPILER && PORT__BITS == 64
    _fseeki64(f->filepointer, start, SEEK_SET);
#else
    fseek(f->filepointer, start, SEEK_SET);
#endif
    return nlines;
}

int16_t PFgetint16(PF_s* f) {
    int16_t res = 0;
    if (f && f->filepointer)
        PFread(&res, sizeof(int16_t), 1, f);
    return res;
}
int32_t PFgetint32(PF_s* f) {
    int32_t res = 0;
    if (f && f->filepointer)
        PFread(&res, sizeof(int32_t), 1, f);
    return res;
}
int64_t PFgetint64(PF_s* f) {
    int64_t res = 0;
    if (f && f->filepointer)
        PFread(&res, sizeof(int64_t), 1, f);
    return res;
}
float PFgetfloat(PF_s* f) {
    float res = 0.0;
    if (f && f->filepointer)
        PFread(&res, sizeof(float), 1, f);
    return res;
}
double PFgetdouble(PF_s* f) {
    double res = 0.0;
    if (f && f->filepointer)
        PFread(&res, sizeof(double), 1, f);
    return res;
}

int PFpollc(PF_s* f){
    if(f->cursor<f->cbs)return f->curbuffer[f->cursor];
    if(f->nbs == 0) return -1;
    return f->nextbuffer[0];
}
int16_t PFpollint16(PF_s* f) {
    int16_t res = 0;
    if (f && f->filepointer)
        PFpoll(&res, sizeof(int16_t), 1, f);
    return res;
}
int32_t PFpollint32(PF_s* f) {
    int32_t res = 0;
    if (f && f->filepointer)
        PFpoll(&res, sizeof(int32_t), 1, f);
    return res;
}
int64_t PFpollint64(PF_s* f) {
    int64_t res = 0;
    if (f && f->filepointer)
        PFpoll(&res, sizeof(int64_t), 1, f);
    return res;
}
float PFpollfloat(PF_s* f) {
    float res = 0.0;
    if (f && f->filepointer)
        PFpoll(&res, sizeof(float), 1, f);
    return res;
}
double PFpolldouble(PF_s* f) {
    double res = 0.0;
    if (f && f->filepointer)
        PFpoll(&res, sizeof(double), 1, f);
    return res;
}

size_t local__PFnextlinestart(char* buffer, size_t start, size_t bufsize){
    if (buffer[start] == 0)return bufsize;
    if (buffer[start] == '\n')return start + 1;
    start++;
    while(start<bufsize){
        if(buffer[start]==0)return bufsize;
        if(buffer[start]=='\n')return start+1;
        if(buffer[start-1]=='\r')return start;
        start++;
    }
    return start;
}

char* PFreadline(PF_s* f) {
    size_t totlinelen;
    size_t offset;
    size_t readbytes;
    size_t bufpos;
    size_t nls;
    int eob_reached;
    char* line;
    if (!f || !f->filepointer) return NULL;
    if (f->cursor == f->cbs && f->nbs == 0)return NULL;
    bufpos = (size_t)f->cursor;
    totlinelen = 0;
    nls = local__PFnextlinestart(f->curbuffer,f->cursor,f->cbs);
    readbytes = nls - f->cursor;
    line = dgn_malloc(readbytes + 1);
    memcpy(line, f->curbuffer + f->cursor, readbytes);
    totlinelen = readbytes;
    while (f->cbs == f->tbs && nls==f->cbs) {
        local__PFnext(f);
        if (totlinelen == 0 || line[totlinelen-1] != '\n') {
            nls = local__PFnextlinestart(f->curbuffer, 0, f->cbs);
            readbytes = nls;
            offset = totlinelen;
            totlinelen += readbytes;
            line = dgn_realloc(line, totlinelen + 1);
            memcpy(line + offset, f->curbuffer, readbytes);
        }
        else {
            nls = 0;
        }
    }
    totlinelen--;
    while(line[totlinelen]=='\n' || line[totlinelen]=='\r')totlinelen--;
    line[totlinelen+1] = 0;
    f->cursor = nls;
    return line;
}
char* PFpollline(PF_s* f) {
    size_t totlinelen;
    size_t offset;
    size_t readbytes;
    size_t bufpos;
    size_t nls;
    char* line;
    if (!f || !f->filepointer) return NULL;
    if (f->cursor == f->cbs && f->nbs == 0)return NULL;
    bufpos = (size_t)f->cursor;
    totlinelen = 0;
    nls = local__PFnextlinestart(f->curbuffer,f->cursor,f->cbs);
    readbytes = nls - f->cursor;
    totlinelen = readbytes;
    if(f->cbs == f->tbs && nls==f->cbs) {
        nls = local__PFnextlinestart(f->nextbuffer,0,f->nbs);
        readbytes = nls;
        offset = totlinelen;
        totlinelen += readbytes;
        line = dgn_malloc(totlinelen + 1);
        memcpy(line + offset, f->nextbuffer, readbytes);
    }
    else {
        line = dgn_malloc(readbytes + 1);
    }
    memcpy(line, f->curbuffer + f->cursor, readbytes);
    while(line[totlinelen]=='\n' || line[totlinelen]=='\r')totlinelen--;
    line[totlinelen+1] = 0;
    return line;
}

errcode_t PFputint16(PF_s* f, int16_t input) {
    if (PFwrite(&input, sizeof(int16_t), 1, f) != sizeof(int16_t))return 1;
    return 0;
}
errcode_t PFputint32(PF_s* f, int32_t input) {
    if (PFwrite(&input, sizeof(int32_t), 1, f) != sizeof(int32_t))return 1;
    return 0;
}
errcode_t PFputint64(PF_s* f, int64_t input) {
    if (PFwrite(&input, sizeof(int64_t), 1, f) != sizeof(int64_t))return 1;
    return 0;
}
errcode_t PFputfloat(PF_s* f, float input) {
    if (PFwrite(&input, sizeof(float), 1, f) != sizeof(float))return 1;
    return 0;
}
errcode_t PFputdouble(PF_s* f, double input) {
    if (PFwrite(&input, sizeof(double), 1, f) != sizeof(double))return 1;
    return 0;
}

