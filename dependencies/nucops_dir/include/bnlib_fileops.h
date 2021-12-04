#ifndef DEF_FILESYS
#define DEF_FILESYS

/*#include "bnlib_params.h"*/
#include <stdint.h>


#ifndef PFBUFFER
#define PFBUFFER    0x100000
#endif

#ifndef PFSTDIOENABLE
#define PFSTDIOENABLE 1
#endif

typedef int errcode_t;
typedef struct portablefile_s PF_s;

/* adapted std calls*/
errcode_t PFopen(PF_s** f, const char* path, const char* mode);
errcode_t PFclose(PF_s* f);
size_t PFread(void* buffer, size_t element_size, size_t element_count, PF_s* f);
size_t PFwrite(const void* buffer, size_t element_size, size_t element_count, PF_s* f);
long long PFtell(PF_s* f);
errcode_t PFseek(PF_s* f, long long offset, int origin);
errcode_t PFrewind(PF_s* f);
int PFgetc(PF_s* f);
errcode_t PFputc(int c, PF_s* f);
errcode_t PFprintf(PF_s* f, const char* fmt, ...);
errcode_t PFvprintf(PF_s* f, const char* fmt, va_list argp);
/* extra functions */
errcode_t PFopen_spec(PF_s** f, const char* path, const char* mode, size_t buffersize);
size_t PFpoll(void* buffer, size_t element_size, size_t element_count, PF_s* f);
errcode_t PFprintx(PF_s* f, char c, size_t ntimes);
int PFseemslikeASCIItext(PF_s* f, size_t bytes_to_test);
char* PFgetfullpathptr(PF_s* f);
char* PFgetbasenameptr(PF_s* f);
FILE* PFgetFILE(PF_s* f);
size_t PFnumlines(PF_s* f, int countempty);

int16_t PFgetint16(PF_s* f);
int32_t PFgetint32(PF_s* f);
int64_t PFgetint64(PF_s* f);
float PFgetfloat(PF_s* f);
double PFgetdouble(PF_s* f);
int PFpollc(PF_s* f);
int16_t PFpollint16(PF_s* f);
int32_t PFpollint32(PF_s* f);
int64_t PFpollint64(PF_s* f);
float PFpollfloat(PF_s* f);
double PFpolldouble(PF_s* f);
/* note: PFreadline removes \n and \r at the end of lines */
char* PFreadline(PF_s* f);
char* PFpollline(PF_s* f);

errcode_t PFputint16(PF_s* f, int16_t input);
errcode_t PFputint32(PF_s* f, int32_t input);
errcode_t PFputint64(PF_s* f, int64_t input);
errcode_t PFputfloat(PF_s* f, float input);
errcode_t PFputdouble(PF_s* f, double input);


#endif

