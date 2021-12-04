/*
MIT License

Copyright (c) 2019 Gleb Goussarov

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <time.h>
#include "textparsing.h"
#include "osportstd.h"

int os_fileexists(const char* path) {
#ifdef _WIN32
    if (_access(path, 0) != -1)
        return 1;
#else
    if (access(path, F_OK) != -1)
        return 1;
#endif
    return 0;
}
char* os_rmdirname(char* path) {
    size_t i;
    i = strlen(path);
    while (i > 0 && path[i - 1] != '/' && path[i - 1] != '\\') {
        i--;
    }
    return path + i;
}
char* os_extractdirname(char* path) {
    char* result;
    size_t i;
    i = strlen(path);
    while (i > 0 && path[i - 1] != '/' && path[i - 1] != '\\') {
        i--;
    }
    if (i == 0) result = "";
    else {
        result =(char*)  malloc(i + 1);
        memcpy(result, path, i);
        result[i] = 0;
    }
    return result;
}
#ifdef _WIN32
/* helper functions for windows-specific prcoess output parsing */
size_t _win32_ReadFromPipe(HANDLE g_hChildStd_OUT_Rd, char* result, size_t buffersize, int* nullflag, int wait) {
    DWORD readbytes, bytesavailable, bltm;
    char tmpc;
    int tmp;
    readbytes = 0;
    OVERLAPPED sto = { 0 };

    if (!wait)
        PeekNamedPipe(g_hChildStd_OUT_Rd, &tmpc, 1, &readbytes, &bytesavailable, &bltm);
    else
        bytesavailable = 1;
    if (bytesavailable == 0 || !ReadFile(g_hChildStd_OUT_Rd, result, (int32_t)buffersize - 1, &readbytes, &sto)) {
        *nullflag = 1;
    }
    else {
        *nullflag = 0;
    }
    tmp = GetLastError();
    result[readbytes] = 0;
    return readbytes;
}
HANDLE _win32_CreateChildProcess(char* prog, int argc, char** argv,HANDLE g_hChildStd_IN_Rd, HANDLE g_hChildStd_IN_Wr, HANDLE g_hChildStd_OUT_Rd, HANDLE g_hChildStd_OUT_Wr) {
    PROCESS_INFORMATION procinfo;
    STARTUPINFO startinfo;
    SECURITY_ATTRIBUTES secattr;
    HANDLE hresult;
    char eof = EOF;

    hresult = INVALID_HANDLE_VALUE;

    char* cmdline;
    cmdline = text_join(argv, "\" \"", "\"","\"", argc);
    /* note : for the sake of clarity I avoid using microsoft macros when a serviceable ANSI-C alternative exists */
    memset(&procinfo, 0, sizeof(PROCESS_INFORMATION));
    memset(&startinfo, 0, sizeof(STARTUPINFO));
    startinfo.cb = sizeof(STARTUPINFO);
    startinfo.hStdError = GetStdHandle(STD_ERROR_HANDLE);
    startinfo.hStdOutput = g_hChildStd_OUT_Wr;
    startinfo.hStdInput = g_hChildStd_IN_Rd;
    startinfo.dwFlags |= STARTF_USESTDHANDLES;
    memset(&secattr, 0, sizeof(SECURITY_ATTRIBUTES));
    secattr.bInheritHandle = 1;
    if (!CreateProcess(prog, cmdline, &secattr, &secattr, 1, 0, NULL, NULL, &startinfo, &procinfo)) {
        fprintf(stderr, "Process did not succeed: %s\n", cmdline);
    }
    else {
        CloseHandle(g_hChildStd_OUT_Wr);
        CloseHandle(g_hChildStd_IN_Rd);
        CloseHandle(procinfo.hThread);
        WaitForInputIdle(procinfo.hProcess, INFINITE);
        hresult = procinfo.hProcess;
    }
    return hresult;
}
#endif

char* os_stdoutfromexec(char* program, int argc, char** argv, size_t* outlen, size_t bufsize) {
    char* buffer;
    char* result;
    size_t totlen;
    size_t curlen;
    size_t oldlen;
#ifdef _WIN32
    SECURITY_ATTRIBUTES sattr;
    int nf;

    HANDLE g_hChildStd_IN_Rd = NULL;
    HANDLE g_hChildStd_IN_Wr = NULL;
    HANDLE g_hChildStd_OUT_Rd = NULL;
    HANDLE g_hChildStd_OUT_Wr = NULL;
    HANDLE hprochandle;
    int32_t exitstatus;

    sattr.nLength = sizeof(SECURITY_ATTRIBUTES);
    sattr.bInheritHandle = TRUE;
    sattr.lpSecurityDescriptor = NULL;
    /*
    Initialize the pipes and make sure the child does not inherit write/read access
    to the wrong end of each pipe. Then launch the external process.
    */
    if (!CreatePipe(&g_hChildStd_OUT_Rd, &g_hChildStd_OUT_Wr, &sattr, 0))
        fprintf(stderr, "Could not open Pipe to child process STDOUT\n");
    if (!SetHandleInformation(g_hChildStd_OUT_Rd, HANDLE_FLAG_INHERIT, 0))
        fprintf(stderr, "Could not prevent inheritance of child STDOUT read access\n");
    if (!CreatePipe(&g_hChildStd_IN_Rd, &g_hChildStd_IN_Wr, &sattr, 0))
        fprintf(stderr, "Could not open Pipe to child process STDIN\n");
    if (!SetHandleInformation(g_hChildStd_IN_Wr, HANDLE_FLAG_INHERIT, 0))
        fprintf(stderr, "Could not prevent inheritance of child STDIN write access\n");
    hprochandle = _win32_CreateChildProcess(program, argc, argv, g_hChildStd_IN_Rd, g_hChildStd_IN_Wr, g_hChildStd_OUT_Rd, g_hChildStd_OUT_Wr);

    nf = 0;
    totlen = oldlen = 0;
    buffer = (char*)calloc(bufsize, 1);
    curlen = _win32_ReadFromPipe(g_hChildStd_OUT_Rd, buffer, bufsize, &nf, 0);
    result = (char*)malloc(curlen+1);
    GetExitCodeProcess(hprochandle, &exitstatus);
    /*
    nf can be raised if the process does not generate its stdout output quickly enough.
    Also, some programs will not termina before their ouput has been parsed.
    This explains the conditions below
    */
    while (!nf || exitstatus == STILL_ACTIVE) {
        curlen = _win32_ReadFromPipe(g_hChildStd_OUT_Rd, buffer, bufsize, &nf, 0);
        if (curlen > 0) {
            totlen += curlen;
            result = (char*)realloc(result, totlen + 1);
            memcpy(result + oldlen, buffer, curlen);
            oldlen = totlen;
        }
        GetExitCodeProcess(hprochandle, &exitstatus);
        if (nf && exitstatus != STILL_ACTIVE)Sleep(1);
    }
    curlen = _win32_ReadFromPipe(g_hChildStd_OUT_Rd, buffer, bufsize, &nf, 0);
    if (curlen > 0) {
        totlen += curlen;
        result = (char*)realloc(result, totlen + 1);
        memcpy(result + oldlen, buffer, curlen);
    }
    free(buffer);
    result[totlen] = 0;
#else
    int exitstatus;
    int pipefd[6];
    pid_t cpid;
    char** nargv;
    result = NULL;
    if (argc < 1) return NULL;

    if (pipe(pipefd) == -1) {
        fprintf(stderr, "could not create stdout pipe to %s:\n", program);
        return NULL;
    }
    if (pipe(pipefd+2) == -1) {
	close(pipefd[0]); close(pipefd[1]);
        fprintf(stderr, "could not create stderr pipe to %s\n", program);
        return NULL;
    }
    if (pipe(pipefd+4) == -1) {
	close(pipefd[0]); close(pipefd[1]);
	close(pipefd[2]); close(pipefd[3]);
        fprintf(stderr, "could not create stdin pipe to %s\n", program);
        return NULL;
    }
    /* write pipe is not used */
    cpid = fork();
    if (cpid == -1) {
        fprintf(stderr, "could not fork process\n");
    }
    else if (cpid != 0) {
        /* parent process */
        close(pipefd[1]);
        close(pipefd[3]);
        close(pipefd[4]);

        totlen = oldlen = 0;
        buffer = (char*)calloc(bufsize, 1);
        curlen = 1000;
        curlen = read(pipefd[0], buffer, bufsize - 1);
        while ( curlen > 0) {
            totlen += curlen;
            result = (char*)realloc(result, totlen + 1);
            memcpy(result + oldlen, buffer, curlen);
            result[totlen] = 0;
            oldlen = totlen;
            curlen = read(pipefd[0], buffer, bufsize - 1);
            if(waitpid(cpid, &exitstatus, WNOHANG))break;
        }
        while ( curlen > 0) {
            totlen += curlen;
            result = (char*)realloc(result, totlen + 1);
            memcpy(result + oldlen, buffer, curlen);
            result[totlen] = 0;
            oldlen = totlen;
            curlen = read(pipefd[0], buffer, bufsize - 1);
        }
        if(!result || strlen(result)<1) {
            fprintf(stderr,"program produced no output\n");
        }
        free(buffer);
        close(pipefd[0]);
        close(pipefd[2]);
	close(pipefd[5]);
        wait(0);
    }
    else {
        /* child process */
        close(pipefd[0]);
	close(pipefd[2]);
        close(pipefd[5]);
        while ((dup2(pipefd[1], STDOUT_FILENO) == -1) && (errno == EINTR)) {}
        while ((dup2(pipefd[3], STDERR_FILENO) == -1) && (errno == EINTR)) {}
        while ((dup2(pipefd[4], STDIN_FILENO) == -1) && (errno == EINTR)) {}

        nargv=malloc(sizeof(char*)*(argc+1));
        memcpy(nargv,argv,sizeof(char*)*(argc));
        nargv[argc]=NULL;
        exitstatus = execv(nargv[0], nargv);
	if(exitstatus){
            fprintf(stdout,"program exited with status [%d], Arguments (%d) were:\n",errno, argc);
            curlen = 0;
/*
            while(nargv[curlen] != NULL) {
                fprintf(stdout,"%s ",nargv[curlen++]);
            }
            fprintf(stdout,"\n");*/
        }
        close(pipefd[1]);
        close(pipefd[3]);
	close(pipefd[4]);
        free(nargv);
        _exit((char)errno);
    }
#endif
    *outlen = totlen;
    return result;
}

void* memcpyalloc(void* data, size_t datasize) {
    void* result;
    result = malloc(datasize);
    memcpy(result, data, datasize);
    return result;
}


unimem PFreadunimem16(PF_t* f) {
    uint16_t tmpsz;
    unimem output;
    tmpsz = PFgetint16(f);
    output = (unimem)malloc((size_t)tmpsz);
    PFread(((char*)output) + 1, 1, tmpsz, f);
    return output;
}
unimem PFreadunimem32(PF_t* f) {
    uint32_t tmpsz;
    unimem output;
    tmpsz = PFgetint32(f);
    output = (unimem)malloc((size_t)tmpsz);
    PFread(((char*)output) + 1, 1, tmpsz, f);
    return output;
}
unimem PFreadunimem64(PF_t* f) {
    uint64_t tmpsz;
    unimem output;
    tmpsz = PFgetint64(f);
    output = (unimem)malloc((size_t)tmpsz);
    PFread(((char*)output) + 1, 1, (size_t)tmpsz, f);
    return output;
}