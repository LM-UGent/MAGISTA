#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "nucops_winsplit.h"
#include "nucops_select.h"
#include "nucops_concat.h"
#include "nucops_mutate.h"
#include "nucops_fq2fa.h"
#include "nucops_shorten.h"
#include "nucops_seqstats.h"
#include "nucops_summaries.h"
#include "nucops_renameseqs.h"

#ifdef _DEBUG
int main(int argc, char** argv) {
    /*
    char* fake_argv[] = { "seqsummary","-i","C:\\Users\\ggoussar\\Documents\\Visual Studio 2015\\Projects\\nucops2\\nucops2\\reads.few.fastq","-o","reads.Q20.summary.tsv"};
    int fake_argc;
    fake_argc = 5;
    nucops_seqsummary(fake_argc, fake_argv);
    */
    char* fake_argv[] = { "renameseqs","-i","C:\\Users\\ggoussar\\Documents\\Visual Studio 2015\\Projects\\nucops2\\nucops2\\reads.few.fastq","-q","none" };
    int fake_argc;
    fake_argc = 5;
    nucops_renameseqs(fake_argc, fake_argv);
    system("PAUSE");
    return 0;
}
#endif
