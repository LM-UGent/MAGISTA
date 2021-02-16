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

int main(int argc, char** argv) {

    if (argc < 2) {
        fprintf(stderr, "Usage: %s <command> [options]\n", argv[0]);
        fprintf(stderr, "Available commands are: 'concat', 'select', 'winsplit', 'fastasummary'\n\n");
    }
    else if (strcmp(argv[1], "fastasummary") == 0) return nucops_fastasummary(argc - 1, argv + 1);
    else if (strcmp(argv[1], "concat") == 0) return nucops_concat(argc - 1, argv + 1);
    else if (strcmp(argv[1], "select") == 0) return nucops_select(argc - 1, argv + 1);
    else if (strcmp(argv[1], "winsplit") == 0) return nucops_winsplit(argc - 1, argv + 1);
    else if (strcmp(argv[1], "mutate") == 0) return nucops_mutate(argc - 1, argv + 1);
    else if (strcmp(argv[1], "fq2fa") == 0) return nucops_fq2fa(argc - 1, argv + 1);
    else if (strcmp(argv[1], "GC")==0) return nucops_GC(argc - 1, argv + 1);
    else if (strcmp(argv[1], "shorten")==0) return nucops_shorten(argc - 1, argv + 1);
    else {
        fprintf(stderr, "Unrecognized command '%s'\n\n", argv[1]);
    }

    return 1;
}

