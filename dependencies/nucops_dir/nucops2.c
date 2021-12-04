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
#include "nucops_compliance.h"
#include "nucops_samcov.h"

#ifndef _DEBUG
int main(int argc, char** argv) {
    char* alt_args[2];
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <command> [options]\n", argv[0]);
        fprintf(stderr, "Available commands are: 'complies', 'concat', 'select', 'winsplit', 'fastasummary', 'fq2fa', 'seqsummary', 'samcov'\n\n");
        fprintf(stderr, "Alternative Usage: %s <filename>\n", argv[0]);
        fprintf(stderr, "Equivalent to: <%s> fastasummary <filename>\n", argv[0]);
        fprintf(stderr, "Note: if the filename matches a command name, this will not work.\n\n");
    }
    else if (strcmp(argv[1], "fastasummary") == 0) return nucops_fastasummary(argc - 1, argv + 1);
    else if (strcmp(argv[1], "seqsummary") == 0) return nucops_seqsummary(argc - 1, argv + 1);
    else if (strcmp(argv[1], "concat") == 0) return nucops_concat(argc - 1, argv + 1);
    else if (strcmp(argv[1], "select") == 0) return nucops_select(argc - 1, argv + 1);
    else if (strcmp(argv[1], "splitfile") == 0) return nucops_splitfile(argc - 1, argv + 1);
    else if (strcmp(argv[1], "winsplit") == 0) return nucops_winsplit(argc - 1, argv + 1);
    else if (strcmp(argv[1], "mutate") == 0) return nucops_mutate(argc - 1, argv + 1);
    else if (strcmp(argv[1], "fq2fa") == 0) return nucops_fq2fa(argc - 1, argv + 1);
    else if (strcmp(argv[1], "GC")==0) return nucops_GC(argc - 1, argv + 1);
    else if (strcmp(argv[1], "shorten")==0) return nucops_shorten(argc - 1, argv + 1);
    else if (strcmp(argv[1], "renameseqs") == 0) return nucops_renameseqs(argc - 1, argv + 1);
    else if (strcmp(argv[1], "complies") == 0) return nucops_compliance(argc - 1, argv + 1);
    else if (strcmp(argv[1], "samcov") == 0) return nucops_samcov(argc - 1, argv + 1);
    else {
        fprintf(stderr, "Unrecognized command '%s' treated as file name.\n", argv[1]);
        alt_args[0] = "fastasummary";
        alt_args[1] = argv[1];
        nucops_fastasummary(2, alt_args);
#ifdef _WIN32
        system("PAUSE");
#endif
    }

    return 1;
}
#endif

