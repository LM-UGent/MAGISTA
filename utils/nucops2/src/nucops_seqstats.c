#include <stdio.h>
#include "osportstd.h"
#include "sequence_base.h"

int nucops_GC(int argc, char** argv){
  size_t A,C,G,T,N;
  size_t a,c,g,t,n;
  size_t nb;
  size_t ncontigs, i;
  contig_s** contigs;
  if(argc != 2 || argv[1][0]=='-'){
    fprintf(stderr,"Usage: nucops2 GC <filename>\n");
  }
  else{
    contigs = contigs_from_fasta(argv[1], &ncontigs);
    A = C = G = T = N = 0;
    for(i=0;i<ncontigs;i++){
      contig_countACGTN(contigs[i],&a,&c,&g,&t,&n);
      A += a;
      C += c;
      G += g;
      T += t;
      N += n;
    }
    nb = A+C+G+T+N;
    if(nb>0){
      fprintf(stdout,"%f\n",((double)(G+C))/((double)nb));
    }
    else{
      fprintf(stdout,"%f\n",0.0);
    }
  }
  return 0;
}
