#ifndef BLOSUM_H
#define BLOSUM_H


/*=======================================================================
    (C) Copyright 1991, Fred Hutchinson Cancer Research Center
    motifj.h  Header file for PROTOMAT programs
    NOTE for Silicon Graphics users:  The type of scores in
    struct score should be changed from char to int to get correct
    processing (but not for SUN!)
------------------------------------------------------------------------
     6/29/90  J. Henikoff
>>>>>>>>>>>>>>>>>>>>>>>>>Blocks 8.0<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
     8/1/94  Changed MAXFREQ from 450 to 700 for PS01033
            NOTE: PS00028 wants 899 dups for 98 sequences => MAXFREQ = 997
            Changed VERSION from 7 to 8
     2/1/95  Increased MAX_LENGTH from 4000 to 5500
    7/11/95  Changed MOTAUTO to MOTAUTO4 & added MOTAUTO3
>>>>>>>>>>>>>>>>>>>>>>>>>Blocks 9.0<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    12/5/95  Changed MAXSEQS from 700 to 715 for PS50011
             Changed MAXFREQ from 700 to 715 also
    7/23/97  Added database type for Proclass PCFam
   10/23/97  Added struct pb_counts for pb_weights()
   10/28/97  Changed MOTAUTO4 from 2 to 3 and MOTAUTO3 from 12 to 6
    1/18/98  Changed MAXSEQS & MAXFREQ from 715 to 980
    2/27/98  Changed MAXFREQ from 980 to 1000 for PS00022
    5/20/15  Changed MAXFREQ from 1000 to 16583 for
=========================================================================*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "math.h"
#include "configurations.h"
#include "block.h"

struct block {                          /* Block structure */
   char ac[10];
   int nseq, width, strength, nclus;
   int aa[MAXSEQS][MAX_MERGE_WIDTH];    /* aas for seq */
   double weight[MAXSEQS];              /* seq weights found in block */
   int cluster[MAXSEQS];                /* cluster # for seq */
   int ncluster[MAXSEQS];               /* #seqs in same cluster */
   double totdiag, totoffd, wtot;
} Block;


class blosum
{
    double Counts[AAS][AAS];
    unsigned long TotPairs, TotSeqs, TotAas, TotWidth;
    double FTotPairs, FTotWeight;
    unsigned long TotBlk;     /* # blocks contributing to matrix */
    double TotClump;        /* # clumps contributing to matrix */
    unsigned long TotSeg;     /* # seqs contributing to matrix */
    unsigned long AaPairs[AAS];   /* Pairs for an AA */
    double FAaPairs[AAS], PBParameter;
    unsigned long AaFreq[AAS];    /* AA frequency */
    int MinStr, MaxStr, Cluster;
    char datfile[FNAMELEN];
    int iscale;
public:
    blosum();
    int  read_dat(FILE* fdat);
    void fill_block(FILE* fdat);
    void cluster_seqs();
    void count_block();
    void count_weight();
    void count_position();
    void count_cluster();
    void read_parameters(FILE* fdat, int argc,char *argv[]);


    void calculate_matrix();

    int  aachar_to_num(char ch);
    char *num_to_aachar(int num);

};

#endif // BLOSUM_H
