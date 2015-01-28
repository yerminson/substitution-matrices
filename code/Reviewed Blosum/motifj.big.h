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
=========================================================================*/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "math.h"

#define VERSION            8    /*  motifj version number */
#define YES                1
#define NO                 0
#define ESC               27
#define CR                13
#define LF                10
#define UMIN(x, y)    ( (x<y) ? x : y)   /* UNIX min macro */
#define UMAX(x, y)    ( (x>y) ? x : y)   /* UNIX max macro */
/*   INDEX & INDEXCOL compute the sequential indices for the lower half of an
     nxn symmetric matrix given row & column coordinates.  Lower half has
     n(n-1)/2 entries; col=0,n-2 and row=col+1,n-1; col has n-col-1 rows.
     Index runs from 0 to n(n-1)/2 - 1 down columns with
     (index=0)==(col=0,row=1) and (index=n(n-1)/2-1)==(col=n-2,row=n-1).  */
#define INDEXCOL(n, col)    ( col*n - (col*(col+1))/2 )
#define INDEX(n, col, row)  ( col*n - (col*(col+3))/2 - 1 + row )

#define randomize()       srand((unsigned)time(NULL))  /* Seed rand() */

#define MAX_DISTANCE      24  /* Max spacing between aminos of motif */
#define MIN_DISTANCE      1     /* Min distance specification */
#define MAXSEQS       1000000   /* Max number of sequences to be analyzed */
#define MINSEQS           2     /* Min number of sequences to be analyzed */
#define MAXFREQ           1000  /* Max occurences of motif in all seqs */
#define MAX_LENGTH    5500  /* Max length of each sequence */
#define MIN_DOMAIN_WIDTH  10  /* Minimum width */
#define MAX_DOMAIN_WIDTH  55  /* Maximum width */
#define MAX_MERGE_WIDTH   55  /* Max. width of merged blocks */
#define RELEVANT_MOTIFS   50  /* Only top scoring motifs are retained */
#define MAX_MOTIFS    100 /* Buffer motifs before discarding */
#define MINSCORE          1     /* Min block trimming column score (0-2500)*/
#define CLTHRES            80   /* Clustering identity percentage (0-100)*/
#define DROPSCORE         -10   /* Default std devs *10 for dropping block */
#define MOTAUTO4            3   /* max. # motifs for run type 4  */
#define MOTAUTO3      6   /* min. # motifs for run type 3 */
#define MAXBLK             15   /* max # blocks for shotgun assembly */
#define MAXTITLE     75   /* max sequence title length */

#define PROTEIN_SUBDIRECTORY "pros/"  /* Subdirectory containing proteins */
#define PROTEIN_EXTENSION    ".pro"    /* Extension for all protein files */
#define READ                 "r"       /* Code to read disk files */
#define SNAMELEN              11       /* Max length of sequence name */
#define IDLEN                 10       /* Max length of db id */
#define FNAMELEN              80       /* Max length of file name */
#define MAXLINE               240      /* Max line length for ASCII file */
#define MATSIZE               21       /* Scoring matrix dimension */
#define HIGHPASS              4        /* Default high pass filter value */

#define MAXDB 6       /* Max. # database formats */
#define GB 0        /* GenBank type */
#define PIR 1       /* PIR type */
#define EMBL 2        /* EMBL type */
#define UNI 3       /* UNIVERSAL type */
#define VMS 4       /* PIR/VMS type */
#define PROC 5        /* Proclass PCFam file */

#define round(x) (x<0?ceil((x)-0.5):floor((x)+0.5))



/* Declare new data types */
typedef unsigned char *aa_type[20][20][MAX_DISTANCE];

/* Structure to store information about each motif: */
/*  NOTE: integer fields defined as unsigned char to save space;
    they must not exceed 255 in value */
struct motif_struct {
  unsigned char aa1, aa2, aa3, distance1, distance2;
  /*  freq is the number of sequences with this motif */
  int freq, dups;
  /*   seq_no[freq] lists the sequence numbers that have the motif,
       pos[freq] lists the offset of the motif in the corresponding
       sequences; so pos[x] is the offset into sequence # seq_no[x], 
       NOT into sequence # x */
  int seq_no[MAXFREQ], pos[MAXFREQ];
  int score, scores[MAX_DOMAIN_WIDTH], domain, mots;
  char group, sub_group;
  };

/* Structure to store information about groups for motif map routine: */
struct group_struct {
  int group_no, sub_no, position;
  };

/*-------------------------------------------------------------------*/
/*  merged_motif is an array of motifs.  Each entry is one or more   */
/*    motifs.  Each can be thought of as a block of sequences aligned*/
/*    around all the motifs.                                         */
/*-------------------------------------------------------------------*/
struct merged_motif {
  int dropped;      /* YES if dropped */
  char aa[3];     /* Amino acid motif */
  int nmotif;     /* number of motifs, >= 0  */
          /* 0 => merged (inactive) */
  int nident;     /* number of identities */
          /* ..occurs >=SIGNIF in a col. */
  int max_score;      /* max. score of merged motifs */
  int domain;     /* displacement of merged motif blocks */
          /* block width = domain+1 */
  int distance;     /* width of motifs within block */
  int dups;     /* total # dups in all seqs */
  int loffset;      /* 1st position of motif within block*/
  int leftpos[MAXSEQS];   /* leftmost position of motif   */
          /*... for each sequence*/
  int cluster[MAXSEQS];   /* cluster number for each seq */
  int scores[MAX_MERGE_WIDTH];  /* column scores */
  int t_loffset;      /* 1st position of motif within */
          /*  trimmed block */
  int t_domain;     /* trimmed block width=t_domain+1*/
  int t_score;      /* score over trimmed block */
  int position[MAXSEQS];    /* position of block within each seq*/
  int maxpos;     /* max position in any seq */
  int minpos;     /* min position in any seq */
  int in_degree;      /* in-degree of block in DAG */
  int out_degree;     /* out-degree of block in DAG */
};

/*------------------------------------------------------------------------*/
/*   sequences contains all information about the sequences.              */
/*------------------------------------------------------------------------*/
struct sequences {
  int num;    /* number of sequences */
  int totlen;   /* total length of all sequences */
  int *len;   /* lengths of each sequence */
  int *offlen;    /* offset to start of each sequence */
  char *name;   /* 10 char name of each sequence */
  char *seq;    /* sequence bases */
};

/*------------------------------------------------------------------------*/
/*   aux_seq is a list of blocks ordered left to right within a sequence  */
/*  Each sequence has the same number of blocks, but the blocks may be
     arranged in different orders in each sequence.                       */
/*------------------------------------------------------------------------*/
struct aux_seq {
  int block[RELEVANT_MOTIFS]; /* index of block in each position*/
};

struct temp {     /* temporary structure for sorting */
  int value;
  int index;
  int flag;
};

/*-----------------------------------------------------------------------*/
/*    Structure for pairs of sequences.                                  */
/*     pair should be allocated as an array, & the number of the         */
/*     sequences forming the pair inferred from the array index.         */
/*-----------------------------------------------------------------------*/
struct pair {
  int score;    /* # of identities within trimmed block */
  int cluster;    /* cluster # for this pair */
};

/*-----------------------------------------------------------------------*/
/*  block_list is a list of blocks in an order that is consistent among
    all sequences;  it implies a partial multiple alignment.
    The next_block list is a list of all blocks, some of which may
    overlap.  The doubly linked next_best/prev_best list is a subset
    consisting of the best non-overlapping blocks in the list.           */
/*-----------------------------------------------------------------------*/
struct block_list {
  int b;      /* index of block (merged_motif) */
  int minprev;    /* min. distance from previous block */
  int maxprev;    /* max. distance from previous block */
  struct block_list *next_block;  /* all blocks in the list  */
  struct block_list *next_best;   /* best blocks in the list */
  struct block_list *prev_best;
};

/*------------------------------------------------------------------------*/
/*  path is a list of all possible paths through the blocks in all seqs.
     The first_block list includes all blocks, including those that
     possible overlap.  The first_best list includes the best sub-path
     of non-overlapping blocks.                                           */
/*------------------------------------------------------------------------*/
struct path {
  int nblocks;    /* # of blocks in path */
  int nbest;    /* # of blocks in best sub-path */
  int naas;   /* # of AAs in best sub-path */
  unsigned long totscore; /* sum of scores of blocks in best sub-path*/
  int totmotif;   /* # motifs in best sub-path*/
  int totident;   /* # conserved residues in best sub_path */
  int nseqs;    /* # seqs in best sub-path */
  int seqs[MAXSEQS];  /* YES if path holds for sequence */
  struct block_list *first_block;  /* first block in path */
  struct block_list *first_best;   /* first block in best sub-path*/

  struct path *next_path;
};

/*----------------------------------------------------------------------*/
/*  Ajacency matrix representation of distances between blocks in
    all sequences.                                                      */
/*----------------------------------------------------------------------*/
struct matrix {
  int npos;   /* # of seqs with positive diff in cell */
  int maxdiff;    /* maximum positive difference in cell  */
/*  int dist[MAXSEQS];   distance from row to col for seq s */
  int mark;   /* all-purpose flag */
};

struct follow_data {
  struct path *path;
  int mark[RELEVANT_MOTIFS][RELEVANT_MOTIFS];
};

/* -------- Scoring matrix structure ----------------------------------*/
struct score {
  char scores[MATSIZE][MATSIZE];  /* valid range -127 to +128 */
  int highpass;     /* high pass filter value */
};

/*--------Sequence database format structure -------------------------*/
struct db_info {                             /*  Database format info */
   char *type, *start, *desc, *seq, *end;
   int title_offset, seq_offset;
};
/*-Structure to split up a file name of the form <directory>\<name>.ext -*/
struct split_name {
  int dir_len, file_len, name_len;
};
/*------ Structure to hold the contents of a .lis or .lst file ------*/
struct db_id {
   char entry[IDLEN+1];   /* sequence name */
   char pir[IDLEN+1];   /* PIR entry name */
   char ps[2];      /* PS type=T, F or P */
   int len;     /* sequence length */
   int frag;      /* YES if seq is a fragment */
   int lst;     /* seq in .lst file */
   int found;     /* seq found in database */
   int block;     /* seq found in block */
   int search;      /* used by excluded.c => use seq for search*/
   int rank;      /* used by matodat.c, fastodat.c */
   int score;     /* used by matodat.c, fastodat.c */
   double pvalue;   /* P-value */
   struct db_id *next;
   struct db_id *prior;
};
/*  For pb_weights(), from BLIMPS  */
struct pb_counts {    /* for pb_weights() */
  double diffaas;   /* # of different aas in position */
  double naas[MATSIZE];   /* # of each aa in position */
};
