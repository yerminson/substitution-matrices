#ifndef CONFIGURATIONS_H
#define CONFIGURATIONS_H



#define VERSION            8    /*  motifj version number */
#define YES                1
#define NO                 0
#define ESC               27
#define CR                13
#define LF                10
#define UMIN(x, y)    ( (x<y) ? x : y)   /* UNIX min macro */
#define UMAX(x, y)    ( (x>y) ? x : y)   /* UNIX max macro */

/*
  INDEX & INDEXCOL compute the sequ ential indices for the lower half of an
  nxn symmetric matrix given row & column coordinates.  Lower half has
  n(n-1)/2 entries; col=0,n-2 and row=col+1,n-1; col has n-col-1 rows.
  Index runs from 0 to n(n-1)/2 - 1 down columns with
  (index=0)==(col=0,row=1) and (index=n(n-1)/2-1)==(col=n-2,row=n-1).
*/
#define INDEXCOL(n, col)    ( col*n - (col*(col+1))/2 )
#define INDEX(n, col, row)  ( col*n - (col*(col+3)) /2 - 1 + row )

#define randomize()       srand((unsigned)time(NULL))  /* Seed rand() */

#define MAX_DISTANCE      24  /* Max spacing between aminos of motif */
#define MIN_DISTANCE      1     /* Min distance specification */
#define MAXSEQS       16583   /* Max number of sequences to be analyzed */
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


#define AAS 20
#define MINSTR 0        /* Min. strength to count */
#define MAXSTR 9999       /* Max. strength to count */
#define PBTOTAL 4.0       /* Total column weight for position-based */
// IMPORTANT CODING CHANGE: Note that in motifj.big.h, you must change
// MAX_MERGE_WIDTH to be the maximum block width for your database.

/*-----------------------------------------------------------------------*/
/*    Structure for pairs of sequences.                                  */
/*     pair should be allocated as an array, & the number of the         */
/*     sequences forming the pair inferred from the array index.         */
/*-----------------------------------------------------------------------*/
struct pair {
  int score;    /* # of identities within trimmed block */
  int cluster;    /* cluster # for this pair */
};


class configurations
{


public:
  configurations();
};

#endif // CONFIGURATIONS_H
