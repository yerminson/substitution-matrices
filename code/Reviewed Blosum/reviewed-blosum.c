/*=======================================================================
    BLOSUM: (C) Copyright 1992, Fred Hutchinson Cancer Research Center
    Program to compute "blosum" scoring matrix from blocks database.

    blosum blocks.dat minstr maxstr cluster scale

    blocks.dat = a blocks database
    minstr = minimum strength of blocks to count (0-9999)
    maxstr = maximum strength of blocks to count (0-9999)
    cluster = re-clustering threshold (n or e or 0-100)
        If "n", does not cluster sequences.
        If "w", uses weights as found in blocks.
        If "pn", computes position-based weights and uses
            1/n as the parameter.
        If "e", uses clusters as found in blocks (clusters of
        sequences are separated by a blank line). If 0-100,
        clusters sequences based on that percentage of pair-
        wise identities.
    scale = scale for matrix;
        0: Let program decide based on entropy
        n: 1/n bits

   Count all pairwise substitutions in all blocks in blocks.dat
   with strength between minstr and maxstr. Pairs are only counted
   for sequences in different clusters.

    Output includes:
        Standard output - various reports, including raw frequencies.
        blosumx.qij file - target frequencies (odds ratios), x=clustering %
        blosumx.sij file - scoring matrix in bits
        blosumx.iij file - scoring matrix rounded to integers and scaled
                             according to relative entropy
--------------------------------------------------------------------------
      9/23/91  J. Henikoff
>>>>>>>>>>>>>>>>>>>> Blocks 5.0 <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      9/18/92  Changed scaling to variable bit fractions as in NCBI pam.c.
               Changed sij[][] from int to double.
               Changed positive matrix to +minsij instead of +8.
      9/24/92  Added printout of floating point sij in bit units.
      9/26/92  Write .qij, .sij and .iij output files.
               Compute B, Z and X entries as weighted averages.
     10/17/92  Added comments. Removed Risler routine.
      12/1/92  Added B,Z,X columns to positive matrix.
      1/14/93  Replaced MAXNSEQ with MAXSEQS from motifj.h.
      1/18/93  Fixed problem with BZ value.
       6/7/93  Added scale option. 
      9/22/93  Added B, Z and X to .sij file
      9/29/93  Fixed problem reading blocks db with sequence weights
>>>>>>>>>>>>>>>>>>>>>  Blocks 7.x <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      1/27/94  Added option to make a single matrix using sequence weights
      4/18/94  Normalize existing sequence weights (divide by total).
      4/20/94  Changed count to be product of sequence weights.
      4/23/94  Added option for position-based weight series.
      4/27/94  Separated total pairs from total weight (FTotPairs, FTotWeight).
      4/28/94  Modified position-based weight option.
      5/ 2/94  For position-based weight, assume total column weight is 4.0
               instead of 1.0 to get a wider range of matrices.
      7/26/94   Skip over tabs in blocks (fill_block)
>>>>>>>>>>>>>>>>>>>>> Blocks 8 <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      8/29/94  Report average number of clusters per block, Block.nclus
      8/31/94   Report weighted average of clusters per block, TotClump
==========================================================================*/

#include "motifj.big.h"
#include <math.h>

#define AAS 20
#define MINSTR 0        /* Min. strength to count */
#define MAXSTR 9999       /* Max. strength to count */
#define PBTOTAL 4.0       /* Total column weight for position-based */
// IMPORTANT CODING CHANGE: Note that in motifj.big.h, you must change
// MAX_MERGE_WIDTH to be the maximum block width for your database. 

/*---------------------Functions----------------------------------*/
int read_dat();
void fill_block();
void count_block();
void count_weight();
void count_cluster();
void count_position();
void cluster_seqs();
int aachar_to_num();
char *num_to_aachar();
/*---------------------Functions in motmisc.o ----------------------*/


struct block {                          /* Block structure */
   char ac[10];
   int nseq, width, strength, nclus;
   int aa[MAXSEQS][MAX_MERGE_WIDTH];    /* aas for seq */
   double weight[MAXSEQS];              /* seq weights found in block */
   int cluster[MAXSEQS];                /* cluster # for seq */
   int ncluster[MAXSEQS];               /* #seqs in same cluster */
   double totdiag, totoffd, wtot;
} Block;


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

/*=====================================================================*/
void
main(argc, argv)
  int argc;
  char *argv[];
{
  FILE *fdat, *fout;
  char datfile[FNAMELEN], outfile[FNAMELEN], ctemp[6];
  unsigned long totaas;
  int totblk, row, col, tij[AAS][AAS], maxsij, minsij, sumsij;
  int iscale, is, itemp;
  double sij[AAS + 3][AAS + 3], totpairs, totdiag, totoffd;
  double dtemp, dtemp1, dtemp2, x, xx;
  double ftotpairs, ftotdiag, ftotoffd;
  double s, fij, fifj, oij, entropy, expected;
  double sumpi, sumqij;

  printf("BLOSUM: (C) Copyright 1992, Fred Hutchinson Cancer");
  printf(" Research Center\n");

   /*------------Arg 1, blocks database name ----------------*/
  if (argc > 1)
    strcpy(datfile, argv[1]);
  else {
    printf("Enter name of blocks database:\n");
    gets(datfile);
  }
  if ((fdat = fopen(datfile, "r")) == NULL) {
    printf("Cannot open %s\n", datfile);
    exit(-1);
  } else
    printf("Reading %s\n", datfile);
   /*----------Arg 2, minimum block strength----------------------*/
  if (argc > 2)
    MinStr = atoi(argv[2]);
  else {
    MinStr = MINSTR;
    printf("Enter minimum block strength [%d]: ", MinStr);
    gets(ctemp);
    if (strlen(ctemp))
      MinStr = atoi(ctemp);
  }
   /*-----------Arg 3, Maximum block strength----------------------*/
  if (argc > 3)
    MaxStr = atoi(argv[3]);
  else {
    MaxStr = MAXSTR;
    printf("Enter maximum block strength [%d]: ", MaxStr);
    gets(ctemp);
    if (strlen(ctemp))
      MaxStr = atoi(ctemp);
  }
  printf("Minimum block strength=%d, Maximum block strength=%d\n", MinStr, MaxStr);

   /*----------Arg 4, Clustering ------------------------------------*/
  Cluster = -1;       /* No clustering */
  if (argc > 4)
    strcpy(ctemp, argv[4]);
  else {
    printf("Enter n for no clustering\n");
    printf("   or a number between 0 and 100 for percent identity clustering\n");
    printf("   or e for existing (implied) clusters\n");
    printf("   or w for existing sequence weights\n");
    printf("   or pn for position-based weights, PB weight = 1/n\n");
    printf("Enter clustering identity percentage or n/e/w/pn [e]: ");
    gets(ctemp);
  }
  /*
     e => -1, n=> -2, w=> -3 p=> -4 
   */
  Cluster = -1;
  if (strlen(ctemp))
    if (ctemp[0] == 'n' || ctemp[0] == 'N')
      Cluster = -2;
    else if (ctemp[0] == 'w' || ctemp[0] == 'W')
      Cluster = -3;
    else if (ctemp[0] == 'p' || ctemp[0] == 'P') {
      Cluster = -4;
      itemp = atoi(ctemp + 1);
      if (itemp > 0)
        PBParameter = sqrt((double) 1.0 / itemp);
      else
        PBParameter = 1.0;
    } else if (ctemp[0] != 'e' && ctemp[0] != 'E')
      Cluster = atoi(ctemp);
  if (Cluster > 100)
    Cluster = -1;
  if (Cluster == -1)
    printf("Existing clustering will be used\n");
  else if (Cluster == -2)
    printf("No clustering will be used\n");
  else if (Cluster == -3)
    printf("Sequence weights will be used\n");
  else if (Cluster == -4) {
    printf("Position-based weights will be used");
    printf(" with parameter = %.2f\n", PBParameter);
  } else
    printf("Re-clustering percentage = %d\n", Cluster);

   /*-----------Arg 5, Scale --------------------------------------*/
  if (argc > 5)
    iscale = atoi(argv[5]);
  else {
    iscale = 0;
    printf("Enter scale n for 1/n bits [%d]: ", iscale);
    gets(ctemp);
    if (strlen(ctemp))
      iscale = atoi(ctemp);
  }
  if (iscale < 0)
    iscale = 0;
  if (iscale > 100)
    iscale = 100;
  if (iscale == 0)
    printf("Scale based on relative entropy\n");
  else
    printf("Requested scale = 1/%d bits\n", iscale);

   /*----------------Initialize--------------------------------------*/
  TotPairs = TotSeqs = TotWidth = TotAas = (unsigned long) 0;
  TotBlk = TotSeg = (unsigned long) 0;
  TotClump = FTotPairs = FTotWeight = 0.0;
  sumsij = 0;
  minsij = 999;
  maxsij = -999;
  for (row = 0; row < AAS; row++) {
    AaPairs[row] = (unsigned long) 0;
    AaFreq[row] = (unsigned long) 0;
    FAaPairs[row] = (double) 0.0;
    for (col = 0; col < AAS; col++) {
      Counts[row][col] = sij[row][col] = (double) 0.0;
    }
  }

   /*---------------Read the blocks ----------------------------------*/
  totblk = read_dat(fdat);
  fclose(fdat);
  printf("\n%d blocks processed", totblk);
  printf(", %ld blocks contributed pairs to matrix\n", TotBlk);
  if (Cluster >= 0) {
    printf(" %lf clumps contributed pairs to matrix (%f)\n", TotClump, TotClump / FTotWeight);
    printf(" %ld segments contributed pairs to matrix (%f)\n", TotSeg, (float) TotSeg / TotBlk);
  }
  printf("%f total pairs, %f total weight\n", FTotPairs, FTotWeight);
  printf(" %ld total sequences, ", TotSeqs);
  printf("%ld total columns, ", TotWidth);
  printf("%ld total AAs\n\n", TotAas);
/*----------------------------------------------------------------------*/
/*-------------------The frequency matrix ------------------------------*/
  printf("Frequencies = fij pairs (off-diagonals = 2*fij):\n");
  for (col = 0; col < AAS; col++)
    printf("    %1s  ", num_to_aachar(col));
  printf("\n");
  totdiag = totoffd = totpairs = (double) 0.0;
  ftotdiag = ftotoffd = ftotpairs = (double) 0.0;
  for (row = 0; row < AAS; row++) {
    for (col = 0; col <= row; col++) {
      if (row == col) {
        printf("%9.2f ", Counts[row][col]);
        totdiag += Counts[row][col];
      } else {
        printf("%9.2f ", (Counts[row][col] + Counts[col][row]));
        totoffd += Counts[row][col] + Counts[col][row];
      }
    }
    printf("\n");
  }
  totpairs = totdiag + totoffd;
  ftotpairs = ftotdiag + ftotoffd;
  for (col = 0; col < AAS; col++)
    printf("   %1s  ", num_to_aachar(col));
  printf("\n");
/*--------------Print the target frequencies, qij ----------------------*/
  if (Cluster >= 0)
    sprintf(outfile, "blosum%d.qij", Cluster);
  else if (Cluster == -4)
    sprintf(outfile, "blosump%d.qij", (int) (1. / PBParameter));
  else if (Cluster == -3)
    sprintf(outfile, "blosumw.qij");
  else if (Cluster == -2)
    sprintf(outfile, "blosumn.qij");
  else if (Cluster == -1)
    sprintf(outfile, "blosume.qij");
  if ((fout = fopen(outfile, "wt")) == NULL)
    fout = stdout;
  printf("\nTarget Probabilities=qij in %s\n", outfile);
  fprintf(fout, "#  BLOSUM Clustered Target Frequencies=qij\n");
  fprintf(fout, "#  Blocks Database = %s\n", datfile);
  if (Cluster >= 0)
    fprintf(fout, "#  Cluster Percentage: >= %d\n", Cluster);
  else if (Cluster == -4)
    fprintf(fout, "#  Position-based Clustering Parameter: <= %.4f\n", PBParameter);
  else if (Cluster == -3)
    fprintf(fout, "#  Explicit Sequence Weights Used\n");
  else if (Cluster == -2)
    fprintf(fout, "#  No Sequence Clustering Used\n");
  else if (Cluster == -1)
    fprintf(fout, "#  Existing Clusters Used\n");
  for (col = 0; col < AAS; col++)
    fprintf(fout, "   %1s   ", num_to_aachar(col));
  fprintf(fout, "\n");
  sumqij = 0.0;
  for (row = 0; row < AAS; row++) {
    for (col = 0; col <= row; col++) {
      if (col == row) {
        fij = Counts[row][col];
        sumqij += fij / totpairs;
      } else {
        fij = (Counts[row][col] + Counts[col][row]) / 2.0;
        sumqij += 2.0 * fij / totpairs;
      }
      fprintf(fout, "%.4f ", (double) fij / totpairs);
    }
    fprintf(fout, "\n");
  }
  fclose(fout);
  /*
     for (col=0; col<AAS; col++) printf(" %1s ", num_to_aachar(col)); 
   */
  printf("sumqij = %f\n", sumqij);
/*----------------Marginal frequencies -------------------------------*/
  printf("\nMarginal Probabilities by AA = p(i,*):\n");
  for (col = 0; col < AAS; col++) {
    for (row = 0; row < AAS; row++)
      if (col == row)
        FAaPairs[row] += Counts[row][col];
      else
        FAaPairs[row] += (Counts[row][col] + Counts[col][row]) / 2.0;
  }
  sumpi = 0.0;
  for (col = 0; col < AAS; col++)
    printf("   %1s  ", num_to_aachar(col));
  printf("\n");
  for (col = 0; col < AAS; col++) {
    printf("%.3f ", FAaPairs[col] / totpairs);  /* p(col,*) */
    sumpi += FAaPairs[col] / totpairs;
  }
  printf("\nsumpi=%.3f", sumpi);

  totaas = (unsigned long) 0;
  printf("\nAA Frequencies = ai:\n");
  for (col = 0; col < AAS; col++) {
    printf("%8ld ", AaFreq[col]);
    totaas += AaFreq[col];
  }
/*---------------Amino acid frequencies ----------------------------*/
  printf("\nAA Probabilities = fi:\n");
  sumpi = 0.0;
  for (col = 0; col < AAS; col++)
    printf("   %1s  ", num_to_aachar(col));
  printf("\n");
  for (col = 0; col < AAS; col++) {
    printf("%.3f ", (double) AaFreq[col] / totaas); /* fi */
    sumpi += (double) AaFreq[col] / totaas;
  }
  printf("\nsumpi=%.3f", sumpi);

  printf("\n\ntotpairs=%.3f, FTotWeight=%.3f", totpairs, FTotWeight);
  printf("\n totdiag=%.3f, totoffd=%.3f", totdiag, totoffd);
  printf(" totaas=%ld, TotAas=%ld\n", totaas, TotAas);

/*-----------------The bit matrix--------------------------------*/
  entropy = expected = 0.0;
  ftotpairs = totpairs;
  for (row = 0; row < AAS; row++) {
    /*
       fij = pair frequencies = f(i,j) 
     */
    /*
       q(i,j) = fij/ftotpairs 
     */
    /*
       fifj = ftotpairs*p(i,*) * ftotpairs*p(*,j) 
     */
    /*
       p(i,*)*p(*,j) = fifj/(ftotpairs*ftotpairs) 
     */
    /*
       oij = odds ratio = q(i,j)/(p(i,*) * p(*,j)) 
     */
    /*
       s = score in bits = s(i,j) 
     */
    for (col = 0; col <= row; col++) {
      if (row == col)
        fij = (double) Counts[row][col];
      else
        fij = (double) (Counts[row][col] + Counts[col][row]) / 2.0;
      /*
         FAaPairs[i] = f(i,*) 
       */
      fifj = (double) FAaPairs[row] * FAaPairs[col];
      if (fifj > 0.000000001)
        oij = ftotpairs * fij / fifj;
      else
        oij = 0.0;
   /*--- Log odds ratio ----*/
      if (oij > 0.000000001)
        s = log((double) oij);
      else
        s = -20.0;    /* minus infinity */
   /*--- Round off for the log base 10 matrix ----*/
      dtemp = 10.0 * s / log(10.0);
                  /*-- log base 10 */
      tij[row][col] = round(dtemp);
   /*--- Log base 2 (bit) matrix ---*/
      s /= log(2.0);    /* log base 2 = bits */
   /*--- compute entropy & expected value in bits */
      if (row == col) {
        entropy += (double) s *fij / ftotpairs;
        expected += (double) s *fifj / (ftotpairs * ftotpairs);
      } else {
        entropy += (double) 2.0 *s * fij / ftotpairs;
        expected += (double) 2.0 *s * fifj / (ftotpairs * ftotpairs);
      }
      sij[row][col] = s;
    }
  }
  printf("\nEntropy=%.4f bits, expected=%.4f bits\n", entropy, expected);

/*----------------Fill in symmetric part of matrix----------------*/
  for (row = 0; row < AAS; row++)
    for (col = 0; col <= row; col++)
      if (row != col)
        sij[col][row] = sij[row][col];

/*------- Compute the B, Z and X columns ------------------------*/
   /*--------------Compute values for B (20) as the weighted average of
    N (2) and D (3) ----------------*/
   /*--------------Compute values for Z (21) as the weighted average of
    Q (5) and E (6) ----------------*/
  for (col = 0; col < AAS; col++) {
    dtemp = (double) (AaFreq[2] * sij[col][2] + AaFreq[3] * sij[col][3]);
    dtemp /= (AaFreq[2] + AaFreq[3]);
    sij[col][20] = sij[20][col] = dtemp;
    dtemp = (double) (AaFreq[5] * sij[col][5] + AaFreq[6] * sij[col][6]);
    dtemp /= (AaFreq[5] + AaFreq[6]);
    sij[col][21] = sij[21][col] = dtemp;
  }
   /*----- value for BB is weighted average of all N x D combos ----*/
  dtemp = (double) AaFreq[2] * AaFreq[2] * sij[2][2];
  dtemp += (double) AaFreq[2] * AaFreq[3] * sij[2][3];
  dtemp += (double) AaFreq[3] * AaFreq[2] * sij[3][2];
  dtemp += (double) AaFreq[3] * AaFreq[3] * sij[3][3];
  dtemp1 = (double) AaFreq[2] + AaFreq[3];
  dtemp /= (dtemp1 * dtemp1);
  sij[20][20] = dtemp;
  // IMPORTANT CODING CHANGE, implemented by Henikoff... fixed
  // an overflow issue for calculating BZ (ZB) values.
   /*--- score for ZB is weighted average of all (N,D)x(Q,E) combos------*/
  dtemp = (double) AaFreq[2] * AaFreq[5] * sij[2][5];
  dtemp += (double) AaFreq[2] * AaFreq[6] * sij[2][6];
  dtemp += (double) AaFreq[3] * AaFreq[5] * sij[3][5];
  dtemp += (double) AaFreq[3] * AaFreq[6] * sij[3][6];
  dtemp2 = (double) AaFreq[2] + AaFreq[3];
  dtemp2 *= (double) AaFreq[5] + AaFreq[6];
  dtemp /= dtemp2;
  dtemp1 = (double) AaFreq[5] * AaFreq[2] * sij[5][2];
  dtemp1 += (double) AaFreq[5] * AaFreq[3] * sij[5][3];
  dtemp1 += (double) AaFreq[6] * AaFreq[2] * sij[6][2];
  dtemp1 += (double) AaFreq[6] * AaFreq[3] * sij[6][3];
  dtemp1 /= dtemp2;
  dtemp *= dtemp1;
  sij[21][20] = sij[20][21] = dtemp;
   /*---- value for ZZ is weighted average of all Q x E combos -----*/
  dtemp = (double) AaFreq[5] * AaFreq[5] * sij[5][5];
  dtemp += (double) AaFreq[5] * AaFreq[6] * sij[5][6];
  dtemp += (double) AaFreq[6] * AaFreq[5] * sij[6][5];
  dtemp += (double) AaFreq[6] * AaFreq[6] * sij[6][6];
  dtemp1 = (double) AaFreq[5] + AaFreq[6];
  dtemp /= (dtemp1 * dtemp1);
  sij[21][21] = dtemp;
   /*-----------Compute values for unknown entry X (22) as the average
    of row scores weighted by frequency ----------------------*/
  xx = dtemp = 0.0;
  for (row = 0; row < AAS; row++) {
    x = 0.0;
    for (col = 0; col < AAS; col++) {
      x += (double) AaFreq[col] * sij[row][col];
      xx += (double) AaFreq[row] * AaFreq[col] * sij[row][col];
      dtemp += (double) AaFreq[row] * AaFreq[col];
    }
    x /= (double) totaas;
    sij[row][22] = sij[22][row] = x;
  }
  dtemp1 = xx / dtemp;
  sij[22][22] = dtemp1;
   /*--------Now fill in (X) x (B,Z) ------------------------------*/
  x = xx = 0.0;
  for (row = 0; row < AAS; row++) {
    // IMPORTANT CODING CHANGES: Two, actually. First, the sij
    // column index in xx should match the AaFreq index.  They
    // initially transposed them.
    // Second, they still haven't fixed the possibility of
    // an overflow error.  They need to typecase each
    // separate multiplication grouping due to order of
    // operations.
    x += (double) AaFreq[row] * AaFreq[2] * sij[row][2] +
      (double) AaFreq[row] * AaFreq[3] * sij[row][3];
    xx += (double) AaFreq[row] * AaFreq[5] * sij[row][5] +
      (double) AaFreq[row] * AaFreq[6] * sij[row][6];
    // x += (double) AaFreq[row]*AaFreq[2]*sij[row][2] +
    // AaFreq[row]*AaFreq[3]*sij[row][3];
    // xx += (double) AaFreq[row]*AaFreq[5]*sij[row][6] +
    // AaFreq[row]*AaFreq[5]*sij[row][6];
  }
  x /= (double) (AaFreq[2] + AaFreq[3]);
  x /= (double) totaas;
  xx /= (double) (AaFreq[5] + AaFreq[6]);
  xx /= (double) totaas;
  sij[20][22] = sij[22][20] = x;
  sij[21][22] = sij[22][21] = xx;

/*----Print the sij matrix in bit units ----------------*/
  if (Cluster >= 0)
    sprintf(outfile, "blosum%d.sij", Cluster);
  else if (Cluster == -4)
    sprintf(outfile, "blosump%d.sij", (int) (1. / PBParameter));
  else if (Cluster == -3)
    sprintf(outfile, "blosumw.sij");
  else if (Cluster == -2)
    sprintf(outfile, "blosumn.sij");
  else if (Cluster == -1)
    sprintf(outfile, "blosume.sij");
  if ((fout = fopen(outfile, "wt")) == NULL)
    fout = stdout;
  printf("\nScoring matrix in bit units=sij in %s\n", outfile);
  fprintf(fout, "#  BLOSUM Clustered Scoring Matrix in Bit Units=sij\n");
  fprintf(fout, "#  Blocks Database = %s\n", datfile);
  if (Cluster >= 0)
    fprintf(fout, "#  Cluster Percentage: >= %d\n", Cluster);
  else if (Cluster == -4)
    fprintf(fout, "#  Position-based Clustering Parameter: <= %.4f\n", PBParameter);
  else if (Cluster == -3)
    fprintf(fout, "#  Explicit Sequence Weights Used\n");
  else if (Cluster == -2)
    fprintf(fout, "#  No Sequence Clustering Used\n");
  else if (Cluster == -1)
    fprintf(fout, "#  Existing Clusters Used\n");
  fprintf(fout, "#  Entropy = % 8.4f, Expected = % 8.4f\n", entropy, expected);
  for (col = 0; col < AAS; col++)
    fprintf(fout, "   %1s     ", num_to_aachar(col));
  fprintf(fout, "B   Z   X\n");
  fprintf(fout, "\n");
  for (row = 0; row < AAS + 3; row++) {
    for (col = 0; col <= row; col++)
      fprintf(fout, "% 8.4f ", sij[row][col]);
    fprintf(fout, "\n");
  }
  fclose(fout);

/*-------------Determine the scale based on entropy-------------------*/
  if (iscale == 0) {
    dtemp = 2.0 / sqrt(entropy);
    iscale = round(dtemp);
    if (iscale < 2)
      iscale = 2;
  }

/*-------------Print the integer matrix-------------------------------*/
  if (Cluster >= 0)
    sprintf(outfile, "blosum%d.iij", Cluster);
  else if (Cluster == -4)
    sprintf(outfile, "blosump%d.iij", (int) (1. / PBParameter));
  else if (Cluster == -3)
    sprintf(outfile, "blosumw.iij", Cluster);
  else if (Cluster == -2)
    sprintf(outfile, "blosumn.iij");
  else if (Cluster == -1)
    sprintf(outfile, "blosume.iij");
  if ((fout = fopen(outfile, "wt")) == NULL)
    fout = stdout;
  printf("\nInteger scoring matrix in 1/%d bit units in %s\n", iscale, outfile);
  fprintf(fout, "#  BLOSUM Clustered Scoring Matrix in 1/%d Bit Units\n", iscale);
  fprintf(fout, "#  Blocks Database = %s\n", datfile);
  if (Cluster >= 0)
    fprintf(fout, "#  Cluster Percentage: >= %d\n", Cluster);
  else if (Cluster == -4)
    fprintf(fout, "#  Position-based Clustering Parameter: <= %.4f\n", PBParameter);
  else if (Cluster == -3)
    fprintf(fout, "#  Explicit Sequence Weights Used\n");
  else if (Cluster == -2)
    fprintf(fout, "#  No Sequence Clustering Used\n");
  else if (Cluster == -1)
    fprintf(fout, "#  Existing Clusters Used\n");
  fprintf(fout, "#  Entropy = % 8.4f, Expected = % 8.4f\n", entropy, expected);
  for (col = 0; col < AAS; col++)
    fprintf(fout, " %1s  ", num_to_aachar(col));
  fprintf(fout, "B   Z   X\n");
  for (row = 0; row < AAS + 3; row++) {
    for (col = 0; col <= row; col++) {
      s = (double) sij[row][col] * iscale;
      is = round(s);
      fprintf(fout, "%3d ", is);
      sumsij += is;
      if (is < minsij)
        minsij = is;
      if (is > maxsij)
        maxsij = is;
    }
    fprintf(fout, "\n");
  }
  fclose(fout);
  /*
     for (col=0; col<AAS; col++) printf(" %1s ", num_to_aachar(col)); 
   */
  printf("\nMaximum=%d, minimum=%d, sum=%d\n", maxsij, minsij, sumsij);

/*--------------------Add the high pass value = minsij--------------------*/
  printf("\nPositive matrix (+%d):\n", 0 - minsij);
  for (col = 0; col < AAS; col++)
    printf(" %1s  ", num_to_aachar(col));
  printf("B   Z   X\n");
  for (row = 0; row < AAS + 3; row++) {
    for (col = 0; col <= row; col++) {
      s = (double) sij[row][col] * iscale;
      is = round(s) - minsij;
      printf("%3d ", is);
    }
    printf("\n");
  }
  for (col = 0; col < AAS; col++)
    printf(" %1s  ", num_to_aachar(col));
  printf("B   Z   X\n");
/*-----------------The log base 10 matrix ----------------------------*/
  printf("\n10 times log base 10 matrix:\n");
  for (col = 0; col < AAS; col++)
    printf(" %1s  ", num_to_aachar(col));
  printf("\n");
  for (row = 0; row < AAS; row++) {
    for (col = 0; col <= row; col++)
      printf("%3d ", tij[row][col]);
    printf("\n");
  }
  for (col = 0; col < AAS; col++)
    printf(" %1s  ", num_to_aachar(col));
  printf("\n");
/*----------------------------------------------------------------------*/
  exit(0);
}               /* end of main */
/*=====================================================================
      Read the blocks database
=======================================================================*/
int
read_dat(fdat)
  FILE *fdat;
{
  int totblk;
  char line[MAXLINE], *ptr, *ptr1;

  totblk = 0;
   /*SPSS*/
    /*
       printf("\nBlock diagonal off-diagonal total strength\n"); 
     */
    while (!feof(fdat) && fgets(line, MAXLINE, fdat) != NULL) {
    if (strncmp(line, "AC   ", 5) == 0) {
      strncpy(Block.ac, line + 5, 8);
      Block.ac[9] = '\0';
    } else if (strncmp(line, "BL   ", 5) == 0) {
      Block.strength = 0;
      ptr = strstr(line, "strength=");
      if (ptr != NULL) {
        ptr1 = strtok(ptr, "=");
        ptr1 = strtok(NULL, "\n\t\r");
        Block.strength = atoi(ptr1);

      } else
        Block.strength = MinStr;  /* use block if no strength field */
      if (Block.strength >= MinStr && Block.strength <= MaxStr) {

        fprintf(stderr, "Block %d:%s\n", totblk, Block.ac);
        fflush(stderr);

        fill_block(fdat);

        if (Cluster >= 0)
          cluster_seqs(); /* re-cluster sequences */
        if (Cluster == -2)
          count_block();  /* no clustering */
        else if (Cluster == -3)
          count_weight();
        else if (Cluster == -4)
          count_position();
        else
          count_cluster();  /* clustering */
        totblk++;
      }
    }
  }

  return (totblk);
}               /* end of read_dat */
/*====================================================================
     Fill up the block structure
=======================================================================*/
void
fill_block(fdat)
  FILE *fdat;
{
  int done, i, li, n, cluster, ncluster;
  char line[MAXLINE], *ptr;

  Block.nseq = Block.width = 0;
  Block.totdiag = Block.totoffd = Block.wtot = (double) 0;
  cluster = ncluster = 0;
  done = NO;
  while (!done && !feof(fdat) && fgets(line, MAXLINE, fdat) != NULL) {
    /*fprintf(stderr, "%d\t%s", Block.nseq, line);*/
    /*fflush(stderr);*/
    if (strlen(line) == 1) {  /* blank line => new cluster */
      /*
         Set #seqs in cluster to seqs in previous cluster 
       */
      if (ncluster > 0)
        for (n = 0; n < Block.nseq; n++)
          if (Block.cluster[n] == cluster)
            Block.ncluster[n] = ncluster;
      cluster++;
      ncluster = 0;
    } else if (strlen(line) > 1) {
      if (strncmp(line, "//", 2) == 0)
        done = YES;
      else if (strlen(line) > 20) {
        li = 0;     /* skip over sequence name & offset */
        while (line[li] != ')')
          li++;
        li++;
        while (line[li] == ' ')
          li++;   /* skip over spaces */
        /*
           assuming no spaces within block ! 
         */
        ptr = strtok(line + li, " \r\n\t");
        for (i = 0; i < strlen(ptr); i++) {
          Block.aa[Block.nseq][i] = aachar_to_num(ptr[i]);
          if (Block.aa[Block.nseq][i] >= 0 && Block.aa[Block.nseq][i] < AAS) {
            AaFreq[Block.aa[Block.nseq][i]]++;
            TotAas++;
          }
          /*
             else { printf("\nBad char: %s %c", Block.ac, ptr[i]); } 
           */
        }
        li += strlen(ptr);
        ptr = strtok(NULL, " \r\n\t");
        if (ptr != NULL) {
          Block.weight[Block.nseq] = atof(ptr);
          Block.wtot += Block.weight[Block.nseq];
        } else
          Block.weight[Block.nseq] = 0;
        Block.cluster[Block.nseq] = cluster;
        ncluster++;   /* # seqs in current cluster */
        Block.width = i;
        Block.nseq++;
      }
    }
  }
  /*
     Compute weights for the last cluster 
   */
  if (ncluster > 0)
    for (n = 0; n < Block.nseq; n++)
      if (Block.cluster[n] == cluster)
        Block.ncluster[n] = ncluster;
  TotSeqs += Block.nseq;
  TotWidth += Block.width;
}               /* end of fill_block */
/*===================================================================*/
void
count_block()
{
  int seq1, seq2, col;

  for (col = 0; col < Block.width; col++)
    for (seq1 = 0; seq1 < Block.nseq; seq1++)
      for (seq2 = seq1 + 1; seq2 < Block.nseq; seq2++)
        if (Block.aa[seq1][col] >= 0 && Block.aa[seq1][col] < AAS &&
          Block.aa[seq2][col] >= 0 && Block.aa[seq2][col] < AAS) {
          Counts[Block.aa[seq1][col]][Block.aa[seq2][col]] += 1.0;
          FTotPairs += 1.0;
          FTotWeight += 1.0;
          if (Block.aa[seq1][col] == Block.aa[seq2][col])
            Block.totdiag += 1.0;
          else
            Block.totoffd += 1.0;
        }
  /*
     Output for SPSS job to count pairs by strength 
   */
  TotBlk++;
  /*
     printf("%s %9.2f %9.2f %9.2f %d\n", Block.ac, Block.totdiag, Block.totoffd,
     Block.totoffd+Block.totdiag, Block.strength); 
   */

}               /* end of count_block */
/*===================================================================*/
void
count_position()
{
  int seq1, seq2, col, aa;
  double weight, pb[AAS], diffaas, naas[AAS];

  /*
     printf("\n%s nseq=%d width=%d\n", Block.ac, Block.nseq, Block.width);
   */
  for (col = 0; col < Block.width; col++) {
    /*
       compute the position-based weight & clump for each aa in col 
     */
    diffaas = 0;
    for (aa = 0; aa < AAS; aa++) {
      naas[aa] = pb[aa] = 0.0;
    }

    /*
       count how many of each AA in the column 
     */
    for (seq1 = 0; seq1 < Block.nseq; seq1++)
      if (Block.aa[seq1][col] >= 0 && Block.aa[seq1][col] < AAS) {
        naas[Block.aa[seq1][col]] += 1;
        Block.cluster[seq1] = seq1; /* initial values */
        Block.ncluster[seq1] = 1;
      }

    /*
       number of different types of AAs in col 
     */
    for (aa = 0; aa < AAS; aa++)
      if (naas[aa] > 0.0)
        diffaas++;

    /*
       now compute the position-based weight for each AA 
     */
    for (aa = 0; aa < AAS; aa++)
      if (diffaas > 0.0 && naas[aa] > 0.0) {
        pb[aa] = sqrt(1.0 / (diffaas * naas[aa]));
        /*
           printf("%.4f ", pb[aa]);
         */
        /*
           clump rows now by assigning same cluster 
         */
        /*
           now the cluster number is the sequence number, so assign a cluster number always 
           larger than that 
         */
        if (pb[aa] < PBParameter)
          for (seq1 = 0; seq1 < Block.nseq; seq1++)
            if (Block.aa[seq1][col] == aa) {
              Block.cluster[seq1] = MAXSEQS + aa;
              Block.ncluster[seq1] = naas[aa];
            }
      }
    /*
       printf("\n");
     */

    /*-------------------count between column clumps now ----------*/
    for (seq1 = 0; seq1 < Block.nseq; seq1++)
      for (seq2 = seq1 + 1; seq2 < Block.nseq; seq2++)
        if (Block.aa[seq1][col] >= 0 && Block.aa[seq1][col] < AAS &&
          Block.aa[seq2][col] >= 0 && Block.aa[seq2][col] < AAS &&
          Block.cluster[seq1] != Block.cluster[seq2]) {
          weight = (double) 1.0 / Block.ncluster[seq2];
          weight *= (double) 1.0 / Block.ncluster[seq1];
          Counts[Block.aa[seq1][col]][Block.aa[seq2][col]] += weight;
          FTotPairs += 1.0;
          FTotWeight += weight;
          if (Block.aa[seq1][col] == Block.aa[seq2][col])
            Block.totdiag += weight;
          else
            Block.totoffd += weight;
        }
  }             /* end of col */

  /*
     SPSS 
   */
  if (Block.totoffd > 0.0 || Block.totdiag > 0.0) {
    TotBlk++;
    /*
       printf("%s %9.2f %9.2f %9.2f %d\n", Block.ac, Block.totdiag, Block.totoffd,
       Block.totdiag+Block.totoffd, Block.strength); 
     */
  }
}               /* end of count_position */
/*===================================================================*/
void
count_weight()
{
  int seq1, seq2, col;
  double weight;

  for (col = 0; col < Block.width; col++)
    for (seq1 = 0; seq1 < Block.nseq; seq1++)
      for (seq2 = seq1 + 1; seq2 < Block.nseq; seq2++)
        if (Block.aa[seq1][col] >= 0 && Block.aa[seq1][col] < AAS &&
          Block.aa[seq2][col] >= 0 && Block.aa[seq2][col] < AAS) {
          weight = Block.nseq * Block.weight[seq2] / Block.wtot;
          weight *= Block.nseq * Block.weight[seq1] / Block.wtot;
          Counts[Block.aa[seq1][col]][Block.aa[seq2][col]] += weight;
          FTotPairs += 1.0;
          FTotWeight += weight;
          if (Block.aa[seq1][col] == Block.aa[seq2][col])
            Block.totdiag += weight;
          else
            Block.totoffd += weight;
        }
  /*
     SPSS 
   */
  if (Block.totoffd > 0.0 || Block.totdiag > 0.0) {
    TotBlk++;
    /*
       printf("%s %9.2f %9.2f %9.2f %d\n", Block.ac, Block.totdiag, Block.totoffd,
       Block.totdiag+Block.totoffd, Block.strength); 
     */
  }
}               /* end of count_weight */
/*===================================================================*/
void
count_cluster()
{
  int seq1, seq2, col;
  double weight;

  for (col = 0; col < Block.width; col++)
    for (seq1 = 0; seq1 < Block.nseq; seq1++)
      for (seq2 = seq1 + 1; seq2 < Block.nseq; seq2++)
        if (Block.aa[seq1][col] >= 0 && Block.aa[seq1][col] < AAS &&
          Block.aa[seq2][col] >= 0 && Block.aa[seq2][col] < AAS &&
          Block.cluster[seq1] != Block.cluster[seq2]) {
          // IMPORTANT CODING CHANGE, implemented by Henikoff in blimps
          // code.  This will make the weights proper.
          weight = (double) 1.0 / Block.ncluster[seq2];
          weight *= (double) 1.0 / Block.ncluster[seq1];


          Counts[Block.aa[seq1][col]][Block.aa[seq2][col]] += weight;
          FTotPairs += 1.0;
          FTotWeight += weight;
          if (Block.aa[seq1][col] == Block.aa[seq2][col])
            Block.totdiag += weight;
          else
            Block.totoffd += weight;
        }
  /*
     SPSS 
   */
  if (Block.totoffd > 0.0 || Block.totdiag > 0.0) {
    TotBlk++;
    TotClump += (Block.totdiag + Block.totoffd) * Block.nclus;
    TotSeg += Block.nseq;
    /*
       printf("%s %9.2f %9.2f %9.2f %d\n", Block.ac, Block.totdiag, Block.totoffd,
       Block.totdiag+Block.totoffd, Block.strength); 
     */
  }
}
char *num_to_aachar(num)
int num;
{
  switch (num) {
    case 0: return("A");
    case 1: return("R");
    case 2: return("N");
    case 3: return("D");
    case 4: return("C");
    case 5: return("Q");
    case 6: return("E");
    case 7: return("G");
    case 8: return("H");
    case 9: return("I");
    case 10: return("L");
    case 11: return("K");
    case 12: return("M");
    case 13: return("F");
    case 14: return("P");
    case 15: return("S");
    case 16: return("T");
    case 17: return("W");
    case 18: return("Y");
    case 19: return("V");
    case 20: return("X");
    case -1: return(".");
    default: return("*");   /* Should never happen */
    }
}
/*======================================================================*/
/* Amino acid to number                  */
/*  B is changed to D, Z is changed to E, O and J are changed to X  */
int aachar_to_num(ch)
char ch;
{
  switch (ch) {
    case 'A': return(0);
    case 'R': return(1);
    case 'N': return(2);
    case 'D': return(3);
    case 'B': return(3);
    case 'C': return(4);
    case 'Q': return(5);
    case 'E': return(6);
    case 'Z': return(6);
    case 'G': return(7);
    case 'H': return(8);
    case 'I': return(9);
    case 'L': return(10);
    case 'K': return(11);
    case 'M': return(12);
    case 'F': return(13);
    case 'P': return(14);
    case 'S': return(15);
    case 'T': return(16);
    case 'W': return(17);
    case 'Y': return(18);
    case 'V': return(19);
    case 'J': return(20);
    case 'O': return(20);
    case 'X': return(20);
    case '.': return(-1);
    default: return(-1);
    }
}


               /* end of count_cluster */
/*======================================================================*/
/*
   Cluster sequences in a block based on the number of 
 */
/*
   identities within the block. Sets Block.cluster & Block.ncluster 
 */
/*
   Sets Block.nclus = total number of clusters 
 */
/*
   1. Compute number of identities for each possible pair of seqs. 
 */
/*
   Results stored in lower half of matrix (pairs).  
 */
/*
   2. Use clustering threshold % of # of AAs in trimmed block.  
 */
/*
   3. Cluster recursively by traversing cols, rows of matrix.  
 */
/*
   UNIX NOTE: Program aborts when running under UNIX at free(pairs), so use the fixed size
   declaration pairs & remove the malloc() & free() calls when compiling for UNIX 
 */
/*======================================================================*/
void
cluster_seqs()
{
  int clus, npair, threshold, s1, s2, l1, l2, px, i, i1, i2;
  int nclus[MAXSEQS], minclus, oldclus;
  struct pair *pairs;
  /*
     UNIX struct pair pairs[MAXSEQS*(MAXSEQS-1)/2]; 
   */

  npair = Block.nseq * (Block.nseq - 1) / 2;
  pairs = (struct pair *) malloc(npair * sizeof(struct pair));
  if (pairs == NULL) {
    printf("\ncluster_seqs: Unable to allocate pair structure!\n");
    exit(-1);
  }
  threshold = (int) (Cluster * (Block.width)) / 100;

  /*
     Compute scores for all possible pairs of sequences 
   */
  for (s1 = 0; s1 < Block.nseq - 1; s1++) { /* col = 0, n-2 */
    l1 = 0;
    for (s2 = s1 + 1; s2 < Block.nseq; s2++) {  /* row = col+1, n-1 */
      l2 = 0;
      px = INDEX(Block.nseq, s1, s2);
      pairs[px].score = 0;
      pairs[px].cluster = -1;
      for (i = 0; i <= Block.width; i++) {
        i1 = l1 + i;
        i2 = l2 + i;
        if (i1 >= 0 && i1 < Block.width &&
          i2 >= 0 && i2 < Block.width && Block.aa[s1][i1] == Block.aa[s2][i2])
          pairs[px].score += 1;
      }
    }           /* end of s2 */
  }             /* end of s1 */

  /*
     Print scores 
   */
  /*
     printf("\nThreshold=%d", threshold); for (s2=1; s2<Block.nseq; s2++) { printf ("\n"); for
     (s1=0; s1<s2; s1++) { px = INDEX(Block.nseq, s1, s2); printf(" %.3d", pairs[px].score); } } 
   */

/*-------Cluster if score exceeds threshold by scanning cols (s1) */
  for (s1 = 0; s1 < Block.nseq; s1++) {
    Block.cluster[s1] = -1; /* clear out old values */
    Block.ncluster[s1] = 1;
    nclus[s1] = 0;
  }
  clus = 0;         /* cluster number */
  for (s1 = 0; s1 < Block.nseq - 1; s1++) /* col = 0, n-2 */
    for (s2 = s1 + 1; s2 < Block.nseq; s2++) {  /* row = col+1, n-1 */
      px = INDEX(Block.nseq, s1, s2);
      if (pairs[px].score >= threshold) { /* cluster this pair */
        if (Block.cluster[s1] < 0) {  /* s1 not yet clustered */
          if (Block.cluster[s2] < 0) {  /* new cluster */
            Block.cluster[s1] = clus++;
            Block.cluster[s2] = Block.cluster[s1];
          } else    /* use s2's cluster */
            Block.cluster[s1] = Block.cluster[s2];
        }
        /*
           use s1's cluster if it has one and s2 doesn't 
         */
        else if (Block.cluster[s1] >= 0 && Block.cluster[s2] < 0)
          Block.cluster[s2] = Block.cluster[s1];
        /*
           merge the two clusters into the lower number 
         */
        else if (Block.cluster[s1] >= 0 && Block.cluster[s2] >= 0) {
          minclus = Block.cluster[s1];
          oldclus = Block.cluster[s2];
          if (Block.cluster[s2] < Block.cluster[s1]) {
            minclus = Block.cluster[s2];
            oldclus = Block.cluster[s1];
          }
          for (i1 = 0; i1 < Block.nseq; i1++)
            if (Block.cluster[i1] == oldclus)
              Block.cluster[i1] = minclus;
        }
      }         /* end of if pairs */
    }           /* end of s2 */

   /*---  Set Block.ncluster, get rid of negative cluster numbers --*/
  // IMPORTANT CODING CHANGE: oK, not really important, but labelled as such for
  // consistency... attempt to avoid negative cluster numbers had failed.
  for (s1 = 0; s1 < Block.nseq; s1++)
    if (Block.cluster[s1] >= 0)
      nclus[Block.cluster[s1]]++;
  // if (Block.ncluster[s1] >= 0) nclus[Block.cluster[s1]]++;
  for (s1 = 0; s1 < Block.nseq; s1++) {
    if (Block.cluster[s1] < 0) {
      Block.cluster[s1] = clus++;
      Block.ncluster[s1] = 1;
    } else
      Block.ncluster[s1] = nclus[Block.cluster[s1]];

  }
  /*
     Count the total number of clusters and put in Block.nclus, the numbers in Block.ncluster[]
     are arbitrary 
   */
  Block.nclus = clus;     /* number of clumps */
  free(pairs);
}               /* end of cluster_seqs */