#include <math.h>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <vector>

using namespace std;

// Number of aminoacids 
#define AAS 20
// Min. strength to count 
#define MINSTR 0
// Max. strength to count 
#define MAXSTR 9999
// Definition of YES word 
#define YES                1
// Definition of NO word 
#define NO                 0

/*
 * INDEX & INDEXCOL compute the sequential indices for the lower half of
 * an nxn symmetric matrix given row & column coordinates.  Lower half
 * has n(n-1)/2 entries; col=0,n-2 and row=col+1,n-1; col has n-col-1
 * rows.
 * Index runs from 0 to n(n-1)/2 - 1 down columns with
 * (index=0)==(col=0,row=1) and (index=n(n-1)/2-1)==(col=n-2,row=n-1).
 */
#define INDEX(n, col, row)  ( col*n - (col*(col+3))/2 - 1 + row )
// Max number of sequences to be analyzed
#define MAXSEQS 16583
// Max. width of merged blocks
#define MAX_MERGE_WIDTH 55
// Max length of file name
#define FNAMELEN 80
// Max line length for ASCII file
#define MAXLINE 240

/*
 * Structure for pairs of sequences
 * pair should be allocated as an array, & the number of the
 * sequences forming the pair inferred from the array index
 */
struct pairScoreCluster {
  // # of identities within trimmed block
  int score;
  // cluster # for this pair
  int cluster;
};

// Block structure 
struct block { 
  char ac[11];
  int nseq, width, strength, nclus;
  // Aminoacids for sequence 
  int aa[MAXSEQS][MAX_MERGE_WIDTH];
  // seq weights found in block 
  double weight[MAXSEQS];
  // cluster # for seq
  int cluster[MAXSEQS];
  // #seqs in same cluster
  int ncluster[MAXSEQS];
  double totdiag, totoffd, wtot;
} Block;

int read_dat(FILE *fdat);
void fill_block(FILE *fdat);
void count_block();
void count_weight();
void count_cluster();
void count_position();
void cluster_seqs();
int aa_to_num(string aa);
string num_to_aa(int num);

// Values after counting the different aminoacids
vector<vector <double> > Counts(AAS);

// Total for Pairs , Sequences , Total Aminoacids and the total of width
unsigned long TotPairs, TotSeqs, TotAas, TotWidth;

double FTotPairs, FTotWeight;

// # blocks contributing to matrix
unsigned long TotBlk;

// # clumps contributing to matrix
double TotClump;

// # seqs contributing to matrix 
unsigned long TotSeg;

// Pairs for an AA
unsigned long AaPairs[AAS];

double FAaPairs[AAS];

double PBParameter;

// AA frequency
unsigned long AaFreq[AAS];

int MinStr, MaxStr, Cluster;

void initialize(){
  for(int i=0;i<AAS;i++){
    Counts[i].resize(AAS);
  }
}

int main( int argc, char *argv[]) {
  // Input file and output file
  FILE *fdat, *fout;

  // Name of the database file (Input File)
  string datfile;
  // Name for the ouput file
  char outfile[80];
  // Temporal variable to save values received from Console
  string ctemp;

  // Total of Aminoacids
  unsigned long totaas;

  vector<vector<int> >tij(AAS);

  for(int i=0;i<AAS;i++){
    tij[i].resize(AAS);
  }
  
  int totblk, row, col, maxsij, minsij, sumsij;

  int iscale, is, itemp;

  vector<vector<double> >sij(AAS+3);

  for(int i=0;i<AAS;i++){
    sij[i].resize(AAS);
  }
  

  double totpairs, totdiag, totoffd;
  double dtemp;
  double ftotpairs, ftotdiag, ftotoffd;
  double s, fij, fifj, oij, entropy, expected;
  double sumpi, sumqij;

  // Blocks Database Name
  cin >> datfile;

  if ((fdat = fopen(datfile.c_str(), "r")) == NULL) {
    printf("Cannot open %s\n", datfile.c_str());
    exit(-1);
  } else{
    printf("Reading %s\n", datfile.c_str());
  }

  // Minimum block strength

  MinStr = MINSTR;
  printf("Enter minimum block strength [%d]: ", MinStr);
  cin >> ctemp;;
  if (ctemp.length())
    MinStr = atoi(ctemp.c_str());

  // Maximum block strength
  MaxStr = MAXSTR;
  printf("Enter maximum block strength [%d]: ", MaxStr);
  cin >> ctemp;
  if (ctemp.length())
    MaxStr = atoi(ctemp.c_str());

  printf("Minimum block strength=%d, Maximum block strength=%d\n", MinStr, MaxStr);

  // Clustering
  Cluster = -1;

  printf("Enter n for no clustering\n");
  printf("   or a number between 0 and 100 for percent identity clustering\n");
  printf("   or e for existing (implied) clusters\n");
  printf("   or w for existing sequence weights\n");
  printf("   or pn for position-based weights, PB weight = 1/n\n");
  printf("Enter clustering identity percentage or n/e/w/pn [e]: ");
  cin >> ctemp;

  // e => -1, n=> -2, w=> -3 p=> -4

  Cluster = -1;

  if (ctemp.length()){
    if (ctemp[0] == 'n' || ctemp[0] == 'N')
      Cluster = -2;
    else if (ctemp[0] == 'w' || ctemp[0] == 'W')
      Cluster = -3;
    else if (ctemp[0] == 'p' || ctemp[0] == 'P') {
      Cluster = -4;
      itemp = atoi(ctemp.substr(1,1).c_str());
      if (itemp > 0)
        PBParameter = sqrt((double) 1.0 / itemp);
      else
        PBParameter = 1.0;
    } else if (ctemp[0] != 'e' && ctemp[0] != 'E')
      Cluster = atoi(ctemp.c_str());
  }
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

  // Scale
  iscale = 0;
  printf("Enter scale n for 1/n bits [%d]: ", iscale);
    cin >> ctemp;
    if (ctemp.length())
      iscale = atoi(ctemp.c_str());

  if (iscale < 0)
    iscale = 0;
  if (iscale > 100)
    iscale = 100;
  if (iscale == 0)
    printf("Scale based on relative entropy\n");
  else
    printf("Requested scale = 1/%d bits\n", iscale);
  
  // Initialize
  initialize();
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

  // Read Blocks
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

  // The frequency matrix fij
  printf("Frequencies = fij pairs (off-diagonals = 2*fij):\n");

  for (col = 0; col < AAS; col++)
    printf("    %1s  ", num_to_aa(col).c_str());
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
    printf("    %1s  ", num_to_aa(col).c_str());
  printf("\n");

  // Print the targe frequencies qij
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
  fprintf(fout, "#  Blocks Database = %s\n", datfile.c_str());

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
    fprintf(fout, "   %1s   ", num_to_aa(col).c_str());
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

  printf("sumqij = %f\n", sumqij);
  // Marginal frequencies
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
    printf("    %1s  ", num_to_aa(col).c_str());
  printf("\n");

  for (col = 0; col < AAS; col++) {
    printf("%.3f ", FAaPairs[col] / totpairs); /* p(col,*) */
    sumpi += FAaPairs[col] / totpairs;
  }

  printf("\nsumpi=%.3f", sumpi);

  totaas = (unsigned long) 0;

  printf("\nAA Frequencies = ai:\n");

  for (col = 0; col < AAS; col++) {
    printf("%8ld ", AaFreq[col]);
    totaas += AaFreq[col];
  }
  // Aminoacid frequencies
  printf("\nAA Probabilities = fi:\n");

  sumpi = 0.0;

  for (col = 0; col < AAS; col++)
    printf("    %1s  ", num_to_aa(col).c_str());
  printf("\n");

  for (col = 0; col < AAS; col++) {
    printf("%.3f ", (double) AaFreq[col] / totaas);	/* fi */
    sumpi += (double) AaFreq[col] / totaas;
  }
  printf("\nsumpi=%.3f", sumpi);

  printf("\n\ntotpairs=%.3f, FTotWeight=%.3f", totpairs, FTotWeight);
  printf("\n totdiag=%.3f, totoffd=%.3f", totdiag, totoffd);
  printf(" totaas=%ld, TotAas=%ld\n", totaas, TotAas);

  // The bit matrix
  entropy = expected = 0.0;
  ftotpairs = totpairs;
  for (row = 0; row < AAS; row++) {
    for (col = 0; col <= row; col++) {
      if (row == col)
        fij = (double) Counts[row][col];
      else
        fij = (double) (Counts[row][col] + Counts[col][row]) / 2.0;
      // FAaPairs[i] = f(i,*)
      fifj = (double) FAaPairs[row] * FAaPairs[col];
      if (fifj > 0.000000001)
        oij = ftotpairs * fij / fifj;
      else
        oij = 0.0;
      // Log odds ratio 
      if (oij > 0.000000001)
        s = log((double) oij);
      else
        s = -20.0;    // minus infinity
      // Round off for the log base 10 matrix
      dtemp = 10.0 * s / log(10.0);
      // log base 10
      tij[row][col] = round(dtemp);
      // Log base 2 (bit) matrix
      s /= log(2.0);    // log base 2 = bits 
      // compute entropy & expected value in bits
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


  // Fill in symmetric part of matrix
  for (row = 0; row < AAS; row++)
    for (col = 0; col <= row; col++)
      if (row != col)
        sij[col][row] = sij[row][col];

  // Print the sij matrix in bit units
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
  fprintf(fout, "#  Blocks Database = %s\n", datfile.c_str());

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
    fprintf(fout, "   %1s     ", num_to_aa(col).c_str());
  fprintf(fout, "\n");

  for (row = 0; row < AAS; row++) {
    for (col = 0; col <= row; col++)
      fprintf(fout, "% 8.4f ", sij[row][col]);
    fprintf(fout, "\n");
  }

  fclose(fout);

  // Determine the scale based on entropy
  if (iscale == 0) {
    dtemp = 2.0 / sqrt(entropy);
  
    iscale = round(dtemp);
    if (iscale < 2)
      iscale = 2;
  }

  // Print the integer matrix
  if (Cluster >= 0)
    sprintf(outfile, "blosum%d.iij", Cluster);
  else if (Cluster == -4)
    sprintf(outfile, "blosump%d.iij", (int) (1. / PBParameter));
  else if (Cluster == -3)
    sprintf(outfile, "blosumw.iij");
  else if (Cluster == -2)
    sprintf(outfile, "blosumn.iij");
  else if (Cluster == -1)
    sprintf(outfile, "blosume.iij");
  if ((fout = fopen(outfile, "wt")) == NULL)
    fout = stdout;
  
  printf("\nInteger scoring matrix in 1/%d bit units in %s\n", iscale, outfile);

  fprintf(fout, "#  BLOSUM Clustered Scoring Matrix in 1/%d Bit Units\n", iscale);
  fprintf(fout, "#  Blocks Database = %s\n", datfile.c_str());

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
    fprintf(fout, " %1s  ", num_to_aa(col).c_str());
  fprintf(fout, "\n");

  for (row = 0; row < AAS; row++) {
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

  // Add the high pass value = minsij
  printf("\nPositive matrix (+%d):\n", 0 - minsij);
  for (col = 0; col < AAS; col++)
    printf("    %1s  ", num_to_aa(col).c_str());
  printf("\n");
  
  for (row = 0; row < AAS; row++) {
    for (col = 0; col <= row; col++) {
      s = (double) sij[row][col] * iscale;
      is = round(s) - minsij;
      printf("%3d ", is);
    }
    printf("\n");
  }

  for (col = 0; col < AAS; col++)
    printf("    %1s  ", num_to_aa(col).c_str());
  printf("\n");


  // The log base 10 matrix
  sprintf(outfile, "blosumLog10.iij");
  if ((fout = fopen(outfile, "wt")) == NULL)
    fout = stdout;
  fprintf(fout,"\n10 times log base 10 matrix:\n");
  for (col = 0; col < AAS; col++)
    fprintf(fout," %1s  ", num_to_aa(col).c_str());
  fprintf(fout,"\n");

  for (row = 0; row < AAS; row++) {
    for (col = 0; col <= row; col++)
      fprintf(fout,"%3d ", tij[row][col]);
    fprintf(fout,"\n");
  }

  for (col = 0; col < AAS; col++)
    fprintf(fout," %1s  ", num_to_aa(col).c_str());
  fprintf(fout,"\n");
  fclose(fout);

  exit(0);
}/* end of main */

// Read Blocks Database
int read_dat(FILE *fdat) {
  int totblk;

  // Temporal variables used for reading the file
  char line[MAXLINE], *ptr, *ptr1;

  // Initializate the variable totblk to 0
  totblk = 0;

  while (!feof(fdat) && fgets(line, MAXLINE, fdat) != NULL) {
    if (strncmp(line, "AC   ", 5) == 0) {
      strncpy(Block.ac, line + 5, 9);
      Block.ac[10] = '\0';
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
} /* end of read_dat */

// Fill up the block structure
void fill_block(FILE *fdat){
  int done, i, li, n, cluster, ncluster;
  char line[MAXLINE], *ptr;

  Block.nseq = Block.width = 0;
  Block.totdiag = Block.totoffd = Block.wtot = (double) 0;
  cluster = ncluster = 0;
  done = NO;

  while (!done && !feof(fdat) && fgets(line, MAXLINE, fdat) != NULL) {
    // blank line => new cluster
    if (strlen(line) == 1) {  
      // Set #seqs in cluster to seqs in previous cluster
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
        // assuming no spaces within block !
        ptr = strtok(line + li, " \r\n\t");
        
        for (i = 0; i < (int)strlen(ptr); i++) {
          string test(ptr);
          test = test.substr(i,1);

          Block.aa[Block.nseq][i] = aa_to_num(test);
          if (Block.aa[Block.nseq][i] >= 0 && Block.aa[Block.nseq][i] < AAS) {
            AaFreq[Block.aa[Block.nseq][i]]++;
            TotAas++;
          }
        }
        li += strlen(ptr);
        ptr = strtok(NULL, " \r\n\t");
        if (ptr != NULL) {
          Block.weight[Block.nseq] = atof(ptr);
          Block.wtot += Block.weight[Block.nseq];
        } else
          Block.weight[Block.nseq] = 0;
        Block.cluster[Block.nseq] = cluster;
        // # seqs in current cluster
        ncluster++;  
        Block.width = i;
        Block.nseq++;
      }
    }
  }

  // Compute weights for the last cluster
  if (ncluster > 0)
    for (n = 0; n < Block.nseq; n++)
      if (Block.cluster[n] == cluster)
        Block.ncluster[n] = ncluster;
  TotSeqs += Block.nseq;
  TotWidth += Block.width;
} /* end of fill_block */

void count_block() {
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

  TotBlk++;
}               /* end of count_block */

void count_position() {
  int seq1, seq2, col, aa;
  double weight, pb[AAS], diffaas, naas[AAS];

  for (col = 0; col < Block.width; col++) {
    // compute the position-based weight & clump for each aa in col
    diffaas = 0;
    for (aa = 0; aa < AAS; aa++) {
      naas[aa] = pb[aa] = 0.0;
    }

    // count how many of each AA in the column
    for (seq1 = 0; seq1 < Block.nseq; seq1++)
      if (Block.aa[seq1][col] >= 0 && Block.aa[seq1][col] < AAS) {
        naas[Block.aa[seq1][col]] += 1;
        Block.cluster[seq1] = seq1; /* initial values */
        Block.ncluster[seq1] = 1;
      }

    // number of different types of AAs in col
    
    for (aa = 0; aa < AAS; aa++)
      if (naas[aa] > 0.0)
        diffaas++;

    // now compute the position-based weight for each AA
    for (aa = 0; aa < AAS; aa++)
      if (diffaas > 0.0 && naas[aa] > 0.0) {
        pb[aa] = sqrt(1.0 / (diffaas * naas[aa]));

        // clump rows now by assigning same cluster
        /*
         * now the cluster number is the sequence number, so assign a cluster number always
         * larger than that
        */
        if (pb[aa] < PBParameter)
          for (seq1 = 0; seq1 < Block.nseq; seq1++)
            if (Block.aa[seq1][col] == aa) {
              Block.cluster[seq1] = MAXSEQS + aa;
              Block.ncluster[seq1] = naas[aa];
            }
      }

    // count between column clumps now 
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
  }


  // SPSS
  if (Block.totoffd > 0.0 || Block.totdiag > 0.0) {
    TotBlk++;
  }
} /* end of count_position */

void count_weight() {
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
  // SPSS
  if (Block.totoffd > 0.0 || Block.totdiag > 0.0) {
    TotBlk++;
  }
} /* end of count_weight */

void count_cluster() {
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

  // SPSS
  if (Block.totoffd > 0.0 || Block.totdiag > 0.0) {
    TotBlk++;
    TotClump += (Block.totdiag + Block.totoffd) * Block.nclus;
    TotSeg += Block.nseq;
  }
} /* end of count_cluster */

/*
 * Cluster sequences in a block based on the number of
 * identities within the block. Sets Block.cluster & Block.ncluster
 * Sets Block.nclus = total number of clusters
 * 
 * 1. Compute number of identities for each possible pair of seqs.
 * Results stored in lower half of matrix (pairs).
 * 
 * 2. Use clustering threshold % of # of AAs in trimmed block.
 * 
 * 3. Cluster recursively by traversing cols, rows of matrix.
 * 
 * UNIX NOTE: Program aborts when running under UNIX at free(pairs), so use the fixed size
 * declaration pairs & remove the malloc() & free() calls when compiling for UNIX
*/

void cluster_seqs() {
  int  s1, s2;
  int clus, npair, threshold, l1, l2, px, i, i1, i2;
  int nclus[MAXSEQS], minclus, oldclus;
  struct pairScoreCluster *pairs;

  // UNIX struct pair pairs[MAXSEQS*(MAXSEQS-1)/2];

  npair = Block.nseq * (Block.nseq - 1) / 2;
  pairs = (struct pairScoreCluster *) malloc(npair * sizeof(struct pairScoreCluster));
  if (pairs == NULL) {
    printf("\ncluster_seqs: Unable to allocate pair structure!\n");
    exit(-1);
  }

  threshold = (int) (Cluster * (Block.width)) / 100;

  // Compute scores for all possible pairs of sequences

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
    } /* end of s2 */
  } /* end of s1 */

  // Cluster if score exceeds threshold by scanning cols (s1)
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

        // use s1's cluster if it has one and s2 doesn't
        
        else if (Block.cluster[s1] >= 0 && Block.cluster[s2] < 0)
          Block.cluster[s2] = Block.cluster[s1];

        // merge the two clusters into the lower number

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
      } /* end of if pairs */
    } /* end of s2 */

  // Set Block.ncluster, get rid of negative cluster numbers 
  // IMPORTANT CODING CHANGE: oK, not really important, but labelled as such for
  // consistency... attempt to avoid negative cluster numbers had failed.
  for (s1 = 0; s1 < Block.nseq; s1++)
    if (Block.cluster[s1] >= 0)
      nclus[Block.cluster[s1]]++;
  for (s1 = 0; s1 < Block.nseq; s1++) {
    if (Block.cluster[s1] < 0) {
      Block.cluster[s1] = clus++;
      Block.ncluster[s1] = 1;
    } else
      Block.ncluster[s1] = nclus[Block.cluster[s1]];

  }
  
  /*
   * Count the total number of clusters and put in Block.nclus, the numbers in Block.ncluster[]
   * are arbitrary
  */
  // number of clumps
  Block.nclus = clus;     
  free(pairs);
} /* end of cluster_seqs */
int aa_to_num(string aa){ 
  if ( aa == "A") 
    return 0;
  else if( aa == "R")  
    return 1;
  else if( aa == "N")  
    return 2;
  else if( aa == "D")  
    return 3;
  else if( aa == "C")  
    return 4;
  else if( aa == "Q")  
    return 5;
  else if( aa == "E")  
    return 6;
  else if( aa == "G")  
    return 7;
  else if( aa == "H")  
    return 8;
  else if( aa == "I")  
    return 9;
  else if( aa == "L")  
    return 10;
  else if( aa == "K")  
    return 11;
  else if( aa == "M")  
    return 12;
  else if( aa == "F")  
    return 13;
  else if( aa == "P")  
    return 14;
  else if( aa == "S")  
    return 15;
  else if( aa == "T")  
    return 16;
  else if( aa == "W")  
    return 17;
  else if( aa == "Y")  
    return 18;
  else if( aa == "V")  
    return 19;
  return 0;
} 

string num_to_aa(int num) { 
  if ( num == 0) 
    return "A";
  else if( num == 1) 
    return "R";
  else if( num == 2) 
    return "N";
  else if( num == 3) 
    return "D";
  else if( num == 4) 
    return "C";
  else if( num == 5) 
    return "Q";
  else if( num == 6) 
    return "E";
  else if( num == 7) 
    return "G";
  else if( num == 8) 
    return "H";
  else if( num == 9) 
    return "I";
  else if( num == 10) 
    return "L";
  else if( num == 11) 
    return "K";
  else if( num == 12) 
    return "M";
  else if( num == 13) 
    return "F";
  else if( num == 14) 
    return "P";
  else if( num == 15) 
    return "S";
  else if( num == 16) 
    return "T";
  else if( num == 17) 
    return "W";
  else if( num == 18) 
    return "Y";
  else if( num == 19) 
    return "V";
  return "A";
} 

