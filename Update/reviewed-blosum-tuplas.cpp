#include <math.h>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <vector>

using namespace std;

// Number of aminoacids 
#define AAS 400
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
        int cont = 0;
        for (i = 0; i < (int)strlen(ptr); i+=2) {
          string test(ptr);
          test = test.substr(i,2);
          if(test.length() < 2)
            continue;
          Block.aa[Block.nseq][cont] = aa_to_num(test);
          if (Block.aa[Block.nseq][cont] >= 0 && Block.aa[Block.nseq][cont] < AAS) {
            AaFreq[Block.aa[Block.nseq][cont]]++;
            TotAas++;
          }
          cont++;
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
        Block.width = cont;
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
  if ( aa == "AA") 
    return 0;
  else if( aa == "AR")  
    return 1;
  else if( aa == "AN")  
    return 2;
  else if( aa == "AD")  
    return 3;
  else if( aa == "AC")  
    return 4;
  else if( aa == "AQ")  
    return 5;
  else if( aa == "AE")  
    return 6;
  else if( aa == "AG")  
    return 7;
  else if( aa == "AH")  
    return 8;
  else if( aa == "AI")  
    return 9;
  else if( aa == "AL")  
    return 10;
  else if( aa == "AK")  
    return 11;
  else if( aa == "AM")  
    return 12;
  else if( aa == "AF")  
    return 13;
  else if( aa == "AP")  
    return 14;
  else if( aa == "AS")  
    return 15;
  else if( aa == "AT")  
    return 16;
  else if( aa == "AW")  
    return 17;
  else if( aa == "AY")  
    return 18;
  else if( aa == "AV")  
    return 19;
  else if( aa == "RA")  
    return 20;
  else if( aa == "RR")  
    return 21;
  else if( aa == "RN")  
    return 22;
  else if( aa == "RD")  
    return 23;
  else if( aa == "RC")  
    return 24;
  else if( aa == "RQ")  
    return 25;
  else if( aa == "RE")  
    return 26;
  else if( aa == "RG")  
    return 27;
  else if( aa == "RH")  
    return 28;
  else if( aa == "RI")  
    return 29;
  else if( aa == "RL")  
    return 30;
  else if( aa == "RK")  
    return 31;
  else if( aa == "RM")  
    return 32;
  else if( aa == "RF")  
    return 33;
  else if( aa == "RP")  
    return 34;
  else if( aa == "RS")  
    return 35;
  else if( aa == "RT")  
    return 36;
  else if( aa == "RW")  
    return 37;
  else if( aa == "RY")  
    return 38;
  else if( aa == "RV")  
    return 39;
  else if( aa == "NA")  
    return 40;
  else if( aa == "NR")  
    return 41;
  else if( aa == "NN")  
    return 42;
  else if( aa == "ND")  
    return 43;
  else if( aa == "NC")  
    return 44;
  else if( aa == "NQ")  
    return 45;
  else if( aa == "NE")  
    return 46;
  else if( aa == "NG")  
    return 47;
  else if( aa == "NH")  
    return 48;
  else if( aa == "NI")  
    return 49;
  else if( aa == "NL")  
    return 50;
  else if( aa == "NK")  
    return 51;
  else if( aa == "NM")  
    return 52;
  else if( aa == "NF")  
    return 53;
  else if( aa == "NP")  
    return 54;
  else if( aa == "NS")  
    return 55;
  else if( aa == "NT")  
    return 56;
  else if( aa == "NW")  
    return 57;
  else if( aa == "NY")  
    return 58;
  else if( aa == "NV")  
    return 59;
  else if( aa == "DA")  
    return 60;
  else if( aa == "DR")  
    return 61;
  else if( aa == "DN")  
    return 62;
  else if( aa == "DD")  
    return 63;
  else if( aa == "DC")  
    return 64;
  else if( aa == "DQ")  
    return 65;
  else if( aa == "DE")  
    return 66;
  else if( aa == "DG")  
    return 67;
  else if( aa == "DH")  
    return 68;
  else if( aa == "DI")  
    return 69;
  else if( aa == "DL")  
    return 70;
  else if( aa == "DK")  
    return 71;
  else if( aa == "DM")  
    return 72;
  else if( aa == "DF")  
    return 73;
  else if( aa == "DP")  
    return 74;
  else if( aa == "DS")  
    return 75;
  else if( aa == "DT")  
    return 76;
  else if( aa == "DW")  
    return 77;
  else if( aa == "DY")  
    return 78;
  else if( aa == "DV")  
    return 79;
  else if( aa == "CA")  
    return 80;
  else if( aa == "CR")  
    return 81;
  else if( aa == "CN")  
    return 82;
  else if( aa == "CD")  
    return 83;
  else if( aa == "CC")  
    return 84;
  else if( aa == "CQ")  
    return 85;
  else if( aa == "CE")  
    return 86;
  else if( aa == "CG")  
    return 87;
  else if( aa == "CH")  
    return 88;
  else if( aa == "CI")  
    return 89;
  else if( aa == "CL")  
    return 90;
  else if( aa == "CK")  
    return 91;
  else if( aa == "CM")  
    return 92;
  else if( aa == "CF")  
    return 93;
  else if( aa == "CP")  
    return 94;
  else if( aa == "CS")  
    return 95;
  else if( aa == "CT")  
    return 96;
  else if( aa == "CW")  
    return 97;
  else if( aa == "CY")  
    return 98;
  else if( aa == "CV")  
    return 99;
  else if( aa == "QA")  
    return 100;
  else if( aa == "QR")  
    return 101;
  else if( aa == "QN")  
    return 102;
  else if( aa == "QD")  
    return 103;
  else if( aa == "QC")  
    return 104;
  else if( aa == "QQ")  
    return 105;
  else if( aa == "QE")  
    return 106;
  else if( aa == "QG")  
    return 107;
  else if( aa == "QH")  
    return 108;
  else if( aa == "QI")  
    return 109;
  else if( aa == "QL")  
    return 110;
  else if( aa == "QK")  
    return 111;
  else if( aa == "QM")  
    return 112;
  else if( aa == "QF")  
    return 113;
  else if( aa == "QP")  
    return 114;
  else if( aa == "QS")  
    return 115;
  else if( aa == "QT")  
    return 116;
  else if( aa == "QW")  
    return 117;
  else if( aa == "QY")  
    return 118;
  else if( aa == "QV")  
    return 119;
  else if( aa == "EA")  
    return 120;
  else if( aa == "ER")  
    return 121;
  else if( aa == "EN")  
    return 122;
  else if( aa == "ED")  
    return 123;
  else if( aa == "EC")  
    return 124;
  else if( aa == "EQ")  
    return 125;
  else if( aa == "EE")  
    return 126;
  else if( aa == "EG")  
    return 127;
  else if( aa == "EH")  
    return 128;
  else if( aa == "EI")  
    return 129;
  else if( aa == "EL")  
    return 130;
  else if( aa == "EK")  
    return 131;
  else if( aa == "EM")  
    return 132;
  else if( aa == "EF")  
    return 133;
  else if( aa == "EP")  
    return 134;
  else if( aa == "ES")  
    return 135;
  else if( aa == "ET")  
    return 136;
  else if( aa == "EW")  
    return 137;
  else if( aa == "EY")  
    return 138;
  else if( aa == "EV")  
    return 139;
  else if( aa == "GA")  
    return 140;
  else if( aa == "GR")  
    return 141;
  else if( aa == "GN")  
    return 142;
  else if( aa == "GD")  
    return 143;
  else if( aa == "GC")  
    return 144;
  else if( aa == "GQ")  
    return 145;
  else if( aa == "GE")  
    return 146;
  else if( aa == "GG")  
    return 147;
  else if( aa == "GH")  
    return 148;
  else if( aa == "GI")  
    return 149;
  else if( aa == "GL")  
    return 150;
  else if( aa == "GK")  
    return 151;
  else if( aa == "GM")  
    return 152;
  else if( aa == "GF")  
    return 153;
  else if( aa == "GP")  
    return 154;
  else if( aa == "GS")  
    return 155;
  else if( aa == "GT")  
    return 156;
  else if( aa == "GW")  
    return 157;
  else if( aa == "GY")  
    return 158;
  else if( aa == "GV")  
    return 159;
  else if( aa == "HA")  
    return 160;
  else if( aa == "HR")  
    return 161;
  else if( aa == "HN")  
    return 162;
  else if( aa == "HD")  
    return 163;
  else if( aa == "HC")  
    return 164;
  else if( aa == "HQ")  
    return 165;
  else if( aa == "HE")  
    return 166;
  else if( aa == "HG")  
    return 167;
  else if( aa == "HH")  
    return 168;
  else if( aa == "HI")  
    return 169;
  else if( aa == "HL")  
    return 170;
  else if( aa == "HK")  
    return 171;
  else if( aa == "HM")  
    return 172;
  else if( aa == "HF")  
    return 173;
  else if( aa == "HP")  
    return 174;
  else if( aa == "HS")  
    return 175;
  else if( aa == "HT")  
    return 176;
  else if( aa == "HW")  
    return 177;
  else if( aa == "HY")  
    return 178;
  else if( aa == "HV")  
    return 179;
  else if( aa == "IA")  
    return 180;
  else if( aa == "IR")  
    return 181;
  else if( aa == "IN")  
    return 182;
  else if( aa == "ID")  
    return 183;
  else if( aa == "IC")  
    return 184;
  else if( aa == "IQ")  
    return 185;
  else if( aa == "IE")  
    return 186;
  else if( aa == "IG")  
    return 187;
  else if( aa == "IH")  
    return 188;
  else if( aa == "II")  
    return 189;
  else if( aa == "IL")  
    return 190;
  else if( aa == "IK")  
    return 191;
  else if( aa == "IM")  
    return 192;
  else if( aa == "IF")  
    return 193;
  else if( aa == "IP")  
    return 194;
  else if( aa == "IS")  
    return 195;
  else if( aa == "IT")  
    return 196;
  else if( aa == "IW")  
    return 197;
  else if( aa == "IY")  
    return 198;
  else if( aa == "IV")  
    return 199;
  else if( aa == "LA")  
    return 200;
  else if( aa == "LR")  
    return 201;
  else if( aa == "LN")  
    return 202;
  else if( aa == "LD")  
    return 203;
  else if( aa == "LC")  
    return 204;
  else if( aa == "LQ")  
    return 205;
  else if( aa == "LE")  
    return 206;
  else if( aa == "LG")  
    return 207;
  else if( aa == "LH")  
    return 208;
  else if( aa == "LI")  
    return 209;
  else if( aa == "LL")  
    return 210;
  else if( aa == "LK")  
    return 211;
  else if( aa == "LM")  
    return 212;
  else if( aa == "LF")  
    return 213;
  else if( aa == "LP")  
    return 214;
  else if( aa == "LS")  
    return 215;
  else if( aa == "LT")  
    return 216;
  else if( aa == "LW")  
    return 217;
  else if( aa == "LY")  
    return 218;
  else if( aa == "LV")  
    return 219;
  else if( aa == "KA")  
    return 220;
  else if( aa == "KR")  
    return 221;
  else if( aa == "KN")  
    return 222;
  else if( aa == "KD")  
    return 223;
  else if( aa == "KC")  
    return 224;
  else if( aa == "KQ")  
    return 225;
  else if( aa == "KE")  
    return 226;
  else if( aa == "KG")  
    return 227;
  else if( aa == "KH")  
    return 228;
  else if( aa == "KI")  
    return 229;
  else if( aa == "KL")  
    return 230;
  else if( aa == "KK")  
    return 231;
  else if( aa == "KM")  
    return 232;
  else if( aa == "KF")  
    return 233;
  else if( aa == "KP")  
    return 234;
  else if( aa == "KS")  
    return 235;
  else if( aa == "KT")  
    return 236;
  else if( aa == "KW")  
    return 237;
  else if( aa == "KY")  
    return 238;
  else if( aa == "KV")  
    return 239;
  else if( aa == "MA")  
    return 240;
  else if( aa == "MR")  
    return 241;
  else if( aa == "MN")  
    return 242;
  else if( aa == "MD")  
    return 243;
  else if( aa == "MC")  
    return 244;
  else if( aa == "MQ")  
    return 245;
  else if( aa == "ME")  
    return 246;
  else if( aa == "MG")  
    return 247;
  else if( aa == "MH")  
    return 248;
  else if( aa == "MI")  
    return 249;
  else if( aa == "ML")  
    return 250;
  else if( aa == "MK")  
    return 251;
  else if( aa == "MM")  
    return 252;
  else if( aa == "MF")  
    return 253;
  else if( aa == "MP")  
    return 254;
  else if( aa == "MS")  
    return 255;
  else if( aa == "MT")  
    return 256;
  else if( aa == "MW")  
    return 257;
  else if( aa == "MY")  
    return 258;
  else if( aa == "MV")  
    return 259;
  else if( aa == "FA")  
    return 260;
  else if( aa == "FR")  
    return 261;
  else if( aa == "FN")  
    return 262;
  else if( aa == "FD")  
    return 263;
  else if( aa == "FC")  
    return 264;
  else if( aa == "FQ")  
    return 265;
  else if( aa == "FE")  
    return 266;
  else if( aa == "FG")  
    return 267;
  else if( aa == "FH")  
    return 268;
  else if( aa == "FI")  
    return 269;
  else if( aa == "FL")  
    return 270;
  else if( aa == "FK")  
    return 271;
  else if( aa == "FM")  
    return 272;
  else if( aa == "FF")  
    return 273;
  else if( aa == "FP")  
    return 274;
  else if( aa == "FS")  
    return 275;
  else if( aa == "FT")  
    return 276;
  else if( aa == "FW")  
    return 277;
  else if( aa == "FY")  
    return 278;
  else if( aa == "FV")  
    return 279;
  else if( aa == "PA")  
    return 280;
  else if( aa == "PR")  
    return 281;
  else if( aa == "PN")  
    return 282;
  else if( aa == "PD")  
    return 283;
  else if( aa == "PC")  
    return 284;
  else if( aa == "PQ")  
    return 285;
  else if( aa == "PE")  
    return 286;
  else if( aa == "PG")  
    return 287;
  else if( aa == "PH")  
    return 288;
  else if( aa == "PI")  
    return 289;
  else if( aa == "PL")  
    return 290;
  else if( aa == "PK")  
    return 291;
  else if( aa == "PM")  
    return 292;
  else if( aa == "PF")  
    return 293;
  else if( aa == "PP")  
    return 294;
  else if( aa == "PS")  
    return 295;
  else if( aa == "PT")  
    return 296;
  else if( aa == "PW")  
    return 297;
  else if( aa == "PY")  
    return 298;
  else if( aa == "PV")  
    return 299;
  else if( aa == "SA")  
    return 300;
  else if( aa == "SR")  
    return 301;
  else if( aa == "SN")  
    return 302;
  else if( aa == "SD")  
    return 303;
  else if( aa == "SC")  
    return 304;
  else if( aa == "SQ")  
    return 305;
  else if( aa == "SE")  
    return 306;
  else if( aa == "SG")  
    return 307;
  else if( aa == "SH")  
    return 308;
  else if( aa == "SI")  
    return 309;
  else if( aa == "SL")  
    return 310;
  else if( aa == "SK")  
    return 311;
  else if( aa == "SM")  
    return 312;
  else if( aa == "SF")  
    return 313;
  else if( aa == "SP")  
    return 314;
  else if( aa == "SS")  
    return 315;
  else if( aa == "ST")  
    return 316;
  else if( aa == "SW")  
    return 317;
  else if( aa == "SY")  
    return 318;
  else if( aa == "SV")  
    return 319;
  else if( aa == "TA")  
    return 320;
  else if( aa == "TR")  
    return 321;
  else if( aa == "TN")  
    return 322;
  else if( aa == "TD")  
    return 323;
  else if( aa == "TC")  
    return 324;
  else if( aa == "TQ")  
    return 325;
  else if( aa == "TE")  
    return 326;
  else if( aa == "TG")  
    return 327;
  else if( aa == "TH")  
    return 328;
  else if( aa == "TI")  
    return 329;
  else if( aa == "TL")  
    return 330;
  else if( aa == "TK")  
    return 331;
  else if( aa == "TM")  
    return 332;
  else if( aa == "TF")  
    return 333;
  else if( aa == "TP")  
    return 334;
  else if( aa == "TS")  
    return 335;
  else if( aa == "TT")  
    return 336;
  else if( aa == "TW")  
    return 337;
  else if( aa == "TY")  
    return 338;
  else if( aa == "TV")  
    return 339;
  else if( aa == "WA")  
    return 340;
  else if( aa == "WR")  
    return 341;
  else if( aa == "WN")  
    return 342;
  else if( aa == "WD")  
    return 343;
  else if( aa == "WC")  
    return 344;
  else if( aa == "WQ")  
    return 345;
  else if( aa == "WE")  
    return 346;
  else if( aa == "WG")  
    return 347;
  else if( aa == "WH")  
    return 348;
  else if( aa == "WI")  
    return 349;
  else if( aa == "WL")  
    return 350;
  else if( aa == "WK")  
    return 351;
  else if( aa == "WM")  
    return 352;
  else if( aa == "WF")  
    return 353;
  else if( aa == "WP")  
    return 354;
  else if( aa == "WS")  
    return 355;
  else if( aa == "WT")  
    return 356;
  else if( aa == "WW")  
    return 357;
  else if( aa == "WY")  
    return 358;
  else if( aa == "WV")  
    return 359;
  else if( aa == "YA")  
    return 360;
  else if( aa == "YR")  
    return 361;
  else if( aa == "YN")  
    return 362;
  else if( aa == "YD")  
    return 363;
  else if( aa == "YC")  
    return 364;
  else if( aa == "YQ")  
    return 365;
  else if( aa == "YE")  
    return 366;
  else if( aa == "YG")  
    return 367;
  else if( aa == "YH")  
    return 368;
  else if( aa == "YI")  
    return 369;
  else if( aa == "YL")  
    return 370;
  else if( aa == "YK")  
    return 371;
  else if( aa == "YM")  
    return 372;
  else if( aa == "YF")  
    return 373;
  else if( aa == "YP")  
    return 374;
  else if( aa == "YS")  
    return 375;
  else if( aa == "YT")  
    return 376;
  else if( aa == "YW")  
    return 377;
  else if( aa == "YY")  
    return 378;
  else if( aa == "YV")  
    return 379;
  else if( aa == "VA")  
    return 380;
  else if( aa == "VR")  
    return 381;
  else if( aa == "VN")  
    return 382;
  else if( aa == "VD")  
    return 383;
  else if( aa == "VC")  
    return 384;
  else if( aa == "VQ")  
    return 385;
  else if( aa == "VE")  
    return 386;
  else if( aa == "VG")  
    return 387;
  else if( aa == "VH")  
    return 388;
  else if( aa == "VI")  
    return 389;
  else if( aa == "VL")  
    return 390;
  else if( aa == "VK")  
    return 391;
  else if( aa == "VM")  
    return 392;
  else if( aa == "VF")  
    return 393;
  else if( aa == "VP")  
    return 394;
  else if( aa == "VS")  
    return 395;
  else if( aa == "VT")  
    return 396;
  else if( aa == "VW")  
    return 397;
  else if( aa == "VY")  
    return 398;
  else if( aa == "VV")  
    return 399;
  return 0;
} 

string num_to_aa(int num) { 
  if ( num == 0) 
    return "AA";
  else if( num == 1) 
    return "AR";
  else if( num == 2) 
    return "AN";
  else if( num == 3) 
    return "AD";
  else if( num == 4) 
    return "AC";
  else if( num == 5) 
    return "AQ";
  else if( num == 6) 
    return "AE";
  else if( num == 7) 
    return "AG";
  else if( num == 8) 
    return "AH";
  else if( num == 9) 
    return "AI";
  else if( num == 10) 
    return "AL";
  else if( num == 11) 
    return "AK";
  else if( num == 12) 
    return "AM";
  else if( num == 13) 
    return "AF";
  else if( num == 14) 
    return "AP";
  else if( num == 15) 
    return "AS";
  else if( num == 16) 
    return "AT";
  else if( num == 17) 
    return "AW";
  else if( num == 18) 
    return "AY";
  else if( num == 19) 
    return "AV";
  else if( num == 20) 
    return "RA";
  else if( num == 21) 
    return "RR";
  else if( num == 22) 
    return "RN";
  else if( num == 23) 
    return "RD";
  else if( num == 24) 
    return "RC";
  else if( num == 25) 
    return "RQ";
  else if( num == 26) 
    return "RE";
  else if( num == 27) 
    return "RG";
  else if( num == 28) 
    return "RH";
  else if( num == 29) 
    return "RI";
  else if( num == 30) 
    return "RL";
  else if( num == 31) 
    return "RK";
  else if( num == 32) 
    return "RM";
  else if( num == 33) 
    return "RF";
  else if( num == 34) 
    return "RP";
  else if( num == 35) 
    return "RS";
  else if( num == 36) 
    return "RT";
  else if( num == 37) 
    return "RW";
  else if( num == 38) 
    return "RY";
  else if( num == 39) 
    return "RV";
  else if( num == 40) 
    return "NA";
  else if( num == 41) 
    return "NR";
  else if( num == 42) 
    return "NN";
  else if( num == 43) 
    return "ND";
  else if( num == 44) 
    return "NC";
  else if( num == 45) 
    return "NQ";
  else if( num == 46) 
    return "NE";
  else if( num == 47) 
    return "NG";
  else if( num == 48) 
    return "NH";
  else if( num == 49) 
    return "NI";
  else if( num == 50) 
    return "NL";
  else if( num == 51) 
    return "NK";
  else if( num == 52) 
    return "NM";
  else if( num == 53) 
    return "NF";
  else if( num == 54) 
    return "NP";
  else if( num == 55) 
    return "NS";
  else if( num == 56) 
    return "NT";
  else if( num == 57) 
    return "NW";
  else if( num == 58) 
    return "NY";
  else if( num == 59) 
    return "NV";
  else if( num == 60) 
    return "DA";
  else if( num == 61) 
    return "DR";
  else if( num == 62) 
    return "DN";
  else if( num == 63) 
    return "DD";
  else if( num == 64) 
    return "DC";
  else if( num == 65) 
    return "DQ";
  else if( num == 66) 
    return "DE";
  else if( num == 67) 
    return "DG";
  else if( num == 68) 
    return "DH";
  else if( num == 69) 
    return "DI";
  else if( num == 70) 
    return "DL";
  else if( num == 71) 
    return "DK";
  else if( num == 72) 
    return "DM";
  else if( num == 73) 
    return "DF";
  else if( num == 74) 
    return "DP";
  else if( num == 75) 
    return "DS";
  else if( num == 76) 
    return "DT";
  else if( num == 77) 
    return "DW";
  else if( num == 78) 
    return "DY";
  else if( num == 79) 
    return "DV";
  else if( num == 80) 
    return "CA";
  else if( num == 81) 
    return "CR";
  else if( num == 82) 
    return "CN";
  else if( num == 83) 
    return "CD";
  else if( num == 84) 
    return "CC";
  else if( num == 85) 
    return "CQ";
  else if( num == 86) 
    return "CE";
  else if( num == 87) 
    return "CG";
  else if( num == 88) 
    return "CH";
  else if( num == 89) 
    return "CI";
  else if( num == 90) 
    return "CL";
  else if( num == 91) 
    return "CK";
  else if( num == 92) 
    return "CM";
  else if( num == 93) 
    return "CF";
  else if( num == 94) 
    return "CP";
  else if( num == 95) 
    return "CS";
  else if( num == 96) 
    return "CT";
  else if( num == 97) 
    return "CW";
  else if( num == 98) 
    return "CY";
  else if( num == 99) 
    return "CV";
  else if( num == 100) 
    return "QA";
  else if( num == 101) 
    return "QR";
  else if( num == 102) 
    return "QN";
  else if( num == 103) 
    return "QD";
  else if( num == 104) 
    return "QC";
  else if( num == 105) 
    return "QQ";
  else if( num == 106) 
    return "QE";
  else if( num == 107) 
    return "QG";
  else if( num == 108) 
    return "QH";
  else if( num == 109) 
    return "QI";
  else if( num == 110) 
    return "QL";
  else if( num == 111) 
    return "QK";
  else if( num == 112) 
    return "QM";
  else if( num == 113) 
    return "QF";
  else if( num == 114) 
    return "QP";
  else if( num == 115) 
    return "QS";
  else if( num == 116) 
    return "QT";
  else if( num == 117) 
    return "QW";
  else if( num == 118) 
    return "QY";
  else if( num == 119) 
    return "QV";
  else if( num == 120) 
    return "EA";
  else if( num == 121) 
    return "ER";
  else if( num == 122) 
    return "EN";
  else if( num == 123) 
    return "ED";
  else if( num == 124) 
    return "EC";
  else if( num == 125) 
    return "EQ";
  else if( num == 126) 
    return "EE";
  else if( num == 127) 
    return "EG";
  else if( num == 128) 
    return "EH";
  else if( num == 129) 
    return "EI";
  else if( num == 130) 
    return "EL";
  else if( num == 131) 
    return "EK";
  else if( num == 132) 
    return "EM";
  else if( num == 133) 
    return "EF";
  else if( num == 134) 
    return "EP";
  else if( num == 135) 
    return "ES";
  else if( num == 136) 
    return "ET";
  else if( num == 137) 
    return "EW";
  else if( num == 138) 
    return "EY";
  else if( num == 139) 
    return "EV";
  else if( num == 140) 
    return "GA";
  else if( num == 141) 
    return "GR";
  else if( num == 142) 
    return "GN";
  else if( num == 143) 
    return "GD";
  else if( num == 144) 
    return "GC";
  else if( num == 145) 
    return "GQ";
  else if( num == 146) 
    return "GE";
  else if( num == 147) 
    return "GG";
  else if( num == 148) 
    return "GH";
  else if( num == 149) 
    return "GI";
  else if( num == 150) 
    return "GL";
  else if( num == 151) 
    return "GK";
  else if( num == 152) 
    return "GM";
  else if( num == 153) 
    return "GF";
  else if( num == 154) 
    return "GP";
  else if( num == 155) 
    return "GS";
  else if( num == 156) 
    return "GT";
  else if( num == 157) 
    return "GW";
  else if( num == 158) 
    return "GY";
  else if( num == 159) 
    return "GV";
  else if( num == 160) 
    return "HA";
  else if( num == 161) 
    return "HR";
  else if( num == 162) 
    return "HN";
  else if( num == 163) 
    return "HD";
  else if( num == 164) 
    return "HC";
  else if( num == 165) 
    return "HQ";
  else if( num == 166) 
    return "HE";
  else if( num == 167) 
    return "HG";
  else if( num == 168) 
    return "HH";
  else if( num == 169) 
    return "HI";
  else if( num == 170) 
    return "HL";
  else if( num == 171) 
    return "HK";
  else if( num == 172) 
    return "HM";
  else if( num == 173) 
    return "HF";
  else if( num == 174) 
    return "HP";
  else if( num == 175) 
    return "HS";
  else if( num == 176) 
    return "HT";
  else if( num == 177) 
    return "HW";
  else if( num == 178) 
    return "HY";
  else if( num == 179) 
    return "HV";
  else if( num == 180) 
    return "IA";
  else if( num == 181) 
    return "IR";
  else if( num == 182) 
    return "IN";
  else if( num == 183) 
    return "ID";
  else if( num == 184) 
    return "IC";
  else if( num == 185) 
    return "IQ";
  else if( num == 186) 
    return "IE";
  else if( num == 187) 
    return "IG";
  else if( num == 188) 
    return "IH";
  else if( num == 189) 
    return "II";
  else if( num == 190) 
    return "IL";
  else if( num == 191) 
    return "IK";
  else if( num == 192) 
    return "IM";
  else if( num == 193) 
    return "IF";
  else if( num == 194) 
    return "IP";
  else if( num == 195) 
    return "IS";
  else if( num == 196) 
    return "IT";
  else if( num == 197) 
    return "IW";
  else if( num == 198) 
    return "IY";
  else if( num == 199) 
    return "IV";
  else if( num == 200) 
    return "LA";
  else if( num == 201) 
    return "LR";
  else if( num == 202) 
    return "LN";
  else if( num == 203) 
    return "LD";
  else if( num == 204) 
    return "LC";
  else if( num == 205) 
    return "LQ";
  else if( num == 206) 
    return "LE";
  else if( num == 207) 
    return "LG";
  else if( num == 208) 
    return "LH";
  else if( num == 209) 
    return "LI";
  else if( num == 210) 
    return "LL";
  else if( num == 211) 
    return "LK";
  else if( num == 212) 
    return "LM";
  else if( num == 213) 
    return "LF";
  else if( num == 214) 
    return "LP";
  else if( num == 215) 
    return "LS";
  else if( num == 216) 
    return "LT";
  else if( num == 217) 
    return "LW";
  else if( num == 218) 
    return "LY";
  else if( num == 219) 
    return "LV";
  else if( num == 220) 
    return "KA";
  else if( num == 221) 
    return "KR";
  else if( num == 222) 
    return "KN";
  else if( num == 223) 
    return "KD";
  else if( num == 224) 
    return "KC";
  else if( num == 225) 
    return "KQ";
  else if( num == 226) 
    return "KE";
  else if( num == 227) 
    return "KG";
  else if( num == 228) 
    return "KH";
  else if( num == 229) 
    return "KI";
  else if( num == 230) 
    return "KL";
  else if( num == 231) 
    return "KK";
  else if( num == 232) 
    return "KM";
  else if( num == 233) 
    return "KF";
  else if( num == 234) 
    return "KP";
  else if( num == 235) 
    return "KS";
  else if( num == 236) 
    return "KT";
  else if( num == 237) 
    return "KW";
  else if( num == 238) 
    return "KY";
  else if( num == 239) 
    return "KV";
  else if( num == 240) 
    return "MA";
  else if( num == 241) 
    return "MR";
  else if( num == 242) 
    return "MN";
  else if( num == 243) 
    return "MD";
  else if( num == 244) 
    return "MC";
  else if( num == 245) 
    return "MQ";
  else if( num == 246) 
    return "ME";
  else if( num == 247) 
    return "MG";
  else if( num == 248) 
    return "MH";
  else if( num == 249) 
    return "MI";
  else if( num == 250) 
    return "ML";
  else if( num == 251) 
    return "MK";
  else if( num == 252) 
    return "MM";
  else if( num == 253) 
    return "MF";
  else if( num == 254) 
    return "MP";
  else if( num == 255) 
    return "MS";
  else if( num == 256) 
    return "MT";
  else if( num == 257) 
    return "MW";
  else if( num == 258) 
    return "MY";
  else if( num == 259) 
    return "MV";
  else if( num == 260) 
    return "FA";
  else if( num == 261) 
    return "FR";
  else if( num == 262) 
    return "FN";
  else if( num == 263) 
    return "FD";
  else if( num == 264) 
    return "FC";
  else if( num == 265) 
    return "FQ";
  else if( num == 266) 
    return "FE";
  else if( num == 267) 
    return "FG";
  else if( num == 268) 
    return "FH";
  else if( num == 269) 
    return "FI";
  else if( num == 270) 
    return "FL";
  else if( num == 271) 
    return "FK";
  else if( num == 272) 
    return "FM";
  else if( num == 273) 
    return "FF";
  else if( num == 274) 
    return "FP";
  else if( num == 275) 
    return "FS";
  else if( num == 276) 
    return "FT";
  else if( num == 277) 
    return "FW";
  else if( num == 278) 
    return "FY";
  else if( num == 279) 
    return "FV";
  else if( num == 280) 
    return "PA";
  else if( num == 281) 
    return "PR";
  else if( num == 282) 
    return "PN";
  else if( num == 283) 
    return "PD";
  else if( num == 284) 
    return "PC";
  else if( num == 285) 
    return "PQ";
  else if( num == 286) 
    return "PE";
  else if( num == 287) 
    return "PG";
  else if( num == 288) 
    return "PH";
  else if( num == 289) 
    return "PI";
  else if( num == 290) 
    return "PL";
  else if( num == 291) 
    return "PK";
  else if( num == 292) 
    return "PM";
  else if( num == 293) 
    return "PF";
  else if( num == 294) 
    return "PP";
  else if( num == 295) 
    return "PS";
  else if( num == 296) 
    return "PT";
  else if( num == 297) 
    return "PW";
  else if( num == 298) 
    return "PY";
  else if( num == 299) 
    return "PV";
  else if( num == 300) 
    return "SA";
  else if( num == 301) 
    return "SR";
  else if( num == 302) 
    return "SN";
  else if( num == 303) 
    return "SD";
  else if( num == 304) 
    return "SC";
  else if( num == 305) 
    return "SQ";
  else if( num == 306) 
    return "SE";
  else if( num == 307) 
    return "SG";
  else if( num == 308) 
    return "SH";
  else if( num == 309) 
    return "SI";
  else if( num == 310) 
    return "SL";
  else if( num == 311) 
    return "SK";
  else if( num == 312) 
    return "SM";
  else if( num == 313) 
    return "SF";
  else if( num == 314) 
    return "SP";
  else if( num == 315) 
    return "SS";
  else if( num == 316) 
    return "ST";
  else if( num == 317) 
    return "SW";
  else if( num == 318) 
    return "SY";
  else if( num == 319) 
    return "SV";
  else if( num == 320) 
    return "TA";
  else if( num == 321) 
    return "TR";
  else if( num == 322) 
    return "TN";
  else if( num == 323) 
    return "TD";
  else if( num == 324) 
    return "TC";
  else if( num == 325) 
    return "TQ";
  else if( num == 326) 
    return "TE";
  else if( num == 327) 
    return "TG";
  else if( num == 328) 
    return "TH";
  else if( num == 329) 
    return "TI";
  else if( num == 330) 
    return "TL";
  else if( num == 331) 
    return "TK";
  else if( num == 332) 
    return "TM";
  else if( num == 333) 
    return "TF";
  else if( num == 334) 
    return "TP";
  else if( num == 335) 
    return "TS";
  else if( num == 336) 
    return "TT";
  else if( num == 337) 
    return "TW";
  else if( num == 338) 
    return "TY";
  else if( num == 339) 
    return "TV";
  else if( num == 340) 
    return "WA";
  else if( num == 341) 
    return "WR";
  else if( num == 342) 
    return "WN";
  else if( num == 343) 
    return "WD";
  else if( num == 344) 
    return "WC";
  else if( num == 345) 
    return "WQ";
  else if( num == 346) 
    return "WE";
  else if( num == 347) 
    return "WG";
  else if( num == 348) 
    return "WH";
  else if( num == 349) 
    return "WI";
  else if( num == 350) 
    return "WL";
  else if( num == 351) 
    return "WK";
  else if( num == 352) 
    return "WM";
  else if( num == 353) 
    return "WF";
  else if( num == 354) 
    return "WP";
  else if( num == 355) 
    return "WS";
  else if( num == 356) 
    return "WT";
  else if( num == 357) 
    return "WW";
  else if( num == 358) 
    return "WY";
  else if( num == 359) 
    return "WV";
  else if( num == 360) 
    return "YA";
  else if( num == 361) 
    return "YR";
  else if( num == 362) 
    return "YN";
  else if( num == 363) 
    return "YD";
  else if( num == 364) 
    return "YC";
  else if( num == 365) 
    return "YQ";
  else if( num == 366) 
    return "YE";
  else if( num == 367) 
    return "YG";
  else if( num == 368) 
    return "YH";
  else if( num == 369) 
    return "YI";
  else if( num == 370) 
    return "YL";
  else if( num == 371) 
    return "YK";
  else if( num == 372) 
    return "YM";
  else if( num == 373) 
    return "YF";
  else if( num == 374) 
    return "YP";
  else if( num == 375) 
    return "YS";
  else if( num == 376) 
    return "YT";
  else if( num == 377) 
    return "YW";
  else if( num == 378) 
    return "YY";
  else if( num == 379) 
    return "YV";
  else if( num == 380) 
    return "VA";
  else if( num == 381) 
    return "VR";
  else if( num == 382) 
    return "VN";
  else if( num == 383) 
    return "VD";
  else if( num == 384) 
    return "VC";
  else if( num == 385) 
    return "VQ";
  else if( num == 386) 
    return "VE";
  else if( num == 387) 
    return "VG";
  else if( num == 388) 
    return "VH";
  else if( num == 389) 
    return "VI";
  else if( num == 390) 
    return "VL";
  else if( num == 391) 
    return "VK";
  else if( num == 392) 
    return "VM";
  else if( num == 393) 
    return "VF";
  else if( num == 394) 
    return "VP";
  else if( num == 395) 
    return "VS";
  else if( num == 396) 
    return "VT";
  else if( num == 397) 
    return "VW";
  else if( num == 398) 
    return "VY";
  else if( num == 399) 
    return "VV";
  return "AA";
} 
