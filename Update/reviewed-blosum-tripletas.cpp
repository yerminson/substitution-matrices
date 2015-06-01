#include <math.h>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <vector>

using namespace std;

// Number of aminoacids 
#define AAS 8000
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
        for (i = 0; i < (int)strlen(ptr); i+=3) {
          string test(ptr);
          test = test.substr(i,3);
          if(test.length() < 3)
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
  if ( aa == "AAA") 
    return 0;
  else if( aa == "AAR")  
    return 1;
  else if( aa == "AAN")  
    return 2;
  else if( aa == "AAD")  
    return 3;
  else if( aa == "AAC")  
    return 4;
  else if( aa == "AAQ")  
    return 5;
  else if( aa == "AAE")  
    return 6;
  else if( aa == "AAG")  
    return 7;
  else if( aa == "AAH")  
    return 8;
  else if( aa == "AAI")  
    return 9;
  else if( aa == "AAL")  
    return 10;
  else if( aa == "AAK")  
    return 11;
  else if( aa == "AAM")  
    return 12;
  else if( aa == "AAF")  
    return 13;
  else if( aa == "AAP")  
    return 14;
  else if( aa == "AAS")  
    return 15;
  else if( aa == "AAT")  
    return 16;
  else if( aa == "AAW")  
    return 17;
  else if( aa == "AAY")  
    return 18;
  else if( aa == "AAV")  
    return 19;
  else if( aa == "ARA")  
    return 20;
  else if( aa == "ARR")  
    return 21;
  else if( aa == "ARN")  
    return 22;
  else if( aa == "ARD")  
    return 23;
  else if( aa == "ARC")  
    return 24;
  else if( aa == "ARQ")  
    return 25;
  else if( aa == "ARE")  
    return 26;
  else if( aa == "ARG")  
    return 27;
  else if( aa == "ARH")  
    return 28;
  else if( aa == "ARI")  
    return 29;
  else if( aa == "ARL")  
    return 30;
  else if( aa == "ARK")  
    return 31;
  else if( aa == "ARM")  
    return 32;
  else if( aa == "ARF")  
    return 33;
  else if( aa == "ARP")  
    return 34;
  else if( aa == "ARS")  
    return 35;
  else if( aa == "ART")  
    return 36;
  else if( aa == "ARW")  
    return 37;
  else if( aa == "ARY")  
    return 38;
  else if( aa == "ARV")  
    return 39;
  else if( aa == "ANA")  
    return 40;
  else if( aa == "ANR")  
    return 41;
  else if( aa == "ANN")  
    return 42;
  else if( aa == "AND")  
    return 43;
  else if( aa == "ANC")  
    return 44;
  else if( aa == "ANQ")  
    return 45;
  else if( aa == "ANE")  
    return 46;
  else if( aa == "ANG")  
    return 47;
  else if( aa == "ANH")  
    return 48;
  else if( aa == "ANI")  
    return 49;
  else if( aa == "ANL")  
    return 50;
  else if( aa == "ANK")  
    return 51;
  else if( aa == "ANM")  
    return 52;
  else if( aa == "ANF")  
    return 53;
  else if( aa == "ANP")  
    return 54;
  else if( aa == "ANS")  
    return 55;
  else if( aa == "ANT")  
    return 56;
  else if( aa == "ANW")  
    return 57;
  else if( aa == "ANY")  
    return 58;
  else if( aa == "ANV")  
    return 59;
  else if( aa == "ADA")  
    return 60;
  else if( aa == "ADR")  
    return 61;
  else if( aa == "ADN")  
    return 62;
  else if( aa == "ADD")  
    return 63;
  else if( aa == "ADC")  
    return 64;
  else if( aa == "ADQ")  
    return 65;
  else if( aa == "ADE")  
    return 66;
  else if( aa == "ADG")  
    return 67;
  else if( aa == "ADH")  
    return 68;
  else if( aa == "ADI")  
    return 69;
  else if( aa == "ADL")  
    return 70;
  else if( aa == "ADK")  
    return 71;
  else if( aa == "ADM")  
    return 72;
  else if( aa == "ADF")  
    return 73;
  else if( aa == "ADP")  
    return 74;
  else if( aa == "ADS")  
    return 75;
  else if( aa == "ADT")  
    return 76;
  else if( aa == "ADW")  
    return 77;
  else if( aa == "ADY")  
    return 78;
  else if( aa == "ADV")  
    return 79;
  else if( aa == "ACA")  
    return 80;
  else if( aa == "ACR")  
    return 81;
  else if( aa == "ACN")  
    return 82;
  else if( aa == "ACD")  
    return 83;
  else if( aa == "ACC")  
    return 84;
  else if( aa == "ACQ")  
    return 85;
  else if( aa == "ACE")  
    return 86;
  else if( aa == "ACG")  
    return 87;
  else if( aa == "ACH")  
    return 88;
  else if( aa == "ACI")  
    return 89;
  else if( aa == "ACL")  
    return 90;
  else if( aa == "ACK")  
    return 91;
  else if( aa == "ACM")  
    return 92;
  else if( aa == "ACF")  
    return 93;
  else if( aa == "ACP")  
    return 94;
  else if( aa == "ACS")  
    return 95;
  else if( aa == "ACT")  
    return 96;
  else if( aa == "ACW")  
    return 97;
  else if( aa == "ACY")  
    return 98;
  else if( aa == "ACV")  
    return 99;
  else if( aa == "AQA")  
    return 100;
  else if( aa == "AQR")  
    return 101;
  else if( aa == "AQN")  
    return 102;
  else if( aa == "AQD")  
    return 103;
  else if( aa == "AQC")  
    return 104;
  else if( aa == "AQQ")  
    return 105;
  else if( aa == "AQE")  
    return 106;
  else if( aa == "AQG")  
    return 107;
  else if( aa == "AQH")  
    return 108;
  else if( aa == "AQI")  
    return 109;
  else if( aa == "AQL")  
    return 110;
  else if( aa == "AQK")  
    return 111;
  else if( aa == "AQM")  
    return 112;
  else if( aa == "AQF")  
    return 113;
  else if( aa == "AQP")  
    return 114;
  else if( aa == "AQS")  
    return 115;
  else if( aa == "AQT")  
    return 116;
  else if( aa == "AQW")  
    return 117;
  else if( aa == "AQY")  
    return 118;
  else if( aa == "AQV")  
    return 119;
  else if( aa == "AEA")  
    return 120;
  else if( aa == "AER")  
    return 121;
  else if( aa == "AEN")  
    return 122;
  else if( aa == "AED")  
    return 123;
  else if( aa == "AEC")  
    return 124;
  else if( aa == "AEQ")  
    return 125;
  else if( aa == "AEE")  
    return 126;
  else if( aa == "AEG")  
    return 127;
  else if( aa == "AEH")  
    return 128;
  else if( aa == "AEI")  
    return 129;
  else if( aa == "AEL")  
    return 130;
  else if( aa == "AEK")  
    return 131;
  else if( aa == "AEM")  
    return 132;
  else if( aa == "AEF")  
    return 133;
  else if( aa == "AEP")  
    return 134;
  else if( aa == "AES")  
    return 135;
  else if( aa == "AET")  
    return 136;
  else if( aa == "AEW")  
    return 137;
  else if( aa == "AEY")  
    return 138;
  else if( aa == "AEV")  
    return 139;
  else if( aa == "AGA")  
    return 140;
  else if( aa == "AGR")  
    return 141;
  else if( aa == "AGN")  
    return 142;
  else if( aa == "AGD")  
    return 143;
  else if( aa == "AGC")  
    return 144;
  else if( aa == "AGQ")  
    return 145;
  else if( aa == "AGE")  
    return 146;
  else if( aa == "AGG")  
    return 147;
  else if( aa == "AGH")  
    return 148;
  else if( aa == "AGI")  
    return 149;
  else if( aa == "AGL")  
    return 150;
  else if( aa == "AGK")  
    return 151;
  else if( aa == "AGM")  
    return 152;
  else if( aa == "AGF")  
    return 153;
  else if( aa == "AGP")  
    return 154;
  else if( aa == "AGS")  
    return 155;
  else if( aa == "AGT")  
    return 156;
  else if( aa == "AGW")  
    return 157;
  else if( aa == "AGY")  
    return 158;
  else if( aa == "AGV")  
    return 159;
  else if( aa == "AHA")  
    return 160;
  else if( aa == "AHR")  
    return 161;
  else if( aa == "AHN")  
    return 162;
  else if( aa == "AHD")  
    return 163;
  else if( aa == "AHC")  
    return 164;
  else if( aa == "AHQ")  
    return 165;
  else if( aa == "AHE")  
    return 166;
  else if( aa == "AHG")  
    return 167;
  else if( aa == "AHH")  
    return 168;
  else if( aa == "AHI")  
    return 169;
  else if( aa == "AHL")  
    return 170;
  else if( aa == "AHK")  
    return 171;
  else if( aa == "AHM")  
    return 172;
  else if( aa == "AHF")  
    return 173;
  else if( aa == "AHP")  
    return 174;
  else if( aa == "AHS")  
    return 175;
  else if( aa == "AHT")  
    return 176;
  else if( aa == "AHW")  
    return 177;
  else if( aa == "AHY")  
    return 178;
  else if( aa == "AHV")  
    return 179;
  else if( aa == "AIA")  
    return 180;
  else if( aa == "AIR")  
    return 181;
  else if( aa == "AIN")  
    return 182;
  else if( aa == "AID")  
    return 183;
  else if( aa == "AIC")  
    return 184;
  else if( aa == "AIQ")  
    return 185;
  else if( aa == "AIE")  
    return 186;
  else if( aa == "AIG")  
    return 187;
  else if( aa == "AIH")  
    return 188;
  else if( aa == "AII")  
    return 189;
  else if( aa == "AIL")  
    return 190;
  else if( aa == "AIK")  
    return 191;
  else if( aa == "AIM")  
    return 192;
  else if( aa == "AIF")  
    return 193;
  else if( aa == "AIP")  
    return 194;
  else if( aa == "AIS")  
    return 195;
  else if( aa == "AIT")  
    return 196;
  else if( aa == "AIW")  
    return 197;
  else if( aa == "AIY")  
    return 198;
  else if( aa == "AIV")  
    return 199;
  else if( aa == "ALA")  
    return 200;
  else if( aa == "ALR")  
    return 201;
  else if( aa == "ALN")  
    return 202;
  else if( aa == "ALD")  
    return 203;
  else if( aa == "ALC")  
    return 204;
  else if( aa == "ALQ")  
    return 205;
  else if( aa == "ALE")  
    return 206;
  else if( aa == "ALG")  
    return 207;
  else if( aa == "ALH")  
    return 208;
  else if( aa == "ALI")  
    return 209;
  else if( aa == "ALL")  
    return 210;
  else if( aa == "ALK")  
    return 211;
  else if( aa == "ALM")  
    return 212;
  else if( aa == "ALF")  
    return 213;
  else if( aa == "ALP")  
    return 214;
  else if( aa == "ALS")  
    return 215;
  else if( aa == "ALT")  
    return 216;
  else if( aa == "ALW")  
    return 217;
  else if( aa == "ALY")  
    return 218;
  else if( aa == "ALV")  
    return 219;
  else if( aa == "AKA")  
    return 220;
  else if( aa == "AKR")  
    return 221;
  else if( aa == "AKN")  
    return 222;
  else if( aa == "AKD")  
    return 223;
  else if( aa == "AKC")  
    return 224;
  else if( aa == "AKQ")  
    return 225;
  else if( aa == "AKE")  
    return 226;
  else if( aa == "AKG")  
    return 227;
  else if( aa == "AKH")  
    return 228;
  else if( aa == "AKI")  
    return 229;
  else if( aa == "AKL")  
    return 230;
  else if( aa == "AKK")  
    return 231;
  else if( aa == "AKM")  
    return 232;
  else if( aa == "AKF")  
    return 233;
  else if( aa == "AKP")  
    return 234;
  else if( aa == "AKS")  
    return 235;
  else if( aa == "AKT")  
    return 236;
  else if( aa == "AKW")  
    return 237;
  else if( aa == "AKY")  
    return 238;
  else if( aa == "AKV")  
    return 239;
  else if( aa == "AMA")  
    return 240;
  else if( aa == "AMR")  
    return 241;
  else if( aa == "AMN")  
    return 242;
  else if( aa == "AMD")  
    return 243;
  else if( aa == "AMC")  
    return 244;
  else if( aa == "AMQ")  
    return 245;
  else if( aa == "AME")  
    return 246;
  else if( aa == "AMG")  
    return 247;
  else if( aa == "AMH")  
    return 248;
  else if( aa == "AMI")  
    return 249;
  else if( aa == "AML")  
    return 250;
  else if( aa == "AMK")  
    return 251;
  else if( aa == "AMM")  
    return 252;
  else if( aa == "AMF")  
    return 253;
  else if( aa == "AMP")  
    return 254;
  else if( aa == "AMS")  
    return 255;
  else if( aa == "AMT")  
    return 256;
  else if( aa == "AMW")  
    return 257;
  else if( aa == "AMY")  
    return 258;
  else if( aa == "AMV")  
    return 259;
  else if( aa == "AFA")  
    return 260;
  else if( aa == "AFR")  
    return 261;
  else if( aa == "AFN")  
    return 262;
  else if( aa == "AFD")  
    return 263;
  else if( aa == "AFC")  
    return 264;
  else if( aa == "AFQ")  
    return 265;
  else if( aa == "AFE")  
    return 266;
  else if( aa == "AFG")  
    return 267;
  else if( aa == "AFH")  
    return 268;
  else if( aa == "AFI")  
    return 269;
  else if( aa == "AFL")  
    return 270;
  else if( aa == "AFK")  
    return 271;
  else if( aa == "AFM")  
    return 272;
  else if( aa == "AFF")  
    return 273;
  else if( aa == "AFP")  
    return 274;
  else if( aa == "AFS")  
    return 275;
  else if( aa == "AFT")  
    return 276;
  else if( aa == "AFW")  
    return 277;
  else if( aa == "AFY")  
    return 278;
  else if( aa == "AFV")  
    return 279;
  else if( aa == "APA")  
    return 280;
  else if( aa == "APR")  
    return 281;
  else if( aa == "APN")  
    return 282;
  else if( aa == "APD")  
    return 283;
  else if( aa == "APC")  
    return 284;
  else if( aa == "APQ")  
    return 285;
  else if( aa == "APE")  
    return 286;
  else if( aa == "APG")  
    return 287;
  else if( aa == "APH")  
    return 288;
  else if( aa == "API")  
    return 289;
  else if( aa == "APL")  
    return 290;
  else if( aa == "APK")  
    return 291;
  else if( aa == "APM")  
    return 292;
  else if( aa == "APF")  
    return 293;
  else if( aa == "APP")  
    return 294;
  else if( aa == "APS")  
    return 295;
  else if( aa == "APT")  
    return 296;
  else if( aa == "APW")  
    return 297;
  else if( aa == "APY")  
    return 298;
  else if( aa == "APV")  
    return 299;
  else if( aa == "ASA")  
    return 300;
  else if( aa == "ASR")  
    return 301;
  else if( aa == "ASN")  
    return 302;
  else if( aa == "ASD")  
    return 303;
  else if( aa == "ASC")  
    return 304;
  else if( aa == "ASQ")  
    return 305;
  else if( aa == "ASE")  
    return 306;
  else if( aa == "ASG")  
    return 307;
  else if( aa == "ASH")  
    return 308;
  else if( aa == "ASI")  
    return 309;
  else if( aa == "ASL")  
    return 310;
  else if( aa == "ASK")  
    return 311;
  else if( aa == "ASM")  
    return 312;
  else if( aa == "ASF")  
    return 313;
  else if( aa == "ASP")  
    return 314;
  else if( aa == "ASS")  
    return 315;
  else if( aa == "AST")  
    return 316;
  else if( aa == "ASW")  
    return 317;
  else if( aa == "ASY")  
    return 318;
  else if( aa == "ASV")  
    return 319;
  else if( aa == "ATA")  
    return 320;
  else if( aa == "ATR")  
    return 321;
  else if( aa == "ATN")  
    return 322;
  else if( aa == "ATD")  
    return 323;
  else if( aa == "ATC")  
    return 324;
  else if( aa == "ATQ")  
    return 325;
  else if( aa == "ATE")  
    return 326;
  else if( aa == "ATG")  
    return 327;
  else if( aa == "ATH")  
    return 328;
  else if( aa == "ATI")  
    return 329;
  else if( aa == "ATL")  
    return 330;
  else if( aa == "ATK")  
    return 331;
  else if( aa == "ATM")  
    return 332;
  else if( aa == "ATF")  
    return 333;
  else if( aa == "ATP")  
    return 334;
  else if( aa == "ATS")  
    return 335;
  else if( aa == "ATT")  
    return 336;
  else if( aa == "ATW")  
    return 337;
  else if( aa == "ATY")  
    return 338;
  else if( aa == "ATV")  
    return 339;
  else if( aa == "AWA")  
    return 340;
  else if( aa == "AWR")  
    return 341;
  else if( aa == "AWN")  
    return 342;
  else if( aa == "AWD")  
    return 343;
  else if( aa == "AWC")  
    return 344;
  else if( aa == "AWQ")  
    return 345;
  else if( aa == "AWE")  
    return 346;
  else if( aa == "AWG")  
    return 347;
  else if( aa == "AWH")  
    return 348;
  else if( aa == "AWI")  
    return 349;
  else if( aa == "AWL")  
    return 350;
  else if( aa == "AWK")  
    return 351;
  else if( aa == "AWM")  
    return 352;
  else if( aa == "AWF")  
    return 353;
  else if( aa == "AWP")  
    return 354;
  else if( aa == "AWS")  
    return 355;
  else if( aa == "AWT")  
    return 356;
  else if( aa == "AWW")  
    return 357;
  else if( aa == "AWY")  
    return 358;
  else if( aa == "AWV")  
    return 359;
  else if( aa == "AYA")  
    return 360;
  else if( aa == "AYR")  
    return 361;
  else if( aa == "AYN")  
    return 362;
  else if( aa == "AYD")  
    return 363;
  else if( aa == "AYC")  
    return 364;
  else if( aa == "AYQ")  
    return 365;
  else if( aa == "AYE")  
    return 366;
  else if( aa == "AYG")  
    return 367;
  else if( aa == "AYH")  
    return 368;
  else if( aa == "AYI")  
    return 369;
  else if( aa == "AYL")  
    return 370;
  else if( aa == "AYK")  
    return 371;
  else if( aa == "AYM")  
    return 372;
  else if( aa == "AYF")  
    return 373;
  else if( aa == "AYP")  
    return 374;
  else if( aa == "AYS")  
    return 375;
  else if( aa == "AYT")  
    return 376;
  else if( aa == "AYW")  
    return 377;
  else if( aa == "AYY")  
    return 378;
  else if( aa == "AYV")  
    return 379;
  else if( aa == "AVA")  
    return 380;
  else if( aa == "AVR")  
    return 381;
  else if( aa == "AVN")  
    return 382;
  else if( aa == "AVD")  
    return 383;
  else if( aa == "AVC")  
    return 384;
  else if( aa == "AVQ")  
    return 385;
  else if( aa == "AVE")  
    return 386;
  else if( aa == "AVG")  
    return 387;
  else if( aa == "AVH")  
    return 388;
  else if( aa == "AVI")  
    return 389;
  else if( aa == "AVL")  
    return 390;
  else if( aa == "AVK")  
    return 391;
  else if( aa == "AVM")  
    return 392;
  else if( aa == "AVF")  
    return 393;
  else if( aa == "AVP")  
    return 394;
  else if( aa == "AVS")  
    return 395;
  else if( aa == "AVT")  
    return 396;
  else if( aa == "AVW")  
    return 397;
  else if( aa == "AVY")  
    return 398;
  else if( aa == "AVV")  
    return 399;
  else if( aa == "RAA")  
    return 400;
  else if( aa == "RAR")  
    return 401;
  else if( aa == "RAN")  
    return 402;
  else if( aa == "RAD")  
    return 403;
  else if( aa == "RAC")  
    return 404;
  else if( aa == "RAQ")  
    return 405;
  else if( aa == "RAE")  
    return 406;
  else if( aa == "RAG")  
    return 407;
  else if( aa == "RAH")  
    return 408;
  else if( aa == "RAI")  
    return 409;
  else if( aa == "RAL")  
    return 410;
  else if( aa == "RAK")  
    return 411;
  else if( aa == "RAM")  
    return 412;
  else if( aa == "RAF")  
    return 413;
  else if( aa == "RAP")  
    return 414;
  else if( aa == "RAS")  
    return 415;
  else if( aa == "RAT")  
    return 416;
  else if( aa == "RAW")  
    return 417;
  else if( aa == "RAY")  
    return 418;
  else if( aa == "RAV")  
    return 419;
  else if( aa == "RRA")  
    return 420;
  else if( aa == "RRR")  
    return 421;
  else if( aa == "RRN")  
    return 422;
  else if( aa == "RRD")  
    return 423;
  else if( aa == "RRC")  
    return 424;
  else if( aa == "RRQ")  
    return 425;
  else if( aa == "RRE")  
    return 426;
  else if( aa == "RRG")  
    return 427;
  else if( aa == "RRH")  
    return 428;
  else if( aa == "RRI")  
    return 429;
  else if( aa == "RRL")  
    return 430;
  else if( aa == "RRK")  
    return 431;
  else if( aa == "RRM")  
    return 432;
  else if( aa == "RRF")  
    return 433;
  else if( aa == "RRP")  
    return 434;
  else if( aa == "RRS")  
    return 435;
  else if( aa == "RRT")  
    return 436;
  else if( aa == "RRW")  
    return 437;
  else if( aa == "RRY")  
    return 438;
  else if( aa == "RRV")  
    return 439;
  else if( aa == "RNA")  
    return 440;
  else if( aa == "RNR")  
    return 441;
  else if( aa == "RNN")  
    return 442;
  else if( aa == "RND")  
    return 443;
  else if( aa == "RNC")  
    return 444;
  else if( aa == "RNQ")  
    return 445;
  else if( aa == "RNE")  
    return 446;
  else if( aa == "RNG")  
    return 447;
  else if( aa == "RNH")  
    return 448;
  else if( aa == "RNI")  
    return 449;
  else if( aa == "RNL")  
    return 450;
  else if( aa == "RNK")  
    return 451;
  else if( aa == "RNM")  
    return 452;
  else if( aa == "RNF")  
    return 453;
  else if( aa == "RNP")  
    return 454;
  else if( aa == "RNS")  
    return 455;
  else if( aa == "RNT")  
    return 456;
  else if( aa == "RNW")  
    return 457;
  else if( aa == "RNY")  
    return 458;
  else if( aa == "RNV")  
    return 459;
  else if( aa == "RDA")  
    return 460;
  else if( aa == "RDR")  
    return 461;
  else if( aa == "RDN")  
    return 462;
  else if( aa == "RDD")  
    return 463;
  else if( aa == "RDC")  
    return 464;
  else if( aa == "RDQ")  
    return 465;
  else if( aa == "RDE")  
    return 466;
  else if( aa == "RDG")  
    return 467;
  else if( aa == "RDH")  
    return 468;
  else if( aa == "RDI")  
    return 469;
  else if( aa == "RDL")  
    return 470;
  else if( aa == "RDK")  
    return 471;
  else if( aa == "RDM")  
    return 472;
  else if( aa == "RDF")  
    return 473;
  else if( aa == "RDP")  
    return 474;
  else if( aa == "RDS")  
    return 475;
  else if( aa == "RDT")  
    return 476;
  else if( aa == "RDW")  
    return 477;
  else if( aa == "RDY")  
    return 478;
  else if( aa == "RDV")  
    return 479;
  else if( aa == "RCA")  
    return 480;
  else if( aa == "RCR")  
    return 481;
  else if( aa == "RCN")  
    return 482;
  else if( aa == "RCD")  
    return 483;
  else if( aa == "RCC")  
    return 484;
  else if( aa == "RCQ")  
    return 485;
  else if( aa == "RCE")  
    return 486;
  else if( aa == "RCG")  
    return 487;
  else if( aa == "RCH")  
    return 488;
  else if( aa == "RCI")  
    return 489;
  else if( aa == "RCL")  
    return 490;
  else if( aa == "RCK")  
    return 491;
  else if( aa == "RCM")  
    return 492;
  else if( aa == "RCF")  
    return 493;
  else if( aa == "RCP")  
    return 494;
  else if( aa == "RCS")  
    return 495;
  else if( aa == "RCT")  
    return 496;
  else if( aa == "RCW")  
    return 497;
  else if( aa == "RCY")  
    return 498;
  else if( aa == "RCV")  
    return 499;
  else if( aa == "RQA")  
    return 500;
  else if( aa == "RQR")  
    return 501;
  else if( aa == "RQN")  
    return 502;
  else if( aa == "RQD")  
    return 503;
  else if( aa == "RQC")  
    return 504;
  else if( aa == "RQQ")  
    return 505;
  else if( aa == "RQE")  
    return 506;
  else if( aa == "RQG")  
    return 507;
  else if( aa == "RQH")  
    return 508;
  else if( aa == "RQI")  
    return 509;
  else if( aa == "RQL")  
    return 510;
  else if( aa == "RQK")  
    return 511;
  else if( aa == "RQM")  
    return 512;
  else if( aa == "RQF")  
    return 513;
  else if( aa == "RQP")  
    return 514;
  else if( aa == "RQS")  
    return 515;
  else if( aa == "RQT")  
    return 516;
  else if( aa == "RQW")  
    return 517;
  else if( aa == "RQY")  
    return 518;
  else if( aa == "RQV")  
    return 519;
  else if( aa == "REA")  
    return 520;
  else if( aa == "RER")  
    return 521;
  else if( aa == "REN")  
    return 522;
  else if( aa == "RED")  
    return 523;
  else if( aa == "REC")  
    return 524;
  else if( aa == "REQ")  
    return 525;
  else if( aa == "REE")  
    return 526;
  else if( aa == "REG")  
    return 527;
  else if( aa == "REH")  
    return 528;
  else if( aa == "REI")  
    return 529;
  else if( aa == "REL")  
    return 530;
  else if( aa == "REK")  
    return 531;
  else if( aa == "REM")  
    return 532;
  else if( aa == "REF")  
    return 533;
  else if( aa == "REP")  
    return 534;
  else if( aa == "RES")  
    return 535;
  else if( aa == "RET")  
    return 536;
  else if( aa == "REW")  
    return 537;
  else if( aa == "REY")  
    return 538;
  else if( aa == "REV")  
    return 539;
  else if( aa == "RGA")  
    return 540;
  else if( aa == "RGR")  
    return 541;
  else if( aa == "RGN")  
    return 542;
  else if( aa == "RGD")  
    return 543;
  else if( aa == "RGC")  
    return 544;
  else if( aa == "RGQ")  
    return 545;
  else if( aa == "RGE")  
    return 546;
  else if( aa == "RGG")  
    return 547;
  else if( aa == "RGH")  
    return 548;
  else if( aa == "RGI")  
    return 549;
  else if( aa == "RGL")  
    return 550;
  else if( aa == "RGK")  
    return 551;
  else if( aa == "RGM")  
    return 552;
  else if( aa == "RGF")  
    return 553;
  else if( aa == "RGP")  
    return 554;
  else if( aa == "RGS")  
    return 555;
  else if( aa == "RGT")  
    return 556;
  else if( aa == "RGW")  
    return 557;
  else if( aa == "RGY")  
    return 558;
  else if( aa == "RGV")  
    return 559;
  else if( aa == "RHA")  
    return 560;
  else if( aa == "RHR")  
    return 561;
  else if( aa == "RHN")  
    return 562;
  else if( aa == "RHD")  
    return 563;
  else if( aa == "RHC")  
    return 564;
  else if( aa == "RHQ")  
    return 565;
  else if( aa == "RHE")  
    return 566;
  else if( aa == "RHG")  
    return 567;
  else if( aa == "RHH")  
    return 568;
  else if( aa == "RHI")  
    return 569;
  else if( aa == "RHL")  
    return 570;
  else if( aa == "RHK")  
    return 571;
  else if( aa == "RHM")  
    return 572;
  else if( aa == "RHF")  
    return 573;
  else if( aa == "RHP")  
    return 574;
  else if( aa == "RHS")  
    return 575;
  else if( aa == "RHT")  
    return 576;
  else if( aa == "RHW")  
    return 577;
  else if( aa == "RHY")  
    return 578;
  else if( aa == "RHV")  
    return 579;
  else if( aa == "RIA")  
    return 580;
  else if( aa == "RIR")  
    return 581;
  else if( aa == "RIN")  
    return 582;
  else if( aa == "RID")  
    return 583;
  else if( aa == "RIC")  
    return 584;
  else if( aa == "RIQ")  
    return 585;
  else if( aa == "RIE")  
    return 586;
  else if( aa == "RIG")  
    return 587;
  else if( aa == "RIH")  
    return 588;
  else if( aa == "RII")  
    return 589;
  else if( aa == "RIL")  
    return 590;
  else if( aa == "RIK")  
    return 591;
  else if( aa == "RIM")  
    return 592;
  else if( aa == "RIF")  
    return 593;
  else if( aa == "RIP")  
    return 594;
  else if( aa == "RIS")  
    return 595;
  else if( aa == "RIT")  
    return 596;
  else if( aa == "RIW")  
    return 597;
  else if( aa == "RIY")  
    return 598;
  else if( aa == "RIV")  
    return 599;
  else if( aa == "RLA")  
    return 600;
  else if( aa == "RLR")  
    return 601;
  else if( aa == "RLN")  
    return 602;
  else if( aa == "RLD")  
    return 603;
  else if( aa == "RLC")  
    return 604;
  else if( aa == "RLQ")  
    return 605;
  else if( aa == "RLE")  
    return 606;
  else if( aa == "RLG")  
    return 607;
  else if( aa == "RLH")  
    return 608;
  else if( aa == "RLI")  
    return 609;
  else if( aa == "RLL")  
    return 610;
  else if( aa == "RLK")  
    return 611;
  else if( aa == "RLM")  
    return 612;
  else if( aa == "RLF")  
    return 613;
  else if( aa == "RLP")  
    return 614;
  else if( aa == "RLS")  
    return 615;
  else if( aa == "RLT")  
    return 616;
  else if( aa == "RLW")  
    return 617;
  else if( aa == "RLY")  
    return 618;
  else if( aa == "RLV")  
    return 619;
  else if( aa == "RKA")  
    return 620;
  else if( aa == "RKR")  
    return 621;
  else if( aa == "RKN")  
    return 622;
  else if( aa == "RKD")  
    return 623;
  else if( aa == "RKC")  
    return 624;
  else if( aa == "RKQ")  
    return 625;
  else if( aa == "RKE")  
    return 626;
  else if( aa == "RKG")  
    return 627;
  else if( aa == "RKH")  
    return 628;
  else if( aa == "RKI")  
    return 629;
  else if( aa == "RKL")  
    return 630;
  else if( aa == "RKK")  
    return 631;
  else if( aa == "RKM")  
    return 632;
  else if( aa == "RKF")  
    return 633;
  else if( aa == "RKP")  
    return 634;
  else if( aa == "RKS")  
    return 635;
  else if( aa == "RKT")  
    return 636;
  else if( aa == "RKW")  
    return 637;
  else if( aa == "RKY")  
    return 638;
  else if( aa == "RKV")  
    return 639;
  else if( aa == "RMA")  
    return 640;
  else if( aa == "RMR")  
    return 641;
  else if( aa == "RMN")  
    return 642;
  else if( aa == "RMD")  
    return 643;
  else if( aa == "RMC")  
    return 644;
  else if( aa == "RMQ")  
    return 645;
  else if( aa == "RME")  
    return 646;
  else if( aa == "RMG")  
    return 647;
  else if( aa == "RMH")  
    return 648;
  else if( aa == "RMI")  
    return 649;
  else if( aa == "RML")  
    return 650;
  else if( aa == "RMK")  
    return 651;
  else if( aa == "RMM")  
    return 652;
  else if( aa == "RMF")  
    return 653;
  else if( aa == "RMP")  
    return 654;
  else if( aa == "RMS")  
    return 655;
  else if( aa == "RMT")  
    return 656;
  else if( aa == "RMW")  
    return 657;
  else if( aa == "RMY")  
    return 658;
  else if( aa == "RMV")  
    return 659;
  else if( aa == "RFA")  
    return 660;
  else if( aa == "RFR")  
    return 661;
  else if( aa == "RFN")  
    return 662;
  else if( aa == "RFD")  
    return 663;
  else if( aa == "RFC")  
    return 664;
  else if( aa == "RFQ")  
    return 665;
  else if( aa == "RFE")  
    return 666;
  else if( aa == "RFG")  
    return 667;
  else if( aa == "RFH")  
    return 668;
  else if( aa == "RFI")  
    return 669;
  else if( aa == "RFL")  
    return 670;
  else if( aa == "RFK")  
    return 671;
  else if( aa == "RFM")  
    return 672;
  else if( aa == "RFF")  
    return 673;
  else if( aa == "RFP")  
    return 674;
  else if( aa == "RFS")  
    return 675;
  else if( aa == "RFT")  
    return 676;
  else if( aa == "RFW")  
    return 677;
  else if( aa == "RFY")  
    return 678;
  else if( aa == "RFV")  
    return 679;
  else if( aa == "RPA")  
    return 680;
  else if( aa == "RPR")  
    return 681;
  else if( aa == "RPN")  
    return 682;
  else if( aa == "RPD")  
    return 683;
  else if( aa == "RPC")  
    return 684;
  else if( aa == "RPQ")  
    return 685;
  else if( aa == "RPE")  
    return 686;
  else if( aa == "RPG")  
    return 687;
  else if( aa == "RPH")  
    return 688;
  else if( aa == "RPI")  
    return 689;
  else if( aa == "RPL")  
    return 690;
  else if( aa == "RPK")  
    return 691;
  else if( aa == "RPM")  
    return 692;
  else if( aa == "RPF")  
    return 693;
  else if( aa == "RPP")  
    return 694;
  else if( aa == "RPS")  
    return 695;
  else if( aa == "RPT")  
    return 696;
  else if( aa == "RPW")  
    return 697;
  else if( aa == "RPY")  
    return 698;
  else if( aa == "RPV")  
    return 699;
  else if( aa == "RSA")  
    return 700;
  else if( aa == "RSR")  
    return 701;
  else if( aa == "RSN")  
    return 702;
  else if( aa == "RSD")  
    return 703;
  else if( aa == "RSC")  
    return 704;
  else if( aa == "RSQ")  
    return 705;
  else if( aa == "RSE")  
    return 706;
  else if( aa == "RSG")  
    return 707;
  else if( aa == "RSH")  
    return 708;
  else if( aa == "RSI")  
    return 709;
  else if( aa == "RSL")  
    return 710;
  else if( aa == "RSK")  
    return 711;
  else if( aa == "RSM")  
    return 712;
  else if( aa == "RSF")  
    return 713;
  else if( aa == "RSP")  
    return 714;
  else if( aa == "RSS")  
    return 715;
  else if( aa == "RST")  
    return 716;
  else if( aa == "RSW")  
    return 717;
  else if( aa == "RSY")  
    return 718;
  else if( aa == "RSV")  
    return 719;
  else if( aa == "RTA")  
    return 720;
  else if( aa == "RTR")  
    return 721;
  else if( aa == "RTN")  
    return 722;
  else if( aa == "RTD")  
    return 723;
  else if( aa == "RTC")  
    return 724;
  else if( aa == "RTQ")  
    return 725;
  else if( aa == "RTE")  
    return 726;
  else if( aa == "RTG")  
    return 727;
  else if( aa == "RTH")  
    return 728;
  else if( aa == "RTI")  
    return 729;
  else if( aa == "RTL")  
    return 730;
  else if( aa == "RTK")  
    return 731;
  else if( aa == "RTM")  
    return 732;
  else if( aa == "RTF")  
    return 733;
  else if( aa == "RTP")  
    return 734;
  else if( aa == "RTS")  
    return 735;
  else if( aa == "RTT")  
    return 736;
  else if( aa == "RTW")  
    return 737;
  else if( aa == "RTY")  
    return 738;
  else if( aa == "RTV")  
    return 739;
  else if( aa == "RWA")  
    return 740;
  else if( aa == "RWR")  
    return 741;
  else if( aa == "RWN")  
    return 742;
  else if( aa == "RWD")  
    return 743;
  else if( aa == "RWC")  
    return 744;
  else if( aa == "RWQ")  
    return 745;
  else if( aa == "RWE")  
    return 746;
  else if( aa == "RWG")  
    return 747;
  else if( aa == "RWH")  
    return 748;
  else if( aa == "RWI")  
    return 749;
  else if( aa == "RWL")  
    return 750;
  else if( aa == "RWK")  
    return 751;
  else if( aa == "RWM")  
    return 752;
  else if( aa == "RWF")  
    return 753;
  else if( aa == "RWP")  
    return 754;
  else if( aa == "RWS")  
    return 755;
  else if( aa == "RWT")  
    return 756;
  else if( aa == "RWW")  
    return 757;
  else if( aa == "RWY")  
    return 758;
  else if( aa == "RWV")  
    return 759;
  else if( aa == "RYA")  
    return 760;
  else if( aa == "RYR")  
    return 761;
  else if( aa == "RYN")  
    return 762;
  else if( aa == "RYD")  
    return 763;
  else if( aa == "RYC")  
    return 764;
  else if( aa == "RYQ")  
    return 765;
  else if( aa == "RYE")  
    return 766;
  else if( aa == "RYG")  
    return 767;
  else if( aa == "RYH")  
    return 768;
  else if( aa == "RYI")  
    return 769;
  else if( aa == "RYL")  
    return 770;
  else if( aa == "RYK")  
    return 771;
  else if( aa == "RYM")  
    return 772;
  else if( aa == "RYF")  
    return 773;
  else if( aa == "RYP")  
    return 774;
  else if( aa == "RYS")  
    return 775;
  else if( aa == "RYT")  
    return 776;
  else if( aa == "RYW")  
    return 777;
  else if( aa == "RYY")  
    return 778;
  else if( aa == "RYV")  
    return 779;
  else if( aa == "RVA")  
    return 780;
  else if( aa == "RVR")  
    return 781;
  else if( aa == "RVN")  
    return 782;
  else if( aa == "RVD")  
    return 783;
  else if( aa == "RVC")  
    return 784;
  else if( aa == "RVQ")  
    return 785;
  else if( aa == "RVE")  
    return 786;
  else if( aa == "RVG")  
    return 787;
  else if( aa == "RVH")  
    return 788;
  else if( aa == "RVI")  
    return 789;
  else if( aa == "RVL")  
    return 790;
  else if( aa == "RVK")  
    return 791;
  else if( aa == "RVM")  
    return 792;
  else if( aa == "RVF")  
    return 793;
  else if( aa == "RVP")  
    return 794;
  else if( aa == "RVS")  
    return 795;
  else if( aa == "RVT")  
    return 796;
  else if( aa == "RVW")  
    return 797;
  else if( aa == "RVY")  
    return 798;
  else if( aa == "RVV")  
    return 799;
  else if( aa == "NAA")  
    return 800;
  else if( aa == "NAR")  
    return 801;
  else if( aa == "NAN")  
    return 802;
  else if( aa == "NAD")  
    return 803;
  else if( aa == "NAC")  
    return 804;
  else if( aa == "NAQ")  
    return 805;
  else if( aa == "NAE")  
    return 806;
  else if( aa == "NAG")  
    return 807;
  else if( aa == "NAH")  
    return 808;
  else if( aa == "NAI")  
    return 809;
  else if( aa == "NAL")  
    return 810;
  else if( aa == "NAK")  
    return 811;
  else if( aa == "NAM")  
    return 812;
  else if( aa == "NAF")  
    return 813;
  else if( aa == "NAP")  
    return 814;
  else if( aa == "NAS")  
    return 815;
  else if( aa == "NAT")  
    return 816;
  else if( aa == "NAW")  
    return 817;
  else if( aa == "NAY")  
    return 818;
  else if( aa == "NAV")  
    return 819;
  else if( aa == "NRA")  
    return 820;
  else if( aa == "NRR")  
    return 821;
  else if( aa == "NRN")  
    return 822;
  else if( aa == "NRD")  
    return 823;
  else if( aa == "NRC")  
    return 824;
  else if( aa == "NRQ")  
    return 825;
  else if( aa == "NRE")  
    return 826;
  else if( aa == "NRG")  
    return 827;
  else if( aa == "NRH")  
    return 828;
  else if( aa == "NRI")  
    return 829;
  else if( aa == "NRL")  
    return 830;
  else if( aa == "NRK")  
    return 831;
  else if( aa == "NRM")  
    return 832;
  else if( aa == "NRF")  
    return 833;
  else if( aa == "NRP")  
    return 834;
  else if( aa == "NRS")  
    return 835;
  else if( aa == "NRT")  
    return 836;
  else if( aa == "NRW")  
    return 837;
  else if( aa == "NRY")  
    return 838;
  else if( aa == "NRV")  
    return 839;
  else if( aa == "NNA")  
    return 840;
  else if( aa == "NNR")  
    return 841;
  else if( aa == "NNN")  
    return 842;
  else if( aa == "NND")  
    return 843;
  else if( aa == "NNC")  
    return 844;
  else if( aa == "NNQ")  
    return 845;
  else if( aa == "NNE")  
    return 846;
  else if( aa == "NNG")  
    return 847;
  else if( aa == "NNH")  
    return 848;
  else if( aa == "NNI")  
    return 849;
  else if( aa == "NNL")  
    return 850;
  else if( aa == "NNK")  
    return 851;
  else if( aa == "NNM")  
    return 852;
  else if( aa == "NNF")  
    return 853;
  else if( aa == "NNP")  
    return 854;
  else if( aa == "NNS")  
    return 855;
  else if( aa == "NNT")  
    return 856;
  else if( aa == "NNW")  
    return 857;
  else if( aa == "NNY")  
    return 858;
  else if( aa == "NNV")  
    return 859;
  else if( aa == "NDA")  
    return 860;
  else if( aa == "NDR")  
    return 861;
  else if( aa == "NDN")  
    return 862;
  else if( aa == "NDD")  
    return 863;
  else if( aa == "NDC")  
    return 864;
  else if( aa == "NDQ")  
    return 865;
  else if( aa == "NDE")  
    return 866;
  else if( aa == "NDG")  
    return 867;
  else if( aa == "NDH")  
    return 868;
  else if( aa == "NDI")  
    return 869;
  else if( aa == "NDL")  
    return 870;
  else if( aa == "NDK")  
    return 871;
  else if( aa == "NDM")  
    return 872;
  else if( aa == "NDF")  
    return 873;
  else if( aa == "NDP")  
    return 874;
  else if( aa == "NDS")  
    return 875;
  else if( aa == "NDT")  
    return 876;
  else if( aa == "NDW")  
    return 877;
  else if( aa == "NDY")  
    return 878;
  else if( aa == "NDV")  
    return 879;
  else if( aa == "NCA")  
    return 880;
  else if( aa == "NCR")  
    return 881;
  else if( aa == "NCN")  
    return 882;
  else if( aa == "NCD")  
    return 883;
  else if( aa == "NCC")  
    return 884;
  else if( aa == "NCQ")  
    return 885;
  else if( aa == "NCE")  
    return 886;
  else if( aa == "NCG")  
    return 887;
  else if( aa == "NCH")  
    return 888;
  else if( aa == "NCI")  
    return 889;
  else if( aa == "NCL")  
    return 890;
  else if( aa == "NCK")  
    return 891;
  else if( aa == "NCM")  
    return 892;
  else if( aa == "NCF")  
    return 893;
  else if( aa == "NCP")  
    return 894;
  else if( aa == "NCS")  
    return 895;
  else if( aa == "NCT")  
    return 896;
  else if( aa == "NCW")  
    return 897;
  else if( aa == "NCY")  
    return 898;
  else if( aa == "NCV")  
    return 899;
  else if( aa == "NQA")  
    return 900;
  else if( aa == "NQR")  
    return 901;
  else if( aa == "NQN")  
    return 902;
  else if( aa == "NQD")  
    return 903;
  else if( aa == "NQC")  
    return 904;
  else if( aa == "NQQ")  
    return 905;
  else if( aa == "NQE")  
    return 906;
  else if( aa == "NQG")  
    return 907;
  else if( aa == "NQH")  
    return 908;
  else if( aa == "NQI")  
    return 909;
  else if( aa == "NQL")  
    return 910;
  else if( aa == "NQK")  
    return 911;
  else if( aa == "NQM")  
    return 912;
  else if( aa == "NQF")  
    return 913;
  else if( aa == "NQP")  
    return 914;
  else if( aa == "NQS")  
    return 915;
  else if( aa == "NQT")  
    return 916;
  else if( aa == "NQW")  
    return 917;
  else if( aa == "NQY")  
    return 918;
  else if( aa == "NQV")  
    return 919;
  else if( aa == "NEA")  
    return 920;
  else if( aa == "NER")  
    return 921;
  else if( aa == "NEN")  
    return 922;
  else if( aa == "NED")  
    return 923;
  else if( aa == "NEC")  
    return 924;
  else if( aa == "NEQ")  
    return 925;
  else if( aa == "NEE")  
    return 926;
  else if( aa == "NEG")  
    return 927;
  else if( aa == "NEH")  
    return 928;
  else if( aa == "NEI")  
    return 929;
  else if( aa == "NEL")  
    return 930;
  else if( aa == "NEK")  
    return 931;
  else if( aa == "NEM")  
    return 932;
  else if( aa == "NEF")  
    return 933;
  else if( aa == "NEP")  
    return 934;
  else if( aa == "NES")  
    return 935;
  else if( aa == "NET")  
    return 936;
  else if( aa == "NEW")  
    return 937;
  else if( aa == "NEY")  
    return 938;
  else if( aa == "NEV")  
    return 939;
  else if( aa == "NGA")  
    return 940;
  else if( aa == "NGR")  
    return 941;
  else if( aa == "NGN")  
    return 942;
  else if( aa == "NGD")  
    return 943;
  else if( aa == "NGC")  
    return 944;
  else if( aa == "NGQ")  
    return 945;
  else if( aa == "NGE")  
    return 946;
  else if( aa == "NGG")  
    return 947;
  else if( aa == "NGH")  
    return 948;
  else if( aa == "NGI")  
    return 949;
  else if( aa == "NGL")  
    return 950;
  else if( aa == "NGK")  
    return 951;
  else if( aa == "NGM")  
    return 952;
  else if( aa == "NGF")  
    return 953;
  else if( aa == "NGP")  
    return 954;
  else if( aa == "NGS")  
    return 955;
  else if( aa == "NGT")  
    return 956;
  else if( aa == "NGW")  
    return 957;
  else if( aa == "NGY")  
    return 958;
  else if( aa == "NGV")  
    return 959;
  else if( aa == "NHA")  
    return 960;
  else if( aa == "NHR")  
    return 961;
  else if( aa == "NHN")  
    return 962;
  else if( aa == "NHD")  
    return 963;
  else if( aa == "NHC")  
    return 964;
  else if( aa == "NHQ")  
    return 965;
  else if( aa == "NHE")  
    return 966;
  else if( aa == "NHG")  
    return 967;
  else if( aa == "NHH")  
    return 968;
  else if( aa == "NHI")  
    return 969;
  else if( aa == "NHL")  
    return 970;
  else if( aa == "NHK")  
    return 971;
  else if( aa == "NHM")  
    return 972;
  else if( aa == "NHF")  
    return 973;
  else if( aa == "NHP")  
    return 974;
  else if( aa == "NHS")  
    return 975;
  else if( aa == "NHT")  
    return 976;
  else if( aa == "NHW")  
    return 977;
  else if( aa == "NHY")  
    return 978;
  else if( aa == "NHV")  
    return 979;
  else if( aa == "NIA")  
    return 980;
  else if( aa == "NIR")  
    return 981;
  else if( aa == "NIN")  
    return 982;
  else if( aa == "NID")  
    return 983;
  else if( aa == "NIC")  
    return 984;
  else if( aa == "NIQ")  
    return 985;
  else if( aa == "NIE")  
    return 986;
  else if( aa == "NIG")  
    return 987;
  else if( aa == "NIH")  
    return 988;
  else if( aa == "NII")  
    return 989;
  else if( aa == "NIL")  
    return 990;
  else if( aa == "NIK")  
    return 991;
  else if( aa == "NIM")  
    return 992;
  else if( aa == "NIF")  
    return 993;
  else if( aa == "NIP")  
    return 994;
  else if( aa == "NIS")  
    return 995;
  else if( aa == "NIT")  
    return 996;
  else if( aa == "NIW")  
    return 997;
  else if( aa == "NIY")  
    return 998;
  else if( aa == "NIV")  
    return 999;
  else if( aa == "NLA")  
    return 1000;
  else if( aa == "NLR")  
    return 1001;
  else if( aa == "NLN")  
    return 1002;
  else if( aa == "NLD")  
    return 1003;
  else if( aa == "NLC")  
    return 1004;
  else if( aa == "NLQ")  
    return 1005;
  else if( aa == "NLE")  
    return 1006;
  else if( aa == "NLG")  
    return 1007;
  else if( aa == "NLH")  
    return 1008;
  else if( aa == "NLI")  
    return 1009;
  else if( aa == "NLL")  
    return 1010;
  else if( aa == "NLK")  
    return 1011;
  else if( aa == "NLM")  
    return 1012;
  else if( aa == "NLF")  
    return 1013;
  else if( aa == "NLP")  
    return 1014;
  else if( aa == "NLS")  
    return 1015;
  else if( aa == "NLT")  
    return 1016;
  else if( aa == "NLW")  
    return 1017;
  else if( aa == "NLY")  
    return 1018;
  else if( aa == "NLV")  
    return 1019;
  else if( aa == "NKA")  
    return 1020;
  else if( aa == "NKR")  
    return 1021;
  else if( aa == "NKN")  
    return 1022;
  else if( aa == "NKD")  
    return 1023;
  else if( aa == "NKC")  
    return 1024;
  else if( aa == "NKQ")  
    return 1025;
  else if( aa == "NKE")  
    return 1026;
  else if( aa == "NKG")  
    return 1027;
  else if( aa == "NKH")  
    return 1028;
  else if( aa == "NKI")  
    return 1029;
  else if( aa == "NKL")  
    return 1030;
  else if( aa == "NKK")  
    return 1031;
  else if( aa == "NKM")  
    return 1032;
  else if( aa == "NKF")  
    return 1033;
  else if( aa == "NKP")  
    return 1034;
  else if( aa == "NKS")  
    return 1035;
  else if( aa == "NKT")  
    return 1036;
  else if( aa == "NKW")  
    return 1037;
  else if( aa == "NKY")  
    return 1038;
  else if( aa == "NKV")  
    return 1039;
  else if( aa == "NMA")  
    return 1040;
  else if( aa == "NMR")  
    return 1041;
  else if( aa == "NMN")  
    return 1042;
  else if( aa == "NMD")  
    return 1043;
  else if( aa == "NMC")  
    return 1044;
  else if( aa == "NMQ")  
    return 1045;
  else if( aa == "NME")  
    return 1046;
  else if( aa == "NMG")  
    return 1047;
  else if( aa == "NMH")  
    return 1048;
  else if( aa == "NMI")  
    return 1049;
  else if( aa == "NML")  
    return 1050;
  else if( aa == "NMK")  
    return 1051;
  else if( aa == "NMM")  
    return 1052;
  else if( aa == "NMF")  
    return 1053;
  else if( aa == "NMP")  
    return 1054;
  else if( aa == "NMS")  
    return 1055;
  else if( aa == "NMT")  
    return 1056;
  else if( aa == "NMW")  
    return 1057;
  else if( aa == "NMY")  
    return 1058;
  else if( aa == "NMV")  
    return 1059;
  else if( aa == "NFA")  
    return 1060;
  else if( aa == "NFR")  
    return 1061;
  else if( aa == "NFN")  
    return 1062;
  else if( aa == "NFD")  
    return 1063;
  else if( aa == "NFC")  
    return 1064;
  else if( aa == "NFQ")  
    return 1065;
  else if( aa == "NFE")  
    return 1066;
  else if( aa == "NFG")  
    return 1067;
  else if( aa == "NFH")  
    return 1068;
  else if( aa == "NFI")  
    return 1069;
  else if( aa == "NFL")  
    return 1070;
  else if( aa == "NFK")  
    return 1071;
  else if( aa == "NFM")  
    return 1072;
  else if( aa == "NFF")  
    return 1073;
  else if( aa == "NFP")  
    return 1074;
  else if( aa == "NFS")  
    return 1075;
  else if( aa == "NFT")  
    return 1076;
  else if( aa == "NFW")  
    return 1077;
  else if( aa == "NFY")  
    return 1078;
  else if( aa == "NFV")  
    return 1079;
  else if( aa == "NPA")  
    return 1080;
  else if( aa == "NPR")  
    return 1081;
  else if( aa == "NPN")  
    return 1082;
  else if( aa == "NPD")  
    return 1083;
  else if( aa == "NPC")  
    return 1084;
  else if( aa == "NPQ")  
    return 1085;
  else if( aa == "NPE")  
    return 1086;
  else if( aa == "NPG")  
    return 1087;
  else if( aa == "NPH")  
    return 1088;
  else if( aa == "NPI")  
    return 1089;
  else if( aa == "NPL")  
    return 1090;
  else if( aa == "NPK")  
    return 1091;
  else if( aa == "NPM")  
    return 1092;
  else if( aa == "NPF")  
    return 1093;
  else if( aa == "NPP")  
    return 1094;
  else if( aa == "NPS")  
    return 1095;
  else if( aa == "NPT")  
    return 1096;
  else if( aa == "NPW")  
    return 1097;
  else if( aa == "NPY")  
    return 1098;
  else if( aa == "NPV")  
    return 1099;
  else if( aa == "NSA")  
    return 1100;
  else if( aa == "NSR")  
    return 1101;
  else if( aa == "NSN")  
    return 1102;
  else if( aa == "NSD")  
    return 1103;
  else if( aa == "NSC")  
    return 1104;
  else if( aa == "NSQ")  
    return 1105;
  else if( aa == "NSE")  
    return 1106;
  else if( aa == "NSG")  
    return 1107;
  else if( aa == "NSH")  
    return 1108;
  else if( aa == "NSI")  
    return 1109;
  else if( aa == "NSL")  
    return 1110;
  else if( aa == "NSK")  
    return 1111;
  else if( aa == "NSM")  
    return 1112;
  else if( aa == "NSF")  
    return 1113;
  else if( aa == "NSP")  
    return 1114;
  else if( aa == "NSS")  
    return 1115;
  else if( aa == "NST")  
    return 1116;
  else if( aa == "NSW")  
    return 1117;
  else if( aa == "NSY")  
    return 1118;
  else if( aa == "NSV")  
    return 1119;
  else if( aa == "NTA")  
    return 1120;
  else if( aa == "NTR")  
    return 1121;
  else if( aa == "NTN")  
    return 1122;
  else if( aa == "NTD")  
    return 1123;
  else if( aa == "NTC")  
    return 1124;
  else if( aa == "NTQ")  
    return 1125;
  else if( aa == "NTE")  
    return 1126;
  else if( aa == "NTG")  
    return 1127;
  else if( aa == "NTH")  
    return 1128;
  else if( aa == "NTI")  
    return 1129;
  else if( aa == "NTL")  
    return 1130;
  else if( aa == "NTK")  
    return 1131;
  else if( aa == "NTM")  
    return 1132;
  else if( aa == "NTF")  
    return 1133;
  else if( aa == "NTP")  
    return 1134;
  else if( aa == "NTS")  
    return 1135;
  else if( aa == "NTT")  
    return 1136;
  else if( aa == "NTW")  
    return 1137;
  else if( aa == "NTY")  
    return 1138;
  else if( aa == "NTV")  
    return 1139;
  else if( aa == "NWA")  
    return 1140;
  else if( aa == "NWR")  
    return 1141;
  else if( aa == "NWN")  
    return 1142;
  else if( aa == "NWD")  
    return 1143;
  else if( aa == "NWC")  
    return 1144;
  else if( aa == "NWQ")  
    return 1145;
  else if( aa == "NWE")  
    return 1146;
  else if( aa == "NWG")  
    return 1147;
  else if( aa == "NWH")  
    return 1148;
  else if( aa == "NWI")  
    return 1149;
  else if( aa == "NWL")  
    return 1150;
  else if( aa == "NWK")  
    return 1151;
  else if( aa == "NWM")  
    return 1152;
  else if( aa == "NWF")  
    return 1153;
  else if( aa == "NWP")  
    return 1154;
  else if( aa == "NWS")  
    return 1155;
  else if( aa == "NWT")  
    return 1156;
  else if( aa == "NWW")  
    return 1157;
  else if( aa == "NWY")  
    return 1158;
  else if( aa == "NWV")  
    return 1159;
  else if( aa == "NYA")  
    return 1160;
  else if( aa == "NYR")  
    return 1161;
  else if( aa == "NYN")  
    return 1162;
  else if( aa == "NYD")  
    return 1163;
  else if( aa == "NYC")  
    return 1164;
  else if( aa == "NYQ")  
    return 1165;
  else if( aa == "NYE")  
    return 1166;
  else if( aa == "NYG")  
    return 1167;
  else if( aa == "NYH")  
    return 1168;
  else if( aa == "NYI")  
    return 1169;
  else if( aa == "NYL")  
    return 1170;
  else if( aa == "NYK")  
    return 1171;
  else if( aa == "NYM")  
    return 1172;
  else if( aa == "NYF")  
    return 1173;
  else if( aa == "NYP")  
    return 1174;
  else if( aa == "NYS")  
    return 1175;
  else if( aa == "NYT")  
    return 1176;
  else if( aa == "NYW")  
    return 1177;
  else if( aa == "NYY")  
    return 1178;
  else if( aa == "NYV")  
    return 1179;
  else if( aa == "NVA")  
    return 1180;
  else if( aa == "NVR")  
    return 1181;
  else if( aa == "NVN")  
    return 1182;
  else if( aa == "NVD")  
    return 1183;
  else if( aa == "NVC")  
    return 1184;
  else if( aa == "NVQ")  
    return 1185;
  else if( aa == "NVE")  
    return 1186;
  else if( aa == "NVG")  
    return 1187;
  else if( aa == "NVH")  
    return 1188;
  else if( aa == "NVI")  
    return 1189;
  else if( aa == "NVL")  
    return 1190;
  else if( aa == "NVK")  
    return 1191;
  else if( aa == "NVM")  
    return 1192;
  else if( aa == "NVF")  
    return 1193;
  else if( aa == "NVP")  
    return 1194;
  else if( aa == "NVS")  
    return 1195;
  else if( aa == "NVT")  
    return 1196;
  else if( aa == "NVW")  
    return 1197;
  else if( aa == "NVY")  
    return 1198;
  else if( aa == "NVV")  
    return 1199;
  else if( aa == "DAA")  
    return 1200;
  else if( aa == "DAR")  
    return 1201;
  else if( aa == "DAN")  
    return 1202;
  else if( aa == "DAD")  
    return 1203;
  else if( aa == "DAC")  
    return 1204;
  else if( aa == "DAQ")  
    return 1205;
  else if( aa == "DAE")  
    return 1206;
  else if( aa == "DAG")  
    return 1207;
  else if( aa == "DAH")  
    return 1208;
  else if( aa == "DAI")  
    return 1209;
  else if( aa == "DAL")  
    return 1210;
  else if( aa == "DAK")  
    return 1211;
  else if( aa == "DAM")  
    return 1212;
  else if( aa == "DAF")  
    return 1213;
  else if( aa == "DAP")  
    return 1214;
  else if( aa == "DAS")  
    return 1215;
  else if( aa == "DAT")  
    return 1216;
  else if( aa == "DAW")  
    return 1217;
  else if( aa == "DAY")  
    return 1218;
  else if( aa == "DAV")  
    return 1219;
  else if( aa == "DRA")  
    return 1220;
  else if( aa == "DRR")  
    return 1221;
  else if( aa == "DRN")  
    return 1222;
  else if( aa == "DRD")  
    return 1223;
  else if( aa == "DRC")  
    return 1224;
  else if( aa == "DRQ")  
    return 1225;
  else if( aa == "DRE")  
    return 1226;
  else if( aa == "DRG")  
    return 1227;
  else if( aa == "DRH")  
    return 1228;
  else if( aa == "DRI")  
    return 1229;
  else if( aa == "DRL")  
    return 1230;
  else if( aa == "DRK")  
    return 1231;
  else if( aa == "DRM")  
    return 1232;
  else if( aa == "DRF")  
    return 1233;
  else if( aa == "DRP")  
    return 1234;
  else if( aa == "DRS")  
    return 1235;
  else if( aa == "DRT")  
    return 1236;
  else if( aa == "DRW")  
    return 1237;
  else if( aa == "DRY")  
    return 1238;
  else if( aa == "DRV")  
    return 1239;
  else if( aa == "DNA")  
    return 1240;
  else if( aa == "DNR")  
    return 1241;
  else if( aa == "DNN")  
    return 1242;
  else if( aa == "DND")  
    return 1243;
  else if( aa == "DNC")  
    return 1244;
  else if( aa == "DNQ")  
    return 1245;
  else if( aa == "DNE")  
    return 1246;
  else if( aa == "DNG")  
    return 1247;
  else if( aa == "DNH")  
    return 1248;
  else if( aa == "DNI")  
    return 1249;
  else if( aa == "DNL")  
    return 1250;
  else if( aa == "DNK")  
    return 1251;
  else if( aa == "DNM")  
    return 1252;
  else if( aa == "DNF")  
    return 1253;
  else if( aa == "DNP")  
    return 1254;
  else if( aa == "DNS")  
    return 1255;
  else if( aa == "DNT")  
    return 1256;
  else if( aa == "DNW")  
    return 1257;
  else if( aa == "DNY")  
    return 1258;
  else if( aa == "DNV")  
    return 1259;
  else if( aa == "DDA")  
    return 1260;
  else if( aa == "DDR")  
    return 1261;
  else if( aa == "DDN")  
    return 1262;
  else if( aa == "DDD")  
    return 1263;
  else if( aa == "DDC")  
    return 1264;
  else if( aa == "DDQ")  
    return 1265;
  else if( aa == "DDE")  
    return 1266;
  else if( aa == "DDG")  
    return 1267;
  else if( aa == "DDH")  
    return 1268;
  else if( aa == "DDI")  
    return 1269;
  else if( aa == "DDL")  
    return 1270;
  else if( aa == "DDK")  
    return 1271;
  else if( aa == "DDM")  
    return 1272;
  else if( aa == "DDF")  
    return 1273;
  else if( aa == "DDP")  
    return 1274;
  else if( aa == "DDS")  
    return 1275;
  else if( aa == "DDT")  
    return 1276;
  else if( aa == "DDW")  
    return 1277;
  else if( aa == "DDY")  
    return 1278;
  else if( aa == "DDV")  
    return 1279;
  else if( aa == "DCA")  
    return 1280;
  else if( aa == "DCR")  
    return 1281;
  else if( aa == "DCN")  
    return 1282;
  else if( aa == "DCD")  
    return 1283;
  else if( aa == "DCC")  
    return 1284;
  else if( aa == "DCQ")  
    return 1285;
  else if( aa == "DCE")  
    return 1286;
  else if( aa == "DCG")  
    return 1287;
  else if( aa == "DCH")  
    return 1288;
  else if( aa == "DCI")  
    return 1289;
  else if( aa == "DCL")  
    return 1290;
  else if( aa == "DCK")  
    return 1291;
  else if( aa == "DCM")  
    return 1292;
  else if( aa == "DCF")  
    return 1293;
  else if( aa == "DCP")  
    return 1294;
  else if( aa == "DCS")  
    return 1295;
  else if( aa == "DCT")  
    return 1296;
  else if( aa == "DCW")  
    return 1297;
  else if( aa == "DCY")  
    return 1298;
  else if( aa == "DCV")  
    return 1299;
  else if( aa == "DQA")  
    return 1300;
  else if( aa == "DQR")  
    return 1301;
  else if( aa == "DQN")  
    return 1302;
  else if( aa == "DQD")  
    return 1303;
  else if( aa == "DQC")  
    return 1304;
  else if( aa == "DQQ")  
    return 1305;
  else if( aa == "DQE")  
    return 1306;
  else if( aa == "DQG")  
    return 1307;
  else if( aa == "DQH")  
    return 1308;
  else if( aa == "DQI")  
    return 1309;
  else if( aa == "DQL")  
    return 1310;
  else if( aa == "DQK")  
    return 1311;
  else if( aa == "DQM")  
    return 1312;
  else if( aa == "DQF")  
    return 1313;
  else if( aa == "DQP")  
    return 1314;
  else if( aa == "DQS")  
    return 1315;
  else if( aa == "DQT")  
    return 1316;
  else if( aa == "DQW")  
    return 1317;
  else if( aa == "DQY")  
    return 1318;
  else if( aa == "DQV")  
    return 1319;
  else if( aa == "DEA")  
    return 1320;
  else if( aa == "DER")  
    return 1321;
  else if( aa == "DEN")  
    return 1322;
  else if( aa == "DED")  
    return 1323;
  else if( aa == "DEC")  
    return 1324;
  else if( aa == "DEQ")  
    return 1325;
  else if( aa == "DEE")  
    return 1326;
  else if( aa == "DEG")  
    return 1327;
  else if( aa == "DEH")  
    return 1328;
  else if( aa == "DEI")  
    return 1329;
  else if( aa == "DEL")  
    return 1330;
  else if( aa == "DEK")  
    return 1331;
  else if( aa == "DEM")  
    return 1332;
  else if( aa == "DEF")  
    return 1333;
  else if( aa == "DEP")  
    return 1334;
  else if( aa == "DES")  
    return 1335;
  else if( aa == "DET")  
    return 1336;
  else if( aa == "DEW")  
    return 1337;
  else if( aa == "DEY")  
    return 1338;
  else if( aa == "DEV")  
    return 1339;
  else if( aa == "DGA")  
    return 1340;
  else if( aa == "DGR")  
    return 1341;
  else if( aa == "DGN")  
    return 1342;
  else if( aa == "DGD")  
    return 1343;
  else if( aa == "DGC")  
    return 1344;
  else if( aa == "DGQ")  
    return 1345;
  else if( aa == "DGE")  
    return 1346;
  else if( aa == "DGG")  
    return 1347;
  else if( aa == "DGH")  
    return 1348;
  else if( aa == "DGI")  
    return 1349;
  else if( aa == "DGL")  
    return 1350;
  else if( aa == "DGK")  
    return 1351;
  else if( aa == "DGM")  
    return 1352;
  else if( aa == "DGF")  
    return 1353;
  else if( aa == "DGP")  
    return 1354;
  else if( aa == "DGS")  
    return 1355;
  else if( aa == "DGT")  
    return 1356;
  else if( aa == "DGW")  
    return 1357;
  else if( aa == "DGY")  
    return 1358;
  else if( aa == "DGV")  
    return 1359;
  else if( aa == "DHA")  
    return 1360;
  else if( aa == "DHR")  
    return 1361;
  else if( aa == "DHN")  
    return 1362;
  else if( aa == "DHD")  
    return 1363;
  else if( aa == "DHC")  
    return 1364;
  else if( aa == "DHQ")  
    return 1365;
  else if( aa == "DHE")  
    return 1366;
  else if( aa == "DHG")  
    return 1367;
  else if( aa == "DHH")  
    return 1368;
  else if( aa == "DHI")  
    return 1369;
  else if( aa == "DHL")  
    return 1370;
  else if( aa == "DHK")  
    return 1371;
  else if( aa == "DHM")  
    return 1372;
  else if( aa == "DHF")  
    return 1373;
  else if( aa == "DHP")  
    return 1374;
  else if( aa == "DHS")  
    return 1375;
  else if( aa == "DHT")  
    return 1376;
  else if( aa == "DHW")  
    return 1377;
  else if( aa == "DHY")  
    return 1378;
  else if( aa == "DHV")  
    return 1379;
  else if( aa == "DIA")  
    return 1380;
  else if( aa == "DIR")  
    return 1381;
  else if( aa == "DIN")  
    return 1382;
  else if( aa == "DID")  
    return 1383;
  else if( aa == "DIC")  
    return 1384;
  else if( aa == "DIQ")  
    return 1385;
  else if( aa == "DIE")  
    return 1386;
  else if( aa == "DIG")  
    return 1387;
  else if( aa == "DIH")  
    return 1388;
  else if( aa == "DII")  
    return 1389;
  else if( aa == "DIL")  
    return 1390;
  else if( aa == "DIK")  
    return 1391;
  else if( aa == "DIM")  
    return 1392;
  else if( aa == "DIF")  
    return 1393;
  else if( aa == "DIP")  
    return 1394;
  else if( aa == "DIS")  
    return 1395;
  else if( aa == "DIT")  
    return 1396;
  else if( aa == "DIW")  
    return 1397;
  else if( aa == "DIY")  
    return 1398;
  else if( aa == "DIV")  
    return 1399;
  else if( aa == "DLA")  
    return 1400;
  else if( aa == "DLR")  
    return 1401;
  else if( aa == "DLN")  
    return 1402;
  else if( aa == "DLD")  
    return 1403;
  else if( aa == "DLC")  
    return 1404;
  else if( aa == "DLQ")  
    return 1405;
  else if( aa == "DLE")  
    return 1406;
  else if( aa == "DLG")  
    return 1407;
  else if( aa == "DLH")  
    return 1408;
  else if( aa == "DLI")  
    return 1409;
  else if( aa == "DLL")  
    return 1410;
  else if( aa == "DLK")  
    return 1411;
  else if( aa == "DLM")  
    return 1412;
  else if( aa == "DLF")  
    return 1413;
  else if( aa == "DLP")  
    return 1414;
  else if( aa == "DLS")  
    return 1415;
  else if( aa == "DLT")  
    return 1416;
  else if( aa == "DLW")  
    return 1417;
  else if( aa == "DLY")  
    return 1418;
  else if( aa == "DLV")  
    return 1419;
  else if( aa == "DKA")  
    return 1420;
  else if( aa == "DKR")  
    return 1421;
  else if( aa == "DKN")  
    return 1422;
  else if( aa == "DKD")  
    return 1423;
  else if( aa == "DKC")  
    return 1424;
  else if( aa == "DKQ")  
    return 1425;
  else if( aa == "DKE")  
    return 1426;
  else if( aa == "DKG")  
    return 1427;
  else if( aa == "DKH")  
    return 1428;
  else if( aa == "DKI")  
    return 1429;
  else if( aa == "DKL")  
    return 1430;
  else if( aa == "DKK")  
    return 1431;
  else if( aa == "DKM")  
    return 1432;
  else if( aa == "DKF")  
    return 1433;
  else if( aa == "DKP")  
    return 1434;
  else if( aa == "DKS")  
    return 1435;
  else if( aa == "DKT")  
    return 1436;
  else if( aa == "DKW")  
    return 1437;
  else if( aa == "DKY")  
    return 1438;
  else if( aa == "DKV")  
    return 1439;
  else if( aa == "DMA")  
    return 1440;
  else if( aa == "DMR")  
    return 1441;
  else if( aa == "DMN")  
    return 1442;
  else if( aa == "DMD")  
    return 1443;
  else if( aa == "DMC")  
    return 1444;
  else if( aa == "DMQ")  
    return 1445;
  else if( aa == "DME")  
    return 1446;
  else if( aa == "DMG")  
    return 1447;
  else if( aa == "DMH")  
    return 1448;
  else if( aa == "DMI")  
    return 1449;
  else if( aa == "DML")  
    return 1450;
  else if( aa == "DMK")  
    return 1451;
  else if( aa == "DMM")  
    return 1452;
  else if( aa == "DMF")  
    return 1453;
  else if( aa == "DMP")  
    return 1454;
  else if( aa == "DMS")  
    return 1455;
  else if( aa == "DMT")  
    return 1456;
  else if( aa == "DMW")  
    return 1457;
  else if( aa == "DMY")  
    return 1458;
  else if( aa == "DMV")  
    return 1459;
  else if( aa == "DFA")  
    return 1460;
  else if( aa == "DFR")  
    return 1461;
  else if( aa == "DFN")  
    return 1462;
  else if( aa == "DFD")  
    return 1463;
  else if( aa == "DFC")  
    return 1464;
  else if( aa == "DFQ")  
    return 1465;
  else if( aa == "DFE")  
    return 1466;
  else if( aa == "DFG")  
    return 1467;
  else if( aa == "DFH")  
    return 1468;
  else if( aa == "DFI")  
    return 1469;
  else if( aa == "DFL")  
    return 1470;
  else if( aa == "DFK")  
    return 1471;
  else if( aa == "DFM")  
    return 1472;
  else if( aa == "DFF")  
    return 1473;
  else if( aa == "DFP")  
    return 1474;
  else if( aa == "DFS")  
    return 1475;
  else if( aa == "DFT")  
    return 1476;
  else if( aa == "DFW")  
    return 1477;
  else if( aa == "DFY")  
    return 1478;
  else if( aa == "DFV")  
    return 1479;
  else if( aa == "DPA")  
    return 1480;
  else if( aa == "DPR")  
    return 1481;
  else if( aa == "DPN")  
    return 1482;
  else if( aa == "DPD")  
    return 1483;
  else if( aa == "DPC")  
    return 1484;
  else if( aa == "DPQ")  
    return 1485;
  else if( aa == "DPE")  
    return 1486;
  else if( aa == "DPG")  
    return 1487;
  else if( aa == "DPH")  
    return 1488;
  else if( aa == "DPI")  
    return 1489;
  else if( aa == "DPL")  
    return 1490;
  else if( aa == "DPK")  
    return 1491;
  else if( aa == "DPM")  
    return 1492;
  else if( aa == "DPF")  
    return 1493;
  else if( aa == "DPP")  
    return 1494;
  else if( aa == "DPS")  
    return 1495;
  else if( aa == "DPT")  
    return 1496;
  else if( aa == "DPW")  
    return 1497;
  else if( aa == "DPY")  
    return 1498;
  else if( aa == "DPV")  
    return 1499;
  else if( aa == "DSA")  
    return 1500;
  else if( aa == "DSR")  
    return 1501;
  else if( aa == "DSN")  
    return 1502;
  else if( aa == "DSD")  
    return 1503;
  else if( aa == "DSC")  
    return 1504;
  else if( aa == "DSQ")  
    return 1505;
  else if( aa == "DSE")  
    return 1506;
  else if( aa == "DSG")  
    return 1507;
  else if( aa == "DSH")  
    return 1508;
  else if( aa == "DSI")  
    return 1509;
  else if( aa == "DSL")  
    return 1510;
  else if( aa == "DSK")  
    return 1511;
  else if( aa == "DSM")  
    return 1512;
  else if( aa == "DSF")  
    return 1513;
  else if( aa == "DSP")  
    return 1514;
  else if( aa == "DSS")  
    return 1515;
  else if( aa == "DST")  
    return 1516;
  else if( aa == "DSW")  
    return 1517;
  else if( aa == "DSY")  
    return 1518;
  else if( aa == "DSV")  
    return 1519;
  else if( aa == "DTA")  
    return 1520;
  else if( aa == "DTR")  
    return 1521;
  else if( aa == "DTN")  
    return 1522;
  else if( aa == "DTD")  
    return 1523;
  else if( aa == "DTC")  
    return 1524;
  else if( aa == "DTQ")  
    return 1525;
  else if( aa == "DTE")  
    return 1526;
  else if( aa == "DTG")  
    return 1527;
  else if( aa == "DTH")  
    return 1528;
  else if( aa == "DTI")  
    return 1529;
  else if( aa == "DTL")  
    return 1530;
  else if( aa == "DTK")  
    return 1531;
  else if( aa == "DTM")  
    return 1532;
  else if( aa == "DTF")  
    return 1533;
  else if( aa == "DTP")  
    return 1534;
  else if( aa == "DTS")  
    return 1535;
  else if( aa == "DTT")  
    return 1536;
  else if( aa == "DTW")  
    return 1537;
  else if( aa == "DTY")  
    return 1538;
  else if( aa == "DTV")  
    return 1539;
  else if( aa == "DWA")  
    return 1540;
  else if( aa == "DWR")  
    return 1541;
  else if( aa == "DWN")  
    return 1542;
  else if( aa == "DWD")  
    return 1543;
  else if( aa == "DWC")  
    return 1544;
  else if( aa == "DWQ")  
    return 1545;
  else if( aa == "DWE")  
    return 1546;
  else if( aa == "DWG")  
    return 1547;
  else if( aa == "DWH")  
    return 1548;
  else if( aa == "DWI")  
    return 1549;
  else if( aa == "DWL")  
    return 1550;
  else if( aa == "DWK")  
    return 1551;
  else if( aa == "DWM")  
    return 1552;
  else if( aa == "DWF")  
    return 1553;
  else if( aa == "DWP")  
    return 1554;
  else if( aa == "DWS")  
    return 1555;
  else if( aa == "DWT")  
    return 1556;
  else if( aa == "DWW")  
    return 1557;
  else if( aa == "DWY")  
    return 1558;
  else if( aa == "DWV")  
    return 1559;
  else if( aa == "DYA")  
    return 1560;
  else if( aa == "DYR")  
    return 1561;
  else if( aa == "DYN")  
    return 1562;
  else if( aa == "DYD")  
    return 1563;
  else if( aa == "DYC")  
    return 1564;
  else if( aa == "DYQ")  
    return 1565;
  else if( aa == "DYE")  
    return 1566;
  else if( aa == "DYG")  
    return 1567;
  else if( aa == "DYH")  
    return 1568;
  else if( aa == "DYI")  
    return 1569;
  else if( aa == "DYL")  
    return 1570;
  else if( aa == "DYK")  
    return 1571;
  else if( aa == "DYM")  
    return 1572;
  else if( aa == "DYF")  
    return 1573;
  else if( aa == "DYP")  
    return 1574;
  else if( aa == "DYS")  
    return 1575;
  else if( aa == "DYT")  
    return 1576;
  else if( aa == "DYW")  
    return 1577;
  else if( aa == "DYY")  
    return 1578;
  else if( aa == "DYV")  
    return 1579;
  else if( aa == "DVA")  
    return 1580;
  else if( aa == "DVR")  
    return 1581;
  else if( aa == "DVN")  
    return 1582;
  else if( aa == "DVD")  
    return 1583;
  else if( aa == "DVC")  
    return 1584;
  else if( aa == "DVQ")  
    return 1585;
  else if( aa == "DVE")  
    return 1586;
  else if( aa == "DVG")  
    return 1587;
  else if( aa == "DVH")  
    return 1588;
  else if( aa == "DVI")  
    return 1589;
  else if( aa == "DVL")  
    return 1590;
  else if( aa == "DVK")  
    return 1591;
  else if( aa == "DVM")  
    return 1592;
  else if( aa == "DVF")  
    return 1593;
  else if( aa == "DVP")  
    return 1594;
  else if( aa == "DVS")  
    return 1595;
  else if( aa == "DVT")  
    return 1596;
  else if( aa == "DVW")  
    return 1597;
  else if( aa == "DVY")  
    return 1598;
  else if( aa == "DVV")  
    return 1599;
  else if( aa == "CAA")  
    return 1600;
  else if( aa == "CAR")  
    return 1601;
  else if( aa == "CAN")  
    return 1602;
  else if( aa == "CAD")  
    return 1603;
  else if( aa == "CAC")  
    return 1604;
  else if( aa == "CAQ")  
    return 1605;
  else if( aa == "CAE")  
    return 1606;
  else if( aa == "CAG")  
    return 1607;
  else if( aa == "CAH")  
    return 1608;
  else if( aa == "CAI")  
    return 1609;
  else if( aa == "CAL")  
    return 1610;
  else if( aa == "CAK")  
    return 1611;
  else if( aa == "CAM")  
    return 1612;
  else if( aa == "CAF")  
    return 1613;
  else if( aa == "CAP")  
    return 1614;
  else if( aa == "CAS")  
    return 1615;
  else if( aa == "CAT")  
    return 1616;
  else if( aa == "CAW")  
    return 1617;
  else if( aa == "CAY")  
    return 1618;
  else if( aa == "CAV")  
    return 1619;
  else if( aa == "CRA")  
    return 1620;
  else if( aa == "CRR")  
    return 1621;
  else if( aa == "CRN")  
    return 1622;
  else if( aa == "CRD")  
    return 1623;
  else if( aa == "CRC")  
    return 1624;
  else if( aa == "CRQ")  
    return 1625;
  else if( aa == "CRE")  
    return 1626;
  else if( aa == "CRG")  
    return 1627;
  else if( aa == "CRH")  
    return 1628;
  else if( aa == "CRI")  
    return 1629;
  else if( aa == "CRL")  
    return 1630;
  else if( aa == "CRK")  
    return 1631;
  else if( aa == "CRM")  
    return 1632;
  else if( aa == "CRF")  
    return 1633;
  else if( aa == "CRP")  
    return 1634;
  else if( aa == "CRS")  
    return 1635;
  else if( aa == "CRT")  
    return 1636;
  else if( aa == "CRW")  
    return 1637;
  else if( aa == "CRY")  
    return 1638;
  else if( aa == "CRV")  
    return 1639;
  else if( aa == "CNA")  
    return 1640;
  else if( aa == "CNR")  
    return 1641;
  else if( aa == "CNN")  
    return 1642;
  else if( aa == "CND")  
    return 1643;
  else if( aa == "CNC")  
    return 1644;
  else if( aa == "CNQ")  
    return 1645;
  else if( aa == "CNE")  
    return 1646;
  else if( aa == "CNG")  
    return 1647;
  else if( aa == "CNH")  
    return 1648;
  else if( aa == "CNI")  
    return 1649;
  else if( aa == "CNL")  
    return 1650;
  else if( aa == "CNK")  
    return 1651;
  else if( aa == "CNM")  
    return 1652;
  else if( aa == "CNF")  
    return 1653;
  else if( aa == "CNP")  
    return 1654;
  else if( aa == "CNS")  
    return 1655;
  else if( aa == "CNT")  
    return 1656;
  else if( aa == "CNW")  
    return 1657;
  else if( aa == "CNY")  
    return 1658;
  else if( aa == "CNV")  
    return 1659;
  else if( aa == "CDA")  
    return 1660;
  else if( aa == "CDR")  
    return 1661;
  else if( aa == "CDN")  
    return 1662;
  else if( aa == "CDD")  
    return 1663;
  else if( aa == "CDC")  
    return 1664;
  else if( aa == "CDQ")  
    return 1665;
  else if( aa == "CDE")  
    return 1666;
  else if( aa == "CDG")  
    return 1667;
  else if( aa == "CDH")  
    return 1668;
  else if( aa == "CDI")  
    return 1669;
  else if( aa == "CDL")  
    return 1670;
  else if( aa == "CDK")  
    return 1671;
  else if( aa == "CDM")  
    return 1672;
  else if( aa == "CDF")  
    return 1673;
  else if( aa == "CDP")  
    return 1674;
  else if( aa == "CDS")  
    return 1675;
  else if( aa == "CDT")  
    return 1676;
  else if( aa == "CDW")  
    return 1677;
  else if( aa == "CDY")  
    return 1678;
  else if( aa == "CDV")  
    return 1679;
  else if( aa == "CCA")  
    return 1680;
  else if( aa == "CCR")  
    return 1681;
  else if( aa == "CCN")  
    return 1682;
  else if( aa == "CCD")  
    return 1683;
  else if( aa == "CCC")  
    return 1684;
  else if( aa == "CCQ")  
    return 1685;
  else if( aa == "CCE")  
    return 1686;
  else if( aa == "CCG")  
    return 1687;
  else if( aa == "CCH")  
    return 1688;
  else if( aa == "CCI")  
    return 1689;
  else if( aa == "CCL")  
    return 1690;
  else if( aa == "CCK")  
    return 1691;
  else if( aa == "CCM")  
    return 1692;
  else if( aa == "CCF")  
    return 1693;
  else if( aa == "CCP")  
    return 1694;
  else if( aa == "CCS")  
    return 1695;
  else if( aa == "CCT")  
    return 1696;
  else if( aa == "CCW")  
    return 1697;
  else if( aa == "CCY")  
    return 1698;
  else if( aa == "CCV")  
    return 1699;
  else if( aa == "CQA")  
    return 1700;
  else if( aa == "CQR")  
    return 1701;
  else if( aa == "CQN")  
    return 1702;
  else if( aa == "CQD")  
    return 1703;
  else if( aa == "CQC")  
    return 1704;
  else if( aa == "CQQ")  
    return 1705;
  else if( aa == "CQE")  
    return 1706;
  else if( aa == "CQG")  
    return 1707;
  else if( aa == "CQH")  
    return 1708;
  else if( aa == "CQI")  
    return 1709;
  else if( aa == "CQL")  
    return 1710;
  else if( aa == "CQK")  
    return 1711;
  else if( aa == "CQM")  
    return 1712;
  else if( aa == "CQF")  
    return 1713;
  else if( aa == "CQP")  
    return 1714;
  else if( aa == "CQS")  
    return 1715;
  else if( aa == "CQT")  
    return 1716;
  else if( aa == "CQW")  
    return 1717;
  else if( aa == "CQY")  
    return 1718;
  else if( aa == "CQV")  
    return 1719;
  else if( aa == "CEA")  
    return 1720;
  else if( aa == "CER")  
    return 1721;
  else if( aa == "CEN")  
    return 1722;
  else if( aa == "CED")  
    return 1723;
  else if( aa == "CEC")  
    return 1724;
  else if( aa == "CEQ")  
    return 1725;
  else if( aa == "CEE")  
    return 1726;
  else if( aa == "CEG")  
    return 1727;
  else if( aa == "CEH")  
    return 1728;
  else if( aa == "CEI")  
    return 1729;
  else if( aa == "CEL")  
    return 1730;
  else if( aa == "CEK")  
    return 1731;
  else if( aa == "CEM")  
    return 1732;
  else if( aa == "CEF")  
    return 1733;
  else if( aa == "CEP")  
    return 1734;
  else if( aa == "CES")  
    return 1735;
  else if( aa == "CET")  
    return 1736;
  else if( aa == "CEW")  
    return 1737;
  else if( aa == "CEY")  
    return 1738;
  else if( aa == "CEV")  
    return 1739;
  else if( aa == "CGA")  
    return 1740;
  else if( aa == "CGR")  
    return 1741;
  else if( aa == "CGN")  
    return 1742;
  else if( aa == "CGD")  
    return 1743;
  else if( aa == "CGC")  
    return 1744;
  else if( aa == "CGQ")  
    return 1745;
  else if( aa == "CGE")  
    return 1746;
  else if( aa == "CGG")  
    return 1747;
  else if( aa == "CGH")  
    return 1748;
  else if( aa == "CGI")  
    return 1749;
  else if( aa == "CGL")  
    return 1750;
  else if( aa == "CGK")  
    return 1751;
  else if( aa == "CGM")  
    return 1752;
  else if( aa == "CGF")  
    return 1753;
  else if( aa == "CGP")  
    return 1754;
  else if( aa == "CGS")  
    return 1755;
  else if( aa == "CGT")  
    return 1756;
  else if( aa == "CGW")  
    return 1757;
  else if( aa == "CGY")  
    return 1758;
  else if( aa == "CGV")  
    return 1759;
  else if( aa == "CHA")  
    return 1760;
  else if( aa == "CHR")  
    return 1761;
  else if( aa == "CHN")  
    return 1762;
  else if( aa == "CHD")  
    return 1763;
  else if( aa == "CHC")  
    return 1764;
  else if( aa == "CHQ")  
    return 1765;
  else if( aa == "CHE")  
    return 1766;
  else if( aa == "CHG")  
    return 1767;
  else if( aa == "CHH")  
    return 1768;
  else if( aa == "CHI")  
    return 1769;
  else if( aa == "CHL")  
    return 1770;
  else if( aa == "CHK")  
    return 1771;
  else if( aa == "CHM")  
    return 1772;
  else if( aa == "CHF")  
    return 1773;
  else if( aa == "CHP")  
    return 1774;
  else if( aa == "CHS")  
    return 1775;
  else if( aa == "CHT")  
    return 1776;
  else if( aa == "CHW")  
    return 1777;
  else if( aa == "CHY")  
    return 1778;
  else if( aa == "CHV")  
    return 1779;
  else if( aa == "CIA")  
    return 1780;
  else if( aa == "CIR")  
    return 1781;
  else if( aa == "CIN")  
    return 1782;
  else if( aa == "CID")  
    return 1783;
  else if( aa == "CIC")  
    return 1784;
  else if( aa == "CIQ")  
    return 1785;
  else if( aa == "CIE")  
    return 1786;
  else if( aa == "CIG")  
    return 1787;
  else if( aa == "CIH")  
    return 1788;
  else if( aa == "CII")  
    return 1789;
  else if( aa == "CIL")  
    return 1790;
  else if( aa == "CIK")  
    return 1791;
  else if( aa == "CIM")  
    return 1792;
  else if( aa == "CIF")  
    return 1793;
  else if( aa == "CIP")  
    return 1794;
  else if( aa == "CIS")  
    return 1795;
  else if( aa == "CIT")  
    return 1796;
  else if( aa == "CIW")  
    return 1797;
  else if( aa == "CIY")  
    return 1798;
  else if( aa == "CIV")  
    return 1799;
  else if( aa == "CLA")  
    return 1800;
  else if( aa == "CLR")  
    return 1801;
  else if( aa == "CLN")  
    return 1802;
  else if( aa == "CLD")  
    return 1803;
  else if( aa == "CLC")  
    return 1804;
  else if( aa == "CLQ")  
    return 1805;
  else if( aa == "CLE")  
    return 1806;
  else if( aa == "CLG")  
    return 1807;
  else if( aa == "CLH")  
    return 1808;
  else if( aa == "CLI")  
    return 1809;
  else if( aa == "CLL")  
    return 1810;
  else if( aa == "CLK")  
    return 1811;
  else if( aa == "CLM")  
    return 1812;
  else if( aa == "CLF")  
    return 1813;
  else if( aa == "CLP")  
    return 1814;
  else if( aa == "CLS")  
    return 1815;
  else if( aa == "CLT")  
    return 1816;
  else if( aa == "CLW")  
    return 1817;
  else if( aa == "CLY")  
    return 1818;
  else if( aa == "CLV")  
    return 1819;
  else if( aa == "CKA")  
    return 1820;
  else if( aa == "CKR")  
    return 1821;
  else if( aa == "CKN")  
    return 1822;
  else if( aa == "CKD")  
    return 1823;
  else if( aa == "CKC")  
    return 1824;
  else if( aa == "CKQ")  
    return 1825;
  else if( aa == "CKE")  
    return 1826;
  else if( aa == "CKG")  
    return 1827;
  else if( aa == "CKH")  
    return 1828;
  else if( aa == "CKI")  
    return 1829;
  else if( aa == "CKL")  
    return 1830;
  else if( aa == "CKK")  
    return 1831;
  else if( aa == "CKM")  
    return 1832;
  else if( aa == "CKF")  
    return 1833;
  else if( aa == "CKP")  
    return 1834;
  else if( aa == "CKS")  
    return 1835;
  else if( aa == "CKT")  
    return 1836;
  else if( aa == "CKW")  
    return 1837;
  else if( aa == "CKY")  
    return 1838;
  else if( aa == "CKV")  
    return 1839;
  else if( aa == "CMA")  
    return 1840;
  else if( aa == "CMR")  
    return 1841;
  else if( aa == "CMN")  
    return 1842;
  else if( aa == "CMD")  
    return 1843;
  else if( aa == "CMC")  
    return 1844;
  else if( aa == "CMQ")  
    return 1845;
  else if( aa == "CME")  
    return 1846;
  else if( aa == "CMG")  
    return 1847;
  else if( aa == "CMH")  
    return 1848;
  else if( aa == "CMI")  
    return 1849;
  else if( aa == "CML")  
    return 1850;
  else if( aa == "CMK")  
    return 1851;
  else if( aa == "CMM")  
    return 1852;
  else if( aa == "CMF")  
    return 1853;
  else if( aa == "CMP")  
    return 1854;
  else if( aa == "CMS")  
    return 1855;
  else if( aa == "CMT")  
    return 1856;
  else if( aa == "CMW")  
    return 1857;
  else if( aa == "CMY")  
    return 1858;
  else if( aa == "CMV")  
    return 1859;
  else if( aa == "CFA")  
    return 1860;
  else if( aa == "CFR")  
    return 1861;
  else if( aa == "CFN")  
    return 1862;
  else if( aa == "CFD")  
    return 1863;
  else if( aa == "CFC")  
    return 1864;
  else if( aa == "CFQ")  
    return 1865;
  else if( aa == "CFE")  
    return 1866;
  else if( aa == "CFG")  
    return 1867;
  else if( aa == "CFH")  
    return 1868;
  else if( aa == "CFI")  
    return 1869;
  else if( aa == "CFL")  
    return 1870;
  else if( aa == "CFK")  
    return 1871;
  else if( aa == "CFM")  
    return 1872;
  else if( aa == "CFF")  
    return 1873;
  else if( aa == "CFP")  
    return 1874;
  else if( aa == "CFS")  
    return 1875;
  else if( aa == "CFT")  
    return 1876;
  else if( aa == "CFW")  
    return 1877;
  else if( aa == "CFY")  
    return 1878;
  else if( aa == "CFV")  
    return 1879;
  else if( aa == "CPA")  
    return 1880;
  else if( aa == "CPR")  
    return 1881;
  else if( aa == "CPN")  
    return 1882;
  else if( aa == "CPD")  
    return 1883;
  else if( aa == "CPC")  
    return 1884;
  else if( aa == "CPQ")  
    return 1885;
  else if( aa == "CPE")  
    return 1886;
  else if( aa == "CPG")  
    return 1887;
  else if( aa == "CPH")  
    return 1888;
  else if( aa == "CPI")  
    return 1889;
  else if( aa == "CPL")  
    return 1890;
  else if( aa == "CPK")  
    return 1891;
  else if( aa == "CPM")  
    return 1892;
  else if( aa == "CPF")  
    return 1893;
  else if( aa == "CPP")  
    return 1894;
  else if( aa == "CPS")  
    return 1895;
  else if( aa == "CPT")  
    return 1896;
  else if( aa == "CPW")  
    return 1897;
  else if( aa == "CPY")  
    return 1898;
  else if( aa == "CPV")  
    return 1899;
  else if( aa == "CSA")  
    return 1900;
  else if( aa == "CSR")  
    return 1901;
  else if( aa == "CSN")  
    return 1902;
  else if( aa == "CSD")  
    return 1903;
  else if( aa == "CSC")  
    return 1904;
  else if( aa == "CSQ")  
    return 1905;
  else if( aa == "CSE")  
    return 1906;
  else if( aa == "CSG")  
    return 1907;
  else if( aa == "CSH")  
    return 1908;
  else if( aa == "CSI")  
    return 1909;
  else if( aa == "CSL")  
    return 1910;
  else if( aa == "CSK")  
    return 1911;
  else if( aa == "CSM")  
    return 1912;
  else if( aa == "CSF")  
    return 1913;
  else if( aa == "CSP")  
    return 1914;
  else if( aa == "CSS")  
    return 1915;
  else if( aa == "CST")  
    return 1916;
  else if( aa == "CSW")  
    return 1917;
  else if( aa == "CSY")  
    return 1918;
  else if( aa == "CSV")  
    return 1919;
  else if( aa == "CTA")  
    return 1920;
  else if( aa == "CTR")  
    return 1921;
  else if( aa == "CTN")  
    return 1922;
  else if( aa == "CTD")  
    return 1923;
  else if( aa == "CTC")  
    return 1924;
  else if( aa == "CTQ")  
    return 1925;
  else if( aa == "CTE")  
    return 1926;
  else if( aa == "CTG")  
    return 1927;
  else if( aa == "CTH")  
    return 1928;
  else if( aa == "CTI")  
    return 1929;
  else if( aa == "CTL")  
    return 1930;
  else if( aa == "CTK")  
    return 1931;
  else if( aa == "CTM")  
    return 1932;
  else if( aa == "CTF")  
    return 1933;
  else if( aa == "CTP")  
    return 1934;
  else if( aa == "CTS")  
    return 1935;
  else if( aa == "CTT")  
    return 1936;
  else if( aa == "CTW")  
    return 1937;
  else if( aa == "CTY")  
    return 1938;
  else if( aa == "CTV")  
    return 1939;
  else if( aa == "CWA")  
    return 1940;
  else if( aa == "CWR")  
    return 1941;
  else if( aa == "CWN")  
    return 1942;
  else if( aa == "CWD")  
    return 1943;
  else if( aa == "CWC")  
    return 1944;
  else if( aa == "CWQ")  
    return 1945;
  else if( aa == "CWE")  
    return 1946;
  else if( aa == "CWG")  
    return 1947;
  else if( aa == "CWH")  
    return 1948;
  else if( aa == "CWI")  
    return 1949;
  else if( aa == "CWL")  
    return 1950;
  else if( aa == "CWK")  
    return 1951;
  else if( aa == "CWM")  
    return 1952;
  else if( aa == "CWF")  
    return 1953;
  else if( aa == "CWP")  
    return 1954;
  else if( aa == "CWS")  
    return 1955;
  else if( aa == "CWT")  
    return 1956;
  else if( aa == "CWW")  
    return 1957;
  else if( aa == "CWY")  
    return 1958;
  else if( aa == "CWV")  
    return 1959;
  else if( aa == "CYA")  
    return 1960;
  else if( aa == "CYR")  
    return 1961;
  else if( aa == "CYN")  
    return 1962;
  else if( aa == "CYD")  
    return 1963;
  else if( aa == "CYC")  
    return 1964;
  else if( aa == "CYQ")  
    return 1965;
  else if( aa == "CYE")  
    return 1966;
  else if( aa == "CYG")  
    return 1967;
  else if( aa == "CYH")  
    return 1968;
  else if( aa == "CYI")  
    return 1969;
  else if( aa == "CYL")  
    return 1970;
  else if( aa == "CYK")  
    return 1971;
  else if( aa == "CYM")  
    return 1972;
  else if( aa == "CYF")  
    return 1973;
  else if( aa == "CYP")  
    return 1974;
  else if( aa == "CYS")  
    return 1975;
  else if( aa == "CYT")  
    return 1976;
  else if( aa == "CYW")  
    return 1977;
  else if( aa == "CYY")  
    return 1978;
  else if( aa == "CYV")  
    return 1979;
  else if( aa == "CVA")  
    return 1980;
  else if( aa == "CVR")  
    return 1981;
  else if( aa == "CVN")  
    return 1982;
  else if( aa == "CVD")  
    return 1983;
  else if( aa == "CVC")  
    return 1984;
  else if( aa == "CVQ")  
    return 1985;
  else if( aa == "CVE")  
    return 1986;
  else if( aa == "CVG")  
    return 1987;
  else if( aa == "CVH")  
    return 1988;
  else if( aa == "CVI")  
    return 1989;
  else if( aa == "CVL")  
    return 1990;
  else if( aa == "CVK")  
    return 1991;
  else if( aa == "CVM")  
    return 1992;
  else if( aa == "CVF")  
    return 1993;
  else if( aa == "CVP")  
    return 1994;
  else if( aa == "CVS")  
    return 1995;
  else if( aa == "CVT")  
    return 1996;
  else if( aa == "CVW")  
    return 1997;
  else if( aa == "CVY")  
    return 1998;
  else if( aa == "CVV")  
    return 1999;
  else if( aa == "QAA")  
    return 2000;
  else if( aa == "QAR")  
    return 2001;
  else if( aa == "QAN")  
    return 2002;
  else if( aa == "QAD")  
    return 2003;
  else if( aa == "QAC")  
    return 2004;
  else if( aa == "QAQ")  
    return 2005;
  else if( aa == "QAE")  
    return 2006;
  else if( aa == "QAG")  
    return 2007;
  else if( aa == "QAH")  
    return 2008;
  else if( aa == "QAI")  
    return 2009;
  else if( aa == "QAL")  
    return 2010;
  else if( aa == "QAK")  
    return 2011;
  else if( aa == "QAM")  
    return 2012;
  else if( aa == "QAF")  
    return 2013;
  else if( aa == "QAP")  
    return 2014;
  else if( aa == "QAS")  
    return 2015;
  else if( aa == "QAT")  
    return 2016;
  else if( aa == "QAW")  
    return 2017;
  else if( aa == "QAY")  
    return 2018;
  else if( aa == "QAV")  
    return 2019;
  else if( aa == "QRA")  
    return 2020;
  else if( aa == "QRR")  
    return 2021;
  else if( aa == "QRN")  
    return 2022;
  else if( aa == "QRD")  
    return 2023;
  else if( aa == "QRC")  
    return 2024;
  else if( aa == "QRQ")  
    return 2025;
  else if( aa == "QRE")  
    return 2026;
  else if( aa == "QRG")  
    return 2027;
  else if( aa == "QRH")  
    return 2028;
  else if( aa == "QRI")  
    return 2029;
  else if( aa == "QRL")  
    return 2030;
  else if( aa == "QRK")  
    return 2031;
  else if( aa == "QRM")  
    return 2032;
  else if( aa == "QRF")  
    return 2033;
  else if( aa == "QRP")  
    return 2034;
  else if( aa == "QRS")  
    return 2035;
  else if( aa == "QRT")  
    return 2036;
  else if( aa == "QRW")  
    return 2037;
  else if( aa == "QRY")  
    return 2038;
  else if( aa == "QRV")  
    return 2039;
  else if( aa == "QNA")  
    return 2040;
  else if( aa == "QNR")  
    return 2041;
  else if( aa == "QNN")  
    return 2042;
  else if( aa == "QND")  
    return 2043;
  else if( aa == "QNC")  
    return 2044;
  else if( aa == "QNQ")  
    return 2045;
  else if( aa == "QNE")  
    return 2046;
  else if( aa == "QNG")  
    return 2047;
  else if( aa == "QNH")  
    return 2048;
  else if( aa == "QNI")  
    return 2049;
  else if( aa == "QNL")  
    return 2050;
  else if( aa == "QNK")  
    return 2051;
  else if( aa == "QNM")  
    return 2052;
  else if( aa == "QNF")  
    return 2053;
  else if( aa == "QNP")  
    return 2054;
  else if( aa == "QNS")  
    return 2055;
  else if( aa == "QNT")  
    return 2056;
  else if( aa == "QNW")  
    return 2057;
  else if( aa == "QNY")  
    return 2058;
  else if( aa == "QNV")  
    return 2059;
  else if( aa == "QDA")  
    return 2060;
  else if( aa == "QDR")  
    return 2061;
  else if( aa == "QDN")  
    return 2062;
  else if( aa == "QDD")  
    return 2063;
  else if( aa == "QDC")  
    return 2064;
  else if( aa == "QDQ")  
    return 2065;
  else if( aa == "QDE")  
    return 2066;
  else if( aa == "QDG")  
    return 2067;
  else if( aa == "QDH")  
    return 2068;
  else if( aa == "QDI")  
    return 2069;
  else if( aa == "QDL")  
    return 2070;
  else if( aa == "QDK")  
    return 2071;
  else if( aa == "QDM")  
    return 2072;
  else if( aa == "QDF")  
    return 2073;
  else if( aa == "QDP")  
    return 2074;
  else if( aa == "QDS")  
    return 2075;
  else if( aa == "QDT")  
    return 2076;
  else if( aa == "QDW")  
    return 2077;
  else if( aa == "QDY")  
    return 2078;
  else if( aa == "QDV")  
    return 2079;
  else if( aa == "QCA")  
    return 2080;
  else if( aa == "QCR")  
    return 2081;
  else if( aa == "QCN")  
    return 2082;
  else if( aa == "QCD")  
    return 2083;
  else if( aa == "QCC")  
    return 2084;
  else if( aa == "QCQ")  
    return 2085;
  else if( aa == "QCE")  
    return 2086;
  else if( aa == "QCG")  
    return 2087;
  else if( aa == "QCH")  
    return 2088;
  else if( aa == "QCI")  
    return 2089;
  else if( aa == "QCL")  
    return 2090;
  else if( aa == "QCK")  
    return 2091;
  else if( aa == "QCM")  
    return 2092;
  else if( aa == "QCF")  
    return 2093;
  else if( aa == "QCP")  
    return 2094;
  else if( aa == "QCS")  
    return 2095;
  else if( aa == "QCT")  
    return 2096;
  else if( aa == "QCW")  
    return 2097;
  else if( aa == "QCY")  
    return 2098;
  else if( aa == "QCV")  
    return 2099;
  else if( aa == "QQA")  
    return 2100;
  else if( aa == "QQR")  
    return 2101;
  else if( aa == "QQN")  
    return 2102;
  else if( aa == "QQD")  
    return 2103;
  else if( aa == "QQC")  
    return 2104;
  else if( aa == "QQQ")  
    return 2105;
  else if( aa == "QQE")  
    return 2106;
  else if( aa == "QQG")  
    return 2107;
  else if( aa == "QQH")  
    return 2108;
  else if( aa == "QQI")  
    return 2109;
  else if( aa == "QQL")  
    return 2110;
  else if( aa == "QQK")  
    return 2111;
  else if( aa == "QQM")  
    return 2112;
  else if( aa == "QQF")  
    return 2113;
  else if( aa == "QQP")  
    return 2114;
  else if( aa == "QQS")  
    return 2115;
  else if( aa == "QQT")  
    return 2116;
  else if( aa == "QQW")  
    return 2117;
  else if( aa == "QQY")  
    return 2118;
  else if( aa == "QQV")  
    return 2119;
  else if( aa == "QEA")  
    return 2120;
  else if( aa == "QER")  
    return 2121;
  else if( aa == "QEN")  
    return 2122;
  else if( aa == "QED")  
    return 2123;
  else if( aa == "QEC")  
    return 2124;
  else if( aa == "QEQ")  
    return 2125;
  else if( aa == "QEE")  
    return 2126;
  else if( aa == "QEG")  
    return 2127;
  else if( aa == "QEH")  
    return 2128;
  else if( aa == "QEI")  
    return 2129;
  else if( aa == "QEL")  
    return 2130;
  else if( aa == "QEK")  
    return 2131;
  else if( aa == "QEM")  
    return 2132;
  else if( aa == "QEF")  
    return 2133;
  else if( aa == "QEP")  
    return 2134;
  else if( aa == "QES")  
    return 2135;
  else if( aa == "QET")  
    return 2136;
  else if( aa == "QEW")  
    return 2137;
  else if( aa == "QEY")  
    return 2138;
  else if( aa == "QEV")  
    return 2139;
  else if( aa == "QGA")  
    return 2140;
  else if( aa == "QGR")  
    return 2141;
  else if( aa == "QGN")  
    return 2142;
  else if( aa == "QGD")  
    return 2143;
  else if( aa == "QGC")  
    return 2144;
  else if( aa == "QGQ")  
    return 2145;
  else if( aa == "QGE")  
    return 2146;
  else if( aa == "QGG")  
    return 2147;
  else if( aa == "QGH")  
    return 2148;
  else if( aa == "QGI")  
    return 2149;
  else if( aa == "QGL")  
    return 2150;
  else if( aa == "QGK")  
    return 2151;
  else if( aa == "QGM")  
    return 2152;
  else if( aa == "QGF")  
    return 2153;
  else if( aa == "QGP")  
    return 2154;
  else if( aa == "QGS")  
    return 2155;
  else if( aa == "QGT")  
    return 2156;
  else if( aa == "QGW")  
    return 2157;
  else if( aa == "QGY")  
    return 2158;
  else if( aa == "QGV")  
    return 2159;
  else if( aa == "QHA")  
    return 2160;
  else if( aa == "QHR")  
    return 2161;
  else if( aa == "QHN")  
    return 2162;
  else if( aa == "QHD")  
    return 2163;
  else if( aa == "QHC")  
    return 2164;
  else if( aa == "QHQ")  
    return 2165;
  else if( aa == "QHE")  
    return 2166;
  else if( aa == "QHG")  
    return 2167;
  else if( aa == "QHH")  
    return 2168;
  else if( aa == "QHI")  
    return 2169;
  else if( aa == "QHL")  
    return 2170;
  else if( aa == "QHK")  
    return 2171;
  else if( aa == "QHM")  
    return 2172;
  else if( aa == "QHF")  
    return 2173;
  else if( aa == "QHP")  
    return 2174;
  else if( aa == "QHS")  
    return 2175;
  else if( aa == "QHT")  
    return 2176;
  else if( aa == "QHW")  
    return 2177;
  else if( aa == "QHY")  
    return 2178;
  else if( aa == "QHV")  
    return 2179;
  else if( aa == "QIA")  
    return 2180;
  else if( aa == "QIR")  
    return 2181;
  else if( aa == "QIN")  
    return 2182;
  else if( aa == "QID")  
    return 2183;
  else if( aa == "QIC")  
    return 2184;
  else if( aa == "QIQ")  
    return 2185;
  else if( aa == "QIE")  
    return 2186;
  else if( aa == "QIG")  
    return 2187;
  else if( aa == "QIH")  
    return 2188;
  else if( aa == "QII")  
    return 2189;
  else if( aa == "QIL")  
    return 2190;
  else if( aa == "QIK")  
    return 2191;
  else if( aa == "QIM")  
    return 2192;
  else if( aa == "QIF")  
    return 2193;
  else if( aa == "QIP")  
    return 2194;
  else if( aa == "QIS")  
    return 2195;
  else if( aa == "QIT")  
    return 2196;
  else if( aa == "QIW")  
    return 2197;
  else if( aa == "QIY")  
    return 2198;
  else if( aa == "QIV")  
    return 2199;
  else if( aa == "QLA")  
    return 2200;
  else if( aa == "QLR")  
    return 2201;
  else if( aa == "QLN")  
    return 2202;
  else if( aa == "QLD")  
    return 2203;
  else if( aa == "QLC")  
    return 2204;
  else if( aa == "QLQ")  
    return 2205;
  else if( aa == "QLE")  
    return 2206;
  else if( aa == "QLG")  
    return 2207;
  else if( aa == "QLH")  
    return 2208;
  else if( aa == "QLI")  
    return 2209;
  else if( aa == "QLL")  
    return 2210;
  else if( aa == "QLK")  
    return 2211;
  else if( aa == "QLM")  
    return 2212;
  else if( aa == "QLF")  
    return 2213;
  else if( aa == "QLP")  
    return 2214;
  else if( aa == "QLS")  
    return 2215;
  else if( aa == "QLT")  
    return 2216;
  else if( aa == "QLW")  
    return 2217;
  else if( aa == "QLY")  
    return 2218;
  else if( aa == "QLV")  
    return 2219;
  else if( aa == "QKA")  
    return 2220;
  else if( aa == "QKR")  
    return 2221;
  else if( aa == "QKN")  
    return 2222;
  else if( aa == "QKD")  
    return 2223;
  else if( aa == "QKC")  
    return 2224;
  else if( aa == "QKQ")  
    return 2225;
  else if( aa == "QKE")  
    return 2226;
  else if( aa == "QKG")  
    return 2227;
  else if( aa == "QKH")  
    return 2228;
  else if( aa == "QKI")  
    return 2229;
  else if( aa == "QKL")  
    return 2230;
  else if( aa == "QKK")  
    return 2231;
  else if( aa == "QKM")  
    return 2232;
  else if( aa == "QKF")  
    return 2233;
  else if( aa == "QKP")  
    return 2234;
  else if( aa == "QKS")  
    return 2235;
  else if( aa == "QKT")  
    return 2236;
  else if( aa == "QKW")  
    return 2237;
  else if( aa == "QKY")  
    return 2238;
  else if( aa == "QKV")  
    return 2239;
  else if( aa == "QMA")  
    return 2240;
  else if( aa == "QMR")  
    return 2241;
  else if( aa == "QMN")  
    return 2242;
  else if( aa == "QMD")  
    return 2243;
  else if( aa == "QMC")  
    return 2244;
  else if( aa == "QMQ")  
    return 2245;
  else if( aa == "QME")  
    return 2246;
  else if( aa == "QMG")  
    return 2247;
  else if( aa == "QMH")  
    return 2248;
  else if( aa == "QMI")  
    return 2249;
  else if( aa == "QML")  
    return 2250;
  else if( aa == "QMK")  
    return 2251;
  else if( aa == "QMM")  
    return 2252;
  else if( aa == "QMF")  
    return 2253;
  else if( aa == "QMP")  
    return 2254;
  else if( aa == "QMS")  
    return 2255;
  else if( aa == "QMT")  
    return 2256;
  else if( aa == "QMW")  
    return 2257;
  else if( aa == "QMY")  
    return 2258;
  else if( aa == "QMV")  
    return 2259;
  else if( aa == "QFA")  
    return 2260;
  else if( aa == "QFR")  
    return 2261;
  else if( aa == "QFN")  
    return 2262;
  else if( aa == "QFD")  
    return 2263;
  else if( aa == "QFC")  
    return 2264;
  else if( aa == "QFQ")  
    return 2265;
  else if( aa == "QFE")  
    return 2266;
  else if( aa == "QFG")  
    return 2267;
  else if( aa == "QFH")  
    return 2268;
  else if( aa == "QFI")  
    return 2269;
  else if( aa == "QFL")  
    return 2270;
  else if( aa == "QFK")  
    return 2271;
  else if( aa == "QFM")  
    return 2272;
  else if( aa == "QFF")  
    return 2273;
  else if( aa == "QFP")  
    return 2274;
  else if( aa == "QFS")  
    return 2275;
  else if( aa == "QFT")  
    return 2276;
  else if( aa == "QFW")  
    return 2277;
  else if( aa == "QFY")  
    return 2278;
  else if( aa == "QFV")  
    return 2279;
  else if( aa == "QPA")  
    return 2280;
  else if( aa == "QPR")  
    return 2281;
  else if( aa == "QPN")  
    return 2282;
  else if( aa == "QPD")  
    return 2283;
  else if( aa == "QPC")  
    return 2284;
  else if( aa == "QPQ")  
    return 2285;
  else if( aa == "QPE")  
    return 2286;
  else if( aa == "QPG")  
    return 2287;
  else if( aa == "QPH")  
    return 2288;
  else if( aa == "QPI")  
    return 2289;
  else if( aa == "QPL")  
    return 2290;
  else if( aa == "QPK")  
    return 2291;
  else if( aa == "QPM")  
    return 2292;
  else if( aa == "QPF")  
    return 2293;
  else if( aa == "QPP")  
    return 2294;
  else if( aa == "QPS")  
    return 2295;
  else if( aa == "QPT")  
    return 2296;
  else if( aa == "QPW")  
    return 2297;
  else if( aa == "QPY")  
    return 2298;
  else if( aa == "QPV")  
    return 2299;
  else if( aa == "QSA")  
    return 2300;
  else if( aa == "QSR")  
    return 2301;
  else if( aa == "QSN")  
    return 2302;
  else if( aa == "QSD")  
    return 2303;
  else if( aa == "QSC")  
    return 2304;
  else if( aa == "QSQ")  
    return 2305;
  else if( aa == "QSE")  
    return 2306;
  else if( aa == "QSG")  
    return 2307;
  else if( aa == "QSH")  
    return 2308;
  else if( aa == "QSI")  
    return 2309;
  else if( aa == "QSL")  
    return 2310;
  else if( aa == "QSK")  
    return 2311;
  else if( aa == "QSM")  
    return 2312;
  else if( aa == "QSF")  
    return 2313;
  else if( aa == "QSP")  
    return 2314;
  else if( aa == "QSS")  
    return 2315;
  else if( aa == "QST")  
    return 2316;
  else if( aa == "QSW")  
    return 2317;
  else if( aa == "QSY")  
    return 2318;
  else if( aa == "QSV")  
    return 2319;
  else if( aa == "QTA")  
    return 2320;
  else if( aa == "QTR")  
    return 2321;
  else if( aa == "QTN")  
    return 2322;
  else if( aa == "QTD")  
    return 2323;
  else if( aa == "QTC")  
    return 2324;
  else if( aa == "QTQ")  
    return 2325;
  else if( aa == "QTE")  
    return 2326;
  else if( aa == "QTG")  
    return 2327;
  else if( aa == "QTH")  
    return 2328;
  else if( aa == "QTI")  
    return 2329;
  else if( aa == "QTL")  
    return 2330;
  else if( aa == "QTK")  
    return 2331;
  else if( aa == "QTM")  
    return 2332;
  else if( aa == "QTF")  
    return 2333;
  else if( aa == "QTP")  
    return 2334;
  else if( aa == "QTS")  
    return 2335;
  else if( aa == "QTT")  
    return 2336;
  else if( aa == "QTW")  
    return 2337;
  else if( aa == "QTY")  
    return 2338;
  else if( aa == "QTV")  
    return 2339;
  else if( aa == "QWA")  
    return 2340;
  else if( aa == "QWR")  
    return 2341;
  else if( aa == "QWN")  
    return 2342;
  else if( aa == "QWD")  
    return 2343;
  else if( aa == "QWC")  
    return 2344;
  else if( aa == "QWQ")  
    return 2345;
  else if( aa == "QWE")  
    return 2346;
  else if( aa == "QWG")  
    return 2347;
  else if( aa == "QWH")  
    return 2348;
  else if( aa == "QWI")  
    return 2349;
  else if( aa == "QWL")  
    return 2350;
  else if( aa == "QWK")  
    return 2351;
  else if( aa == "QWM")  
    return 2352;
  else if( aa == "QWF")  
    return 2353;
  else if( aa == "QWP")  
    return 2354;
  else if( aa == "QWS")  
    return 2355;
  else if( aa == "QWT")  
    return 2356;
  else if( aa == "QWW")  
    return 2357;
  else if( aa == "QWY")  
    return 2358;
  else if( aa == "QWV")  
    return 2359;
  else if( aa == "QYA")  
    return 2360;
  else if( aa == "QYR")  
    return 2361;
  else if( aa == "QYN")  
    return 2362;
  else if( aa == "QYD")  
    return 2363;
  else if( aa == "QYC")  
    return 2364;
  else if( aa == "QYQ")  
    return 2365;
  else if( aa == "QYE")  
    return 2366;
  else if( aa == "QYG")  
    return 2367;
  else if( aa == "QYH")  
    return 2368;
  else if( aa == "QYI")  
    return 2369;
  else if( aa == "QYL")  
    return 2370;
  else if( aa == "QYK")  
    return 2371;
  else if( aa == "QYM")  
    return 2372;
  else if( aa == "QYF")  
    return 2373;
  else if( aa == "QYP")  
    return 2374;
  else if( aa == "QYS")  
    return 2375;
  else if( aa == "QYT")  
    return 2376;
  else if( aa == "QYW")  
    return 2377;
  else if( aa == "QYY")  
    return 2378;
  else if( aa == "QYV")  
    return 2379;
  else if( aa == "QVA")  
    return 2380;
  else if( aa == "QVR")  
    return 2381;
  else if( aa == "QVN")  
    return 2382;
  else if( aa == "QVD")  
    return 2383;
  else if( aa == "QVC")  
    return 2384;
  else if( aa == "QVQ")  
    return 2385;
  else if( aa == "QVE")  
    return 2386;
  else if( aa == "QVG")  
    return 2387;
  else if( aa == "QVH")  
    return 2388;
  else if( aa == "QVI")  
    return 2389;
  else if( aa == "QVL")  
    return 2390;
  else if( aa == "QVK")  
    return 2391;
  else if( aa == "QVM")  
    return 2392;
  else if( aa == "QVF")  
    return 2393;
  else if( aa == "QVP")  
    return 2394;
  else if( aa == "QVS")  
    return 2395;
  else if( aa == "QVT")  
    return 2396;
  else if( aa == "QVW")  
    return 2397;
  else if( aa == "QVY")  
    return 2398;
  else if( aa == "QVV")  
    return 2399;
  else if( aa == "EAA")  
    return 2400;
  else if( aa == "EAR")  
    return 2401;
  else if( aa == "EAN")  
    return 2402;
  else if( aa == "EAD")  
    return 2403;
  else if( aa == "EAC")  
    return 2404;
  else if( aa == "EAQ")  
    return 2405;
  else if( aa == "EAE")  
    return 2406;
  else if( aa == "EAG")  
    return 2407;
  else if( aa == "EAH")  
    return 2408;
  else if( aa == "EAI")  
    return 2409;
  else if( aa == "EAL")  
    return 2410;
  else if( aa == "EAK")  
    return 2411;
  else if( aa == "EAM")  
    return 2412;
  else if( aa == "EAF")  
    return 2413;
  else if( aa == "EAP")  
    return 2414;
  else if( aa == "EAS")  
    return 2415;
  else if( aa == "EAT")  
    return 2416;
  else if( aa == "EAW")  
    return 2417;
  else if( aa == "EAY")  
    return 2418;
  else if( aa == "EAV")  
    return 2419;
  else if( aa == "ERA")  
    return 2420;
  else if( aa == "ERR")  
    return 2421;
  else if( aa == "ERN")  
    return 2422;
  else if( aa == "ERD")  
    return 2423;
  else if( aa == "ERC")  
    return 2424;
  else if( aa == "ERQ")  
    return 2425;
  else if( aa == "ERE")  
    return 2426;
  else if( aa == "ERG")  
    return 2427;
  else if( aa == "ERH")  
    return 2428;
  else if( aa == "ERI")  
    return 2429;
  else if( aa == "ERL")  
    return 2430;
  else if( aa == "ERK")  
    return 2431;
  else if( aa == "ERM")  
    return 2432;
  else if( aa == "ERF")  
    return 2433;
  else if( aa == "ERP")  
    return 2434;
  else if( aa == "ERS")  
    return 2435;
  else if( aa == "ERT")  
    return 2436;
  else if( aa == "ERW")  
    return 2437;
  else if( aa == "ERY")  
    return 2438;
  else if( aa == "ERV")  
    return 2439;
  else if( aa == "ENA")  
    return 2440;
  else if( aa == "ENR")  
    return 2441;
  else if( aa == "ENN")  
    return 2442;
  else if( aa == "END")  
    return 2443;
  else if( aa == "ENC")  
    return 2444;
  else if( aa == "ENQ")  
    return 2445;
  else if( aa == "ENE")  
    return 2446;
  else if( aa == "ENG")  
    return 2447;
  else if( aa == "ENH")  
    return 2448;
  else if( aa == "ENI")  
    return 2449;
  else if( aa == "ENL")  
    return 2450;
  else if( aa == "ENK")  
    return 2451;
  else if( aa == "ENM")  
    return 2452;
  else if( aa == "ENF")  
    return 2453;
  else if( aa == "ENP")  
    return 2454;
  else if( aa == "ENS")  
    return 2455;
  else if( aa == "ENT")  
    return 2456;
  else if( aa == "ENW")  
    return 2457;
  else if( aa == "ENY")  
    return 2458;
  else if( aa == "ENV")  
    return 2459;
  else if( aa == "EDA")  
    return 2460;
  else if( aa == "EDR")  
    return 2461;
  else if( aa == "EDN")  
    return 2462;
  else if( aa == "EDD")  
    return 2463;
  else if( aa == "EDC")  
    return 2464;
  else if( aa == "EDQ")  
    return 2465;
  else if( aa == "EDE")  
    return 2466;
  else if( aa == "EDG")  
    return 2467;
  else if( aa == "EDH")  
    return 2468;
  else if( aa == "EDI")  
    return 2469;
  else if( aa == "EDL")  
    return 2470;
  else if( aa == "EDK")  
    return 2471;
  else if( aa == "EDM")  
    return 2472;
  else if( aa == "EDF")  
    return 2473;
  else if( aa == "EDP")  
    return 2474;
  else if( aa == "EDS")  
    return 2475;
  else if( aa == "EDT")  
    return 2476;
  else if( aa == "EDW")  
    return 2477;
  else if( aa == "EDY")  
    return 2478;
  else if( aa == "EDV")  
    return 2479;
  else if( aa == "ECA")  
    return 2480;
  else if( aa == "ECR")  
    return 2481;
  else if( aa == "ECN")  
    return 2482;
  else if( aa == "ECD")  
    return 2483;
  else if( aa == "ECC")  
    return 2484;
  else if( aa == "ECQ")  
    return 2485;
  else if( aa == "ECE")  
    return 2486;
  else if( aa == "ECG")  
    return 2487;
  else if( aa == "ECH")  
    return 2488;
  else if( aa == "ECI")  
    return 2489;
  else if( aa == "ECL")  
    return 2490;
  else if( aa == "ECK")  
    return 2491;
  else if( aa == "ECM")  
    return 2492;
  else if( aa == "ECF")  
    return 2493;
  else if( aa == "ECP")  
    return 2494;
  else if( aa == "ECS")  
    return 2495;
  else if( aa == "ECT")  
    return 2496;
  else if( aa == "ECW")  
    return 2497;
  else if( aa == "ECY")  
    return 2498;
  else if( aa == "ECV")  
    return 2499;
  else if( aa == "EQA")  
    return 2500;
  else if( aa == "EQR")  
    return 2501;
  else if( aa == "EQN")  
    return 2502;
  else if( aa == "EQD")  
    return 2503;
  else if( aa == "EQC")  
    return 2504;
  else if( aa == "EQQ")  
    return 2505;
  else if( aa == "EQE")  
    return 2506;
  else if( aa == "EQG")  
    return 2507;
  else if( aa == "EQH")  
    return 2508;
  else if( aa == "EQI")  
    return 2509;
  else if( aa == "EQL")  
    return 2510;
  else if( aa == "EQK")  
    return 2511;
  else if( aa == "EQM")  
    return 2512;
  else if( aa == "EQF")  
    return 2513;
  else if( aa == "EQP")  
    return 2514;
  else if( aa == "EQS")  
    return 2515;
  else if( aa == "EQT")  
    return 2516;
  else if( aa == "EQW")  
    return 2517;
  else if( aa == "EQY")  
    return 2518;
  else if( aa == "EQV")  
    return 2519;
  else if( aa == "EEA")  
    return 2520;
  else if( aa == "EER")  
    return 2521;
  else if( aa == "EEN")  
    return 2522;
  else if( aa == "EED")  
    return 2523;
  else if( aa == "EEC")  
    return 2524;
  else if( aa == "EEQ")  
    return 2525;
  else if( aa == "EEE")  
    return 2526;
  else if( aa == "EEG")  
    return 2527;
  else if( aa == "EEH")  
    return 2528;
  else if( aa == "EEI")  
    return 2529;
  else if( aa == "EEL")  
    return 2530;
  else if( aa == "EEK")  
    return 2531;
  else if( aa == "EEM")  
    return 2532;
  else if( aa == "EEF")  
    return 2533;
  else if( aa == "EEP")  
    return 2534;
  else if( aa == "EES")  
    return 2535;
  else if( aa == "EET")  
    return 2536;
  else if( aa == "EEW")  
    return 2537;
  else if( aa == "EEY")  
    return 2538;
  else if( aa == "EEV")  
    return 2539;
  else if( aa == "EGA")  
    return 2540;
  else if( aa == "EGR")  
    return 2541;
  else if( aa == "EGN")  
    return 2542;
  else if( aa == "EGD")  
    return 2543;
  else if( aa == "EGC")  
    return 2544;
  else if( aa == "EGQ")  
    return 2545;
  else if( aa == "EGE")  
    return 2546;
  else if( aa == "EGG")  
    return 2547;
  else if( aa == "EGH")  
    return 2548;
  else if( aa == "EGI")  
    return 2549;
  else if( aa == "EGL")  
    return 2550;
  else if( aa == "EGK")  
    return 2551;
  else if( aa == "EGM")  
    return 2552;
  else if( aa == "EGF")  
    return 2553;
  else if( aa == "EGP")  
    return 2554;
  else if( aa == "EGS")  
    return 2555;
  else if( aa == "EGT")  
    return 2556;
  else if( aa == "EGW")  
    return 2557;
  else if( aa == "EGY")  
    return 2558;
  else if( aa == "EGV")  
    return 2559;
  else if( aa == "EHA")  
    return 2560;
  else if( aa == "EHR")  
    return 2561;
  else if( aa == "EHN")  
    return 2562;
  else if( aa == "EHD")  
    return 2563;
  else if( aa == "EHC")  
    return 2564;
  else if( aa == "EHQ")  
    return 2565;
  else if( aa == "EHE")  
    return 2566;
  else if( aa == "EHG")  
    return 2567;
  else if( aa == "EHH")  
    return 2568;
  else if( aa == "EHI")  
    return 2569;
  else if( aa == "EHL")  
    return 2570;
  else if( aa == "EHK")  
    return 2571;
  else if( aa == "EHM")  
    return 2572;
  else if( aa == "EHF")  
    return 2573;
  else if( aa == "EHP")  
    return 2574;
  else if( aa == "EHS")  
    return 2575;
  else if( aa == "EHT")  
    return 2576;
  else if( aa == "EHW")  
    return 2577;
  else if( aa == "EHY")  
    return 2578;
  else if( aa == "EHV")  
    return 2579;
  else if( aa == "EIA")  
    return 2580;
  else if( aa == "EIR")  
    return 2581;
  else if( aa == "EIN")  
    return 2582;
  else if( aa == "EID")  
    return 2583;
  else if( aa == "EIC")  
    return 2584;
  else if( aa == "EIQ")  
    return 2585;
  else if( aa == "EIE")  
    return 2586;
  else if( aa == "EIG")  
    return 2587;
  else if( aa == "EIH")  
    return 2588;
  else if( aa == "EII")  
    return 2589;
  else if( aa == "EIL")  
    return 2590;
  else if( aa == "EIK")  
    return 2591;
  else if( aa == "EIM")  
    return 2592;
  else if( aa == "EIF")  
    return 2593;
  else if( aa == "EIP")  
    return 2594;
  else if( aa == "EIS")  
    return 2595;
  else if( aa == "EIT")  
    return 2596;
  else if( aa == "EIW")  
    return 2597;
  else if( aa == "EIY")  
    return 2598;
  else if( aa == "EIV")  
    return 2599;
  else if( aa == "ELA")  
    return 2600;
  else if( aa == "ELR")  
    return 2601;
  else if( aa == "ELN")  
    return 2602;
  else if( aa == "ELD")  
    return 2603;
  else if( aa == "ELC")  
    return 2604;
  else if( aa == "ELQ")  
    return 2605;
  else if( aa == "ELE")  
    return 2606;
  else if( aa == "ELG")  
    return 2607;
  else if( aa == "ELH")  
    return 2608;
  else if( aa == "ELI")  
    return 2609;
  else if( aa == "ELL")  
    return 2610;
  else if( aa == "ELK")  
    return 2611;
  else if( aa == "ELM")  
    return 2612;
  else if( aa == "ELF")  
    return 2613;
  else if( aa == "ELP")  
    return 2614;
  else if( aa == "ELS")  
    return 2615;
  else if( aa == "ELT")  
    return 2616;
  else if( aa == "ELW")  
    return 2617;
  else if( aa == "ELY")  
    return 2618;
  else if( aa == "ELV")  
    return 2619;
  else if( aa == "EKA")  
    return 2620;
  else if( aa == "EKR")  
    return 2621;
  else if( aa == "EKN")  
    return 2622;
  else if( aa == "EKD")  
    return 2623;
  else if( aa == "EKC")  
    return 2624;
  else if( aa == "EKQ")  
    return 2625;
  else if( aa == "EKE")  
    return 2626;
  else if( aa == "EKG")  
    return 2627;
  else if( aa == "EKH")  
    return 2628;
  else if( aa == "EKI")  
    return 2629;
  else if( aa == "EKL")  
    return 2630;
  else if( aa == "EKK")  
    return 2631;
  else if( aa == "EKM")  
    return 2632;
  else if( aa == "EKF")  
    return 2633;
  else if( aa == "EKP")  
    return 2634;
  else if( aa == "EKS")  
    return 2635;
  else if( aa == "EKT")  
    return 2636;
  else if( aa == "EKW")  
    return 2637;
  else if( aa == "EKY")  
    return 2638;
  else if( aa == "EKV")  
    return 2639;
  else if( aa == "EMA")  
    return 2640;
  else if( aa == "EMR")  
    return 2641;
  else if( aa == "EMN")  
    return 2642;
  else if( aa == "EMD")  
    return 2643;
  else if( aa == "EMC")  
    return 2644;
  else if( aa == "EMQ")  
    return 2645;
  else if( aa == "EME")  
    return 2646;
  else if( aa == "EMG")  
    return 2647;
  else if( aa == "EMH")  
    return 2648;
  else if( aa == "EMI")  
    return 2649;
  else if( aa == "EML")  
    return 2650;
  else if( aa == "EMK")  
    return 2651;
  else if( aa == "EMM")  
    return 2652;
  else if( aa == "EMF")  
    return 2653;
  else if( aa == "EMP")  
    return 2654;
  else if( aa == "EMS")  
    return 2655;
  else if( aa == "EMT")  
    return 2656;
  else if( aa == "EMW")  
    return 2657;
  else if( aa == "EMY")  
    return 2658;
  else if( aa == "EMV")  
    return 2659;
  else if( aa == "EFA")  
    return 2660;
  else if( aa == "EFR")  
    return 2661;
  else if( aa == "EFN")  
    return 2662;
  else if( aa == "EFD")  
    return 2663;
  else if( aa == "EFC")  
    return 2664;
  else if( aa == "EFQ")  
    return 2665;
  else if( aa == "EFE")  
    return 2666;
  else if( aa == "EFG")  
    return 2667;
  else if( aa == "EFH")  
    return 2668;
  else if( aa == "EFI")  
    return 2669;
  else if( aa == "EFL")  
    return 2670;
  else if( aa == "EFK")  
    return 2671;
  else if( aa == "EFM")  
    return 2672;
  else if( aa == "EFF")  
    return 2673;
  else if( aa == "EFP")  
    return 2674;
  else if( aa == "EFS")  
    return 2675;
  else if( aa == "EFT")  
    return 2676;
  else if( aa == "EFW")  
    return 2677;
  else if( aa == "EFY")  
    return 2678;
  else if( aa == "EFV")  
    return 2679;
  else if( aa == "EPA")  
    return 2680;
  else if( aa == "EPR")  
    return 2681;
  else if( aa == "EPN")  
    return 2682;
  else if( aa == "EPD")  
    return 2683;
  else if( aa == "EPC")  
    return 2684;
  else if( aa == "EPQ")  
    return 2685;
  else if( aa == "EPE")  
    return 2686;
  else if( aa == "EPG")  
    return 2687;
  else if( aa == "EPH")  
    return 2688;
  else if( aa == "EPI")  
    return 2689;
  else if( aa == "EPL")  
    return 2690;
  else if( aa == "EPK")  
    return 2691;
  else if( aa == "EPM")  
    return 2692;
  else if( aa == "EPF")  
    return 2693;
  else if( aa == "EPP")  
    return 2694;
  else if( aa == "EPS")  
    return 2695;
  else if( aa == "EPT")  
    return 2696;
  else if( aa == "EPW")  
    return 2697;
  else if( aa == "EPY")  
    return 2698;
  else if( aa == "EPV")  
    return 2699;
  else if( aa == "ESA")  
    return 2700;
  else if( aa == "ESR")  
    return 2701;
  else if( aa == "ESN")  
    return 2702;
  else if( aa == "ESD")  
    return 2703;
  else if( aa == "ESC")  
    return 2704;
  else if( aa == "ESQ")  
    return 2705;
  else if( aa == "ESE")  
    return 2706;
  else if( aa == "ESG")  
    return 2707;
  else if( aa == "ESH")  
    return 2708;
  else if( aa == "ESI")  
    return 2709;
  else if( aa == "ESL")  
    return 2710;
  else if( aa == "ESK")  
    return 2711;
  else if( aa == "ESM")  
    return 2712;
  else if( aa == "ESF")  
    return 2713;
  else if( aa == "ESP")  
    return 2714;
  else if( aa == "ESS")  
    return 2715;
  else if( aa == "EST")  
    return 2716;
  else if( aa == "ESW")  
    return 2717;
  else if( aa == "ESY")  
    return 2718;
  else if( aa == "ESV")  
    return 2719;
  else if( aa == "ETA")  
    return 2720;
  else if( aa == "ETR")  
    return 2721;
  else if( aa == "ETN")  
    return 2722;
  else if( aa == "ETD")  
    return 2723;
  else if( aa == "ETC")  
    return 2724;
  else if( aa == "ETQ")  
    return 2725;
  else if( aa == "ETE")  
    return 2726;
  else if( aa == "ETG")  
    return 2727;
  else if( aa == "ETH")  
    return 2728;
  else if( aa == "ETI")  
    return 2729;
  else if( aa == "ETL")  
    return 2730;
  else if( aa == "ETK")  
    return 2731;
  else if( aa == "ETM")  
    return 2732;
  else if( aa == "ETF")  
    return 2733;
  else if( aa == "ETP")  
    return 2734;
  else if( aa == "ETS")  
    return 2735;
  else if( aa == "ETT")  
    return 2736;
  else if( aa == "ETW")  
    return 2737;
  else if( aa == "ETY")  
    return 2738;
  else if( aa == "ETV")  
    return 2739;
  else if( aa == "EWA")  
    return 2740;
  else if( aa == "EWR")  
    return 2741;
  else if( aa == "EWN")  
    return 2742;
  else if( aa == "EWD")  
    return 2743;
  else if( aa == "EWC")  
    return 2744;
  else if( aa == "EWQ")  
    return 2745;
  else if( aa == "EWE")  
    return 2746;
  else if( aa == "EWG")  
    return 2747;
  else if( aa == "EWH")  
    return 2748;
  else if( aa == "EWI")  
    return 2749;
  else if( aa == "EWL")  
    return 2750;
  else if( aa == "EWK")  
    return 2751;
  else if( aa == "EWM")  
    return 2752;
  else if( aa == "EWF")  
    return 2753;
  else if( aa == "EWP")  
    return 2754;
  else if( aa == "EWS")  
    return 2755;
  else if( aa == "EWT")  
    return 2756;
  else if( aa == "EWW")  
    return 2757;
  else if( aa == "EWY")  
    return 2758;
  else if( aa == "EWV")  
    return 2759;
  else if( aa == "EYA")  
    return 2760;
  else if( aa == "EYR")  
    return 2761;
  else if( aa == "EYN")  
    return 2762;
  else if( aa == "EYD")  
    return 2763;
  else if( aa == "EYC")  
    return 2764;
  else if( aa == "EYQ")  
    return 2765;
  else if( aa == "EYE")  
    return 2766;
  else if( aa == "EYG")  
    return 2767;
  else if( aa == "EYH")  
    return 2768;
  else if( aa == "EYI")  
    return 2769;
  else if( aa == "EYL")  
    return 2770;
  else if( aa == "EYK")  
    return 2771;
  else if( aa == "EYM")  
    return 2772;
  else if( aa == "EYF")  
    return 2773;
  else if( aa == "EYP")  
    return 2774;
  else if( aa == "EYS")  
    return 2775;
  else if( aa == "EYT")  
    return 2776;
  else if( aa == "EYW")  
    return 2777;
  else if( aa == "EYY")  
    return 2778;
  else if( aa == "EYV")  
    return 2779;
  else if( aa == "EVA")  
    return 2780;
  else if( aa == "EVR")  
    return 2781;
  else if( aa == "EVN")  
    return 2782;
  else if( aa == "EVD")  
    return 2783;
  else if( aa == "EVC")  
    return 2784;
  else if( aa == "EVQ")  
    return 2785;
  else if( aa == "EVE")  
    return 2786;
  else if( aa == "EVG")  
    return 2787;
  else if( aa == "EVH")  
    return 2788;
  else if( aa == "EVI")  
    return 2789;
  else if( aa == "EVL")  
    return 2790;
  else if( aa == "EVK")  
    return 2791;
  else if( aa == "EVM")  
    return 2792;
  else if( aa == "EVF")  
    return 2793;
  else if( aa == "EVP")  
    return 2794;
  else if( aa == "EVS")  
    return 2795;
  else if( aa == "EVT")  
    return 2796;
  else if( aa == "EVW")  
    return 2797;
  else if( aa == "EVY")  
    return 2798;
  else if( aa == "EVV")  
    return 2799;
  else if( aa == "GAA")  
    return 2800;
  else if( aa == "GAR")  
    return 2801;
  else if( aa == "GAN")  
    return 2802;
  else if( aa == "GAD")  
    return 2803;
  else if( aa == "GAC")  
    return 2804;
  else if( aa == "GAQ")  
    return 2805;
  else if( aa == "GAE")  
    return 2806;
  else if( aa == "GAG")  
    return 2807;
  else if( aa == "GAH")  
    return 2808;
  else if( aa == "GAI")  
    return 2809;
  else if( aa == "GAL")  
    return 2810;
  else if( aa == "GAK")  
    return 2811;
  else if( aa == "GAM")  
    return 2812;
  else if( aa == "GAF")  
    return 2813;
  else if( aa == "GAP")  
    return 2814;
  else if( aa == "GAS")  
    return 2815;
  else if( aa == "GAT")  
    return 2816;
  else if( aa == "GAW")  
    return 2817;
  else if( aa == "GAY")  
    return 2818;
  else if( aa == "GAV")  
    return 2819;
  else if( aa == "GRA")  
    return 2820;
  else if( aa == "GRR")  
    return 2821;
  else if( aa == "GRN")  
    return 2822;
  else if( aa == "GRD")  
    return 2823;
  else if( aa == "GRC")  
    return 2824;
  else if( aa == "GRQ")  
    return 2825;
  else if( aa == "GRE")  
    return 2826;
  else if( aa == "GRG")  
    return 2827;
  else if( aa == "GRH")  
    return 2828;
  else if( aa == "GRI")  
    return 2829;
  else if( aa == "GRL")  
    return 2830;
  else if( aa == "GRK")  
    return 2831;
  else if( aa == "GRM")  
    return 2832;
  else if( aa == "GRF")  
    return 2833;
  else if( aa == "GRP")  
    return 2834;
  else if( aa == "GRS")  
    return 2835;
  else if( aa == "GRT")  
    return 2836;
  else if( aa == "GRW")  
    return 2837;
  else if( aa == "GRY")  
    return 2838;
  else if( aa == "GRV")  
    return 2839;
  else if( aa == "GNA")  
    return 2840;
  else if( aa == "GNR")  
    return 2841;
  else if( aa == "GNN")  
    return 2842;
  else if( aa == "GND")  
    return 2843;
  else if( aa == "GNC")  
    return 2844;
  else if( aa == "GNQ")  
    return 2845;
  else if( aa == "GNE")  
    return 2846;
  else if( aa == "GNG")  
    return 2847;
  else if( aa == "GNH")  
    return 2848;
  else if( aa == "GNI")  
    return 2849;
  else if( aa == "GNL")  
    return 2850;
  else if( aa == "GNK")  
    return 2851;
  else if( aa == "GNM")  
    return 2852;
  else if( aa == "GNF")  
    return 2853;
  else if( aa == "GNP")  
    return 2854;
  else if( aa == "GNS")  
    return 2855;
  else if( aa == "GNT")  
    return 2856;
  else if( aa == "GNW")  
    return 2857;
  else if( aa == "GNY")  
    return 2858;
  else if( aa == "GNV")  
    return 2859;
  else if( aa == "GDA")  
    return 2860;
  else if( aa == "GDR")  
    return 2861;
  else if( aa == "GDN")  
    return 2862;
  else if( aa == "GDD")  
    return 2863;
  else if( aa == "GDC")  
    return 2864;
  else if( aa == "GDQ")  
    return 2865;
  else if( aa == "GDE")  
    return 2866;
  else if( aa == "GDG")  
    return 2867;
  else if( aa == "GDH")  
    return 2868;
  else if( aa == "GDI")  
    return 2869;
  else if( aa == "GDL")  
    return 2870;
  else if( aa == "GDK")  
    return 2871;
  else if( aa == "GDM")  
    return 2872;
  else if( aa == "GDF")  
    return 2873;
  else if( aa == "GDP")  
    return 2874;
  else if( aa == "GDS")  
    return 2875;
  else if( aa == "GDT")  
    return 2876;
  else if( aa == "GDW")  
    return 2877;
  else if( aa == "GDY")  
    return 2878;
  else if( aa == "GDV")  
    return 2879;
  else if( aa == "GCA")  
    return 2880;
  else if( aa == "GCR")  
    return 2881;
  else if( aa == "GCN")  
    return 2882;
  else if( aa == "GCD")  
    return 2883;
  else if( aa == "GCC")  
    return 2884;
  else if( aa == "GCQ")  
    return 2885;
  else if( aa == "GCE")  
    return 2886;
  else if( aa == "GCG")  
    return 2887;
  else if( aa == "GCH")  
    return 2888;
  else if( aa == "GCI")  
    return 2889;
  else if( aa == "GCL")  
    return 2890;
  else if( aa == "GCK")  
    return 2891;
  else if( aa == "GCM")  
    return 2892;
  else if( aa == "GCF")  
    return 2893;
  else if( aa == "GCP")  
    return 2894;
  else if( aa == "GCS")  
    return 2895;
  else if( aa == "GCT")  
    return 2896;
  else if( aa == "GCW")  
    return 2897;
  else if( aa == "GCY")  
    return 2898;
  else if( aa == "GCV")  
    return 2899;
  else if( aa == "GQA")  
    return 2900;
  else if( aa == "GQR")  
    return 2901;
  else if( aa == "GQN")  
    return 2902;
  else if( aa == "GQD")  
    return 2903;
  else if( aa == "GQC")  
    return 2904;
  else if( aa == "GQQ")  
    return 2905;
  else if( aa == "GQE")  
    return 2906;
  else if( aa == "GQG")  
    return 2907;
  else if( aa == "GQH")  
    return 2908;
  else if( aa == "GQI")  
    return 2909;
  else if( aa == "GQL")  
    return 2910;
  else if( aa == "GQK")  
    return 2911;
  else if( aa == "GQM")  
    return 2912;
  else if( aa == "GQF")  
    return 2913;
  else if( aa == "GQP")  
    return 2914;
  else if( aa == "GQS")  
    return 2915;
  else if( aa == "GQT")  
    return 2916;
  else if( aa == "GQW")  
    return 2917;
  else if( aa == "GQY")  
    return 2918;
  else if( aa == "GQV")  
    return 2919;
  else if( aa == "GEA")  
    return 2920;
  else if( aa == "GER")  
    return 2921;
  else if( aa == "GEN")  
    return 2922;
  else if( aa == "GED")  
    return 2923;
  else if( aa == "GEC")  
    return 2924;
  else if( aa == "GEQ")  
    return 2925;
  else if( aa == "GEE")  
    return 2926;
  else if( aa == "GEG")  
    return 2927;
  else if( aa == "GEH")  
    return 2928;
  else if( aa == "GEI")  
    return 2929;
  else if( aa == "GEL")  
    return 2930;
  else if( aa == "GEK")  
    return 2931;
  else if( aa == "GEM")  
    return 2932;
  else if( aa == "GEF")  
    return 2933;
  else if( aa == "GEP")  
    return 2934;
  else if( aa == "GES")  
    return 2935;
  else if( aa == "GET")  
    return 2936;
  else if( aa == "GEW")  
    return 2937;
  else if( aa == "GEY")  
    return 2938;
  else if( aa == "GEV")  
    return 2939;
  else if( aa == "GGA")  
    return 2940;
  else if( aa == "GGR")  
    return 2941;
  else if( aa == "GGN")  
    return 2942;
  else if( aa == "GGD")  
    return 2943;
  else if( aa == "GGC")  
    return 2944;
  else if( aa == "GGQ")  
    return 2945;
  else if( aa == "GGE")  
    return 2946;
  else if( aa == "GGG")  
    return 2947;
  else if( aa == "GGH")  
    return 2948;
  else if( aa == "GGI")  
    return 2949;
  else if( aa == "GGL")  
    return 2950;
  else if( aa == "GGK")  
    return 2951;
  else if( aa == "GGM")  
    return 2952;
  else if( aa == "GGF")  
    return 2953;
  else if( aa == "GGP")  
    return 2954;
  else if( aa == "GGS")  
    return 2955;
  else if( aa == "GGT")  
    return 2956;
  else if( aa == "GGW")  
    return 2957;
  else if( aa == "GGY")  
    return 2958;
  else if( aa == "GGV")  
    return 2959;
  else if( aa == "GHA")  
    return 2960;
  else if( aa == "GHR")  
    return 2961;
  else if( aa == "GHN")  
    return 2962;
  else if( aa == "GHD")  
    return 2963;
  else if( aa == "GHC")  
    return 2964;
  else if( aa == "GHQ")  
    return 2965;
  else if( aa == "GHE")  
    return 2966;
  else if( aa == "GHG")  
    return 2967;
  else if( aa == "GHH")  
    return 2968;
  else if( aa == "GHI")  
    return 2969;
  else if( aa == "GHL")  
    return 2970;
  else if( aa == "GHK")  
    return 2971;
  else if( aa == "GHM")  
    return 2972;
  else if( aa == "GHF")  
    return 2973;
  else if( aa == "GHP")  
    return 2974;
  else if( aa == "GHS")  
    return 2975;
  else if( aa == "GHT")  
    return 2976;
  else if( aa == "GHW")  
    return 2977;
  else if( aa == "GHY")  
    return 2978;
  else if( aa == "GHV")  
    return 2979;
  else if( aa == "GIA")  
    return 2980;
  else if( aa == "GIR")  
    return 2981;
  else if( aa == "GIN")  
    return 2982;
  else if( aa == "GID")  
    return 2983;
  else if( aa == "GIC")  
    return 2984;
  else if( aa == "GIQ")  
    return 2985;
  else if( aa == "GIE")  
    return 2986;
  else if( aa == "GIG")  
    return 2987;
  else if( aa == "GIH")  
    return 2988;
  else if( aa == "GII")  
    return 2989;
  else if( aa == "GIL")  
    return 2990;
  else if( aa == "GIK")  
    return 2991;
  else if( aa == "GIM")  
    return 2992;
  else if( aa == "GIF")  
    return 2993;
  else if( aa == "GIP")  
    return 2994;
  else if( aa == "GIS")  
    return 2995;
  else if( aa == "GIT")  
    return 2996;
  else if( aa == "GIW")  
    return 2997;
  else if( aa == "GIY")  
    return 2998;
  else if( aa == "GIV")  
    return 2999;
  else if( aa == "GLA")  
    return 3000;
  else if( aa == "GLR")  
    return 3001;
  else if( aa == "GLN")  
    return 3002;
  else if( aa == "GLD")  
    return 3003;
  else if( aa == "GLC")  
    return 3004;
  else if( aa == "GLQ")  
    return 3005;
  else if( aa == "GLE")  
    return 3006;
  else if( aa == "GLG")  
    return 3007;
  else if( aa == "GLH")  
    return 3008;
  else if( aa == "GLI")  
    return 3009;
  else if( aa == "GLL")  
    return 3010;
  else if( aa == "GLK")  
    return 3011;
  else if( aa == "GLM")  
    return 3012;
  else if( aa == "GLF")  
    return 3013;
  else if( aa == "GLP")  
    return 3014;
  else if( aa == "GLS")  
    return 3015;
  else if( aa == "GLT")  
    return 3016;
  else if( aa == "GLW")  
    return 3017;
  else if( aa == "GLY")  
    return 3018;
  else if( aa == "GLV")  
    return 3019;
  else if( aa == "GKA")  
    return 3020;
  else if( aa == "GKR")  
    return 3021;
  else if( aa == "GKN")  
    return 3022;
  else if( aa == "GKD")  
    return 3023;
  else if( aa == "GKC")  
    return 3024;
  else if( aa == "GKQ")  
    return 3025;
  else if( aa == "GKE")  
    return 3026;
  else if( aa == "GKG")  
    return 3027;
  else if( aa == "GKH")  
    return 3028;
  else if( aa == "GKI")  
    return 3029;
  else if( aa == "GKL")  
    return 3030;
  else if( aa == "GKK")  
    return 3031;
  else if( aa == "GKM")  
    return 3032;
  else if( aa == "GKF")  
    return 3033;
  else if( aa == "GKP")  
    return 3034;
  else if( aa == "GKS")  
    return 3035;
  else if( aa == "GKT")  
    return 3036;
  else if( aa == "GKW")  
    return 3037;
  else if( aa == "GKY")  
    return 3038;
  else if( aa == "GKV")  
    return 3039;
  else if( aa == "GMA")  
    return 3040;
  else if( aa == "GMR")  
    return 3041;
  else if( aa == "GMN")  
    return 3042;
  else if( aa == "GMD")  
    return 3043;
  else if( aa == "GMC")  
    return 3044;
  else if( aa == "GMQ")  
    return 3045;
  else if( aa == "GME")  
    return 3046;
  else if( aa == "GMG")  
    return 3047;
  else if( aa == "GMH")  
    return 3048;
  else if( aa == "GMI")  
    return 3049;
  else if( aa == "GML")  
    return 3050;
  else if( aa == "GMK")  
    return 3051;
  else if( aa == "GMM")  
    return 3052;
  else if( aa == "GMF")  
    return 3053;
  else if( aa == "GMP")  
    return 3054;
  else if( aa == "GMS")  
    return 3055;
  else if( aa == "GMT")  
    return 3056;
  else if( aa == "GMW")  
    return 3057;
  else if( aa == "GMY")  
    return 3058;
  else if( aa == "GMV")  
    return 3059;
  else if( aa == "GFA")  
    return 3060;
  else if( aa == "GFR")  
    return 3061;
  else if( aa == "GFN")  
    return 3062;
  else if( aa == "GFD")  
    return 3063;
  else if( aa == "GFC")  
    return 3064;
  else if( aa == "GFQ")  
    return 3065;
  else if( aa == "GFE")  
    return 3066;
  else if( aa == "GFG")  
    return 3067;
  else if( aa == "GFH")  
    return 3068;
  else if( aa == "GFI")  
    return 3069;
  else if( aa == "GFL")  
    return 3070;
  else if( aa == "GFK")  
    return 3071;
  else if( aa == "GFM")  
    return 3072;
  else if( aa == "GFF")  
    return 3073;
  else if( aa == "GFP")  
    return 3074;
  else if( aa == "GFS")  
    return 3075;
  else if( aa == "GFT")  
    return 3076;
  else if( aa == "GFW")  
    return 3077;
  else if( aa == "GFY")  
    return 3078;
  else if( aa == "GFV")  
    return 3079;
  else if( aa == "GPA")  
    return 3080;
  else if( aa == "GPR")  
    return 3081;
  else if( aa == "GPN")  
    return 3082;
  else if( aa == "GPD")  
    return 3083;
  else if( aa == "GPC")  
    return 3084;
  else if( aa == "GPQ")  
    return 3085;
  else if( aa == "GPE")  
    return 3086;
  else if( aa == "GPG")  
    return 3087;
  else if( aa == "GPH")  
    return 3088;
  else if( aa == "GPI")  
    return 3089;
  else if( aa == "GPL")  
    return 3090;
  else if( aa == "GPK")  
    return 3091;
  else if( aa == "GPM")  
    return 3092;
  else if( aa == "GPF")  
    return 3093;
  else if( aa == "GPP")  
    return 3094;
  else if( aa == "GPS")  
    return 3095;
  else if( aa == "GPT")  
    return 3096;
  else if( aa == "GPW")  
    return 3097;
  else if( aa == "GPY")  
    return 3098;
  else if( aa == "GPV")  
    return 3099;
  else if( aa == "GSA")  
    return 3100;
  else if( aa == "GSR")  
    return 3101;
  else if( aa == "GSN")  
    return 3102;
  else if( aa == "GSD")  
    return 3103;
  else if( aa == "GSC")  
    return 3104;
  else if( aa == "GSQ")  
    return 3105;
  else if( aa == "GSE")  
    return 3106;
  else if( aa == "GSG")  
    return 3107;
  else if( aa == "GSH")  
    return 3108;
  else if( aa == "GSI")  
    return 3109;
  else if( aa == "GSL")  
    return 3110;
  else if( aa == "GSK")  
    return 3111;
  else if( aa == "GSM")  
    return 3112;
  else if( aa == "GSF")  
    return 3113;
  else if( aa == "GSP")  
    return 3114;
  else if( aa == "GSS")  
    return 3115;
  else if( aa == "GST")  
    return 3116;
  else if( aa == "GSW")  
    return 3117;
  else if( aa == "GSY")  
    return 3118;
  else if( aa == "GSV")  
    return 3119;
  else if( aa == "GTA")  
    return 3120;
  else if( aa == "GTR")  
    return 3121;
  else if( aa == "GTN")  
    return 3122;
  else if( aa == "GTD")  
    return 3123;
  else if( aa == "GTC")  
    return 3124;
  else if( aa == "GTQ")  
    return 3125;
  else if( aa == "GTE")  
    return 3126;
  else if( aa == "GTG")  
    return 3127;
  else if( aa == "GTH")  
    return 3128;
  else if( aa == "GTI")  
    return 3129;
  else if( aa == "GTL")  
    return 3130;
  else if( aa == "GTK")  
    return 3131;
  else if( aa == "GTM")  
    return 3132;
  else if( aa == "GTF")  
    return 3133;
  else if( aa == "GTP")  
    return 3134;
  else if( aa == "GTS")  
    return 3135;
  else if( aa == "GTT")  
    return 3136;
  else if( aa == "GTW")  
    return 3137;
  else if( aa == "GTY")  
    return 3138;
  else if( aa == "GTV")  
    return 3139;
  else if( aa == "GWA")  
    return 3140;
  else if( aa == "GWR")  
    return 3141;
  else if( aa == "GWN")  
    return 3142;
  else if( aa == "GWD")  
    return 3143;
  else if( aa == "GWC")  
    return 3144;
  else if( aa == "GWQ")  
    return 3145;
  else if( aa == "GWE")  
    return 3146;
  else if( aa == "GWG")  
    return 3147;
  else if( aa == "GWH")  
    return 3148;
  else if( aa == "GWI")  
    return 3149;
  else if( aa == "GWL")  
    return 3150;
  else if( aa == "GWK")  
    return 3151;
  else if( aa == "GWM")  
    return 3152;
  else if( aa == "GWF")  
    return 3153;
  else if( aa == "GWP")  
    return 3154;
  else if( aa == "GWS")  
    return 3155;
  else if( aa == "GWT")  
    return 3156;
  else if( aa == "GWW")  
    return 3157;
  else if( aa == "GWY")  
    return 3158;
  else if( aa == "GWV")  
    return 3159;
  else if( aa == "GYA")  
    return 3160;
  else if( aa == "GYR")  
    return 3161;
  else if( aa == "GYN")  
    return 3162;
  else if( aa == "GYD")  
    return 3163;
  else if( aa == "GYC")  
    return 3164;
  else if( aa == "GYQ")  
    return 3165;
  else if( aa == "GYE")  
    return 3166;
  else if( aa == "GYG")  
    return 3167;
  else if( aa == "GYH")  
    return 3168;
  else if( aa == "GYI")  
    return 3169;
  else if( aa == "GYL")  
    return 3170;
  else if( aa == "GYK")  
    return 3171;
  else if( aa == "GYM")  
    return 3172;
  else if( aa == "GYF")  
    return 3173;
  else if( aa == "GYP")  
    return 3174;
  else if( aa == "GYS")  
    return 3175;
  else if( aa == "GYT")  
    return 3176;
  else if( aa == "GYW")  
    return 3177;
  else if( aa == "GYY")  
    return 3178;
  else if( aa == "GYV")  
    return 3179;
  else if( aa == "GVA")  
    return 3180;
  else if( aa == "GVR")  
    return 3181;
  else if( aa == "GVN")  
    return 3182;
  else if( aa == "GVD")  
    return 3183;
  else if( aa == "GVC")  
    return 3184;
  else if( aa == "GVQ")  
    return 3185;
  else if( aa == "GVE")  
    return 3186;
  else if( aa == "GVG")  
    return 3187;
  else if( aa == "GVH")  
    return 3188;
  else if( aa == "GVI")  
    return 3189;
  else if( aa == "GVL")  
    return 3190;
  else if( aa == "GVK")  
    return 3191;
  else if( aa == "GVM")  
    return 3192;
  else if( aa == "GVF")  
    return 3193;
  else if( aa == "GVP")  
    return 3194;
  else if( aa == "GVS")  
    return 3195;
  else if( aa == "GVT")  
    return 3196;
  else if( aa == "GVW")  
    return 3197;
  else if( aa == "GVY")  
    return 3198;
  else if( aa == "GVV")  
    return 3199;
  else if( aa == "HAA")  
    return 3200;
  else if( aa == "HAR")  
    return 3201;
  else if( aa == "HAN")  
    return 3202;
  else if( aa == "HAD")  
    return 3203;
  else if( aa == "HAC")  
    return 3204;
  else if( aa == "HAQ")  
    return 3205;
  else if( aa == "HAE")  
    return 3206;
  else if( aa == "HAG")  
    return 3207;
  else if( aa == "HAH")  
    return 3208;
  else if( aa == "HAI")  
    return 3209;
  else if( aa == "HAL")  
    return 3210;
  else if( aa == "HAK")  
    return 3211;
  else if( aa == "HAM")  
    return 3212;
  else if( aa == "HAF")  
    return 3213;
  else if( aa == "HAP")  
    return 3214;
  else if( aa == "HAS")  
    return 3215;
  else if( aa == "HAT")  
    return 3216;
  else if( aa == "HAW")  
    return 3217;
  else if( aa == "HAY")  
    return 3218;
  else if( aa == "HAV")  
    return 3219;
  else if( aa == "HRA")  
    return 3220;
  else if( aa == "HRR")  
    return 3221;
  else if( aa == "HRN")  
    return 3222;
  else if( aa == "HRD")  
    return 3223;
  else if( aa == "HRC")  
    return 3224;
  else if( aa == "HRQ")  
    return 3225;
  else if( aa == "HRE")  
    return 3226;
  else if( aa == "HRG")  
    return 3227;
  else if( aa == "HRH")  
    return 3228;
  else if( aa == "HRI")  
    return 3229;
  else if( aa == "HRL")  
    return 3230;
  else if( aa == "HRK")  
    return 3231;
  else if( aa == "HRM")  
    return 3232;
  else if( aa == "HRF")  
    return 3233;
  else if( aa == "HRP")  
    return 3234;
  else if( aa == "HRS")  
    return 3235;
  else if( aa == "HRT")  
    return 3236;
  else if( aa == "HRW")  
    return 3237;
  else if( aa == "HRY")  
    return 3238;
  else if( aa == "HRV")  
    return 3239;
  else if( aa == "HNA")  
    return 3240;
  else if( aa == "HNR")  
    return 3241;
  else if( aa == "HNN")  
    return 3242;
  else if( aa == "HND")  
    return 3243;
  else if( aa == "HNC")  
    return 3244;
  else if( aa == "HNQ")  
    return 3245;
  else if( aa == "HNE")  
    return 3246;
  else if( aa == "HNG")  
    return 3247;
  else if( aa == "HNH")  
    return 3248;
  else if( aa == "HNI")  
    return 3249;
  else if( aa == "HNL")  
    return 3250;
  else if( aa == "HNK")  
    return 3251;
  else if( aa == "HNM")  
    return 3252;
  else if( aa == "HNF")  
    return 3253;
  else if( aa == "HNP")  
    return 3254;
  else if( aa == "HNS")  
    return 3255;
  else if( aa == "HNT")  
    return 3256;
  else if( aa == "HNW")  
    return 3257;
  else if( aa == "HNY")  
    return 3258;
  else if( aa == "HNV")  
    return 3259;
  else if( aa == "HDA")  
    return 3260;
  else if( aa == "HDR")  
    return 3261;
  else if( aa == "HDN")  
    return 3262;
  else if( aa == "HDD")  
    return 3263;
  else if( aa == "HDC")  
    return 3264;
  else if( aa == "HDQ")  
    return 3265;
  else if( aa == "HDE")  
    return 3266;
  else if( aa == "HDG")  
    return 3267;
  else if( aa == "HDH")  
    return 3268;
  else if( aa == "HDI")  
    return 3269;
  else if( aa == "HDL")  
    return 3270;
  else if( aa == "HDK")  
    return 3271;
  else if( aa == "HDM")  
    return 3272;
  else if( aa == "HDF")  
    return 3273;
  else if( aa == "HDP")  
    return 3274;
  else if( aa == "HDS")  
    return 3275;
  else if( aa == "HDT")  
    return 3276;
  else if( aa == "HDW")  
    return 3277;
  else if( aa == "HDY")  
    return 3278;
  else if( aa == "HDV")  
    return 3279;
  else if( aa == "HCA")  
    return 3280;
  else if( aa == "HCR")  
    return 3281;
  else if( aa == "HCN")  
    return 3282;
  else if( aa == "HCD")  
    return 3283;
  else if( aa == "HCC")  
    return 3284;
  else if( aa == "HCQ")  
    return 3285;
  else if( aa == "HCE")  
    return 3286;
  else if( aa == "HCG")  
    return 3287;
  else if( aa == "HCH")  
    return 3288;
  else if( aa == "HCI")  
    return 3289;
  else if( aa == "HCL")  
    return 3290;
  else if( aa == "HCK")  
    return 3291;
  else if( aa == "HCM")  
    return 3292;
  else if( aa == "HCF")  
    return 3293;
  else if( aa == "HCP")  
    return 3294;
  else if( aa == "HCS")  
    return 3295;
  else if( aa == "HCT")  
    return 3296;
  else if( aa == "HCW")  
    return 3297;
  else if( aa == "HCY")  
    return 3298;
  else if( aa == "HCV")  
    return 3299;
  else if( aa == "HQA")  
    return 3300;
  else if( aa == "HQR")  
    return 3301;
  else if( aa == "HQN")  
    return 3302;
  else if( aa == "HQD")  
    return 3303;
  else if( aa == "HQC")  
    return 3304;
  else if( aa == "HQQ")  
    return 3305;
  else if( aa == "HQE")  
    return 3306;
  else if( aa == "HQG")  
    return 3307;
  else if( aa == "HQH")  
    return 3308;
  else if( aa == "HQI")  
    return 3309;
  else if( aa == "HQL")  
    return 3310;
  else if( aa == "HQK")  
    return 3311;
  else if( aa == "HQM")  
    return 3312;
  else if( aa == "HQF")  
    return 3313;
  else if( aa == "HQP")  
    return 3314;
  else if( aa == "HQS")  
    return 3315;
  else if( aa == "HQT")  
    return 3316;
  else if( aa == "HQW")  
    return 3317;
  else if( aa == "HQY")  
    return 3318;
  else if( aa == "HQV")  
    return 3319;
  else if( aa == "HEA")  
    return 3320;
  else if( aa == "HER")  
    return 3321;
  else if( aa == "HEN")  
    return 3322;
  else if( aa == "HED")  
    return 3323;
  else if( aa == "HEC")  
    return 3324;
  else if( aa == "HEQ")  
    return 3325;
  else if( aa == "HEE")  
    return 3326;
  else if( aa == "HEG")  
    return 3327;
  else if( aa == "HEH")  
    return 3328;
  else if( aa == "HEI")  
    return 3329;
  else if( aa == "HEL")  
    return 3330;
  else if( aa == "HEK")  
    return 3331;
  else if( aa == "HEM")  
    return 3332;
  else if( aa == "HEF")  
    return 3333;
  else if( aa == "HEP")  
    return 3334;
  else if( aa == "HES")  
    return 3335;
  else if( aa == "HET")  
    return 3336;
  else if( aa == "HEW")  
    return 3337;
  else if( aa == "HEY")  
    return 3338;
  else if( aa == "HEV")  
    return 3339;
  else if( aa == "HGA")  
    return 3340;
  else if( aa == "HGR")  
    return 3341;
  else if( aa == "HGN")  
    return 3342;
  else if( aa == "HGD")  
    return 3343;
  else if( aa == "HGC")  
    return 3344;
  else if( aa == "HGQ")  
    return 3345;
  else if( aa == "HGE")  
    return 3346;
  else if( aa == "HGG")  
    return 3347;
  else if( aa == "HGH")  
    return 3348;
  else if( aa == "HGI")  
    return 3349;
  else if( aa == "HGL")  
    return 3350;
  else if( aa == "HGK")  
    return 3351;
  else if( aa == "HGM")  
    return 3352;
  else if( aa == "HGF")  
    return 3353;
  else if( aa == "HGP")  
    return 3354;
  else if( aa == "HGS")  
    return 3355;
  else if( aa == "HGT")  
    return 3356;
  else if( aa == "HGW")  
    return 3357;
  else if( aa == "HGY")  
    return 3358;
  else if( aa == "HGV")  
    return 3359;
  else if( aa == "HHA")  
    return 3360;
  else if( aa == "HHR")  
    return 3361;
  else if( aa == "HHN")  
    return 3362;
  else if( aa == "HHD")  
    return 3363;
  else if( aa == "HHC")  
    return 3364;
  else if( aa == "HHQ")  
    return 3365;
  else if( aa == "HHE")  
    return 3366;
  else if( aa == "HHG")  
    return 3367;
  else if( aa == "HHH")  
    return 3368;
  else if( aa == "HHI")  
    return 3369;
  else if( aa == "HHL")  
    return 3370;
  else if( aa == "HHK")  
    return 3371;
  else if( aa == "HHM")  
    return 3372;
  else if( aa == "HHF")  
    return 3373;
  else if( aa == "HHP")  
    return 3374;
  else if( aa == "HHS")  
    return 3375;
  else if( aa == "HHT")  
    return 3376;
  else if( aa == "HHW")  
    return 3377;
  else if( aa == "HHY")  
    return 3378;
  else if( aa == "HHV")  
    return 3379;
  else if( aa == "HIA")  
    return 3380;
  else if( aa == "HIR")  
    return 3381;
  else if( aa == "HIN")  
    return 3382;
  else if( aa == "HID")  
    return 3383;
  else if( aa == "HIC")  
    return 3384;
  else if( aa == "HIQ")  
    return 3385;
  else if( aa == "HIE")  
    return 3386;
  else if( aa == "HIG")  
    return 3387;
  else if( aa == "HIH")  
    return 3388;
  else if( aa == "HII")  
    return 3389;
  else if( aa == "HIL")  
    return 3390;
  else if( aa == "HIK")  
    return 3391;
  else if( aa == "HIM")  
    return 3392;
  else if( aa == "HIF")  
    return 3393;
  else if( aa == "HIP")  
    return 3394;
  else if( aa == "HIS")  
    return 3395;
  else if( aa == "HIT")  
    return 3396;
  else if( aa == "HIW")  
    return 3397;
  else if( aa == "HIY")  
    return 3398;
  else if( aa == "HIV")  
    return 3399;
  else if( aa == "HLA")  
    return 3400;
  else if( aa == "HLR")  
    return 3401;
  else if( aa == "HLN")  
    return 3402;
  else if( aa == "HLD")  
    return 3403;
  else if( aa == "HLC")  
    return 3404;
  else if( aa == "HLQ")  
    return 3405;
  else if( aa == "HLE")  
    return 3406;
  else if( aa == "HLG")  
    return 3407;
  else if( aa == "HLH")  
    return 3408;
  else if( aa == "HLI")  
    return 3409;
  else if( aa == "HLL")  
    return 3410;
  else if( aa == "HLK")  
    return 3411;
  else if( aa == "HLM")  
    return 3412;
  else if( aa == "HLF")  
    return 3413;
  else if( aa == "HLP")  
    return 3414;
  else if( aa == "HLS")  
    return 3415;
  else if( aa == "HLT")  
    return 3416;
  else if( aa == "HLW")  
    return 3417;
  else if( aa == "HLY")  
    return 3418;
  else if( aa == "HLV")  
    return 3419;
  else if( aa == "HKA")  
    return 3420;
  else if( aa == "HKR")  
    return 3421;
  else if( aa == "HKN")  
    return 3422;
  else if( aa == "HKD")  
    return 3423;
  else if( aa == "HKC")  
    return 3424;
  else if( aa == "HKQ")  
    return 3425;
  else if( aa == "HKE")  
    return 3426;
  else if( aa == "HKG")  
    return 3427;
  else if( aa == "HKH")  
    return 3428;
  else if( aa == "HKI")  
    return 3429;
  else if( aa == "HKL")  
    return 3430;
  else if( aa == "HKK")  
    return 3431;
  else if( aa == "HKM")  
    return 3432;
  else if( aa == "HKF")  
    return 3433;
  else if( aa == "HKP")  
    return 3434;
  else if( aa == "HKS")  
    return 3435;
  else if( aa == "HKT")  
    return 3436;
  else if( aa == "HKW")  
    return 3437;
  else if( aa == "HKY")  
    return 3438;
  else if( aa == "HKV")  
    return 3439;
  else if( aa == "HMA")  
    return 3440;
  else if( aa == "HMR")  
    return 3441;
  else if( aa == "HMN")  
    return 3442;
  else if( aa == "HMD")  
    return 3443;
  else if( aa == "HMC")  
    return 3444;
  else if( aa == "HMQ")  
    return 3445;
  else if( aa == "HME")  
    return 3446;
  else if( aa == "HMG")  
    return 3447;
  else if( aa == "HMH")  
    return 3448;
  else if( aa == "HMI")  
    return 3449;
  else if( aa == "HML")  
    return 3450;
  else if( aa == "HMK")  
    return 3451;
  else if( aa == "HMM")  
    return 3452;
  else if( aa == "HMF")  
    return 3453;
  else if( aa == "HMP")  
    return 3454;
  else if( aa == "HMS")  
    return 3455;
  else if( aa == "HMT")  
    return 3456;
  else if( aa == "HMW")  
    return 3457;
  else if( aa == "HMY")  
    return 3458;
  else if( aa == "HMV")  
    return 3459;
  else if( aa == "HFA")  
    return 3460;
  else if( aa == "HFR")  
    return 3461;
  else if( aa == "HFN")  
    return 3462;
  else if( aa == "HFD")  
    return 3463;
  else if( aa == "HFC")  
    return 3464;
  else if( aa == "HFQ")  
    return 3465;
  else if( aa == "HFE")  
    return 3466;
  else if( aa == "HFG")  
    return 3467;
  else if( aa == "HFH")  
    return 3468;
  else if( aa == "HFI")  
    return 3469;
  else if( aa == "HFL")  
    return 3470;
  else if( aa == "HFK")  
    return 3471;
  else if( aa == "HFM")  
    return 3472;
  else if( aa == "HFF")  
    return 3473;
  else if( aa == "HFP")  
    return 3474;
  else if( aa == "HFS")  
    return 3475;
  else if( aa == "HFT")  
    return 3476;
  else if( aa == "HFW")  
    return 3477;
  else if( aa == "HFY")  
    return 3478;
  else if( aa == "HFV")  
    return 3479;
  else if( aa == "HPA")  
    return 3480;
  else if( aa == "HPR")  
    return 3481;
  else if( aa == "HPN")  
    return 3482;
  else if( aa == "HPD")  
    return 3483;
  else if( aa == "HPC")  
    return 3484;
  else if( aa == "HPQ")  
    return 3485;
  else if( aa == "HPE")  
    return 3486;
  else if( aa == "HPG")  
    return 3487;
  else if( aa == "HPH")  
    return 3488;
  else if( aa == "HPI")  
    return 3489;
  else if( aa == "HPL")  
    return 3490;
  else if( aa == "HPK")  
    return 3491;
  else if( aa == "HPM")  
    return 3492;
  else if( aa == "HPF")  
    return 3493;
  else if( aa == "HPP")  
    return 3494;
  else if( aa == "HPS")  
    return 3495;
  else if( aa == "HPT")  
    return 3496;
  else if( aa == "HPW")  
    return 3497;
  else if( aa == "HPY")  
    return 3498;
  else if( aa == "HPV")  
    return 3499;
  else if( aa == "HSA")  
    return 3500;
  else if( aa == "HSR")  
    return 3501;
  else if( aa == "HSN")  
    return 3502;
  else if( aa == "HSD")  
    return 3503;
  else if( aa == "HSC")  
    return 3504;
  else if( aa == "HSQ")  
    return 3505;
  else if( aa == "HSE")  
    return 3506;
  else if( aa == "HSG")  
    return 3507;
  else if( aa == "HSH")  
    return 3508;
  else if( aa == "HSI")  
    return 3509;
  else if( aa == "HSL")  
    return 3510;
  else if( aa == "HSK")  
    return 3511;
  else if( aa == "HSM")  
    return 3512;
  else if( aa == "HSF")  
    return 3513;
  else if( aa == "HSP")  
    return 3514;
  else if( aa == "HSS")  
    return 3515;
  else if( aa == "HST")  
    return 3516;
  else if( aa == "HSW")  
    return 3517;
  else if( aa == "HSY")  
    return 3518;
  else if( aa == "HSV")  
    return 3519;
  else if( aa == "HTA")  
    return 3520;
  else if( aa == "HTR")  
    return 3521;
  else if( aa == "HTN")  
    return 3522;
  else if( aa == "HTD")  
    return 3523;
  else if( aa == "HTC")  
    return 3524;
  else if( aa == "HTQ")  
    return 3525;
  else if( aa == "HTE")  
    return 3526;
  else if( aa == "HTG")  
    return 3527;
  else if( aa == "HTH")  
    return 3528;
  else if( aa == "HTI")  
    return 3529;
  else if( aa == "HTL")  
    return 3530;
  else if( aa == "HTK")  
    return 3531;
  else if( aa == "HTM")  
    return 3532;
  else if( aa == "HTF")  
    return 3533;
  else if( aa == "HTP")  
    return 3534;
  else if( aa == "HTS")  
    return 3535;
  else if( aa == "HTT")  
    return 3536;
  else if( aa == "HTW")  
    return 3537;
  else if( aa == "HTY")  
    return 3538;
  else if( aa == "HTV")  
    return 3539;
  else if( aa == "HWA")  
    return 3540;
  else if( aa == "HWR")  
    return 3541;
  else if( aa == "HWN")  
    return 3542;
  else if( aa == "HWD")  
    return 3543;
  else if( aa == "HWC")  
    return 3544;
  else if( aa == "HWQ")  
    return 3545;
  else if( aa == "HWE")  
    return 3546;
  else if( aa == "HWG")  
    return 3547;
  else if( aa == "HWH")  
    return 3548;
  else if( aa == "HWI")  
    return 3549;
  else if( aa == "HWL")  
    return 3550;
  else if( aa == "HWK")  
    return 3551;
  else if( aa == "HWM")  
    return 3552;
  else if( aa == "HWF")  
    return 3553;
  else if( aa == "HWP")  
    return 3554;
  else if( aa == "HWS")  
    return 3555;
  else if( aa == "HWT")  
    return 3556;
  else if( aa == "HWW")  
    return 3557;
  else if( aa == "HWY")  
    return 3558;
  else if( aa == "HWV")  
    return 3559;
  else if( aa == "HYA")  
    return 3560;
  else if( aa == "HYR")  
    return 3561;
  else if( aa == "HYN")  
    return 3562;
  else if( aa == "HYD")  
    return 3563;
  else if( aa == "HYC")  
    return 3564;
  else if( aa == "HYQ")  
    return 3565;
  else if( aa == "HYE")  
    return 3566;
  else if( aa == "HYG")  
    return 3567;
  else if( aa == "HYH")  
    return 3568;
  else if( aa == "HYI")  
    return 3569;
  else if( aa == "HYL")  
    return 3570;
  else if( aa == "HYK")  
    return 3571;
  else if( aa == "HYM")  
    return 3572;
  else if( aa == "HYF")  
    return 3573;
  else if( aa == "HYP")  
    return 3574;
  else if( aa == "HYS")  
    return 3575;
  else if( aa == "HYT")  
    return 3576;
  else if( aa == "HYW")  
    return 3577;
  else if( aa == "HYY")  
    return 3578;
  else if( aa == "HYV")  
    return 3579;
  else if( aa == "HVA")  
    return 3580;
  else if( aa == "HVR")  
    return 3581;
  else if( aa == "HVN")  
    return 3582;
  else if( aa == "HVD")  
    return 3583;
  else if( aa == "HVC")  
    return 3584;
  else if( aa == "HVQ")  
    return 3585;
  else if( aa == "HVE")  
    return 3586;
  else if( aa == "HVG")  
    return 3587;
  else if( aa == "HVH")  
    return 3588;
  else if( aa == "HVI")  
    return 3589;
  else if( aa == "HVL")  
    return 3590;
  else if( aa == "HVK")  
    return 3591;
  else if( aa == "HVM")  
    return 3592;
  else if( aa == "HVF")  
    return 3593;
  else if( aa == "HVP")  
    return 3594;
  else if( aa == "HVS")  
    return 3595;
  else if( aa == "HVT")  
    return 3596;
  else if( aa == "HVW")  
    return 3597;
  else if( aa == "HVY")  
    return 3598;
  else if( aa == "HVV")  
    return 3599;
  else if( aa == "IAA")  
    return 3600;
  else if( aa == "IAR")  
    return 3601;
  else if( aa == "IAN")  
    return 3602;
  else if( aa == "IAD")  
    return 3603;
  else if( aa == "IAC")  
    return 3604;
  else if( aa == "IAQ")  
    return 3605;
  else if( aa == "IAE")  
    return 3606;
  else if( aa == "IAG")  
    return 3607;
  else if( aa == "IAH")  
    return 3608;
  else if( aa == "IAI")  
    return 3609;
  else if( aa == "IAL")  
    return 3610;
  else if( aa == "IAK")  
    return 3611;
  else if( aa == "IAM")  
    return 3612;
  else if( aa == "IAF")  
    return 3613;
  else if( aa == "IAP")  
    return 3614;
  else if( aa == "IAS")  
    return 3615;
  else if( aa == "IAT")  
    return 3616;
  else if( aa == "IAW")  
    return 3617;
  else if( aa == "IAY")  
    return 3618;
  else if( aa == "IAV")  
    return 3619;
  else if( aa == "IRA")  
    return 3620;
  else if( aa == "IRR")  
    return 3621;
  else if( aa == "IRN")  
    return 3622;
  else if( aa == "IRD")  
    return 3623;
  else if( aa == "IRC")  
    return 3624;
  else if( aa == "IRQ")  
    return 3625;
  else if( aa == "IRE")  
    return 3626;
  else if( aa == "IRG")  
    return 3627;
  else if( aa == "IRH")  
    return 3628;
  else if( aa == "IRI")  
    return 3629;
  else if( aa == "IRL")  
    return 3630;
  else if( aa == "IRK")  
    return 3631;
  else if( aa == "IRM")  
    return 3632;
  else if( aa == "IRF")  
    return 3633;
  else if( aa == "IRP")  
    return 3634;
  else if( aa == "IRS")  
    return 3635;
  else if( aa == "IRT")  
    return 3636;
  else if( aa == "IRW")  
    return 3637;
  else if( aa == "IRY")  
    return 3638;
  else if( aa == "IRV")  
    return 3639;
  else if( aa == "INA")  
    return 3640;
  else if( aa == "INR")  
    return 3641;
  else if( aa == "INN")  
    return 3642;
  else if( aa == "IND")  
    return 3643;
  else if( aa == "INC")  
    return 3644;
  else if( aa == "INQ")  
    return 3645;
  else if( aa == "INE")  
    return 3646;
  else if( aa == "ING")  
    return 3647;
  else if( aa == "INH")  
    return 3648;
  else if( aa == "INI")  
    return 3649;
  else if( aa == "INL")  
    return 3650;
  else if( aa == "INK")  
    return 3651;
  else if( aa == "INM")  
    return 3652;
  else if( aa == "INF")  
    return 3653;
  else if( aa == "INP")  
    return 3654;
  else if( aa == "INS")  
    return 3655;
  else if( aa == "INT")  
    return 3656;
  else if( aa == "INW")  
    return 3657;
  else if( aa == "INY")  
    return 3658;
  else if( aa == "INV")  
    return 3659;
  else if( aa == "IDA")  
    return 3660;
  else if( aa == "IDR")  
    return 3661;
  else if( aa == "IDN")  
    return 3662;
  else if( aa == "IDD")  
    return 3663;
  else if( aa == "IDC")  
    return 3664;
  else if( aa == "IDQ")  
    return 3665;
  else if( aa == "IDE")  
    return 3666;
  else if( aa == "IDG")  
    return 3667;
  else if( aa == "IDH")  
    return 3668;
  else if( aa == "IDI")  
    return 3669;
  else if( aa == "IDL")  
    return 3670;
  else if( aa == "IDK")  
    return 3671;
  else if( aa == "IDM")  
    return 3672;
  else if( aa == "IDF")  
    return 3673;
  else if( aa == "IDP")  
    return 3674;
  else if( aa == "IDS")  
    return 3675;
  else if( aa == "IDT")  
    return 3676;
  else if( aa == "IDW")  
    return 3677;
  else if( aa == "IDY")  
    return 3678;
  else if( aa == "IDV")  
    return 3679;
  else if( aa == "ICA")  
    return 3680;
  else if( aa == "ICR")  
    return 3681;
  else if( aa == "ICN")  
    return 3682;
  else if( aa == "ICD")  
    return 3683;
  else if( aa == "ICC")  
    return 3684;
  else if( aa == "ICQ")  
    return 3685;
  else if( aa == "ICE")  
    return 3686;
  else if( aa == "ICG")  
    return 3687;
  else if( aa == "ICH")  
    return 3688;
  else if( aa == "ICI")  
    return 3689;
  else if( aa == "ICL")  
    return 3690;
  else if( aa == "ICK")  
    return 3691;
  else if( aa == "ICM")  
    return 3692;
  else if( aa == "ICF")  
    return 3693;
  else if( aa == "ICP")  
    return 3694;
  else if( aa == "ICS")  
    return 3695;
  else if( aa == "ICT")  
    return 3696;
  else if( aa == "ICW")  
    return 3697;
  else if( aa == "ICY")  
    return 3698;
  else if( aa == "ICV")  
    return 3699;
  else if( aa == "IQA")  
    return 3700;
  else if( aa == "IQR")  
    return 3701;
  else if( aa == "IQN")  
    return 3702;
  else if( aa == "IQD")  
    return 3703;
  else if( aa == "IQC")  
    return 3704;
  else if( aa == "IQQ")  
    return 3705;
  else if( aa == "IQE")  
    return 3706;
  else if( aa == "IQG")  
    return 3707;
  else if( aa == "IQH")  
    return 3708;
  else if( aa == "IQI")  
    return 3709;
  else if( aa == "IQL")  
    return 3710;
  else if( aa == "IQK")  
    return 3711;
  else if( aa == "IQM")  
    return 3712;
  else if( aa == "IQF")  
    return 3713;
  else if( aa == "IQP")  
    return 3714;
  else if( aa == "IQS")  
    return 3715;
  else if( aa == "IQT")  
    return 3716;
  else if( aa == "IQW")  
    return 3717;
  else if( aa == "IQY")  
    return 3718;
  else if( aa == "IQV")  
    return 3719;
  else if( aa == "IEA")  
    return 3720;
  else if( aa == "IER")  
    return 3721;
  else if( aa == "IEN")  
    return 3722;
  else if( aa == "IED")  
    return 3723;
  else if( aa == "IEC")  
    return 3724;
  else if( aa == "IEQ")  
    return 3725;
  else if( aa == "IEE")  
    return 3726;
  else if( aa == "IEG")  
    return 3727;
  else if( aa == "IEH")  
    return 3728;
  else if( aa == "IEI")  
    return 3729;
  else if( aa == "IEL")  
    return 3730;
  else if( aa == "IEK")  
    return 3731;
  else if( aa == "IEM")  
    return 3732;
  else if( aa == "IEF")  
    return 3733;
  else if( aa == "IEP")  
    return 3734;
  else if( aa == "IES")  
    return 3735;
  else if( aa == "IET")  
    return 3736;
  else if( aa == "IEW")  
    return 3737;
  else if( aa == "IEY")  
    return 3738;
  else if( aa == "IEV")  
    return 3739;
  else if( aa == "IGA")  
    return 3740;
  else if( aa == "IGR")  
    return 3741;
  else if( aa == "IGN")  
    return 3742;
  else if( aa == "IGD")  
    return 3743;
  else if( aa == "IGC")  
    return 3744;
  else if( aa == "IGQ")  
    return 3745;
  else if( aa == "IGE")  
    return 3746;
  else if( aa == "IGG")  
    return 3747;
  else if( aa == "IGH")  
    return 3748;
  else if( aa == "IGI")  
    return 3749;
  else if( aa == "IGL")  
    return 3750;
  else if( aa == "IGK")  
    return 3751;
  else if( aa == "IGM")  
    return 3752;
  else if( aa == "IGF")  
    return 3753;
  else if( aa == "IGP")  
    return 3754;
  else if( aa == "IGS")  
    return 3755;
  else if( aa == "IGT")  
    return 3756;
  else if( aa == "IGW")  
    return 3757;
  else if( aa == "IGY")  
    return 3758;
  else if( aa == "IGV")  
    return 3759;
  else if( aa == "IHA")  
    return 3760;
  else if( aa == "IHR")  
    return 3761;
  else if( aa == "IHN")  
    return 3762;
  else if( aa == "IHD")  
    return 3763;
  else if( aa == "IHC")  
    return 3764;
  else if( aa == "IHQ")  
    return 3765;
  else if( aa == "IHE")  
    return 3766;
  else if( aa == "IHG")  
    return 3767;
  else if( aa == "IHH")  
    return 3768;
  else if( aa == "IHI")  
    return 3769;
  else if( aa == "IHL")  
    return 3770;
  else if( aa == "IHK")  
    return 3771;
  else if( aa == "IHM")  
    return 3772;
  else if( aa == "IHF")  
    return 3773;
  else if( aa == "IHP")  
    return 3774;
  else if( aa == "IHS")  
    return 3775;
  else if( aa == "IHT")  
    return 3776;
  else if( aa == "IHW")  
    return 3777;
  else if( aa == "IHY")  
    return 3778;
  else if( aa == "IHV")  
    return 3779;
  else if( aa == "IIA")  
    return 3780;
  else if( aa == "IIR")  
    return 3781;
  else if( aa == "IIN")  
    return 3782;
  else if( aa == "IID")  
    return 3783;
  else if( aa == "IIC")  
    return 3784;
  else if( aa == "IIQ")  
    return 3785;
  else if( aa == "IIE")  
    return 3786;
  else if( aa == "IIG")  
    return 3787;
  else if( aa == "IIH")  
    return 3788;
  else if( aa == "III")  
    return 3789;
  else if( aa == "IIL")  
    return 3790;
  else if( aa == "IIK")  
    return 3791;
  else if( aa == "IIM")  
    return 3792;
  else if( aa == "IIF")  
    return 3793;
  else if( aa == "IIP")  
    return 3794;
  else if( aa == "IIS")  
    return 3795;
  else if( aa == "IIT")  
    return 3796;
  else if( aa == "IIW")  
    return 3797;
  else if( aa == "IIY")  
    return 3798;
  else if( aa == "IIV")  
    return 3799;
  else if( aa == "ILA")  
    return 3800;
  else if( aa == "ILR")  
    return 3801;
  else if( aa == "ILN")  
    return 3802;
  else if( aa == "ILD")  
    return 3803;
  else if( aa == "ILC")  
    return 3804;
  else if( aa == "ILQ")  
    return 3805;
  else if( aa == "ILE")  
    return 3806;
  else if( aa == "ILG")  
    return 3807;
  else if( aa == "ILH")  
    return 3808;
  else if( aa == "ILI")  
    return 3809;
  else if( aa == "ILL")  
    return 3810;
  else if( aa == "ILK")  
    return 3811;
  else if( aa == "ILM")  
    return 3812;
  else if( aa == "ILF")  
    return 3813;
  else if( aa == "ILP")  
    return 3814;
  else if( aa == "ILS")  
    return 3815;
  else if( aa == "ILT")  
    return 3816;
  else if( aa == "ILW")  
    return 3817;
  else if( aa == "ILY")  
    return 3818;
  else if( aa == "ILV")  
    return 3819;
  else if( aa == "IKA")  
    return 3820;
  else if( aa == "IKR")  
    return 3821;
  else if( aa == "IKN")  
    return 3822;
  else if( aa == "IKD")  
    return 3823;
  else if( aa == "IKC")  
    return 3824;
  else if( aa == "IKQ")  
    return 3825;
  else if( aa == "IKE")  
    return 3826;
  else if( aa == "IKG")  
    return 3827;
  else if( aa == "IKH")  
    return 3828;
  else if( aa == "IKI")  
    return 3829;
  else if( aa == "IKL")  
    return 3830;
  else if( aa == "IKK")  
    return 3831;
  else if( aa == "IKM")  
    return 3832;
  else if( aa == "IKF")  
    return 3833;
  else if( aa == "IKP")  
    return 3834;
  else if( aa == "IKS")  
    return 3835;
  else if( aa == "IKT")  
    return 3836;
  else if( aa == "IKW")  
    return 3837;
  else if( aa == "IKY")  
    return 3838;
  else if( aa == "IKV")  
    return 3839;
  else if( aa == "IMA")  
    return 3840;
  else if( aa == "IMR")  
    return 3841;
  else if( aa == "IMN")  
    return 3842;
  else if( aa == "IMD")  
    return 3843;
  else if( aa == "IMC")  
    return 3844;
  else if( aa == "IMQ")  
    return 3845;
  else if( aa == "IME")  
    return 3846;
  else if( aa == "IMG")  
    return 3847;
  else if( aa == "IMH")  
    return 3848;
  else if( aa == "IMI")  
    return 3849;
  else if( aa == "IML")  
    return 3850;
  else if( aa == "IMK")  
    return 3851;
  else if( aa == "IMM")  
    return 3852;
  else if( aa == "IMF")  
    return 3853;
  else if( aa == "IMP")  
    return 3854;
  else if( aa == "IMS")  
    return 3855;
  else if( aa == "IMT")  
    return 3856;
  else if( aa == "IMW")  
    return 3857;
  else if( aa == "IMY")  
    return 3858;
  else if( aa == "IMV")  
    return 3859;
  else if( aa == "IFA")  
    return 3860;
  else if( aa == "IFR")  
    return 3861;
  else if( aa == "IFN")  
    return 3862;
  else if( aa == "IFD")  
    return 3863;
  else if( aa == "IFC")  
    return 3864;
  else if( aa == "IFQ")  
    return 3865;
  else if( aa == "IFE")  
    return 3866;
  else if( aa == "IFG")  
    return 3867;
  else if( aa == "IFH")  
    return 3868;
  else if( aa == "IFI")  
    return 3869;
  else if( aa == "IFL")  
    return 3870;
  else if( aa == "IFK")  
    return 3871;
  else if( aa == "IFM")  
    return 3872;
  else if( aa == "IFF")  
    return 3873;
  else if( aa == "IFP")  
    return 3874;
  else if( aa == "IFS")  
    return 3875;
  else if( aa == "IFT")  
    return 3876;
  else if( aa == "IFW")  
    return 3877;
  else if( aa == "IFY")  
    return 3878;
  else if( aa == "IFV")  
    return 3879;
  else if( aa == "IPA")  
    return 3880;
  else if( aa == "IPR")  
    return 3881;
  else if( aa == "IPN")  
    return 3882;
  else if( aa == "IPD")  
    return 3883;
  else if( aa == "IPC")  
    return 3884;
  else if( aa == "IPQ")  
    return 3885;
  else if( aa == "IPE")  
    return 3886;
  else if( aa == "IPG")  
    return 3887;
  else if( aa == "IPH")  
    return 3888;
  else if( aa == "IPI")  
    return 3889;
  else if( aa == "IPL")  
    return 3890;
  else if( aa == "IPK")  
    return 3891;
  else if( aa == "IPM")  
    return 3892;
  else if( aa == "IPF")  
    return 3893;
  else if( aa == "IPP")  
    return 3894;
  else if( aa == "IPS")  
    return 3895;
  else if( aa == "IPT")  
    return 3896;
  else if( aa == "IPW")  
    return 3897;
  else if( aa == "IPY")  
    return 3898;
  else if( aa == "IPV")  
    return 3899;
  else if( aa == "ISA")  
    return 3900;
  else if( aa == "ISR")  
    return 3901;
  else if( aa == "ISN")  
    return 3902;
  else if( aa == "ISD")  
    return 3903;
  else if( aa == "ISC")  
    return 3904;
  else if( aa == "ISQ")  
    return 3905;
  else if( aa == "ISE")  
    return 3906;
  else if( aa == "ISG")  
    return 3907;
  else if( aa == "ISH")  
    return 3908;
  else if( aa == "ISI")  
    return 3909;
  else if( aa == "ISL")  
    return 3910;
  else if( aa == "ISK")  
    return 3911;
  else if( aa == "ISM")  
    return 3912;
  else if( aa == "ISF")  
    return 3913;
  else if( aa == "ISP")  
    return 3914;
  else if( aa == "ISS")  
    return 3915;
  else if( aa == "IST")  
    return 3916;
  else if( aa == "ISW")  
    return 3917;
  else if( aa == "ISY")  
    return 3918;
  else if( aa == "ISV")  
    return 3919;
  else if( aa == "ITA")  
    return 3920;
  else if( aa == "ITR")  
    return 3921;
  else if( aa == "ITN")  
    return 3922;
  else if( aa == "ITD")  
    return 3923;
  else if( aa == "ITC")  
    return 3924;
  else if( aa == "ITQ")  
    return 3925;
  else if( aa == "ITE")  
    return 3926;
  else if( aa == "ITG")  
    return 3927;
  else if( aa == "ITH")  
    return 3928;
  else if( aa == "ITI")  
    return 3929;
  else if( aa == "ITL")  
    return 3930;
  else if( aa == "ITK")  
    return 3931;
  else if( aa == "ITM")  
    return 3932;
  else if( aa == "ITF")  
    return 3933;
  else if( aa == "ITP")  
    return 3934;
  else if( aa == "ITS")  
    return 3935;
  else if( aa == "ITT")  
    return 3936;
  else if( aa == "ITW")  
    return 3937;
  else if( aa == "ITY")  
    return 3938;
  else if( aa == "ITV")  
    return 3939;
  else if( aa == "IWA")  
    return 3940;
  else if( aa == "IWR")  
    return 3941;
  else if( aa == "IWN")  
    return 3942;
  else if( aa == "IWD")  
    return 3943;
  else if( aa == "IWC")  
    return 3944;
  else if( aa == "IWQ")  
    return 3945;
  else if( aa == "IWE")  
    return 3946;
  else if( aa == "IWG")  
    return 3947;
  else if( aa == "IWH")  
    return 3948;
  else if( aa == "IWI")  
    return 3949;
  else if( aa == "IWL")  
    return 3950;
  else if( aa == "IWK")  
    return 3951;
  else if( aa == "IWM")  
    return 3952;
  else if( aa == "IWF")  
    return 3953;
  else if( aa == "IWP")  
    return 3954;
  else if( aa == "IWS")  
    return 3955;
  else if( aa == "IWT")  
    return 3956;
  else if( aa == "IWW")  
    return 3957;
  else if( aa == "IWY")  
    return 3958;
  else if( aa == "IWV")  
    return 3959;
  else if( aa == "IYA")  
    return 3960;
  else if( aa == "IYR")  
    return 3961;
  else if( aa == "IYN")  
    return 3962;
  else if( aa == "IYD")  
    return 3963;
  else if( aa == "IYC")  
    return 3964;
  else if( aa == "IYQ")  
    return 3965;
  else if( aa == "IYE")  
    return 3966;
  else if( aa == "IYG")  
    return 3967;
  else if( aa == "IYH")  
    return 3968;
  else if( aa == "IYI")  
    return 3969;
  else if( aa == "IYL")  
    return 3970;
  else if( aa == "IYK")  
    return 3971;
  else if( aa == "IYM")  
    return 3972;
  else if( aa == "IYF")  
    return 3973;
  else if( aa == "IYP")  
    return 3974;
  else if( aa == "IYS")  
    return 3975;
  else if( aa == "IYT")  
    return 3976;
  else if( aa == "IYW")  
    return 3977;
  else if( aa == "IYY")  
    return 3978;
  else if( aa == "IYV")  
    return 3979;
  else if( aa == "IVA")  
    return 3980;
  else if( aa == "IVR")  
    return 3981;
  else if( aa == "IVN")  
    return 3982;
  else if( aa == "IVD")  
    return 3983;
  else if( aa == "IVC")  
    return 3984;
  else if( aa == "IVQ")  
    return 3985;
  else if( aa == "IVE")  
    return 3986;
  else if( aa == "IVG")  
    return 3987;
  else if( aa == "IVH")  
    return 3988;
  else if( aa == "IVI")  
    return 3989;
  else if( aa == "IVL")  
    return 3990;
  else if( aa == "IVK")  
    return 3991;
  else if( aa == "IVM")  
    return 3992;
  else if( aa == "IVF")  
    return 3993;
  else if( aa == "IVP")  
    return 3994;
  else if( aa == "IVS")  
    return 3995;
  else if( aa == "IVT")  
    return 3996;
  else if( aa == "IVW")  
    return 3997;
  else if( aa == "IVY")  
    return 3998;
  else if( aa == "IVV")  
    return 3999;
  else if( aa == "LAA")  
    return 4000;
  else if( aa == "LAR")  
    return 4001;
  else if( aa == "LAN")  
    return 4002;
  else if( aa == "LAD")  
    return 4003;
  else if( aa == "LAC")  
    return 4004;
  else if( aa == "LAQ")  
    return 4005;
  else if( aa == "LAE")  
    return 4006;
  else if( aa == "LAG")  
    return 4007;
  else if( aa == "LAH")  
    return 4008;
  else if( aa == "LAI")  
    return 4009;
  else if( aa == "LAL")  
    return 4010;
  else if( aa == "LAK")  
    return 4011;
  else if( aa == "LAM")  
    return 4012;
  else if( aa == "LAF")  
    return 4013;
  else if( aa == "LAP")  
    return 4014;
  else if( aa == "LAS")  
    return 4015;
  else if( aa == "LAT")  
    return 4016;
  else if( aa == "LAW")  
    return 4017;
  else if( aa == "LAY")  
    return 4018;
  else if( aa == "LAV")  
    return 4019;
  else if( aa == "LRA")  
    return 4020;
  else if( aa == "LRR")  
    return 4021;
  else if( aa == "LRN")  
    return 4022;
  else if( aa == "LRD")  
    return 4023;
  else if( aa == "LRC")  
    return 4024;
  else if( aa == "LRQ")  
    return 4025;
  else if( aa == "LRE")  
    return 4026;
  else if( aa == "LRG")  
    return 4027;
  else if( aa == "LRH")  
    return 4028;
  else if( aa == "LRI")  
    return 4029;
  else if( aa == "LRL")  
    return 4030;
  else if( aa == "LRK")  
    return 4031;
  else if( aa == "LRM")  
    return 4032;
  else if( aa == "LRF")  
    return 4033;
  else if( aa == "LRP")  
    return 4034;
  else if( aa == "LRS")  
    return 4035;
  else if( aa == "LRT")  
    return 4036;
  else if( aa == "LRW")  
    return 4037;
  else if( aa == "LRY")  
    return 4038;
  else if( aa == "LRV")  
    return 4039;
  else if( aa == "LNA")  
    return 4040;
  else if( aa == "LNR")  
    return 4041;
  else if( aa == "LNN")  
    return 4042;
  else if( aa == "LND")  
    return 4043;
  else if( aa == "LNC")  
    return 4044;
  else if( aa == "LNQ")  
    return 4045;
  else if( aa == "LNE")  
    return 4046;
  else if( aa == "LNG")  
    return 4047;
  else if( aa == "LNH")  
    return 4048;
  else if( aa == "LNI")  
    return 4049;
  else if( aa == "LNL")  
    return 4050;
  else if( aa == "LNK")  
    return 4051;
  else if( aa == "LNM")  
    return 4052;
  else if( aa == "LNF")  
    return 4053;
  else if( aa == "LNP")  
    return 4054;
  else if( aa == "LNS")  
    return 4055;
  else if( aa == "LNT")  
    return 4056;
  else if( aa == "LNW")  
    return 4057;
  else if( aa == "LNY")  
    return 4058;
  else if( aa == "LNV")  
    return 4059;
  else if( aa == "LDA")  
    return 4060;
  else if( aa == "LDR")  
    return 4061;
  else if( aa == "LDN")  
    return 4062;
  else if( aa == "LDD")  
    return 4063;
  else if( aa == "LDC")  
    return 4064;
  else if( aa == "LDQ")  
    return 4065;
  else if( aa == "LDE")  
    return 4066;
  else if( aa == "LDG")  
    return 4067;
  else if( aa == "LDH")  
    return 4068;
  else if( aa == "LDI")  
    return 4069;
  else if( aa == "LDL")  
    return 4070;
  else if( aa == "LDK")  
    return 4071;
  else if( aa == "LDM")  
    return 4072;
  else if( aa == "LDF")  
    return 4073;
  else if( aa == "LDP")  
    return 4074;
  else if( aa == "LDS")  
    return 4075;
  else if( aa == "LDT")  
    return 4076;
  else if( aa == "LDW")  
    return 4077;
  else if( aa == "LDY")  
    return 4078;
  else if( aa == "LDV")  
    return 4079;
  else if( aa == "LCA")  
    return 4080;
  else if( aa == "LCR")  
    return 4081;
  else if( aa == "LCN")  
    return 4082;
  else if( aa == "LCD")  
    return 4083;
  else if( aa == "LCC")  
    return 4084;
  else if( aa == "LCQ")  
    return 4085;
  else if( aa == "LCE")  
    return 4086;
  else if( aa == "LCG")  
    return 4087;
  else if( aa == "LCH")  
    return 4088;
  else if( aa == "LCI")  
    return 4089;
  else if( aa == "LCL")  
    return 4090;
  else if( aa == "LCK")  
    return 4091;
  else if( aa == "LCM")  
    return 4092;
  else if( aa == "LCF")  
    return 4093;
  else if( aa == "LCP")  
    return 4094;
  else if( aa == "LCS")  
    return 4095;
  else if( aa == "LCT")  
    return 4096;
  else if( aa == "LCW")  
    return 4097;
  else if( aa == "LCY")  
    return 4098;
  else if( aa == "LCV")  
    return 4099;
  else if( aa == "LQA")  
    return 4100;
  else if( aa == "LQR")  
    return 4101;
  else if( aa == "LQN")  
    return 4102;
  else if( aa == "LQD")  
    return 4103;
  else if( aa == "LQC")  
    return 4104;
  else if( aa == "LQQ")  
    return 4105;
  else if( aa == "LQE")  
    return 4106;
  else if( aa == "LQG")  
    return 4107;
  else if( aa == "LQH")  
    return 4108;
  else if( aa == "LQI")  
    return 4109;
  else if( aa == "LQL")  
    return 4110;
  else if( aa == "LQK")  
    return 4111;
  else if( aa == "LQM")  
    return 4112;
  else if( aa == "LQF")  
    return 4113;
  else if( aa == "LQP")  
    return 4114;
  else if( aa == "LQS")  
    return 4115;
  else if( aa == "LQT")  
    return 4116;
  else if( aa == "LQW")  
    return 4117;
  else if( aa == "LQY")  
    return 4118;
  else if( aa == "LQV")  
    return 4119;
  else if( aa == "LEA")  
    return 4120;
  else if( aa == "LER")  
    return 4121;
  else if( aa == "LEN")  
    return 4122;
  else if( aa == "LED")  
    return 4123;
  else if( aa == "LEC")  
    return 4124;
  else if( aa == "LEQ")  
    return 4125;
  else if( aa == "LEE")  
    return 4126;
  else if( aa == "LEG")  
    return 4127;
  else if( aa == "LEH")  
    return 4128;
  else if( aa == "LEI")  
    return 4129;
  else if( aa == "LEL")  
    return 4130;
  else if( aa == "LEK")  
    return 4131;
  else if( aa == "LEM")  
    return 4132;
  else if( aa == "LEF")  
    return 4133;
  else if( aa == "LEP")  
    return 4134;
  else if( aa == "LES")  
    return 4135;
  else if( aa == "LET")  
    return 4136;
  else if( aa == "LEW")  
    return 4137;
  else if( aa == "LEY")  
    return 4138;
  else if( aa == "LEV")  
    return 4139;
  else if( aa == "LGA")  
    return 4140;
  else if( aa == "LGR")  
    return 4141;
  else if( aa == "LGN")  
    return 4142;
  else if( aa == "LGD")  
    return 4143;
  else if( aa == "LGC")  
    return 4144;
  else if( aa == "LGQ")  
    return 4145;
  else if( aa == "LGE")  
    return 4146;
  else if( aa == "LGG")  
    return 4147;
  else if( aa == "LGH")  
    return 4148;
  else if( aa == "LGI")  
    return 4149;
  else if( aa == "LGL")  
    return 4150;
  else if( aa == "LGK")  
    return 4151;
  else if( aa == "LGM")  
    return 4152;
  else if( aa == "LGF")  
    return 4153;
  else if( aa == "LGP")  
    return 4154;
  else if( aa == "LGS")  
    return 4155;
  else if( aa == "LGT")  
    return 4156;
  else if( aa == "LGW")  
    return 4157;
  else if( aa == "LGY")  
    return 4158;
  else if( aa == "LGV")  
    return 4159;
  else if( aa == "LHA")  
    return 4160;
  else if( aa == "LHR")  
    return 4161;
  else if( aa == "LHN")  
    return 4162;
  else if( aa == "LHD")  
    return 4163;
  else if( aa == "LHC")  
    return 4164;
  else if( aa == "LHQ")  
    return 4165;
  else if( aa == "LHE")  
    return 4166;
  else if( aa == "LHG")  
    return 4167;
  else if( aa == "LHH")  
    return 4168;
  else if( aa == "LHI")  
    return 4169;
  else if( aa == "LHL")  
    return 4170;
  else if( aa == "LHK")  
    return 4171;
  else if( aa == "LHM")  
    return 4172;
  else if( aa == "LHF")  
    return 4173;
  else if( aa == "LHP")  
    return 4174;
  else if( aa == "LHS")  
    return 4175;
  else if( aa == "LHT")  
    return 4176;
  else if( aa == "LHW")  
    return 4177;
  else if( aa == "LHY")  
    return 4178;
  else if( aa == "LHV")  
    return 4179;
  else if( aa == "LIA")  
    return 4180;
  else if( aa == "LIR")  
    return 4181;
  else if( aa == "LIN")  
    return 4182;
  else if( aa == "LID")  
    return 4183;
  else if( aa == "LIC")  
    return 4184;
  else if( aa == "LIQ")  
    return 4185;
  else if( aa == "LIE")  
    return 4186;
  else if( aa == "LIG")  
    return 4187;
  else if( aa == "LIH")  
    return 4188;
  else if( aa == "LII")  
    return 4189;
  else if( aa == "LIL")  
    return 4190;
  else if( aa == "LIK")  
    return 4191;
  else if( aa == "LIM")  
    return 4192;
  else if( aa == "LIF")  
    return 4193;
  else if( aa == "LIP")  
    return 4194;
  else if( aa == "LIS")  
    return 4195;
  else if( aa == "LIT")  
    return 4196;
  else if( aa == "LIW")  
    return 4197;
  else if( aa == "LIY")  
    return 4198;
  else if( aa == "LIV")  
    return 4199;
  else if( aa == "LLA")  
    return 4200;
  else if( aa == "LLR")  
    return 4201;
  else if( aa == "LLN")  
    return 4202;
  else if( aa == "LLD")  
    return 4203;
  else if( aa == "LLC")  
    return 4204;
  else if( aa == "LLQ")  
    return 4205;
  else if( aa == "LLE")  
    return 4206;
  else if( aa == "LLG")  
    return 4207;
  else if( aa == "LLH")  
    return 4208;
  else if( aa == "LLI")  
    return 4209;
  else if( aa == "LLL")  
    return 4210;
  else if( aa == "LLK")  
    return 4211;
  else if( aa == "LLM")  
    return 4212;
  else if( aa == "LLF")  
    return 4213;
  else if( aa == "LLP")  
    return 4214;
  else if( aa == "LLS")  
    return 4215;
  else if( aa == "LLT")  
    return 4216;
  else if( aa == "LLW")  
    return 4217;
  else if( aa == "LLY")  
    return 4218;
  else if( aa == "LLV")  
    return 4219;
  else if( aa == "LKA")  
    return 4220;
  else if( aa == "LKR")  
    return 4221;
  else if( aa == "LKN")  
    return 4222;
  else if( aa == "LKD")  
    return 4223;
  else if( aa == "LKC")  
    return 4224;
  else if( aa == "LKQ")  
    return 4225;
  else if( aa == "LKE")  
    return 4226;
  else if( aa == "LKG")  
    return 4227;
  else if( aa == "LKH")  
    return 4228;
  else if( aa == "LKI")  
    return 4229;
  else if( aa == "LKL")  
    return 4230;
  else if( aa == "LKK")  
    return 4231;
  else if( aa == "LKM")  
    return 4232;
  else if( aa == "LKF")  
    return 4233;
  else if( aa == "LKP")  
    return 4234;
  else if( aa == "LKS")  
    return 4235;
  else if( aa == "LKT")  
    return 4236;
  else if( aa == "LKW")  
    return 4237;
  else if( aa == "LKY")  
    return 4238;
  else if( aa == "LKV")  
    return 4239;
  else if( aa == "LMA")  
    return 4240;
  else if( aa == "LMR")  
    return 4241;
  else if( aa == "LMN")  
    return 4242;
  else if( aa == "LMD")  
    return 4243;
  else if( aa == "LMC")  
    return 4244;
  else if( aa == "LMQ")  
    return 4245;
  else if( aa == "LME")  
    return 4246;
  else if( aa == "LMG")  
    return 4247;
  else if( aa == "LMH")  
    return 4248;
  else if( aa == "LMI")  
    return 4249;
  else if( aa == "LML")  
    return 4250;
  else if( aa == "LMK")  
    return 4251;
  else if( aa == "LMM")  
    return 4252;
  else if( aa == "LMF")  
    return 4253;
  else if( aa == "LMP")  
    return 4254;
  else if( aa == "LMS")  
    return 4255;
  else if( aa == "LMT")  
    return 4256;
  else if( aa == "LMW")  
    return 4257;
  else if( aa == "LMY")  
    return 4258;
  else if( aa == "LMV")  
    return 4259;
  else if( aa == "LFA")  
    return 4260;
  else if( aa == "LFR")  
    return 4261;
  else if( aa == "LFN")  
    return 4262;
  else if( aa == "LFD")  
    return 4263;
  else if( aa == "LFC")  
    return 4264;
  else if( aa == "LFQ")  
    return 4265;
  else if( aa == "LFE")  
    return 4266;
  else if( aa == "LFG")  
    return 4267;
  else if( aa == "LFH")  
    return 4268;
  else if( aa == "LFI")  
    return 4269;
  else if( aa == "LFL")  
    return 4270;
  else if( aa == "LFK")  
    return 4271;
  else if( aa == "LFM")  
    return 4272;
  else if( aa == "LFF")  
    return 4273;
  else if( aa == "LFP")  
    return 4274;
  else if( aa == "LFS")  
    return 4275;
  else if( aa == "LFT")  
    return 4276;
  else if( aa == "LFW")  
    return 4277;
  else if( aa == "LFY")  
    return 4278;
  else if( aa == "LFV")  
    return 4279;
  else if( aa == "LPA")  
    return 4280;
  else if( aa == "LPR")  
    return 4281;
  else if( aa == "LPN")  
    return 4282;
  else if( aa == "LPD")  
    return 4283;
  else if( aa == "LPC")  
    return 4284;
  else if( aa == "LPQ")  
    return 4285;
  else if( aa == "LPE")  
    return 4286;
  else if( aa == "LPG")  
    return 4287;
  else if( aa == "LPH")  
    return 4288;
  else if( aa == "LPI")  
    return 4289;
  else if( aa == "LPL")  
    return 4290;
  else if( aa == "LPK")  
    return 4291;
  else if( aa == "LPM")  
    return 4292;
  else if( aa == "LPF")  
    return 4293;
  else if( aa == "LPP")  
    return 4294;
  else if( aa == "LPS")  
    return 4295;
  else if( aa == "LPT")  
    return 4296;
  else if( aa == "LPW")  
    return 4297;
  else if( aa == "LPY")  
    return 4298;
  else if( aa == "LPV")  
    return 4299;
  else if( aa == "LSA")  
    return 4300;
  else if( aa == "LSR")  
    return 4301;
  else if( aa == "LSN")  
    return 4302;
  else if( aa == "LSD")  
    return 4303;
  else if( aa == "LSC")  
    return 4304;
  else if( aa == "LSQ")  
    return 4305;
  else if( aa == "LSE")  
    return 4306;
  else if( aa == "LSG")  
    return 4307;
  else if( aa == "LSH")  
    return 4308;
  else if( aa == "LSI")  
    return 4309;
  else if( aa == "LSL")  
    return 4310;
  else if( aa == "LSK")  
    return 4311;
  else if( aa == "LSM")  
    return 4312;
  else if( aa == "LSF")  
    return 4313;
  else if( aa == "LSP")  
    return 4314;
  else if( aa == "LSS")  
    return 4315;
  else if( aa == "LST")  
    return 4316;
  else if( aa == "LSW")  
    return 4317;
  else if( aa == "LSY")  
    return 4318;
  else if( aa == "LSV")  
    return 4319;
  else if( aa == "LTA")  
    return 4320;
  else if( aa == "LTR")  
    return 4321;
  else if( aa == "LTN")  
    return 4322;
  else if( aa == "LTD")  
    return 4323;
  else if( aa == "LTC")  
    return 4324;
  else if( aa == "LTQ")  
    return 4325;
  else if( aa == "LTE")  
    return 4326;
  else if( aa == "LTG")  
    return 4327;
  else if( aa == "LTH")  
    return 4328;
  else if( aa == "LTI")  
    return 4329;
  else if( aa == "LTL")  
    return 4330;
  else if( aa == "LTK")  
    return 4331;
  else if( aa == "LTM")  
    return 4332;
  else if( aa == "LTF")  
    return 4333;
  else if( aa == "LTP")  
    return 4334;
  else if( aa == "LTS")  
    return 4335;
  else if( aa == "LTT")  
    return 4336;
  else if( aa == "LTW")  
    return 4337;
  else if( aa == "LTY")  
    return 4338;
  else if( aa == "LTV")  
    return 4339;
  else if( aa == "LWA")  
    return 4340;
  else if( aa == "LWR")  
    return 4341;
  else if( aa == "LWN")  
    return 4342;
  else if( aa == "LWD")  
    return 4343;
  else if( aa == "LWC")  
    return 4344;
  else if( aa == "LWQ")  
    return 4345;
  else if( aa == "LWE")  
    return 4346;
  else if( aa == "LWG")  
    return 4347;
  else if( aa == "LWH")  
    return 4348;
  else if( aa == "LWI")  
    return 4349;
  else if( aa == "LWL")  
    return 4350;
  else if( aa == "LWK")  
    return 4351;
  else if( aa == "LWM")  
    return 4352;
  else if( aa == "LWF")  
    return 4353;
  else if( aa == "LWP")  
    return 4354;
  else if( aa == "LWS")  
    return 4355;
  else if( aa == "LWT")  
    return 4356;
  else if( aa == "LWW")  
    return 4357;
  else if( aa == "LWY")  
    return 4358;
  else if( aa == "LWV")  
    return 4359;
  else if( aa == "LYA")  
    return 4360;
  else if( aa == "LYR")  
    return 4361;
  else if( aa == "LYN")  
    return 4362;
  else if( aa == "LYD")  
    return 4363;
  else if( aa == "LYC")  
    return 4364;
  else if( aa == "LYQ")  
    return 4365;
  else if( aa == "LYE")  
    return 4366;
  else if( aa == "LYG")  
    return 4367;
  else if( aa == "LYH")  
    return 4368;
  else if( aa == "LYI")  
    return 4369;
  else if( aa == "LYL")  
    return 4370;
  else if( aa == "LYK")  
    return 4371;
  else if( aa == "LYM")  
    return 4372;
  else if( aa == "LYF")  
    return 4373;
  else if( aa == "LYP")  
    return 4374;
  else if( aa == "LYS")  
    return 4375;
  else if( aa == "LYT")  
    return 4376;
  else if( aa == "LYW")  
    return 4377;
  else if( aa == "LYY")  
    return 4378;
  else if( aa == "LYV")  
    return 4379;
  else if( aa == "LVA")  
    return 4380;
  else if( aa == "LVR")  
    return 4381;
  else if( aa == "LVN")  
    return 4382;
  else if( aa == "LVD")  
    return 4383;
  else if( aa == "LVC")  
    return 4384;
  else if( aa == "LVQ")  
    return 4385;
  else if( aa == "LVE")  
    return 4386;
  else if( aa == "LVG")  
    return 4387;
  else if( aa == "LVH")  
    return 4388;
  else if( aa == "LVI")  
    return 4389;
  else if( aa == "LVL")  
    return 4390;
  else if( aa == "LVK")  
    return 4391;
  else if( aa == "LVM")  
    return 4392;
  else if( aa == "LVF")  
    return 4393;
  else if( aa == "LVP")  
    return 4394;
  else if( aa == "LVS")  
    return 4395;
  else if( aa == "LVT")  
    return 4396;
  else if( aa == "LVW")  
    return 4397;
  else if( aa == "LVY")  
    return 4398;
  else if( aa == "LVV")  
    return 4399;
  else if( aa == "KAA")  
    return 4400;
  else if( aa == "KAR")  
    return 4401;
  else if( aa == "KAN")  
    return 4402;
  else if( aa == "KAD")  
    return 4403;
  else if( aa == "KAC")  
    return 4404;
  else if( aa == "KAQ")  
    return 4405;
  else if( aa == "KAE")  
    return 4406;
  else if( aa == "KAG")  
    return 4407;
  else if( aa == "KAH")  
    return 4408;
  else if( aa == "KAI")  
    return 4409;
  else if( aa == "KAL")  
    return 4410;
  else if( aa == "KAK")  
    return 4411;
  else if( aa == "KAM")  
    return 4412;
  else if( aa == "KAF")  
    return 4413;
  else if( aa == "KAP")  
    return 4414;
  else if( aa == "KAS")  
    return 4415;
  else if( aa == "KAT")  
    return 4416;
  else if( aa == "KAW")  
    return 4417;
  else if( aa == "KAY")  
    return 4418;
  else if( aa == "KAV")  
    return 4419;
  else if( aa == "KRA")  
    return 4420;
  else if( aa == "KRR")  
    return 4421;
  else if( aa == "KRN")  
    return 4422;
  else if( aa == "KRD")  
    return 4423;
  else if( aa == "KRC")  
    return 4424;
  else if( aa == "KRQ")  
    return 4425;
  else if( aa == "KRE")  
    return 4426;
  else if( aa == "KRG")  
    return 4427;
  else if( aa == "KRH")  
    return 4428;
  else if( aa == "KRI")  
    return 4429;
  else if( aa == "KRL")  
    return 4430;
  else if( aa == "KRK")  
    return 4431;
  else if( aa == "KRM")  
    return 4432;
  else if( aa == "KRF")  
    return 4433;
  else if( aa == "KRP")  
    return 4434;
  else if( aa == "KRS")  
    return 4435;
  else if( aa == "KRT")  
    return 4436;
  else if( aa == "KRW")  
    return 4437;
  else if( aa == "KRY")  
    return 4438;
  else if( aa == "KRV")  
    return 4439;
  else if( aa == "KNA")  
    return 4440;
  else if( aa == "KNR")  
    return 4441;
  else if( aa == "KNN")  
    return 4442;
  else if( aa == "KND")  
    return 4443;
  else if( aa == "KNC")  
    return 4444;
  else if( aa == "KNQ")  
    return 4445;
  else if( aa == "KNE")  
    return 4446;
  else if( aa == "KNG")  
    return 4447;
  else if( aa == "KNH")  
    return 4448;
  else if( aa == "KNI")  
    return 4449;
  else if( aa == "KNL")  
    return 4450;
  else if( aa == "KNK")  
    return 4451;
  else if( aa == "KNM")  
    return 4452;
  else if( aa == "KNF")  
    return 4453;
  else if( aa == "KNP")  
    return 4454;
  else if( aa == "KNS")  
    return 4455;
  else if( aa == "KNT")  
    return 4456;
  else if( aa == "KNW")  
    return 4457;
  else if( aa == "KNY")  
    return 4458;
  else if( aa == "KNV")  
    return 4459;
  else if( aa == "KDA")  
    return 4460;
  else if( aa == "KDR")  
    return 4461;
  else if( aa == "KDN")  
    return 4462;
  else if( aa == "KDD")  
    return 4463;
  else if( aa == "KDC")  
    return 4464;
  else if( aa == "KDQ")  
    return 4465;
  else if( aa == "KDE")  
    return 4466;
  else if( aa == "KDG")  
    return 4467;
  else if( aa == "KDH")  
    return 4468;
  else if( aa == "KDI")  
    return 4469;
  else if( aa == "KDL")  
    return 4470;
  else if( aa == "KDK")  
    return 4471;
  else if( aa == "KDM")  
    return 4472;
  else if( aa == "KDF")  
    return 4473;
  else if( aa == "KDP")  
    return 4474;
  else if( aa == "KDS")  
    return 4475;
  else if( aa == "KDT")  
    return 4476;
  else if( aa == "KDW")  
    return 4477;
  else if( aa == "KDY")  
    return 4478;
  else if( aa == "KDV")  
    return 4479;
  else if( aa == "KCA")  
    return 4480;
  else if( aa == "KCR")  
    return 4481;
  else if( aa == "KCN")  
    return 4482;
  else if( aa == "KCD")  
    return 4483;
  else if( aa == "KCC")  
    return 4484;
  else if( aa == "KCQ")  
    return 4485;
  else if( aa == "KCE")  
    return 4486;
  else if( aa == "KCG")  
    return 4487;
  else if( aa == "KCH")  
    return 4488;
  else if( aa == "KCI")  
    return 4489;
  else if( aa == "KCL")  
    return 4490;
  else if( aa == "KCK")  
    return 4491;
  else if( aa == "KCM")  
    return 4492;
  else if( aa == "KCF")  
    return 4493;
  else if( aa == "KCP")  
    return 4494;
  else if( aa == "KCS")  
    return 4495;
  else if( aa == "KCT")  
    return 4496;
  else if( aa == "KCW")  
    return 4497;
  else if( aa == "KCY")  
    return 4498;
  else if( aa == "KCV")  
    return 4499;
  else if( aa == "KQA")  
    return 4500;
  else if( aa == "KQR")  
    return 4501;
  else if( aa == "KQN")  
    return 4502;
  else if( aa == "KQD")  
    return 4503;
  else if( aa == "KQC")  
    return 4504;
  else if( aa == "KQQ")  
    return 4505;
  else if( aa == "KQE")  
    return 4506;
  else if( aa == "KQG")  
    return 4507;
  else if( aa == "KQH")  
    return 4508;
  else if( aa == "KQI")  
    return 4509;
  else if( aa == "KQL")  
    return 4510;
  else if( aa == "KQK")  
    return 4511;
  else if( aa == "KQM")  
    return 4512;
  else if( aa == "KQF")  
    return 4513;
  else if( aa == "KQP")  
    return 4514;
  else if( aa == "KQS")  
    return 4515;
  else if( aa == "KQT")  
    return 4516;
  else if( aa == "KQW")  
    return 4517;
  else if( aa == "KQY")  
    return 4518;
  else if( aa == "KQV")  
    return 4519;
  else if( aa == "KEA")  
    return 4520;
  else if( aa == "KER")  
    return 4521;
  else if( aa == "KEN")  
    return 4522;
  else if( aa == "KED")  
    return 4523;
  else if( aa == "KEC")  
    return 4524;
  else if( aa == "KEQ")  
    return 4525;
  else if( aa == "KEE")  
    return 4526;
  else if( aa == "KEG")  
    return 4527;
  else if( aa == "KEH")  
    return 4528;
  else if( aa == "KEI")  
    return 4529;
  else if( aa == "KEL")  
    return 4530;
  else if( aa == "KEK")  
    return 4531;
  else if( aa == "KEM")  
    return 4532;
  else if( aa == "KEF")  
    return 4533;
  else if( aa == "KEP")  
    return 4534;
  else if( aa == "KES")  
    return 4535;
  else if( aa == "KET")  
    return 4536;
  else if( aa == "KEW")  
    return 4537;
  else if( aa == "KEY")  
    return 4538;
  else if( aa == "KEV")  
    return 4539;
  else if( aa == "KGA")  
    return 4540;
  else if( aa == "KGR")  
    return 4541;
  else if( aa == "KGN")  
    return 4542;
  else if( aa == "KGD")  
    return 4543;
  else if( aa == "KGC")  
    return 4544;
  else if( aa == "KGQ")  
    return 4545;
  else if( aa == "KGE")  
    return 4546;
  else if( aa == "KGG")  
    return 4547;
  else if( aa == "KGH")  
    return 4548;
  else if( aa == "KGI")  
    return 4549;
  else if( aa == "KGL")  
    return 4550;
  else if( aa == "KGK")  
    return 4551;
  else if( aa == "KGM")  
    return 4552;
  else if( aa == "KGF")  
    return 4553;
  else if( aa == "KGP")  
    return 4554;
  else if( aa == "KGS")  
    return 4555;
  else if( aa == "KGT")  
    return 4556;
  else if( aa == "KGW")  
    return 4557;
  else if( aa == "KGY")  
    return 4558;
  else if( aa == "KGV")  
    return 4559;
  else if( aa == "KHA")  
    return 4560;
  else if( aa == "KHR")  
    return 4561;
  else if( aa == "KHN")  
    return 4562;
  else if( aa == "KHD")  
    return 4563;
  else if( aa == "KHC")  
    return 4564;
  else if( aa == "KHQ")  
    return 4565;
  else if( aa == "KHE")  
    return 4566;
  else if( aa == "KHG")  
    return 4567;
  else if( aa == "KHH")  
    return 4568;
  else if( aa == "KHI")  
    return 4569;
  else if( aa == "KHL")  
    return 4570;
  else if( aa == "KHK")  
    return 4571;
  else if( aa == "KHM")  
    return 4572;
  else if( aa == "KHF")  
    return 4573;
  else if( aa == "KHP")  
    return 4574;
  else if( aa == "KHS")  
    return 4575;
  else if( aa == "KHT")  
    return 4576;
  else if( aa == "KHW")  
    return 4577;
  else if( aa == "KHY")  
    return 4578;
  else if( aa == "KHV")  
    return 4579;
  else if( aa == "KIA")  
    return 4580;
  else if( aa == "KIR")  
    return 4581;
  else if( aa == "KIN")  
    return 4582;
  else if( aa == "KID")  
    return 4583;
  else if( aa == "KIC")  
    return 4584;
  else if( aa == "KIQ")  
    return 4585;
  else if( aa == "KIE")  
    return 4586;
  else if( aa == "KIG")  
    return 4587;
  else if( aa == "KIH")  
    return 4588;
  else if( aa == "KII")  
    return 4589;
  else if( aa == "KIL")  
    return 4590;
  else if( aa == "KIK")  
    return 4591;
  else if( aa == "KIM")  
    return 4592;
  else if( aa == "KIF")  
    return 4593;
  else if( aa == "KIP")  
    return 4594;
  else if( aa == "KIS")  
    return 4595;
  else if( aa == "KIT")  
    return 4596;
  else if( aa == "KIW")  
    return 4597;
  else if( aa == "KIY")  
    return 4598;
  else if( aa == "KIV")  
    return 4599;
  else if( aa == "KLA")  
    return 4600;
  else if( aa == "KLR")  
    return 4601;
  else if( aa == "KLN")  
    return 4602;
  else if( aa == "KLD")  
    return 4603;
  else if( aa == "KLC")  
    return 4604;
  else if( aa == "KLQ")  
    return 4605;
  else if( aa == "KLE")  
    return 4606;
  else if( aa == "KLG")  
    return 4607;
  else if( aa == "KLH")  
    return 4608;
  else if( aa == "KLI")  
    return 4609;
  else if( aa == "KLL")  
    return 4610;
  else if( aa == "KLK")  
    return 4611;
  else if( aa == "KLM")  
    return 4612;
  else if( aa == "KLF")  
    return 4613;
  else if( aa == "KLP")  
    return 4614;
  else if( aa == "KLS")  
    return 4615;
  else if( aa == "KLT")  
    return 4616;
  else if( aa == "KLW")  
    return 4617;
  else if( aa == "KLY")  
    return 4618;
  else if( aa == "KLV")  
    return 4619;
  else if( aa == "KKA")  
    return 4620;
  else if( aa == "KKR")  
    return 4621;
  else if( aa == "KKN")  
    return 4622;
  else if( aa == "KKD")  
    return 4623;
  else if( aa == "KKC")  
    return 4624;
  else if( aa == "KKQ")  
    return 4625;
  else if( aa == "KKE")  
    return 4626;
  else if( aa == "KKG")  
    return 4627;
  else if( aa == "KKH")  
    return 4628;
  else if( aa == "KKI")  
    return 4629;
  else if( aa == "KKL")  
    return 4630;
  else if( aa == "KKK")  
    return 4631;
  else if( aa == "KKM")  
    return 4632;
  else if( aa == "KKF")  
    return 4633;
  else if( aa == "KKP")  
    return 4634;
  else if( aa == "KKS")  
    return 4635;
  else if( aa == "KKT")  
    return 4636;
  else if( aa == "KKW")  
    return 4637;
  else if( aa == "KKY")  
    return 4638;
  else if( aa == "KKV")  
    return 4639;
  else if( aa == "KMA")  
    return 4640;
  else if( aa == "KMR")  
    return 4641;
  else if( aa == "KMN")  
    return 4642;
  else if( aa == "KMD")  
    return 4643;
  else if( aa == "KMC")  
    return 4644;
  else if( aa == "KMQ")  
    return 4645;
  else if( aa == "KME")  
    return 4646;
  else if( aa == "KMG")  
    return 4647;
  else if( aa == "KMH")  
    return 4648;
  else if( aa == "KMI")  
    return 4649;
  else if( aa == "KML")  
    return 4650;
  else if( aa == "KMK")  
    return 4651;
  else if( aa == "KMM")  
    return 4652;
  else if( aa == "KMF")  
    return 4653;
  else if( aa == "KMP")  
    return 4654;
  else if( aa == "KMS")  
    return 4655;
  else if( aa == "KMT")  
    return 4656;
  else if( aa == "KMW")  
    return 4657;
  else if( aa == "KMY")  
    return 4658;
  else if( aa == "KMV")  
    return 4659;
  else if( aa == "KFA")  
    return 4660;
  else if( aa == "KFR")  
    return 4661;
  else if( aa == "KFN")  
    return 4662;
  else if( aa == "KFD")  
    return 4663;
  else if( aa == "KFC")  
    return 4664;
  else if( aa == "KFQ")  
    return 4665;
  else if( aa == "KFE")  
    return 4666;
  else if( aa == "KFG")  
    return 4667;
  else if( aa == "KFH")  
    return 4668;
  else if( aa == "KFI")  
    return 4669;
  else if( aa == "KFL")  
    return 4670;
  else if( aa == "KFK")  
    return 4671;
  else if( aa == "KFM")  
    return 4672;
  else if( aa == "KFF")  
    return 4673;
  else if( aa == "KFP")  
    return 4674;
  else if( aa == "KFS")  
    return 4675;
  else if( aa == "KFT")  
    return 4676;
  else if( aa == "KFW")  
    return 4677;
  else if( aa == "KFY")  
    return 4678;
  else if( aa == "KFV")  
    return 4679;
  else if( aa == "KPA")  
    return 4680;
  else if( aa == "KPR")  
    return 4681;
  else if( aa == "KPN")  
    return 4682;
  else if( aa == "KPD")  
    return 4683;
  else if( aa == "KPC")  
    return 4684;
  else if( aa == "KPQ")  
    return 4685;
  else if( aa == "KPE")  
    return 4686;
  else if( aa == "KPG")  
    return 4687;
  else if( aa == "KPH")  
    return 4688;
  else if( aa == "KPI")  
    return 4689;
  else if( aa == "KPL")  
    return 4690;
  else if( aa == "KPK")  
    return 4691;
  else if( aa == "KPM")  
    return 4692;
  else if( aa == "KPF")  
    return 4693;
  else if( aa == "KPP")  
    return 4694;
  else if( aa == "KPS")  
    return 4695;
  else if( aa == "KPT")  
    return 4696;
  else if( aa == "KPW")  
    return 4697;
  else if( aa == "KPY")  
    return 4698;
  else if( aa == "KPV")  
    return 4699;
  else if( aa == "KSA")  
    return 4700;
  else if( aa == "KSR")  
    return 4701;
  else if( aa == "KSN")  
    return 4702;
  else if( aa == "KSD")  
    return 4703;
  else if( aa == "KSC")  
    return 4704;
  else if( aa == "KSQ")  
    return 4705;
  else if( aa == "KSE")  
    return 4706;
  else if( aa == "KSG")  
    return 4707;
  else if( aa == "KSH")  
    return 4708;
  else if( aa == "KSI")  
    return 4709;
  else if( aa == "KSL")  
    return 4710;
  else if( aa == "KSK")  
    return 4711;
  else if( aa == "KSM")  
    return 4712;
  else if( aa == "KSF")  
    return 4713;
  else if( aa == "KSP")  
    return 4714;
  else if( aa == "KSS")  
    return 4715;
  else if( aa == "KST")  
    return 4716;
  else if( aa == "KSW")  
    return 4717;
  else if( aa == "KSY")  
    return 4718;
  else if( aa == "KSV")  
    return 4719;
  else if( aa == "KTA")  
    return 4720;
  else if( aa == "KTR")  
    return 4721;
  else if( aa == "KTN")  
    return 4722;
  else if( aa == "KTD")  
    return 4723;
  else if( aa == "KTC")  
    return 4724;
  else if( aa == "KTQ")  
    return 4725;
  else if( aa == "KTE")  
    return 4726;
  else if( aa == "KTG")  
    return 4727;
  else if( aa == "KTH")  
    return 4728;
  else if( aa == "KTI")  
    return 4729;
  else if( aa == "KTL")  
    return 4730;
  else if( aa == "KTK")  
    return 4731;
  else if( aa == "KTM")  
    return 4732;
  else if( aa == "KTF")  
    return 4733;
  else if( aa == "KTP")  
    return 4734;
  else if( aa == "KTS")  
    return 4735;
  else if( aa == "KTT")  
    return 4736;
  else if( aa == "KTW")  
    return 4737;
  else if( aa == "KTY")  
    return 4738;
  else if( aa == "KTV")  
    return 4739;
  else if( aa == "KWA")  
    return 4740;
  else if( aa == "KWR")  
    return 4741;
  else if( aa == "KWN")  
    return 4742;
  else if( aa == "KWD")  
    return 4743;
  else if( aa == "KWC")  
    return 4744;
  else if( aa == "KWQ")  
    return 4745;
  else if( aa == "KWE")  
    return 4746;
  else if( aa == "KWG")  
    return 4747;
  else if( aa == "KWH")  
    return 4748;
  else if( aa == "KWI")  
    return 4749;
  else if( aa == "KWL")  
    return 4750;
  else if( aa == "KWK")  
    return 4751;
  else if( aa == "KWM")  
    return 4752;
  else if( aa == "KWF")  
    return 4753;
  else if( aa == "KWP")  
    return 4754;
  else if( aa == "KWS")  
    return 4755;
  else if( aa == "KWT")  
    return 4756;
  else if( aa == "KWW")  
    return 4757;
  else if( aa == "KWY")  
    return 4758;
  else if( aa == "KWV")  
    return 4759;
  else if( aa == "KYA")  
    return 4760;
  else if( aa == "KYR")  
    return 4761;
  else if( aa == "KYN")  
    return 4762;
  else if( aa == "KYD")  
    return 4763;
  else if( aa == "KYC")  
    return 4764;
  else if( aa == "KYQ")  
    return 4765;
  else if( aa == "KYE")  
    return 4766;
  else if( aa == "KYG")  
    return 4767;
  else if( aa == "KYH")  
    return 4768;
  else if( aa == "KYI")  
    return 4769;
  else if( aa == "KYL")  
    return 4770;
  else if( aa == "KYK")  
    return 4771;
  else if( aa == "KYM")  
    return 4772;
  else if( aa == "KYF")  
    return 4773;
  else if( aa == "KYP")  
    return 4774;
  else if( aa == "KYS")  
    return 4775;
  else if( aa == "KYT")  
    return 4776;
  else if( aa == "KYW")  
    return 4777;
  else if( aa == "KYY")  
    return 4778;
  else if( aa == "KYV")  
    return 4779;
  else if( aa == "KVA")  
    return 4780;
  else if( aa == "KVR")  
    return 4781;
  else if( aa == "KVN")  
    return 4782;
  else if( aa == "KVD")  
    return 4783;
  else if( aa == "KVC")  
    return 4784;
  else if( aa == "KVQ")  
    return 4785;
  else if( aa == "KVE")  
    return 4786;
  else if( aa == "KVG")  
    return 4787;
  else if( aa == "KVH")  
    return 4788;
  else if( aa == "KVI")  
    return 4789;
  else if( aa == "KVL")  
    return 4790;
  else if( aa == "KVK")  
    return 4791;
  else if( aa == "KVM")  
    return 4792;
  else if( aa == "KVF")  
    return 4793;
  else if( aa == "KVP")  
    return 4794;
  else if( aa == "KVS")  
    return 4795;
  else if( aa == "KVT")  
    return 4796;
  else if( aa == "KVW")  
    return 4797;
  else if( aa == "KVY")  
    return 4798;
  else if( aa == "KVV")  
    return 4799;
  else if( aa == "MAA")  
    return 4800;
  else if( aa == "MAR")  
    return 4801;
  else if( aa == "MAN")  
    return 4802;
  else if( aa == "MAD")  
    return 4803;
  else if( aa == "MAC")  
    return 4804;
  else if( aa == "MAQ")  
    return 4805;
  else if( aa == "MAE")  
    return 4806;
  else if( aa == "MAG")  
    return 4807;
  else if( aa == "MAH")  
    return 4808;
  else if( aa == "MAI")  
    return 4809;
  else if( aa == "MAL")  
    return 4810;
  else if( aa == "MAK")  
    return 4811;
  else if( aa == "MAM")  
    return 4812;
  else if( aa == "MAF")  
    return 4813;
  else if( aa == "MAP")  
    return 4814;
  else if( aa == "MAS")  
    return 4815;
  else if( aa == "MAT")  
    return 4816;
  else if( aa == "MAW")  
    return 4817;
  else if( aa == "MAY")  
    return 4818;
  else if( aa == "MAV")  
    return 4819;
  else if( aa == "MRA")  
    return 4820;
  else if( aa == "MRR")  
    return 4821;
  else if( aa == "MRN")  
    return 4822;
  else if( aa == "MRD")  
    return 4823;
  else if( aa == "MRC")  
    return 4824;
  else if( aa == "MRQ")  
    return 4825;
  else if( aa == "MRE")  
    return 4826;
  else if( aa == "MRG")  
    return 4827;
  else if( aa == "MRH")  
    return 4828;
  else if( aa == "MRI")  
    return 4829;
  else if( aa == "MRL")  
    return 4830;
  else if( aa == "MRK")  
    return 4831;
  else if( aa == "MRM")  
    return 4832;
  else if( aa == "MRF")  
    return 4833;
  else if( aa == "MRP")  
    return 4834;
  else if( aa == "MRS")  
    return 4835;
  else if( aa == "MRT")  
    return 4836;
  else if( aa == "MRW")  
    return 4837;
  else if( aa == "MRY")  
    return 4838;
  else if( aa == "MRV")  
    return 4839;
  else if( aa == "MNA")  
    return 4840;
  else if( aa == "MNR")  
    return 4841;
  else if( aa == "MNN")  
    return 4842;
  else if( aa == "MND")  
    return 4843;
  else if( aa == "MNC")  
    return 4844;
  else if( aa == "MNQ")  
    return 4845;
  else if( aa == "MNE")  
    return 4846;
  else if( aa == "MNG")  
    return 4847;
  else if( aa == "MNH")  
    return 4848;
  else if( aa == "MNI")  
    return 4849;
  else if( aa == "MNL")  
    return 4850;
  else if( aa == "MNK")  
    return 4851;
  else if( aa == "MNM")  
    return 4852;
  else if( aa == "MNF")  
    return 4853;
  else if( aa == "MNP")  
    return 4854;
  else if( aa == "MNS")  
    return 4855;
  else if( aa == "MNT")  
    return 4856;
  else if( aa == "MNW")  
    return 4857;
  else if( aa == "MNY")  
    return 4858;
  else if( aa == "MNV")  
    return 4859;
  else if( aa == "MDA")  
    return 4860;
  else if( aa == "MDR")  
    return 4861;
  else if( aa == "MDN")  
    return 4862;
  else if( aa == "MDD")  
    return 4863;
  else if( aa == "MDC")  
    return 4864;
  else if( aa == "MDQ")  
    return 4865;
  else if( aa == "MDE")  
    return 4866;
  else if( aa == "MDG")  
    return 4867;
  else if( aa == "MDH")  
    return 4868;
  else if( aa == "MDI")  
    return 4869;
  else if( aa == "MDL")  
    return 4870;
  else if( aa == "MDK")  
    return 4871;
  else if( aa == "MDM")  
    return 4872;
  else if( aa == "MDF")  
    return 4873;
  else if( aa == "MDP")  
    return 4874;
  else if( aa == "MDS")  
    return 4875;
  else if( aa == "MDT")  
    return 4876;
  else if( aa == "MDW")  
    return 4877;
  else if( aa == "MDY")  
    return 4878;
  else if( aa == "MDV")  
    return 4879;
  else if( aa == "MCA")  
    return 4880;
  else if( aa == "MCR")  
    return 4881;
  else if( aa == "MCN")  
    return 4882;
  else if( aa == "MCD")  
    return 4883;
  else if( aa == "MCC")  
    return 4884;
  else if( aa == "MCQ")  
    return 4885;
  else if( aa == "MCE")  
    return 4886;
  else if( aa == "MCG")  
    return 4887;
  else if( aa == "MCH")  
    return 4888;
  else if( aa == "MCI")  
    return 4889;
  else if( aa == "MCL")  
    return 4890;
  else if( aa == "MCK")  
    return 4891;
  else if( aa == "MCM")  
    return 4892;
  else if( aa == "MCF")  
    return 4893;
  else if( aa == "MCP")  
    return 4894;
  else if( aa == "MCS")  
    return 4895;
  else if( aa == "MCT")  
    return 4896;
  else if( aa == "MCW")  
    return 4897;
  else if( aa == "MCY")  
    return 4898;
  else if( aa == "MCV")  
    return 4899;
  else if( aa == "MQA")  
    return 4900;
  else if( aa == "MQR")  
    return 4901;
  else if( aa == "MQN")  
    return 4902;
  else if( aa == "MQD")  
    return 4903;
  else if( aa == "MQC")  
    return 4904;
  else if( aa == "MQQ")  
    return 4905;
  else if( aa == "MQE")  
    return 4906;
  else if( aa == "MQG")  
    return 4907;
  else if( aa == "MQH")  
    return 4908;
  else if( aa == "MQI")  
    return 4909;
  else if( aa == "MQL")  
    return 4910;
  else if( aa == "MQK")  
    return 4911;
  else if( aa == "MQM")  
    return 4912;
  else if( aa == "MQF")  
    return 4913;
  else if( aa == "MQP")  
    return 4914;
  else if( aa == "MQS")  
    return 4915;
  else if( aa == "MQT")  
    return 4916;
  else if( aa == "MQW")  
    return 4917;
  else if( aa == "MQY")  
    return 4918;
  else if( aa == "MQV")  
    return 4919;
  else if( aa == "MEA")  
    return 4920;
  else if( aa == "MER")  
    return 4921;
  else if( aa == "MEN")  
    return 4922;
  else if( aa == "MED")  
    return 4923;
  else if( aa == "MEC")  
    return 4924;
  else if( aa == "MEQ")  
    return 4925;
  else if( aa == "MEE")  
    return 4926;
  else if( aa == "MEG")  
    return 4927;
  else if( aa == "MEH")  
    return 4928;
  else if( aa == "MEI")  
    return 4929;
  else if( aa == "MEL")  
    return 4930;
  else if( aa == "MEK")  
    return 4931;
  else if( aa == "MEM")  
    return 4932;
  else if( aa == "MEF")  
    return 4933;
  else if( aa == "MEP")  
    return 4934;
  else if( aa == "MES")  
    return 4935;
  else if( aa == "MET")  
    return 4936;
  else if( aa == "MEW")  
    return 4937;
  else if( aa == "MEY")  
    return 4938;
  else if( aa == "MEV")  
    return 4939;
  else if( aa == "MGA")  
    return 4940;
  else if( aa == "MGR")  
    return 4941;
  else if( aa == "MGN")  
    return 4942;
  else if( aa == "MGD")  
    return 4943;
  else if( aa == "MGC")  
    return 4944;
  else if( aa == "MGQ")  
    return 4945;
  else if( aa == "MGE")  
    return 4946;
  else if( aa == "MGG")  
    return 4947;
  else if( aa == "MGH")  
    return 4948;
  else if( aa == "MGI")  
    return 4949;
  else if( aa == "MGL")  
    return 4950;
  else if( aa == "MGK")  
    return 4951;
  else if( aa == "MGM")  
    return 4952;
  else if( aa == "MGF")  
    return 4953;
  else if( aa == "MGP")  
    return 4954;
  else if( aa == "MGS")  
    return 4955;
  else if( aa == "MGT")  
    return 4956;
  else if( aa == "MGW")  
    return 4957;
  else if( aa == "MGY")  
    return 4958;
  else if( aa == "MGV")  
    return 4959;
  else if( aa == "MHA")  
    return 4960;
  else if( aa == "MHR")  
    return 4961;
  else if( aa == "MHN")  
    return 4962;
  else if( aa == "MHD")  
    return 4963;
  else if( aa == "MHC")  
    return 4964;
  else if( aa == "MHQ")  
    return 4965;
  else if( aa == "MHE")  
    return 4966;
  else if( aa == "MHG")  
    return 4967;
  else if( aa == "MHH")  
    return 4968;
  else if( aa == "MHI")  
    return 4969;
  else if( aa == "MHL")  
    return 4970;
  else if( aa == "MHK")  
    return 4971;
  else if( aa == "MHM")  
    return 4972;
  else if( aa == "MHF")  
    return 4973;
  else if( aa == "MHP")  
    return 4974;
  else if( aa == "MHS")  
    return 4975;
  else if( aa == "MHT")  
    return 4976;
  else if( aa == "MHW")  
    return 4977;
  else if( aa == "MHY")  
    return 4978;
  else if( aa == "MHV")  
    return 4979;
  else if( aa == "MIA")  
    return 4980;
  else if( aa == "MIR")  
    return 4981;
  else if( aa == "MIN")  
    return 4982;
  else if( aa == "MID")  
    return 4983;
  else if( aa == "MIC")  
    return 4984;
  else if( aa == "MIQ")  
    return 4985;
  else if( aa == "MIE")  
    return 4986;
  else if( aa == "MIG")  
    return 4987;
  else if( aa == "MIH")  
    return 4988;
  else if( aa == "MII")  
    return 4989;
  else if( aa == "MIL")  
    return 4990;
  else if( aa == "MIK")  
    return 4991;
  else if( aa == "MIM")  
    return 4992;
  else if( aa == "MIF")  
    return 4993;
  else if( aa == "MIP")  
    return 4994;
  else if( aa == "MIS")  
    return 4995;
  else if( aa == "MIT")  
    return 4996;
  else if( aa == "MIW")  
    return 4997;
  else if( aa == "MIY")  
    return 4998;
  else if( aa == "MIV")  
    return 4999;
  else if( aa == "MLA")  
    return 5000;
  else if( aa == "MLR")  
    return 5001;
  else if( aa == "MLN")  
    return 5002;
  else if( aa == "MLD")  
    return 5003;
  else if( aa == "MLC")  
    return 5004;
  else if( aa == "MLQ")  
    return 5005;
  else if( aa == "MLE")  
    return 5006;
  else if( aa == "MLG")  
    return 5007;
  else if( aa == "MLH")  
    return 5008;
  else if( aa == "MLI")  
    return 5009;
  else if( aa == "MLL")  
    return 5010;
  else if( aa == "MLK")  
    return 5011;
  else if( aa == "MLM")  
    return 5012;
  else if( aa == "MLF")  
    return 5013;
  else if( aa == "MLP")  
    return 5014;
  else if( aa == "MLS")  
    return 5015;
  else if( aa == "MLT")  
    return 5016;
  else if( aa == "MLW")  
    return 5017;
  else if( aa == "MLY")  
    return 5018;
  else if( aa == "MLV")  
    return 5019;
  else if( aa == "MKA")  
    return 5020;
  else if( aa == "MKR")  
    return 5021;
  else if( aa == "MKN")  
    return 5022;
  else if( aa == "MKD")  
    return 5023;
  else if( aa == "MKC")  
    return 5024;
  else if( aa == "MKQ")  
    return 5025;
  else if( aa == "MKE")  
    return 5026;
  else if( aa == "MKG")  
    return 5027;
  else if( aa == "MKH")  
    return 5028;
  else if( aa == "MKI")  
    return 5029;
  else if( aa == "MKL")  
    return 5030;
  else if( aa == "MKK")  
    return 5031;
  else if( aa == "MKM")  
    return 5032;
  else if( aa == "MKF")  
    return 5033;
  else if( aa == "MKP")  
    return 5034;
  else if( aa == "MKS")  
    return 5035;
  else if( aa == "MKT")  
    return 5036;
  else if( aa == "MKW")  
    return 5037;
  else if( aa == "MKY")  
    return 5038;
  else if( aa == "MKV")  
    return 5039;
  else if( aa == "MMA")  
    return 5040;
  else if( aa == "MMR")  
    return 5041;
  else if( aa == "MMN")  
    return 5042;
  else if( aa == "MMD")  
    return 5043;
  else if( aa == "MMC")  
    return 5044;
  else if( aa == "MMQ")  
    return 5045;
  else if( aa == "MME")  
    return 5046;
  else if( aa == "MMG")  
    return 5047;
  else if( aa == "MMH")  
    return 5048;
  else if( aa == "MMI")  
    return 5049;
  else if( aa == "MML")  
    return 5050;
  else if( aa == "MMK")  
    return 5051;
  else if( aa == "MMM")  
    return 5052;
  else if( aa == "MMF")  
    return 5053;
  else if( aa == "MMP")  
    return 5054;
  else if( aa == "MMS")  
    return 5055;
  else if( aa == "MMT")  
    return 5056;
  else if( aa == "MMW")  
    return 5057;
  else if( aa == "MMY")  
    return 5058;
  else if( aa == "MMV")  
    return 5059;
  else if( aa == "MFA")  
    return 5060;
  else if( aa == "MFR")  
    return 5061;
  else if( aa == "MFN")  
    return 5062;
  else if( aa == "MFD")  
    return 5063;
  else if( aa == "MFC")  
    return 5064;
  else if( aa == "MFQ")  
    return 5065;
  else if( aa == "MFE")  
    return 5066;
  else if( aa == "MFG")  
    return 5067;
  else if( aa == "MFH")  
    return 5068;
  else if( aa == "MFI")  
    return 5069;
  else if( aa == "MFL")  
    return 5070;
  else if( aa == "MFK")  
    return 5071;
  else if( aa == "MFM")  
    return 5072;
  else if( aa == "MFF")  
    return 5073;
  else if( aa == "MFP")  
    return 5074;
  else if( aa == "MFS")  
    return 5075;
  else if( aa == "MFT")  
    return 5076;
  else if( aa == "MFW")  
    return 5077;
  else if( aa == "MFY")  
    return 5078;
  else if( aa == "MFV")  
    return 5079;
  else if( aa == "MPA")  
    return 5080;
  else if( aa == "MPR")  
    return 5081;
  else if( aa == "MPN")  
    return 5082;
  else if( aa == "MPD")  
    return 5083;
  else if( aa == "MPC")  
    return 5084;
  else if( aa == "MPQ")  
    return 5085;
  else if( aa == "MPE")  
    return 5086;
  else if( aa == "MPG")  
    return 5087;
  else if( aa == "MPH")  
    return 5088;
  else if( aa == "MPI")  
    return 5089;
  else if( aa == "MPL")  
    return 5090;
  else if( aa == "MPK")  
    return 5091;
  else if( aa == "MPM")  
    return 5092;
  else if( aa == "MPF")  
    return 5093;
  else if( aa == "MPP")  
    return 5094;
  else if( aa == "MPS")  
    return 5095;
  else if( aa == "MPT")  
    return 5096;
  else if( aa == "MPW")  
    return 5097;
  else if( aa == "MPY")  
    return 5098;
  else if( aa == "MPV")  
    return 5099;
  else if( aa == "MSA")  
    return 5100;
  else if( aa == "MSR")  
    return 5101;
  else if( aa == "MSN")  
    return 5102;
  else if( aa == "MSD")  
    return 5103;
  else if( aa == "MSC")  
    return 5104;
  else if( aa == "MSQ")  
    return 5105;
  else if( aa == "MSE")  
    return 5106;
  else if( aa == "MSG")  
    return 5107;
  else if( aa == "MSH")  
    return 5108;
  else if( aa == "MSI")  
    return 5109;
  else if( aa == "MSL")  
    return 5110;
  else if( aa == "MSK")  
    return 5111;
  else if( aa == "MSM")  
    return 5112;
  else if( aa == "MSF")  
    return 5113;
  else if( aa == "MSP")  
    return 5114;
  else if( aa == "MSS")  
    return 5115;
  else if( aa == "MST")  
    return 5116;
  else if( aa == "MSW")  
    return 5117;
  else if( aa == "MSY")  
    return 5118;
  else if( aa == "MSV")  
    return 5119;
  else if( aa == "MTA")  
    return 5120;
  else if( aa == "MTR")  
    return 5121;
  else if( aa == "MTN")  
    return 5122;
  else if( aa == "MTD")  
    return 5123;
  else if( aa == "MTC")  
    return 5124;
  else if( aa == "MTQ")  
    return 5125;
  else if( aa == "MTE")  
    return 5126;
  else if( aa == "MTG")  
    return 5127;
  else if( aa == "MTH")  
    return 5128;
  else if( aa == "MTI")  
    return 5129;
  else if( aa == "MTL")  
    return 5130;
  else if( aa == "MTK")  
    return 5131;
  else if( aa == "MTM")  
    return 5132;
  else if( aa == "MTF")  
    return 5133;
  else if( aa == "MTP")  
    return 5134;
  else if( aa == "MTS")  
    return 5135;
  else if( aa == "MTT")  
    return 5136;
  else if( aa == "MTW")  
    return 5137;
  else if( aa == "MTY")  
    return 5138;
  else if( aa == "MTV")  
    return 5139;
  else if( aa == "MWA")  
    return 5140;
  else if( aa == "MWR")  
    return 5141;
  else if( aa == "MWN")  
    return 5142;
  else if( aa == "MWD")  
    return 5143;
  else if( aa == "MWC")  
    return 5144;
  else if( aa == "MWQ")  
    return 5145;
  else if( aa == "MWE")  
    return 5146;
  else if( aa == "MWG")  
    return 5147;
  else if( aa == "MWH")  
    return 5148;
  else if( aa == "MWI")  
    return 5149;
  else if( aa == "MWL")  
    return 5150;
  else if( aa == "MWK")  
    return 5151;
  else if( aa == "MWM")  
    return 5152;
  else if( aa == "MWF")  
    return 5153;
  else if( aa == "MWP")  
    return 5154;
  else if( aa == "MWS")  
    return 5155;
  else if( aa == "MWT")  
    return 5156;
  else if( aa == "MWW")  
    return 5157;
  else if( aa == "MWY")  
    return 5158;
  else if( aa == "MWV")  
    return 5159;
  else if( aa == "MYA")  
    return 5160;
  else if( aa == "MYR")  
    return 5161;
  else if( aa == "MYN")  
    return 5162;
  else if( aa == "MYD")  
    return 5163;
  else if( aa == "MYC")  
    return 5164;
  else if( aa == "MYQ")  
    return 5165;
  else if( aa == "MYE")  
    return 5166;
  else if( aa == "MYG")  
    return 5167;
  else if( aa == "MYH")  
    return 5168;
  else if( aa == "MYI")  
    return 5169;
  else if( aa == "MYL")  
    return 5170;
  else if( aa == "MYK")  
    return 5171;
  else if( aa == "MYM")  
    return 5172;
  else if( aa == "MYF")  
    return 5173;
  else if( aa == "MYP")  
    return 5174;
  else if( aa == "MYS")  
    return 5175;
  else if( aa == "MYT")  
    return 5176;
  else if( aa == "MYW")  
    return 5177;
  else if( aa == "MYY")  
    return 5178;
  else if( aa == "MYV")  
    return 5179;
  else if( aa == "MVA")  
    return 5180;
  else if( aa == "MVR")  
    return 5181;
  else if( aa == "MVN")  
    return 5182;
  else if( aa == "MVD")  
    return 5183;
  else if( aa == "MVC")  
    return 5184;
  else if( aa == "MVQ")  
    return 5185;
  else if( aa == "MVE")  
    return 5186;
  else if( aa == "MVG")  
    return 5187;
  else if( aa == "MVH")  
    return 5188;
  else if( aa == "MVI")  
    return 5189;
  else if( aa == "MVL")  
    return 5190;
  else if( aa == "MVK")  
    return 5191;
  else if( aa == "MVM")  
    return 5192;
  else if( aa == "MVF")  
    return 5193;
  else if( aa == "MVP")  
    return 5194;
  else if( aa == "MVS")  
    return 5195;
  else if( aa == "MVT")  
    return 5196;
  else if( aa == "MVW")  
    return 5197;
  else if( aa == "MVY")  
    return 5198;
  else if( aa == "MVV")  
    return 5199;
  else if( aa == "FAA")  
    return 5200;
  else if( aa == "FAR")  
    return 5201;
  else if( aa == "FAN")  
    return 5202;
  else if( aa == "FAD")  
    return 5203;
  else if( aa == "FAC")  
    return 5204;
  else if( aa == "FAQ")  
    return 5205;
  else if( aa == "FAE")  
    return 5206;
  else if( aa == "FAG")  
    return 5207;
  else if( aa == "FAH")  
    return 5208;
  else if( aa == "FAI")  
    return 5209;
  else if( aa == "FAL")  
    return 5210;
  else if( aa == "FAK")  
    return 5211;
  else if( aa == "FAM")  
    return 5212;
  else if( aa == "FAF")  
    return 5213;
  else if( aa == "FAP")  
    return 5214;
  else if( aa == "FAS")  
    return 5215;
  else if( aa == "FAT")  
    return 5216;
  else if( aa == "FAW")  
    return 5217;
  else if( aa == "FAY")  
    return 5218;
  else if( aa == "FAV")  
    return 5219;
  else if( aa == "FRA")  
    return 5220;
  else if( aa == "FRR")  
    return 5221;
  else if( aa == "FRN")  
    return 5222;
  else if( aa == "FRD")  
    return 5223;
  else if( aa == "FRC")  
    return 5224;
  else if( aa == "FRQ")  
    return 5225;
  else if( aa == "FRE")  
    return 5226;
  else if( aa == "FRG")  
    return 5227;
  else if( aa == "FRH")  
    return 5228;
  else if( aa == "FRI")  
    return 5229;
  else if( aa == "FRL")  
    return 5230;
  else if( aa == "FRK")  
    return 5231;
  else if( aa == "FRM")  
    return 5232;
  else if( aa == "FRF")  
    return 5233;
  else if( aa == "FRP")  
    return 5234;
  else if( aa == "FRS")  
    return 5235;
  else if( aa == "FRT")  
    return 5236;
  else if( aa == "FRW")  
    return 5237;
  else if( aa == "FRY")  
    return 5238;
  else if( aa == "FRV")  
    return 5239;
  else if( aa == "FNA")  
    return 5240;
  else if( aa == "FNR")  
    return 5241;
  else if( aa == "FNN")  
    return 5242;
  else if( aa == "FND")  
    return 5243;
  else if( aa == "FNC")  
    return 5244;
  else if( aa == "FNQ")  
    return 5245;
  else if( aa == "FNE")  
    return 5246;
  else if( aa == "FNG")  
    return 5247;
  else if( aa == "FNH")  
    return 5248;
  else if( aa == "FNI")  
    return 5249;
  else if( aa == "FNL")  
    return 5250;
  else if( aa == "FNK")  
    return 5251;
  else if( aa == "FNM")  
    return 5252;
  else if( aa == "FNF")  
    return 5253;
  else if( aa == "FNP")  
    return 5254;
  else if( aa == "FNS")  
    return 5255;
  else if( aa == "FNT")  
    return 5256;
  else if( aa == "FNW")  
    return 5257;
  else if( aa == "FNY")  
    return 5258;
  else if( aa == "FNV")  
    return 5259;
  else if( aa == "FDA")  
    return 5260;
  else if( aa == "FDR")  
    return 5261;
  else if( aa == "FDN")  
    return 5262;
  else if( aa == "FDD")  
    return 5263;
  else if( aa == "FDC")  
    return 5264;
  else if( aa == "FDQ")  
    return 5265;
  else if( aa == "FDE")  
    return 5266;
  else if( aa == "FDG")  
    return 5267;
  else if( aa == "FDH")  
    return 5268;
  else if( aa == "FDI")  
    return 5269;
  else if( aa == "FDL")  
    return 5270;
  else if( aa == "FDK")  
    return 5271;
  else if( aa == "FDM")  
    return 5272;
  else if( aa == "FDF")  
    return 5273;
  else if( aa == "FDP")  
    return 5274;
  else if( aa == "FDS")  
    return 5275;
  else if( aa == "FDT")  
    return 5276;
  else if( aa == "FDW")  
    return 5277;
  else if( aa == "FDY")  
    return 5278;
  else if( aa == "FDV")  
    return 5279;
  else if( aa == "FCA")  
    return 5280;
  else if( aa == "FCR")  
    return 5281;
  else if( aa == "FCN")  
    return 5282;
  else if( aa == "FCD")  
    return 5283;
  else if( aa == "FCC")  
    return 5284;
  else if( aa == "FCQ")  
    return 5285;
  else if( aa == "FCE")  
    return 5286;
  else if( aa == "FCG")  
    return 5287;
  else if( aa == "FCH")  
    return 5288;
  else if( aa == "FCI")  
    return 5289;
  else if( aa == "FCL")  
    return 5290;
  else if( aa == "FCK")  
    return 5291;
  else if( aa == "FCM")  
    return 5292;
  else if( aa == "FCF")  
    return 5293;
  else if( aa == "FCP")  
    return 5294;
  else if( aa == "FCS")  
    return 5295;
  else if( aa == "FCT")  
    return 5296;
  else if( aa == "FCW")  
    return 5297;
  else if( aa == "FCY")  
    return 5298;
  else if( aa == "FCV")  
    return 5299;
  else if( aa == "FQA")  
    return 5300;
  else if( aa == "FQR")  
    return 5301;
  else if( aa == "FQN")  
    return 5302;
  else if( aa == "FQD")  
    return 5303;
  else if( aa == "FQC")  
    return 5304;
  else if( aa == "FQQ")  
    return 5305;
  else if( aa == "FQE")  
    return 5306;
  else if( aa == "FQG")  
    return 5307;
  else if( aa == "FQH")  
    return 5308;
  else if( aa == "FQI")  
    return 5309;
  else if( aa == "FQL")  
    return 5310;
  else if( aa == "FQK")  
    return 5311;
  else if( aa == "FQM")  
    return 5312;
  else if( aa == "FQF")  
    return 5313;
  else if( aa == "FQP")  
    return 5314;
  else if( aa == "FQS")  
    return 5315;
  else if( aa == "FQT")  
    return 5316;
  else if( aa == "FQW")  
    return 5317;
  else if( aa == "FQY")  
    return 5318;
  else if( aa == "FQV")  
    return 5319;
  else if( aa == "FEA")  
    return 5320;
  else if( aa == "FER")  
    return 5321;
  else if( aa == "FEN")  
    return 5322;
  else if( aa == "FED")  
    return 5323;
  else if( aa == "FEC")  
    return 5324;
  else if( aa == "FEQ")  
    return 5325;
  else if( aa == "FEE")  
    return 5326;
  else if( aa == "FEG")  
    return 5327;
  else if( aa == "FEH")  
    return 5328;
  else if( aa == "FEI")  
    return 5329;
  else if( aa == "FEL")  
    return 5330;
  else if( aa == "FEK")  
    return 5331;
  else if( aa == "FEM")  
    return 5332;
  else if( aa == "FEF")  
    return 5333;
  else if( aa == "FEP")  
    return 5334;
  else if( aa == "FES")  
    return 5335;
  else if( aa == "FET")  
    return 5336;
  else if( aa == "FEW")  
    return 5337;
  else if( aa == "FEY")  
    return 5338;
  else if( aa == "FEV")  
    return 5339;
  else if( aa == "FGA")  
    return 5340;
  else if( aa == "FGR")  
    return 5341;
  else if( aa == "FGN")  
    return 5342;
  else if( aa == "FGD")  
    return 5343;
  else if( aa == "FGC")  
    return 5344;
  else if( aa == "FGQ")  
    return 5345;
  else if( aa == "FGE")  
    return 5346;
  else if( aa == "FGG")  
    return 5347;
  else if( aa == "FGH")  
    return 5348;
  else if( aa == "FGI")  
    return 5349;
  else if( aa == "FGL")  
    return 5350;
  else if( aa == "FGK")  
    return 5351;
  else if( aa == "FGM")  
    return 5352;
  else if( aa == "FGF")  
    return 5353;
  else if( aa == "FGP")  
    return 5354;
  else if( aa == "FGS")  
    return 5355;
  else if( aa == "FGT")  
    return 5356;
  else if( aa == "FGW")  
    return 5357;
  else if( aa == "FGY")  
    return 5358;
  else if( aa == "FGV")  
    return 5359;
  else if( aa == "FHA")  
    return 5360;
  else if( aa == "FHR")  
    return 5361;
  else if( aa == "FHN")  
    return 5362;
  else if( aa == "FHD")  
    return 5363;
  else if( aa == "FHC")  
    return 5364;
  else if( aa == "FHQ")  
    return 5365;
  else if( aa == "FHE")  
    return 5366;
  else if( aa == "FHG")  
    return 5367;
  else if( aa == "FHH")  
    return 5368;
  else if( aa == "FHI")  
    return 5369;
  else if( aa == "FHL")  
    return 5370;
  else if( aa == "FHK")  
    return 5371;
  else if( aa == "FHM")  
    return 5372;
  else if( aa == "FHF")  
    return 5373;
  else if( aa == "FHP")  
    return 5374;
  else if( aa == "FHS")  
    return 5375;
  else if( aa == "FHT")  
    return 5376;
  else if( aa == "FHW")  
    return 5377;
  else if( aa == "FHY")  
    return 5378;
  else if( aa == "FHV")  
    return 5379;
  else if( aa == "FIA")  
    return 5380;
  else if( aa == "FIR")  
    return 5381;
  else if( aa == "FIN")  
    return 5382;
  else if( aa == "FID")  
    return 5383;
  else if( aa == "FIC")  
    return 5384;
  else if( aa == "FIQ")  
    return 5385;
  else if( aa == "FIE")  
    return 5386;
  else if( aa == "FIG")  
    return 5387;
  else if( aa == "FIH")  
    return 5388;
  else if( aa == "FII")  
    return 5389;
  else if( aa == "FIL")  
    return 5390;
  else if( aa == "FIK")  
    return 5391;
  else if( aa == "FIM")  
    return 5392;
  else if( aa == "FIF")  
    return 5393;
  else if( aa == "FIP")  
    return 5394;
  else if( aa == "FIS")  
    return 5395;
  else if( aa == "FIT")  
    return 5396;
  else if( aa == "FIW")  
    return 5397;
  else if( aa == "FIY")  
    return 5398;
  else if( aa == "FIV")  
    return 5399;
  else if( aa == "FLA")  
    return 5400;
  else if( aa == "FLR")  
    return 5401;
  else if( aa == "FLN")  
    return 5402;
  else if( aa == "FLD")  
    return 5403;
  else if( aa == "FLC")  
    return 5404;
  else if( aa == "FLQ")  
    return 5405;
  else if( aa == "FLE")  
    return 5406;
  else if( aa == "FLG")  
    return 5407;
  else if( aa == "FLH")  
    return 5408;
  else if( aa == "FLI")  
    return 5409;
  else if( aa == "FLL")  
    return 5410;
  else if( aa == "FLK")  
    return 5411;
  else if( aa == "FLM")  
    return 5412;
  else if( aa == "FLF")  
    return 5413;
  else if( aa == "FLP")  
    return 5414;
  else if( aa == "FLS")  
    return 5415;
  else if( aa == "FLT")  
    return 5416;
  else if( aa == "FLW")  
    return 5417;
  else if( aa == "FLY")  
    return 5418;
  else if( aa == "FLV")  
    return 5419;
  else if( aa == "FKA")  
    return 5420;
  else if( aa == "FKR")  
    return 5421;
  else if( aa == "FKN")  
    return 5422;
  else if( aa == "FKD")  
    return 5423;
  else if( aa == "FKC")  
    return 5424;
  else if( aa == "FKQ")  
    return 5425;
  else if( aa == "FKE")  
    return 5426;
  else if( aa == "FKG")  
    return 5427;
  else if( aa == "FKH")  
    return 5428;
  else if( aa == "FKI")  
    return 5429;
  else if( aa == "FKL")  
    return 5430;
  else if( aa == "FKK")  
    return 5431;
  else if( aa == "FKM")  
    return 5432;
  else if( aa == "FKF")  
    return 5433;
  else if( aa == "FKP")  
    return 5434;
  else if( aa == "FKS")  
    return 5435;
  else if( aa == "FKT")  
    return 5436;
  else if( aa == "FKW")  
    return 5437;
  else if( aa == "FKY")  
    return 5438;
  else if( aa == "FKV")  
    return 5439;
  else if( aa == "FMA")  
    return 5440;
  else if( aa == "FMR")  
    return 5441;
  else if( aa == "FMN")  
    return 5442;
  else if( aa == "FMD")  
    return 5443;
  else if( aa == "FMC")  
    return 5444;
  else if( aa == "FMQ")  
    return 5445;
  else if( aa == "FME")  
    return 5446;
  else if( aa == "FMG")  
    return 5447;
  else if( aa == "FMH")  
    return 5448;
  else if( aa == "FMI")  
    return 5449;
  else if( aa == "FML")  
    return 5450;
  else if( aa == "FMK")  
    return 5451;
  else if( aa == "FMM")  
    return 5452;
  else if( aa == "FMF")  
    return 5453;
  else if( aa == "FMP")  
    return 5454;
  else if( aa == "FMS")  
    return 5455;
  else if( aa == "FMT")  
    return 5456;
  else if( aa == "FMW")  
    return 5457;
  else if( aa == "FMY")  
    return 5458;
  else if( aa == "FMV")  
    return 5459;
  else if( aa == "FFA")  
    return 5460;
  else if( aa == "FFR")  
    return 5461;
  else if( aa == "FFN")  
    return 5462;
  else if( aa == "FFD")  
    return 5463;
  else if( aa == "FFC")  
    return 5464;
  else if( aa == "FFQ")  
    return 5465;
  else if( aa == "FFE")  
    return 5466;
  else if( aa == "FFG")  
    return 5467;
  else if( aa == "FFH")  
    return 5468;
  else if( aa == "FFI")  
    return 5469;
  else if( aa == "FFL")  
    return 5470;
  else if( aa == "FFK")  
    return 5471;
  else if( aa == "FFM")  
    return 5472;
  else if( aa == "FFF")  
    return 5473;
  else if( aa == "FFP")  
    return 5474;
  else if( aa == "FFS")  
    return 5475;
  else if( aa == "FFT")  
    return 5476;
  else if( aa == "FFW")  
    return 5477;
  else if( aa == "FFY")  
    return 5478;
  else if( aa == "FFV")  
    return 5479;
  else if( aa == "FPA")  
    return 5480;
  else if( aa == "FPR")  
    return 5481;
  else if( aa == "FPN")  
    return 5482;
  else if( aa == "FPD")  
    return 5483;
  else if( aa == "FPC")  
    return 5484;
  else if( aa == "FPQ")  
    return 5485;
  else if( aa == "FPE")  
    return 5486;
  else if( aa == "FPG")  
    return 5487;
  else if( aa == "FPH")  
    return 5488;
  else if( aa == "FPI")  
    return 5489;
  else if( aa == "FPL")  
    return 5490;
  else if( aa == "FPK")  
    return 5491;
  else if( aa == "FPM")  
    return 5492;
  else if( aa == "FPF")  
    return 5493;
  else if( aa == "FPP")  
    return 5494;
  else if( aa == "FPS")  
    return 5495;
  else if( aa == "FPT")  
    return 5496;
  else if( aa == "FPW")  
    return 5497;
  else if( aa == "FPY")  
    return 5498;
  else if( aa == "FPV")  
    return 5499;
  else if( aa == "FSA")  
    return 5500;
  else if( aa == "FSR")  
    return 5501;
  else if( aa == "FSN")  
    return 5502;
  else if( aa == "FSD")  
    return 5503;
  else if( aa == "FSC")  
    return 5504;
  else if( aa == "FSQ")  
    return 5505;
  else if( aa == "FSE")  
    return 5506;
  else if( aa == "FSG")  
    return 5507;
  else if( aa == "FSH")  
    return 5508;
  else if( aa == "FSI")  
    return 5509;
  else if( aa == "FSL")  
    return 5510;
  else if( aa == "FSK")  
    return 5511;
  else if( aa == "FSM")  
    return 5512;
  else if( aa == "FSF")  
    return 5513;
  else if( aa == "FSP")  
    return 5514;
  else if( aa == "FSS")  
    return 5515;
  else if( aa == "FST")  
    return 5516;
  else if( aa == "FSW")  
    return 5517;
  else if( aa == "FSY")  
    return 5518;
  else if( aa == "FSV")  
    return 5519;
  else if( aa == "FTA")  
    return 5520;
  else if( aa == "FTR")  
    return 5521;
  else if( aa == "FTN")  
    return 5522;
  else if( aa == "FTD")  
    return 5523;
  else if( aa == "FTC")  
    return 5524;
  else if( aa == "FTQ")  
    return 5525;
  else if( aa == "FTE")  
    return 5526;
  else if( aa == "FTG")  
    return 5527;
  else if( aa == "FTH")  
    return 5528;
  else if( aa == "FTI")  
    return 5529;
  else if( aa == "FTL")  
    return 5530;
  else if( aa == "FTK")  
    return 5531;
  else if( aa == "FTM")  
    return 5532;
  else if( aa == "FTF")  
    return 5533;
  else if( aa == "FTP")  
    return 5534;
  else if( aa == "FTS")  
    return 5535;
  else if( aa == "FTT")  
    return 5536;
  else if( aa == "FTW")  
    return 5537;
  else if( aa == "FTY")  
    return 5538;
  else if( aa == "FTV")  
    return 5539;
  else if( aa == "FWA")  
    return 5540;
  else if( aa == "FWR")  
    return 5541;
  else if( aa == "FWN")  
    return 5542;
  else if( aa == "FWD")  
    return 5543;
  else if( aa == "FWC")  
    return 5544;
  else if( aa == "FWQ")  
    return 5545;
  else if( aa == "FWE")  
    return 5546;
  else if( aa == "FWG")  
    return 5547;
  else if( aa == "FWH")  
    return 5548;
  else if( aa == "FWI")  
    return 5549;
  else if( aa == "FWL")  
    return 5550;
  else if( aa == "FWK")  
    return 5551;
  else if( aa == "FWM")  
    return 5552;
  else if( aa == "FWF")  
    return 5553;
  else if( aa == "FWP")  
    return 5554;
  else if( aa == "FWS")  
    return 5555;
  else if( aa == "FWT")  
    return 5556;
  else if( aa == "FWW")  
    return 5557;
  else if( aa == "FWY")  
    return 5558;
  else if( aa == "FWV")  
    return 5559;
  else if( aa == "FYA")  
    return 5560;
  else if( aa == "FYR")  
    return 5561;
  else if( aa == "FYN")  
    return 5562;
  else if( aa == "FYD")  
    return 5563;
  else if( aa == "FYC")  
    return 5564;
  else if( aa == "FYQ")  
    return 5565;
  else if( aa == "FYE")  
    return 5566;
  else if( aa == "FYG")  
    return 5567;
  else if( aa == "FYH")  
    return 5568;
  else if( aa == "FYI")  
    return 5569;
  else if( aa == "FYL")  
    return 5570;
  else if( aa == "FYK")  
    return 5571;
  else if( aa == "FYM")  
    return 5572;
  else if( aa == "FYF")  
    return 5573;
  else if( aa == "FYP")  
    return 5574;
  else if( aa == "FYS")  
    return 5575;
  else if( aa == "FYT")  
    return 5576;
  else if( aa == "FYW")  
    return 5577;
  else if( aa == "FYY")  
    return 5578;
  else if( aa == "FYV")  
    return 5579;
  else if( aa == "FVA")  
    return 5580;
  else if( aa == "FVR")  
    return 5581;
  else if( aa == "FVN")  
    return 5582;
  else if( aa == "FVD")  
    return 5583;
  else if( aa == "FVC")  
    return 5584;
  else if( aa == "FVQ")  
    return 5585;
  else if( aa == "FVE")  
    return 5586;
  else if( aa == "FVG")  
    return 5587;
  else if( aa == "FVH")  
    return 5588;
  else if( aa == "FVI")  
    return 5589;
  else if( aa == "FVL")  
    return 5590;
  else if( aa == "FVK")  
    return 5591;
  else if( aa == "FVM")  
    return 5592;
  else if( aa == "FVF")  
    return 5593;
  else if( aa == "FVP")  
    return 5594;
  else if( aa == "FVS")  
    return 5595;
  else if( aa == "FVT")  
    return 5596;
  else if( aa == "FVW")  
    return 5597;
  else if( aa == "FVY")  
    return 5598;
  else if( aa == "FVV")  
    return 5599;
  else if( aa == "PAA")  
    return 5600;
  else if( aa == "PAR")  
    return 5601;
  else if( aa == "PAN")  
    return 5602;
  else if( aa == "PAD")  
    return 5603;
  else if( aa == "PAC")  
    return 5604;
  else if( aa == "PAQ")  
    return 5605;
  else if( aa == "PAE")  
    return 5606;
  else if( aa == "PAG")  
    return 5607;
  else if( aa == "PAH")  
    return 5608;
  else if( aa == "PAI")  
    return 5609;
  else if( aa == "PAL")  
    return 5610;
  else if( aa == "PAK")  
    return 5611;
  else if( aa == "PAM")  
    return 5612;
  else if( aa == "PAF")  
    return 5613;
  else if( aa == "PAP")  
    return 5614;
  else if( aa == "PAS")  
    return 5615;
  else if( aa == "PAT")  
    return 5616;
  else if( aa == "PAW")  
    return 5617;
  else if( aa == "PAY")  
    return 5618;
  else if( aa == "PAV")  
    return 5619;
  else if( aa == "PRA")  
    return 5620;
  else if( aa == "PRR")  
    return 5621;
  else if( aa == "PRN")  
    return 5622;
  else if( aa == "PRD")  
    return 5623;
  else if( aa == "PRC")  
    return 5624;
  else if( aa == "PRQ")  
    return 5625;
  else if( aa == "PRE")  
    return 5626;
  else if( aa == "PRG")  
    return 5627;
  else if( aa == "PRH")  
    return 5628;
  else if( aa == "PRI")  
    return 5629;
  else if( aa == "PRL")  
    return 5630;
  else if( aa == "PRK")  
    return 5631;
  else if( aa == "PRM")  
    return 5632;
  else if( aa == "PRF")  
    return 5633;
  else if( aa == "PRP")  
    return 5634;
  else if( aa == "PRS")  
    return 5635;
  else if( aa == "PRT")  
    return 5636;
  else if( aa == "PRW")  
    return 5637;
  else if( aa == "PRY")  
    return 5638;
  else if( aa == "PRV")  
    return 5639;
  else if( aa == "PNA")  
    return 5640;
  else if( aa == "PNR")  
    return 5641;
  else if( aa == "PNN")  
    return 5642;
  else if( aa == "PND")  
    return 5643;
  else if( aa == "PNC")  
    return 5644;
  else if( aa == "PNQ")  
    return 5645;
  else if( aa == "PNE")  
    return 5646;
  else if( aa == "PNG")  
    return 5647;
  else if( aa == "PNH")  
    return 5648;
  else if( aa == "PNI")  
    return 5649;
  else if( aa == "PNL")  
    return 5650;
  else if( aa == "PNK")  
    return 5651;
  else if( aa == "PNM")  
    return 5652;
  else if( aa == "PNF")  
    return 5653;
  else if( aa == "PNP")  
    return 5654;
  else if( aa == "PNS")  
    return 5655;
  else if( aa == "PNT")  
    return 5656;
  else if( aa == "PNW")  
    return 5657;
  else if( aa == "PNY")  
    return 5658;
  else if( aa == "PNV")  
    return 5659;
  else if( aa == "PDA")  
    return 5660;
  else if( aa == "PDR")  
    return 5661;
  else if( aa == "PDN")  
    return 5662;
  else if( aa == "PDD")  
    return 5663;
  else if( aa == "PDC")  
    return 5664;
  else if( aa == "PDQ")  
    return 5665;
  else if( aa == "PDE")  
    return 5666;
  else if( aa == "PDG")  
    return 5667;
  else if( aa == "PDH")  
    return 5668;
  else if( aa == "PDI")  
    return 5669;
  else if( aa == "PDL")  
    return 5670;
  else if( aa == "PDK")  
    return 5671;
  else if( aa == "PDM")  
    return 5672;
  else if( aa == "PDF")  
    return 5673;
  else if( aa == "PDP")  
    return 5674;
  else if( aa == "PDS")  
    return 5675;
  else if( aa == "PDT")  
    return 5676;
  else if( aa == "PDW")  
    return 5677;
  else if( aa == "PDY")  
    return 5678;
  else if( aa == "PDV")  
    return 5679;
  else if( aa == "PCA")  
    return 5680;
  else if( aa == "PCR")  
    return 5681;
  else if( aa == "PCN")  
    return 5682;
  else if( aa == "PCD")  
    return 5683;
  else if( aa == "PCC")  
    return 5684;
  else if( aa == "PCQ")  
    return 5685;
  else if( aa == "PCE")  
    return 5686;
  else if( aa == "PCG")  
    return 5687;
  else if( aa == "PCH")  
    return 5688;
  else if( aa == "PCI")  
    return 5689;
  else if( aa == "PCL")  
    return 5690;
  else if( aa == "PCK")  
    return 5691;
  else if( aa == "PCM")  
    return 5692;
  else if( aa == "PCF")  
    return 5693;
  else if( aa == "PCP")  
    return 5694;
  else if( aa == "PCS")  
    return 5695;
  else if( aa == "PCT")  
    return 5696;
  else if( aa == "PCW")  
    return 5697;
  else if( aa == "PCY")  
    return 5698;
  else if( aa == "PCV")  
    return 5699;
  else if( aa == "PQA")  
    return 5700;
  else if( aa == "PQR")  
    return 5701;
  else if( aa == "PQN")  
    return 5702;
  else if( aa == "PQD")  
    return 5703;
  else if( aa == "PQC")  
    return 5704;
  else if( aa == "PQQ")  
    return 5705;
  else if( aa == "PQE")  
    return 5706;
  else if( aa == "PQG")  
    return 5707;
  else if( aa == "PQH")  
    return 5708;
  else if( aa == "PQI")  
    return 5709;
  else if( aa == "PQL")  
    return 5710;
  else if( aa == "PQK")  
    return 5711;
  else if( aa == "PQM")  
    return 5712;
  else if( aa == "PQF")  
    return 5713;
  else if( aa == "PQP")  
    return 5714;
  else if( aa == "PQS")  
    return 5715;
  else if( aa == "PQT")  
    return 5716;
  else if( aa == "PQW")  
    return 5717;
  else if( aa == "PQY")  
    return 5718;
  else if( aa == "PQV")  
    return 5719;
  else if( aa == "PEA")  
    return 5720;
  else if( aa == "PER")  
    return 5721;
  else if( aa == "PEN")  
    return 5722;
  else if( aa == "PED")  
    return 5723;
  else if( aa == "PEC")  
    return 5724;
  else if( aa == "PEQ")  
    return 5725;
  else if( aa == "PEE")  
    return 5726;
  else if( aa == "PEG")  
    return 5727;
  else if( aa == "PEH")  
    return 5728;
  else if( aa == "PEI")  
    return 5729;
  else if( aa == "PEL")  
    return 5730;
  else if( aa == "PEK")  
    return 5731;
  else if( aa == "PEM")  
    return 5732;
  else if( aa == "PEF")  
    return 5733;
  else if( aa == "PEP")  
    return 5734;
  else if( aa == "PES")  
    return 5735;
  else if( aa == "PET")  
    return 5736;
  else if( aa == "PEW")  
    return 5737;
  else if( aa == "PEY")  
    return 5738;
  else if( aa == "PEV")  
    return 5739;
  else if( aa == "PGA")  
    return 5740;
  else if( aa == "PGR")  
    return 5741;
  else if( aa == "PGN")  
    return 5742;
  else if( aa == "PGD")  
    return 5743;
  else if( aa == "PGC")  
    return 5744;
  else if( aa == "PGQ")  
    return 5745;
  else if( aa == "PGE")  
    return 5746;
  else if( aa == "PGG")  
    return 5747;
  else if( aa == "PGH")  
    return 5748;
  else if( aa == "PGI")  
    return 5749;
  else if( aa == "PGL")  
    return 5750;
  else if( aa == "PGK")  
    return 5751;
  else if( aa == "PGM")  
    return 5752;
  else if( aa == "PGF")  
    return 5753;
  else if( aa == "PGP")  
    return 5754;
  else if( aa == "PGS")  
    return 5755;
  else if( aa == "PGT")  
    return 5756;
  else if( aa == "PGW")  
    return 5757;
  else if( aa == "PGY")  
    return 5758;
  else if( aa == "PGV")  
    return 5759;
  else if( aa == "PHA")  
    return 5760;
  else if( aa == "PHR")  
    return 5761;
  else if( aa == "PHN")  
    return 5762;
  else if( aa == "PHD")  
    return 5763;
  else if( aa == "PHC")  
    return 5764;
  else if( aa == "PHQ")  
    return 5765;
  else if( aa == "PHE")  
    return 5766;
  else if( aa == "PHG")  
    return 5767;
  else if( aa == "PHH")  
    return 5768;
  else if( aa == "PHI")  
    return 5769;
  else if( aa == "PHL")  
    return 5770;
  else if( aa == "PHK")  
    return 5771;
  else if( aa == "PHM")  
    return 5772;
  else if( aa == "PHF")  
    return 5773;
  else if( aa == "PHP")  
    return 5774;
  else if( aa == "PHS")  
    return 5775;
  else if( aa == "PHT")  
    return 5776;
  else if( aa == "PHW")  
    return 5777;
  else if( aa == "PHY")  
    return 5778;
  else if( aa == "PHV")  
    return 5779;
  else if( aa == "PIA")  
    return 5780;
  else if( aa == "PIR")  
    return 5781;
  else if( aa == "PIN")  
    return 5782;
  else if( aa == "PID")  
    return 5783;
  else if( aa == "PIC")  
    return 5784;
  else if( aa == "PIQ")  
    return 5785;
  else if( aa == "PIE")  
    return 5786;
  else if( aa == "PIG")  
    return 5787;
  else if( aa == "PIH")  
    return 5788;
  else if( aa == "PII")  
    return 5789;
  else if( aa == "PIL")  
    return 5790;
  else if( aa == "PIK")  
    return 5791;
  else if( aa == "PIM")  
    return 5792;
  else if( aa == "PIF")  
    return 5793;
  else if( aa == "PIP")  
    return 5794;
  else if( aa == "PIS")  
    return 5795;
  else if( aa == "PIT")  
    return 5796;
  else if( aa == "PIW")  
    return 5797;
  else if( aa == "PIY")  
    return 5798;
  else if( aa == "PIV")  
    return 5799;
  else if( aa == "PLA")  
    return 5800;
  else if( aa == "PLR")  
    return 5801;
  else if( aa == "PLN")  
    return 5802;
  else if( aa == "PLD")  
    return 5803;
  else if( aa == "PLC")  
    return 5804;
  else if( aa == "PLQ")  
    return 5805;
  else if( aa == "PLE")  
    return 5806;
  else if( aa == "PLG")  
    return 5807;
  else if( aa == "PLH")  
    return 5808;
  else if( aa == "PLI")  
    return 5809;
  else if( aa == "PLL")  
    return 5810;
  else if( aa == "PLK")  
    return 5811;
  else if( aa == "PLM")  
    return 5812;
  else if( aa == "PLF")  
    return 5813;
  else if( aa == "PLP")  
    return 5814;
  else if( aa == "PLS")  
    return 5815;
  else if( aa == "PLT")  
    return 5816;
  else if( aa == "PLW")  
    return 5817;
  else if( aa == "PLY")  
    return 5818;
  else if( aa == "PLV")  
    return 5819;
  else if( aa == "PKA")  
    return 5820;
  else if( aa == "PKR")  
    return 5821;
  else if( aa == "PKN")  
    return 5822;
  else if( aa == "PKD")  
    return 5823;
  else if( aa == "PKC")  
    return 5824;
  else if( aa == "PKQ")  
    return 5825;
  else if( aa == "PKE")  
    return 5826;
  else if( aa == "PKG")  
    return 5827;
  else if( aa == "PKH")  
    return 5828;
  else if( aa == "PKI")  
    return 5829;
  else if( aa == "PKL")  
    return 5830;
  else if( aa == "PKK")  
    return 5831;
  else if( aa == "PKM")  
    return 5832;
  else if( aa == "PKF")  
    return 5833;
  else if( aa == "PKP")  
    return 5834;
  else if( aa == "PKS")  
    return 5835;
  else if( aa == "PKT")  
    return 5836;
  else if( aa == "PKW")  
    return 5837;
  else if( aa == "PKY")  
    return 5838;
  else if( aa == "PKV")  
    return 5839;
  else if( aa == "PMA")  
    return 5840;
  else if( aa == "PMR")  
    return 5841;
  else if( aa == "PMN")  
    return 5842;
  else if( aa == "PMD")  
    return 5843;
  else if( aa == "PMC")  
    return 5844;
  else if( aa == "PMQ")  
    return 5845;
  else if( aa == "PME")  
    return 5846;
  else if( aa == "PMG")  
    return 5847;
  else if( aa == "PMH")  
    return 5848;
  else if( aa == "PMI")  
    return 5849;
  else if( aa == "PML")  
    return 5850;
  else if( aa == "PMK")  
    return 5851;
  else if( aa == "PMM")  
    return 5852;
  else if( aa == "PMF")  
    return 5853;
  else if( aa == "PMP")  
    return 5854;
  else if( aa == "PMS")  
    return 5855;
  else if( aa == "PMT")  
    return 5856;
  else if( aa == "PMW")  
    return 5857;
  else if( aa == "PMY")  
    return 5858;
  else if( aa == "PMV")  
    return 5859;
  else if( aa == "PFA")  
    return 5860;
  else if( aa == "PFR")  
    return 5861;
  else if( aa == "PFN")  
    return 5862;
  else if( aa == "PFD")  
    return 5863;
  else if( aa == "PFC")  
    return 5864;
  else if( aa == "PFQ")  
    return 5865;
  else if( aa == "PFE")  
    return 5866;
  else if( aa == "PFG")  
    return 5867;
  else if( aa == "PFH")  
    return 5868;
  else if( aa == "PFI")  
    return 5869;
  else if( aa == "PFL")  
    return 5870;
  else if( aa == "PFK")  
    return 5871;
  else if( aa == "PFM")  
    return 5872;
  else if( aa == "PFF")  
    return 5873;
  else if( aa == "PFP")  
    return 5874;
  else if( aa == "PFS")  
    return 5875;
  else if( aa == "PFT")  
    return 5876;
  else if( aa == "PFW")  
    return 5877;
  else if( aa == "PFY")  
    return 5878;
  else if( aa == "PFV")  
    return 5879;
  else if( aa == "PPA")  
    return 5880;
  else if( aa == "PPR")  
    return 5881;
  else if( aa == "PPN")  
    return 5882;
  else if( aa == "PPD")  
    return 5883;
  else if( aa == "PPC")  
    return 5884;
  else if( aa == "PPQ")  
    return 5885;
  else if( aa == "PPE")  
    return 5886;
  else if( aa == "PPG")  
    return 5887;
  else if( aa == "PPH")  
    return 5888;
  else if( aa == "PPI")  
    return 5889;
  else if( aa == "PPL")  
    return 5890;
  else if( aa == "PPK")  
    return 5891;
  else if( aa == "PPM")  
    return 5892;
  else if( aa == "PPF")  
    return 5893;
  else if( aa == "PPP")  
    return 5894;
  else if( aa == "PPS")  
    return 5895;
  else if( aa == "PPT")  
    return 5896;
  else if( aa == "PPW")  
    return 5897;
  else if( aa == "PPY")  
    return 5898;
  else if( aa == "PPV")  
    return 5899;
  else if( aa == "PSA")  
    return 5900;
  else if( aa == "PSR")  
    return 5901;
  else if( aa == "PSN")  
    return 5902;
  else if( aa == "PSD")  
    return 5903;
  else if( aa == "PSC")  
    return 5904;
  else if( aa == "PSQ")  
    return 5905;
  else if( aa == "PSE")  
    return 5906;
  else if( aa == "PSG")  
    return 5907;
  else if( aa == "PSH")  
    return 5908;
  else if( aa == "PSI")  
    return 5909;
  else if( aa == "PSL")  
    return 5910;
  else if( aa == "PSK")  
    return 5911;
  else if( aa == "PSM")  
    return 5912;
  else if( aa == "PSF")  
    return 5913;
  else if( aa == "PSP")  
    return 5914;
  else if( aa == "PSS")  
    return 5915;
  else if( aa == "PST")  
    return 5916;
  else if( aa == "PSW")  
    return 5917;
  else if( aa == "PSY")  
    return 5918;
  else if( aa == "PSV")  
    return 5919;
  else if( aa == "PTA")  
    return 5920;
  else if( aa == "PTR")  
    return 5921;
  else if( aa == "PTN")  
    return 5922;
  else if( aa == "PTD")  
    return 5923;
  else if( aa == "PTC")  
    return 5924;
  else if( aa == "PTQ")  
    return 5925;
  else if( aa == "PTE")  
    return 5926;
  else if( aa == "PTG")  
    return 5927;
  else if( aa == "PTH")  
    return 5928;
  else if( aa == "PTI")  
    return 5929;
  else if( aa == "PTL")  
    return 5930;
  else if( aa == "PTK")  
    return 5931;
  else if( aa == "PTM")  
    return 5932;
  else if( aa == "PTF")  
    return 5933;
  else if( aa == "PTP")  
    return 5934;
  else if( aa == "PTS")  
    return 5935;
  else if( aa == "PTT")  
    return 5936;
  else if( aa == "PTW")  
    return 5937;
  else if( aa == "PTY")  
    return 5938;
  else if( aa == "PTV")  
    return 5939;
  else if( aa == "PWA")  
    return 5940;
  else if( aa == "PWR")  
    return 5941;
  else if( aa == "PWN")  
    return 5942;
  else if( aa == "PWD")  
    return 5943;
  else if( aa == "PWC")  
    return 5944;
  else if( aa == "PWQ")  
    return 5945;
  else if( aa == "PWE")  
    return 5946;
  else if( aa == "PWG")  
    return 5947;
  else if( aa == "PWH")  
    return 5948;
  else if( aa == "PWI")  
    return 5949;
  else if( aa == "PWL")  
    return 5950;
  else if( aa == "PWK")  
    return 5951;
  else if( aa == "PWM")  
    return 5952;
  else if( aa == "PWF")  
    return 5953;
  else if( aa == "PWP")  
    return 5954;
  else if( aa == "PWS")  
    return 5955;
  else if( aa == "PWT")  
    return 5956;
  else if( aa == "PWW")  
    return 5957;
  else if( aa == "PWY")  
    return 5958;
  else if( aa == "PWV")  
    return 5959;
  else if( aa == "PYA")  
    return 5960;
  else if( aa == "PYR")  
    return 5961;
  else if( aa == "PYN")  
    return 5962;
  else if( aa == "PYD")  
    return 5963;
  else if( aa == "PYC")  
    return 5964;
  else if( aa == "PYQ")  
    return 5965;
  else if( aa == "PYE")  
    return 5966;
  else if( aa == "PYG")  
    return 5967;
  else if( aa == "PYH")  
    return 5968;
  else if( aa == "PYI")  
    return 5969;
  else if( aa == "PYL")  
    return 5970;
  else if( aa == "PYK")  
    return 5971;
  else if( aa == "PYM")  
    return 5972;
  else if( aa == "PYF")  
    return 5973;
  else if( aa == "PYP")  
    return 5974;
  else if( aa == "PYS")  
    return 5975;
  else if( aa == "PYT")  
    return 5976;
  else if( aa == "PYW")  
    return 5977;
  else if( aa == "PYY")  
    return 5978;
  else if( aa == "PYV")  
    return 5979;
  else if( aa == "PVA")  
    return 5980;
  else if( aa == "PVR")  
    return 5981;
  else if( aa == "PVN")  
    return 5982;
  else if( aa == "PVD")  
    return 5983;
  else if( aa == "PVC")  
    return 5984;
  else if( aa == "PVQ")  
    return 5985;
  else if( aa == "PVE")  
    return 5986;
  else if( aa == "PVG")  
    return 5987;
  else if( aa == "PVH")  
    return 5988;
  else if( aa == "PVI")  
    return 5989;
  else if( aa == "PVL")  
    return 5990;
  else if( aa == "PVK")  
    return 5991;
  else if( aa == "PVM")  
    return 5992;
  else if( aa == "PVF")  
    return 5993;
  else if( aa == "PVP")  
    return 5994;
  else if( aa == "PVS")  
    return 5995;
  else if( aa == "PVT")  
    return 5996;
  else if( aa == "PVW")  
    return 5997;
  else if( aa == "PVY")  
    return 5998;
  else if( aa == "PVV")  
    return 5999;
  else if( aa == "SAA")  
    return 6000;
  else if( aa == "SAR")  
    return 6001;
  else if( aa == "SAN")  
    return 6002;
  else if( aa == "SAD")  
    return 6003;
  else if( aa == "SAC")  
    return 6004;
  else if( aa == "SAQ")  
    return 6005;
  else if( aa == "SAE")  
    return 6006;
  else if( aa == "SAG")  
    return 6007;
  else if( aa == "SAH")  
    return 6008;
  else if( aa == "SAI")  
    return 6009;
  else if( aa == "SAL")  
    return 6010;
  else if( aa == "SAK")  
    return 6011;
  else if( aa == "SAM")  
    return 6012;
  else if( aa == "SAF")  
    return 6013;
  else if( aa == "SAP")  
    return 6014;
  else if( aa == "SAS")  
    return 6015;
  else if( aa == "SAT")  
    return 6016;
  else if( aa == "SAW")  
    return 6017;
  else if( aa == "SAY")  
    return 6018;
  else if( aa == "SAV")  
    return 6019;
  else if( aa == "SRA")  
    return 6020;
  else if( aa == "SRR")  
    return 6021;
  else if( aa == "SRN")  
    return 6022;
  else if( aa == "SRD")  
    return 6023;
  else if( aa == "SRC")  
    return 6024;
  else if( aa == "SRQ")  
    return 6025;
  else if( aa == "SRE")  
    return 6026;
  else if( aa == "SRG")  
    return 6027;
  else if( aa == "SRH")  
    return 6028;
  else if( aa == "SRI")  
    return 6029;
  else if( aa == "SRL")  
    return 6030;
  else if( aa == "SRK")  
    return 6031;
  else if( aa == "SRM")  
    return 6032;
  else if( aa == "SRF")  
    return 6033;
  else if( aa == "SRP")  
    return 6034;
  else if( aa == "SRS")  
    return 6035;
  else if( aa == "SRT")  
    return 6036;
  else if( aa == "SRW")  
    return 6037;
  else if( aa == "SRY")  
    return 6038;
  else if( aa == "SRV")  
    return 6039;
  else if( aa == "SNA")  
    return 6040;
  else if( aa == "SNR")  
    return 6041;
  else if( aa == "SNN")  
    return 6042;
  else if( aa == "SND")  
    return 6043;
  else if( aa == "SNC")  
    return 6044;
  else if( aa == "SNQ")  
    return 6045;
  else if( aa == "SNE")  
    return 6046;
  else if( aa == "SNG")  
    return 6047;
  else if( aa == "SNH")  
    return 6048;
  else if( aa == "SNI")  
    return 6049;
  else if( aa == "SNL")  
    return 6050;
  else if( aa == "SNK")  
    return 6051;
  else if( aa == "SNM")  
    return 6052;
  else if( aa == "SNF")  
    return 6053;
  else if( aa == "SNP")  
    return 6054;
  else if( aa == "SNS")  
    return 6055;
  else if( aa == "SNT")  
    return 6056;
  else if( aa == "SNW")  
    return 6057;
  else if( aa == "SNY")  
    return 6058;
  else if( aa == "SNV")  
    return 6059;
  else if( aa == "SDA")  
    return 6060;
  else if( aa == "SDR")  
    return 6061;
  else if( aa == "SDN")  
    return 6062;
  else if( aa == "SDD")  
    return 6063;
  else if( aa == "SDC")  
    return 6064;
  else if( aa == "SDQ")  
    return 6065;
  else if( aa == "SDE")  
    return 6066;
  else if( aa == "SDG")  
    return 6067;
  else if( aa == "SDH")  
    return 6068;
  else if( aa == "SDI")  
    return 6069;
  else if( aa == "SDL")  
    return 6070;
  else if( aa == "SDK")  
    return 6071;
  else if( aa == "SDM")  
    return 6072;
  else if( aa == "SDF")  
    return 6073;
  else if( aa == "SDP")  
    return 6074;
  else if( aa == "SDS")  
    return 6075;
  else if( aa == "SDT")  
    return 6076;
  else if( aa == "SDW")  
    return 6077;
  else if( aa == "SDY")  
    return 6078;
  else if( aa == "SDV")  
    return 6079;
  else if( aa == "SCA")  
    return 6080;
  else if( aa == "SCR")  
    return 6081;
  else if( aa == "SCN")  
    return 6082;
  else if( aa == "SCD")  
    return 6083;
  else if( aa == "SCC")  
    return 6084;
  else if( aa == "SCQ")  
    return 6085;
  else if( aa == "SCE")  
    return 6086;
  else if( aa == "SCG")  
    return 6087;
  else if( aa == "SCH")  
    return 6088;
  else if( aa == "SCI")  
    return 6089;
  else if( aa == "SCL")  
    return 6090;
  else if( aa == "SCK")  
    return 6091;
  else if( aa == "SCM")  
    return 6092;
  else if( aa == "SCF")  
    return 6093;
  else if( aa == "SCP")  
    return 6094;
  else if( aa == "SCS")  
    return 6095;
  else if( aa == "SCT")  
    return 6096;
  else if( aa == "SCW")  
    return 6097;
  else if( aa == "SCY")  
    return 6098;
  else if( aa == "SCV")  
    return 6099;
  else if( aa == "SQA")  
    return 6100;
  else if( aa == "SQR")  
    return 6101;
  else if( aa == "SQN")  
    return 6102;
  else if( aa == "SQD")  
    return 6103;
  else if( aa == "SQC")  
    return 6104;
  else if( aa == "SQQ")  
    return 6105;
  else if( aa == "SQE")  
    return 6106;
  else if( aa == "SQG")  
    return 6107;
  else if( aa == "SQH")  
    return 6108;
  else if( aa == "SQI")  
    return 6109;
  else if( aa == "SQL")  
    return 6110;
  else if( aa == "SQK")  
    return 6111;
  else if( aa == "SQM")  
    return 6112;
  else if( aa == "SQF")  
    return 6113;
  else if( aa == "SQP")  
    return 6114;
  else if( aa == "SQS")  
    return 6115;
  else if( aa == "SQT")  
    return 6116;
  else if( aa == "SQW")  
    return 6117;
  else if( aa == "SQY")  
    return 6118;
  else if( aa == "SQV")  
    return 6119;
  else if( aa == "SEA")  
    return 6120;
  else if( aa == "SER")  
    return 6121;
  else if( aa == "SEN")  
    return 6122;
  else if( aa == "SED")  
    return 6123;
  else if( aa == "SEC")  
    return 6124;
  else if( aa == "SEQ")  
    return 6125;
  else if( aa == "SEE")  
    return 6126;
  else if( aa == "SEG")  
    return 6127;
  else if( aa == "SEH")  
    return 6128;
  else if( aa == "SEI")  
    return 6129;
  else if( aa == "SEL")  
    return 6130;
  else if( aa == "SEK")  
    return 6131;
  else if( aa == "SEM")  
    return 6132;
  else if( aa == "SEF")  
    return 6133;
  else if( aa == "SEP")  
    return 6134;
  else if( aa == "SES")  
    return 6135;
  else if( aa == "SET")  
    return 6136;
  else if( aa == "SEW")  
    return 6137;
  else if( aa == "SEY")  
    return 6138;
  else if( aa == "SEV")  
    return 6139;
  else if( aa == "SGA")  
    return 6140;
  else if( aa == "SGR")  
    return 6141;
  else if( aa == "SGN")  
    return 6142;
  else if( aa == "SGD")  
    return 6143;
  else if( aa == "SGC")  
    return 6144;
  else if( aa == "SGQ")  
    return 6145;
  else if( aa == "SGE")  
    return 6146;
  else if( aa == "SGG")  
    return 6147;
  else if( aa == "SGH")  
    return 6148;
  else if( aa == "SGI")  
    return 6149;
  else if( aa == "SGL")  
    return 6150;
  else if( aa == "SGK")  
    return 6151;
  else if( aa == "SGM")  
    return 6152;
  else if( aa == "SGF")  
    return 6153;
  else if( aa == "SGP")  
    return 6154;
  else if( aa == "SGS")  
    return 6155;
  else if( aa == "SGT")  
    return 6156;
  else if( aa == "SGW")  
    return 6157;
  else if( aa == "SGY")  
    return 6158;
  else if( aa == "SGV")  
    return 6159;
  else if( aa == "SHA")  
    return 6160;
  else if( aa == "SHR")  
    return 6161;
  else if( aa == "SHN")  
    return 6162;
  else if( aa == "SHD")  
    return 6163;
  else if( aa == "SHC")  
    return 6164;
  else if( aa == "SHQ")  
    return 6165;
  else if( aa == "SHE")  
    return 6166;
  else if( aa == "SHG")  
    return 6167;
  else if( aa == "SHH")  
    return 6168;
  else if( aa == "SHI")  
    return 6169;
  else if( aa == "SHL")  
    return 6170;
  else if( aa == "SHK")  
    return 6171;
  else if( aa == "SHM")  
    return 6172;
  else if( aa == "SHF")  
    return 6173;
  else if( aa == "SHP")  
    return 6174;
  else if( aa == "SHS")  
    return 6175;
  else if( aa == "SHT")  
    return 6176;
  else if( aa == "SHW")  
    return 6177;
  else if( aa == "SHY")  
    return 6178;
  else if( aa == "SHV")  
    return 6179;
  else if( aa == "SIA")  
    return 6180;
  else if( aa == "SIR")  
    return 6181;
  else if( aa == "SIN")  
    return 6182;
  else if( aa == "SID")  
    return 6183;
  else if( aa == "SIC")  
    return 6184;
  else if( aa == "SIQ")  
    return 6185;
  else if( aa == "SIE")  
    return 6186;
  else if( aa == "SIG")  
    return 6187;
  else if( aa == "SIH")  
    return 6188;
  else if( aa == "SII")  
    return 6189;
  else if( aa == "SIL")  
    return 6190;
  else if( aa == "SIK")  
    return 6191;
  else if( aa == "SIM")  
    return 6192;
  else if( aa == "SIF")  
    return 6193;
  else if( aa == "SIP")  
    return 6194;
  else if( aa == "SIS")  
    return 6195;
  else if( aa == "SIT")  
    return 6196;
  else if( aa == "SIW")  
    return 6197;
  else if( aa == "SIY")  
    return 6198;
  else if( aa == "SIV")  
    return 6199;
  else if( aa == "SLA")  
    return 6200;
  else if( aa == "SLR")  
    return 6201;
  else if( aa == "SLN")  
    return 6202;
  else if( aa == "SLD")  
    return 6203;
  else if( aa == "SLC")  
    return 6204;
  else if( aa == "SLQ")  
    return 6205;
  else if( aa == "SLE")  
    return 6206;
  else if( aa == "SLG")  
    return 6207;
  else if( aa == "SLH")  
    return 6208;
  else if( aa == "SLI")  
    return 6209;
  else if( aa == "SLL")  
    return 6210;
  else if( aa == "SLK")  
    return 6211;
  else if( aa == "SLM")  
    return 6212;
  else if( aa == "SLF")  
    return 6213;
  else if( aa == "SLP")  
    return 6214;
  else if( aa == "SLS")  
    return 6215;
  else if( aa == "SLT")  
    return 6216;
  else if( aa == "SLW")  
    return 6217;
  else if( aa == "SLY")  
    return 6218;
  else if( aa == "SLV")  
    return 6219;
  else if( aa == "SKA")  
    return 6220;
  else if( aa == "SKR")  
    return 6221;
  else if( aa == "SKN")  
    return 6222;
  else if( aa == "SKD")  
    return 6223;
  else if( aa == "SKC")  
    return 6224;
  else if( aa == "SKQ")  
    return 6225;
  else if( aa == "SKE")  
    return 6226;
  else if( aa == "SKG")  
    return 6227;
  else if( aa == "SKH")  
    return 6228;
  else if( aa == "SKI")  
    return 6229;
  else if( aa == "SKL")  
    return 6230;
  else if( aa == "SKK")  
    return 6231;
  else if( aa == "SKM")  
    return 6232;
  else if( aa == "SKF")  
    return 6233;
  else if( aa == "SKP")  
    return 6234;
  else if( aa == "SKS")  
    return 6235;
  else if( aa == "SKT")  
    return 6236;
  else if( aa == "SKW")  
    return 6237;
  else if( aa == "SKY")  
    return 6238;
  else if( aa == "SKV")  
    return 6239;
  else if( aa == "SMA")  
    return 6240;
  else if( aa == "SMR")  
    return 6241;
  else if( aa == "SMN")  
    return 6242;
  else if( aa == "SMD")  
    return 6243;
  else if( aa == "SMC")  
    return 6244;
  else if( aa == "SMQ")  
    return 6245;
  else if( aa == "SME")  
    return 6246;
  else if( aa == "SMG")  
    return 6247;
  else if( aa == "SMH")  
    return 6248;
  else if( aa == "SMI")  
    return 6249;
  else if( aa == "SML")  
    return 6250;
  else if( aa == "SMK")  
    return 6251;
  else if( aa == "SMM")  
    return 6252;
  else if( aa == "SMF")  
    return 6253;
  else if( aa == "SMP")  
    return 6254;
  else if( aa == "SMS")  
    return 6255;
  else if( aa == "SMT")  
    return 6256;
  else if( aa == "SMW")  
    return 6257;
  else if( aa == "SMY")  
    return 6258;
  else if( aa == "SMV")  
    return 6259;
  else if( aa == "SFA")  
    return 6260;
  else if( aa == "SFR")  
    return 6261;
  else if( aa == "SFN")  
    return 6262;
  else if( aa == "SFD")  
    return 6263;
  else if( aa == "SFC")  
    return 6264;
  else if( aa == "SFQ")  
    return 6265;
  else if( aa == "SFE")  
    return 6266;
  else if( aa == "SFG")  
    return 6267;
  else if( aa == "SFH")  
    return 6268;
  else if( aa == "SFI")  
    return 6269;
  else if( aa == "SFL")  
    return 6270;
  else if( aa == "SFK")  
    return 6271;
  else if( aa == "SFM")  
    return 6272;
  else if( aa == "SFF")  
    return 6273;
  else if( aa == "SFP")  
    return 6274;
  else if( aa == "SFS")  
    return 6275;
  else if( aa == "SFT")  
    return 6276;
  else if( aa == "SFW")  
    return 6277;
  else if( aa == "SFY")  
    return 6278;
  else if( aa == "SFV")  
    return 6279;
  else if( aa == "SPA")  
    return 6280;
  else if( aa == "SPR")  
    return 6281;
  else if( aa == "SPN")  
    return 6282;
  else if( aa == "SPD")  
    return 6283;
  else if( aa == "SPC")  
    return 6284;
  else if( aa == "SPQ")  
    return 6285;
  else if( aa == "SPE")  
    return 6286;
  else if( aa == "SPG")  
    return 6287;
  else if( aa == "SPH")  
    return 6288;
  else if( aa == "SPI")  
    return 6289;
  else if( aa == "SPL")  
    return 6290;
  else if( aa == "SPK")  
    return 6291;
  else if( aa == "SPM")  
    return 6292;
  else if( aa == "SPF")  
    return 6293;
  else if( aa == "SPP")  
    return 6294;
  else if( aa == "SPS")  
    return 6295;
  else if( aa == "SPT")  
    return 6296;
  else if( aa == "SPW")  
    return 6297;
  else if( aa == "SPY")  
    return 6298;
  else if( aa == "SPV")  
    return 6299;
  else if( aa == "SSA")  
    return 6300;
  else if( aa == "SSR")  
    return 6301;
  else if( aa == "SSN")  
    return 6302;
  else if( aa == "SSD")  
    return 6303;
  else if( aa == "SSC")  
    return 6304;
  else if( aa == "SSQ")  
    return 6305;
  else if( aa == "SSE")  
    return 6306;
  else if( aa == "SSG")  
    return 6307;
  else if( aa == "SSH")  
    return 6308;
  else if( aa == "SSI")  
    return 6309;
  else if( aa == "SSL")  
    return 6310;
  else if( aa == "SSK")  
    return 6311;
  else if( aa == "SSM")  
    return 6312;
  else if( aa == "SSF")  
    return 6313;
  else if( aa == "SSP")  
    return 6314;
  else if( aa == "SSS")  
    return 6315;
  else if( aa == "SST")  
    return 6316;
  else if( aa == "SSW")  
    return 6317;
  else if( aa == "SSY")  
    return 6318;
  else if( aa == "SSV")  
    return 6319;
  else if( aa == "STA")  
    return 6320;
  else if( aa == "STR")  
    return 6321;
  else if( aa == "STN")  
    return 6322;
  else if( aa == "STD")  
    return 6323;
  else if( aa == "STC")  
    return 6324;
  else if( aa == "STQ")  
    return 6325;
  else if( aa == "STE")  
    return 6326;
  else if( aa == "STG")  
    return 6327;
  else if( aa == "STH")  
    return 6328;
  else if( aa == "STI")  
    return 6329;
  else if( aa == "STL")  
    return 6330;
  else if( aa == "STK")  
    return 6331;
  else if( aa == "STM")  
    return 6332;
  else if( aa == "STF")  
    return 6333;
  else if( aa == "STP")  
    return 6334;
  else if( aa == "STS")  
    return 6335;
  else if( aa == "STT")  
    return 6336;
  else if( aa == "STW")  
    return 6337;
  else if( aa == "STY")  
    return 6338;
  else if( aa == "STV")  
    return 6339;
  else if( aa == "SWA")  
    return 6340;
  else if( aa == "SWR")  
    return 6341;
  else if( aa == "SWN")  
    return 6342;
  else if( aa == "SWD")  
    return 6343;
  else if( aa == "SWC")  
    return 6344;
  else if( aa == "SWQ")  
    return 6345;
  else if( aa == "SWE")  
    return 6346;
  else if( aa == "SWG")  
    return 6347;
  else if( aa == "SWH")  
    return 6348;
  else if( aa == "SWI")  
    return 6349;
  else if( aa == "SWL")  
    return 6350;
  else if( aa == "SWK")  
    return 6351;
  else if( aa == "SWM")  
    return 6352;
  else if( aa == "SWF")  
    return 6353;
  else if( aa == "SWP")  
    return 6354;
  else if( aa == "SWS")  
    return 6355;
  else if( aa == "SWT")  
    return 6356;
  else if( aa == "SWW")  
    return 6357;
  else if( aa == "SWY")  
    return 6358;
  else if( aa == "SWV")  
    return 6359;
  else if( aa == "SYA")  
    return 6360;
  else if( aa == "SYR")  
    return 6361;
  else if( aa == "SYN")  
    return 6362;
  else if( aa == "SYD")  
    return 6363;
  else if( aa == "SYC")  
    return 6364;
  else if( aa == "SYQ")  
    return 6365;
  else if( aa == "SYE")  
    return 6366;
  else if( aa == "SYG")  
    return 6367;
  else if( aa == "SYH")  
    return 6368;
  else if( aa == "SYI")  
    return 6369;
  else if( aa == "SYL")  
    return 6370;
  else if( aa == "SYK")  
    return 6371;
  else if( aa == "SYM")  
    return 6372;
  else if( aa == "SYF")  
    return 6373;
  else if( aa == "SYP")  
    return 6374;
  else if( aa == "SYS")  
    return 6375;
  else if( aa == "SYT")  
    return 6376;
  else if( aa == "SYW")  
    return 6377;
  else if( aa == "SYY")  
    return 6378;
  else if( aa == "SYV")  
    return 6379;
  else if( aa == "SVA")  
    return 6380;
  else if( aa == "SVR")  
    return 6381;
  else if( aa == "SVN")  
    return 6382;
  else if( aa == "SVD")  
    return 6383;
  else if( aa == "SVC")  
    return 6384;
  else if( aa == "SVQ")  
    return 6385;
  else if( aa == "SVE")  
    return 6386;
  else if( aa == "SVG")  
    return 6387;
  else if( aa == "SVH")  
    return 6388;
  else if( aa == "SVI")  
    return 6389;
  else if( aa == "SVL")  
    return 6390;
  else if( aa == "SVK")  
    return 6391;
  else if( aa == "SVM")  
    return 6392;
  else if( aa == "SVF")  
    return 6393;
  else if( aa == "SVP")  
    return 6394;
  else if( aa == "SVS")  
    return 6395;
  else if( aa == "SVT")  
    return 6396;
  else if( aa == "SVW")  
    return 6397;
  else if( aa == "SVY")  
    return 6398;
  else if( aa == "SVV")  
    return 6399;
  else if( aa == "TAA")  
    return 6400;
  else if( aa == "TAR")  
    return 6401;
  else if( aa == "TAN")  
    return 6402;
  else if( aa == "TAD")  
    return 6403;
  else if( aa == "TAC")  
    return 6404;
  else if( aa == "TAQ")  
    return 6405;
  else if( aa == "TAE")  
    return 6406;
  else if( aa == "TAG")  
    return 6407;
  else if( aa == "TAH")  
    return 6408;
  else if( aa == "TAI")  
    return 6409;
  else if( aa == "TAL")  
    return 6410;
  else if( aa == "TAK")  
    return 6411;
  else if( aa == "TAM")  
    return 6412;
  else if( aa == "TAF")  
    return 6413;
  else if( aa == "TAP")  
    return 6414;
  else if( aa == "TAS")  
    return 6415;
  else if( aa == "TAT")  
    return 6416;
  else if( aa == "TAW")  
    return 6417;
  else if( aa == "TAY")  
    return 6418;
  else if( aa == "TAV")  
    return 6419;
  else if( aa == "TRA")  
    return 6420;
  else if( aa == "TRR")  
    return 6421;
  else if( aa == "TRN")  
    return 6422;
  else if( aa == "TRD")  
    return 6423;
  else if( aa == "TRC")  
    return 6424;
  else if( aa == "TRQ")  
    return 6425;
  else if( aa == "TRE")  
    return 6426;
  else if( aa == "TRG")  
    return 6427;
  else if( aa == "TRH")  
    return 6428;
  else if( aa == "TRI")  
    return 6429;
  else if( aa == "TRL")  
    return 6430;
  else if( aa == "TRK")  
    return 6431;
  else if( aa == "TRM")  
    return 6432;
  else if( aa == "TRF")  
    return 6433;
  else if( aa == "TRP")  
    return 6434;
  else if( aa == "TRS")  
    return 6435;
  else if( aa == "TRT")  
    return 6436;
  else if( aa == "TRW")  
    return 6437;
  else if( aa == "TRY")  
    return 6438;
  else if( aa == "TRV")  
    return 6439;
  else if( aa == "TNA")  
    return 6440;
  else if( aa == "TNR")  
    return 6441;
  else if( aa == "TNN")  
    return 6442;
  else if( aa == "TND")  
    return 6443;
  else if( aa == "TNC")  
    return 6444;
  else if( aa == "TNQ")  
    return 6445;
  else if( aa == "TNE")  
    return 6446;
  else if( aa == "TNG")  
    return 6447;
  else if( aa == "TNH")  
    return 6448;
  else if( aa == "TNI")  
    return 6449;
  else if( aa == "TNL")  
    return 6450;
  else if( aa == "TNK")  
    return 6451;
  else if( aa == "TNM")  
    return 6452;
  else if( aa == "TNF")  
    return 6453;
  else if( aa == "TNP")  
    return 6454;
  else if( aa == "TNS")  
    return 6455;
  else if( aa == "TNT")  
    return 6456;
  else if( aa == "TNW")  
    return 6457;
  else if( aa == "TNY")  
    return 6458;
  else if( aa == "TNV")  
    return 6459;
  else if( aa == "TDA")  
    return 6460;
  else if( aa == "TDR")  
    return 6461;
  else if( aa == "TDN")  
    return 6462;
  else if( aa == "TDD")  
    return 6463;
  else if( aa == "TDC")  
    return 6464;
  else if( aa == "TDQ")  
    return 6465;
  else if( aa == "TDE")  
    return 6466;
  else if( aa == "TDG")  
    return 6467;
  else if( aa == "TDH")  
    return 6468;
  else if( aa == "TDI")  
    return 6469;
  else if( aa == "TDL")  
    return 6470;
  else if( aa == "TDK")  
    return 6471;
  else if( aa == "TDM")  
    return 6472;
  else if( aa == "TDF")  
    return 6473;
  else if( aa == "TDP")  
    return 6474;
  else if( aa == "TDS")  
    return 6475;
  else if( aa == "TDT")  
    return 6476;
  else if( aa == "TDW")  
    return 6477;
  else if( aa == "TDY")  
    return 6478;
  else if( aa == "TDV")  
    return 6479;
  else if( aa == "TCA")  
    return 6480;
  else if( aa == "TCR")  
    return 6481;
  else if( aa == "TCN")  
    return 6482;
  else if( aa == "TCD")  
    return 6483;
  else if( aa == "TCC")  
    return 6484;
  else if( aa == "TCQ")  
    return 6485;
  else if( aa == "TCE")  
    return 6486;
  else if( aa == "TCG")  
    return 6487;
  else if( aa == "TCH")  
    return 6488;
  else if( aa == "TCI")  
    return 6489;
  else if( aa == "TCL")  
    return 6490;
  else if( aa == "TCK")  
    return 6491;
  else if( aa == "TCM")  
    return 6492;
  else if( aa == "TCF")  
    return 6493;
  else if( aa == "TCP")  
    return 6494;
  else if( aa == "TCS")  
    return 6495;
  else if( aa == "TCT")  
    return 6496;
  else if( aa == "TCW")  
    return 6497;
  else if( aa == "TCY")  
    return 6498;
  else if( aa == "TCV")  
    return 6499;
  else if( aa == "TQA")  
    return 6500;
  else if( aa == "TQR")  
    return 6501;
  else if( aa == "TQN")  
    return 6502;
  else if( aa == "TQD")  
    return 6503;
  else if( aa == "TQC")  
    return 6504;
  else if( aa == "TQQ")  
    return 6505;
  else if( aa == "TQE")  
    return 6506;
  else if( aa == "TQG")  
    return 6507;
  else if( aa == "TQH")  
    return 6508;
  else if( aa == "TQI")  
    return 6509;
  else if( aa == "TQL")  
    return 6510;
  else if( aa == "TQK")  
    return 6511;
  else if( aa == "TQM")  
    return 6512;
  else if( aa == "TQF")  
    return 6513;
  else if( aa == "TQP")  
    return 6514;
  else if( aa == "TQS")  
    return 6515;
  else if( aa == "TQT")  
    return 6516;
  else if( aa == "TQW")  
    return 6517;
  else if( aa == "TQY")  
    return 6518;
  else if( aa == "TQV")  
    return 6519;
  else if( aa == "TEA")  
    return 6520;
  else if( aa == "TER")  
    return 6521;
  else if( aa == "TEN")  
    return 6522;
  else if( aa == "TED")  
    return 6523;
  else if( aa == "TEC")  
    return 6524;
  else if( aa == "TEQ")  
    return 6525;
  else if( aa == "TEE")  
    return 6526;
  else if( aa == "TEG")  
    return 6527;
  else if( aa == "TEH")  
    return 6528;
  else if( aa == "TEI")  
    return 6529;
  else if( aa == "TEL")  
    return 6530;
  else if( aa == "TEK")  
    return 6531;
  else if( aa == "TEM")  
    return 6532;
  else if( aa == "TEF")  
    return 6533;
  else if( aa == "TEP")  
    return 6534;
  else if( aa == "TES")  
    return 6535;
  else if( aa == "TET")  
    return 6536;
  else if( aa == "TEW")  
    return 6537;
  else if( aa == "TEY")  
    return 6538;
  else if( aa == "TEV")  
    return 6539;
  else if( aa == "TGA")  
    return 6540;
  else if( aa == "TGR")  
    return 6541;
  else if( aa == "TGN")  
    return 6542;
  else if( aa == "TGD")  
    return 6543;
  else if( aa == "TGC")  
    return 6544;
  else if( aa == "TGQ")  
    return 6545;
  else if( aa == "TGE")  
    return 6546;
  else if( aa == "TGG")  
    return 6547;
  else if( aa == "TGH")  
    return 6548;
  else if( aa == "TGI")  
    return 6549;
  else if( aa == "TGL")  
    return 6550;
  else if( aa == "TGK")  
    return 6551;
  else if( aa == "TGM")  
    return 6552;
  else if( aa == "TGF")  
    return 6553;
  else if( aa == "TGP")  
    return 6554;
  else if( aa == "TGS")  
    return 6555;
  else if( aa == "TGT")  
    return 6556;
  else if( aa == "TGW")  
    return 6557;
  else if( aa == "TGY")  
    return 6558;
  else if( aa == "TGV")  
    return 6559;
  else if( aa == "THA")  
    return 6560;
  else if( aa == "THR")  
    return 6561;
  else if( aa == "THN")  
    return 6562;
  else if( aa == "THD")  
    return 6563;
  else if( aa == "THC")  
    return 6564;
  else if( aa == "THQ")  
    return 6565;
  else if( aa == "THE")  
    return 6566;
  else if( aa == "THG")  
    return 6567;
  else if( aa == "THH")  
    return 6568;
  else if( aa == "THI")  
    return 6569;
  else if( aa == "THL")  
    return 6570;
  else if( aa == "THK")  
    return 6571;
  else if( aa == "THM")  
    return 6572;
  else if( aa == "THF")  
    return 6573;
  else if( aa == "THP")  
    return 6574;
  else if( aa == "THS")  
    return 6575;
  else if( aa == "THT")  
    return 6576;
  else if( aa == "THW")  
    return 6577;
  else if( aa == "THY")  
    return 6578;
  else if( aa == "THV")  
    return 6579;
  else if( aa == "TIA")  
    return 6580;
  else if( aa == "TIR")  
    return 6581;
  else if( aa == "TIN")  
    return 6582;
  else if( aa == "TID")  
    return 6583;
  else if( aa == "TIC")  
    return 6584;
  else if( aa == "TIQ")  
    return 6585;
  else if( aa == "TIE")  
    return 6586;
  else if( aa == "TIG")  
    return 6587;
  else if( aa == "TIH")  
    return 6588;
  else if( aa == "TII")  
    return 6589;
  else if( aa == "TIL")  
    return 6590;
  else if( aa == "TIK")  
    return 6591;
  else if( aa == "TIM")  
    return 6592;
  else if( aa == "TIF")  
    return 6593;
  else if( aa == "TIP")  
    return 6594;
  else if( aa == "TIS")  
    return 6595;
  else if( aa == "TIT")  
    return 6596;
  else if( aa == "TIW")  
    return 6597;
  else if( aa == "TIY")  
    return 6598;
  else if( aa == "TIV")  
    return 6599;
  else if( aa == "TLA")  
    return 6600;
  else if( aa == "TLR")  
    return 6601;
  else if( aa == "TLN")  
    return 6602;
  else if( aa == "TLD")  
    return 6603;
  else if( aa == "TLC")  
    return 6604;
  else if( aa == "TLQ")  
    return 6605;
  else if( aa == "TLE")  
    return 6606;
  else if( aa == "TLG")  
    return 6607;
  else if( aa == "TLH")  
    return 6608;
  else if( aa == "TLI")  
    return 6609;
  else if( aa == "TLL")  
    return 6610;
  else if( aa == "TLK")  
    return 6611;
  else if( aa == "TLM")  
    return 6612;
  else if( aa == "TLF")  
    return 6613;
  else if( aa == "TLP")  
    return 6614;
  else if( aa == "TLS")  
    return 6615;
  else if( aa == "TLT")  
    return 6616;
  else if( aa == "TLW")  
    return 6617;
  else if( aa == "TLY")  
    return 6618;
  else if( aa == "TLV")  
    return 6619;
  else if( aa == "TKA")  
    return 6620;
  else if( aa == "TKR")  
    return 6621;
  else if( aa == "TKN")  
    return 6622;
  else if( aa == "TKD")  
    return 6623;
  else if( aa == "TKC")  
    return 6624;
  else if( aa == "TKQ")  
    return 6625;
  else if( aa == "TKE")  
    return 6626;
  else if( aa == "TKG")  
    return 6627;
  else if( aa == "TKH")  
    return 6628;
  else if( aa == "TKI")  
    return 6629;
  else if( aa == "TKL")  
    return 6630;
  else if( aa == "TKK")  
    return 6631;
  else if( aa == "TKM")  
    return 6632;
  else if( aa == "TKF")  
    return 6633;
  else if( aa == "TKP")  
    return 6634;
  else if( aa == "TKS")  
    return 6635;
  else if( aa == "TKT")  
    return 6636;
  else if( aa == "TKW")  
    return 6637;
  else if( aa == "TKY")  
    return 6638;
  else if( aa == "TKV")  
    return 6639;
  else if( aa == "TMA")  
    return 6640;
  else if( aa == "TMR")  
    return 6641;
  else if( aa == "TMN")  
    return 6642;
  else if( aa == "TMD")  
    return 6643;
  else if( aa == "TMC")  
    return 6644;
  else if( aa == "TMQ")  
    return 6645;
  else if( aa == "TME")  
    return 6646;
  else if( aa == "TMG")  
    return 6647;
  else if( aa == "TMH")  
    return 6648;
  else if( aa == "TMI")  
    return 6649;
  else if( aa == "TML")  
    return 6650;
  else if( aa == "TMK")  
    return 6651;
  else if( aa == "TMM")  
    return 6652;
  else if( aa == "TMF")  
    return 6653;
  else if( aa == "TMP")  
    return 6654;
  else if( aa == "TMS")  
    return 6655;
  else if( aa == "TMT")  
    return 6656;
  else if( aa == "TMW")  
    return 6657;
  else if( aa == "TMY")  
    return 6658;
  else if( aa == "TMV")  
    return 6659;
  else if( aa == "TFA")  
    return 6660;
  else if( aa == "TFR")  
    return 6661;
  else if( aa == "TFN")  
    return 6662;
  else if( aa == "TFD")  
    return 6663;
  else if( aa == "TFC")  
    return 6664;
  else if( aa == "TFQ")  
    return 6665;
  else if( aa == "TFE")  
    return 6666;
  else if( aa == "TFG")  
    return 6667;
  else if( aa == "TFH")  
    return 6668;
  else if( aa == "TFI")  
    return 6669;
  else if( aa == "TFL")  
    return 6670;
  else if( aa == "TFK")  
    return 6671;
  else if( aa == "TFM")  
    return 6672;
  else if( aa == "TFF")  
    return 6673;
  else if( aa == "TFP")  
    return 6674;
  else if( aa == "TFS")  
    return 6675;
  else if( aa == "TFT")  
    return 6676;
  else if( aa == "TFW")  
    return 6677;
  else if( aa == "TFY")  
    return 6678;
  else if( aa == "TFV")  
    return 6679;
  else if( aa == "TPA")  
    return 6680;
  else if( aa == "TPR")  
    return 6681;
  else if( aa == "TPN")  
    return 6682;
  else if( aa == "TPD")  
    return 6683;
  else if( aa == "TPC")  
    return 6684;
  else if( aa == "TPQ")  
    return 6685;
  else if( aa == "TPE")  
    return 6686;
  else if( aa == "TPG")  
    return 6687;
  else if( aa == "TPH")  
    return 6688;
  else if( aa == "TPI")  
    return 6689;
  else if( aa == "TPL")  
    return 6690;
  else if( aa == "TPK")  
    return 6691;
  else if( aa == "TPM")  
    return 6692;
  else if( aa == "TPF")  
    return 6693;
  else if( aa == "TPP")  
    return 6694;
  else if( aa == "TPS")  
    return 6695;
  else if( aa == "TPT")  
    return 6696;
  else if( aa == "TPW")  
    return 6697;
  else if( aa == "TPY")  
    return 6698;
  else if( aa == "TPV")  
    return 6699;
  else if( aa == "TSA")  
    return 6700;
  else if( aa == "TSR")  
    return 6701;
  else if( aa == "TSN")  
    return 6702;
  else if( aa == "TSD")  
    return 6703;
  else if( aa == "TSC")  
    return 6704;
  else if( aa == "TSQ")  
    return 6705;
  else if( aa == "TSE")  
    return 6706;
  else if( aa == "TSG")  
    return 6707;
  else if( aa == "TSH")  
    return 6708;
  else if( aa == "TSI")  
    return 6709;
  else if( aa == "TSL")  
    return 6710;
  else if( aa == "TSK")  
    return 6711;
  else if( aa == "TSM")  
    return 6712;
  else if( aa == "TSF")  
    return 6713;
  else if( aa == "TSP")  
    return 6714;
  else if( aa == "TSS")  
    return 6715;
  else if( aa == "TST")  
    return 6716;
  else if( aa == "TSW")  
    return 6717;
  else if( aa == "TSY")  
    return 6718;
  else if( aa == "TSV")  
    return 6719;
  else if( aa == "TTA")  
    return 6720;
  else if( aa == "TTR")  
    return 6721;
  else if( aa == "TTN")  
    return 6722;
  else if( aa == "TTD")  
    return 6723;
  else if( aa == "TTC")  
    return 6724;
  else if( aa == "TTQ")  
    return 6725;
  else if( aa == "TTE")  
    return 6726;
  else if( aa == "TTG")  
    return 6727;
  else if( aa == "TTH")  
    return 6728;
  else if( aa == "TTI")  
    return 6729;
  else if( aa == "TTL")  
    return 6730;
  else if( aa == "TTK")  
    return 6731;
  else if( aa == "TTM")  
    return 6732;
  else if( aa == "TTF")  
    return 6733;
  else if( aa == "TTP")  
    return 6734;
  else if( aa == "TTS")  
    return 6735;
  else if( aa == "TTT")  
    return 6736;
  else if( aa == "TTW")  
    return 6737;
  else if( aa == "TTY")  
    return 6738;
  else if( aa == "TTV")  
    return 6739;
  else if( aa == "TWA")  
    return 6740;
  else if( aa == "TWR")  
    return 6741;
  else if( aa == "TWN")  
    return 6742;
  else if( aa == "TWD")  
    return 6743;
  else if( aa == "TWC")  
    return 6744;
  else if( aa == "TWQ")  
    return 6745;
  else if( aa == "TWE")  
    return 6746;
  else if( aa == "TWG")  
    return 6747;
  else if( aa == "TWH")  
    return 6748;
  else if( aa == "TWI")  
    return 6749;
  else if( aa == "TWL")  
    return 6750;
  else if( aa == "TWK")  
    return 6751;
  else if( aa == "TWM")  
    return 6752;
  else if( aa == "TWF")  
    return 6753;
  else if( aa == "TWP")  
    return 6754;
  else if( aa == "TWS")  
    return 6755;
  else if( aa == "TWT")  
    return 6756;
  else if( aa == "TWW")  
    return 6757;
  else if( aa == "TWY")  
    return 6758;
  else if( aa == "TWV")  
    return 6759;
  else if( aa == "TYA")  
    return 6760;
  else if( aa == "TYR")  
    return 6761;
  else if( aa == "TYN")  
    return 6762;
  else if( aa == "TYD")  
    return 6763;
  else if( aa == "TYC")  
    return 6764;
  else if( aa == "TYQ")  
    return 6765;
  else if( aa == "TYE")  
    return 6766;
  else if( aa == "TYG")  
    return 6767;
  else if( aa == "TYH")  
    return 6768;
  else if( aa == "TYI")  
    return 6769;
  else if( aa == "TYL")  
    return 6770;
  else if( aa == "TYK")  
    return 6771;
  else if( aa == "TYM")  
    return 6772;
  else if( aa == "TYF")  
    return 6773;
  else if( aa == "TYP")  
    return 6774;
  else if( aa == "TYS")  
    return 6775;
  else if( aa == "TYT")  
    return 6776;
  else if( aa == "TYW")  
    return 6777;
  else if( aa == "TYY")  
    return 6778;
  else if( aa == "TYV")  
    return 6779;
  else if( aa == "TVA")  
    return 6780;
  else if( aa == "TVR")  
    return 6781;
  else if( aa == "TVN")  
    return 6782;
  else if( aa == "TVD")  
    return 6783;
  else if( aa == "TVC")  
    return 6784;
  else if( aa == "TVQ")  
    return 6785;
  else if( aa == "TVE")  
    return 6786;
  else if( aa == "TVG")  
    return 6787;
  else if( aa == "TVH")  
    return 6788;
  else if( aa == "TVI")  
    return 6789;
  else if( aa == "TVL")  
    return 6790;
  else if( aa == "TVK")  
    return 6791;
  else if( aa == "TVM")  
    return 6792;
  else if( aa == "TVF")  
    return 6793;
  else if( aa == "TVP")  
    return 6794;
  else if( aa == "TVS")  
    return 6795;
  else if( aa == "TVT")  
    return 6796;
  else if( aa == "TVW")  
    return 6797;
  else if( aa == "TVY")  
    return 6798;
  else if( aa == "TVV")  
    return 6799;
  else if( aa == "WAA")  
    return 6800;
  else if( aa == "WAR")  
    return 6801;
  else if( aa == "WAN")  
    return 6802;
  else if( aa == "WAD")  
    return 6803;
  else if( aa == "WAC")  
    return 6804;
  else if( aa == "WAQ")  
    return 6805;
  else if( aa == "WAE")  
    return 6806;
  else if( aa == "WAG")  
    return 6807;
  else if( aa == "WAH")  
    return 6808;
  else if( aa == "WAI")  
    return 6809;
  else if( aa == "WAL")  
    return 6810;
  else if( aa == "WAK")  
    return 6811;
  else if( aa == "WAM")  
    return 6812;
  else if( aa == "WAF")  
    return 6813;
  else if( aa == "WAP")  
    return 6814;
  else if( aa == "WAS")  
    return 6815;
  else if( aa == "WAT")  
    return 6816;
  else if( aa == "WAW")  
    return 6817;
  else if( aa == "WAY")  
    return 6818;
  else if( aa == "WAV")  
    return 6819;
  else if( aa == "WRA")  
    return 6820;
  else if( aa == "WRR")  
    return 6821;
  else if( aa == "WRN")  
    return 6822;
  else if( aa == "WRD")  
    return 6823;
  else if( aa == "WRC")  
    return 6824;
  else if( aa == "WRQ")  
    return 6825;
  else if( aa == "WRE")  
    return 6826;
  else if( aa == "WRG")  
    return 6827;
  else if( aa == "WRH")  
    return 6828;
  else if( aa == "WRI")  
    return 6829;
  else if( aa == "WRL")  
    return 6830;
  else if( aa == "WRK")  
    return 6831;
  else if( aa == "WRM")  
    return 6832;
  else if( aa == "WRF")  
    return 6833;
  else if( aa == "WRP")  
    return 6834;
  else if( aa == "WRS")  
    return 6835;
  else if( aa == "WRT")  
    return 6836;
  else if( aa == "WRW")  
    return 6837;
  else if( aa == "WRY")  
    return 6838;
  else if( aa == "WRV")  
    return 6839;
  else if( aa == "WNA")  
    return 6840;
  else if( aa == "WNR")  
    return 6841;
  else if( aa == "WNN")  
    return 6842;
  else if( aa == "WND")  
    return 6843;
  else if( aa == "WNC")  
    return 6844;
  else if( aa == "WNQ")  
    return 6845;
  else if( aa == "WNE")  
    return 6846;
  else if( aa == "WNG")  
    return 6847;
  else if( aa == "WNH")  
    return 6848;
  else if( aa == "WNI")  
    return 6849;
  else if( aa == "WNL")  
    return 6850;
  else if( aa == "WNK")  
    return 6851;
  else if( aa == "WNM")  
    return 6852;
  else if( aa == "WNF")  
    return 6853;
  else if( aa == "WNP")  
    return 6854;
  else if( aa == "WNS")  
    return 6855;
  else if( aa == "WNT")  
    return 6856;
  else if( aa == "WNW")  
    return 6857;
  else if( aa == "WNY")  
    return 6858;
  else if( aa == "WNV")  
    return 6859;
  else if( aa == "WDA")  
    return 6860;
  else if( aa == "WDR")  
    return 6861;
  else if( aa == "WDN")  
    return 6862;
  else if( aa == "WDD")  
    return 6863;
  else if( aa == "WDC")  
    return 6864;
  else if( aa == "WDQ")  
    return 6865;
  else if( aa == "WDE")  
    return 6866;
  else if( aa == "WDG")  
    return 6867;
  else if( aa == "WDH")  
    return 6868;
  else if( aa == "WDI")  
    return 6869;
  else if( aa == "WDL")  
    return 6870;
  else if( aa == "WDK")  
    return 6871;
  else if( aa == "WDM")  
    return 6872;
  else if( aa == "WDF")  
    return 6873;
  else if( aa == "WDP")  
    return 6874;
  else if( aa == "WDS")  
    return 6875;
  else if( aa == "WDT")  
    return 6876;
  else if( aa == "WDW")  
    return 6877;
  else if( aa == "WDY")  
    return 6878;
  else if( aa == "WDV")  
    return 6879;
  else if( aa == "WCA")  
    return 6880;
  else if( aa == "WCR")  
    return 6881;
  else if( aa == "WCN")  
    return 6882;
  else if( aa == "WCD")  
    return 6883;
  else if( aa == "WCC")  
    return 6884;
  else if( aa == "WCQ")  
    return 6885;
  else if( aa == "WCE")  
    return 6886;
  else if( aa == "WCG")  
    return 6887;
  else if( aa == "WCH")  
    return 6888;
  else if( aa == "WCI")  
    return 6889;
  else if( aa == "WCL")  
    return 6890;
  else if( aa == "WCK")  
    return 6891;
  else if( aa == "WCM")  
    return 6892;
  else if( aa == "WCF")  
    return 6893;
  else if( aa == "WCP")  
    return 6894;
  else if( aa == "WCS")  
    return 6895;
  else if( aa == "WCT")  
    return 6896;
  else if( aa == "WCW")  
    return 6897;
  else if( aa == "WCY")  
    return 6898;
  else if( aa == "WCV")  
    return 6899;
  else if( aa == "WQA")  
    return 6900;
  else if( aa == "WQR")  
    return 6901;
  else if( aa == "WQN")  
    return 6902;
  else if( aa == "WQD")  
    return 6903;
  else if( aa == "WQC")  
    return 6904;
  else if( aa == "WQQ")  
    return 6905;
  else if( aa == "WQE")  
    return 6906;
  else if( aa == "WQG")  
    return 6907;
  else if( aa == "WQH")  
    return 6908;
  else if( aa == "WQI")  
    return 6909;
  else if( aa == "WQL")  
    return 6910;
  else if( aa == "WQK")  
    return 6911;
  else if( aa == "WQM")  
    return 6912;
  else if( aa == "WQF")  
    return 6913;
  else if( aa == "WQP")  
    return 6914;
  else if( aa == "WQS")  
    return 6915;
  else if( aa == "WQT")  
    return 6916;
  else if( aa == "WQW")  
    return 6917;
  else if( aa == "WQY")  
    return 6918;
  else if( aa == "WQV")  
    return 6919;
  else if( aa == "WEA")  
    return 6920;
  else if( aa == "WER")  
    return 6921;
  else if( aa == "WEN")  
    return 6922;
  else if( aa == "WED")  
    return 6923;
  else if( aa == "WEC")  
    return 6924;
  else if( aa == "WEQ")  
    return 6925;
  else if( aa == "WEE")  
    return 6926;
  else if( aa == "WEG")  
    return 6927;
  else if( aa == "WEH")  
    return 6928;
  else if( aa == "WEI")  
    return 6929;
  else if( aa == "WEL")  
    return 6930;
  else if( aa == "WEK")  
    return 6931;
  else if( aa == "WEM")  
    return 6932;
  else if( aa == "WEF")  
    return 6933;
  else if( aa == "WEP")  
    return 6934;
  else if( aa == "WES")  
    return 6935;
  else if( aa == "WET")  
    return 6936;
  else if( aa == "WEW")  
    return 6937;
  else if( aa == "WEY")  
    return 6938;
  else if( aa == "WEV")  
    return 6939;
  else if( aa == "WGA")  
    return 6940;
  else if( aa == "WGR")  
    return 6941;
  else if( aa == "WGN")  
    return 6942;
  else if( aa == "WGD")  
    return 6943;
  else if( aa == "WGC")  
    return 6944;
  else if( aa == "WGQ")  
    return 6945;
  else if( aa == "WGE")  
    return 6946;
  else if( aa == "WGG")  
    return 6947;
  else if( aa == "WGH")  
    return 6948;
  else if( aa == "WGI")  
    return 6949;
  else if( aa == "WGL")  
    return 6950;
  else if( aa == "WGK")  
    return 6951;
  else if( aa == "WGM")  
    return 6952;
  else if( aa == "WGF")  
    return 6953;
  else if( aa == "WGP")  
    return 6954;
  else if( aa == "WGS")  
    return 6955;
  else if( aa == "WGT")  
    return 6956;
  else if( aa == "WGW")  
    return 6957;
  else if( aa == "WGY")  
    return 6958;
  else if( aa == "WGV")  
    return 6959;
  else if( aa == "WHA")  
    return 6960;
  else if( aa == "WHR")  
    return 6961;
  else if( aa == "WHN")  
    return 6962;
  else if( aa == "WHD")  
    return 6963;
  else if( aa == "WHC")  
    return 6964;
  else if( aa == "WHQ")  
    return 6965;
  else if( aa == "WHE")  
    return 6966;
  else if( aa == "WHG")  
    return 6967;
  else if( aa == "WHH")  
    return 6968;
  else if( aa == "WHI")  
    return 6969;
  else if( aa == "WHL")  
    return 6970;
  else if( aa == "WHK")  
    return 6971;
  else if( aa == "WHM")  
    return 6972;
  else if( aa == "WHF")  
    return 6973;
  else if( aa == "WHP")  
    return 6974;
  else if( aa == "WHS")  
    return 6975;
  else if( aa == "WHT")  
    return 6976;
  else if( aa == "WHW")  
    return 6977;
  else if( aa == "WHY")  
    return 6978;
  else if( aa == "WHV")  
    return 6979;
  else if( aa == "WIA")  
    return 6980;
  else if( aa == "WIR")  
    return 6981;
  else if( aa == "WIN")  
    return 6982;
  else if( aa == "WID")  
    return 6983;
  else if( aa == "WIC")  
    return 6984;
  else if( aa == "WIQ")  
    return 6985;
  else if( aa == "WIE")  
    return 6986;
  else if( aa == "WIG")  
    return 6987;
  else if( aa == "WIH")  
    return 6988;
  else if( aa == "WII")  
    return 6989;
  else if( aa == "WIL")  
    return 6990;
  else if( aa == "WIK")  
    return 6991;
  else if( aa == "WIM")  
    return 6992;
  else if( aa == "WIF")  
    return 6993;
  else if( aa == "WIP")  
    return 6994;
  else if( aa == "WIS")  
    return 6995;
  else if( aa == "WIT")  
    return 6996;
  else if( aa == "WIW")  
    return 6997;
  else if( aa == "WIY")  
    return 6998;
  else if( aa == "WIV")  
    return 6999;
  else if( aa == "WLA")  
    return 7000;
  else if( aa == "WLR")  
    return 7001;
  else if( aa == "WLN")  
    return 7002;
  else if( aa == "WLD")  
    return 7003;
  else if( aa == "WLC")  
    return 7004;
  else if( aa == "WLQ")  
    return 7005;
  else if( aa == "WLE")  
    return 7006;
  else if( aa == "WLG")  
    return 7007;
  else if( aa == "WLH")  
    return 7008;
  else if( aa == "WLI")  
    return 7009;
  else if( aa == "WLL")  
    return 7010;
  else if( aa == "WLK")  
    return 7011;
  else if( aa == "WLM")  
    return 7012;
  else if( aa == "WLF")  
    return 7013;
  else if( aa == "WLP")  
    return 7014;
  else if( aa == "WLS")  
    return 7015;
  else if( aa == "WLT")  
    return 7016;
  else if( aa == "WLW")  
    return 7017;
  else if( aa == "WLY")  
    return 7018;
  else if( aa == "WLV")  
    return 7019;
  else if( aa == "WKA")  
    return 7020;
  else if( aa == "WKR")  
    return 7021;
  else if( aa == "WKN")  
    return 7022;
  else if( aa == "WKD")  
    return 7023;
  else if( aa == "WKC")  
    return 7024;
  else if( aa == "WKQ")  
    return 7025;
  else if( aa == "WKE")  
    return 7026;
  else if( aa == "WKG")  
    return 7027;
  else if( aa == "WKH")  
    return 7028;
  else if( aa == "WKI")  
    return 7029;
  else if( aa == "WKL")  
    return 7030;
  else if( aa == "WKK")  
    return 7031;
  else if( aa == "WKM")  
    return 7032;
  else if( aa == "WKF")  
    return 7033;
  else if( aa == "WKP")  
    return 7034;
  else if( aa == "WKS")  
    return 7035;
  else if( aa == "WKT")  
    return 7036;
  else if( aa == "WKW")  
    return 7037;
  else if( aa == "WKY")  
    return 7038;
  else if( aa == "WKV")  
    return 7039;
  else if( aa == "WMA")  
    return 7040;
  else if( aa == "WMR")  
    return 7041;
  else if( aa == "WMN")  
    return 7042;
  else if( aa == "WMD")  
    return 7043;
  else if( aa == "WMC")  
    return 7044;
  else if( aa == "WMQ")  
    return 7045;
  else if( aa == "WME")  
    return 7046;
  else if( aa == "WMG")  
    return 7047;
  else if( aa == "WMH")  
    return 7048;
  else if( aa == "WMI")  
    return 7049;
  else if( aa == "WML")  
    return 7050;
  else if( aa == "WMK")  
    return 7051;
  else if( aa == "WMM")  
    return 7052;
  else if( aa == "WMF")  
    return 7053;
  else if( aa == "WMP")  
    return 7054;
  else if( aa == "WMS")  
    return 7055;
  else if( aa == "WMT")  
    return 7056;
  else if( aa == "WMW")  
    return 7057;
  else if( aa == "WMY")  
    return 7058;
  else if( aa == "WMV")  
    return 7059;
  else if( aa == "WFA")  
    return 7060;
  else if( aa == "WFR")  
    return 7061;
  else if( aa == "WFN")  
    return 7062;
  else if( aa == "WFD")  
    return 7063;
  else if( aa == "WFC")  
    return 7064;
  else if( aa == "WFQ")  
    return 7065;
  else if( aa == "WFE")  
    return 7066;
  else if( aa == "WFG")  
    return 7067;
  else if( aa == "WFH")  
    return 7068;
  else if( aa == "WFI")  
    return 7069;
  else if( aa == "WFL")  
    return 7070;
  else if( aa == "WFK")  
    return 7071;
  else if( aa == "WFM")  
    return 7072;
  else if( aa == "WFF")  
    return 7073;
  else if( aa == "WFP")  
    return 7074;
  else if( aa == "WFS")  
    return 7075;
  else if( aa == "WFT")  
    return 7076;
  else if( aa == "WFW")  
    return 7077;
  else if( aa == "WFY")  
    return 7078;
  else if( aa == "WFV")  
    return 7079;
  else if( aa == "WPA")  
    return 7080;
  else if( aa == "WPR")  
    return 7081;
  else if( aa == "WPN")  
    return 7082;
  else if( aa == "WPD")  
    return 7083;
  else if( aa == "WPC")  
    return 7084;
  else if( aa == "WPQ")  
    return 7085;
  else if( aa == "WPE")  
    return 7086;
  else if( aa == "WPG")  
    return 7087;
  else if( aa == "WPH")  
    return 7088;
  else if( aa == "WPI")  
    return 7089;
  else if( aa == "WPL")  
    return 7090;
  else if( aa == "WPK")  
    return 7091;
  else if( aa == "WPM")  
    return 7092;
  else if( aa == "WPF")  
    return 7093;
  else if( aa == "WPP")  
    return 7094;
  else if( aa == "WPS")  
    return 7095;
  else if( aa == "WPT")  
    return 7096;
  else if( aa == "WPW")  
    return 7097;
  else if( aa == "WPY")  
    return 7098;
  else if( aa == "WPV")  
    return 7099;
  else if( aa == "WSA")  
    return 7100;
  else if( aa == "WSR")  
    return 7101;
  else if( aa == "WSN")  
    return 7102;
  else if( aa == "WSD")  
    return 7103;
  else if( aa == "WSC")  
    return 7104;
  else if( aa == "WSQ")  
    return 7105;
  else if( aa == "WSE")  
    return 7106;
  else if( aa == "WSG")  
    return 7107;
  else if( aa == "WSH")  
    return 7108;
  else if( aa == "WSI")  
    return 7109;
  else if( aa == "WSL")  
    return 7110;
  else if( aa == "WSK")  
    return 7111;
  else if( aa == "WSM")  
    return 7112;
  else if( aa == "WSF")  
    return 7113;
  else if( aa == "WSP")  
    return 7114;
  else if( aa == "WSS")  
    return 7115;
  else if( aa == "WST")  
    return 7116;
  else if( aa == "WSW")  
    return 7117;
  else if( aa == "WSY")  
    return 7118;
  else if( aa == "WSV")  
    return 7119;
  else if( aa == "WTA")  
    return 7120;
  else if( aa == "WTR")  
    return 7121;
  else if( aa == "WTN")  
    return 7122;
  else if( aa == "WTD")  
    return 7123;
  else if( aa == "WTC")  
    return 7124;
  else if( aa == "WTQ")  
    return 7125;
  else if( aa == "WTE")  
    return 7126;
  else if( aa == "WTG")  
    return 7127;
  else if( aa == "WTH")  
    return 7128;
  else if( aa == "WTI")  
    return 7129;
  else if( aa == "WTL")  
    return 7130;
  else if( aa == "WTK")  
    return 7131;
  else if( aa == "WTM")  
    return 7132;
  else if( aa == "WTF")  
    return 7133;
  else if( aa == "WTP")  
    return 7134;
  else if( aa == "WTS")  
    return 7135;
  else if( aa == "WTT")  
    return 7136;
  else if( aa == "WTW")  
    return 7137;
  else if( aa == "WTY")  
    return 7138;
  else if( aa == "WTV")  
    return 7139;
  else if( aa == "WWA")  
    return 7140;
  else if( aa == "WWR")  
    return 7141;
  else if( aa == "WWN")  
    return 7142;
  else if( aa == "WWD")  
    return 7143;
  else if( aa == "WWC")  
    return 7144;
  else if( aa == "WWQ")  
    return 7145;
  else if( aa == "WWE")  
    return 7146;
  else if( aa == "WWG")  
    return 7147;
  else if( aa == "WWH")  
    return 7148;
  else if( aa == "WWI")  
    return 7149;
  else if( aa == "WWL")  
    return 7150;
  else if( aa == "WWK")  
    return 7151;
  else if( aa == "WWM")  
    return 7152;
  else if( aa == "WWF")  
    return 7153;
  else if( aa == "WWP")  
    return 7154;
  else if( aa == "WWS")  
    return 7155;
  else if( aa == "WWT")  
    return 7156;
  else if( aa == "WWW")  
    return 7157;
  else if( aa == "WWY")  
    return 7158;
  else if( aa == "WWV")  
    return 7159;
  else if( aa == "WYA")  
    return 7160;
  else if( aa == "WYR")  
    return 7161;
  else if( aa == "WYN")  
    return 7162;
  else if( aa == "WYD")  
    return 7163;
  else if( aa == "WYC")  
    return 7164;
  else if( aa == "WYQ")  
    return 7165;
  else if( aa == "WYE")  
    return 7166;
  else if( aa == "WYG")  
    return 7167;
  else if( aa == "WYH")  
    return 7168;
  else if( aa == "WYI")  
    return 7169;
  else if( aa == "WYL")  
    return 7170;
  else if( aa == "WYK")  
    return 7171;
  else if( aa == "WYM")  
    return 7172;
  else if( aa == "WYF")  
    return 7173;
  else if( aa == "WYP")  
    return 7174;
  else if( aa == "WYS")  
    return 7175;
  else if( aa == "WYT")  
    return 7176;
  else if( aa == "WYW")  
    return 7177;
  else if( aa == "WYY")  
    return 7178;
  else if( aa == "WYV")  
    return 7179;
  else if( aa == "WVA")  
    return 7180;
  else if( aa == "WVR")  
    return 7181;
  else if( aa == "WVN")  
    return 7182;
  else if( aa == "WVD")  
    return 7183;
  else if( aa == "WVC")  
    return 7184;
  else if( aa == "WVQ")  
    return 7185;
  else if( aa == "WVE")  
    return 7186;
  else if( aa == "WVG")  
    return 7187;
  else if( aa == "WVH")  
    return 7188;
  else if( aa == "WVI")  
    return 7189;
  else if( aa == "WVL")  
    return 7190;
  else if( aa == "WVK")  
    return 7191;
  else if( aa == "WVM")  
    return 7192;
  else if( aa == "WVF")  
    return 7193;
  else if( aa == "WVP")  
    return 7194;
  else if( aa == "WVS")  
    return 7195;
  else if( aa == "WVT")  
    return 7196;
  else if( aa == "WVW")  
    return 7197;
  else if( aa == "WVY")  
    return 7198;
  else if( aa == "WVV")  
    return 7199;
  else if( aa == "YAA")  
    return 7200;
  else if( aa == "YAR")  
    return 7201;
  else if( aa == "YAN")  
    return 7202;
  else if( aa == "YAD")  
    return 7203;
  else if( aa == "YAC")  
    return 7204;
  else if( aa == "YAQ")  
    return 7205;
  else if( aa == "YAE")  
    return 7206;
  else if( aa == "YAG")  
    return 7207;
  else if( aa == "YAH")  
    return 7208;
  else if( aa == "YAI")  
    return 7209;
  else if( aa == "YAL")  
    return 7210;
  else if( aa == "YAK")  
    return 7211;
  else if( aa == "YAM")  
    return 7212;
  else if( aa == "YAF")  
    return 7213;
  else if( aa == "YAP")  
    return 7214;
  else if( aa == "YAS")  
    return 7215;
  else if( aa == "YAT")  
    return 7216;
  else if( aa == "YAW")  
    return 7217;
  else if( aa == "YAY")  
    return 7218;
  else if( aa == "YAV")  
    return 7219;
  else if( aa == "YRA")  
    return 7220;
  else if( aa == "YRR")  
    return 7221;
  else if( aa == "YRN")  
    return 7222;
  else if( aa == "YRD")  
    return 7223;
  else if( aa == "YRC")  
    return 7224;
  else if( aa == "YRQ")  
    return 7225;
  else if( aa == "YRE")  
    return 7226;
  else if( aa == "YRG")  
    return 7227;
  else if( aa == "YRH")  
    return 7228;
  else if( aa == "YRI")  
    return 7229;
  else if( aa == "YRL")  
    return 7230;
  else if( aa == "YRK")  
    return 7231;
  else if( aa == "YRM")  
    return 7232;
  else if( aa == "YRF")  
    return 7233;
  else if( aa == "YRP")  
    return 7234;
  else if( aa == "YRS")  
    return 7235;
  else if( aa == "YRT")  
    return 7236;
  else if( aa == "YRW")  
    return 7237;
  else if( aa == "YRY")  
    return 7238;
  else if( aa == "YRV")  
    return 7239;
  else if( aa == "YNA")  
    return 7240;
  else if( aa == "YNR")  
    return 7241;
  else if( aa == "YNN")  
    return 7242;
  else if( aa == "YND")  
    return 7243;
  else if( aa == "YNC")  
    return 7244;
  else if( aa == "YNQ")  
    return 7245;
  else if( aa == "YNE")  
    return 7246;
  else if( aa == "YNG")  
    return 7247;
  else if( aa == "YNH")  
    return 7248;
  else if( aa == "YNI")  
    return 7249;
  else if( aa == "YNL")  
    return 7250;
  else if( aa == "YNK")  
    return 7251;
  else if( aa == "YNM")  
    return 7252;
  else if( aa == "YNF")  
    return 7253;
  else if( aa == "YNP")  
    return 7254;
  else if( aa == "YNS")  
    return 7255;
  else if( aa == "YNT")  
    return 7256;
  else if( aa == "YNW")  
    return 7257;
  else if( aa == "YNY")  
    return 7258;
  else if( aa == "YNV")  
    return 7259;
  else if( aa == "YDA")  
    return 7260;
  else if( aa == "YDR")  
    return 7261;
  else if( aa == "YDN")  
    return 7262;
  else if( aa == "YDD")  
    return 7263;
  else if( aa == "YDC")  
    return 7264;
  else if( aa == "YDQ")  
    return 7265;
  else if( aa == "YDE")  
    return 7266;
  else if( aa == "YDG")  
    return 7267;
  else if( aa == "YDH")  
    return 7268;
  else if( aa == "YDI")  
    return 7269;
  else if( aa == "YDL")  
    return 7270;
  else if( aa == "YDK")  
    return 7271;
  else if( aa == "YDM")  
    return 7272;
  else if( aa == "YDF")  
    return 7273;
  else if( aa == "YDP")  
    return 7274;
  else if( aa == "YDS")  
    return 7275;
  else if( aa == "YDT")  
    return 7276;
  else if( aa == "YDW")  
    return 7277;
  else if( aa == "YDY")  
    return 7278;
  else if( aa == "YDV")  
    return 7279;
  else if( aa == "YCA")  
    return 7280;
  else if( aa == "YCR")  
    return 7281;
  else if( aa == "YCN")  
    return 7282;
  else if( aa == "YCD")  
    return 7283;
  else if( aa == "YCC")  
    return 7284;
  else if( aa == "YCQ")  
    return 7285;
  else if( aa == "YCE")  
    return 7286;
  else if( aa == "YCG")  
    return 7287;
  else if( aa == "YCH")  
    return 7288;
  else if( aa == "YCI")  
    return 7289;
  else if( aa == "YCL")  
    return 7290;
  else if( aa == "YCK")  
    return 7291;
  else if( aa == "YCM")  
    return 7292;
  else if( aa == "YCF")  
    return 7293;
  else if( aa == "YCP")  
    return 7294;
  else if( aa == "YCS")  
    return 7295;
  else if( aa == "YCT")  
    return 7296;
  else if( aa == "YCW")  
    return 7297;
  else if( aa == "YCY")  
    return 7298;
  else if( aa == "YCV")  
    return 7299;
  else if( aa == "YQA")  
    return 7300;
  else if( aa == "YQR")  
    return 7301;
  else if( aa == "YQN")  
    return 7302;
  else if( aa == "YQD")  
    return 7303;
  else if( aa == "YQC")  
    return 7304;
  else if( aa == "YQQ")  
    return 7305;
  else if( aa == "YQE")  
    return 7306;
  else if( aa == "YQG")  
    return 7307;
  else if( aa == "YQH")  
    return 7308;
  else if( aa == "YQI")  
    return 7309;
  else if( aa == "YQL")  
    return 7310;
  else if( aa == "YQK")  
    return 7311;
  else if( aa == "YQM")  
    return 7312;
  else if( aa == "YQF")  
    return 7313;
  else if( aa == "YQP")  
    return 7314;
  else if( aa == "YQS")  
    return 7315;
  else if( aa == "YQT")  
    return 7316;
  else if( aa == "YQW")  
    return 7317;
  else if( aa == "YQY")  
    return 7318;
  else if( aa == "YQV")  
    return 7319;
  else if( aa == "YEA")  
    return 7320;
  else if( aa == "YER")  
    return 7321;
  else if( aa == "YEN")  
    return 7322;
  else if( aa == "YED")  
    return 7323;
  else if( aa == "YEC")  
    return 7324;
  else if( aa == "YEQ")  
    return 7325;
  else if( aa == "YEE")  
    return 7326;
  else if( aa == "YEG")  
    return 7327;
  else if( aa == "YEH")  
    return 7328;
  else if( aa == "YEI")  
    return 7329;
  else if( aa == "YEL")  
    return 7330;
  else if( aa == "YEK")  
    return 7331;
  else if( aa == "YEM")  
    return 7332;
  else if( aa == "YEF")  
    return 7333;
  else if( aa == "YEP")  
    return 7334;
  else if( aa == "YES")  
    return 7335;
  else if( aa == "YET")  
    return 7336;
  else if( aa == "YEW")  
    return 7337;
  else if( aa == "YEY")  
    return 7338;
  else if( aa == "YEV")  
    return 7339;
  else if( aa == "YGA")  
    return 7340;
  else if( aa == "YGR")  
    return 7341;
  else if( aa == "YGN")  
    return 7342;
  else if( aa == "YGD")  
    return 7343;
  else if( aa == "YGC")  
    return 7344;
  else if( aa == "YGQ")  
    return 7345;
  else if( aa == "YGE")  
    return 7346;
  else if( aa == "YGG")  
    return 7347;
  else if( aa == "YGH")  
    return 7348;
  else if( aa == "YGI")  
    return 7349;
  else if( aa == "YGL")  
    return 7350;
  else if( aa == "YGK")  
    return 7351;
  else if( aa == "YGM")  
    return 7352;
  else if( aa == "YGF")  
    return 7353;
  else if( aa == "YGP")  
    return 7354;
  else if( aa == "YGS")  
    return 7355;
  else if( aa == "YGT")  
    return 7356;
  else if( aa == "YGW")  
    return 7357;
  else if( aa == "YGY")  
    return 7358;
  else if( aa == "YGV")  
    return 7359;
  else if( aa == "YHA")  
    return 7360;
  else if( aa == "YHR")  
    return 7361;
  else if( aa == "YHN")  
    return 7362;
  else if( aa == "YHD")  
    return 7363;
  else if( aa == "YHC")  
    return 7364;
  else if( aa == "YHQ")  
    return 7365;
  else if( aa == "YHE")  
    return 7366;
  else if( aa == "YHG")  
    return 7367;
  else if( aa == "YHH")  
    return 7368;
  else if( aa == "YHI")  
    return 7369;
  else if( aa == "YHL")  
    return 7370;
  else if( aa == "YHK")  
    return 7371;
  else if( aa == "YHM")  
    return 7372;
  else if( aa == "YHF")  
    return 7373;
  else if( aa == "YHP")  
    return 7374;
  else if( aa == "YHS")  
    return 7375;
  else if( aa == "YHT")  
    return 7376;
  else if( aa == "YHW")  
    return 7377;
  else if( aa == "YHY")  
    return 7378;
  else if( aa == "YHV")  
    return 7379;
  else if( aa == "YIA")  
    return 7380;
  else if( aa == "YIR")  
    return 7381;
  else if( aa == "YIN")  
    return 7382;
  else if( aa == "YID")  
    return 7383;
  else if( aa == "YIC")  
    return 7384;
  else if( aa == "YIQ")  
    return 7385;
  else if( aa == "YIE")  
    return 7386;
  else if( aa == "YIG")  
    return 7387;
  else if( aa == "YIH")  
    return 7388;
  else if( aa == "YII")  
    return 7389;
  else if( aa == "YIL")  
    return 7390;
  else if( aa == "YIK")  
    return 7391;
  else if( aa == "YIM")  
    return 7392;
  else if( aa == "YIF")  
    return 7393;
  else if( aa == "YIP")  
    return 7394;
  else if( aa == "YIS")  
    return 7395;
  else if( aa == "YIT")  
    return 7396;
  else if( aa == "YIW")  
    return 7397;
  else if( aa == "YIY")  
    return 7398;
  else if( aa == "YIV")  
    return 7399;
  else if( aa == "YLA")  
    return 7400;
  else if( aa == "YLR")  
    return 7401;
  else if( aa == "YLN")  
    return 7402;
  else if( aa == "YLD")  
    return 7403;
  else if( aa == "YLC")  
    return 7404;
  else if( aa == "YLQ")  
    return 7405;
  else if( aa == "YLE")  
    return 7406;
  else if( aa == "YLG")  
    return 7407;
  else if( aa == "YLH")  
    return 7408;
  else if( aa == "YLI")  
    return 7409;
  else if( aa == "YLL")  
    return 7410;
  else if( aa == "YLK")  
    return 7411;
  else if( aa == "YLM")  
    return 7412;
  else if( aa == "YLF")  
    return 7413;
  else if( aa == "YLP")  
    return 7414;
  else if( aa == "YLS")  
    return 7415;
  else if( aa == "YLT")  
    return 7416;
  else if( aa == "YLW")  
    return 7417;
  else if( aa == "YLY")  
    return 7418;
  else if( aa == "YLV")  
    return 7419;
  else if( aa == "YKA")  
    return 7420;
  else if( aa == "YKR")  
    return 7421;
  else if( aa == "YKN")  
    return 7422;
  else if( aa == "YKD")  
    return 7423;
  else if( aa == "YKC")  
    return 7424;
  else if( aa == "YKQ")  
    return 7425;
  else if( aa == "YKE")  
    return 7426;
  else if( aa == "YKG")  
    return 7427;
  else if( aa == "YKH")  
    return 7428;
  else if( aa == "YKI")  
    return 7429;
  else if( aa == "YKL")  
    return 7430;
  else if( aa == "YKK")  
    return 7431;
  else if( aa == "YKM")  
    return 7432;
  else if( aa == "YKF")  
    return 7433;
  else if( aa == "YKP")  
    return 7434;
  else if( aa == "YKS")  
    return 7435;
  else if( aa == "YKT")  
    return 7436;
  else if( aa == "YKW")  
    return 7437;
  else if( aa == "YKY")  
    return 7438;
  else if( aa == "YKV")  
    return 7439;
  else if( aa == "YMA")  
    return 7440;
  else if( aa == "YMR")  
    return 7441;
  else if( aa == "YMN")  
    return 7442;
  else if( aa == "YMD")  
    return 7443;
  else if( aa == "YMC")  
    return 7444;
  else if( aa == "YMQ")  
    return 7445;
  else if( aa == "YME")  
    return 7446;
  else if( aa == "YMG")  
    return 7447;
  else if( aa == "YMH")  
    return 7448;
  else if( aa == "YMI")  
    return 7449;
  else if( aa == "YML")  
    return 7450;
  else if( aa == "YMK")  
    return 7451;
  else if( aa == "YMM")  
    return 7452;
  else if( aa == "YMF")  
    return 7453;
  else if( aa == "YMP")  
    return 7454;
  else if( aa == "YMS")  
    return 7455;
  else if( aa == "YMT")  
    return 7456;
  else if( aa == "YMW")  
    return 7457;
  else if( aa == "YMY")  
    return 7458;
  else if( aa == "YMV")  
    return 7459;
  else if( aa == "YFA")  
    return 7460;
  else if( aa == "YFR")  
    return 7461;
  else if( aa == "YFN")  
    return 7462;
  else if( aa == "YFD")  
    return 7463;
  else if( aa == "YFC")  
    return 7464;
  else if( aa == "YFQ")  
    return 7465;
  else if( aa == "YFE")  
    return 7466;
  else if( aa == "YFG")  
    return 7467;
  else if( aa == "YFH")  
    return 7468;
  else if( aa == "YFI")  
    return 7469;
  else if( aa == "YFL")  
    return 7470;
  else if( aa == "YFK")  
    return 7471;
  else if( aa == "YFM")  
    return 7472;
  else if( aa == "YFF")  
    return 7473;
  else if( aa == "YFP")  
    return 7474;
  else if( aa == "YFS")  
    return 7475;
  else if( aa == "YFT")  
    return 7476;
  else if( aa == "YFW")  
    return 7477;
  else if( aa == "YFY")  
    return 7478;
  else if( aa == "YFV")  
    return 7479;
  else if( aa == "YPA")  
    return 7480;
  else if( aa == "YPR")  
    return 7481;
  else if( aa == "YPN")  
    return 7482;
  else if( aa == "YPD")  
    return 7483;
  else if( aa == "YPC")  
    return 7484;
  else if( aa == "YPQ")  
    return 7485;
  else if( aa == "YPE")  
    return 7486;
  else if( aa == "YPG")  
    return 7487;
  else if( aa == "YPH")  
    return 7488;
  else if( aa == "YPI")  
    return 7489;
  else if( aa == "YPL")  
    return 7490;
  else if( aa == "YPK")  
    return 7491;
  else if( aa == "YPM")  
    return 7492;
  else if( aa == "YPF")  
    return 7493;
  else if( aa == "YPP")  
    return 7494;
  else if( aa == "YPS")  
    return 7495;
  else if( aa == "YPT")  
    return 7496;
  else if( aa == "YPW")  
    return 7497;
  else if( aa == "YPY")  
    return 7498;
  else if( aa == "YPV")  
    return 7499;
  else if( aa == "YSA")  
    return 7500;
  else if( aa == "YSR")  
    return 7501;
  else if( aa == "YSN")  
    return 7502;
  else if( aa == "YSD")  
    return 7503;
  else if( aa == "YSC")  
    return 7504;
  else if( aa == "YSQ")  
    return 7505;
  else if( aa == "YSE")  
    return 7506;
  else if( aa == "YSG")  
    return 7507;
  else if( aa == "YSH")  
    return 7508;
  else if( aa == "YSI")  
    return 7509;
  else if( aa == "YSL")  
    return 7510;
  else if( aa == "YSK")  
    return 7511;
  else if( aa == "YSM")  
    return 7512;
  else if( aa == "YSF")  
    return 7513;
  else if( aa == "YSP")  
    return 7514;
  else if( aa == "YSS")  
    return 7515;
  else if( aa == "YST")  
    return 7516;
  else if( aa == "YSW")  
    return 7517;
  else if( aa == "YSY")  
    return 7518;
  else if( aa == "YSV")  
    return 7519;
  else if( aa == "YTA")  
    return 7520;
  else if( aa == "YTR")  
    return 7521;
  else if( aa == "YTN")  
    return 7522;
  else if( aa == "YTD")  
    return 7523;
  else if( aa == "YTC")  
    return 7524;
  else if( aa == "YTQ")  
    return 7525;
  else if( aa == "YTE")  
    return 7526;
  else if( aa == "YTG")  
    return 7527;
  else if( aa == "YTH")  
    return 7528;
  else if( aa == "YTI")  
    return 7529;
  else if( aa == "YTL")  
    return 7530;
  else if( aa == "YTK")  
    return 7531;
  else if( aa == "YTM")  
    return 7532;
  else if( aa == "YTF")  
    return 7533;
  else if( aa == "YTP")  
    return 7534;
  else if( aa == "YTS")  
    return 7535;
  else if( aa == "YTT")  
    return 7536;
  else if( aa == "YTW")  
    return 7537;
  else if( aa == "YTY")  
    return 7538;
  else if( aa == "YTV")  
    return 7539;
  else if( aa == "YWA")  
    return 7540;
  else if( aa == "YWR")  
    return 7541;
  else if( aa == "YWN")  
    return 7542;
  else if( aa == "YWD")  
    return 7543;
  else if( aa == "YWC")  
    return 7544;
  else if( aa == "YWQ")  
    return 7545;
  else if( aa == "YWE")  
    return 7546;
  else if( aa == "YWG")  
    return 7547;
  else if( aa == "YWH")  
    return 7548;
  else if( aa == "YWI")  
    return 7549;
  else if( aa == "YWL")  
    return 7550;
  else if( aa == "YWK")  
    return 7551;
  else if( aa == "YWM")  
    return 7552;
  else if( aa == "YWF")  
    return 7553;
  else if( aa == "YWP")  
    return 7554;
  else if( aa == "YWS")  
    return 7555;
  else if( aa == "YWT")  
    return 7556;
  else if( aa == "YWW")  
    return 7557;
  else if( aa == "YWY")  
    return 7558;
  else if( aa == "YWV")  
    return 7559;
  else if( aa == "YYA")  
    return 7560;
  else if( aa == "YYR")  
    return 7561;
  else if( aa == "YYN")  
    return 7562;
  else if( aa == "YYD")  
    return 7563;
  else if( aa == "YYC")  
    return 7564;
  else if( aa == "YYQ")  
    return 7565;
  else if( aa == "YYE")  
    return 7566;
  else if( aa == "YYG")  
    return 7567;
  else if( aa == "YYH")  
    return 7568;
  else if( aa == "YYI")  
    return 7569;
  else if( aa == "YYL")  
    return 7570;
  else if( aa == "YYK")  
    return 7571;
  else if( aa == "YYM")  
    return 7572;
  else if( aa == "YYF")  
    return 7573;
  else if( aa == "YYP")  
    return 7574;
  else if( aa == "YYS")  
    return 7575;
  else if( aa == "YYT")  
    return 7576;
  else if( aa == "YYW")  
    return 7577;
  else if( aa == "YYY")  
    return 7578;
  else if( aa == "YYV")  
    return 7579;
  else if( aa == "YVA")  
    return 7580;
  else if( aa == "YVR")  
    return 7581;
  else if( aa == "YVN")  
    return 7582;
  else if( aa == "YVD")  
    return 7583;
  else if( aa == "YVC")  
    return 7584;
  else if( aa == "YVQ")  
    return 7585;
  else if( aa == "YVE")  
    return 7586;
  else if( aa == "YVG")  
    return 7587;
  else if( aa == "YVH")  
    return 7588;
  else if( aa == "YVI")  
    return 7589;
  else if( aa == "YVL")  
    return 7590;
  else if( aa == "YVK")  
    return 7591;
  else if( aa == "YVM")  
    return 7592;
  else if( aa == "YVF")  
    return 7593;
  else if( aa == "YVP")  
    return 7594;
  else if( aa == "YVS")  
    return 7595;
  else if( aa == "YVT")  
    return 7596;
  else if( aa == "YVW")  
    return 7597;
  else if( aa == "YVY")  
    return 7598;
  else if( aa == "YVV")  
    return 7599;
  else if( aa == "VAA")  
    return 7600;
  else if( aa == "VAR")  
    return 7601;
  else if( aa == "VAN")  
    return 7602;
  else if( aa == "VAD")  
    return 7603;
  else if( aa == "VAC")  
    return 7604;
  else if( aa == "VAQ")  
    return 7605;
  else if( aa == "VAE")  
    return 7606;
  else if( aa == "VAG")  
    return 7607;
  else if( aa == "VAH")  
    return 7608;
  else if( aa == "VAI")  
    return 7609;
  else if( aa == "VAL")  
    return 7610;
  else if( aa == "VAK")  
    return 7611;
  else if( aa == "VAM")  
    return 7612;
  else if( aa == "VAF")  
    return 7613;
  else if( aa == "VAP")  
    return 7614;
  else if( aa == "VAS")  
    return 7615;
  else if( aa == "VAT")  
    return 7616;
  else if( aa == "VAW")  
    return 7617;
  else if( aa == "VAY")  
    return 7618;
  else if( aa == "VAV")  
    return 7619;
  else if( aa == "VRA")  
    return 7620;
  else if( aa == "VRR")  
    return 7621;
  else if( aa == "VRN")  
    return 7622;
  else if( aa == "VRD")  
    return 7623;
  else if( aa == "VRC")  
    return 7624;
  else if( aa == "VRQ")  
    return 7625;
  else if( aa == "VRE")  
    return 7626;
  else if( aa == "VRG")  
    return 7627;
  else if( aa == "VRH")  
    return 7628;
  else if( aa == "VRI")  
    return 7629;
  else if( aa == "VRL")  
    return 7630;
  else if( aa == "VRK")  
    return 7631;
  else if( aa == "VRM")  
    return 7632;
  else if( aa == "VRF")  
    return 7633;
  else if( aa == "VRP")  
    return 7634;
  else if( aa == "VRS")  
    return 7635;
  else if( aa == "VRT")  
    return 7636;
  else if( aa == "VRW")  
    return 7637;
  else if( aa == "VRY")  
    return 7638;
  else if( aa == "VRV")  
    return 7639;
  else if( aa == "VNA")  
    return 7640;
  else if( aa == "VNR")  
    return 7641;
  else if( aa == "VNN")  
    return 7642;
  else if( aa == "VND")  
    return 7643;
  else if( aa == "VNC")  
    return 7644;
  else if( aa == "VNQ")  
    return 7645;
  else if( aa == "VNE")  
    return 7646;
  else if( aa == "VNG")  
    return 7647;
  else if( aa == "VNH")  
    return 7648;
  else if( aa == "VNI")  
    return 7649;
  else if( aa == "VNL")  
    return 7650;
  else if( aa == "VNK")  
    return 7651;
  else if( aa == "VNM")  
    return 7652;
  else if( aa == "VNF")  
    return 7653;
  else if( aa == "VNP")  
    return 7654;
  else if( aa == "VNS")  
    return 7655;
  else if( aa == "VNT")  
    return 7656;
  else if( aa == "VNW")  
    return 7657;
  else if( aa == "VNY")  
    return 7658;
  else if( aa == "VNV")  
    return 7659;
  else if( aa == "VDA")  
    return 7660;
  else if( aa == "VDR")  
    return 7661;
  else if( aa == "VDN")  
    return 7662;
  else if( aa == "VDD")  
    return 7663;
  else if( aa == "VDC")  
    return 7664;
  else if( aa == "VDQ")  
    return 7665;
  else if( aa == "VDE")  
    return 7666;
  else if( aa == "VDG")  
    return 7667;
  else if( aa == "VDH")  
    return 7668;
  else if( aa == "VDI")  
    return 7669;
  else if( aa == "VDL")  
    return 7670;
  else if( aa == "VDK")  
    return 7671;
  else if( aa == "VDM")  
    return 7672;
  else if( aa == "VDF")  
    return 7673;
  else if( aa == "VDP")  
    return 7674;
  else if( aa == "VDS")  
    return 7675;
  else if( aa == "VDT")  
    return 7676;
  else if( aa == "VDW")  
    return 7677;
  else if( aa == "VDY")  
    return 7678;
  else if( aa == "VDV")  
    return 7679;
  else if( aa == "VCA")  
    return 7680;
  else if( aa == "VCR")  
    return 7681;
  else if( aa == "VCN")  
    return 7682;
  else if( aa == "VCD")  
    return 7683;
  else if( aa == "VCC")  
    return 7684;
  else if( aa == "VCQ")  
    return 7685;
  else if( aa == "VCE")  
    return 7686;
  else if( aa == "VCG")  
    return 7687;
  else if( aa == "VCH")  
    return 7688;
  else if( aa == "VCI")  
    return 7689;
  else if( aa == "VCL")  
    return 7690;
  else if( aa == "VCK")  
    return 7691;
  else if( aa == "VCM")  
    return 7692;
  else if( aa == "VCF")  
    return 7693;
  else if( aa == "VCP")  
    return 7694;
  else if( aa == "VCS")  
    return 7695;
  else if( aa == "VCT")  
    return 7696;
  else if( aa == "VCW")  
    return 7697;
  else if( aa == "VCY")  
    return 7698;
  else if( aa == "VCV")  
    return 7699;
  else if( aa == "VQA")  
    return 7700;
  else if( aa == "VQR")  
    return 7701;
  else if( aa == "VQN")  
    return 7702;
  else if( aa == "VQD")  
    return 7703;
  else if( aa == "VQC")  
    return 7704;
  else if( aa == "VQQ")  
    return 7705;
  else if( aa == "VQE")  
    return 7706;
  else if( aa == "VQG")  
    return 7707;
  else if( aa == "VQH")  
    return 7708;
  else if( aa == "VQI")  
    return 7709;
  else if( aa == "VQL")  
    return 7710;
  else if( aa == "VQK")  
    return 7711;
  else if( aa == "VQM")  
    return 7712;
  else if( aa == "VQF")  
    return 7713;
  else if( aa == "VQP")  
    return 7714;
  else if( aa == "VQS")  
    return 7715;
  else if( aa == "VQT")  
    return 7716;
  else if( aa == "VQW")  
    return 7717;
  else if( aa == "VQY")  
    return 7718;
  else if( aa == "VQV")  
    return 7719;
  else if( aa == "VEA")  
    return 7720;
  else if( aa == "VER")  
    return 7721;
  else if( aa == "VEN")  
    return 7722;
  else if( aa == "VED")  
    return 7723;
  else if( aa == "VEC")  
    return 7724;
  else if( aa == "VEQ")  
    return 7725;
  else if( aa == "VEE")  
    return 7726;
  else if( aa == "VEG")  
    return 7727;
  else if( aa == "VEH")  
    return 7728;
  else if( aa == "VEI")  
    return 7729;
  else if( aa == "VEL")  
    return 7730;
  else if( aa == "VEK")  
    return 7731;
  else if( aa == "VEM")  
    return 7732;
  else if( aa == "VEF")  
    return 7733;
  else if( aa == "VEP")  
    return 7734;
  else if( aa == "VES")  
    return 7735;
  else if( aa == "VET")  
    return 7736;
  else if( aa == "VEW")  
    return 7737;
  else if( aa == "VEY")  
    return 7738;
  else if( aa == "VEV")  
    return 7739;
  else if( aa == "VGA")  
    return 7740;
  else if( aa == "VGR")  
    return 7741;
  else if( aa == "VGN")  
    return 7742;
  else if( aa == "VGD")  
    return 7743;
  else if( aa == "VGC")  
    return 7744;
  else if( aa == "VGQ")  
    return 7745;
  else if( aa == "VGE")  
    return 7746;
  else if( aa == "VGG")  
    return 7747;
  else if( aa == "VGH")  
    return 7748;
  else if( aa == "VGI")  
    return 7749;
  else if( aa == "VGL")  
    return 7750;
  else if( aa == "VGK")  
    return 7751;
  else if( aa == "VGM")  
    return 7752;
  else if( aa == "VGF")  
    return 7753;
  else if( aa == "VGP")  
    return 7754;
  else if( aa == "VGS")  
    return 7755;
  else if( aa == "VGT")  
    return 7756;
  else if( aa == "VGW")  
    return 7757;
  else if( aa == "VGY")  
    return 7758;
  else if( aa == "VGV")  
    return 7759;
  else if( aa == "VHA")  
    return 7760;
  else if( aa == "VHR")  
    return 7761;
  else if( aa == "VHN")  
    return 7762;
  else if( aa == "VHD")  
    return 7763;
  else if( aa == "VHC")  
    return 7764;
  else if( aa == "VHQ")  
    return 7765;
  else if( aa == "VHE")  
    return 7766;
  else if( aa == "VHG")  
    return 7767;
  else if( aa == "VHH")  
    return 7768;
  else if( aa == "VHI")  
    return 7769;
  else if( aa == "VHL")  
    return 7770;
  else if( aa == "VHK")  
    return 7771;
  else if( aa == "VHM")  
    return 7772;
  else if( aa == "VHF")  
    return 7773;
  else if( aa == "VHP")  
    return 7774;
  else if( aa == "VHS")  
    return 7775;
  else if( aa == "VHT")  
    return 7776;
  else if( aa == "VHW")  
    return 7777;
  else if( aa == "VHY")  
    return 7778;
  else if( aa == "VHV")  
    return 7779;
  else if( aa == "VIA")  
    return 7780;
  else if( aa == "VIR")  
    return 7781;
  else if( aa == "VIN")  
    return 7782;
  else if( aa == "VID")  
    return 7783;
  else if( aa == "VIC")  
    return 7784;
  else if( aa == "VIQ")  
    return 7785;
  else if( aa == "VIE")  
    return 7786;
  else if( aa == "VIG")  
    return 7787;
  else if( aa == "VIH")  
    return 7788;
  else if( aa == "VII")  
    return 7789;
  else if( aa == "VIL")  
    return 7790;
  else if( aa == "VIK")  
    return 7791;
  else if( aa == "VIM")  
    return 7792;
  else if( aa == "VIF")  
    return 7793;
  else if( aa == "VIP")  
    return 7794;
  else if( aa == "VIS")  
    return 7795;
  else if( aa == "VIT")  
    return 7796;
  else if( aa == "VIW")  
    return 7797;
  else if( aa == "VIY")  
    return 7798;
  else if( aa == "VIV")  
    return 7799;
  else if( aa == "VLA")  
    return 7800;
  else if( aa == "VLR")  
    return 7801;
  else if( aa == "VLN")  
    return 7802;
  else if( aa == "VLD")  
    return 7803;
  else if( aa == "VLC")  
    return 7804;
  else if( aa == "VLQ")  
    return 7805;
  else if( aa == "VLE")  
    return 7806;
  else if( aa == "VLG")  
    return 7807;
  else if( aa == "VLH")  
    return 7808;
  else if( aa == "VLI")  
    return 7809;
  else if( aa == "VLL")  
    return 7810;
  else if( aa == "VLK")  
    return 7811;
  else if( aa == "VLM")  
    return 7812;
  else if( aa == "VLF")  
    return 7813;
  else if( aa == "VLP")  
    return 7814;
  else if( aa == "VLS")  
    return 7815;
  else if( aa == "VLT")  
    return 7816;
  else if( aa == "VLW")  
    return 7817;
  else if( aa == "VLY")  
    return 7818;
  else if( aa == "VLV")  
    return 7819;
  else if( aa == "VKA")  
    return 7820;
  else if( aa == "VKR")  
    return 7821;
  else if( aa == "VKN")  
    return 7822;
  else if( aa == "VKD")  
    return 7823;
  else if( aa == "VKC")  
    return 7824;
  else if( aa == "VKQ")  
    return 7825;
  else if( aa == "VKE")  
    return 7826;
  else if( aa == "VKG")  
    return 7827;
  else if( aa == "VKH")  
    return 7828;
  else if( aa == "VKI")  
    return 7829;
  else if( aa == "VKL")  
    return 7830;
  else if( aa == "VKK")  
    return 7831;
  else if( aa == "VKM")  
    return 7832;
  else if( aa == "VKF")  
    return 7833;
  else if( aa == "VKP")  
    return 7834;
  else if( aa == "VKS")  
    return 7835;
  else if( aa == "VKT")  
    return 7836;
  else if( aa == "VKW")  
    return 7837;
  else if( aa == "VKY")  
    return 7838;
  else if( aa == "VKV")  
    return 7839;
  else if( aa == "VMA")  
    return 7840;
  else if( aa == "VMR")  
    return 7841;
  else if( aa == "VMN")  
    return 7842;
  else if( aa == "VMD")  
    return 7843;
  else if( aa == "VMC")  
    return 7844;
  else if( aa == "VMQ")  
    return 7845;
  else if( aa == "VME")  
    return 7846;
  else if( aa == "VMG")  
    return 7847;
  else if( aa == "VMH")  
    return 7848;
  else if( aa == "VMI")  
    return 7849;
  else if( aa == "VML")  
    return 7850;
  else if( aa == "VMK")  
    return 7851;
  else if( aa == "VMM")  
    return 7852;
  else if( aa == "VMF")  
    return 7853;
  else if( aa == "VMP")  
    return 7854;
  else if( aa == "VMS")  
    return 7855;
  else if( aa == "VMT")  
    return 7856;
  else if( aa == "VMW")  
    return 7857;
  else if( aa == "VMY")  
    return 7858;
  else if( aa == "VMV")  
    return 7859;
  else if( aa == "VFA")  
    return 7860;
  else if( aa == "VFR")  
    return 7861;
  else if( aa == "VFN")  
    return 7862;
  else if( aa == "VFD")  
    return 7863;
  else if( aa == "VFC")  
    return 7864;
  else if( aa == "VFQ")  
    return 7865;
  else if( aa == "VFE")  
    return 7866;
  else if( aa == "VFG")  
    return 7867;
  else if( aa == "VFH")  
    return 7868;
  else if( aa == "VFI")  
    return 7869;
  else if( aa == "VFL")  
    return 7870;
  else if( aa == "VFK")  
    return 7871;
  else if( aa == "VFM")  
    return 7872;
  else if( aa == "VFF")  
    return 7873;
  else if( aa == "VFP")  
    return 7874;
  else if( aa == "VFS")  
    return 7875;
  else if( aa == "VFT")  
    return 7876;
  else if( aa == "VFW")  
    return 7877;
  else if( aa == "VFY")  
    return 7878;
  else if( aa == "VFV")  
    return 7879;
  else if( aa == "VPA")  
    return 7880;
  else if( aa == "VPR")  
    return 7881;
  else if( aa == "VPN")  
    return 7882;
  else if( aa == "VPD")  
    return 7883;
  else if( aa == "VPC")  
    return 7884;
  else if( aa == "VPQ")  
    return 7885;
  else if( aa == "VPE")  
    return 7886;
  else if( aa == "VPG")  
    return 7887;
  else if( aa == "VPH")  
    return 7888;
  else if( aa == "VPI")  
    return 7889;
  else if( aa == "VPL")  
    return 7890;
  else if( aa == "VPK")  
    return 7891;
  else if( aa == "VPM")  
    return 7892;
  else if( aa == "VPF")  
    return 7893;
  else if( aa == "VPP")  
    return 7894;
  else if( aa == "VPS")  
    return 7895;
  else if( aa == "VPT")  
    return 7896;
  else if( aa == "VPW")  
    return 7897;
  else if( aa == "VPY")  
    return 7898;
  else if( aa == "VPV")  
    return 7899;
  else if( aa == "VSA")  
    return 7900;
  else if( aa == "VSR")  
    return 7901;
  else if( aa == "VSN")  
    return 7902;
  else if( aa == "VSD")  
    return 7903;
  else if( aa == "VSC")  
    return 7904;
  else if( aa == "VSQ")  
    return 7905;
  else if( aa == "VSE")  
    return 7906;
  else if( aa == "VSG")  
    return 7907;
  else if( aa == "VSH")  
    return 7908;
  else if( aa == "VSI")  
    return 7909;
  else if( aa == "VSL")  
    return 7910;
  else if( aa == "VSK")  
    return 7911;
  else if( aa == "VSM")  
    return 7912;
  else if( aa == "VSF")  
    return 7913;
  else if( aa == "VSP")  
    return 7914;
  else if( aa == "VSS")  
    return 7915;
  else if( aa == "VST")  
    return 7916;
  else if( aa == "VSW")  
    return 7917;
  else if( aa == "VSY")  
    return 7918;
  else if( aa == "VSV")  
    return 7919;
  else if( aa == "VTA")  
    return 7920;
  else if( aa == "VTR")  
    return 7921;
  else if( aa == "VTN")  
    return 7922;
  else if( aa == "VTD")  
    return 7923;
  else if( aa == "VTC")  
    return 7924;
  else if( aa == "VTQ")  
    return 7925;
  else if( aa == "VTE")  
    return 7926;
  else if( aa == "VTG")  
    return 7927;
  else if( aa == "VTH")  
    return 7928;
  else if( aa == "VTI")  
    return 7929;
  else if( aa == "VTL")  
    return 7930;
  else if( aa == "VTK")  
    return 7931;
  else if( aa == "VTM")  
    return 7932;
  else if( aa == "VTF")  
    return 7933;
  else if( aa == "VTP")  
    return 7934;
  else if( aa == "VTS")  
    return 7935;
  else if( aa == "VTT")  
    return 7936;
  else if( aa == "VTW")  
    return 7937;
  else if( aa == "VTY")  
    return 7938;
  else if( aa == "VTV")  
    return 7939;
  else if( aa == "VWA")  
    return 7940;
  else if( aa == "VWR")  
    return 7941;
  else if( aa == "VWN")  
    return 7942;
  else if( aa == "VWD")  
    return 7943;
  else if( aa == "VWC")  
    return 7944;
  else if( aa == "VWQ")  
    return 7945;
  else if( aa == "VWE")  
    return 7946;
  else if( aa == "VWG")  
    return 7947;
  else if( aa == "VWH")  
    return 7948;
  else if( aa == "VWI")  
    return 7949;
  else if( aa == "VWL")  
    return 7950;
  else if( aa == "VWK")  
    return 7951;
  else if( aa == "VWM")  
    return 7952;
  else if( aa == "VWF")  
    return 7953;
  else if( aa == "VWP")  
    return 7954;
  else if( aa == "VWS")  
    return 7955;
  else if( aa == "VWT")  
    return 7956;
  else if( aa == "VWW")  
    return 7957;
  else if( aa == "VWY")  
    return 7958;
  else if( aa == "VWV")  
    return 7959;
  else if( aa == "VYA")  
    return 7960;
  else if( aa == "VYR")  
    return 7961;
  else if( aa == "VYN")  
    return 7962;
  else if( aa == "VYD")  
    return 7963;
  else if( aa == "VYC")  
    return 7964;
  else if( aa == "VYQ")  
    return 7965;
  else if( aa == "VYE")  
    return 7966;
  else if( aa == "VYG")  
    return 7967;
  else if( aa == "VYH")  
    return 7968;
  else if( aa == "VYI")  
    return 7969;
  else if( aa == "VYL")  
    return 7970;
  else if( aa == "VYK")  
    return 7971;
  else if( aa == "VYM")  
    return 7972;
  else if( aa == "VYF")  
    return 7973;
  else if( aa == "VYP")  
    return 7974;
  else if( aa == "VYS")  
    return 7975;
  else if( aa == "VYT")  
    return 7976;
  else if( aa == "VYW")  
    return 7977;
  else if( aa == "VYY")  
    return 7978;
  else if( aa == "VYV")  
    return 7979;
  else if( aa == "VVA")  
    return 7980;
  else if( aa == "VVR")  
    return 7981;
  else if( aa == "VVN")  
    return 7982;
  else if( aa == "VVD")  
    return 7983;
  else if( aa == "VVC")  
    return 7984;
  else if( aa == "VVQ")  
    return 7985;
  else if( aa == "VVE")  
    return 7986;
  else if( aa == "VVG")  
    return 7987;
  else if( aa == "VVH")  
    return 7988;
  else if( aa == "VVI")  
    return 7989;
  else if( aa == "VVL")  
    return 7990;
  else if( aa == "VVK")  
    return 7991;
  else if( aa == "VVM")  
    return 7992;
  else if( aa == "VVF")  
    return 7993;
  else if( aa == "VVP")  
    return 7994;
  else if( aa == "VVS")  
    return 7995;
  else if( aa == "VVT")  
    return 7996;
  else if( aa == "VVW")  
    return 7997;
  else if( aa == "VVY")  
    return 7998;
  else if( aa == "VVV")  
    return 7999;
  return 0;
} 

string num_to_aa(int num) { 
  if ( num == 0) 
    return "AAA";
  else if( num == 1) 
    return "AAR";
  else if( num == 2) 
    return "AAN";
  else if( num == 3) 
    return "AAD";
  else if( num == 4) 
    return "AAC";
  else if( num == 5) 
    return "AAQ";
  else if( num == 6) 
    return "AAE";
  else if( num == 7) 
    return "AAG";
  else if( num == 8) 
    return "AAH";
  else if( num == 9) 
    return "AAI";
  else if( num == 10) 
    return "AAL";
  else if( num == 11) 
    return "AAK";
  else if( num == 12) 
    return "AAM";
  else if( num == 13) 
    return "AAF";
  else if( num == 14) 
    return "AAP";
  else if( num == 15) 
    return "AAS";
  else if( num == 16) 
    return "AAT";
  else if( num == 17) 
    return "AAW";
  else if( num == 18) 
    return "AAY";
  else if( num == 19) 
    return "AAV";
  else if( num == 20) 
    return "ARA";
  else if( num == 21) 
    return "ARR";
  else if( num == 22) 
    return "ARN";
  else if( num == 23) 
    return "ARD";
  else if( num == 24) 
    return "ARC";
  else if( num == 25) 
    return "ARQ";
  else if( num == 26) 
    return "ARE";
  else if( num == 27) 
    return "ARG";
  else if( num == 28) 
    return "ARH";
  else if( num == 29) 
    return "ARI";
  else if( num == 30) 
    return "ARL";
  else if( num == 31) 
    return "ARK";
  else if( num == 32) 
    return "ARM";
  else if( num == 33) 
    return "ARF";
  else if( num == 34) 
    return "ARP";
  else if( num == 35) 
    return "ARS";
  else if( num == 36) 
    return "ART";
  else if( num == 37) 
    return "ARW";
  else if( num == 38) 
    return "ARY";
  else if( num == 39) 
    return "ARV";
  else if( num == 40) 
    return "ANA";
  else if( num == 41) 
    return "ANR";
  else if( num == 42) 
    return "ANN";
  else if( num == 43) 
    return "AND";
  else if( num == 44) 
    return "ANC";
  else if( num == 45) 
    return "ANQ";
  else if( num == 46) 
    return "ANE";
  else if( num == 47) 
    return "ANG";
  else if( num == 48) 
    return "ANH";
  else if( num == 49) 
    return "ANI";
  else if( num == 50) 
    return "ANL";
  else if( num == 51) 
    return "ANK";
  else if( num == 52) 
    return "ANM";
  else if( num == 53) 
    return "ANF";
  else if( num == 54) 
    return "ANP";
  else if( num == 55) 
    return "ANS";
  else if( num == 56) 
    return "ANT";
  else if( num == 57) 
    return "ANW";
  else if( num == 58) 
    return "ANY";
  else if( num == 59) 
    return "ANV";
  else if( num == 60) 
    return "ADA";
  else if( num == 61) 
    return "ADR";
  else if( num == 62) 
    return "ADN";
  else if( num == 63) 
    return "ADD";
  else if( num == 64) 
    return "ADC";
  else if( num == 65) 
    return "ADQ";
  else if( num == 66) 
    return "ADE";
  else if( num == 67) 
    return "ADG";
  else if( num == 68) 
    return "ADH";
  else if( num == 69) 
    return "ADI";
  else if( num == 70) 
    return "ADL";
  else if( num == 71) 
    return "ADK";
  else if( num == 72) 
    return "ADM";
  else if( num == 73) 
    return "ADF";
  else if( num == 74) 
    return "ADP";
  else if( num == 75) 
    return "ADS";
  else if( num == 76) 
    return "ADT";
  else if( num == 77) 
    return "ADW";
  else if( num == 78) 
    return "ADY";
  else if( num == 79) 
    return "ADV";
  else if( num == 80) 
    return "ACA";
  else if( num == 81) 
    return "ACR";
  else if( num == 82) 
    return "ACN";
  else if( num == 83) 
    return "ACD";
  else if( num == 84) 
    return "ACC";
  else if( num == 85) 
    return "ACQ";
  else if( num == 86) 
    return "ACE";
  else if( num == 87) 
    return "ACG";
  else if( num == 88) 
    return "ACH";
  else if( num == 89) 
    return "ACI";
  else if( num == 90) 
    return "ACL";
  else if( num == 91) 
    return "ACK";
  else if( num == 92) 
    return "ACM";
  else if( num == 93) 
    return "ACF";
  else if( num == 94) 
    return "ACP";
  else if( num == 95) 
    return "ACS";
  else if( num == 96) 
    return "ACT";
  else if( num == 97) 
    return "ACW";
  else if( num == 98) 
    return "ACY";
  else if( num == 99) 
    return "ACV";
  else if( num == 100) 
    return "AQA";
  else if( num == 101) 
    return "AQR";
  else if( num == 102) 
    return "AQN";
  else if( num == 103) 
    return "AQD";
  else if( num == 104) 
    return "AQC";
  else if( num == 105) 
    return "AQQ";
  else if( num == 106) 
    return "AQE";
  else if( num == 107) 
    return "AQG";
  else if( num == 108) 
    return "AQH";
  else if( num == 109) 
    return "AQI";
  else if( num == 110) 
    return "AQL";
  else if( num == 111) 
    return "AQK";
  else if( num == 112) 
    return "AQM";
  else if( num == 113) 
    return "AQF";
  else if( num == 114) 
    return "AQP";
  else if( num == 115) 
    return "AQS";
  else if( num == 116) 
    return "AQT";
  else if( num == 117) 
    return "AQW";
  else if( num == 118) 
    return "AQY";
  else if( num == 119) 
    return "AQV";
  else if( num == 120) 
    return "AEA";
  else if( num == 121) 
    return "AER";
  else if( num == 122) 
    return "AEN";
  else if( num == 123) 
    return "AED";
  else if( num == 124) 
    return "AEC";
  else if( num == 125) 
    return "AEQ";
  else if( num == 126) 
    return "AEE";
  else if( num == 127) 
    return "AEG";
  else if( num == 128) 
    return "AEH";
  else if( num == 129) 
    return "AEI";
  else if( num == 130) 
    return "AEL";
  else if( num == 131) 
    return "AEK";
  else if( num == 132) 
    return "AEM";
  else if( num == 133) 
    return "AEF";
  else if( num == 134) 
    return "AEP";
  else if( num == 135) 
    return "AES";
  else if( num == 136) 
    return "AET";
  else if( num == 137) 
    return "AEW";
  else if( num == 138) 
    return "AEY";
  else if( num == 139) 
    return "AEV";
  else if( num == 140) 
    return "AGA";
  else if( num == 141) 
    return "AGR";
  else if( num == 142) 
    return "AGN";
  else if( num == 143) 
    return "AGD";
  else if( num == 144) 
    return "AGC";
  else if( num == 145) 
    return "AGQ";
  else if( num == 146) 
    return "AGE";
  else if( num == 147) 
    return "AGG";
  else if( num == 148) 
    return "AGH";
  else if( num == 149) 
    return "AGI";
  else if( num == 150) 
    return "AGL";
  else if( num == 151) 
    return "AGK";
  else if( num == 152) 
    return "AGM";
  else if( num == 153) 
    return "AGF";
  else if( num == 154) 
    return "AGP";
  else if( num == 155) 
    return "AGS";
  else if( num == 156) 
    return "AGT";
  else if( num == 157) 
    return "AGW";
  else if( num == 158) 
    return "AGY";
  else if( num == 159) 
    return "AGV";
  else if( num == 160) 
    return "AHA";
  else if( num == 161) 
    return "AHR";
  else if( num == 162) 
    return "AHN";
  else if( num == 163) 
    return "AHD";
  else if( num == 164) 
    return "AHC";
  else if( num == 165) 
    return "AHQ";
  else if( num == 166) 
    return "AHE";
  else if( num == 167) 
    return "AHG";
  else if( num == 168) 
    return "AHH";
  else if( num == 169) 
    return "AHI";
  else if( num == 170) 
    return "AHL";
  else if( num == 171) 
    return "AHK";
  else if( num == 172) 
    return "AHM";
  else if( num == 173) 
    return "AHF";
  else if( num == 174) 
    return "AHP";
  else if( num == 175) 
    return "AHS";
  else if( num == 176) 
    return "AHT";
  else if( num == 177) 
    return "AHW";
  else if( num == 178) 
    return "AHY";
  else if( num == 179) 
    return "AHV";
  else if( num == 180) 
    return "AIA";
  else if( num == 181) 
    return "AIR";
  else if( num == 182) 
    return "AIN";
  else if( num == 183) 
    return "AID";
  else if( num == 184) 
    return "AIC";
  else if( num == 185) 
    return "AIQ";
  else if( num == 186) 
    return "AIE";
  else if( num == 187) 
    return "AIG";
  else if( num == 188) 
    return "AIH";
  else if( num == 189) 
    return "AII";
  else if( num == 190) 
    return "AIL";
  else if( num == 191) 
    return "AIK";
  else if( num == 192) 
    return "AIM";
  else if( num == 193) 
    return "AIF";
  else if( num == 194) 
    return "AIP";
  else if( num == 195) 
    return "AIS";
  else if( num == 196) 
    return "AIT";
  else if( num == 197) 
    return "AIW";
  else if( num == 198) 
    return "AIY";
  else if( num == 199) 
    return "AIV";
  else if( num == 200) 
    return "ALA";
  else if( num == 201) 
    return "ALR";
  else if( num == 202) 
    return "ALN";
  else if( num == 203) 
    return "ALD";
  else if( num == 204) 
    return "ALC";
  else if( num == 205) 
    return "ALQ";
  else if( num == 206) 
    return "ALE";
  else if( num == 207) 
    return "ALG";
  else if( num == 208) 
    return "ALH";
  else if( num == 209) 
    return "ALI";
  else if( num == 210) 
    return "ALL";
  else if( num == 211) 
    return "ALK";
  else if( num == 212) 
    return "ALM";
  else if( num == 213) 
    return "ALF";
  else if( num == 214) 
    return "ALP";
  else if( num == 215) 
    return "ALS";
  else if( num == 216) 
    return "ALT";
  else if( num == 217) 
    return "ALW";
  else if( num == 218) 
    return "ALY";
  else if( num == 219) 
    return "ALV";
  else if( num == 220) 
    return "AKA";
  else if( num == 221) 
    return "AKR";
  else if( num == 222) 
    return "AKN";
  else if( num == 223) 
    return "AKD";
  else if( num == 224) 
    return "AKC";
  else if( num == 225) 
    return "AKQ";
  else if( num == 226) 
    return "AKE";
  else if( num == 227) 
    return "AKG";
  else if( num == 228) 
    return "AKH";
  else if( num == 229) 
    return "AKI";
  else if( num == 230) 
    return "AKL";
  else if( num == 231) 
    return "AKK";
  else if( num == 232) 
    return "AKM";
  else if( num == 233) 
    return "AKF";
  else if( num == 234) 
    return "AKP";
  else if( num == 235) 
    return "AKS";
  else if( num == 236) 
    return "AKT";
  else if( num == 237) 
    return "AKW";
  else if( num == 238) 
    return "AKY";
  else if( num == 239) 
    return "AKV";
  else if( num == 240) 
    return "AMA";
  else if( num == 241) 
    return "AMR";
  else if( num == 242) 
    return "AMN";
  else if( num == 243) 
    return "AMD";
  else if( num == 244) 
    return "AMC";
  else if( num == 245) 
    return "AMQ";
  else if( num == 246) 
    return "AME";
  else if( num == 247) 
    return "AMG";
  else if( num == 248) 
    return "AMH";
  else if( num == 249) 
    return "AMI";
  else if( num == 250) 
    return "AML";
  else if( num == 251) 
    return "AMK";
  else if( num == 252) 
    return "AMM";
  else if( num == 253) 
    return "AMF";
  else if( num == 254) 
    return "AMP";
  else if( num == 255) 
    return "AMS";
  else if( num == 256) 
    return "AMT";
  else if( num == 257) 
    return "AMW";
  else if( num == 258) 
    return "AMY";
  else if( num == 259) 
    return "AMV";
  else if( num == 260) 
    return "AFA";
  else if( num == 261) 
    return "AFR";
  else if( num == 262) 
    return "AFN";
  else if( num == 263) 
    return "AFD";
  else if( num == 264) 
    return "AFC";
  else if( num == 265) 
    return "AFQ";
  else if( num == 266) 
    return "AFE";
  else if( num == 267) 
    return "AFG";
  else if( num == 268) 
    return "AFH";
  else if( num == 269) 
    return "AFI";
  else if( num == 270) 
    return "AFL";
  else if( num == 271) 
    return "AFK";
  else if( num == 272) 
    return "AFM";
  else if( num == 273) 
    return "AFF";
  else if( num == 274) 
    return "AFP";
  else if( num == 275) 
    return "AFS";
  else if( num == 276) 
    return "AFT";
  else if( num == 277) 
    return "AFW";
  else if( num == 278) 
    return "AFY";
  else if( num == 279) 
    return "AFV";
  else if( num == 280) 
    return "APA";
  else if( num == 281) 
    return "APR";
  else if( num == 282) 
    return "APN";
  else if( num == 283) 
    return "APD";
  else if( num == 284) 
    return "APC";
  else if( num == 285) 
    return "APQ";
  else if( num == 286) 
    return "APE";
  else if( num == 287) 
    return "APG";
  else if( num == 288) 
    return "APH";
  else if( num == 289) 
    return "API";
  else if( num == 290) 
    return "APL";
  else if( num == 291) 
    return "APK";
  else if( num == 292) 
    return "APM";
  else if( num == 293) 
    return "APF";
  else if( num == 294) 
    return "APP";
  else if( num == 295) 
    return "APS";
  else if( num == 296) 
    return "APT";
  else if( num == 297) 
    return "APW";
  else if( num == 298) 
    return "APY";
  else if( num == 299) 
    return "APV";
  else if( num == 300) 
    return "ASA";
  else if( num == 301) 
    return "ASR";
  else if( num == 302) 
    return "ASN";
  else if( num == 303) 
    return "ASD";
  else if( num == 304) 
    return "ASC";
  else if( num == 305) 
    return "ASQ";
  else if( num == 306) 
    return "ASE";
  else if( num == 307) 
    return "ASG";
  else if( num == 308) 
    return "ASH";
  else if( num == 309) 
    return "ASI";
  else if( num == 310) 
    return "ASL";
  else if( num == 311) 
    return "ASK";
  else if( num == 312) 
    return "ASM";
  else if( num == 313) 
    return "ASF";
  else if( num == 314) 
    return "ASP";
  else if( num == 315) 
    return "ASS";
  else if( num == 316) 
    return "AST";
  else if( num == 317) 
    return "ASW";
  else if( num == 318) 
    return "ASY";
  else if( num == 319) 
    return "ASV";
  else if( num == 320) 
    return "ATA";
  else if( num == 321) 
    return "ATR";
  else if( num == 322) 
    return "ATN";
  else if( num == 323) 
    return "ATD";
  else if( num == 324) 
    return "ATC";
  else if( num == 325) 
    return "ATQ";
  else if( num == 326) 
    return "ATE";
  else if( num == 327) 
    return "ATG";
  else if( num == 328) 
    return "ATH";
  else if( num == 329) 
    return "ATI";
  else if( num == 330) 
    return "ATL";
  else if( num == 331) 
    return "ATK";
  else if( num == 332) 
    return "ATM";
  else if( num == 333) 
    return "ATF";
  else if( num == 334) 
    return "ATP";
  else if( num == 335) 
    return "ATS";
  else if( num == 336) 
    return "ATT";
  else if( num == 337) 
    return "ATW";
  else if( num == 338) 
    return "ATY";
  else if( num == 339) 
    return "ATV";
  else if( num == 340) 
    return "AWA";
  else if( num == 341) 
    return "AWR";
  else if( num == 342) 
    return "AWN";
  else if( num == 343) 
    return "AWD";
  else if( num == 344) 
    return "AWC";
  else if( num == 345) 
    return "AWQ";
  else if( num == 346) 
    return "AWE";
  else if( num == 347) 
    return "AWG";
  else if( num == 348) 
    return "AWH";
  else if( num == 349) 
    return "AWI";
  else if( num == 350) 
    return "AWL";
  else if( num == 351) 
    return "AWK";
  else if( num == 352) 
    return "AWM";
  else if( num == 353) 
    return "AWF";
  else if( num == 354) 
    return "AWP";
  else if( num == 355) 
    return "AWS";
  else if( num == 356) 
    return "AWT";
  else if( num == 357) 
    return "AWW";
  else if( num == 358) 
    return "AWY";
  else if( num == 359) 
    return "AWV";
  else if( num == 360) 
    return "AYA";
  else if( num == 361) 
    return "AYR";
  else if( num == 362) 
    return "AYN";
  else if( num == 363) 
    return "AYD";
  else if( num == 364) 
    return "AYC";
  else if( num == 365) 
    return "AYQ";
  else if( num == 366) 
    return "AYE";
  else if( num == 367) 
    return "AYG";
  else if( num == 368) 
    return "AYH";
  else if( num == 369) 
    return "AYI";
  else if( num == 370) 
    return "AYL";
  else if( num == 371) 
    return "AYK";
  else if( num == 372) 
    return "AYM";
  else if( num == 373) 
    return "AYF";
  else if( num == 374) 
    return "AYP";
  else if( num == 375) 
    return "AYS";
  else if( num == 376) 
    return "AYT";
  else if( num == 377) 
    return "AYW";
  else if( num == 378) 
    return "AYY";
  else if( num == 379) 
    return "AYV";
  else if( num == 380) 
    return "AVA";
  else if( num == 381) 
    return "AVR";
  else if( num == 382) 
    return "AVN";
  else if( num == 383) 
    return "AVD";
  else if( num == 384) 
    return "AVC";
  else if( num == 385) 
    return "AVQ";
  else if( num == 386) 
    return "AVE";
  else if( num == 387) 
    return "AVG";
  else if( num == 388) 
    return "AVH";
  else if( num == 389) 
    return "AVI";
  else if( num == 390) 
    return "AVL";
  else if( num == 391) 
    return "AVK";
  else if( num == 392) 
    return "AVM";
  else if( num == 393) 
    return "AVF";
  else if( num == 394) 
    return "AVP";
  else if( num == 395) 
    return "AVS";
  else if( num == 396) 
    return "AVT";
  else if( num == 397) 
    return "AVW";
  else if( num == 398) 
    return "AVY";
  else if( num == 399) 
    return "AVV";
  else if( num == 400) 
    return "RAA";
  else if( num == 401) 
    return "RAR";
  else if( num == 402) 
    return "RAN";
  else if( num == 403) 
    return "RAD";
  else if( num == 404) 
    return "RAC";
  else if( num == 405) 
    return "RAQ";
  else if( num == 406) 
    return "RAE";
  else if( num == 407) 
    return "RAG";
  else if( num == 408) 
    return "RAH";
  else if( num == 409) 
    return "RAI";
  else if( num == 410) 
    return "RAL";
  else if( num == 411) 
    return "RAK";
  else if( num == 412) 
    return "RAM";
  else if( num == 413) 
    return "RAF";
  else if( num == 414) 
    return "RAP";
  else if( num == 415) 
    return "RAS";
  else if( num == 416) 
    return "RAT";
  else if( num == 417) 
    return "RAW";
  else if( num == 418) 
    return "RAY";
  else if( num == 419) 
    return "RAV";
  else if( num == 420) 
    return "RRA";
  else if( num == 421) 
    return "RRR";
  else if( num == 422) 
    return "RRN";
  else if( num == 423) 
    return "RRD";
  else if( num == 424) 
    return "RRC";
  else if( num == 425) 
    return "RRQ";
  else if( num == 426) 
    return "RRE";
  else if( num == 427) 
    return "RRG";
  else if( num == 428) 
    return "RRH";
  else if( num == 429) 
    return "RRI";
  else if( num == 430) 
    return "RRL";
  else if( num == 431) 
    return "RRK";
  else if( num == 432) 
    return "RRM";
  else if( num == 433) 
    return "RRF";
  else if( num == 434) 
    return "RRP";
  else if( num == 435) 
    return "RRS";
  else if( num == 436) 
    return "RRT";
  else if( num == 437) 
    return "RRW";
  else if( num == 438) 
    return "RRY";
  else if( num == 439) 
    return "RRV";
  else if( num == 440) 
    return "RNA";
  else if( num == 441) 
    return "RNR";
  else if( num == 442) 
    return "RNN";
  else if( num == 443) 
    return "RND";
  else if( num == 444) 
    return "RNC";
  else if( num == 445) 
    return "RNQ";
  else if( num == 446) 
    return "RNE";
  else if( num == 447) 
    return "RNG";
  else if( num == 448) 
    return "RNH";
  else if( num == 449) 
    return "RNI";
  else if( num == 450) 
    return "RNL";
  else if( num == 451) 
    return "RNK";
  else if( num == 452) 
    return "RNM";
  else if( num == 453) 
    return "RNF";
  else if( num == 454) 
    return "RNP";
  else if( num == 455) 
    return "RNS";
  else if( num == 456) 
    return "RNT";
  else if( num == 457) 
    return "RNW";
  else if( num == 458) 
    return "RNY";
  else if( num == 459) 
    return "RNV";
  else if( num == 460) 
    return "RDA";
  else if( num == 461) 
    return "RDR";
  else if( num == 462) 
    return "RDN";
  else if( num == 463) 
    return "RDD";
  else if( num == 464) 
    return "RDC";
  else if( num == 465) 
    return "RDQ";
  else if( num == 466) 
    return "RDE";
  else if( num == 467) 
    return "RDG";
  else if( num == 468) 
    return "RDH";
  else if( num == 469) 
    return "RDI";
  else if( num == 470) 
    return "RDL";
  else if( num == 471) 
    return "RDK";
  else if( num == 472) 
    return "RDM";
  else if( num == 473) 
    return "RDF";
  else if( num == 474) 
    return "RDP";
  else if( num == 475) 
    return "RDS";
  else if( num == 476) 
    return "RDT";
  else if( num == 477) 
    return "RDW";
  else if( num == 478) 
    return "RDY";
  else if( num == 479) 
    return "RDV";
  else if( num == 480) 
    return "RCA";
  else if( num == 481) 
    return "RCR";
  else if( num == 482) 
    return "RCN";
  else if( num == 483) 
    return "RCD";
  else if( num == 484) 
    return "RCC";
  else if( num == 485) 
    return "RCQ";
  else if( num == 486) 
    return "RCE";
  else if( num == 487) 
    return "RCG";
  else if( num == 488) 
    return "RCH";
  else if( num == 489) 
    return "RCI";
  else if( num == 490) 
    return "RCL";
  else if( num == 491) 
    return "RCK";
  else if( num == 492) 
    return "RCM";
  else if( num == 493) 
    return "RCF";
  else if( num == 494) 
    return "RCP";
  else if( num == 495) 
    return "RCS";
  else if( num == 496) 
    return "RCT";
  else if( num == 497) 
    return "RCW";
  else if( num == 498) 
    return "RCY";
  else if( num == 499) 
    return "RCV";
  else if( num == 500) 
    return "RQA";
  else if( num == 501) 
    return "RQR";
  else if( num == 502) 
    return "RQN";
  else if( num == 503) 
    return "RQD";
  else if( num == 504) 
    return "RQC";
  else if( num == 505) 
    return "RQQ";
  else if( num == 506) 
    return "RQE";
  else if( num == 507) 
    return "RQG";
  else if( num == 508) 
    return "RQH";
  else if( num == 509) 
    return "RQI";
  else if( num == 510) 
    return "RQL";
  else if( num == 511) 
    return "RQK";
  else if( num == 512) 
    return "RQM";
  else if( num == 513) 
    return "RQF";
  else if( num == 514) 
    return "RQP";
  else if( num == 515) 
    return "RQS";
  else if( num == 516) 
    return "RQT";
  else if( num == 517) 
    return "RQW";
  else if( num == 518) 
    return "RQY";
  else if( num == 519) 
    return "RQV";
  else if( num == 520) 
    return "REA";
  else if( num == 521) 
    return "RER";
  else if( num == 522) 
    return "REN";
  else if( num == 523) 
    return "RED";
  else if( num == 524) 
    return "REC";
  else if( num == 525) 
    return "REQ";
  else if( num == 526) 
    return "REE";
  else if( num == 527) 
    return "REG";
  else if( num == 528) 
    return "REH";
  else if( num == 529) 
    return "REI";
  else if( num == 530) 
    return "REL";
  else if( num == 531) 
    return "REK";
  else if( num == 532) 
    return "REM";
  else if( num == 533) 
    return "REF";
  else if( num == 534) 
    return "REP";
  else if( num == 535) 
    return "RES";
  else if( num == 536) 
    return "RET";
  else if( num == 537) 
    return "REW";
  else if( num == 538) 
    return "REY";
  else if( num == 539) 
    return "REV";
  else if( num == 540) 
    return "RGA";
  else if( num == 541) 
    return "RGR";
  else if( num == 542) 
    return "RGN";
  else if( num == 543) 
    return "RGD";
  else if( num == 544) 
    return "RGC";
  else if( num == 545) 
    return "RGQ";
  else if( num == 546) 
    return "RGE";
  else if( num == 547) 
    return "RGG";
  else if( num == 548) 
    return "RGH";
  else if( num == 549) 
    return "RGI";
  else if( num == 550) 
    return "RGL";
  else if( num == 551) 
    return "RGK";
  else if( num == 552) 
    return "RGM";
  else if( num == 553) 
    return "RGF";
  else if( num == 554) 
    return "RGP";
  else if( num == 555) 
    return "RGS";
  else if( num == 556) 
    return "RGT";
  else if( num == 557) 
    return "RGW";
  else if( num == 558) 
    return "RGY";
  else if( num == 559) 
    return "RGV";
  else if( num == 560) 
    return "RHA";
  else if( num == 561) 
    return "RHR";
  else if( num == 562) 
    return "RHN";
  else if( num == 563) 
    return "RHD";
  else if( num == 564) 
    return "RHC";
  else if( num == 565) 
    return "RHQ";
  else if( num == 566) 
    return "RHE";
  else if( num == 567) 
    return "RHG";
  else if( num == 568) 
    return "RHH";
  else if( num == 569) 
    return "RHI";
  else if( num == 570) 
    return "RHL";
  else if( num == 571) 
    return "RHK";
  else if( num == 572) 
    return "RHM";
  else if( num == 573) 
    return "RHF";
  else if( num == 574) 
    return "RHP";
  else if( num == 575) 
    return "RHS";
  else if( num == 576) 
    return "RHT";
  else if( num == 577) 
    return "RHW";
  else if( num == 578) 
    return "RHY";
  else if( num == 579) 
    return "RHV";
  else if( num == 580) 
    return "RIA";
  else if( num == 581) 
    return "RIR";
  else if( num == 582) 
    return "RIN";
  else if( num == 583) 
    return "RID";
  else if( num == 584) 
    return "RIC";
  else if( num == 585) 
    return "RIQ";
  else if( num == 586) 
    return "RIE";
  else if( num == 587) 
    return "RIG";
  else if( num == 588) 
    return "RIH";
  else if( num == 589) 
    return "RII";
  else if( num == 590) 
    return "RIL";
  else if( num == 591) 
    return "RIK";
  else if( num == 592) 
    return "RIM";
  else if( num == 593) 
    return "RIF";
  else if( num == 594) 
    return "RIP";
  else if( num == 595) 
    return "RIS";
  else if( num == 596) 
    return "RIT";
  else if( num == 597) 
    return "RIW";
  else if( num == 598) 
    return "RIY";
  else if( num == 599) 
    return "RIV";
  else if( num == 600) 
    return "RLA";
  else if( num == 601) 
    return "RLR";
  else if( num == 602) 
    return "RLN";
  else if( num == 603) 
    return "RLD";
  else if( num == 604) 
    return "RLC";
  else if( num == 605) 
    return "RLQ";
  else if( num == 606) 
    return "RLE";
  else if( num == 607) 
    return "RLG";
  else if( num == 608) 
    return "RLH";
  else if( num == 609) 
    return "RLI";
  else if( num == 610) 
    return "RLL";
  else if( num == 611) 
    return "RLK";
  else if( num == 612) 
    return "RLM";
  else if( num == 613) 
    return "RLF";
  else if( num == 614) 
    return "RLP";
  else if( num == 615) 
    return "RLS";
  else if( num == 616) 
    return "RLT";
  else if( num == 617) 
    return "RLW";
  else if( num == 618) 
    return "RLY";
  else if( num == 619) 
    return "RLV";
  else if( num == 620) 
    return "RKA";
  else if( num == 621) 
    return "RKR";
  else if( num == 622) 
    return "RKN";
  else if( num == 623) 
    return "RKD";
  else if( num == 624) 
    return "RKC";
  else if( num == 625) 
    return "RKQ";
  else if( num == 626) 
    return "RKE";
  else if( num == 627) 
    return "RKG";
  else if( num == 628) 
    return "RKH";
  else if( num == 629) 
    return "RKI";
  else if( num == 630) 
    return "RKL";
  else if( num == 631) 
    return "RKK";
  else if( num == 632) 
    return "RKM";
  else if( num == 633) 
    return "RKF";
  else if( num == 634) 
    return "RKP";
  else if( num == 635) 
    return "RKS";
  else if( num == 636) 
    return "RKT";
  else if( num == 637) 
    return "RKW";
  else if( num == 638) 
    return "RKY";
  else if( num == 639) 
    return "RKV";
  else if( num == 640) 
    return "RMA";
  else if( num == 641) 
    return "RMR";
  else if( num == 642) 
    return "RMN";
  else if( num == 643) 
    return "RMD";
  else if( num == 644) 
    return "RMC";
  else if( num == 645) 
    return "RMQ";
  else if( num == 646) 
    return "RME";
  else if( num == 647) 
    return "RMG";
  else if( num == 648) 
    return "RMH";
  else if( num == 649) 
    return "RMI";
  else if( num == 650) 
    return "RML";
  else if( num == 651) 
    return "RMK";
  else if( num == 652) 
    return "RMM";
  else if( num == 653) 
    return "RMF";
  else if( num == 654) 
    return "RMP";
  else if( num == 655) 
    return "RMS";
  else if( num == 656) 
    return "RMT";
  else if( num == 657) 
    return "RMW";
  else if( num == 658) 
    return "RMY";
  else if( num == 659) 
    return "RMV";
  else if( num == 660) 
    return "RFA";
  else if( num == 661) 
    return "RFR";
  else if( num == 662) 
    return "RFN";
  else if( num == 663) 
    return "RFD";
  else if( num == 664) 
    return "RFC";
  else if( num == 665) 
    return "RFQ";
  else if( num == 666) 
    return "RFE";
  else if( num == 667) 
    return "RFG";
  else if( num == 668) 
    return "RFH";
  else if( num == 669) 
    return "RFI";
  else if( num == 670) 
    return "RFL";
  else if( num == 671) 
    return "RFK";
  else if( num == 672) 
    return "RFM";
  else if( num == 673) 
    return "RFF";
  else if( num == 674) 
    return "RFP";
  else if( num == 675) 
    return "RFS";
  else if( num == 676) 
    return "RFT";
  else if( num == 677) 
    return "RFW";
  else if( num == 678) 
    return "RFY";
  else if( num == 679) 
    return "RFV";
  else if( num == 680) 
    return "RPA";
  else if( num == 681) 
    return "RPR";
  else if( num == 682) 
    return "RPN";
  else if( num == 683) 
    return "RPD";
  else if( num == 684) 
    return "RPC";
  else if( num == 685) 
    return "RPQ";
  else if( num == 686) 
    return "RPE";
  else if( num == 687) 
    return "RPG";
  else if( num == 688) 
    return "RPH";
  else if( num == 689) 
    return "RPI";
  else if( num == 690) 
    return "RPL";
  else if( num == 691) 
    return "RPK";
  else if( num == 692) 
    return "RPM";
  else if( num == 693) 
    return "RPF";
  else if( num == 694) 
    return "RPP";
  else if( num == 695) 
    return "RPS";
  else if( num == 696) 
    return "RPT";
  else if( num == 697) 
    return "RPW";
  else if( num == 698) 
    return "RPY";
  else if( num == 699) 
    return "RPV";
  else if( num == 700) 
    return "RSA";
  else if( num == 701) 
    return "RSR";
  else if( num == 702) 
    return "RSN";
  else if( num == 703) 
    return "RSD";
  else if( num == 704) 
    return "RSC";
  else if( num == 705) 
    return "RSQ";
  else if( num == 706) 
    return "RSE";
  else if( num == 707) 
    return "RSG";
  else if( num == 708) 
    return "RSH";
  else if( num == 709) 
    return "RSI";
  else if( num == 710) 
    return "RSL";
  else if( num == 711) 
    return "RSK";
  else if( num == 712) 
    return "RSM";
  else if( num == 713) 
    return "RSF";
  else if( num == 714) 
    return "RSP";
  else if( num == 715) 
    return "RSS";
  else if( num == 716) 
    return "RST";
  else if( num == 717) 
    return "RSW";
  else if( num == 718) 
    return "RSY";
  else if( num == 719) 
    return "RSV";
  else if( num == 720) 
    return "RTA";
  else if( num == 721) 
    return "RTR";
  else if( num == 722) 
    return "RTN";
  else if( num == 723) 
    return "RTD";
  else if( num == 724) 
    return "RTC";
  else if( num == 725) 
    return "RTQ";
  else if( num == 726) 
    return "RTE";
  else if( num == 727) 
    return "RTG";
  else if( num == 728) 
    return "RTH";
  else if( num == 729) 
    return "RTI";
  else if( num == 730) 
    return "RTL";
  else if( num == 731) 
    return "RTK";
  else if( num == 732) 
    return "RTM";
  else if( num == 733) 
    return "RTF";
  else if( num == 734) 
    return "RTP";
  else if( num == 735) 
    return "RTS";
  else if( num == 736) 
    return "RTT";
  else if( num == 737) 
    return "RTW";
  else if( num == 738) 
    return "RTY";
  else if( num == 739) 
    return "RTV";
  else if( num == 740) 
    return "RWA";
  else if( num == 741) 
    return "RWR";
  else if( num == 742) 
    return "RWN";
  else if( num == 743) 
    return "RWD";
  else if( num == 744) 
    return "RWC";
  else if( num == 745) 
    return "RWQ";
  else if( num == 746) 
    return "RWE";
  else if( num == 747) 
    return "RWG";
  else if( num == 748) 
    return "RWH";
  else if( num == 749) 
    return "RWI";
  else if( num == 750) 
    return "RWL";
  else if( num == 751) 
    return "RWK";
  else if( num == 752) 
    return "RWM";
  else if( num == 753) 
    return "RWF";
  else if( num == 754) 
    return "RWP";
  else if( num == 755) 
    return "RWS";
  else if( num == 756) 
    return "RWT";
  else if( num == 757) 
    return "RWW";
  else if( num == 758) 
    return "RWY";
  else if( num == 759) 
    return "RWV";
  else if( num == 760) 
    return "RYA";
  else if( num == 761) 
    return "RYR";
  else if( num == 762) 
    return "RYN";
  else if( num == 763) 
    return "RYD";
  else if( num == 764) 
    return "RYC";
  else if( num == 765) 
    return "RYQ";
  else if( num == 766) 
    return "RYE";
  else if( num == 767) 
    return "RYG";
  else if( num == 768) 
    return "RYH";
  else if( num == 769) 
    return "RYI";
  else if( num == 770) 
    return "RYL";
  else if( num == 771) 
    return "RYK";
  else if( num == 772) 
    return "RYM";
  else if( num == 773) 
    return "RYF";
  else if( num == 774) 
    return "RYP";
  else if( num == 775) 
    return "RYS";
  else if( num == 776) 
    return "RYT";
  else if( num == 777) 
    return "RYW";
  else if( num == 778) 
    return "RYY";
  else if( num == 779) 
    return "RYV";
  else if( num == 780) 
    return "RVA";
  else if( num == 781) 
    return "RVR";
  else if( num == 782) 
    return "RVN";
  else if( num == 783) 
    return "RVD";
  else if( num == 784) 
    return "RVC";
  else if( num == 785) 
    return "RVQ";
  else if( num == 786) 
    return "RVE";
  else if( num == 787) 
    return "RVG";
  else if( num == 788) 
    return "RVH";
  else if( num == 789) 
    return "RVI";
  else if( num == 790) 
    return "RVL";
  else if( num == 791) 
    return "RVK";
  else if( num == 792) 
    return "RVM";
  else if( num == 793) 
    return "RVF";
  else if( num == 794) 
    return "RVP";
  else if( num == 795) 
    return "RVS";
  else if( num == 796) 
    return "RVT";
  else if( num == 797) 
    return "RVW";
  else if( num == 798) 
    return "RVY";
  else if( num == 799) 
    return "RVV";
  else if( num == 800) 
    return "NAA";
  else if( num == 801) 
    return "NAR";
  else if( num == 802) 
    return "NAN";
  else if( num == 803) 
    return "NAD";
  else if( num == 804) 
    return "NAC";
  else if( num == 805) 
    return "NAQ";
  else if( num == 806) 
    return "NAE";
  else if( num == 807) 
    return "NAG";
  else if( num == 808) 
    return "NAH";
  else if( num == 809) 
    return "NAI";
  else if( num == 810) 
    return "NAL";
  else if( num == 811) 
    return "NAK";
  else if( num == 812) 
    return "NAM";
  else if( num == 813) 
    return "NAF";
  else if( num == 814) 
    return "NAP";
  else if( num == 815) 
    return "NAS";
  else if( num == 816) 
    return "NAT";
  else if( num == 817) 
    return "NAW";
  else if( num == 818) 
    return "NAY";
  else if( num == 819) 
    return "NAV";
  else if( num == 820) 
    return "NRA";
  else if( num == 821) 
    return "NRR";
  else if( num == 822) 
    return "NRN";
  else if( num == 823) 
    return "NRD";
  else if( num == 824) 
    return "NRC";
  else if( num == 825) 
    return "NRQ";
  else if( num == 826) 
    return "NRE";
  else if( num == 827) 
    return "NRG";
  else if( num == 828) 
    return "NRH";
  else if( num == 829) 
    return "NRI";
  else if( num == 830) 
    return "NRL";
  else if( num == 831) 
    return "NRK";
  else if( num == 832) 
    return "NRM";
  else if( num == 833) 
    return "NRF";
  else if( num == 834) 
    return "NRP";
  else if( num == 835) 
    return "NRS";
  else if( num == 836) 
    return "NRT";
  else if( num == 837) 
    return "NRW";
  else if( num == 838) 
    return "NRY";
  else if( num == 839) 
    return "NRV";
  else if( num == 840) 
    return "NNA";
  else if( num == 841) 
    return "NNR";
  else if( num == 842) 
    return "NNN";
  else if( num == 843) 
    return "NND";
  else if( num == 844) 
    return "NNC";
  else if( num == 845) 
    return "NNQ";
  else if( num == 846) 
    return "NNE";
  else if( num == 847) 
    return "NNG";
  else if( num == 848) 
    return "NNH";
  else if( num == 849) 
    return "NNI";
  else if( num == 850) 
    return "NNL";
  else if( num == 851) 
    return "NNK";
  else if( num == 852) 
    return "NNM";
  else if( num == 853) 
    return "NNF";
  else if( num == 854) 
    return "NNP";
  else if( num == 855) 
    return "NNS";
  else if( num == 856) 
    return "NNT";
  else if( num == 857) 
    return "NNW";
  else if( num == 858) 
    return "NNY";
  else if( num == 859) 
    return "NNV";
  else if( num == 860) 
    return "NDA";
  else if( num == 861) 
    return "NDR";
  else if( num == 862) 
    return "NDN";
  else if( num == 863) 
    return "NDD";
  else if( num == 864) 
    return "NDC";
  else if( num == 865) 
    return "NDQ";
  else if( num == 866) 
    return "NDE";
  else if( num == 867) 
    return "NDG";
  else if( num == 868) 
    return "NDH";
  else if( num == 869) 
    return "NDI";
  else if( num == 870) 
    return "NDL";
  else if( num == 871) 
    return "NDK";
  else if( num == 872) 
    return "NDM";
  else if( num == 873) 
    return "NDF";
  else if( num == 874) 
    return "NDP";
  else if( num == 875) 
    return "NDS";
  else if( num == 876) 
    return "NDT";
  else if( num == 877) 
    return "NDW";
  else if( num == 878) 
    return "NDY";
  else if( num == 879) 
    return "NDV";
  else if( num == 880) 
    return "NCA";
  else if( num == 881) 
    return "NCR";
  else if( num == 882) 
    return "NCN";
  else if( num == 883) 
    return "NCD";
  else if( num == 884) 
    return "NCC";
  else if( num == 885) 
    return "NCQ";
  else if( num == 886) 
    return "NCE";
  else if( num == 887) 
    return "NCG";
  else if( num == 888) 
    return "NCH";
  else if( num == 889) 
    return "NCI";
  else if( num == 890) 
    return "NCL";
  else if( num == 891) 
    return "NCK";
  else if( num == 892) 
    return "NCM";
  else if( num == 893) 
    return "NCF";
  else if( num == 894) 
    return "NCP";
  else if( num == 895) 
    return "NCS";
  else if( num == 896) 
    return "NCT";
  else if( num == 897) 
    return "NCW";
  else if( num == 898) 
    return "NCY";
  else if( num == 899) 
    return "NCV";
  else if( num == 900) 
    return "NQA";
  else if( num == 901) 
    return "NQR";
  else if( num == 902) 
    return "NQN";
  else if( num == 903) 
    return "NQD";
  else if( num == 904) 
    return "NQC";
  else if( num == 905) 
    return "NQQ";
  else if( num == 906) 
    return "NQE";
  else if( num == 907) 
    return "NQG";
  else if( num == 908) 
    return "NQH";
  else if( num == 909) 
    return "NQI";
  else if( num == 910) 
    return "NQL";
  else if( num == 911) 
    return "NQK";
  else if( num == 912) 
    return "NQM";
  else if( num == 913) 
    return "NQF";
  else if( num == 914) 
    return "NQP";
  else if( num == 915) 
    return "NQS";
  else if( num == 916) 
    return "NQT";
  else if( num == 917) 
    return "NQW";
  else if( num == 918) 
    return "NQY";
  else if( num == 919) 
    return "NQV";
  else if( num == 920) 
    return "NEA";
  else if( num == 921) 
    return "NER";
  else if( num == 922) 
    return "NEN";
  else if( num == 923) 
    return "NED";
  else if( num == 924) 
    return "NEC";
  else if( num == 925) 
    return "NEQ";
  else if( num == 926) 
    return "NEE";
  else if( num == 927) 
    return "NEG";
  else if( num == 928) 
    return "NEH";
  else if( num == 929) 
    return "NEI";
  else if( num == 930) 
    return "NEL";
  else if( num == 931) 
    return "NEK";
  else if( num == 932) 
    return "NEM";
  else if( num == 933) 
    return "NEF";
  else if( num == 934) 
    return "NEP";
  else if( num == 935) 
    return "NES";
  else if( num == 936) 
    return "NET";
  else if( num == 937) 
    return "NEW";
  else if( num == 938) 
    return "NEY";
  else if( num == 939) 
    return "NEV";
  else if( num == 940) 
    return "NGA";
  else if( num == 941) 
    return "NGR";
  else if( num == 942) 
    return "NGN";
  else if( num == 943) 
    return "NGD";
  else if( num == 944) 
    return "NGC";
  else if( num == 945) 
    return "NGQ";
  else if( num == 946) 
    return "NGE";
  else if( num == 947) 
    return "NGG";
  else if( num == 948) 
    return "NGH";
  else if( num == 949) 
    return "NGI";
  else if( num == 950) 
    return "NGL";
  else if( num == 951) 
    return "NGK";
  else if( num == 952) 
    return "NGM";
  else if( num == 953) 
    return "NGF";
  else if( num == 954) 
    return "NGP";
  else if( num == 955) 
    return "NGS";
  else if( num == 956) 
    return "NGT";
  else if( num == 957) 
    return "NGW";
  else if( num == 958) 
    return "NGY";
  else if( num == 959) 
    return "NGV";
  else if( num == 960) 
    return "NHA";
  else if( num == 961) 
    return "NHR";
  else if( num == 962) 
    return "NHN";
  else if( num == 963) 
    return "NHD";
  else if( num == 964) 
    return "NHC";
  else if( num == 965) 
    return "NHQ";
  else if( num == 966) 
    return "NHE";
  else if( num == 967) 
    return "NHG";
  else if( num == 968) 
    return "NHH";
  else if( num == 969) 
    return "NHI";
  else if( num == 970) 
    return "NHL";
  else if( num == 971) 
    return "NHK";
  else if( num == 972) 
    return "NHM";
  else if( num == 973) 
    return "NHF";
  else if( num == 974) 
    return "NHP";
  else if( num == 975) 
    return "NHS";
  else if( num == 976) 
    return "NHT";
  else if( num == 977) 
    return "NHW";
  else if( num == 978) 
    return "NHY";
  else if( num == 979) 
    return "NHV";
  else if( num == 980) 
    return "NIA";
  else if( num == 981) 
    return "NIR";
  else if( num == 982) 
    return "NIN";
  else if( num == 983) 
    return "NID";
  else if( num == 984) 
    return "NIC";
  else if( num == 985) 
    return "NIQ";
  else if( num == 986) 
    return "NIE";
  else if( num == 987) 
    return "NIG";
  else if( num == 988) 
    return "NIH";
  else if( num == 989) 
    return "NII";
  else if( num == 990) 
    return "NIL";
  else if( num == 991) 
    return "NIK";
  else if( num == 992) 
    return "NIM";
  else if( num == 993) 
    return "NIF";
  else if( num == 994) 
    return "NIP";
  else if( num == 995) 
    return "NIS";
  else if( num == 996) 
    return "NIT";
  else if( num == 997) 
    return "NIW";
  else if( num == 998) 
    return "NIY";
  else if( num == 999) 
    return "NIV";
  else if( num == 1000) 
    return "NLA";
  else if( num == 1001) 
    return "NLR";
  else if( num == 1002) 
    return "NLN";
  else if( num == 1003) 
    return "NLD";
  else if( num == 1004) 
    return "NLC";
  else if( num == 1005) 
    return "NLQ";
  else if( num == 1006) 
    return "NLE";
  else if( num == 1007) 
    return "NLG";
  else if( num == 1008) 
    return "NLH";
  else if( num == 1009) 
    return "NLI";
  else if( num == 1010) 
    return "NLL";
  else if( num == 1011) 
    return "NLK";
  else if( num == 1012) 
    return "NLM";
  else if( num == 1013) 
    return "NLF";
  else if( num == 1014) 
    return "NLP";
  else if( num == 1015) 
    return "NLS";
  else if( num == 1016) 
    return "NLT";
  else if( num == 1017) 
    return "NLW";
  else if( num == 1018) 
    return "NLY";
  else if( num == 1019) 
    return "NLV";
  else if( num == 1020) 
    return "NKA";
  else if( num == 1021) 
    return "NKR";
  else if( num == 1022) 
    return "NKN";
  else if( num == 1023) 
    return "NKD";
  else if( num == 1024) 
    return "NKC";
  else if( num == 1025) 
    return "NKQ";
  else if( num == 1026) 
    return "NKE";
  else if( num == 1027) 
    return "NKG";
  else if( num == 1028) 
    return "NKH";
  else if( num == 1029) 
    return "NKI";
  else if( num == 1030) 
    return "NKL";
  else if( num == 1031) 
    return "NKK";
  else if( num == 1032) 
    return "NKM";
  else if( num == 1033) 
    return "NKF";
  else if( num == 1034) 
    return "NKP";
  else if( num == 1035) 
    return "NKS";
  else if( num == 1036) 
    return "NKT";
  else if( num == 1037) 
    return "NKW";
  else if( num == 1038) 
    return "NKY";
  else if( num == 1039) 
    return "NKV";
  else if( num == 1040) 
    return "NMA";
  else if( num == 1041) 
    return "NMR";
  else if( num == 1042) 
    return "NMN";
  else if( num == 1043) 
    return "NMD";
  else if( num == 1044) 
    return "NMC";
  else if( num == 1045) 
    return "NMQ";
  else if( num == 1046) 
    return "NME";
  else if( num == 1047) 
    return "NMG";
  else if( num == 1048) 
    return "NMH";
  else if( num == 1049) 
    return "NMI";
  else if( num == 1050) 
    return "NML";
  else if( num == 1051) 
    return "NMK";
  else if( num == 1052) 
    return "NMM";
  else if( num == 1053) 
    return "NMF";
  else if( num == 1054) 
    return "NMP";
  else if( num == 1055) 
    return "NMS";
  else if( num == 1056) 
    return "NMT";
  else if( num == 1057) 
    return "NMW";
  else if( num == 1058) 
    return "NMY";
  else if( num == 1059) 
    return "NMV";
  else if( num == 1060) 
    return "NFA";
  else if( num == 1061) 
    return "NFR";
  else if( num == 1062) 
    return "NFN";
  else if( num == 1063) 
    return "NFD";
  else if( num == 1064) 
    return "NFC";
  else if( num == 1065) 
    return "NFQ";
  else if( num == 1066) 
    return "NFE";
  else if( num == 1067) 
    return "NFG";
  else if( num == 1068) 
    return "NFH";
  else if( num == 1069) 
    return "NFI";
  else if( num == 1070) 
    return "NFL";
  else if( num == 1071) 
    return "NFK";
  else if( num == 1072) 
    return "NFM";
  else if( num == 1073) 
    return "NFF";
  else if( num == 1074) 
    return "NFP";
  else if( num == 1075) 
    return "NFS";
  else if( num == 1076) 
    return "NFT";
  else if( num == 1077) 
    return "NFW";
  else if( num == 1078) 
    return "NFY";
  else if( num == 1079) 
    return "NFV";
  else if( num == 1080) 
    return "NPA";
  else if( num == 1081) 
    return "NPR";
  else if( num == 1082) 
    return "NPN";
  else if( num == 1083) 
    return "NPD";
  else if( num == 1084) 
    return "NPC";
  else if( num == 1085) 
    return "NPQ";
  else if( num == 1086) 
    return "NPE";
  else if( num == 1087) 
    return "NPG";
  else if( num == 1088) 
    return "NPH";
  else if( num == 1089) 
    return "NPI";
  else if( num == 1090) 
    return "NPL";
  else if( num == 1091) 
    return "NPK";
  else if( num == 1092) 
    return "NPM";
  else if( num == 1093) 
    return "NPF";
  else if( num == 1094) 
    return "NPP";
  else if( num == 1095) 
    return "NPS";
  else if( num == 1096) 
    return "NPT";
  else if( num == 1097) 
    return "NPW";
  else if( num == 1098) 
    return "NPY";
  else if( num == 1099) 
    return "NPV";
  else if( num == 1100) 
    return "NSA";
  else if( num == 1101) 
    return "NSR";
  else if( num == 1102) 
    return "NSN";
  else if( num == 1103) 
    return "NSD";
  else if( num == 1104) 
    return "NSC";
  else if( num == 1105) 
    return "NSQ";
  else if( num == 1106) 
    return "NSE";
  else if( num == 1107) 
    return "NSG";
  else if( num == 1108) 
    return "NSH";
  else if( num == 1109) 
    return "NSI";
  else if( num == 1110) 
    return "NSL";
  else if( num == 1111) 
    return "NSK";
  else if( num == 1112) 
    return "NSM";
  else if( num == 1113) 
    return "NSF";
  else if( num == 1114) 
    return "NSP";
  else if( num == 1115) 
    return "NSS";
  else if( num == 1116) 
    return "NST";
  else if( num == 1117) 
    return "NSW";
  else if( num == 1118) 
    return "NSY";
  else if( num == 1119) 
    return "NSV";
  else if( num == 1120) 
    return "NTA";
  else if( num == 1121) 
    return "NTR";
  else if( num == 1122) 
    return "NTN";
  else if( num == 1123) 
    return "NTD";
  else if( num == 1124) 
    return "NTC";
  else if( num == 1125) 
    return "NTQ";
  else if( num == 1126) 
    return "NTE";
  else if( num == 1127) 
    return "NTG";
  else if( num == 1128) 
    return "NTH";
  else if( num == 1129) 
    return "NTI";
  else if( num == 1130) 
    return "NTL";
  else if( num == 1131) 
    return "NTK";
  else if( num == 1132) 
    return "NTM";
  else if( num == 1133) 
    return "NTF";
  else if( num == 1134) 
    return "NTP";
  else if( num == 1135) 
    return "NTS";
  else if( num == 1136) 
    return "NTT";
  else if( num == 1137) 
    return "NTW";
  else if( num == 1138) 
    return "NTY";
  else if( num == 1139) 
    return "NTV";
  else if( num == 1140) 
    return "NWA";
  else if( num == 1141) 
    return "NWR";
  else if( num == 1142) 
    return "NWN";
  else if( num == 1143) 
    return "NWD";
  else if( num == 1144) 
    return "NWC";
  else if( num == 1145) 
    return "NWQ";
  else if( num == 1146) 
    return "NWE";
  else if( num == 1147) 
    return "NWG";
  else if( num == 1148) 
    return "NWH";
  else if( num == 1149) 
    return "NWI";
  else if( num == 1150) 
    return "NWL";
  else if( num == 1151) 
    return "NWK";
  else if( num == 1152) 
    return "NWM";
  else if( num == 1153) 
    return "NWF";
  else if( num == 1154) 
    return "NWP";
  else if( num == 1155) 
    return "NWS";
  else if( num == 1156) 
    return "NWT";
  else if( num == 1157) 
    return "NWW";
  else if( num == 1158) 
    return "NWY";
  else if( num == 1159) 
    return "NWV";
  else if( num == 1160) 
    return "NYA";
  else if( num == 1161) 
    return "NYR";
  else if( num == 1162) 
    return "NYN";
  else if( num == 1163) 
    return "NYD";
  else if( num == 1164) 
    return "NYC";
  else if( num == 1165) 
    return "NYQ";
  else if( num == 1166) 
    return "NYE";
  else if( num == 1167) 
    return "NYG";
  else if( num == 1168) 
    return "NYH";
  else if( num == 1169) 
    return "NYI";
  else if( num == 1170) 
    return "NYL";
  else if( num == 1171) 
    return "NYK";
  else if( num == 1172) 
    return "NYM";
  else if( num == 1173) 
    return "NYF";
  else if( num == 1174) 
    return "NYP";
  else if( num == 1175) 
    return "NYS";
  else if( num == 1176) 
    return "NYT";
  else if( num == 1177) 
    return "NYW";
  else if( num == 1178) 
    return "NYY";
  else if( num == 1179) 
    return "NYV";
  else if( num == 1180) 
    return "NVA";
  else if( num == 1181) 
    return "NVR";
  else if( num == 1182) 
    return "NVN";
  else if( num == 1183) 
    return "NVD";
  else if( num == 1184) 
    return "NVC";
  else if( num == 1185) 
    return "NVQ";
  else if( num == 1186) 
    return "NVE";
  else if( num == 1187) 
    return "NVG";
  else if( num == 1188) 
    return "NVH";
  else if( num == 1189) 
    return "NVI";
  else if( num == 1190) 
    return "NVL";
  else if( num == 1191) 
    return "NVK";
  else if( num == 1192) 
    return "NVM";
  else if( num == 1193) 
    return "NVF";
  else if( num == 1194) 
    return "NVP";
  else if( num == 1195) 
    return "NVS";
  else if( num == 1196) 
    return "NVT";
  else if( num == 1197) 
    return "NVW";
  else if( num == 1198) 
    return "NVY";
  else if( num == 1199) 
    return "NVV";
  else if( num == 1200) 
    return "DAA";
  else if( num == 1201) 
    return "DAR";
  else if( num == 1202) 
    return "DAN";
  else if( num == 1203) 
    return "DAD";
  else if( num == 1204) 
    return "DAC";
  else if( num == 1205) 
    return "DAQ";
  else if( num == 1206) 
    return "DAE";
  else if( num == 1207) 
    return "DAG";
  else if( num == 1208) 
    return "DAH";
  else if( num == 1209) 
    return "DAI";
  else if( num == 1210) 
    return "DAL";
  else if( num == 1211) 
    return "DAK";
  else if( num == 1212) 
    return "DAM";
  else if( num == 1213) 
    return "DAF";
  else if( num == 1214) 
    return "DAP";
  else if( num == 1215) 
    return "DAS";
  else if( num == 1216) 
    return "DAT";
  else if( num == 1217) 
    return "DAW";
  else if( num == 1218) 
    return "DAY";
  else if( num == 1219) 
    return "DAV";
  else if( num == 1220) 
    return "DRA";
  else if( num == 1221) 
    return "DRR";
  else if( num == 1222) 
    return "DRN";
  else if( num == 1223) 
    return "DRD";
  else if( num == 1224) 
    return "DRC";
  else if( num == 1225) 
    return "DRQ";
  else if( num == 1226) 
    return "DRE";
  else if( num == 1227) 
    return "DRG";
  else if( num == 1228) 
    return "DRH";
  else if( num == 1229) 
    return "DRI";
  else if( num == 1230) 
    return "DRL";
  else if( num == 1231) 
    return "DRK";
  else if( num == 1232) 
    return "DRM";
  else if( num == 1233) 
    return "DRF";
  else if( num == 1234) 
    return "DRP";
  else if( num == 1235) 
    return "DRS";
  else if( num == 1236) 
    return "DRT";
  else if( num == 1237) 
    return "DRW";
  else if( num == 1238) 
    return "DRY";
  else if( num == 1239) 
    return "DRV";
  else if( num == 1240) 
    return "DNA";
  else if( num == 1241) 
    return "DNR";
  else if( num == 1242) 
    return "DNN";
  else if( num == 1243) 
    return "DND";
  else if( num == 1244) 
    return "DNC";
  else if( num == 1245) 
    return "DNQ";
  else if( num == 1246) 
    return "DNE";
  else if( num == 1247) 
    return "DNG";
  else if( num == 1248) 
    return "DNH";
  else if( num == 1249) 
    return "DNI";
  else if( num == 1250) 
    return "DNL";
  else if( num == 1251) 
    return "DNK";
  else if( num == 1252) 
    return "DNM";
  else if( num == 1253) 
    return "DNF";
  else if( num == 1254) 
    return "DNP";
  else if( num == 1255) 
    return "DNS";
  else if( num == 1256) 
    return "DNT";
  else if( num == 1257) 
    return "DNW";
  else if( num == 1258) 
    return "DNY";
  else if( num == 1259) 
    return "DNV";
  else if( num == 1260) 
    return "DDA";
  else if( num == 1261) 
    return "DDR";
  else if( num == 1262) 
    return "DDN";
  else if( num == 1263) 
    return "DDD";
  else if( num == 1264) 
    return "DDC";
  else if( num == 1265) 
    return "DDQ";
  else if( num == 1266) 
    return "DDE";
  else if( num == 1267) 
    return "DDG";
  else if( num == 1268) 
    return "DDH";
  else if( num == 1269) 
    return "DDI";
  else if( num == 1270) 
    return "DDL";
  else if( num == 1271) 
    return "DDK";
  else if( num == 1272) 
    return "DDM";
  else if( num == 1273) 
    return "DDF";
  else if( num == 1274) 
    return "DDP";
  else if( num == 1275) 
    return "DDS";
  else if( num == 1276) 
    return "DDT";
  else if( num == 1277) 
    return "DDW";
  else if( num == 1278) 
    return "DDY";
  else if( num == 1279) 
    return "DDV";
  else if( num == 1280) 
    return "DCA";
  else if( num == 1281) 
    return "DCR";
  else if( num == 1282) 
    return "DCN";
  else if( num == 1283) 
    return "DCD";
  else if( num == 1284) 
    return "DCC";
  else if( num == 1285) 
    return "DCQ";
  else if( num == 1286) 
    return "DCE";
  else if( num == 1287) 
    return "DCG";
  else if( num == 1288) 
    return "DCH";
  else if( num == 1289) 
    return "DCI";
  else if( num == 1290) 
    return "DCL";
  else if( num == 1291) 
    return "DCK";
  else if( num == 1292) 
    return "DCM";
  else if( num == 1293) 
    return "DCF";
  else if( num == 1294) 
    return "DCP";
  else if( num == 1295) 
    return "DCS";
  else if( num == 1296) 
    return "DCT";
  else if( num == 1297) 
    return "DCW";
  else if( num == 1298) 
    return "DCY";
  else if( num == 1299) 
    return "DCV";
  else if( num == 1300) 
    return "DQA";
  else if( num == 1301) 
    return "DQR";
  else if( num == 1302) 
    return "DQN";
  else if( num == 1303) 
    return "DQD";
  else if( num == 1304) 
    return "DQC";
  else if( num == 1305) 
    return "DQQ";
  else if( num == 1306) 
    return "DQE";
  else if( num == 1307) 
    return "DQG";
  else if( num == 1308) 
    return "DQH";
  else if( num == 1309) 
    return "DQI";
  else if( num == 1310) 
    return "DQL";
  else if( num == 1311) 
    return "DQK";
  else if( num == 1312) 
    return "DQM";
  else if( num == 1313) 
    return "DQF";
  else if( num == 1314) 
    return "DQP";
  else if( num == 1315) 
    return "DQS";
  else if( num == 1316) 
    return "DQT";
  else if( num == 1317) 
    return "DQW";
  else if( num == 1318) 
    return "DQY";
  else if( num == 1319) 
    return "DQV";
  else if( num == 1320) 
    return "DEA";
  else if( num == 1321) 
    return "DER";
  else if( num == 1322) 
    return "DEN";
  else if( num == 1323) 
    return "DED";
  else if( num == 1324) 
    return "DEC";
  else if( num == 1325) 
    return "DEQ";
  else if( num == 1326) 
    return "DEE";
  else if( num == 1327) 
    return "DEG";
  else if( num == 1328) 
    return "DEH";
  else if( num == 1329) 
    return "DEI";
  else if( num == 1330) 
    return "DEL";
  else if( num == 1331) 
    return "DEK";
  else if( num == 1332) 
    return "DEM";
  else if( num == 1333) 
    return "DEF";
  else if( num == 1334) 
    return "DEP";
  else if( num == 1335) 
    return "DES";
  else if( num == 1336) 
    return "DET";
  else if( num == 1337) 
    return "DEW";
  else if( num == 1338) 
    return "DEY";
  else if( num == 1339) 
    return "DEV";
  else if( num == 1340) 
    return "DGA";
  else if( num == 1341) 
    return "DGR";
  else if( num == 1342) 
    return "DGN";
  else if( num == 1343) 
    return "DGD";
  else if( num == 1344) 
    return "DGC";
  else if( num == 1345) 
    return "DGQ";
  else if( num == 1346) 
    return "DGE";
  else if( num == 1347) 
    return "DGG";
  else if( num == 1348) 
    return "DGH";
  else if( num == 1349) 
    return "DGI";
  else if( num == 1350) 
    return "DGL";
  else if( num == 1351) 
    return "DGK";
  else if( num == 1352) 
    return "DGM";
  else if( num == 1353) 
    return "DGF";
  else if( num == 1354) 
    return "DGP";
  else if( num == 1355) 
    return "DGS";
  else if( num == 1356) 
    return "DGT";
  else if( num == 1357) 
    return "DGW";
  else if( num == 1358) 
    return "DGY";
  else if( num == 1359) 
    return "DGV";
  else if( num == 1360) 
    return "DHA";
  else if( num == 1361) 
    return "DHR";
  else if( num == 1362) 
    return "DHN";
  else if( num == 1363) 
    return "DHD";
  else if( num == 1364) 
    return "DHC";
  else if( num == 1365) 
    return "DHQ";
  else if( num == 1366) 
    return "DHE";
  else if( num == 1367) 
    return "DHG";
  else if( num == 1368) 
    return "DHH";
  else if( num == 1369) 
    return "DHI";
  else if( num == 1370) 
    return "DHL";
  else if( num == 1371) 
    return "DHK";
  else if( num == 1372) 
    return "DHM";
  else if( num == 1373) 
    return "DHF";
  else if( num == 1374) 
    return "DHP";
  else if( num == 1375) 
    return "DHS";
  else if( num == 1376) 
    return "DHT";
  else if( num == 1377) 
    return "DHW";
  else if( num == 1378) 
    return "DHY";
  else if( num == 1379) 
    return "DHV";
  else if( num == 1380) 
    return "DIA";
  else if( num == 1381) 
    return "DIR";
  else if( num == 1382) 
    return "DIN";
  else if( num == 1383) 
    return "DID";
  else if( num == 1384) 
    return "DIC";
  else if( num == 1385) 
    return "DIQ";
  else if( num == 1386) 
    return "DIE";
  else if( num == 1387) 
    return "DIG";
  else if( num == 1388) 
    return "DIH";
  else if( num == 1389) 
    return "DII";
  else if( num == 1390) 
    return "DIL";
  else if( num == 1391) 
    return "DIK";
  else if( num == 1392) 
    return "DIM";
  else if( num == 1393) 
    return "DIF";
  else if( num == 1394) 
    return "DIP";
  else if( num == 1395) 
    return "DIS";
  else if( num == 1396) 
    return "DIT";
  else if( num == 1397) 
    return "DIW";
  else if( num == 1398) 
    return "DIY";
  else if( num == 1399) 
    return "DIV";
  else if( num == 1400) 
    return "DLA";
  else if( num == 1401) 
    return "DLR";
  else if( num == 1402) 
    return "DLN";
  else if( num == 1403) 
    return "DLD";
  else if( num == 1404) 
    return "DLC";
  else if( num == 1405) 
    return "DLQ";
  else if( num == 1406) 
    return "DLE";
  else if( num == 1407) 
    return "DLG";
  else if( num == 1408) 
    return "DLH";
  else if( num == 1409) 
    return "DLI";
  else if( num == 1410) 
    return "DLL";
  else if( num == 1411) 
    return "DLK";
  else if( num == 1412) 
    return "DLM";
  else if( num == 1413) 
    return "DLF";
  else if( num == 1414) 
    return "DLP";
  else if( num == 1415) 
    return "DLS";
  else if( num == 1416) 
    return "DLT";
  else if( num == 1417) 
    return "DLW";
  else if( num == 1418) 
    return "DLY";
  else if( num == 1419) 
    return "DLV";
  else if( num == 1420) 
    return "DKA";
  else if( num == 1421) 
    return "DKR";
  else if( num == 1422) 
    return "DKN";
  else if( num == 1423) 
    return "DKD";
  else if( num == 1424) 
    return "DKC";
  else if( num == 1425) 
    return "DKQ";
  else if( num == 1426) 
    return "DKE";
  else if( num == 1427) 
    return "DKG";
  else if( num == 1428) 
    return "DKH";
  else if( num == 1429) 
    return "DKI";
  else if( num == 1430) 
    return "DKL";
  else if( num == 1431) 
    return "DKK";
  else if( num == 1432) 
    return "DKM";
  else if( num == 1433) 
    return "DKF";
  else if( num == 1434) 
    return "DKP";
  else if( num == 1435) 
    return "DKS";
  else if( num == 1436) 
    return "DKT";
  else if( num == 1437) 
    return "DKW";
  else if( num == 1438) 
    return "DKY";
  else if( num == 1439) 
    return "DKV";
  else if( num == 1440) 
    return "DMA";
  else if( num == 1441) 
    return "DMR";
  else if( num == 1442) 
    return "DMN";
  else if( num == 1443) 
    return "DMD";
  else if( num == 1444) 
    return "DMC";
  else if( num == 1445) 
    return "DMQ";
  else if( num == 1446) 
    return "DME";
  else if( num == 1447) 
    return "DMG";
  else if( num == 1448) 
    return "DMH";
  else if( num == 1449) 
    return "DMI";
  else if( num == 1450) 
    return "DML";
  else if( num == 1451) 
    return "DMK";
  else if( num == 1452) 
    return "DMM";
  else if( num == 1453) 
    return "DMF";
  else if( num == 1454) 
    return "DMP";
  else if( num == 1455) 
    return "DMS";
  else if( num == 1456) 
    return "DMT";
  else if( num == 1457) 
    return "DMW";
  else if( num == 1458) 
    return "DMY";
  else if( num == 1459) 
    return "DMV";
  else if( num == 1460) 
    return "DFA";
  else if( num == 1461) 
    return "DFR";
  else if( num == 1462) 
    return "DFN";
  else if( num == 1463) 
    return "DFD";
  else if( num == 1464) 
    return "DFC";
  else if( num == 1465) 
    return "DFQ";
  else if( num == 1466) 
    return "DFE";
  else if( num == 1467) 
    return "DFG";
  else if( num == 1468) 
    return "DFH";
  else if( num == 1469) 
    return "DFI";
  else if( num == 1470) 
    return "DFL";
  else if( num == 1471) 
    return "DFK";
  else if( num == 1472) 
    return "DFM";
  else if( num == 1473) 
    return "DFF";
  else if( num == 1474) 
    return "DFP";
  else if( num == 1475) 
    return "DFS";
  else if( num == 1476) 
    return "DFT";
  else if( num == 1477) 
    return "DFW";
  else if( num == 1478) 
    return "DFY";
  else if( num == 1479) 
    return "DFV";
  else if( num == 1480) 
    return "DPA";
  else if( num == 1481) 
    return "DPR";
  else if( num == 1482) 
    return "DPN";
  else if( num == 1483) 
    return "DPD";
  else if( num == 1484) 
    return "DPC";
  else if( num == 1485) 
    return "DPQ";
  else if( num == 1486) 
    return "DPE";
  else if( num == 1487) 
    return "DPG";
  else if( num == 1488) 
    return "DPH";
  else if( num == 1489) 
    return "DPI";
  else if( num == 1490) 
    return "DPL";
  else if( num == 1491) 
    return "DPK";
  else if( num == 1492) 
    return "DPM";
  else if( num == 1493) 
    return "DPF";
  else if( num == 1494) 
    return "DPP";
  else if( num == 1495) 
    return "DPS";
  else if( num == 1496) 
    return "DPT";
  else if( num == 1497) 
    return "DPW";
  else if( num == 1498) 
    return "DPY";
  else if( num == 1499) 
    return "DPV";
  else if( num == 1500) 
    return "DSA";
  else if( num == 1501) 
    return "DSR";
  else if( num == 1502) 
    return "DSN";
  else if( num == 1503) 
    return "DSD";
  else if( num == 1504) 
    return "DSC";
  else if( num == 1505) 
    return "DSQ";
  else if( num == 1506) 
    return "DSE";
  else if( num == 1507) 
    return "DSG";
  else if( num == 1508) 
    return "DSH";
  else if( num == 1509) 
    return "DSI";
  else if( num == 1510) 
    return "DSL";
  else if( num == 1511) 
    return "DSK";
  else if( num == 1512) 
    return "DSM";
  else if( num == 1513) 
    return "DSF";
  else if( num == 1514) 
    return "DSP";
  else if( num == 1515) 
    return "DSS";
  else if( num == 1516) 
    return "DST";
  else if( num == 1517) 
    return "DSW";
  else if( num == 1518) 
    return "DSY";
  else if( num == 1519) 
    return "DSV";
  else if( num == 1520) 
    return "DTA";
  else if( num == 1521) 
    return "DTR";
  else if( num == 1522) 
    return "DTN";
  else if( num == 1523) 
    return "DTD";
  else if( num == 1524) 
    return "DTC";
  else if( num == 1525) 
    return "DTQ";
  else if( num == 1526) 
    return "DTE";
  else if( num == 1527) 
    return "DTG";
  else if( num == 1528) 
    return "DTH";
  else if( num == 1529) 
    return "DTI";
  else if( num == 1530) 
    return "DTL";
  else if( num == 1531) 
    return "DTK";
  else if( num == 1532) 
    return "DTM";
  else if( num == 1533) 
    return "DTF";
  else if( num == 1534) 
    return "DTP";
  else if( num == 1535) 
    return "DTS";
  else if( num == 1536) 
    return "DTT";
  else if( num == 1537) 
    return "DTW";
  else if( num == 1538) 
    return "DTY";
  else if( num == 1539) 
    return "DTV";
  else if( num == 1540) 
    return "DWA";
  else if( num == 1541) 
    return "DWR";
  else if( num == 1542) 
    return "DWN";
  else if( num == 1543) 
    return "DWD";
  else if( num == 1544) 
    return "DWC";
  else if( num == 1545) 
    return "DWQ";
  else if( num == 1546) 
    return "DWE";
  else if( num == 1547) 
    return "DWG";
  else if( num == 1548) 
    return "DWH";
  else if( num == 1549) 
    return "DWI";
  else if( num == 1550) 
    return "DWL";
  else if( num == 1551) 
    return "DWK";
  else if( num == 1552) 
    return "DWM";
  else if( num == 1553) 
    return "DWF";
  else if( num == 1554) 
    return "DWP";
  else if( num == 1555) 
    return "DWS";
  else if( num == 1556) 
    return "DWT";
  else if( num == 1557) 
    return "DWW";
  else if( num == 1558) 
    return "DWY";
  else if( num == 1559) 
    return "DWV";
  else if( num == 1560) 
    return "DYA";
  else if( num == 1561) 
    return "DYR";
  else if( num == 1562) 
    return "DYN";
  else if( num == 1563) 
    return "DYD";
  else if( num == 1564) 
    return "DYC";
  else if( num == 1565) 
    return "DYQ";
  else if( num == 1566) 
    return "DYE";
  else if( num == 1567) 
    return "DYG";
  else if( num == 1568) 
    return "DYH";
  else if( num == 1569) 
    return "DYI";
  else if( num == 1570) 
    return "DYL";
  else if( num == 1571) 
    return "DYK";
  else if( num == 1572) 
    return "DYM";
  else if( num == 1573) 
    return "DYF";
  else if( num == 1574) 
    return "DYP";
  else if( num == 1575) 
    return "DYS";
  else if( num == 1576) 
    return "DYT";
  else if( num == 1577) 
    return "DYW";
  else if( num == 1578) 
    return "DYY";
  else if( num == 1579) 
    return "DYV";
  else if( num == 1580) 
    return "DVA";
  else if( num == 1581) 
    return "DVR";
  else if( num == 1582) 
    return "DVN";
  else if( num == 1583) 
    return "DVD";
  else if( num == 1584) 
    return "DVC";
  else if( num == 1585) 
    return "DVQ";
  else if( num == 1586) 
    return "DVE";
  else if( num == 1587) 
    return "DVG";
  else if( num == 1588) 
    return "DVH";
  else if( num == 1589) 
    return "DVI";
  else if( num == 1590) 
    return "DVL";
  else if( num == 1591) 
    return "DVK";
  else if( num == 1592) 
    return "DVM";
  else if( num == 1593) 
    return "DVF";
  else if( num == 1594) 
    return "DVP";
  else if( num == 1595) 
    return "DVS";
  else if( num == 1596) 
    return "DVT";
  else if( num == 1597) 
    return "DVW";
  else if( num == 1598) 
    return "DVY";
  else if( num == 1599) 
    return "DVV";
  else if( num == 1600) 
    return "CAA";
  else if( num == 1601) 
    return "CAR";
  else if( num == 1602) 
    return "CAN";
  else if( num == 1603) 
    return "CAD";
  else if( num == 1604) 
    return "CAC";
  else if( num == 1605) 
    return "CAQ";
  else if( num == 1606) 
    return "CAE";
  else if( num == 1607) 
    return "CAG";
  else if( num == 1608) 
    return "CAH";
  else if( num == 1609) 
    return "CAI";
  else if( num == 1610) 
    return "CAL";
  else if( num == 1611) 
    return "CAK";
  else if( num == 1612) 
    return "CAM";
  else if( num == 1613) 
    return "CAF";
  else if( num == 1614) 
    return "CAP";
  else if( num == 1615) 
    return "CAS";
  else if( num == 1616) 
    return "CAT";
  else if( num == 1617) 
    return "CAW";
  else if( num == 1618) 
    return "CAY";
  else if( num == 1619) 
    return "CAV";
  else if( num == 1620) 
    return "CRA";
  else if( num == 1621) 
    return "CRR";
  else if( num == 1622) 
    return "CRN";
  else if( num == 1623) 
    return "CRD";
  else if( num == 1624) 
    return "CRC";
  else if( num == 1625) 
    return "CRQ";
  else if( num == 1626) 
    return "CRE";
  else if( num == 1627) 
    return "CRG";
  else if( num == 1628) 
    return "CRH";
  else if( num == 1629) 
    return "CRI";
  else if( num == 1630) 
    return "CRL";
  else if( num == 1631) 
    return "CRK";
  else if( num == 1632) 
    return "CRM";
  else if( num == 1633) 
    return "CRF";
  else if( num == 1634) 
    return "CRP";
  else if( num == 1635) 
    return "CRS";
  else if( num == 1636) 
    return "CRT";
  else if( num == 1637) 
    return "CRW";
  else if( num == 1638) 
    return "CRY";
  else if( num == 1639) 
    return "CRV";
  else if( num == 1640) 
    return "CNA";
  else if( num == 1641) 
    return "CNR";
  else if( num == 1642) 
    return "CNN";
  else if( num == 1643) 
    return "CND";
  else if( num == 1644) 
    return "CNC";
  else if( num == 1645) 
    return "CNQ";
  else if( num == 1646) 
    return "CNE";
  else if( num == 1647) 
    return "CNG";
  else if( num == 1648) 
    return "CNH";
  else if( num == 1649) 
    return "CNI";
  else if( num == 1650) 
    return "CNL";
  else if( num == 1651) 
    return "CNK";
  else if( num == 1652) 
    return "CNM";
  else if( num == 1653) 
    return "CNF";
  else if( num == 1654) 
    return "CNP";
  else if( num == 1655) 
    return "CNS";
  else if( num == 1656) 
    return "CNT";
  else if( num == 1657) 
    return "CNW";
  else if( num == 1658) 
    return "CNY";
  else if( num == 1659) 
    return "CNV";
  else if( num == 1660) 
    return "CDA";
  else if( num == 1661) 
    return "CDR";
  else if( num == 1662) 
    return "CDN";
  else if( num == 1663) 
    return "CDD";
  else if( num == 1664) 
    return "CDC";
  else if( num == 1665) 
    return "CDQ";
  else if( num == 1666) 
    return "CDE";
  else if( num == 1667) 
    return "CDG";
  else if( num == 1668) 
    return "CDH";
  else if( num == 1669) 
    return "CDI";
  else if( num == 1670) 
    return "CDL";
  else if( num == 1671) 
    return "CDK";
  else if( num == 1672) 
    return "CDM";
  else if( num == 1673) 
    return "CDF";
  else if( num == 1674) 
    return "CDP";
  else if( num == 1675) 
    return "CDS";
  else if( num == 1676) 
    return "CDT";
  else if( num == 1677) 
    return "CDW";
  else if( num == 1678) 
    return "CDY";
  else if( num == 1679) 
    return "CDV";
  else if( num == 1680) 
    return "CCA";
  else if( num == 1681) 
    return "CCR";
  else if( num == 1682) 
    return "CCN";
  else if( num == 1683) 
    return "CCD";
  else if( num == 1684) 
    return "CCC";
  else if( num == 1685) 
    return "CCQ";
  else if( num == 1686) 
    return "CCE";
  else if( num == 1687) 
    return "CCG";
  else if( num == 1688) 
    return "CCH";
  else if( num == 1689) 
    return "CCI";
  else if( num == 1690) 
    return "CCL";
  else if( num == 1691) 
    return "CCK";
  else if( num == 1692) 
    return "CCM";
  else if( num == 1693) 
    return "CCF";
  else if( num == 1694) 
    return "CCP";
  else if( num == 1695) 
    return "CCS";
  else if( num == 1696) 
    return "CCT";
  else if( num == 1697) 
    return "CCW";
  else if( num == 1698) 
    return "CCY";
  else if( num == 1699) 
    return "CCV";
  else if( num == 1700) 
    return "CQA";
  else if( num == 1701) 
    return "CQR";
  else if( num == 1702) 
    return "CQN";
  else if( num == 1703) 
    return "CQD";
  else if( num == 1704) 
    return "CQC";
  else if( num == 1705) 
    return "CQQ";
  else if( num == 1706) 
    return "CQE";
  else if( num == 1707) 
    return "CQG";
  else if( num == 1708) 
    return "CQH";
  else if( num == 1709) 
    return "CQI";
  else if( num == 1710) 
    return "CQL";
  else if( num == 1711) 
    return "CQK";
  else if( num == 1712) 
    return "CQM";
  else if( num == 1713) 
    return "CQF";
  else if( num == 1714) 
    return "CQP";
  else if( num == 1715) 
    return "CQS";
  else if( num == 1716) 
    return "CQT";
  else if( num == 1717) 
    return "CQW";
  else if( num == 1718) 
    return "CQY";
  else if( num == 1719) 
    return "CQV";
  else if( num == 1720) 
    return "CEA";
  else if( num == 1721) 
    return "CER";
  else if( num == 1722) 
    return "CEN";
  else if( num == 1723) 
    return "CED";
  else if( num == 1724) 
    return "CEC";
  else if( num == 1725) 
    return "CEQ";
  else if( num == 1726) 
    return "CEE";
  else if( num == 1727) 
    return "CEG";
  else if( num == 1728) 
    return "CEH";
  else if( num == 1729) 
    return "CEI";
  else if( num == 1730) 
    return "CEL";
  else if( num == 1731) 
    return "CEK";
  else if( num == 1732) 
    return "CEM";
  else if( num == 1733) 
    return "CEF";
  else if( num == 1734) 
    return "CEP";
  else if( num == 1735) 
    return "CES";
  else if( num == 1736) 
    return "CET";
  else if( num == 1737) 
    return "CEW";
  else if( num == 1738) 
    return "CEY";
  else if( num == 1739) 
    return "CEV";
  else if( num == 1740) 
    return "CGA";
  else if( num == 1741) 
    return "CGR";
  else if( num == 1742) 
    return "CGN";
  else if( num == 1743) 
    return "CGD";
  else if( num == 1744) 
    return "CGC";
  else if( num == 1745) 
    return "CGQ";
  else if( num == 1746) 
    return "CGE";
  else if( num == 1747) 
    return "CGG";
  else if( num == 1748) 
    return "CGH";
  else if( num == 1749) 
    return "CGI";
  else if( num == 1750) 
    return "CGL";
  else if( num == 1751) 
    return "CGK";
  else if( num == 1752) 
    return "CGM";
  else if( num == 1753) 
    return "CGF";
  else if( num == 1754) 
    return "CGP";
  else if( num == 1755) 
    return "CGS";
  else if( num == 1756) 
    return "CGT";
  else if( num == 1757) 
    return "CGW";
  else if( num == 1758) 
    return "CGY";
  else if( num == 1759) 
    return "CGV";
  else if( num == 1760) 
    return "CHA";
  else if( num == 1761) 
    return "CHR";
  else if( num == 1762) 
    return "CHN";
  else if( num == 1763) 
    return "CHD";
  else if( num == 1764) 
    return "CHC";
  else if( num == 1765) 
    return "CHQ";
  else if( num == 1766) 
    return "CHE";
  else if( num == 1767) 
    return "CHG";
  else if( num == 1768) 
    return "CHH";
  else if( num == 1769) 
    return "CHI";
  else if( num == 1770) 
    return "CHL";
  else if( num == 1771) 
    return "CHK";
  else if( num == 1772) 
    return "CHM";
  else if( num == 1773) 
    return "CHF";
  else if( num == 1774) 
    return "CHP";
  else if( num == 1775) 
    return "CHS";
  else if( num == 1776) 
    return "CHT";
  else if( num == 1777) 
    return "CHW";
  else if( num == 1778) 
    return "CHY";
  else if( num == 1779) 
    return "CHV";
  else if( num == 1780) 
    return "CIA";
  else if( num == 1781) 
    return "CIR";
  else if( num == 1782) 
    return "CIN";
  else if( num == 1783) 
    return "CID";
  else if( num == 1784) 
    return "CIC";
  else if( num == 1785) 
    return "CIQ";
  else if( num == 1786) 
    return "CIE";
  else if( num == 1787) 
    return "CIG";
  else if( num == 1788) 
    return "CIH";
  else if( num == 1789) 
    return "CII";
  else if( num == 1790) 
    return "CIL";
  else if( num == 1791) 
    return "CIK";
  else if( num == 1792) 
    return "CIM";
  else if( num == 1793) 
    return "CIF";
  else if( num == 1794) 
    return "CIP";
  else if( num == 1795) 
    return "CIS";
  else if( num == 1796) 
    return "CIT";
  else if( num == 1797) 
    return "CIW";
  else if( num == 1798) 
    return "CIY";
  else if( num == 1799) 
    return "CIV";
  else if( num == 1800) 
    return "CLA";
  else if( num == 1801) 
    return "CLR";
  else if( num == 1802) 
    return "CLN";
  else if( num == 1803) 
    return "CLD";
  else if( num == 1804) 
    return "CLC";
  else if( num == 1805) 
    return "CLQ";
  else if( num == 1806) 
    return "CLE";
  else if( num == 1807) 
    return "CLG";
  else if( num == 1808) 
    return "CLH";
  else if( num == 1809) 
    return "CLI";
  else if( num == 1810) 
    return "CLL";
  else if( num == 1811) 
    return "CLK";
  else if( num == 1812) 
    return "CLM";
  else if( num == 1813) 
    return "CLF";
  else if( num == 1814) 
    return "CLP";
  else if( num == 1815) 
    return "CLS";
  else if( num == 1816) 
    return "CLT";
  else if( num == 1817) 
    return "CLW";
  else if( num == 1818) 
    return "CLY";
  else if( num == 1819) 
    return "CLV";
  else if( num == 1820) 
    return "CKA";
  else if( num == 1821) 
    return "CKR";
  else if( num == 1822) 
    return "CKN";
  else if( num == 1823) 
    return "CKD";
  else if( num == 1824) 
    return "CKC";
  else if( num == 1825) 
    return "CKQ";
  else if( num == 1826) 
    return "CKE";
  else if( num == 1827) 
    return "CKG";
  else if( num == 1828) 
    return "CKH";
  else if( num == 1829) 
    return "CKI";
  else if( num == 1830) 
    return "CKL";
  else if( num == 1831) 
    return "CKK";
  else if( num == 1832) 
    return "CKM";
  else if( num == 1833) 
    return "CKF";
  else if( num == 1834) 
    return "CKP";
  else if( num == 1835) 
    return "CKS";
  else if( num == 1836) 
    return "CKT";
  else if( num == 1837) 
    return "CKW";
  else if( num == 1838) 
    return "CKY";
  else if( num == 1839) 
    return "CKV";
  else if( num == 1840) 
    return "CMA";
  else if( num == 1841) 
    return "CMR";
  else if( num == 1842) 
    return "CMN";
  else if( num == 1843) 
    return "CMD";
  else if( num == 1844) 
    return "CMC";
  else if( num == 1845) 
    return "CMQ";
  else if( num == 1846) 
    return "CME";
  else if( num == 1847) 
    return "CMG";
  else if( num == 1848) 
    return "CMH";
  else if( num == 1849) 
    return "CMI";
  else if( num == 1850) 
    return "CML";
  else if( num == 1851) 
    return "CMK";
  else if( num == 1852) 
    return "CMM";
  else if( num == 1853) 
    return "CMF";
  else if( num == 1854) 
    return "CMP";
  else if( num == 1855) 
    return "CMS";
  else if( num == 1856) 
    return "CMT";
  else if( num == 1857) 
    return "CMW";
  else if( num == 1858) 
    return "CMY";
  else if( num == 1859) 
    return "CMV";
  else if( num == 1860) 
    return "CFA";
  else if( num == 1861) 
    return "CFR";
  else if( num == 1862) 
    return "CFN";
  else if( num == 1863) 
    return "CFD";
  else if( num == 1864) 
    return "CFC";
  else if( num == 1865) 
    return "CFQ";
  else if( num == 1866) 
    return "CFE";
  else if( num == 1867) 
    return "CFG";
  else if( num == 1868) 
    return "CFH";
  else if( num == 1869) 
    return "CFI";
  else if( num == 1870) 
    return "CFL";
  else if( num == 1871) 
    return "CFK";
  else if( num == 1872) 
    return "CFM";
  else if( num == 1873) 
    return "CFF";
  else if( num == 1874) 
    return "CFP";
  else if( num == 1875) 
    return "CFS";
  else if( num == 1876) 
    return "CFT";
  else if( num == 1877) 
    return "CFW";
  else if( num == 1878) 
    return "CFY";
  else if( num == 1879) 
    return "CFV";
  else if( num == 1880) 
    return "CPA";
  else if( num == 1881) 
    return "CPR";
  else if( num == 1882) 
    return "CPN";
  else if( num == 1883) 
    return "CPD";
  else if( num == 1884) 
    return "CPC";
  else if( num == 1885) 
    return "CPQ";
  else if( num == 1886) 
    return "CPE";
  else if( num == 1887) 
    return "CPG";
  else if( num == 1888) 
    return "CPH";
  else if( num == 1889) 
    return "CPI";
  else if( num == 1890) 
    return "CPL";
  else if( num == 1891) 
    return "CPK";
  else if( num == 1892) 
    return "CPM";
  else if( num == 1893) 
    return "CPF";
  else if( num == 1894) 
    return "CPP";
  else if( num == 1895) 
    return "CPS";
  else if( num == 1896) 
    return "CPT";
  else if( num == 1897) 
    return "CPW";
  else if( num == 1898) 
    return "CPY";
  else if( num == 1899) 
    return "CPV";
  else if( num == 1900) 
    return "CSA";
  else if( num == 1901) 
    return "CSR";
  else if( num == 1902) 
    return "CSN";
  else if( num == 1903) 
    return "CSD";
  else if( num == 1904) 
    return "CSC";
  else if( num == 1905) 
    return "CSQ";
  else if( num == 1906) 
    return "CSE";
  else if( num == 1907) 
    return "CSG";
  else if( num == 1908) 
    return "CSH";
  else if( num == 1909) 
    return "CSI";
  else if( num == 1910) 
    return "CSL";
  else if( num == 1911) 
    return "CSK";
  else if( num == 1912) 
    return "CSM";
  else if( num == 1913) 
    return "CSF";
  else if( num == 1914) 
    return "CSP";
  else if( num == 1915) 
    return "CSS";
  else if( num == 1916) 
    return "CST";
  else if( num == 1917) 
    return "CSW";
  else if( num == 1918) 
    return "CSY";
  else if( num == 1919) 
    return "CSV";
  else if( num == 1920) 
    return "CTA";
  else if( num == 1921) 
    return "CTR";
  else if( num == 1922) 
    return "CTN";
  else if( num == 1923) 
    return "CTD";
  else if( num == 1924) 
    return "CTC";
  else if( num == 1925) 
    return "CTQ";
  else if( num == 1926) 
    return "CTE";
  else if( num == 1927) 
    return "CTG";
  else if( num == 1928) 
    return "CTH";
  else if( num == 1929) 
    return "CTI";
  else if( num == 1930) 
    return "CTL";
  else if( num == 1931) 
    return "CTK";
  else if( num == 1932) 
    return "CTM";
  else if( num == 1933) 
    return "CTF";
  else if( num == 1934) 
    return "CTP";
  else if( num == 1935) 
    return "CTS";
  else if( num == 1936) 
    return "CTT";
  else if( num == 1937) 
    return "CTW";
  else if( num == 1938) 
    return "CTY";
  else if( num == 1939) 
    return "CTV";
  else if( num == 1940) 
    return "CWA";
  else if( num == 1941) 
    return "CWR";
  else if( num == 1942) 
    return "CWN";
  else if( num == 1943) 
    return "CWD";
  else if( num == 1944) 
    return "CWC";
  else if( num == 1945) 
    return "CWQ";
  else if( num == 1946) 
    return "CWE";
  else if( num == 1947) 
    return "CWG";
  else if( num == 1948) 
    return "CWH";
  else if( num == 1949) 
    return "CWI";
  else if( num == 1950) 
    return "CWL";
  else if( num == 1951) 
    return "CWK";
  else if( num == 1952) 
    return "CWM";
  else if( num == 1953) 
    return "CWF";
  else if( num == 1954) 
    return "CWP";
  else if( num == 1955) 
    return "CWS";
  else if( num == 1956) 
    return "CWT";
  else if( num == 1957) 
    return "CWW";
  else if( num == 1958) 
    return "CWY";
  else if( num == 1959) 
    return "CWV";
  else if( num == 1960) 
    return "CYA";
  else if( num == 1961) 
    return "CYR";
  else if( num == 1962) 
    return "CYN";
  else if( num == 1963) 
    return "CYD";
  else if( num == 1964) 
    return "CYC";
  else if( num == 1965) 
    return "CYQ";
  else if( num == 1966) 
    return "CYE";
  else if( num == 1967) 
    return "CYG";
  else if( num == 1968) 
    return "CYH";
  else if( num == 1969) 
    return "CYI";
  else if( num == 1970) 
    return "CYL";
  else if( num == 1971) 
    return "CYK";
  else if( num == 1972) 
    return "CYM";
  else if( num == 1973) 
    return "CYF";
  else if( num == 1974) 
    return "CYP";
  else if( num == 1975) 
    return "CYS";
  else if( num == 1976) 
    return "CYT";
  else if( num == 1977) 
    return "CYW";
  else if( num == 1978) 
    return "CYY";
  else if( num == 1979) 
    return "CYV";
  else if( num == 1980) 
    return "CVA";
  else if( num == 1981) 
    return "CVR";
  else if( num == 1982) 
    return "CVN";
  else if( num == 1983) 
    return "CVD";
  else if( num == 1984) 
    return "CVC";
  else if( num == 1985) 
    return "CVQ";
  else if( num == 1986) 
    return "CVE";
  else if( num == 1987) 
    return "CVG";
  else if( num == 1988) 
    return "CVH";
  else if( num == 1989) 
    return "CVI";
  else if( num == 1990) 
    return "CVL";
  else if( num == 1991) 
    return "CVK";
  else if( num == 1992) 
    return "CVM";
  else if( num == 1993) 
    return "CVF";
  else if( num == 1994) 
    return "CVP";
  else if( num == 1995) 
    return "CVS";
  else if( num == 1996) 
    return "CVT";
  else if( num == 1997) 
    return "CVW";
  else if( num == 1998) 
    return "CVY";
  else if( num == 1999) 
    return "CVV";
  else if( num == 2000) 
    return "QAA";
  else if( num == 2001) 
    return "QAR";
  else if( num == 2002) 
    return "QAN";
  else if( num == 2003) 
    return "QAD";
  else if( num == 2004) 
    return "QAC";
  else if( num == 2005) 
    return "QAQ";
  else if( num == 2006) 
    return "QAE";
  else if( num == 2007) 
    return "QAG";
  else if( num == 2008) 
    return "QAH";
  else if( num == 2009) 
    return "QAI";
  else if( num == 2010) 
    return "QAL";
  else if( num == 2011) 
    return "QAK";
  else if( num == 2012) 
    return "QAM";
  else if( num == 2013) 
    return "QAF";
  else if( num == 2014) 
    return "QAP";
  else if( num == 2015) 
    return "QAS";
  else if( num == 2016) 
    return "QAT";
  else if( num == 2017) 
    return "QAW";
  else if( num == 2018) 
    return "QAY";
  else if( num == 2019) 
    return "QAV";
  else if( num == 2020) 
    return "QRA";
  else if( num == 2021) 
    return "QRR";
  else if( num == 2022) 
    return "QRN";
  else if( num == 2023) 
    return "QRD";
  else if( num == 2024) 
    return "QRC";
  else if( num == 2025) 
    return "QRQ";
  else if( num == 2026) 
    return "QRE";
  else if( num == 2027) 
    return "QRG";
  else if( num == 2028) 
    return "QRH";
  else if( num == 2029) 
    return "QRI";
  else if( num == 2030) 
    return "QRL";
  else if( num == 2031) 
    return "QRK";
  else if( num == 2032) 
    return "QRM";
  else if( num == 2033) 
    return "QRF";
  else if( num == 2034) 
    return "QRP";
  else if( num == 2035) 
    return "QRS";
  else if( num == 2036) 
    return "QRT";
  else if( num == 2037) 
    return "QRW";
  else if( num == 2038) 
    return "QRY";
  else if( num == 2039) 
    return "QRV";
  else if( num == 2040) 
    return "QNA";
  else if( num == 2041) 
    return "QNR";
  else if( num == 2042) 
    return "QNN";
  else if( num == 2043) 
    return "QND";
  else if( num == 2044) 
    return "QNC";
  else if( num == 2045) 
    return "QNQ";
  else if( num == 2046) 
    return "QNE";
  else if( num == 2047) 
    return "QNG";
  else if( num == 2048) 
    return "QNH";
  else if( num == 2049) 
    return "QNI";
  else if( num == 2050) 
    return "QNL";
  else if( num == 2051) 
    return "QNK";
  else if( num == 2052) 
    return "QNM";
  else if( num == 2053) 
    return "QNF";
  else if( num == 2054) 
    return "QNP";
  else if( num == 2055) 
    return "QNS";
  else if( num == 2056) 
    return "QNT";
  else if( num == 2057) 
    return "QNW";
  else if( num == 2058) 
    return "QNY";
  else if( num == 2059) 
    return "QNV";
  else if( num == 2060) 
    return "QDA";
  else if( num == 2061) 
    return "QDR";
  else if( num == 2062) 
    return "QDN";
  else if( num == 2063) 
    return "QDD";
  else if( num == 2064) 
    return "QDC";
  else if( num == 2065) 
    return "QDQ";
  else if( num == 2066) 
    return "QDE";
  else if( num == 2067) 
    return "QDG";
  else if( num == 2068) 
    return "QDH";
  else if( num == 2069) 
    return "QDI";
  else if( num == 2070) 
    return "QDL";
  else if( num == 2071) 
    return "QDK";
  else if( num == 2072) 
    return "QDM";
  else if( num == 2073) 
    return "QDF";
  else if( num == 2074) 
    return "QDP";
  else if( num == 2075) 
    return "QDS";
  else if( num == 2076) 
    return "QDT";
  else if( num == 2077) 
    return "QDW";
  else if( num == 2078) 
    return "QDY";
  else if( num == 2079) 
    return "QDV";
  else if( num == 2080) 
    return "QCA";
  else if( num == 2081) 
    return "QCR";
  else if( num == 2082) 
    return "QCN";
  else if( num == 2083) 
    return "QCD";
  else if( num == 2084) 
    return "QCC";
  else if( num == 2085) 
    return "QCQ";
  else if( num == 2086) 
    return "QCE";
  else if( num == 2087) 
    return "QCG";
  else if( num == 2088) 
    return "QCH";
  else if( num == 2089) 
    return "QCI";
  else if( num == 2090) 
    return "QCL";
  else if( num == 2091) 
    return "QCK";
  else if( num == 2092) 
    return "QCM";
  else if( num == 2093) 
    return "QCF";
  else if( num == 2094) 
    return "QCP";
  else if( num == 2095) 
    return "QCS";
  else if( num == 2096) 
    return "QCT";
  else if( num == 2097) 
    return "QCW";
  else if( num == 2098) 
    return "QCY";
  else if( num == 2099) 
    return "QCV";
  else if( num == 2100) 
    return "QQA";
  else if( num == 2101) 
    return "QQR";
  else if( num == 2102) 
    return "QQN";
  else if( num == 2103) 
    return "QQD";
  else if( num == 2104) 
    return "QQC";
  else if( num == 2105) 
    return "QQQ";
  else if( num == 2106) 
    return "QQE";
  else if( num == 2107) 
    return "QQG";
  else if( num == 2108) 
    return "QQH";
  else if( num == 2109) 
    return "QQI";
  else if( num == 2110) 
    return "QQL";
  else if( num == 2111) 
    return "QQK";
  else if( num == 2112) 
    return "QQM";
  else if( num == 2113) 
    return "QQF";
  else if( num == 2114) 
    return "QQP";
  else if( num == 2115) 
    return "QQS";
  else if( num == 2116) 
    return "QQT";
  else if( num == 2117) 
    return "QQW";
  else if( num == 2118) 
    return "QQY";
  else if( num == 2119) 
    return "QQV";
  else if( num == 2120) 
    return "QEA";
  else if( num == 2121) 
    return "QER";
  else if( num == 2122) 
    return "QEN";
  else if( num == 2123) 
    return "QED";
  else if( num == 2124) 
    return "QEC";
  else if( num == 2125) 
    return "QEQ";
  else if( num == 2126) 
    return "QEE";
  else if( num == 2127) 
    return "QEG";
  else if( num == 2128) 
    return "QEH";
  else if( num == 2129) 
    return "QEI";
  else if( num == 2130) 
    return "QEL";
  else if( num == 2131) 
    return "QEK";
  else if( num == 2132) 
    return "QEM";
  else if( num == 2133) 
    return "QEF";
  else if( num == 2134) 
    return "QEP";
  else if( num == 2135) 
    return "QES";
  else if( num == 2136) 
    return "QET";
  else if( num == 2137) 
    return "QEW";
  else if( num == 2138) 
    return "QEY";
  else if( num == 2139) 
    return "QEV";
  else if( num == 2140) 
    return "QGA";
  else if( num == 2141) 
    return "QGR";
  else if( num == 2142) 
    return "QGN";
  else if( num == 2143) 
    return "QGD";
  else if( num == 2144) 
    return "QGC";
  else if( num == 2145) 
    return "QGQ";
  else if( num == 2146) 
    return "QGE";
  else if( num == 2147) 
    return "QGG";
  else if( num == 2148) 
    return "QGH";
  else if( num == 2149) 
    return "QGI";
  else if( num == 2150) 
    return "QGL";
  else if( num == 2151) 
    return "QGK";
  else if( num == 2152) 
    return "QGM";
  else if( num == 2153) 
    return "QGF";
  else if( num == 2154) 
    return "QGP";
  else if( num == 2155) 
    return "QGS";
  else if( num == 2156) 
    return "QGT";
  else if( num == 2157) 
    return "QGW";
  else if( num == 2158) 
    return "QGY";
  else if( num == 2159) 
    return "QGV";
  else if( num == 2160) 
    return "QHA";
  else if( num == 2161) 
    return "QHR";
  else if( num == 2162) 
    return "QHN";
  else if( num == 2163) 
    return "QHD";
  else if( num == 2164) 
    return "QHC";
  else if( num == 2165) 
    return "QHQ";
  else if( num == 2166) 
    return "QHE";
  else if( num == 2167) 
    return "QHG";
  else if( num == 2168) 
    return "QHH";
  else if( num == 2169) 
    return "QHI";
  else if( num == 2170) 
    return "QHL";
  else if( num == 2171) 
    return "QHK";
  else if( num == 2172) 
    return "QHM";
  else if( num == 2173) 
    return "QHF";
  else if( num == 2174) 
    return "QHP";
  else if( num == 2175) 
    return "QHS";
  else if( num == 2176) 
    return "QHT";
  else if( num == 2177) 
    return "QHW";
  else if( num == 2178) 
    return "QHY";
  else if( num == 2179) 
    return "QHV";
  else if( num == 2180) 
    return "QIA";
  else if( num == 2181) 
    return "QIR";
  else if( num == 2182) 
    return "QIN";
  else if( num == 2183) 
    return "QID";
  else if( num == 2184) 
    return "QIC";
  else if( num == 2185) 
    return "QIQ";
  else if( num == 2186) 
    return "QIE";
  else if( num == 2187) 
    return "QIG";
  else if( num == 2188) 
    return "QIH";
  else if( num == 2189) 
    return "QII";
  else if( num == 2190) 
    return "QIL";
  else if( num == 2191) 
    return "QIK";
  else if( num == 2192) 
    return "QIM";
  else if( num == 2193) 
    return "QIF";
  else if( num == 2194) 
    return "QIP";
  else if( num == 2195) 
    return "QIS";
  else if( num == 2196) 
    return "QIT";
  else if( num == 2197) 
    return "QIW";
  else if( num == 2198) 
    return "QIY";
  else if( num == 2199) 
    return "QIV";
  else if( num == 2200) 
    return "QLA";
  else if( num == 2201) 
    return "QLR";
  else if( num == 2202) 
    return "QLN";
  else if( num == 2203) 
    return "QLD";
  else if( num == 2204) 
    return "QLC";
  else if( num == 2205) 
    return "QLQ";
  else if( num == 2206) 
    return "QLE";
  else if( num == 2207) 
    return "QLG";
  else if( num == 2208) 
    return "QLH";
  else if( num == 2209) 
    return "QLI";
  else if( num == 2210) 
    return "QLL";
  else if( num == 2211) 
    return "QLK";
  else if( num == 2212) 
    return "QLM";
  else if( num == 2213) 
    return "QLF";
  else if( num == 2214) 
    return "QLP";
  else if( num == 2215) 
    return "QLS";
  else if( num == 2216) 
    return "QLT";
  else if( num == 2217) 
    return "QLW";
  else if( num == 2218) 
    return "QLY";
  else if( num == 2219) 
    return "QLV";
  else if( num == 2220) 
    return "QKA";
  else if( num == 2221) 
    return "QKR";
  else if( num == 2222) 
    return "QKN";
  else if( num == 2223) 
    return "QKD";
  else if( num == 2224) 
    return "QKC";
  else if( num == 2225) 
    return "QKQ";
  else if( num == 2226) 
    return "QKE";
  else if( num == 2227) 
    return "QKG";
  else if( num == 2228) 
    return "QKH";
  else if( num == 2229) 
    return "QKI";
  else if( num == 2230) 
    return "QKL";
  else if( num == 2231) 
    return "QKK";
  else if( num == 2232) 
    return "QKM";
  else if( num == 2233) 
    return "QKF";
  else if( num == 2234) 
    return "QKP";
  else if( num == 2235) 
    return "QKS";
  else if( num == 2236) 
    return "QKT";
  else if( num == 2237) 
    return "QKW";
  else if( num == 2238) 
    return "QKY";
  else if( num == 2239) 
    return "QKV";
  else if( num == 2240) 
    return "QMA";
  else if( num == 2241) 
    return "QMR";
  else if( num == 2242) 
    return "QMN";
  else if( num == 2243) 
    return "QMD";
  else if( num == 2244) 
    return "QMC";
  else if( num == 2245) 
    return "QMQ";
  else if( num == 2246) 
    return "QME";
  else if( num == 2247) 
    return "QMG";
  else if( num == 2248) 
    return "QMH";
  else if( num == 2249) 
    return "QMI";
  else if( num == 2250) 
    return "QML";
  else if( num == 2251) 
    return "QMK";
  else if( num == 2252) 
    return "QMM";
  else if( num == 2253) 
    return "QMF";
  else if( num == 2254) 
    return "QMP";
  else if( num == 2255) 
    return "QMS";
  else if( num == 2256) 
    return "QMT";
  else if( num == 2257) 
    return "QMW";
  else if( num == 2258) 
    return "QMY";
  else if( num == 2259) 
    return "QMV";
  else if( num == 2260) 
    return "QFA";
  else if( num == 2261) 
    return "QFR";
  else if( num == 2262) 
    return "QFN";
  else if( num == 2263) 
    return "QFD";
  else if( num == 2264) 
    return "QFC";
  else if( num == 2265) 
    return "QFQ";
  else if( num == 2266) 
    return "QFE";
  else if( num == 2267) 
    return "QFG";
  else if( num == 2268) 
    return "QFH";
  else if( num == 2269) 
    return "QFI";
  else if( num == 2270) 
    return "QFL";
  else if( num == 2271) 
    return "QFK";
  else if( num == 2272) 
    return "QFM";
  else if( num == 2273) 
    return "QFF";
  else if( num == 2274) 
    return "QFP";
  else if( num == 2275) 
    return "QFS";
  else if( num == 2276) 
    return "QFT";
  else if( num == 2277) 
    return "QFW";
  else if( num == 2278) 
    return "QFY";
  else if( num == 2279) 
    return "QFV";
  else if( num == 2280) 
    return "QPA";
  else if( num == 2281) 
    return "QPR";
  else if( num == 2282) 
    return "QPN";
  else if( num == 2283) 
    return "QPD";
  else if( num == 2284) 
    return "QPC";
  else if( num == 2285) 
    return "QPQ";
  else if( num == 2286) 
    return "QPE";
  else if( num == 2287) 
    return "QPG";
  else if( num == 2288) 
    return "QPH";
  else if( num == 2289) 
    return "QPI";
  else if( num == 2290) 
    return "QPL";
  else if( num == 2291) 
    return "QPK";
  else if( num == 2292) 
    return "QPM";
  else if( num == 2293) 
    return "QPF";
  else if( num == 2294) 
    return "QPP";
  else if( num == 2295) 
    return "QPS";
  else if( num == 2296) 
    return "QPT";
  else if( num == 2297) 
    return "QPW";
  else if( num == 2298) 
    return "QPY";
  else if( num == 2299) 
    return "QPV";
  else if( num == 2300) 
    return "QSA";
  else if( num == 2301) 
    return "QSR";
  else if( num == 2302) 
    return "QSN";
  else if( num == 2303) 
    return "QSD";
  else if( num == 2304) 
    return "QSC";
  else if( num == 2305) 
    return "QSQ";
  else if( num == 2306) 
    return "QSE";
  else if( num == 2307) 
    return "QSG";
  else if( num == 2308) 
    return "QSH";
  else if( num == 2309) 
    return "QSI";
  else if( num == 2310) 
    return "QSL";
  else if( num == 2311) 
    return "QSK";
  else if( num == 2312) 
    return "QSM";
  else if( num == 2313) 
    return "QSF";
  else if( num == 2314) 
    return "QSP";
  else if( num == 2315) 
    return "QSS";
  else if( num == 2316) 
    return "QST";
  else if( num == 2317) 
    return "QSW";
  else if( num == 2318) 
    return "QSY";
  else if( num == 2319) 
    return "QSV";
  else if( num == 2320) 
    return "QTA";
  else if( num == 2321) 
    return "QTR";
  else if( num == 2322) 
    return "QTN";
  else if( num == 2323) 
    return "QTD";
  else if( num == 2324) 
    return "QTC";
  else if( num == 2325) 
    return "QTQ";
  else if( num == 2326) 
    return "QTE";
  else if( num == 2327) 
    return "QTG";
  else if( num == 2328) 
    return "QTH";
  else if( num == 2329) 
    return "QTI";
  else if( num == 2330) 
    return "QTL";
  else if( num == 2331) 
    return "QTK";
  else if( num == 2332) 
    return "QTM";
  else if( num == 2333) 
    return "QTF";
  else if( num == 2334) 
    return "QTP";
  else if( num == 2335) 
    return "QTS";
  else if( num == 2336) 
    return "QTT";
  else if( num == 2337) 
    return "QTW";
  else if( num == 2338) 
    return "QTY";
  else if( num == 2339) 
    return "QTV";
  else if( num == 2340) 
    return "QWA";
  else if( num == 2341) 
    return "QWR";
  else if( num == 2342) 
    return "QWN";
  else if( num == 2343) 
    return "QWD";
  else if( num == 2344) 
    return "QWC";
  else if( num == 2345) 
    return "QWQ";
  else if( num == 2346) 
    return "QWE";
  else if( num == 2347) 
    return "QWG";
  else if( num == 2348) 
    return "QWH";
  else if( num == 2349) 
    return "QWI";
  else if( num == 2350) 
    return "QWL";
  else if( num == 2351) 
    return "QWK";
  else if( num == 2352) 
    return "QWM";
  else if( num == 2353) 
    return "QWF";
  else if( num == 2354) 
    return "QWP";
  else if( num == 2355) 
    return "QWS";
  else if( num == 2356) 
    return "QWT";
  else if( num == 2357) 
    return "QWW";
  else if( num == 2358) 
    return "QWY";
  else if( num == 2359) 
    return "QWV";
  else if( num == 2360) 
    return "QYA";
  else if( num == 2361) 
    return "QYR";
  else if( num == 2362) 
    return "QYN";
  else if( num == 2363) 
    return "QYD";
  else if( num == 2364) 
    return "QYC";
  else if( num == 2365) 
    return "QYQ";
  else if( num == 2366) 
    return "QYE";
  else if( num == 2367) 
    return "QYG";
  else if( num == 2368) 
    return "QYH";
  else if( num == 2369) 
    return "QYI";
  else if( num == 2370) 
    return "QYL";
  else if( num == 2371) 
    return "QYK";
  else if( num == 2372) 
    return "QYM";
  else if( num == 2373) 
    return "QYF";
  else if( num == 2374) 
    return "QYP";
  else if( num == 2375) 
    return "QYS";
  else if( num == 2376) 
    return "QYT";
  else if( num == 2377) 
    return "QYW";
  else if( num == 2378) 
    return "QYY";
  else if( num == 2379) 
    return "QYV";
  else if( num == 2380) 
    return "QVA";
  else if( num == 2381) 
    return "QVR";
  else if( num == 2382) 
    return "QVN";
  else if( num == 2383) 
    return "QVD";
  else if( num == 2384) 
    return "QVC";
  else if( num == 2385) 
    return "QVQ";
  else if( num == 2386) 
    return "QVE";
  else if( num == 2387) 
    return "QVG";
  else if( num == 2388) 
    return "QVH";
  else if( num == 2389) 
    return "QVI";
  else if( num == 2390) 
    return "QVL";
  else if( num == 2391) 
    return "QVK";
  else if( num == 2392) 
    return "QVM";
  else if( num == 2393) 
    return "QVF";
  else if( num == 2394) 
    return "QVP";
  else if( num == 2395) 
    return "QVS";
  else if( num == 2396) 
    return "QVT";
  else if( num == 2397) 
    return "QVW";
  else if( num == 2398) 
    return "QVY";
  else if( num == 2399) 
    return "QVV";
  else if( num == 2400) 
    return "EAA";
  else if( num == 2401) 
    return "EAR";
  else if( num == 2402) 
    return "EAN";
  else if( num == 2403) 
    return "EAD";
  else if( num == 2404) 
    return "EAC";
  else if( num == 2405) 
    return "EAQ";
  else if( num == 2406) 
    return "EAE";
  else if( num == 2407) 
    return "EAG";
  else if( num == 2408) 
    return "EAH";
  else if( num == 2409) 
    return "EAI";
  else if( num == 2410) 
    return "EAL";
  else if( num == 2411) 
    return "EAK";
  else if( num == 2412) 
    return "EAM";
  else if( num == 2413) 
    return "EAF";
  else if( num == 2414) 
    return "EAP";
  else if( num == 2415) 
    return "EAS";
  else if( num == 2416) 
    return "EAT";
  else if( num == 2417) 
    return "EAW";
  else if( num == 2418) 
    return "EAY";
  else if( num == 2419) 
    return "EAV";
  else if( num == 2420) 
    return "ERA";
  else if( num == 2421) 
    return "ERR";
  else if( num == 2422) 
    return "ERN";
  else if( num == 2423) 
    return "ERD";
  else if( num == 2424) 
    return "ERC";
  else if( num == 2425) 
    return "ERQ";
  else if( num == 2426) 
    return "ERE";
  else if( num == 2427) 
    return "ERG";
  else if( num == 2428) 
    return "ERH";
  else if( num == 2429) 
    return "ERI";
  else if( num == 2430) 
    return "ERL";
  else if( num == 2431) 
    return "ERK";
  else if( num == 2432) 
    return "ERM";
  else if( num == 2433) 
    return "ERF";
  else if( num == 2434) 
    return "ERP";
  else if( num == 2435) 
    return "ERS";
  else if( num == 2436) 
    return "ERT";
  else if( num == 2437) 
    return "ERW";
  else if( num == 2438) 
    return "ERY";
  else if( num == 2439) 
    return "ERV";
  else if( num == 2440) 
    return "ENA";
  else if( num == 2441) 
    return "ENR";
  else if( num == 2442) 
    return "ENN";
  else if( num == 2443) 
    return "END";
  else if( num == 2444) 
    return "ENC";
  else if( num == 2445) 
    return "ENQ";
  else if( num == 2446) 
    return "ENE";
  else if( num == 2447) 
    return "ENG";
  else if( num == 2448) 
    return "ENH";
  else if( num == 2449) 
    return "ENI";
  else if( num == 2450) 
    return "ENL";
  else if( num == 2451) 
    return "ENK";
  else if( num == 2452) 
    return "ENM";
  else if( num == 2453) 
    return "ENF";
  else if( num == 2454) 
    return "ENP";
  else if( num == 2455) 
    return "ENS";
  else if( num == 2456) 
    return "ENT";
  else if( num == 2457) 
    return "ENW";
  else if( num == 2458) 
    return "ENY";
  else if( num == 2459) 
    return "ENV";
  else if( num == 2460) 
    return "EDA";
  else if( num == 2461) 
    return "EDR";
  else if( num == 2462) 
    return "EDN";
  else if( num == 2463) 
    return "EDD";
  else if( num == 2464) 
    return "EDC";
  else if( num == 2465) 
    return "EDQ";
  else if( num == 2466) 
    return "EDE";
  else if( num == 2467) 
    return "EDG";
  else if( num == 2468) 
    return "EDH";
  else if( num == 2469) 
    return "EDI";
  else if( num == 2470) 
    return "EDL";
  else if( num == 2471) 
    return "EDK";
  else if( num == 2472) 
    return "EDM";
  else if( num == 2473) 
    return "EDF";
  else if( num == 2474) 
    return "EDP";
  else if( num == 2475) 
    return "EDS";
  else if( num == 2476) 
    return "EDT";
  else if( num == 2477) 
    return "EDW";
  else if( num == 2478) 
    return "EDY";
  else if( num == 2479) 
    return "EDV";
  else if( num == 2480) 
    return "ECA";
  else if( num == 2481) 
    return "ECR";
  else if( num == 2482) 
    return "ECN";
  else if( num == 2483) 
    return "ECD";
  else if( num == 2484) 
    return "ECC";
  else if( num == 2485) 
    return "ECQ";
  else if( num == 2486) 
    return "ECE";
  else if( num == 2487) 
    return "ECG";
  else if( num == 2488) 
    return "ECH";
  else if( num == 2489) 
    return "ECI";
  else if( num == 2490) 
    return "ECL";
  else if( num == 2491) 
    return "ECK";
  else if( num == 2492) 
    return "ECM";
  else if( num == 2493) 
    return "ECF";
  else if( num == 2494) 
    return "ECP";
  else if( num == 2495) 
    return "ECS";
  else if( num == 2496) 
    return "ECT";
  else if( num == 2497) 
    return "ECW";
  else if( num == 2498) 
    return "ECY";
  else if( num == 2499) 
    return "ECV";
  else if( num == 2500) 
    return "EQA";
  else if( num == 2501) 
    return "EQR";
  else if( num == 2502) 
    return "EQN";
  else if( num == 2503) 
    return "EQD";
  else if( num == 2504) 
    return "EQC";
  else if( num == 2505) 
    return "EQQ";
  else if( num == 2506) 
    return "EQE";
  else if( num == 2507) 
    return "EQG";
  else if( num == 2508) 
    return "EQH";
  else if( num == 2509) 
    return "EQI";
  else if( num == 2510) 
    return "EQL";
  else if( num == 2511) 
    return "EQK";
  else if( num == 2512) 
    return "EQM";
  else if( num == 2513) 
    return "EQF";
  else if( num == 2514) 
    return "EQP";
  else if( num == 2515) 
    return "EQS";
  else if( num == 2516) 
    return "EQT";
  else if( num == 2517) 
    return "EQW";
  else if( num == 2518) 
    return "EQY";
  else if( num == 2519) 
    return "EQV";
  else if( num == 2520) 
    return "EEA";
  else if( num == 2521) 
    return "EER";
  else if( num == 2522) 
    return "EEN";
  else if( num == 2523) 
    return "EED";
  else if( num == 2524) 
    return "EEC";
  else if( num == 2525) 
    return "EEQ";
  else if( num == 2526) 
    return "EEE";
  else if( num == 2527) 
    return "EEG";
  else if( num == 2528) 
    return "EEH";
  else if( num == 2529) 
    return "EEI";
  else if( num == 2530) 
    return "EEL";
  else if( num == 2531) 
    return "EEK";
  else if( num == 2532) 
    return "EEM";
  else if( num == 2533) 
    return "EEF";
  else if( num == 2534) 
    return "EEP";
  else if( num == 2535) 
    return "EES";
  else if( num == 2536) 
    return "EET";
  else if( num == 2537) 
    return "EEW";
  else if( num == 2538) 
    return "EEY";
  else if( num == 2539) 
    return "EEV";
  else if( num == 2540) 
    return "EGA";
  else if( num == 2541) 
    return "EGR";
  else if( num == 2542) 
    return "EGN";
  else if( num == 2543) 
    return "EGD";
  else if( num == 2544) 
    return "EGC";
  else if( num == 2545) 
    return "EGQ";
  else if( num == 2546) 
    return "EGE";
  else if( num == 2547) 
    return "EGG";
  else if( num == 2548) 
    return "EGH";
  else if( num == 2549) 
    return "EGI";
  else if( num == 2550) 
    return "EGL";
  else if( num == 2551) 
    return "EGK";
  else if( num == 2552) 
    return "EGM";
  else if( num == 2553) 
    return "EGF";
  else if( num == 2554) 
    return "EGP";
  else if( num == 2555) 
    return "EGS";
  else if( num == 2556) 
    return "EGT";
  else if( num == 2557) 
    return "EGW";
  else if( num == 2558) 
    return "EGY";
  else if( num == 2559) 
    return "EGV";
  else if( num == 2560) 
    return "EHA";
  else if( num == 2561) 
    return "EHR";
  else if( num == 2562) 
    return "EHN";
  else if( num == 2563) 
    return "EHD";
  else if( num == 2564) 
    return "EHC";
  else if( num == 2565) 
    return "EHQ";
  else if( num == 2566) 
    return "EHE";
  else if( num == 2567) 
    return "EHG";
  else if( num == 2568) 
    return "EHH";
  else if( num == 2569) 
    return "EHI";
  else if( num == 2570) 
    return "EHL";
  else if( num == 2571) 
    return "EHK";
  else if( num == 2572) 
    return "EHM";
  else if( num == 2573) 
    return "EHF";
  else if( num == 2574) 
    return "EHP";
  else if( num == 2575) 
    return "EHS";
  else if( num == 2576) 
    return "EHT";
  else if( num == 2577) 
    return "EHW";
  else if( num == 2578) 
    return "EHY";
  else if( num == 2579) 
    return "EHV";
  else if( num == 2580) 
    return "EIA";
  else if( num == 2581) 
    return "EIR";
  else if( num == 2582) 
    return "EIN";
  else if( num == 2583) 
    return "EID";
  else if( num == 2584) 
    return "EIC";
  else if( num == 2585) 
    return "EIQ";
  else if( num == 2586) 
    return "EIE";
  else if( num == 2587) 
    return "EIG";
  else if( num == 2588) 
    return "EIH";
  else if( num == 2589) 
    return "EII";
  else if( num == 2590) 
    return "EIL";
  else if( num == 2591) 
    return "EIK";
  else if( num == 2592) 
    return "EIM";
  else if( num == 2593) 
    return "EIF";
  else if( num == 2594) 
    return "EIP";
  else if( num == 2595) 
    return "EIS";
  else if( num == 2596) 
    return "EIT";
  else if( num == 2597) 
    return "EIW";
  else if( num == 2598) 
    return "EIY";
  else if( num == 2599) 
    return "EIV";
  else if( num == 2600) 
    return "ELA";
  else if( num == 2601) 
    return "ELR";
  else if( num == 2602) 
    return "ELN";
  else if( num == 2603) 
    return "ELD";
  else if( num == 2604) 
    return "ELC";
  else if( num == 2605) 
    return "ELQ";
  else if( num == 2606) 
    return "ELE";
  else if( num == 2607) 
    return "ELG";
  else if( num == 2608) 
    return "ELH";
  else if( num == 2609) 
    return "ELI";
  else if( num == 2610) 
    return "ELL";
  else if( num == 2611) 
    return "ELK";
  else if( num == 2612) 
    return "ELM";
  else if( num == 2613) 
    return "ELF";
  else if( num == 2614) 
    return "ELP";
  else if( num == 2615) 
    return "ELS";
  else if( num == 2616) 
    return "ELT";
  else if( num == 2617) 
    return "ELW";
  else if( num == 2618) 
    return "ELY";
  else if( num == 2619) 
    return "ELV";
  else if( num == 2620) 
    return "EKA";
  else if( num == 2621) 
    return "EKR";
  else if( num == 2622) 
    return "EKN";
  else if( num == 2623) 
    return "EKD";
  else if( num == 2624) 
    return "EKC";
  else if( num == 2625) 
    return "EKQ";
  else if( num == 2626) 
    return "EKE";
  else if( num == 2627) 
    return "EKG";
  else if( num == 2628) 
    return "EKH";
  else if( num == 2629) 
    return "EKI";
  else if( num == 2630) 
    return "EKL";
  else if( num == 2631) 
    return "EKK";
  else if( num == 2632) 
    return "EKM";
  else if( num == 2633) 
    return "EKF";
  else if( num == 2634) 
    return "EKP";
  else if( num == 2635) 
    return "EKS";
  else if( num == 2636) 
    return "EKT";
  else if( num == 2637) 
    return "EKW";
  else if( num == 2638) 
    return "EKY";
  else if( num == 2639) 
    return "EKV";
  else if( num == 2640) 
    return "EMA";
  else if( num == 2641) 
    return "EMR";
  else if( num == 2642) 
    return "EMN";
  else if( num == 2643) 
    return "EMD";
  else if( num == 2644) 
    return "EMC";
  else if( num == 2645) 
    return "EMQ";
  else if( num == 2646) 
    return "EME";
  else if( num == 2647) 
    return "EMG";
  else if( num == 2648) 
    return "EMH";
  else if( num == 2649) 
    return "EMI";
  else if( num == 2650) 
    return "EML";
  else if( num == 2651) 
    return "EMK";
  else if( num == 2652) 
    return "EMM";
  else if( num == 2653) 
    return "EMF";
  else if( num == 2654) 
    return "EMP";
  else if( num == 2655) 
    return "EMS";
  else if( num == 2656) 
    return "EMT";
  else if( num == 2657) 
    return "EMW";
  else if( num == 2658) 
    return "EMY";
  else if( num == 2659) 
    return "EMV";
  else if( num == 2660) 
    return "EFA";
  else if( num == 2661) 
    return "EFR";
  else if( num == 2662) 
    return "EFN";
  else if( num == 2663) 
    return "EFD";
  else if( num == 2664) 
    return "EFC";
  else if( num == 2665) 
    return "EFQ";
  else if( num == 2666) 
    return "EFE";
  else if( num == 2667) 
    return "EFG";
  else if( num == 2668) 
    return "EFH";
  else if( num == 2669) 
    return "EFI";
  else if( num == 2670) 
    return "EFL";
  else if( num == 2671) 
    return "EFK";
  else if( num == 2672) 
    return "EFM";
  else if( num == 2673) 
    return "EFF";
  else if( num == 2674) 
    return "EFP";
  else if( num == 2675) 
    return "EFS";
  else if( num == 2676) 
    return "EFT";
  else if( num == 2677) 
    return "EFW";
  else if( num == 2678) 
    return "EFY";
  else if( num == 2679) 
    return "EFV";
  else if( num == 2680) 
    return "EPA";
  else if( num == 2681) 
    return "EPR";
  else if( num == 2682) 
    return "EPN";
  else if( num == 2683) 
    return "EPD";
  else if( num == 2684) 
    return "EPC";
  else if( num == 2685) 
    return "EPQ";
  else if( num == 2686) 
    return "EPE";
  else if( num == 2687) 
    return "EPG";
  else if( num == 2688) 
    return "EPH";
  else if( num == 2689) 
    return "EPI";
  else if( num == 2690) 
    return "EPL";
  else if( num == 2691) 
    return "EPK";
  else if( num == 2692) 
    return "EPM";
  else if( num == 2693) 
    return "EPF";
  else if( num == 2694) 
    return "EPP";
  else if( num == 2695) 
    return "EPS";
  else if( num == 2696) 
    return "EPT";
  else if( num == 2697) 
    return "EPW";
  else if( num == 2698) 
    return "EPY";
  else if( num == 2699) 
    return "EPV";
  else if( num == 2700) 
    return "ESA";
  else if( num == 2701) 
    return "ESR";
  else if( num == 2702) 
    return "ESN";
  else if( num == 2703) 
    return "ESD";
  else if( num == 2704) 
    return "ESC";
  else if( num == 2705) 
    return "ESQ";
  else if( num == 2706) 
    return "ESE";
  else if( num == 2707) 
    return "ESG";
  else if( num == 2708) 
    return "ESH";
  else if( num == 2709) 
    return "ESI";
  else if( num == 2710) 
    return "ESL";
  else if( num == 2711) 
    return "ESK";
  else if( num == 2712) 
    return "ESM";
  else if( num == 2713) 
    return "ESF";
  else if( num == 2714) 
    return "ESP";
  else if( num == 2715) 
    return "ESS";
  else if( num == 2716) 
    return "EST";
  else if( num == 2717) 
    return "ESW";
  else if( num == 2718) 
    return "ESY";
  else if( num == 2719) 
    return "ESV";
  else if( num == 2720) 
    return "ETA";
  else if( num == 2721) 
    return "ETR";
  else if( num == 2722) 
    return "ETN";
  else if( num == 2723) 
    return "ETD";
  else if( num == 2724) 
    return "ETC";
  else if( num == 2725) 
    return "ETQ";
  else if( num == 2726) 
    return "ETE";
  else if( num == 2727) 
    return "ETG";
  else if( num == 2728) 
    return "ETH";
  else if( num == 2729) 
    return "ETI";
  else if( num == 2730) 
    return "ETL";
  else if( num == 2731) 
    return "ETK";
  else if( num == 2732) 
    return "ETM";
  else if( num == 2733) 
    return "ETF";
  else if( num == 2734) 
    return "ETP";
  else if( num == 2735) 
    return "ETS";
  else if( num == 2736) 
    return "ETT";
  else if( num == 2737) 
    return "ETW";
  else if( num == 2738) 
    return "ETY";
  else if( num == 2739) 
    return "ETV";
  else if( num == 2740) 
    return "EWA";
  else if( num == 2741) 
    return "EWR";
  else if( num == 2742) 
    return "EWN";
  else if( num == 2743) 
    return "EWD";
  else if( num == 2744) 
    return "EWC";
  else if( num == 2745) 
    return "EWQ";
  else if( num == 2746) 
    return "EWE";
  else if( num == 2747) 
    return "EWG";
  else if( num == 2748) 
    return "EWH";
  else if( num == 2749) 
    return "EWI";
  else if( num == 2750) 
    return "EWL";
  else if( num == 2751) 
    return "EWK";
  else if( num == 2752) 
    return "EWM";
  else if( num == 2753) 
    return "EWF";
  else if( num == 2754) 
    return "EWP";
  else if( num == 2755) 
    return "EWS";
  else if( num == 2756) 
    return "EWT";
  else if( num == 2757) 
    return "EWW";
  else if( num == 2758) 
    return "EWY";
  else if( num == 2759) 
    return "EWV";
  else if( num == 2760) 
    return "EYA";
  else if( num == 2761) 
    return "EYR";
  else if( num == 2762) 
    return "EYN";
  else if( num == 2763) 
    return "EYD";
  else if( num == 2764) 
    return "EYC";
  else if( num == 2765) 
    return "EYQ";
  else if( num == 2766) 
    return "EYE";
  else if( num == 2767) 
    return "EYG";
  else if( num == 2768) 
    return "EYH";
  else if( num == 2769) 
    return "EYI";
  else if( num == 2770) 
    return "EYL";
  else if( num == 2771) 
    return "EYK";
  else if( num == 2772) 
    return "EYM";
  else if( num == 2773) 
    return "EYF";
  else if( num == 2774) 
    return "EYP";
  else if( num == 2775) 
    return "EYS";
  else if( num == 2776) 
    return "EYT";
  else if( num == 2777) 
    return "EYW";
  else if( num == 2778) 
    return "EYY";
  else if( num == 2779) 
    return "EYV";
  else if( num == 2780) 
    return "EVA";
  else if( num == 2781) 
    return "EVR";
  else if( num == 2782) 
    return "EVN";
  else if( num == 2783) 
    return "EVD";
  else if( num == 2784) 
    return "EVC";
  else if( num == 2785) 
    return "EVQ";
  else if( num == 2786) 
    return "EVE";
  else if( num == 2787) 
    return "EVG";
  else if( num == 2788) 
    return "EVH";
  else if( num == 2789) 
    return "EVI";
  else if( num == 2790) 
    return "EVL";
  else if( num == 2791) 
    return "EVK";
  else if( num == 2792) 
    return "EVM";
  else if( num == 2793) 
    return "EVF";
  else if( num == 2794) 
    return "EVP";
  else if( num == 2795) 
    return "EVS";
  else if( num == 2796) 
    return "EVT";
  else if( num == 2797) 
    return "EVW";
  else if( num == 2798) 
    return "EVY";
  else if( num == 2799) 
    return "EVV";
  else if( num == 2800) 
    return "GAA";
  else if( num == 2801) 
    return "GAR";
  else if( num == 2802) 
    return "GAN";
  else if( num == 2803) 
    return "GAD";
  else if( num == 2804) 
    return "GAC";
  else if( num == 2805) 
    return "GAQ";
  else if( num == 2806) 
    return "GAE";
  else if( num == 2807) 
    return "GAG";
  else if( num == 2808) 
    return "GAH";
  else if( num == 2809) 
    return "GAI";
  else if( num == 2810) 
    return "GAL";
  else if( num == 2811) 
    return "GAK";
  else if( num == 2812) 
    return "GAM";
  else if( num == 2813) 
    return "GAF";
  else if( num == 2814) 
    return "GAP";
  else if( num == 2815) 
    return "GAS";
  else if( num == 2816) 
    return "GAT";
  else if( num == 2817) 
    return "GAW";
  else if( num == 2818) 
    return "GAY";
  else if( num == 2819) 
    return "GAV";
  else if( num == 2820) 
    return "GRA";
  else if( num == 2821) 
    return "GRR";
  else if( num == 2822) 
    return "GRN";
  else if( num == 2823) 
    return "GRD";
  else if( num == 2824) 
    return "GRC";
  else if( num == 2825) 
    return "GRQ";
  else if( num == 2826) 
    return "GRE";
  else if( num == 2827) 
    return "GRG";
  else if( num == 2828) 
    return "GRH";
  else if( num == 2829) 
    return "GRI";
  else if( num == 2830) 
    return "GRL";
  else if( num == 2831) 
    return "GRK";
  else if( num == 2832) 
    return "GRM";
  else if( num == 2833) 
    return "GRF";
  else if( num == 2834) 
    return "GRP";
  else if( num == 2835) 
    return "GRS";
  else if( num == 2836) 
    return "GRT";
  else if( num == 2837) 
    return "GRW";
  else if( num == 2838) 
    return "GRY";
  else if( num == 2839) 
    return "GRV";
  else if( num == 2840) 
    return "GNA";
  else if( num == 2841) 
    return "GNR";
  else if( num == 2842) 
    return "GNN";
  else if( num == 2843) 
    return "GND";
  else if( num == 2844) 
    return "GNC";
  else if( num == 2845) 
    return "GNQ";
  else if( num == 2846) 
    return "GNE";
  else if( num == 2847) 
    return "GNG";
  else if( num == 2848) 
    return "GNH";
  else if( num == 2849) 
    return "GNI";
  else if( num == 2850) 
    return "GNL";
  else if( num == 2851) 
    return "GNK";
  else if( num == 2852) 
    return "GNM";
  else if( num == 2853) 
    return "GNF";
  else if( num == 2854) 
    return "GNP";
  else if( num == 2855) 
    return "GNS";
  else if( num == 2856) 
    return "GNT";
  else if( num == 2857) 
    return "GNW";
  else if( num == 2858) 
    return "GNY";
  else if( num == 2859) 
    return "GNV";
  else if( num == 2860) 
    return "GDA";
  else if( num == 2861) 
    return "GDR";
  else if( num == 2862) 
    return "GDN";
  else if( num == 2863) 
    return "GDD";
  else if( num == 2864) 
    return "GDC";
  else if( num == 2865) 
    return "GDQ";
  else if( num == 2866) 
    return "GDE";
  else if( num == 2867) 
    return "GDG";
  else if( num == 2868) 
    return "GDH";
  else if( num == 2869) 
    return "GDI";
  else if( num == 2870) 
    return "GDL";
  else if( num == 2871) 
    return "GDK";
  else if( num == 2872) 
    return "GDM";
  else if( num == 2873) 
    return "GDF";
  else if( num == 2874) 
    return "GDP";
  else if( num == 2875) 
    return "GDS";
  else if( num == 2876) 
    return "GDT";
  else if( num == 2877) 
    return "GDW";
  else if( num == 2878) 
    return "GDY";
  else if( num == 2879) 
    return "GDV";
  else if( num == 2880) 
    return "GCA";
  else if( num == 2881) 
    return "GCR";
  else if( num == 2882) 
    return "GCN";
  else if( num == 2883) 
    return "GCD";
  else if( num == 2884) 
    return "GCC";
  else if( num == 2885) 
    return "GCQ";
  else if( num == 2886) 
    return "GCE";
  else if( num == 2887) 
    return "GCG";
  else if( num == 2888) 
    return "GCH";
  else if( num == 2889) 
    return "GCI";
  else if( num == 2890) 
    return "GCL";
  else if( num == 2891) 
    return "GCK";
  else if( num == 2892) 
    return "GCM";
  else if( num == 2893) 
    return "GCF";
  else if( num == 2894) 
    return "GCP";
  else if( num == 2895) 
    return "GCS";
  else if( num == 2896) 
    return "GCT";
  else if( num == 2897) 
    return "GCW";
  else if( num == 2898) 
    return "GCY";
  else if( num == 2899) 
    return "GCV";
  else if( num == 2900) 
    return "GQA";
  else if( num == 2901) 
    return "GQR";
  else if( num == 2902) 
    return "GQN";
  else if( num == 2903) 
    return "GQD";
  else if( num == 2904) 
    return "GQC";
  else if( num == 2905) 
    return "GQQ";
  else if( num == 2906) 
    return "GQE";
  else if( num == 2907) 
    return "GQG";
  else if( num == 2908) 
    return "GQH";
  else if( num == 2909) 
    return "GQI";
  else if( num == 2910) 
    return "GQL";
  else if( num == 2911) 
    return "GQK";
  else if( num == 2912) 
    return "GQM";
  else if( num == 2913) 
    return "GQF";
  else if( num == 2914) 
    return "GQP";
  else if( num == 2915) 
    return "GQS";
  else if( num == 2916) 
    return "GQT";
  else if( num == 2917) 
    return "GQW";
  else if( num == 2918) 
    return "GQY";
  else if( num == 2919) 
    return "GQV";
  else if( num == 2920) 
    return "GEA";
  else if( num == 2921) 
    return "GER";
  else if( num == 2922) 
    return "GEN";
  else if( num == 2923) 
    return "GED";
  else if( num == 2924) 
    return "GEC";
  else if( num == 2925) 
    return "GEQ";
  else if( num == 2926) 
    return "GEE";
  else if( num == 2927) 
    return "GEG";
  else if( num == 2928) 
    return "GEH";
  else if( num == 2929) 
    return "GEI";
  else if( num == 2930) 
    return "GEL";
  else if( num == 2931) 
    return "GEK";
  else if( num == 2932) 
    return "GEM";
  else if( num == 2933) 
    return "GEF";
  else if( num == 2934) 
    return "GEP";
  else if( num == 2935) 
    return "GES";
  else if( num == 2936) 
    return "GET";
  else if( num == 2937) 
    return "GEW";
  else if( num == 2938) 
    return "GEY";
  else if( num == 2939) 
    return "GEV";
  else if( num == 2940) 
    return "GGA";
  else if( num == 2941) 
    return "GGR";
  else if( num == 2942) 
    return "GGN";
  else if( num == 2943) 
    return "GGD";
  else if( num == 2944) 
    return "GGC";
  else if( num == 2945) 
    return "GGQ";
  else if( num == 2946) 
    return "GGE";
  else if( num == 2947) 
    return "GGG";
  else if( num == 2948) 
    return "GGH";
  else if( num == 2949) 
    return "GGI";
  else if( num == 2950) 
    return "GGL";
  else if( num == 2951) 
    return "GGK";
  else if( num == 2952) 
    return "GGM";
  else if( num == 2953) 
    return "GGF";
  else if( num == 2954) 
    return "GGP";
  else if( num == 2955) 
    return "GGS";
  else if( num == 2956) 
    return "GGT";
  else if( num == 2957) 
    return "GGW";
  else if( num == 2958) 
    return "GGY";
  else if( num == 2959) 
    return "GGV";
  else if( num == 2960) 
    return "GHA";
  else if( num == 2961) 
    return "GHR";
  else if( num == 2962) 
    return "GHN";
  else if( num == 2963) 
    return "GHD";
  else if( num == 2964) 
    return "GHC";
  else if( num == 2965) 
    return "GHQ";
  else if( num == 2966) 
    return "GHE";
  else if( num == 2967) 
    return "GHG";
  else if( num == 2968) 
    return "GHH";
  else if( num == 2969) 
    return "GHI";
  else if( num == 2970) 
    return "GHL";
  else if( num == 2971) 
    return "GHK";
  else if( num == 2972) 
    return "GHM";
  else if( num == 2973) 
    return "GHF";
  else if( num == 2974) 
    return "GHP";
  else if( num == 2975) 
    return "GHS";
  else if( num == 2976) 
    return "GHT";
  else if( num == 2977) 
    return "GHW";
  else if( num == 2978) 
    return "GHY";
  else if( num == 2979) 
    return "GHV";
  else if( num == 2980) 
    return "GIA";
  else if( num == 2981) 
    return "GIR";
  else if( num == 2982) 
    return "GIN";
  else if( num == 2983) 
    return "GID";
  else if( num == 2984) 
    return "GIC";
  else if( num == 2985) 
    return "GIQ";
  else if( num == 2986) 
    return "GIE";
  else if( num == 2987) 
    return "GIG";
  else if( num == 2988) 
    return "GIH";
  else if( num == 2989) 
    return "GII";
  else if( num == 2990) 
    return "GIL";
  else if( num == 2991) 
    return "GIK";
  else if( num == 2992) 
    return "GIM";
  else if( num == 2993) 
    return "GIF";
  else if( num == 2994) 
    return "GIP";
  else if( num == 2995) 
    return "GIS";
  else if( num == 2996) 
    return "GIT";
  else if( num == 2997) 
    return "GIW";
  else if( num == 2998) 
    return "GIY";
  else if( num == 2999) 
    return "GIV";
  else if( num == 3000) 
    return "GLA";
  else if( num == 3001) 
    return "GLR";
  else if( num == 3002) 
    return "GLN";
  else if( num == 3003) 
    return "GLD";
  else if( num == 3004) 
    return "GLC";
  else if( num == 3005) 
    return "GLQ";
  else if( num == 3006) 
    return "GLE";
  else if( num == 3007) 
    return "GLG";
  else if( num == 3008) 
    return "GLH";
  else if( num == 3009) 
    return "GLI";
  else if( num == 3010) 
    return "GLL";
  else if( num == 3011) 
    return "GLK";
  else if( num == 3012) 
    return "GLM";
  else if( num == 3013) 
    return "GLF";
  else if( num == 3014) 
    return "GLP";
  else if( num == 3015) 
    return "GLS";
  else if( num == 3016) 
    return "GLT";
  else if( num == 3017) 
    return "GLW";
  else if( num == 3018) 
    return "GLY";
  else if( num == 3019) 
    return "GLV";
  else if( num == 3020) 
    return "GKA";
  else if( num == 3021) 
    return "GKR";
  else if( num == 3022) 
    return "GKN";
  else if( num == 3023) 
    return "GKD";
  else if( num == 3024) 
    return "GKC";
  else if( num == 3025) 
    return "GKQ";
  else if( num == 3026) 
    return "GKE";
  else if( num == 3027) 
    return "GKG";
  else if( num == 3028) 
    return "GKH";
  else if( num == 3029) 
    return "GKI";
  else if( num == 3030) 
    return "GKL";
  else if( num == 3031) 
    return "GKK";
  else if( num == 3032) 
    return "GKM";
  else if( num == 3033) 
    return "GKF";
  else if( num == 3034) 
    return "GKP";
  else if( num == 3035) 
    return "GKS";
  else if( num == 3036) 
    return "GKT";
  else if( num == 3037) 
    return "GKW";
  else if( num == 3038) 
    return "GKY";
  else if( num == 3039) 
    return "GKV";
  else if( num == 3040) 
    return "GMA";
  else if( num == 3041) 
    return "GMR";
  else if( num == 3042) 
    return "GMN";
  else if( num == 3043) 
    return "GMD";
  else if( num == 3044) 
    return "GMC";
  else if( num == 3045) 
    return "GMQ";
  else if( num == 3046) 
    return "GME";
  else if( num == 3047) 
    return "GMG";
  else if( num == 3048) 
    return "GMH";
  else if( num == 3049) 
    return "GMI";
  else if( num == 3050) 
    return "GML";
  else if( num == 3051) 
    return "GMK";
  else if( num == 3052) 
    return "GMM";
  else if( num == 3053) 
    return "GMF";
  else if( num == 3054) 
    return "GMP";
  else if( num == 3055) 
    return "GMS";
  else if( num == 3056) 
    return "GMT";
  else if( num == 3057) 
    return "GMW";
  else if( num == 3058) 
    return "GMY";
  else if( num == 3059) 
    return "GMV";
  else if( num == 3060) 
    return "GFA";
  else if( num == 3061) 
    return "GFR";
  else if( num == 3062) 
    return "GFN";
  else if( num == 3063) 
    return "GFD";
  else if( num == 3064) 
    return "GFC";
  else if( num == 3065) 
    return "GFQ";
  else if( num == 3066) 
    return "GFE";
  else if( num == 3067) 
    return "GFG";
  else if( num == 3068) 
    return "GFH";
  else if( num == 3069) 
    return "GFI";
  else if( num == 3070) 
    return "GFL";
  else if( num == 3071) 
    return "GFK";
  else if( num == 3072) 
    return "GFM";
  else if( num == 3073) 
    return "GFF";
  else if( num == 3074) 
    return "GFP";
  else if( num == 3075) 
    return "GFS";
  else if( num == 3076) 
    return "GFT";
  else if( num == 3077) 
    return "GFW";
  else if( num == 3078) 
    return "GFY";
  else if( num == 3079) 
    return "GFV";
  else if( num == 3080) 
    return "GPA";
  else if( num == 3081) 
    return "GPR";
  else if( num == 3082) 
    return "GPN";
  else if( num == 3083) 
    return "GPD";
  else if( num == 3084) 
    return "GPC";
  else if( num == 3085) 
    return "GPQ";
  else if( num == 3086) 
    return "GPE";
  else if( num == 3087) 
    return "GPG";
  else if( num == 3088) 
    return "GPH";
  else if( num == 3089) 
    return "GPI";
  else if( num == 3090) 
    return "GPL";
  else if( num == 3091) 
    return "GPK";
  else if( num == 3092) 
    return "GPM";
  else if( num == 3093) 
    return "GPF";
  else if( num == 3094) 
    return "GPP";
  else if( num == 3095) 
    return "GPS";
  else if( num == 3096) 
    return "GPT";
  else if( num == 3097) 
    return "GPW";
  else if( num == 3098) 
    return "GPY";
  else if( num == 3099) 
    return "GPV";
  else if( num == 3100) 
    return "GSA";
  else if( num == 3101) 
    return "GSR";
  else if( num == 3102) 
    return "GSN";
  else if( num == 3103) 
    return "GSD";
  else if( num == 3104) 
    return "GSC";
  else if( num == 3105) 
    return "GSQ";
  else if( num == 3106) 
    return "GSE";
  else if( num == 3107) 
    return "GSG";
  else if( num == 3108) 
    return "GSH";
  else if( num == 3109) 
    return "GSI";
  else if( num == 3110) 
    return "GSL";
  else if( num == 3111) 
    return "GSK";
  else if( num == 3112) 
    return "GSM";
  else if( num == 3113) 
    return "GSF";
  else if( num == 3114) 
    return "GSP";
  else if( num == 3115) 
    return "GSS";
  else if( num == 3116) 
    return "GST";
  else if( num == 3117) 
    return "GSW";
  else if( num == 3118) 
    return "GSY";
  else if( num == 3119) 
    return "GSV";
  else if( num == 3120) 
    return "GTA";
  else if( num == 3121) 
    return "GTR";
  else if( num == 3122) 
    return "GTN";
  else if( num == 3123) 
    return "GTD";
  else if( num == 3124) 
    return "GTC";
  else if( num == 3125) 
    return "GTQ";
  else if( num == 3126) 
    return "GTE";
  else if( num == 3127) 
    return "GTG";
  else if( num == 3128) 
    return "GTH";
  else if( num == 3129) 
    return "GTI";
  else if( num == 3130) 
    return "GTL";
  else if( num == 3131) 
    return "GTK";
  else if( num == 3132) 
    return "GTM";
  else if( num == 3133) 
    return "GTF";
  else if( num == 3134) 
    return "GTP";
  else if( num == 3135) 
    return "GTS";
  else if( num == 3136) 
    return "GTT";
  else if( num == 3137) 
    return "GTW";
  else if( num == 3138) 
    return "GTY";
  else if( num == 3139) 
    return "GTV";
  else if( num == 3140) 
    return "GWA";
  else if( num == 3141) 
    return "GWR";
  else if( num == 3142) 
    return "GWN";
  else if( num == 3143) 
    return "GWD";
  else if( num == 3144) 
    return "GWC";
  else if( num == 3145) 
    return "GWQ";
  else if( num == 3146) 
    return "GWE";
  else if( num == 3147) 
    return "GWG";
  else if( num == 3148) 
    return "GWH";
  else if( num == 3149) 
    return "GWI";
  else if( num == 3150) 
    return "GWL";
  else if( num == 3151) 
    return "GWK";
  else if( num == 3152) 
    return "GWM";
  else if( num == 3153) 
    return "GWF";
  else if( num == 3154) 
    return "GWP";
  else if( num == 3155) 
    return "GWS";
  else if( num == 3156) 
    return "GWT";
  else if( num == 3157) 
    return "GWW";
  else if( num == 3158) 
    return "GWY";
  else if( num == 3159) 
    return "GWV";
  else if( num == 3160) 
    return "GYA";
  else if( num == 3161) 
    return "GYR";
  else if( num == 3162) 
    return "GYN";
  else if( num == 3163) 
    return "GYD";
  else if( num == 3164) 
    return "GYC";
  else if( num == 3165) 
    return "GYQ";
  else if( num == 3166) 
    return "GYE";
  else if( num == 3167) 
    return "GYG";
  else if( num == 3168) 
    return "GYH";
  else if( num == 3169) 
    return "GYI";
  else if( num == 3170) 
    return "GYL";
  else if( num == 3171) 
    return "GYK";
  else if( num == 3172) 
    return "GYM";
  else if( num == 3173) 
    return "GYF";
  else if( num == 3174) 
    return "GYP";
  else if( num == 3175) 
    return "GYS";
  else if( num == 3176) 
    return "GYT";
  else if( num == 3177) 
    return "GYW";
  else if( num == 3178) 
    return "GYY";
  else if( num == 3179) 
    return "GYV";
  else if( num == 3180) 
    return "GVA";
  else if( num == 3181) 
    return "GVR";
  else if( num == 3182) 
    return "GVN";
  else if( num == 3183) 
    return "GVD";
  else if( num == 3184) 
    return "GVC";
  else if( num == 3185) 
    return "GVQ";
  else if( num == 3186) 
    return "GVE";
  else if( num == 3187) 
    return "GVG";
  else if( num == 3188) 
    return "GVH";
  else if( num == 3189) 
    return "GVI";
  else if( num == 3190) 
    return "GVL";
  else if( num == 3191) 
    return "GVK";
  else if( num == 3192) 
    return "GVM";
  else if( num == 3193) 
    return "GVF";
  else if( num == 3194) 
    return "GVP";
  else if( num == 3195) 
    return "GVS";
  else if( num == 3196) 
    return "GVT";
  else if( num == 3197) 
    return "GVW";
  else if( num == 3198) 
    return "GVY";
  else if( num == 3199) 
    return "GVV";
  else if( num == 3200) 
    return "HAA";
  else if( num == 3201) 
    return "HAR";
  else if( num == 3202) 
    return "HAN";
  else if( num == 3203) 
    return "HAD";
  else if( num == 3204) 
    return "HAC";
  else if( num == 3205) 
    return "HAQ";
  else if( num == 3206) 
    return "HAE";
  else if( num == 3207) 
    return "HAG";
  else if( num == 3208) 
    return "HAH";
  else if( num == 3209) 
    return "HAI";
  else if( num == 3210) 
    return "HAL";
  else if( num == 3211) 
    return "HAK";
  else if( num == 3212) 
    return "HAM";
  else if( num == 3213) 
    return "HAF";
  else if( num == 3214) 
    return "HAP";
  else if( num == 3215) 
    return "HAS";
  else if( num == 3216) 
    return "HAT";
  else if( num == 3217) 
    return "HAW";
  else if( num == 3218) 
    return "HAY";
  else if( num == 3219) 
    return "HAV";
  else if( num == 3220) 
    return "HRA";
  else if( num == 3221) 
    return "HRR";
  else if( num == 3222) 
    return "HRN";
  else if( num == 3223) 
    return "HRD";
  else if( num == 3224) 
    return "HRC";
  else if( num == 3225) 
    return "HRQ";
  else if( num == 3226) 
    return "HRE";
  else if( num == 3227) 
    return "HRG";
  else if( num == 3228) 
    return "HRH";
  else if( num == 3229) 
    return "HRI";
  else if( num == 3230) 
    return "HRL";
  else if( num == 3231) 
    return "HRK";
  else if( num == 3232) 
    return "HRM";
  else if( num == 3233) 
    return "HRF";
  else if( num == 3234) 
    return "HRP";
  else if( num == 3235) 
    return "HRS";
  else if( num == 3236) 
    return "HRT";
  else if( num == 3237) 
    return "HRW";
  else if( num == 3238) 
    return "HRY";
  else if( num == 3239) 
    return "HRV";
  else if( num == 3240) 
    return "HNA";
  else if( num == 3241) 
    return "HNR";
  else if( num == 3242) 
    return "HNN";
  else if( num == 3243) 
    return "HND";
  else if( num == 3244) 
    return "HNC";
  else if( num == 3245) 
    return "HNQ";
  else if( num == 3246) 
    return "HNE";
  else if( num == 3247) 
    return "HNG";
  else if( num == 3248) 
    return "HNH";
  else if( num == 3249) 
    return "HNI";
  else if( num == 3250) 
    return "HNL";
  else if( num == 3251) 
    return "HNK";
  else if( num == 3252) 
    return "HNM";
  else if( num == 3253) 
    return "HNF";
  else if( num == 3254) 
    return "HNP";
  else if( num == 3255) 
    return "HNS";
  else if( num == 3256) 
    return "HNT";
  else if( num == 3257) 
    return "HNW";
  else if( num == 3258) 
    return "HNY";
  else if( num == 3259) 
    return "HNV";
  else if( num == 3260) 
    return "HDA";
  else if( num == 3261) 
    return "HDR";
  else if( num == 3262) 
    return "HDN";
  else if( num == 3263) 
    return "HDD";
  else if( num == 3264) 
    return "HDC";
  else if( num == 3265) 
    return "HDQ";
  else if( num == 3266) 
    return "HDE";
  else if( num == 3267) 
    return "HDG";
  else if( num == 3268) 
    return "HDH";
  else if( num == 3269) 
    return "HDI";
  else if( num == 3270) 
    return "HDL";
  else if( num == 3271) 
    return "HDK";
  else if( num == 3272) 
    return "HDM";
  else if( num == 3273) 
    return "HDF";
  else if( num == 3274) 
    return "HDP";
  else if( num == 3275) 
    return "HDS";
  else if( num == 3276) 
    return "HDT";
  else if( num == 3277) 
    return "HDW";
  else if( num == 3278) 
    return "HDY";
  else if( num == 3279) 
    return "HDV";
  else if( num == 3280) 
    return "HCA";
  else if( num == 3281) 
    return "HCR";
  else if( num == 3282) 
    return "HCN";
  else if( num == 3283) 
    return "HCD";
  else if( num == 3284) 
    return "HCC";
  else if( num == 3285) 
    return "HCQ";
  else if( num == 3286) 
    return "HCE";
  else if( num == 3287) 
    return "HCG";
  else if( num == 3288) 
    return "HCH";
  else if( num == 3289) 
    return "HCI";
  else if( num == 3290) 
    return "HCL";
  else if( num == 3291) 
    return "HCK";
  else if( num == 3292) 
    return "HCM";
  else if( num == 3293) 
    return "HCF";
  else if( num == 3294) 
    return "HCP";
  else if( num == 3295) 
    return "HCS";
  else if( num == 3296) 
    return "HCT";
  else if( num == 3297) 
    return "HCW";
  else if( num == 3298) 
    return "HCY";
  else if( num == 3299) 
    return "HCV";
  else if( num == 3300) 
    return "HQA";
  else if( num == 3301) 
    return "HQR";
  else if( num == 3302) 
    return "HQN";
  else if( num == 3303) 
    return "HQD";
  else if( num == 3304) 
    return "HQC";
  else if( num == 3305) 
    return "HQQ";
  else if( num == 3306) 
    return "HQE";
  else if( num == 3307) 
    return "HQG";
  else if( num == 3308) 
    return "HQH";
  else if( num == 3309) 
    return "HQI";
  else if( num == 3310) 
    return "HQL";
  else if( num == 3311) 
    return "HQK";
  else if( num == 3312) 
    return "HQM";
  else if( num == 3313) 
    return "HQF";
  else if( num == 3314) 
    return "HQP";
  else if( num == 3315) 
    return "HQS";
  else if( num == 3316) 
    return "HQT";
  else if( num == 3317) 
    return "HQW";
  else if( num == 3318) 
    return "HQY";
  else if( num == 3319) 
    return "HQV";
  else if( num == 3320) 
    return "HEA";
  else if( num == 3321) 
    return "HER";
  else if( num == 3322) 
    return "HEN";
  else if( num == 3323) 
    return "HED";
  else if( num == 3324) 
    return "HEC";
  else if( num == 3325) 
    return "HEQ";
  else if( num == 3326) 
    return "HEE";
  else if( num == 3327) 
    return "HEG";
  else if( num == 3328) 
    return "HEH";
  else if( num == 3329) 
    return "HEI";
  else if( num == 3330) 
    return "HEL";
  else if( num == 3331) 
    return "HEK";
  else if( num == 3332) 
    return "HEM";
  else if( num == 3333) 
    return "HEF";
  else if( num == 3334) 
    return "HEP";
  else if( num == 3335) 
    return "HES";
  else if( num == 3336) 
    return "HET";
  else if( num == 3337) 
    return "HEW";
  else if( num == 3338) 
    return "HEY";
  else if( num == 3339) 
    return "HEV";
  else if( num == 3340) 
    return "HGA";
  else if( num == 3341) 
    return "HGR";
  else if( num == 3342) 
    return "HGN";
  else if( num == 3343) 
    return "HGD";
  else if( num == 3344) 
    return "HGC";
  else if( num == 3345) 
    return "HGQ";
  else if( num == 3346) 
    return "HGE";
  else if( num == 3347) 
    return "HGG";
  else if( num == 3348) 
    return "HGH";
  else if( num == 3349) 
    return "HGI";
  else if( num == 3350) 
    return "HGL";
  else if( num == 3351) 
    return "HGK";
  else if( num == 3352) 
    return "HGM";
  else if( num == 3353) 
    return "HGF";
  else if( num == 3354) 
    return "HGP";
  else if( num == 3355) 
    return "HGS";
  else if( num == 3356) 
    return "HGT";
  else if( num == 3357) 
    return "HGW";
  else if( num == 3358) 
    return "HGY";
  else if( num == 3359) 
    return "HGV";
  else if( num == 3360) 
    return "HHA";
  else if( num == 3361) 
    return "HHR";
  else if( num == 3362) 
    return "HHN";
  else if( num == 3363) 
    return "HHD";
  else if( num == 3364) 
    return "HHC";
  else if( num == 3365) 
    return "HHQ";
  else if( num == 3366) 
    return "HHE";
  else if( num == 3367) 
    return "HHG";
  else if( num == 3368) 
    return "HHH";
  else if( num == 3369) 
    return "HHI";
  else if( num == 3370) 
    return "HHL";
  else if( num == 3371) 
    return "HHK";
  else if( num == 3372) 
    return "HHM";
  else if( num == 3373) 
    return "HHF";
  else if( num == 3374) 
    return "HHP";
  else if( num == 3375) 
    return "HHS";
  else if( num == 3376) 
    return "HHT";
  else if( num == 3377) 
    return "HHW";
  else if( num == 3378) 
    return "HHY";
  else if( num == 3379) 
    return "HHV";
  else if( num == 3380) 
    return "HIA";
  else if( num == 3381) 
    return "HIR";
  else if( num == 3382) 
    return "HIN";
  else if( num == 3383) 
    return "HID";
  else if( num == 3384) 
    return "HIC";
  else if( num == 3385) 
    return "HIQ";
  else if( num == 3386) 
    return "HIE";
  else if( num == 3387) 
    return "HIG";
  else if( num == 3388) 
    return "HIH";
  else if( num == 3389) 
    return "HII";
  else if( num == 3390) 
    return "HIL";
  else if( num == 3391) 
    return "HIK";
  else if( num == 3392) 
    return "HIM";
  else if( num == 3393) 
    return "HIF";
  else if( num == 3394) 
    return "HIP";
  else if( num == 3395) 
    return "HIS";
  else if( num == 3396) 
    return "HIT";
  else if( num == 3397) 
    return "HIW";
  else if( num == 3398) 
    return "HIY";
  else if( num == 3399) 
    return "HIV";
  else if( num == 3400) 
    return "HLA";
  else if( num == 3401) 
    return "HLR";
  else if( num == 3402) 
    return "HLN";
  else if( num == 3403) 
    return "HLD";
  else if( num == 3404) 
    return "HLC";
  else if( num == 3405) 
    return "HLQ";
  else if( num == 3406) 
    return "HLE";
  else if( num == 3407) 
    return "HLG";
  else if( num == 3408) 
    return "HLH";
  else if( num == 3409) 
    return "HLI";
  else if( num == 3410) 
    return "HLL";
  else if( num == 3411) 
    return "HLK";
  else if( num == 3412) 
    return "HLM";
  else if( num == 3413) 
    return "HLF";
  else if( num == 3414) 
    return "HLP";
  else if( num == 3415) 
    return "HLS";
  else if( num == 3416) 
    return "HLT";
  else if( num == 3417) 
    return "HLW";
  else if( num == 3418) 
    return "HLY";
  else if( num == 3419) 
    return "HLV";
  else if( num == 3420) 
    return "HKA";
  else if( num == 3421) 
    return "HKR";
  else if( num == 3422) 
    return "HKN";
  else if( num == 3423) 
    return "HKD";
  else if( num == 3424) 
    return "HKC";
  else if( num == 3425) 
    return "HKQ";
  else if( num == 3426) 
    return "HKE";
  else if( num == 3427) 
    return "HKG";
  else if( num == 3428) 
    return "HKH";
  else if( num == 3429) 
    return "HKI";
  else if( num == 3430) 
    return "HKL";
  else if( num == 3431) 
    return "HKK";
  else if( num == 3432) 
    return "HKM";
  else if( num == 3433) 
    return "HKF";
  else if( num == 3434) 
    return "HKP";
  else if( num == 3435) 
    return "HKS";
  else if( num == 3436) 
    return "HKT";
  else if( num == 3437) 
    return "HKW";
  else if( num == 3438) 
    return "HKY";
  else if( num == 3439) 
    return "HKV";
  else if( num == 3440) 
    return "HMA";
  else if( num == 3441) 
    return "HMR";
  else if( num == 3442) 
    return "HMN";
  else if( num == 3443) 
    return "HMD";
  else if( num == 3444) 
    return "HMC";
  else if( num == 3445) 
    return "HMQ";
  else if( num == 3446) 
    return "HME";
  else if( num == 3447) 
    return "HMG";
  else if( num == 3448) 
    return "HMH";
  else if( num == 3449) 
    return "HMI";
  else if( num == 3450) 
    return "HML";
  else if( num == 3451) 
    return "HMK";
  else if( num == 3452) 
    return "HMM";
  else if( num == 3453) 
    return "HMF";
  else if( num == 3454) 
    return "HMP";
  else if( num == 3455) 
    return "HMS";
  else if( num == 3456) 
    return "HMT";
  else if( num == 3457) 
    return "HMW";
  else if( num == 3458) 
    return "HMY";
  else if( num == 3459) 
    return "HMV";
  else if( num == 3460) 
    return "HFA";
  else if( num == 3461) 
    return "HFR";
  else if( num == 3462) 
    return "HFN";
  else if( num == 3463) 
    return "HFD";
  else if( num == 3464) 
    return "HFC";
  else if( num == 3465) 
    return "HFQ";
  else if( num == 3466) 
    return "HFE";
  else if( num == 3467) 
    return "HFG";
  else if( num == 3468) 
    return "HFH";
  else if( num == 3469) 
    return "HFI";
  else if( num == 3470) 
    return "HFL";
  else if( num == 3471) 
    return "HFK";
  else if( num == 3472) 
    return "HFM";
  else if( num == 3473) 
    return "HFF";
  else if( num == 3474) 
    return "HFP";
  else if( num == 3475) 
    return "HFS";
  else if( num == 3476) 
    return "HFT";
  else if( num == 3477) 
    return "HFW";
  else if( num == 3478) 
    return "HFY";
  else if( num == 3479) 
    return "HFV";
  else if( num == 3480) 
    return "HPA";
  else if( num == 3481) 
    return "HPR";
  else if( num == 3482) 
    return "HPN";
  else if( num == 3483) 
    return "HPD";
  else if( num == 3484) 
    return "HPC";
  else if( num == 3485) 
    return "HPQ";
  else if( num == 3486) 
    return "HPE";
  else if( num == 3487) 
    return "HPG";
  else if( num == 3488) 
    return "HPH";
  else if( num == 3489) 
    return "HPI";
  else if( num == 3490) 
    return "HPL";
  else if( num == 3491) 
    return "HPK";
  else if( num == 3492) 
    return "HPM";
  else if( num == 3493) 
    return "HPF";
  else if( num == 3494) 
    return "HPP";
  else if( num == 3495) 
    return "HPS";
  else if( num == 3496) 
    return "HPT";
  else if( num == 3497) 
    return "HPW";
  else if( num == 3498) 
    return "HPY";
  else if( num == 3499) 
    return "HPV";
  else if( num == 3500) 
    return "HSA";
  else if( num == 3501) 
    return "HSR";
  else if( num == 3502) 
    return "HSN";
  else if( num == 3503) 
    return "HSD";
  else if( num == 3504) 
    return "HSC";
  else if( num == 3505) 
    return "HSQ";
  else if( num == 3506) 
    return "HSE";
  else if( num == 3507) 
    return "HSG";
  else if( num == 3508) 
    return "HSH";
  else if( num == 3509) 
    return "HSI";
  else if( num == 3510) 
    return "HSL";
  else if( num == 3511) 
    return "HSK";
  else if( num == 3512) 
    return "HSM";
  else if( num == 3513) 
    return "HSF";
  else if( num == 3514) 
    return "HSP";
  else if( num == 3515) 
    return "HSS";
  else if( num == 3516) 
    return "HST";
  else if( num == 3517) 
    return "HSW";
  else if( num == 3518) 
    return "HSY";
  else if( num == 3519) 
    return "HSV";
  else if( num == 3520) 
    return "HTA";
  else if( num == 3521) 
    return "HTR";
  else if( num == 3522) 
    return "HTN";
  else if( num == 3523) 
    return "HTD";
  else if( num == 3524) 
    return "HTC";
  else if( num == 3525) 
    return "HTQ";
  else if( num == 3526) 
    return "HTE";
  else if( num == 3527) 
    return "HTG";
  else if( num == 3528) 
    return "HTH";
  else if( num == 3529) 
    return "HTI";
  else if( num == 3530) 
    return "HTL";
  else if( num == 3531) 
    return "HTK";
  else if( num == 3532) 
    return "HTM";
  else if( num == 3533) 
    return "HTF";
  else if( num == 3534) 
    return "HTP";
  else if( num == 3535) 
    return "HTS";
  else if( num == 3536) 
    return "HTT";
  else if( num == 3537) 
    return "HTW";
  else if( num == 3538) 
    return "HTY";
  else if( num == 3539) 
    return "HTV";
  else if( num == 3540) 
    return "HWA";
  else if( num == 3541) 
    return "HWR";
  else if( num == 3542) 
    return "HWN";
  else if( num == 3543) 
    return "HWD";
  else if( num == 3544) 
    return "HWC";
  else if( num == 3545) 
    return "HWQ";
  else if( num == 3546) 
    return "HWE";
  else if( num == 3547) 
    return "HWG";
  else if( num == 3548) 
    return "HWH";
  else if( num == 3549) 
    return "HWI";
  else if( num == 3550) 
    return "HWL";
  else if( num == 3551) 
    return "HWK";
  else if( num == 3552) 
    return "HWM";
  else if( num == 3553) 
    return "HWF";
  else if( num == 3554) 
    return "HWP";
  else if( num == 3555) 
    return "HWS";
  else if( num == 3556) 
    return "HWT";
  else if( num == 3557) 
    return "HWW";
  else if( num == 3558) 
    return "HWY";
  else if( num == 3559) 
    return "HWV";
  else if( num == 3560) 
    return "HYA";
  else if( num == 3561) 
    return "HYR";
  else if( num == 3562) 
    return "HYN";
  else if( num == 3563) 
    return "HYD";
  else if( num == 3564) 
    return "HYC";
  else if( num == 3565) 
    return "HYQ";
  else if( num == 3566) 
    return "HYE";
  else if( num == 3567) 
    return "HYG";
  else if( num == 3568) 
    return "HYH";
  else if( num == 3569) 
    return "HYI";
  else if( num == 3570) 
    return "HYL";
  else if( num == 3571) 
    return "HYK";
  else if( num == 3572) 
    return "HYM";
  else if( num == 3573) 
    return "HYF";
  else if( num == 3574) 
    return "HYP";
  else if( num == 3575) 
    return "HYS";
  else if( num == 3576) 
    return "HYT";
  else if( num == 3577) 
    return "HYW";
  else if( num == 3578) 
    return "HYY";
  else if( num == 3579) 
    return "HYV";
  else if( num == 3580) 
    return "HVA";
  else if( num == 3581) 
    return "HVR";
  else if( num == 3582) 
    return "HVN";
  else if( num == 3583) 
    return "HVD";
  else if( num == 3584) 
    return "HVC";
  else if( num == 3585) 
    return "HVQ";
  else if( num == 3586) 
    return "HVE";
  else if( num == 3587) 
    return "HVG";
  else if( num == 3588) 
    return "HVH";
  else if( num == 3589) 
    return "HVI";
  else if( num == 3590) 
    return "HVL";
  else if( num == 3591) 
    return "HVK";
  else if( num == 3592) 
    return "HVM";
  else if( num == 3593) 
    return "HVF";
  else if( num == 3594) 
    return "HVP";
  else if( num == 3595) 
    return "HVS";
  else if( num == 3596) 
    return "HVT";
  else if( num == 3597) 
    return "HVW";
  else if( num == 3598) 
    return "HVY";
  else if( num == 3599) 
    return "HVV";
  else if( num == 3600) 
    return "IAA";
  else if( num == 3601) 
    return "IAR";
  else if( num == 3602) 
    return "IAN";
  else if( num == 3603) 
    return "IAD";
  else if( num == 3604) 
    return "IAC";
  else if( num == 3605) 
    return "IAQ";
  else if( num == 3606) 
    return "IAE";
  else if( num == 3607) 
    return "IAG";
  else if( num == 3608) 
    return "IAH";
  else if( num == 3609) 
    return "IAI";
  else if( num == 3610) 
    return "IAL";
  else if( num == 3611) 
    return "IAK";
  else if( num == 3612) 
    return "IAM";
  else if( num == 3613) 
    return "IAF";
  else if( num == 3614) 
    return "IAP";
  else if( num == 3615) 
    return "IAS";
  else if( num == 3616) 
    return "IAT";
  else if( num == 3617) 
    return "IAW";
  else if( num == 3618) 
    return "IAY";
  else if( num == 3619) 
    return "IAV";
  else if( num == 3620) 
    return "IRA";
  else if( num == 3621) 
    return "IRR";
  else if( num == 3622) 
    return "IRN";
  else if( num == 3623) 
    return "IRD";
  else if( num == 3624) 
    return "IRC";
  else if( num == 3625) 
    return "IRQ";
  else if( num == 3626) 
    return "IRE";
  else if( num == 3627) 
    return "IRG";
  else if( num == 3628) 
    return "IRH";
  else if( num == 3629) 
    return "IRI";
  else if( num == 3630) 
    return "IRL";
  else if( num == 3631) 
    return "IRK";
  else if( num == 3632) 
    return "IRM";
  else if( num == 3633) 
    return "IRF";
  else if( num == 3634) 
    return "IRP";
  else if( num == 3635) 
    return "IRS";
  else if( num == 3636) 
    return "IRT";
  else if( num == 3637) 
    return "IRW";
  else if( num == 3638) 
    return "IRY";
  else if( num == 3639) 
    return "IRV";
  else if( num == 3640) 
    return "INA";
  else if( num == 3641) 
    return "INR";
  else if( num == 3642) 
    return "INN";
  else if( num == 3643) 
    return "IND";
  else if( num == 3644) 
    return "INC";
  else if( num == 3645) 
    return "INQ";
  else if( num == 3646) 
    return "INE";
  else if( num == 3647) 
    return "ING";
  else if( num == 3648) 
    return "INH";
  else if( num == 3649) 
    return "INI";
  else if( num == 3650) 
    return "INL";
  else if( num == 3651) 
    return "INK";
  else if( num == 3652) 
    return "INM";
  else if( num == 3653) 
    return "INF";
  else if( num == 3654) 
    return "INP";
  else if( num == 3655) 
    return "INS";
  else if( num == 3656) 
    return "INT";
  else if( num == 3657) 
    return "INW";
  else if( num == 3658) 
    return "INY";
  else if( num == 3659) 
    return "INV";
  else if( num == 3660) 
    return "IDA";
  else if( num == 3661) 
    return "IDR";
  else if( num == 3662) 
    return "IDN";
  else if( num == 3663) 
    return "IDD";
  else if( num == 3664) 
    return "IDC";
  else if( num == 3665) 
    return "IDQ";
  else if( num == 3666) 
    return "IDE";
  else if( num == 3667) 
    return "IDG";
  else if( num == 3668) 
    return "IDH";
  else if( num == 3669) 
    return "IDI";
  else if( num == 3670) 
    return "IDL";
  else if( num == 3671) 
    return "IDK";
  else if( num == 3672) 
    return "IDM";
  else if( num == 3673) 
    return "IDF";
  else if( num == 3674) 
    return "IDP";
  else if( num == 3675) 
    return "IDS";
  else if( num == 3676) 
    return "IDT";
  else if( num == 3677) 
    return "IDW";
  else if( num == 3678) 
    return "IDY";
  else if( num == 3679) 
    return "IDV";
  else if( num == 3680) 
    return "ICA";
  else if( num == 3681) 
    return "ICR";
  else if( num == 3682) 
    return "ICN";
  else if( num == 3683) 
    return "ICD";
  else if( num == 3684) 
    return "ICC";
  else if( num == 3685) 
    return "ICQ";
  else if( num == 3686) 
    return "ICE";
  else if( num == 3687) 
    return "ICG";
  else if( num == 3688) 
    return "ICH";
  else if( num == 3689) 
    return "ICI";
  else if( num == 3690) 
    return "ICL";
  else if( num == 3691) 
    return "ICK";
  else if( num == 3692) 
    return "ICM";
  else if( num == 3693) 
    return "ICF";
  else if( num == 3694) 
    return "ICP";
  else if( num == 3695) 
    return "ICS";
  else if( num == 3696) 
    return "ICT";
  else if( num == 3697) 
    return "ICW";
  else if( num == 3698) 
    return "ICY";
  else if( num == 3699) 
    return "ICV";
  else if( num == 3700) 
    return "IQA";
  else if( num == 3701) 
    return "IQR";
  else if( num == 3702) 
    return "IQN";
  else if( num == 3703) 
    return "IQD";
  else if( num == 3704) 
    return "IQC";
  else if( num == 3705) 
    return "IQQ";
  else if( num == 3706) 
    return "IQE";
  else if( num == 3707) 
    return "IQG";
  else if( num == 3708) 
    return "IQH";
  else if( num == 3709) 
    return "IQI";
  else if( num == 3710) 
    return "IQL";
  else if( num == 3711) 
    return "IQK";
  else if( num == 3712) 
    return "IQM";
  else if( num == 3713) 
    return "IQF";
  else if( num == 3714) 
    return "IQP";
  else if( num == 3715) 
    return "IQS";
  else if( num == 3716) 
    return "IQT";
  else if( num == 3717) 
    return "IQW";
  else if( num == 3718) 
    return "IQY";
  else if( num == 3719) 
    return "IQV";
  else if( num == 3720) 
    return "IEA";
  else if( num == 3721) 
    return "IER";
  else if( num == 3722) 
    return "IEN";
  else if( num == 3723) 
    return "IED";
  else if( num == 3724) 
    return "IEC";
  else if( num == 3725) 
    return "IEQ";
  else if( num == 3726) 
    return "IEE";
  else if( num == 3727) 
    return "IEG";
  else if( num == 3728) 
    return "IEH";
  else if( num == 3729) 
    return "IEI";
  else if( num == 3730) 
    return "IEL";
  else if( num == 3731) 
    return "IEK";
  else if( num == 3732) 
    return "IEM";
  else if( num == 3733) 
    return "IEF";
  else if( num == 3734) 
    return "IEP";
  else if( num == 3735) 
    return "IES";
  else if( num == 3736) 
    return "IET";
  else if( num == 3737) 
    return "IEW";
  else if( num == 3738) 
    return "IEY";
  else if( num == 3739) 
    return "IEV";
  else if( num == 3740) 
    return "IGA";
  else if( num == 3741) 
    return "IGR";
  else if( num == 3742) 
    return "IGN";
  else if( num == 3743) 
    return "IGD";
  else if( num == 3744) 
    return "IGC";
  else if( num == 3745) 
    return "IGQ";
  else if( num == 3746) 
    return "IGE";
  else if( num == 3747) 
    return "IGG";
  else if( num == 3748) 
    return "IGH";
  else if( num == 3749) 
    return "IGI";
  else if( num == 3750) 
    return "IGL";
  else if( num == 3751) 
    return "IGK";
  else if( num == 3752) 
    return "IGM";
  else if( num == 3753) 
    return "IGF";
  else if( num == 3754) 
    return "IGP";
  else if( num == 3755) 
    return "IGS";
  else if( num == 3756) 
    return "IGT";
  else if( num == 3757) 
    return "IGW";
  else if( num == 3758) 
    return "IGY";
  else if( num == 3759) 
    return "IGV";
  else if( num == 3760) 
    return "IHA";
  else if( num == 3761) 
    return "IHR";
  else if( num == 3762) 
    return "IHN";
  else if( num == 3763) 
    return "IHD";
  else if( num == 3764) 
    return "IHC";
  else if( num == 3765) 
    return "IHQ";
  else if( num == 3766) 
    return "IHE";
  else if( num == 3767) 
    return "IHG";
  else if( num == 3768) 
    return "IHH";
  else if( num == 3769) 
    return "IHI";
  else if( num == 3770) 
    return "IHL";
  else if( num == 3771) 
    return "IHK";
  else if( num == 3772) 
    return "IHM";
  else if( num == 3773) 
    return "IHF";
  else if( num == 3774) 
    return "IHP";
  else if( num == 3775) 
    return "IHS";
  else if( num == 3776) 
    return "IHT";
  else if( num == 3777) 
    return "IHW";
  else if( num == 3778) 
    return "IHY";
  else if( num == 3779) 
    return "IHV";
  else if( num == 3780) 
    return "IIA";
  else if( num == 3781) 
    return "IIR";
  else if( num == 3782) 
    return "IIN";
  else if( num == 3783) 
    return "IID";
  else if( num == 3784) 
    return "IIC";
  else if( num == 3785) 
    return "IIQ";
  else if( num == 3786) 
    return "IIE";
  else if( num == 3787) 
    return "IIG";
  else if( num == 3788) 
    return "IIH";
  else if( num == 3789) 
    return "III";
  else if( num == 3790) 
    return "IIL";
  else if( num == 3791) 
    return "IIK";
  else if( num == 3792) 
    return "IIM";
  else if( num == 3793) 
    return "IIF";
  else if( num == 3794) 
    return "IIP";
  else if( num == 3795) 
    return "IIS";
  else if( num == 3796) 
    return "IIT";
  else if( num == 3797) 
    return "IIW";
  else if( num == 3798) 
    return "IIY";
  else if( num == 3799) 
    return "IIV";
  else if( num == 3800) 
    return "ILA";
  else if( num == 3801) 
    return "ILR";
  else if( num == 3802) 
    return "ILN";
  else if( num == 3803) 
    return "ILD";
  else if( num == 3804) 
    return "ILC";
  else if( num == 3805) 
    return "ILQ";
  else if( num == 3806) 
    return "ILE";
  else if( num == 3807) 
    return "ILG";
  else if( num == 3808) 
    return "ILH";
  else if( num == 3809) 
    return "ILI";
  else if( num == 3810) 
    return "ILL";
  else if( num == 3811) 
    return "ILK";
  else if( num == 3812) 
    return "ILM";
  else if( num == 3813) 
    return "ILF";
  else if( num == 3814) 
    return "ILP";
  else if( num == 3815) 
    return "ILS";
  else if( num == 3816) 
    return "ILT";
  else if( num == 3817) 
    return "ILW";
  else if( num == 3818) 
    return "ILY";
  else if( num == 3819) 
    return "ILV";
  else if( num == 3820) 
    return "IKA";
  else if( num == 3821) 
    return "IKR";
  else if( num == 3822) 
    return "IKN";
  else if( num == 3823) 
    return "IKD";
  else if( num == 3824) 
    return "IKC";
  else if( num == 3825) 
    return "IKQ";
  else if( num == 3826) 
    return "IKE";
  else if( num == 3827) 
    return "IKG";
  else if( num == 3828) 
    return "IKH";
  else if( num == 3829) 
    return "IKI";
  else if( num == 3830) 
    return "IKL";
  else if( num == 3831) 
    return "IKK";
  else if( num == 3832) 
    return "IKM";
  else if( num == 3833) 
    return "IKF";
  else if( num == 3834) 
    return "IKP";
  else if( num == 3835) 
    return "IKS";
  else if( num == 3836) 
    return "IKT";
  else if( num == 3837) 
    return "IKW";
  else if( num == 3838) 
    return "IKY";
  else if( num == 3839) 
    return "IKV";
  else if( num == 3840) 
    return "IMA";
  else if( num == 3841) 
    return "IMR";
  else if( num == 3842) 
    return "IMN";
  else if( num == 3843) 
    return "IMD";
  else if( num == 3844) 
    return "IMC";
  else if( num == 3845) 
    return "IMQ";
  else if( num == 3846) 
    return "IME";
  else if( num == 3847) 
    return "IMG";
  else if( num == 3848) 
    return "IMH";
  else if( num == 3849) 
    return "IMI";
  else if( num == 3850) 
    return "IML";
  else if( num == 3851) 
    return "IMK";
  else if( num == 3852) 
    return "IMM";
  else if( num == 3853) 
    return "IMF";
  else if( num == 3854) 
    return "IMP";
  else if( num == 3855) 
    return "IMS";
  else if( num == 3856) 
    return "IMT";
  else if( num == 3857) 
    return "IMW";
  else if( num == 3858) 
    return "IMY";
  else if( num == 3859) 
    return "IMV";
  else if( num == 3860) 
    return "IFA";
  else if( num == 3861) 
    return "IFR";
  else if( num == 3862) 
    return "IFN";
  else if( num == 3863) 
    return "IFD";
  else if( num == 3864) 
    return "IFC";
  else if( num == 3865) 
    return "IFQ";
  else if( num == 3866) 
    return "IFE";
  else if( num == 3867) 
    return "IFG";
  else if( num == 3868) 
    return "IFH";
  else if( num == 3869) 
    return "IFI";
  else if( num == 3870) 
    return "IFL";
  else if( num == 3871) 
    return "IFK";
  else if( num == 3872) 
    return "IFM";
  else if( num == 3873) 
    return "IFF";
  else if( num == 3874) 
    return "IFP";
  else if( num == 3875) 
    return "IFS";
  else if( num == 3876) 
    return "IFT";
  else if( num == 3877) 
    return "IFW";
  else if( num == 3878) 
    return "IFY";
  else if( num == 3879) 
    return "IFV";
  else if( num == 3880) 
    return "IPA";
  else if( num == 3881) 
    return "IPR";
  else if( num == 3882) 
    return "IPN";
  else if( num == 3883) 
    return "IPD";
  else if( num == 3884) 
    return "IPC";
  else if( num == 3885) 
    return "IPQ";
  else if( num == 3886) 
    return "IPE";
  else if( num == 3887) 
    return "IPG";
  else if( num == 3888) 
    return "IPH";
  else if( num == 3889) 
    return "IPI";
  else if( num == 3890) 
    return "IPL";
  else if( num == 3891) 
    return "IPK";
  else if( num == 3892) 
    return "IPM";
  else if( num == 3893) 
    return "IPF";
  else if( num == 3894) 
    return "IPP";
  else if( num == 3895) 
    return "IPS";
  else if( num == 3896) 
    return "IPT";
  else if( num == 3897) 
    return "IPW";
  else if( num == 3898) 
    return "IPY";
  else if( num == 3899) 
    return "IPV";
  else if( num == 3900) 
    return "ISA";
  else if( num == 3901) 
    return "ISR";
  else if( num == 3902) 
    return "ISN";
  else if( num == 3903) 
    return "ISD";
  else if( num == 3904) 
    return "ISC";
  else if( num == 3905) 
    return "ISQ";
  else if( num == 3906) 
    return "ISE";
  else if( num == 3907) 
    return "ISG";
  else if( num == 3908) 
    return "ISH";
  else if( num == 3909) 
    return "ISI";
  else if( num == 3910) 
    return "ISL";
  else if( num == 3911) 
    return "ISK";
  else if( num == 3912) 
    return "ISM";
  else if( num == 3913) 
    return "ISF";
  else if( num == 3914) 
    return "ISP";
  else if( num == 3915) 
    return "ISS";
  else if( num == 3916) 
    return "IST";
  else if( num == 3917) 
    return "ISW";
  else if( num == 3918) 
    return "ISY";
  else if( num == 3919) 
    return "ISV";
  else if( num == 3920) 
    return "ITA";
  else if( num == 3921) 
    return "ITR";
  else if( num == 3922) 
    return "ITN";
  else if( num == 3923) 
    return "ITD";
  else if( num == 3924) 
    return "ITC";
  else if( num == 3925) 
    return "ITQ";
  else if( num == 3926) 
    return "ITE";
  else if( num == 3927) 
    return "ITG";
  else if( num == 3928) 
    return "ITH";
  else if( num == 3929) 
    return "ITI";
  else if( num == 3930) 
    return "ITL";
  else if( num == 3931) 
    return "ITK";
  else if( num == 3932) 
    return "ITM";
  else if( num == 3933) 
    return "ITF";
  else if( num == 3934) 
    return "ITP";
  else if( num == 3935) 
    return "ITS";
  else if( num == 3936) 
    return "ITT";
  else if( num == 3937) 
    return "ITW";
  else if( num == 3938) 
    return "ITY";
  else if( num == 3939) 
    return "ITV";
  else if( num == 3940) 
    return "IWA";
  else if( num == 3941) 
    return "IWR";
  else if( num == 3942) 
    return "IWN";
  else if( num == 3943) 
    return "IWD";
  else if( num == 3944) 
    return "IWC";
  else if( num == 3945) 
    return "IWQ";
  else if( num == 3946) 
    return "IWE";
  else if( num == 3947) 
    return "IWG";
  else if( num == 3948) 
    return "IWH";
  else if( num == 3949) 
    return "IWI";
  else if( num == 3950) 
    return "IWL";
  else if( num == 3951) 
    return "IWK";
  else if( num == 3952) 
    return "IWM";
  else if( num == 3953) 
    return "IWF";
  else if( num == 3954) 
    return "IWP";
  else if( num == 3955) 
    return "IWS";
  else if( num == 3956) 
    return "IWT";
  else if( num == 3957) 
    return "IWW";
  else if( num == 3958) 
    return "IWY";
  else if( num == 3959) 
    return "IWV";
  else if( num == 3960) 
    return "IYA";
  else if( num == 3961) 
    return "IYR";
  else if( num == 3962) 
    return "IYN";
  else if( num == 3963) 
    return "IYD";
  else if( num == 3964) 
    return "IYC";
  else if( num == 3965) 
    return "IYQ";
  else if( num == 3966) 
    return "IYE";
  else if( num == 3967) 
    return "IYG";
  else if( num == 3968) 
    return "IYH";
  else if( num == 3969) 
    return "IYI";
  else if( num == 3970) 
    return "IYL";
  else if( num == 3971) 
    return "IYK";
  else if( num == 3972) 
    return "IYM";
  else if( num == 3973) 
    return "IYF";
  else if( num == 3974) 
    return "IYP";
  else if( num == 3975) 
    return "IYS";
  else if( num == 3976) 
    return "IYT";
  else if( num == 3977) 
    return "IYW";
  else if( num == 3978) 
    return "IYY";
  else if( num == 3979) 
    return "IYV";
  else if( num == 3980) 
    return "IVA";
  else if( num == 3981) 
    return "IVR";
  else if( num == 3982) 
    return "IVN";
  else if( num == 3983) 
    return "IVD";
  else if( num == 3984) 
    return "IVC";
  else if( num == 3985) 
    return "IVQ";
  else if( num == 3986) 
    return "IVE";
  else if( num == 3987) 
    return "IVG";
  else if( num == 3988) 
    return "IVH";
  else if( num == 3989) 
    return "IVI";
  else if( num == 3990) 
    return "IVL";
  else if( num == 3991) 
    return "IVK";
  else if( num == 3992) 
    return "IVM";
  else if( num == 3993) 
    return "IVF";
  else if( num == 3994) 
    return "IVP";
  else if( num == 3995) 
    return "IVS";
  else if( num == 3996) 
    return "IVT";
  else if( num == 3997) 
    return "IVW";
  else if( num == 3998) 
    return "IVY";
  else if( num == 3999) 
    return "IVV";
  else if( num == 4000) 
    return "LAA";
  else if( num == 4001) 
    return "LAR";
  else if( num == 4002) 
    return "LAN";
  else if( num == 4003) 
    return "LAD";
  else if( num == 4004) 
    return "LAC";
  else if( num == 4005) 
    return "LAQ";
  else if( num == 4006) 
    return "LAE";
  else if( num == 4007) 
    return "LAG";
  else if( num == 4008) 
    return "LAH";
  else if( num == 4009) 
    return "LAI";
  else if( num == 4010) 
    return "LAL";
  else if( num == 4011) 
    return "LAK";
  else if( num == 4012) 
    return "LAM";
  else if( num == 4013) 
    return "LAF";
  else if( num == 4014) 
    return "LAP";
  else if( num == 4015) 
    return "LAS";
  else if( num == 4016) 
    return "LAT";
  else if( num == 4017) 
    return "LAW";
  else if( num == 4018) 
    return "LAY";
  else if( num == 4019) 
    return "LAV";
  else if( num == 4020) 
    return "LRA";
  else if( num == 4021) 
    return "LRR";
  else if( num == 4022) 
    return "LRN";
  else if( num == 4023) 
    return "LRD";
  else if( num == 4024) 
    return "LRC";
  else if( num == 4025) 
    return "LRQ";
  else if( num == 4026) 
    return "LRE";
  else if( num == 4027) 
    return "LRG";
  else if( num == 4028) 
    return "LRH";
  else if( num == 4029) 
    return "LRI";
  else if( num == 4030) 
    return "LRL";
  else if( num == 4031) 
    return "LRK";
  else if( num == 4032) 
    return "LRM";
  else if( num == 4033) 
    return "LRF";
  else if( num == 4034) 
    return "LRP";
  else if( num == 4035) 
    return "LRS";
  else if( num == 4036) 
    return "LRT";
  else if( num == 4037) 
    return "LRW";
  else if( num == 4038) 
    return "LRY";
  else if( num == 4039) 
    return "LRV";
  else if( num == 4040) 
    return "LNA";
  else if( num == 4041) 
    return "LNR";
  else if( num == 4042) 
    return "LNN";
  else if( num == 4043) 
    return "LND";
  else if( num == 4044) 
    return "LNC";
  else if( num == 4045) 
    return "LNQ";
  else if( num == 4046) 
    return "LNE";
  else if( num == 4047) 
    return "LNG";
  else if( num == 4048) 
    return "LNH";
  else if( num == 4049) 
    return "LNI";
  else if( num == 4050) 
    return "LNL";
  else if( num == 4051) 
    return "LNK";
  else if( num == 4052) 
    return "LNM";
  else if( num == 4053) 
    return "LNF";
  else if( num == 4054) 
    return "LNP";
  else if( num == 4055) 
    return "LNS";
  else if( num == 4056) 
    return "LNT";
  else if( num == 4057) 
    return "LNW";
  else if( num == 4058) 
    return "LNY";
  else if( num == 4059) 
    return "LNV";
  else if( num == 4060) 
    return "LDA";
  else if( num == 4061) 
    return "LDR";
  else if( num == 4062) 
    return "LDN";
  else if( num == 4063) 
    return "LDD";
  else if( num == 4064) 
    return "LDC";
  else if( num == 4065) 
    return "LDQ";
  else if( num == 4066) 
    return "LDE";
  else if( num == 4067) 
    return "LDG";
  else if( num == 4068) 
    return "LDH";
  else if( num == 4069) 
    return "LDI";
  else if( num == 4070) 
    return "LDL";
  else if( num == 4071) 
    return "LDK";
  else if( num == 4072) 
    return "LDM";
  else if( num == 4073) 
    return "LDF";
  else if( num == 4074) 
    return "LDP";
  else if( num == 4075) 
    return "LDS";
  else if( num == 4076) 
    return "LDT";
  else if( num == 4077) 
    return "LDW";
  else if( num == 4078) 
    return "LDY";
  else if( num == 4079) 
    return "LDV";
  else if( num == 4080) 
    return "LCA";
  else if( num == 4081) 
    return "LCR";
  else if( num == 4082) 
    return "LCN";
  else if( num == 4083) 
    return "LCD";
  else if( num == 4084) 
    return "LCC";
  else if( num == 4085) 
    return "LCQ";
  else if( num == 4086) 
    return "LCE";
  else if( num == 4087) 
    return "LCG";
  else if( num == 4088) 
    return "LCH";
  else if( num == 4089) 
    return "LCI";
  else if( num == 4090) 
    return "LCL";
  else if( num == 4091) 
    return "LCK";
  else if( num == 4092) 
    return "LCM";
  else if( num == 4093) 
    return "LCF";
  else if( num == 4094) 
    return "LCP";
  else if( num == 4095) 
    return "LCS";
  else if( num == 4096) 
    return "LCT";
  else if( num == 4097) 
    return "LCW";
  else if( num == 4098) 
    return "LCY";
  else if( num == 4099) 
    return "LCV";
  else if( num == 4100) 
    return "LQA";
  else if( num == 4101) 
    return "LQR";
  else if( num == 4102) 
    return "LQN";
  else if( num == 4103) 
    return "LQD";
  else if( num == 4104) 
    return "LQC";
  else if( num == 4105) 
    return "LQQ";
  else if( num == 4106) 
    return "LQE";
  else if( num == 4107) 
    return "LQG";
  else if( num == 4108) 
    return "LQH";
  else if( num == 4109) 
    return "LQI";
  else if( num == 4110) 
    return "LQL";
  else if( num == 4111) 
    return "LQK";
  else if( num == 4112) 
    return "LQM";
  else if( num == 4113) 
    return "LQF";
  else if( num == 4114) 
    return "LQP";
  else if( num == 4115) 
    return "LQS";
  else if( num == 4116) 
    return "LQT";
  else if( num == 4117) 
    return "LQW";
  else if( num == 4118) 
    return "LQY";
  else if( num == 4119) 
    return "LQV";
  else if( num == 4120) 
    return "LEA";
  else if( num == 4121) 
    return "LER";
  else if( num == 4122) 
    return "LEN";
  else if( num == 4123) 
    return "LED";
  else if( num == 4124) 
    return "LEC";
  else if( num == 4125) 
    return "LEQ";
  else if( num == 4126) 
    return "LEE";
  else if( num == 4127) 
    return "LEG";
  else if( num == 4128) 
    return "LEH";
  else if( num == 4129) 
    return "LEI";
  else if( num == 4130) 
    return "LEL";
  else if( num == 4131) 
    return "LEK";
  else if( num == 4132) 
    return "LEM";
  else if( num == 4133) 
    return "LEF";
  else if( num == 4134) 
    return "LEP";
  else if( num == 4135) 
    return "LES";
  else if( num == 4136) 
    return "LET";
  else if( num == 4137) 
    return "LEW";
  else if( num == 4138) 
    return "LEY";
  else if( num == 4139) 
    return "LEV";
  else if( num == 4140) 
    return "LGA";
  else if( num == 4141) 
    return "LGR";
  else if( num == 4142) 
    return "LGN";
  else if( num == 4143) 
    return "LGD";
  else if( num == 4144) 
    return "LGC";
  else if( num == 4145) 
    return "LGQ";
  else if( num == 4146) 
    return "LGE";
  else if( num == 4147) 
    return "LGG";
  else if( num == 4148) 
    return "LGH";
  else if( num == 4149) 
    return "LGI";
  else if( num == 4150) 
    return "LGL";
  else if( num == 4151) 
    return "LGK";
  else if( num == 4152) 
    return "LGM";
  else if( num == 4153) 
    return "LGF";
  else if( num == 4154) 
    return "LGP";
  else if( num == 4155) 
    return "LGS";
  else if( num == 4156) 
    return "LGT";
  else if( num == 4157) 
    return "LGW";
  else if( num == 4158) 
    return "LGY";
  else if( num == 4159) 
    return "LGV";
  else if( num == 4160) 
    return "LHA";
  else if( num == 4161) 
    return "LHR";
  else if( num == 4162) 
    return "LHN";
  else if( num == 4163) 
    return "LHD";
  else if( num == 4164) 
    return "LHC";
  else if( num == 4165) 
    return "LHQ";
  else if( num == 4166) 
    return "LHE";
  else if( num == 4167) 
    return "LHG";
  else if( num == 4168) 
    return "LHH";
  else if( num == 4169) 
    return "LHI";
  else if( num == 4170) 
    return "LHL";
  else if( num == 4171) 
    return "LHK";
  else if( num == 4172) 
    return "LHM";
  else if( num == 4173) 
    return "LHF";
  else if( num == 4174) 
    return "LHP";
  else if( num == 4175) 
    return "LHS";
  else if( num == 4176) 
    return "LHT";
  else if( num == 4177) 
    return "LHW";
  else if( num == 4178) 
    return "LHY";
  else if( num == 4179) 
    return "LHV";
  else if( num == 4180) 
    return "LIA";
  else if( num == 4181) 
    return "LIR";
  else if( num == 4182) 
    return "LIN";
  else if( num == 4183) 
    return "LID";
  else if( num == 4184) 
    return "LIC";
  else if( num == 4185) 
    return "LIQ";
  else if( num == 4186) 
    return "LIE";
  else if( num == 4187) 
    return "LIG";
  else if( num == 4188) 
    return "LIH";
  else if( num == 4189) 
    return "LII";
  else if( num == 4190) 
    return "LIL";
  else if( num == 4191) 
    return "LIK";
  else if( num == 4192) 
    return "LIM";
  else if( num == 4193) 
    return "LIF";
  else if( num == 4194) 
    return "LIP";
  else if( num == 4195) 
    return "LIS";
  else if( num == 4196) 
    return "LIT";
  else if( num == 4197) 
    return "LIW";
  else if( num == 4198) 
    return "LIY";
  else if( num == 4199) 
    return "LIV";
  else if( num == 4200) 
    return "LLA";
  else if( num == 4201) 
    return "LLR";
  else if( num == 4202) 
    return "LLN";
  else if( num == 4203) 
    return "LLD";
  else if( num == 4204) 
    return "LLC";
  else if( num == 4205) 
    return "LLQ";
  else if( num == 4206) 
    return "LLE";
  else if( num == 4207) 
    return "LLG";
  else if( num == 4208) 
    return "LLH";
  else if( num == 4209) 
    return "LLI";
  else if( num == 4210) 
    return "LLL";
  else if( num == 4211) 
    return "LLK";
  else if( num == 4212) 
    return "LLM";
  else if( num == 4213) 
    return "LLF";
  else if( num == 4214) 
    return "LLP";
  else if( num == 4215) 
    return "LLS";
  else if( num == 4216) 
    return "LLT";
  else if( num == 4217) 
    return "LLW";
  else if( num == 4218) 
    return "LLY";
  else if( num == 4219) 
    return "LLV";
  else if( num == 4220) 
    return "LKA";
  else if( num == 4221) 
    return "LKR";
  else if( num == 4222) 
    return "LKN";
  else if( num == 4223) 
    return "LKD";
  else if( num == 4224) 
    return "LKC";
  else if( num == 4225) 
    return "LKQ";
  else if( num == 4226) 
    return "LKE";
  else if( num == 4227) 
    return "LKG";
  else if( num == 4228) 
    return "LKH";
  else if( num == 4229) 
    return "LKI";
  else if( num == 4230) 
    return "LKL";
  else if( num == 4231) 
    return "LKK";
  else if( num == 4232) 
    return "LKM";
  else if( num == 4233) 
    return "LKF";
  else if( num == 4234) 
    return "LKP";
  else if( num == 4235) 
    return "LKS";
  else if( num == 4236) 
    return "LKT";
  else if( num == 4237) 
    return "LKW";
  else if( num == 4238) 
    return "LKY";
  else if( num == 4239) 
    return "LKV";
  else if( num == 4240) 
    return "LMA";
  else if( num == 4241) 
    return "LMR";
  else if( num == 4242) 
    return "LMN";
  else if( num == 4243) 
    return "LMD";
  else if( num == 4244) 
    return "LMC";
  else if( num == 4245) 
    return "LMQ";
  else if( num == 4246) 
    return "LME";
  else if( num == 4247) 
    return "LMG";
  else if( num == 4248) 
    return "LMH";
  else if( num == 4249) 
    return "LMI";
  else if( num == 4250) 
    return "LML";
  else if( num == 4251) 
    return "LMK";
  else if( num == 4252) 
    return "LMM";
  else if( num == 4253) 
    return "LMF";
  else if( num == 4254) 
    return "LMP";
  else if( num == 4255) 
    return "LMS";
  else if( num == 4256) 
    return "LMT";
  else if( num == 4257) 
    return "LMW";
  else if( num == 4258) 
    return "LMY";
  else if( num == 4259) 
    return "LMV";
  else if( num == 4260) 
    return "LFA";
  else if( num == 4261) 
    return "LFR";
  else if( num == 4262) 
    return "LFN";
  else if( num == 4263) 
    return "LFD";
  else if( num == 4264) 
    return "LFC";
  else if( num == 4265) 
    return "LFQ";
  else if( num == 4266) 
    return "LFE";
  else if( num == 4267) 
    return "LFG";
  else if( num == 4268) 
    return "LFH";
  else if( num == 4269) 
    return "LFI";
  else if( num == 4270) 
    return "LFL";
  else if( num == 4271) 
    return "LFK";
  else if( num == 4272) 
    return "LFM";
  else if( num == 4273) 
    return "LFF";
  else if( num == 4274) 
    return "LFP";
  else if( num == 4275) 
    return "LFS";
  else if( num == 4276) 
    return "LFT";
  else if( num == 4277) 
    return "LFW";
  else if( num == 4278) 
    return "LFY";
  else if( num == 4279) 
    return "LFV";
  else if( num == 4280) 
    return "LPA";
  else if( num == 4281) 
    return "LPR";
  else if( num == 4282) 
    return "LPN";
  else if( num == 4283) 
    return "LPD";
  else if( num == 4284) 
    return "LPC";
  else if( num == 4285) 
    return "LPQ";
  else if( num == 4286) 
    return "LPE";
  else if( num == 4287) 
    return "LPG";
  else if( num == 4288) 
    return "LPH";
  else if( num == 4289) 
    return "LPI";
  else if( num == 4290) 
    return "LPL";
  else if( num == 4291) 
    return "LPK";
  else if( num == 4292) 
    return "LPM";
  else if( num == 4293) 
    return "LPF";
  else if( num == 4294) 
    return "LPP";
  else if( num == 4295) 
    return "LPS";
  else if( num == 4296) 
    return "LPT";
  else if( num == 4297) 
    return "LPW";
  else if( num == 4298) 
    return "LPY";
  else if( num == 4299) 
    return "LPV";
  else if( num == 4300) 
    return "LSA";
  else if( num == 4301) 
    return "LSR";
  else if( num == 4302) 
    return "LSN";
  else if( num == 4303) 
    return "LSD";
  else if( num == 4304) 
    return "LSC";
  else if( num == 4305) 
    return "LSQ";
  else if( num == 4306) 
    return "LSE";
  else if( num == 4307) 
    return "LSG";
  else if( num == 4308) 
    return "LSH";
  else if( num == 4309) 
    return "LSI";
  else if( num == 4310) 
    return "LSL";
  else if( num == 4311) 
    return "LSK";
  else if( num == 4312) 
    return "LSM";
  else if( num == 4313) 
    return "LSF";
  else if( num == 4314) 
    return "LSP";
  else if( num == 4315) 
    return "LSS";
  else if( num == 4316) 
    return "LST";
  else if( num == 4317) 
    return "LSW";
  else if( num == 4318) 
    return "LSY";
  else if( num == 4319) 
    return "LSV";
  else if( num == 4320) 
    return "LTA";
  else if( num == 4321) 
    return "LTR";
  else if( num == 4322) 
    return "LTN";
  else if( num == 4323) 
    return "LTD";
  else if( num == 4324) 
    return "LTC";
  else if( num == 4325) 
    return "LTQ";
  else if( num == 4326) 
    return "LTE";
  else if( num == 4327) 
    return "LTG";
  else if( num == 4328) 
    return "LTH";
  else if( num == 4329) 
    return "LTI";
  else if( num == 4330) 
    return "LTL";
  else if( num == 4331) 
    return "LTK";
  else if( num == 4332) 
    return "LTM";
  else if( num == 4333) 
    return "LTF";
  else if( num == 4334) 
    return "LTP";
  else if( num == 4335) 
    return "LTS";
  else if( num == 4336) 
    return "LTT";
  else if( num == 4337) 
    return "LTW";
  else if( num == 4338) 
    return "LTY";
  else if( num == 4339) 
    return "LTV";
  else if( num == 4340) 
    return "LWA";
  else if( num == 4341) 
    return "LWR";
  else if( num == 4342) 
    return "LWN";
  else if( num == 4343) 
    return "LWD";
  else if( num == 4344) 
    return "LWC";
  else if( num == 4345) 
    return "LWQ";
  else if( num == 4346) 
    return "LWE";
  else if( num == 4347) 
    return "LWG";
  else if( num == 4348) 
    return "LWH";
  else if( num == 4349) 
    return "LWI";
  else if( num == 4350) 
    return "LWL";
  else if( num == 4351) 
    return "LWK";
  else if( num == 4352) 
    return "LWM";
  else if( num == 4353) 
    return "LWF";
  else if( num == 4354) 
    return "LWP";
  else if( num == 4355) 
    return "LWS";
  else if( num == 4356) 
    return "LWT";
  else if( num == 4357) 
    return "LWW";
  else if( num == 4358) 
    return "LWY";
  else if( num == 4359) 
    return "LWV";
  else if( num == 4360) 
    return "LYA";
  else if( num == 4361) 
    return "LYR";
  else if( num == 4362) 
    return "LYN";
  else if( num == 4363) 
    return "LYD";
  else if( num == 4364) 
    return "LYC";
  else if( num == 4365) 
    return "LYQ";
  else if( num == 4366) 
    return "LYE";
  else if( num == 4367) 
    return "LYG";
  else if( num == 4368) 
    return "LYH";
  else if( num == 4369) 
    return "LYI";
  else if( num == 4370) 
    return "LYL";
  else if( num == 4371) 
    return "LYK";
  else if( num == 4372) 
    return "LYM";
  else if( num == 4373) 
    return "LYF";
  else if( num == 4374) 
    return "LYP";
  else if( num == 4375) 
    return "LYS";
  else if( num == 4376) 
    return "LYT";
  else if( num == 4377) 
    return "LYW";
  else if( num == 4378) 
    return "LYY";
  else if( num == 4379) 
    return "LYV";
  else if( num == 4380) 
    return "LVA";
  else if( num == 4381) 
    return "LVR";
  else if( num == 4382) 
    return "LVN";
  else if( num == 4383) 
    return "LVD";
  else if( num == 4384) 
    return "LVC";
  else if( num == 4385) 
    return "LVQ";
  else if( num == 4386) 
    return "LVE";
  else if( num == 4387) 
    return "LVG";
  else if( num == 4388) 
    return "LVH";
  else if( num == 4389) 
    return "LVI";
  else if( num == 4390) 
    return "LVL";
  else if( num == 4391) 
    return "LVK";
  else if( num == 4392) 
    return "LVM";
  else if( num == 4393) 
    return "LVF";
  else if( num == 4394) 
    return "LVP";
  else if( num == 4395) 
    return "LVS";
  else if( num == 4396) 
    return "LVT";
  else if( num == 4397) 
    return "LVW";
  else if( num == 4398) 
    return "LVY";
  else if( num == 4399) 
    return "LVV";
  else if( num == 4400) 
    return "KAA";
  else if( num == 4401) 
    return "KAR";
  else if( num == 4402) 
    return "KAN";
  else if( num == 4403) 
    return "KAD";
  else if( num == 4404) 
    return "KAC";
  else if( num == 4405) 
    return "KAQ";
  else if( num == 4406) 
    return "KAE";
  else if( num == 4407) 
    return "KAG";
  else if( num == 4408) 
    return "KAH";
  else if( num == 4409) 
    return "KAI";
  else if( num == 4410) 
    return "KAL";
  else if( num == 4411) 
    return "KAK";
  else if( num == 4412) 
    return "KAM";
  else if( num == 4413) 
    return "KAF";
  else if( num == 4414) 
    return "KAP";
  else if( num == 4415) 
    return "KAS";
  else if( num == 4416) 
    return "KAT";
  else if( num == 4417) 
    return "KAW";
  else if( num == 4418) 
    return "KAY";
  else if( num == 4419) 
    return "KAV";
  else if( num == 4420) 
    return "KRA";
  else if( num == 4421) 
    return "KRR";
  else if( num == 4422) 
    return "KRN";
  else if( num == 4423) 
    return "KRD";
  else if( num == 4424) 
    return "KRC";
  else if( num == 4425) 
    return "KRQ";
  else if( num == 4426) 
    return "KRE";
  else if( num == 4427) 
    return "KRG";
  else if( num == 4428) 
    return "KRH";
  else if( num == 4429) 
    return "KRI";
  else if( num == 4430) 
    return "KRL";
  else if( num == 4431) 
    return "KRK";
  else if( num == 4432) 
    return "KRM";
  else if( num == 4433) 
    return "KRF";
  else if( num == 4434) 
    return "KRP";
  else if( num == 4435) 
    return "KRS";
  else if( num == 4436) 
    return "KRT";
  else if( num == 4437) 
    return "KRW";
  else if( num == 4438) 
    return "KRY";
  else if( num == 4439) 
    return "KRV";
  else if( num == 4440) 
    return "KNA";
  else if( num == 4441) 
    return "KNR";
  else if( num == 4442) 
    return "KNN";
  else if( num == 4443) 
    return "KND";
  else if( num == 4444) 
    return "KNC";
  else if( num == 4445) 
    return "KNQ";
  else if( num == 4446) 
    return "KNE";
  else if( num == 4447) 
    return "KNG";
  else if( num == 4448) 
    return "KNH";
  else if( num == 4449) 
    return "KNI";
  else if( num == 4450) 
    return "KNL";
  else if( num == 4451) 
    return "KNK";
  else if( num == 4452) 
    return "KNM";
  else if( num == 4453) 
    return "KNF";
  else if( num == 4454) 
    return "KNP";
  else if( num == 4455) 
    return "KNS";
  else if( num == 4456) 
    return "KNT";
  else if( num == 4457) 
    return "KNW";
  else if( num == 4458) 
    return "KNY";
  else if( num == 4459) 
    return "KNV";
  else if( num == 4460) 
    return "KDA";
  else if( num == 4461) 
    return "KDR";
  else if( num == 4462) 
    return "KDN";
  else if( num == 4463) 
    return "KDD";
  else if( num == 4464) 
    return "KDC";
  else if( num == 4465) 
    return "KDQ";
  else if( num == 4466) 
    return "KDE";
  else if( num == 4467) 
    return "KDG";
  else if( num == 4468) 
    return "KDH";
  else if( num == 4469) 
    return "KDI";
  else if( num == 4470) 
    return "KDL";
  else if( num == 4471) 
    return "KDK";
  else if( num == 4472) 
    return "KDM";
  else if( num == 4473) 
    return "KDF";
  else if( num == 4474) 
    return "KDP";
  else if( num == 4475) 
    return "KDS";
  else if( num == 4476) 
    return "KDT";
  else if( num == 4477) 
    return "KDW";
  else if( num == 4478) 
    return "KDY";
  else if( num == 4479) 
    return "KDV";
  else if( num == 4480) 
    return "KCA";
  else if( num == 4481) 
    return "KCR";
  else if( num == 4482) 
    return "KCN";
  else if( num == 4483) 
    return "KCD";
  else if( num == 4484) 
    return "KCC";
  else if( num == 4485) 
    return "KCQ";
  else if( num == 4486) 
    return "KCE";
  else if( num == 4487) 
    return "KCG";
  else if( num == 4488) 
    return "KCH";
  else if( num == 4489) 
    return "KCI";
  else if( num == 4490) 
    return "KCL";
  else if( num == 4491) 
    return "KCK";
  else if( num == 4492) 
    return "KCM";
  else if( num == 4493) 
    return "KCF";
  else if( num == 4494) 
    return "KCP";
  else if( num == 4495) 
    return "KCS";
  else if( num == 4496) 
    return "KCT";
  else if( num == 4497) 
    return "KCW";
  else if( num == 4498) 
    return "KCY";
  else if( num == 4499) 
    return "KCV";
  else if( num == 4500) 
    return "KQA";
  else if( num == 4501) 
    return "KQR";
  else if( num == 4502) 
    return "KQN";
  else if( num == 4503) 
    return "KQD";
  else if( num == 4504) 
    return "KQC";
  else if( num == 4505) 
    return "KQQ";
  else if( num == 4506) 
    return "KQE";
  else if( num == 4507) 
    return "KQG";
  else if( num == 4508) 
    return "KQH";
  else if( num == 4509) 
    return "KQI";
  else if( num == 4510) 
    return "KQL";
  else if( num == 4511) 
    return "KQK";
  else if( num == 4512) 
    return "KQM";
  else if( num == 4513) 
    return "KQF";
  else if( num == 4514) 
    return "KQP";
  else if( num == 4515) 
    return "KQS";
  else if( num == 4516) 
    return "KQT";
  else if( num == 4517) 
    return "KQW";
  else if( num == 4518) 
    return "KQY";
  else if( num == 4519) 
    return "KQV";
  else if( num == 4520) 
    return "KEA";
  else if( num == 4521) 
    return "KER";
  else if( num == 4522) 
    return "KEN";
  else if( num == 4523) 
    return "KED";
  else if( num == 4524) 
    return "KEC";
  else if( num == 4525) 
    return "KEQ";
  else if( num == 4526) 
    return "KEE";
  else if( num == 4527) 
    return "KEG";
  else if( num == 4528) 
    return "KEH";
  else if( num == 4529) 
    return "KEI";
  else if( num == 4530) 
    return "KEL";
  else if( num == 4531) 
    return "KEK";
  else if( num == 4532) 
    return "KEM";
  else if( num == 4533) 
    return "KEF";
  else if( num == 4534) 
    return "KEP";
  else if( num == 4535) 
    return "KES";
  else if( num == 4536) 
    return "KET";
  else if( num == 4537) 
    return "KEW";
  else if( num == 4538) 
    return "KEY";
  else if( num == 4539) 
    return "KEV";
  else if( num == 4540) 
    return "KGA";
  else if( num == 4541) 
    return "KGR";
  else if( num == 4542) 
    return "KGN";
  else if( num == 4543) 
    return "KGD";
  else if( num == 4544) 
    return "KGC";
  else if( num == 4545) 
    return "KGQ";
  else if( num == 4546) 
    return "KGE";
  else if( num == 4547) 
    return "KGG";
  else if( num == 4548) 
    return "KGH";
  else if( num == 4549) 
    return "KGI";
  else if( num == 4550) 
    return "KGL";
  else if( num == 4551) 
    return "KGK";
  else if( num == 4552) 
    return "KGM";
  else if( num == 4553) 
    return "KGF";
  else if( num == 4554) 
    return "KGP";
  else if( num == 4555) 
    return "KGS";
  else if( num == 4556) 
    return "KGT";
  else if( num == 4557) 
    return "KGW";
  else if( num == 4558) 
    return "KGY";
  else if( num == 4559) 
    return "KGV";
  else if( num == 4560) 
    return "KHA";
  else if( num == 4561) 
    return "KHR";
  else if( num == 4562) 
    return "KHN";
  else if( num == 4563) 
    return "KHD";
  else if( num == 4564) 
    return "KHC";
  else if( num == 4565) 
    return "KHQ";
  else if( num == 4566) 
    return "KHE";
  else if( num == 4567) 
    return "KHG";
  else if( num == 4568) 
    return "KHH";
  else if( num == 4569) 
    return "KHI";
  else if( num == 4570) 
    return "KHL";
  else if( num == 4571) 
    return "KHK";
  else if( num == 4572) 
    return "KHM";
  else if( num == 4573) 
    return "KHF";
  else if( num == 4574) 
    return "KHP";
  else if( num == 4575) 
    return "KHS";
  else if( num == 4576) 
    return "KHT";
  else if( num == 4577) 
    return "KHW";
  else if( num == 4578) 
    return "KHY";
  else if( num == 4579) 
    return "KHV";
  else if( num == 4580) 
    return "KIA";
  else if( num == 4581) 
    return "KIR";
  else if( num == 4582) 
    return "KIN";
  else if( num == 4583) 
    return "KID";
  else if( num == 4584) 
    return "KIC";
  else if( num == 4585) 
    return "KIQ";
  else if( num == 4586) 
    return "KIE";
  else if( num == 4587) 
    return "KIG";
  else if( num == 4588) 
    return "KIH";
  else if( num == 4589) 
    return "KII";
  else if( num == 4590) 
    return "KIL";
  else if( num == 4591) 
    return "KIK";
  else if( num == 4592) 
    return "KIM";
  else if( num == 4593) 
    return "KIF";
  else if( num == 4594) 
    return "KIP";
  else if( num == 4595) 
    return "KIS";
  else if( num == 4596) 
    return "KIT";
  else if( num == 4597) 
    return "KIW";
  else if( num == 4598) 
    return "KIY";
  else if( num == 4599) 
    return "KIV";
  else if( num == 4600) 
    return "KLA";
  else if( num == 4601) 
    return "KLR";
  else if( num == 4602) 
    return "KLN";
  else if( num == 4603) 
    return "KLD";
  else if( num == 4604) 
    return "KLC";
  else if( num == 4605) 
    return "KLQ";
  else if( num == 4606) 
    return "KLE";
  else if( num == 4607) 
    return "KLG";
  else if( num == 4608) 
    return "KLH";
  else if( num == 4609) 
    return "KLI";
  else if( num == 4610) 
    return "KLL";
  else if( num == 4611) 
    return "KLK";
  else if( num == 4612) 
    return "KLM";
  else if( num == 4613) 
    return "KLF";
  else if( num == 4614) 
    return "KLP";
  else if( num == 4615) 
    return "KLS";
  else if( num == 4616) 
    return "KLT";
  else if( num == 4617) 
    return "KLW";
  else if( num == 4618) 
    return "KLY";
  else if( num == 4619) 
    return "KLV";
  else if( num == 4620) 
    return "KKA";
  else if( num == 4621) 
    return "KKR";
  else if( num == 4622) 
    return "KKN";
  else if( num == 4623) 
    return "KKD";
  else if( num == 4624) 
    return "KKC";
  else if( num == 4625) 
    return "KKQ";
  else if( num == 4626) 
    return "KKE";
  else if( num == 4627) 
    return "KKG";
  else if( num == 4628) 
    return "KKH";
  else if( num == 4629) 
    return "KKI";
  else if( num == 4630) 
    return "KKL";
  else if( num == 4631) 
    return "KKK";
  else if( num == 4632) 
    return "KKM";
  else if( num == 4633) 
    return "KKF";
  else if( num == 4634) 
    return "KKP";
  else if( num == 4635) 
    return "KKS";
  else if( num == 4636) 
    return "KKT";
  else if( num == 4637) 
    return "KKW";
  else if( num == 4638) 
    return "KKY";
  else if( num == 4639) 
    return "KKV";
  else if( num == 4640) 
    return "KMA";
  else if( num == 4641) 
    return "KMR";
  else if( num == 4642) 
    return "KMN";
  else if( num == 4643) 
    return "KMD";
  else if( num == 4644) 
    return "KMC";
  else if( num == 4645) 
    return "KMQ";
  else if( num == 4646) 
    return "KME";
  else if( num == 4647) 
    return "KMG";
  else if( num == 4648) 
    return "KMH";
  else if( num == 4649) 
    return "KMI";
  else if( num == 4650) 
    return "KML";
  else if( num == 4651) 
    return "KMK";
  else if( num == 4652) 
    return "KMM";
  else if( num == 4653) 
    return "KMF";
  else if( num == 4654) 
    return "KMP";
  else if( num == 4655) 
    return "KMS";
  else if( num == 4656) 
    return "KMT";
  else if( num == 4657) 
    return "KMW";
  else if( num == 4658) 
    return "KMY";
  else if( num == 4659) 
    return "KMV";
  else if( num == 4660) 
    return "KFA";
  else if( num == 4661) 
    return "KFR";
  else if( num == 4662) 
    return "KFN";
  else if( num == 4663) 
    return "KFD";
  else if( num == 4664) 
    return "KFC";
  else if( num == 4665) 
    return "KFQ";
  else if( num == 4666) 
    return "KFE";
  else if( num == 4667) 
    return "KFG";
  else if( num == 4668) 
    return "KFH";
  else if( num == 4669) 
    return "KFI";
  else if( num == 4670) 
    return "KFL";
  else if( num == 4671) 
    return "KFK";
  else if( num == 4672) 
    return "KFM";
  else if( num == 4673) 
    return "KFF";
  else if( num == 4674) 
    return "KFP";
  else if( num == 4675) 
    return "KFS";
  else if( num == 4676) 
    return "KFT";
  else if( num == 4677) 
    return "KFW";
  else if( num == 4678) 
    return "KFY";
  else if( num == 4679) 
    return "KFV";
  else if( num == 4680) 
    return "KPA";
  else if( num == 4681) 
    return "KPR";
  else if( num == 4682) 
    return "KPN";
  else if( num == 4683) 
    return "KPD";
  else if( num == 4684) 
    return "KPC";
  else if( num == 4685) 
    return "KPQ";
  else if( num == 4686) 
    return "KPE";
  else if( num == 4687) 
    return "KPG";
  else if( num == 4688) 
    return "KPH";
  else if( num == 4689) 
    return "KPI";
  else if( num == 4690) 
    return "KPL";
  else if( num == 4691) 
    return "KPK";
  else if( num == 4692) 
    return "KPM";
  else if( num == 4693) 
    return "KPF";
  else if( num == 4694) 
    return "KPP";
  else if( num == 4695) 
    return "KPS";
  else if( num == 4696) 
    return "KPT";
  else if( num == 4697) 
    return "KPW";
  else if( num == 4698) 
    return "KPY";
  else if( num == 4699) 
    return "KPV";
  else if( num == 4700) 
    return "KSA";
  else if( num == 4701) 
    return "KSR";
  else if( num == 4702) 
    return "KSN";
  else if( num == 4703) 
    return "KSD";
  else if( num == 4704) 
    return "KSC";
  else if( num == 4705) 
    return "KSQ";
  else if( num == 4706) 
    return "KSE";
  else if( num == 4707) 
    return "KSG";
  else if( num == 4708) 
    return "KSH";
  else if( num == 4709) 
    return "KSI";
  else if( num == 4710) 
    return "KSL";
  else if( num == 4711) 
    return "KSK";
  else if( num == 4712) 
    return "KSM";
  else if( num == 4713) 
    return "KSF";
  else if( num == 4714) 
    return "KSP";
  else if( num == 4715) 
    return "KSS";
  else if( num == 4716) 
    return "KST";
  else if( num == 4717) 
    return "KSW";
  else if( num == 4718) 
    return "KSY";
  else if( num == 4719) 
    return "KSV";
  else if( num == 4720) 
    return "KTA";
  else if( num == 4721) 
    return "KTR";
  else if( num == 4722) 
    return "KTN";
  else if( num == 4723) 
    return "KTD";
  else if( num == 4724) 
    return "KTC";
  else if( num == 4725) 
    return "KTQ";
  else if( num == 4726) 
    return "KTE";
  else if( num == 4727) 
    return "KTG";
  else if( num == 4728) 
    return "KTH";
  else if( num == 4729) 
    return "KTI";
  else if( num == 4730) 
    return "KTL";
  else if( num == 4731) 
    return "KTK";
  else if( num == 4732) 
    return "KTM";
  else if( num == 4733) 
    return "KTF";
  else if( num == 4734) 
    return "KTP";
  else if( num == 4735) 
    return "KTS";
  else if( num == 4736) 
    return "KTT";
  else if( num == 4737) 
    return "KTW";
  else if( num == 4738) 
    return "KTY";
  else if( num == 4739) 
    return "KTV";
  else if( num == 4740) 
    return "KWA";
  else if( num == 4741) 
    return "KWR";
  else if( num == 4742) 
    return "KWN";
  else if( num == 4743) 
    return "KWD";
  else if( num == 4744) 
    return "KWC";
  else if( num == 4745) 
    return "KWQ";
  else if( num == 4746) 
    return "KWE";
  else if( num == 4747) 
    return "KWG";
  else if( num == 4748) 
    return "KWH";
  else if( num == 4749) 
    return "KWI";
  else if( num == 4750) 
    return "KWL";
  else if( num == 4751) 
    return "KWK";
  else if( num == 4752) 
    return "KWM";
  else if( num == 4753) 
    return "KWF";
  else if( num == 4754) 
    return "KWP";
  else if( num == 4755) 
    return "KWS";
  else if( num == 4756) 
    return "KWT";
  else if( num == 4757) 
    return "KWW";
  else if( num == 4758) 
    return "KWY";
  else if( num == 4759) 
    return "KWV";
  else if( num == 4760) 
    return "KYA";
  else if( num == 4761) 
    return "KYR";
  else if( num == 4762) 
    return "KYN";
  else if( num == 4763) 
    return "KYD";
  else if( num == 4764) 
    return "KYC";
  else if( num == 4765) 
    return "KYQ";
  else if( num == 4766) 
    return "KYE";
  else if( num == 4767) 
    return "KYG";
  else if( num == 4768) 
    return "KYH";
  else if( num == 4769) 
    return "KYI";
  else if( num == 4770) 
    return "KYL";
  else if( num == 4771) 
    return "KYK";
  else if( num == 4772) 
    return "KYM";
  else if( num == 4773) 
    return "KYF";
  else if( num == 4774) 
    return "KYP";
  else if( num == 4775) 
    return "KYS";
  else if( num == 4776) 
    return "KYT";
  else if( num == 4777) 
    return "KYW";
  else if( num == 4778) 
    return "KYY";
  else if( num == 4779) 
    return "KYV";
  else if( num == 4780) 
    return "KVA";
  else if( num == 4781) 
    return "KVR";
  else if( num == 4782) 
    return "KVN";
  else if( num == 4783) 
    return "KVD";
  else if( num == 4784) 
    return "KVC";
  else if( num == 4785) 
    return "KVQ";
  else if( num == 4786) 
    return "KVE";
  else if( num == 4787) 
    return "KVG";
  else if( num == 4788) 
    return "KVH";
  else if( num == 4789) 
    return "KVI";
  else if( num == 4790) 
    return "KVL";
  else if( num == 4791) 
    return "KVK";
  else if( num == 4792) 
    return "KVM";
  else if( num == 4793) 
    return "KVF";
  else if( num == 4794) 
    return "KVP";
  else if( num == 4795) 
    return "KVS";
  else if( num == 4796) 
    return "KVT";
  else if( num == 4797) 
    return "KVW";
  else if( num == 4798) 
    return "KVY";
  else if( num == 4799) 
    return "KVV";
  else if( num == 4800) 
    return "MAA";
  else if( num == 4801) 
    return "MAR";
  else if( num == 4802) 
    return "MAN";
  else if( num == 4803) 
    return "MAD";
  else if( num == 4804) 
    return "MAC";
  else if( num == 4805) 
    return "MAQ";
  else if( num == 4806) 
    return "MAE";
  else if( num == 4807) 
    return "MAG";
  else if( num == 4808) 
    return "MAH";
  else if( num == 4809) 
    return "MAI";
  else if( num == 4810) 
    return "MAL";
  else if( num == 4811) 
    return "MAK";
  else if( num == 4812) 
    return "MAM";
  else if( num == 4813) 
    return "MAF";
  else if( num == 4814) 
    return "MAP";
  else if( num == 4815) 
    return "MAS";
  else if( num == 4816) 
    return "MAT";
  else if( num == 4817) 
    return "MAW";
  else if( num == 4818) 
    return "MAY";
  else if( num == 4819) 
    return "MAV";
  else if( num == 4820) 
    return "MRA";
  else if( num == 4821) 
    return "MRR";
  else if( num == 4822) 
    return "MRN";
  else if( num == 4823) 
    return "MRD";
  else if( num == 4824) 
    return "MRC";
  else if( num == 4825) 
    return "MRQ";
  else if( num == 4826) 
    return "MRE";
  else if( num == 4827) 
    return "MRG";
  else if( num == 4828) 
    return "MRH";
  else if( num == 4829) 
    return "MRI";
  else if( num == 4830) 
    return "MRL";
  else if( num == 4831) 
    return "MRK";
  else if( num == 4832) 
    return "MRM";
  else if( num == 4833) 
    return "MRF";
  else if( num == 4834) 
    return "MRP";
  else if( num == 4835) 
    return "MRS";
  else if( num == 4836) 
    return "MRT";
  else if( num == 4837) 
    return "MRW";
  else if( num == 4838) 
    return "MRY";
  else if( num == 4839) 
    return "MRV";
  else if( num == 4840) 
    return "MNA";
  else if( num == 4841) 
    return "MNR";
  else if( num == 4842) 
    return "MNN";
  else if( num == 4843) 
    return "MND";
  else if( num == 4844) 
    return "MNC";
  else if( num == 4845) 
    return "MNQ";
  else if( num == 4846) 
    return "MNE";
  else if( num == 4847) 
    return "MNG";
  else if( num == 4848) 
    return "MNH";
  else if( num == 4849) 
    return "MNI";
  else if( num == 4850) 
    return "MNL";
  else if( num == 4851) 
    return "MNK";
  else if( num == 4852) 
    return "MNM";
  else if( num == 4853) 
    return "MNF";
  else if( num == 4854) 
    return "MNP";
  else if( num == 4855) 
    return "MNS";
  else if( num == 4856) 
    return "MNT";
  else if( num == 4857) 
    return "MNW";
  else if( num == 4858) 
    return "MNY";
  else if( num == 4859) 
    return "MNV";
  else if( num == 4860) 
    return "MDA";
  else if( num == 4861) 
    return "MDR";
  else if( num == 4862) 
    return "MDN";
  else if( num == 4863) 
    return "MDD";
  else if( num == 4864) 
    return "MDC";
  else if( num == 4865) 
    return "MDQ";
  else if( num == 4866) 
    return "MDE";
  else if( num == 4867) 
    return "MDG";
  else if( num == 4868) 
    return "MDH";
  else if( num == 4869) 
    return "MDI";
  else if( num == 4870) 
    return "MDL";
  else if( num == 4871) 
    return "MDK";
  else if( num == 4872) 
    return "MDM";
  else if( num == 4873) 
    return "MDF";
  else if( num == 4874) 
    return "MDP";
  else if( num == 4875) 
    return "MDS";
  else if( num == 4876) 
    return "MDT";
  else if( num == 4877) 
    return "MDW";
  else if( num == 4878) 
    return "MDY";
  else if( num == 4879) 
    return "MDV";
  else if( num == 4880) 
    return "MCA";
  else if( num == 4881) 
    return "MCR";
  else if( num == 4882) 
    return "MCN";
  else if( num == 4883) 
    return "MCD";
  else if( num == 4884) 
    return "MCC";
  else if( num == 4885) 
    return "MCQ";
  else if( num == 4886) 
    return "MCE";
  else if( num == 4887) 
    return "MCG";
  else if( num == 4888) 
    return "MCH";
  else if( num == 4889) 
    return "MCI";
  else if( num == 4890) 
    return "MCL";
  else if( num == 4891) 
    return "MCK";
  else if( num == 4892) 
    return "MCM";
  else if( num == 4893) 
    return "MCF";
  else if( num == 4894) 
    return "MCP";
  else if( num == 4895) 
    return "MCS";
  else if( num == 4896) 
    return "MCT";
  else if( num == 4897) 
    return "MCW";
  else if( num == 4898) 
    return "MCY";
  else if( num == 4899) 
    return "MCV";
  else if( num == 4900) 
    return "MQA";
  else if( num == 4901) 
    return "MQR";
  else if( num == 4902) 
    return "MQN";
  else if( num == 4903) 
    return "MQD";
  else if( num == 4904) 
    return "MQC";
  else if( num == 4905) 
    return "MQQ";
  else if( num == 4906) 
    return "MQE";
  else if( num == 4907) 
    return "MQG";
  else if( num == 4908) 
    return "MQH";
  else if( num == 4909) 
    return "MQI";
  else if( num == 4910) 
    return "MQL";
  else if( num == 4911) 
    return "MQK";
  else if( num == 4912) 
    return "MQM";
  else if( num == 4913) 
    return "MQF";
  else if( num == 4914) 
    return "MQP";
  else if( num == 4915) 
    return "MQS";
  else if( num == 4916) 
    return "MQT";
  else if( num == 4917) 
    return "MQW";
  else if( num == 4918) 
    return "MQY";
  else if( num == 4919) 
    return "MQV";
  else if( num == 4920) 
    return "MEA";
  else if( num == 4921) 
    return "MER";
  else if( num == 4922) 
    return "MEN";
  else if( num == 4923) 
    return "MED";
  else if( num == 4924) 
    return "MEC";
  else if( num == 4925) 
    return "MEQ";
  else if( num == 4926) 
    return "MEE";
  else if( num == 4927) 
    return "MEG";
  else if( num == 4928) 
    return "MEH";
  else if( num == 4929) 
    return "MEI";
  else if( num == 4930) 
    return "MEL";
  else if( num == 4931) 
    return "MEK";
  else if( num == 4932) 
    return "MEM";
  else if( num == 4933) 
    return "MEF";
  else if( num == 4934) 
    return "MEP";
  else if( num == 4935) 
    return "MES";
  else if( num == 4936) 
    return "MET";
  else if( num == 4937) 
    return "MEW";
  else if( num == 4938) 
    return "MEY";
  else if( num == 4939) 
    return "MEV";
  else if( num == 4940) 
    return "MGA";
  else if( num == 4941) 
    return "MGR";
  else if( num == 4942) 
    return "MGN";
  else if( num == 4943) 
    return "MGD";
  else if( num == 4944) 
    return "MGC";
  else if( num == 4945) 
    return "MGQ";
  else if( num == 4946) 
    return "MGE";
  else if( num == 4947) 
    return "MGG";
  else if( num == 4948) 
    return "MGH";
  else if( num == 4949) 
    return "MGI";
  else if( num == 4950) 
    return "MGL";
  else if( num == 4951) 
    return "MGK";
  else if( num == 4952) 
    return "MGM";
  else if( num == 4953) 
    return "MGF";
  else if( num == 4954) 
    return "MGP";
  else if( num == 4955) 
    return "MGS";
  else if( num == 4956) 
    return "MGT";
  else if( num == 4957) 
    return "MGW";
  else if( num == 4958) 
    return "MGY";
  else if( num == 4959) 
    return "MGV";
  else if( num == 4960) 
    return "MHA";
  else if( num == 4961) 
    return "MHR";
  else if( num == 4962) 
    return "MHN";
  else if( num == 4963) 
    return "MHD";
  else if( num == 4964) 
    return "MHC";
  else if( num == 4965) 
    return "MHQ";
  else if( num == 4966) 
    return "MHE";
  else if( num == 4967) 
    return "MHG";
  else if( num == 4968) 
    return "MHH";
  else if( num == 4969) 
    return "MHI";
  else if( num == 4970) 
    return "MHL";
  else if( num == 4971) 
    return "MHK";
  else if( num == 4972) 
    return "MHM";
  else if( num == 4973) 
    return "MHF";
  else if( num == 4974) 
    return "MHP";
  else if( num == 4975) 
    return "MHS";
  else if( num == 4976) 
    return "MHT";
  else if( num == 4977) 
    return "MHW";
  else if( num == 4978) 
    return "MHY";
  else if( num == 4979) 
    return "MHV";
  else if( num == 4980) 
    return "MIA";
  else if( num == 4981) 
    return "MIR";
  else if( num == 4982) 
    return "MIN";
  else if( num == 4983) 
    return "MID";
  else if( num == 4984) 
    return "MIC";
  else if( num == 4985) 
    return "MIQ";
  else if( num == 4986) 
    return "MIE";
  else if( num == 4987) 
    return "MIG";
  else if( num == 4988) 
    return "MIH";
  else if( num == 4989) 
    return "MII";
  else if( num == 4990) 
    return "MIL";
  else if( num == 4991) 
    return "MIK";
  else if( num == 4992) 
    return "MIM";
  else if( num == 4993) 
    return "MIF";
  else if( num == 4994) 
    return "MIP";
  else if( num == 4995) 
    return "MIS";
  else if( num == 4996) 
    return "MIT";
  else if( num == 4997) 
    return "MIW";
  else if( num == 4998) 
    return "MIY";
  else if( num == 4999) 
    return "MIV";
  else if( num == 5000) 
    return "MLA";
  else if( num == 5001) 
    return "MLR";
  else if( num == 5002) 
    return "MLN";
  else if( num == 5003) 
    return "MLD";
  else if( num == 5004) 
    return "MLC";
  else if( num == 5005) 
    return "MLQ";
  else if( num == 5006) 
    return "MLE";
  else if( num == 5007) 
    return "MLG";
  else if( num == 5008) 
    return "MLH";
  else if( num == 5009) 
    return "MLI";
  else if( num == 5010) 
    return "MLL";
  else if( num == 5011) 
    return "MLK";
  else if( num == 5012) 
    return "MLM";
  else if( num == 5013) 
    return "MLF";
  else if( num == 5014) 
    return "MLP";
  else if( num == 5015) 
    return "MLS";
  else if( num == 5016) 
    return "MLT";
  else if( num == 5017) 
    return "MLW";
  else if( num == 5018) 
    return "MLY";
  else if( num == 5019) 
    return "MLV";
  else if( num == 5020) 
    return "MKA";
  else if( num == 5021) 
    return "MKR";
  else if( num == 5022) 
    return "MKN";
  else if( num == 5023) 
    return "MKD";
  else if( num == 5024) 
    return "MKC";
  else if( num == 5025) 
    return "MKQ";
  else if( num == 5026) 
    return "MKE";
  else if( num == 5027) 
    return "MKG";
  else if( num == 5028) 
    return "MKH";
  else if( num == 5029) 
    return "MKI";
  else if( num == 5030) 
    return "MKL";
  else if( num == 5031) 
    return "MKK";
  else if( num == 5032) 
    return "MKM";
  else if( num == 5033) 
    return "MKF";
  else if( num == 5034) 
    return "MKP";
  else if( num == 5035) 
    return "MKS";
  else if( num == 5036) 
    return "MKT";
  else if( num == 5037) 
    return "MKW";
  else if( num == 5038) 
    return "MKY";
  else if( num == 5039) 
    return "MKV";
  else if( num == 5040) 
    return "MMA";
  else if( num == 5041) 
    return "MMR";
  else if( num == 5042) 
    return "MMN";
  else if( num == 5043) 
    return "MMD";
  else if( num == 5044) 
    return "MMC";
  else if( num == 5045) 
    return "MMQ";
  else if( num == 5046) 
    return "MME";
  else if( num == 5047) 
    return "MMG";
  else if( num == 5048) 
    return "MMH";
  else if( num == 5049) 
    return "MMI";
  else if( num == 5050) 
    return "MML";
  else if( num == 5051) 
    return "MMK";
  else if( num == 5052) 
    return "MMM";
  else if( num == 5053) 
    return "MMF";
  else if( num == 5054) 
    return "MMP";
  else if( num == 5055) 
    return "MMS";
  else if( num == 5056) 
    return "MMT";
  else if( num == 5057) 
    return "MMW";
  else if( num == 5058) 
    return "MMY";
  else if( num == 5059) 
    return "MMV";
  else if( num == 5060) 
    return "MFA";
  else if( num == 5061) 
    return "MFR";
  else if( num == 5062) 
    return "MFN";
  else if( num == 5063) 
    return "MFD";
  else if( num == 5064) 
    return "MFC";
  else if( num == 5065) 
    return "MFQ";
  else if( num == 5066) 
    return "MFE";
  else if( num == 5067) 
    return "MFG";
  else if( num == 5068) 
    return "MFH";
  else if( num == 5069) 
    return "MFI";
  else if( num == 5070) 
    return "MFL";
  else if( num == 5071) 
    return "MFK";
  else if( num == 5072) 
    return "MFM";
  else if( num == 5073) 
    return "MFF";
  else if( num == 5074) 
    return "MFP";
  else if( num == 5075) 
    return "MFS";
  else if( num == 5076) 
    return "MFT";
  else if( num == 5077) 
    return "MFW";
  else if( num == 5078) 
    return "MFY";
  else if( num == 5079) 
    return "MFV";
  else if( num == 5080) 
    return "MPA";
  else if( num == 5081) 
    return "MPR";
  else if( num == 5082) 
    return "MPN";
  else if( num == 5083) 
    return "MPD";
  else if( num == 5084) 
    return "MPC";
  else if( num == 5085) 
    return "MPQ";
  else if( num == 5086) 
    return "MPE";
  else if( num == 5087) 
    return "MPG";
  else if( num == 5088) 
    return "MPH";
  else if( num == 5089) 
    return "MPI";
  else if( num == 5090) 
    return "MPL";
  else if( num == 5091) 
    return "MPK";
  else if( num == 5092) 
    return "MPM";
  else if( num == 5093) 
    return "MPF";
  else if( num == 5094) 
    return "MPP";
  else if( num == 5095) 
    return "MPS";
  else if( num == 5096) 
    return "MPT";
  else if( num == 5097) 
    return "MPW";
  else if( num == 5098) 
    return "MPY";
  else if( num == 5099) 
    return "MPV";
  else if( num == 5100) 
    return "MSA";
  else if( num == 5101) 
    return "MSR";
  else if( num == 5102) 
    return "MSN";
  else if( num == 5103) 
    return "MSD";
  else if( num == 5104) 
    return "MSC";
  else if( num == 5105) 
    return "MSQ";
  else if( num == 5106) 
    return "MSE";
  else if( num == 5107) 
    return "MSG";
  else if( num == 5108) 
    return "MSH";
  else if( num == 5109) 
    return "MSI";
  else if( num == 5110) 
    return "MSL";
  else if( num == 5111) 
    return "MSK";
  else if( num == 5112) 
    return "MSM";
  else if( num == 5113) 
    return "MSF";
  else if( num == 5114) 
    return "MSP";
  else if( num == 5115) 
    return "MSS";
  else if( num == 5116) 
    return "MST";
  else if( num == 5117) 
    return "MSW";
  else if( num == 5118) 
    return "MSY";
  else if( num == 5119) 
    return "MSV";
  else if( num == 5120) 
    return "MTA";
  else if( num == 5121) 
    return "MTR";
  else if( num == 5122) 
    return "MTN";
  else if( num == 5123) 
    return "MTD";
  else if( num == 5124) 
    return "MTC";
  else if( num == 5125) 
    return "MTQ";
  else if( num == 5126) 
    return "MTE";
  else if( num == 5127) 
    return "MTG";
  else if( num == 5128) 
    return "MTH";
  else if( num == 5129) 
    return "MTI";
  else if( num == 5130) 
    return "MTL";
  else if( num == 5131) 
    return "MTK";
  else if( num == 5132) 
    return "MTM";
  else if( num == 5133) 
    return "MTF";
  else if( num == 5134) 
    return "MTP";
  else if( num == 5135) 
    return "MTS";
  else if( num == 5136) 
    return "MTT";
  else if( num == 5137) 
    return "MTW";
  else if( num == 5138) 
    return "MTY";
  else if( num == 5139) 
    return "MTV";
  else if( num == 5140) 
    return "MWA";
  else if( num == 5141) 
    return "MWR";
  else if( num == 5142) 
    return "MWN";
  else if( num == 5143) 
    return "MWD";
  else if( num == 5144) 
    return "MWC";
  else if( num == 5145) 
    return "MWQ";
  else if( num == 5146) 
    return "MWE";
  else if( num == 5147) 
    return "MWG";
  else if( num == 5148) 
    return "MWH";
  else if( num == 5149) 
    return "MWI";
  else if( num == 5150) 
    return "MWL";
  else if( num == 5151) 
    return "MWK";
  else if( num == 5152) 
    return "MWM";
  else if( num == 5153) 
    return "MWF";
  else if( num == 5154) 
    return "MWP";
  else if( num == 5155) 
    return "MWS";
  else if( num == 5156) 
    return "MWT";
  else if( num == 5157) 
    return "MWW";
  else if( num == 5158) 
    return "MWY";
  else if( num == 5159) 
    return "MWV";
  else if( num == 5160) 
    return "MYA";
  else if( num == 5161) 
    return "MYR";
  else if( num == 5162) 
    return "MYN";
  else if( num == 5163) 
    return "MYD";
  else if( num == 5164) 
    return "MYC";
  else if( num == 5165) 
    return "MYQ";
  else if( num == 5166) 
    return "MYE";
  else if( num == 5167) 
    return "MYG";
  else if( num == 5168) 
    return "MYH";
  else if( num == 5169) 
    return "MYI";
  else if( num == 5170) 
    return "MYL";
  else if( num == 5171) 
    return "MYK";
  else if( num == 5172) 
    return "MYM";
  else if( num == 5173) 
    return "MYF";
  else if( num == 5174) 
    return "MYP";
  else if( num == 5175) 
    return "MYS";
  else if( num == 5176) 
    return "MYT";
  else if( num == 5177) 
    return "MYW";
  else if( num == 5178) 
    return "MYY";
  else if( num == 5179) 
    return "MYV";
  else if( num == 5180) 
    return "MVA";
  else if( num == 5181) 
    return "MVR";
  else if( num == 5182) 
    return "MVN";
  else if( num == 5183) 
    return "MVD";
  else if( num == 5184) 
    return "MVC";
  else if( num == 5185) 
    return "MVQ";
  else if( num == 5186) 
    return "MVE";
  else if( num == 5187) 
    return "MVG";
  else if( num == 5188) 
    return "MVH";
  else if( num == 5189) 
    return "MVI";
  else if( num == 5190) 
    return "MVL";
  else if( num == 5191) 
    return "MVK";
  else if( num == 5192) 
    return "MVM";
  else if( num == 5193) 
    return "MVF";
  else if( num == 5194) 
    return "MVP";
  else if( num == 5195) 
    return "MVS";
  else if( num == 5196) 
    return "MVT";
  else if( num == 5197) 
    return "MVW";
  else if( num == 5198) 
    return "MVY";
  else if( num == 5199) 
    return "MVV";
  else if( num == 5200) 
    return "FAA";
  else if( num == 5201) 
    return "FAR";
  else if( num == 5202) 
    return "FAN";
  else if( num == 5203) 
    return "FAD";
  else if( num == 5204) 
    return "FAC";
  else if( num == 5205) 
    return "FAQ";
  else if( num == 5206) 
    return "FAE";
  else if( num == 5207) 
    return "FAG";
  else if( num == 5208) 
    return "FAH";
  else if( num == 5209) 
    return "FAI";
  else if( num == 5210) 
    return "FAL";
  else if( num == 5211) 
    return "FAK";
  else if( num == 5212) 
    return "FAM";
  else if( num == 5213) 
    return "FAF";
  else if( num == 5214) 
    return "FAP";
  else if( num == 5215) 
    return "FAS";
  else if( num == 5216) 
    return "FAT";
  else if( num == 5217) 
    return "FAW";
  else if( num == 5218) 
    return "FAY";
  else if( num == 5219) 
    return "FAV";
  else if( num == 5220) 
    return "FRA";
  else if( num == 5221) 
    return "FRR";
  else if( num == 5222) 
    return "FRN";
  else if( num == 5223) 
    return "FRD";
  else if( num == 5224) 
    return "FRC";
  else if( num == 5225) 
    return "FRQ";
  else if( num == 5226) 
    return "FRE";
  else if( num == 5227) 
    return "FRG";
  else if( num == 5228) 
    return "FRH";
  else if( num == 5229) 
    return "FRI";
  else if( num == 5230) 
    return "FRL";
  else if( num == 5231) 
    return "FRK";
  else if( num == 5232) 
    return "FRM";
  else if( num == 5233) 
    return "FRF";
  else if( num == 5234) 
    return "FRP";
  else if( num == 5235) 
    return "FRS";
  else if( num == 5236) 
    return "FRT";
  else if( num == 5237) 
    return "FRW";
  else if( num == 5238) 
    return "FRY";
  else if( num == 5239) 
    return "FRV";
  else if( num == 5240) 
    return "FNA";
  else if( num == 5241) 
    return "FNR";
  else if( num == 5242) 
    return "FNN";
  else if( num == 5243) 
    return "FND";
  else if( num == 5244) 
    return "FNC";
  else if( num == 5245) 
    return "FNQ";
  else if( num == 5246) 
    return "FNE";
  else if( num == 5247) 
    return "FNG";
  else if( num == 5248) 
    return "FNH";
  else if( num == 5249) 
    return "FNI";
  else if( num == 5250) 
    return "FNL";
  else if( num == 5251) 
    return "FNK";
  else if( num == 5252) 
    return "FNM";
  else if( num == 5253) 
    return "FNF";
  else if( num == 5254) 
    return "FNP";
  else if( num == 5255) 
    return "FNS";
  else if( num == 5256) 
    return "FNT";
  else if( num == 5257) 
    return "FNW";
  else if( num == 5258) 
    return "FNY";
  else if( num == 5259) 
    return "FNV";
  else if( num == 5260) 
    return "FDA";
  else if( num == 5261) 
    return "FDR";
  else if( num == 5262) 
    return "FDN";
  else if( num == 5263) 
    return "FDD";
  else if( num == 5264) 
    return "FDC";
  else if( num == 5265) 
    return "FDQ";
  else if( num == 5266) 
    return "FDE";
  else if( num == 5267) 
    return "FDG";
  else if( num == 5268) 
    return "FDH";
  else if( num == 5269) 
    return "FDI";
  else if( num == 5270) 
    return "FDL";
  else if( num == 5271) 
    return "FDK";
  else if( num == 5272) 
    return "FDM";
  else if( num == 5273) 
    return "FDF";
  else if( num == 5274) 
    return "FDP";
  else if( num == 5275) 
    return "FDS";
  else if( num == 5276) 
    return "FDT";
  else if( num == 5277) 
    return "FDW";
  else if( num == 5278) 
    return "FDY";
  else if( num == 5279) 
    return "FDV";
  else if( num == 5280) 
    return "FCA";
  else if( num == 5281) 
    return "FCR";
  else if( num == 5282) 
    return "FCN";
  else if( num == 5283) 
    return "FCD";
  else if( num == 5284) 
    return "FCC";
  else if( num == 5285) 
    return "FCQ";
  else if( num == 5286) 
    return "FCE";
  else if( num == 5287) 
    return "FCG";
  else if( num == 5288) 
    return "FCH";
  else if( num == 5289) 
    return "FCI";
  else if( num == 5290) 
    return "FCL";
  else if( num == 5291) 
    return "FCK";
  else if( num == 5292) 
    return "FCM";
  else if( num == 5293) 
    return "FCF";
  else if( num == 5294) 
    return "FCP";
  else if( num == 5295) 
    return "FCS";
  else if( num == 5296) 
    return "FCT";
  else if( num == 5297) 
    return "FCW";
  else if( num == 5298) 
    return "FCY";
  else if( num == 5299) 
    return "FCV";
  else if( num == 5300) 
    return "FQA";
  else if( num == 5301) 
    return "FQR";
  else if( num == 5302) 
    return "FQN";
  else if( num == 5303) 
    return "FQD";
  else if( num == 5304) 
    return "FQC";
  else if( num == 5305) 
    return "FQQ";
  else if( num == 5306) 
    return "FQE";
  else if( num == 5307) 
    return "FQG";
  else if( num == 5308) 
    return "FQH";
  else if( num == 5309) 
    return "FQI";
  else if( num == 5310) 
    return "FQL";
  else if( num == 5311) 
    return "FQK";
  else if( num == 5312) 
    return "FQM";
  else if( num == 5313) 
    return "FQF";
  else if( num == 5314) 
    return "FQP";
  else if( num == 5315) 
    return "FQS";
  else if( num == 5316) 
    return "FQT";
  else if( num == 5317) 
    return "FQW";
  else if( num == 5318) 
    return "FQY";
  else if( num == 5319) 
    return "FQV";
  else if( num == 5320) 
    return "FEA";
  else if( num == 5321) 
    return "FER";
  else if( num == 5322) 
    return "FEN";
  else if( num == 5323) 
    return "FED";
  else if( num == 5324) 
    return "FEC";
  else if( num == 5325) 
    return "FEQ";
  else if( num == 5326) 
    return "FEE";
  else if( num == 5327) 
    return "FEG";
  else if( num == 5328) 
    return "FEH";
  else if( num == 5329) 
    return "FEI";
  else if( num == 5330) 
    return "FEL";
  else if( num == 5331) 
    return "FEK";
  else if( num == 5332) 
    return "FEM";
  else if( num == 5333) 
    return "FEF";
  else if( num == 5334) 
    return "FEP";
  else if( num == 5335) 
    return "FES";
  else if( num == 5336) 
    return "FET";
  else if( num == 5337) 
    return "FEW";
  else if( num == 5338) 
    return "FEY";
  else if( num == 5339) 
    return "FEV";
  else if( num == 5340) 
    return "FGA";
  else if( num == 5341) 
    return "FGR";
  else if( num == 5342) 
    return "FGN";
  else if( num == 5343) 
    return "FGD";
  else if( num == 5344) 
    return "FGC";
  else if( num == 5345) 
    return "FGQ";
  else if( num == 5346) 
    return "FGE";
  else if( num == 5347) 
    return "FGG";
  else if( num == 5348) 
    return "FGH";
  else if( num == 5349) 
    return "FGI";
  else if( num == 5350) 
    return "FGL";
  else if( num == 5351) 
    return "FGK";
  else if( num == 5352) 
    return "FGM";
  else if( num == 5353) 
    return "FGF";
  else if( num == 5354) 
    return "FGP";
  else if( num == 5355) 
    return "FGS";
  else if( num == 5356) 
    return "FGT";
  else if( num == 5357) 
    return "FGW";
  else if( num == 5358) 
    return "FGY";
  else if( num == 5359) 
    return "FGV";
  else if( num == 5360) 
    return "FHA";
  else if( num == 5361) 
    return "FHR";
  else if( num == 5362) 
    return "FHN";
  else if( num == 5363) 
    return "FHD";
  else if( num == 5364) 
    return "FHC";
  else if( num == 5365) 
    return "FHQ";
  else if( num == 5366) 
    return "FHE";
  else if( num == 5367) 
    return "FHG";
  else if( num == 5368) 
    return "FHH";
  else if( num == 5369) 
    return "FHI";
  else if( num == 5370) 
    return "FHL";
  else if( num == 5371) 
    return "FHK";
  else if( num == 5372) 
    return "FHM";
  else if( num == 5373) 
    return "FHF";
  else if( num == 5374) 
    return "FHP";
  else if( num == 5375) 
    return "FHS";
  else if( num == 5376) 
    return "FHT";
  else if( num == 5377) 
    return "FHW";
  else if( num == 5378) 
    return "FHY";
  else if( num == 5379) 
    return "FHV";
  else if( num == 5380) 
    return "FIA";
  else if( num == 5381) 
    return "FIR";
  else if( num == 5382) 
    return "FIN";
  else if( num == 5383) 
    return "FID";
  else if( num == 5384) 
    return "FIC";
  else if( num == 5385) 
    return "FIQ";
  else if( num == 5386) 
    return "FIE";
  else if( num == 5387) 
    return "FIG";
  else if( num == 5388) 
    return "FIH";
  else if( num == 5389) 
    return "FII";
  else if( num == 5390) 
    return "FIL";
  else if( num == 5391) 
    return "FIK";
  else if( num == 5392) 
    return "FIM";
  else if( num == 5393) 
    return "FIF";
  else if( num == 5394) 
    return "FIP";
  else if( num == 5395) 
    return "FIS";
  else if( num == 5396) 
    return "FIT";
  else if( num == 5397) 
    return "FIW";
  else if( num == 5398) 
    return "FIY";
  else if( num == 5399) 
    return "FIV";
  else if( num == 5400) 
    return "FLA";
  else if( num == 5401) 
    return "FLR";
  else if( num == 5402) 
    return "FLN";
  else if( num == 5403) 
    return "FLD";
  else if( num == 5404) 
    return "FLC";
  else if( num == 5405) 
    return "FLQ";
  else if( num == 5406) 
    return "FLE";
  else if( num == 5407) 
    return "FLG";
  else if( num == 5408) 
    return "FLH";
  else if( num == 5409) 
    return "FLI";
  else if( num == 5410) 
    return "FLL";
  else if( num == 5411) 
    return "FLK";
  else if( num == 5412) 
    return "FLM";
  else if( num == 5413) 
    return "FLF";
  else if( num == 5414) 
    return "FLP";
  else if( num == 5415) 
    return "FLS";
  else if( num == 5416) 
    return "FLT";
  else if( num == 5417) 
    return "FLW";
  else if( num == 5418) 
    return "FLY";
  else if( num == 5419) 
    return "FLV";
  else if( num == 5420) 
    return "FKA";
  else if( num == 5421) 
    return "FKR";
  else if( num == 5422) 
    return "FKN";
  else if( num == 5423) 
    return "FKD";
  else if( num == 5424) 
    return "FKC";
  else if( num == 5425) 
    return "FKQ";
  else if( num == 5426) 
    return "FKE";
  else if( num == 5427) 
    return "FKG";
  else if( num == 5428) 
    return "FKH";
  else if( num == 5429) 
    return "FKI";
  else if( num == 5430) 
    return "FKL";
  else if( num == 5431) 
    return "FKK";
  else if( num == 5432) 
    return "FKM";
  else if( num == 5433) 
    return "FKF";
  else if( num == 5434) 
    return "FKP";
  else if( num == 5435) 
    return "FKS";
  else if( num == 5436) 
    return "FKT";
  else if( num == 5437) 
    return "FKW";
  else if( num == 5438) 
    return "FKY";
  else if( num == 5439) 
    return "FKV";
  else if( num == 5440) 
    return "FMA";
  else if( num == 5441) 
    return "FMR";
  else if( num == 5442) 
    return "FMN";
  else if( num == 5443) 
    return "FMD";
  else if( num == 5444) 
    return "FMC";
  else if( num == 5445) 
    return "FMQ";
  else if( num == 5446) 
    return "FME";
  else if( num == 5447) 
    return "FMG";
  else if( num == 5448) 
    return "FMH";
  else if( num == 5449) 
    return "FMI";
  else if( num == 5450) 
    return "FML";
  else if( num == 5451) 
    return "FMK";
  else if( num == 5452) 
    return "FMM";
  else if( num == 5453) 
    return "FMF";
  else if( num == 5454) 
    return "FMP";
  else if( num == 5455) 
    return "FMS";
  else if( num == 5456) 
    return "FMT";
  else if( num == 5457) 
    return "FMW";
  else if( num == 5458) 
    return "FMY";
  else if( num == 5459) 
    return "FMV";
  else if( num == 5460) 
    return "FFA";
  else if( num == 5461) 
    return "FFR";
  else if( num == 5462) 
    return "FFN";
  else if( num == 5463) 
    return "FFD";
  else if( num == 5464) 
    return "FFC";
  else if( num == 5465) 
    return "FFQ";
  else if( num == 5466) 
    return "FFE";
  else if( num == 5467) 
    return "FFG";
  else if( num == 5468) 
    return "FFH";
  else if( num == 5469) 
    return "FFI";
  else if( num == 5470) 
    return "FFL";
  else if( num == 5471) 
    return "FFK";
  else if( num == 5472) 
    return "FFM";
  else if( num == 5473) 
    return "FFF";
  else if( num == 5474) 
    return "FFP";
  else if( num == 5475) 
    return "FFS";
  else if( num == 5476) 
    return "FFT";
  else if( num == 5477) 
    return "FFW";
  else if( num == 5478) 
    return "FFY";
  else if( num == 5479) 
    return "FFV";
  else if( num == 5480) 
    return "FPA";
  else if( num == 5481) 
    return "FPR";
  else if( num == 5482) 
    return "FPN";
  else if( num == 5483) 
    return "FPD";
  else if( num == 5484) 
    return "FPC";
  else if( num == 5485) 
    return "FPQ";
  else if( num == 5486) 
    return "FPE";
  else if( num == 5487) 
    return "FPG";
  else if( num == 5488) 
    return "FPH";
  else if( num == 5489) 
    return "FPI";
  else if( num == 5490) 
    return "FPL";
  else if( num == 5491) 
    return "FPK";
  else if( num == 5492) 
    return "FPM";
  else if( num == 5493) 
    return "FPF";
  else if( num == 5494) 
    return "FPP";
  else if( num == 5495) 
    return "FPS";
  else if( num == 5496) 
    return "FPT";
  else if( num == 5497) 
    return "FPW";
  else if( num == 5498) 
    return "FPY";
  else if( num == 5499) 
    return "FPV";
  else if( num == 5500) 
    return "FSA";
  else if( num == 5501) 
    return "FSR";
  else if( num == 5502) 
    return "FSN";
  else if( num == 5503) 
    return "FSD";
  else if( num == 5504) 
    return "FSC";
  else if( num == 5505) 
    return "FSQ";
  else if( num == 5506) 
    return "FSE";
  else if( num == 5507) 
    return "FSG";
  else if( num == 5508) 
    return "FSH";
  else if( num == 5509) 
    return "FSI";
  else if( num == 5510) 
    return "FSL";
  else if( num == 5511) 
    return "FSK";
  else if( num == 5512) 
    return "FSM";
  else if( num == 5513) 
    return "FSF";
  else if( num == 5514) 
    return "FSP";
  else if( num == 5515) 
    return "FSS";
  else if( num == 5516) 
    return "FST";
  else if( num == 5517) 
    return "FSW";
  else if( num == 5518) 
    return "FSY";
  else if( num == 5519) 
    return "FSV";
  else if( num == 5520) 
    return "FTA";
  else if( num == 5521) 
    return "FTR";
  else if( num == 5522) 
    return "FTN";
  else if( num == 5523) 
    return "FTD";
  else if( num == 5524) 
    return "FTC";
  else if( num == 5525) 
    return "FTQ";
  else if( num == 5526) 
    return "FTE";
  else if( num == 5527) 
    return "FTG";
  else if( num == 5528) 
    return "FTH";
  else if( num == 5529) 
    return "FTI";
  else if( num == 5530) 
    return "FTL";
  else if( num == 5531) 
    return "FTK";
  else if( num == 5532) 
    return "FTM";
  else if( num == 5533) 
    return "FTF";
  else if( num == 5534) 
    return "FTP";
  else if( num == 5535) 
    return "FTS";
  else if( num == 5536) 
    return "FTT";
  else if( num == 5537) 
    return "FTW";
  else if( num == 5538) 
    return "FTY";
  else if( num == 5539) 
    return "FTV";
  else if( num == 5540) 
    return "FWA";
  else if( num == 5541) 
    return "FWR";
  else if( num == 5542) 
    return "FWN";
  else if( num == 5543) 
    return "FWD";
  else if( num == 5544) 
    return "FWC";
  else if( num == 5545) 
    return "FWQ";
  else if( num == 5546) 
    return "FWE";
  else if( num == 5547) 
    return "FWG";
  else if( num == 5548) 
    return "FWH";
  else if( num == 5549) 
    return "FWI";
  else if( num == 5550) 
    return "FWL";
  else if( num == 5551) 
    return "FWK";
  else if( num == 5552) 
    return "FWM";
  else if( num == 5553) 
    return "FWF";
  else if( num == 5554) 
    return "FWP";
  else if( num == 5555) 
    return "FWS";
  else if( num == 5556) 
    return "FWT";
  else if( num == 5557) 
    return "FWW";
  else if( num == 5558) 
    return "FWY";
  else if( num == 5559) 
    return "FWV";
  else if( num == 5560) 
    return "FYA";
  else if( num == 5561) 
    return "FYR";
  else if( num == 5562) 
    return "FYN";
  else if( num == 5563) 
    return "FYD";
  else if( num == 5564) 
    return "FYC";
  else if( num == 5565) 
    return "FYQ";
  else if( num == 5566) 
    return "FYE";
  else if( num == 5567) 
    return "FYG";
  else if( num == 5568) 
    return "FYH";
  else if( num == 5569) 
    return "FYI";
  else if( num == 5570) 
    return "FYL";
  else if( num == 5571) 
    return "FYK";
  else if( num == 5572) 
    return "FYM";
  else if( num == 5573) 
    return "FYF";
  else if( num == 5574) 
    return "FYP";
  else if( num == 5575) 
    return "FYS";
  else if( num == 5576) 
    return "FYT";
  else if( num == 5577) 
    return "FYW";
  else if( num == 5578) 
    return "FYY";
  else if( num == 5579) 
    return "FYV";
  else if( num == 5580) 
    return "FVA";
  else if( num == 5581) 
    return "FVR";
  else if( num == 5582) 
    return "FVN";
  else if( num == 5583) 
    return "FVD";
  else if( num == 5584) 
    return "FVC";
  else if( num == 5585) 
    return "FVQ";
  else if( num == 5586) 
    return "FVE";
  else if( num == 5587) 
    return "FVG";
  else if( num == 5588) 
    return "FVH";
  else if( num == 5589) 
    return "FVI";
  else if( num == 5590) 
    return "FVL";
  else if( num == 5591) 
    return "FVK";
  else if( num == 5592) 
    return "FVM";
  else if( num == 5593) 
    return "FVF";
  else if( num == 5594) 
    return "FVP";
  else if( num == 5595) 
    return "FVS";
  else if( num == 5596) 
    return "FVT";
  else if( num == 5597) 
    return "FVW";
  else if( num == 5598) 
    return "FVY";
  else if( num == 5599) 
    return "FVV";
  else if( num == 5600) 
    return "PAA";
  else if( num == 5601) 
    return "PAR";
  else if( num == 5602) 
    return "PAN";
  else if( num == 5603) 
    return "PAD";
  else if( num == 5604) 
    return "PAC";
  else if( num == 5605) 
    return "PAQ";
  else if( num == 5606) 
    return "PAE";
  else if( num == 5607) 
    return "PAG";
  else if( num == 5608) 
    return "PAH";
  else if( num == 5609) 
    return "PAI";
  else if( num == 5610) 
    return "PAL";
  else if( num == 5611) 
    return "PAK";
  else if( num == 5612) 
    return "PAM";
  else if( num == 5613) 
    return "PAF";
  else if( num == 5614) 
    return "PAP";
  else if( num == 5615) 
    return "PAS";
  else if( num == 5616) 
    return "PAT";
  else if( num == 5617) 
    return "PAW";
  else if( num == 5618) 
    return "PAY";
  else if( num == 5619) 
    return "PAV";
  else if( num == 5620) 
    return "PRA";
  else if( num == 5621) 
    return "PRR";
  else if( num == 5622) 
    return "PRN";
  else if( num == 5623) 
    return "PRD";
  else if( num == 5624) 
    return "PRC";
  else if( num == 5625) 
    return "PRQ";
  else if( num == 5626) 
    return "PRE";
  else if( num == 5627) 
    return "PRG";
  else if( num == 5628) 
    return "PRH";
  else if( num == 5629) 
    return "PRI";
  else if( num == 5630) 
    return "PRL";
  else if( num == 5631) 
    return "PRK";
  else if( num == 5632) 
    return "PRM";
  else if( num == 5633) 
    return "PRF";
  else if( num == 5634) 
    return "PRP";
  else if( num == 5635) 
    return "PRS";
  else if( num == 5636) 
    return "PRT";
  else if( num == 5637) 
    return "PRW";
  else if( num == 5638) 
    return "PRY";
  else if( num == 5639) 
    return "PRV";
  else if( num == 5640) 
    return "PNA";
  else if( num == 5641) 
    return "PNR";
  else if( num == 5642) 
    return "PNN";
  else if( num == 5643) 
    return "PND";
  else if( num == 5644) 
    return "PNC";
  else if( num == 5645) 
    return "PNQ";
  else if( num == 5646) 
    return "PNE";
  else if( num == 5647) 
    return "PNG";
  else if( num == 5648) 
    return "PNH";
  else if( num == 5649) 
    return "PNI";
  else if( num == 5650) 
    return "PNL";
  else if( num == 5651) 
    return "PNK";
  else if( num == 5652) 
    return "PNM";
  else if( num == 5653) 
    return "PNF";
  else if( num == 5654) 
    return "PNP";
  else if( num == 5655) 
    return "PNS";
  else if( num == 5656) 
    return "PNT";
  else if( num == 5657) 
    return "PNW";
  else if( num == 5658) 
    return "PNY";
  else if( num == 5659) 
    return "PNV";
  else if( num == 5660) 
    return "PDA";
  else if( num == 5661) 
    return "PDR";
  else if( num == 5662) 
    return "PDN";
  else if( num == 5663) 
    return "PDD";
  else if( num == 5664) 
    return "PDC";
  else if( num == 5665) 
    return "PDQ";
  else if( num == 5666) 
    return "PDE";
  else if( num == 5667) 
    return "PDG";
  else if( num == 5668) 
    return "PDH";
  else if( num == 5669) 
    return "PDI";
  else if( num == 5670) 
    return "PDL";
  else if( num == 5671) 
    return "PDK";
  else if( num == 5672) 
    return "PDM";
  else if( num == 5673) 
    return "PDF";
  else if( num == 5674) 
    return "PDP";
  else if( num == 5675) 
    return "PDS";
  else if( num == 5676) 
    return "PDT";
  else if( num == 5677) 
    return "PDW";
  else if( num == 5678) 
    return "PDY";
  else if( num == 5679) 
    return "PDV";
  else if( num == 5680) 
    return "PCA";
  else if( num == 5681) 
    return "PCR";
  else if( num == 5682) 
    return "PCN";
  else if( num == 5683) 
    return "PCD";
  else if( num == 5684) 
    return "PCC";
  else if( num == 5685) 
    return "PCQ";
  else if( num == 5686) 
    return "PCE";
  else if( num == 5687) 
    return "PCG";
  else if( num == 5688) 
    return "PCH";
  else if( num == 5689) 
    return "PCI";
  else if( num == 5690) 
    return "PCL";
  else if( num == 5691) 
    return "PCK";
  else if( num == 5692) 
    return "PCM";
  else if( num == 5693) 
    return "PCF";
  else if( num == 5694) 
    return "PCP";
  else if( num == 5695) 
    return "PCS";
  else if( num == 5696) 
    return "PCT";
  else if( num == 5697) 
    return "PCW";
  else if( num == 5698) 
    return "PCY";
  else if( num == 5699) 
    return "PCV";
  else if( num == 5700) 
    return "PQA";
  else if( num == 5701) 
    return "PQR";
  else if( num == 5702) 
    return "PQN";
  else if( num == 5703) 
    return "PQD";
  else if( num == 5704) 
    return "PQC";
  else if( num == 5705) 
    return "PQQ";
  else if( num == 5706) 
    return "PQE";
  else if( num == 5707) 
    return "PQG";
  else if( num == 5708) 
    return "PQH";
  else if( num == 5709) 
    return "PQI";
  else if( num == 5710) 
    return "PQL";
  else if( num == 5711) 
    return "PQK";
  else if( num == 5712) 
    return "PQM";
  else if( num == 5713) 
    return "PQF";
  else if( num == 5714) 
    return "PQP";
  else if( num == 5715) 
    return "PQS";
  else if( num == 5716) 
    return "PQT";
  else if( num == 5717) 
    return "PQW";
  else if( num == 5718) 
    return "PQY";
  else if( num == 5719) 
    return "PQV";
  else if( num == 5720) 
    return "PEA";
  else if( num == 5721) 
    return "PER";
  else if( num == 5722) 
    return "PEN";
  else if( num == 5723) 
    return "PED";
  else if( num == 5724) 
    return "PEC";
  else if( num == 5725) 
    return "PEQ";
  else if( num == 5726) 
    return "PEE";
  else if( num == 5727) 
    return "PEG";
  else if( num == 5728) 
    return "PEH";
  else if( num == 5729) 
    return "PEI";
  else if( num == 5730) 
    return "PEL";
  else if( num == 5731) 
    return "PEK";
  else if( num == 5732) 
    return "PEM";
  else if( num == 5733) 
    return "PEF";
  else if( num == 5734) 
    return "PEP";
  else if( num == 5735) 
    return "PES";
  else if( num == 5736) 
    return "PET";
  else if( num == 5737) 
    return "PEW";
  else if( num == 5738) 
    return "PEY";
  else if( num == 5739) 
    return "PEV";
  else if( num == 5740) 
    return "PGA";
  else if( num == 5741) 
    return "PGR";
  else if( num == 5742) 
    return "PGN";
  else if( num == 5743) 
    return "PGD";
  else if( num == 5744) 
    return "PGC";
  else if( num == 5745) 
    return "PGQ";
  else if( num == 5746) 
    return "PGE";
  else if( num == 5747) 
    return "PGG";
  else if( num == 5748) 
    return "PGH";
  else if( num == 5749) 
    return "PGI";
  else if( num == 5750) 
    return "PGL";
  else if( num == 5751) 
    return "PGK";
  else if( num == 5752) 
    return "PGM";
  else if( num == 5753) 
    return "PGF";
  else if( num == 5754) 
    return "PGP";
  else if( num == 5755) 
    return "PGS";
  else if( num == 5756) 
    return "PGT";
  else if( num == 5757) 
    return "PGW";
  else if( num == 5758) 
    return "PGY";
  else if( num == 5759) 
    return "PGV";
  else if( num == 5760) 
    return "PHA";
  else if( num == 5761) 
    return "PHR";
  else if( num == 5762) 
    return "PHN";
  else if( num == 5763) 
    return "PHD";
  else if( num == 5764) 
    return "PHC";
  else if( num == 5765) 
    return "PHQ";
  else if( num == 5766) 
    return "PHE";
  else if( num == 5767) 
    return "PHG";
  else if( num == 5768) 
    return "PHH";
  else if( num == 5769) 
    return "PHI";
  else if( num == 5770) 
    return "PHL";
  else if( num == 5771) 
    return "PHK";
  else if( num == 5772) 
    return "PHM";
  else if( num == 5773) 
    return "PHF";
  else if( num == 5774) 
    return "PHP";
  else if( num == 5775) 
    return "PHS";
  else if( num == 5776) 
    return "PHT";
  else if( num == 5777) 
    return "PHW";
  else if( num == 5778) 
    return "PHY";
  else if( num == 5779) 
    return "PHV";
  else if( num == 5780) 
    return "PIA";
  else if( num == 5781) 
    return "PIR";
  else if( num == 5782) 
    return "PIN";
  else if( num == 5783) 
    return "PID";
  else if( num == 5784) 
    return "PIC";
  else if( num == 5785) 
    return "PIQ";
  else if( num == 5786) 
    return "PIE";
  else if( num == 5787) 
    return "PIG";
  else if( num == 5788) 
    return "PIH";
  else if( num == 5789) 
    return "PII";
  else if( num == 5790) 
    return "PIL";
  else if( num == 5791) 
    return "PIK";
  else if( num == 5792) 
    return "PIM";
  else if( num == 5793) 
    return "PIF";
  else if( num == 5794) 
    return "PIP";
  else if( num == 5795) 
    return "PIS";
  else if( num == 5796) 
    return "PIT";
  else if( num == 5797) 
    return "PIW";
  else if( num == 5798) 
    return "PIY";
  else if( num == 5799) 
    return "PIV";
  else if( num == 5800) 
    return "PLA";
  else if( num == 5801) 
    return "PLR";
  else if( num == 5802) 
    return "PLN";
  else if( num == 5803) 
    return "PLD";
  else if( num == 5804) 
    return "PLC";
  else if( num == 5805) 
    return "PLQ";
  else if( num == 5806) 
    return "PLE";
  else if( num == 5807) 
    return "PLG";
  else if( num == 5808) 
    return "PLH";
  else if( num == 5809) 
    return "PLI";
  else if( num == 5810) 
    return "PLL";
  else if( num == 5811) 
    return "PLK";
  else if( num == 5812) 
    return "PLM";
  else if( num == 5813) 
    return "PLF";
  else if( num == 5814) 
    return "PLP";
  else if( num == 5815) 
    return "PLS";
  else if( num == 5816) 
    return "PLT";
  else if( num == 5817) 
    return "PLW";
  else if( num == 5818) 
    return "PLY";
  else if( num == 5819) 
    return "PLV";
  else if( num == 5820) 
    return "PKA";
  else if( num == 5821) 
    return "PKR";
  else if( num == 5822) 
    return "PKN";
  else if( num == 5823) 
    return "PKD";
  else if( num == 5824) 
    return "PKC";
  else if( num == 5825) 
    return "PKQ";
  else if( num == 5826) 
    return "PKE";
  else if( num == 5827) 
    return "PKG";
  else if( num == 5828) 
    return "PKH";
  else if( num == 5829) 
    return "PKI";
  else if( num == 5830) 
    return "PKL";
  else if( num == 5831) 
    return "PKK";
  else if( num == 5832) 
    return "PKM";
  else if( num == 5833) 
    return "PKF";
  else if( num == 5834) 
    return "PKP";
  else if( num == 5835) 
    return "PKS";
  else if( num == 5836) 
    return "PKT";
  else if( num == 5837) 
    return "PKW";
  else if( num == 5838) 
    return "PKY";
  else if( num == 5839) 
    return "PKV";
  else if( num == 5840) 
    return "PMA";
  else if( num == 5841) 
    return "PMR";
  else if( num == 5842) 
    return "PMN";
  else if( num == 5843) 
    return "PMD";
  else if( num == 5844) 
    return "PMC";
  else if( num == 5845) 
    return "PMQ";
  else if( num == 5846) 
    return "PME";
  else if( num == 5847) 
    return "PMG";
  else if( num == 5848) 
    return "PMH";
  else if( num == 5849) 
    return "PMI";
  else if( num == 5850) 
    return "PML";
  else if( num == 5851) 
    return "PMK";
  else if( num == 5852) 
    return "PMM";
  else if( num == 5853) 
    return "PMF";
  else if( num == 5854) 
    return "PMP";
  else if( num == 5855) 
    return "PMS";
  else if( num == 5856) 
    return "PMT";
  else if( num == 5857) 
    return "PMW";
  else if( num == 5858) 
    return "PMY";
  else if( num == 5859) 
    return "PMV";
  else if( num == 5860) 
    return "PFA";
  else if( num == 5861) 
    return "PFR";
  else if( num == 5862) 
    return "PFN";
  else if( num == 5863) 
    return "PFD";
  else if( num == 5864) 
    return "PFC";
  else if( num == 5865) 
    return "PFQ";
  else if( num == 5866) 
    return "PFE";
  else if( num == 5867) 
    return "PFG";
  else if( num == 5868) 
    return "PFH";
  else if( num == 5869) 
    return "PFI";
  else if( num == 5870) 
    return "PFL";
  else if( num == 5871) 
    return "PFK";
  else if( num == 5872) 
    return "PFM";
  else if( num == 5873) 
    return "PFF";
  else if( num == 5874) 
    return "PFP";
  else if( num == 5875) 
    return "PFS";
  else if( num == 5876) 
    return "PFT";
  else if( num == 5877) 
    return "PFW";
  else if( num == 5878) 
    return "PFY";
  else if( num == 5879) 
    return "PFV";
  else if( num == 5880) 
    return "PPA";
  else if( num == 5881) 
    return "PPR";
  else if( num == 5882) 
    return "PPN";
  else if( num == 5883) 
    return "PPD";
  else if( num == 5884) 
    return "PPC";
  else if( num == 5885) 
    return "PPQ";
  else if( num == 5886) 
    return "PPE";
  else if( num == 5887) 
    return "PPG";
  else if( num == 5888) 
    return "PPH";
  else if( num == 5889) 
    return "PPI";
  else if( num == 5890) 
    return "PPL";
  else if( num == 5891) 
    return "PPK";
  else if( num == 5892) 
    return "PPM";
  else if( num == 5893) 
    return "PPF";
  else if( num == 5894) 
    return "PPP";
  else if( num == 5895) 
    return "PPS";
  else if( num == 5896) 
    return "PPT";
  else if( num == 5897) 
    return "PPW";
  else if( num == 5898) 
    return "PPY";
  else if( num == 5899) 
    return "PPV";
  else if( num == 5900) 
    return "PSA";
  else if( num == 5901) 
    return "PSR";
  else if( num == 5902) 
    return "PSN";
  else if( num == 5903) 
    return "PSD";
  else if( num == 5904) 
    return "PSC";
  else if( num == 5905) 
    return "PSQ";
  else if( num == 5906) 
    return "PSE";
  else if( num == 5907) 
    return "PSG";
  else if( num == 5908) 
    return "PSH";
  else if( num == 5909) 
    return "PSI";
  else if( num == 5910) 
    return "PSL";
  else if( num == 5911) 
    return "PSK";
  else if( num == 5912) 
    return "PSM";
  else if( num == 5913) 
    return "PSF";
  else if( num == 5914) 
    return "PSP";
  else if( num == 5915) 
    return "PSS";
  else if( num == 5916) 
    return "PST";
  else if( num == 5917) 
    return "PSW";
  else if( num == 5918) 
    return "PSY";
  else if( num == 5919) 
    return "PSV";
  else if( num == 5920) 
    return "PTA";
  else if( num == 5921) 
    return "PTR";
  else if( num == 5922) 
    return "PTN";
  else if( num == 5923) 
    return "PTD";
  else if( num == 5924) 
    return "PTC";
  else if( num == 5925) 
    return "PTQ";
  else if( num == 5926) 
    return "PTE";
  else if( num == 5927) 
    return "PTG";
  else if( num == 5928) 
    return "PTH";
  else if( num == 5929) 
    return "PTI";
  else if( num == 5930) 
    return "PTL";
  else if( num == 5931) 
    return "PTK";
  else if( num == 5932) 
    return "PTM";
  else if( num == 5933) 
    return "PTF";
  else if( num == 5934) 
    return "PTP";
  else if( num == 5935) 
    return "PTS";
  else if( num == 5936) 
    return "PTT";
  else if( num == 5937) 
    return "PTW";
  else if( num == 5938) 
    return "PTY";
  else if( num == 5939) 
    return "PTV";
  else if( num == 5940) 
    return "PWA";
  else if( num == 5941) 
    return "PWR";
  else if( num == 5942) 
    return "PWN";
  else if( num == 5943) 
    return "PWD";
  else if( num == 5944) 
    return "PWC";
  else if( num == 5945) 
    return "PWQ";
  else if( num == 5946) 
    return "PWE";
  else if( num == 5947) 
    return "PWG";
  else if( num == 5948) 
    return "PWH";
  else if( num == 5949) 
    return "PWI";
  else if( num == 5950) 
    return "PWL";
  else if( num == 5951) 
    return "PWK";
  else if( num == 5952) 
    return "PWM";
  else if( num == 5953) 
    return "PWF";
  else if( num == 5954) 
    return "PWP";
  else if( num == 5955) 
    return "PWS";
  else if( num == 5956) 
    return "PWT";
  else if( num == 5957) 
    return "PWW";
  else if( num == 5958) 
    return "PWY";
  else if( num == 5959) 
    return "PWV";
  else if( num == 5960) 
    return "PYA";
  else if( num == 5961) 
    return "PYR";
  else if( num == 5962) 
    return "PYN";
  else if( num == 5963) 
    return "PYD";
  else if( num == 5964) 
    return "PYC";
  else if( num == 5965) 
    return "PYQ";
  else if( num == 5966) 
    return "PYE";
  else if( num == 5967) 
    return "PYG";
  else if( num == 5968) 
    return "PYH";
  else if( num == 5969) 
    return "PYI";
  else if( num == 5970) 
    return "PYL";
  else if( num == 5971) 
    return "PYK";
  else if( num == 5972) 
    return "PYM";
  else if( num == 5973) 
    return "PYF";
  else if( num == 5974) 
    return "PYP";
  else if( num == 5975) 
    return "PYS";
  else if( num == 5976) 
    return "PYT";
  else if( num == 5977) 
    return "PYW";
  else if( num == 5978) 
    return "PYY";
  else if( num == 5979) 
    return "PYV";
  else if( num == 5980) 
    return "PVA";
  else if( num == 5981) 
    return "PVR";
  else if( num == 5982) 
    return "PVN";
  else if( num == 5983) 
    return "PVD";
  else if( num == 5984) 
    return "PVC";
  else if( num == 5985) 
    return "PVQ";
  else if( num == 5986) 
    return "PVE";
  else if( num == 5987) 
    return "PVG";
  else if( num == 5988) 
    return "PVH";
  else if( num == 5989) 
    return "PVI";
  else if( num == 5990) 
    return "PVL";
  else if( num == 5991) 
    return "PVK";
  else if( num == 5992) 
    return "PVM";
  else if( num == 5993) 
    return "PVF";
  else if( num == 5994) 
    return "PVP";
  else if( num == 5995) 
    return "PVS";
  else if( num == 5996) 
    return "PVT";
  else if( num == 5997) 
    return "PVW";
  else if( num == 5998) 
    return "PVY";
  else if( num == 5999) 
    return "PVV";
  else if( num == 6000) 
    return "SAA";
  else if( num == 6001) 
    return "SAR";
  else if( num == 6002) 
    return "SAN";
  else if( num == 6003) 
    return "SAD";
  else if( num == 6004) 
    return "SAC";
  else if( num == 6005) 
    return "SAQ";
  else if( num == 6006) 
    return "SAE";
  else if( num == 6007) 
    return "SAG";
  else if( num == 6008) 
    return "SAH";
  else if( num == 6009) 
    return "SAI";
  else if( num == 6010) 
    return "SAL";
  else if( num == 6011) 
    return "SAK";
  else if( num == 6012) 
    return "SAM";
  else if( num == 6013) 
    return "SAF";
  else if( num == 6014) 
    return "SAP";
  else if( num == 6015) 
    return "SAS";
  else if( num == 6016) 
    return "SAT";
  else if( num == 6017) 
    return "SAW";
  else if( num == 6018) 
    return "SAY";
  else if( num == 6019) 
    return "SAV";
  else if( num == 6020) 
    return "SRA";
  else if( num == 6021) 
    return "SRR";
  else if( num == 6022) 
    return "SRN";
  else if( num == 6023) 
    return "SRD";
  else if( num == 6024) 
    return "SRC";
  else if( num == 6025) 
    return "SRQ";
  else if( num == 6026) 
    return "SRE";
  else if( num == 6027) 
    return "SRG";
  else if( num == 6028) 
    return "SRH";
  else if( num == 6029) 
    return "SRI";
  else if( num == 6030) 
    return "SRL";
  else if( num == 6031) 
    return "SRK";
  else if( num == 6032) 
    return "SRM";
  else if( num == 6033) 
    return "SRF";
  else if( num == 6034) 
    return "SRP";
  else if( num == 6035) 
    return "SRS";
  else if( num == 6036) 
    return "SRT";
  else if( num == 6037) 
    return "SRW";
  else if( num == 6038) 
    return "SRY";
  else if( num == 6039) 
    return "SRV";
  else if( num == 6040) 
    return "SNA";
  else if( num == 6041) 
    return "SNR";
  else if( num == 6042) 
    return "SNN";
  else if( num == 6043) 
    return "SND";
  else if( num == 6044) 
    return "SNC";
  else if( num == 6045) 
    return "SNQ";
  else if( num == 6046) 
    return "SNE";
  else if( num == 6047) 
    return "SNG";
  else if( num == 6048) 
    return "SNH";
  else if( num == 6049) 
    return "SNI";
  else if( num == 6050) 
    return "SNL";
  else if( num == 6051) 
    return "SNK";
  else if( num == 6052) 
    return "SNM";
  else if( num == 6053) 
    return "SNF";
  else if( num == 6054) 
    return "SNP";
  else if( num == 6055) 
    return "SNS";
  else if( num == 6056) 
    return "SNT";
  else if( num == 6057) 
    return "SNW";
  else if( num == 6058) 
    return "SNY";
  else if( num == 6059) 
    return "SNV";
  else if( num == 6060) 
    return "SDA";
  else if( num == 6061) 
    return "SDR";
  else if( num == 6062) 
    return "SDN";
  else if( num == 6063) 
    return "SDD";
  else if( num == 6064) 
    return "SDC";
  else if( num == 6065) 
    return "SDQ";
  else if( num == 6066) 
    return "SDE";
  else if( num == 6067) 
    return "SDG";
  else if( num == 6068) 
    return "SDH";
  else if( num == 6069) 
    return "SDI";
  else if( num == 6070) 
    return "SDL";
  else if( num == 6071) 
    return "SDK";
  else if( num == 6072) 
    return "SDM";
  else if( num == 6073) 
    return "SDF";
  else if( num == 6074) 
    return "SDP";
  else if( num == 6075) 
    return "SDS";
  else if( num == 6076) 
    return "SDT";
  else if( num == 6077) 
    return "SDW";
  else if( num == 6078) 
    return "SDY";
  else if( num == 6079) 
    return "SDV";
  else if( num == 6080) 
    return "SCA";
  else if( num == 6081) 
    return "SCR";
  else if( num == 6082) 
    return "SCN";
  else if( num == 6083) 
    return "SCD";
  else if( num == 6084) 
    return "SCC";
  else if( num == 6085) 
    return "SCQ";
  else if( num == 6086) 
    return "SCE";
  else if( num == 6087) 
    return "SCG";
  else if( num == 6088) 
    return "SCH";
  else if( num == 6089) 
    return "SCI";
  else if( num == 6090) 
    return "SCL";
  else if( num == 6091) 
    return "SCK";
  else if( num == 6092) 
    return "SCM";
  else if( num == 6093) 
    return "SCF";
  else if( num == 6094) 
    return "SCP";
  else if( num == 6095) 
    return "SCS";
  else if( num == 6096) 
    return "SCT";
  else if( num == 6097) 
    return "SCW";
  else if( num == 6098) 
    return "SCY";
  else if( num == 6099) 
    return "SCV";
  else if( num == 6100) 
    return "SQA";
  else if( num == 6101) 
    return "SQR";
  else if( num == 6102) 
    return "SQN";
  else if( num == 6103) 
    return "SQD";
  else if( num == 6104) 
    return "SQC";
  else if( num == 6105) 
    return "SQQ";
  else if( num == 6106) 
    return "SQE";
  else if( num == 6107) 
    return "SQG";
  else if( num == 6108) 
    return "SQH";
  else if( num == 6109) 
    return "SQI";
  else if( num == 6110) 
    return "SQL";
  else if( num == 6111) 
    return "SQK";
  else if( num == 6112) 
    return "SQM";
  else if( num == 6113) 
    return "SQF";
  else if( num == 6114) 
    return "SQP";
  else if( num == 6115) 
    return "SQS";
  else if( num == 6116) 
    return "SQT";
  else if( num == 6117) 
    return "SQW";
  else if( num == 6118) 
    return "SQY";
  else if( num == 6119) 
    return "SQV";
  else if( num == 6120) 
    return "SEA";
  else if( num == 6121) 
    return "SER";
  else if( num == 6122) 
    return "SEN";
  else if( num == 6123) 
    return "SED";
  else if( num == 6124) 
    return "SEC";
  else if( num == 6125) 
    return "SEQ";
  else if( num == 6126) 
    return "SEE";
  else if( num == 6127) 
    return "SEG";
  else if( num == 6128) 
    return "SEH";
  else if( num == 6129) 
    return "SEI";
  else if( num == 6130) 
    return "SEL";
  else if( num == 6131) 
    return "SEK";
  else if( num == 6132) 
    return "SEM";
  else if( num == 6133) 
    return "SEF";
  else if( num == 6134) 
    return "SEP";
  else if( num == 6135) 
    return "SES";
  else if( num == 6136) 
    return "SET";
  else if( num == 6137) 
    return "SEW";
  else if( num == 6138) 
    return "SEY";
  else if( num == 6139) 
    return "SEV";
  else if( num == 6140) 
    return "SGA";
  else if( num == 6141) 
    return "SGR";
  else if( num == 6142) 
    return "SGN";
  else if( num == 6143) 
    return "SGD";
  else if( num == 6144) 
    return "SGC";
  else if( num == 6145) 
    return "SGQ";
  else if( num == 6146) 
    return "SGE";
  else if( num == 6147) 
    return "SGG";
  else if( num == 6148) 
    return "SGH";
  else if( num == 6149) 
    return "SGI";
  else if( num == 6150) 
    return "SGL";
  else if( num == 6151) 
    return "SGK";
  else if( num == 6152) 
    return "SGM";
  else if( num == 6153) 
    return "SGF";
  else if( num == 6154) 
    return "SGP";
  else if( num == 6155) 
    return "SGS";
  else if( num == 6156) 
    return "SGT";
  else if( num == 6157) 
    return "SGW";
  else if( num == 6158) 
    return "SGY";
  else if( num == 6159) 
    return "SGV";
  else if( num == 6160) 
    return "SHA";
  else if( num == 6161) 
    return "SHR";
  else if( num == 6162) 
    return "SHN";
  else if( num == 6163) 
    return "SHD";
  else if( num == 6164) 
    return "SHC";
  else if( num == 6165) 
    return "SHQ";
  else if( num == 6166) 
    return "SHE";
  else if( num == 6167) 
    return "SHG";
  else if( num == 6168) 
    return "SHH";
  else if( num == 6169) 
    return "SHI";
  else if( num == 6170) 
    return "SHL";
  else if( num == 6171) 
    return "SHK";
  else if( num == 6172) 
    return "SHM";
  else if( num == 6173) 
    return "SHF";
  else if( num == 6174) 
    return "SHP";
  else if( num == 6175) 
    return "SHS";
  else if( num == 6176) 
    return "SHT";
  else if( num == 6177) 
    return "SHW";
  else if( num == 6178) 
    return "SHY";
  else if( num == 6179) 
    return "SHV";
  else if( num == 6180) 
    return "SIA";
  else if( num == 6181) 
    return "SIR";
  else if( num == 6182) 
    return "SIN";
  else if( num == 6183) 
    return "SID";
  else if( num == 6184) 
    return "SIC";
  else if( num == 6185) 
    return "SIQ";
  else if( num == 6186) 
    return "SIE";
  else if( num == 6187) 
    return "SIG";
  else if( num == 6188) 
    return "SIH";
  else if( num == 6189) 
    return "SII";
  else if( num == 6190) 
    return "SIL";
  else if( num == 6191) 
    return "SIK";
  else if( num == 6192) 
    return "SIM";
  else if( num == 6193) 
    return "SIF";
  else if( num == 6194) 
    return "SIP";
  else if( num == 6195) 
    return "SIS";
  else if( num == 6196) 
    return "SIT";
  else if( num == 6197) 
    return "SIW";
  else if( num == 6198) 
    return "SIY";
  else if( num == 6199) 
    return "SIV";
  else if( num == 6200) 
    return "SLA";
  else if( num == 6201) 
    return "SLR";
  else if( num == 6202) 
    return "SLN";
  else if( num == 6203) 
    return "SLD";
  else if( num == 6204) 
    return "SLC";
  else if( num == 6205) 
    return "SLQ";
  else if( num == 6206) 
    return "SLE";
  else if( num == 6207) 
    return "SLG";
  else if( num == 6208) 
    return "SLH";
  else if( num == 6209) 
    return "SLI";
  else if( num == 6210) 
    return "SLL";
  else if( num == 6211) 
    return "SLK";
  else if( num == 6212) 
    return "SLM";
  else if( num == 6213) 
    return "SLF";
  else if( num == 6214) 
    return "SLP";
  else if( num == 6215) 
    return "SLS";
  else if( num == 6216) 
    return "SLT";
  else if( num == 6217) 
    return "SLW";
  else if( num == 6218) 
    return "SLY";
  else if( num == 6219) 
    return "SLV";
  else if( num == 6220) 
    return "SKA";
  else if( num == 6221) 
    return "SKR";
  else if( num == 6222) 
    return "SKN";
  else if( num == 6223) 
    return "SKD";
  else if( num == 6224) 
    return "SKC";
  else if( num == 6225) 
    return "SKQ";
  else if( num == 6226) 
    return "SKE";
  else if( num == 6227) 
    return "SKG";
  else if( num == 6228) 
    return "SKH";
  else if( num == 6229) 
    return "SKI";
  else if( num == 6230) 
    return "SKL";
  else if( num == 6231) 
    return "SKK";
  else if( num == 6232) 
    return "SKM";
  else if( num == 6233) 
    return "SKF";
  else if( num == 6234) 
    return "SKP";
  else if( num == 6235) 
    return "SKS";
  else if( num == 6236) 
    return "SKT";
  else if( num == 6237) 
    return "SKW";
  else if( num == 6238) 
    return "SKY";
  else if( num == 6239) 
    return "SKV";
  else if( num == 6240) 
    return "SMA";
  else if( num == 6241) 
    return "SMR";
  else if( num == 6242) 
    return "SMN";
  else if( num == 6243) 
    return "SMD";
  else if( num == 6244) 
    return "SMC";
  else if( num == 6245) 
    return "SMQ";
  else if( num == 6246) 
    return "SME";
  else if( num == 6247) 
    return "SMG";
  else if( num == 6248) 
    return "SMH";
  else if( num == 6249) 
    return "SMI";
  else if( num == 6250) 
    return "SML";
  else if( num == 6251) 
    return "SMK";
  else if( num == 6252) 
    return "SMM";
  else if( num == 6253) 
    return "SMF";
  else if( num == 6254) 
    return "SMP";
  else if( num == 6255) 
    return "SMS";
  else if( num == 6256) 
    return "SMT";
  else if( num == 6257) 
    return "SMW";
  else if( num == 6258) 
    return "SMY";
  else if( num == 6259) 
    return "SMV";
  else if( num == 6260) 
    return "SFA";
  else if( num == 6261) 
    return "SFR";
  else if( num == 6262) 
    return "SFN";
  else if( num == 6263) 
    return "SFD";
  else if( num == 6264) 
    return "SFC";
  else if( num == 6265) 
    return "SFQ";
  else if( num == 6266) 
    return "SFE";
  else if( num == 6267) 
    return "SFG";
  else if( num == 6268) 
    return "SFH";
  else if( num == 6269) 
    return "SFI";
  else if( num == 6270) 
    return "SFL";
  else if( num == 6271) 
    return "SFK";
  else if( num == 6272) 
    return "SFM";
  else if( num == 6273) 
    return "SFF";
  else if( num == 6274) 
    return "SFP";
  else if( num == 6275) 
    return "SFS";
  else if( num == 6276) 
    return "SFT";
  else if( num == 6277) 
    return "SFW";
  else if( num == 6278) 
    return "SFY";
  else if( num == 6279) 
    return "SFV";
  else if( num == 6280) 
    return "SPA";
  else if( num == 6281) 
    return "SPR";
  else if( num == 6282) 
    return "SPN";
  else if( num == 6283) 
    return "SPD";
  else if( num == 6284) 
    return "SPC";
  else if( num == 6285) 
    return "SPQ";
  else if( num == 6286) 
    return "SPE";
  else if( num == 6287) 
    return "SPG";
  else if( num == 6288) 
    return "SPH";
  else if( num == 6289) 
    return "SPI";
  else if( num == 6290) 
    return "SPL";
  else if( num == 6291) 
    return "SPK";
  else if( num == 6292) 
    return "SPM";
  else if( num == 6293) 
    return "SPF";
  else if( num == 6294) 
    return "SPP";
  else if( num == 6295) 
    return "SPS";
  else if( num == 6296) 
    return "SPT";
  else if( num == 6297) 
    return "SPW";
  else if( num == 6298) 
    return "SPY";
  else if( num == 6299) 
    return "SPV";
  else if( num == 6300) 
    return "SSA";
  else if( num == 6301) 
    return "SSR";
  else if( num == 6302) 
    return "SSN";
  else if( num == 6303) 
    return "SSD";
  else if( num == 6304) 
    return "SSC";
  else if( num == 6305) 
    return "SSQ";
  else if( num == 6306) 
    return "SSE";
  else if( num == 6307) 
    return "SSG";
  else if( num == 6308) 
    return "SSH";
  else if( num == 6309) 
    return "SSI";
  else if( num == 6310) 
    return "SSL";
  else if( num == 6311) 
    return "SSK";
  else if( num == 6312) 
    return "SSM";
  else if( num == 6313) 
    return "SSF";
  else if( num == 6314) 
    return "SSP";
  else if( num == 6315) 
    return "SSS";
  else if( num == 6316) 
    return "SST";
  else if( num == 6317) 
    return "SSW";
  else if( num == 6318) 
    return "SSY";
  else if( num == 6319) 
    return "SSV";
  else if( num == 6320) 
    return "STA";
  else if( num == 6321) 
    return "STR";
  else if( num == 6322) 
    return "STN";
  else if( num == 6323) 
    return "STD";
  else if( num == 6324) 
    return "STC";
  else if( num == 6325) 
    return "STQ";
  else if( num == 6326) 
    return "STE";
  else if( num == 6327) 
    return "STG";
  else if( num == 6328) 
    return "STH";
  else if( num == 6329) 
    return "STI";
  else if( num == 6330) 
    return "STL";
  else if( num == 6331) 
    return "STK";
  else if( num == 6332) 
    return "STM";
  else if( num == 6333) 
    return "STF";
  else if( num == 6334) 
    return "STP";
  else if( num == 6335) 
    return "STS";
  else if( num == 6336) 
    return "STT";
  else if( num == 6337) 
    return "STW";
  else if( num == 6338) 
    return "STY";
  else if( num == 6339) 
    return "STV";
  else if( num == 6340) 
    return "SWA";
  else if( num == 6341) 
    return "SWR";
  else if( num == 6342) 
    return "SWN";
  else if( num == 6343) 
    return "SWD";
  else if( num == 6344) 
    return "SWC";
  else if( num == 6345) 
    return "SWQ";
  else if( num == 6346) 
    return "SWE";
  else if( num == 6347) 
    return "SWG";
  else if( num == 6348) 
    return "SWH";
  else if( num == 6349) 
    return "SWI";
  else if( num == 6350) 
    return "SWL";
  else if( num == 6351) 
    return "SWK";
  else if( num == 6352) 
    return "SWM";
  else if( num == 6353) 
    return "SWF";
  else if( num == 6354) 
    return "SWP";
  else if( num == 6355) 
    return "SWS";
  else if( num == 6356) 
    return "SWT";
  else if( num == 6357) 
    return "SWW";
  else if( num == 6358) 
    return "SWY";
  else if( num == 6359) 
    return "SWV";
  else if( num == 6360) 
    return "SYA";
  else if( num == 6361) 
    return "SYR";
  else if( num == 6362) 
    return "SYN";
  else if( num == 6363) 
    return "SYD";
  else if( num == 6364) 
    return "SYC";
  else if( num == 6365) 
    return "SYQ";
  else if( num == 6366) 
    return "SYE";
  else if( num == 6367) 
    return "SYG";
  else if( num == 6368) 
    return "SYH";
  else if( num == 6369) 
    return "SYI";
  else if( num == 6370) 
    return "SYL";
  else if( num == 6371) 
    return "SYK";
  else if( num == 6372) 
    return "SYM";
  else if( num == 6373) 
    return "SYF";
  else if( num == 6374) 
    return "SYP";
  else if( num == 6375) 
    return "SYS";
  else if( num == 6376) 
    return "SYT";
  else if( num == 6377) 
    return "SYW";
  else if( num == 6378) 
    return "SYY";
  else if( num == 6379) 
    return "SYV";
  else if( num == 6380) 
    return "SVA";
  else if( num == 6381) 
    return "SVR";
  else if( num == 6382) 
    return "SVN";
  else if( num == 6383) 
    return "SVD";
  else if( num == 6384) 
    return "SVC";
  else if( num == 6385) 
    return "SVQ";
  else if( num == 6386) 
    return "SVE";
  else if( num == 6387) 
    return "SVG";
  else if( num == 6388) 
    return "SVH";
  else if( num == 6389) 
    return "SVI";
  else if( num == 6390) 
    return "SVL";
  else if( num == 6391) 
    return "SVK";
  else if( num == 6392) 
    return "SVM";
  else if( num == 6393) 
    return "SVF";
  else if( num == 6394) 
    return "SVP";
  else if( num == 6395) 
    return "SVS";
  else if( num == 6396) 
    return "SVT";
  else if( num == 6397) 
    return "SVW";
  else if( num == 6398) 
    return "SVY";
  else if( num == 6399) 
    return "SVV";
  else if( num == 6400) 
    return "TAA";
  else if( num == 6401) 
    return "TAR";
  else if( num == 6402) 
    return "TAN";
  else if( num == 6403) 
    return "TAD";
  else if( num == 6404) 
    return "TAC";
  else if( num == 6405) 
    return "TAQ";
  else if( num == 6406) 
    return "TAE";
  else if( num == 6407) 
    return "TAG";
  else if( num == 6408) 
    return "TAH";
  else if( num == 6409) 
    return "TAI";
  else if( num == 6410) 
    return "TAL";
  else if( num == 6411) 
    return "TAK";
  else if( num == 6412) 
    return "TAM";
  else if( num == 6413) 
    return "TAF";
  else if( num == 6414) 
    return "TAP";
  else if( num == 6415) 
    return "TAS";
  else if( num == 6416) 
    return "TAT";
  else if( num == 6417) 
    return "TAW";
  else if( num == 6418) 
    return "TAY";
  else if( num == 6419) 
    return "TAV";
  else if( num == 6420) 
    return "TRA";
  else if( num == 6421) 
    return "TRR";
  else if( num == 6422) 
    return "TRN";
  else if( num == 6423) 
    return "TRD";
  else if( num == 6424) 
    return "TRC";
  else if( num == 6425) 
    return "TRQ";
  else if( num == 6426) 
    return "TRE";
  else if( num == 6427) 
    return "TRG";
  else if( num == 6428) 
    return "TRH";
  else if( num == 6429) 
    return "TRI";
  else if( num == 6430) 
    return "TRL";
  else if( num == 6431) 
    return "TRK";
  else if( num == 6432) 
    return "TRM";
  else if( num == 6433) 
    return "TRF";
  else if( num == 6434) 
    return "TRP";
  else if( num == 6435) 
    return "TRS";
  else if( num == 6436) 
    return "TRT";
  else if( num == 6437) 
    return "TRW";
  else if( num == 6438) 
    return "TRY";
  else if( num == 6439) 
    return "TRV";
  else if( num == 6440) 
    return "TNA";
  else if( num == 6441) 
    return "TNR";
  else if( num == 6442) 
    return "TNN";
  else if( num == 6443) 
    return "TND";
  else if( num == 6444) 
    return "TNC";
  else if( num == 6445) 
    return "TNQ";
  else if( num == 6446) 
    return "TNE";
  else if( num == 6447) 
    return "TNG";
  else if( num == 6448) 
    return "TNH";
  else if( num == 6449) 
    return "TNI";
  else if( num == 6450) 
    return "TNL";
  else if( num == 6451) 
    return "TNK";
  else if( num == 6452) 
    return "TNM";
  else if( num == 6453) 
    return "TNF";
  else if( num == 6454) 
    return "TNP";
  else if( num == 6455) 
    return "TNS";
  else if( num == 6456) 
    return "TNT";
  else if( num == 6457) 
    return "TNW";
  else if( num == 6458) 
    return "TNY";
  else if( num == 6459) 
    return "TNV";
  else if( num == 6460) 
    return "TDA";
  else if( num == 6461) 
    return "TDR";
  else if( num == 6462) 
    return "TDN";
  else if( num == 6463) 
    return "TDD";
  else if( num == 6464) 
    return "TDC";
  else if( num == 6465) 
    return "TDQ";
  else if( num == 6466) 
    return "TDE";
  else if( num == 6467) 
    return "TDG";
  else if( num == 6468) 
    return "TDH";
  else if( num == 6469) 
    return "TDI";
  else if( num == 6470) 
    return "TDL";
  else if( num == 6471) 
    return "TDK";
  else if( num == 6472) 
    return "TDM";
  else if( num == 6473) 
    return "TDF";
  else if( num == 6474) 
    return "TDP";
  else if( num == 6475) 
    return "TDS";
  else if( num == 6476) 
    return "TDT";
  else if( num == 6477) 
    return "TDW";
  else if( num == 6478) 
    return "TDY";
  else if( num == 6479) 
    return "TDV";
  else if( num == 6480) 
    return "TCA";
  else if( num == 6481) 
    return "TCR";
  else if( num == 6482) 
    return "TCN";
  else if( num == 6483) 
    return "TCD";
  else if( num == 6484) 
    return "TCC";
  else if( num == 6485) 
    return "TCQ";
  else if( num == 6486) 
    return "TCE";
  else if( num == 6487) 
    return "TCG";
  else if( num == 6488) 
    return "TCH";
  else if( num == 6489) 
    return "TCI";
  else if( num == 6490) 
    return "TCL";
  else if( num == 6491) 
    return "TCK";
  else if( num == 6492) 
    return "TCM";
  else if( num == 6493) 
    return "TCF";
  else if( num == 6494) 
    return "TCP";
  else if( num == 6495) 
    return "TCS";
  else if( num == 6496) 
    return "TCT";
  else if( num == 6497) 
    return "TCW";
  else if( num == 6498) 
    return "TCY";
  else if( num == 6499) 
    return "TCV";
  else if( num == 6500) 
    return "TQA";
  else if( num == 6501) 
    return "TQR";
  else if( num == 6502) 
    return "TQN";
  else if( num == 6503) 
    return "TQD";
  else if( num == 6504) 
    return "TQC";
  else if( num == 6505) 
    return "TQQ";
  else if( num == 6506) 
    return "TQE";
  else if( num == 6507) 
    return "TQG";
  else if( num == 6508) 
    return "TQH";
  else if( num == 6509) 
    return "TQI";
  else if( num == 6510) 
    return "TQL";
  else if( num == 6511) 
    return "TQK";
  else if( num == 6512) 
    return "TQM";
  else if( num == 6513) 
    return "TQF";
  else if( num == 6514) 
    return "TQP";
  else if( num == 6515) 
    return "TQS";
  else if( num == 6516) 
    return "TQT";
  else if( num == 6517) 
    return "TQW";
  else if( num == 6518) 
    return "TQY";
  else if( num == 6519) 
    return "TQV";
  else if( num == 6520) 
    return "TEA";
  else if( num == 6521) 
    return "TER";
  else if( num == 6522) 
    return "TEN";
  else if( num == 6523) 
    return "TED";
  else if( num == 6524) 
    return "TEC";
  else if( num == 6525) 
    return "TEQ";
  else if( num == 6526) 
    return "TEE";
  else if( num == 6527) 
    return "TEG";
  else if( num == 6528) 
    return "TEH";
  else if( num == 6529) 
    return "TEI";
  else if( num == 6530) 
    return "TEL";
  else if( num == 6531) 
    return "TEK";
  else if( num == 6532) 
    return "TEM";
  else if( num == 6533) 
    return "TEF";
  else if( num == 6534) 
    return "TEP";
  else if( num == 6535) 
    return "TES";
  else if( num == 6536) 
    return "TET";
  else if( num == 6537) 
    return "TEW";
  else if( num == 6538) 
    return "TEY";
  else if( num == 6539) 
    return "TEV";
  else if( num == 6540) 
    return "TGA";
  else if( num == 6541) 
    return "TGR";
  else if( num == 6542) 
    return "TGN";
  else if( num == 6543) 
    return "TGD";
  else if( num == 6544) 
    return "TGC";
  else if( num == 6545) 
    return "TGQ";
  else if( num == 6546) 
    return "TGE";
  else if( num == 6547) 
    return "TGG";
  else if( num == 6548) 
    return "TGH";
  else if( num == 6549) 
    return "TGI";
  else if( num == 6550) 
    return "TGL";
  else if( num == 6551) 
    return "TGK";
  else if( num == 6552) 
    return "TGM";
  else if( num == 6553) 
    return "TGF";
  else if( num == 6554) 
    return "TGP";
  else if( num == 6555) 
    return "TGS";
  else if( num == 6556) 
    return "TGT";
  else if( num == 6557) 
    return "TGW";
  else if( num == 6558) 
    return "TGY";
  else if( num == 6559) 
    return "TGV";
  else if( num == 6560) 
    return "THA";
  else if( num == 6561) 
    return "THR";
  else if( num == 6562) 
    return "THN";
  else if( num == 6563) 
    return "THD";
  else if( num == 6564) 
    return "THC";
  else if( num == 6565) 
    return "THQ";
  else if( num == 6566) 
    return "THE";
  else if( num == 6567) 
    return "THG";
  else if( num == 6568) 
    return "THH";
  else if( num == 6569) 
    return "THI";
  else if( num == 6570) 
    return "THL";
  else if( num == 6571) 
    return "THK";
  else if( num == 6572) 
    return "THM";
  else if( num == 6573) 
    return "THF";
  else if( num == 6574) 
    return "THP";
  else if( num == 6575) 
    return "THS";
  else if( num == 6576) 
    return "THT";
  else if( num == 6577) 
    return "THW";
  else if( num == 6578) 
    return "THY";
  else if( num == 6579) 
    return "THV";
  else if( num == 6580) 
    return "TIA";
  else if( num == 6581) 
    return "TIR";
  else if( num == 6582) 
    return "TIN";
  else if( num == 6583) 
    return "TID";
  else if( num == 6584) 
    return "TIC";
  else if( num == 6585) 
    return "TIQ";
  else if( num == 6586) 
    return "TIE";
  else if( num == 6587) 
    return "TIG";
  else if( num == 6588) 
    return "TIH";
  else if( num == 6589) 
    return "TII";
  else if( num == 6590) 
    return "TIL";
  else if( num == 6591) 
    return "TIK";
  else if( num == 6592) 
    return "TIM";
  else if( num == 6593) 
    return "TIF";
  else if( num == 6594) 
    return "TIP";
  else if( num == 6595) 
    return "TIS";
  else if( num == 6596) 
    return "TIT";
  else if( num == 6597) 
    return "TIW";
  else if( num == 6598) 
    return "TIY";
  else if( num == 6599) 
    return "TIV";
  else if( num == 6600) 
    return "TLA";
  else if( num == 6601) 
    return "TLR";
  else if( num == 6602) 
    return "TLN";
  else if( num == 6603) 
    return "TLD";
  else if( num == 6604) 
    return "TLC";
  else if( num == 6605) 
    return "TLQ";
  else if( num == 6606) 
    return "TLE";
  else if( num == 6607) 
    return "TLG";
  else if( num == 6608) 
    return "TLH";
  else if( num == 6609) 
    return "TLI";
  else if( num == 6610) 
    return "TLL";
  else if( num == 6611) 
    return "TLK";
  else if( num == 6612) 
    return "TLM";
  else if( num == 6613) 
    return "TLF";
  else if( num == 6614) 
    return "TLP";
  else if( num == 6615) 
    return "TLS";
  else if( num == 6616) 
    return "TLT";
  else if( num == 6617) 
    return "TLW";
  else if( num == 6618) 
    return "TLY";
  else if( num == 6619) 
    return "TLV";
  else if( num == 6620) 
    return "TKA";
  else if( num == 6621) 
    return "TKR";
  else if( num == 6622) 
    return "TKN";
  else if( num == 6623) 
    return "TKD";
  else if( num == 6624) 
    return "TKC";
  else if( num == 6625) 
    return "TKQ";
  else if( num == 6626) 
    return "TKE";
  else if( num == 6627) 
    return "TKG";
  else if( num == 6628) 
    return "TKH";
  else if( num == 6629) 
    return "TKI";
  else if( num == 6630) 
    return "TKL";
  else if( num == 6631) 
    return "TKK";
  else if( num == 6632) 
    return "TKM";
  else if( num == 6633) 
    return "TKF";
  else if( num == 6634) 
    return "TKP";
  else if( num == 6635) 
    return "TKS";
  else if( num == 6636) 
    return "TKT";
  else if( num == 6637) 
    return "TKW";
  else if( num == 6638) 
    return "TKY";
  else if( num == 6639) 
    return "TKV";
  else if( num == 6640) 
    return "TMA";
  else if( num == 6641) 
    return "TMR";
  else if( num == 6642) 
    return "TMN";
  else if( num == 6643) 
    return "TMD";
  else if( num == 6644) 
    return "TMC";
  else if( num == 6645) 
    return "TMQ";
  else if( num == 6646) 
    return "TME";
  else if( num == 6647) 
    return "TMG";
  else if( num == 6648) 
    return "TMH";
  else if( num == 6649) 
    return "TMI";
  else if( num == 6650) 
    return "TML";
  else if( num == 6651) 
    return "TMK";
  else if( num == 6652) 
    return "TMM";
  else if( num == 6653) 
    return "TMF";
  else if( num == 6654) 
    return "TMP";
  else if( num == 6655) 
    return "TMS";
  else if( num == 6656) 
    return "TMT";
  else if( num == 6657) 
    return "TMW";
  else if( num == 6658) 
    return "TMY";
  else if( num == 6659) 
    return "TMV";
  else if( num == 6660) 
    return "TFA";
  else if( num == 6661) 
    return "TFR";
  else if( num == 6662) 
    return "TFN";
  else if( num == 6663) 
    return "TFD";
  else if( num == 6664) 
    return "TFC";
  else if( num == 6665) 
    return "TFQ";
  else if( num == 6666) 
    return "TFE";
  else if( num == 6667) 
    return "TFG";
  else if( num == 6668) 
    return "TFH";
  else if( num == 6669) 
    return "TFI";
  else if( num == 6670) 
    return "TFL";
  else if( num == 6671) 
    return "TFK";
  else if( num == 6672) 
    return "TFM";
  else if( num == 6673) 
    return "TFF";
  else if( num == 6674) 
    return "TFP";
  else if( num == 6675) 
    return "TFS";
  else if( num == 6676) 
    return "TFT";
  else if( num == 6677) 
    return "TFW";
  else if( num == 6678) 
    return "TFY";
  else if( num == 6679) 
    return "TFV";
  else if( num == 6680) 
    return "TPA";
  else if( num == 6681) 
    return "TPR";
  else if( num == 6682) 
    return "TPN";
  else if( num == 6683) 
    return "TPD";
  else if( num == 6684) 
    return "TPC";
  else if( num == 6685) 
    return "TPQ";
  else if( num == 6686) 
    return "TPE";
  else if( num == 6687) 
    return "TPG";
  else if( num == 6688) 
    return "TPH";
  else if( num == 6689) 
    return "TPI";
  else if( num == 6690) 
    return "TPL";
  else if( num == 6691) 
    return "TPK";
  else if( num == 6692) 
    return "TPM";
  else if( num == 6693) 
    return "TPF";
  else if( num == 6694) 
    return "TPP";
  else if( num == 6695) 
    return "TPS";
  else if( num == 6696) 
    return "TPT";
  else if( num == 6697) 
    return "TPW";
  else if( num == 6698) 
    return "TPY";
  else if( num == 6699) 
    return "TPV";
  else if( num == 6700) 
    return "TSA";
  else if( num == 6701) 
    return "TSR";
  else if( num == 6702) 
    return "TSN";
  else if( num == 6703) 
    return "TSD";
  else if( num == 6704) 
    return "TSC";
  else if( num == 6705) 
    return "TSQ";
  else if( num == 6706) 
    return "TSE";
  else if( num == 6707) 
    return "TSG";
  else if( num == 6708) 
    return "TSH";
  else if( num == 6709) 
    return "TSI";
  else if( num == 6710) 
    return "TSL";
  else if( num == 6711) 
    return "TSK";
  else if( num == 6712) 
    return "TSM";
  else if( num == 6713) 
    return "TSF";
  else if( num == 6714) 
    return "TSP";
  else if( num == 6715) 
    return "TSS";
  else if( num == 6716) 
    return "TST";
  else if( num == 6717) 
    return "TSW";
  else if( num == 6718) 
    return "TSY";
  else if( num == 6719) 
    return "TSV";
  else if( num == 6720) 
    return "TTA";
  else if( num == 6721) 
    return "TTR";
  else if( num == 6722) 
    return "TTN";
  else if( num == 6723) 
    return "TTD";
  else if( num == 6724) 
    return "TTC";
  else if( num == 6725) 
    return "TTQ";
  else if( num == 6726) 
    return "TTE";
  else if( num == 6727) 
    return "TTG";
  else if( num == 6728) 
    return "TTH";
  else if( num == 6729) 
    return "TTI";
  else if( num == 6730) 
    return "TTL";
  else if( num == 6731) 
    return "TTK";
  else if( num == 6732) 
    return "TTM";
  else if( num == 6733) 
    return "TTF";
  else if( num == 6734) 
    return "TTP";
  else if( num == 6735) 
    return "TTS";
  else if( num == 6736) 
    return "TTT";
  else if( num == 6737) 
    return "TTW";
  else if( num == 6738) 
    return "TTY";
  else if( num == 6739) 
    return "TTV";
  else if( num == 6740) 
    return "TWA";
  else if( num == 6741) 
    return "TWR";
  else if( num == 6742) 
    return "TWN";
  else if( num == 6743) 
    return "TWD";
  else if( num == 6744) 
    return "TWC";
  else if( num == 6745) 
    return "TWQ";
  else if( num == 6746) 
    return "TWE";
  else if( num == 6747) 
    return "TWG";
  else if( num == 6748) 
    return "TWH";
  else if( num == 6749) 
    return "TWI";
  else if( num == 6750) 
    return "TWL";
  else if( num == 6751) 
    return "TWK";
  else if( num == 6752) 
    return "TWM";
  else if( num == 6753) 
    return "TWF";
  else if( num == 6754) 
    return "TWP";
  else if( num == 6755) 
    return "TWS";
  else if( num == 6756) 
    return "TWT";
  else if( num == 6757) 
    return "TWW";
  else if( num == 6758) 
    return "TWY";
  else if( num == 6759) 
    return "TWV";
  else if( num == 6760) 
    return "TYA";
  else if( num == 6761) 
    return "TYR";
  else if( num == 6762) 
    return "TYN";
  else if( num == 6763) 
    return "TYD";
  else if( num == 6764) 
    return "TYC";
  else if( num == 6765) 
    return "TYQ";
  else if( num == 6766) 
    return "TYE";
  else if( num == 6767) 
    return "TYG";
  else if( num == 6768) 
    return "TYH";
  else if( num == 6769) 
    return "TYI";
  else if( num == 6770) 
    return "TYL";
  else if( num == 6771) 
    return "TYK";
  else if( num == 6772) 
    return "TYM";
  else if( num == 6773) 
    return "TYF";
  else if( num == 6774) 
    return "TYP";
  else if( num == 6775) 
    return "TYS";
  else if( num == 6776) 
    return "TYT";
  else if( num == 6777) 
    return "TYW";
  else if( num == 6778) 
    return "TYY";
  else if( num == 6779) 
    return "TYV";
  else if( num == 6780) 
    return "TVA";
  else if( num == 6781) 
    return "TVR";
  else if( num == 6782) 
    return "TVN";
  else if( num == 6783) 
    return "TVD";
  else if( num == 6784) 
    return "TVC";
  else if( num == 6785) 
    return "TVQ";
  else if( num == 6786) 
    return "TVE";
  else if( num == 6787) 
    return "TVG";
  else if( num == 6788) 
    return "TVH";
  else if( num == 6789) 
    return "TVI";
  else if( num == 6790) 
    return "TVL";
  else if( num == 6791) 
    return "TVK";
  else if( num == 6792) 
    return "TVM";
  else if( num == 6793) 
    return "TVF";
  else if( num == 6794) 
    return "TVP";
  else if( num == 6795) 
    return "TVS";
  else if( num == 6796) 
    return "TVT";
  else if( num == 6797) 
    return "TVW";
  else if( num == 6798) 
    return "TVY";
  else if( num == 6799) 
    return "TVV";
  else if( num == 6800) 
    return "WAA";
  else if( num == 6801) 
    return "WAR";
  else if( num == 6802) 
    return "WAN";
  else if( num == 6803) 
    return "WAD";
  else if( num == 6804) 
    return "WAC";
  else if( num == 6805) 
    return "WAQ";
  else if( num == 6806) 
    return "WAE";
  else if( num == 6807) 
    return "WAG";
  else if( num == 6808) 
    return "WAH";
  else if( num == 6809) 
    return "WAI";
  else if( num == 6810) 
    return "WAL";
  else if( num == 6811) 
    return "WAK";
  else if( num == 6812) 
    return "WAM";
  else if( num == 6813) 
    return "WAF";
  else if( num == 6814) 
    return "WAP";
  else if( num == 6815) 
    return "WAS";
  else if( num == 6816) 
    return "WAT";
  else if( num == 6817) 
    return "WAW";
  else if( num == 6818) 
    return "WAY";
  else if( num == 6819) 
    return "WAV";
  else if( num == 6820) 
    return "WRA";
  else if( num == 6821) 
    return "WRR";
  else if( num == 6822) 
    return "WRN";
  else if( num == 6823) 
    return "WRD";
  else if( num == 6824) 
    return "WRC";
  else if( num == 6825) 
    return "WRQ";
  else if( num == 6826) 
    return "WRE";
  else if( num == 6827) 
    return "WRG";
  else if( num == 6828) 
    return "WRH";
  else if( num == 6829) 
    return "WRI";
  else if( num == 6830) 
    return "WRL";
  else if( num == 6831) 
    return "WRK";
  else if( num == 6832) 
    return "WRM";
  else if( num == 6833) 
    return "WRF";
  else if( num == 6834) 
    return "WRP";
  else if( num == 6835) 
    return "WRS";
  else if( num == 6836) 
    return "WRT";
  else if( num == 6837) 
    return "WRW";
  else if( num == 6838) 
    return "WRY";
  else if( num == 6839) 
    return "WRV";
  else if( num == 6840) 
    return "WNA";
  else if( num == 6841) 
    return "WNR";
  else if( num == 6842) 
    return "WNN";
  else if( num == 6843) 
    return "WND";
  else if( num == 6844) 
    return "WNC";
  else if( num == 6845) 
    return "WNQ";
  else if( num == 6846) 
    return "WNE";
  else if( num == 6847) 
    return "WNG";
  else if( num == 6848) 
    return "WNH";
  else if( num == 6849) 
    return "WNI";
  else if( num == 6850) 
    return "WNL";
  else if( num == 6851) 
    return "WNK";
  else if( num == 6852) 
    return "WNM";
  else if( num == 6853) 
    return "WNF";
  else if( num == 6854) 
    return "WNP";
  else if( num == 6855) 
    return "WNS";
  else if( num == 6856) 
    return "WNT";
  else if( num == 6857) 
    return "WNW";
  else if( num == 6858) 
    return "WNY";
  else if( num == 6859) 
    return "WNV";
  else if( num == 6860) 
    return "WDA";
  else if( num == 6861) 
    return "WDR";
  else if( num == 6862) 
    return "WDN";
  else if( num == 6863) 
    return "WDD";
  else if( num == 6864) 
    return "WDC";
  else if( num == 6865) 
    return "WDQ";
  else if( num == 6866) 
    return "WDE";
  else if( num == 6867) 
    return "WDG";
  else if( num == 6868) 
    return "WDH";
  else if( num == 6869) 
    return "WDI";
  else if( num == 6870) 
    return "WDL";
  else if( num == 6871) 
    return "WDK";
  else if( num == 6872) 
    return "WDM";
  else if( num == 6873) 
    return "WDF";
  else if( num == 6874) 
    return "WDP";
  else if( num == 6875) 
    return "WDS";
  else if( num == 6876) 
    return "WDT";
  else if( num == 6877) 
    return "WDW";
  else if( num == 6878) 
    return "WDY";
  else if( num == 6879) 
    return "WDV";
  else if( num == 6880) 
    return "WCA";
  else if( num == 6881) 
    return "WCR";
  else if( num == 6882) 
    return "WCN";
  else if( num == 6883) 
    return "WCD";
  else if( num == 6884) 
    return "WCC";
  else if( num == 6885) 
    return "WCQ";
  else if( num == 6886) 
    return "WCE";
  else if( num == 6887) 
    return "WCG";
  else if( num == 6888) 
    return "WCH";
  else if( num == 6889) 
    return "WCI";
  else if( num == 6890) 
    return "WCL";
  else if( num == 6891) 
    return "WCK";
  else if( num == 6892) 
    return "WCM";
  else if( num == 6893) 
    return "WCF";
  else if( num == 6894) 
    return "WCP";
  else if( num == 6895) 
    return "WCS";
  else if( num == 6896) 
    return "WCT";
  else if( num == 6897) 
    return "WCW";
  else if( num == 6898) 
    return "WCY";
  else if( num == 6899) 
    return "WCV";
  else if( num == 6900) 
    return "WQA";
  else if( num == 6901) 
    return "WQR";
  else if( num == 6902) 
    return "WQN";
  else if( num == 6903) 
    return "WQD";
  else if( num == 6904) 
    return "WQC";
  else if( num == 6905) 
    return "WQQ";
  else if( num == 6906) 
    return "WQE";
  else if( num == 6907) 
    return "WQG";
  else if( num == 6908) 
    return "WQH";
  else if( num == 6909) 
    return "WQI";
  else if( num == 6910) 
    return "WQL";
  else if( num == 6911) 
    return "WQK";
  else if( num == 6912) 
    return "WQM";
  else if( num == 6913) 
    return "WQF";
  else if( num == 6914) 
    return "WQP";
  else if( num == 6915) 
    return "WQS";
  else if( num == 6916) 
    return "WQT";
  else if( num == 6917) 
    return "WQW";
  else if( num == 6918) 
    return "WQY";
  else if( num == 6919) 
    return "WQV";
  else if( num == 6920) 
    return "WEA";
  else if( num == 6921) 
    return "WER";
  else if( num == 6922) 
    return "WEN";
  else if( num == 6923) 
    return "WED";
  else if( num == 6924) 
    return "WEC";
  else if( num == 6925) 
    return "WEQ";
  else if( num == 6926) 
    return "WEE";
  else if( num == 6927) 
    return "WEG";
  else if( num == 6928) 
    return "WEH";
  else if( num == 6929) 
    return "WEI";
  else if( num == 6930) 
    return "WEL";
  else if( num == 6931) 
    return "WEK";
  else if( num == 6932) 
    return "WEM";
  else if( num == 6933) 
    return "WEF";
  else if( num == 6934) 
    return "WEP";
  else if( num == 6935) 
    return "WES";
  else if( num == 6936) 
    return "WET";
  else if( num == 6937) 
    return "WEW";
  else if( num == 6938) 
    return "WEY";
  else if( num == 6939) 
    return "WEV";
  else if( num == 6940) 
    return "WGA";
  else if( num == 6941) 
    return "WGR";
  else if( num == 6942) 
    return "WGN";
  else if( num == 6943) 
    return "WGD";
  else if( num == 6944) 
    return "WGC";
  else if( num == 6945) 
    return "WGQ";
  else if( num == 6946) 
    return "WGE";
  else if( num == 6947) 
    return "WGG";
  else if( num == 6948) 
    return "WGH";
  else if( num == 6949) 
    return "WGI";
  else if( num == 6950) 
    return "WGL";
  else if( num == 6951) 
    return "WGK";
  else if( num == 6952) 
    return "WGM";
  else if( num == 6953) 
    return "WGF";
  else if( num == 6954) 
    return "WGP";
  else if( num == 6955) 
    return "WGS";
  else if( num == 6956) 
    return "WGT";
  else if( num == 6957) 
    return "WGW";
  else if( num == 6958) 
    return "WGY";
  else if( num == 6959) 
    return "WGV";
  else if( num == 6960) 
    return "WHA";
  else if( num == 6961) 
    return "WHR";
  else if( num == 6962) 
    return "WHN";
  else if( num == 6963) 
    return "WHD";
  else if( num == 6964) 
    return "WHC";
  else if( num == 6965) 
    return "WHQ";
  else if( num == 6966) 
    return "WHE";
  else if( num == 6967) 
    return "WHG";
  else if( num == 6968) 
    return "WHH";
  else if( num == 6969) 
    return "WHI";
  else if( num == 6970) 
    return "WHL";
  else if( num == 6971) 
    return "WHK";
  else if( num == 6972) 
    return "WHM";
  else if( num == 6973) 
    return "WHF";
  else if( num == 6974) 
    return "WHP";
  else if( num == 6975) 
    return "WHS";
  else if( num == 6976) 
    return "WHT";
  else if( num == 6977) 
    return "WHW";
  else if( num == 6978) 
    return "WHY";
  else if( num == 6979) 
    return "WHV";
  else if( num == 6980) 
    return "WIA";
  else if( num == 6981) 
    return "WIR";
  else if( num == 6982) 
    return "WIN";
  else if( num == 6983) 
    return "WID";
  else if( num == 6984) 
    return "WIC";
  else if( num == 6985) 
    return "WIQ";
  else if( num == 6986) 
    return "WIE";
  else if( num == 6987) 
    return "WIG";
  else if( num == 6988) 
    return "WIH";
  else if( num == 6989) 
    return "WII";
  else if( num == 6990) 
    return "WIL";
  else if( num == 6991) 
    return "WIK";
  else if( num == 6992) 
    return "WIM";
  else if( num == 6993) 
    return "WIF";
  else if( num == 6994) 
    return "WIP";
  else if( num == 6995) 
    return "WIS";
  else if( num == 6996) 
    return "WIT";
  else if( num == 6997) 
    return "WIW";
  else if( num == 6998) 
    return "WIY";
  else if( num == 6999) 
    return "WIV";
  else if( num == 7000) 
    return "WLA";
  else if( num == 7001) 
    return "WLR";
  else if( num == 7002) 
    return "WLN";
  else if( num == 7003) 
    return "WLD";
  else if( num == 7004) 
    return "WLC";
  else if( num == 7005) 
    return "WLQ";
  else if( num == 7006) 
    return "WLE";
  else if( num == 7007) 
    return "WLG";
  else if( num == 7008) 
    return "WLH";
  else if( num == 7009) 
    return "WLI";
  else if( num == 7010) 
    return "WLL";
  else if( num == 7011) 
    return "WLK";
  else if( num == 7012) 
    return "WLM";
  else if( num == 7013) 
    return "WLF";
  else if( num == 7014) 
    return "WLP";
  else if( num == 7015) 
    return "WLS";
  else if( num == 7016) 
    return "WLT";
  else if( num == 7017) 
    return "WLW";
  else if( num == 7018) 
    return "WLY";
  else if( num == 7019) 
    return "WLV";
  else if( num == 7020) 
    return "WKA";
  else if( num == 7021) 
    return "WKR";
  else if( num == 7022) 
    return "WKN";
  else if( num == 7023) 
    return "WKD";
  else if( num == 7024) 
    return "WKC";
  else if( num == 7025) 
    return "WKQ";
  else if( num == 7026) 
    return "WKE";
  else if( num == 7027) 
    return "WKG";
  else if( num == 7028) 
    return "WKH";
  else if( num == 7029) 
    return "WKI";
  else if( num == 7030) 
    return "WKL";
  else if( num == 7031) 
    return "WKK";
  else if( num == 7032) 
    return "WKM";
  else if( num == 7033) 
    return "WKF";
  else if( num == 7034) 
    return "WKP";
  else if( num == 7035) 
    return "WKS";
  else if( num == 7036) 
    return "WKT";
  else if( num == 7037) 
    return "WKW";
  else if( num == 7038) 
    return "WKY";
  else if( num == 7039) 
    return "WKV";
  else if( num == 7040) 
    return "WMA";
  else if( num == 7041) 
    return "WMR";
  else if( num == 7042) 
    return "WMN";
  else if( num == 7043) 
    return "WMD";
  else if( num == 7044) 
    return "WMC";
  else if( num == 7045) 
    return "WMQ";
  else if( num == 7046) 
    return "WME";
  else if( num == 7047) 
    return "WMG";
  else if( num == 7048) 
    return "WMH";
  else if( num == 7049) 
    return "WMI";
  else if( num == 7050) 
    return "WML";
  else if( num == 7051) 
    return "WMK";
  else if( num == 7052) 
    return "WMM";
  else if( num == 7053) 
    return "WMF";
  else if( num == 7054) 
    return "WMP";
  else if( num == 7055) 
    return "WMS";
  else if( num == 7056) 
    return "WMT";
  else if( num == 7057) 
    return "WMW";
  else if( num == 7058) 
    return "WMY";
  else if( num == 7059) 
    return "WMV";
  else if( num == 7060) 
    return "WFA";
  else if( num == 7061) 
    return "WFR";
  else if( num == 7062) 
    return "WFN";
  else if( num == 7063) 
    return "WFD";
  else if( num == 7064) 
    return "WFC";
  else if( num == 7065) 
    return "WFQ";
  else if( num == 7066) 
    return "WFE";
  else if( num == 7067) 
    return "WFG";
  else if( num == 7068) 
    return "WFH";
  else if( num == 7069) 
    return "WFI";
  else if( num == 7070) 
    return "WFL";
  else if( num == 7071) 
    return "WFK";
  else if( num == 7072) 
    return "WFM";
  else if( num == 7073) 
    return "WFF";
  else if( num == 7074) 
    return "WFP";
  else if( num == 7075) 
    return "WFS";
  else if( num == 7076) 
    return "WFT";
  else if( num == 7077) 
    return "WFW";
  else if( num == 7078) 
    return "WFY";
  else if( num == 7079) 
    return "WFV";
  else if( num == 7080) 
    return "WPA";
  else if( num == 7081) 
    return "WPR";
  else if( num == 7082) 
    return "WPN";
  else if( num == 7083) 
    return "WPD";
  else if( num == 7084) 
    return "WPC";
  else if( num == 7085) 
    return "WPQ";
  else if( num == 7086) 
    return "WPE";
  else if( num == 7087) 
    return "WPG";
  else if( num == 7088) 
    return "WPH";
  else if( num == 7089) 
    return "WPI";
  else if( num == 7090) 
    return "WPL";
  else if( num == 7091) 
    return "WPK";
  else if( num == 7092) 
    return "WPM";
  else if( num == 7093) 
    return "WPF";
  else if( num == 7094) 
    return "WPP";
  else if( num == 7095) 
    return "WPS";
  else if( num == 7096) 
    return "WPT";
  else if( num == 7097) 
    return "WPW";
  else if( num == 7098) 
    return "WPY";
  else if( num == 7099) 
    return "WPV";
  else if( num == 7100) 
    return "WSA";
  else if( num == 7101) 
    return "WSR";
  else if( num == 7102) 
    return "WSN";
  else if( num == 7103) 
    return "WSD";
  else if( num == 7104) 
    return "WSC";
  else if( num == 7105) 
    return "WSQ";
  else if( num == 7106) 
    return "WSE";
  else if( num == 7107) 
    return "WSG";
  else if( num == 7108) 
    return "WSH";
  else if( num == 7109) 
    return "WSI";
  else if( num == 7110) 
    return "WSL";
  else if( num == 7111) 
    return "WSK";
  else if( num == 7112) 
    return "WSM";
  else if( num == 7113) 
    return "WSF";
  else if( num == 7114) 
    return "WSP";
  else if( num == 7115) 
    return "WSS";
  else if( num == 7116) 
    return "WST";
  else if( num == 7117) 
    return "WSW";
  else if( num == 7118) 
    return "WSY";
  else if( num == 7119) 
    return "WSV";
  else if( num == 7120) 
    return "WTA";
  else if( num == 7121) 
    return "WTR";
  else if( num == 7122) 
    return "WTN";
  else if( num == 7123) 
    return "WTD";
  else if( num == 7124) 
    return "WTC";
  else if( num == 7125) 
    return "WTQ";
  else if( num == 7126) 
    return "WTE";
  else if( num == 7127) 
    return "WTG";
  else if( num == 7128) 
    return "WTH";
  else if( num == 7129) 
    return "WTI";
  else if( num == 7130) 
    return "WTL";
  else if( num == 7131) 
    return "WTK";
  else if( num == 7132) 
    return "WTM";
  else if( num == 7133) 
    return "WTF";
  else if( num == 7134) 
    return "WTP";
  else if( num == 7135) 
    return "WTS";
  else if( num == 7136) 
    return "WTT";
  else if( num == 7137) 
    return "WTW";
  else if( num == 7138) 
    return "WTY";
  else if( num == 7139) 
    return "WTV";
  else if( num == 7140) 
    return "WWA";
  else if( num == 7141) 
    return "WWR";
  else if( num == 7142) 
    return "WWN";
  else if( num == 7143) 
    return "WWD";
  else if( num == 7144) 
    return "WWC";
  else if( num == 7145) 
    return "WWQ";
  else if( num == 7146) 
    return "WWE";
  else if( num == 7147) 
    return "WWG";
  else if( num == 7148) 
    return "WWH";
  else if( num == 7149) 
    return "WWI";
  else if( num == 7150) 
    return "WWL";
  else if( num == 7151) 
    return "WWK";
  else if( num == 7152) 
    return "WWM";
  else if( num == 7153) 
    return "WWF";
  else if( num == 7154) 
    return "WWP";
  else if( num == 7155) 
    return "WWS";
  else if( num == 7156) 
    return "WWT";
  else if( num == 7157) 
    return "WWW";
  else if( num == 7158) 
    return "WWY";
  else if( num == 7159) 
    return "WWV";
  else if( num == 7160) 
    return "WYA";
  else if( num == 7161) 
    return "WYR";
  else if( num == 7162) 
    return "WYN";
  else if( num == 7163) 
    return "WYD";
  else if( num == 7164) 
    return "WYC";
  else if( num == 7165) 
    return "WYQ";
  else if( num == 7166) 
    return "WYE";
  else if( num == 7167) 
    return "WYG";
  else if( num == 7168) 
    return "WYH";
  else if( num == 7169) 
    return "WYI";
  else if( num == 7170) 
    return "WYL";
  else if( num == 7171) 
    return "WYK";
  else if( num == 7172) 
    return "WYM";
  else if( num == 7173) 
    return "WYF";
  else if( num == 7174) 
    return "WYP";
  else if( num == 7175) 
    return "WYS";
  else if( num == 7176) 
    return "WYT";
  else if( num == 7177) 
    return "WYW";
  else if( num == 7178) 
    return "WYY";
  else if( num == 7179) 
    return "WYV";
  else if( num == 7180) 
    return "WVA";
  else if( num == 7181) 
    return "WVR";
  else if( num == 7182) 
    return "WVN";
  else if( num == 7183) 
    return "WVD";
  else if( num == 7184) 
    return "WVC";
  else if( num == 7185) 
    return "WVQ";
  else if( num == 7186) 
    return "WVE";
  else if( num == 7187) 
    return "WVG";
  else if( num == 7188) 
    return "WVH";
  else if( num == 7189) 
    return "WVI";
  else if( num == 7190) 
    return "WVL";
  else if( num == 7191) 
    return "WVK";
  else if( num == 7192) 
    return "WVM";
  else if( num == 7193) 
    return "WVF";
  else if( num == 7194) 
    return "WVP";
  else if( num == 7195) 
    return "WVS";
  else if( num == 7196) 
    return "WVT";
  else if( num == 7197) 
    return "WVW";
  else if( num == 7198) 
    return "WVY";
  else if( num == 7199) 
    return "WVV";
  else if( num == 7200) 
    return "YAA";
  else if( num == 7201) 
    return "YAR";
  else if( num == 7202) 
    return "YAN";
  else if( num == 7203) 
    return "YAD";
  else if( num == 7204) 
    return "YAC";
  else if( num == 7205) 
    return "YAQ";
  else if( num == 7206) 
    return "YAE";
  else if( num == 7207) 
    return "YAG";
  else if( num == 7208) 
    return "YAH";
  else if( num == 7209) 
    return "YAI";
  else if( num == 7210) 
    return "YAL";
  else if( num == 7211) 
    return "YAK";
  else if( num == 7212) 
    return "YAM";
  else if( num == 7213) 
    return "YAF";
  else if( num == 7214) 
    return "YAP";
  else if( num == 7215) 
    return "YAS";
  else if( num == 7216) 
    return "YAT";
  else if( num == 7217) 
    return "YAW";
  else if( num == 7218) 
    return "YAY";
  else if( num == 7219) 
    return "YAV";
  else if( num == 7220) 
    return "YRA";
  else if( num == 7221) 
    return "YRR";
  else if( num == 7222) 
    return "YRN";
  else if( num == 7223) 
    return "YRD";
  else if( num == 7224) 
    return "YRC";
  else if( num == 7225) 
    return "YRQ";
  else if( num == 7226) 
    return "YRE";
  else if( num == 7227) 
    return "YRG";
  else if( num == 7228) 
    return "YRH";
  else if( num == 7229) 
    return "YRI";
  else if( num == 7230) 
    return "YRL";
  else if( num == 7231) 
    return "YRK";
  else if( num == 7232) 
    return "YRM";
  else if( num == 7233) 
    return "YRF";
  else if( num == 7234) 
    return "YRP";
  else if( num == 7235) 
    return "YRS";
  else if( num == 7236) 
    return "YRT";
  else if( num == 7237) 
    return "YRW";
  else if( num == 7238) 
    return "YRY";
  else if( num == 7239) 
    return "YRV";
  else if( num == 7240) 
    return "YNA";
  else if( num == 7241) 
    return "YNR";
  else if( num == 7242) 
    return "YNN";
  else if( num == 7243) 
    return "YND";
  else if( num == 7244) 
    return "YNC";
  else if( num == 7245) 
    return "YNQ";
  else if( num == 7246) 
    return "YNE";
  else if( num == 7247) 
    return "YNG";
  else if( num == 7248) 
    return "YNH";
  else if( num == 7249) 
    return "YNI";
  else if( num == 7250) 
    return "YNL";
  else if( num == 7251) 
    return "YNK";
  else if( num == 7252) 
    return "YNM";
  else if( num == 7253) 
    return "YNF";
  else if( num == 7254) 
    return "YNP";
  else if( num == 7255) 
    return "YNS";
  else if( num == 7256) 
    return "YNT";
  else if( num == 7257) 
    return "YNW";
  else if( num == 7258) 
    return "YNY";
  else if( num == 7259) 
    return "YNV";
  else if( num == 7260) 
    return "YDA";
  else if( num == 7261) 
    return "YDR";
  else if( num == 7262) 
    return "YDN";
  else if( num == 7263) 
    return "YDD";
  else if( num == 7264) 
    return "YDC";
  else if( num == 7265) 
    return "YDQ";
  else if( num == 7266) 
    return "YDE";
  else if( num == 7267) 
    return "YDG";
  else if( num == 7268) 
    return "YDH";
  else if( num == 7269) 
    return "YDI";
  else if( num == 7270) 
    return "YDL";
  else if( num == 7271) 
    return "YDK";
  else if( num == 7272) 
    return "YDM";
  else if( num == 7273) 
    return "YDF";
  else if( num == 7274) 
    return "YDP";
  else if( num == 7275) 
    return "YDS";
  else if( num == 7276) 
    return "YDT";
  else if( num == 7277) 
    return "YDW";
  else if( num == 7278) 
    return "YDY";
  else if( num == 7279) 
    return "YDV";
  else if( num == 7280) 
    return "YCA";
  else if( num == 7281) 
    return "YCR";
  else if( num == 7282) 
    return "YCN";
  else if( num == 7283) 
    return "YCD";
  else if( num == 7284) 
    return "YCC";
  else if( num == 7285) 
    return "YCQ";
  else if( num == 7286) 
    return "YCE";
  else if( num == 7287) 
    return "YCG";
  else if( num == 7288) 
    return "YCH";
  else if( num == 7289) 
    return "YCI";
  else if( num == 7290) 
    return "YCL";
  else if( num == 7291) 
    return "YCK";
  else if( num == 7292) 
    return "YCM";
  else if( num == 7293) 
    return "YCF";
  else if( num == 7294) 
    return "YCP";
  else if( num == 7295) 
    return "YCS";
  else if( num == 7296) 
    return "YCT";
  else if( num == 7297) 
    return "YCW";
  else if( num == 7298) 
    return "YCY";
  else if( num == 7299) 
    return "YCV";
  else if( num == 7300) 
    return "YQA";
  else if( num == 7301) 
    return "YQR";
  else if( num == 7302) 
    return "YQN";
  else if( num == 7303) 
    return "YQD";
  else if( num == 7304) 
    return "YQC";
  else if( num == 7305) 
    return "YQQ";
  else if( num == 7306) 
    return "YQE";
  else if( num == 7307) 
    return "YQG";
  else if( num == 7308) 
    return "YQH";
  else if( num == 7309) 
    return "YQI";
  else if( num == 7310) 
    return "YQL";
  else if( num == 7311) 
    return "YQK";
  else if( num == 7312) 
    return "YQM";
  else if( num == 7313) 
    return "YQF";
  else if( num == 7314) 
    return "YQP";
  else if( num == 7315) 
    return "YQS";
  else if( num == 7316) 
    return "YQT";
  else if( num == 7317) 
    return "YQW";
  else if( num == 7318) 
    return "YQY";
  else if( num == 7319) 
    return "YQV";
  else if( num == 7320) 
    return "YEA";
  else if( num == 7321) 
    return "YER";
  else if( num == 7322) 
    return "YEN";
  else if( num == 7323) 
    return "YED";
  else if( num == 7324) 
    return "YEC";
  else if( num == 7325) 
    return "YEQ";
  else if( num == 7326) 
    return "YEE";
  else if( num == 7327) 
    return "YEG";
  else if( num == 7328) 
    return "YEH";
  else if( num == 7329) 
    return "YEI";
  else if( num == 7330) 
    return "YEL";
  else if( num == 7331) 
    return "YEK";
  else if( num == 7332) 
    return "YEM";
  else if( num == 7333) 
    return "YEF";
  else if( num == 7334) 
    return "YEP";
  else if( num == 7335) 
    return "YES";
  else if( num == 7336) 
    return "YET";
  else if( num == 7337) 
    return "YEW";
  else if( num == 7338) 
    return "YEY";
  else if( num == 7339) 
    return "YEV";
  else if( num == 7340) 
    return "YGA";
  else if( num == 7341) 
    return "YGR";
  else if( num == 7342) 
    return "YGN";
  else if( num == 7343) 
    return "YGD";
  else if( num == 7344) 
    return "YGC";
  else if( num == 7345) 
    return "YGQ";
  else if( num == 7346) 
    return "YGE";
  else if( num == 7347) 
    return "YGG";
  else if( num == 7348) 
    return "YGH";
  else if( num == 7349) 
    return "YGI";
  else if( num == 7350) 
    return "YGL";
  else if( num == 7351) 
    return "YGK";
  else if( num == 7352) 
    return "YGM";
  else if( num == 7353) 
    return "YGF";
  else if( num == 7354) 
    return "YGP";
  else if( num == 7355) 
    return "YGS";
  else if( num == 7356) 
    return "YGT";
  else if( num == 7357) 
    return "YGW";
  else if( num == 7358) 
    return "YGY";
  else if( num == 7359) 
    return "YGV";
  else if( num == 7360) 
    return "YHA";
  else if( num == 7361) 
    return "YHR";
  else if( num == 7362) 
    return "YHN";
  else if( num == 7363) 
    return "YHD";
  else if( num == 7364) 
    return "YHC";
  else if( num == 7365) 
    return "YHQ";
  else if( num == 7366) 
    return "YHE";
  else if( num == 7367) 
    return "YHG";
  else if( num == 7368) 
    return "YHH";
  else if( num == 7369) 
    return "YHI";
  else if( num == 7370) 
    return "YHL";
  else if( num == 7371) 
    return "YHK";
  else if( num == 7372) 
    return "YHM";
  else if( num == 7373) 
    return "YHF";
  else if( num == 7374) 
    return "YHP";
  else if( num == 7375) 
    return "YHS";
  else if( num == 7376) 
    return "YHT";
  else if( num == 7377) 
    return "YHW";
  else if( num == 7378) 
    return "YHY";
  else if( num == 7379) 
    return "YHV";
  else if( num == 7380) 
    return "YIA";
  else if( num == 7381) 
    return "YIR";
  else if( num == 7382) 
    return "YIN";
  else if( num == 7383) 
    return "YID";
  else if( num == 7384) 
    return "YIC";
  else if( num == 7385) 
    return "YIQ";
  else if( num == 7386) 
    return "YIE";
  else if( num == 7387) 
    return "YIG";
  else if( num == 7388) 
    return "YIH";
  else if( num == 7389) 
    return "YII";
  else if( num == 7390) 
    return "YIL";
  else if( num == 7391) 
    return "YIK";
  else if( num == 7392) 
    return "YIM";
  else if( num == 7393) 
    return "YIF";
  else if( num == 7394) 
    return "YIP";
  else if( num == 7395) 
    return "YIS";
  else if( num == 7396) 
    return "YIT";
  else if( num == 7397) 
    return "YIW";
  else if( num == 7398) 
    return "YIY";
  else if( num == 7399) 
    return "YIV";
  else if( num == 7400) 
    return "YLA";
  else if( num == 7401) 
    return "YLR";
  else if( num == 7402) 
    return "YLN";
  else if( num == 7403) 
    return "YLD";
  else if( num == 7404) 
    return "YLC";
  else if( num == 7405) 
    return "YLQ";
  else if( num == 7406) 
    return "YLE";
  else if( num == 7407) 
    return "YLG";
  else if( num == 7408) 
    return "YLH";
  else if( num == 7409) 
    return "YLI";
  else if( num == 7410) 
    return "YLL";
  else if( num == 7411) 
    return "YLK";
  else if( num == 7412) 
    return "YLM";
  else if( num == 7413) 
    return "YLF";
  else if( num == 7414) 
    return "YLP";
  else if( num == 7415) 
    return "YLS";
  else if( num == 7416) 
    return "YLT";
  else if( num == 7417) 
    return "YLW";
  else if( num == 7418) 
    return "YLY";
  else if( num == 7419) 
    return "YLV";
  else if( num == 7420) 
    return "YKA";
  else if( num == 7421) 
    return "YKR";
  else if( num == 7422) 
    return "YKN";
  else if( num == 7423) 
    return "YKD";
  else if( num == 7424) 
    return "YKC";
  else if( num == 7425) 
    return "YKQ";
  else if( num == 7426) 
    return "YKE";
  else if( num == 7427) 
    return "YKG";
  else if( num == 7428) 
    return "YKH";
  else if( num == 7429) 
    return "YKI";
  else if( num == 7430) 
    return "YKL";
  else if( num == 7431) 
    return "YKK";
  else if( num == 7432) 
    return "YKM";
  else if( num == 7433) 
    return "YKF";
  else if( num == 7434) 
    return "YKP";
  else if( num == 7435) 
    return "YKS";
  else if( num == 7436) 
    return "YKT";
  else if( num == 7437) 
    return "YKW";
  else if( num == 7438) 
    return "YKY";
  else if( num == 7439) 
    return "YKV";
  else if( num == 7440) 
    return "YMA";
  else if( num == 7441) 
    return "YMR";
  else if( num == 7442) 
    return "YMN";
  else if( num == 7443) 
    return "YMD";
  else if( num == 7444) 
    return "YMC";
  else if( num == 7445) 
    return "YMQ";
  else if( num == 7446) 
    return "YME";
  else if( num == 7447) 
    return "YMG";
  else if( num == 7448) 
    return "YMH";
  else if( num == 7449) 
    return "YMI";
  else if( num == 7450) 
    return "YML";
  else if( num == 7451) 
    return "YMK";
  else if( num == 7452) 
    return "YMM";
  else if( num == 7453) 
    return "YMF";
  else if( num == 7454) 
    return "YMP";
  else if( num == 7455) 
    return "YMS";
  else if( num == 7456) 
    return "YMT";
  else if( num == 7457) 
    return "YMW";
  else if( num == 7458) 
    return "YMY";
  else if( num == 7459) 
    return "YMV";
  else if( num == 7460) 
    return "YFA";
  else if( num == 7461) 
    return "YFR";
  else if( num == 7462) 
    return "YFN";
  else if( num == 7463) 
    return "YFD";
  else if( num == 7464) 
    return "YFC";
  else if( num == 7465) 
    return "YFQ";
  else if( num == 7466) 
    return "YFE";
  else if( num == 7467) 
    return "YFG";
  else if( num == 7468) 
    return "YFH";
  else if( num == 7469) 
    return "YFI";
  else if( num == 7470) 
    return "YFL";
  else if( num == 7471) 
    return "YFK";
  else if( num == 7472) 
    return "YFM";
  else if( num == 7473) 
    return "YFF";
  else if( num == 7474) 
    return "YFP";
  else if( num == 7475) 
    return "YFS";
  else if( num == 7476) 
    return "YFT";
  else if( num == 7477) 
    return "YFW";
  else if( num == 7478) 
    return "YFY";
  else if( num == 7479) 
    return "YFV";
  else if( num == 7480) 
    return "YPA";
  else if( num == 7481) 
    return "YPR";
  else if( num == 7482) 
    return "YPN";
  else if( num == 7483) 
    return "YPD";
  else if( num == 7484) 
    return "YPC";
  else if( num == 7485) 
    return "YPQ";
  else if( num == 7486) 
    return "YPE";
  else if( num == 7487) 
    return "YPG";
  else if( num == 7488) 
    return "YPH";
  else if( num == 7489) 
    return "YPI";
  else if( num == 7490) 
    return "YPL";
  else if( num == 7491) 
    return "YPK";
  else if( num == 7492) 
    return "YPM";
  else if( num == 7493) 
    return "YPF";
  else if( num == 7494) 
    return "YPP";
  else if( num == 7495) 
    return "YPS";
  else if( num == 7496) 
    return "YPT";
  else if( num == 7497) 
    return "YPW";
  else if( num == 7498) 
    return "YPY";
  else if( num == 7499) 
    return "YPV";
  else if( num == 7500) 
    return "YSA";
  else if( num == 7501) 
    return "YSR";
  else if( num == 7502) 
    return "YSN";
  else if( num == 7503) 
    return "YSD";
  else if( num == 7504) 
    return "YSC";
  else if( num == 7505) 
    return "YSQ";
  else if( num == 7506) 
    return "YSE";
  else if( num == 7507) 
    return "YSG";
  else if( num == 7508) 
    return "YSH";
  else if( num == 7509) 
    return "YSI";
  else if( num == 7510) 
    return "YSL";
  else if( num == 7511) 
    return "YSK";
  else if( num == 7512) 
    return "YSM";
  else if( num == 7513) 
    return "YSF";
  else if( num == 7514) 
    return "YSP";
  else if( num == 7515) 
    return "YSS";
  else if( num == 7516) 
    return "YST";
  else if( num == 7517) 
    return "YSW";
  else if( num == 7518) 
    return "YSY";
  else if( num == 7519) 
    return "YSV";
  else if( num == 7520) 
    return "YTA";
  else if( num == 7521) 
    return "YTR";
  else if( num == 7522) 
    return "YTN";
  else if( num == 7523) 
    return "YTD";
  else if( num == 7524) 
    return "YTC";
  else if( num == 7525) 
    return "YTQ";
  else if( num == 7526) 
    return "YTE";
  else if( num == 7527) 
    return "YTG";
  else if( num == 7528) 
    return "YTH";
  else if( num == 7529) 
    return "YTI";
  else if( num == 7530) 
    return "YTL";
  else if( num == 7531) 
    return "YTK";
  else if( num == 7532) 
    return "YTM";
  else if( num == 7533) 
    return "YTF";
  else if( num == 7534) 
    return "YTP";
  else if( num == 7535) 
    return "YTS";
  else if( num == 7536) 
    return "YTT";
  else if( num == 7537) 
    return "YTW";
  else if( num == 7538) 
    return "YTY";
  else if( num == 7539) 
    return "YTV";
  else if( num == 7540) 
    return "YWA";
  else if( num == 7541) 
    return "YWR";
  else if( num == 7542) 
    return "YWN";
  else if( num == 7543) 
    return "YWD";
  else if( num == 7544) 
    return "YWC";
  else if( num == 7545) 
    return "YWQ";
  else if( num == 7546) 
    return "YWE";
  else if( num == 7547) 
    return "YWG";
  else if( num == 7548) 
    return "YWH";
  else if( num == 7549) 
    return "YWI";
  else if( num == 7550) 
    return "YWL";
  else if( num == 7551) 
    return "YWK";
  else if( num == 7552) 
    return "YWM";
  else if( num == 7553) 
    return "YWF";
  else if( num == 7554) 
    return "YWP";
  else if( num == 7555) 
    return "YWS";
  else if( num == 7556) 
    return "YWT";
  else if( num == 7557) 
    return "YWW";
  else if( num == 7558) 
    return "YWY";
  else if( num == 7559) 
    return "YWV";
  else if( num == 7560) 
    return "YYA";
  else if( num == 7561) 
    return "YYR";
  else if( num == 7562) 
    return "YYN";
  else if( num == 7563) 
    return "YYD";
  else if( num == 7564) 
    return "YYC";
  else if( num == 7565) 
    return "YYQ";
  else if( num == 7566) 
    return "YYE";
  else if( num == 7567) 
    return "YYG";
  else if( num == 7568) 
    return "YYH";
  else if( num == 7569) 
    return "YYI";
  else if( num == 7570) 
    return "YYL";
  else if( num == 7571) 
    return "YYK";
  else if( num == 7572) 
    return "YYM";
  else if( num == 7573) 
    return "YYF";
  else if( num == 7574) 
    return "YYP";
  else if( num == 7575) 
    return "YYS";
  else if( num == 7576) 
    return "YYT";
  else if( num == 7577) 
    return "YYW";
  else if( num == 7578) 
    return "YYY";
  else if( num == 7579) 
    return "YYV";
  else if( num == 7580) 
    return "YVA";
  else if( num == 7581) 
    return "YVR";
  else if( num == 7582) 
    return "YVN";
  else if( num == 7583) 
    return "YVD";
  else if( num == 7584) 
    return "YVC";
  else if( num == 7585) 
    return "YVQ";
  else if( num == 7586) 
    return "YVE";
  else if( num == 7587) 
    return "YVG";
  else if( num == 7588) 
    return "YVH";
  else if( num == 7589) 
    return "YVI";
  else if( num == 7590) 
    return "YVL";
  else if( num == 7591) 
    return "YVK";
  else if( num == 7592) 
    return "YVM";
  else if( num == 7593) 
    return "YVF";
  else if( num == 7594) 
    return "YVP";
  else if( num == 7595) 
    return "YVS";
  else if( num == 7596) 
    return "YVT";
  else if( num == 7597) 
    return "YVW";
  else if( num == 7598) 
    return "YVY";
  else if( num == 7599) 
    return "YVV";
  else if( num == 7600) 
    return "VAA";
  else if( num == 7601) 
    return "VAR";
  else if( num == 7602) 
    return "VAN";
  else if( num == 7603) 
    return "VAD";
  else if( num == 7604) 
    return "VAC";
  else if( num == 7605) 
    return "VAQ";
  else if( num == 7606) 
    return "VAE";
  else if( num == 7607) 
    return "VAG";
  else if( num == 7608) 
    return "VAH";
  else if( num == 7609) 
    return "VAI";
  else if( num == 7610) 
    return "VAL";
  else if( num == 7611) 
    return "VAK";
  else if( num == 7612) 
    return "VAM";
  else if( num == 7613) 
    return "VAF";
  else if( num == 7614) 
    return "VAP";
  else if( num == 7615) 
    return "VAS";
  else if( num == 7616) 
    return "VAT";
  else if( num == 7617) 
    return "VAW";
  else if( num == 7618) 
    return "VAY";
  else if( num == 7619) 
    return "VAV";
  else if( num == 7620) 
    return "VRA";
  else if( num == 7621) 
    return "VRR";
  else if( num == 7622) 
    return "VRN";
  else if( num == 7623) 
    return "VRD";
  else if( num == 7624) 
    return "VRC";
  else if( num == 7625) 
    return "VRQ";
  else if( num == 7626) 
    return "VRE";
  else if( num == 7627) 
    return "VRG";
  else if( num == 7628) 
    return "VRH";
  else if( num == 7629) 
    return "VRI";
  else if( num == 7630) 
    return "VRL";
  else if( num == 7631) 
    return "VRK";
  else if( num == 7632) 
    return "VRM";
  else if( num == 7633) 
    return "VRF";
  else if( num == 7634) 
    return "VRP";
  else if( num == 7635) 
    return "VRS";
  else if( num == 7636) 
    return "VRT";
  else if( num == 7637) 
    return "VRW";
  else if( num == 7638) 
    return "VRY";
  else if( num == 7639) 
    return "VRV";
  else if( num == 7640) 
    return "VNA";
  else if( num == 7641) 
    return "VNR";
  else if( num == 7642) 
    return "VNN";
  else if( num == 7643) 
    return "VND";
  else if( num == 7644) 
    return "VNC";
  else if( num == 7645) 
    return "VNQ";
  else if( num == 7646) 
    return "VNE";
  else if( num == 7647) 
    return "VNG";
  else if( num == 7648) 
    return "VNH";
  else if( num == 7649) 
    return "VNI";
  else if( num == 7650) 
    return "VNL";
  else if( num == 7651) 
    return "VNK";
  else if( num == 7652) 
    return "VNM";
  else if( num == 7653) 
    return "VNF";
  else if( num == 7654) 
    return "VNP";
  else if( num == 7655) 
    return "VNS";
  else if( num == 7656) 
    return "VNT";
  else if( num == 7657) 
    return "VNW";
  else if( num == 7658) 
    return "VNY";
  else if( num == 7659) 
    return "VNV";
  else if( num == 7660) 
    return "VDA";
  else if( num == 7661) 
    return "VDR";
  else if( num == 7662) 
    return "VDN";
  else if( num == 7663) 
    return "VDD";
  else if( num == 7664) 
    return "VDC";
  else if( num == 7665) 
    return "VDQ";
  else if( num == 7666) 
    return "VDE";
  else if( num == 7667) 
    return "VDG";
  else if( num == 7668) 
    return "VDH";
  else if( num == 7669) 
    return "VDI";
  else if( num == 7670) 
    return "VDL";
  else if( num == 7671) 
    return "VDK";
  else if( num == 7672) 
    return "VDM";
  else if( num == 7673) 
    return "VDF";
  else if( num == 7674) 
    return "VDP";
  else if( num == 7675) 
    return "VDS";
  else if( num == 7676) 
    return "VDT";
  else if( num == 7677) 
    return "VDW";
  else if( num == 7678) 
    return "VDY";
  else if( num == 7679) 
    return "VDV";
  else if( num == 7680) 
    return "VCA";
  else if( num == 7681) 
    return "VCR";
  else if( num == 7682) 
    return "VCN";
  else if( num == 7683) 
    return "VCD";
  else if( num == 7684) 
    return "VCC";
  else if( num == 7685) 
    return "VCQ";
  else if( num == 7686) 
    return "VCE";
  else if( num == 7687) 
    return "VCG";
  else if( num == 7688) 
    return "VCH";
  else if( num == 7689) 
    return "VCI";
  else if( num == 7690) 
    return "VCL";
  else if( num == 7691) 
    return "VCK";
  else if( num == 7692) 
    return "VCM";
  else if( num == 7693) 
    return "VCF";
  else if( num == 7694) 
    return "VCP";
  else if( num == 7695) 
    return "VCS";
  else if( num == 7696) 
    return "VCT";
  else if( num == 7697) 
    return "VCW";
  else if( num == 7698) 
    return "VCY";
  else if( num == 7699) 
    return "VCV";
  else if( num == 7700) 
    return "VQA";
  else if( num == 7701) 
    return "VQR";
  else if( num == 7702) 
    return "VQN";
  else if( num == 7703) 
    return "VQD";
  else if( num == 7704) 
    return "VQC";
  else if( num == 7705) 
    return "VQQ";
  else if( num == 7706) 
    return "VQE";
  else if( num == 7707) 
    return "VQG";
  else if( num == 7708) 
    return "VQH";
  else if( num == 7709) 
    return "VQI";
  else if( num == 7710) 
    return "VQL";
  else if( num == 7711) 
    return "VQK";
  else if( num == 7712) 
    return "VQM";
  else if( num == 7713) 
    return "VQF";
  else if( num == 7714) 
    return "VQP";
  else if( num == 7715) 
    return "VQS";
  else if( num == 7716) 
    return "VQT";
  else if( num == 7717) 
    return "VQW";
  else if( num == 7718) 
    return "VQY";
  else if( num == 7719) 
    return "VQV";
  else if( num == 7720) 
    return "VEA";
  else if( num == 7721) 
    return "VER";
  else if( num == 7722) 
    return "VEN";
  else if( num == 7723) 
    return "VED";
  else if( num == 7724) 
    return "VEC";
  else if( num == 7725) 
    return "VEQ";
  else if( num == 7726) 
    return "VEE";
  else if( num == 7727) 
    return "VEG";
  else if( num == 7728) 
    return "VEH";
  else if( num == 7729) 
    return "VEI";
  else if( num == 7730) 
    return "VEL";
  else if( num == 7731) 
    return "VEK";
  else if( num == 7732) 
    return "VEM";
  else if( num == 7733) 
    return "VEF";
  else if( num == 7734) 
    return "VEP";
  else if( num == 7735) 
    return "VES";
  else if( num == 7736) 
    return "VET";
  else if( num == 7737) 
    return "VEW";
  else if( num == 7738) 
    return "VEY";
  else if( num == 7739) 
    return "VEV";
  else if( num == 7740) 
    return "VGA";
  else if( num == 7741) 
    return "VGR";
  else if( num == 7742) 
    return "VGN";
  else if( num == 7743) 
    return "VGD";
  else if( num == 7744) 
    return "VGC";
  else if( num == 7745) 
    return "VGQ";
  else if( num == 7746) 
    return "VGE";
  else if( num == 7747) 
    return "VGG";
  else if( num == 7748) 
    return "VGH";
  else if( num == 7749) 
    return "VGI";
  else if( num == 7750) 
    return "VGL";
  else if( num == 7751) 
    return "VGK";
  else if( num == 7752) 
    return "VGM";
  else if( num == 7753) 
    return "VGF";
  else if( num == 7754) 
    return "VGP";
  else if( num == 7755) 
    return "VGS";
  else if( num == 7756) 
    return "VGT";
  else if( num == 7757) 
    return "VGW";
  else if( num == 7758) 
    return "VGY";
  else if( num == 7759) 
    return "VGV";
  else if( num == 7760) 
    return "VHA";
  else if( num == 7761) 
    return "VHR";
  else if( num == 7762) 
    return "VHN";
  else if( num == 7763) 
    return "VHD";
  else if( num == 7764) 
    return "VHC";
  else if( num == 7765) 
    return "VHQ";
  else if( num == 7766) 
    return "VHE";
  else if( num == 7767) 
    return "VHG";
  else if( num == 7768) 
    return "VHH";
  else if( num == 7769) 
    return "VHI";
  else if( num == 7770) 
    return "VHL";
  else if( num == 7771) 
    return "VHK";
  else if( num == 7772) 
    return "VHM";
  else if( num == 7773) 
    return "VHF";
  else if( num == 7774) 
    return "VHP";
  else if( num == 7775) 
    return "VHS";
  else if( num == 7776) 
    return "VHT";
  else if( num == 7777) 
    return "VHW";
  else if( num == 7778) 
    return "VHY";
  else if( num == 7779) 
    return "VHV";
  else if( num == 7780) 
    return "VIA";
  else if( num == 7781) 
    return "VIR";
  else if( num == 7782) 
    return "VIN";
  else if( num == 7783) 
    return "VID";
  else if( num == 7784) 
    return "VIC";
  else if( num == 7785) 
    return "VIQ";
  else if( num == 7786) 
    return "VIE";
  else if( num == 7787) 
    return "VIG";
  else if( num == 7788) 
    return "VIH";
  else if( num == 7789) 
    return "VII";
  else if( num == 7790) 
    return "VIL";
  else if( num == 7791) 
    return "VIK";
  else if( num == 7792) 
    return "VIM";
  else if( num == 7793) 
    return "VIF";
  else if( num == 7794) 
    return "VIP";
  else if( num == 7795) 
    return "VIS";
  else if( num == 7796) 
    return "VIT";
  else if( num == 7797) 
    return "VIW";
  else if( num == 7798) 
    return "VIY";
  else if( num == 7799) 
    return "VIV";
  else if( num == 7800) 
    return "VLA";
  else if( num == 7801) 
    return "VLR";
  else if( num == 7802) 
    return "VLN";
  else if( num == 7803) 
    return "VLD";
  else if( num == 7804) 
    return "VLC";
  else if( num == 7805) 
    return "VLQ";
  else if( num == 7806) 
    return "VLE";
  else if( num == 7807) 
    return "VLG";
  else if( num == 7808) 
    return "VLH";
  else if( num == 7809) 
    return "VLI";
  else if( num == 7810) 
    return "VLL";
  else if( num == 7811) 
    return "VLK";
  else if( num == 7812) 
    return "VLM";
  else if( num == 7813) 
    return "VLF";
  else if( num == 7814) 
    return "VLP";
  else if( num == 7815) 
    return "VLS";
  else if( num == 7816) 
    return "VLT";
  else if( num == 7817) 
    return "VLW";
  else if( num == 7818) 
    return "VLY";
  else if( num == 7819) 
    return "VLV";
  else if( num == 7820) 
    return "VKA";
  else if( num == 7821) 
    return "VKR";
  else if( num == 7822) 
    return "VKN";
  else if( num == 7823) 
    return "VKD";
  else if( num == 7824) 
    return "VKC";
  else if( num == 7825) 
    return "VKQ";
  else if( num == 7826) 
    return "VKE";
  else if( num == 7827) 
    return "VKG";
  else if( num == 7828) 
    return "VKH";
  else if( num == 7829) 
    return "VKI";
  else if( num == 7830) 
    return "VKL";
  else if( num == 7831) 
    return "VKK";
  else if( num == 7832) 
    return "VKM";
  else if( num == 7833) 
    return "VKF";
  else if( num == 7834) 
    return "VKP";
  else if( num == 7835) 
    return "VKS";
  else if( num == 7836) 
    return "VKT";
  else if( num == 7837) 
    return "VKW";
  else if( num == 7838) 
    return "VKY";
  else if( num == 7839) 
    return "VKV";
  else if( num == 7840) 
    return "VMA";
  else if( num == 7841) 
    return "VMR";
  else if( num == 7842) 
    return "VMN";
  else if( num == 7843) 
    return "VMD";
  else if( num == 7844) 
    return "VMC";
  else if( num == 7845) 
    return "VMQ";
  else if( num == 7846) 
    return "VME";
  else if( num == 7847) 
    return "VMG";
  else if( num == 7848) 
    return "VMH";
  else if( num == 7849) 
    return "VMI";
  else if( num == 7850) 
    return "VML";
  else if( num == 7851) 
    return "VMK";
  else if( num == 7852) 
    return "VMM";
  else if( num == 7853) 
    return "VMF";
  else if( num == 7854) 
    return "VMP";
  else if( num == 7855) 
    return "VMS";
  else if( num == 7856) 
    return "VMT";
  else if( num == 7857) 
    return "VMW";
  else if( num == 7858) 
    return "VMY";
  else if( num == 7859) 
    return "VMV";
  else if( num == 7860) 
    return "VFA";
  else if( num == 7861) 
    return "VFR";
  else if( num == 7862) 
    return "VFN";
  else if( num == 7863) 
    return "VFD";
  else if( num == 7864) 
    return "VFC";
  else if( num == 7865) 
    return "VFQ";
  else if( num == 7866) 
    return "VFE";
  else if( num == 7867) 
    return "VFG";
  else if( num == 7868) 
    return "VFH";
  else if( num == 7869) 
    return "VFI";
  else if( num == 7870) 
    return "VFL";
  else if( num == 7871) 
    return "VFK";
  else if( num == 7872) 
    return "VFM";
  else if( num == 7873) 
    return "VFF";
  else if( num == 7874) 
    return "VFP";
  else if( num == 7875) 
    return "VFS";
  else if( num == 7876) 
    return "VFT";
  else if( num == 7877) 
    return "VFW";
  else if( num == 7878) 
    return "VFY";
  else if( num == 7879) 
    return "VFV";
  else if( num == 7880) 
    return "VPA";
  else if( num == 7881) 
    return "VPR";
  else if( num == 7882) 
    return "VPN";
  else if( num == 7883) 
    return "VPD";
  else if( num == 7884) 
    return "VPC";
  else if( num == 7885) 
    return "VPQ";
  else if( num == 7886) 
    return "VPE";
  else if( num == 7887) 
    return "VPG";
  else if( num == 7888) 
    return "VPH";
  else if( num == 7889) 
    return "VPI";
  else if( num == 7890) 
    return "VPL";
  else if( num == 7891) 
    return "VPK";
  else if( num == 7892) 
    return "VPM";
  else if( num == 7893) 
    return "VPF";
  else if( num == 7894) 
    return "VPP";
  else if( num == 7895) 
    return "VPS";
  else if( num == 7896) 
    return "VPT";
  else if( num == 7897) 
    return "VPW";
  else if( num == 7898) 
    return "VPY";
  else if( num == 7899) 
    return "VPV";
  else if( num == 7900) 
    return "VSA";
  else if( num == 7901) 
    return "VSR";
  else if( num == 7902) 
    return "VSN";
  else if( num == 7903) 
    return "VSD";
  else if( num == 7904) 
    return "VSC";
  else if( num == 7905) 
    return "VSQ";
  else if( num == 7906) 
    return "VSE";
  else if( num == 7907) 
    return "VSG";
  else if( num == 7908) 
    return "VSH";
  else if( num == 7909) 
    return "VSI";
  else if( num == 7910) 
    return "VSL";
  else if( num == 7911) 
    return "VSK";
  else if( num == 7912) 
    return "VSM";
  else if( num == 7913) 
    return "VSF";
  else if( num == 7914) 
    return "VSP";
  else if( num == 7915) 
    return "VSS";
  else if( num == 7916) 
    return "VST";
  else if( num == 7917) 
    return "VSW";
  else if( num == 7918) 
    return "VSY";
  else if( num == 7919) 
    return "VSV";
  else if( num == 7920) 
    return "VTA";
  else if( num == 7921) 
    return "VTR";
  else if( num == 7922) 
    return "VTN";
  else if( num == 7923) 
    return "VTD";
  else if( num == 7924) 
    return "VTC";
  else if( num == 7925) 
    return "VTQ";
  else if( num == 7926) 
    return "VTE";
  else if( num == 7927) 
    return "VTG";
  else if( num == 7928) 
    return "VTH";
  else if( num == 7929) 
    return "VTI";
  else if( num == 7930) 
    return "VTL";
  else if( num == 7931) 
    return "VTK";
  else if( num == 7932) 
    return "VTM";
  else if( num == 7933) 
    return "VTF";
  else if( num == 7934) 
    return "VTP";
  else if( num == 7935) 
    return "VTS";
  else if( num == 7936) 
    return "VTT";
  else if( num == 7937) 
    return "VTW";
  else if( num == 7938) 
    return "VTY";
  else if( num == 7939) 
    return "VTV";
  else if( num == 7940) 
    return "VWA";
  else if( num == 7941) 
    return "VWR";
  else if( num == 7942) 
    return "VWN";
  else if( num == 7943) 
    return "VWD";
  else if( num == 7944) 
    return "VWC";
  else if( num == 7945) 
    return "VWQ";
  else if( num == 7946) 
    return "VWE";
  else if( num == 7947) 
    return "VWG";
  else if( num == 7948) 
    return "VWH";
  else if( num == 7949) 
    return "VWI";
  else if( num == 7950) 
    return "VWL";
  else if( num == 7951) 
    return "VWK";
  else if( num == 7952) 
    return "VWM";
  else if( num == 7953) 
    return "VWF";
  else if( num == 7954) 
    return "VWP";
  else if( num == 7955) 
    return "VWS";
  else if( num == 7956) 
    return "VWT";
  else if( num == 7957) 
    return "VWW";
  else if( num == 7958) 
    return "VWY";
  else if( num == 7959) 
    return "VWV";
  else if( num == 7960) 
    return "VYA";
  else if( num == 7961) 
    return "VYR";
  else if( num == 7962) 
    return "VYN";
  else if( num == 7963) 
    return "VYD";
  else if( num == 7964) 
    return "VYC";
  else if( num == 7965) 
    return "VYQ";
  else if( num == 7966) 
    return "VYE";
  else if( num == 7967) 
    return "VYG";
  else if( num == 7968) 
    return "VYH";
  else if( num == 7969) 
    return "VYI";
  else if( num == 7970) 
    return "VYL";
  else if( num == 7971) 
    return "VYK";
  else if( num == 7972) 
    return "VYM";
  else if( num == 7973) 
    return "VYF";
  else if( num == 7974) 
    return "VYP";
  else if( num == 7975) 
    return "VYS";
  else if( num == 7976) 
    return "VYT";
  else if( num == 7977) 
    return "VYW";
  else if( num == 7978) 
    return "VYY";
  else if( num == 7979) 
    return "VYV";
  else if( num == 7980) 
    return "VVA";
  else if( num == 7981) 
    return "VVR";
  else if( num == 7982) 
    return "VVN";
  else if( num == 7983) 
    return "VVD";
  else if( num == 7984) 
    return "VVC";
  else if( num == 7985) 
    return "VVQ";
  else if( num == 7986) 
    return "VVE";
  else if( num == 7987) 
    return "VVG";
  else if( num == 7988) 
    return "VVH";
  else if( num == 7989) 
    return "VVI";
  else if( num == 7990) 
    return "VVL";
  else if( num == 7991) 
    return "VVK";
  else if( num == 7992) 
    return "VVM";
  else if( num == 7993) 
    return "VVF";
  else if( num == 7994) 
    return "VVP";
  else if( num == 7995) 
    return "VVS";
  else if( num == 7996) 
    return "VVT";
  else if( num == 7997) 
    return "VVW";
  else if( num == 7998) 
    return "VVY";
  else if( num == 7999) 
    return "VVV";
  return "AAA";
} 
