#ifndef BLOCK_H
#define BLOCK_H

#include "configurations.h"

class Block
{

  char ac[10];
  int nseq, width, strength, nclus;
  int aa[MAXSEQS][MAX_MERGE_WIDTH];    /* aas for seq */
  double weight[MAXSEQS];              /* seq weights found in block */
  int cluster[MAXSEQS];                /* cluster # for seq */
  int ncluster[MAXSEQS];               /* #seqs in same cluster */
  double totdiag, totoffd, wtot;

public:
    Block();

    //  Getters
    char* getAc();
    int getNseq();
    int getWidth();
    int getStrength();
    int getNclus();
    int** getAa();
    double* getWeight();
    int* getCluster();
    int* getNcluster();
    double getTotdiag();
    double getTotoffd();
    double getWtot();

    //  Setters
    void setAc(char* ac);
    void setNSeq(int nseq);
    void setWidth(int width);
    void setStrength(int strength);
    void setNclus(int nclus);
    void setAa(int** aa);
    void setWeight(double* weight);
    void setCluster(int* cluster);
    void setNcluster(int* ncluster);
    void setTotdiag(double totdiag);
    void setTotoffd(double totoffd);
    void setWtot(double wtot);


};

#endif // BLOCK_H
