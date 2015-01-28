/*=======================================================================
(C) Copyright 1991, Fred Hutchinson Cancer Research Center
     motmisc.c    Miscellaneous PROTOMAT routines.
-------------------------------------------------------------------------
   6/20/91  J. Henikoff
  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 6.0  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   2/3/93   Changed offset for PIR format in init_dbs().
            Changed extract_seqs() to work when database is not in sorted
      order (true for PIR). Added check_entry().
   2/18/93  Changed extract_seqs() to add PIR entry name to db_id structure
      from DR line when db type is EMBL.
            Added initialization of pir field to makedbid().
            Added pir field to get_ids().
  7/19/93   Use X instead of J for all routines: X is don't care/placeholder
            amino acid (num_to_aachar() ).
  7/19/93   Limit sequence title to MAXTITLE characters (extract_seqs()).
  7/21/93   Fixed GENBANK problem (lower case in extract_seqs()).
  7/31/93   Added strnjcmp() routine for systems without strnicmp().
 >>>>>>>>>>>>>>>>>>>>>>>>>>>> 6.2 <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  10/27/93  Changed extract_seqs() to only extract sequences that are not
            fragments or Prosite false positives.
  11/ 4/93  Changed extract_seqs() to check frag flag for extracting fragments
  11/ 6/93  Changed extract_seqs() to optionally put all sequences in one file.
            Got rid of extra new line after title.
  11/ 7/93  Changed default matrix from bl60 to bl62.
=========================================================================*/
/*DOS #include <ctype.h>
#include <io.h>
#include <dir.h>
#include <process.h>
DOS*/
#include <sys/types.h>
#include <sys/dir.h>
#include "motifj.big.h"

/*---- Global scoring matrix , order is :
      A R N D C Q E G H I L K M F P S T W Y V X   -----------*/

/* BLOSUM60 matrix with each cell offset by 8 to make all scores
   non-negative */
int bl60_highpass = 8;
char bl60_matrix[21][21]={
12, 7, 7, 6, 7, 7, 7, 8, 7, 7, 7, 7, 7, 6, 8, 9, 8, 5, 6, 8, 8, /*A*/
 7,13, 8, 7, 5, 9, 8, 6, 8, 5, 6,10, 7, 6, 6, 7, 7, 6, 6, 6, 8, /*R*/
 7, 8,14, 9, 6, 8, 8, 8, 8, 5, 5, 8, 6, 5, 6, 9, 8, 6, 6, 5, 8, /*N*/
 6, 7, 9,14, 4, 8, 9, 7, 7, 5, 5, 7, 5, 5, 6, 8, 7, 4, 6, 5, 8, /*D*/
 7, 5, 6, 4,17, 6, 5, 6, 5, 7, 7, 5, 6, 6, 6, 7, 7, 5, 5, 7, 8, /*C*/
 7, 9, 8, 8, 6,13,10, 6, 9, 6, 6, 9, 8, 5, 7, 8, 7, 6, 6, 6, 8, /*Q*/
 7, 8, 8, 9, 5,10,13, 6, 8, 5, 6, 9, 6, 5, 7, 8, 7, 6, 6, 6, 8,
 8, 6, 8, 7, 6, 6, 6,14, 6, 5, 5, 7, 5, 5, 6, 8, 7, 6, 5, 5, 8,
 7, 8, 8, 7, 5, 9, 8, 6,15, 5, 5, 7, 6, 7, 6, 7, 6, 6, 9, 5, 8,
 7, 5, 5, 5, 7, 6, 5, 5, 5,12,10, 6, 9, 8, 5, 6, 7, 6, 7,11, 8,
 7, 6, 5, 5, 7, 6, 6, 5, 5,10,12, 6,10, 8, 6, 6, 7, 6, 7, 9, 8,
 7,10, 8, 7, 5, 9, 9, 7, 7, 6, 6,12, 7, 5, 7, 8, 7, 5, 6, 6, 8,
 7, 7, 6, 5, 6, 8, 6, 5, 6, 9,10, 7,14, 8, 6, 6, 7, 6, 7, 9, 8,
 6, 6, 5, 5, 6, 5, 5, 5, 7, 8, 8, 5, 8,14, 5, 6, 6, 9,11, 7, 8,
 8, 6, 6, 6, 6, 7, 7, 6, 6, 5, 6, 7, 6, 5,15, 7, 7, 5, 6, 6, 8,
 9, 7, 9, 8, 7, 8, 8, 8, 7, 6, 6, 8, 6, 6, 7,12, 9, 6, 6, 6, 8,
 8, 7, 8, 7, 7, 7, 7, 7, 6, 7, 7, 7, 7, 6, 7, 9,12, 6, 6, 8, 8,
 5, 6, 6, 4, 5, 6, 6, 6, 6, 6, 6, 5, 6, 9, 5, 6, 6,18,10, 6, 8,
 6, 6, 6, 6, 5, 6, 6, 5, 9, 7, 7, 6, 7,11, 6, 6, 6,10,15, 7, 8,
 8, 6, 5, 5, 7, 6, 6, 5, 5,11, 9, 6, 9, 7, 6, 6, 8, 6, 7,12, 8,
 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 0  /* X */
};

/* BLOSUM62 matrix with each cell offset by 4 to make all scores
   non-negative */
int bl62_highpass = 4;
/* A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  X  */
char bl62_matrix[21][21]={
   8, 3, 2, 2, 4, 3, 3, 4, 2, 3, 3, 3, 3, 2, 3, 5, 4, 1, 2, 4, 0, 
   3, 9, 4, 2, 1, 5, 4, 2, 4, 1, 2, 6, 3, 1, 2, 3, 3, 1, 2, 1, 0, 
   2, 4,10, 5, 1, 4, 4, 4, 5, 1, 1, 4, 2, 1, 2, 5, 4, 0, 2, 1, 0, 
   2, 2, 5,10, 1, 4, 6, 3, 3, 1, 0, 3, 1, 1, 3, 4, 3, 0, 1, 1, 0, 
   4, 1, 1, 1,13, 1, 0, 1, 1, 3, 3, 1, 3, 2, 1, 3, 3, 2, 2, 3, 0, 
   3, 5, 4, 4, 1, 9, 6, 2, 4, 1, 2, 5, 4, 1, 3, 4, 3, 2, 3, 2, 0, 
   3, 4, 4, 6, 0, 6, 9, 2, 4, 1, 1, 5, 2, 1, 3, 4, 3, 1, 2, 2, 0, 
   4, 2, 4, 3, 1, 2, 2,10, 2, 0, 0, 2, 1, 1, 2, 4, 2, 2, 1, 1, 0, 
   2, 4, 5, 3, 1, 4, 4, 2,12, 1, 1, 3, 2, 3, 2, 3, 2, 2, 6, 1, 0, 
   3, 1, 1, 1, 3, 1, 1, 0, 1, 8, 6, 1, 5, 4, 1, 2, 3, 1, 3, 7, 0, 
   3, 2, 1, 0, 3, 2, 1, 0, 1, 6, 8, 2, 6, 4, 1, 2, 3, 2, 3, 5, 0, 
   3, 6, 4, 3, 1, 5, 5, 2, 3, 1, 2, 9, 3, 1, 3, 4, 3, 1, 2, 2, 0, 
   3, 3, 2, 1, 3, 4, 2, 1, 2, 5, 6, 3, 9, 4, 2, 3, 3, 3, 3, 5, 0, 
   2, 1, 1, 1, 2, 1, 1, 1, 3, 4, 4, 1, 4,10, 0, 2, 2, 5, 7, 3, 0, 
   3, 2, 2, 3, 1, 3, 3, 2, 2, 1, 1, 3, 2, 0,11, 3, 3, 0, 1, 2, 0, 
   5, 3, 5, 4, 3, 4, 4, 4, 3, 2, 2, 4, 3, 2, 3, 8, 5, 1, 2, 2, 0, 
   4, 3, 4, 3, 3, 3, 3, 2, 2, 3, 3, 3, 3, 2, 3, 5, 9, 2, 2, 4, 0, 
   1, 1, 0, 0, 2, 2, 1, 2, 2, 1, 2, 1, 3, 5, 0, 1, 2,15, 6, 1, 0, 
   2, 2, 2, 1, 2, 3, 2, 1, 6, 3, 3, 2, 3, 7, 1, 2, 2, 6,11, 3, 0, 
   4, 1, 1, 1, 3, 2, 2, 1, 1, 7, 5, 2, 5, 3, 2, 2, 4, 1, 3, 8, 0, 
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0  
};
/*=======================================================================*/
/* Number to amino acid                  */
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
/*=====================================================================*/
/* Use internal numerical representation to print amino acid: */
void pr_num_to_aa(num)
char num;
{
  switch (num) {
    case 0: printf("A"); break;
    case 1: printf("R"); break;
    case 2: printf("N"); break;
    case 3: printf("D"); break;
    case 4: printf("C"); break;
    case 5: printf("Q"); break;
    case 6: printf("E"); break;
    case 7: printf("G"); break;
    case 8: printf("H"); break;
    case 9: printf("I"); break;
    case 10: printf("L"); break;
    case 11: printf("K"); break;
    case 12: printf("M"); break;
    case 13: printf("F"); break;
    case 14: printf("P"); break;
    case 15: printf("S"); break;
    case 16: printf("T"); break;
    case 17: printf("W"); break;
    case 18: printf("Y"); break;
    case 19: printf("V"); break;
    case 20: printf("."); break;
    case -1: printf("."); break;
    default: printf("*");   /* Should never happen */
    }
}
/*======================================================================*/
void pr_num_to_aa_space(c)
char c;
{
  pr_num_to_aa(c);
  printf(" ");
}
/*=======================================================================
      getscore reads a file containing a scoring matrix and
      loads it into Score[MATSIZE][MATSIZE].  Assumes alphabet for file is
      listed on first non-blank line.
=========================================================================*/
void getscore(matrix)
struct score *matrix;
{
   char filename[FNAMELEN], line[MAXLINE], chigh[6], *ptr;
   FILE *fin, *fstp;
   int alpha[MATSIZE+10], nrows, ncols, row, col, i;

   if ((fstp = fopen("protomat.stp", "rt")) == NULL)
   {
      fin = NULL;
      strcpy(filename, "def");
   }
   else
   {
      line[0] = filename[0] = '\0';
      while(fgets(line, sizeof(line), fstp) != NULL)
      {
    if (strncmp(line, "SCORE", 5) == 0)
    {
       ptr = strtok(line, " ,\t\n");
       if (ptr != NULL)
       {
    ptr = strtok(NULL, " ,\t\n");
    if (ptr != NULL) strcpy(filename, ptr);
       }
       if ((fin = fopen(filename, "rt")) == NULL)
       {
  printf("Could not open %s, using default BLOSUM scoring matrix\n",
          filename);
    strcpy(filename, "def");
       }
    }
    else if (strncmp(line, "HIGH", 4) == 0)
    {
       ptr = strtok(line, " ,\t\n");
       if (ptr != NULL)
       {
    ptr = strtok(NULL, " ,\t\n");
    if (ptr != NULL) strcpy(chigh, ptr);
    matrix->highpass = atoi(chigh);
       }
    }
      }
      fclose(fstp);
   }

/*----------Read file until first non-blank line --------------*/
   if (fin != NULL)
   {
      printf("\nUsing scoring matrix from %s\n", filename);
      line[0] = '\0';
      while (strlen(line) < 1 && fgets(line, sizeof(line), fin) != NULL)
      ;
/*------See if the first line has characters on it ------------*/
      for (col=0; col < 30; col++) alpha[col] = -1;
      if (strstr(line, "A") != NULL)  /* This line has characters */
      {
   row = 0; /* # of alphabetic characters on the line */
   for (i=0; i<strlen(line); i++)
   {
      col = aachar_to_num(line[i]);
      if (col >= 0)
      {
         alpha[row] = col;
         row++;
      }
      else if (isalpha(line[i])) row++; /* skip over other alpha */
   }
      }
/*-------Get the data values now ------------*/
      for (row=0; row<MATSIZE; row++)
  for (col=0; col<MATSIZE; col++)
     matrix->scores[row][col] = -99;    /* Null value */
      nrows = 0;
      line[0] = '\0';
      while (fgets(line, sizeof(line), fin) != NULL)
      {
   if (strlen(line) > 1)
   {
      if (alpha[nrows] >= 0 && alpha[nrows] < MATSIZE)
      {
         row = alpha[nrows]; ncols = 0;
         ptr = strtok(line, " ,\t\n");
         while (ptr != NULL)
         {
      if (strspn(ptr, "+-0123456789") == strlen(ptr))
      {
         col = alpha[ncols];
         if (col >= 0 && col < MATSIZE)
      matrix->scores[row][col] = atoi(ptr);
         ncols++;
      }
      ptr = strtok(NULL, " ,\t\n");
         }
      }
      nrows++;
   }
      }

/*-------If some entries are still missing, assume symmetry ---------*/
      for (row=0; row<MATSIZE; row++)
      {
  for (col=0; col<MATSIZE; col++)
  {
     if (matrix->scores[row][col] == -99)
      matrix->scores[row][col] = matrix->scores[col][row];
/*     printf("%2d ", matrix->scores[row][col]); */
  }
/*  printf("\n"); */
      }
      fclose(fin);
   }
   else    /*   no input file  */
   {
      printf("\nUsing BLOSUM62 scoring matrix.\n");
      matrix->highpass = bl62_highpass;
      for (row=0; row<MATSIZE; row++)
      {
   for (col=0; col<MATSIZE; col++)
   {
      matrix->scores[row][col] = bl62_matrix[row][col];
/*      printf("%2d ", matrix->scores[row][col]);*/
   }
/*   printf("\n");*/
      }
   }
   matrix->highpass *= 100;
/*   printf("HighPass = %d", matrix->highpass);*/
}   /* end of getscore() */
/*======================================================================*/
void init_dbs(dbs)
struct db_info *dbs[];
{
   int i;

   for (i=0; i<MAXDB; i++)
   {
      dbs[i] = (struct db_info *) malloc(sizeof(struct db_info));
      dbs[i]->type = (char *) malloc(10 * sizeof(char));
      dbs[i]->start = (char *) malloc(12 * sizeof(char));
      dbs[i]->desc = (char *) malloc(12 * sizeof(char));
      dbs[i]->seq = (char *)  malloc(12 * sizeof(char));
      dbs[i]->end = (char *)  malloc(6 * sizeof(char));
   }

   dbs[GB]->type = "GENBANK";
   dbs[GB]->start = "LOCUS";
   dbs[GB]->desc = "DEFINITION";
   dbs[GB]->seq = "ORIGIN";
   dbs[GB]->end = "//";
   dbs[GB]->title_offset = 12;
   dbs[GB]->seq_offset = 10;

   dbs[PIR]->type = "PIR";
   dbs[PIR]->start = "ENTRY";
   dbs[PIR]->desc = "TITLE";
   dbs[PIR]->seq = "SEQUENCE";
   dbs[PIR]->end = "///";
   dbs[PIR]->title_offset = 16;
   dbs[PIR]->seq_offset = 8;

   dbs[EMBL]->type = "EMBL";
   dbs[EMBL]->start = "ID";
   dbs[EMBL]->desc = "DE";
   dbs[EMBL]->seq = "SQ";
   dbs[EMBL]->end = "//";
   dbs[EMBL]->title_offset = 5;
   dbs[EMBL]->seq_offset = 5;

   dbs[UNI]->type = "UNI";
   dbs[UNI]->start = ">";
   dbs[UNI]->desc = ">";
   dbs[UNI]->seq = "";
   dbs[UNI]->end = "*";
   dbs[UNI]->title_offset = 1;
   dbs[UNI]->seq_offset = 0;

   dbs[VMS]->type = "VMS";
   dbs[VMS]->start = ">";   /* first line */
   dbs[VMS]->desc = "";     /* second line */
   dbs[VMS]->seq = "";      /* third line */
   dbs[VMS]->end = "*";
   dbs[VMS]->title_offset = 4;    /* first line only */
   dbs[VMS]->seq_offset = 0;

}   /*  end of init_dbs */
/*======================================================================
      type_dbs() determines what type a database is from the allowable
      types
========================================================================*/
int type_dbs(fin, dbs)
FILE *fin;
struct db_info *dbs[];
{
   int db, i;
   char line[MAXLINE];

   db = -1;
/*---------  Figure out what type of input file it is ------------------*/
   while (db < 0 && fgets(line, sizeof(line), fin) != NULL)
      for (i=0; i<MAXDB; i++)
   if (strncmp(line, dbs[i]->start, strlen(dbs[i]->start)) == 0)
      db = i;
   if (db == VMS && line[3] != ';') /* start=='>' is ambiguous */
  db = UNI;
   if (db < 0 || db >= MAXDB) db = -1;  /* can't tell what it is */

   rewind (fin);
   return(db);
}  /* end of type_dbs */
/*====================================================================*/
/*  This is Kernighan & Ritchie's ASCII to integer conversion (p. 58) */
/*====================================================================*/
int kr_atoi(s)
char s[];
{
   int i, n, sign;

   for (i=0; s[i]==' ' || s[i]=='\n' || s[i]=='\t'; i++)
    ;
   sign = 1;
   if (s[i] == '+' || s[i] == '-')
      sign = (s[i++]=='+') ? 1 : -1;
   for (n=0; s[i] >= '0' && s[i] <= '9'; i++)
      n = 10 * n + s[i] - '0';
   return(sign * n);
}  /* end of kr_atoi */
/*=====================================================================*/
/*  This is Kernighan & Ritchie's integer to ASCII conversion (p. 60) */
/*====================================================================*/
void kr_itoa(n, s)
int n;
char s[];
{
   int c, i, j, sign;

   sign = n;
   if (sign < 0) n = -n;
   i = 0;
   do {
      s[i++] = n % 10 + '0';
   }  while ( (n /= 10) > 0);
   if (sign < 0) s[i++] = '-';
   s[i] = '\0';
   for (i=0, j=strlen(s)-1; i<j; i++, j--)
   {
      c = s[i];
      s[i] = s[j];
      s[j] = c;
   }
}   /*  end of kr_itoa */
/*=====================================================================
     Locate the directory, file name and extension in a file name
========================================================================*/
struct split_name *split_names(filename)
char *filename;
{
   struct split_name *new;
   int i, ext_len;

   ext_len = 0;
   new = (struct split_name *) malloc(sizeof(struct split_name));
   new->dir_len=new->file_len=new->name_len=0;
   i = strlen(filename);
   /*-------  Read the file name backwards ---------------*/
   while (i>=0 && (!ext_len || !new->dir_len))
   {
      /*---  Last period in string => file extension ----*/
      if (filename[i] == '.') ext_len = strlen(filename)-i;
      /*--- Last slash in string => directory -----------*/
      if (filename[i] == '/' && new->dir_len == 0) new->dir_len = i+1;
      /*--- Last colon and no slash after it => DOS directory -----*/
/*      if (filename[i] == ':' && new->dir_len == 0) new->dir_len = i+1; */
      i--;
   }
   new->file_len = strlen(filename)-new->dir_len;
   new->name_len = new->file_len - ext_len;

   return(new);
}
/*========================================================================
  dir_dos(): DOS code to get name of directory from line and
       create it if necessary
===========================================================================*/
/*char *dir_dos(line)
char *line;
{
   char tname[FNAMELEN], mem[MAXLINE], pros[FNAMELEN], *ptr;
   char filename[FNAMELEN];
   int test;
   FILE *ftmp;

   pros[0] = '\0';
   if (line[0] != '>' &&
    (strstr(line, "\\") != NULL || strstr(line, ":") != NULL) )
   {
         ptr = strtok(line, "\n\r");
         strcpy(pros, ptr);
   if (pros[strlen(pros)-1] != '\\') strcat(pros, "\\");
   }
*/
/*-------------------Create the directory ---------------------------------*/
/*
   if (strlen(pros))
   {
      tmpnam(tname);
      strcpy(filename, pros);
      strcat(filename, tname);
      if ( (ftmp=fopen(filename, "w"))== NULL)
      {
   strcpy(tname, pros);     
   tname[strlen(pros)-1] = '\0';
   sprintf(mem, "md %s", tname);
   test = system(mem);
   if (test == 0) printf("\nCreated directory %s", tname);
   else
   {
      printf("\nUnable to create directory %s", tname);
      printf("\nProtein files will be placed in current directory");
      pros[0] = '\0';
   }
      }
      else
      {
   fclose(ftmp);
   unlink(filename);
      }
   }
   return(pros);
} 
*/
/*=======================================================================
  dir_unix(): UNIX code to get name of directory from line and
        create it if necessary
==========================================================================*/
char *dir_unix(line)
char *line;
{
   char tname[FNAMELEN], mem[MAXLINE], pros[FNAMELEN], *ptr;
   int test;
   DIR *dp;

   pros[0] = '\0';
   if (line[0] != '>' && (strstr(line, "/") != NULL != NULL) )
   {
         ptr = strtok(line, "\n\r");
         strcpy(pros, ptr);
   if (pros[strlen(pros)-1] != '/') strcat(pros, "/");
   }
/*-----------------Create the directory ---------------------------------*/
   if (strlen(pros))
   {
      strcpy(tname, pros);
      tname[strlen(pros)-1] = '\0';
      sprintf(mem, "mkdir %s", tname);
      if ((dp=opendir(tname))==NULL)
      {
         test=system(mem);
         if (test == 0) printf("\nCreated directory %s\n", tname);
         else
         {
            printf("\nUnable to create directory %s", tname);
            printf("\nProtein files will be placed in current directory\n");
            pros[0] = '\0';
         }
      }
   }
   return(pros);
}    /* end of dir_unix */
/*======================================================================
     Create & intialize a db_id structure
========================================================================*/
struct db_id *makedbid()
{
   struct db_id *new;

   new = (struct db_id *) malloc (sizeof(struct db_id));
   new->entry[0] = '\0';
   new->pir[0] = '\0';
   new->ps[0] = '\0';
   new->len = 0;
   new->rank = new->score = 0;
   new->lst = NO;
   new->found = NO;
   new->block = NO;
   new->frag = NO;
   new->search = NO;
   new->pvalue = (double) 0.0;
   new->next = NULL;
   new->prior = NULL;
   return(new);
}  /*  end of makedbid */

/*======================================================================
     get_ids() reads a .lis or .lst file & inserts the sequences
       found in it into a list sorted by sequence name.
========================================================================*/
int get_ids(flis, ids)
FILE *flis;
struct db_id *ids;
{
   char line[MAXLINE], ctemp[10], *ptr;
   struct db_id *id, *last, *new;
   int len, nids = 0;

   while(!feof(flis) && fgets(line, MAXLINE, flis) != NULL &&
   strlen(line) > 2)
   {    /* skip over title or directory lines */
      if (line[0] != '>' && strstr(line, "/") == NULL &&
    strstr(line,"\\") == NULL && strstr(line, ":") == NULL)
      {
   nids += 1;
/*-----  Copy up to the first space or carriage return ------------*/
   len = strcspn(line, " \t\r\n");
   if (len > IDLEN) len = IDLEN;  /* No id should be longer that this*/
   new = makedbid();
   strncpy(new->entry, line, len);  new->entry[len] = '\0';
   last = ids;  id = ids->next;
/*-------- Get any other information on the .lis file line --------*/
   if (strstr(line, "FRAGMENT") != NULL) new->frag = YES;
   if (strstr(line, "LST") != NULL) new->lst = YES;
   ptr = strstr(line, "PS=");
   if (ptr != NULL)
      strncpy(new->ps, ptr+3, 1); new->ps[1] = '\0';
   ptr = strstr(line, "LENGTH=");
   if (ptr != NULL)
   {
      len = strcspn(ptr+7, " \t\r\n");
      if (len > 0)
      {
         strncpy(ctemp, ptr+7, len); ctemp[len] = '\0';
         new->len = atoi(ctemp);
      }
   }
   ptr = strstr(line, "PIR=");
   if (ptr != NULL)
   {
      len = strcspn(ptr+4, " \t\r\n");
      if (len > 0)
      {
         strncpy(new->pir, ptr+4, len); new->pir[len] = '\0';
      }
   }
/*------  Insert id into a sorted list ----------------------------*/
   while (id != NULL && id->entry != NULL &&
       strcmp(id->entry, new->entry) < 0)
   {
      last = id;
      id = id->next;
   }
   new->prior = last;
   new->next = id;
   last->next = new;
   if (id != NULL) id->prior = new;
   new = NULL;
      }
   }
   return(nids);
}  /*  end of get_ids */
/*======================================================================
     Extracts the sequences from file fin if they appear in the
     sorted list ids.
=======================================================================*/
int extract_seqs(nids, dbs, fin, ids, pros, fout, frag)
int nids, frag;
struct db_info *dbs[MAXDB];
FILE *fin, *fout;
struct db_id *ids;
char *pros;
{
   struct db_id *check_entry();
   struct db_id *id;
   int nseq, i, db, out, here, nout;
   char line[MAXLINE], title[MAXLINE], temp[MAXLINE], *ptr;
   char foutname[FNAMELEN], entry[IDLEN+1];

   nseq = nout = 0;
   if (fout == NULL) here = YES;  /* create output files here */
   else here = NO;    
   db = type_dbs(fin, dbs);
   if (db < 0 || db >= MAXDB)
   {
      printf("\nCannot determine type of input file");
      return(-1);
   }
   printf("\nProcessing input file as %s", dbs[db]->type);
   if (db==VMS)
   {
      printf("\n WARNING: Titles are sometimes truncated in this format;");
      printf("\n          I may not be able to distinguish fragments.");
   }


   /*--------Get the first record  -------------------------------*/
   fgets(line, MAXLINE, fin);
   /*---------Read and process the rest of the file -----------*/
   while (!feof(fin) && nseq < nids)
   {
      if (strncmp(line, dbs[db]->start, strlen(dbs[db]->start)) == 0)
      {
   /*----- Make VMS work like UNI by getting the 2nd title line --*/
   if (db == VMS )  /* get the 2nd line, too */
   {
     ptr=strtok(line, "\n\r");    /* get rid of CRLF */
     strcat(line, " ");
     if (fgets(temp, MAXLINE, fin) != NULL) strcat(line, temp);
   }

   /*--------Check to see if this sequence is in the list --------*/
   strncpy(entry, &line[dbs[db]->title_offset], IDLEN+1);
   entry[IDLEN+1] = '\0';
   id = check_entry(ids, entry);
   if (id != NULL)    /*  want this entry */
   {
      nseq++;  id->found = YES; id->len = 0;
      /*----  Check to see if this sequence is a fragment ----- */
      if (strstr(line, "FRAGMENT") != NULL ||
       strstr(line, "fragment") != NULL) id->frag = 1;
           /*----- Don't create an output file if the sequence is a
                   fragment or a Prosite "possible" ----------------*/
            out = YES;
            if (!frag &&
                (id->frag || (strcmp(id->ps, "P") == 0))) out = NO;
     /*------ Create the output file ---------*/
     /*------    NOTE: since id->entry is 10 chars
       and DOS file names are only 8, DOS file name may be truncated
       and therefore not unique --------*/
            if (out)
            {
               nout++;
         printf("\n%d. Entry %s found...", nseq, id->entry);
               if (here)    /* separate file */
               {  
            strcpy(foutname, pros);
            strcpy(temp, id->entry);  temp[SNAMELEN-1] = '\0';
            strcat(foutname, temp);
            if (db == GB) strcat(foutname, ".dna");
            else strcat(foutname, ".pro");
            printf("Creating %s", foutname);
            /*---Open file: should check whether it already exists... */
            if ( (fout = fopen(foutname, "w+t")) == NULL)
            {
         printf("\nCannot open %s\n", foutname);
           return(-1);
            }
               }
            }
            else
            {
         printf("\n%d. Entry %s found...but not extracted",
                       nseq, id->entry);
               printf(" because FRAGMENT or PROSITE=P");
            }

      /*------- Set up the title, write it out if UNI/VMS ----*/
      if (db == UNI || db == VMS)
      {
       strcpy(temp, &line[dbs[db]->title_offset]);
                   if (strlen(temp) > MAXTITLE)
                   { temp[MAXTITLE-1] = '\n'; temp[MAXTITLE] = '\0'; }
                   if (out) fprintf(fout, ">%s", temp);
      }
      else    /*  EMBL, GENBANK, PIR */
      {
      ptr = strtok(&line[dbs[db]->title_offset], " ");
      strcpy(title, ptr);
      }

      /*------ Process other records for this entry -------*/
      while (!feof(fin) && fgets(line, MAXLINE, fin) != NULL &&
     strncmp(line, dbs[db]->end, strlen(dbs[db]->end)) != 0 &&
     strncmp(line, dbs[db]->start, strlen(dbs[db]->start)) != 0)
      {
         /* ---For UNI & VMS, just write it back out & count length- */
         if (db == UNI || db == VMS)
         {
         if (out) fputs(line, fout);
         for (i=0; i<strlen(line); i++)
                       if (isalpha(line[i])) id->len = id->len + 1;
         }
         else    /* EMBL, GENBANK, PIR */
         {
       /*-----There might be more title to get ---------*/
       if (strncmp(line,dbs[db]->desc,strlen(dbs[db]->desc))==0 &&
       strlen(title) + strlen(line) < MAXLINE )
       {
      strcpy(temp, title); strcat(temp, " ");
      strcat(temp, &line[dbs[db]->title_offset]);
      ptr = strtok(temp, "\n\r");   /* remove CRLF */
      strcpy(title, temp);
       }
                   /*--- Kludge to get PIR entry name from EMBL files ---*/
                   else if (db == EMBL && strncmp(line, "DR   PIR", 8) == 0)
                   {
                       ptr = strtok(line+9, ";");
                       if (ptr != NULL)
                       {
                          ptr = strtok(NULL, ".\n\r");
                          if (ptr != NULL) 
                          {
                             if (ptr[0] == ' ') strcpy(id->pir, ptr+1);
                             else strcpy(id->pir, ptr);
                          }
                       }
                   }
       /*------ Process the sequence ------------------------*/
       if (strncmp(line,dbs[db]->seq,strlen(dbs[db]->seq))==0)
       {
       /*- Check again to see if this sequence is a fragment -*/
          if (strstr(title, "FRAGMENT") != NULL ||
          strstr(title, "fragment") != NULL) id->frag = 1;
          /*-- Write out the title ---*/
                      strcpy(temp, title);
                      if (strlen(temp) > MAXTITLE) temp[MAXTITLE] = '\0';
          if (out) fprintf(fout, ">%s\n", temp);  /*print title*/
          /*---  Write out the sequence & count its length --*/
          while(!feof(fin) &&
         fgets(line, MAXLINE, fin) != NULL &&
        strncmp(line,dbs[db]->end,strlen(dbs[db]->end))!=0)
          {
                          if (strpbrk(&line[dbs[db]->seq_offset],
                              "ARNDCQEGHILKMFPSTWYVBZXactg") != NULL)
                          {
                 if (out) fprintf(fout, &line[dbs[db]->seq_offset]);
                             for (i=dbs[db]->seq_offset; i<strlen(line); i++)
                               if (isalpha(line[i])) id->len = id->len + 1;
                          }
          }  /*  end of sequence lines*/
      }  /* end of if ->seq */
         }   /*  end of if not UNI and not VMS */
      }  /*  end of records for entry */
      if (out && here) fclose(fout);

      /*---- Extra processing for VMS -------*/
      if (db == VMS &&
    strncmp(line, dbs[db]->start, strlen(dbs[db]->start)) == 0)
      {                           /* get the 2nd line, too */
       ptr = strtok(line, "\n\r");        /* remove CRLF */
       strcat(line, " ");
       if (fgets(temp, MAXLINE, fin) != NULL) strcat(line, temp);
      }
   }  /*  end of found entry */
   else fgets(line, MAXLINE, fin);
      }   /*  end of if ->start */
      else fgets(line, MAXLINE, fin);
   }  /*  end of db file */

   return(nout);
}  /* end of extract_seqs */
/*====================================================================
    See if a an entry name is in a list of entries
======================================================================*/
struct db_id *check_entry(ids, entry)
struct db_id *ids;
char *entry;
{
   struct db_id *id;

   id = ids->next;
   while (id != NULL)
   {
      if (strncmp(id->entry, entry, strlen(id->entry)) == 0 &&
    id->found == NO)
        return(id);
      id = id->next;
   }
   return(NULL);
}  /* end of check_entry */

/*==================================================================
  strnjcmp.c
  Case in-sensitive version of strncmp for alphabetic characters
    A-Z are decimal 065-090
    a-z are decimal 097-122  (A+32, etc.)
=====================================================================*/
int strnjcmp(st1, st2, nchar)
char *st1, *st2;
int nchar;
{
   int i;
   char up1[2], up2[2];

   up1[1] = up2[1] = '\0';
   /*--------Convert st1 and st2 to upper case for comparison ---------*/
   for (i=0; i<nchar; i++)
   {
      up1[0] = st1[i]; up2[0] = st2[i];
      if (up1[0] >= 97 && up1[0] <= 122) up1[0] = up1[0] - 32;
      if (up2[0] >= 97 && up2[0] <= 122) up2[0] = up2[0] - 32;
      if (up1[0] != up2[0]) return(up1[0] - up2[0]);
   }
   return(0);
}