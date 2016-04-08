/*=======================================================================
(C) Copyright 1999, Fred Hutchinson Cancer Research Center
     protomat.c
     motmisc.c    Miscellaneous PROTOMAT routines.
	minus those routines duplicated in blimps
-------------------------------------------------------------------------
   2/11/99  J. Henikoff
   6/25/00  Updates to get_ids() and check_entry() for multi-part seq ids
=========================================================================*/
#include <sys/types.h>
#include <dirent.h>
#include <ctype.h>
/*	blimps library header, has stdio, etc.  */
#include <protomat.h>

/*---- Global scoring matrix , order is :
		  A R N D C Q E G H I L K M F P S T W Y V X   -----------*/

/* BLOSUM60 matrix with each cell offset by 8 to make all scores
   non-negative */
int bl60_highpass = 8;
char bl60_matrix[21][21]={
 {12, 7, 7, 6, 7, 7, 7, 8, 7, 7, 7, 7, 7, 6, 8, 9, 8, 5, 6, 8, 8}, /*A*/
 { 7,13, 8, 7, 5, 9, 8, 6, 8, 5, 6,10, 7, 6, 6, 7, 7, 6, 6, 6, 8}, /*R*/
 { 7, 8,14, 9, 6, 8, 8, 8, 8, 5, 5, 8, 6, 5, 6, 9, 8, 6, 6, 5, 8}, /*N*/
 { 6, 7, 9,14, 4, 8, 9, 7, 7, 5, 5, 7, 5, 5, 6, 8, 7, 4, 6, 5, 8}, /*D*/
 {7, 5, 6, 4,17, 6, 5, 6, 5, 7, 7, 5, 6, 6, 6, 7, 7, 5, 5, 7, 8}, /*C*/
 {7, 9, 8, 8, 6,13,10, 6, 9, 6, 6, 9, 8, 5, 7, 8, 7, 6, 6, 6, 8}, /*Q*/
 {7, 8, 8, 9, 5,10,13, 6, 8, 5, 6, 9, 6, 5, 7, 8, 7, 6, 6, 6, 8},
 {8, 6, 8, 7, 6, 6, 6,14, 6, 5, 5, 7, 5, 5, 6, 8, 7, 6, 5, 5, 8},
 {7, 8, 8, 7, 5, 9, 8, 6,15, 5, 5, 7, 6, 7, 6, 7, 6, 6, 9, 5, 8},
 {7, 5, 5, 5, 7, 6, 5, 5, 5,12,10, 6, 9, 8, 5, 6, 7, 6, 7,11, 8},
 {7, 6, 5, 5, 7, 6, 6, 5, 5,10,12, 6,10, 8, 6, 6, 7, 6, 7, 9, 8},
 {7,10, 8, 7, 5, 9, 9, 7, 7, 6, 6,12, 7, 5, 7, 8, 7, 5, 6, 6, 8},
 {7, 7, 6, 5, 6, 8, 6, 5, 6, 9,10, 7,14, 8, 6, 6, 7, 6, 7, 9, 8},
 {6, 6, 5, 5, 6, 5, 5, 5, 7, 8, 8, 5, 8,14, 5, 6, 6, 9,11, 7, 8},
 {8, 6, 6, 6, 6, 7, 7, 6, 6, 5, 6, 7, 6, 5,15, 7, 7, 5, 6, 6, 8},
 {9, 7, 9, 8, 7, 8, 8, 8, 7, 6, 6, 8, 6, 6, 7,12, 9, 6, 6, 6, 8},
 {8, 7, 8, 7, 7, 7, 7, 7, 6, 7, 7, 7, 7, 6, 7, 9,12, 6, 6, 8, 8},
 {5, 6, 6, 4, 5, 6, 6, 6, 6, 6, 6, 5, 6, 9, 5, 6, 6,18,10, 6, 8},
 {6, 6, 6, 6, 5, 6, 6, 5, 9, 7, 7, 6, 7,11, 6, 6, 6,10,15, 7, 8},
 {8, 6, 5, 5, 7, 6, 6, 5, 5,11, 9, 6, 9, 7, 6, 6, 8, 6, 7,12, 8},
 {8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 0}  /* X */
};

/* BLOSUM62 matrix with each cell offset by 4 to make all scores
   non-negative */
int bl62_highpass = 4;
/* A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  X  */
char bl62_matrix[21][21]={
{8, 3, 2, 2, 4, 3, 3, 4, 2, 3, 3, 3, 3, 2, 3, 5, 4, 1, 2, 4, 0}, 
{3, 9, 4, 2, 1, 5, 4, 2, 4, 1, 2, 6, 3, 1, 2, 3, 3, 1, 2, 1, 0}, 
{2, 4,10, 5, 1, 4, 4, 4, 5, 1, 1, 4, 2, 1, 2, 5, 4, 0, 2, 1, 0}, 
{2, 2, 5,10, 1, 4, 6, 3, 3, 1, 0, 3, 1, 1, 3, 4, 3, 0, 1, 1, 0}, 
{4, 1, 1, 1,13, 1, 0, 1, 1, 3, 3, 1, 3, 2, 1, 3, 3, 2, 2, 3, 0}, 
{3, 5, 4, 4, 1, 9, 6, 2, 4, 1, 2, 5, 4, 1, 3, 4, 3, 2, 3, 2, 0}, 
{3, 4, 4, 6, 0, 6, 9, 2, 4, 1, 1, 5, 2, 1, 3, 4, 3, 1, 2, 2, 0}, 
{4, 2, 4, 3, 1, 2, 2,10, 2, 0, 0, 2, 1, 1, 2, 4, 2, 2, 1, 1, 0}, 
{2, 4, 5, 3, 1, 4, 4, 2,12, 1, 1, 3, 2, 3, 2, 3, 2, 2, 6, 1, 0}, 
{3, 1, 1, 1, 3, 1, 1, 0, 1, 8, 6, 1, 5, 4, 1, 2, 3, 1, 3, 7, 0}, 
{3, 2, 1, 0, 3, 2, 1, 0, 1, 6, 8, 2, 6, 4, 1, 2, 3, 2, 3, 5, 0}, 
{3, 6, 4, 3, 1, 5, 5, 2, 3, 1, 2, 9, 3, 1, 3, 4, 3, 1, 2, 2, 0}, 
{3, 3, 2, 1, 3, 4, 2, 1, 2, 5, 6, 3, 9, 4, 2, 3, 3, 3, 3, 5, 0}, 
{2, 1, 1, 1, 2, 1, 1, 1, 3, 4, 4, 1, 4,10, 0, 2, 2, 5, 7, 3, 0}, 
{3, 2, 2, 3, 1, 3, 3, 2, 2, 1, 1, 3, 2, 0,11, 3, 3, 0, 1, 2, 0}, 
{5, 3, 5, 4, 3, 4, 4, 4, 3, 2, 2, 4, 3, 2, 3, 8, 5, 1, 2, 2, 0}, 
{4, 3, 4, 3, 3, 3, 3, 2, 2, 3, 3, 3, 3, 2, 3, 5, 9, 2, 2, 4, 0}, 
{1, 1, 0, 0, 2, 2, 1, 2, 2, 1, 2, 1, 3, 5, 0, 1, 2,15, 6, 1, 0}, 
{2, 2, 2, 1, 2, 3, 2, 1, 6, 3, 3, 2, 3, 7, 1, 2, 2, 6,11, 3, 0}, 
{4, 1, 1, 1, 3, 2, 2, 1, 1, 7, 5, 2, 5, 3, 2, 2, 4, 1, 3, 8, 0}, 
{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}  
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
    default: return("*");		/* Should never happen */
    }
}
/*======================================================================*/
/* Amino acid to number                  */
/*  B is changed to D, Z is changed to E, O and J are changed to X  */
int aachar_to_num(ch)
char ch;
{
  if (ch >= 97 && ch <= 122) ch = ch - 32;	/* Make letters upper case */
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
    default: printf("*");		/* Should never happen */
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
   FILE *fin=NULL, *fstp;
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
      if (strstr(line, "A") != NULL)	/* This line has characters */
      {
	 row = 0;	/* # of alphabetic characters on the line */
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
	   matrix->scores[row][col] = -99;		/* Null value */
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
/*	   printf("%2d ", matrix->scores[row][col]); */
	}
/*	printf("\n"); */
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
/*	    printf("%2d ", matrix->scores[row][col]);*/
	 }
/*	 printf("\n");*/
      }
   }
   matrix->highpass *= 100;
/*   printf("HighPass = %d", matrix->highpass);*/
}   /* end of getscore() */
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
   /* pros changed to static local variable which, who's address can be 
    * returned */
   char tname[FNAMELEN], mem[MAXLINE], *ptr;
   static char pros[FNAMELEN];
   int test;
   DIR *dp;

   pros[0] = '\0';
   if (line[0] != '>' && (strstr(line, "/") != NULL) )
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
   new->entry[0] = new->full_entry[0] = '\0';
   new->ps[0] = '\0';
   new->info[0] = '\0';
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
   char line[MAXLINE], ctemp[20], *ptr;
   struct db_id *id, *last, *new;
   int len, nids, i, bar;

   nids = 0;
   while(!feof(flis) && fgets(line, MAXLINE, flis) != NULL)
   {		/* skip over title or directory lines */
      if (strlen(line) && 
          line[0] != '>' && strstr(line, "/") == NULL &&
	  strstr(line,"\\") == NULL && strstr(line, ":") == NULL)
      {
	 nids += 1;
	 new = makedbid();
/*-----  Copy up to the first space or carriage return ------------*/
	 len = strcspn(line, " \t\r\n");
	 if (len > 2*SNAMELEN) len = 2*SNAMELEN;
	 strncpy(new->full_entry, line, len);  new->full_entry[len] = '\0';
	 if (len > SNAMELEN) len = SNAMELEN;
	 strncpy(new->entry, line, len);  new->entry[len] = '\0';
	 /* Only want max of two parts in id->entry   */
         bar=NO; i=0;
	 for (i=0; i<len; i++)
         {
		if (new->entry[i] == '|')
		{
		   if (!bar) { bar=YES;   }
		   else { new->entry[i] = '\0'; i = len;  }
		}
         }
	 last = ids;  id = ids->next;
/*-------- Get any other information on the .lis file line --------*/
	 if (strstr(line, "FRAGMENT") != NULL) new->frag = YES;
	 if (strstr(line, "BLOCK") != NULL) new->block = YES;
	 if (strstr(line, "LST") != NULL) new->lst = YES;
	 ptr = strstr(line, "PS=");
	 if (ptr != NULL)
         {
	    strncpy(new->ps, ptr+3, 1); new->ps[1] = '\0';
            if (strcmp(new->ps, "P") == 0) new->frag = YES;
         }
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
/*====================================================================
    See if a an entry name is in a list of entries
======================================================================*/
struct db_id *check_entry(ids, entry)
struct db_id *ids;
char *entry;
{
   struct db_id *id;
   int i, tot, bar, len;
   char *ptr, idtemp[SNAMELEN*6 + 6], dbtemp[SNAMELEN*6 + 6];
   char ctemp[SNAMELEN*6 + 6];
   char idlist[10][SNAMELEN*6 + 6];

   id = ids->next;
   while (id != NULL)
   {
      if (id->found == NO)
      {
         if (strcmp(id->entry, "-") == 0)  /* just take the first sequence */
         {
             return(id);
         }
         else
         {
            /*	Store the id names because can't use strtok
		for two things simultaneously    */
            strcpy(idtemp, id->entry);
            tot = 0;
            ptr = strtok(idtemp, "| \0\r\n\t");
            while(ptr != NULL)
            {
	      if (tot < 10) strcpy(idlist[tot++], ptr);
              ptr = strtok(NULL, "| \0\r\n\t");
            }
	    len = strlen(entry);
	    if (len > SNAMELEN*6) len = SNAMELEN*6;
            strncpy(dbtemp, entry, len); dbtemp[len] = '\0';
            ptr = strtok(dbtemp, "| \0\r\n\t");
            while(ptr != NULL)
            {
               for (i=0; i< tot; i++)
               {
	         if (strcmp(ptr, idlist[i]) == 0) 
                 {
                   strncpy(ctemp, entry, len); ctemp[len] = '\0'; 
	           /* Only want max of two parts in id->entry   */
                   bar=NO; i=0;
	           for (i=0; i<strlen(ctemp); i++)
                   {
		          if (ctemp[i] == '|')
		          {
		             if (!bar) { bar=YES;   }
		             else { ctemp[i] = '\0'; i = len;  }
		          }
                   }
                    strcpy(id->full_entry, ctemp);
                    return(id);
                 }
               }
               ptr = strtok(NULL, "| \0\r\n\t");
            }
         }
      }
      id = id->next;
   }
   return(NULL);
}  /* end of check_entry */
