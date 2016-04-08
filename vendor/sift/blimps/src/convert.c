/* (C) Copyright 1993-2001, Fred Hutchinson Cancer Research Center */
/* Use, modification or distribution of these programs is subject to */
/* the terms of the non-commercial licensing agreement in license.h. */

/* convert.c: functions for different methods of converting a block into a */
/*            matrix */
/* Written by: Bill Alford */
/* Change log information is at the end of the file. */


/*	system headers not in global.h */
#include <assert.h>
#include <math.h>
/*	blimps library headers */
#include <global.h>
#include <options.h>
#include <blocks.h>	/* includes sequences.h, output.h */
#include <matrix.h>	/* includes pattern.h */
#include <residues.h>
#include <frequency.h>
#include <protomat.h>	/* for MAXLINE */
#include <convert.h>
/*	headers in current directory */
#include "blimps.h"


/*
 * Exported variables and data structures
 */
double RTot;
struct float_qij *Qij;

/*
 * Local variables and data structures
 */


/*
 * Function definitions
 */


/*
 * --OLD-- Sequence weighting methods.  Various methods of giving different
 * weights to the sequences.
 */

static double *clustered_weights(/*block*/);
static double *pre_weighted_sequences(/*block*/);



/*
 * --OLD-- Matrix construction methods.  Build matricies from blocks and
 * sequence weights.
 */

static void basic_matrix_construction(/*block, seq_weight, matrix*/);







/*
 * Conversion methods
 *   void default_conversion_method(block, matrix)
 *   void original_conversion_method_cleaned_up(block, matrix)
 *   void original_conversion_method(block, matrix)
 *   void pre_weighted_conversion_method(block, matrix)
 *   void altschul_data_dependent_conversion_method(block, matrix)
 *   void gribskov_conversion_method(block, matrix)
 */




/*
 * block_to_matrix
 *   converts a block into a matrix (possition specific matrix?) according
 *   to the rule specified in BlockToMatrixConversionMethod.  The block field
 *   of the matrix is set in this function.
 *   Parameters:
 *     Block *block: the block to convert
 *   Return codes: Returns the pointer to the new Matrix.
 *   Error codes:
 */

Matrix *block_to_matrix(block, conversion_method)
     Block *block;
     int conversion_method;
{
  Matrix *matrix;
  char *tmp;

  /* get new matrix */
  matrix = new_matrix(block->width);

  /* initialize the pattern */
  matrix->patterns = NULL;

  /* copy the relevant block information into the matrix */
  matrix->block = block;

  strncpy(Buffer, block->id, SMALL_BUFF_LENGTH);
  /* NOTE: Chance to goof by replacing the wrong "BLOCK" */
  tmp = strstr(Buffer, "BLOCK");
  if (tmp != NULL) {
    strncpy(strstr(Buffer, "BLOCK"), "MATRIX\0", 7);
    strncpy(matrix->id, Buffer, SMALL_BUFF_LENGTH);
  }
  else {
    strncpy(matrix->id,
	    strncat(Buffer, "; MATRIX", SMALL_BUFF_LENGTH - strlen(Buffer)),
	    SMALL_BUFF_LENGTH);
  }

  strcpy(matrix->ac, block->ac);
  strncpy(matrix->de, block->de, DESC_WIDTH);
  strncpy(matrix->ma, block->bl, SMALL_BUFF_LENGTH);
  /* Just leave it BLxxxxx */
  strcpy(matrix->number, block->number);
  strncpy(matrix->motif, block->motif, 20);
  /* This 20 is from the size in the struct */

  matrix->width = block->width;
  matrix->max_length = matrix->width;
  matrix->percentile = block->percentile;
  matrix->strength = block->strength;
  matrix->num_sequences = block->num_sequences;


  switch (conversion_method) {
  case 0:   /* seq wts from clumps, odds ratios */
    original_conversion_method_cleaned_up(block, matrix);
    break;
  case 1:   /* seq wts from clumps, odds ratios */
    original_conversion_method(block, matrix);
    break;
  case 2:   /* seq wts from block else pb_weights(), odds ratios */
    pre_weighted_conversion_method(block, matrix);
    break;

  /*  DEFAULT is case 3 (scale=0)  */
  /* cases 3-40 use seq wts from block else pb_weights()  */
  /* cases 3-6, 10, 20 use position-specific pseudo counts */
  /*  require frequency[] = default.amino.frq, RTot and  */
  /*  Qij = default.qij which must be initialized by calling program */
  case 3:   /* Altschul's data-dependent method of computing a PSSM */
	    /* log_e(odds ratios)=nats, positive values */
    altschul_data_dependent_conversion_method(block, matrix, 0);
    break;
  case 4:   /* log_2(odds ratios)=bits, signed integers */
    altschul_data_dependent_conversion_method(block, matrix, 1);
    break;
  case 5:   /* log_2(odds ratios)/2= half bits, signed integers */
    altschul_data_dependent_conversion_method(block, matrix, 2);
    break;
  case 6:   /* log_2(odds ratios)/3= third bits, signed integers */
    altschul_data_dependent_conversion_method(block, matrix, 3);
    break;
  case 9:   /* pseudo-counts close to zero */
	    /* log_e(odds ratios)=nats, positive values */
    altschul_data_dependent_conversion_method(block, matrix, 9);
    break;

  case 10:   /* odds ratios */
    altschul_data_dependent_conversion_method(block, matrix, 10);
    break;
  case 20:   /* (count[aa]+pseudo-count[aa])/(totcount+totpseudo) */
    altschul_data_dependent_conversion_method(block, matrix, 20);
    break;
  case 21:   /* count[aa]/totcount */
    altschul_data_dependent_conversion_method(block, matrix, 21);
    break;

  case 30:   /* Gribskov's average score method; loads default.sij */
    gribskov_conversion_method(block, matrix);
    break;

  case 40:   /* SIFT method; requires default.rank, default.diri   */
	     /* which are loaded here */
    SIFT_conversion_method(block, matrix);
    break;

  default:  /* the default case */

    sprintf(ErrorBuffer,
	    "Invalid block to matrix conversion method specified, %d.",
	    conversion_method);
    ErrorReport(WARNING_ERR_LVL);

    sprintf(ErrorBuffer,     /* ^^^^----------------vvvvvvvvvvvvvvvvvvv */
	    "Using the default conversion method of Altschul's data-dependent method.\n");
    ErrorReport(WARNING_ERR_LVL);
    altschul_data_dependent_conversion_method(block, matrix, 0);
    break;

  } /* end switch of conversion types */

  /* return the matrix */
  return matrix;
}



/*
 * original_conversion_method_cleaned_up
 *   The original conversion method.  This is done by weighted average of the
 *   clusters.  This follows the method in patmat but has been written in
 *   a more legible manner and has the following changes:
 *     no more reliance on flag values for B, Z, and X
 *     X, '-', '*', & non-code scores are read straight from the frequencies
 *     the frequencies for B and Z are ignored, when a B or Z is encountered
 *       it is partitioned between D & N or E & Q.
 *     the matrix scores for B and Z are computed from the matrix scores of
 *       D & N and E & Q.
 *   Parameters:
 *     Block *block:   the block to be converted
 *     Matrix *matrix: where the resulting matrix will be put
 *   Error codes: none
 */

void original_conversion_method_cleaned_up(block, matrix)
     Block *block;
     Matrix *matrix;
{
  double *seq_weight;		/* the contribution of this sequence to the */
				/* block.  = 1/(num seq in cluster) */

  /* get the weights of the sequences */
  seq_weight = clustered_weights(block);

  basic_matrix_construction(block, seq_weight, matrix);

  free(seq_weight);
}



/*
 * pre_weighted_conversion_method
 *   This conversion method uses the pre-weighted sequence scores and
 *   uses the basic matrix construction method.  This method has the
 *   following properties:
 *     X, '-', '*', & non-code scores are read straight from the frequencies
 *     the frequencies for B and Z are ignored, when a B or Z is encountered
 *       it is partitioned between D & N or E & Q.
 *     the matrix scores for B and Z are computed from the matrix scores of
 *       D & N and E & Q.
 *   Parameters:
 *     Block *block:   the block to be converted
 *     Matrix *matrix: where the resulting matrix will be put
 *   Error codes: none
 */

void pre_weighted_conversion_method(block, matrix)
     Block *block;
     Matrix *matrix;
{
  double *seq_weight;		/* the contribution of this sequence to the */
				/* block.  = the value given in the block */

  /* get the weights of the sequences */
  seq_weight = pre_weighted_sequences(block);

  if (seq_weight != NULL) {
    basic_matrix_construction(block, seq_weight, matrix);

    free(seq_weight);
  }
  else {
    /* remember that seq_weight has been freed already if NULL is returned */
    sprintf(ErrorBuffer,
	    "All weights in the block were less than or equal to zero.");
    ErrorReport(WARNING_ERR_LVL);
    sprintf(ErrorBuffer,
	    "Computing position-based sequence weights.\n");
    ErrorReport(WARNING_ERR_LVL);

    pb_weights(block);
    seq_weight = pre_weighted_sequences(block);

    if (seq_weight != NULL) {
       basic_matrix_construction(block, seq_weight, matrix);
       free(seq_weight);
    }


  }

} /* end of pre_weighted_convertion_method */




/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

/*
 * altschul data dependent conversion method aux. structures/defines
 */

/* these are from JGH's makealts, I am not sure where these are from exctly. */
/* I think they are the actual values from the aabet.h structures. */
/* AAS and AASALL are not matched up with any of the aabet.h defines. */
#define AASALL 26
/* 0- 1A 2R 3N 4D 5C 6Q 7E 8G 9H 10I 11L 12K 13M 14F 15P 16S 17T 18W 19Y
   20V 21B 22Z 23X 24* 25J,O,U */

#define MAXWIDTH 400			/* Max. block width */

struct working {	/* Working information for one column */
  double cnt[AASALL];		/* Sequence-weighted counts */
  double totcnt;
  double raw[AASALL];		/* unweighted counts */
  double totraw;
  double reg[AASALL];		/*  pseudo counts */
  double totreg;
  double probn[MAXDIRI];        /* for dirichlet */
  double probj[MAXDIRI];        /* for dirichlet */

};

struct work_pssm {              /* Working PSSM area */
  double** value;   /* double value[MAXWIDTH][AASALL]*/
  double* sum ;   /* double sum[MAXWIDTH] */
};

struct rank_cell {
        int aa;
        double value;
};

/*
 * altschul data dependent conversion method aux. functions
 */

/*==================================================================*/
static struct working *make_col()
{
  struct working *col;
  int aa;

  CheckMem(
	   col = (struct working *) malloc(sizeof(struct working))
	   );

  col->totcnt = col->totreg = 0.0;
  for (aa=0; aa < AASALL; aa++) {
    col->cnt[aa] = col->reg[aa] = 0.0;
  }

  return col;
}  /* end of make_col */


/*=====================================================================*/
static struct work_pssm *make_work_pssm(length)
int length;
{
  struct work_pssm *pssm;
  int pos, aa;
  double* double_pointer;

  CheckMem(
           pssm = (struct work_pssm *) malloc(sizeof(struct work_pssm))
           );

        double_pointer = (double*) calloc (length * AASALL, sizeof (double));
        pssm->value = (double **) calloc (length, sizeof (double*));
        pssm->sum = (double *) calloc (length, sizeof (double));

  for (pos = 0; pos < length; pos++) {
        pssm->value[pos] = &(double_pointer[pos * AASALL]);
  }

  for (pos = 0; pos < length; pos++) {
    pssm->sum[pos] = 0.0;
    for (aa=0; aa < AASALL; aa++) {
      pssm->value[pos][aa] = 0.0;
    }
  }

  return pssm;
}  /* end of make_work_pssm */
/*=====================================================================*/
void free_work_pssm(pssm)
struct work_pssm *pssm;
{
  free(pssm->value[0]);
  free(pssm->value);
  free(pssm->sum);
  free(pssm);
}  /* end of free_work_pssm */


/*======================================================================*/
static void counts(block, col, pos)
     Block *block;
     struct working *col;
     int pos;
{
  int seq, aa, aa1;

  col->totcnt = col->totraw = col->totreg = 0.0;
  for (aa = 0; aa < AASALL; aa++) {
    col->cnt[aa] = col->raw[aa] = col->reg[aa] = 0.0;
  }

  /*  Only count the real 20 aas, combine B(21) with D & Z(22) with E  */
  for (seq = 0; seq < block->num_sequences; seq++) {
    aa = block->residues[seq][pos];
    if (aa == 21) aa = 4;			/* combine B with D */
    if (aa == 22) aa = 7;			/* combine Z with E */
    if (aa >= 1 && aa < AAS) {
      col->cnt[aa] += block->sequences[seq].weight;
      col->totcnt += block->sequences[seq].weight;
      col->raw[aa] += 1.0;
      col->totraw += 1.0;
    }
    else {
/*
      sprintf(ErrorBuffer,
	      "Uncounted \"residue\" character for %s: %c\n",
	      block->number, aa_btoa[block->residues[seq][pos]]);
      ErrorReport(INFO_ERR_LVL);
*/
      /* If not one of the basic aas, divide the count among them */
      for (aa1 = 1; aa1 < AAS; aa1++)
      {
         col->cnt[aa1] += (block->sequences[seq].weight / 20.0);
         col->raw[aa1] += (1.0 / 20.0) ;
      }
      col->totcnt += block->sequences[seq].weight;
      col->totraw += 1.0;
    }
  }
}  /* end of counts */


/*=====================================================================*/
static int count_residues(col)
     struct working *col;
{
  int aa, nr;

  nr = 0;
  for (aa = 1; aa < AAS; aa++) {
    if (col->cnt[aa] > 0.0) nr++;
  }

  return nr;
}  /*  end of count_residues */


/*=======================================================================*/
static void pseudo_alts(col, qij, epsilon)
     struct working *col;
     struct float_qij *qij;
     double epsilon;
{
  int aa, row;

  /*---------- get the pseudo counts -------------------------------*/
  for (aa=1; aa < AAS; aa++) {
    col->reg[aa] = 0.0;
    for (row = 1; row < AAS; row++) {
      col->reg[aa] += (col->cnt[row] * qij->value[aa][row] / qij->marg[row]);
    }
    col->reg[aa] *= epsilon;
    if (col->totcnt > 0.0) col->reg[aa] /= col->totcnt;
    col->totreg += col->reg[aa];
  }
}  /* end of pseudo_alts */


/*========================================================================
	Computes scores for B, Z, X, -(gap) and *(stop) in a column
	of a PSSM using other scores in the column
	frequency[] is a global array created by load_frequencies()
	why is it an argument here?
==========================================================================*/
/* SOME DUPLICATE CODE FROM basic_matrix_construction(). */
static void compute_BZX(frequency, matrix, col)
     double *frequency;
     Matrix *matrix;
     int col;
{
  int aa;
  double dmean, dmin;
  double part_D;		/* the partition of D for B. */
				/* = freq[D] / ( freq[D] + freq[N] ) */
  double part_N;		/* the partition of N for B. */
				/* = freq[N] / ( freq[D] + freq[N] ) */
  double part_E;		/* the partition of E for Z. */
				/* = freq[E] / ( freq[E] + freq[Q] ) */
  double part_Q;		/* the partition of Q for Z. */
				/* = freq[Q] / ( freq[E] + freq[Q] ) */
  /*
   * find the partitions of D, N, E, and Q for B and Z
   */
  part_D = frequency[aa_atob['D']] /
    ( frequency[aa_atob['D']] + frequency[aa_atob['N']] );
  part_N = frequency[aa_atob['N']] /
    ( frequency[aa_atob['D']] + frequency[aa_atob['N']] );
  part_E = frequency[aa_atob['E']] /
    ( frequency[aa_atob['E']] + frequency[aa_atob['Q']] );
  part_Q = frequency[aa_atob['Q']] /
    ( frequency[aa_atob['E']] + frequency[aa_atob['Q']] );

    /* fill in the matrix for B, Z, X, gap, stop and non */
    matrix->weights[aa_atob['B']][col] =
      (part_D * matrix->weights[aa_atob['D']][col] +
       part_N * matrix->weights[aa_atob['N']][col]);
    matrix->weights[aa_atob['Z']][col] =
      (part_E * matrix->weights[aa_atob['E']][col] +
       part_Q * matrix->weights[aa_atob['Q']][col]);

    /*   X or unk gets the weighted average score; - and * get the min score */
    dmin = 999.99;
    dmean = 0.0;
    for (aa=1; aa<20; aa++)
    {
       dmean += frequency[aa] * matrix->weights[aa][col];
       if (matrix->weights[aa][col] < dmin) dmin = matrix->weights[aa][col];
    }
    matrix->weights[aa_atob['X']][col] = dmean;
    matrix->weights[aa_atob[25]][col] = dmean;		/* unknown res */
    matrix->weights[aa_atob['-']][col] = dmin;
    if (dmin > 0.0) matrix->weights[aa_atob['*']][col] = 0.0;
    else            matrix->weights[aa_atob['*']][col] = dmin;

}  /* end of compute_BZX */


/*=========================================================================
      Adds negative minval to give all positive matrix,
      then multiplies by 99/maxval to give scores ranging from 0 to 99
      NOTE: Not 0 to 100 because "output_matrix" routine might not leave
	    enough space.
===========================================================================*/
static void positive_matrix(freqs, pssm, matrix)
     double *freqs;
     struct work_pssm *pssm;
     Matrix *matrix;
{
  int pos, aa;
  double factor, maxval, minval, dtemp;

  minval = 9999.9;
  maxval = -9999.9;
  for (pos = 0; pos < matrix->width; pos++) {
    for (aa=1; aa < AAS; aa++) {
      if (pssm->value[pos][aa] < minval) minval = pssm->value[pos][aa];
      if (pssm->value[pos][aa] > maxval) maxval = pssm->value[pos][aa];
    }
  }

  if (minval < 0.0) {
    factor = 99.0 / (maxval - minval);
  }
  else {
    factor = 99.0 / maxval;
  }
  if (factor < 1.0) {
    factor = 1.0;
  }
  for (pos = 0; pos < matrix->width; pos++) {
    for (aa=1; aa < AAS; aa++) {
      if (minval < 0.0) {
	dtemp = factor * (pssm->value[pos][aa] - minval);
      }
      else {
	dtemp = factor * pssm->value[pos][aa];
      }
      matrix->weights[aa][pos] = (MatType) dtemp;
    }
    compute_BZX(freqs, matrix, pos);
  }  /*  end of for pos */
}  /* end of positive_matrix */


/*==========================================================================
     Uses Altschul's method of getting pseudo-counts with a qij matrix,
	scale=0 => DEFAULT nats, positive integers
	scale=1 => bits, signed integers
	scale=2 => half bits, signed integers
	scale=3 => third bits, signed integers
	scale=9 => few pseudo counts, bits, positive integers
	scale=10 => odds ratios
	scale=20 => counts + pseudo counts
	scale=21 => counts
===========================================================================*/
static void make_alts(block, matrix, freqs, qij, RTot, scale)
     Block *block;
     Matrix *matrix;
     double *freqs;
     struct float_qij *qij;
     double RTot;		/* Total R for Altschul */
     int scale;
{
  double factor, dtemp, epsilon;
  int pos, aa;
  struct working *col;
  struct work_pssm *pssm;

  factor = 1.0;
  if (scale > 0 && scale < 9) factor = (double) scale / log(2.0);
  col = make_col();
  pssm = make_work_pssm(block->width);

  /*--------------  Do one position at a time -------------------*/
  for (pos = 0; pos < block->width; pos++)
  {
    /*-------- count the number of each aa in this position ------------*/
    counts(block, col, pos);

    /*-------- determine total number of pseudo-counts in column ------*/
    /*  9 => compute minimal number of pseudo-counts; RTot should
	be a large number in this case, such as 1000  */
    if (scale == 9)
/*
    {  epsilon = (double) col->totcnt / RTot ; }
*/
    {  epsilon = (double) 1.0 / RTot ; }
    else
    {  epsilon = (double) RTot * count_residues(col);  }

    /*---------- get the pseudo counts -------------------------------*/
    pseudo_alts(col, qij, epsilon);

    /*---------   Fill in the matrix entries --------------------*/
    pssm->sum[pos] = 0.0;
    for (aa=1; aa < AAS; aa++)
    {
       /*    Count proportions; returned if scale == 21   */
       if (scale == 21)
       {
          pssm->value[pos][aa] = col->cnt[aa];
          if ( col->totcnt > 0.0)
             pssm->value[pos][aa] /= col->totcnt;
       }

       /*    Count+pseudo proportions; returned if scale == 20   */
       if (scale <= 20)
       {
          pssm->value[pos][aa] = col->cnt[aa] + col->reg[aa];
          if ( (col->totcnt + col->totreg) > 0.0)
             pssm->value[pos][aa] /= (col->totcnt + col->totreg);
       }

       /*     Odds ratios; returned if scale == 10   */
       if (scale < 20 && freqs[aa] > 0.0)
          pssm->value[pos][aa] /= freqs[aa];

       /*   take the log of the odds ratio  */
       if (scale < 10 && pssm->value[pos][aa] > 0.0)
          pssm->value[pos][aa] = log(pssm->value[pos][aa]);

       pssm->sum[pos] += pssm->value[pos][aa];

       /*  scale the matrix  */
       dtemp = factor * pssm->value[pos][aa];
       matrix->weights[aa][pos] = dtemp;

     }  /* end of aa */

     compute_BZX(freqs, matrix, pos);
  }  /*  end of for pos */

  /*------ Now make the final scores; make log scores non-neg    */
  /*       by subtracting the min. value for all positions ----- */
  if (scale==0 || scale==9) positive_matrix(freqs, pssm, matrix);

  free(col);
  free_work_pssm(pssm);
}  /* end of make_alts */

/*==========================================================================
     Uses Gribskov's method
     Assumes a non-negative substition matrix (min. value = 0)
     Sets 20 aas plus B and Z from the matrix; X, etc to zero = min. value
     Index for B is 21, for Z is 22
===========================================================================*/
void make_gribs(block, matrix, subst)
Block *block;
Matrix *matrix;
struct float_qij *subst;
{
   double sum, temp[AASALL];
   int pos, aa, aa1;
   struct working *col;

   col = make_col();

   for (pos = 0; pos < block->width; pos++)
   {
      /*-------- count the number of each aa in this position ------------*/
      counts(block, col, pos);

      sum = 0.0;
      for (aa=0; aa <= 22; aa++)
      {
         temp[aa] = 0.0;
         for (aa1=1; aa1 <= 22 ; aa1++)
         {
            temp[aa] += col->cnt[aa1] * subst->value[aa][aa1];
         }
         temp[aa] /= col->totcnt;
         sum += temp[aa];
      }
      for (aa=1; aa <= 22; aa++)
      {
/*		ends up rounding to zero!
         dtemp = temp[aa];
         matrix->weights[aa][pos] = round(dtemp);
*/
         matrix->weights[aa][pos] = temp[aa];
      }
      /*  Set X, * and - to zero */
      matrix->weights[0][pos] = matrix->weights[23][pos] = 0.0;
      matrix->weights[24][pos] = matrix->weights[25][pos] = 0.0;
   }  /*  end of for pos */

   free(col);
}   /*  end of make_gribs */

/*======================================================================*/
struct float_qij *load_qij( FILE *fin)
{
  char line[LARGE_BUFF_LENGTH], *ptr;
  double total;
  int alpha[AAS], nrows, ncols, row, col, i;
  struct float_qij *new;

  CheckMem (
	    new = (struct float_qij *) malloc(sizeof(struct float_qij))
	    );

  /*----------Read file until first non-blank line --------------*/
  /* Skip comments at beginning of file - 1st char = #, > or ;   */
  line[0] = '\0';
  while (((int) strlen(line) < 1 ||
           line[0]=='#' || line[0]=='>' || line[0]==';')
	 && fgets(line, sizeof(line), fin) != NULL)
    ;

  /*------See if the first line has characters on it ------------*/
  for (col=0; col < AAS; col++) alpha[col] = -1;
  if (strstr(line, "A") != NULL) {	/* This line has characters */
    row = 0;	/* # of alphabetic characters on the line */
    for (i=0; i< (int) strlen(line); i++) {
      col = -1;
      col = aa_atob[(int)line[i]];
      if (col >= 0 && col < AAS) {
	alpha[row] = col;
	row++;
      }
      else if (isalpha(line[i])) {
	row++;
      }
    }
  }

  /*-------Get the data values now ------------*/
  for (row=0; row<AAS; row++) {
    for (col=0; col<AAS; col++) {
      new->value[row][col] = -1.0;		/* Null value */
    }
  }

  nrows = 0;
  line[0] = '\0';
  while (fgets(line, sizeof(line), fin) != NULL) {
    if ((int) strlen(line) > 1 && nrows < AAS) {
      if (alpha[nrows] >= 0 && alpha[nrows] < AAS) {
	row = alpha[nrows]; ncols = 0;
	ptr = strtok(line, " ,\n");
	while (ptr != NULL) {
	  if (strspn(ptr, ".+-0123456789") == strlen(ptr)) {
	    col = alpha[ncols];
	    if (col >= 0 && col < AAS) {
	      new->value[row][col] = (double) atof(ptr);
	    }
	    ncols++;
	  }
	  ptr = strtok(NULL, " ,\n");
	}
      }
      nrows++;
    }
  }

  /*-------If some entries are still missing, assume symmetry ---------*/
  for (row=0; row<AAS; row++) {
    for (col=0; col<AAS; col++) {
      if (new->value[row][col] < 0.0) {
	new->value[row][col] = new->value[col][row];
      }
    }
  }

  /*-------compute the marginal probabilities ---------*/
  total = 0.0;
  for (row=1; row<AAS; row++) {
    new->marg[row] = 0.0;
    for (col=1; col<AAS; col++) {
      new->marg[row] += new->value[row][col];
    }
    total += new->marg[row];
  }

  return new;
}  /* end of load_qij */


/*=========================================================================
    Normalize sequence weights to add up to the number of sequences
=========================================================================*/
void normalize(block)
     Block *block;
{
  double sum, factor;
  int seq;

  sum = 0.0;
  for (seq = 0; seq < block->num_sequences; seq++) {
    sum += block->sequences[seq].weight;
  }
  if (sum <= 0.0)   /*  All seqs have zero weight! Add pb weights */
  {
      sprintf(ErrorBuffer,
	    "All weights in the block were less than or equal to zero.");
      ErrorReport(WARNING_ERR_LVL);
      sprintf(ErrorBuffer,
	    "Computing position-based sequence weights.\n");
      ErrorReport(WARNING_ERR_LVL);

      pb_weights(block);
      for (seq = 0; seq < block->num_sequences; seq++) {
        sum += block->sequences[seq].weight;
      }
  }
  if (sum > 0.0) factor = block->num_sequences / sum;
  else factor = 1.0;
  for (seq = 0; seq < block->num_sequences; seq++) {
    block->sequences[seq].weight *= factor;
  }
}  /*  end of normalize */


/*
 * altschul_data_dependent_conversion_method
 *   Expects global variables frequency[], Qij[][] and RTot
 *   Parameters:
 *     Block *block:   the block to be converted
 *     Matrix *matrix: where the resulting matrix will be put
	scale = 0	log odds ratios, positive values
		1	bits, log odds ratios
		2	half bits, log odds ratios
         	3	third bits, log odds ratios
		9	few pseudo counts
		10	odds ratios
		20	(counts+pseudo counts)/tot
		21	counts
 *   Error codes: Checks MAXWIDTH, global variables
 */
void altschul_data_dependent_conversion_method(block, matrix, scale)
     Block *block;
     Matrix *matrix;
     int scale;
{
  int itemp;

  if (block->width > MAXWIDTH)
  {
      itemp = MAXWIDTH;
      sprintf(ErrorBuffer, "convert: Block is too wide, unable to continue (max=%d).\n",
		itemp);
      ErrorReport(FATAL_ERR_LVL);
  }
  if (Qij == NULL) {
      sprintf(ErrorBuffer, "Qij matrix missing, unable to continue.\n");
      ErrorReport(FATAL_ERR_LVL);
    }

  /* check to see if the block has sequence weights */

  normalize(block);   /*  Make weights sum to number of sequences */
  make_alts(block, matrix, frequency, Qij, RTot, scale);

} /* end of altschul_data_dependent_conversion_method  */

/*
 * gribskov_dependent_conversion_method
 *   Expects global variables frequency[]
 *   Loads default.sij
 *   Parameters:
 *     Block *block:   the block to be converted
 *     Matrix *matrix: where the resulting matrix will be put
 *   Error codes: none
 */

void
gribskov_conversion_method(block, matrix)
     Block *block;
     Matrix *matrix;
{
  char sijname[SMALL_BUFF_LENGTH], *blimps_dir;
  struct float_qij *sij_matrix;
  FILE *fp=NULL;
  /* Is this an error?  fp is never opened */

  /*  load the substitution matrix */
  blimps_dir = getenv("BLIMPS_DIR");
  if (blimps_dir != NULL)
  {
     sprintf(sijname, "%s/docs/default.sij", blimps_dir);
  }
  else
  {
     sprintf(sijname, "default.sij");
  }

  sij_matrix = load_qij(fp);
  fclose(fp);
  if (sij_matrix == NULL)
  {
      sprintf(ErrorBuffer,
	      "gribskov_conversion_method: default.sij matrix missing, Cannot continue.\n");
      ErrorReport(FATAL_ERR_LVL);
  }

  /* check to see if the block has sequence weights */
  normalize(block);   /*  Make weights sum to number of sequences */
  make_gribs(block, matrix, sij_matrix);

}

/* end of gribskov_dependent_conversion_method
 * Sequence weighting methods.  Various methods of giving different
 * weights to the sequences.
 */

/*
 * static double *clustered_weights(block)
 * static double *pre_weighted_sequences(block)
 * void pb_weights(block)
 */


/*
 * clustered_weights
 *   This sequence weighing depends on the clustering of the block.
 *   The weight of a sequence is one over the number of sequences in
 *   that cluster.
 *   Parameters:
 *     Block* block: the block that the sequences are from.
 *   Return codes: a pointer to the sequence weights array.
 *   Error codes:
 */

static double *clustered_weights(block)
     Block* block;
{
  double *seq_weight;		/* the contribution of this sequence to the */
				/* block.  = 1/(num seq in cluster) */

  int clust;			/* a cluster counter */
  int seq;			/* a sequence counter */

  int seq_num;			/* keep's track of where to start when */
				/* filling seq_weight*/

  double weight;		/* a temporary holder for the weight */


  /*
   * create and fill the seq_weight array.
   */
  /* allocate space */
  CheckMem(
	   seq_weight = (double *) calloc(block->num_sequences, sizeof(double))
	   );

  /* fill the seq_weight array */
  seq_num = 0;
  for (clust=0; clust<block->num_clusters; clust++) {
    weight = 1.0/(double)block->clusters[clust].num_sequences;
    /* fill each sequence in the cluster */
    for (seq=seq_num; seq<seq_num+block->clusters[clust].num_sequences;
	 seq++) {
      seq_weight[seq] = weight;
    }
    seq_num += block->clusters[clust].num_sequences;
  }

  return seq_weight;
} /* end of clustered_weights */

/*
 * pre_weighted_sequences
 *   This sequence weighing was done beforehand.  The weight of the
 *   sequence was read from the block.
 *   Parameters:
 *     Block* block: the block that the sequences are from.
 *   Return codes: a pointer to the sequence weights array.
 *   Error codes: NULL if all of the weights in the block are <= zero.
 */

static double *pre_weighted_sequences(block)
     Block* block;
{
  double *seq_weight;		/* the contribution of this sequence to the */
				/* block. */

  int seq;			/* a sequence counter */

  int num_seqs;			/* keep's track of total sequences */

  double max_weight;

  /*
   * create and fill the seq_weight array.
   */
  /* allocate space */
  CheckMem(
	   seq_weight = (double *) calloc(block->num_sequences, sizeof(double))
	   );

  num_seqs = block->num_sequences;

  /* fill each sequence weight from the block */
  max_weight = 0.0;
  for (seq=0; seq<num_seqs; seq++) {
    seq_weight[seq] = block->sequences[seq].weight;
    if (seq_weight[seq] > max_weight) {
      max_weight = seq_weight[seq];
    }
  }

  if (max_weight <= 0) {
    /* all weights were zero */
    free(seq_weight);

    return NULL;  /* the error is handled in the calling function */
  }

  return seq_weight;
}  /* end of pre_weighted_sequences  */

/*=======================================================================
      Compute Steve's position-based sequence weights
	NOTE:  Characters - (0), X (23) and * (24) are ignored
	       Characters 1-22 are the 20 basic aas plus B & Z
	MATRIX_AA_WIDTH = 26, defined in matrix.h
========================================================================*/
void pb_weights(block)
Block *block;
{
   struct pb_counts *pb;
   //double factor, dtemp;
   double dtemp;
   int seq, pos, aa, width;

   width = block->width;
   CheckMem(
	pb = (struct pb_counts *) malloc(width*sizeof(struct pb_counts))
    );

   for (pos = 0; pos < width; pos++)
   {
      pb[pos].diffaas = 0.0;
      for (aa = 0; aa < MATRIX_AA_WIDTH; aa++)
         pb[pos].naas[aa] = (double) 0.0;
   }

   for (pos = 0; pos < width; pos++)
      for (seq = 0; seq < block->num_sequences; seq++)
      {
         if (block->residues[seq][pos] >= 1 &&
             block->residues[seq][pos] <= 22)
         {
            pb[pos].naas[block->residues[seq][pos]] += 1;
         }
         else
         {
            sprintf(ErrorBuffer, "pb_weights:%d ignored for %s\n",
                 block->residues[seq][pos], block->number);
	    ErrorReport(INFO_ERR_LVL);
         }
      }

   // factor = 1.0;
   for (pos = 0; pos < width; pos++)
   {
      for (aa = 1; aa <= 22; aa++)
      {
         if (pb[pos].naas[aa] > 0.0)
         {
            pb[pos].diffaas += 1;	/* # of different types of aas in pos */
         }
      }
   }

   for (seq = 0; seq < block->num_sequences; seq++)
   {
      block->sequences[seq].weight = 0.0;
      for (pos = 0; pos < width; pos++)
      {
         aa = block->residues[seq][pos];
         dtemp = pb[pos].diffaas * pb[pos].naas[aa];
         if (dtemp > 0.0)
            block->sequences[seq].weight += 1.0 / dtemp;
      }
   }

   free(pb);
}  /* end of pb_weights */




/*
 * Matrix construction methods.  Build matricies from blocks and
 * sequence weights.
 */

/*
 * static void basic_matrix_construction(block, seq_weight, matrix)
 *
 */


/*
 * basic_matrix_construction
 *   This is the most general matrix construction method.  Each
 *   occurence of a residue in a sequence contributes the sequence
 *   weight of the sequence divided by the frequency of the residue to
 *   the total weight of the residue in the matrix at that column.
 *   This contribution is scaled by the total of the other residue
 *   contributions for the column of the matrix.  The following are
 *   exceptions to the method:
 *     X, '-', '*', & non-code scores are read straight from the frequencies
 *     the frequencies for B and Z are ignored, when a B or Z is encountered
 *       it is partitioned between D & N or E & Q.
 *     the matrix scores for B and Z are computed from the matrix scores of
 *       D & N and E & Q.
 *   Parameters:
 *     Block *block:       the block to be converted
 *     double *seq_weight: the weights of each sequence.
 *     Matrix *matrix:     where the resulting matrix will be put
 *   Error codes: none
 */

static void basic_matrix_construction(block, seq_weight, matrix)
     Block *block;
     double *seq_weight;
     Matrix *matrix;
{
  int seq;			/* a sequence counter */
  int col;			/* a column counter */
  int aa;			/* an amino acid counter */

  int num_sequences;		/* the number of sequences in the block */

  Residue res;			/* a residue.  Used for keeping hold of */
				/* block->residues[seq][col] */

  double total;			/* the sum of the seq_weight/freq[] for each */
				/* column (every amino acid).  This is used */
				/* to scale each matrix entry for a column */
				/* to a percentage */

  double count[MATRIX_AA_WIDTH]; /* the sum of the weights for each amino */
				 /* acid in the column */

  double part_D;		/* the partition of D for B. */
				/* = freq[D] / ( freq[D] + freq[N] ) */
  double part_N;		/* the partition of N for B. */
				/* = freq[N] / ( freq[D] + freq[N] ) */
  double part_E;		/* the partition of E for Z. */
				/* = freq[E] / ( freq[E] + freq[Q] ) */
  double part_Q;		/* the partition of Q for Z. */
				/* = freq[Q] / ( freq[E] + freq[Q] ) */


  /*
   * find the partitions of D, N, E, and Q for B and Z
   */
  part_D = frequency[aa_atob['D']] /
    ( frequency[aa_atob['D']] + frequency[aa_atob['N']] );
  part_N = frequency[aa_atob['N']] /
    ( frequency[aa_atob['D']] + frequency[aa_atob['N']] );
  part_E = frequency[aa_atob['E']] /
    ( frequency[aa_atob['E']] + frequency[aa_atob['Q']] );
  part_Q = frequency[aa_atob['Q']] /
    ( frequency[aa_atob['E']] + frequency[aa_atob['Q']] );


  /*
   * make the matrix
   */
  /* for every column in the block, count up the countribution of each aa, */
  /* and fill the matrix column */
  for (col=0; col<block->width; col++) {
    /* reset the counts array and the total */
    for (aa=0; aa<MATRIX_AA_WIDTH; aa++) {
      count[aa] = 0.0;
    }
    total = 0.0;

    /* for every sequence in the block build up the count and total data */
    num_sequences = block->num_sequences;
    for (seq=0; seq<num_sequences; seq++) {
      res = block->residues[seq][col];
      /* if not B, Z, X, gap, stop and non, do the regular */
      if ((res != aa_atob['B']) &&
	  (res != aa_atob['Z']) &&
	  (res != aa_atob['X']) &&
	  (res != aa_atob['-']) && /* gap */
	  (res != aa_atob['*']) && /* stop */
	  (res != AAID_NAR)) {     /* non */
	if (frequency[res] != 0.0) { /* try to protect from div by zero */
	  count[res] += seq_weight[seq] / frequency[res];
	  total      += seq_weight[seq] / frequency[res];
	}
      }
      /* otherwise, if it is B, partition B between D and N */
      else if (res == aa_atob['B']) {
	if (frequency[aa_atob['D']] != 0.0) { /* try to protect div by zero */
	  count[aa_atob['D']] +=
	    (part_D * seq_weight[seq]) / frequency[aa_atob['D']];
	  total +=
	    (part_D * seq_weight[seq]) / frequency[aa_atob['D']];
	}
	if (frequency[aa_atob['N']] != 0.0) { /* try to protect div by zero */
	  count[aa_atob['N']] +=
	    (part_N * seq_weight[seq]) / frequency[aa_atob['N']];
	  total +=
	    (part_N * seq_weight[seq]) / frequency[aa_atob['N']];
	}
      }
      /* otherwise, if it is Z, partition Z between E and Q */
      else if (res == aa_atob['Z']) {
	if (frequency[aa_atob['E']] != 0.0) { /* try to protect div by zero */
	  count[aa_atob['E']] +=
	    (part_E * seq_weight[seq]) / frequency[aa_atob['E']];
	  total +=
	    (part_E * seq_weight[seq]) / frequency[aa_atob['E']];
	}
	if (frequency[aa_atob['Q']] != 0.0) { /* try to protect div by zero */
	  count[aa_atob['Q']] +=
	    (part_Q * seq_weight[seq]) / frequency[aa_atob['Q']];
	  total +=
	    (part_Q * seq_weight[seq]) / frequency[aa_atob['Q']];
	}
      }
      /* otherwise, if it is X, gap, stop or non, don't calculate */
    }

    /* for every amino acid, fill in the matrix for this column, unless it */
    /* is B, Z, X, gap, stop or non */
    for (aa=0; aa<MATRIX_AA_WIDTH; aa++) {
      if ((aa != aa_atob['B']) &&
	  (aa != aa_atob['Z']) &&
	  (aa != aa_atob['X']) &&
	  (aa != aa_atob['-']) && /* gap */
	  (aa != aa_atob['*']) && /* stop */
	  (aa != AAID_NAR)) {	  /* non */
	if (total != 0.0) {	/* try to protect from div by zero */
	  matrix->weights[aa][col] = count[aa] * 100.0 / total;
	}
	else {
	  matrix->weights[aa][col] = 0.0;
	}
      }
    }

    /* fill in the matrix for B, Z, X, gap, stop and non */
    matrix->weights[aa_atob['B']][col] =
      (part_D * matrix->weights[aa_atob['D']][col] +
       part_N * matrix->weights[aa_atob['N']][col]);
    matrix->weights[aa_atob['Z']][col] =
      (part_E * matrix->weights[aa_atob['E']][col] +
       part_Q * matrix->weights[aa_atob['Q']][col]);
    matrix->weights[aa_atob['X']][col] = frequency[aa_atob['X']];
    matrix->weights[aa_atob['-']][col] = frequency[aa_atob['-']];
    matrix->weights[aa_atob['*']][col] = frequency[aa_atob['*']];
    matrix->weights[AAID_NAR][col]     = frequency[AAID_NAR];

  } /* end for each column */

} /* end of basic_matrix_construction */









/*
 * The original PATMAT method.  Hardly ever done anymore.
 */

/*
 * original_conversion_method
 *   The original conversion method.  This is done by weighted average of the
 *   clusters.  This follows the method in patmat.
 *   Parameters:
 *     Block *block:   the block to be converted
 *     Matrix *matrix: where the resulting matrix will be put
 *   Error codes: none
 */

void original_conversion_method(block, matrix)
     Block *block;
     Matrix *matrix;
{
  int aa, pos, clust, seq;
  int num_clusters, num_sequences, width;
  Residue *sequence;
  double *divisor;

  double tmp_double;

  double *number_occurred[MATRIX_AA_WIDTH]; /* [amino acid][position] */



  /* allocate the space for the number occured 2-d array ([amino acid][pos]) */
  CheckMem(
	 number_occurred[0] = (double *) calloc(matrix->width*MATRIX_AA_WIDTH,
						sizeof(double))
	   );

  /* setup the array of pointers */
  for (aa=0; aa < MATRIX_AA_WIDTH; aa++) {
    number_occurred[aa] = number_occurred[0] + aa*matrix->width;
  }

  /* get the space for the divisor array */
  CheckMem(
	   divisor = (double *) calloc(block->width, sizeof(double))
	   );



  width = block->width;

  /*
   * get the number of residues at each position.  clusters use
   * fractional values.
   */

  /* for each cluster */
  num_clusters = block->num_clusters;
  for (clust=0; clust < num_clusters; clust++) {

    /* for each sequence in the cluster */
    num_sequences = block->clusters[clust].num_sequences;
    for (seq=0; seq < num_sequences; seq++) {

      /* for each position in the sequence */
      sequence = block->clusters[clust].sequences[seq].sequence;
      for (pos=0; pos < width; pos++) {

	/* increase the number of amino acids seen at this position */
	/* in the number_occured matrix */
	number_occurred[sequence[pos]][pos] +=
	  (double) 1.0/(double)num_sequences; /* type casting to be safe */

      } /* end for each position */
    } /* end for each sequence */
  } /* end for each cluster */

  /*
   * calculate the divisor for later use.
   * divisor[i] = summation over all amino acids of
   *                ( num_occured[aa][i] / frequency[aa] )
   */

  /* for each position */
  for (pos=0; pos < width; pos++) {

    /* for all amino acids */
    for (aa=0; aa < MATRIX_AA_WIDTH; aa++) {

      if (frequency[aa] > 0.0) {
	if (number_occurred[aa][pos] > 0.0) {
	  divisor[pos] += number_occurred[aa][pos] / frequency[aa];
	}
      }

    } /* end for all amino acids */
  } /* end for each position */


  /*
   * calculate the matrix weights.
   *                              ( num_occur[aa][i] / freq[aa] ) * 100
   * in general: weights[aa][i] = -------------------------------------
   *                                     divisor[i]
   */

  /* for each amino acid */
  for (aa=0; aa < MATRIX_AA_WIDTH; aa++) {

    /* if there is a valid weight (frequency is not -1.0) */
    if ((int)frequency[aa] != -1) { /* typecast to be safe */

      /* for each position */
      for (pos=0; pos < width; pos++) {

	if (frequency[aa] < -1.0) {
	  tmp_double = frequency[aa];
	}
	else if (number_occurred[aa][pos] > 0) {
	  tmp_double = ( (number_occurred[aa][pos] / frequency[aa]) * 100 )
	               / divisor[pos];
	}
	else {
	  tmp_double = 0.0;
	}

	matrix->weights[aa][pos] = tmp_double;

      } /* end for each position */
    } /* end if valid weight */
  } /* end for each amino acid */

  /* free temporarily allocated space */
  free(divisor);
  free(number_occurred[0]);

}  /* end of original_conversion_method  */

/*===================================================================*/
/*  Pauline's additions to convert.c for SIFT */

void SIFT_conversion_method(block, pssm)
     Block *block;
     Matrix *pssm;
{
     int diri_option, gap_option, exp_option, subtract_option;

	diri_option = gap_option = exp_option = TRUE;
	subtract_option = FALSE;

      if (Qij == NULL) {
          sprintf(ErrorBuffer, "Qij matrix missing, unable to continue.\n");
          ErrorReport(FATAL_ERR_LVL);
        }
        normalize  (block);
        SIFT_pssm(block, pssm, frequency, Qij, diri_option, gap_option,
                          exp_option, subtract_option);
}  /* end of SIFT_conversion_method */

/*==========================================================================
10/20/00 SIFT pseudocounts
13-dirichlet component
option:  diri_pseudocounts
	  TRUE: use 13-component Dirichlet
	  FALSE use BLOSUM62 qij's
option:  gap: TRUE -allow everything
	      FALSE - just look at 20 amino acids
option:  exp (m)
	 m=0 # of amino acids at pos (default)
	 m=1  similiarity scale (more conservative)

option:	subtract_threshold TRUE : scores -= SIFT_TOLERANCE
			   FALSE: leave as original scores
===========================================================================*/
void SIFT_pssm(block, matrix, freqs, qij,
		diri_pseudocounts, gap_option,
		exp_option, subtract_threshold)
Block *block;
Matrix *matrix;
double *freqs;
struct float_qij *qij;
int diri_pseudocounts, gap_option, exp_option, subtract_threshold;
{
  FILE *rfp;
  double dtemp, epsilon;
  int pos, aa,  itemp;
  struct working *col;
  struct work_pssm *pssm;
  int original_aa;
  struct diri* diric;
  int max_aa;
  double min_freq;
  int div_by_max;
  struct float_qij *rank_matrix;
  char *blimps_dir, diriname[SMALL_BUFF_LENGTH], rankname[SMALL_BUFF_LENGTH];

  blimps_dir = getenv("BLIMPS_DIR");
   if (blimps_dir != NULL)
   {
      sprintf(diriname, "%s/docs/default.diri", blimps_dir);
      sprintf(rankname, "%s/docs/default.rank", blimps_dir);
/*>>>>> Need to check that rankname & diriname actually exists
        for WWW servers would be better to check in current
        directory first  <<<<<<*/
   }
   else
   {
      sprintf(diriname, "default.diri");
      sprintf(rankname, "default.rank");
   }

  diric = load_diri (diriname);
  rank_matrix = NULL;

  if (exp_option == 1)
  {
     if ( (rfp = fopen(rankname, "r") ) == NULL)
     {
         sprintf(ErrorBuffer, "SIFT_pssm(): Cannot open %s\n", rankname);
         ErrorReport(FATAL_ERR_LVL);
     }
     else
     {
        rank_matrix = load_qij(rfp);
        fclose(rfp);
     }
  }

  div_by_max = TRUE;
  col = make_col();
  pssm = make_work_pssm(block->width);

  /*--------------  Do one position at a time -------------------*/
  for (pos = 0; pos < block->width; pos++)
  {
    /*-------- count the number of each aa in this position ------------*/
    if (gap_option) {
	counts (block,col, pos);
    } else {
	counts_nogaps(block, col, pos);
   }
    /*-------- determine total number of pseudo-counts in column ------*/
    epsilon = 0;
    itemp = count_residues (col);
    original_aa = block->residues[0][pos];

    if (exp_option == 1) {
	/*printf ("pos %d :  ", pos); */
	dtemp = similarity_dependent_scale(col, rank_matrix, original_aa);
	if (itemp == 1) {
		epsilon = 0;
	} else {
		epsilon = exp (dtemp);
	}
    } else {
	if (itemp == 1) { epsilon = 0; }
        else if (itemp > 7) { epsilon = 1000; }
        else { epsilon = exp( (double) itemp) ; }

   }

    /*---------- get the pseudo counts -------------------------------*/
    if (diri_pseudocounts) {
	pseudo_diri (col, diric, epsilon );
    } else {
	pseudo_alts(col, qij, epsilon);
    }

    /*---------   Fill in the matrix entries --------------------*/
    pssm->sum[pos] = 0.0;
    for (aa=1; aa < AAS; aa++)
    {

	pssm->value[pos][aa] = col->cnt[aa] + epsilon * col->reg[aa];
        if ( (col->totcnt + col->totreg) > 0.0)
             pssm->value[pos][aa] /= (col->totcnt + epsilon);


    } /* end of aa, pssm->values filled for a given pos */

   if (div_by_max) {
        max_aa = find_max_aa_pssm (pssm, pos);
	min_freq = pssm->value[pos][max_aa];
       for (aa = 1; aa < AAS; aa++) {
		pssm->value[pos][aa] /= min_freq;
	}
   }

   if (subtract_threshold) {
	   for (aa = 1; aa < AAS; aa++) {
	        pssm->value[pos][aa] -= SIFT_TOLERANCE;
   	}
   }

       pssm->sum[pos] += pssm->value[pos][aa];

    	for (aa = 1; aa < AAS; aa++) {
	  matrix->weights[aa][pos] = pssm->value[pos][aa];
	}
	original_aa = block->residues[0][pos];
        if ( (matrix->weights[original_aa][pos] < SIFT_TOLERANCE
                         && (!subtract_threshold)) ||
              (matrix->weights[original_aa][pos] < 0.0 && subtract_threshold) )
        {
	   sprintf (ErrorBuffer, "SIFT_pssm(): Amino acid %c at pos %d in your original sequence was not allowed by the prediction.\n",
		aa_btoa[original_aa], pos);
           ErrorReport(WARNING_ERR_LVL);
	}

	} /* end of for pos */

  if (diri_pseudocounts) {
	free (diric);
  }

  free(col);
  free_work_pssm(pssm);

}  /* end of SIFT_pssm */

/*======================================================================
   Dirichlet input file order (0-19) is  ARNDCQEGHILKMFPSTWYV
   Blimps order (0-24) is               -ARNDCQEGHILKMFPSTWYVBZX*
   Store Dirichlet alphas in positions 1-20 to match blimps
=======================================================================*/
struct diri *load_diri (filename)
char filename[LARGE_BUFF_LENGTH];
{
   FILE* fin;
   int i, aa, numc, type;
   char line[MAXLINE], *ptr;
   struct diri *diric;
   double denom;
   double background_frequency[AAS];

   if ((fin = fopen (filename, "r")) == NULL)
   {
      sprintf (ErrorBuffer, "dirichlet(): Cannot open dirichlet file %s\n",
			 filename);
      ErrorReport(FATAL_ERR_LVL);
   }

   diric = (struct diri *) malloc(sizeof(struct diri));
   if (diric == NULL)
   {
      sprintf (ErrorBuffer, "dirichlet(): OUT OF MEMORY\n");
      ErrorReport(FATAL_ERR_LVL);
   }
   numc = 0;
   while (numc < MAXDIRI && !feof(fin) && fgets(line, MAXLINE, fin) != NULL)
   {
      type = 0;
      if (strstr(line, "Mixture=") != NULL) type = 1;
      if (strstr(line, "Alpha=") != NULL)  type = 2;
      if (type > 0)
      {
         ptr = strtok(line, "="); ptr = strtok(NULL, "\t\n\r ");
         switch (type)
         {
            case 1:
               diric->q[numc] = atof(ptr);
               break;
            case 2:
               diric->altot[numc] = atof(ptr);
               aa = 1;
               while (aa < AAS && ptr != NULL)
               {
                  ptr = strtok(NULL, "\t\n\r ");
                  diric->alpha[numc][aa++] = atof(ptr);
               }
               numc++;
               break;
            default:
               break;
         }
      }
   }  /* end of while */
   diric->ncomp = numc;
	for (i = 0; i < numc; i++) {
		for (aa = 1; aa < AAS; aa++) {
			diric->alpha_normalized[i][aa] = diric->alpha[i][aa]/diric->altot[i];
		}
	}

	for (aa = 1; aa < AAS; aa++) {
		denom = 0.0;
		for (i = 0; i < numc; i++) {
			denom += (diric->q[i] * diric->alpha[i][aa] /
							diric->altot[i]);
		}
		background_frequency[aa] = denom;
	}
	for (i =0; i < numc; i++) {
	/*	printf ("Component %d ratio of aa relative to background \n", i); */
		for (aa = 1; aa < AAS; aa++) {
			diric->frequency_to_background_ratio[i][aa] =
	diric->alpha[i][aa]/(diric->altot[i] * background_frequency[aa]);
		}
	}
   fclose(fin);
   return(diric);
}   /* end of load_diri */

/*========================================================================
==========================================================================*/
/* only count valid amino acids, no considerationto gaps */
void counts_nogaps (block, col, pos)
Block *block;
struct working *col;
int pos;
{
 int seq, aa;

  col->totcnt = col->totraw = col->totreg = 0.0;
  for (aa = 0; aa < AASALL; aa++) {
    col->cnt[aa] = col->raw[aa] = col->reg[aa] = 0.0;
  }

  /*  Only count the real 20 aas, combine B(21) with D & Z(22) with E  */
  for (seq = 0; seq < block->num_sequences; seq++) {
    aa = block->residues[seq][pos];
    if (aa >= 1 && aa < AAS) {
      col->cnt[aa] += block->sequences[seq].weight;
      col->totcnt += block->sequences[seq].weight;
      col->raw[aa] += 1.0;
      col->totraw += 1.0;
    }
  } /* end of for */
} /* end of counts_nogaps */


/*========================================================================
==========================================================================*/
double similarity_dependent_scale (col, rank_matrix, original_aa)
struct working *col;
struct float_qij *rank_matrix;
int original_aa;
{
  int aa, rank; double sum;
   int n;
  int max_aa;

  max_aa = find_max_aa_col(col);
  original_aa = max_aa; /* change to max aa 10/24/00 */

  n = 0;
  sum = 0.0;
  for (aa = 1; aa <= 20 ; aa++) {
        if (col->cnt[aa] > 0.0) {
                n++;
              rank = (int) rank_matrix->value[original_aa][aa];
              sum += rank * col->cnt[aa]/col->totcnt;
/*		printf ("\taa %c rank %d weight %.3f", aa_btoa[aa], rank, col->cnt[aa]/col->totcnt); */
        }
   }
/*   printf (" sim score %.3f\n", sum); */
  return sum;

} /* end of similarity_dependent_scale */


/*=====================================================================*/
/*========================================================================
==========================================================================*/
void pseudo_diri(col, diric, epsilon)
struct working *col;
struct diri *diric;
double epsilon;
{
   int j, aa;
   double denom, dtemp;
   double total;

      /*-----------   compute equation (3), Prob(n|j) ------------  */
      for (j = 0; j < diric->ncomp; j++)
      {
         col->probn[j] = lgamma(col->totcnt + 1.0) + lgamma(diric->altot[j]);
         col->probn[j] -= lgamma(col->totcnt + diric->altot[j]);

	/*   Note range on aa varies; Diric values only for 1-20 */
         for (aa = 1; aa < AAS; aa++)
         {
         /*  ni = cnt[i] */
	    if (col->cnt[aa] >= 0.0)
            {
               dtemp = lgamma(col->cnt[aa] + diric->alpha[j][aa]);
/*  printf ("aa %c cnt %.3f alpha %.3f dtemp %.3f\n", aa_btoa[aa], col->cnt[aa], diric->alpha[j][aa], dtemp); */
		dtemp -= lgamma(col->cnt[aa] + 1.0);
	   	dtemp -= lgamma(diric->alpha[j][aa]);
		col->probn[j] += dtemp;
/*               printf ("col->probn is now %.3f prob is %.3f \n\n", col->probn[j], exp (col->probn[j])) ; */
	     } /* end of if > = 0.0 */
         } /* end of for amino acids */
      } /* end of for all j components */

      /*------ compute sum qk * p(n|k) using logs & exponents ----------*/
      denom = log(diric->q[0]) + col->probn[0];
      for (j = 1; j < diric->ncomp; j++)
      {
         dtemp = log(diric->q[j]) + col->probn[j];
         denom = add_logs(denom, dtemp);
      }
/* printf ("%.3f denom number of comp %d\n", denom, diric->ncomp); */

      /*   compute equation (3), Prob(j|n)  */
      for (j = 0; j < diric->ncomp; j++)
      {
         col->probj[j] = log(diric->q[j]) + col->probn[j] - denom;
/*       printf("j=%d probn[j]=%f probj[j]=%f\n",
             j, exp(col->probn[j]), exp(col->probj[j])); */

      }

      /* ----- compute equation (4), ni + bi, bi = Prob(j|n) * alpha[j][i]  */

	 for (aa = 1; aa < AAS; aa++)
      {
         for (j = 0; j < diric->ncomp; j++) {
		col->reg[aa] +=
				(exp(col->probj[j]) * diric->alpha[j][aa]);
			/*  * epsilon); */
            }
	 col->totreg += col->reg[aa];
	}
/*	printf (" total prob. %.2f\n", col->totreg); */
/* scale dirichlet to probabilities */
	total = 0.0;
	for (aa = 1; aa < AAS; aa ++) {
		col->reg[aa] /= col->totreg;
		total += col->reg[aa];
	}
/*	printf ("%.2f total\n", total); */

/*	printf ("%c %.3f", aa_btoa[aa], col->reg[aa]);
	if (col->cnt[aa] > 0.0) {printf ("*\n"); } else { printf ("\n");} */

}  /* end of pseudo_diri */

/*========================================================================
==========================================================================*/
int find_max_aa_col (col)
struct working *col;
{
	int aa, max_aa;
	double max = 0.0;

        max_aa = -1;
	for (aa = 1; aa < AAS; aa++) {
		if (col->cnt[aa] > max) {
			max = col->cnt[aa];
			max_aa = aa;
		}
	}
	return max_aa;

} /* end of find_max_aa_col */

/*========================================================================
==========================================================================*/
int find_max_aa_pssm (pssm, pos)
struct work_pssm *pssm;
int pos;
{
	int aa, max_aa;
	double max;

        max_aa = -1;
	max = 0.0;
	for (aa = 1; aa < AAS; aa++) {
		if (pssm->value[pos][aa] > max) {
			max_aa = aa;
			max = pssm->value[pos][aa];
		}
	}
	return max_aa;
} /* end of find_max_aa_pssm */


/*===============================================================
      Returns log(e^x + e^y)
================================================================*/
double add_logs(lx, ly)
double lx, ly;
{
   if (lx > ly)  return(lx + log(1.0 + exp(ly - lx)));
   else          return(ly + log(1.0 + exp(lx - ly)));
}  /*  end of add_logs */


/*=========================================================================*/
/* Change log information follows. */
/*
 Changes since version 3.5:
  8/27/03 Fix problem in make_alts() with scale=10,20,21,30
 Changes since version 3.4:
 12/29/00 Added Pauline's routines for SIFT (SIFT_prediction, etc.)
	  Modified struct working, struct work_pssm and make_work_pssm()
	  (width is now dynamic) Added free_work_pssm()
  6/16/00 Added type 22 (pseudo-counts nearly zero) to block_to_matrix()
 Changes since version 3.3.2:
  5/24/00 Increased MAXWIDTH from 100 to 400, but need to do away with
	  this restriction!
 Changes since version 3.3:
  9/ 4/99 Added type 21 (counts only) to block_to_matrix()
 Changes since version 3.2.5:
  2/22/99  Removed block from positive_matrix()
 Changes since version 3.2.4:
 12/12/98  compute_BZX(): Use 0 for * column instead of min value, unless
	   min value is negative.
 Changes since version 3.2.3:
  8/11/98  convert.c: MAXWIDTH increased, error check added.
  4/ 8/98  Avoid division by zero in pseudo_alts() & block_to_matrix().
           Divide counts of aas other than the basic 20 among those 20.
 Changes since version 3.1:
  1/30/97  Changed no weights error message from SERIOUS to WARNING.
 10/10/96  Fixed bug in pb_weights() when X (#23) in sequence.
  9/25/96  Changed compute_BZX() to compute average score for X.
*/
