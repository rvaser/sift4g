/* (C) Copyright 2000, Fred Hutchinson Cancer Research Center */
/* Use, modification or distribution of these programs is subject to */
/* the terms of the non-commercial licensing agreement in license.h. */

/* convert.c: functions for different methods of converting a block into a */
/*            matrix */
/* Written by: Bill Alford */
/* Change log information is at the end of the file. */
/* PN_convert.c copied from convert.c to allow for large block widths */

/*07-13-00 changed in counts() so that if an X is not observed, counts is NOT
evenly distributed throughout.  For gaps, counts is evenly distributed
throughout.
*/
/*	system headers not in global.h */
#ifndef _PN_CONVERT_C_
#define _PN_CONVERT_C_

#include <assert.h>
#include <math.h>

#include "PN_convert.h"
#include "queue.h"
#include "SortList.h"
#include "constants.h"

#define TOLERANCE_PROB_THRESHOLD .05
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
 * Sequence weighting methods.  Various methods of giving different
 * weights to the sequences.
 */

static double *clustered_weights(/*block*/);
static double *pre_weighted_sequences(/*block*/);
void pb_weights(/*block*/);

/*
 * Matrix construction methods.  Build matricies from blocks and
 * sequence weights.
 */

static void basic_matrix_construction(/*block, seq_weight, matrix*/);

Matrix* copy_block_info_to_matrix (Block* block)
{
  Matrix* matrix;
  char* tmp;
  char Buffertemp[LARGE_BUFF_LENGTH];
  /* Aug. 18, changed Buffer (a global variable) to local Buffertemp */

  /* get new matrix */
  matrix = new_matrix(block->width);

  /* initialize the pattern */
  matrix->patterns = NULL;

  /* copy the relevant block information into the matrix */
  matrix->block = block;

  strncpy(Buffertemp, block->id, SMALL_BUFF_LENGTH);
  /* NOTE: Chance to goof by replacing the wrong "BLOCK" */
  tmp = strstr(Buffertemp, "BLOCK");
  if (tmp != NULL) {
    strncpy(strstr(Buffertemp, "BLOCK"), "MATRIX\0", 7);
    strncpy(matrix->id, Buffertemp, SMALL_BUFF_LENGTH);
  }
  else {
    strncpy(matrix->id,
            strncat(Buffertemp, "; MATRIX", SMALL_BUFF_LENGTH - strlen(Buffertemp)),
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

  return matrix;

} /* end of copy_block_info_to_matrix */


Matrix *PN_block_to_matrix(block, conversion_method)
     Block *block;
     int conversion_method;
{
  Matrix *matrix;
  char *tmp;
  matrix =  copy_block_info_to_matrix (block);
  switch (conversion_method) {
  case 2: /* use the weights from the block file */
    pre_weighted_conversion_method (block, matrix);
    break;
  case 3:   /* Altschul's data-dependent method of computing a PSSM */
	    /* from a block */
    PN_altschul_data_dependent_conversion_method(block, matrix, 0);
    break;
  case 5:   /* Version of case 3 that scales the matrix in half bits*/
    PN_altschul_data_dependent_conversion_method(block, matrix, 2);
    break;
  case 6:   /* Version of case 3 that scales the matrix in third bits*/
    PN_altschul_data_dependent_conversion_method(block, matrix, 3);
    break;
  case 10:   /* Version of case 3 that returns the odds ratios */
    PN_altschul_data_dependent_conversion_method(block, matrix, 10);
    break;
  case 20:   /* Version of case 3 that returns the counts+pseudo-counts */
    PN_altschul_data_dependent_conversion_method(block, matrix, 20);
    break;
  case 21:   /* Version of case 3 that returns the counts */
    PN_altschul_data_dependent_conversion_method(block, matrix, 21);
    break;
  case 22:  /* instead of adding RTot * diffaas*/
	    /* as in option 20, add RTot * (diffass - 1), thus not adding*/
	    /* any pseudocounts to completely conserved columns */
	    /* qij's as pseudocounts */
    PN_altschul_data_dependent_conversion_method (block, matrix, 22);
    break;
  case 23: /* qij's as pseudocounts (pseudoalts) */
	   /* exp. pseudocounts*/
    PN_altschul_data_dependent_conversion_method (block, matrix, 23);
    break;
  case 24:
    PN_altschul_data_dependent_conversion_method (block, matrix, 24);
    break;
  case 25:
    PN_altschul_data_dependent_conversion_method (block, matrix, 25);
    break;
  case 26:
    PN_altschul_data_dependent_conversion_method (block, matrix, 26);
    break;
  case 27:
    /* dirichlet pseudocounts, exp weighting, div by max to score */
    PN_altschul_data_dependent_conversion_method (block, matrix, 27);
    break;
  case 28:
    /* assign a scoring matrix such as BLOSUM62's sij's to matrix.
	take first sequence as the standard */
    scoring_matrix_profile (block, matrix);
    break;
  case 29:
        /*probabilities estimated by dirichlet,  no scaling to condiitional prob*/
	PN_altschul_data_dependent_conversion_method (block, matrix, 29);
	break;
  case 30:
	/* 27 with cutoff PROB_INTOLERANCE_THRESHOLD so that scores >= 0.0
	are tolerant */
	PN_altschul_data_dependent_conversion_method (block, matrix, 30);
	break;
  case 31:
	qij_matrix_profile (block, matrix, Qij);
	break;
   default: /* the default case */
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

/** AUGUST 3, 2000 -- make PN_block_to_matrix_RTot so you can pass in RTot &
set the number of pseudocounts */

Matrix *PN_block_to_matrix_RTot_par(block, conversion_method, RTotPar)
     Block *block;
     int conversion_method;
     double RTotPar;
{
  Matrix *matrix;
  char *tmp;

  matrix =  copy_block_info_to_matrix (block);

  if (conversion_method == 22 || conversion_method == 23 || conversion_method == 25) {
	RTot = RTotPar;
       PN_altschul_data_dependent_conversion_method (block, matrix, conversion_method);
        RTot = 5.0; /* set it back to default */
  } else {
   	 sprintf(ErrorBuffer,
            "Invalid block to matrix conversion method specified, %d.",
            conversion_method);
	    ErrorReport(WARNING_ERR_LVL);

	 sprintf(ErrorBuffer,     /* ^^^^----------------vvvvvvvvvvvvvvvvvvv */                    "Using the default conversion method of Altschul's data-dependent method.\n");

	 ErrorReport(WARNING_ERR_LVL);
    	altschul_data_dependent_conversion_method(block, matrix, 0);
  }

  /* return the matrix */
  return matrix;


} /* end PN_block_to_matrix_Rtotpar */



/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/
/* ######  PN modifications here, for allocating width memory ######*/
/*
 * altschul data dependent conversion method aux. structures/defines
 */

/*
 * altschul data dependent conversion method aux. functions
 */

/*==================================================================*/
struct working *make_col()
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
/* removed static */
struct work_pssm *PN_make_pssm(int length)
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
}  /* end of PN_make_pssm */

/* only count valid amino acids, no considerationto gaps */
void counts_no_gaps (Block* block, struct working* col, int pos)
{
 int seq, aa, aa1;

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
} /* end of subroutine */

/*======================================================================*/
void counts(block, col, pos)
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
    else { /* not a basic residue */
/*
      sprintf(ErrorBuffer,
	      "Uncounted \"residue\" character for %s: %c\n",
	      block->number, aa_btoa[block->residues[seq][pos]]);
      ErrorReport(INFO_ERR_LVL);
*/
      /* If not one of the basic aas, divide the count among them gaps*/
  if (aa_btoa[aa] != 'X') {
	for (aa1 = 1; aa1 < AAS; aa1++)
      {
         col->cnt[aa1] += (block->sequences[seq].weight / 20.0);
         col->raw[aa1] += (1.0 / 20.0) ;
      }
      col->totcnt += block->sequences[seq].weight;
      col->totraw += 1.0;
   } /* end of if not X */
   }
  }
}  /* end of counts */


/*=====================================================================*/
int count_residues(col)
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
static void pseudo_alts(col, qij, beta)
     struct working *col;
     struct float_qij *qij;
     double beta;
{
  int aa, row;

   col->totreg = 0.0;

  /*---------- get the pseudo counts -------------------------------*/
  for (aa=1; aa < AAS; aa++) {
    col->reg[aa] = 0.0;
    for (row = 1; row < AAS; row++) {
      col->reg[aa] += (col->cnt[row] * qij->value[aa][row] / qij->marg[row]);
    }
    col->reg[aa] *= beta;
    if (col->totcnt > 0.0) col->reg[aa] /= col->totcnt;
    col->totreg += col->reg[aa];
  }
}  /* end of pseudo_alts */

/*===============Pauline's range ***********************************/

int
find_min_aa_in_pssm (struct working* col, struct work_pssm* pssm, const int pos)
{
	int aa, min_aa = 0;
	double min_value;

	min_value = 1000;

	for (aa = 0; aa < AAS; aa++) {
		if (col->cnt[aa] > 0.0 && pssm->value[pos][aa] < min_value) {
			min_value = pssm->value[pos][aa];
			min_aa = aa;
		}
	}

	return min_aa;

} /* end of find_min_aa_in_pssm */

void range_on_pssm (struct working* col, struct work_pssm* pssm, const int pos)
{
        int min_aa, aa;
	double min;
	double threshold;

	threshold = 0.9;

	min_aa = find_min_aa_in_pssm(col, pssm, pos);
        min = pssm->value[pos][min_aa] * threshold;
        for (aa = 0; aa < AAS; aa++) {
                pssm->value[pos][aa] -= min;
	}


}

/* assign scoring matrix's values (such as BLOSUM62) to matrix, using
first sequence in block as the standard */
void
scoring_matrix_profile (Block* block, Matrix* matrix)
{
	// int sij[24][24];
	FILE* matrixfp;
	char matrix_file[LARGE_BUFF_LENGTH];
	int pos, original_aa, aa;

	/* read matrix */
    // char* blimps_dir = getenv ("BLIMPS_DIR");
	// sprintf (matrix_file, "%s/docs/blosum62.bla.new", blimps_dir);

/*	strcpy (matrix_file, "/tuna/blocks_5.0/blosum/blosum62.bla.new");*/

   // if ( (matrixfp = fopen (matrix_file, "r")) == NULL) {
    //    fprintf (errorfp, "can't read the matrix %s\n", matrix_file);
    //    exit (-1);
   //}
    //    read_substitution_scoring_matrix (matrixfp, sij);

	//fclose (matrixfp);
    int (*sij)[24] = default_blosum62;

	for (pos = 0; pos < block->width; pos++) {
		original_aa = block->residues[0][pos];
		for (aa = 1; aa < AAS; aa++) {
			matrix->weights[aa][pos] = sij[original_aa -1][aa -1];
		}
	}

} /* end of scoring_matrix_profile */

void
qij_matrix_profile (Block* block, Matrix* matrix, struct float_qij *qij)
{
	int pos, original_aa, aa;

	for (pos = 0; pos < block->width; pos++) {
		original_aa = block->residues[0][pos];
		for (aa = 1; aa < AAS; aa++) {
			matrix->weights[aa][pos] = qij->value[aa][original_aa];
			matrix->weights[aa][pos] /= qij->value[original_aa][original_aa];
					/* Henikoff qij->marg[original_aa] ; */
			matrix->weights[aa][pos] -= TOLERANCE_PROB_THRESHOLD;
		}
	}

} /* end of qij_matrix_profile */


static void Pauline_pseudo_alts (col, qij, epsilon)
	struct working* col;
     struct float_qij *qij;
     double epsilon;
{
  int aa, row;

	col->totreg = 0.0;
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
}  /* end of Pauline_pseudo_alts */

/*=====================================================================*/
static void log_odds(pos, freqs, col, pssm)
     int pos;
     double *freqs;
     struct working *col;
     struct work_pssm *pssm;
{
  int aa;

  pssm->sum[pos] = 0.0;
  for (aa=1; aa < AAS; aa++)
  {
    pssm->value[pos][aa] = col->cnt[aa] + col->reg[aa];
    pssm->value[pos][aa] /= (col->totcnt + col->totreg);

    /*     compute the odds ratio   */
    if (freqs[aa] > 0.0) {
      pssm->value[pos][aa] /= freqs[aa];
    }

    /*   take the log of the odds ratio  */
    if (pssm->value[pos][aa] > 0.0) {
      pssm->value[pos][aa] = log(pssm->value[pos][aa]);
    }

    pssm->sum[pos] += pssm->value[pos][aa];
  }
}  /* end of log_odds */


/*========================================================================
	Computes scores for B, Z, X, -(gap) and *(stop) in a column
	of a PSSM using other scores in the column
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

void
initialize_cutoff_freqs (double* newfreqs, double* oldfreqs)
{
	int aa;

/* old version. changed July 2, 2001 to make it linux compatible
	assert (oldfreqs[MATRIX_AA_WIDTH] != NULL);
	assert (newfreqs[MATRIX_AA_WIDTH] != NULL);
*/

	assert (oldfreqs != NULL);
	assert (newfreqs != NULL);

	for (aa = 0; aa < MATRIX_AA_WIDTH; aa++) {
		newfreqs[aa] = 0.05; /* oldfreqs[aa]; */
	}
}


/* pseudocounts from Altschul's 1997 NAR gapped blast and PSI-BLAST paper */
void
psiblast_alts (Block* block, Matrix* matrix, double* freqs, struct float_qij* qij)
{
	double beta;
	int pos, alpha, aa, original_aa;
	double min_freq;
	struct working *col;
        struct work_pssm *pssm;

	beta = 10.0;
        col = make_col();
	pssm = PN_make_pssm(block->width);

 /*--------------  Do one position at a time -------------------*/
  for (pos = 0; pos < block->width; pos++)
  {
    /*-------- count the number of each aa in this position ------------*/
    counts(block, col, pos);
    alpha = count_residues(col);
    alpha--;

    pseudo_alts(col, qij, beta);
    /* col->reg contains beta * gi's */

/*---------   Fill in the matrix entries --------------------*/
   pssm->sum[pos] = 0.0;
    	for (aa=1; aa < AAS; aa++)
    	{
       	   pssm->value[pos][aa] = (alpha*col->cnt[aa]/col->totcnt)
				 + col->reg[aa];
       		 if ( (alpha + beta) > 0.0)
             	pssm->value[pos][aa] /= (alpha + beta);
    	} /* end of aa, pssm->values filled for a given pos */

	/* turn into conditional probability on orignal amino acid*/
        original_aa = block->residues[0][pos];
        min_freq = pssm->value[pos][original_aa];

   	for (aa = 1; aa < AAS; aa++) {
       		pssm->value[pos][aa] /= min_freq;
       		pssm->value[pos][aa] -= 0.05;
       		pssm->sum[pos] += pssm->value[pos][aa];
       		matrix->weights[aa][pos] = pssm->value[pos][aa];
   	}  /* end of aa */
  } /* end of for pos */

 free(col);
 free(pssm->value[0]);
 free(pssm->value); free(pssm->sum);
 free(pssm);

}

int
find_max_aa_in_col (struct working* col)
{
	int aa, max_aa = 1;
	double max = 0.0;

	for (aa = 1; aa < AAS; aa++) {
		if (col->cnt[aa] > max) {
			max = col->cnt[aa];
			max_aa = aa;
		}
	}
	return max_aa;

} /* end of find_max_aa_in_col */

int
find_max_aa_in_pssm (struct work_pssm* pssm , int pos)
{
	int aa, max_aa = 1;
	double max;

	max = 0.0;
	for (aa = 1; aa < AAS; aa++) {
		if (pssm->value[pos][aa] > max) {
			max_aa = aa;
			max = pssm->value[pos][aa];
		}
	}
	return max_aa;
}

/*==========================================================================
     Uses Altschul's method of getting pseudo-counts with a qij matrix,
===========================================================================*/
static void PN_make_alts(block, matrix, freqs, qij, RTot, scale)
     Block *block;
     Matrix *matrix;
     double *freqs;
     struct float_qij *qij;
     double RTot;               /* Total R for Altschul */
     int scale;
{
  double avg_score, factor, dtemp, dtemp2, epsilon, original_score;
  int pos, aa,  itemp, normal_aa;
  struct working *col;
  struct work_pssm *pssm;
  int odds; int log_odds;
  int original_aa;
  struct dirichlet* diric;
  int min_aa, max_aa;
  double min_freq = 1, no_of_gaps;
  double cutoff_freq[MATRIX_AA_WIDTH];
  int rank_matrix[AAS][AAS];
  double number_of_observed_aas;
  int diri_pseudo;
  char diri_file[LARGE_BUFF_LENGTH];

	if (scale == 30) {
		SIFT_alts (block, matrix, freqs, qij, 1,0, 1, 0);
                return;
	}
  diri_pseudo = FALSE;
  factor = 1.0;
  if (scale > 0 && scale < 10) factor = (double) scale / log(2.0);
  col = make_col();
  pssm = PN_make_pssm(block->width);

  // char* blimps_dir = getenv ("BLIMPS_DIR");
  // printf ("blimps dir is %s\n", blimps_dir);
  // sprintf (diri_file , "%s/docs/default.diri", blimps_dir);
  // printf  ("diri file is %s\n", diri_file);

/*  diric = load_diri_file ("/howard2/pauline/src/merge-opt.13compnew"); */
  // diric = load_diri_file (diri_file);
  diric = &default_dirichlet;
  construct_rank_matrix ( rank_matrix);
printf ("diri file loaded\n");
  /*--------------  Do one position at a time -------------------*/
  for (pos = 0; pos < block->width; pos++)
  {
    /*-------- count the number of each aa in this position ------------*/
    counts(block, col, pos);
    /*-------- determine total number of pseudo-counts in column ------*/
    epsilon = 0;

    itemp = count_residues (col);
    if (scale == 27 || scale == 29 || scale == 30) {
	if (itemp ==1 ) {epsilon = 0; }
	else {
		original_aa = block->residues[0][pos];
                dtemp = similarity_dependent_scale_0 ( col, rank_matrix, original_aa);
		epsilon = exp (dtemp);
       }

    } else if (scale == 20) {
        epsilon = (double) RTot * count_residues(col) ;
    } else if (scale == 22 || scale == 23 || scale == 25 || scale == 26) {
        if (itemp == 1) {
                epsilon = 0;
        } else {
                epsilon = (double) RTot * itemp;
        }
   }

    /*---------- get the pseudo counts -------------------------------*/
    if (scale == 27 || scale == 29 || scale == 30 || scale == 26) {
                diri_pseudo = TRUE;
		pseudo_diric (col, diric, epsilon );
    } else {
        pseudo_alts(col, qij, epsilon);
    }
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
       if (scale <= 20 || scale == 22 || scale == 25 || scale == 23
                || scale == 26)

       {
          pssm->value[pos][aa] = col->cnt[aa] + col->reg[aa];
        if ( (col->totcnt + col->totreg) > 0.0)
             pssm->value[pos][aa] /= (col->totcnt + col->totreg);
       }
	if (diri_pseudo)  { /* diri has frequencies instead of actual counts*/
		pssm->value[pos][aa] = col->cnt[aa] + (epsilon * col->reg[aa]);
		if ((col->totcnt + col->totreg) > 0.0) {
			pssm->value[pos][aa] /= (col->totcnt + epsilon);
		} else {
			fprintf (stderr, "ERROR!  Dividing by 0.\n");
		}
	}

    } /* end of aa, pssm->values filled for a given pos */

   if (scale == 27 || scale == 26 || scale == 30) {
        max_aa = find_max_aa_in_pssm (pssm, pos);
        /* original_aa = block->residues[0][pos]; */
        min_freq = pssm->value[pos][max_aa];
/*      printf ("divide pos %d by %.2f\n", pos, min_freq); */
    } /* end of if scale == 27 */

   for (aa = 1; aa < AAS; aa++) {
        if (scale == 27) {
                pssm->value[pos][aa] /= min_freq;
        } else if (scale == 30 || scale == 26) {
                pssm->value[pos][aa] /= min_freq;
                pssm->value[pos][aa] /= freqs[aa];
        	pssm->value[pos][aa] = log (pssm->value[pos][aa]);
		pssm->value[pos][aa] += 0.05; /* to compensate forthreshold */
	}

       pssm->sum[pos] += pssm->value[pos][aa];

       /*  scale the matrix  */
       if (scale)
       {
          dtemp = factor * pssm->value[pos][aa];
          matrix->weights[aa][pos] = dtemp;
        }

     }  /* end of aa */

     /*compute_BZX(frequency, matrix, pos);*/

        original_aa = block->residues[0][pos];
        if (matrix->weights[original_aa][pos] < 0.0) {
                fprintf (errorfp, "WARNING!!! Amino acid %c at pos %d in your original sequence was not allowed by the prediction.\n", aa_btoa[original_aa], pos);
        }

        } /* end of for pos */
  /*------ Now make the final scores; make log scores non-neg    */
  /*       by subtracting the min. value for all positions ----- */
  if (!scale) positive_matrix(freqs, pssm, matrix);

printf ("left make _alts\n");
  // free (diric);
  free(col);
  free_struct_work_pssm (pssm );

}  /* end of PN_make_alts */

/*==========================================================================
10/20/99 SIFT pseudocounts
13-dirichlet component
option:  diri_pseudocounts
	  TRUE: use 13-component Dirichlet
	  FALSE use BLOSUM62 qij's
option:  gap: TRUE -allow everything
	      FALSE - just look at 20 amino acids
option:  exp (m)
	 m=0 # of amino acids at pos (default)
	 m=1  similiarity scale (more conservative)

option:	subtract_threshold TRUE : scores -= TOLERANCE_PROB_THRESHOLD
			   FALSE: leave as original scores
===========================================================================*/

Matrix*
SIFT_prediction (Block* block, int diri_option, int gap_option, int exp_option, int subtract_option)
{
	Matrix* pssm;

	pssm =  copy_block_info_to_matrix (block);

        if (Qij == NULL) {
                printf ("qij matrix missng\n");
                exit (-1);
        }
        normalize  (block);
        SIFT_alts(block, pssm, frequency, Qij, diri_option, gap_option,
                          exp_option, subtract_option);

/*	printf ("SIFT predictions\n");
	print_matrix (pssm); */
	return pssm;
}

double
similarity_dependent_scale_0_from_block (Block* block, int pos)
{
	double dtemp;
	struct working* col;
	int original_aa;
	int rank_matrix[AAS][AAS];

	construct_rank_matrix ( rank_matrix);

	col = make_col();
	counts_no_gaps(block, col, pos);
	original_aa = block->residues[0][pos];
	dtemp = similarity_dependent_scale_0(col, rank_matrix, original_aa);
	free (col);
	return dtemp;

}
void calculate_counts (Block* block)
{
	int pos, aa;
	struct working* col;
	col = make_col();

	for (pos = 0; pos < block->width; pos++)
	{
		counts_no_gaps (block, col, pos);
		for (aa = 1; aa < AAS; aa++) {
			printf ("pos %d aa %c counts %.2f weight %.2f\n",
			pos, aa_btoa[aa], col->cnt[aa], col->totcnt);
		}
	}
}

void SIFT_alts(Block* block, Matrix* matrix, double* freqs,
		struct float_qij* qij,
		int diri_pseudocounts, int gap_option,
		int exp_option, int subtract_threshold)
{
  double dtemp, dtemp2, epsilon;
  int pos, aa,  itemp;
  struct working *col;
  struct work_pssm *pssm;
  int original_aa;
  struct dirichlet* diric;
  int min_aa, max_aa;
  double min_freq, no_of_gaps;
  double number_of_observed_aas;
  int div_by_max;
  int rank_matrix[AAS][AAS];
  FILE* rfp;
  char diri_file[LARGE_BUFF_LENGTH];

  div_by_max = TRUE;

  col = make_col();
  pssm = PN_make_pssm(block->width);

  // char* blimps_dir = getenv ("BLIMPS_DIR");
  // sprintf (diri_file , "%s/docs/default.diri", blimps_dir);

  // diric = load_diri_file (diri_file);
  diric = &default_dirichlet;

  if (exp_option == 1) {

	construct_rank_matrix ( rank_matrix);
   }

  /*--------------  Do one position at a time -------------------*/
  for (pos = 0; pos < block->width; pos++)
  {
    /*-------- count the number of each aa in this position ------------*/
    if (gap_option) {
	counts (block,col, pos);
    } else {
	counts_no_gaps(block, col, pos);
   }

    /*-------- determine total number of pseudo-counts in column ------*/
    epsilon = 0;
    itemp = count_residues (col);
    original_aa = block->residues[0][pos];

    if (exp_option == 1) {
	/*printf ("pos %d :  ", pos); */
	dtemp = similarity_dependent_scale_0(col, rank_matrix, original_aa);
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
	pseudo_diric (col, diric, epsilon );
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
        max_aa = find_max_aa_in_pssm (pssm, pos);
/*printf ("%c",aa_btoa[max_aa]) ; */
	min_freq = pssm->value[pos][max_aa];
       for (aa = 1; aa < AAS; aa++) {
		pssm->value[pos][aa] /= min_freq;
	}
   }

   if (subtract_threshold) {
	   for (aa = 1; aa < AAS; aa++) {
	        pssm->value[pos][aa] -= TOLERANCE_PROB_THRESHOLD;
   	}
   }

       pssm->sum[pos] += pssm->value[pos][aa];

    	for (aa = 1; aa < AAS; aa++) {
	  matrix->weights[aa][pos] = pssm->value[pos][aa];
	}
/*	if (gap_option) {
		assign_gap_positive_scores (matrix, pos);
	} */
	original_aa = block->residues[0][pos];
	if (  (matrix->weights[original_aa][pos] < TOLERANCE_PROB_THRESHOLD
			 && (!subtract_threshold))
		||
	      (matrix->weights[original_aa][pos] < 0.0 && subtract_threshold)
	) {
		fprintf (errorfp, "WARNING!!! Amino acid %c at pos %d in your original sequence had score %.3f and was not allowed by the prediction.<BR>\n", aa_btoa[original_aa], pos + 1, matrix->weights[original_aa][pos]);
	}
	} /* end of for pos */

  if (diri_pseudocounts) {
	// free (diric);
  }

  free(col);
  free(pssm->value[0]);
  free(pssm->value); free(pssm->sum);
  free(pssm);

}  /* end of PN_make_alts */

void SIFT_alts_test(Block* block, Matrix* matrix, double* freqs,
                struct float_qij* qij,
                int diri_pseudocounts, int gap_option,
                int exp_option, int subtract_threshold)
{
  double dtemp, dtemp2, epsilon;
  int pos, aa,  itemp;
  struct working *col;
  struct work_pssm *pssm;
  int original_aa;
  struct dirichlet* diric;
  int min_aa, max_aa;
  double min_freq, no_of_gaps;
  double number_of_observed_aas;
  int div_by_max;
  int rank_matrix[AAS][AAS];
  FILE* rfp;
  char diri_file[LARGE_BUFF_LENGTH];
  char matrix_file[LARGE_BUFF_LENGTH];

  div_by_max = TRUE;


  col = make_col();
  pssm = PN_make_pssm(block->width);
  diri_pseudocounts = TRUE;

  // char* blimps_dir = getenv ("BLIMPS_DIR");
  // sprintf (diri_file , "%s/docs/default.diri", blimps_dir);

  // diric = load_diri_file (diri_file);
  diric = &default_dirichlet;
                /* 13 component Dirichlet 01/02/00  */

  if (exp_option == 1) {

        construct_rank_matrix ( rank_matrix);
   }
  // sprintf (matrix_file, "%s/docs/default.qij", blimps_dir);
  init_frq_qij_for_matrix (matrix_file);

  /*--------------  Do one position at a time -------------------*/
  for (pos = 0; pos < block->width; pos++)
  {
    /*-------- count the number of each aa in this position ------------*/
    if (gap_option) {
        counts (block,col, pos);
    } else {
        counts_no_gaps(block, col, pos);
   }
    /*-------- determine total number of pseudo-counts in column ------*/
    epsilon = 0;
    itemp = count_residues (col);
    original_aa = block->residues[0][pos];

    if (exp_option == 1) {
        /*printf ("pos %d :  ", pos); */
        dtemp = similarity_dependent_scale_0(col, rank_matrix, original_aa);
        if (itemp == 1) {
                epsilon = 0;
        } else {
                epsilon = dtemp * 100;
        }
    } else {
        if (itemp == 1) { epsilon = 0; }
        else if (itemp > 7) { epsilon = 1000; }
        else { epsilon = exp( (double) itemp) ; }

   }
    /*---------- get the pseudo counts -------------------------------*/
    if (diri_pseudocounts) {
        pseudo_diric (col, diric, epsilon );
    } else {
        pseudo_alts(col, qij, epsilon);
    }

    /*---------   Fill in the matrix entries --------------------*/
    pssm->sum[pos] = 0.0;
    for (aa=1; aa < AAS; aa++)
    {

        if (diri_pseudocounts) {
		pssm->value[pos][aa] = col->cnt[aa] + epsilon* col->reg[aa];
        	if ( (col->totcnt + col->totreg) > 0.0)
             		pssm->value[pos][aa] /= (col->totcnt + epsilon);
	} else {
		pssm->value[pos][aa] = col->cnt[aa] + col->reg[aa];
		if ( (col->totcnt + col->totreg) > 0.0)
                        pssm->value[pos][aa] /= (col->totcnt + col->totreg);
	}

    } /* end of aa, pssm->values filled for a given pos */

   if (div_by_max) {
        max_aa = find_max_aa_in_pssm (pssm, pos);
/*printf ("%c",aa_btoa[max_aa]) ; */
        min_freq = pssm->value[pos][original_aa];
       for (aa = 1; aa < AAS; aa++) {
                pssm->value[pos][aa] = min_freq - pssm->value[pos][aa];
        }
   }

   if (subtract_threshold) {
           for (aa = 1; aa < AAS; aa++) {
                pssm->value[pos][aa] -= TOLERANCE_PROB_THRESHOLD;
        }
   }

       pssm->sum[pos] += pssm->value[pos][aa];

        for (aa = 1; aa < AAS; aa++) {
          matrix->weights[aa][pos] = pssm->value[pos][aa];
        }
/*      if (gap_option) {
                assign_gap_positive_scores (matrix, pos);
        } */
        original_aa = block->residues[0][pos];
        if (  (matrix->weights[original_aa][pos] < TOLERANCE_PROB_THRESHOLD
                         && (!subtract_threshold))
                ||
              (matrix->weights[original_aa][pos] < 0.0 && subtract_threshold)

        ) {
            fprintf (errorfp, "WARNING!!! Amino acid %c at pos %d in your original sequence had score %.3f and was not allowed by the prediction. <BR>\n", aa_btoa[original_aa], pos, matrix->weights[original_aa][pos]);
        }

        } /* end of for pos */

  if (diri_pseudocounts) {
    // free (diric);
  }
  init_frq_qij(); /* return back to normal */
  free(col);
  free(pssm->value[0]);
  free(pssm->value); free(pssm->sum);
  free(pssm);

}  /* end of SIFT_alts_test */

double
count_gaps (Matrix* matrix, int pos)
{
        double gaps;
        int seq;

        assert (matrix->block != NULL);
        gaps = 0;
        for (seq = 0; seq < matrix->num_sequences; seq++) {
                if (matrix->block->residues[seq][pos] == aa_atob['-'] ) {
                        gaps += matrix->block->sequences[seq].weight;
                }
        }
	return gaps;
}

void
gap_pseudocounts ( struct working* col, double no_of_gaps)
{
	double pseudocount;
	int aa;

	pseudocount = exp (no_of_gaps);
	for (aa = 1; aa < AAS; aa++) {
		col->reg[aa] += pseudocount;
		col->totreg += pseudocount;
	}

}

void
ratio_mutated_to_normal (struct working* col, int standard_aa)
{
	int aa;

	for (aa = 1; aa < AAS; aa++) {
		col->reg[aa] /= col->reg[standard_aa];
	}

} /* end ratio_mutated_to_normal */

/***********end of PN_mixed_alts ***************/

void
assign_gap_positive_scores (Matrix* matrix, int pos)
{
	int gaps;
	int seq;

	assert (matrix->block != NULL);
/* 	printf ("entering gap_positive scores\n"); */
	gaps = 0;
	for (seq = 0; seq < matrix->num_sequences; seq++) {
		if (matrix->block->residues[seq][pos] == aa_atob['-'] ) {
			gaps++;
		}
	}
	if (gaps == 0) { /* printf ("no gaps\n"); */return; }
	else { /* printf ("gaps %d\n", gaps) ; */
			assign_positive_scores_to_column (matrix, pos, gaps); }

} /* end of assign_gap_positive_scores */

void
assign_positive_scores_to_column (Matrix* matrix, int pos, int gaps)
{
	double min; int seq;
	int aa;

	min = 0.0;
	for (aa = 1; aa < AAS; aa++) {
		matrix->weights[aa][pos] = 0.05;
	}
} /* end of positive_column */

/*=========================================================================
    Normalize sequence weights to add up to the number of sequences
=========================================================================*/
/*static void normalize(block) deleted because made extern in BLIMPS 3.5 */

/*
 * altschul_data_dependent_conversion_method
 *   Expects global variables frequency[], Qij[][] and RTot
 *   Parameters:
 *     Block *block:   the block to be converted
 *     Matrix *matrix: where the resulting matrix will be put
 *   Error codes: none
 */


void
PN_altschul_data_dependent_conversion_method(block, matrix, scale)
     Block *block;
     Matrix *matrix;
     int scale;
{
  assert (Qij != NULL);
  assert (frequency != NULL);

  if (Qij == NULL) {
      sprintf(ErrorBuffer, "Qij matrix missing, unable to continue.\n");
      ErrorReport(FATAL_ERR_LVL);
    }

  /* check to see if the block has sequence weights */
printf ("enter normalize block\n");
  normalize(block);   /*  Make weights sum to number of sequences */
printf ("enter PN_make_alts\n");
  PN_make_alts(block, matrix, frequency, Qij, RTot, scale);
printf ("exit PN_make_alts\n");
}

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

void
multiply_pssm_by_100 (Matrix* pssm)
{
	int pos, aa;

	for (pos = 0; pos < pssm->width; pos++) {
		for (aa = 1; aa < AAS; aa++) {
			pssm->weights[aa][pos] *= 100;
		}
	}
}

Sequence*
get_consensus (Matrix* matrix)
{
	Residue res;
	Sequence* seq;
	int pos;

	seq = (Sequence *) malloc (sizeof (Sequence));
	seq->sequence = (Residue *) calloc (matrix->width, sizeof (Residue) );
	strcpy (seq->name, "CONSENSUS");
	seq->info[0] = '\0';
	seq->position = 0;
	seq->length = matrix->width; seq->max_length=matrix->width;
	seq->type = AA_SEQ; seq->weight = 0;

	for (pos = 0; pos < matrix->width; pos++) {
		res = find_max_aa_in_pos (matrix, pos);
		seq->sequence[pos] = res;
		printf ("res %c pos %d\n", aa_btoa[seq->sequence[pos]], pos);
	}
	return seq;

} /* end get_consensus */

Residue
find_max_aa_in_pos (Matrix* matrix, int pos)
{
        int aa;
        Residue max_aa;
        double max = 0.0;

	max_aa = aa_atob['X']; /* defaultis X */
        for (aa = 1; aa < AAS; aa++) {
                if (matrix->weights[aa][pos] > max) {
                        max = matrix->weights[aa][pos];
                        max_aa = aa;
                }
        }
        return max_aa;

} /* end of find_max_aa_in_pos */

/* choose a residue in proporition to its weight */
Residue
sampling_residue (MatType *weights[MATRIX_AA_WIDTH], int pos)
{
	double total;
	int aa;
	double col[MATRIX_AA_WIDTH];
	int random_aa ;
	double cum_prob[MATRIX_AA_WIDTH];

	for (aa = 0; aa < MATRIX_AA_WIDTH; aa++) {
		col[aa] = 0.0;
		if (( aa != aa_atob['B']) &&
		    (aa != aa_atob['Z']) &&
          	    (aa != aa_atob['X']) &&
          	    (aa != aa_atob['-']) && /* gap */
          	    (aa != aa_atob['*']) && /* stop */
          	    (aa != AAID_NAR)) {     /* non */
		col[aa] = weights[aa][pos];
		total += col[aa];
	/*	printf ("%.3f total %.3f, " , col[aa], total); */
		}
	}
/*	if (total == 0.0) { is this causing me to get res '-' ?
		return (aa_atob['X']);
	} else { */
		cumulative_prob (col, cum_prob);
		random_aa = random_pick (cum_prob);
		return (random_aa); /* returning an int -- is this defined
					in the binary representation of aa */
/*	} */
}

void
cumulative_prob (double col[MATRIX_AA_WIDTH], double cum_prob[MATRIX_AA_WIDTH])
{
	int i;
	double index;


	index = 0;

	for (i = 0; i <MATRIX_AA_WIDTH; i++) {
		cum_prob[i] = index + col[i];
		index = cum_prob[i];
	}

}


/* returns a state , given a cumulative probability vector P with
MATRIX_AA_WIDTH states ranging from 0 to 1 */
int
random_pick (double* P)
{
        double random;
        int i = 0;

/*	printf ("random pick: "); */
        random =   (double) rand() /  (double) RAND_MAX * 100 ;
        while (random > P[i] && i < MATRIX_AA_WIDTH) {
/*                printf ("| %.3f, %c |", P[i], aa_btoa[i]); */
		i++;
        }
/*	printf ("picked %c with rand %.3f prob %.3f\n", aa_btoa[i], random, P[i]); */
	if (i == MATRIX_AA_WIDTH) { /* no aa in this column, return an X */
		return (aa_atob['X']);
	} else {
        	return (i);
	}
}

void output_matrix_sPN (matrix, omfp, style)
     Matrix *matrix;
     FILE *omfp;
     int style;
{
  char c;
  int l;

  if (style == INT_OUTPUT) {
    fprintf(omfp, "ID   %s\n", matrix->id);
    fprintf(omfp, "AC   %s\n", matrix->ac);
    fprintf(omfp, "DE   %s\n", matrix->de);
    fprintf(omfp, "MA   %s\n", matrix->ma);

    fprintf(omfp, " ");
    for (c='A'; c <= 'Z'; c++) {
      if ((c != 'J') && (c != 'O') && (c != 'U')) {
        fprintf(omfp, " %c  ", c);
      }
    }
    fprintf(omfp, " *  ");
    fprintf(omfp, " -\n");

    for (l=0; l<matrix->width; l++) {
      for (c='A'; c <= 'Z'; c++) {
        if ((c != 'J') && (c != 'O') && (c != 'U')) {
          fprintf(omfp, "%3d ", (int) round(matrix->weights[aa_atob[c]][l]));
        }
      }
      fprintf(omfp, "%3d ", (int) round(matrix->weights[aa_atob['*']][l]));
      fprintf(omfp, "%3d\n", (int) round(matrix->weights[aa_atob['-']][l]));
    }

    fprintf(omfp, "//\n");
  }
  else if (style == FLOAT_OUTPUT) {
    fprintf(omfp, "ID   %s\n", matrix->id);
    fprintf(omfp, "AC   %s\n", matrix->ac);
    fprintf(omfp, "DE   %s\n", matrix->de);
    fprintf(omfp, "MA   %s\n", matrix->ma);

    fprintf(omfp, "pos ");
    for (c='A'; c < 'Z'; c++) {
      if ((c != 'J') && (c != 'O') && (c != 'U') && (c != 'B') && (c != 'X') ) {
        fprintf(omfp, "  %c  ", c);
      }
    }
   /* fprintf(omfp, " *  ");
    fprintf(omfp, " -\n");
*/
    fprintf (omfp, "\n");
    for (l=0; l<matrix->width; l++) {
      fprintf (omfp, "%2d ", l + matrix->block->sequences[0].position);
      for (c='A'; c < 'Z'; c++) {
        if ((c != 'J') && (c != 'O') && (c != 'U') && (c != 'B') && (c != 'X') ) {
          if (matrix->weights[aa_atob[c]][l] > 0) {
		fprintf(omfp, "  %.2f", matrix->weights[aa_atob[c]][l]);
          } else {
		fprintf(omfp, " %.2f", matrix->weights[aa_atob[c]][l]);
	  } /* end of print if matrix wegiths > 0 */
        }
      }
/*      fprintf(omfp, " %.1f ", matrix->weights[aa_atob['*']][l]);
      fprintf(omfp, " %.1f\n", matrix->weights[aa_atob['-']][l]); */
      fprintf (omfp, "\n");
    }

    fprintf(omfp, "//\n");
  }
  else { /* unknown */
    sprintf(ErrorBuffer,
            "Unknown output type: %d, using integer output\n", style);
    ErrorReport(WARNING_ERR_LVL);
    output_matrix_s(matrix, omfp, INT_OUTPUT);
  }
}  /*  end of output_matrix_s() */

void read_substitution_scoring_matrix (fin, scores)
FILE *fin;
int scores[24][24];
{
   char line[132], *ptr;
   int alpha[24], nrows, ncols, row, col, i;

/*----------Read file until first non-blank line --------------*/
/* Skip comments at beginning of file - 1st char = #, > or ;   */
   line[0] = '\0';
   while ((strlen(line) < (size_t) 1 ||
           line[0]=='#' || line[0]=='>' || line[0]==';')
          && fgets(line, sizeof(line), fin) != NULL)
            ;
/*------See if the first line has characters on it ------------*/
   for (col=0; col < 24; col++) alpha[col] = -1;
   if (strstr(line, "A") != NULL)       /* This line has characters */
   {
      row = 0;  /* # of alphabetic characters on the line */
      for (i=0; i < (int) strlen(line); i++)
      {
         col = -1;
         if (line[i] == 'A') col = 0;
         if (line[i] == 'R') col = 1;
         if (line[i] == 'N') col = 2;
         if (line[i] == 'D') col = 3;
         if (line[i] == 'C') col = 4;
         if (line[i] == 'Q') col = 5;
         if (line[i] == 'E') col = 6;
         if (line[i] == 'G') col = 7;
         if (line[i] == 'H') col = 8;
         if (line[i] == 'I') col = 9;
         if (line[i] == 'L') col = 10;
         if (line[i] == 'K') col = 11;
         if (line[i] == 'M') col = 12;
         if (line[i] == 'F') col = 13;
         if (line[i] == 'P') col = 14;
         if (line[i] == 'S') col = 15;
         if (line[i] == 'T') col = 16;
         if (line[i] == 'W') col = 17;
         if (line[i] == 'Y') col = 18;
         if (line[i] == 'V') col = 19;
         if (line[i] == 'B') col = 20;
         if (line[i] == 'Z') col = 21;
         if (line[i]=='O' || line[i]=='X') col = 22;
         if (line[i]=='J' || line[i]=='*') col = 23;
         if (col >= 0)
         {
            alpha[row] = col;
            row++;
         }
         else if (isalpha(line[i])) row++;
      }
   }
/*-------Get the data values now ------------*/
   for (row=0; row<24; row++)
     for (col=0; col< 24; col++)
        scores[row][col] = -999;                /* Null value */
   nrows = 0;
   line[0] = '\0';
   while (fgets(line, sizeof(line), fin) != NULL)
   {
      if (strlen(line) > (size_t) 1 && nrows < 24)
      {
         if (alpha[nrows] >= 0 && alpha[nrows] < 24)
         {
            row = alpha[nrows]; ncols = 0;
            ptr = strtok(line, " ,\n");
            while (ptr != NULL)
            {
               if (strspn(ptr, "+-0123456789") == strlen(ptr))
               {
                  col = alpha[ncols];
                  if (col >= 0 && col < 24)
                     scores[row][col] = atoi(ptr);
                  ncols++;
               }
               ptr = strtok(NULL, " ,\n");
            }
         }
         nrows++;
      }
   }

/*-------If some entries are still missing, assume symmetry ---------*/
   for (row=0; row<24; row++)
   {
     for (col=0; col<24; col++)
     {
        if (scores[row][col] == -999) scores[row][col] = scores[col][row];
        if (row==20 && scores[row][col]==-999)  /*  B -> D */
           scores[row][col] = scores[3][col];
        if (row==21 && scores[row][col]==-999)  /* Z -> E */
           scores[row][col] = scores[6][col];
        if (col==20 && scores[row][col]==-999)  /*  B -> D */
           scores[row][col] = scores[row][3];
        if (col==21 && scores[row][col]==-999)  /* Z -> E */
           scores[row][col] = scores[row][6];
     }
   }
} /* end of read_substitution_scoring_matrix */

/*=====================================================================*/
/*   This routine initializes the blimps global variables
        Qij, RTot and frequency[]
Input Parameters .qij file
		 .out file (Jorja's format, to get frq)
*/
void init_frq_qij_for_matrix(char qijname[LARGE_BUFF_LENGTH])
{
   FILE *fqij;

   if (Qij != NULL) {
		free (Qij);
		Qij = NULL;
   }
   if ( (fqij=fopen(qijname, "r")) != NULL)
   { Qij = load_qij(fqij); fclose(fqij);  }
   else {printf ("could not open %s qij file\n", qijname); }
   RTot = LOCAL_QIJ_RTOT;
}  /* end of frq_qij */

/*===============================================================
      Returns log(e^x + e^y)
================================================================*/
double addlogs(double lx, double ly)
{
   if (lx > ly)  return(lx + log(1.0 + exp(ly - lx)));
   else          return(ly + log(1.0 + exp(lx - ly)));
}  /*  end of addlogs */


/*=====================================================================*/
void pseudo_diric(struct working* col, struct dirichlet* diric, double epsilon)
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
         denom = addlogs(denom, dtemp);
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

}  /* end of pseudo_diric */

/* aa observed in block that has small pseudocount */
int
find_min_aa_in_col_reg (struct working* col)
{
        int aa, min_aa = 1;
        double min;

        min = 10000;

        for (aa = 1; aa < AAS; aa++) {
                if (col->cnt[aa] > 0.0 && col->reg[aa] < min) {
                        min_aa = aa;
                        min = col->reg[aa];
                }
        }
        return min_aa;
}

char
find_min_aa_in_col (struct working* col)
{
	int aa, min_aa = 1;
	double min;

	min = 1000;

	for (aa = 1; aa < AAS; aa++) {
		if (col->cnt[aa] > 0.0 && col->cnt[aa] < min) {
			min_aa = aa;
			min = col->cnt[aa];
		}
	}
	return (aa_btoa[min_aa]);
}


/*======================================================================
   Dirichlet input file order (0-19) is  ARNDCQEGHILKMFPSTWYV
   Blimps order (0-24) is               -ARNDCQEGHILKMFPSTWYVBZX*
   Store Dirichlet alphas in positions 1-20 to match blimps
=======================================================================*/
struct dirichlet *load_diri_file (char filename[LARGE_BUFF_LENGTH])
{
	FILE* fin;
	int i, aa, numc, type;
        char line[MAXLINE], *ptr;
        struct dirichlet *diric;
	double denom;
	double background_frequency[21];
    printf ("filename is %s\n", filename);
	if ((fin = fopen (filename, "r")) == NULL) {
		fprintf (errorfp, "Cannot open dirichlet file %s\n", filename);
		exit (-1);
	}

   diric = (struct dirichlet *) malloc(sizeof(struct dirichlet));
   if (diric == NULL)
   {
      fprintf(errorfp, "\nOUT OF MEMORY\n");
      exit(-1);
   }
   numc = 0;
   while (numc < MAXCOMP && !feof(fin) && fgets(line, MAXLINE, fin) != NULL)
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



void
cutoff_target_freq (double* cutoff_freq, int original_aa, struct float_qij *qij)
{

	int aa;

	for (aa = 1; aa <= 20; aa++)  {
		cutoff_freq[aa] = qij->value[aa][original_aa] / qij->marg[original_aa];
	}
} /* end of cutoff_target_freq */


/*=====================================================================*/
/*   This routine initializes the blimps global variables
        Qij, RTot and frequency[]       */
void init_frq_qij()
{
   // FILE *fqij;
   // char qijname[80], frqname[80];
   // char* blimps_dir = getenv ("BLIMPS_DIR");
   // if (blimps_dir != NULL)
   // {
   //   sprintf(frqname, "%s/docs/default.amino.frq", blimps_dir);
   //   sprintf(qijname, "%s/docs/default.qij", blimps_dir);
   // }
   // else
   // {
   //   sprintf(frqname, "default.amino.frq");
   //   sprintf(qijname, "default.qij");
   // }
   // load_frequencies(frqname);           /* creates frequency[]  */
   get_default_frequencies();
   // Qij = NULL;
   // if ( (fqij=fopen(qijname, "r")) != NULL)
   // { Qij = load_qij(fqij); fclose(fqij);  }
   Qij = &default_qij;
   RTot = LOCAL_QIJ_RTOT;
}  /* end of frq_qij */

double
similarity_dependent_scale_0 (struct working* col,
                                int rank_matrix[AAS][AAS],
                                int original_aa)
{

  int aa, rank; double sum;
        double total;
   int n;
  int max_aa;

  max_aa = find_max_aa_in_col(col);
  original_aa = max_aa; /* change to max aa 10/24/00 */

  n = 0;
  sum = 0.0;
  for (aa = 1; aa <= 20 ; aa++) {
        if (col->cnt[aa] > 0.0) {
                n++;
              rank = rank_matrix[original_aa][aa];
              sum += rank * col->cnt[aa]/col->totcnt;
/*		printf ("\taa %c rank %d weight %.3f", aa_btoa[aa], rank, col->cnt[aa]/col->totcnt); */
        }
   }
/*   printf (" sim score %.3f\n", sum); */
  return sum;

}

void
construct_rank_matrix (int rank_matrix[AAS][AAS])
{
        struct rank_cell aalist[AAS];
        int original_aa;
	struct float_qij* rank_Qij;
	// FILE* qijfp;
	int aa;
	// char filename[80];

    // char* blimps_dir = getenv ("BLIMPS_DIR");
	// sprintf (filename, "%s/docs/default.sij", blimps_dir);

	// if ((qijfp = fopen (filename, "r")) == NULL) {
	//	fprintf (stderr, "Can't opeen blosum62.sij file in construct_rank_matrix looked in %s\n", blimps_dir);
	//	exit (-1);
	//}

	//rank_Qij = load_qij (qijfp);
    //fclose (qijfp);
    rank_Qij = &default_sij;

    for (original_aa = 0; original_aa < AAS; original_aa++) {
        copy_values_to_ranklist (aalist, original_aa, rank_Qij);
        sort_ranks (aalist);
        assign_ranks (rank_matrix, aalist, original_aa);
    }
	// free (rank_Qij);
}

void
assign_ranks (int rank_matrix[AAS][AAS], struct rank_cell aalist[AAS],
                int original_aa)
{
        int rank, aa;

        for (rank = 0; rank < AAS; rank++) {
                aa = aalist[rank].aa;
                rank_matrix[original_aa][aa] = rank + 1;
			/* otherwise rank starts from 0
                             to AAS-1 rather than 1 to AAS */
        }
}

void
copy_values_to_ranklist (struct rank_cell aalist[AAS], int original_aa,
			struct float_qij* rankQij) {

        int aa;

	/* 01/03/00 rankQij has -1 for gap, which will mess up
	ranks.  Assign it a really negative value so it will have
	the lowest rank */

	aalist[0].aa = 0;
	aalist[0].value = -50;

        for (aa = 1; aa < AAS; aa++) {
                aalist[aa].aa = aa;
                aalist[aa].value = rankQij->value[original_aa][aa];
        }
}

void sort_ranks (struct rank_cell aalist[AAS])
{
        int i, j;
        struct rank_cell tmp;

        for (i = 1; i < AAS; i++) {
                tmp.aa = aalist[i].aa;
                tmp.value = aalist[i].value;
                for (j = i; j > 0 && aalist[j-1].value < tmp.value; j--) {
                        aalist[j].value = aalist[j-1].value;
                        aalist[j].aa = aalist[j-1].aa;
                }
                aalist[j].value = tmp.value;
                aalist[j].aa = tmp.aa;
        }
}

void
megaprior (Block* block, struct dirichlet* diric)
{
        double s, k, n, b;
        int i, pos, aa;
        struct working *col;

        b = n = 0.0;
        k = 10;

        for (i = 0; i < diric->ncomp; i++) {
                b += diric->altot[i];
        }

        col = make_col();

        for (pos = 0; pos < block->width; pos++)
        {
                counts(block, col, pos);
                n += count_residues(col);
        }

        s = k * n / b;

        for (i = 0; i < diric->ncomp; i++) {
                for (aa = 1; aa < AAS; aa++) {
                        diric->alpha[i][aa] *= s;
                }
                diric->altot[i] *= s;
        }

}

/* counts the number of real amino acids at a given position */
double
number_of_real_aminoacids (Block* block, int pos)
{
	int seq;
	int count;

	count = 0;
	for (seq = 0; seq < block->num_sequences; seq++) {
		if (block->residues[seq][pos] >= 1 &&
			block->residues[seq][pos] <=20) {
			count++;
		}
	}
	return count;
} /* end of number_of_real_aminoacids */

double*
calculate_basic_aa_fraction (Block* block)
{
	double* fraction_stored;
	double percent_obsaa_div_seq;
	int pos;

	fraction_stored = (double *) calloc (block->width, sizeof (double));
	for (pos = 0; pos < block->width; pos++) {
		percent_obsaa_div_seq = (double) number_of_real_aminoacids (block, pos);
		percent_obsaa_div_seq /= block->num_sequences;
		fraction_stored[pos] = percent_obsaa_div_seq;
	}
	return fraction_stored;
}

int*
calculate_basic_aa_stored (Block* block)
{
        int* residues_stored = NULL;
        int pos;

        residues_stored = (int *) calloc (block->width, sizeof (int));
        for (pos = 0; pos < block->width; pos++) {
                residues_stored[pos] = number_of_real_aminoacids (block, pos);
        }
        return residues_stored;
}

void
free_struct_work_pssm (struct work_pssm* pssm)
{
	free(pssm->value[0]); /* holy cow, Ithink this works. This points
                               to space allocated in double pointer*/
	free(pssm->value); /* this points to an array of block length .
		  each cell points to values allocated by double pointer*/
	free(pssm->sum);
	free(pssm);
	pssm = NULL;

}

#endif
