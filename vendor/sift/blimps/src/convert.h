
/* Use, modification or distribution of these programs is subject to */
/* the terms of the non-commercial licensing agreement in license.h. */

/* convert.h: functions for different methods of converting a block into a */
/*            matrix */
/* Written by: Bill Alford */
/* Change log information is at the end of the file. */

#ifndef CONVERT_H_
#define CONVERT_H_

/*    Defaults for Altschul's substitution probability pseudo-count
	conversion method (3)  */
#define LOCAL_QIJ_FILE "default.qij"
#define LOCAL_QIJ_RTOT 5.0
#define AAS 21

/*
 * Exported variables and data structures
 */
extern double RTot;		/* total number of pseudo-counts */

struct float_qij {		/* qij matrix and marginal frequencies */
  double value[AAS][AAS];
  double marg[AAS];
};
extern struct float_qij *Qij;

struct pb_counts {		/* for pb_weights() */
  double diffaas;		/* # of different aas in position */
  double naas[MATRIX_AA_WIDTH]; /* # of each aa in position */
};

/* Pauline Ng's definitions (renamed) */

/* Dirichlet stuff */
#define MAXDIRI 40 /* max # of Dirichlet components (was MAXCOMP) */
/* Dirichlet mixtures  (was struct dirichlet) */
struct diri {   
   int ncomp;           /* number of components */
   double q[MAXDIRI];
   double altot[MAXDIRI];
   double alpha[MAXDIRI][AAS];
   double alpha_normalized[MAXDIRI][AAS]; /* Pauline, relative freq. of amino
                                        acid to component (normalized against
                                        component length*/
   double frequency_to_background_ratio[MAXDIRI][AAS];
                                /* frequency to background ratio as calculated
                                in TableII of Sjolander's 1996 CABIOS paper*/
};

#define SIFT_TOLERANCE .05	/* was TOLERANCE_PROB_THRESHOLD */

/*
 * Exported functions
 */
extern struct float_qij *load_qij( FILE *fin);
extern void pb_weights();
extern void normalize();


/*
 * block_to_matrix
 *   converts a block into a matrix (possition specific matrix?) according
 *   to the rule specified in BlockToMatrixConversionType.  The block field
 *   of the matrix is set in this function.
 *   Parameters:
 *     Block *block: the block to convert
 *   Return codes: Returns the pointer to the new Matrix.
 *   Error codes:
 */

extern Matrix *block_to_matrix();







/*
 * original_conversion_method
 *   The original conversion method.  This is done by weighted average of the 
 *   clusters.  This follows the method in patmat.
 *   Parameters:
 *     Block *block:   the block to be converted
 *     Matrix *matrix: where the resulting matrix will be put
 *   Error codes: none
 */

extern void original_conversion_method();


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
extern void original_conversion_method_cleaned_up();

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
extern void pre_weighted_conversion_method();

extern void pb_weights();

extern void altschul_data_dependent_conversion_method();
extern void gribskov_conversion_method();

/* Pauline's routines (renamed) */
extern void SIFT_conversion_method();
extern void SIFT_pssm();			/* SIFT_alts() */
extern struct diri *load_diri();		/* load_diri_file() */
extern void counts_nogaps();			/* counts_no_gaps */
extern double similarity_dependent_scale();	/* ..._0() */
extern void pseudo_diri();			/* pseudo_diric() */
extern int find_max_aa_col();			/* ..._in_col() */
extern int find_max_aa_pssm();			/* ..._in_pssm() */
extern double add_logs();			/* add_logs() */


#endif /*  CONVERT_H_ */

/* Change log information follows. 
    Changes since version 3.4:
Made normalize() external
12/29/00 Added Pauline's SIFT stuff.
    Changes for version 3.0.1:
    5/23/96 Added gribskov_conversion_method()
    Changes for version 3.0.0:
    9/16/95 Added pb_weights() routine.  JGH
 *  9/9/95 Added LOCAL_QIJ_FILE & LOCAL_QIJ_RTOT   JGH
*/
