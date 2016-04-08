/* (C) Copyright 1993, Fred Hutchinson Cancer Research Center */
/* Use, modification or distribution of these programs is subject to */
/* the terms of the non-commercial licensing agreement in license.h. */

/* convert.h: functions for different methods of converting a block into a */
/*            matrix */
/* Written by: Bill Alford */
/* Change log information is at the end of the file. */

/* #################Pauline's stuff ###################*/
/* this version of PN_Convert.h is modified to allow for variable
PSSM width */
#ifndef PN_CONVERT_H_
#define PN_CONVERT_H_

#include <stdio.h>

#include "blimps/blocksprogs.h"

/* structures */

#define MAXCOMP 40 /* max # of Dirichlet components */

extern FILE* errorfp;

/* Dirichlet stuff */
struct dirichlet {      /* Dirichlet mixture information */
   int ncomp;           /* number of components */
   double q[MAXCOMP];
   double altot[MAXCOMP];
   double alpha[MAXCOMP][AAS];
   double alpha_normalized[MAXCOMP][AAS]; /* Pauline, relative freq. of amino
					acid to component (normalized against
					component length*/
   double frequency_to_background_ratio[MAXCOMP][AAS];
				/* frequency to background ratio as calculated
				in TableII of Sjolander's 1996 CABIOS paper*/
};

/* these are from JGH's makealts, I am not sure where these are from exctly. */
/* I think they are the actual values from the aabet.h structures. */
/* AAS and AASALL are not matched up with any of the aabet.h defines. */
#define AASALL 26
/* 0- 1A 2R 3N 4D 5C 6Q 7E 8G 9H 10I 11L 12K 13M 14F 15P 16S 17T 18W 19Y
   20V 21B 22Z 23X 24* 25J,O,U */


struct working {        /* Working information for one column */
  double cnt[AASALL];           /* Sequence-weighted counts */
  double totcnt;
  double raw[AASALL];           /* unweighted counts */
  double totraw;
  double reg[AASALL];           /*  pseudo counts */
  double totreg;
  double probn[MAXCOMP];        /* for dirichlet */
  double probj[MAXCOMP];        /* for dirichlet */
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
 * block_to_matrix
 *   converts a block into a matrix (possition specific matrix?) according
 *   to the rule specified in BlockToMatrixConversionType.  The block field
 *   of the matrix is set in this function.
 *   Parameters:
 *     Block *block: the block to convert
 *   Return codes: Returns the pointer to the new Matrix.
 *   Error codes:
 */

extern Matrix *PN_block_to_matrix();

extern void PN_altschul_data_dependent_conversion_method();

void PN_mixed_alts();

extern void multiply_pssm_by_100 (Matrix* pssm);

Sequence* get_consensus (Matrix* matrix);

Residue sampling_residue (MatType *weights[MATRIX_AA_WIDTH], int pos);

void cumulative_prob (double col[MATRIX_AA_WIDTH], double cum_prob[MATRIX_AA_WIDTH]);

int random_pick (double* P);

void assign_gap_positive_scores (Matrix* matrix, int pos);

void assign_positive_scores_to_column (Matrix* matrix, int pos, int gaps);

void read_substitution_scoring_matrix ();

void init_frq_qij_for_matrix (char qijname[LARGE_BUFF_LENGTH]);

double addlogs(double lx, double ly);

void pseudo_diric(struct working* col, struct dirichlet* diric, double epsilon);

struct dirichlet *load_diri_file (char filename[LARGE_BUFF_LENGTH]);

int allowed_diri (int component, int aa);

char find_min_aa_in_diri (struct working* col, double diri_array[MAXCOMP][AAS], int comp);

char find_min_aa_in_col (struct working* col);

void add_diri_pseudo(struct working* col, struct dirichlet* diric );

void ratio_mutated_to_normal (struct working* col, int standard_aa);

double min_aa_with_positive_zero_matrix_score (int reference_aa,
						struct working* col,
                                                         int sij[24][24],
                                         struct float_qij* qij, int* min_aa);

void sphere_on_marginals (struct working* col,
                        struct float_qij* qij,int sij[24][24],
                        double epsilon, int original_aa);
void sphere_pseudocounts (Block* block, Matrix* matrix, double* freqs, struct float_qij* qij);

/*static void normalize(); # compile error 2/14/2002??? didn't get this
before, comment it out */

void megaprior (Block* block, struct dirichlet* diric);

void inverse_megaprior (struct working* col, struct dirichlet* diric);

void information_dependent_cutoff (double* cutoff_freq, struct working* col, int num_sequences);

void initialize_cutoff_freqs (double* newfreqs, double* oldfreqs);

void scoring_matrix_profile (Block* block, Matrix* matrix);

void information_dependent_scale_1 (struct work_pssm *pssm, int pos);

void information_scale_0 (struct work_pssm *pssm, int pos, struct working* col, int num_sequences);

void information_scale_1 (struct work_pssm *pssm, int pos, struct working* col, int num_sequences);

void cutoff_target_freq (double* cutoff_freq, int original_aa, struct float_qij *qij);

double avg_score_of_observed_aas (struct work_pssm *pssm, struct working* col, int pos);

double similarity_dependent_scale_0 ( struct working* col,
				int rank_matrix[AAS][AAS],
				int original_aa);

void construct_rank_matrix (int rank_matrix[AAS][AAS]);

void assign_ranks (int rank_matrix[AAS][AAS], struct rank_cell aalist[AAS], int original_aa);

void copy_values_to_ranklist (struct rank_cell aalist[AAS], int original_aa,
			struct float_qij* rankQij);

void sort_ranks (struct rank_cell aalist[AAS]);

double standard_deviation (struct work_pssm *pssm, struct working* col, int pos);

struct working *make_col();

void counts();

int count_residues(col);

void psiblast_alts (Block* block, Matrix* matrix, double* freqs, struct float_qij* qij);

int find_max_aa_in_pssm (struct work_pssm* pssm, int pos);

void init_frq_qij();

double count_gaps (Matrix* matrix, int pos);

void gap_pseudocounts ( struct working* col, double no_of_gaps);

void qij_matrix_profile (Block* block, Matrix* matrix, struct float_qij *qij);

static void SIFT_alts(Block* block, Matrix* matrix, double* freqs,
                struct float_qij* qij,
                int diri_pseudocounts, int gap_option,
                int exp_option, int subtract_threshold);

Residue find_max_aa_in_pos (Matrix* matrix, int pos);

void free_struct_work_pssm (struct work_pssm* pssm);


void SIFT_alts_test(Block* block, Matrix* matrix, double* freqs,
                struct float_qij* qij,
                int diri_pseudocounts, int gap_option,
                int exp_option, int subtract_threshold);

#endif /*  CONVERT_H_ */
