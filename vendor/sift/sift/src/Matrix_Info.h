/* (C) Copyright 2000, Fred Hutchinson Cancer Research Center */
/* Use, modification or distribution of these programs is subject to */
/* the terms of the non-commercial licensing agreement in license.h. */

#ifndef MATRIX_INFO_H_
#define MATRIX_INFO_H_

#define MATRIX_LENGTH_INCREASE_SIZE 80 /* the default length of the matrix */
                                       /* and the size to increase by if */
                                       /* more room is needed (doubtful) */


/* matrix with experimental info - for SNPS project 4/15/00

   like Jorja's matrix.h except additional field to store whether the
   substitution can be observed

5/26/00 changed calculate_information_R so can calculate information for
	PSSMs.  But X's can't be counted in correction factor for PSSMS.

07-12-00 removed (int) round of matrix-> weights when evaluating score so that frq. < .5 are
counted as ob
 */

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "blimps/blocksprogs.h"

#define TRUE 1
#define FALSE 0
#define AAs 20

#define allowed 0
#define less_severe 3
#define more_severe 4
#define not_allowed 1
#define no_data 2
#define NO_DATA 2000

struct aanode {
        Residue aa;
        struct aanode* next;
};
typedef struct aanode* AAnode;

extern FILE* errorfp;

/* int allowed= 0;
int less_severe = 3;
int more_severe = 4;
int not_allowed= 1;
int no_data=2;
int NO_DATA = 2000;
*/
/* I don't know why I can't build an array of this type */
/*enum InfoType { allowed, not_allowed, no_data}; */

struct matrix_with_info {

	Matrix* block_matrix;

	int *info[MATRIX_AA_WIDTH]; /* weight matrix containing the */
					/* information/experimental data */
	double* logo_heights[MATRIX_AA_WIDTH]; /* contains heights of each
						aa according to LOGO */
	double* R; /* contains total information at position as defined
			by Schneider */

	double* avg_R_per_residue; /* avg_R_per_residue : R / # of res,
				     hopefully reflects
					diverse, low info. columns */
};

typedef struct matrix_with_info Matrix_Info;


Matrix_Info* initialize_matrix_info (Matrix* matrix);

void free_Matrix_Info (Matrix_Info* matrix_info);

void read_experimental_info (Matrix_Info* matrix, FILE* infofp,
			int beg, int end, const fpos_t start_of_file);

void output_Matrix_Info (Matrix_Info* info_matrix, FILE* omfp, int option, int severity_option,
			 int style);
/* option 1: for PSSMs with negative scores.  counts 0 as neutral,
	     i.e. tolerated

   option 2: for PSSMs with + and 0 scores, such as matrices containing only
	     counts

   style : FLOAT_OUTPUT or INT_OUTPUT
*/

void output_Matrix_Info_float (Matrix_Info* info_matrix, FILE* omfp, int option, int severity_option);

/*======================================================================
calculate_information_R (Matrix_Info* matrix)

calculates the total information at position in R array and fills in
height of aa at position pos height = f (aa, pos) * R (pos) in
2-dim. array logo_heights
=======================================================================*/

double calculate_information_R (Matrix_Info* matrix, int error_correction_option);

void read_matrix_20_aa_body();

void normalize_matrix (Matrix* matrix);

void allow_min_and_above (Matrix* matrix);

void convert_experimental_info_to_values (Matrix_Info* info_matrix);

double R_at_pos (Matrix* matrix, int pos);

Residue* read_polymorphism_data (FILE* infofp, Sequence* seq);


/* May 5, 2002 allow multiple polymorphisms at a single position */

AAnode* read_multiple_polymorphism_data (FILE* infofp, Sequence* seq);

void free_AAnode (AAnode aanode);

void print_multiple_polymorphism_data (AAnode* aanodes, int length);

void print_stat_info_on_pssm (Matrix_Info* info_matrix);

int index_for_char (char c);

#endif
