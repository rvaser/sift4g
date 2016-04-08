/* (C) Copyright 1993, Fred Hutchinson Cancer Research Center */
/* Use, modification or distribution of these programs is subject to */
/* the terms of the non-commercial licensing agreement in license.h. */

/* matrix.h: Matrix datatype and matrix manipulation functions */
/* Written by: Bill Alford */
/* Change log information is at the end of the file. */

#ifndef MATRIX_H_
#define MATRIX_H_

#include "pattern.h"

/*
 * Exported variables and data structures
 */

typedef double MatType;		/* so we can easily change to double if */
				/* needed.  NOTE:  read_matrix_body will */
				/* need to be changed to reflect any changes */
				/* to the type of MatType. */
				/* and convert.c it rounds values to an int */
				/* before storing */
				/* WARNING!!! in convert.c there is a LOT of */
				/* rounding! (and I mean a LOT). */

#define NUMBER_WIDTH 20
#define DESC_WIDTH   80

#define MATRIX_AA_WIDTH     26	/* 20 aa, 3 ambig, 1 stop, 1 gap, and */
				/* 1 for non-codes J, O, U. */
				/* the +2 to account for the gap and */
                                /* codes that are not defined ( - and */
				/* J, O, U ) */

struct matrix_struct {
  Block *block;			/* the block that this matrix is from */
  char id[SMALL_BUFF_LENGTH];	/* Matrix ID string */
  char ac[SMALL_BUFF_LENGTH];	/* Matrix AC string */
  char de[DESC_WIDTH];		/* Matrix DE string */
  char ma[SMALL_BUFF_LENGTH];	/* Matrix MA string */
  char number[NUMBER_WIDTH];	/* the number of this matrix */
  char motif[20];		/* motif */
  int width;			/* the length of the sequence/matrix */
  int num_sequences;		/* the number of sequences that were in the */
				/* generating data. */
  int percentile;		/* the block/matrix percentile */
  int strength;			/* the block/matrix strength */
  int max_length;		/* max length when reading in */
  Pattern **patterns;		/* the patterns for this PSSM, an array
				   of ptrs to Patterns, NULL */
  int undefined;		/* for future use */
  void *undefined_ptr;		/* for future use */
  MatType *weights[MATRIX_AA_WIDTH]; /* the weight matrix */
};
typedef struct matrix_struct Matrix;



/*
 * Exported functions
 */

/*
 * read_a_matrix
 *   reads a matrix from the data base and returns a pointer to the new
 *   matrix data structure
 *   Parameters:
 *     FILE *mfp: a pointer to the database file with the matrix.
 *   Error codes: NULL if a matrix was not read
 */

extern Matrix *read_a_matrix ();


/*
 * matrix_comparison
 *   Compares two matricies.   It compares by the value in the block->id
 *   field if it exists.
 *   Parameters:
 *     MatrixListEntry a, b: the entries to compare
 *   Return codes:  a return value < 0 if a < b, a return value = 0 if a == b,
 *                  and a return value > 0 if a > b
 *   Error codes: none
 */

extern int matrix_comparison();


/*
 * new_matrix
 *   allocates space for a matrix of length len, sets up the data structure
 *   and returns a pointer to the Matrix.
 *   Parameters:
 *     int len: the length of the sequence/matrix
 *   Error codes:
 */

extern Matrix *new_matrix();


/*
 * free_matrix
 *   Deletes the matrix and the sub elements.
 *   Parameters:
 *     Matrix *matrix: the matrix to free
 *   Return code: none
 *   Error code: none
 */

extern void free_matrix();


/*
 * print_matrix
 *   Prints a Matrix data structure.  Primarily for debugging purposes.
 *   Parameters:
 *     Matrix *matrix:  the matrix to print
 *   Error Codes: none
 */

extern void print_matrix();


/*
 * output_matrix
 *   Outputs a matrix data structure to the given file.
 *   Parameters:
 *     Matrix *matrix: the matrix to output
 *     FILE   *omfp:   the output matrix file pointer
 *   Return code: none
 *   Error code: none
 */

extern void output_matrix();


/*
 * output_matrix_s
 *   Outputs a matrix data structure to the given file with the
 *   specified style of data.
 *   Parameters:
 *     Matrix *matrix: the matrix to output
 *     FILE   *omfp:   the output matrix file pointer
 *     int    style:   the kind of output (INT_OUTPUT, FLOAT_OUTPUT)
 *   Return code: none
 *   Error code: none
 */

extern void output_matrix_s();

/*
 * output_matrix_st
 *   Outputs a matrix data structure to the given file with the
 *   specified style of data. Optionally only outputs ACTG.
 *   Parameters:
 *     Matrix *matrix: the matrix to output
 *     FILE   *omfp:   the output matrix file pointer
 *     int    style:   the kind of output (INT_OUTPUT, FLOAT_OUTPUT)
 *     int    type:    the type of matrix (AA_SEQ, NA_SEQ)
 *   Return code: none
 *   Error code: none
 */

extern void output_matrix_st();




#endif /*  MATRIX_H_ */

/* Change log information follows. */
/*
  Changes since 3.2.5:
  2/11/99  Renamed Matrix->length to Matrix->width.
  Changes since 3.1:
  2/14/97  Added output_matrix_st().
  Changes since 3.0.0:
  3/25/96   Changed MatType from int to double.   JGH
*/
