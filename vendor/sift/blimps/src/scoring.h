/* (C) Copyright 1993-9, Fred Hutchinson Cancer Research Center */
/* Use, modification or distribution of these programs is subject to */
/* the terms of the non-commercial licensing agreement in license.h. */

/* scoring.h: functions for different methods of scoring a sequence against */
/*            a matrix */
/* Written by: Bill Alford */
/* Change log information is at the end of the file. */

#ifndef SCORING_H_
#define SCORING_H_

/*
 * Exported variables and data structures
 */

#define MIN_SAVE_SCORE 945	/* If SV option used, never save lower score*/
#define MAX_HIST_SCORE  6000	/* the top of the histogram */
/* MAX_HIST_SCORE is 6000 = 60 x 100 = max block width x max weight */
#define HIST_BOX_SIZE   20	/* must evenly divide MAX_HIST_SCORE */
#define NUM_HIST_COLS   ( MAX_HIST_SCORE / HIST_BOX_SIZE )

extern double Alignments_Done;
extern double Scores_Done;
extern int histogram[];
extern Boolean DoHistogram;


/*
 * Exported functions
 */

/*
 * default_scoring_method
 *   The default scoring method.  Basically a score is the sum of the scores 
 *   in the matrix at each position-amino acid pair for one alignment.
 *   Parameters:
 *     Sequence *seq: the sequence to score
 *     Matrix *matrix: the matrix used for scoring
 *     int frame: the frame we are looking at (0 = aa seq.,
 *                1, 2, 3, -1, -2, -3 = na seqs.)
 *     Boolean repeats: enter repeats into the score_list?
 *     int search_type: the type of search we are doing.  If it is a 
 *                      SEARCH_TYPE_MATRIX search (block vs sequences), 
 *                      do NOT divide by the block percentile.  If it is a 
 *                      SEARCH_TYPE_BLOCK search (sequence vs blocks), DO 
 *                      divide by the block percentile. 
 *   Error codes: none
 */

extern void default_scoring_method();


/* 
 * print_histogram
 *   ouputs the histogram array from scoring.[ch] to the output stream ofp.
 *   Parameters:
 *     FILE *ofp: the output file pointer
 *   Return codes: none
 *   Error codes: none
 */

extern void print_histogram();




#endif /*  SCORING_H_ */

/* Change log information follows. */
/*
   Changes since version 3.2.5:
    5/25/99 Renamed scores_done to Scores_Done & alignemnts_done
	    to Alignments_Done & changed them from unsigned long to double
    2/ 3/99 Added MIN_SAVE_SCORE
   10/19/94 Increased MAX_HIST_SCORE from 2000 to 6000 = 60 x 100 =
	    max block width x max weight
 */

