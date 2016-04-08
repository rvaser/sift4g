/* (C) Copyright 1993-7, Fred Hutchinson Cancer Research Center */
/* Use, modification or distribution of these programs is subject to */
/* the terms of the non-commercial licensing agreement in license.h. */

/* scores.h: score structure manipulation and scoring routines */
/* Written by: Bill Alford */
/* Change log information is at the end of the file. */

#ifndef SCORES_H_
#define SCORES_H_

#include <matrix.h>

/*
 * Exported variables and data structures
 */

struct score_struct {
  int score;			/* the score of the matrix, used as the key */
  int strength;			/* the block/matrix strength */
  int position;			/* the position in the sequence where we */
				/* started comparing */
  int seq_length;		/* the length of the sequence scored */
  int frame;			/* translated nucleic acid frame (1, 2, 3, */
				/* -1, -2, -3). a frame of zero is an amino */
				/* acid sequence */
  char matrix_number[NUMBER_WIDTH]; /* the matrix/block number */
  char sequence_number[NUMBER_WIDTH]; /* the sequence number/name */
  char matrix_desc[DESC_WIDTH];	/* the matrix/block description */
  char sequence_desc[DESC_WIDTH]; /* the sequence description */
  char consensus[DESC_WIDTH];	/* the "consensus" sequence */
};
typedef struct score_struct Score;


/*
 * Exported functions
 */

/* 
 * score_and_enter
 *   computes the score of the sequence against the matrix.  If repeats is
 *   TRUE, all of the pairwise scores are entered into the score_list, 
 *   otherwise only the best score is entered.
 *   The scoring algorithm is selected by setting SequenceMatrixScoringMethod
 *   appropriately.  
 *   Note: the sequence must be an AA_SEQ.
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
 *   Return Codes: TRUE if the sequence is an AA_SEQ, FALSE otherwise.  The
 *                 score is inserted into the scores list.
 *   Error Codes: returns FALSE if the sequence was not an AA_SEQ and does not 
 *                score the sequence.
 */

extern int score_and_enter();


/*
 * score_comparison
 *   Compares two score records.   It compares by the value in the score
 *   field.  
 *   Parameters:
 *     ScoreListEntry a, b: the entries to compare
 *   Return codes:  a return value > 0 if a > b, a return value = 0 if a == b,
 *                  and a return value < 0 if a < b
 *   Error codes: none
 */

extern int score_comparison();


/*
 * neg_score_comparison
 *   Returns the negative of score_comparison.  Effectively the opposite
 *   logic of score_comparison. 
 *   Parameters:
 *     ScoreListEntry a, b: the entries to compare
 *   Return codes:  a return value < 0 if a < b, a return value = 0 if a == b,
 *                  and a return value > 0 if a > b
 *   Error codes: none
 */
extern int neg_score_comparison();


/* 
 * free_score
 *   Deletes the score list entry and the sub elements.
 *   Parameters: 
 *     Score *score: the score to free
 *   Return code: none
 *   Error code: none
 */

extern void free_score();


/* 
 * output_scores
 *   Prints out the score list in the final format.
 *   Parameters:
 *     int num_to_report: the number of scores to report, if <1 reports all.
 *     FILE *ofp:         the file to output scores to.
 *   Error codes: none
 */

extern void output_scores();



#endif /*  SCORES_H_ */

/* Change log information follows. */
/* 
 *
 * Revision 0.110  1993/08/10  02:50:53  billa
 * System version update.  Added convert.[ch], scoring.[ch] and version.[ch].
 *
 * Revision 0.106  1993/08/10  02:26:37  billa
 * Moved the code in score_and_enter that did the scoring of a sequence and
 * a matrix into the function default_scoring_method in scoring.[ch].
 * Also moved create_score_list_entry from scores.c to scoring.c.
 *
 * Revision 0.101  1993/08/04  19:47:14  billa
 * Reorganized file structure and changed includes.
 * Moved to scores.[ch] from lists.[ch]:
 *   score_comparison, neg_score_comparison, free_score,
 *   create_score_list_entry
 * Moved to scores.[ch] from matrix.[ch]:
 *   score_and_enter
 *
 * Revision 0.100  1993/08/03  17:47:24  billa
 * Creation.  System version 0.100.  Pre file structure reorganization.
 *
 */

