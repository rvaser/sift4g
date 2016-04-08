/* (C) Copyright 1993-7, Fred Hutchinson Cancer Research Center */
/* Use, modification or distribution of these programs is subject to */
/* the terms of the non-commercial licensing agreement in license.h. */

/* scores.c: score structure manipulation and scoring routines */
/* Written by: Bill Alford */
/* Change log information is at the end of the file. */

/*	system headers not in global.h */
#include <math.h>
/*	blimps library headers */
#include <global.h>
#include <blocks.h>	/* includes sequences.h, output.h */
#include <matrix.h>	/* includes pattern.h */
/*	headers in current directory */
#include "scores.h"
#include "blimps.h"
#include "scoring.h"
#include "lists.h"


/*
 * Exported variables and data structures
 */

/*
 * Local variables and data structures
 */
static int output_formatted_score();
/* Never Used */
/* static void print_score(); */
/* static void print_score_summary(); */
static void score_and_enter_switch();
static void select_with_pattern();



/*
 * Function definitions
 */



/* 
 * score_and_enter
 *   Calls the appropriate function for scoring.  If it is to select patterns
 *   the function select_with_pattern() is called, otherwise it calls the 
 *   appropriate scoring method in score_and_enter_switch().
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

int score_and_enter(seq, matrix, frame, repeats, search_type)
     Sequence *seq;
     Matrix *matrix;
     int frame;
     Boolean repeats;
     int search_type;
{
  /* check for the correct sequence type */
  if (seq->type != AA_SEQ) {
    sprintf(ErrorBuffer, 
	    "score_and_enter(): Tried to score a sequence that is not an");
    ErrorReport(PROGRAM_ERR_LVL);
    sprintf(ErrorBuffer,
	    "                   amino acid sequence.  Not scoring the");
    ErrorReport(PROGRAM_ERR_LVL);
    sprintf(ErrorBuffer,
	    "                   sequence %s\n", seq->name);
    ErrorReport(PROGRAM_ERR_LVL);
    return FALSE;
  }

  if (UsePatterns) {
    /* select the patterns to score and score them */
    select_with_pattern(seq, matrix, frame, repeats, search_type);
  }
  else {
    /* we are scorint the whole sequence */
    score_and_enter_switch(seq, matrix, frame, repeats, search_type, FALSE);
  }

  return TRUE;
  
}


/*
 * score_and_enter_switch
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
 *   Return Codes: none
 *   Error Codes: none
 */

static void
score_and_enter_switch(seq, matrix, frame, repeats, search_type, 
		       pattern_search)
     Sequence *seq;
     Matrix *matrix;
     int frame;
     Boolean repeats;
     int search_type;
     Boolean pattern_search;
{
  switch (SequenceMatrixScoringMethod) {
  case 0:
  default:

    if (SequenceMatrixScoringMethod != 0) {
      sprintf(ErrorBuffer, 
	      "Invalid matrix-sequence scoring method specified, %d.",
	      SequenceMatrixScoringMethod);
      ErrorReport(WARNING_ERR_LVL);
      sprintf(ErrorBuffer,
	      "Setting the method to the default value and using the");
      ErrorReport(WARNING_ERR_LVL);
      sprintf(ErrorBuffer,
	      "default scoring method.\n");
      ErrorReport(WARNING_ERR_LVL);
      SequenceMatrixScoringMethod = 0;
    }

    default_scoring_method(seq, matrix, frame, repeats, search_type, 
			   pattern_search);

    break;

  } /* end switch of SequenceMatrixScoringMethod */
}



static void
select_with_pattern(seq, matrix, frame, repeats, search_type)
     Sequence *seq;
     Matrix *matrix;
     int frame;
     Boolean repeats;
     int search_type;
{
  int x, y;
  int seq_length = seq->length;
  int mat_length = matrix->width -1;
  int seq_pos_scored_length = (mat_length*2) + seq_length;
  Boolean seq_pos_scored[1000];

  Pattern **patterns;
  int seq_compare_start;	/* the place to start, in the seq_pos_scored */
				/* reference frame */

/*printf("mat_length: %d \tseq_length: %d\n", mat_length, seq_length);*/

  /* build the patterns for the matrix if not done so already */
  scan_patterns(matrix);

  if ((matrix->patterns == NULL) ||
      (matrix->patterns[0] == NULL)) {
    sprintf(ErrorBuffer,
	    "No pattern for matrix %s, Scoring the entire matrix\n",
	    matrix->number);
    ErrorReport(WARNING_ERR_LVL);

    score_and_enter_switch(seq, matrix, frame, repeats, search_type, FALSE);
    return;
  }

  for (x=0; x<(mat_length+seq_length); x++) {
    /* Fill the array with falses */
    seq_pos_scored[x] = FALSE;
  }

  patterns = matrix->patterns;

  for (x=0; x<seq_pos_scored_length; x++) {
    for (y=0; patterns[y] != NULL; y++) {

      /* some precalculations of reused stuff */
      seq_compare_start = x - patterns[y]->beg_offset;

      /* check that it is within range of the array */
      if ((seq_compare_start >= 0) &&
	  (seq_compare_start < seq_pos_scored_length)) {
	/* check that it hasn't been scored */
	if (seq_pos_scored[seq_compare_start] == FALSE) {
	  /* check if this is scorable */ 
	  if (pattern_matches(seq, x-mat_length, patterns[y])) {
	    
	    /* shift it back to the sequences' reference point */
	    seq->undefined = seq_compare_start - mat_length;
	    /*printf("Possition is: %d  \tScore pos is: %d\n", x-mat_length,seq->undefined);*/
	    seq_pos_scored[seq_compare_start] = TRUE;
	    score_and_enter_switch(seq, matrix, frame, repeats, 
				   search_type, TRUE);
	  }
	}
      }
    }
  }
}


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

int score_comparison(a, b)
     Score *a, *b;
{
  return (a->score - b->score);	
}

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
int neg_score_comparison(a, b)
     Score *a, *b;
{
  return (- (score_comparison(a, b))); /* reversed because we want the */
				       /* highest score first */
}

/* 
 * free_score
 *   Deletes the score list entry and the sub elements.
 *   Parameters: 
 *     Score *score: the score to free
 *   Return code: none
 *   Error code: none
 */

void free_score(score) 
     Score *score;
{
  free(score);
}

/*
 * Score printing.
 */

/* 
 * output_scores
 *   Prints out the score list in the final format.
 *   Parameters:
 *     int num_to_report: the number of scores to report, if <1 reports all.
 *     FILE *ofp:         the file to output scores to
 *   Error codes: none
 */

void output_scores(num_to_report, ofp)
     int num_to_report;
     FILE *ofp;
{
  /* make the PrintScores list */
  DoForSL(Scores, enter_score_into_print_scores, NULL);
  
  DoForSL(PrintScores, output_formatted_score, ofp);
                              /* pass the print function and the output file */
}


/* 
 * output_formatted_score
 *   Prints a score and related data to the specified file.  This function
 *   controls what the final output for the scores will look like.
 *   Parameters:
 *     Score *score:  the score to print
 *     FILE *ofp: the output file pointer.  May be stdout.
 *   Return code: SL_CONTINUE, always.  So that when going through the list
 *                it will always behave the same.  See the function 
 *                description of DoForSL in the skiplist functions.doc and 
 *                the defines in skiplist.h.
 *   Error codes: none
 */

static int output_formatted_score(score, ofp)
     Score *score;
     FILE *ofp;
{
  if (SearchType == SEARCH_TYPE_BLOCK) { /* seq vs blocks */
    fprintf(ofp, 
	    "%-11.11s %-50.50s      %4d   %4d %2d %6d %s\n",
	    score->matrix_number, /* accession number */
	    score->matrix_desc,	/* description, de or info field trimed */
	    score->strength,	/* strength */
	    score->score,	/* score */
	    score->frame,	/* reading frame */
	    score->position,	/* position in the sequence */
	    score->consensus);	/* the consensus string */

  }
  else if (SearchType == SEARCH_TYPE_MATRIX) { /* block vs sequences */
    fprintf(ofp, 
	    "%-20.20s %-62.62s %4d %2d %5d %6d %s\n",
	    score->sequence_number, /* accession number */
	    score->sequence_desc, /* description, de or info field trimed */
	    score->score,	/* score */
	    score->frame,	/* reading frame */
	    score->position,	/* position in the sequence */
	    score->seq_length,  /* length of the sequence */
	    score->consensus);		/* the consensus string */
    
  }
  else if (SearchType == SEARCH_TYPE_UNKNOWN) {	/* blocks vs sequences */
    fprintf(ofp, 
	    "%-8s vs %-7s   %-20s... & %-20s... %4d %4d %2d %6d %6d %s\n",
	    score->matrix_number, /* accession number */
	    score->sequence_number, /* accession number */
	    score->matrix_desc,	/* description, de or info field trimed */
	    score->sequence_desc, /* description, de or info field trimed */
	    score->strength,	/* strength */
	    score->score,	/* score */
	    score->frame,	/* reading frame */
	    score->position,	/* position in the sequence */
	    score->seq_length,  /* length of the sequence */
	    score->consensus);		/* the consensus string */

  }
  else {
    fprintf(ofp, 
	    "%-8s %4u %-45s      %4d   %4d %2d %6d %s\n",
	    "--------",		/* accession number */
	    0,			/* record number, I don't use */
	    "---------------------------------------------", 
				/* description, de or info field trimed */
	    score->strength,	/* strength */
	    score->score,	/* score */
	    score->frame,	/* reading frame */
	    score->position,	/* position in the sequence */
	    score->consensus);	/* the consensus string */
    
  }

  return SL_CONTINUE;
}


/*
 * print_score
 *   Prints a Score data structure.  Primarily for debugging purposes.
 *   Parameters: 
 *     Score *score:  the score to print
 *   Error Codes: none
 */
/*
static void print_score(score)
     Score *score;
{
  printf("--- score ---\n");

  printf("\tScore: %d\n", score->score);
  
  printf("\tPosition: %d\n", score->position);

  printf("\tFrame: %d\n", score->frame);

  printf("\tMat Num: %s\n", score->matrix_number);

  printf("\tSeq Num: %s\n", score->sequence_number);

  printf("\tMat Desc: %s\n", score->matrix_desc);

  printf("\tSeq Desc: %s\n", score->sequence_desc);

  printf("\tConsensus: %s\n", score->consensus);

  printf("--- ----- ---\n");
}
*/

/*
 * print_score_summary
 *   Prints a Score data structure in a compact format.  Primarily for 
 *   debugging purposes.
 *   Parameters: 
 *     Score *score:  the score to print
 *   Error Codes: none
 */
/*
static void print_score_summary(score)
     Score *score;
{
  printf("S:% 6d\tP:% 6d\tF:% 6d\n", score->score, score->position, score->frame);

}
*/



/* Change log information follows. */
/* 
 Changes since 3.3.1:
 10/ 2/99  Added more includes for standard C libraries.
 Changes since 3.2.5:
 4/29/99 output_formatted_score() Don't print the 0 before description.

 Changes since 3.2.3:
 5/13/98 Print 1st 20 instead of 1st 11 chars of sequence names.
 *
 *
 */

