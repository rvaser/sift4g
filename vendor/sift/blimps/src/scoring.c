/* (C) Copyright 1993-9, Fred Hutchinson Cancer Research Center */
/* Use, modification or distribution of these programs is subject to */
/* the terms of the non-commercial licensing agreement in license.h. */

/* scoring.c: functions for different methods of scoring a sequence against */
/*            a matrix */
/* Written by: Bill Alford */
/* Change log information is at the end of the file. */

/*	system headers not in global.h  */
#include <math.h>
/*	blimps library headers */
#include <global.h>
#include <blocks.h>	/* includes sequences.h, output.h */
#include <matrix.h>	/* includes pattern.h */
#include <residues.h>
/*	headers in current directory */
#include "scoring.h"
#include "blimps.h"
#include "scores.h"
#include "lists.h"
/* #include <files.h>	*/	/* for DoHistogram */



/*
 * Exported variables and data structures
 */

double Alignments_Done = 0.0;
double Scores_Done = 0.0;
int histogram[NUM_HIST_COLS];  /* according to ansii C these initialize as 0 */
Boolean DoHistogram = FALSE;

/*
 * Local variables and data structures
 */

static void enter_score();
static Score *make_score();


/*
 * Function definitions
 */

/*************************************************************************
 *************************************************************************
 *                                                                       *
 *                    I M P O R T A N T   N O T E ! ! !                  *
 *                                                                       *
 *      All scoring methods need to check if UsePatterns is set.         *
 *      If it is set, the scoring method must only score the first       *
 *      full alignment.  eg pos. 1- pos. 1, 2-2, 3-3, etc...             *
 *                                                                       *
 *************************************************************************
 *************************************************************************
 */
/*
compile block ! ! ! @ # $ % ^ &* * ( ) _ + @^ *#&^$ 
stuck...
Need to decide if the repeats is needed for the scoring of patterns
  If not, continue with the current stuff, otherwise, back to the chalkboard

Note, currently planning on using the optional integer field in the sequence
structure to carry the information on where the pattern matching is to start.
There will also be another parameter to tell if there is pattern matching
involved (note, this invalidates the above comment note).
*/
/* most recent approach is to tell the scoring method if it is a */
/* patterns search.  If it is then we treat it as a score all repeats. */
/* The method to try for not scoring repeats is this.  Make the */
/* scoring return info about the score if it was a patterns search. */
/* This way, since the patterns score is only a single score, the */
/* calling function can keep track of the scores and enter only the */
/* best one into the scores list.  Not that it will need to keep track */
/* of the other information needed also. */


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
 *     Boolean pattern_search: whether or not to use the pattern search stuff
 *   Error codes: none
 */

void default_scoring_method(seq, matrix, frame, repeats, search_type,
			    pattern_search) 
     Sequence *seq;
     Matrix *matrix;
     int frame;
     Boolean repeats;
     int search_type;
     Boolean pattern_search;
{
  int aa, mat_pos, seq_pos;
  int matrix_length, seq_length;
  int scan_pos;
  double max_aa_score, max_score;
  double seq_score=0, max_seq_score=0;
  int best_position=0;
  Boolean seen_pattern;

  /* get the lengths, gets rid of the indirrection */
  matrix_length = matrix->width;
  seq_length    = seq->length;
  
  /* Make sure these lengths are > 0 */
  if (matrix_length <= 0) {
    sprintf(ErrorBuffer,
	    "default_scoring_method(): matrix %s has length <= 0",
		matrix->number);
    ErrorReport(PROGRAM_ERR_LVL);
    sprintf(ErrorBuffer,
	    "default_scoring_method(): Not scoring the sequence.\n");
    ErrorReport(PROGRAM_ERR_LVL);
    return;
  }
  if (seq_length <= 0) {
    sprintf(ErrorBuffer,
      "default_scoring_method(): sequence %s has length (amino acids) <= 0",
		seq->name);
    ErrorReport(PROGRAM_ERR_LVL);
    sprintf(ErrorBuffer,
	    "default_scoring_method(): Not scoring the sequence.\n");
    ErrorReport(PROGRAM_ERR_LVL);
    return;
  }

  /* compute the maximum possible score for this matrix */
  /* this is used in standardizing the scores in a BLOCK search with a  */
  /* block without a percentile ranking (only calulate max_score if this */
  /* is the case) */
  max_score = 0.0;
  if ((search_type == SEARCH_TYPE_BLOCK) &&
      (matrix->percentile <= 0)) {
    for (mat_pos=0; mat_pos < matrix_length; mat_pos++) {
      max_aa_score = 0.0;
      for (aa=0; aa < MATRIX_AA_WIDTH; aa++) {
	if (matrix->weights[aa][mat_pos] > max_aa_score) {
	  max_aa_score = matrix->weights[aa][mat_pos];
	}
      }
      max_score = max_score + max_aa_score;
    }
  }

  

  /* 
   * score the matrix vs the sequence for each pairwise line up 
   */
  max_seq_score = 0.0;
  
  /* for every alignment of the matrix and sequence score the alignment. */
  /* there are (seq_length + matrix_length  - 1) different alignments. */
  /* Note: seq_pos is the relative position of the 1st matrix column to the */
  /*       1st sequence column */
  /* Note: the indexing is shifted to make the calculation of scan_pos */
  /*       (see below) easier/faster */
  seen_pattern = FALSE;
  for (seq_pos= -matrix_length+1; seq_pos < seq_length; seq_pos++) 
  {
    if (pattern_search) {
      /* This is a pattern!  Only score the first full alignment. */
      if (seen_pattern) {
	break; /* break out of the for loop.  I know it is a hack, but this *
		* does the least damage to the other logic. */
      }
      else {
	seen_pattern = TRUE; /* stop after this one time through the loop */
	seq_pos = seq->undefined;
      }
    }
    

    /* 
     * score this alignment
     */
    seq_score = 0.0;

    /* for each position in the matrix add the weight to the total score */
    /* Note: mat_pos is the current column of the matrix */
    for (mat_pos=0; mat_pos < matrix_length; mat_pos++) 
    {
      /* make sure the sequence and the matrix overlap at this point */
      /* Note: scan_pos is where the current matrix column is in the */
      /*       sequence */
      scan_pos = seq_pos + mat_pos;
      if ((scan_pos >= 0) && (scan_pos < seq_length)) { /* if in the seq */
	seq_score += matrix->weights[seq->sequence[scan_pos]][mat_pos];
      }
      else {		/* not in the sequence, score as gap */
	seq_score += matrix->weights[aa_atob['-']][mat_pos]; 
      }
      
    } /* end of for mat_pos: score this alignment */
    Alignments_Done += 1.0;
    
    /*
     * adjust the score if it is a SEARCH_TYPE_BLOCK search 
     */
    if (search_type == SEARCH_TYPE_BLOCK) {
      if (matrix->percentile > 0) {
	seq_score = (int) ((seq_score/matrix->percentile ) 
			   * 1000 );
      }
      else {
	seq_score = (int) (((seq_score/max_score) 
			    * log10((double) matrix->width))
			   * 1000 );
      }
    } /* end adjust score */

    if (seq_score > max_seq_score) {
	max_seq_score = seq_score;
	best_position = seq_pos;
    }

    /*   Save all scores for repeats */   
    if (repeats && (!SavedScoresFlag || seq_score >= MIN_SAVE_SCORE) )
    {
         enter_score((int)seq_score, seq_pos, frame, matrix, seq);
    }
  } /* end of for seq_pos */ 
    
  /*
   * if repeats were not allowed then save the best score seen
   */
 if (!repeats && (!SavedScoresFlag || seq_score >= MIN_SAVE_SCORE) )
 {
    enter_score((int)max_seq_score, best_position, frame, matrix, seq);
 }
  
}  /*  end of default_scoring_method */


/*
 * enter_score
 *   takes the data for a score entry and enters it into the list if it is
 *   ok to enter it. Note that not entering the score if it will be deleted 
 *   anyway saves time with the allocating and deallocating of score records.
 *   Parameters:
 *     int score:          the score
 *     int position:       the position
 *     int frame:          the frame we are looking at 
 *                             (0 = aa seq., 1, 2, 3, -1, -2, -3 na seqs.)
 *     Matrix *matrix:     the matrix
 *     Sequence *sequence: the sequence
 *   Error codes: none
 */

static void enter_score(score, position, frame, 
			matrix, sequence)
     int score;
     int position;
     int frame;
     Matrix *matrix;
     Sequence *sequence;
{  
  int result;

  /* increase the number of scores done */
  Scores_Done += 1.0;

  /* enter into the histogram array if keeping track */
  if (DoHistogram) {
    histogram[min((max(0, (int)score/HIST_BOX_SIZE)), (NUM_HIST_COLS-1))]++;
  }

  /* if score is worth entering into the list, enter it */
  /*   if the number of scores in the list is not yet to the max, we */
  /*   still need to enter the score, no matter what */
  if (!(score < MinScoreOfList) || 
      ((NumberToReport >= 0) && 
       (NumInSL(Scores) < NumberToReport))) {
    /* save score into score_list */
    result = 
      insert_in_score_list(make_score(score, position, frame, 
				      matrix, sequence));
    if (result == FALSE) {
      sprintf(ErrorBuffer, 
	      "enter_score(): Error placing matrix score into list\n");
      ErrorReport(PROGRAM_ERR_LVL);
    }
  }
}



/*
 * make_score
 *   allocates space for the structure that contains the score, sequence
 *   position matrix, and sequence information to be added to the Scores list.
 *   Parameters:
 *     int score:          the score
 *     int position:       the position
 *     int frame:          the frame we are looking at 
 *                             (0 = aa seq., 1, 2, 3, -1, -2, -3 na seqs.)
 *     Matrix *matrix:     the matrix
 *     Sequence *sequence: the sequence
 *   Return codes: returns a ScoreListEntry to be used later.
 *   Error codes:
 */

static Score *make_score(score, position, frame, matrix, sequence)
     int score;
     int position;
     int frame;
     Matrix *matrix;
     Sequence *sequence;
{  
  int block_length;		/* the length of the block */
  int seq_length;		/* the length of the sequence */
  int block_pos;		/* the position in the block */
  int scan_pos;			/* where we are related to block[0] */
  int seq_loop;			/* looping index through all the sequences */
				/* in the block*/
  char *consensus;
  Residue block_res, seq_res;	/* temporary places for the residues that */
				/* are being compared */
  Block *block;

  Score *new_score;

  CheckMem(
	   new_score = (Score *) malloc(sizeof(struct score_struct))
	   );

  new_score->score = score;
  new_score->strength = matrix->strength;
  new_score->position = position;
  new_score->seq_length = sequence->length;
  new_score->frame = frame;
  strncpy(new_score->sequence_number, sequence->name, NUMBER_WIDTH);
  strncpy(new_score->sequence_desc, sequence->info, DESC_WIDTH);
  if (matrix->block != NULL) {
    /* there is a block to use */
    strncpy(new_score->matrix_number, matrix->block->number, NUMBER_WIDTH);
    strncpy(new_score->matrix_desc, matrix->block->de, DESC_WIDTH);
  }
  else {			
    /* no block to use */
    strncpy(new_score->matrix_number, matrix->number, NUMBER_WIDTH);
    strncpy(new_score->matrix_desc, matrix->de, DESC_WIDTH);
  }

  /* if there is a block to look at and use to make the consensus sequence */
  /* then make the consensus sequence */
  /* otherwise, just use the sequence segment as the consensus    JGH*/
  consensus = new_score->consensus;
  seq_length = sequence->length;
  
  if (matrix->block != NULL) {
    block = matrix->block;

    block_length = block->width;
    
    /* for each position in the block figure out what case to use in the */
    /* consensus sequence to be printed out */
    /* Note: block_pos is the current column of the block */
    for (block_pos=0; block_pos < min(block_length, DESC_WIDTH); block_pos++) {

      /* make sure the sequence and the block overlap at this point */
      /* Note: scan_pos is where the current block column is in the */
      /*       sequence */
      scan_pos = position + block_pos;
      if (scan_pos < 0) {	/* before the sequence */
	consensus[block_pos] = ' ';
      }
      else if ((scan_pos >= 0) && (scan_pos < seq_length)) { /* in the seq */
	/*
	 * check for a consensus 
	 */
	consensus[block_pos] = tolower(aa_btoa[sequence->sequence[scan_pos]]);
	/* for every sequence in the block, check the residue.  if there */
	/* is one that is the same as the sequence, then change the case */
	seq_res = sequence->sequence[scan_pos];
	for (seq_loop=0; seq_loop< block->num_sequences; seq_loop++) {
	  block_res = block->sequences[seq_loop].sequence[block_pos];
	  if (seq_res == block_res) {
	    consensus[block_pos] = toupper(aa_btoa[seq_res]);
	    break;
	  }
	  /* if it is an ambiguity code, change the case */
	  else if (((seq_res == aa_atob['B']) &&
		    ((block_res == aa_atob['D']) || 
		     (block_res == aa_atob['N']))) ||
		   ((seq_res == aa_atob['Z']) &&
		    ((block_res == aa_atob['E']) || 
		     (block_res == aa_atob['Q']))) ||
		   ((block_res == aa_atob['B']) &&
		    ((seq_res == aa_atob['D']) || 
		     (seq_res == aa_atob['N']))) ||
		   ((block_res == aa_atob['Z']) &&
		    ((seq_res == aa_atob['E']) || 
		     (seq_res == aa_atob['Q'])))) {
	    consensus[block_pos] = toupper(aa_btoa[seq_res]);
	    break;
	  }
	}
	/* make sure that all X's are lower case */
	if (seq_res == aa_atob['X']) {
	  consensus[block_pos] = tolower(aa_btoa[seq_res]);
	}
      }
      else {			/* after the sequence */
	break;			/* don't need to keep going */
      }
    }
    consensus[block_pos] = '\0';
  }
  else { /* no block to make a consensus sequence from */
 
    block_length = matrix->width;
     
    /* for each position in the block figure out what case to use in the */
    /* consensus sequence to be printed out */
    /* Note: block_pos is the current column of the block */
    for (block_pos=0; block_pos < min(block_length, DESC_WIDTH); block_pos++) {
 
      /* make sure the sequence and the block overlap at this point */
      /* Note: scan_pos is where the current block column is in the */
      /*       sequence */
      scan_pos = position + block_pos;
      if (scan_pos < 0) {	/* before the sequence */
 	consensus[block_pos] = ' ';
      }
      else if ((scan_pos >= 0) && (scan_pos < seq_length)) { /* in the seq */
 	/*
 	 * copy the sequence as the "consensus"
 	 */
 	consensus[block_pos] = tolower(aa_btoa[sequence->sequence[scan_pos]]);
 	/* for every sequence in the block, check the residue.  if there */
 	/* is one that is the same as the sequence, then change the case */
      }
      else {			/* after the sequence */
 	break;			/* don't need to keep going */
      }
    }
    consensus[block_pos] = '\0';
  }

  return new_score;
}




/* 
 * print_histogram
 *   ouputs the histogram array from scoring.[ch] to the output stream ofp.
 *   Parameters:
 *     FILE *ofp: the output file pointer
 *   Return codes: none
 *   Error codes: none
 */

void print_histogram(ofp)
     FILE *ofp;
{
  int i, maxi;

  /*  Find highest non-empty cell   JGH */
  maxi = NUM_HIST_COLS - 1;
  while (histogram[maxi] == 0) {
    maxi--;
  }

  fprintf(ofp, "Histogram:\n");
/*  fprintf(ofp, "         < %-5d : %d\n", HIST_BOX_SIZE, histogram[0]);*/

  for (i=0; i<=maxi; i++) {
    fprintf(ofp, "   %5d - %-5d : %d\n", i*HIST_BOX_SIZE, 
	    (i+1)*HIST_BOX_SIZE-1,
	    histogram[i]);
  }

  fprintf(ofp, "   %5d -       : %d\n", 
	  (i)*HIST_BOX_SIZE, 
	  0);
  fprintf(ofp, "\n");
}





/* Change log information follows. */
/*
  Changes since version 3.3.1:
  10/ 2/99 Added includes for more standard C libraries.
  Changes since version 3.2.5:
  5/25/99 Scores_Done & Alignments_Done changed to double.
  2/3/99  default_scoring_method():
          If SavedScoresFlag, don't save scores < MIN_SAVE_SCORE
  Changes since version 3.1:
  1/16/97  Clarified error message. JGH
  Changes since version 3.0.0:
  3/25/96  Print out sequence name if error reading db.  JGH
*/
