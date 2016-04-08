/* (C) Copyright 1993-7, Fred Hutchinson Cancer Research Center */
/* Use, modification or distribution of these programs is subject to */
/* the terms of the non-commercial licensing agreement in license.h. */

/* lists.c: definitions and functions for various ordered lists */
/* Written by: Bill Alford */
/* Change log information is at the end of the file. */

/*	system headers not in global.h  */
/*	blimps library headers not in global.h */
#include <global.h>
#include <blocks.h>	/* includes sequences.h  */
#ifndef NO_MATRIX
#include <matrix.h>	/* includes pattern.h */
#endif
/*	headers in current directory */
#include "lists.h"
#include "blimps.h"
#ifndef NO_SCORES
#include "scores.h"
#endif


/*
 * Exported variables and data structures
 */

ScoreList    Scores;
ScoreList    PrintScores;
MatrixList   Matrices;
BlockList    Blocks;
SequenceList Sequences;

int MinScoreOfList = -1;	/* the minimum score of the score list */


/*
 * Local variables and data structures
 */


/*
 * Function definitions
 */


/* 
 * initialize_lists
 *   initializes the lists.
 *   Parameters: none
 *   Error codes:
 */

void initialize_lists()
{
#ifndef NO_SCORES
  Scores      = NewSL(score_comparison,     free_score, ALLOW_DUPLICATES);
  PrintScores = NewSL(neg_score_comparison, free_score, ALLOW_DUPLICATES);
#endif
#ifndef NO_MATRIX
  Matrices   = NewSL(matrix_comparison,    NULL, ALLOW_DUPLICATES);
#endif
  Blocks      = NewSL(block_comparison,     NULL, ALLOW_DUPLICATES);
  Sequences   = NewSL(sequence_comparison,  NULL, ALLOW_DUPLICATES);
}




/*
 * Score List functions
 *
 */
#ifndef NO_SCORES

/*
 * insert_in_score_list
 *   Inserts the score entry into the score list and makes sure that the 
 *   number of elements n the list are fewer than the NumberToReport if there
 *   is a limit.
 *   Parameters:
 *     Score score: the score entry.
 *   Return codes: the same as InsertSL.
 *   Error codes: none
 */

int insert_in_score_list(score)
     Score *score;
{
  int ret_value;

  ret_value = InsertSL(Scores, score);
  if (NumberToReport >= 0) {
    DoForSL(Scores, limit_Scores_list_size, NULL);
  }
  return ret_value;
}


/* 
 * limit_Scores_list_size
 *   Used in insert_in_score_list to limit the size of the score list.
 *   Parameters: 
 *     Score *score: the score list entry
 *     void *arg: unused, here to match the procedure call in DoForSL
 *   Return codes: SL_DELETE if the list is too big, SL_QUIT otherwise.
 *   Error codes: none
 */

int limit_Scores_list_size(score, arg)
     Score *score;
     void *arg;  /* not used */
{
/*>>> Seems to cause memory problems?
  if (NumInSL(Scores) > NumberToReport &&
      SavedScoresFlag && score->score > 1000 && NumberToReport < 4990)
  {		NumberToReport += 10;    }
*/

  if (NumInSL(Scores) > NumberToReport) {
    return SL_DELETE;
  }
  else {
    MinScoreOfList = score->score;
    return SL_QUIT; /* done limiting size, stop the DoForSL */
  }
}


/* 
 * enter_score_into_print_scores
 *   Places the score from Scores into PrintScores.  PrintScores is a reverse
 *   ordered list so that the printing comes out with the largest score
 *   first.
 *   Parameters:
 *     Score *score: the score list entry
 *     void *arg: unused, here to match the procedure call in DoForSL
 *   Return code: SL_CONTINUE always.
 *   Error codes: none
 */

int enter_score_into_print_scores(score, arg)
     Score *score;
     void *arg;
{
  InsertSL(PrintScores, score);
  return SL_CONTINUE;
}

#endif








/* Change log information follows. */
/* 
 Changes since version 3.2.5:
 3/3/99 Corrected spelling of Matrices.
        Note: SavedScoresFlag changes to limit_Scores_list_size() were removed,
	      see comments there.
 */
