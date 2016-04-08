/* (C) Copyright 1993-7, Fred Hutchinson Cancer Research Center */
/* Use, modification or distribution of these programs is subject to */
/* the terms of the non-commercial licensing agreement in license.h. */

/* lists.h: definitions and functions for various ordered lists */
/* Written by: Bill Alford */
/* Change log information is at the end of the file. */

#ifndef LISTS_H_
#define LISTS_H_

#include <skiplist.h>

/*
 * Exported variables and data structures
 */


typedef SkipList ScoreList;
typedef SkipList MatrixList;
typedef SkipList BlockList;
typedef SkipList SequenceList;


extern ScoreList    Scores;
extern ScoreList    PrintScores;
extern MatrixList   Matrices;
extern BlockList    Blocks;
extern SequenceList Sequences;
extern Boolean SavedScoresFlag;

extern int MinScoreOfList;	/* the minimum score of the score list */
				/* didn't make a function to do this */
				/* because it would be worse than just */
				/* having a global variable */



/*
 * Exported functions
 */


/* 
 * initialize_lists
 *   initializes the lists.
 *   Parameters: none
 *   Error codes:
 */

extern void initialize_lists();

#ifndef NO_SCORES
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

extern int enter_score_into_print_scores();

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

extern int insert_in_score_list();

/* 
 * limit_Scores_list_size
 *   Used in insert_in_score_list to limit the size of the score list.
 *   Parameters: 
 *     Score *score: the score list entry
 *     void *arg: unused, here to match the procedure call in DoForSL
 *   Return codes: SL_DELETE if the list is too big, SL_QUIT otherwise.
 *   Error codes: none
 */

int limit_Scores_list_size();

#endif



#endif /*  LISTS_H_ */

/* Change log information follows. */
/* 
  Changes since version 3.2.5:
  1/23/99 Added SavedScoresFlag, used by lists.c, blimps.c, scoring.c, config.c
*/
