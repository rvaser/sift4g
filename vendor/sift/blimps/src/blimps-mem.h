/* (C) Copyright 1993-7, Fred Hutchinson Cancer Research Center */
/* Use, modification or distribution of these programs is subject to */
/* the terms of the non-commercial licensing agreement in license.h. */

/* blimps_mem.h: Memory management function for blimps. */
/* Written by: Bill Alford */
/* Change log information is at the end of the file. */

#ifndef BLIMPS_MEM_H_
#define BLIMPS_MEM_H_

/*
 * Exported variables and data structures
 */

/*
 * Exported macros
 */

/*
 * Exported functions
 */

/*
 * blimps_reclaim_space
 *   Tries to reclaim some of the allocated memory.  It tries to
 *   reclaim the memory from the various lists that are used by
 *   decreasing the elements that can be stored in the lists.
 *   NOTE: This function assumes that the searching has already begun.
 *         If that is not the case it will report that it was unable
 *         to get space because the list size and the NumberToReport
 *         will likely be zero.  This doesn't really matter much
 *         because if we have to reclaim space that early in the run, 
 *         there are problems.
 *   NOTE: Because of the SavedNodes list in the skiplist package many
 *         calls to blimps_reclaim_space may be needed before the free lists
 *         are filled and actual memory is freed.  This is not that
 *         big of a problem, considering that a balance will be
 *         reached with the nodes in the free lists.  The only time
 *         this would cause the Scores list size to reduce to zero is
 *         if there is so little memory that there can be only about
 *         20 scores in the list.  That amount of memory is about 100
 *         Kbytes. 
 *   Parameters: none
 *   Return codes: TRUE if it was able to get free some space, FALSE
 *                 if not.
 *   Error codes:  FALSE if it unable to get space.
 */


extern Boolean blimps_reclaim_space();


#endif /*  BLIMPS_MEM_H_ */




/* Change log information follows. */
/* 
 *
 * Revision 2.1017  1994/04/27  02:29:29  billa
 * Removed RECLAIM_SPACE define.
 *
 * Revision 2.1011  1994/04/26  02:39:09  billa
 * Removed BLIMPS dependent variables from global.[ch] and put in blimps.[ch].
 *
 * Revision 1.1000  1993/08/27  16:56:34  billa
 * Creation.  Added the function reclaim_space() and the macro CheckMem(A)
 * to handle the situation when there is no more memory to allocate.
 *
 */

