/* (C) Copyright 1993-7, Fred Hutchinson Cancer Research Center */
/* Use, modification or distribution of these programs is subject to */
/* the terms of the non-commercial licensing agreement in license.h. */

/* memory.h: Memory management functions and macros. */
/* Written by: Bill Alford */
/* Change log information is at the end of the file. */

#ifndef MEMORY_H_
#define MEMORY_H_

/*
 * Exported variables and data structures
 */

#define REPORT_DECREASE_AMOUNT  1  /* amount to decrease NumberToReport by */
				   /* if there is not enough room */
				   /* NOTE: MUST be greater than zero */


/*
 * Exported macros
 */

/********************************************************************
 * This macro is supposed to enclose malloc, calloc, and realloc 
 * assignments.  For example:                             
 * CheckMem(                                              
 *          new_block = (Block *) malloc(sizeof(Block))   
 *          );                                            
 * 
 * Note that if using around a realloc you should be assigning the new 
 * memory pointer to a temporary variable so that you do not loose the 
 * original pointer to memory if the CheckMem macro has to loop.
 * For example:                          
 *  Residue *tmp_ptr;  
 *  ...  
 *  CheckMem(                                                        
 *           tmp_ptr = (Residue *) realloc(seq->sequence,            
 *                                         seq->max_length *   
 *                                         sizeof(Residue))    
 *           );                                                      
 *  seq->sequence = tmp_ptr;
 * 
 ********************************************************************/

#ifndef NO_RECLAIM
#define CheckMem(A)		           \
{                                          \
  while(((A) == NULL) &&                   \
	((reclaim_space()) ?               \
	 TRUE :                            \
	 (sprintf(ErrorBuffer,             \
		 "Unable to recover enough memory to continue.  Aborting.\n"),\
	  ErrorReport(FATAL_ERR_LVL),      \
	  FALSE)                           \
	 )                                 \
	);                                 \
}
#else
#define CheckMem(A)                        \
{                                          \
  while(((A) == NULL) &&                   \
        (sprintf(ErrorBuffer,              \
                 "Unable to allocate memory.  Aborting.\n"),\
         ErrorReport(FATAL_ERR_LVL),       \
         FALSE)                            \
        );                                 \
}
#endif /* NO_RECLAIM */


/*
 * Exported functions
 */

/*
 * reclaim_space
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
 *         calls to reclaim_space may be needed before the free lists
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

#ifndef NO_RECLAIM
extern Boolean reclaim_space();
#endif /* NO_RECLAIM */


/*
 * init_reclaim_space
 *   Sets up the reclaim_space function to call the passed function.
 *   Parameters:
 *     Boolean (*rec_func)(): The reclaiming function.
 *   Return codes: none.
 *   Error codes: none.
 */

extern void init_reclaim_space();


#endif /*  MEMORY_H_ */




/* Change log information follows. */
/* 
 * Revision 2.1021  1994/05/09  20:22:18  billa
 * Made the reclaim_space() function call the programmer supplied memory
 * recovery function rather than using the blimps specific code.
 *
 * Revision 2.1017  1994/04/27  02:29:29  billa
 * Removed RECLAIM_SPACE define.
 *
 * Revision 1.1000  1993/08/27  16:56:34  billa
 * Creation.  Added the function reclaim_space() and the macro CheckMem(A)
 * to handle the situation when there is no more memory to allocate.
 *
 */



