/* (C) Copyright 1993-9, Fred Hutchinson Cancer Research Center */
/* Use, modification or distribution of these programs is subject to */
/* the terms of the non-commercial licensing agreement in license.h. */

/* memory.c: Memory management functions. */
/* Written by: Bill Alford */
/* Change log information is at the end of the file. */

/*	system headers not in global.h  */
/*	blimps library headers not in global.h */
#include <global.h>

/*
 * Exported variables and data structures
 */

/*
 * Local variables and data structures
 */

static Boolean (*RecFunc)() = NULL;

/*
 * Function definitions
 */

/*
 * reclaim_space
 *   Tries to reclaim some of the allocated memory.  It tries to
 *   reclaim the memory the program used using the function specified in 
 *   init_reclaim_space().
 *   NOTE: The reclaiming function must have the same return and error codes
 *         as below!
 *   Parameters: none
 *   Return codes: TRUE if it was able to get free some space, FALSE
 *                 if not.
 *   Error codes:  FALSE if it unable to get space.
 */

Boolean reclaim_space()
{
  if (RecFunc == NULL) {
    sprintf(ErrorBuffer, "Reclaim space function is not defined.  No memory will be reclaimed.\n");
    ErrorReport(WARNING_ERR_LVL);
    return FALSE;
  }
  else {
    return (*RecFunc)();
  }
}
  


/*
 * init_reclaim_space
 *   Sets up the reclaim_space function to call the passed function.
 *   Parameters:
 *     Boolean (*rec_func)(): The reclaiming function.
 *   Return codes: none.
 *   Error codes: none.
 */

void init_reclaim_space(rec_func) 
     Boolean (*rec_func)();
{
  if (rec_func == NULL) {
    sprintf(ErrorBuffer, "Reclaim space function is not defined.  No memory will be reclaimed.\n");
    ErrorReport(WARNING_ERR_LVL);
  }
  RecFunc = rec_func;
}





/* Change log information follows. */
/* 
 * Revision 2.2010  1995/07/28  23:47:14  billa
 * Added new convert method and pattern matching and minor updates merged.
 *
 */

