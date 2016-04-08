/* (C) Copyright 1993-7, Fred Hutchinson Cancer Research Center */
/* Use, modification or distribution of these programs is subject to */
/* the terms of the non-commercial licensing agreement in license.h. */

/* blimps-mem.c: Memory management function for blimps. */
/* Written by: Bill Alford */
/* Change log information is at the end of the file. */


/*	system headers not in global.h */
/*	blimps library headers not in global.h */
#include <global.h>
/*	headers in current directory */
#include "blimps-mem.h"
#include "blimps.h"
#include "lists.h"


/*
 * Exported variables and data structures
 */

/*
 * Local variables and data structures
 */

static int CallsWhenToReduce  = 10;
static int NumberOfTimesCalled = 0;
#define DECREASE_AMOUNT (CallsWhenToReduce*.25 + 1) /* 10 7 5 3 2 1 0 */
#define REPORT_DECREASE_AMOUNT  1  /* amount to decrease NumberToReport by */
				   /* if there is not enough room */
				   /* NOTE: MUST be greater than zero */



/*
 * Function definitions
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

Boolean blimps_reclaim_space()
{
  /* note for later.  make a skiplist function to turn on severe space */
  /* saving.  In the function it will empty out the saved empty list */
  /* nodes and set a variable that will turn off the saving of the */
  /* nodes.  This should free up a _small_ amount of space.  Might not */
  /* be worth doing.  This skiplist function would be called from here */
  /* when the space is very tight. */

  int list_size;
  int prev_report_size;
  int ret_val;

  /* increase the count of the number of times this has been called */
  if (!(NumberOfTimesCalled<0)) {
    NumberOfTimesCalled++;
  }

  /* get the size of the Scores list for future use */
  list_size = NumInSL(Scores);
  prev_report_size = NumberToReport;

  if ((NumberToReport == 0) &&
      (list_size == 0)) {
    /* There is nothing we can do to save memory.  NumberToReport is */
    /* non-zero once the runs start.  If the list_size is zero, we deleted */
    /* everything we could. */
    sprintf(ErrorBuffer,
	    "Unable to recover memory.\n");
    ErrorReport(SERIOUS_ERR_LVL);
    return FALSE;
  }
  else if (NumberToReport < 0) { /* trying to save all, but can't */
    /* reset the NumberToReport */
    if (SearchType == SEARCH_TYPE_MATRIX) { /* blk vs sequences */
      NumberToReport = MAX_DEFAULT_MATRIX_SEARCH_NUM_TO_REPORT;
    }
    else if (SearchType == SEARCH_TYPE_BLOCK) {	/* seq vs blocks */
      NumberToReport = MAX_DEFAULT_BLOCK_SEARCH_NUM_TO_REPORT;
    }
    else {
      NumberToReport = max(MAX_DEFAULT_BLOCK_SEARCH_NUM_TO_REPORT, 
			   MAX_DEFAULT_MATRIX_SEARCH_NUM_TO_REPORT);
    }

    NumberToReport = list_size - REPORT_DECREASE_AMOUNT;
    if ( NumberToReport < 0 ) {
      NumberToReport = 1;
    }

    /* report the error */
    sprintf(ErrorBuffer,
	    "Decreasing the NumberToReport from saving all (%d currently) to %d.",
	    list_size, NumberToReport);
    ErrorReport(INFO_ERR_LVL);
    sprintf(ErrorBuffer,
	    "Ran out of memory (1), fewer scores will be reported.\n");
    ErrorReport(WARNING_ERR_LVL);

    /* free up some space with the new NumberToReport */
    DoForSL(Scores, limit_Scores_list_size, NULL);

    return TRUE;
  }
  else if (list_size <= NumberToReport) { /* can not report as many */
					  /* scores as we wanted to */
    /* reset the NumberToReport */
    NumberToReport = list_size - REPORT_DECREASE_AMOUNT;
    if ( NumberToReport < 0 ) {
      /* if we have not totally decreased the saved nodes to none, lower */
      /* the amount of nodes to save and report one, otherwise report zero */
      if (NumberOfTimesCalled >= 0) { /* NumberOfTimesCalled <0 if there is */
				      /* no more decreasing possible */
	sprintf(ErrorBuffer,
		"Decreasing the number of saved nodes in the list.");
	ErrorReport(INFO_ERR_LVL);
	ret_val = LowerSavedNodesLevel();
	if (ret_val < 0) {
	  sprintf(ErrorBuffer,
		  "The saved nodes level of the list has been reduced to %d.", 
		  ret_val);
	  ErrorReport(INFO_ERR_LVL);
	  sprintf(ErrorBuffer,
		  "The list saved nodes cannot be reduced any more.\n");
	  ErrorReport(WARNING_ERR_LVL);
	  NumberOfTimesCalled = -1; /* set up to stop calling if nothing */
				    /* else can be done */
	}
	else {
	  sprintf(ErrorBuffer,
		  "The saved nodes level of the list has been reduced to %d.\n", 
		  ret_val);
	  ErrorReport(INFO_ERR_LVL);
	}
	NumberToReport = 1;
      }
      else {			/* Can't decrease the saved nodes anymore */
	NumberToReport = 0;	/* Don't want to overshoot zero and get into */
				/* an infinite loop */
				/* NOTE: This will allow it to keep going */
				/*       but no scores will be reported */
      }
    }

    /* if this function has been called too often, reduce the level of */
    /* saving in the skiplists */
    if ((NumberOfTimesCalled >= 0) &&
	(NumberOfTimesCalled > CallsWhenToReduce)) {
      NumberOfTimesCalled = 0;
      CallsWhenToReduce = CallsWhenToReduce - DECREASE_AMOUNT;
      sprintf(ErrorBuffer,
	      "Decreasing the number of saved nodes in the list.");
      ErrorReport(INFO_ERR_LVL);
      ret_val = LowerSavedNodesLevel();
      if (ret_val < 0) {
	sprintf(ErrorBuffer,
		"The saved nodes level of the list has been reduced to %d.", 
		ret_val);
	ErrorReport(INFO_ERR_LVL);
	sprintf(ErrorBuffer,
		"The list saved nodes cannot be reduced any more.\n");
	ErrorReport(WARNING_ERR_LVL);
	NumberOfTimesCalled = -1; /* set up to stop calling if nothing */
				  /* else can be done */
      }
      else {
	sprintf(ErrorBuffer,
		"The saved nodes level of the list has been reduced to %d.\n", 
		ret_val);
	ErrorReport(INFO_ERR_LVL);
      }
    }

    /* report the error */
    sprintf(ErrorBuffer,
	    "Decreasing the NumberToReport from %d to %d.",
	    prev_report_size, NumberToReport);
    ErrorReport(INFO_ERR_LVL);
    sprintf(ErrorBuffer,
	    "Ran out of memory (2), fewer scores will be reported.\n");
    ErrorReport(WARNING_ERR_LVL);

    /* free up some space with the new NumberToReport */
    DoForSL(Scores, limit_Scores_list_size, NULL);

    return TRUE;
  }
  else {			/* list_size > NumberToReport, need to */
				/* decrease the size of the list. */
    /* Should never get to this point.  The list size should always be */
    /* kept <= the NumberToReport. */

    /* report the error */
    sprintf(ErrorBuffer,
	    "blimps_reclaim_space(): NumInSL(Scores) is greater than NumberToReport.");
    ErrorReport(PROGRAM_ERR_LVL);
    sprintf(ErrorBuffer,
	    "                 NumInSL(Scores) is supposed to be maintained");
    ErrorReport(PROGRAM_ERR_LVL);
    sprintf(ErrorBuffer,
	    "                 less than or equal to the NumberToReport.");
    ErrorReport(PROGRAM_ERR_LVL);
    sprintf(ErrorBuffer,
	    "                 NumberToReport = %d, NumInSL(Scores) = %d.",
	    NumberToReport, NumInSL(Scores));
    ErrorReport(PROGRAM_ERR_LVL);
    sprintf(ErrorBuffer,
	    "                 Ran out of memory, removing the excess scores");
    ErrorReport(PROGRAM_ERR_LVL);
    sprintf(ErrorBuffer,
	    "                 in the Scores list.\n");
    ErrorReport(PROGRAM_ERR_LVL);

    /* free up some space with the new NumberToReport */
    DoForSL(Scores, limit_Scores_list_size, NULL);

    return TRUE;
  }
}







/* Change log information follows. */
/* 
 * Revision 1.1000  1993/08/27  17:00:11  billa
 * Creation.  Added the function reclaim_space() and the macro CheckMem(A)
 * to handle the situation when there is no more memory to allocate.
 *
 */

