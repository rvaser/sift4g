/* (C) Copyright 1993-7, Fred Hutchinson Cancer Research Center */
/* Use, modification or distribution of these programs is subject to */
/* the terms of the non-commercial licensing agreement in license.h. */

/* errors.c: error reporting functions */
/* Written by: Bill Alford */
/* Change log information is at the end of the file. */


/*	system headers not in global.h  */
/*	blimps library headers */
#include <global.h>

/*
 * Exported variables and data structures
 */

char ErrorBuffer[LARGE_BUFF_LENGTH];
int  ErrorLevelReport;


/*
 * Exported variables and data structures
 */

static char ErrorFile[SMALL_BUFF_LENGTH];

/*
 * Function definitions
 */

/*
 * set_error_file_name
 *   Sets the name of the error output file to use.
 *   Parameters:
 *     char *error_file: the name of the file.
 *   Error codes: None
 */

void set_error_file_name(error_file)
     char *error_file;
{
  strncpy(ErrorFile, error_file, SMALL_BUFF_LENGTH);
}

/*
 * ErrorReport
 *   ErrorReport prints out the error message in the ErrorBuffer to stderr and
 *   the ErrorFile if it can be opened.
 *   Note: The error message must be placed in ErrorBuffer before ErrorReport
 *         is called.
 *   Parameters:
 *     int err_level:  The error level reporting at.  If the error level is
 *                     >= FATAL_ERR_LVL, exit() is called with the error level.
 *   Error codes: None
 */

void ErrorReport(err_level)
     int err_level;
{
  FILE *efp; 
  Boolean dont_skip_error_file_report;
  
  /* check to see if we report */
  if ((err_level >= ErrorLevelReport) || (err_level == FATAL_ERR_LVL)) {

    /* open the error file for appending */
    efp = fopen(ErrorFile, "a");
    
    if (efp == NULL) {
      dont_skip_error_file_report = FALSE;
    }
    else {
      dont_skip_error_file_report = TRUE;
    }
  
    switch (err_level) {
    case INFO_ERR_LVL:
      fprintf(stderr, "Information: %s\n", ErrorBuffer);
      if (dont_skip_error_file_report) {
	fprintf(efp, "Information: %s\n", ErrorBuffer);
      }
      break;
    case WARNING_ERR_LVL:
      fprintf(stderr, "Warning: %s\n", ErrorBuffer);
      if (dont_skip_error_file_report) {
	fprintf(efp, "Warning: %s\n", ErrorBuffer);
      }
      break;
    case SERIOUS_ERR_LVL:
      fprintf(stderr, "Serious: %s\n", ErrorBuffer);
      if (dont_skip_error_file_report) {
	fprintf(efp, "Serious: %s\n", ErrorBuffer);
      }
      break;
    case PROGRAM_ERR_LVL:
      fprintf(stderr, "Program Error: %s\n", ErrorBuffer);
      if (dont_skip_error_file_report) {
	fprintf(efp, "Program Error: %s\n", ErrorBuffer);
      }
      break;
    case FATAL_ERR_LVL:
    default:
      fprintf(stderr, "Fatal Error: %s\n", ErrorBuffer);
      if (dont_skip_error_file_report) {
	fprintf(efp, "Fatal Error: %s\n", ErrorBuffer);
      }
      if (dont_skip_error_file_report) {
	fclose(efp);		/* close the file */
      }
      
      /* exit(????); choose a good error level to return */
      exit(FATAL_ERR_LVL);	/* closes all the open files */
    }
    
    ErrorBuffer[0] = '\0';	/* clear the string incase the caller does */
				/* not setup ErrorBuffer correctly for the */
				/* next call */
    if (dont_skip_error_file_report) {
      fclose(efp);		/* close the file */
    }
  }
}



/*
 * ABRT_signal_handler
 *   Catches the SIGABRT signal and announces the probable cause that was
 *   seen in testing.  It appears that when there is very little memory
 *   left free() has a hard time deallocating memory.  I'm guessing that it
 *   is when it is trying to put the free memory block into a list, the
 *   error that occurs is: "getfreehdr: out of memory".  The error seems to
 *   have only occured in the function call free() after many small blocks
 *   of memory have been freed.  It probably cannot handle as many free 
 *   blocks of memory as it is getting.
 *   NOTE: The memory recover functions "should" work.  The only
 *     problem appears to be in the function free().  The function is
 *     always called with a valid pointer to a busy block, allocated
 *     previously by a malloc function, when the error occurs.  And for
 *     about 10 calls before the error the pointer is OK too.  The
 *     problem I'm having with free() may be particular to this
 *     development system (SunOS Release 4.1.3) or I might have done
 *     something that causes this side effect.
 *   NOTE: The function recover_memory() in memory.c handles the recovery
 *     of memory when there is no more allocatable memory.  See the macro
 *     CheckMem(A) in memory.h to see how the function is called.
 *   Parameters: none
 *   Error codes: none, aborts the program.
 */

void ABRT_signal_handler()
{
  sprintf(ErrorBuffer,
	  "Caught the SIGABRT (6) signal error.");
  ErrorReport(PROGRAM_ERR_LVL);
  sprintf(ErrorBuffer,
	  "In testing, this signal occured when all the available");
  ErrorReport(PROGRAM_ERR_LVL);
  sprintf(ErrorBuffer,
	  "memory was filled and when trying to free memory to continue,");
  ErrorReport(PROGRAM_ERR_LVL);
  sprintf(ErrorBuffer,
	  "by using free(), after many small blocks had been freed.");
  ErrorReport(PROGRAM_ERR_LVL);
  sprintf(ErrorBuffer,
	  "If that is the case, try setting the number of scores to");
  ErrorReport(PROGRAM_ERR_LVL);
  sprintf(ErrorBuffer,
	  "report to a lower value than at the point memory was filled. ");
  ErrorReport(PROGRAM_ERR_LVL);
  sprintf(ErrorBuffer,
	  "Sorry for the inconvienience, but the memory recovery");
  ErrorReport(PROGRAM_ERR_LVL);
  sprintf(ErrorBuffer,
	  "functions \"should\" have worked.  ");
  ErrorReport(PROGRAM_ERR_LVL);
  sprintf(ErrorBuffer,
	  "See the errors.c source code for programming details.\n");
  ErrorReport(PROGRAM_ERR_LVL);

  sprintf(ErrorBuffer,
	  "Aborting due to the SIGABRT signal.\n");
  ErrorReport(FATAL_ERR_LVL);
}


/* Change log information follows. */
/* 
 * Revision 1.1  1993/08/17  04:18:16  billa
 * Added INFO_ERR_LVL for low level messages/errors and ErrorLevelReport
 * to control the level of error reporting.
 *
 * Revision 0.51  1993/07/23  17:40:18  billa
 * Changed the exit return number to FATAL_ERR_LVL.  Was empty.
 *
 * Revision 0.20  1993/06/23  16:46:49  billa
 * creation, error reporting to stderr and to the error file from
 * ErrorReport()
 *
 */



