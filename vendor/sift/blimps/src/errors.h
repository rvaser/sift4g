/* (C) Copyright 1993-7, Fred Hutchinson Cancer Research Center */
/* Use, modification or distribution of these programs is subject to */
/* the terms of the non-commercial licensing agreement in license.h. */

/* errors.h: error reporting functions */
/* Written by: Bill Alford */
/* Change log information is at the end of the file. */

#ifndef ERRORS_H_
#define ERRORS_H_

/*
 * Exported variables and data structures
 */

#define INFO_ERR_LVL    1	/* General information. */
#define WARNING_ERR_LVL 2	/* Probable user error. */
				/* Generally recoverable. */
#define SERIOUS_ERR_LVL 3	/* Probable user error. */
				/* Difficult to recover from. */
#define PROGRAM_ERR_LVL 4	/* Programing error, problem should have */
				/* been caught before it arose. */
				/* Recoverability varies. */
#define FATAL_ERR_LVL   5	/* Various causes, user or programmer. */
				/* Unable to recover from. */
				/* Exits program with FATAL_ERR_LVL value. */


extern char ErrorBuffer[LARGE_BUFF_LENGTH];
extern int  ErrorLevelReport;


/*
 * Exported functions
 */

/*
 * set_error_file_name
 *   Sets the name of the error output file to use.
 *   Parameters:
 *     char *error_file: the name of the file.
 *   Error codes: None
 */

extern void set_error_file_name();


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

extern void ErrorReport();


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

extern void ABRT_signal_handler();

#endif /*  ERRORS_H_ */

/* Change log information follows. */
/*
 * Revision 0.20  1993/06/23  16:46:49  billa
 * creation, error reporting to stderr and to the error file from
 * ErrorReport()
 *
 */



