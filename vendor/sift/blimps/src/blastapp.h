/* (C) Copyright 1993, Fred Hutchinson Cancer Research Center */
/* Use, modification or distribution of these programs is subject to */
/* the terms of the non-commercial licensing agreement in license.h. */

/* blastapp.h: set up defines for aabet.h, alphabet.h, and gcode.h */
/* Modified by: Bill Alford */
/* Change log information is at the end of the file. */

/**************************************************************************
*                                                                         *
*              National Center for Biotechnology Information              *
*       Bldg. 38A, NIH,  8600 Rockville Pike,  Bethesda, MD  20894        *
*                                                                         *
*                             COPYRIGHT NOTICE                            *
*                                                                         *
* This software/database is categorized as "United States Government      *
* Work" under the terms of the United States Copyright Act.  It was       *
* produced as part of the author's official duties as a Government        *
* employee and thus can not be copyrighted.  This software/database is    *
* freely available to the public for use without a copyright notice.      *
* Restrictions can not be placed on its present or future use.            *
*                                                                         *
* Although all reasonable efforts have been taken to ensure the accuracy  *
* and reliability of the software and data, the National Library of       *
* Medicine (NLM) and the U.S. Government do not and can not warrant the   *
* performance or results that may be obtained by using this software,     *
* data, or derivative works thereof.  The NLM and the U.S. Government     *
* disclaim any and all warranties, expressed or implied, as to the        *
* performance, merchantability or fitness for any particular purpose or   *
* use.                                                                    *
*                                                                         *
* In any work or product derived from this material, proper attribution   *
* of the author(s) as the source of the software or data would be         *
* appreciated.                                                            *
*                                                                         *

Karlin, Samuel and Stephen F. Altschul (1990).  Methods for assessing
the statistical significance of molecular sequence features by using
general scoring schemes.  Proc. Natl. Acad. Sci. USA 87:2264-2268.

Altschul, Stephen F., Warren Gish, Webb Miller, Eugene W. Myers, and
David J. Lipman (1990).  A basic local alignment search tool.
J. Mol. Biol. 215:403-410.
**************************************************************************/
/* (For a better look, set tabstops every 4 columns) */
#ifndef __BLAST_H__
#define __BLAST_H__
/*
#include <signal.h>
#include <gish.h>
*/
#include <limits.h> /* defines CHAR_BIT for aabet.h and ntbet.h */

/* EXTERN should be defined in the main program module only */
#ifndef EXTERN
#define EXTERN	extern
#ifndef INITIAL
#define INITIAL(x)
#endif
#else
#ifndef INITIAL
#define INITIAL(x) = x
#endif
#ifndef INIT
#define INIT
#endif
#endif

/* much deleted here */

/* added */
#define PROTO(x) () 


#endif /* __BLAST_H__ */




/* Change log information follows.
 * $Log: blastapp.h,v $
 * Revision 1.1  2011-05-24 16:32:31  gsims
 * *** empty log message ***
 *
 * Revision 1.2  2011-02-23 04:43:46  gsims
 * *** empty log message ***
 *
 * Revision 1.1.1.1  2011-02-22 22:16:02  gsims
 * Initial Import of sift 4.0.4
 *
 * Revision 2.2010  1995/07/28  23:47:14  billa
 * Added new convert method and pattern matching and minor updates merged.
 *
 * Revision 2.2001  1994/05/18  19:17:57  billa
 * Removed global.c and moved Buffer variable to strutil.[ch].
 *
 * Revision 2.2000  1994/05/10  18:39:13  billa
 * System version 2.2 A [2.2000]
 * Goofed on dev. version.
 *
 * Revision 2.1100  1994/05/10  18:27:07  billa
 * System version 2.2 A [2.1100].
 * Created blimps library and fixed up minor bugs.  No major changes to the
 * program.
 *
 * Revision 2.1025  1994/05/10  00:51:47  billa
 * Added SCCS patterns to trick the what command into printing out the RCS
 * identifying patterns.
 *
 * Revision 2.1020  1994/04/27  21:59:11  billa
 * After library creation.  In testing phase.
 *
 * Revision 2.1015  1994/04/26  03:15:34  billa
 * Pre-library checkin.
 *
 * Revision 2.1011  1994/04/26  03:14:16  billa
 * Pre-library checkin.
 *
 * Revision 2.1010  1994/04/26  00:28:18  billa
 * Cleaned up warnings from gcc switches.
 *
 * Revision 2.1000  1994/03/09  22:23:59  billa
 * System version 2.1 A [2.1000]
 * Added reading and scoring of precomputed site specific scoring matricies.
 * Reconfigured files.[ch] into files.[ch] and config.[ch].
 *
 * Revision 2.0  1993/12/06  19:31:18  billa
 * System version 2.0 A [2.0000].
 * Changed source to remove some warnings.  Added sequence weights to the
 * blocks.  Added a conversion method to use the weights.  Made the method
 * the default.  Added some ifndef's in preparation of creating a library
 * of functions.  Added the ability to choose the genetic code to use for
 * translation.  Added the ability to read a flat file as a sequence if the
 * other types sequences fail to match.  Added some defines for strdup().
 *
 * Revision 1.1120  1993/12/02  03:44:22  billa
 * Development version 1.1120.  Right after copyright/license info added.
 *
 * Revision 1.1101  1993/12/02  03:37:12  billa
 * Added the copyright notice to the top of the file.
 *
 * Revision 1.1100  1993/09/09  17:23:19  billa
 * System version 1.1 B [1.1100].
 *
 * Revision 1.1000  1993/08/27  17:01:31  billa
 * System version 1.1 A, development version 1.1000.
 * Added the function reclaim_space() and the macro CheckMem(A) to handle
 * the situation when there is no more memory to allocate.
 *
 * Revision 1.0  1993/08/14  00:50:52  billa
 * System version 1.00
 *
 * Revision 0.110  1993/08/10  02:50:53  billa
 * System version update.  Added convert.[ch], scoring.[ch] and version.[ch].
 *
 * Revision 0.105  1993/08/04  19:53:40  billa
 * System version 0.105.  Post file structure reorganization.
 *
 * Revision 0.100  1993/08/03  00:48:46  billa
 * System version 0.100.  Pre file structure reorganization.
 *
 * Revision 0.11  1993/07/16  22:16:00  billa
 * Added RCS information.  Added PROTO definition.
 * 
 */






