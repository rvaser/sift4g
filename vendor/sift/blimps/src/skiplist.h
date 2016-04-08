/* (C) Copyright 1993, Fred Hutchinson Cancer Research Center */
/* Use, modification or distribution of these programs is subject to */
/* the terms of the non-commercial licensing agreement in license.h. */

/* skiplist.h: This header file contains the definitions for use with the */
/*             generic SkipList package. */
/* Change log information is at the end of the file. */
/* most of this file is untouched from the distribution I started from, Bill */

/* This header file contains the definitions for use with the generic
 * SkipList package.
 *
 *      -- THIS CODE IS UNDER THE GNU COPYLEFT --
 *
 *    Dominic Giampaolo (nick@maxine.wpi.edu) 
 */

#ifndef SKIPLIST_H
#define SKIPLIST_H


/* RAND_MAX should be defined if you are using an ANSI compiler system,
 * but alas it isn't always.  You should define it to be the correct
 * value for whatever your library rand() function returns.
 *
 * Under unix (mach, bsd, etc), that's 2^31 - 1.  On my Amiga at home 
 * it's 2^15 - 1.  It would be wise to verify what your compiler uses
 * for RAND_MAX (the maximum value returned from rand()) because otherwise
 * the code will _not_ work.
 */
#ifndef RAND_MAX
#define RAND_MAX (0x7fffffff)
#endif


#define ALLOW_DUPLICATES  1   /* allow or disallow duplicates in a list */
#define NO_DUPLICATES	  0
#define DUPLICATE_ITEM	 -1   /* ret val from InsertSL if dups not allowed */


/* typedef's */
typedef struct SLNodeStruct *SLNode;

struct SLNodeStruct
{
  void	 *key;
  SLNode  forward[1]; /* variable sized array of forward pointers */
};

typedef struct _SkipList
{
  struct SLNodeStruct *header;	   /* pointer to header */

  int  (*compare)();
  void (*freeitem)();

  int flags;
  int level;			   /* max index+1 of the forward array */

  int count;                       /* number of elements in the list */
} *SkipList;



/* protos */
SkipList   NewSL();
void	   FreeSL();
int	   InsertSL();
int	   DeleteSL();		/* delete the matching key. */
				/* Only deletes the first one if duplicates. */
void	  *SearchSL();
void	   DoForSL();
void       DoForRangeSL();

int        NumInSL();

void      *Nth(); /* returns the Nth element in the list. */
		  /* Counting starts at 0. */

int LowerSavedNodesLevel();	/* function to lower the amount of nodes */
				/* saved in the SavedNodes list.  Returns */
				/* the new level of saving. */

/* These defines are to be used as return values from the function
 * you pass to DoForSL().  They can be or'ed together to do multiple
 * things (like delete a node and then quit going through the list).
 */
#define  SL_CONTINUE  0x00
#define  SL_DELETE    0x01
#define  SL_QUIT      0x02

#endif	/* SKIPLIST_H */


/* added by Bill Alford
 * Change log information follows. 
 * $Log: skiplist.h,v $
 * Revision 1.3  2011-09-28 22:10:55  gsims
 * *** empty log message ***
 *
 * Revision 1.2  2011-06-03 18:33:01  gsims
 * *** empty log message ***
 *
 * Revision 1.1  2011-05-24 16:32:32  gsims
 * *** empty log message ***
 *
 * Revision 1.1.1.1  2011-02-22 22:16:02  gsims
 * Initial Import of sift 4.0.4
 *
 * Revision 2.2010  1995/07/28  23:47:14  billa
 * Added new convert method and pattern matching and minor updates merged.
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
 * Revision 0.51  1993/07/02  00:55:15  billa
 * un-ANSIfied the prototypes.  Removed the argument declarations to be
 * backward compatible.
 *
 * Revision 0.50  1993/07/02  00:36:31  billa
 * Entry into the system.  Addition of skiplists to the system.
 *
 */






