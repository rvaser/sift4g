/* (C) Copyright 1993, Fred Hutchinson Cancer Research Center */
/* Use, modification or distribution of these programs is subject to */
/* the terms of the non-commercial licensing agreement in license.h. */

/* skiplist.c: This file contains a heavily hacked and generalized version */
/*             of the Example skiplist code distributed on mimsy.cs.umd.edu. */
/* Change log information is at the end of the file. */
/* most of this file is untouched from the distribution I started from, Bill */

/* This file contains a heavily hacked and generalized version of the
   Example skiplist code distributed on mimsy.cs.umd.edu.

 Here is a short excerpt from the original comment :

	 Example of Skip List source code for C :

	   Skip Lists are a probabilistic alternative to balanced trees,
	 as described in the June 1990 issue of CACM and were invented by
	 William Pugh in 1987.

 These are my additions :

     This file contains my (Dominic Giampaolo's) heavily hacked version
   of skip lists.  These work on any arbitrary data by using callback
   functions which you supply (at list creation time) to do the data
   comparisons.  You could instantly use this package to implement a
   symbol table for a compiler which would be blazingly fast and
   require zippo effort on your part.

     I've changed the function names (not to protect the innocent, but
   to make them easier to read :) and changed the data structures a bit.
   I've ansi'fied the code, added prototypes, and changed all those ugly
   do/while's to for loops.  I also removed the dependance on those silly
   NIL items at the end of the list (it just checks for regular NULL
   pointers instead).  Additionally, the code is more easily reentrant now,
   and doesn't depend on any global variables.  The code is quite a bit
   different looking than it originally was, but the underlying algorithims
   (of course) remain unchanged.

              -- THIS CODE IS UNDER THE GNU COPYLEFT --


      Dominic Giampaolo (nick@maxine.wpi.edu)
      Aug, 1991
*/
/*	system headers not in global.h*/
#include <assert.h>
#include <time.h>
/*	blimps library headers */
#include <global.h>
#include <skiplist.h>

/* define's */
#ifndef TRUE
#define TRUE  1
#endif
#ifndef FALSE
#define FALSE 0
#endif

#define MaxNumberOfLevels 16
#define MaxLevel (MaxNumberOfLevels-1)
/* #define NewNodeOfLevel(x) (SLNode)malloc(sizeof(struct SLNodeStruct)+(x)*sizeof(SLNode *)) */ /* commented out by Bill Alford */


/* private proto */
static int RandomLevelSL();

/* added by Bill Alford */
static void FreeSLNode(); /* node, level */
static SLNode NewNodeOfLevel(); /* level */


/* save structure */
#define SLNODES_TO_SAVE_PER_LEVEL_AT_START 10
int SLNodesToSave = SLNODES_TO_SAVE_PER_LEVEL_AT_START;

struct {
  int num_saved;
  int num_to_save;
  SLNode saved_nodes;
} SavedNodes[MaxNumberOfLevels];

int InitializedSaveList = FALSE;


/* functions added by Bill Alford */

static SLNode NewNodeOfLevel(level) 
     int level;
{
  register SLNode tmp;

  level--; /* change indexing to start at 0 */

  if (!SavedNodes[level].num_saved) {
    /* This prints out the level of the node being created. */
    /* It is useful for debugging. */
    /*
    fprintf(stderr, "making node level: %d\n", level);
    */
    CheckMem(
	     tmp = (SLNode) malloc(sizeof(struct SLNodeStruct) +
				   (level+1)*sizeof(SLNode *))
	     );
    return tmp;	       /* needed to use the tmp to work around the macro */
  }
  else {
    tmp = SavedNodes[level].saved_nodes;
    SavedNodes[level].saved_nodes = tmp->forward[0];
    SavedNodes[level].num_saved--;
    return tmp;
  }
}

static void FreeSLNode(node, level)
     SLNode node;
     int level;
{
  level--; /* change indexing to start at 0 */

  if (SavedNodes[level].num_saved >= SavedNodes[level].num_to_save) {
    free(node);
  }
  else {
    /* save the node */
    /* NOTE: This behaves like a stack.  It might be better to use another */
    /*       style of save list, maybe like FIFO.  This might help */
    /*       _slightly_ in the fragmentation problem. */
    node->forward[0] = SavedNodes[level].saved_nodes;
    SavedNodes[level].saved_nodes = node;
    SavedNodes[level].num_saved++;
  }
  /* This prints out the number of nodes in each list. */
  /* It is useful for debugging. */
  /*
  {
    int i;
    for (i=0; i<MaxNumberOfLevels; i++) {
      fprintf(stderr, "%2d ", SavedNodes[i].num_saved);
    }
    fprintf(stderr, "\n");
  }
  */
}

static void SetSavedNodesLevel(nodes_to_save)
     int nodes_to_save;
{
  int i;

  for (i=0; i<MaxNumberOfLevels; i++) {
    SavedNodes[i].num_to_save = 
      (i==0 ? 
       nodes_to_save :
       (SavedNodes[i-1].num_to_save -
	((int) ((double)(SavedNodes[i-1].num_to_save*i) * .25 + 0.5))));
  }
}

int LowerSavedNodesLevel()
{
  SLNode tmp;
  int i;

  if (SavedNodes[0].num_to_save <= 0) {
    return -1;			/* unable to lower the save level any more */
  }
  SetSavedNodesLevel(SavedNodes[0].num_to_save - 1);
  for (i=0; i<MaxNumberOfLevels; i++) {
    while (SavedNodes[i].num_saved > SavedNodes[i].num_to_save) {
      tmp = SavedNodes[i].saved_nodes;
      SavedNodes[i].saved_nodes = tmp->forward[0];
      SavedNodes[i].num_saved--;
      free(tmp);
    }
  }
  return SavedNodes[0].num_to_save; /* the new save level */
}

/* end added by Bill Alford */

/* functions */
SkipList  NewSL(compare, freeitem, flags)
     int (*compare)();
     void (*freeitem)();
     int flags;
{
  SkipList l;
  int i;

  /* added by Bill Alford */ 
  /* initialize the saved nodes list */
  if (!InitializedSaveList) {
    for (i=0; i<MaxNumberOfLevels; i++) {
      SavedNodes[i].num_saved = 0;
      SavedNodes[i].saved_nodes = NULL;
      SetSavedNodesLevel(SLNodesToSave);
    }
    InitializedSaveList = TRUE;
  }
  /* end added by Bill Alford */

  if (compare == NULL)    /* need at least a compare function... */
    return NULL;

  CheckMem(
	   l = (SkipList)malloc(sizeof(struct _SkipList))
	   );
  if (l == NULL)
    return NULL;

  l->level = 1;
  l->header = NewNodeOfLevel(MaxNumberOfLevels);
  if (l->header == NULL)
    { free(l); return NULL; }

  for(i=0; i < MaxNumberOfLevels; i++)
    l->header->forward[i] = NULL;
  l->header->key = NULL;		 /* just to be sure */

  srand(time(NULL) | 0x01);   /* seed with an odd number */

  l->compare	 = compare;
  l->freeitem	 = freeitem;
  l->flags	 = flags;

  l->count = 0;

  return(l);
}


void FreeSL(l)
     SkipList l;
{
  register SLNode p,q;
  void (*freeitem)() = l->freeitem;

  if (l == NULL || l->header == NULL)
    return;

  p = l->header;	   /* free header node first, because it doesn't */
  q = p->forward[0];	   /* have a real key to it			 */
  free(p);
  p = q;

  while (p != NULL)
   {
     q = p->forward[0];
     if (freeitem)
       (*freeitem)(p->key);
     free(p);
     p = q;
   }

  free(l);
}




/*
 *   This RandomLevelSL function generates a very good representation of
 *   p=.25 (or p=.5, etc).  The number of nodes of each level works out
 *   to be very very close to what they should be.  I didn't check it
 *   statistically,  but on large data sets, I imagine it's +/- 5% of what
 *   it should be.  This P value is good for lists of up to 64K elements.
 *   
 *   For more info about the P value, see the papers by Mr. Pugh (available
 *   in postscript from mimsy.umd.edu).
 */
#define P_50   (RAND_MAX / 2)     /* p value of .50   */
#define P_25   (RAND_MAX / 4)     /* p value of .25   */
#define P_125  (RAND_MAX / 8)     /* p value of .125  */

static int RandomLevelSL(l)
     SkipList l;
{
  register int level = 0;

  while(rand() < P_25)
    {
      level++;
    }

  return (level > MaxLevel ? MaxLevel : level);
}


int InsertSL(l, key)
     SkipList l;
     void *key;
{
  register int i,k;
  SLNode update[MaxNumberOfLevels];
  register SLNode p,q=0;
  int (*compare)() = l->compare;

  p = l->header;
  
  for(k = l->level-1; k >= 0; k--)
   {
     while((q = p->forward[k]) && (*compare)(q->key, key) < 0)
	p = q;

     update[k] = p;
   }

  if ((l->flags & ALLOW_DUPLICATES) == FALSE) /* if no duplicates allowed */
   if (q && (*compare)(q->key, key) == 0)     /* item is a duplicate */
     {
       return DUPLICATE_ITEM;
     }


  k = RandomLevelSL(l);
  if (k >= l->level)
    {
      l->level++;
      k = l->level - 1;
      update[k] = l->header;
    }

  q = NewNodeOfLevel(k+1);

  if (q == NULL)
    return FALSE;

  l->count++;             /* update the number of nodes in the list */
  
  q->key = key;
  for(i=0; i < k; i++)
    q->forward[i] = NULL;

  for(; k >= 0; k--)
    {
      p = update[k];
      q->forward[k] = p->forward[k];
      p->forward[k] = q;
    }

  return TRUE;
}



int DeleteSL(l, key)
     SkipList l;
     void *key;
{
  register int k,m;
  SLNode update[MaxNumberOfLevels];
  register SLNode p,q;
  int  (*compare)()  = l->compare;
  void (*freeitem)() = l->freeitem;

  p = l->header;

  for(k=l->level-1; k >= 0; k--)
   {
     while((q = p->forward[k]) && (*compare)(q->key, key) < 0)
	p = q;

     update[k] = p;
   }
  q = p->forward[0];

  if (q && (*compare)(q->key, key) == 0)
    {
      m = l->level - 1;
      for(k=0; k <= m; k++)
	{
	  p = update[k];
	  if (p == NULL || p->forward[k] != q)
	    break;
	  p->forward[k] = q->forward[k];
	}

      l->count--;
      
      if (freeitem)
	(*freeitem)(q->key);
      
      /* added by Bill Alford */
      FreeSLNode(q, k);
      /* end added by Bill Alford */

      m = l->level - 1;
      while(l->header->forward[m] == NULL && m > 0)
	m--;

      l->level = m + 1;
      return TRUE;
    }
  else
    return FALSE;
}



void *SearchSL(l, key)
     SkipList l; 
     void *key;
{
  register int k;
  register SLNode p=0,q=0;
  int (*compare)() = l->compare;

  p = l->header;
  
  for(k=l->level-1; k >= 0; k--)
   {
     while((q = p->forward[k]) && (*compare)(q->key, key) < 0)
	p = q;
   }

  if (q == NULL || (*compare)(q->key, key) != 0)
    return NULL;

  return q->key;
}


void DoForSL(l, function, arg)
     SkipList l;
     int (*function)();
     void *arg;
{
  register SLNode p,q, fix;
  register int k, ret;
  SLNode save[MaxNumberOfLevels], who[MaxNumberOfLevels];
  void (*freeitem)() = l->freeitem;
  

  if (l == NULL || l->header == NULL || function == NULL)
    return;

  /* this is to make compilers happy.  Since the important values are copied */
  /* into the array below */
  for (k=0; k < MaxNumberOfLevels; k++)	{	/*JGH*/
    save[k] = who[k] = NULL; 
  }

  p = l->header;	   /* skip header node because it isn't a real node */

  /* Save the initial header info
   */
  for(k=0; k < l->level; k++)
   {
     save[k] = p->forward[k];
     who[k]  = p;
   }

  p = p->forward[0];      /* skip to the first data node */

  while (p != NULL)
   {
     q = p->forward[0];
     ret = (*function)(p->key, arg);

     if (ret & SL_DELETE)
      {
	k = 0;
	while(save[k] == p)
	 {
	   fix = who[k];
	   fix->forward[k] = p->forward[k];
	   save[k] = p->forward[k];
	   k++;
	 }

	l->count--;         /* decrement the count of items */
	
	if (freeitem)
	  (*freeitem)(p->key, arg);
	/* added by Bill Alford */
	FreeSLNode(p, k);
	/* end added by Bill Alford */
      }
     else
      {
	k = 0;
	while(save[k] == p)
	 {
	   save[k] = p->forward[k];
	   who[k]  = p;
	   k++;
	 }
      }
     
     if (ret & SL_QUIT)
       break;

     p = q;    /* advance to the next one */
   }
}


void DoForRangeSL(l, key, compare, func, arg)
     SkipList l;
     void *key;
     int (*compare)();
     int (*func)();
     void *arg;
{
  register int k;
  SLNode update[MaxNumberOfLevels];
  register SLNode p,q;
  void (*freeitem)() = l->freeitem;
  int  ret;

  p = l->header;

  for(k=l->level-1; k >= 0; k--)
   {
     while((q = p->forward[k]) && (*compare)(q->key, key) < 0)
	p = q;

     update[k] = p;
   }
  p = p->forward[0];

  if (p == NULL || (*compare)(p->key, key) != 0)   /* then nothing matched */
    return;

  do
   {
     q = p->forward[0];            /* save next pointer */
     ret = (*func)(p->key, arg);
     
     if (ret & SL_DELETE)
      {
	for(k=0; k < l->level && update[k] && update[k]->forward[k] == p; k++)
	  update[k]->forward[k] = p->forward[k];

	l->count--;         /* decrement the count of items */
	
	if (freeitem)
	  (*freeitem)(p->key, arg);
	/* added by Bill Alford */
	FreeSLNode(p, k);
	/* end added by Bill Alford */
      }
     
     if (ret & SL_QUIT)
       break;

     p = q;    /* advance to the next one */
   }
  while(p != NULL && (*compare)(p->key, key) == 0);

}


int NumInSL(l)
     SkipList l;
{
  return l->count;
}


/* 
 * Nth
 *
 * returns the Nth element in the list.  Counting starts at 0.
 *
 */

void *
Nth(l, n)
     SkipList l;
     int n;
{
  int x;
  SLNode p;

  /* range check */
  if ((NumInSL(l) < n) || (n < 0)) {
    return NULL;
  }
  
  p = l->header; /* skip header node because it isn't a real node */
  p = p->forward[0];      /* skip to the first data node */

  for (x=0; x<n; x++) {
    p = p->forward[0];
  }
  
  assert(p->key != NULL);
  return p->key;
}


/* added by Bill Alford */
/* Change log information follows. */
/*
   06/23/94 Added Nth()   BJA
   10/4/94 Initialized pointers in DoForSLC()   jgh 

 * $Log: sl.c,v $
 * Revision 1.2  2011-06-03 18:33:01  gsims
 * *** empty log message ***
 *
 * Revision 1.1  2011-05-24 16:32:30  gsims
 * *** empty log message ***
 *
 * Revision 1.1.1.1  2011-02-22 22:16:01  gsims
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
 * Revision 0.52  1993/07/16  20:10:05  billa
 * Added a save list for deleted nodes so that the lists will not thrash
 * memory mallocing and freeing memory.
 *
 * Revision 0.51  1993/07/02  01:02:28  billa
 * un-ANSIfied the function argument decarations to be backward compatible.
 *
 * Revision 0.50  1993/07/02  00:36:31  billa
 * Entry into the system.  Addition of skiplists to the system.
 *
 */






