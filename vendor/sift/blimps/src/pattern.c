/* (C) Copyright 1993, Fred Hutchinson Cancer Research Center */
/* Use, modification or distribution of these programs is subject to */
/* the terms of the non-commercial licensing agreement in license.h. */

/* pattern.c: */
/* Written by: Bill Alford */
/* Change log information is at the end of the file. */

/*	system headers not in global.h   */
#include <assert.h>
/*	blimps library headers  */
#include <global.h>
#include <files.h>
#include <blocks.h>	/* includes sequences.h, output.h */
#include <matrix.h>	/* includes pattern.h */
#include <residues.h>
#include <skiplist.h>


/*
 * Exported variables and data structures
 */

Boolean UsePatterns = FALSE;

/*
 * Local variables and data structures
 */

SkipList Patterns;
Boolean ListInitialized = FALSE;

#define PATTERN_ARRAY_START 30
#define PATTERN_ARRAY_INC 20

/*
 * Function definitions
 */

static Boolean headers_match();
static Boolean search_for_pattern();
static Pattern *scan_pattern();
static void find_residues();
static PatternResidue *parse_residue();
/* Never Used */
/* static void print_patterns(); */
/* static void print_pattern();  */
/* static void print_residue(); */
/* Never used */
/* static int print_residue_for_list(); */



/*
 * scan_patterns
 *
 * scans in the pattern entries from the pattern file for the given
 * matrix and puts them in the matrix.
 *
 */

void
scan_patterns(matrix)
     Matrix *matrix;
{
  int max_patterns_for_array = 0;
  int num_patterns_seen = 0;
  Boolean cont;
  Pattern **patterns;
  Pattern **tmp_ptr;
  SkipList residue_list;

  char *pat;

  FILE *pat_file;

  if (! UsePatterns) {
    matrix->patterns = NULL;
    return;
  }

  /* is the pattern already here? */
  if (matrix->patterns != NULL) {
    /* if so, then do nothing */
    return;
  }
  else {
    /* otherwise, we need to read in the patterns */

    /* allocate space for the array of pointers to patterns */
    CheckMem(
	     patterns = (Pattern **) malloc(sizeof(Pattern *) * PATTERN_ARRAY_START)
	     );

    max_patterns_for_array = PATTERN_ARRAY_START;
    num_patterns_seen = 0;

    /* prepare the residue list for storing residues. */
    /* this is done outside of the call to scan_pattern() so as to cut */
    /* down on initialization/cleanup time */
    /* Note that Patterns is global to save some startup/cleanup time */
    if (! ListInitialized) {
      Patterns = NewSL(residue_compare_function, NULL, ALLOW_DUPLICATES);
      ListInitialized = TRUE;
    }
    residue_list = Patterns;

/*printf("Searching for: %s\n", matrix->number);*/

    /* find the pattern entry that matches the PSSM */
    if (! headers_match(matrix)) {
/*printf("Pattern not found\n");*/
      rewind_file(PATTERN_FILES);

      if (! search_for_pattern(matrix)) {
	sprintf(ErrorBuffer,
		"Unable to find pattern for PSSM %s\n", matrix->number);
	ErrorReport(WARNING_ERR_LVL);

	matrix->patterns = NULL;
	return;
      }
    }

    /* file should be positioned somewhere in the pattern header */
    cont = FALSE;
    pat_file = get_file(PATTERN_FILES);
    /* there should be something there regardless */
    /* probably should mention bad file format */
    assert(pat_file != NULL);

    /* scan until after the PA line */
    while (fgets(Buffer, LARGE_BUFF_LENGTH, pat_file) &&
	   !((Buffer[0] == 'P') && (Buffer[1] == 'A')));

    /* read in each pattern line until "//" is seen */

    cont = TRUE;    do {
      /* read in a line */
      fgets(Buffer, LARGE_BUFF_LENGTH, pat_file);

      /* if it is "//" then we are done */
      if ((Buffer[0] == '/') && (Buffer[1] == '/')) {
	cont = FALSE;
      }
      else {
	/* shift to the beginning of the pattern */
	pat = get_token(Buffer);

	/* put a null at the end of the pattern */
	remove_trailing_whitespace(pat);

	/* increase the size if need be */
	/* (adding one early to allow room for the NULL ptr */
	if (num_patterns_seen >= (max_patterns_for_array-1)) {
	  max_patterns_for_array += PATTERN_ARRAY_INC;
	  CheckMem(
		   tmp_ptr = (Pattern **) realloc(patterns,
						  max_patterns_for_array *
						  sizeof(Pattern *))
		   );
	  patterns = tmp_ptr;
	}

	/* make null for safety */
	patterns[num_patterns_seen] = NULL;
	patterns[num_patterns_seen+1] = NULL;

	/* scan in this pattern and place it in the array */
	patterns[num_patterns_seen] = scan_pattern(pat, residue_list);

	/* only move ahead if something was scanned */
	if (patterns[num_patterns_seen] != NULL) {
	  num_patterns_seen++;
	}
      }
    } while (cont);

    /* make the last position NULL */
    patterns[num_patterns_seen] = NULL;

    /* the patterns have been read in and we can update the matrix */
    matrix->patterns = patterns;

/*print_patterns(matrix->patterns);*/

    assert(NumInSL(residue_list) == 0);

  }
}

Boolean
pattern_matches(seq, compare_start, pattern)
     Sequence *seq;
     int compare_start;
     Pattern *pattern;
{
  Residue *sequence;
  int seq_length;
  PatternResidue *p;
  Boolean ok=FALSE;
  int x, num, offset;

  sequence = seq->sequence;
  seq_length = seq->length;


  /* optimize for the common case of one residue */
  if (pattern->num_residues == 1) {
    /* get first residue */
    p = pattern->pat;

    /* check the residue(s) */

    offset = p->offset;

    /* Is this an error?  ok is never set */

    /* only check if it is ok if it is within the sequence */
    if (((compare_start+offset) < seq_length) &&
	((compare_start+offset) >= 0)) {

      num = p->num_residues;

      /* optimize for the common case by only checking the one */
      if (num == 1) {
	/* this saves an assignment and an additional compare and a
	   loop back*/
	if (aa_btoa[sequence[compare_start+offset]] ==
	    toupper(p->residues[0])) {
	  /* can return true here since there is only one residue */
	  return TRUE;
	}
      }
      else {
	for (x=0; !ok && x<num; x++) {
	  if (aa_btoa[sequence[compare_start+offset]] ==
	      toupper(p->residues[x])) {
	    /* can return true here since there is only one residue */
	    return TRUE;
	  }
	}
      }

      /* if we have reached this point nothing matched */
      return FALSE;
    }
  }
  /* otherwise do the general case */
  else {

    /* get first residue */
    p = pattern->pat;

    /* while the residue is there (the pointer != NULL) */
    while (p != NULL) {
      /* check the residue(s) */

      offset = p->offset;

      /* only check if it is ok if it is within the sequence */
      if (((compare_start+offset) < seq_length) &&
	  ((compare_start+offset) >= 0)) {
	ok = FALSE;

	num = p->num_residues;

	/* optimize for the common case by only checking the one */
	if (num == 1) {
	  /* this saves an assignment and an additional compare and a
	     loop back*/
	  if (aa_btoa[sequence[compare_start+offset]] ==
	      toupper(p->residues[0])) {
	    ok = TRUE;
	  }
	}
	else {
	  for (x=0; !ok && x<num; x++) {
	    if (aa_btoa[sequence[compare_start+offset]] ==
		toupper(p->residues[x])) {
	      ok = TRUE;
	    }
	  }
	}

	/* if it/they didn't match, return FALSE */
	if (! ok) {
	  return FALSE;
	}
      }

      /* get the next residue */
      p = p->next;
    }
  }
  return TRUE;
}


/*************************************************************************
 *************************************************************************
 *
 *             Local Functions
 *
 *************************************************************************
 *************************************************************************
 */



/*
 * residue_compare_function
 *
 * compare two residue structures.  This is for the skiplist ordering.
 * The shorter residues are placed first, then they are sorted by
 * alphabetical order.
 *
 */

int
residue_compare_function(res1, res2)
     PatternResidue *res1;
     PatternResidue *res2;
{
  if (res1->num_residues != res2->num_residues) {
    return res1->num_residues - res2->num_residues;
  }
  else {
    return strcmp(res1->residues, res2->residues);
  }
  assert(0);
}


/*
 * headers_match
 *
 * scan the pattern file from the current point and if the matrix and
 * header match (have the same base name) return true.
 *
 */

static Boolean
headers_match(matrix)
     Matrix *matrix;
{
  char *tmp;
  FILE *pat_file;


  /* foreach of the pattern files */
  while ((pat_file = get_file(PATTERN_FILES)) != NULL) {
    while (fgets(Buffer, LARGE_BUFF_LENGTH, pat_file) && !feof(pat_file)) {
      /* if we are past the header, stop and return FALSE */
      if ((Buffer[0] == '/') && (Buffer[1] == '/')) {
	return FALSE;
      }
      /* otherwise, if it is the AC line, check */
      else if ((Buffer[0] == 'A') && (Buffer[1] == 'C')) {
	/* get the number of the pattern.  Be sure to truncate the ';' */
	tmp = get_token(Buffer);
	tmp = get_token(NULL);
	if (tmp[strlen(tmp)-1] == ';') {
	  tmp[strlen(tmp)-1] = '\0';
	}

	/* must skip over the BL/MA stuff (so add two) */
	return (strcmp(tmp+2, (matrix->number)+2) == 0);
      }
      /* otherwise, skip it and go to the next line */
    }
  }
  /* Should never see this */
  return FALSE;
}

/*
 * search_for_pattern
 *
 * scan from the current file position in the pattern file searching
 * for the pattern.  If none found, return FALSE, otherwise return
 * TRUE.
 *
 */

static Boolean
search_for_pattern(matrix)
     Matrix *matrix;
{
  FILE *pat_file;
  char *tmp;

  /* foreach of the pattern files */
  while ((pat_file = get_file(PATTERN_FILES)) != NULL) {
    while (fgets(Buffer, LARGE_BUFF_LENGTH, pat_file) && !feof(pat_file)) {
      if ((Buffer[0] == 'A') && (Buffer[1] == 'C')) {
	/* get the number of the pattern.  Be sure to truncate the ';' */
	tmp = get_token(Buffer);
	tmp = get_token(NULL);
	if (tmp[strlen(tmp)-1] == ';') {
	  tmp[strlen(tmp)-1] = '\0';
	}

	if ((strcmp(tmp+2, (matrix->number)+2) == 0)) {
	  return TRUE;
	}
      }
    }
  }

  /* reached the end of the file, so return FALSE */
  return FALSE;
}





/*
 * scan_pattern
 *
 * scan the current pattern (the string is only stuff in the pattern)
 * and return a pointer to the pattern data structure
 *
 * NULL is returned if there are no residues for the pattern
 *
 */

static Pattern *
scan_pattern(pat, residue_list)
     char *pat;
     SkipList residue_list; /* declared externally for reuse to save time */
				/* it should start and end empty */
{
  Pattern *pattern;
  PatternResidue *p, *s;
  int start;

  assert(NumInSL(residue_list) == 0);

  CheckMem(
	   pattern = (Pattern *) malloc(sizeof(Pattern))
	   );

  /* get all the residues in the pattern and place them in the ordered list */
  find_residues(pat, residue_list);
  /* the list is now full of residues */

  /* prints the contents of the ordered residues list */
  /* DoForSL(residue_list, print_residue_for_list, NULL); */

  if (NumInSL(residue_list) == 0) {
    free(pattern);
    return NULL;
  }

  pattern->num_residues = NumInSL(residue_list);

  p = Nth(residue_list, 0);
  DeleteSL(residue_list, p);

  start = p->offset;
  p->offset = 0;
  pattern->pat = p;
  pattern->beg_offset = start;

  while (NumInSL(residue_list) != 0) {
    s = Nth(residue_list, 0);
    DeleteSL(residue_list, s);
    s->offset = s->offset - start;
    p->next = s;
    p = s;
  }

  p->next = NULL; /* just being safe */

  assert(NumInSL(residue_list) == 0);

  return pattern;
}



static void
find_residues(pat, residue_list)
     char *pat;
     SkipList residue_list;	/* to fill with the residues */
{
  int offset = 0;
  PatternResidue *patres;

  while (pat[0] != '\0') {
    if ((pat[0] != 'x') && (pat[0] != 'X')) {
      patres = parse_residue(offset, &pat);
      InsertSL(residue_list, patres);
    }
    offset++;
    pat++;
  }
}


/*
 *
 * parse_residue
 *
 * creates a pattern residue to be added to the pattern.  The residue
 * at the current possition is the one used.  THIS DOES NOT CHECK THAT
 * THE RESIDUE IS A VALID RESIDUE (eg 'X').
 *
 */

static PatternResidue *
parse_residue(offset, pat)
     int offset;
     char **pat;		/* this is so that the possition in the */
				/* pattern string is updated */
{
  PatternResidue *pr;
  int num = 0;
  char *tmp;

  if (**pat != '[') {		/* it is only one residue */

    /* declare space for the pattern residue and the character string */
    CheckMem(
	     pr = (PatternResidue *) malloc(sizeof(PatternResidue))
	     );
    CheckMem(
	     pr->residues = (char *) malloc(sizeof(char)*2)
	     );

    /* initialize the structure */
    pr->offset = offset;
    pr->num_residues = 1;
    pr->next = NULL;

    /* set the residues string */
    pr->residues[0] = **pat;
    pr->residues[1] = '\0';

    /* don't move the move the pattern pointer/marker ahead since it
       is done once in the calling function. */

    return pr;
  }
  else {			/* it is more than one residue */
    num = 0;
    tmp = *pat;
    tmp++;			/* move over the '[' */
    while (*tmp != ']') {
      num++;
      tmp++;
    }
    *tmp = '\0';		/* make place a null to end the residues */

    /* declare space for the pattern residue and the character string */
    CheckMem(
	     pr = (PatternResidue *) malloc(sizeof(PatternResidue))
	     );
    CheckMem(
	     pr->residues = (char *) malloc(sizeof(char)*(num+1))
	     );

    /* initialize the structure */
    pr->offset = offset;
    pr->num_residues = num;
    pr->next = NULL;

    /* set the residues string */
    strcpy(pr->residues, (*pat)+1);
    pr->residues[num] = '\0';

    /* move the pattern pointer/marker ahead, all but one.  It is
       moved ahead one position automatically */
    *pat = *pat + 1 + num;

    return pr;
  }

  assert(0);			/* shouldn't get here */

}

/* never used */
/*
static void
print_patterns(Pattern **patterns)
{
  int x;
  for (x = 0; patterns[x] != NULL; x++) {
    printf("%d: ", x);
    print_pattern(patterns[x]);
  }
}

static void
print_pattern(pattern)
     Pattern *pattern;
{
  if (pattern == NULL) {
    printf("(null)\n");
  }

  printf("%d: ", pattern->beg_offset);
  print_residue(pattern->pat);
}

*/
/*
void
static print_residue(pat)
     PatternResidue *pat;
{
  if (pat == NULL) {
    printf(";\n");
    return;
  }
  printf("(%d, [%s]) ", pat->offset, pat->residues);
  print_residue(pat->next);
}
*/
/*
int
static print_residue_for_list(pat, arg)
     PatternResidue *pat;
     void *arg;
{
  print_residue(pat);
  return SL_CONTINUE;
}
*/

/* Change log information follows.
 * $Log: pattern.c,v $
 * Revision 1.2  2011-09-28 22:10:55  gsims
 * *** empty log message ***
 *
 * Revision 1.1  2011-05-24 16:32:28  gsims
 * *** empty log message ***
 *
 * Revision 1.1.1.1  2011-02-22 22:16:00  gsims
 * Initial Import of sift 4.0.4
 *
 * Revision 2.2010  1995/07/28  23:47:14  billa
 * Added new convert method and pattern matching and minor updates merged.
 *
 */
