/* (C) Copyright 1993, Fred Hutchinson Cancer Research Center */
/* Use, modification or distribution of these programs is subject to */
/* the terms of the non-commercial licensing agreement in license.h. */

/* pattern.h:  */
/* Written by: Bill Alford */
/* Change log information is at the end of the file. */

#ifndef PATTERN_H_
#define PATTERN_H_

/*
 * Exported variables and data structures
 */

extern Boolean UsePatterns;  /* initially FALSE */

typedef struct pattern_residue_struct PatternResidue;
struct pattern_residue_struct {
  int offset;			/* the offset from the current position,
				   the current position is the previous 
				   PatternResidue or the starting location */
  int num_residues;		/* the number of residues */
  char *residues;		/* the residues */
  PatternResidue *next; /* the next pattern res. in the list */
};


struct pattern_struct {
  int beg_offset;		/* the offset from here to the beginning of
				   the sequence */
  int num_residues;		/* the number of sub residues */
  PatternResidue *pat;		/* the list of residues in the pattern */
};
typedef struct pattern_struct Pattern;




extern void scan_patterns();

extern Boolean pattern_matches();

extern int residue_compare_function();

#endif /*  PATTERN_H_ */

/* Change log information follows.
 * $Log: pattern.h,v $
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
 */

