/* (C) Copyright 1993, Fred Hutchinson Cancer Research Center */
/* Use, modification or distribution of these programs is subject to */
/* the terms of the non-commercial licensing agreement in license.h. */

/* alphabet.h: Definitions for the nucleotide and amino acid alphabets in */
/*             aabet.h and ntbet.h */
/* Modified by: Bill Alford */
/* Change log information is at the end of the file. */

#ifndef _ALPHABET_
#define _ALPHABET_

/*
	For better performance on most platforms, ALPHASIZE_MAX and ALPHAVAL_MAX
	should have values of 2**N and (2**N - 1), respectively, for some N.
	Choose an N that is just large enough for the job, to avoid needlessly
	thrashing in small CPU caches.
*/
#define ALPHASIZE_MAX	128	/* Max. no. of letters permitted in an alphabet */
#define ALPHAVAL_MAX	127	/* Max. numerical value permitted for a letter */


/* Ambiguous residues are defined using this structure */
typedef struct degen {
		char	residue;	/* ASCII letter for this residue */
		int		ndegen;		/* no. of residues this ASCII letter matches */
		char	*list;		/* null-terminated list of matching letters */
	} Degen, DegenPtr;

	/* Use DEGENLIST macro to more easily manage data initialization */
#define DEGENLIST(Chr, List)	{ Chr, (sizeof(List)-1), (List) }

#endif /* !_ALPHABET_ */




/* Change log information follows. */
/*
 * Revision 2.2010  1995/07/28  23:47:14  billa
 * 
 */

