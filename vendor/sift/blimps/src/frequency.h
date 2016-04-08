/* (C) Copyright 1993, Fred Hutchinson Cancer Research Center */
/* Use, modification or distribution of these programs is subject to */
/* the terms of the non-commercial licensing agreement in license.h. */

/* frequency.h: hardcoded amino acid frequency information for matricies */
/* Written by: Bill Alford */
/* Change log information is at the end of the file. */

#ifndef FREQUENCY_H_
#define FREQUENCY_H_

/* the default files with the amino acid frequencies */
#define LOCAL_AMINO_FREQUENCY_FILE   "default.amino.frq" 
#define LOCAL_CODON_FREQUENCY_FILE   "default.codon.frq" 


/*
 * Exported variables and data structures
 */


extern double frequency[MATRIX_AA_WIDTH];

extern double Codon_Usage[64];


/*
 * Exported functions
 */

/*
 * load_frequencies
 *   reads in the frequencies if there is a frequency file.
 *   Parameters:
 *     char *fname: The frequency file filename
 *   Return codes: TRUE if there was a frequency file, FALSE if not.
 *   Creates the global array frequency[]
 */

extern Boolean load_frequencies();


/*
 *  load_codons
 *    reads in the codon usages if there is a codon usage file.
 *    Parameters:
 *      FILE *ffp:  The codon usage file.
 *    Creates the global array Codon_Usage[]
 */

extern Boolean load_codons();

/*
 *  frq_qij
 *    Calls load_frequencies() and load_Qij to initialize the
 *    global variables frequency[], Qij, RTot
 */
extern void frq_qij();

#endif /*  FREQUENCY_H_ */
/* Change log information follows. */
/*
   Changes since 3.1:
   2/14/97 Added Codon_Usage[] and load_codons()
 * Revision 2.2011  1995/07/29  02:34:27  billa
 * Made load_frequencies more portable.
 *
 *
 */
