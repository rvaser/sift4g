/* (C) Copyright 1993-7, Fred Hutchinson Cancer Research Center */
/* Use, modification or distribution of these programs is subject to */
/* the terms of the non-commercial licensing agreement in license.h. */

/* frequency.c: amino acid frequency information  */
/* Written by: Bill Alford */
/* Change log information is at the end of the file. */

/*	system headers not in global.h */
/*	blimps library headers */
#include <global.h>
#include <blocks.h>	/* includes sequences.h for NA_SEQ */
#include <matrix.h>	/* includes pattern.h */
#include <frequency.h>
#include <convert.h>	/* For Qij stuff */
#include <files.h>
#include <residues.h>

/*
 * Exported variables and data structures
 */

double frequency[MATRIX_AA_WIDTH];
double Codon_Usage[64];


/*
 * Local variables and data structures
 */




/*
 * Function definitions
 */

static double read_freq();

/*
 * Frequency file related functions
 *
 *   Boolean load_aa_frequencies()
 *   static double read_freq(ffp)
 *
 */

/*
 * load_frequencies
 *   reads in the frequencies if there is a frequency file.
 *   puts them in global array frequency[]
 *   Parameters:
 *     char *fname: The frequency file filename
 *   Return codes: TRUE if there was a frequency file and the data was 
 *                 sucessfully read in, FALSE if not.
 *   Error codes: none
 */

Boolean load_frequencies(fname)
     char *fname;		/* the frequency file filename */
{
  int i;
  FILE *ffp;			/* the frequency file pointer */
  double J_weight, O_weight, U_weight;

  ffp = fopen(fname, "r");

  if (ffp == NULL) {
    sprintf(ErrorBuffer, "load_frequencies: Unable to open frequency file: %s", 
	    fname);
    ErrorReport(WARNING_ERR_LVL);
/*
    sprintf(ErrorBuffer, 
	    "load_frequencies: Setting all frequencies to one.");
    ErrorReport(SERIOUS_ERR_LVL);
*/
    for(i=0; i<MATRIX_AA_WIDTH; i++) {
      frequency[i] = (double) 1.0;
    }      
    return FALSE;
  }

#define NAR '|'   /* a char that get converted to AAID_NAR for the codes so */
		  /* that when reading in the frequencies entries for J, O, */
		  /* an U do not overwrite the entry for X.  Also if the */
		  /* idea of have non-code alpha characters being treated as */
		  /* X ever changes, there shouldn't be a change needed */
		  /* (assuming that J, O, U are weighted the same). */

  /* old style format */
  frequency[aa_atob['-']] = read_freq(ffp);
  frequency[aa_atob['A']] = read_freq(ffp);
  frequency[aa_atob['B']] = read_freq(ffp);
  frequency[aa_atob['C']] = read_freq(ffp);
  frequency[aa_atob['D']] = read_freq(ffp);
  frequency[aa_atob['E']] = read_freq(ffp);
  frequency[aa_atob['F']] = read_freq(ffp);
  frequency[aa_atob['G']] = read_freq(ffp);
  frequency[aa_atob['H']] = read_freq(ffp);
  frequency[aa_atob['I']] = read_freq(ffp);

  J_weight = read_freq(ffp);
  frequency[aa_atob[NAR]] = J_weight; /* 'J' */

  frequency[aa_atob['K']] = read_freq(ffp);
  frequency[aa_atob['L']] = read_freq(ffp);
  frequency[aa_atob['M']] = read_freq(ffp);
  frequency[aa_atob['N']] = read_freq(ffp);

  O_weight = read_freq(ffp);
  frequency[aa_atob[NAR]] = O_weight; /* 'O' */

  frequency[aa_atob['P']] = read_freq(ffp);
  frequency[aa_atob['Q']] = read_freq(ffp);
  frequency[aa_atob['R']] = read_freq(ffp);
  frequency[aa_atob['S']] = read_freq(ffp);
  frequency[aa_atob['T']] = read_freq(ffp);

  U_weight = read_freq(ffp);
  frequency[aa_atob[NAR]] = U_weight; /* 'U' */

  frequency[aa_atob['V']] = read_freq(ffp);
  frequency[aa_atob['W']] = read_freq(ffp);
  frequency[aa_atob['X']] = read_freq(ffp);
  frequency[aa_atob['Y']] = read_freq(ffp);
  frequency[aa_atob['Z']] = read_freq(ffp);
  frequency[aa_atob['*']] = read_freq(ffp);

#undef NAR


  /* check the frequency array for zero values */
  for (i=0; i<MATRIX_AA_WIDTH; i++) {
    if (frequency[i] <= 0.0) {
      /* announce warning */
      if (i != AAID_NAR) {
	sprintf(ErrorBuffer, 
		"load_frequencies: Read a frequency of value of %f for: %c\n", 
		frequency[i], 
		aa_btoa[i]);
      }
      else {
	sprintf(ErrorBuffer, 
		"load_frequencies: Read a frequency of value of %f for: not a residue\n",
		frequency[i]);
      }
      ErrorReport(INFO_ERR_LVL);
    }
  }


  /* Check to make sure the non-aa codes are the same.  These values will */
  /* be used for filling the position of the AAID_NAR spot, the last one */
  /* seen will be used. */
  if ( J_weight != O_weight ) {
    sprintf(ErrorBuffer, 
	    "load_frequencies: The non-amino acid codes J and O weights are different");
    ErrorReport(INFO_ERR_LVL);
  }
  if ( O_weight != U_weight ) {
    sprintf(ErrorBuffer, 
	    "load_frequencies: The non-amino acid codes O and U weights are different");
    ErrorReport(INFO_ERR_LVL);
  }
  if ( U_weight != J_weight ) {
    sprintf(ErrorBuffer, 
	    "load_frequencies: The non-amino acid codes U and J weights are different");
    ErrorReport(INFO_ERR_LVL);
  }
  if (( J_weight != O_weight ) ||
      ( O_weight != U_weight ) ||
      ( U_weight != J_weight )) {
    sprintf(ErrorBuffer, 
	    "load_frequencies: The value for the non-amino acid will be the last weight seen.");
    ErrorReport(INFO_ERR_LVL);
    sprintf(ErrorBuffer, 
	    "load_frequencies: Currently the program ignores non-residue characters so this");
    ErrorReport(INFO_ERR_LVL);
    sprintf(ErrorBuffer,
	    "field has no effect.\n");
    ErrorReport(INFO_ERR_LVL);
  }

  return TRUE;
}


/*
 * read_freq
 *   Reads a frequency value from the frequency file.  A value is the first
 *   number on a line.  If a line does not have a number first, the line is
 *   skipped.
 *   Parameters:
 *     FILE *ffp: the frequency file pointer.
 *   Return codes: The next frequency in the file.  The default value of zero
 *                 if there was an error in the reading.
 *   Error codes: none.  Can not determine between reading a zero and 
 *                returning it or if there was no number read.
 */

static double read_freq(ffp)
     FILE *ffp;
{
  double freq;

  freq = 0.0;

  while (fgets(Buffer, LARGE_BUFF_LENGTH, ffp) &&
	 sscanf(Buffer, "%lg", &freq) != 1);

  return freq;
}


/*
 * load_codons
 *   reads in the codons if there is a codon file
 *   Parameters:
 *     fcod codon file pointer
 *   Return codes: TRUE if there was a codon file and the data was 
 *                 sucessfully read in, FALSE if not.
 *   Error codes: none
 *   Loads array Codon_Usage[]
 */

Boolean load_codons(ffp)
	FILE *ffp;
{
  int i;

  if (ffp == NULL) {
    sprintf(ErrorBuffer, "Unable to open codon usage file");
    ErrorReport(SERIOUS_ERR_LVL);
    sprintf(ErrorBuffer, "load_codons: Setting all codon usages to 1.\n");
    ErrorReport(SERIOUS_ERR_LVL);
    for(i=0; i<64; i++) Codon_Usage[i] = (double) 1.0;
    return FALSE;
  }

  /*---- NOTE:  Assumes the codon usage file is in the correct order!  
     from AAA to UUU in ACGU order - see default.codon.use------*/
  for (i=0; i<64; i++) Codon_Usage[i] = read_freq(ffp);

  /* check the array for zero values */
  for (i=0; i<64; i++) {
    if (Codon_Usage[i] <= 0.0) {
      /* announce warning */
	sprintf(ErrorBuffer, 
		"load_codons: Read a codon usage value of %f for: %c\n", 
		Codon_Usage[i], aa_btoa[i]);
        ErrorReport(INFO_ERR_LVL);
    }
  }

  return TRUE;
} /* end of load_codons */

/*========================================================================
 This routine initializes the blimps global variables Qij, RTot 
	and frequency[].  If BLIMPS_DIR is set, assumes files are there,
	otherwise tries current directory.
==========================================================================*/
void frq_qij()
{
   FILE *fqij;
   char *blimps_dir, qijname[SMALL_BUFF_LENGTH], frqname[SMALL_BUFF_LENGTH];

   blimps_dir = getenv("BLIMPS_DIR");
   if (blimps_dir != NULL) 
   {
      sprintf(frqname, "%s/docs/default.amino.frq", blimps_dir);
      sprintf(qijname, "%s/docs/default.qij", blimps_dir);
/*>>>>> Need to check that frqname & qijname actually exist;
	for WWW servers would be better to check in current
	directory first  <<<<<<*/
   }
   else
   {
      sprintf(frqname, "default.amino.frq");
      sprintf(qijname, "default.qij");
   }
   load_frequencies(frqname);		/* creates frequency[]  */
   Qij = NULL;
   if ( (fqij=fopen(qijname, "r")) != NULL)
   { Qij = load_qij(fqij); fclose(fqij);  }
   RTot = LOCAL_QIJ_RTOT;
/*
   Matrix *matrix = block_to_matrix(Block *block, 3);
*/
}  /* end of frq_qij */

/* Change log information follows. */
/* 
 * Changes since 3.5:
   7/29/02  Modified error messages
 * Changes since 3.2:
   2/11/99  Added frq_qij()  JGH
 * Changes since 3.1:
   2/14/97  Added load_codons()  JGH
 */

