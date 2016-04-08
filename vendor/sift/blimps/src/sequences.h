/* (C) Copyright 1993, Fred Hutchinson Cancer Research Center */
/* Use, modification or distribution of these programs is subject to */
/* the terms of the non-commercial licensing agreement in license.h. */

/* sequences.h: sequence manipulation functions */
/* Written by: Bill Alford */
/* Change log information is at the end of the file. */

#ifndef SEQUENCES_H_
#define SEQUENCES_H_

#include "strutil.h"

/*
 * Exported variables and data structures
 */
#define AA_SEQ 0	/* Amino acid sequence designator */
#define NA_SEQ 1	/* Nucleic acid sequence designator */
#define UNKNOWN_SEQ 2	/* Unknown sequence designator */

#define NA_FREQ 0.85	/* the percent of 'A', 'C', 'G'. & 'T' in */
			/* the sequence to consider it a NA_SEQ */

/* the Residue type */
typedef unsigned char Residue;

struct sequence_struct {
  char name[SMALL_BUFF_LENGTH];	   /* the name of the sequence */
  char info[SMALL_BUFF_LENGTH];	   /* the info line of the sequence */
  int position;			   /* the position of the subsequence (if */
				   /* is one) in the sequence */
  int length;			   /* the length of the sequence */
  int max_length;		   /* the maximum length of memory allocated */
  int type;			   /* type of sequence: AA_SEQ or NA_SEQ */
  double weight;		   /* the weight of the sequence if the */
				   /* sequence is in a block */
  int undefined;		   /* for future use */
  double undefined_dbl;		   /* for future use */
  void *undefined_ptr;		   /* for future use */
  Residue *sequence;		   /* the sequence of either protein or DNA */
};
typedef struct sequence_struct Sequence;




/* Sequence Database types */
/*--------Sequence database format structure -------------------------*/
#define MAXDB   10		/* Max. # database formats */
#define GB       0		/* GenBank type */
#define PIR      1		/* PIR type */
#define EMBL     2		/* EMBL type.  SwissProt is the same. */
#define GCG      3		/* plain GCG type */
#define GCG_GB   4		/* GCG with GenBank header */
#define GCG_PIR  5		/* GCG with PIR header */
#define GCG_EMBL 6		/* GCG with EMBL/SwissProt header */
#define FASTA    7		/* FASTA type */
#define UNI      8  		/* UNIVERSAL type */
#define FLAT     9		/* FLAT, no title, just sequence */
/*#define VMS 3*/ 		/* PIR/VMS type */
/* NOTE: VMS must be the the last type and UNI must be just before it.  See */
/*       the code in type_dbs for the reason */
/* NOTE: Fasta and Universa (and VMS) formats are very similar.  The */
/*       only difference between Fasta and Universa is that there is a */
/*       '*' at the end of the Universa format.  The reason for the */
/*       Universa format is that there was a time when patmat needed */
/*       the '*' end-of-sequence delimiter.  There is no longer a */
/*       need with blimps, but the universa format handling is being */
/*       kept around for backward compatibility. */
/* NOTE: The GCG formats consists of a header file of the format of the */
/*       originating database and the sequence protion in the GCG */
/*       format.  This obviously can lead to confusion, especially */
/*       since the characters in the second line of the GCG sequence */
/*       format will signal that the sequence is an amino acid */
/*       sequence (bad when the GCG sequence is from a DNA database). */
/* NOTE: The plain GCG format is a kludge.  There is code to handle */
/*       it, but since it is one of the formats that does not have an */
/*       end delimiter there is the possibility that there will be */
/*       junk added to the sequence if there is more than one */
/*       sequence in the file. */

struct db_info {		/* Sequence Database format info */
   char *type;
   char *start;
   char *desc;
   char *seq;
   char *end;
   int title_offset;
   int seq_offset;
};

extern struct db_info DbInfo[MAXDB];



/*
 * Exported functions
 */

/*
 * read_a_sequence
 *   reads a sequence from the data base and returns a pointer to the new
 *   Sequence data structure
 *   Parameters:
 *     FILE *sfp: the sequence file pointer
 *     int db:    the database type of the file
 *     int type:  the type of sequences in the file
 *   Error codes: Returns NULL if there was no sequence to read.
 */

extern Sequence *read_a_sequence();

/*
 * read_sequence
 *   read_sequence reads from the string the sequence of the given type into
 *   the passed in Sequence data structure
 *   Parameters:
 *     Sequence *seq: the Sequence data structure to put the data in
 *     int seq_type:  the type of sequence (amino acid or nucleic acid)
 *     int start_pos: the position in the sequence to start entering the data
 *     char *str:     the string that has the sequence (NULL terminated)
 *   Note: The entire string up to the NULL character is considered the
 *         sequence.
 *   Preconditions: seq.length must have a the value of the size allocated to
 *                  the residues
 *   Postconditions: the residues from start_pos to either seq.length or
 *                   return value + start_pos in the seq.sequence have been
 *                   updated to the residues passed in by the string
 *                   seq.type is set to seq_type
 *   Return codes: returns the number of residues added to the sequence
 *   Error codes: returns the negative value of the number of residues left
 *                over after the sequence has been filled.
 *   Note: A return value of zero is means that no residues were
 *         added to the sequence.
 */

extern int read_sequence();

/*
 * sequence_type
 *   tries to guess the type of the sequence passed in the string.
 *   Note: if the sequence can not be matched to an NA_SEQ, it falls through
 *         to being a AA_SEQ.
 *   Parameters:
 *     char *str: the sequence string, no more than SMALL_BUFF_LENGTH long
 *   Return Codes: Returns either AA_SEQ or NA_SEQ depending on
 *                 if the sequence is an Amino Acid or Nucleic Acid.
 *   Error Codes: none
 */

extern int sequence_type();


/*
 * sequence_comparison
 *   Compares two sequences.   It compares by the value in the name
 *   field if it exists.
 *   Parameters:
 *     SequenceListEntry a, b: the entries to compare
 *   Return codes:  a return value < 0 if a < b, a return value = 0 if a == b,
 *                  and a return value > 0 if a > b
 *   Error codes: none
 */

extern int sequence_comparison();

/*
 * free_sequence
 *   Deletes the sequence and the sub elements.
 *   Parameters:
 *     Sequence *sequence: the sequence to free
 *   Return code: none
 *   Error code: none
 */

extern void free_sequence();



/*
 * translate_sequence
 *   Given a frame (1, 2, 3, -1, -2, or -3), translates the DNA sequence
 *   (NA_SEQ) into an amino acid sequence.
 *   Parameters:
 *     Sequence *seq: the DNA sequence
 *     int frame:     the frame to read
 *     unsigned char gcode[64]:    the genetic code (forward).
 *                                   for frames 1, 2, 3
 *     unsigned char revgcode[64]: the genetic code (reverse complement)
 *                                   for frames -1, -2, -3
 *   Return code: a pointer to the new translated sequence
 *   Error codes: returns the original sequence if the sequence is not a
 *                NA_SEQ or the translation frame is incorrect.
 */

extern Sequence *translate_sequence();

/*  untranslate an amino acid sequence to a degenerate DNA sequence  */
extern Sequence *untranslate_sequence();


/* NOTE: These should only be used in config.c when determining the */
/*       database types.  These should be used to set the values in the */
/*       file lists. */

/*
 * type_dbs
 *   Figures out the type of the database to be read.
 *   This is from motomisc.c used in the universa conversion program.
 *   Parameters:
 *     FILE *fin:            the database file
 *     struct db_info dbs[]: the information of the structure of the db's
 *   Return codes: the type of database the file is.
 *   Error code: the database number -1 is returned if the type of database
 *               cannot be figured out
 */
extern int type_dbs();


/*
 * seq_type_dbs
 *   Finds the type or the database given.
 *   Parameters:
 *     FILE *fin:            the database file
 *     struct db_info dbs[]: the information of the structure of the db's
 *     int db:               the number of this database (indexes into dbs)
 *   Return codes: Returns the type of sequence of the first entry in the
 *                 database.  This assumes that the entire database is of this
 *                 type.
 *   Error codes: none
 */

extern int seq_type_dbs();

/*
 * print_sequence
 *   Prints a Sequence data structure.  Primarily for debugging purposes.
 *   Parameters:
 *     Sequence *seq:  the sequence to print
 *   Error Codes: none
 */

extern void print_sequence();

/*
 * output_sequence
 *   Outputs a Sequence data structure to the given file.
 *   Parameters:
 *     Sequence *seq:  the sequence to print
 *     FILE *osfp:     the output sequence file pointer
 *   Error Codes: none
 */

extern void output_sequence();


/*
 * resize_sequence
 *   Adds more memory for the storage of the residues of a sequence.
 *   Parameter:
 *     Sequence *seq: the sequence to resize
 *   Error codes: none
 */

extern void resize_sequence();

#endif /*  SEQUENCES_H_ */

/* Change log information follows. */
/*
   Changes since version 3.3.2:
   12/15/99 Added double undefined_dbl field.
   Changes since version 3.2:
   11/24/97 Made resize_sequence() external.
   Changes since version 3.1:
   1/20/97  Added untranslate_sequence()
   Changes since version 2.x:
   3/20/96  Lowered NA_FREQ from 0.90 to 0.85 after user complaints. JGH
 *
 */
