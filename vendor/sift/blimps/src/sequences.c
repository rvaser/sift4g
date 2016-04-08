/* (C) Copyright 1993, Fred Hutchinson Cancer Research Center */
/* Use, modification or distribution of these programs is subject to */
/* the terms of the non-commercial licensing agreement in license.h. */

/* sequences.c: sequence manipulation functions */
/* Written by: Bill Alford */
/* Change log information is at the end of the file. */

/*	system headers not in global.h  */
/*	blimps library headers */
#include <global.h>
#include <sequences.h>
#include <residues.h>

/*
 * Exported variables and data structures
 */


struct db_info DbInfo[MAXDB] =
{
  {
    "GenBank",			/* type */
    "LOCUS",			/* start */
    "DEFINITION",		/* desc */
    "ORIGIN",			/* seq */
    "//",			/* end */
    12,				/* title_offset */
    10 },			/* seq_offset */
  {   /* pir has numbers in the first line following the sequence signif. */
    "PIR",			/* type */
    "ENTRY",			/* start */
    "TITLE",			/* desc */
    "SEQUENCE",			/* seq */
    "///",			/* end */
    16,				/* title_offset */
    8 },			/* seq_offset */
  {   /* EMBL and SwissProt have the same id's */
    "EMBL/SwissProt",		/* type */
    "ID   ",			/* start */
    "DE   ",			/* desc */
    "SQ",			/* seq */
    "//",			/* end */
    5,				/* title_offset */
    5 },			/* seq_offset */
  {   /* GCG has non-sequence info in the first three lines following the */
      /* sequence signifier */
      /* NOTE: these delimiters are used to help control the flow of */
      /*       the sequence algorithms.  If they are changed be very */
      /*       careful to make sure the changes work. */
    "GCG",			/* type */
    "",				/* start */
    "",				/* desc */
    "",				/* seq */
    "",				/* end */
    0,				/* title_offset */
    10 },			/* seq_offset */
  {   /* GCG has non-sequence info in the first three lines following the */
      /* sequence signifier */
    "GCG-GenBank",		/* type */
    "LOCUS",			/* start */
    "DEFINITION",		/* desc */
    "ORIGIN",			/* seq */
    "LOCUS",			/* end */
    12,				/* title_offset */
    10 },			/* seq_offset */
  {   /* GCG has non-sequence info in the first three lines following the */
      /* sequence signifier */
      /* This is not what the GCG PIR entries look like.  They have another */
      /* format that would be very difficult to parse with the current */
      /* algorithm.  As it is a GCG PIR sequence will be recognized as plain */
      /* GCG. */
    "GCG-PIR",			/* type */
    "ENTRY",			/* start */
    "TITLE",			/* desc */
    "SEQUENCE",			/* seq */
    "ENTRY",			/* end */
    16,				/* title_offset */
    10 },			/* seq_offset */
  {   /* GCG has non-sequence info in the first three lines following the */
      /* sequence signifier */
      /* EMBL and SwissProt have the same id's */
    "GCG-EMBL",			/* type */
    "ID",			/* start */
    "DE",			/* desc */
    "SQ",			/* seq */
    "ID",			/* end */
    5,				/* title_offset */
    10 },			/* seq_offset */
  {
    "Fasta",			/* type */
    ">",			/* start */
    ">",			/* desc */
    "",				/* seq */
    ">",			/* end */
    1,				/* title_offset */
    0 },			/* seq_offset */
  {
    "Universa",			/* type */
    ">",			/* start */
    ">",			/* desc */
    "",				/* seq */
    "*",			/* end */
    1,				/* title_offset */
    0 },			/* seq_offset */
  {
    "Flat-file, no title",	/* type */
    "",				/* start */
    "",				/* desc */
    "",				/* seq */
    "",				/* end */
    0,				/* title_offset */
    0 }/*,*/			/* seq_offset */
/*  {   Uncomment the lines in function type_dbs */
/*    "VMS",			*//* type */
/*    ">",			*//* start */	   /* first line */
/*    "",		        *//* desc */	   /* second line */
/*    "",			*//* seq */	   /* third line */
/*    "*",			*//* end */
/*    4,			*//* title_offset */ /* first line only */
/*    0 },			*//* seq_offset */
};




/*
 * Local variables and data structures
 */

#define SEQUENCE_RESIDUE_INITIAL_SIZE 1000
#define SEQUENCE_RESIDUE_INCREASE_SIZE 1000  /* default number of residues to */
				   /* allocate if the number of residues */
				   /* in the sequence is not specified */


/*
 * Function definitions
 */



#ifndef NO_SEQUENCE_READING

/*
 * Sequence file related functions
 *
 *   Sequence *read_a_sequence()
 *   int sequence_type(str)
 *   int read_sequence(seq, seq_type, start_pos, str)
 *   void free_sequence(seq)
 *   void process_sequence_db_line(seq, buff)
 *   void resize_sequence(seq)
 *
 */

static int read_sequence_header();
static void process_sequence_db_line();
void resize_sequence();


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

Sequence *read_a_sequence(sfp, db, type)
     FILE *sfp;			/* the sequence file pointer */
     int db;			/* the database type of the file*/
     int type;			/* the type of sequences in the file */
{
  Sequence *new_sequence;

  char *ptr;

  /*>>>BIG PROBLEM here if sfp is rewound between calls to read_a_sequence()
	because lbuff[] isn't cleared!!! <<<*/
  static char lbuff[LARGE_BUFF_LENGTH];	/* local Buffer to save the info */
					/* from the last line of the */
					/* previous sequence.  FASTA and GCG */
					/* formats have no special separator */
					/* between sequences other than the */
					/* beginning key */

/*  sfp = get_file(SEQUENCE_FILES);*/

  if (sfp == NULL) {
    /* no more data to read into sequences */
    new_sequence = NULL;
    return new_sequence;
  }

  /* get the database type.  Note: shouldn't directly access SequenceFiles */
/*  db = get_sequence_db_db_type(); */

  /* allocate space for the new Sequence */
  CheckMem(
	   new_sequence = (Sequence *) malloc(sizeof(Sequence))
	   );


  /*
   * Read the header
   */

/*>>> reset static lbuff[] first time sfp is read, otherwise has contents
      from the last time read_a_sequence() was called and sfp may have been
      rewound in the interim, or may even be a different file, but if
      it's initialized here only reads every other fasta record.
      The only time lbuff[] is initialized is when read_sequence_header()
	hits eof, so calling programs must read to end of file before
	rewinding!
  lbuff[0] = '\0';
<<<*/
  if (read_sequence_header(sfp, db, new_sequence, lbuff) == FALSE) {
    free(new_sequence);
    return NULL;
  }


  /*
   * Read the sequence
   */

  /* set the position to zero since this _is_ the sequence */
  /* and not a subsequence like in a block */
  new_sequence->position = 0;

  /* set the max_length to an initial value and allocate space */
  new_sequence->max_length = SEQUENCE_RESIDUE_INITIAL_SIZE;
  CheckMem(
	  new_sequence->sequence = (Residue *) calloc(new_sequence->max_length,
					      sizeof(Residue))
	   );

  /* set the length to zero, will increase as we read in the residues */
  new_sequence->length = 0;
  new_sequence->undefined = 0;
  new_sequence->undefined_dbl = 0.0;
  new_sequence->undefined_ptr = NULL;

  /* set the type of the sequence */
  new_sequence->type = type;  /* get_sequence_db_seq_type(); */

  do {
    /* If there is is a blank line or the line is shorter than the */
    /* sequence offset skip it so that the data from the previous line */
    /* is not read again.  This can happen by indexing over the end of */
    /* line ('\0') character. */
    if (!blank_line(lbuff) && ((int)strlen(lbuff) > DbInfo[db].seq_offset)) {
      process_sequence_db_line(new_sequence, &(lbuff[DbInfo[db].seq_offset]));
    }
  } while (!feof(sfp) &&
	   fgets(lbuff, LARGE_BUFF_LENGTH, sfp) != NULL &&
	   ((strncmp(lbuff, DbInfo[db].end, strlen(DbInfo[db].end)) != 0) ||
	    (db == FLAT) ||
	    ((db == GCG) &&
	     !(((ptr = strstr(lbuff, "Length:")) != NULL) &&
	       /*((ptr = strstr(ptr, "Type:")) != NULL) &&*/
	       ((ptr = strstr(ptr, "Check:")) != NULL) &&
	       ((ptr = strstr(ptr, "..")) != NULL)))));

  /* return the new sequence */
  return new_sequence;

} /* end read_a_sequence */



/*
 * read_sequence_header
 *   Reads the header information for the sequence.
 *   Parameters:
 *     FILE *sfp:    the sequence file pointer.
 *     int db:       the database type
 *     Sequence *new_sequence: the sequence to put the data in.
 *     char lbuff[]: the saved Buffer for the last line read.
 *   Error Codes: FALSE if the sequence could not be read, TRUE otherwise
>>>>If fasta title line is > LARGE_BUFF_LENGTH treats excess as sequence!<<<
 */

static int read_sequence_header(sfp, db, new_sequence, lbuff)
     FILE *sfp;			/* the sequence file pointer */
     int db;			/* the database type */
     Sequence *new_sequence;
     char lbuff[LARGE_BUFF_LENGTH];
{
  char *ptr;
  char title[LARGE_BUFF_LENGTH], temp[LARGE_BUFF_LENGTH];

  static Boolean first_time = TRUE;

  strcpy(new_sequence->name, "Unknown");
  strcpy(new_sequence->info, "Unknown");

  if (db != FLAT) {	/* if this is the flat file, there is no title */

    temp[0] = '\0';
    title[0] = '\0';

    /* move to the start of the sequence */
    /* NOTE: the order of the while loop is important.  In FASTA and GCG */
    /*       format the last line read when reading the previous sequence */
    /*       will have the title information of this sequence.  In order */
    /*       not to overwrite the information, it is necessary to fail */
    /*       before getting the next line. */
    if (first_time) {
      fgets(lbuff, LARGE_BUFF_LENGTH, sfp);
      lbuff = eat_whitespace(lbuff);
      first_time = FALSE;
    }
    while (strncmp(lbuff, DbInfo[db].start, strlen(DbInfo[db].start)) != 0 &&
	   fgets(lbuff, LARGE_BUFF_LENGTH, sfp) != NULL );

    /* if we reached the end of the file, then throw away the sequence */
    if (feof(sfp)) {
      return FALSE;		/* it did not read a sequence */
    }


    /* ------Working on the DbInfo[db].start line now, get ->name and also
             ->info if DbInfo[db].desc == .start             --------------*/
    /* if this is a simple GCG format, move to the info part.  the way */
    /* the simple GCG dbs entry is set up, it will work with this algorithm. */
    if (db == GCG)
    {
      while(!(((ptr = strstr(lbuff, "Length:")) != NULL) &&
	      /*((ptr = strstr(ptr, "Type:")) != NULL) &&*/
	      ((ptr = strstr(ptr, "Check:")) != NULL) &&
	      ((ptr = strstr(ptr, "..")) != NULL)) &&
	    (fgets(lbuff, LARGE_BUFF_LENGTH, sfp) != NULL));
      /* if we reached the end of the file, then throw away the sequence */
      if (feof(sfp)) {
	return FALSE;		/* it did not read a sequence */
      }
    } /* end of GCG format */

    /* get the name of the sequence. */
    strcpy(temp, lbuff);
    ptr = get_token(&temp[DbInfo[db].title_offset]);
    if (ptr == NULL) { strcpy(new_sequence->name, "UNKNOWN"); }
    else
    {
       if (strlen(ptr) > SMALL_BUFF_LENGTH)
       {  ptr[SMALL_BUFF_LENGTH] = '\0'; }
       strcpy(new_sequence->name, ptr);
    }
    temp[0] = '\0';

    /* check to see if the start and desc lines are the same, eg for
       fasta and universa formats. If so, use the rest of this line
       for ->info */
    if (strcmp(DbInfo[db].start, DbInfo[db].desc) == 0)
    {
     strcpy(temp, &lbuff[DbInfo[db].title_offset + strlen(new_sequence->name)]);
      remove_trailing_whitespace(temp);         /* get rid of CRLF */
      sprintf(title, "%s ", temp);
    }

    /*-----------------Done with DbInfo[db].start now ------------*/
    /* read up to the sequence and if see desc. information, save */
    /* NOTE: most of this is from universa.c */
    while(fgets(lbuff, LARGE_BUFF_LENGTH, sfp) != NULL &&
	  strncmp(lbuff, DbInfo[db].seq, strlen(DbInfo[db].seq)) != 0)
    {
      if (strncmp(lbuff, DbInfo[db].desc, strlen(DbInfo[db].desc)) == 0 &&
	  ((int)(strlen(title) + strlen(lbuff) +1) < LARGE_BUFF_LENGTH) )
      {
	if (title[0] != '\0') {
	  strncpy(temp, title, LARGE_BUFF_LENGTH);
	  strncat(temp, " ", LARGE_BUFF_LENGTH-strlen(temp));
	}
	strncat(temp, &lbuff[DbInfo[db].title_offset],
		LARGE_BUFF_LENGTH-strlen(temp));
	remove_trailing_whitespace(temp);         /* get rid of CRLF */
	strncpy(title, temp, LARGE_BUFF_LENGTH);
      }
    }   /*  end of DbInfo[db].desc */

    /* save the info that has been collected */
    strncpy(new_sequence->info, title, SMALL_BUFF_LENGTH);

    /*-----------------Done with DbInfo[db].desc now ------------*/
    /* if this type of database does not start right off with the sequence */
    /* data (like most databases), read the next line to be on top of the */
    /* sequence data. */
/*>>>> for fasta & universa, if lbuff doesn't have end-of-line,
then read more until first eol, seq starts after first eol. Sometimes
the title line is longer than LARGE_BUFF_LENGTH <<*/
    if (strcmp(DbInfo[db].seq, "") != 0 && !feof(sfp))
    {
      /* read in the next line, get to the sequence */
      fgets(lbuff, LARGE_BUFF_LENGTH, sfp);
    }

    /* if PIR, skip the top row of numbers */
    if (db == PIR && !feof(sfp))
    {
      fgets(lbuff, LARGE_BUFF_LENGTH, sfp);
    }
    /* if one of the GCG types (but not the simple GCG), */
    /* skip the top three rows */
    if (((db == GCG_GB) || (db == GCG_PIR) || ( db == GCG_EMBL)) &&
	!feof(sfp))
    {
      fgets(lbuff, LARGE_BUFF_LENGTH, sfp);
      fgets(lbuff, LARGE_BUFF_LENGTH, sfp);
      fgets(lbuff, LARGE_BUFF_LENGTH, sfp);
    }

  }  /* end of if the format is not FLAT (got the title and description */

  /* if we reached the end of the file, then throw away the sequence */
  if (feof(sfp)) {
    return FALSE;		/* it did not read a sequence */
  }

  return TRUE;			/* read in the header of the sequence  */
}  /*  end of read_sequence_header  */





/*
 * process_sequence_db_line
 *   Passed one line of the sequence portion of the database.  Reads the
 *   characters and converts them to the apropriate internal representation.
 *   Removes the excess spaces in the middle and ends.
 *   NOTE: this function adds the residues onto the end of the sequence
 *         passed to it.
 *   Parameters:
 *     Sequence *seq: the sequence to put the residues into.
 *     char *buff:    the string with the characters to residues to
 *                      add to the sequence
 *   Error codes: none
 */

static void process_sequence_db_line(seq, buff)
     Sequence *seq;
     char *buff;
{
  int start_pos;
  int estimated_length;
  int saved_length;
  int num_entered;

  estimated_length = strlen(buff);
  saved_length = seq->length;

  /* if there are more residues than there is space for, resize. */
  /* Note: The buff may have characters that are not residues, this is */
  /*       ok, there will just be a little extra space */
  while ((saved_length + estimated_length) > seq->max_length) {
    resize_sequence(seq);
  }

  /* temporarily increase the room to put the sequence into */
  start_pos = seq->length;
  seq->length = seq->max_length;

  num_entered = read_sequence(seq, seq->type, start_pos, buff);

/*  if (num_entered < 0) {  shouldn't happen */
    /* Error, wasn't enough room for all the sequences */
    /* resize to get enough room for the extra residues } */

  seq->length = saved_length + num_entered;

}
#endif /* NO_SEQUENCE_READING */


/*
 * resize_sequence
 *   Adds more memory for the storage of the residues of a sequence.
 *   Parameter:
 *     Sequence *seq: the sequence to resize
 *   Error codes: none
 */

void resize_sequence(seq)
     Sequence *seq;
{
  Residue *tmp_ptr;		/* don't want to ruin pointer in case there */
				/* is a repeated realloc call */

  /* announce info */
  sprintf(ErrorBuffer,
	  "Ran out of room for sequence %s.  Allocating more space.\n",
	  seq->name);
  ErrorReport(INFO_ERR_LVL);

  seq->max_length += SEQUENCE_RESIDUE_INCREASE_SIZE;
  CheckMem(
	   tmp_ptr = (Residue *) realloc(seq->sequence,
					 seq->max_length *
					 sizeof(Residue))
	   );
  seq->sequence = tmp_ptr;
}



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

int sequence_type(str)
     char *str;			/* no more than SMALL_BUFF_LENGTH long */
{
  int length;
  int num_na;
  char *str_pntr;
  char upper_str[SMALL_BUFF_LENGTH], *c_upper_str;
  char residue;

  /* Sets of residues:
   *
   *                               AA-amb = -ABCDEFGHI KLMN PQRST VWXYZ*
   *                              !AA-amb =           J    O     U
   *
   *                               NA-amb = -ABCD  GH  K MN   RST VW Y
   *                              !NA-amb =      EF  IJ L  OPQ   U  X Z
   *
   *    AA-amb intersect NA-amb           = -ABCD  GH  K MN   RST VW Y
   * !(AA-amb intersect NA-amb) & !AA-amb
   *                  = !NA-amb - !AA-amb =      EF  I  L   PQ      X Z
   *                             the rest =           J    O     U
   *
   *                            AA-wo-amb = -A CDEFGHI KLMN PQRST VW Y *
   *                           !AA-wo-amb =   B       J    O     U  X Z
   *
   *                            NA-wo-amb = -A C   G      N     T
   *                           !NA-wo-amb =   B DEF HIJKLM OPQRS UVWXYZ
   *
   */

  /* NOTE, U could be in a nucleotide sequence, if someone uses an RNA seq. */


  for (c_upper_str=upper_str; *str != '\0'; str++) {
    if (aa_atob[(int)*str] <= AAID_MAX) { /* get rid of bad characters.  It is ok */
                                     /* to use aa_atob since all of the NT */
				     /* residue characters are included in */
				     /* that set. */
      *c_upper_str = toupper(*str);
      c_upper_str++;
    }
  }

  *c_upper_str = '\0';

  length = strlen(upper_str);

  /* checking first assuming using ambiguity codes */
  /* if see any in (!NA-amb - !AA-amb) this sequence is an AA_SEQ */
  for (str_pntr=upper_str; *str_pntr; str_pntr++) {
    residue = *str_pntr;
    if ((residue == 'E') ||
	(residue == 'F') ||
	(residue == 'I') ||
	(residue == 'L') ||
	(residue == 'P') ||	/* not using 'X', 'Z', or '*' since someone */
	(residue == 'Q')) {	/* might use in a NA_SEQ */
      return AA_SEQ;
    }
  }

  /* check frequency of 'A', 'C', 'G', 'T', 'U' & 'N'  if > 85% (NA_FREQ), */
  /* assume NA_SEQ */
  num_na = 0;
  for (str_pntr=upper_str; *str_pntr; str_pntr++) {
    residue = *str_pntr;
    if ((residue == 'A') ||
	(residue == 'C') ||
	(residue == 'G') ||
	(residue == 'T') ||
	(residue == 'U') ||
	(residue == 'N')) {
      num_na++;
    }
  }
  if ( (double)((double)num_na/(double)length) > NA_FREQ ) {
    return NA_SEQ;
  }
  else {
    return AA_SEQ;
  }

}


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

int read_sequence(seq, seq_type, start_pos, str)
     Sequence *seq;
     int seq_type;
     int start_pos;
     char *str;
{
  Residue *sequence;


/*   return 0; */


  int room_left;
  int num_seen;
  int num_left;

  int return_code;

  int stemp;
  char c;

  sequence = &(seq->sequence[start_pos]);
  room_left = seq->length - start_pos;

  num_seen = 0;

  if (str == NULL) { return num_seen; }

  c = *str;

  /*  NOTE:  the stemp business was added to side-step a problem
      with aa_atob[c] that I don't understand, crops up if there
      are numbers in the sequence; for instance, if c = '4' then
      stemp is a huge negative number ...
  */
  if (seq_type == AA_SEQ) {
   /* while ( (c != NULL) && (c != '\0') && (num_seen < room_left) ) { */
    while ( (c != 0) && (c != '\0') && (num_seen < room_left) ) {
      stemp = aa_atob[(int)c];
      if (stemp < 0 || stemp > AAID_MAX) stemp = AAID_MAX + 1;
      if ( stemp <= AAID_MAX ) {
	*sequence = stemp;
	sequence++;
	num_seen++;
      }
      str++;
      c = *str;
    }
  }
  else if (seq_type == NA_SEQ) {
    while ( (c != '\0') && (num_seen < room_left) ) {
      if ( nt_atob[(int)c] <= NUCID_MAX ) {
	*sequence = nt_atob[(int)c];
	sequence++;
	num_seen++;
      }
      str++;
      c = *str;
    }
  }
  else {
    sprintf(ErrorBuffer,
	    "read_sequence(): Unknown sequence type,");
    ErrorReport(PROGRAM_ERR_LVL);
    sprintf(ErrorBuffer,
	    "                 asumming amino acid sequence\n");
    ErrorReport(PROGRAM_ERR_LVL);
    while ( (c != '\0') && (num_seen < room_left) ) {
      if ( aa_atob[(int)c] <= AAID_MAX ) {
	*sequence = aa_atob[(int)c];
	sequence++;
	num_seen++;
      }
      str++;
      c = *str;
    }
  }


  if (c != '\0') {		/* num_seen >= room_left */
    /* count the number that were left */
    num_left = 0;
    if (seq_type == NA_SEQ) {
      while (c != '\0') {
	if (nt_atob[(int)c] <= NUCID_MAX ) {
	  num_left++;
	}
	str++;
	c = *str;
      }
    }
    else {			/* if not an NA_SEQ assuming an AA_SEQ */
      while (c != '\0') {
	if (aa_atob[(int)c] <= AAID_MAX ) {
	  num_left++;
	}
	str++;
	c = *str;
      }
    }

    if (num_left == 0) {	/* read in all residues so just return */
      return_code = num_seen;	/* the number seen */
    }
    else {
      return_code = - num_left;
    }

  }
  else {		      /* c == '\0', reached the end of the residues */
    return_code = num_seen;
  }

  return return_code;

}


/*
 * free_sequence
 *   Deletes the sequence and the sub elements.
 *   Parameters:
 *     Sequence *sequence: the sequence to free
 *   Return code: none
 *   Error code: none
 */

void free_sequence(sequence)
     Sequence *sequence;
{
  free(sequence->sequence);
  free(sequence);
}


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

Sequence *translate_sequence(seq, frame, gcode, revgcode)
     Sequence *seq;
     int frame;
     unsigned char gcode[64];
     unsigned char revgcode[64];        /* genetic codes */
{
  Sequence *new_seq;
  Residue *res;
  Residue *new_res;

  int i, new_length;

  new_length = -1;

  /* make sure it is a NA_SEQ */
  if ( seq->type != NA_SEQ ) {
    sprintf(ErrorBuffer,
	    "translate_sequence(): Not a nucleic acid sequence, not translating.");
    ErrorReport(PROGRAM_ERR_LVL);
    sprintf(ErrorBuffer,
	    "                      Returning the original sequence.\n");
    ErrorReport(PROGRAM_ERR_LVL);
    return seq;
  }

  /* make sure it is the right kind of frame */
  if ( (frame == 0) || (frame < -3) || (frame > 3) ) {
    sprintf(ErrorBuffer,
	    "translate_sequence(): Unknown translation frame: %d.",
	    frame);
    ErrorReport(PROGRAM_ERR_LVL);
    sprintf(ErrorBuffer,
	    "                      Returning the original sequence.\n");
    ErrorReport(PROGRAM_ERR_LVL);
    return seq;
  }

  /* allocate space for new sequence */
  CheckMem(
	   new_seq = (Sequence *) malloc(sizeof(Sequence))
	   );

  /* allocate space for the new residues */
  CheckMem(
	   new_res = (Residue *) calloc((seq->length+2)/3,
					sizeof(Residue))
	   );
  new_seq->max_length = (seq->length+2)/3;

  /* translate the codons */
  res = seq->sequence;
  if (frame > 0) {		/* 1, 2, or 3 */
    frame--;
    for (i=frame, new_length=0; i<= seq->length - 3; i+=3, new_length++) {
      new_res[new_length] = (Residue) codon2aa(gcode,
				     res[i],
				     res[i+1],
				     res[i+2]);
    }
    /* do the last 1, or 2 residues */
    if ((seq->length - i) == 1) {
      new_res[new_length] = (Residue) codon2aa(gcode,
				     res[i],
				     nt_atob[UNKNOWN_NT_CHR],
				     nt_atob[UNKNOWN_NT_CHR]);
      new_length++;
    }
    else if ((seq->length - i) == 2) {
      new_res[new_length] = (Residue) codon2aa(gcode,
				     res[i],
				     res[i+1],
				     nt_atob[UNKNOWN_NT_CHR]);
      new_length++;
    }
  }

  else if (frame < 0) {
    frame = frame*-1;
    frame--;
    for (i=(seq->length-frame-1), new_length=0; i>=2; i-=3, new_length++) {
      new_res[new_length] = (Residue) codon2aa(revgcode,
				     res[i],
				     res[i-1],
				     res[i-2]);
    }
    /* do the last 1 or 2 residues */
    if (i == 0) {
      new_res[new_length] = (Residue) codon2aa(revgcode,
				     res[i],
				     nt_atob[UNKNOWN_NT_CHR],
				     nt_atob[UNKNOWN_NT_CHR]);
      new_length++;
    }
    else if (i == 1) {
      new_res[new_length] = (Residue) codon2aa(revgcode,
				     res[i],
				     res[i-1],
				     nt_atob[UNKNOWN_NT_CHR]);
      new_length++;
    }
  }

  /* copy the sequence info into the new sequence */
  strncpy(new_seq->name, seq->name, SMALL_BUFF_LENGTH);
  strncpy(new_seq->info, seq->info, SMALL_BUFF_LENGTH);

  new_seq->position = seq->position;
  new_seq->length   = new_length;
  /* new_seq->max_length set when residues were allocated */

  new_seq->type = AA_SEQ;

  new_seq->sequence = new_res;

  return new_seq;

}


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

int sequence_comparison(a, b)
     Sequence *a, *b;
{
  char *sa, *sb;

  if (a->name != NULL) {
    sa = a->name;
  }
  else {
    if (b->name != NULL) {
      return -1;		/* NULL < something */
    }
    else {
      return 0;			/* NULL = NULL */
    }
  }

  if (b->name != NULL) {
    sb = b->name;
  }
  else {
    return 1;			/* something > NULL */
  }

  return (int) (strcmp(sa, sb));
}





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
/* NOTE: DO NOT RETURN FLAT AS THE DATABASE, WANT TO NOTE THAT IT IS NOT */
/*       A GOOD FORMAT */

/* originally from motomisc.c.  used in universa.c */
/*======================================================================
      type_dbs() determines what type a database is from the allowable
      types
========================================================================*/
int type_dbs(fin, dbs)
FILE *fin;
struct db_info dbs[];
{
   int db, i;
   char line[LARGE_BUFF_LENGTH];
   char *ptr;

   db = -1;

   /* Figure out what type of input file it is */
   while (db < 0 && fgets(line, sizeof(line), fin) != NULL) {
     for (i=0; i<MAXDB; i++) {
       /* don't check against GCG or FLAT, they match everything.  Don't */
       /* check against GCG-* because we will check later. */
       if ((i!=GCG) && (i!=GCG_GB) && (i!=GCG_PIR) && (i!=GCG_EMBL) &&
	   (i!=FLAT)) {
	 if (strncmp(line, dbs[i].start, strlen(dbs[i].start)) == 0) {
	   db = i;
	 }
       }
     }
     /* if this line didn't match anything, check if it is plain GCG */
     if ((db < 0) &&
	 ((ptr = strstr(line, "Length:")) != NULL) &&
	 /*((ptr = strstr(ptr, "Type:")) != NULL) &&*/
	 ((ptr = strstr(ptr, "Check:")) != NULL) &&
	 ((ptr = strstr(ptr, "..")) != NULL)) {
       db = GCG;
     }
   }


   /* decide if the type is a GCG type */
   if ((db == GB) || (db == PIR) || ( db == EMBL)) {
     /* announce about situation */
     sprintf(ErrorBuffer,
	     "Deciding if this is a GCG sequence database format of %s.",
	     dbs[db].type);
     ErrorReport(INFO_ERR_LVL);

     /* skip to the beginning of the sequence */
     while (fgets(line, LARGE_BUFF_LENGTH, fin) &&
	    !((strncmp(line, dbs[GB].seq, strlen(dbs[GB].seq)) == 0) ||
	      (strncmp(line, dbs[PIR].seq, strlen(dbs[PIR].seq)) == 0) ||
	      (strncmp(line, dbs[EMBL].seq, strlen(dbs[EMBL].seq)) == 0)));

     /* get the line after the sequence that would have GCG info to */
     /* check against */
     fgets(line, LARGE_BUFF_LENGTH, fin);
     fgets(line, LARGE_BUFF_LENGTH, fin);

     /* check for the GCG identifiers.  They must be in order. */
     if (((ptr = strstr(line, "Length:")) != NULL) &&
	 /*((ptr = strstr(ptr, "Type:")) != NULL) &&*/
	 ((ptr = strstr(ptr, "Check:")) != NULL) &&
	 ((ptr = strstr(ptr, "..")) != NULL)) {
       /* Finish announcement */
       sprintf(ErrorBuffer,
	       "Decided it is a GCG database format of %s.",
	       dbs[db].type);
       ErrorReport(INFO_ERR_LVL);

       /* found the identifiers, it is a GCG format. */
       switch (db) {
       case GB:
	 db = GCG_GB;
	 break;
       case PIR:
	 db = GCG_PIR;
	 break;
       case EMBL:
	 db = GCG_EMBL;
	 break;
       default: /* error */
	 sprintf(ErrorBuffer,
		 "type_dbs(): Unknown database format for GCG type.");
	 ErrorReport(PROGRAM_ERR_LVL);
	 sprintf(ErrorBuffer,
		 "            Not all of the database formats were handled.\n");
	 ErrorReport(PROGRAM_ERR_LVL);
	 break;
       }
     }
     else {
       sprintf(ErrorBuffer,
	       "Decided it is not a GCG database format of %s.",
	       dbs[db].type);
       ErrorReport(INFO_ERR_LVL);
     }
   }


   /* NOTE: when VMS is added, the check for VMS type must be done first */
   /*       since the sequence is all on one line in VMS format and the next */
   /*       line is the start of the next sequence, which looks like FASTA */
   /*       format */
   /*  Don't separate UNI and FASTA */
   if ( db == UNI ) db = FASTA;

   if (db < 0 || db >= MAXDB) {
     db = -1;			/* can't tell what it is */
     sprintf(ErrorBuffer,
	     "Could not figure out database format.\n");
   }
   else {
     sprintf(ErrorBuffer,
	     "Decided on %s database format.\n", dbs[db].type);
   }
   ErrorReport(INFO_ERR_LVL);

   rewind (fin);

   return(db);
}  /* end of type_dbs */


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

int seq_type_dbs(fin, dbs, db, seqtype)
     FILE *fin;
     struct db_info dbs[];
     int db;
     int seqtype;
{
  char seq_check_buf[SMALL_BUFF_LENGTH];
  int check_length;
  char *ptr;

  /* skip the title to the sequence */
  if (db >= 0 && db != GCG)
  {
    while(fgets(Buffer, LARGE_BUFF_LENGTH, fin) != NULL &&
	  strncmp(Buffer, dbs[db].seq, strlen(dbs[db].seq)) != 0);
  }
  else if (db != FLAT) { /* don't do anything if it is just sequence */
    while(fgets(Buffer, LARGE_BUFF_LENGTH, fin) != NULL &&
	  !(((ptr = strstr(Buffer, "Length:")) != NULL) &&
	    /*((ptr = strstr(ptr, "Type:")) != NULL) &&*/
	    ((ptr = strstr(ptr, "Check:")) != NULL) &&
	    ((ptr = strstr(ptr, "..")) != NULL)));
    fgets(Buffer, LARGE_BUFF_LENGTH, fin);
  }

  /* determine the sequence type */
  /* if PIR, skip the top row of numbers */
  if (db == PIR) {
    fgets(Buffer, LARGE_BUFF_LENGTH, fin);
  }
  /* if one of the GCG types, skip the top three rows */
  if ((db == GCG_GB) || (db == GCG_PIR) || ( db == GCG_EMBL)) {
    fgets(Buffer, LARGE_BUFF_LENGTH, fin);
    fgets(Buffer, LARGE_BUFF_LENGTH, fin);
    fgets(Buffer, LARGE_BUFF_LENGTH, fin);
  }

  check_length = 0;	   /* how much has been filled in the seq_check_buf */

  /* fill the seq_check_buf with sequence characters */
  /* db types with null end require special handling */
  if (db >= 0)
  {
     while(!feof(fin) && fgets(Buffer, LARGE_BUFF_LENGTH, fin) != NULL &&
	((strncmp(Buffer, dbs[db].end, strlen(dbs[db].end)) != 0) ||
         (db == FLAT) ||
	 ((db == GCG) &&
	  !(((ptr = strstr(Buffer, "Length:")) != NULL) &&
	    /*((ptr = strstr(ptr, "Type:")) != NULL) &&*/
	    ((ptr = strstr(ptr, "Check:")) != NULL) &&
	    ((ptr = strstr(ptr, "..")) != NULL))))        &&
	(check_length < SMALL_BUFF_LENGTH-1))
     {	/* -1 for the '\0' */
       strncpy(seq_check_buf + check_length,
	    Buffer+dbs[db].seq_offset,
	    SMALL_BUFF_LENGTH - check_length - 1);
       seq_check_buf[SMALL_BUFF_LENGTH-1] = '\0';
       check_length = strlen(seq_check_buf);
     }
  }

  rewind(fin);

  /*----  If the TY config file option was specified, use it ------*/
  if (seqtype == AA_SEQ || seqtype == NA_SEQ)
     return seqtype;
  else
     return sequence_type(seq_check_buf); /* get the sequence type */

} /* end of seq_type_dbs */



/*
 * Sequence printing.
 */

/*
 * print_sequence
 *   Prints a Sequence data structure.  Primarily for debugging purposes.
 *   Parameters:
 *     Sequence *seq:  the sequence to print
 *   Error Codes: none
 */

void print_sequence(seq)
     Sequence *seq;
{
  int k;

  printf(">--- sequence ---\n");

  printf("%s\n", seq->info);
  printf(">%s \t", seq->name);
  printf("( %d:", seq->position);
  printf("%d, ",  seq->length);

  if (seq->type == AA_SEQ) {
    printf("AA) \t");
printf("\n");
    for(k=0; k<seq->length; k++) {
      printf("%c", aa_btoa[seq->sequence[k]]);
    }
  }
  else if (seq->type == NA_SEQ) {
    printf("NA) \t");
printf("\n");
    for(k=0; k<seq->length; k++) {
      printf("%c", nt_btoa[seq->sequence[k]]);
    }
  }
  else {
    printf("\?\?) \t");
printf("\n");
    for(k=0; k<seq->length; k++) {
      printf("%c", aa_btoa[seq->sequence[k]]);
    }
  }

  printf("\n");
  printf(">--- -------- ---\n");
}


/*
 * output_sequence
 *   Outputs a Sequence data structure to the given file.
 *   Parameters:
 *     Sequence *seq:  the sequence to print
 *     FILE *osfp:     the output sequence file pointer
 *   Error Codes: none
 */

void output_sequence(seq, osfp)
     Sequence *seq;
     FILE *osfp;
{
  int k;


  fprintf(osfp, ">%s  %s\n", seq->name, seq->info);

  if (seq->type == AA_SEQ) {
    for(k=0; k<seq->length; k++) {
      fprintf(osfp, "%c", aa_btoa[seq->sequence[k]]);
      if ((k+1)%60 == 0) {
	fprintf(osfp, "\n");
      }
    }
  }
  else if (seq->type == NA_SEQ) {
    for(k=0; k<seq->length; k++) {
      fprintf(osfp, "%c", nt_btoa[seq->sequence[k]]);
      if ((k+1)%60 == 0) {
	fprintf(osfp, "\n");
      }
    }
  }
  else {
    for(k=0; k<seq->length; k++) {
      if ((k+1)%60 == 0) {
	fprintf(osfp, "\n");
      }
      fprintf(osfp, "%c", aa_btoa[seq->sequence[k]]);
    }
  }

  fprintf(osfp, "\n");
}




/*
 * untranslate_sequence
 *   Translates an animo acid sequence into a  DNA sequence of degenerate
 *    codons
 *   Parameters:
 *     Sequence *seq: the amino acid sequence sequence
 *   Return code: a pointer to the new untranslated sequence
 *   Error codes: returns the original sequence if the sequence is not a
 *                AA_SEQ
 */

Sequence *untranslate_sequence(seq)
     Sequence *seq;
{
  Sequence *new_seq;
  Residue *res;
  Residue *new_res;
  int i;


  /* make sure it is a AA_SEQ */
  if ( seq->type != AA_SEQ ) {
    sprintf(ErrorBuffer,
    "untranslate_sequence(): Not an amino acid sequence, not untranslating.");
    ErrorReport(PROGRAM_ERR_LVL);
    sprintf(ErrorBuffer,
	    "                      Returning the original sequence.\n");
    ErrorReport(PROGRAM_ERR_LVL);
    return seq;
  }

  /* allocate space for new sequence */
  CheckMem(
	   new_seq = (Sequence *) malloc(sizeof(Sequence))
	   );

  /* allocate space for the new residues */
  CheckMem(
	   new_res = (Residue *) calloc((seq->length)*3,
					sizeof(Residue))
	   );
  new_seq->length = new_seq->max_length = (seq->length) * 3;
  new_seq->type = NA_SEQ;
  new_seq->position = seq->position;
  strcpy(new_seq->name, seq->name);
  strcpy(new_seq->info, seq->info);

  /* untranslate the amino acids */
  res = seq->sequence;
  for (i=0; i < seq->length; i++)
  {
     aa2codon(res[i], &new_res[i*3], &new_res[i*3+1], &new_res[i*3+2]);
  }

  /* copy the sequence info into the new sequence */
  strncpy(new_seq->name, seq->name, SMALL_BUFF_LENGTH);
  strncpy(new_seq->info, seq->info, SMALL_BUFF_LENGTH);

  new_seq->sequence = new_res;

  return new_seq;

}




/* Change log information follows. */
/*
  >Changes since version 3.6:
   4/15/04 Protect against aborts due to nulls strings
  >Changes since version 3.3.2:
  12/27/99 seq_types_db()  Check for db >= 0
  >Changes since version 3.2.5:
   5/22/99  Increased SEQUENCE_RESIDUE_INCREASE_SIZE from 500 to 1000.
	    Fixed unallocated memory problem in translate_sequence()
   4/22/99  Changed db->start for EMBL from "ID" to "ID   ", etc.
  >Changes since version 3.2.2:
  12/27/97  Increased SEQUENCE_RESIDUE_INCREASE_SIZE from 350 to 500.
            Added SEQUENCE_RESIDUE_INITIAL_SIZE = 1000.
  12/20/97  Don't duplicate ->name in ->info field.
  11/24/97  Made resize_sequence() external.
  10/25/97  read_sequence_header(): Use eat_whitespace() to get rid of
	    leading whitespace in sequence header (WWW problem).
  >Changes since version 3.2.1:
   7/12/97  Changed seq_type_dbs() to work correctly when dbs[db].end is null.
   7/ 7/97  Discarded old "universa" format where sequences have * at end
            because blimps required the * to be on a separate line. Now
	    this format is treated as fasta and any asterisk is scored as
	    as stop codon, which gets the minimum PSSM score for the 20
	    major AAs.
  >Changes since version 3.2:
   4/4/97   Initialized sequence->name & ->info in read_sequence_header().
  >Changes since version 3.1:
   1/20/97  Added untranslate_sequence().
  >Changes since version 3.0.0:
   4/24/96  Added seqtype to seq_type_dbs().
            Treat N as a valid DNA character in sequence_type().
   5/16/96  Problem with aa_atob[] in read_sequence() if sequence
	    has numbers in it ...
*/
