/* (C) Copyright 1993, Fred Hutchinson Cancer Research Center */
/* Use, modification or distribution of these programs is subject to */
/* the terms of the non-commercial licensing agreement in license.h. */

/* config.c: reads the configuration file. */
/* Written by: Bill Alford */
/* Change log information is at the end of the file. */


/*	system headers not in global.h  */
/*	blimps library headers */
#include <global.h>
#include <files.h>
#include <sequences.h>
#include <pattern.h>
#include <options.h> 
/*	headers in current directory */
#include "config.h"
#include "blimps.h"
#include "scoring.h"		/* for DoHistogram */
#include "lists.h"


/* To prevent ansi compiler warnings */

int strcasecmp(const char *s1, const char *s2);


/*
 * Exported variables and data structures
 */
 
/*
 * Local variables and data structures
 */
 
static Boolean BL_field_seen = FALSE;
static Boolean MA_field_seen = FALSE;

static void read_config_file_end();


/*
 * Function definitions
 */
 
/* 
 * Configuration file related functions
 *
 *   int read_config_file (filename) 
 *   void read_config_file_end()
 *   
 */

/* 
 * read config file functions
 *   These are functions for each posible config file entry.  The first
 *   two characters of the function name are usually the config file
 *   field.
 *   Parameters: See the individual functions code.  Most do not have
 *               parameters and rely on the global buffer and get_token.
 *   Return codes: Most are none.  See the individual functions code. 
 *   Error codes:  Most are none.  See the individual functions code. 
 */

static void unknown_conf_key(buf_token)
     char *buf_token;
{
  sprintf(ErrorBuffer, 
	  "Unknown configuration file key: %s\n", buf_token);
  ErrorReport(INFO_ERR_LVL);
}

static void BL_block_filename()
{
  char *buf_token;

  BL_field_seen = TRUE;
  if (MA_field_seen) {
    sprintf(ErrorBuffer,
	    "MA field already used.  Cannot have both the BL and the MA");
    ErrorReport(WARNING_ERR_LVL);
    sprintf(ErrorBuffer,
	    "fields.  Ignoring the BL line.\n");
    ErrorReport(WARNING_ERR_LVL);
  }

  buf_token = get_token(NULL);
  insert_file(buf_token, BLOCK_FILES);
}

static void PA_pattern_filename()
{
  char *buf_token;

  sprintf(ErrorBuffer,
	  "Using the PAttern switch causes the REpeats switch to be ignored.");
  ErrorReport(WARNING_ERR_LVL);
  sprintf(ErrorBuffer,
	  "Alignments reported will include repeats.\n");
  ErrorReport(WARNING_ERR_LVL);

  UsePatterns = TRUE; 
  buf_token = get_token(NULL);
/*printf("Saw pattern filename: %s\n", buf_token);*/
  insert_file(buf_token, PATTERN_FILES); 
/*printf("%d pattern files\n", number_of_files(PATTERN_FILES));*/
}

static void CO_conversion_method()
{
  char *buf_token;

  /* the scoring method to use */
  buf_token = get_token(NULL);
  sscanf(buf_token, "%d", &BlockToMatrixConversionMethod);
}

static void DB_db_filename()
{
  char *buf_token;

  buf_token = get_token(NULL);
  insert_file(buf_token, DATABASE_FILES);
}

static void ER_error_level_to_report()
{
  char *buf_token;

  /* the reporting level to do */
  buf_token = get_token(NULL);
  sscanf(buf_token, "%d", &ErrorLevelReport);
}

static void EX_export_matr_filename()
{
  char *buf_token;

  buf_token = get_token(NULL);
  strncpy(ExportMatrixFile, buf_token, SMALL_BUFF_LENGTH);
}

static void FR_freq_filename()
{
  char *buf_token;

  buf_token = get_token(NULL);
  insert_file(buf_token, FREQUENCY_FILE);
}

static void GE_genetic_code()
{
  char *buf_token;

  /* the genetic code to use */
  buf_token = get_token(NULL);
  sscanf(buf_token, "%d", &GeneticCodeInitializer);
}



/* strcasecmp is not a standard ANSI library function it needs
 * a function prototype */

static void HI_histogram()
{
  char *buf_token;
	
  buf_token = get_token(NULL);
  if ((strcasecmp(buf_token, "TRUE") == 0) ||
      (strcasecmp(buf_token, "YES") == 0) ||
      (strcasecmp(buf_token, "T") == 0) ||
      (strcasecmp(buf_token, "Y") == 0) ||
      ((isdigit(buf_token[0])) && (strcmp(buf_token, "0") != 0))) {
    DoHistogram = TRUE;
  }
  else {
    DoHistogram = FALSE;
  }
}

static void MA_Matrix_db_filename()
{
  char *buf_token;

  MA_field_seen = TRUE;
  if (BL_field_seen) {
    sprintf(ErrorBuffer,
	    "BL field already used.  Cannot have both the MA and the BL");
    ErrorReport(WARNING_ERR_LVL);
    sprintf(ErrorBuffer,
	    "fields.  Ignoring the MA line.\n");
    ErrorReport(WARNING_ERR_LVL);
  }

  buf_token = get_token(NULL);
  insert_file(buf_token, BLOCK_FILES);
}

static void NU_number_to_report()
{
  char *buf_token;

  /* the number of results to save/report */
  buf_token = get_token(NULL);
  sscanf(buf_token, "%d", &NumberToReport);
}

static void OP_options_for_algs()
{
  insert_into_options(NULL);
}

static void OU_output_filename()
{
  char *buf_token;

  buf_token = get_token(NULL);
  strncpy(OutputFile, buf_token, SMALL_BUFF_LENGTH);
}

static void RE_repeats_allowed()
{
  char *buf_token;
	
  buf_token = get_token(NULL);
  if ((strcasecmp(buf_token, "TRUE") == 0) ||
      (strcasecmp(buf_token, "YES") == 0) ||
      (strcasecmp(buf_token, "T") == 0) ||
      (strcasecmp(buf_token, "Y") == 0) ||
      ((isdigit(buf_token[0])) && (strcmp(buf_token, "0") != 0))) {
    RepeatsAllowed = TRUE;
  }
  else {
    RepeatsAllowed = FALSE;
  }
}

static void SC_scoring_method()
{
  char *buf_token;

  /* the scoring method to use */
  buf_token = get_token(NULL);
  sscanf(buf_token, "%d", &SequenceMatrixScoringMethod);
}

static void SE_search_type()
{
  char *buf_token;

  sprintf(ErrorBuffer,
	  "Old style configuration file.");
  ErrorReport(INFO_ERR_LVL);
  sprintf(ErrorBuffer,
	  "New style does not use the SE field.\n");
  ErrorReport(INFO_ERR_LVL);

  buf_token = get_token(NULL);
  if (!strcasecmp(buf_token, "MATRIX")) {
    SearchType = SEARCH_TYPE_MATRIX;
  }
  else if (!strcasecmp(buf_token, "BLOCK")) {
    SearchType = SEARCH_TYPE_BLOCK;
  }
  else {
    SearchType = SEARCH_TYPE_UNKNOWN;
    sprintf(ErrorBuffer, 
	    "Unknown search type: %s\n", buf_token);
    ErrorReport(WARNING_ERR_LVL);
  }
}

static void SQ_sequence_filename()
{
  char *buf_token;
  
  buf_token = get_token(NULL);
  insert_file(buf_token, SEQUENCE_FILES);
}

static void ST_strands_to_search()
{
  char *buf_token;

  /* number of strands to search: 1 (forward), -1 (reverse), 2 (both) */
  buf_token = get_token(NULL);
  sscanf(buf_token, "%d", &StrandsToSearch);
  if (StrandsToSearch < 0) StrandsToSearch = -1;
  if (StrandsToSearch == 0 || StrandsToSearch > 2) StrandsToSearch = 2;
}

static void SV_saved_scores()
{
  char *buf_token;
	
  buf_token = get_token(NULL);
  if ((strcasecmp(buf_token, "TRUE") == 0) ||
      (strcasecmp(buf_token, "YES") == 0) ||
      (strcasecmp(buf_token, "T") == 0) ||
      (strcasecmp(buf_token, "Y") == 0) ||
      ((isdigit(buf_token[0])) && (strcmp(buf_token, "0") != 0))) {
    SavedScoresFlag = TRUE;
  }
  else {
    SavedScoresFlag = FALSE;
  }
}

static void TY_sequence_type()
{
  char *buf_token;
  
  /* DNA => force DNA, AA => force AA  */
  buf_token = get_token(NULL);
  if (!strcasecmp(buf_token, "DNA"))   { SequenceType = NA_SEQ; }
  else 
  {
     if (!strcasecmp(buf_token, "AA")) { SequenceType = AA_SEQ; }
     else                              { SequenceType = UNKNOWN_SEQ; }
  }
}



/*
 * init_config_vars
 *   Sets the initial and default values of the configuration file related
 *   variables.
 *   Parameters: none
 *   Error codes: none
 */

static void init_config_vars()
{

  /* initialize the file lists */
  BlockFiles.fp    = NULL;
  DatabaseFiles.fp = NULL;
  FrequencyFile.fp = NULL;
  SequenceFiles.fp = NULL;

  BlockFiles.db_type    = -1;
  DatabaseFiles.db_type = -1;
  FrequencyFile.db_type = -1;
  SequenceFiles.db_type = -1;
  
  BlockFiles.seq_type    = UNKNOWN_SEQ;
  DatabaseFiles.seq_type = UNKNOWN_SEQ;
  FrequencyFile.seq_type = UNKNOWN_SEQ;
  SequenceFiles.seq_type = UNKNOWN_SEQ;
  
  BlockFiles.num_files    = 0;
  DatabaseFiles.num_files = 0;
  FrequencyFile.num_files = 0;
  SequenceFiles.num_files = 0;
  
  BlockFiles.max_files    = 0;
  DatabaseFiles.max_files = 0;
  FrequencyFile.max_files = 0;
  SequenceFiles.max_files = 0;
  
  BlockFiles.cur_file    = 0;
  DatabaseFiles.cur_file = 0;
  FrequencyFile.cur_file = 0;
  SequenceFiles.cur_file = 0;

  /* initialize some of the other variables to their default values */
  GeneticCodeInitializer = 0;		/* the standard genetic code */

  StrandsToSearch = 2;		/* == 2 if want to search both strands */
  RepeatsAllowed  = TRUE;	/*  RE  */
  SavedScoresFlag = FALSE;	/*  SV  */
  NumberToReport  = 0;		/* <0 means all, 0 means judge, */
				/* >0 means use that number */
  SearchType      = SEARCH_TYPE_UNSET;
  SequenceType	= UNKNOWN_SEQ;		/* sequence type */
  DoHistogram = FALSE;

  BlockToMatrixConversionMethod = 3; /* default method is pseudo-counts */
  SequenceMatrixScoringMethod   = 0; /* default method is zero */

  ErrorLevelReport = WARNING_ERR_LVL; /* report all errors */

}


/*
 * read_config_file
 *   reads the configuration file given by the filename
 *   Parameters:
 *     char * filename  : the configuration filename
 *   Return codes: 1 (TRUE) if read file sucessfully
 *                 0 (FALSE) if had problems, but without a major error
 *   Error codes:
 */

int read_config_file(filename) 
     char *filename;
{
  FILE *cfp;  /* the configuration file pointer */
  char *buf_token;
  int str_length;


  /* initialize the variabls and set the defaults */
  init_config_vars();


  /* create the error file by removing the config file extension and adding *
   * the suffix .err .  If there is no extension, just add .err to the name */
  /* NOTE: this is written to be un*x specific */
  
  strncpy(Buffer, filename, LARGE_BUFF_LENGTH);
  str_length = strlen(Buffer);
  str_length--;
  while ((Buffer[str_length] != '.') &&
	 (Buffer[str_length] != '/')) {
    str_length--;
  }
  if (Buffer[str_length] == '.') { /* there is an extension */
    strncpy(&(Buffer[str_length]), ".err\0", 
	    LARGE_BUFF_LENGTH - str_length);
  }
  else {			/* there was no extension */
    strncpy(&(Buffer[strlen(Buffer)]), ".err\0", 
	    LARGE_BUFF_LENGTH - strlen(Buffer));
  }
  set_error_file_name(Buffer);

  /* open the configuration file for reading */
  cfp = fopen(filename, "r");

  if (cfp == NULL) {
    sprintf(ErrorBuffer, 
	    "Configuration file \"%s\" not found.  Exiting\n", 
	    filename);
    ErrorReport(FATAL_ERR_LVL);	/* exit with an error */
  }


  /* read in the parameters */
  /* loop until EOF or '//' (detected in the switch) */
  while (!feof(cfp) && fgets(Buffer, LARGE_BUFF_LENGTH, cfp)) {
    if (!blank_line(Buffer)) { /* don't want to look if this is a blank line */

      buf_token = get_token(Buffer);

      /* huge switch statement of chars to determine which config entry */
      switch ( toupper(buf_token[0]) ) {
      case 'B':
	switch ( toupper(buf_token[1]) ) {
	case 'L': /* BL: Block filename */
	  BL_block_filename();
	  break;
	default:
	  unknown_conf_key(buf_token);
	  break;
	}
	break;
      case 'C':
	switch ( toupper(buf_token[1]) ) {
	case 'O': /* CO: Conversion method for converting a block to a */
		  /*     matrix */
	  CO_conversion_method();
	  break;
	default:
	  unknown_conf_key(buf_token);
	  break;
	}
	break;
      case 'D':
	switch ( toupper(buf_token[1]) ) {
	case 'B': /* DB: Database filename */
	  DB_db_filename();
	  break;
	default:
	  unknown_conf_key(buf_token);
	  break;
	}
	break;
      case 'E':
	switch ( toupper(buf_token[1]) ) {
	case 'R': /* ER: Error level to report */
	  ER_error_level_to_report();
	  break;
	case 'X': /* EX: Export Matrix filename */
	  EX_export_matr_filename();
	  break;
	default:
	  unknown_conf_key(buf_token);
	  break;
	}
	break;
      case 'F':
	switch ( toupper(buf_token[1]) ) {
	case 'R': /* FR: Codon Frequency Filename */
	  FR_freq_filename();
	  break;
	default:
	  unknown_conf_key(buf_token);
	  break;
	}
	break;
      case 'G':
	switch ( toupper(buf_token[1]) ) {
	case 'E': /* GE: Genetic code */
	  GE_genetic_code();
	  break;
	default:
	  unknown_conf_key(buf_token);
	  break;
	}
	break;
      case 'H':
	switch ( toupper(buf_token[1]) ) {
	case 'I': /* HI: Histogram */
	  HI_histogram();
	  break;
	default:
	  unknown_conf_key(buf_token);
	  break;
	}
	break;
      case 'M':
	switch ( toupper(buf_token[1]) ) {
	case 'A': /* MA: Matrix database filename*/
	  MA_Matrix_db_filename();
	  break;
	default:
	  unknown_conf_key(buf_token);
	  break;
	}
	break;
      case 'N':
	switch ( toupper(buf_token[1]) ) {
	case 'U':	/* NU: Number to report */
	  NU_number_to_report();
	  break;
	default:
	  unknown_conf_key(buf_token);
	  break;
	}
	break;
      case 'O':
	switch ( toupper(buf_token[1]) ) {
	case 'P': /* OP: Options for algorithms */
	  OP_options_for_algs();
	  break;
	case 'U': /* OU: Output filename */
	  OU_output_filename();
	  break;
	default:
	  unknown_conf_key(buf_token);
	  break;
	}
	break;
      case 'P':
	switch ( toupper(buf_token[1]) ) {
	case 'A': /* PA: The pattern file(s) */
	  PA_pattern_filename();
	  break;
	case '4': /* P4: Parameter, for backward compatibility */
	  sprintf(ErrorBuffer,
		  "Old style configuration file.");
	  ErrorReport(INFO_ERR_LVL);
	  sprintf(ErrorBuffer,
		  "New style uses the ST field instead of P4.\n");
	  ErrorReport(INFO_ERR_LVL);
	  ST_strands_to_search();
	  break;
	case '5': /* P5: Parameter, for backward compatibility */
	  sprintf(ErrorBuffer,
		  "Old style configuration file.");
	  ErrorReport(INFO_ERR_LVL);
	  sprintf(ErrorBuffer,
		  "New style uses the RE field instead of P5.\n");
	  ErrorReport(INFO_ERR_LVL);
	  RE_repeats_allowed();
	  break;
	case '9': /* P9: Parameter, for backward compatibility */
	  sprintf(ErrorBuffer,
		  "Old style configuration file.");
	  ErrorReport(INFO_ERR_LVL);
	  sprintf(ErrorBuffer,
		  "New style uses the NU field instead of P9.\n");
	  ErrorReport(INFO_ERR_LVL);
	  NU_number_to_report();
	  break;
	default:
	  unknown_conf_key(buf_token);
	  break;
	}
	break;
      case 'R':
	switch ( toupper(buf_token[1]) ) {
	case 'E':	/* RE: Repeats allowed */
	  RE_repeats_allowed();
	  break;
	default:
	  unknown_conf_key(buf_token);
	  break;
	}
	break;
      case 'S':
	switch ( toupper(buf_token[1]) ) {
	case 'C': /* SC: Scoring method */
	  SC_scoring_method();
	  break;
	case 'E': /* SE: Search type */
	  SE_search_type();
	  break;
	case 'Q': /* SQ: Sequence filename */
	  SQ_sequence_filename();
	  break;
	case 'T': /* ST: Strands to search */
	  ST_strands_to_search();
	  break;
	case 'V': /* SV: Sequence filename */
	  SV_saved_scores();
	  break;
	default:
	  unknown_conf_key(buf_token);
	  break;
	}
	break;
      case 'T':
	switch ( toupper(buf_token[1]) ) {
	case 'Y': /* TY: Sequence type  */
	  TY_sequence_type();
	  break;
	default:
	  unknown_conf_key(buf_token);
	  break;
	}
      case '/':
	switch ( toupper(buf_token[1]) ) {
	case '/': /* //: end of configuration file */
	  /* close the configuration file */
	  fclose(cfp);
	  
	  read_config_file_end();
	  
	  /* return */
	  return TRUE;
	  break;
	default:
	  unknown_conf_key(buf_token);
	  break;
	}
	break;
      default:
	unknown_conf_key(buf_token);
	break;
      } /* end switch */
    } /* end if not a blank line */
  } /* end while loop */

  /* close the configuration file */
  fclose(cfp);

  read_config_file_end();

  /* return */
  return FALSE; /* didn't read a "//", hit the end of the file */
}


/*
 * read_config_file_end
 *   Does the final clean up of reading the configuration file.  Makes
 *   sure the database files are in the right places.
 *   Parameters: none
 *   Error codes: none
 */

static void read_config_file_end()
{
  FILE *tfp;			/* temp file pointer */

  /* move the database to the correct list */

  /* NOTE: In the future there will be DB, SQ, BL, and matrix fields.  The */
  /*       if statements will be very bad.  Mayb have a number and each time */
  /*       one of the fields is seen, set a flag in the number.  Then have */
  /*       a case statement based off of that number. */
  /* SOLUTION: The matricies are stored in the blocks files until after this */
  /*           phase.  After we figure out the search type then determine if */
  /*           the files are blocks or matricies.  Since there can be only */
  /*           blocks or matricies, not both, we can do this. */

  switch (SearchType) {
  case SEARCH_TYPE_UNSET :

    /* determine the type of search. */
    /* SQ vs DB = BLOCK; BL vs DB = MATRIX; SQ vs BL = UNKNOWN; */

    /* if there are database files, should be ok */
    if (number_of_files(DATABASE_FILES) != 0) {

      /* if no sequence files, probably BL vs DB = SEARCH_TYPE_MATRIX */
      if (number_of_files(SEQUENCE_FILES) == 0) {

	/* make sure there is a block to use, */
	/* if not, there is no way to search */
	if (number_of_files(BLOCK_FILES) != 0) {
	  SearchType = SEARCH_TYPE_MATRIX;
	  SequenceFiles = DatabaseFiles;
	}
	else {
	  sprintf(ErrorBuffer,
		  "No block file (BL) or sequence file (SQ) given, aborting.\n");
	  ErrorReport(FATAL_ERR_LVL);
	}

      }

      /* since there is sequence files, if no block files, */
      /* probably SQ vs DB = SEARCH_TYPE_BLOCK */
      else if (number_of_files(BLOCK_FILES) == 0) {

	/* make sure there is a sequence file to use, */
	/* if not, there is no way to search */
	if (number_of_files(SEQUENCE_FILES) != 0) {
	  SearchType = SEARCH_TYPE_BLOCK;
	  BlockFiles = DatabaseFiles;
	}
	else {
	  /* won't be reached, left incase of re-arrangement */
	  sprintf(ErrorBuffer,
		  "No sequence file (SQ) or block file (BL) given, aborting.\n");
	  ErrorReport(FATAL_ERR_LVL);
	}

      }

      /* since there is all three kinds of files (SQ, BL and DB), guess that */
      /* the search will be UNKNOWN */
      else {
	sprintf(ErrorBuffer,
		"Configuration file has SQ, BL, and DB entries.");
	ErrorReport(WARNING_ERR_LVL);
	sprintf(ErrorBuffer,
		"Unable to determine the type of search.");
	ErrorReport(WARNING_ERR_LVL);
	sprintf(ErrorBuffer, 
		"Assuming SQ vs BL search.  Setting SearchType to UNKNOWN.\n");
	ErrorReport(WARNING_ERR_LVL);
	SearchType = SEARCH_TYPE_UNKNOWN;
      }

    }

    /* since there is no database file (no DB line) if there are both */
    /* sequence and block files we can try an UNKNOWN search */
    else {

      /* if either one of the sequence files or block files are missing, */
      /* the search cannot be done */
      if ((number_of_files(SEQUENCE_FILES) == 0) ||
	  (number_of_files(BLOCK_FILES) == 0)) {
	/* missing one of BL and SQ */
	sprintf(ErrorBuffer,
		"Missing either the BL or SQ field, aborting.\n");
	ErrorReport(FATAL_ERR_LVL);
      }

      /* both sequence and block files are present (but no database files), */
      /* we can try an UNKNOWN search */
      else {
	/* has both SQ and BL, make SEARCH_TYPE_UNKNOWN */
	sprintf(ErrorBuffer,
		"Configuration file has SQ and BL entries.");
	ErrorReport(WARNING_ERR_LVL);
	sprintf(ErrorBuffer,
		"Unable to determine the type of search.");
	ErrorReport(WARNING_ERR_LVL);
	sprintf(ErrorBuffer, 
		"Assuming SQ vs BL search.  Setting SearchType to UNKNOWN.\n");
	ErrorReport(WARNING_ERR_LVL);
	SearchType = SEARCH_TYPE_UNKNOWN;
      }

    } /* end if there are database files */

    break;
  case SEARCH_TYPE_MATRIX :
    /* blocks vs databases (sequences) */
    /* check if there are files in the database files */
    if (number_of_files(DATABASE_FILES) == 0) {
      if (number_of_files(SEQUENCE_FILES) == 0) {
	sprintf(ErrorBuffer,
		"No sequence files given for the sequence database, aborting.\n");
	ErrorReport(FATAL_ERR_LVL);
      }
      else {
	sprintf(ErrorBuffer,
		"No DB field given, assuming data is in the SQ field.\n");
	ErrorReport(WARNING_ERR_LVL);
      }
    }
    else {
      if (number_of_files(BLOCK_FILES) == 0) {
	sprintf(ErrorBuffer,
		"No block file given (BL field), aborting.\n");
	ErrorReport(FATAL_ERR_LVL);
      }
      else {
	SequenceFiles = DatabaseFiles;
      }
    }
    break;
  case SEARCH_TYPE_BLOCK :
    /* sequences vs databases (blocks) */
    /* check if there are files in the database files */
    if (number_of_files(DATABASE_FILES) == 0) {
      if (number_of_files(BLOCK_FILES) == 0) {
	sprintf(ErrorBuffer,
		"No block files given for the block database, aborting.\n");
	ErrorReport(FATAL_ERR_LVL);
      }
      else {
	sprintf(ErrorBuffer,
		"No DB field given, assuming data is in the BL field.\n");
	ErrorReport(WARNING_ERR_LVL);
      }
    }
    else {
      if (number_of_files(SEQUENCE_FILES) == 0) {
	sprintf(ErrorBuffer,
		"No sequence file given (SQ field), aborting.\n");
	ErrorReport(FATAL_ERR_LVL);
      }
      else {
	BlockFiles = DatabaseFiles;
      }
    }
    break;
  case SEARCH_TYPE_UNKNOWN :
    /* unknown search type, assuming sequences vs blocks is OK */
    sprintf(ErrorBuffer, 
	    "Unknown search type given, assuming SQ vs BL search.\n");
    ErrorReport(WARNING_ERR_LVL);
    break;
  default:
    /* BIG ERROR, should never get here.  Memory corrupted */
    sprintf(ErrorBuffer, 
	    "read_config_file_end(): SearchType variable corrupted.");
    ErrorReport(PROGRAM_ERR_LVL);
    sprintf(ErrorBuffer,
	    "                        Should never have reached this part of the program.");
    ErrorReport(PROGRAM_ERR_LVL);
    sprintf(ErrorBuffer,
	    "                        Possible memory problems.");
    ErrorReport(PROGRAM_ERR_LVL);
    sprintf(ErrorBuffer,
	    "                        The rest of the run may be invalid.");
    ErrorReport(PROGRAM_ERR_LVL);
    sprintf(ErrorBuffer, 
	    "                        Assuming SQ vs BL search.  Setting SearchType to UNKNOWN.\n");
    ErrorReport(PROGRAM_ERR_LVL);
    SearchType = SEARCH_TYPE_UNKNOWN;
    break;
  }

  /* determine if there are blocks or matricies in the BlockFiles */
  SiteSpecificScoringMatrixType = SSSM_BLOCK;
  if (MA_field_seen) {		/* the user says it is matricies */
    MatrixFiles = BlockFiles;
    SiteSpecificScoringMatrixType = SSSM_PRECOMP_MAT;
  }
  else if (! BL_field_seen) {	/* the user didn't say it's blocks so check */
    rewind_file(BLOCK_FILES);
    tfp = get_file(BLOCK_FILES);
    if (! (tfp == NULL)) {
      /* there is data in the files, so check the type */
      /* search for the MA line in a matrix */
      while (fgets(Buffer, LARGE_BUFF_LENGTH, tfp) &&
	     !(((Buffer[0] == 'M') && (Buffer[1] == 'A')) ||
	       ((Buffer[0] == 'B') && (Buffer[1] == 'L'))) &&
	     !feof(tfp));
      if ((Buffer[0] == 'M') && (Buffer[1] == 'A')) {
	MatrixFiles = BlockFiles;
	SiteSpecificScoringMatrixType = SSSM_PRECOMP_MAT;
      } /* else leave it as a blocks database */
    }
    /* if no data at all let it default to blocks and let the rest handle it */
    /* The system needs these to be rewound to start */
    rewind_file(BLOCK_FILES);	
    /* being safe, the fp was copied from BlockFiles */
    rewind_file(MATRIX_FILES);
  }
    


  /* determine sequence database db_type and seq_type */
  /* NOTE: doing direct accessing because need to set various fields */
  SequenceFiles.fp = fopen(SequenceFiles.file_names[0], "r");
  
  if (SequenceFiles.fp != NULL) {
    SequenceFiles.db_type = type_dbs(SequenceFiles.fp, DbInfo);
    if (SequenceFiles.db_type < 0 || SequenceFiles.db_type >= MAXDB) {
      /* unknown database type */
      sprintf(ErrorBuffer, 
	      "Unknown sequence format.  Known formats are:");
      ErrorReport(WARNING_ERR_LVL);
      sprintf(ErrorBuffer,
	      "Fasta, GenBank, PIR, Swiss Prot, and GCG.");
      ErrorReport(WARNING_ERR_LVL);
      sprintf(ErrorBuffer,
	      "Assuming there is no title and just sequence.\n");
      ErrorReport(WARNING_ERR_LVL);
      SequenceFiles.db_type = FLAT;
    }
    SequenceFiles.seq_type = seq_type_dbs(SequenceFiles.fp, DbInfo, 
			  SequenceFiles.db_type, SequenceType);
  }
  else {
    sprintf(ErrorBuffer, "Unable to read sequence file: %s ",
	    SequenceFiles.file_names[0]);
    ErrorReport(SERIOUS_ERR_LVL);
    sprintf(ErrorBuffer, "to determine the type of database\n");
    ErrorReport(SERIOUS_ERR_LVL);
  }

  fclose(SequenceFiles.fp);
  SequenceFiles.fp = NULL;	/* NOTE: make sure to set to NULL.  This */
				/* field is checked in various file */
				/* functions and is assumed to be NULL at */
				/* the first time called */
}

/* Change log information follows. 
  Changes since version 3.2.5:
   1/23/99 Added SV option to activate SavedScores
           Added lists.h
  Changes since version 3.2.4:
  12/12/98 Modified TY to differentiate between values DNA & AA.
           Modified ST to accept value = -1 => minus strand only.
  Changes since version 3.0.0:
   4/22/96 Added TYpe option to specify sequence type as DNA or protein.
  Changes for version 3.0.0:
   9/16/95 Changed default conversion method from 2 to 3, line 342.  JGH
*/
