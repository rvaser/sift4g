/* (C) Copyright 1993, Fred Hutchinson Cancer Research Center */
/* Use, modification or distribution of these programs is subject to */
/* the terms of the non-commercial licensing agreement in license.h. */

/* files.c: access to configuration and database files */
/* Written by: Bill Alford */
/* Change log information is at the end of the file. */


/*	system headers not in global.h  */
/*	blimps library headers */
#include <global.h>
#include <files.h>
#include <sequences.h>
#include <string.h>

/*
 * Exported variables and data structures definitions
 */

char ExportMatrixFile[SMALL_BUFF_LENGTH];
char OutputFile[SMALL_BUFF_LENGTH];
char ErrorFile[SMALL_BUFF_LENGTH];

FileList BlockFiles;
FileList MatrixFiles;
FileList DatabaseFiles;
FileList FrequencyFile;		/* Note: should be only one file in the list */
FileList SequenceFiles; 
FileList PatternFiles; 


/*
 * Local variables and data structures
 */


#define FILE_LIST_INCREASE_SIZE 20





/*
 * Function definitions
 */

static FileList *get_file_list();




/*
 * Basic File routines
 *
 *   static FileList *get_file_list(file_group)
 *   void insert_file(file_name, file_group)
 *   FILE *get_file(file_group)
 *   void rewind_file(file_group)
 *   char *get_current_file_name(file_group) 
 *   char *get_file_name(file_num, file_group) 
 *   int number_of_files(file_group) 
 *   void close_file(file_group)
 *   int type_dbs(fin,dbs)
 *   int get_sequence_db_seq_type()
 *   int get_sequence_db_db_type()
 *
 */


/* 
 * get_file_list
 *   returns a pointer to the file list of the requested file group.
 *   Parameters:
 *     int file_group: the group to get a file pointer for.
 *   Retrun codes: the pointer to the file list
 *   Error codes:  NULL if the requested file group is unknown
 */

static FileList *get_file_list(file_group)
     int file_group;
{
  FileList *this_list;

  switch (file_group) {
  case BLOCK_FILES :
    this_list = &(BlockFiles);
    break;
  case MATRIX_FILES :
    this_list = &(MatrixFiles);
    break;
  case DATABASE_FILES :
    this_list = &(DatabaseFiles);
    break;
  case FREQUENCY_FILE :		/* Note: should be only one in the list */
    this_list = &(FrequencyFile);
    break;
  case SEQUENCE_FILES :
    this_list = &(SequenceFiles);
    break;
  case PATTERN_FILES :
    this_list = &(PatternFiles);
    break;
  default:
    sprintf(ErrorBuffer, 
	    "get_file_list(): Unknown file group %d.", file_group);
    ErrorReport(PROGRAM_ERR_LVL);
    this_list = NULL;
  } /* end switch */

  return this_list;
}


/*
 * insert_file
 *   insert_file inserts the file_name into the file_group list
 *   Parameters:
 *     char *file_name: the name of the file
 *     int file_group: the group to get a file pointer for.
 *   Return codes: none
 *   Error codes:
 */

void insert_file(file_name, file_group)
     char *file_name;
     int file_group;
{
  FileList *this_list;

  this_list = get_file_list(file_group);

  if (this_list == NULL) {
    sprintf(ErrorBuffer, 
	    "insert_file(): Bad file list, not entering filename \"%s\"",
	    file_name);
    ErrorReport(PROGRAM_ERR_LVL);
    sprintf(ErrorBuffer,
	    "               into a file list\n");
    ErrorReport(PROGRAM_ERR_LVL);
  }

  /* if this is the first time (this_list->num_files == 0) allocate space */
  if (this_list->num_files == 0) {
    CheckMem(
	     this_list->file_names = (char **) calloc(FILE_LIST_INCREASE_SIZE,
						      sizeof(char *))
	     );
    this_list->max_files = FILE_LIST_INCREASE_SIZE;
  }

  /* increase file number count */
  this_list->num_files++;

  /* if this is a frequency file list and there is already a file, */
  /* warn that extras will be ignored */
  if ((file_group == FREQUENCY_FILE) && (this_list->num_files > 1)) {
    sprintf(ErrorBuffer, "Already have a frequency file: %s",
	   this_list->file_names[0]);
    ErrorReport(WARNING_ERR_LVL);
    sprintf(ErrorBuffer, "The file %s will be ignored.\n", file_name);
    ErrorReport(WARNING_ERR_LVL);
  }

  /* if there is no room, allocate more space */
  if (this_list->num_files > this_list->max_files) {
    char **tmp_ptr;		/* don't want to ruin pointer in case there */
				/* is a repeated realloc call */
    
    CheckMem(
	     tmp_ptr = (char **) realloc(this_list->file_names,
					 (this_list->max_files +
					  FILE_LIST_INCREASE_SIZE) *
					 sizeof(char *))
	     );
    this_list->file_names = tmp_ptr;
    this_list->max_files += FILE_LIST_INCREASE_SIZE;
  }

/*#ifndef NO_STRDUP
  this_list->file_names[this_list->num_files - 1] = strdup(file_name);
#else*/
/* gsims: strdup isn't an -ansi function */
  CheckMem(
          this_list->file_names[this_list->num_files - 1] = (char *)
            calloc(strlen(file_name) + 1, sizeof(char))	/* +1 for '\0' */
          );
  strncpy(this_list->file_names[this_list->num_files - 1], file_name, 
         strlen(file_name));
  this_list->file_names[this_list->num_files - 1][strlen(file_name)+1] = '\0'; 
/*#endif*/

}

/*
 * get_file
 *   get_file retrieves the file pointer of the requested group.
 *   Parameters:
 *     int file_group: the group to get a file pointer for.
 *   Return codes: Returns the file pointer to the requested file group
 *   Error codes:
 */

FILE *get_file(file_group) 
     int file_group;
{
  FileList *this_list;
  register char c;

  this_list = get_file_list(file_group);

  if (this_list == NULL) {
    sprintf(ErrorBuffer, "get_file(): Bad file list, unable to open a file\n");
    ErrorReport(PROGRAM_ERR_LVL);
    return NULL;
  } /* end switch */
  
  if (this_list->fp == NULL) {  /* initially set to NULL in declaration */
    this_list->fp = fopen(this_list->file_names[this_list->cur_file], "r");
    /* the sequence info was set at the end of the reading of the */
    /* configuration file and when the file group is rewound.  At this */
    /* time db_type and seq_type are set */
  }
  if (this_list->fp == NULL) {
    /* report error */
    switch (file_group) {
    case BLOCK_FILES :
      sprintf(ErrorBuffer, "Unable to open block file: %s\n", 
	      this_list->file_names[this_list->cur_file]);
      ErrorReport(SERIOUS_ERR_LVL);
      return NULL;
      break;
    case MATRIX_FILES :
      sprintf(ErrorBuffer, "Unable to open matrix file: %s\n", 
	      this_list->file_names[this_list->cur_file]);
      ErrorReport(SERIOUS_ERR_LVL);
      return NULL;
      break;
    case DATABASE_FILES :
      sprintf(ErrorBuffer, "Unable to open database file: %s\n", 
	      this_list->file_names[this_list->cur_file]);
      ErrorReport(SERIOUS_ERR_LVL);
      return NULL;
      break;
    case FREQUENCY_FILE :	/* Note: should be only one in the list */
      sprintf(ErrorBuffer, "Unable to open frequency file: %s\n", 
	      this_list->file_names[this_list->cur_file]);
      ErrorReport(SERIOUS_ERR_LVL);
      return NULL;
      break;
    case SEQUENCE_FILES :
      sprintf(ErrorBuffer, "Unable to read sequence file: %s\n", 
	      this_list->file_names[this_list->cur_file]);
      ErrorReport(SERIOUS_ERR_LVL);
      return NULL;
      break;
    case PATTERN_FILES :
      sprintf(ErrorBuffer, "Unable to read pattern file: %s\n", 
	      this_list->file_names[this_list->cur_file]);
      ErrorReport(SERIOUS_ERR_LVL);
      return NULL;
      break;
    default:
      sprintf(ErrorBuffer, "get_file(): Unknown file group");
      ErrorReport(PROGRAM_ERR_LVL);
      sprintf(ErrorBuffer, "            Unable to open a file\n");
      ErrorReport(PROGRAM_ERR_LVL);
    }
  }
  
  /* HACK!!! since fgets stops reading at a newline, it doesn't see the EOF */
  /* and does not set it.  This is a way to make sure that the EOF is raised */
  c = getc(this_list->fp);
  ungetc(c, this_list->fp);

  if (feof(this_list->fp)) {
    /* if more files to look at */
    if (this_list->cur_file < this_list->num_files -1) { /* -1 because */
			/* cur_files starts at zero, num_files starts at one */
      this_list->cur_file++;
      fclose(this_list->fp); 
      this_list->fp = fopen(this_list->file_names[this_list->cur_file], "r");

  /*-------------------------------------------------------------------*/
  /* determine sequence database db_type and seq_type */
  /* Code take from config.c   */
  /* NOTE: doing direct accessing because need to set various fields */
  
      if (file_group == SEQUENCE_FILES)
      {
  if (this_list->fp != NULL) {
    this_list->db_type = type_dbs(this_list->fp, DbInfo);
    if (this_list->db_type < 0 || this_list->db_type >= MAXDB) {
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
      this_list->db_type = FLAT;
    }
    this_list->seq_type = seq_type_dbs(this_list->fp, DbInfo, 
			  this_list->db_type, this_list);
  }
  else {
    sprintf(ErrorBuffer, "Unable to read sequence file: %s ",
	    this_list->file_names[0]);
    ErrorReport(SERIOUS_ERR_LVL);
    sprintf(ErrorBuffer, "to determine the type of database\n");
    ErrorReport(SERIOUS_ERR_LVL);
  }
     }  /* end of SEQUENCE_FILES  */
  /*-------------------------------------------------------------------*/


      if (this_list->fp == NULL) {
	/* report error */
	switch (file_group) {
	case BLOCK_FILES :
	  sprintf(ErrorBuffer, "Unable to open block file: %s\n", 
		  this_list->file_names[this_list->cur_file]);
	  ErrorReport(SERIOUS_ERR_LVL);
	  return NULL;
	  break;
	case MATRIX_FILES :
	  sprintf(ErrorBuffer, "Unable to open matrix file: %s\n", 
		  this_list->file_names[this_list->cur_file]);
	  ErrorReport(SERIOUS_ERR_LVL);
	  return NULL;
	  break;
	case DATABASE_FILES :
	  sprintf(ErrorBuffer, "Unable to open database file: %s\n", 
		  this_list->file_names[this_list->cur_file]);
	  ErrorReport(SERIOUS_ERR_LVL);
	  return NULL;
	  break;
	case FREQUENCY_FILE :	/* Note: should be only one in the list */
	  sprintf(ErrorBuffer, "Unable to open frequency file: %s\n", 
		  this_list->file_names[this_list->cur_file]);
	  ErrorReport(SERIOUS_ERR_LVL);
	  return NULL;
	  break;
	case SEQUENCE_FILES :
	  sprintf(ErrorBuffer, "Unable to read sequence file: %s", 
		  this_list->file_names[this_list->cur_file]);
	  ErrorReport(SERIOUS_ERR_LVL);
	  return NULL;
	  break;
	case PATTERN_FILES :
	  sprintf(ErrorBuffer, "Unable to read pattern file: %s", 
		  this_list->file_names[this_list->cur_file]);
	  ErrorReport(SERIOUS_ERR_LVL);
	  return NULL;
	  break;
	default:
	  sprintf(ErrorBuffer, "get_file(): Unknown file group");
	  ErrorReport(PROGRAM_ERR_LVL);
	  sprintf(ErrorBuffer, "            Unable to open a file\n");
	  ErrorReport(PROGRAM_ERR_LVL);
	}
      }
    }
    else { /* no more files, close and exit */
      fclose(this_list->fp); 
      return NULL;
    }
  }
  return (this_list->fp);
  
}


/*
 * rewind_file
 *   rewind_file sets/resets the file pointer of the specified group back to
 *   the first file, first character.
 *   Parameters:
 *     int file_group: the group rewind.
 *   Return codes: none
 *   Error codes:
 */

void rewind_file(file_group) 
     int file_group;
{
  FileList *this_list;

  this_list = get_file_list(file_group);

  if (this_list == NULL) {
    sprintf(ErrorBuffer, 
	    "rewind_file(): Bad file group, unable to rewind file list\n");
    ErrorReport(PROGRAM_ERR_LVL);
  }
  
  if (this_list->fp != NULL) {
    fclose(this_list->fp); 
  }

  /* set the current file to the first file */
  this_list->cur_file = 0;

  /* If want to check the sequence database type and sequence */
  /* type for each file, put the code here.  The database type and */
  /* the sequence type are set at the end of reading the */
  /* configuration file. */ 

  /* set the file pointer to null so that the next call to get_file will */
  /* start at the beginning */
  this_list->fp = NULL;

}


/*
 * get_current_file_name
 *   Returns a pointer to the string of the current file name. 
 *   Parameters:
 *     int file_group: the file group for the name
 *   Return codes: a string of the filename
 *   Error codes:  NULL if there is no filename
 */

char *get_current_file_name(file_group) 
     int file_group;
{
  FileList *this_list;

  this_list = get_file_list(file_group);

  if (this_list == NULL) {
    sprintf(ErrorBuffer, 
	    "get_current_file_name(): Bad file group, unable to return a file name\n");
    ErrorReport(PROGRAM_ERR_LVL);
    return NULL;
  }
  
  if (this_list->num_files <= 0) {
    sprintf(ErrorBuffer, 
	    "get_current_file_name(): No files in the file group");
    ErrorReport(PROGRAM_ERR_LVL);
    sprintf(ErrorBuffer, 
	    "                         Unable to return a file name\n");
    ErrorReport(PROGRAM_ERR_LVL);
    return NULL;
  }
  
  /* return a pointer to the current file name */
  return (char *) this_list->file_names[this_list->cur_file];

}


/*
 * get_file_name
 *   Returns a pointer to the string of the requested file name.  
 *   Parameters:
 *     int file_num:   which file to return (0..num_files);
 *     int file_group: the file group for the name
 *   Return codes: a string of the filename
 *   Error codes:  NULL if there is no filename
 */

char *get_file_name(file_num, file_group) 
     int file_num;
     int file_group;
{
  FileList *this_list;

  this_list = get_file_list(file_group);

  if (this_list == NULL) {
    sprintf(ErrorBuffer, 
	    "get_file_name(): Bad file group, unable to return a file name\n");
    ErrorReport(PROGRAM_ERR_LVL);
    return NULL;
  }
  
  if (this_list->num_files <= 0) {
    sprintf(ErrorBuffer, "get_file_name(): No files in the file group");
    ErrorReport(PROGRAM_ERR_LVL);
    sprintf(ErrorBuffer, "                 Unable to return a file name\n");
    ErrorReport(PROGRAM_ERR_LVL);
    return NULL;
  }
  
  if (this_list->num_files < file_num) {
    sprintf(ErrorBuffer, 
	    "get_file_name(): Requested file out of the range of the files");
    ErrorReport(PROGRAM_ERR_LVL);
    sprintf(ErrorBuffer,
	    "                 in the file group");
    ErrorReport(PROGRAM_ERR_LVL);
    sprintf(ErrorBuffer, 
	    "                 Unable to return a file name\n");
    ErrorReport(PROGRAM_ERR_LVL);
    return NULL;
  }

  /* return a copy of the current file name */
  return (char *) this_list->file_names[file_num];

}


/*
 * number_of_files
 *   Returns the number of files in the file group.
 *   Parameters: 
 *     int file_group: the group rewind.
 *   Return codes: the number of files in the file group.
 *   Error codes:  -1 if the file group is unknown
 */

int number_of_files(file_group) 
     int file_group;
{
  FileList *this_list;

  this_list = get_file_list(file_group);

  if (this_list == NULL) {
    sprintf(ErrorBuffer, 
	    "number_of_files(): Bad file group, unable to report number of files\n");
    ErrorReport(PROGRAM_ERR_LVL);
    return -1;
  }

  return this_list->num_files;
}


/*
 * close_file
 *   closes the selected file group.
 *   Parameters: 
 *     int file_group: the group rewind.
 *   Return codes: none
 *   Error codes:
 *   Note: right now this function is identical to rewind_file() except it
 *         does not set the current file to the beginning.
 */

void close_file(file_group) 
     int file_group;
{
  FileList *this_list;

  this_list = get_file_list(file_group);

  if (this_list == NULL) {
    sprintf(ErrorBuffer, 
	    "close_file() : Bad file group, unable to rewind file list\n");
    ErrorReport(PROGRAM_ERR_LVL);
  }
  
  if (this_list->fp != NULL) {
    fclose(this_list->fp); 
  }
  
  /* set the file pointer to null so that the next call to get_file will */
  /* start at the beginning */
  this_list->fp = NULL;

}



/* NOTE: it might be good in the future to make these into generic */
/*       functions that are passed the type of file structure.  There */
/*       might become situations where there are sequence and database */
/*       types for blocks. */

/*
 * get_sequence_db_seq_type
 *   Returns the type of sequences of the databases.  
 *   NOTE: this assumes that all the sequences in a database file are
 *         of the same type, and that all of the database files are of
 *         the same type.
 *   NOTE: SequenceFiles.db_type is set only at the end of reading the 
 *         configuration file.
 *   Parameters: none
 *   Return codes: the different sequence types
 *   Error codes:  none
 */

int get_sequence_db_seq_type()
{
  return SequenceFiles.seq_type;    
}

/*
 * get_sequence_db_db_type
 *   Returns the kind of the database of the seqeuence database file 
 *   currently being read.  
 *   NOTE: SequenceFiles.db_type is set each time a new sequence file is
 *         opened.
 *   Parameters: none
 *   Return codes: the different sequence db types
 *   Error codes:  none
 */

int get_sequence_db_db_type()
{
  return SequenceFiles.db_type;    
}


/* Change log information follows. */
/*  Changes since version 3.0.0:
 7/11/96  Set SequenceFiles.db_type & .seq_type in get_file().
*/
