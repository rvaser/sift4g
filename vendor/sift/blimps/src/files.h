/* (C) Copyright 1993, Fred Hutchinson Cancer Research Center */
/* Use, modification or distribution of these programs is subject to */
/* the terms of the non-commercial licensing agreement in license.h. */

/* files.h: access to configuration and database files */
/* Written by: Bill Alford */
/* Change log information is at the end of the file. */

#ifndef FILES_H_
#define FILES_H_

/*
 * Exported variables and data structures
 */

#define BLOCK_FILES    0
#define MATRIX_FILES   1
#define DATABASE_FILES 2
#define FREQUENCY_FILE 3
#define SEQUENCE_FILES 4
#define PATTERN_FILES 5



extern char ExportMatrixFile[SMALL_BUFF_LENGTH];
extern char OutputFile[SMALL_BUFF_LENGTH];

/* NOTE: The file lists are here only so that config.c has access to */
/*       the file_list struct.  When reading in the config file it is */
/*       necessary to handle the file lists directly.  The code should */
/*       be rewritten so that this is not necessary. */

struct file_list {
  FILE *fp;			/* pointer to the currently open file */
  int db_type;			/* the type of database this is */
  int seq_type;			/* the sequence type of the database, if */
				/* this is a seq. database. */
  int num_files;		/* the number of files in the list */
  int max_files;		/* the current maximum number of files */
				/* allowed in the list.  Can grow. */
  int cur_file;			/* the current position in the file list */
  char **file_names;		/* the array of strings of the file names */
};
typedef struct file_list FileList;

extern FileList BlockFiles;
extern FileList MatrixFiles;
extern FileList DatabaseFiles;
extern FileList FrequencyFile;	/* Note: should be only one in the list */
extern FileList SequenceFiles; 
extern FileList PatternFiles; 




/*
 * Exported functions
 */


/*
 * insert_file
 *   insert_file inserts the file_name into the file_group list
 *   Parameters:
 *     char *file_name: the name of the file
 *     int file_group: the group to get a file pointer for.
 *   Return codes: none
 *   Error codes:
 */

extern void insert_file();


/*
 * get_file
 *   get_file retrieves the file pointer of the requested group.
 *   Parameters:
 *     int file_group: the group to get a file pointer for.
 *   Return codes: Returns the file pointer to the requested file group
 *   Error codes:
 */

extern FILE *get_file();


/*
 * rewind_file
 *   rewind_file sets/resets the file pointer of the specified group back to
 *   the first file, first character.
 *   Parameters:
 *     int file_group: the group rewind.
 *   Return codes: none
 *   Error codes:
 */

extern void rewind_file();


/*
 * get_current_file_name
 *   Returns a pointer to the string of the current file name.
 *   Parameters:
 *     int file_group: the group rewind.
 *   Return codes: a string of the filename
 *   Error codes:  NULL if there is no filename
 */

extern char *get_current_file_name();


/*
 * get_file_name
 *   Returns a string of the requested file name.  (this is a copy)
 *   Parameters:
 *     int file_num:   which file to return (0..num_files);
 *     int file_group: the file group for the name
 *   Return codes: a string of the filename
 *   Error codes:  NULL if there is no filename
 */

char *get_file_name();


/*
 * number_of_files
 *   Returns the number of files in the file group.
 *   Parameters: 
 *     int file_group: the group rewind.
 *   Return codes: the number of files in the file group.
 *   Error codes:  -1 if the file group is unknown
 */

extern int number_of_files();

extern void file_assign();

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

extern int get_sequence_db_seq_type();


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

extern int get_sequence_db_db_type();




/*
 * Other functions
 */


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

extern void close_file();




#endif /*  FILES_H_ */

/* Change log information follows. */
/* Changes since version 3.0.0:
*/
