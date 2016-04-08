/* (C) Copyright 1993-2000, Fred Hutchinson Cancer Research Center */
/* Use, modification or distribution of these programs is subject to */
/* the terms of the non-commercial licensing agreement in license.h. */

/* matrix.c: matrix manipulation functions */
/* Written by: Bill Alford */
/* Change log information is at the end of the file. */

/*	system headers not in global.h   */
/*	blimps library headers */
#include <global.h>
#include <blocks.h>	/* includes sequences.h, output.h */
#include <matrix.h>	/* includes pattern.h */
#include <residues.h>

/*
 * Exported variables and data structures
 */

/*
 * Local variables and data structures
 */

#define MATRIX_LENGTH_INCREASE_SIZE 80 /* the default length of the matrix */
				       /* and the size to increase by if */
				       /* more room is needed (doubtful) */

/*
 * Function definitions
 */

static Boolean read_matrix_header();
static void read_matrix_body();
static void resize_matrix();

/*
 * read_a_matrix
 *   reads a matrix from the data base and returns a pointer to the new
 *   matrix data structure
 *   Parameters:  
 *     FILE *mfp: a pointer to the database file with the matrix.
 *   Error codes: NULL if a matrix was not read
 */

Matrix *read_a_matrix (mfp)
     FILE *mfp;			/* matrix file pointer */
{
  Matrix *matrix;


  /* get the matrix file pointer */
/*  mfp = get_file(MATRIX_FILES);*/
  if (mfp == NULL) {
    /* no more data to read into matricies */
    matrix = NULL;
    return matrix;
  }

  /* allocate space for a new matrix */
  CheckMem(
	   matrix = (Matrix *) malloc(sizeof(Matrix))
	   );

  /* set the max size of the length */
  matrix->max_length = 0;

  matrix->block = NULL;	/* there is no block */

  /* read header */
  if (read_matrix_header(mfp, matrix) == FALSE) {
    free(matrix);
    return NULL;
  }

  /* read body */
  read_matrix_body(mfp, matrix); 

  /* set the patterns to NULL */
  matrix->patterns = NULL;

  /* return matrix */
  return matrix;
}  /* end of read_a_matrix() */


/*
 * read_matrix_header
 *   Reads the header information for the matrix.
 *   Parameters:
 *     FILE *mfp:      the matrix file pointer. 
 *     Matrix *matrix: the matrix to put the data in.
 *   Error Codes: FALSE if the matrix could not be read, TRUE otherwise
 */

static Boolean read_matrix_header(mfp, matrix)
     FILE *mfp;
     Matrix *matrix;
{
  char *buf2, *buf3;

  /* scan for either ID, AC, DE, MA */
  fgets(Buffer, LARGE_BUFF_LENGTH, mfp);

  while ( !(((Buffer[0] == 'I') && (Buffer[1] == 'D')) ||  /* ID */
	    ((Buffer[0] == 'A') && (Buffer[1] == 'C')) ||  /* AC */
	    ((Buffer[0] == 'D') && (Buffer[1] == 'E')) ||  /* DE */
	    ((Buffer[0] == 'M') && (Buffer[1] == 'A')))    /* MA */
	 && !feof(mfp) ) { 
    fgets(Buffer, LARGE_BUFF_LENGTH, mfp);
  }

  /* if we reached the end of the file, then there was no matrix */
  if (feof(mfp)) {
    return FALSE;			/* it did not read the matrix */
  }

  /* read ID if it exists */
  if ((Buffer[0] == 'I') && (Buffer[1] == 'D')) {
    buf2 = &Buffer[2];		/* Eat the ID at the beginining */
    buf2 = eat_whitespace(buf2); /* get a pointer to the string after ID */
    remove_trailing_whitespace(buf2); /* remove the \n at the end of the line */
    strncpy(matrix->id, buf2, SMALL_BUFF_LENGTH); /* copy the string into the matrix entry */
    fgets(Buffer, LARGE_BUFF_LENGTH, mfp);	/* read the next line */
  }
  else { /* no ID field, raise error */
    matrix->id[0] = '\0';	/* make sure string is empty */
    sprintf(ErrorBuffer, "Error in matrix file format.  No ID line.\n");
    ErrorReport(WARNING_ERR_LVL);
  }

  /* read AC if it exists */
  if ((Buffer[0] == 'A') && (Buffer[1] == 'C')) {
    buf2 = &Buffer[2];		/* Eat the AC at the beginining */
    buf2 = eat_whitespace(buf2); /* get a pointer to the string after AC */
    remove_trailing_whitespace(buf2); /* remove the \n at the end of the line */
    strncpy(matrix->ac, buf2, SMALL_BUFF_LENGTH); /* copy the string into the matrix entry */
    buf2 = get_token(buf2);	/* get the number */
    if (buf2[strlen(buf2)-1] == ';') {
      /* remove the ';' at the end of the number */
      buf2[strlen(buf2)-1] = '\0'; 
    }
    strncpy(matrix->number, buf2, NUMBER_WIDTH); /* copy the string into the matrix entry */
    fgets(Buffer, LARGE_BUFF_LENGTH, mfp);	/* read the next line */
  }
  else { /* no AC field, raise error */
    matrix->ac[0] = '\0';	/* make sure string is empty */
    matrix->number[0] = '\0';	/* make sure string is empty */
    sprintf(ErrorBuffer, "Error in matrix file format.  No AC line.\n");
    ErrorReport(WARNING_ERR_LVL);
  }

  /* read DE if it exists */
  if ((Buffer[0] == 'D') && (Buffer[1] == 'E')) {
    buf2 = &Buffer[2];		/* Eat the DE at the beginining */
    buf2 = eat_whitespace(buf2); /* get a pointer to the string after DE */
    remove_trailing_whitespace(buf2); /* remove the \n at the end of the line */
    strncpy(matrix->de, buf2, DESC_WIDTH); /* copy the string into the matrix entry */
    fgets(Buffer, LARGE_BUFF_LENGTH, mfp); /* read the next line */
  }
  else { /* no DE field, raise error */
    matrix->de[0] = '\0';	/* make sure string is empty */
    sprintf(ErrorBuffer, "Error in matrix file format.  No DE line.\n");
    ErrorReport(WARNING_ERR_LVL);
  }
  
  /* read MA if it does not exist raise an error */
  if ((Buffer[0] == 'M') && (Buffer[1] == 'A')) {
    buf2 = &Buffer[2];		/* Eat the MA at the beginining */
    buf2 = eat_whitespace(buf2); /* get a pointer to the string after MA */
    remove_trailing_whitespace(buf2); /* remove the \n at the end of */
					     /* the line */
    strncpy(matrix->ma, buf2, SMALL_BUFF_LENGTH);	/* copy the string into the matrix entry */
    
    /* scan and process the MA line */
    /* find the motif */
    buf3 = strstr(Buffer, "motif");
    if (buf3!=NULL) {
      buf3 -= 4;		/* move back to the motif type */
      sscanf(buf3, "%s", matrix->motif);
    }
    else {
      matrix->motif[0] = '\0';
    }

    /* find the width */
    buf3 = strstr(Buffer, "width=");
    if (buf3!=NULL) {
      sscanf(buf3, "width=%d;", &(matrix->width));
    }
    else {
      /* report an error, it is serious if the width is not set */
      sprintf(ErrorBuffer, "No width field for matrix %s", matrix->number);
      ErrorReport(SERIOUS_ERR_LVL);
      sprintf(ErrorBuffer, "Setting width to zero\n");
      ErrorReport(SERIOUS_ERR_LVL);
      matrix->width = 0;
    }

    /* find the seqs */
    buf3 = strstr(Buffer, "seqs=");
    if (buf3!=NULL) {
      sscanf(buf3, "seqs=%d", &(matrix->num_sequences));
      if (matrix->num_sequences < 0) { /* incase the number is there but is <1 */
	matrix->num_sequences = 0;
      }
    }
    else {
      matrix->num_sequences = 0;
    }
    
    /* find the 99.5% */
    buf3 = strstr(Buffer, "99.5%=");
    if (buf3!=NULL) {
      sscanf(buf3, "99.5%%=%d;", &(matrix->percentile));
    }
    else {
      matrix->percentile = 0;
    }

    /* find the strength */
    buf3 = strstr(Buffer, "strength=");
    if (buf3!=NULL) {
      sscanf(buf3, "strength=%d", &(matrix->strength));
    }
    else {
      matrix->strength = 0;
    }
  }
  else { /* no MA field, raise error */
    matrix->ma[0] = '\0';	/* make sure the MA string is empty */
    matrix->motif[0] = '\0';	/* make sure the motif string is empty */
    matrix->width = 0;		/* make sure the width is small */
    matrix->percentile = 0;	
    matrix->strength = 0;	
    matrix->num_sequences = 0;	
    matrix->max_length = 0;	
    sprintf(ErrorBuffer, "Error in matrix file format.  No MA line.");
    ErrorReport(SERIOUS_ERR_LVL);
    sprintf(ErrorBuffer, "Attempting to set values to be able to continue.");
    ErrorReport(SERIOUS_ERR_LVL);
    sprintf(ErrorBuffer, "The first column of the matrix will be missed.\n");
    ErrorReport(SERIOUS_ERR_LVL);
  }

  return TRUE;
}

/*
 * read_matrix_body
 *   Reads the body information of the matrix.
 *   Parameters:
 *     FILE *mfp:      the matrix file pointer. 
 *     Matrix *matrix: the matrix to put the data in.
 *   Error Codes: needed
 */

static void read_matrix_body(mfp, matrix)
     FILE *mfp;
     Matrix *matrix;
{
  int aa, len;
  float dtemp;
  char c;
  char *num_buf;
  
  if (matrix->width <= 0) {
    matrix->max_length = MATRIX_LENGTH_INCREASE_SIZE;
  }
  else {
    matrix->max_length = matrix->width;
  }
  
  /* allocate space for the weights arrays */
    CheckMem(
	matrix->weights[0] = (MatType *) calloc(matrix->width*MATRIX_AA_WIDTH,
						sizeof(MatType))
	   );

  /* setup the array of pointers, remember weights is an array of arrays */
  for (aa=0; aa<MATRIX_AA_WIDTH; aa++) {
    matrix->weights[aa] = matrix->weights[0] + aa*matrix->width;
  }

  /*
   * read in the matrix values
   */
  
  len = 0;

  /* skip the letters */
  fgets(Buffer, LARGE_BUFF_LENGTH, mfp);

  /* read in the numbers */
  while (fgets(Buffer, LARGE_BUFF_LENGTH, mfp) &&
	 !blank_line(Buffer) &&
	 !((Buffer[0] == '/') && (Buffer[1] == '/'))) {

    if (len > matrix->max_length) {
      resize_matrix(matrix);
    }

    num_buf = get_token(Buffer);
    /*  NOTE:  Assumes matrix is in integer format !  */
    for (c='A'; c<='Z'; c++) {
      if ((c != 'J') && (c != 'O') && (c != 'U')) {
	sscanf(num_buf, "%f", &dtemp);
        matrix->weights[aa_atob[(int)c]][len] = (MatType) dtemp;
	num_buf = get_token(NULL);
	/* do not try to catch long lines by reading in more.  Do this */
	/* only if code has been added to take care of the possibility of */
	/* splitting a number in half. */
	/* no more tokens left, so get the next line */
	/*
	if (num_buf == NULL) {	
	  fgets(Buffer, LARGE_BUFF_LENGTH, mfp);
	  num_buf = get_token(Buffer);
	}
        */
      }
    }
    c = '*';
    sscanf(num_buf, "%f", &dtemp);
    matrix->weights[aa_atob[(int)c]][len] = (MatType) dtemp;
    num_buf = get_token(NULL);
    /* do not try to catch long lines by reading in more.  Do this */
    /* only if code has been added to take care of the possibility of */
    /* splitting a number in half. */
    /* no more tokens left, so get the next line */
    /*
    if (num_buf == NULL) {	
      fgets(Buffer, LARGE_BUFF_LENGTH, mfp);
      num_buf = get_token(Buffer);
    }
    */
    c = '-';
    sscanf(num_buf, "%f", &dtemp);
    matrix->weights[aa_atob[(int)c]][len] = (MatType) dtemp;

    len++;
  }
}  /* end of read_matrix_body */


/*
 * new_matrix
 *   allocates space for a matrix of length len, sets up the data structure
 *   and returns a pointer to the Matrix.
 *   Parameters:
 *     int len: the length of the sequence/matrix
 *   Error codes:
 */

Matrix *new_matrix(len) 
     int len;
{
  int aa;
  Matrix *matrix;

  /* get the space for the Matrix data structure */
  CheckMem(
	   matrix = (Matrix *) malloc(sizeof(Matrix))
	   );

  /* initialize */
  matrix->block = NULL;
  matrix->id[0] = matrix->ac[0] = matrix->de[0] = '\0';
  sprintf(matrix->ma, "width=%d;", len);
  matrix->number[0] = matrix->motif[0] = '\0';
  matrix->width = matrix->max_length = len;
  matrix->num_sequences = 0;
  matrix->percentile = matrix->strength = matrix->undefined = 0;
  matrix->patterns = NULL;
  matrix->undefined_ptr = NULL;

  /* allocate the space for the matrix weights */
  CheckMem(
	matrix->weights[0] = (MatType *) calloc(matrix->width*MATRIX_AA_WIDTH,
						sizeof(MatType))
	   );

  /* setup the array of pointers, remember weights is an array of arrays */
  for (aa=0; aa<MATRIX_AA_WIDTH; aa++) {
    matrix->weights[aa] = matrix->weights[0] + aa*matrix->width;
  }

  return matrix;
}  /* end of new_matrix */


/* 
 * resize_matrix
 *   Increases the memory for the storage of the matrix weights in a matrix.
 *   Parameter:
 *     Matrix *matrix: the matrix to resize
 *   Error codes: none
 */

static void resize_matrix(matrix)
     Matrix *matrix;
{
  /* remember to allocate new space and then to copy over. */
  /* matrix.weights is an array of arrays */
  fprintf(stderr, "resize_matrix() not finished.  Exiting.\n");
  exit(1000);
}


  
/* 
 * free_matrix
 *   Deletes the matrix and the sub elements.
 *   Parameters: 
 *     Matrix *matrix: the matrix to free
 *   Return code: none
 *   Error code: none
 */

void free_matrix(matrix) 
     Matrix *matrix;
{
  if (matrix->patterns != NULL) {
    free(matrix->patterns);
  }
  free(matrix->weights[0]);
  free(matrix);
}


/*
 * matrix_comparison
 *   Compares two matricies.   It compares by the value in the block->id
 *   field if it exists.
 *   Parameters:
 *     MatrixListEntry a, b: the entries to compare
 *   Return codes:  a return value < 0 if a < b, a return value = 0 if a == b,
 *                  and a return value > 0 if a > b
 *   Error codes: none
 */

int matrix_comparison(a, b)
     Matrix *a, *b;
{
  return (int) (strcmp(a->number, b->number));
}




/*
 * Matrix printing.
 */

/*
 * print_matrix
 *   Prints a Matrix data structure.  Primarily for debugging purposes.
 *   Parameters: 
 *     Matrix *matrix:  the matrix to print
 *   Error Codes: none
 */

void print_matrix(matrix)
     Matrix *matrix;
{
  int pos;
  int low, high;
  char c;

#define MATRIX_PRINT_WIDTH 18

  high = MATRIX_PRINT_WIDTH;
  for (low=0; 
       low<high && low<matrix->width; 
       low = (high+=MATRIX_PRINT_WIDTH) - MATRIX_PRINT_WIDTH) {



  printf("\n");

  if (matrix->width > 99) {
    printf("  |");
    for (pos=low; pos<high &&  pos<matrix->width; pos++) {
      if (pos > 99) {
	printf("% 4d", pos%100);
      }
      else {
	printf("    ");
      }
    }
  }

  printf("\n");
  if (matrix->width > 9) {
    printf("  |");
    for (pos=low; pos<high &&  pos<matrix->width; pos++) {
      if (pos > 9) {
	printf("% 4d", (pos-((pos/100)*100)) / 10);
      }
      else {
	printf("    ");
      }
    }
  }
  
  printf("\n");
  printf("  |");
  for (pos=low; pos<high &&  pos<matrix->width; pos++) {
    printf("% 4d", (pos-((pos/10)*10)));
  }
  
  printf("\n");
  printf("--+");
  for (pos=low; pos<high &&  pos<matrix->width; pos++) {
    printf("----");
  }

  printf("\n");
  
 
  for (c='A'; c<='Z'; c++) {
    printf("%c |", c);
    for (pos=low; pos<high &&  pos<matrix->width; pos++) {
      printf("% 6.4f", matrix->weights[aa_atob[(int)c]][pos]);
    }
    printf("\n");
  }
  c = '*';
  printf("%c |", c);
  for (pos=low; pos<high &&  pos<matrix->width; pos++) {
    printf("% 6.4f", matrix->weights[aa_atob[(int)c]][pos]);
  }
  printf("\n");
  c = '-';
  printf("%c |", c);
  for (pos=low; pos<high &&  pos<matrix->width; pos++) {
    printf("% 6.4f", matrix->weights[aa_atob[(int)c]][pos]);
  }
  printf("\n");

  }

}

/* 
 * output_matrix
 *   Outputs a matrix data structure to the given file.
 *   Parameters: 
 *     Matrix *matrix: the matrix to output
 *     FILE   *omfp:   the output matrix file pointer
 *   Return code: none
 *   Error code: none
 */

void output_matrix(matrix, omfp)
     Matrix *matrix;
     FILE *omfp;
{
  output_matrix_st(matrix, omfp, INT_OUTPUT, AA_SEQ);
}


/* 
 * output_matrix_s
 *   Outputs a matrix data structure to the given file with the 
 *   specified style of data.
 *   Parameters: 
 *     Matrix *matrix: the matrix to output
 *     FILE   *omfp:   the output matrix file pointer
 *     int    style:   the kind of output (INT_OUTPUT, FLOAT_OUTPUT)
 *   Return code: none
 *   Error code: none
 */

void output_matrix_s(matrix, omfp, style)
     Matrix *matrix;
     FILE *omfp;
     int style;
{
  char c;
  int l;

  if (style == INT_OUTPUT) {
    fprintf(omfp, "ID   %s\n", matrix->id);
    fprintf(omfp, "AC   %s\n", matrix->ac);
    fprintf(omfp, "DE   %s\n", matrix->de);
    fprintf(omfp, "MA   %s\n", matrix->ma);
    
    fprintf(omfp, " ");
    for (c='A'; c <= 'Z'; c++) {
      if ((c != 'J') && (c != 'O') && (c != 'U')) {
	fprintf(omfp, " %c  ", c);
      }
    }
    fprintf(omfp, " *  ");
    fprintf(omfp, " -\n");
    
    for (l=0; l<matrix->width; l++) {
      for (c='A'; c <= 'Z'; c++) {
	if ((c != 'J') && (c != 'O') && (c != 'U')) {
	  fprintf(omfp, "%3d ", (int) round(matrix->weights[aa_atob[(int)c]][l]));
	}
      }
      fprintf(omfp, "%3d ", (int) round(matrix->weights[aa_atob['*']][l]));
      fprintf(omfp, "%3d\n", (int) round(matrix->weights[aa_atob['-']][l]));
    }  
    
    fprintf(omfp, "//\n");
  }
  else if (style == FLOAT_OUTPUT) {
    fprintf(omfp, "ID   %s\n", matrix->id);
    fprintf(omfp, "AC   %s\n", matrix->ac);
    fprintf(omfp, "DE   %s\n", matrix->de);
    fprintf(omfp, "MA   %s\n", matrix->ma);
    
    fprintf(omfp, " ");
    for (c='A'; c <= 'Z'; c++) {
      if ((c != 'J') && (c != 'O') && (c != 'U')) {
	fprintf(omfp, " %c  ", c);
      }
    }
    fprintf(omfp, " *  ");
    fprintf(omfp, " -\n");
    
    for (l=0; l<matrix->width; l++) {
      for (c='A'; c <= 'Z'; c++) {
	if ((c != 'J') && (c != 'O') && (c != 'U')) {
	  fprintf(omfp, "% 6.4f ", matrix->weights[aa_atob[(int)c]][l]);
	}
      }
      fprintf(omfp, "% 6.4f ", matrix->weights[aa_atob['*']][l]);
      fprintf(omfp, "% 6.4f\n", matrix->weights[aa_atob['-']][l]);
    }  
    
    fprintf(omfp, "//\n");
  } 
  else { /* unknown */
    sprintf(ErrorBuffer, 
	    "Unknown output type: %d, using integer output\n", style);
    ErrorReport(WARNING_ERR_LVL);
    output_matrix_s(matrix, omfp, INT_OUTPUT);
  }
}  /*  end of output_matrix_s() */


/* 
 * output_matrix_st
 *   Outputs a matrix data structure to the given file with the 
 *   specified style of data. NA_SEQ => only ouput ACTG columns.
 *   Parameters: 
 *     Matrix *matrix: the matrix to output
 *     FILE   *omfp:   the output matrix file pointer
 *     int    style:   the kind of output (INT_OUTPUT, FLOAT_OUTPUT)
 *     int    type:    matrix type (AA_SEQ, NA_SEQ)
 *   Return code: none
 *   Error code: none
 */

void output_matrix_st(matrix, omfp, style, type)
     Matrix *matrix;
     FILE *omfp;
     int style, type;
{
  char c;
  int l;

  if ( (style != INT_OUTPUT) && (style != FLOAT_OUTPUT) )
  {
    sprintf(ErrorBuffer, 
	    "Unknown output style: %d, using integer output\n", style);
    ErrorReport(WARNING_ERR_LVL);
    style = INT_OUTPUT;
  }
  if ( (type != AA_SEQ) && (type != NA_SEQ) )
  {
    sprintf(ErrorBuffer, 
	    "Unknown sequence type: %d, using amino acid\n", type);
    ErrorReport(WARNING_ERR_LVL);
    type = AA_SEQ;
  }

  fprintf(omfp, "ID   %s\n", matrix->id);
  fprintf(omfp, "AC   %s\n", matrix->ac);
  fprintf(omfp, "DE   %s\n", matrix->de);
  fprintf(omfp, "MA   %s\n", matrix->ma);
  
  if (type == AA_SEQ)
  {
    for (c='A'; c <= 'Z'; c++) {
      if ((c != 'J') && (c != 'O') && (c != 'U'))
      {
         if (style == FLOAT_OUTPUT) fprintf(omfp, "    %c     ", c);
         else fprintf(omfp, "  %c ", c);
      }
    }
    if (style == FLOAT_OUTPUT) fprintf(omfp, "    *         -\n");
    else                       fprintf(omfp, "  *   -\n");
  }
  else
  {
     if (style == FLOAT_OUTPUT)
        fprintf(omfp, "    A         C         G         T\n");
     else
        fprintf(omfp, "  A   C   G   T\n");
  }

  if (style == INT_OUTPUT) {
    for (l=0; l<matrix->width; l++) {
      if (type == AA_SEQ)
      {
         for (c='A'; c <= 'Z'; c++) {
	   if ((c != 'J') && (c != 'O') && (c != 'U')) {
	     fprintf(omfp, "%3d ", (int) round(matrix->weights[aa_atob[(int)c]][l]));
	   }
         }
         fprintf(omfp, "%3d ", (int) round(matrix->weights[aa_atob['*']][l]));
         fprintf(omfp, "%3d\n", (int) round(matrix->weights[aa_atob['-']][l]));
      }
      else
      {
         for (c='A'; c <= 'T'; c++) {
	   if ((c == 'A') || (c == 'C') || (c == 'G') || (c == 'T')) {
	     fprintf(omfp, "%3d ", (int) round(matrix->weights[aa_atob[(int)c]][l]));
           }
         }
         fprintf(omfp, "\n");
      }  
    } /*  end of for l */
  }   /* end of INT_OUTPUT */

  else if (style == FLOAT_OUTPUT) {
    for (l=0; l<matrix->width; l++) {
      if (type == AA_SEQ)
      {
         for (c='A'; c <= 'Z'; c++) {
      	   if ((c != 'J') && (c != 'O') && (c != 'U')) {
	     fprintf(omfp, "% 9.4f ", matrix->weights[aa_atob[(int)c]][l]);
	   }
         }
         fprintf(omfp, "% 9.4f ", matrix->weights[aa_atob['*']][l]);
         fprintf(omfp, "% 9.4f\n", matrix->weights[aa_atob['-']][l]);
      }
      else
      {
         for (c='A'; c <= 'T'; c++) {
	   if ((c == 'A') || (c == 'C') || (c == 'G') || (c == 'T')) {
	     fprintf(omfp, "% 9.4f ", matrix->weights[aa_atob[(int)c]][l]);
           }
         }
         fprintf(omfp, "\n");
      }
    }  
    
  }  /*  end of FLOAT */ 
  fprintf(omfp, "//\n");
}  /*  end of output_matrix_st */

/* Change log information follows. */
/*
  Changes since version 3.3.2:
   6/ 7/00 Put width in matrix->ma in new_matrix()
  Changes since version 3.1:
   2/14/97 Added output_matrix_st()
  11/18/96 Changed new_matrix() to initialize fields.  JGH
  Changes since version 3.0.0:
  4/11/96  Changed read_matrix_body() to sscanf with %f instead of %d. JGH
  3/25/96  Changed MatType for matrix->weights from int to double. JGH
             Still assumes input matrix is integer & still exports integers.
*/

