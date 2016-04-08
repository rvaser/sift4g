/* (C) Copyright 1993, Fred Hutchinson Cancer Research Center */

/* Use, modification or distribution of these programs is subject to */
/* the terms of the non-commercial licensing agreement in license.h. */

/* blocks.c: basic block operations on the Blocks data type */
/* Written by: Bill Alford */
/* Change log information is at the end of the file. */

/*	system headers not in global.h */
/*	blimps library headers */
#include <global.h>
#include <blocks.h>
#include <residues.h>

/*
 * Exported variables and data structures
 */

/*
 * Local variables and data structures
 */

#define BLOCK_SEQUENCE_INCREASE_SIZE 20 /* default number of sequences to */
				   /* allocate if the number of sequences */
				   /* in the block is not specified */
#define BLOCK_CLUSTER_INCREASE_SIZE 2 /* default number of sequences to */
				   /* allocate if the number of sequences */
				   /* in the block is not specified */


#define DEFAULT_BLOCK_SEQUENCE_WEIGHT 0.0 /* the default weight of each */
					  /* sequence if there is no weight */
					  /* defined in the block */

/*
 * Function definitions
 */

/* 
 * Block file related functions
 *
 *   Block *read_a_block ()
 *   Boolean read_block_header(bfp, block)
 *   int read_block_body(bfp, block)
 *   int next_cluster(bfp, block, num_clusters_seen, num_sequences_seen)
 *   void resize_block_sequences(block)
 *   void resize_block_clusters(block)
 *   void free_block(block)
 *
 */

/*
static Boolean read_block_header();
static void read_block_body();
static void next_cluster();
static void resize_block_sequences();
static void resize_block_clusters();
*/


/*
 * read_a_block
 *   reads a block from the data base and returns a pointer to the new
 *   block data structure
 *   Parameters:  
 *     FILE *bfp: the file pointer the the blocks database/file
 *   Error codes: NULL if a block was not read
 */

static int BlockBufLength;
static char *BlockBuffer;

Block *read_a_block (bfp)
     FILE *bfp;			/* block file pointer */
{
  Block *new_block;


  /* get the block file pointer */
/*  bfp = get_file(BLOCK_FILES);*/
  if (bfp == NULL) {
    /* no more data to read into blocks */
    new_block = NULL;
    return new_block;
  }

  /* allocate space for a new block */
  CheckMem(
	   new_block = (Block *) malloc(sizeof(Block))
	   );

  /* set the max sizes of the clusters and sequences */
  new_block->max_clusters = 0;
  new_block->max_sequences = 0;
  new_block->undefined = 0;
  new_block->undefined_dbl = 0.0;
  new_block->undefined_ptr = NULL;

  /* read header */
  if (read_block_header(bfp, new_block) == FALSE) {
    free(new_block);
    return NULL;
  }

  /* setup to handle extra wide blocks */
  if (new_block->width + 30 < EXTRA_LARGE_BUFF) {
    BlockBufLength = EXTRA_LARGE_BUFF;
    BlockBuffer = Buffer;
  }
  else {
    BlockBufLength = new_block->width + 30;
    CheckMem(
	     BlockBuffer = malloc(sizeof(char)*BlockBufLength)
	     );
  }

  /* read body */
  read_block_body(bfp, new_block); 

  /* free the block buffer if it was malloc'd */
  if (BlockBuffer != Buffer) {
    free(BlockBuffer);
  }

  /* return block */
  return new_block;
}


/*
 * read_block_header
 *   Reads the header information for the block.
 ID   <id>; BLOCK
 AC   <ac>;
 DE   <de>
 BL   <block info>

 *   Parameters:
 *     FILE *bfp:    the block file pointer. 
 *     Block *block: the block to put the data in.
 *   Error Codes: FALSE if the block could not be read, TRUE otherwise
 */

Boolean read_block_header(bfp, block)
     FILE *bfp;
     Block *block;
{
  char *buf2, *buf3;
  Boolean bflag;			/* found a real block */
  Boolean idflag, acflag, deflag, blflag;
  int i, done;

  bflag = idflag = acflag = deflag = blflag = FALSE;

  block->bl[0] = '\0';	/* make sure the BL string is empty */
  block->motif[0] = '\0';	/* make sure the motif string is empty */
  block->width = 0; /* make sure the width is small */
  block->percentile = 0;	
  block->strength = 0;	
  block->max_sequences = 0;	
  block->num_sequences = 0;	
  block->max_clusters = 0;	
  block->num_clusters = 0;	
  block->min_prev = block->max_prev = 0;

  /* scan for first ID, AC, DE, BL */
  fgets(Buffer, LARGE_BUFF_LENGTH, bfp);

  /*  need to check further, can be fooled (see mablock)  */
  while (  (strncmp(Buffer, "ID   ", 5)) != 0 &&
           (strncmp(Buffer, "AC   ", 5)) != 0 &&
           (strncmp(Buffer, "DE   ", 5)) != 0 &&
           (strncmp(Buffer, "BL   ", 5)) != 0 
	 && !feof(bfp) ) { 
    fgets(Buffer, LARGE_BUFF_LENGTH, bfp);
  }

  /* if we reached the end of the file, then there was no block */
  if (feof(bfp)) { return FALSE; }

  /* read ID if it exists */
  if ( strncmp(Buffer, "ID   ", 5) == 0) 
  {
     if (Buffer[5] != ' ')
     {
       if (strstr(Buffer, "; MATRIX") == NULL)
       {
          idflag = bflag = TRUE;
          buf2 = &Buffer[2];		/* Eat the ID at the beginining */
          buf2 = eat_whitespace(buf2); /* a pointer to the string after ID */
          remove_trailing_whitespace(buf2); 
          strncpy(block->id, buf2, SMALL_BUFF_LENGTH);
       }
       else    /*  this looks like a pssm */
       {
          sprintf(ErrorBuffer, "MATRIX on ID line; not a BLOCK");
          ErrorReport(SERIOUS_ERR_LVL);
          return(bflag);
       }
     }
     else 
     { /* invalid ID field, raise error */
       strcpy(block->id, "none; BLOCK");
       sprintf(ErrorBuffer, "Error in block file format.  Invalid ID line:");
       ErrorReport(WARNING_ERR_LVL);
       sprintf(ErrorBuffer, "%s\n", Buffer);
       ErrorReport(WARNING_ERR_LVL);
     }

     fgets(Buffer, LARGE_BUFF_LENGTH, bfp);	/* read the next line */
     while (strncmp(Buffer, "CC   ", 5) == 0)	/* skip over comment lines */
     { fgets(Buffer, LARGE_BUFF_LENGTH, bfp); }
  }
  else 
  { /* no valid ID field, raise error */
    strcpy(block->id, "none; BLOCK");
    sprintf(ErrorBuffer, "Error in block file format.  No ID line.\n");
    ErrorReport(WARNING_ERR_LVL);
  }

  /* read AC if it exists */
  if ( strncmp(Buffer, "AC   ", 5) == 0)
  {
     if (Buffer[5] != ' ')
     {
       acflag = bflag = TRUE;
       buf2 = &Buffer[2];		/* Eat the AC at the beginining */
       buf2 = eat_whitespace(buf2); /* get a pointer to the string after AC */
       remove_trailing_whitespace(buf2); /* remove the \n */
       strncpy(block->ac, buf2, SMALL_BUFF_LENGTH);
       buf2 = get_token(buf2);	/* get the number */
       if (buf2[strlen(buf2)-1] == ';') 
       {
         buf2[strlen(buf2)-1] = '\0';      /* remove the ';'*/
       }
       strncpy(block->number, buf2, SMALL_BUFF_LENGTH);
       strcpy(block->family, block->number);
       i = strlen(block->family); 
       if (i > MAXAC) { block->family[MAXAC] = '\0'; i = MAXAC; }
       done = FALSE;
       while (!done && i >= MINAC)
       {
          if (isalpha(block->family[i])){ block->family[i] = '\0'; done=TRUE; }
          else  { i--;  }
       }
       /* ; distance from previous block=(min,max) */
       strcpy(Buffer, block->ac);		/* Buffer's been destroyed */
       buf3 = strstr(Buffer, "block=(");
       if (buf3!=NULL) 
       {
         sscanf(buf3, "block=(%d,%d);", &(block->min_prev), &(block->max_prev));
       }
     }
     else
     {
       strcpy(block->ac, "none; distance from previous block=( , )");
       block->number[0] = '\0';	/* make sure string is empty */
       block->family[0] = '\0';	/* make sure string is empty */
       sprintf(ErrorBuffer, "Error in block file format.  Invalid AC line:");
       ErrorReport(WARNING_ERR_LVL);
       sprintf(ErrorBuffer, "%s\n", Buffer);
       ErrorReport(WARNING_ERR_LVL);
     }
     fgets(Buffer, LARGE_BUFF_LENGTH, bfp);	/* read the next line */
     while (strncmp(Buffer, "CC   ", 5) == 0)	/* skip over comment lines */
     { fgets(Buffer, LARGE_BUFF_LENGTH, bfp); }
  }
  else { /* no AC field, raise error */
    strcpy(block->ac, "none; distance from previous block=( , )");
    block->number[0] = '\0';	/* make sure string is empty */
    sprintf(ErrorBuffer, "Error in block file format.  No AC line.\n");
    ErrorReport(WARNING_ERR_LVL);
  }

  /* read DE if it exists */
  /* Be more strict about accepting DE to start a block if haven't seen
	ID or AC line yet */
  if (strncmp(Buffer, "DE   ", 5) == 0)
  {
     if ( Buffer[5] != ' ' || idflag || acflag ) 
     {
       deflag = bflag = TRUE;
       buf2 = &Buffer[2];		/* Eat the DE at the beginining */
       buf2 = eat_whitespace(buf2); /* get a pointer to the string after DE */
       remove_trailing_whitespace(buf2); /* remove the \n  */
       strncpy(block->de, buf2, SMALL_BUFF_LENGTH);
     }
     else 
     {
       strcpy(block->de, "none");
       sprintf(ErrorBuffer, "Error in block file format.  Invalid DE line:");
       ErrorReport(WARNING_ERR_LVL);
       sprintf(ErrorBuffer, "%s\n", Buffer);
       ErrorReport(WARNING_ERR_LVL);
     }
     fgets(Buffer, LARGE_BUFF_LENGTH, bfp);	/* read the next line */
     while (strncmp(Buffer, "CC   ", 5) == 0)	/* skip over comment lines */
     { fgets(Buffer, LARGE_BUFF_LENGTH, bfp); }
  }
  else
  {
    strcpy(block->de, "none");
    sprintf(ErrorBuffer, "Error in block file format.  No DE line.\n");
    ErrorReport(WARNING_ERR_LVL);
  }
  
  /* read BL if it does not exist raise an error */
  /* Blocks from MEME have only BL line, no ID|AC|DE   */
  /* Be more strict about accepting BL to start a block if haven't seen
	ID or AC line yet */
  if (strncmp(Buffer, "BL   ", 5) == 0)
  {
     if ( Buffer[5] != ' ' || idflag || acflag )
     {
       blflag = bflag = TRUE;
       buf2 = &Buffer[2];		/* Eat the BL at the beginining */
       buf2 = eat_whitespace(buf2); /* get a pointer to the string after BL */
       remove_trailing_whitespace(buf2); /* remove the \n  */
       strncpy(block->bl, buf2, SMALL_BUFF_LENGTH);
    
       /* scan and process the BL line */
       /* find the motif */
       buf3 = strstr(Buffer, "motif");
       if (buf3!=NULL) 
       {
         buf3 -= 4;		/* move back to the motif type */
         sscanf(buf3, "%s", block->motif);
       }

       /* find the width */
       buf3 = strstr(Buffer, "width=");
       if (buf3!=NULL) 
       { sscanf(buf3, "width=%d;", &(block->width)); }
       else 
       {
         /* report an error, it is serious if the width is not set */
         sprintf(ErrorBuffer, "No width field for block %s", block->number);
         ErrorReport(SERIOUS_ERR_LVL);
         sprintf(ErrorBuffer, "Setting width to zero\n");
         ErrorReport(SERIOUS_ERR_LVL);
         block->width = 0;
       }

       /* find the seqs */
       buf3 = strstr(Buffer, "seqs=");
       if (buf3!=NULL) 
       {
         sscanf(buf3, "seqs=%d", &(block->num_sequences));
         if (block->num_sequences < 0) { block->num_sequences = 0; }
       }
    
       /* find the 99.5% */
       buf3 = strstr(Buffer, "99.5%=");
       if (buf3!=NULL) { sscanf(buf3, "99.5%%=%d;", &(block->percentile)); }

       /* find the strength */
       buf3 = strstr(Buffer, "strength=");
       if (buf3!=NULL) { sscanf(buf3, "strength=%d", &(block->strength)); }

     }  /* end of valid BL */
     else 
     { 
       sprintf(ErrorBuffer, "Error in block file format.  Invalid BL line:");
       ErrorReport(SERIOUS_ERR_LVL);
       sprintf(ErrorBuffer, "%s\n", Buffer);
       ErrorReport(SERIOUS_ERR_LVL);
     }
   }    /* end of "BL   " */
   else
   {
       sprintf(ErrorBuffer, "Error in block file format.  No BL line.\n");
       ErrorReport(SERIOUS_ERR_LVL);
   }

  return(bflag);
}  /* end of read_block_header  */


/*
 * read_block_body
 *   Reads the body information of the block.
 *   Parameters:
 *     FILE *bfp:    the block file pointer. 
 *     Block *block: the block to put the data in.
 *   Error Codes:
 *   Programming note: I'm trying to malloc as large of chunks as
 *     possible hoping that there will be a slight speed increase in
 *     the program due to locality of reference (less page swapping).
 */

void read_block_body(bfp, block)
     FILE *bfp;
     Block *block;
{
  int i;
  
  Sequence *sequence_pointer;
  Residue  *residue_pointer;
  
  if (block->num_sequences <= 0) {
    block->max_sequences = BLOCK_SEQUENCE_INCREASE_SIZE;
  }
  else {
    block->max_sequences = block->num_sequences;
  }

  /* allocate space for all the Sequences of the block sequence array */
  CheckMem(
	   sequence_pointer = (Sequence *) calloc(block->max_sequences,
						  sizeof(Sequence))
	   );

  /* initialize the block sequences array */
  block->sequences = sequence_pointer;
  /* allocate space for all residues */
  CheckMem(
	   residue_pointer = (Residue *) calloc(block->width *
						block->max_sequences,
						sizeof(Residue))
	   );

  /* initialize the residues 2-d array */
  CheckMem(
	   block->residues = 
	   (Residue **) calloc(block->max_sequences,
			       sizeof(Residue *))
	   );

/* testing: */ /*
  for(i=0; i<(block->width *block->max_sequences); i++) {
    residue_pointer[i] = i;
  }
*/

  /* initialize sequences and residues */
  for(i=0; i<block->max_sequences; i++) {
    sequence_pointer[i].length = block->width;
    sequence_pointer[i].sequence = 
      &(residue_pointer[i * block->width]);
    block->residues[i] = &(residue_pointer[i * block->width]);
  }


/* testing: */ /*
  for(i=0;i<block->max_sequences;i++) {
    printf("%d\t",block->sequences[i].sequence[0]);
  }
  printf("\n");
*/

  /* read the clusters */
  /* num clusters seen = 0, num seqs seen = 0 */
  next_cluster(bfp, block, 0, 0); 
}  /* end of read_block_body */


/*
 * next_cluster
 *   next_cluster reads in the next set of sequences for a cluster.
 *   When the end of block, "//", is seen the space for the clusters
 *   is allocated.  This procedure calls itself recursively so when
 *   the memory is allocated all of the clusters have been seen and
 *   the data is stored within each clusters procedure.
 *   Parameters:
 *     FILE *bfp: the block file pointer
 *     Block *block: the block of the clusters.  Where the sequences are.
 *     int num_clusters_seen: the number of clusters seen so far
 *     int num_sequences_seen: the number of sequences seen so far
 *   Precondition: There is no space allocated for the clusters.  Attempts
 *                 to access will result in ambiguous behavior
 *   Postcondition: The space for has been allocated.  The clusters after
 *                  the current number seen are correctly set up.
 *   Error Codes:
 */

void next_cluster(bfp, block, num_clusters_seen, num_sequences_seen)
     FILE *bfp;
     Block *block;
     int num_clusters_seen;
     int num_sequences_seen;
{
  int i, eob;
  int low_sequence;
  int return_value;

  char *buf, *name_buf;

  low_sequence = num_sequences_seen;	/* this clusters starting Sequence */
  eob = FALSE;		/*  end of block flag */
  fgets(BlockBuffer, BlockBufLength, bfp); /* read a line */
  if (feof(bfp)) eob = TRUE;
  if (BlockBuffer[0] == '/' && BlockBuffer[1] == '/') eob = TRUE;

  /* while not a blank line and not the end of the block read the sequences */
  /* also be sure not at next block (ID   ) or end of file */
  while (!blank_line(BlockBuffer) && !eob)
  {
    /* get the name, getting early for use in error reporting if needed */
    name_buf = get_token(BlockBuffer);
  
    /* it is a sequence, so make sure there is enough space for the sequence */
    if (block->max_sequences < num_sequences_seen +1) { /* indexing starts */
  	    /* at 1 for block->max_sequences and 0 for num_sequences_seen  */
      /* more sequences than space */

      /* announce error */
      sprintf(ErrorBuffer, 
	      "No space allocated for sequence %s in block %s.",
	      name_buf, block->number);
      ErrorReport(INFO_ERR_LVL);
      sprintf(ErrorBuffer,
	      "Allocating more space\n");
      ErrorReport(INFO_ERR_LVL);
      
      /* get more space */
      resize_block_sequences(block);
    }

    /* enter the sequence name into the data struct */
    strncpy(block->sequences[num_sequences_seen].name, name_buf, SMALL_BUFF_LENGTH);

    buf = get_token(NULL);		/* eat the "(" or "(#)" */
    if (buf != NULL && buf[1] != '\0') {
      /* this is mainly the old style format */
      sscanf(buf, "(%d)", &(block->sequences[num_sequences_seen].position));
    }
    else {     
      buf = get_token(NULL);		/* get the part with the position */
      if (buf != NULL)
      {  sscanf(buf, "%d)", &(block->sequences[num_sequences_seen].position));}
    }
    
    buf = get_token(NULL);		/* skip over position to the seq. */
    

    /* read the sequence from the string into the Sequence */

    return_value = read_sequence(&(block->sequences[num_sequences_seen]), 
				 AA_SEQ, 0, buf); 
				/* the Sequence, sequence type, starting */
				/* position, string with the sequence */
    /* check to see if there was an error */
    
    if (return_value < 0) {
      /* Error, more residues in the sequence than expected */
      sprintf(ErrorBuffer, "Error reading sequence %s in block %s,", 
	      block->sequences[num_sequences_seen].name,
	      block->number);
      ErrorReport(WARNING_ERR_LVL);
      sprintf(ErrorBuffer, 
	      "%d more residues in the sequence than the expected %d.",
	      - return_value,
	      block->sequences[num_sequences_seen].length);
      ErrorReport(WARNING_ERR_LVL);
      /*  reset the lengths & reread the sequence */
      block->width -= return_value;
      block->sequences[num_sequences_seen].length -= return_value;
      resize_block_sequences(block);
      sprintf(ErrorBuffer, 
	      "Resetting block width to %d.\n",
	      block->sequences[num_sequences_seen].length);
      ErrorReport(WARNING_ERR_LVL);
      return_value = read_sequence(&(block->sequences[num_sequences_seen]), 
				 AA_SEQ, 0, buf); 
    }
    else if (return_value < block->sequences[num_sequences_seen].length) {
      /* Error, not enough residues for the sequence, filling with blanks */
      sprintf(ErrorBuffer, "Error reading sequence %s in block %s,", 
	      block->sequences[num_sequences_seen].name,
	      block->number);
      ErrorReport(WARNING_ERR_LVL);
      sprintf(ErrorBuffer, 
	      "not enough residues to fill the sequence.");
      ErrorReport(WARNING_ERR_LVL);
      sprintf(ErrorBuffer, 
	      "Filling the rest of the sequence with blanks.\n");
      ErrorReport(WARNING_ERR_LVL);

      /* filling with blanks */
      if (block->sequences[num_sequences_seen].type == AA_SEQ) {
	for(i=return_value;i<block->sequences[num_sequences_seen].length;i++) {
	  block->sequences[num_sequences_seen].sequence[i] = aa_atob['-'];
	}
      }
      else if (block->sequences[num_sequences_seen].type == NA_SEQ) {
	for(i=return_value;i<block->sequences[num_sequences_seen].length;i++) {
	  block->sequences[num_sequences_seen].sequence[i] = nt_atob['-'];
	}
      }
      else {
	sprintf(ErrorBuffer, 
		"next_cluster(): Unknown sequence type, assuming an amino acid sequence.\n");
	ErrorReport(PROGRAM_ERR_LVL);
	for(i=return_value;i<block->sequences[num_sequences_seen].length;i++) {
	  block->sequences[num_sequences_seen].sequence[i] = aa_atob['-'];
	}
      }
    }
    else if (return_value > block->sequences[num_sequences_seen].length) {
      /* BIG Error, an undocumented return value from read_sequence */
      sprintf(ErrorBuffer, 
	      "next_cluster(): Error reading sequence %s in block %s,", 
	      block->sequences[num_sequences_seen].name,
	      block->number);
      ErrorReport(PROGRAM_ERR_LVL);
      sprintf(ErrorBuffer, 
	      "                Undocumented return value, %d, from read_sequence().\n",
	      return_value);
      ErrorReport(PROGRAM_ERR_LVL);
    } /* else, ret_value == seq.length, OK */

    /* if there is a sequence weight, read it */
    buf = get_token(NULL);
    if (buf != NULL) {
      /* there is a weight */
      sscanf(buf, "%lg", &block->sequences[num_sequences_seen].weight);
    }
    else {
      /* there was no weight */
      block->sequences[num_sequences_seen].weight = 
	DEFAULT_BLOCK_SEQUENCE_WEIGHT;
    }
    
    /* update sequence count */
    num_sequences_seen++;

    /* get the next line */
    fgets(BlockBuffer, BlockBufLength, bfp);	
    if (feof(bfp)) eob = TRUE;
    if (BlockBuffer[0] == '/' && BlockBuffer[1] == '/') eob = TRUE;
    /*  Assumes no extra sequence will be named "ID     " !!!  */
    if (num_sequences_seen >= block->num_sequences &&
        strncmp(BlockBuffer, "ID   ", 5) != 0)
    {
       eob = TRUE;
       /*  should also backup bfp one record; see ftell()/fseek()  */
    }
  } /* end of while */

  /* if at the end of block: */
  if (eob) { 

    block->num_clusters = num_clusters_seen + 1; /* reached end saw this one */
    block->max_clusters = block->num_clusters;  

    /* allocate space for all the clusters */
    CheckMem(
	     block->clusters = (Cluster *) calloc(block->max_clusters,
						  sizeof(Cluster))
	     );


    /* report size difference */
    if (block->num_sequences != num_sequences_seen) {
      sprintf(ErrorBuffer, 
	      "Number of sequences in block %s read from file as %d but is %d\n",
	      block->number,
	      block->num_sequences,
	      num_sequences_seen);
      ErrorReport(INFO_ERR_LVL);
      block->num_sequences = num_sequences_seen;
    }

  }
  else {
    /* read in the next cluster */
    next_cluster(bfp, block, num_clusters_seen+1, num_sequences_seen);
  }

  /* assign the sequences to the cluster, space has now been allocated */
  block->clusters[num_clusters_seen].num_sequences = 
    num_sequences_seen-low_sequence;
  
  block->clusters[num_clusters_seen].sequences = 
    &(block->sequences[low_sequence]);
}  /*  end of next_cluster */


/* 
 * resize_block_sequences
 *   Increases the memory for the storage of the sequences of a block.
 *   Parameter:
 *     Block *block: the block to resize
 *   Error codes: none
 */

void resize_block_sequences(Block *block)
{
  int i, j, s;
  void *tmp_ptr;		/* don't want to ruin pointer in case there */
				/* is a repeated realloc call */

  /* announce info */
  sprintf(ErrorBuffer, 
	  "Allocating more space for sequences in block %s.\n",
	  block->number);
  ErrorReport(INFO_ERR_LVL);

  
  /* get more space */
  block->max_sequences += BLOCK_SEQUENCE_INCREASE_SIZE;
  
  CheckMem(
	   tmp_ptr = 
	   (Sequence *) realloc(block->sequences,
				sizeof(Sequence) *
				block->max_sequences)
	   );
  block->sequences = (Sequence *) tmp_ptr;
  
  CheckMem(
	   tmp_ptr = 
	   (Residue *) realloc(block->sequences[0].sequence,
			       sizeof(Residue) *
			       block->width *
			       block->max_sequences)
	   );
  block->sequences[0].sequence = (Residue *) tmp_ptr;

  CheckMem(
	   tmp_ptr = 
	   (Residue **) realloc(block->residues,
				sizeof(Residue *) *
				block->max_sequences)
	   );
  block->residues = (Residue **) tmp_ptr;
  
  /* re-initialize sequences and residues, same data, just different place */
  for(i=0; i<block->max_sequences; i++) {
    block->sequences[i].length = block->width;
    block->sequences[i].sequence = 
      &(block->sequences[0].sequence[i * block->width]);
    block->residues[i] = 
      &(block->sequences[0].sequence[i * block->width]);
  }
  /*  Need to re-initialize clusters; point to first seq,
      assumes sequences are in cluster order  */
  s = 0;	/* sequence number */
  for (j=0; j<block->num_clusters; j++)
  {
      block->clusters[j].sequences = &(block->sequences[s]);
      for (i=0; i<block->clusters[j].num_sequences; i++)
              s++;
  }

}


/* 
 * resize_block_clusters
 *   Increases the memory for the storage of the clusters of a block.
 *   Parameter:
 *     Block *block: the block to resize
 *   Error codes: none
 */

void resize_block_clusters(block)
     Block *block;
{
  Cluster *tmp_ptr;		/* don't want to ruin pointer in case there */
				/* is a repeated realloc call */

  /* announce info */
  sprintf(ErrorBuffer, 
	  "Allocating more space for clusters in block %s.\n",
	  block->number);
  ErrorReport(INFO_ERR_LVL);

  
  block->max_clusters += BLOCK_CLUSTER_INCREASE_SIZE;

  CheckMem(
	   tmp_ptr = (Cluster *) realloc(block->clusters,
					 block->max_clusters *
					 sizeof(Cluster))
	   );
  block->clusters = tmp_ptr;
}


/* 
 * free_block
 *   Frees the block and the sub elements.
 *   Parameters: 
 *     Score *score: the score to free
 *   Return code: none
 *   Error code: none
 */

void free_block(block) 
     Block *block;
{
  free(block->residues);
  free(block->clusters);
  free(block->sequences[0].sequence);
  free(block->sequences);
  free(block);
}


/*
 * block_comparison
 *   Compares two blocks.   It compares by the value in the id
 *   field if it exists.
 *   Parameters:
 *     BlockListEntry a, b: the entries to compare
 *   Return codes:  a return value < 0 if a < b, a return value = 0 if a == b,
 *                  and a return value > 0 if a > b
 *   Error codes: none
 */

int block_comparison(a, b)
     Block *a, *b;
{
  char *sa, *sb;

  if (a->id != NULL) {    
    sa = a->id;
  }
  else {
    if (b->id != NULL) {  
      return -1;		/* NULL < something */
    }
    else {
      return 0;			/* NULL = NULL */
    }
  }
  
  if (b->id != NULL) {    
    sb = b->id;
  }
  else {
    return 1;			/* something > NULL */
  }

  return (int) (strcmp(sa, sb));
}



/*
 * Block printing.
 */

/*
 * print_block
 *   Prints a block data structure.  Primarily for debugging purposes.
 *   Parameters: 
 *     Block *block:  the block to print
 *   Error Codes: none
 */

void print_block(block)
     Block *block;
{
  int i,j,k;

  printf("--- block ---\n");
  printf("ID\t%s\n", block->id);
  printf("AC\t%s\n", block->ac);
  printf("DE\t%s\n", block->de);
  printf("BL\t%s\n", block->bl);
  printf("motif:\t%s\n", block->motif);
  printf("width=%d; 99.5%%=%d; strength=%d;\n", 
	 block->width,
	 block->percentile,
	 block->strength);
  printf("num sequences: %d\tnum clusters: %d\n",
	 block->num_sequences,
	 block->num_clusters);

  /* print clusters */
  for (i=0; i<block->num_clusters; i++) {
    for (j=0; j<block->clusters[i].num_sequences; j++) {
/*      print_sequence(&(block->clusters[i].sequences[j])); */
      printf("%s \t( %d:%d, ", 
	     block->clusters[i].sequences[j].name,
	     block->clusters[i].sequences[j].position,
	     block->clusters[i].sequences[j].length);
      if (block->clusters[i].sequences[j].type == AA_SEQ) {
	printf("AA) \t");
	for(k=0; k<block->clusters[i].sequences[j].length; k++) {
	  printf("%c", aa_btoa[block->clusters[i].sequences[j].sequence[k]]);
	}
      }
      else if (block->clusters[i].sequences[j].type == NA_SEQ) {
	printf("NA) \t");
	for(k=0; k<block->clusters[i].sequences[j].length; k++) {
	  printf("%c", nt_btoa[block->clusters[i].sequences[j].sequence[k]]);
	}
      }    
      else {
	printf("\?\?) \t");
	for(k=0; k<block->clusters[i].sequences[j].length; k++) {
	  printf("%c", aa_btoa[block->clusters[i].sequences[j].sequence[k]]);
	}
      }
      printf("\n");
    } /* end of cluster */
    printf("\n");
  }  

  /* print sequences */
  printf("--- ----- ---\n");
  
}


/*
 * ouput_block
 *   Outputs a block data structure to the given file. 
 *   Parameters: 
 *     Block *block:  the block to print
 *     FILE  *obfp:   the ouput block file pointer
 *   Error Codes: none
 */

void output_block(block, obfp)
     Block *block;
     FILE *obfp;
{
  output_block_s(block, obfp, INT_OUTPUT);
}


/*
 * ouput_block_s
 *   Outputs a block data structure to the given file with the 
 *   specified style of data.
 *   Parameters: 
 *     Block *block:  the block to print
 *     FILE  *obfp:   the ouput block file pointer
 *     int    style:   the kind of output (INT_OUTPUT, FLOAT_OUTPUT)
 *   Error Codes: none
 */

void output_block_s(block, obfp, style)
     Block *block;
     FILE *obfp;
     int style;
{
  int i,j,k, offset, lenb, iend, ic;
  char *ptr, bltemp[132];

  if (!((style == INT_OUTPUT) || (style == FLOAT_OUTPUT))) {
    sprintf(ErrorBuffer, 
	    "Unknown output type: %d, using integer output\n", style);
    ErrorReport(WARNING_ERR_LVL);
    output_block_s(block, obfp, INT_OUTPUT);
    return;
  }

  fprintf(obfp, "ID   %s\n", block->id);
  fprintf(obfp, "AC   %s\n", block->ac);
  fprintf(obfp, "DE   %s\n", block->de);

  /*  Edit BL line: update or add width, seqs, 99.5, strength
      BL    ECA motif; width=40; seqs=34; 99.5%=1833; strength=1412 */
/*  This may wipe out some other stuff on the BL line
  sprintf(block->bl, "%s motif; width=%d; seqs=%d; 99.5%%=%d; strength=%d ",
      block->motif, block->width, block->num_sequences,
      block->percentile, block->strength);
*/
  /*  Rewrite everything after width= if it is found, else append */
  lenb = strlen(block->bl);
  ptr = strstr(block->bl, "width=");
  if (ptr!=NULL) 
  {
     strcpy(bltemp, block->bl);
     block->bl[0] = '\0';
     offset = lenb - strlen(ptr);
     strncat(block->bl, bltemp, offset);	/* preserve first part */
     block->bl[offset] = '\0';
  }
  sprintf(bltemp, " width=%d; seqs=%d; 99.5%%=%d; strength=%d ",
      block->width, block->num_sequences,
      block->percentile, block->strength);
  strcat(block->bl, bltemp);
  fprintf(obfp, "BL   %s", block->bl);

  /* print clusters */
  for (i=0; i<block->num_clusters; i++) {
    fprintf(obfp, "\n");
    for (j=0; j<block->clusters[i].num_sequences; j++) {
      strcpy(bltemp, block->clusters[i].sequences[j].name);
      if (strlen(bltemp) > 20)   /* look for nearest | */
      {
         iend = ic = 20;
         while (ic > 0)
         {
            if (bltemp[ic] == '|') { iend = ic; break; }
            ic--;
         }
         bltemp[iend] = '\0';
      }
      /*  right-justify is %20.20s */
      fprintf(obfp, "%-20s (%4d) ", 
              bltemp,
	      block->clusters[i].sequences[j].position);
      for(k=0; k<block->clusters[i].sequences[j].length; k++) {
	fprintf(obfp, "%c", 
		aa_btoa[block->clusters[i].sequences[j].sequence[k]]);
      }
      if (style == INT_OUTPUT) {
	fprintf(obfp, " %3d\n", 
		(int) round(block->clusters[i].sequences[j].weight));
      }
      else { /* FLOAT_OUTPUT -- only two options allowed with beginning if */
	fprintf(obfp, " %f\n", 
		(double) block->clusters[i].sequences[j].weight);
      }
    } /* end of cluster */
  }  
  /* end of the block */
  fprintf(obfp, "//\n");
  fflush(obfp);
}




/*
 *  new_block(#columns, #rows)
 *     Create a block structure
 */


Block *new_block(ncols, nrows)
int ncols, nrows;
{
   Block *new;
   Residue *residue_pointer;
   int i;

   /* allocate space for a new block */
   CheckMem(
	   new = (Block *) malloc(sizeof(Block))
	   );

   /*  Initialize  */
   new->max_sequences = new->num_sequences = nrows;
   new->width = ncols;
   new->id[0] = new->ac[0] = new->de[0] = '\0'; 
   strcpy(new->id, "id");
   strcpy(new->ac, "ac");
   strcpy(new->de, "de");
   sprintf(new->bl, " ; width=%d; seqs=%d; ", ncols, nrows);
   new->number[0] = new->family[0] = new->motif[0] = '\0';
   new->percentile = new->strength = new->undefined = 0;
   new->undefined_dbl = 0.0;
   new->undefined_ptr = NULL;

   /*  Allocate space for a single "cluster"  */
   new->max_clusters = new->num_clusters = 1;
   CheckMem(
	     new->clusters = (Cluster *) calloc(new->max_clusters,
						  sizeof(Cluster))
	     );


   new->clusters[0].num_sequences = nrows;

  /* allocate space for all the Sequences of the block sequence array */
  CheckMem(
	   new->sequences = (Sequence *) calloc(new->max_sequences,
						  sizeof(Sequence))
	   );

  new->clusters[0].sequences = &(new->sequences[0]);

  /* allocate space for all residues */
  CheckMem(
	   residue_pointer = (Residue *) calloc(new->width *
						new->max_sequences,
						sizeof(Residue))
	   );

  /* initialize the residues 2-d array */
  CheckMem(
	   new->residues = 
	   (Residue **) calloc(new->max_sequences, sizeof(Residue *))
	   );

  /* initialize sequences and residues */
  for(i=0; i<new->max_sequences; i++) {
    new->sequences[i].weight = 1.0;
    new->sequences[i].length = new->width;
    new->sequences[i].sequence = &(residue_pointer[i * new->width]);
    new->residues[i] = new->sequences[i].sequence;
  }

   /* return block */
   return new;
}  /*  end of new_block */

/*
 * read_to_block
 *   positions file at AC
 *   assumes AC is 7 chars; assumes bfp is sorted by AC
 *   Parameters:  
 *     FILE *bfp: the file pointer the the blocks database/file
 *     char *ac:  AC of the block to fine
 *   
 */

Boolean read_to_block(bfp, ac)
     FILE *bfp;			/* block file pointer */
     char *ac;
{
  long idpos;

  idpos = -1;
  while (!(feof(bfp)) && fgets(Buffer, LARGE_BUFF_LENGTH, bfp) != NULL)
  {
     if (strncmp(Buffer, "ID   ", 5) == 0) { idpos = ftell(bfp); }
     else
     {
        if (strncmp(Buffer, "AC   ", 5) == 0) 
        {
           if (strncmp(Buffer+5, ac, 7) > 0)  /* beyond AC */
           {   return(FALSE);  }
           else if (strncmp(Buffer+5, ac, 7) == 0)
           {
              fseek(bfp, idpos, 0);	/* back up to start of block */
              return(TRUE);
           }
        }
     }
  }
  return(FALSE);
}  /*  end of read_to_block  */

/*
 * read_a_block_faster
 *   reads a block from the data base and returns a pointer to the new
 *   block data structure
 *   Parameters:  
 *     FILE *bfp: the file pointer the the blocks database/file
 *   Error codes: NULL if a block was not read
>>>>>still only about 20% faster <<<<<
 */

Block *read_a_block_faster(bfp)
     FILE *bfp;			/* block file pointer */
{
  Block *new_block;
  Residue  *residue_pointer;
  char ctemp[80], ctemp2[20];
  int iseq, pos, eob;

  if (bfp == NULL) { return(NULL); }

  /* allocate space for a new block */
  CheckMem( new_block = (Block *) malloc(sizeof(Block)) );
  new_block->undefined = 0;
  new_block->undefined_dbl = 0.0;
  new_block->undefined_ptr = NULL;

  /* read header */
  if (read_block_header(bfp, new_block) == FALSE) {
    free(new_block);
    return(NULL);
  }

  /* This routine assumes num_sequences is set & correct  */
  if (new_block->num_sequences <= 0)
  {
    read_block_body(bfp, new_block);
    return(new_block);
  }
  
  new_block->max_sequences = new_block->num_sequences;
  /* allocate space for all the Sequences of the block sequence array */
  CheckMem(
	   new_block->sequences = (Sequence *) calloc(new_block->max_sequences,
						  sizeof(Sequence))
	   );

  /* allocate space for all residues */
  CheckMem(
	   residue_pointer = (Residue *) calloc(new_block->width *
						new_block->max_sequences,
						sizeof(Residue))
	   );

   /* initialize the residues 2-d array */
   CheckMem(
	   new_block->residues = 
	   (Residue **) calloc(new_block->max_sequences,
			       sizeof(Residue *))
	   );

   /* initialize sequences and residues */
   for(iseq=0; iseq<new_block->max_sequences; iseq++) {
      new_block->sequences[iseq].length = new_block->width;
      new_block->sequences[iseq].max_length = new_block->width;
      new_block->sequences[iseq].sequence = 
        &(residue_pointer[iseq * new_block->width]);
      new_block->residues[iseq] = &(residue_pointer[iseq * new_block->width]);
    }

    /* allocate space for one cluster */
    new_block->num_clusters = new_block->max_clusters = 1;
    CheckMem(
	     new_block->clusters = (Cluster *) calloc(new_block->max_clusters,
						  sizeof(Cluster))
	     );
   new_block->clusters[0].num_sequences = new_block->num_sequences;
   new_block->clusters[0].sequences = &(new_block->sequences[0]);

  /*  Now get the sequences  */
  eob = FALSE;		/*  end of block flag */
  iseq = 0;
  while (!eob && !(feof(bfp)) && fgets(Buffer, LARGE_BUFF_LENGTH, bfp) != NULL)
  {
     if ( (strncmp(Buffer, "//", 2) == 0)    ||
          (strncmp(Buffer, "ID   ", 5) == 0)    )
                eob = TRUE;
     /*  see if this looks like a sequence line:
         <sp>name<sp>(offset)<sp>sequence<sp>weight  */
     if (iseq < new_block->max_sequences &&
         strlen(Buffer) > 5 && strstr(Buffer, "(") != NULL)
     {
        sscanf(Buffer, "%s (%d) %s %s",
           (char *)&(new_block->sequences[iseq].name),
           &(new_block->sequences[iseq].position),
           (char *)&ctemp,
           (char *)&ctemp2  );

        new_block->sequences[iseq].weight = atof(ctemp2);
        if ( (int) strlen(ctemp) != new_block->width)
        {
          sprintf(ErrorBuffer, "Error in block %s: seq %s\n",
           new_block->number, new_block->sequences[iseq].name);
          ErrorReport(WARNING_ERR_LVL);
        }
        for (pos=0; pos < new_block->width; pos++)
           new_block->sequences[iseq].sequence[pos] = aa_atob[(int)ctemp[pos] ];
        strcpy(new_block->sequences[iseq].info, new_block->sequences[iseq].name);
        iseq++;
     }
  }
  return(new_block);
}  /* end of read_a_block_faster */



/* Change log information follows. */
/*
  Changes since version 3.8:
12/23/06 output_block() increase sequences->name from 18 to 20
  Changes since version 3.6:
 4/15/04 Fix abort in next_cluster() if reading incomplete blocks.
  Changes since version 3.4:
12/23/00 BlockBufLength starts at EXTRA_LARGE_BUFF
  Changes since version 3.3.2:
 2/23/00  output_block() Be sure seq name is < 18 chars
12/17/99  undefined_dbl and family fields added.
  Changes since version 3.2.5:
 5/11/99  Changed read_block_header() to do a better job of recognizing
	  whether the input contains a block.
 2/22/99  Added read_a_block_faster() & read_to_block()
 2/11/99  read_block_header(): get block->min_prev & block->max_prev
 1/28/99  output_block(): sequence name field increased from 10 to 18
  Changes since version 3.2.2:
 1/29/98  Removed all static declarations.
          Re-set cluster pointers in resize_block_sequences()
12/30/97  Changed read_block_header() to do a better job of recognizing
	  whether the input contains a block.
  Changes since version 3.2.1:
 7/10/97  Changed next_cluster() to recover from missing width on BL line.
 7/ 7/97  Modified output_block_s() to rebuild BL line.
  Changes since version 3.2:
 4/16/97  Removed static declaration from resize_block_sequences()
  Changes since version 3.1:
 2/14/97  Added new_block().
          Allow CC comment lines between header lines before BL line.
  Changes since version 3.0.0:
 4/15/96  Fix bug in read_a_block()/next_cluster() if block doesn't end with //.
*/
