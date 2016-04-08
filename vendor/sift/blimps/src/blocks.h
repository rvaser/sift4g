/* (C) Copyright 1993-7, Fred Hutchinson Cancer Research Center */
/* Use, modification or distribution of these programs is subject to */
/* the terms of the non-commercial licensing agreement in license.h. */

/* blocks.h: declaration of blocks data type and basic block operations */
/* Written by: Bill Alford */
/* Change log information is at the end of the file. */

#ifndef BLOCKS_H_
#define BLOCKS_H_

#define MINAC 7			/* min & max possible AC lengths */
#define MAXAC 10

#include "output.h"
#include "sequences.h"

/*
 * Exported variables and data structures
 */


/*
 * Sequence, Cluster, and Block data types
 */


struct cluster_struct {
  int      num_sequences;	   /* number of sequences in the block */
  Sequence *sequences;		   /* the array of Sequences of the cluster */
};
typedef struct cluster_struct Cluster;

struct block_struct {
  char     id[SMALL_BUFF_LENGTH];  /* Block ID string */
  char     ac[SMALL_BUFF_LENGTH];  /* Block AC string */
  char     de[SMALL_BUFF_LENGTH];  /* Block DE string */
  char     bl[SMALL_BUFF_LENGTH];  /* Block BL string */
  char     number[SMALL_BUFF_LENGTH]; /* AC:  block number IPR123456A  */
  char     family[SMALL_BUFF_LENGTH]; /* AC: block family IPR123456 */
  char     motif[20];		   /* BL: motif */
  int      width;		   /* BL: block width */
  int      percentile;		   /* BL: 99.5% score */
  int      strength;		   /* BL: strength (median TP score) */
  int      max_sequences;          /* max number of sequences in the block */
  int      num_sequences;          /* number of sequences in the block */
  int      max_clusters;	   /* max number of clusters in the block */
  int      num_clusters;	   /* number of clusters in the block */
  int      min_prev;               /* AC: min distance from previous block */
  int      max_prev;	           /* AC: max distance from previous block */
  int undefined;		/* for future use */
  double undefined_dbl;		/* for future use */
  void *undefined_ptr;		/* for future use */
  Cluster  *clusters;		   /* the array of Clusters of the block */
  Sequence *sequences;		   /* the array of Sequences of the block */
  Residue  **residues;		   /* the 2-d array of residues [seq][pos] */
};
typedef struct block_struct Block;


/*
 * Exported functions
 */

/*
 * read_a_block
 *   reads a block from the data base and returns a pointer to the new
 *   block data structure
 *   Parameters:
 *     FILE *bfp: the file pointer the the blocks database/file
 *   Error codes: NULL if a block was not read
 */

extern Block *read_a_block ();
extern Block *read_a_block_faster ();
extern Boolean read_to_block();


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

extern int block_comparison();


/*
 * free_block
 *   Frees the block and the sub elements.
 *   Parameters:
 *     Score *score: the score to free
 *   Return code: none
 *   Error code: none
 */
extern void free_block() ;

extern int read_block_header();
extern void read_block_body();
extern void next_cluster();

/*
 * resize_block_sequences
 *   Increases the memory for the storage of the sequences of a block.
 *   Parameter:
 *     Block *block: the block to resize
 *   Error codes: none
 */
/* Check this prototype, it was originally extern,
 * but it is declared as static in the format_block.c code
 * however there is also a duplicate in the code library
 * -lblimps*/
void resize_block_sequences();

extern void resize_block_clusters();

/*
 * print_block
 *   Prints a block data structure.  Primarily for debugging purposes.
 *   Parameters:
 *     Block *block:  the block to print
 *   Error Codes: none
 */

extern void print_block();

/*
 * ouput_block
 *   Outputs a block data structure to the given file.
 *   Parameters:
 *     Block *block:  the block to print
 *     FILE  *obfp:   the ouput block file pointer
 *   Error Codes: none
 */

extern void output_block();

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

extern void output_block_s();

/*
 *  new_block
 *    Creates a new empty Block structrue
 *    Parameters:
 *      ncols:  Number of columns = sequence length
 *      nrows:	Number of rows = number of sequences
 */

extern Block *new_block();






#endif /*  BLOCKS_H_ */

/* Change log information follows. */
/*
 * Changes since 3.3.2:
  12/15/99  Added double undefined_dbl field.
  12/17/99  Added family field, MINAC, MAXAC
 * Changes since 3.2.5:
   2/11/99  Added min_prev & max_prev fields to the blocks structure;
		rename block->sequence_length to block_width.
 * Changes since 3.2.2:
   1/29/98  Added resize_block_clusters(), etc.
 * Changes since 3.2:
   4/16/97  Added resize_block_sequences()
 * Changes since 3.1:
   2/14/97  Added new_block().
 *
 */
