/* (C) Copyright 2000, Fred Hutchinson Cancer Research Center */
/* Use, modification or distribution of these programs is subject to */
/* the terms of the non-commercial licensing agreement in license.h. */

#ifndef _PN_BLOCKS_H_
#define _PN_BLOCKS_H_

#define FALSE 0
#define TRUE 1
#include "stringhash.h"

#include "blimps/blocksprogs.h"

#define MAXSEQ 400 /* Maximum number of sequences */

extern FILE* errorfp;

Block* copy_block (Block* block);

/* add sequences from block 2 to block 1 and reweight */
void add_block2_to_block1 (Block* block1, Block* block2);

/* returns sequence from block */
Sequence* seq_from_block (char seq_id[30], Block* block);

Block* extract_seqs_from_old_block (HashTable seqs_names, Block* block);

/* seqname corresponds to a sequence in original block.  This sequence is added
to block */
void add_seq_to_block (Block* block, Block* originalblock, char seqname[]);

void add_sequence_to_block (Block* block, Sequence* seq, int weight_block);

void remove_last_seq_of_block (Block* block);

/* information / block->length is returned */
double information_per_residue (Block* block);

/* percentage identity is assigned to percentile */
void calculate_percentage_identity (Block* block);

/* puts id names in seqnamehash of those block->sequences that have at least percentage identicalwith sequence 0 (supposed to be the query ) */
void
percentage_identity_with_seq0 (Block* block, HashTable seqnamehash, double percentage );

/* converts all gap characters '-' to 'X' */
void convert_gap_to_X_block (Block* block);

Block* remove_seq0_Xes_from_block (Block* block);

void percentage_identity_with_seq0_seqs (Sequence* seqs[MAXSEQ], int nseqs);

void remove_seqs_percent_identical_to_query (Sequence* seqs[MAXSEQ], int* no_of_seqs,
				double percent_identical);

Block* subblock_of_seqs_with_aa_at_pos (Block* block, int test_position);

double* calculate_info_for_each_pos (Block* block, int error_correction_option);

void print_block_sequences (Block* block, FILE* outfp);

void print_block_ids (Block* block, FILE* fp, int print_first_id);

void copy_weights (Block* block, Block* oldblock);

double calculate_median_info_for_pos_in_block (Block* block, int pos);

double calculate_median_information_of_block (Block* block, int error_correction,
	int suppress_pbweights);

#endif
