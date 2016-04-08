/* (C) Copyright 2000, Fred Hutchinson Cancer Research Center */
/* Use, modification or distribution of these programs is subject to */
/* the terms of the non-commercial licensing agreement in license.h. */

#ifndef _ALIGNMENT_H_
#define _ALIGNMENT_H_

#define MAXLEN 800
#define LINE_LEN 800
#define PHYLIP_WIDTH 10
#include <string.h>
#define true 1
#define false 0

#include "blimps/blocksprogs.h"

#define MAXSEQ 400 /* Maximum number of sequences */

extern FILE* errorfp;

/* This header file contains all subroutines involving multipe alignment
	formats */

/*========================================================================
 Try clustal format, which has multiple lines per sequence with the
 sequence name in the first 16 columns of each line.
 Assumes each set of segments is separated by at least one blank line.
 Assumes each line of a set contains the sequence name, whitespace, residues:
============================================================================*/
int try_clustal(FILE* ifp, Sequence* seqs[MAXSEQ], char desc[SMALL_BUFF_LENGTH]);

/*============================================================================
MSF format:  Comments until alignment begins after "//" line:

//

            1                                                   50
P09254-1    .....KRQED AGYDICVPYN LYLKR..... NEFIKIVLPI IRDWDLQHPS
A37470-1    .TFAPKRDED AGYDIAMPYT AVL....... APGENLHVRL PVAYAADAHA
Q00030-2    DYFAPKRDED AGYDISAQTN ATI....... EPDESYFVEL PIVFSSSNPA
============================================================================*/
int try_msf(FILE* ifp, Sequence *seqs[MAXSEQ], char desc[SMALL_BUFF_LENGTH]);

void change_periods_to_dashes (int nseqs, Sequence* seqs[MAXSEQ]);

void remove_gaps (Sequence* newseqs[MAXSEQ], Sequence* seqs[MAXSEQ], int nseqs, int query_length);
/* removes columns that correspond to -'s in query's sequence.  Query must be first
in the list of sequences.  seqs[0] ) */

void initialize_seqs (Sequence* seqs[MAXSEQ], int nseq, int query_length);

void free_seqs (Sequence* seqs[MAXSEQ], int nseq);

/* ========================================================================

 Process using flat master slave without identities, blunt ends
 Option -m 6 for Psiblast

# of sequences returned

QUERY 1   MKPVTLYDVAEYAGVSYQTVSRVVN---Q--A-S---H--VSAKTREKVEAAMAELNYI- 48
30528 1   MKPVTLYDVAEYAGVSYQTVSRVVN---Q--A-S---H--VSAKTREKVEAAMAELNYI- 48
30529 4   -RTATLEDVARRGRVPADGLRRVLN---R--P-E---V--VSARTREQVIRAMQALHYV- 50
30540 3   ----TIKDVAKMAGVSTTTVSHVIN---K--T-R---F--VAKDTEEAVLSAIKQLNYS- 46

QUERY 49  PN--RVAQQL--AG---KQSLLIGV-A---T-SSLALHAP----S-----QIVAAIKS-- 85
30528 49  PN--RVAQQL--AG---KQSLLIGV-A---T-SSLALHAP----S-----QIVAAIKS-- 85
30529 51  PN--RSAQLL--AG---KAAPSIGL-I---T-ASVTLHAP----S-----QIAAAIKS-- 87
30540 47  PS--AVARSL--KV---NTTKSIGM-I---V-TTSEAPYF----A-----EIIHSVEE-- 83

=======================================================================*/
int flat_master_slave (FILE* fp, Sequence *seqs[MAXSEQ], char desc[SMALL_BUFF_LENGTH]);

/*======================================================================
read_psiblast_header_until_last:
Processes lines until reach the beginning of last alignment.
Returns true if alignment has converged before then, false otherwise

=====================================================================*/

int
read_psiblast_header_until_last(FILE* fp, int max_iterations);

/* read_psiblast_header_until_first
Processes lines until first alignment
*/
void read_psiblast_header_until_first (FILE* fp);

int
get_length (Sequence* seq);

/*======================================================================
   Fix up the sequence names
=======================================================================*/
void fix_names(int nseq, Sequence* seqs[MAXSEQ]);

void output_sequence_nonfasta();

void
extract_seqs (Sequence* newseqs[MAXSEQ], Sequence* oldseqs[MAXSEQ], int nseqs, char names[MAXSEQ][SMALL_BUFF_LENGTH], int new_no_of_seq);
/* copies the sequences contained in the list of names from oldseqs into newseqs.  new_no_of_seq is the number of names in list of names, hence newseqs has new_no_of_seq entries */

void copy_sequence (Sequence* newseq, Sequence* oldseq);

int
compare_alignments (Sequence* align1[MAXSEQ], Sequence* align2[MAXSEQ],
		    int align1_nseq, int align2_nseq);

/*========================================================================
makes a block from a region of sequences. block starts from firstpos
(the position in array, not sequence) and the width is passed in.
if seq_weights_option == TRUE, use Block's weights, otherwise if == FALSE,
calculate
========================================================================*/
Block *make_block(int width, int firstpos, int nseq, Sequence* seqs[MAXSEQ],
		  int use_seq_weights_option);

/* makes a block out of a seqs array.  gaps kept */
Block* block_from_seq_given_pos (Sequence* seqs[MAXSEQ], int nseq, int beg, int end);

int width_with_gaps (Sequence* seq, int beg,int end);

int get_gapped_pos (Sequence* seq, int aa_pos);
/* given amino acid position in ungapped sequence, return index in seq
(which has gaps) that corresponds to that amino acid position */

/* switch QUERY to be in seqs[0] */
void make_query_first (Sequence* seqs[MAXSEQ], int nseqs);

/* subroutine to merge local alignments from psiblast output
If there is an aa in one pos and a - in the corresponding pos in the
other sequence, choses the aa over the -.  There should not be an aa
in both seq1 and seq2 at the same position, this will return an error.
Change seq1.
*/
void merge_res (char* seq1, char* seq2);

/* if option_carve_gaps = true, assign sequence amino acid only if
there is NOT a gap in the query sequence else if option_carve_gaps ==
FALSE, allow gaps in query sequence when processing alignments */

int psiblast_pairwise (FILE* fp, Sequence * seqs[MAXSEQ], int length,
			int option_carve_gaps);

void psiblast_pairwise_more_alignments (FILE* fp, Sequence* subject, int* done,
					int option_carve_gaps);

char* copy_string (char* src);

Sequence*
read_1st_4lines_of_pairwise (FILE* fp, char subject_name[SMALL_BUFF_LENGTH],
                                int length, int option_carve_gaps);

int read_4lines_of_pairwise (FILE* fp, Sequence* subject, int option_carve_gaps);

void convert_gap_to_X (Sequence* seq);

void output_sequence_without_gaps_or_Xes (Sequence* seq, FILE* osfp);

void convert_sequence_with_beg_and_end_gaps_to_X (Sequence* seq);

/* just like Jorja's output_sequence, but doesn't allow a newline if
the length is a multiple of 60 (so that a series of fasta sequences
will not be separated by any newlines */
void output_sequence_clean (Sequence* seq, FILE* outfp);

int get_index_from_seqs (char name[], Sequence* seqs[MAXSEQ],
                        int nseq);

int get_length_nogap_noX (Sequence* seq);

/* returns # of clusters */
int cluster (double clus, Sequence* seqs[MAXSEQ], int nseq, FILE* outfp, FILE* keyoutfp, int option);

void output_sequence_without_beg_and_end_X (Sequence* seq, FILE* osfp);

int seq_not_all_Xes (Sequence* seq);

void sort_seq_by_identity (Sequence* seqs[MAXSEQ], int num_seqs);

#endif
