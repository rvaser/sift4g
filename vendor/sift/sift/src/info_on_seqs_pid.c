/* info_on_seqs.c
1/25/01 gap option changed to false, residue threshold changed to 3
*/
/* Aug 7, 2011 Pauline added percent identity threshold as input */

#include <assert.h>

#include "Alignment.h"
#include "PN_blocks.h"
#include "PN_convert.h"
#include "Matrix_Info.h"
#include "Psiblast.h"
/* #include "Information.c */

#include "info_on_seqs_pid.h"

#define MAXSEQ 400 /* Maximum number of sequences */
#define FALSE 0
#define TRUE 1
#define RESIDUE_THRESHOLD 2
#define NEIGHBOR 3
#define THRESHOLD 0.05
#define ADEQUATE_SEQ_INFO 3.25

void output_predictions (AAnode* polymorph_data, int* residues_stored,
	Matrix_Info* matrix_info, FILE* outfp, Matrix* pssm);

double calculate_info_on_neighbors (int pos, double* info, int* residues_stored,
	int length);

double calculate_scaled_info_on_neighbors (int pos, double* info,
	int* residues_stored, int length);

double mic_with_seq (Block* oldblock, int pos);

void comments_on_info(Block* oldblock, int pos, FILE* outfp, int substitution);

void generate_predictions(Sequence** seqs, int nseqs, FILE* polymorphism_fp,
	int seq_identity, FILE* out_fp) {

    /*
	reduce redundancy with query sequence 01/03/00
	so sequence weighting doesn't include more than 1 sequence
	representing query AND there aren't other really closely
	related sequences that may include the polymorphism
	*/
	remove_seqs_percent_identical_to_query (seqs, &nseqs, (double) seq_identity);

	/* this is an alignment, all sequences should have same length */
	int aa_length = get_length (seqs[0]);

	/* READ POLYMORPHISMS */
	AAnode* polymorph_data = read_multiple_polymorphism_data (polymorphism_fp, seqs[0]);

	Block* block = make_block (aa_length, 0, nseqs, seqs, FALSE);

	double median = calculate_median_information_of_block (block, FALSE, TRUE);

	int* residues_stored = calculate_basic_aa_stored (block);
	struct working* col = make_col();

	int pos;
	for (pos = 0; pos < aa_length; pos++) {
		counts (block, col, pos);
		if (count_residues(col) == 1) {
			// printf ("%d\n", pos+1);
		}
	}
	Matrix* SIFT_pssm = SIFT_prediction (block, TRUE, FALSE, TRUE, FALSE);
	Matrix* counts_pssm = block_to_matrix (block, 2);
	Matrix_Info* matrix_info = initialize_matrix_info (counts_pssm);
	calculate_information_R (matrix_info, TRUE);

	output_predictions(polymorph_data, residues_stored, matrix_info, out_fp, SIFT_pssm);
	/*output_matrix_s (SIFT_pssm, out_fp, FLOAT_OUTPUT); */

	int i;
	for (i = 0; i < aa_length; i++) {
		free_AAnode (polymorph_data[i]);
	}
	free (polymorph_data);

	free_block (block);
	free_matrix (SIFT_pssm);
	free_matrix (counts_pssm);
	free_Matrix_Info (matrix_info);
	free(col);
	free(residues_stored);
}


void
output_predictions (AAnode* polymorph_data, int* residues_stored,
	Matrix_Info* matrix_info, FILE* outfp, Matrix* pssm)
{
	int substitution, original_aa;
	int substitution_exists;
	int pos;
	double median;
	AAnode current;
	Block* block_with_seqs_at_pos;
	double* info_array;
	int block_constructed;

	substitution_exists = FALSE;

	for (pos = 0; pos < pssm->width; pos++) {
		assert (pssm != NULL);
		assert (pssm->block != NULL);
		block_constructed = 0;
		original_aa = pssm->block->residues[0][pos];
	    if (pssm->weights[original_aa][pos] < THRESHOLD) {
			block_constructed = 1;
			block_with_seqs_at_pos = subblock_of_seqs_with_aa_at_pos (pssm->block, pos);
            info_array = calculate_info_for_each_pos (block_with_seqs_at_pos, FALSE);
            median = median_of_array (info_array, pssm->block->width);
			/* there's enough sequence data but aa still predicted
			to be intolerant */
			if (median < ADEQUATE_SEQ_INFO) {
				fprintf (outfp, "WARNING! %c%d not allowed! score: %.2f median: %.2f # of sequence: %d\n",
					aa_btoa[original_aa], pos + 1,
					pssm->weights[original_aa][pos], median, block_with_seqs_at_pos->num_sequences);
            }
		}

		current = polymorph_data[pos]->next; /* nothing in first node*/
		if (current != NULL) {
			/* at least one substitution exists */
			substitution_exists = TRUE;
			if (!block_constructed) {
				block_with_seqs_at_pos = subblock_of_seqs_with_aa_at_pos (pssm->block, pos);
            	info_array =  calculate_info_for_each_pos (block_with_seqs_at_pos, FALSE);
                median = median_of_array (info_array, pssm->block->width);
			}
		while (current != NULL) {
			substitution = current->aa;
			if (residues_stored[pos] >= RESIDUE_THRESHOLD) {

			fprintf (outfp, "%c%d%c\t", aa_btoa[original_aa],pos+1, aa_btoa[substitution]);
				if (pssm->weights[substitution][pos] >= THRESHOLD) {
					fprintf (outfp, "TOLERATED\t%.2f\t", pssm->weights[substitution][pos]);
				} else {
					fprintf (outfp, "DELETERIOUS\t%.2f\t", pssm->weights[substitution][pos]);
				}

				fprintf (outfp, "%.2f\t%d\t%d\n", median,
				block_with_seqs_at_pos->num_sequences,
				pssm->block->num_sequences );
			} else { /* don't have enough aa's in this position to score */
				fprintf (outfp, "%c%d%c\t", aa_btoa[original_aa], pos+1,
				aa_btoa[substitution]);
				fprintf (outfp, "NOT SCORED\n");
			}
			current = current->next;
		} /* end of looking at all subst */

		} /* end of if there was at least one substitution */
		if (block_constructed) {
			free(info_array); free_block (block_with_seqs_at_pos);
			block_constructed = 0;
		}

	} /* end of for pos */

	// printf ("done checking all subst\n");
	if (substitution_exists == FALSE) {
		/* fprintf (outfp, "No substitutions were submitted.<BR>\n"); */
		output_matrix_s (pssm, outfp, FLOAT_OUTPUT);
	}
} /* end of output_predictions procedure */

/* calculating information at a position. no error correction because
want to look at divergence rather than conservation.  so low information
will mean diverged sequences whereas high information due to conservation
or low sequence # means not diverged enough
max information always 4.3
*/

double
mic_with_seq (Block* oldblock, int pos)
{
	/* returns median information content . takes into
	account whether sequences are present at the poistion. */

    double median;
    double* info_array; int prot_length;
	Block* block_with_seqs_at_pos;

	block_with_seqs_at_pos = subblock_of_seqs_with_aa_at_pos (oldblock, pos);
    info_array =  calculate_info_for_each_pos (block_with_seqs_at_pos, FALSE);
    prot_length = oldblock->width;
	median = median_of_array (info_array, prot_length);
    free (info_array);
    free_block (block_with_seqs_at_pos);

} /* end mic_with_seq */

void
comments_on_info(Block* oldblock, int pos, FILE* outfp, int substitution)
{
	double median;
	double* info_array; int prot_length;
	int nseqs_at_pos;
	Block* block_with_seqs_at_pos;
	int i;
	double neighbor_info;
	int* residues_stored;
	Matrix* matrix;

	block_with_seqs_at_pos = subblock_of_seqs_with_aa_at_pos (oldblock, pos);
 	/* matrix = SIFT_prediction (block_with_seqs_at_pos, TRUE, FALSE, TRUE, FALSE);
	printf ("end of block with seqs at pos\n");
	fprintf (outfp, "%d%c %.2f\t", pos, aa_btoa[substitution],
		matrix->weights[substitution][pos]); */

    info_array =  calculate_info_for_each_pos (block_with_seqs_at_pos, FALSE);
	prot_length = oldblock->width;
    /* fprintf (outfp, "%.2f\t", info_array[pos]);
    residues_stored = calculate_basic_aa_stored (oldblock);
	neighbor_info = calculate_info_on_neighbors (pos, info_array, residues_stored, oldblock->width);
    fprintf (outfp, "%.2f\t", neighbor_info); */

    /* max_info = calculate_max_info_possible(pssm->num_sequences); */
    /* fprintf (outfp, "%.2f\t", matrix_info->R[pos]/max_info); */
	/* fprintf (outfp, "%.2f\t", neighbor_info);
    fprintf (outfp, "%d\n", residues_stored[pos]); */

	/* the array will be sorted in call for median_of_array*/
	printf ("looking for median\n");
	median = median_of_array (info_array, prot_length);
	printf ("found mediant %.2f\n", median);
	fprintf (outfp, "%.3f\t%d\t%d\n",
		median, /* information median */
		block_with_seqs_at_pos->num_sequences, /* sequences at that position */
		oldblock->num_sequences /* maximum number of sequences in original block*/
	);
 	free (info_array);
	free_block (block_with_seqs_at_pos);
	assert (oldblock->residues != NULL);
	assert (oldblock->clusters!= NULL);
	assert (oldblock->sequences[0].sequence != NULL);
	assert (oldblock->sequences != NULL);
	assert (oldblock != NULL);
	printf ("finished comments oninfo\n");
}

double
calculate_info_on_neighbors (int pos, double* info, int* residues_stored, int length)
{
	double R;
	int i;

	R= 0.0;
	for (i = pos - NEIGHBOR; i < pos; i++) {
		if (residues_stored[i] >= RESIDUE_THRESHOLD && i >= 0) {
			R += info[i];
		} else {
			return -1000;
		}
	}
	for (i = pos + NEIGHBOR; i > pos; i--) {
		if (residues_stored[i] >= RESIDUE_THRESHOLD && i < length) {
			R += info[i];
		} else {
			return -1000;
		}
	}
	return (R/ (NEIGHBOR *2) );
}

double
calculate_scaled_info_on_neighbors (int pos, double* info, int* residues_stored, int length)
{
    double max_info;
	double R;
	int i;

	R= 0.0;
    for (i = pos - NEIGHBOR; i < pos; i++) {
        if (residues_stored[i] >= RESIDUE_THRESHOLD && i >= 0) {
            max_info = calculate_max_info_possible(residues_stored[i]);
			R += (info[i]/max_info);
        } else {
            return -1000;
        }
    }
    for (i = pos + NEIGHBOR; i > pos; i--) {
        if (residues_stored[i] >= RESIDUE_THRESHOLD && i < length) {
            max_info = calculate_max_info_possible(residues_stored[i]);
			R += (info[i]/max_info);
        } else {
            return -1000;
        }
	}
    return (R / (NEIGHBOR *2));

} /* end of calculate_scaled_info_on_neighbors */

int
number_of_digits (int number)
{
	if (number < 10) {
		return 1;
	} else if (number < 100) {
		return 2;
	} else if (number < 1000) {
		return 3;
	} else if (number < 10000) {
		return 4;
	}
	return 5;
}

void
output_matrix_sPN_web (Matrix* matrix, FILE* omfp, double threshold,
	Sequence* query_seq, double* fraction_stored)
{

	char c; int l;
	int space;
	int dig_number;
	int pos;
    fprintf (omfp, "Each row corresponds to a position in the reference protein.  ");
    fprintf (omfp, "Below each position is the fraction of sequences that contain one of the basic amino acids.  A low fraction indicates the position is either severely gapped or unalignable and has little information.  Expect poor prediction at these positions.<BR>");
	fprintf (omfp, "Each column corresponds to one of the twenty amino acids. <BR>");
    fprintf (omfp, "Each entry contains the score at a particular position (row) for ");
    fprintf (omfp, "an amino acid substitution (column).  Substitutions predicted to");

    fprintf (omfp, " be intolerant are highlighted in red.<BR><BR>\n");
    fprintf (omfp, "<BR>\n");
    fprintf (omfp, "<table cellspacing=0 border=0 width=0 cols=21>\n");

	for (l=0; l<matrix->width; l++) {

	/* every 25 positions, print out amino acids for easier reading*/
	fprintf (omfp, "<tr>");
	if (fmod( (double) l, 25.0) == 0) {
	    fprintf (omfp, "<th>pos</th>\n");
	    for (c='A'; c < 'Z'; c++) {
      		if ((c != 'J') && (c != 'O') && (c != 'U') &&
					(c != 'B') && (c != 'X') ) {
		        fprintf(omfp, "<th>%c</th>\n", c);
      		}
        }
		fprintf (omfp, "</tr>\n");
	}
	space = 24;
	pos = l + matrix->block->sequences[0].position;
	fprintf (omfp, "<tr><th>%2d%c %.2f</th>\n", pos,
		aa_btoa[query_seq->sequence[l]],
		fraction_stored[l]);
    for (c='A'; c < 'Z'; c++) {
    	if ((c != 'J') && (c != 'O') && (c != 'U') && (c != 'B') && (c != 'X') )
		{
        	fprintf (omfp, "<td>\n");
	  		if (matrix->weights[aa_atob[c]][l] < threshold) {
				fprintf (omfp, "<font color=red>%.2f</font>\n", matrix->weights[aa_atob[c]][l]);
	   		} else {
				fprintf (omfp, "%.2f", matrix->weights[aa_atob[c]][l]);
			}
			fprintf (omfp, "</td>\n");
		} /* end of if not amino acid characters*/
      } /* end of for alphabet*/
      fprintf (omfp, "</tr>\n");
    }
	fprintf (omfp, "</table>\n");
} /* end of output_matrix */


/*=====================================================================*/
/* sequences over length threshold */
/* return sequences (in newseqs) that have at least seqs[0]->length* threshold
in alignment. will provide seed for starting blocks. returns the number
of sequences in newseq */

int
sequences_over_length_threshold (Sequence* newseqs[MAXSEQ], Sequence* seqs[MAXSEQ],
	int nseqs, double threshold)
{
	int newseq_index, i, length, query_length;

	query_length = get_length (seqs[0]);
	newseq_index = 0;
	for (i = 0; i < nseqs; i++) {
		length = get_length (seqs[i]);
		if (length >= threshold * query_length) {
			/*newseqs[newseq_index] = copy_sequence (seqs[i]); */
			newseq_index++;
		}
	}

} /* end of sequences_over_length_threshold */
