/* info_on_seqs.c
1/25/01 gap option changed to false, residue threshold changed to 3
*/
/* Aug 7, 2011 Pauline added percent identity threshold as input */

#define EXTERN

#include <assert.h>

#include "Alignment.h"
#include "PN_blocks.h"
#include "PN_convert.h"
#include "Matrix_Info.h"
#include "Psiblast.h"
/* #include "Information.c */

#define MAXSEQ 400 /* Maximum number of sequences */
#define LINE_LEN 800
#define FALSE 0
#define TRUE 1
#define RESIDUE_THRESHOLD 2
#define NEIGHBOR 3
#define THRESHOLD 0.05
#define ADEQUATE_SEQ_INFO 3.25

/* Local routines */

void getargs (int argc, char* argv[], FILE** seqfp, FILE** polymorphfp,
			char outfilename[LARGE_BUFF_LENGTH], int* seq_identity);

void
output_predictions (AAnode* polymorph_data, int* residues_stored,
                    Matrix_Info* matrix_info, FILE* outfp, Matrix* pssm);

double
calculate_info_on_neighbors (int pos, double* info,
                             int* residues_stored, int length);

double calculate_scaled_info_on_neighbors (int pos, double* info,
                                   int* residues_stored,
                                int length);

double mic_with_seq (Block* oldblock, int pos);

FILE* errorfp;
char errorfilename[LARGE_BUFF_LENGTH];

void comments_on_info(Block* oldblock, int pos, FILE* outfp, int substitution);

/* MAIN */

int main
(int argc, char* argv[])
{
	FILE* seqfp; FILE* outfp; FILE* polymorphismfp;
	char outfilename[LARGE_BUFF_LENGTH], currentoutfile[LARGE_BUFF_LENGTH];
	char tempname[LARGE_BUFF_LENGTH];
	char desc[SMALL_BUFF_LENGTH]; int desc_length; char* strptr;
	Sequence *seqs[MAXSEQ];
	int nseqs, aa_length, i, pos;
	Block* block;
	struct working* col;
	int db_type;
	int seq_type;
	Matrix* SIFT_pssm; Matrix* counts_pssm;
	int seq_identity;
	Matrix_Info* matrix_info; /* will store counts pssm and information values */
	AAnode* polymorph_data;
	int* residues_stored; /* basic aa/total number of seq at each pos */
	double median;
	int original_aa;

	ErrorLevelReport = 5;

	init_frq_qij();
	printf ("tell me i've entered\n");
	getargs (argc, argv, &seqfp, &polymorphismfp, outfilename,
							&seq_identity );

        if ((outfp = fopen (outfilename, "w")) == NULL)
        {
                printf ("cannot open file %s \n", outfilename);
                exit (-1);
        }

	nseqs = 0;
     /*-----------------------------------------------------------------*/
      /*   Check next for input file of sequences & assume are aligned */
      db_type = type_dbs(seqfp, DbInfo);
      /*   could set db_type = FLAT if it comes back negative   */
      seq_type = UNKNOWN_SEQ;
      seq_type = seq_type_dbs(seqfp, DbInfo, db_type, seq_type);
      if (seq_type == NA_SEQ)
      {
         fprintf(stderr, "WARNING: Sequences appear to be DNA but will be treated as protein\n");
         seq_type = AA_SEQ;
      }
      rewind(seqfp);
      /*-----------------------------------------------------------------*/
      /*   read fasta sequences into memory                    */
      if (db_type >= 0)
      {
         while ( nseqs < MAXSEQ &&
             (seqs[nseqs] = read_a_sequence(seqfp, db_type, seq_type)) != NULL)
         {
            nseqs++;
         }
      }
      else /*  CLUSTAL or MSF? */
      {
         /* get tail of outfilename, will use this as a description to make
	a mablock file */
	strcpy (tempname, outfilename);
	strptr = strtok (tempname, "/");
	while ((strptr = strtok ((char*) NULL, "/ \t\n\r\0")) != NULL) {
		desc_length = strlen(strptr);
		if (desc_length > SMALL_BUFF_LENGTH) {
			desc_length = SMALL_BUFF_LENGTH;
		}
		strncpy (desc, strptr, desc_length);
		desc[desc_length] = '\0';
	}
	 nseqs = try_clustal(seqfp, seqs, desc);
         if (nseqs <= 0)
         { nseqs = try_msf(seqfp, seqs, desc);
	   change_periods_to_dashes(nseqs, seqs);
          }
          for (i = 0; i < nseqs; i++) {
                convert_sequence_with_beg_and_end_gaps_to_X (seqs[i]);
          }
	}
      fix_names(nseqs, seqs);
	if (nseqs == 0) {
		fprintf (errorfp, "Cannot read sequences.  Check format.\n" );
		exit (-1);
	}
	if (nseqs == MAXSEQ)
	{
		fprintf (stderr, "WARNING: Maximum number of sequences = %d\n", nseqs);
	}
        /*******************
		 reduce redundancy with query sequence 01/03/00
           so sequence weighting doesn't include more than 1 sequence
           representing query AND there aren't other really closely
	   related sequences that may include the polymorphism */
        remove_seqs_percent_identical_to_query (seqs, &nseqs,
					(double) seq_identity);

	aa_length = get_length (seqs[0]);

	/* this is an alignment, all sequences should have same length */

	/* READ POLYMORPHISMS */
	polymorph_data = read_multiple_polymorphism_data (polymorphismfp, seqs[0]);
	if (polymorphismfp != NULL) {
		fclose (polymorphismfp);
	}

	block = make_block (aa_length, 0, nseqs, seqs, FALSE);

	median = calculate_median_information_of_block (block,
				FALSE, TRUE);
	residues_stored = calculate_basic_aa_stored (block);
	col = make_col();
	for (pos = 0; pos < aa_length; pos++) {
		counts (block, col, pos);
		if (count_residues(col) == 1) {
			printf ("%d\n", pos+1);
		}
	}
	SIFT_pssm = SIFT_prediction (block, TRUE, FALSE, TRUE, FALSE);
	counts_pssm = block_to_matrix (block, 2);
	matrix_info = initialize_matrix_info (counts_pssm);
	calculate_information_R (matrix_info, TRUE);

	printf ("about to make predictions\n");
	output_predictions(polymorph_data, residues_stored,
                    matrix_info, outfp, SIFT_pssm);
	/*output_matrix_s (SIFT_pssm, outfp, FLOAT_OUTPUT); */
	printf ("trying to free things here\n");
	fclose (outfp);

        for (i = 0; i < aa_length; i++) {
                free_AAnode (polymorph_data[i]);
        }
	free (polymorph_data);

	free_block (block);
        free_matrix (SIFT_pssm); free_matrix (counts_pssm);
	free_Matrix_Info (matrix_info);
	fclose (seqfp);
	/*free_seqs (seqs, nseqs); sequences have been previously freed */
	fclose (errorfp);
	rm_file (errorfilename);
	exit (0);

} /* end main */

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
			block_with_seqs_at_pos = subblock_of_seqs_with_aa_at_pos
                                                        (pssm->block, pos);
                        info_array =  calculate_info_for_each_pos
                                (block_with_seqs_at_pos, FALSE);
                        median = median_of_array (info_array, pssm->block->width
);
			/* there's enough sequence data but aa still predicted
			to be intolerant */
			if (median < ADEQUATE_SEQ_INFO) {                                			fprintf (outfp, "WARNING! %c%d not allowed! score: %.2f median: %.2f # of sequence: %d\n",
                                                 aa_btoa[original_aa], pos + 1,
				pssm->weights[original_aa][pos], median, block_with_seqs_at_pos->num_sequences	);
                	}
		}

		current = polymorph_data[pos]->next; /* nothing in first node*/
		if (current != NULL) {
			/* at least one substitution exists */
			substitution_exists = TRUE;
			if (!block_constructed) {
			block_with_seqs_at_pos = subblock_of_seqs_with_aa_at_pos
                                                        (pssm->block, pos);
                	info_array =  calculate_info_for_each_pos
                                (block_with_seqs_at_pos, FALSE);
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

	printf ("done checking all subst\n");
	if (substitution_exists == FALSE) {
/*		fprintf (outfp, "No substitutions were submitted.<BR>\n"); */
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
        info_array =  calculate_info_for_each_pos (block_with_seqs_at_pos,
                 FALSE);
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
 /*       matrix = SIFT_prediction (block_with_seqs_at_pos, TRUE, FALSE, TRUE, FALSE);
	printf ("end of block with seqs at pos\n");
         fprintf (outfp, "%d%c %.2f\t", pos, aa_btoa[substitution],
				 matrix->weights[substitution][pos]);
*/
        info_array =  calculate_info_for_each_pos (block_with_seqs_at_pos,
                 FALSE);
	prot_length = oldblock->width;
/*        fprintf (outfp, "%.2f\t", info_array[pos]);
        residues_stored = calculate_basic_aa_stored (oldblock);
	neighbor_info = calculate_info_on_neighbors
                                (pos, info_array, residues_stored, oldblock->width);
       fprintf (outfp, "%.2f\t", neighbor_info);
*/
/*       max_info = calculate_max_info_possible(pssm->num_sequences); */
/*                                fprintf (outfp, "%.2f\t", matrix_info->R[pos]/max_info); */
/*       fprintf (outfp, "%.2f\t", neighbor_info);
       fprintf (outfp, "%d\n", residues_stored[pos]);
*/
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
calculate_info_on_neighbors (int pos, double* info,
			     int* residues_stored, int length)
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
calculate_scaled_info_on_neighbors (int pos, double* info,
				   int* residues_stored,
				int length)
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
        return (R/ (NEIGHBOR *2) );

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
    fprintf (omfp, "an amino acid substitution (column).  Substitutions predicted to"
);
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
/* return sequences (in newseqs ) that have at least seqs[0]->length* threshold
in alignment .  will provide seed for starting blocks. returns the number
of sequences in newseq*/

int
sequences_over_length_threshold
(Sequence* newseqs[MAXSEQ], Sequence* seqs[MAXSEQ], int nseqs, double threshold)
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


void getargs (int argc, char* argv[], FILE** seqfp, FILE** polymorphfp,
			char outfilename[LARGE_BUFF_LENGTH], int* seq_identity )
{
	char seqfilename[LINE_LEN];
	char seq_outfilename[LINE_LEN];
	char substfilename[LARGE_BUFF_LENGTH];

	if (argc < 5)
	{
		printf ("info_on_seqs_pid\n");

	}

	if (argc > 1) strcpy (seqfilename, argv[1]);
	else
	{
		printf ("Enter filename with sequences:\n");
		fgets (seqfilename, LINE_LEN, stdin);
	}
printf ("fawegwa\n");
	if ((*seqfp = fopen (seqfilename, "r")) == NULL)
	{
		printf ("cannot open file %s \n", seqfilename);
		exit (-1);
	}
printf ("eaegrtjkl\n");
	if (argc > 2) strcpy (substfilename, argv[2]);
	else {
		printf ("Enter file with substitutions.\n");
		printf ("format:\n");
		printf ("M1Y\n");
		printf ("S2K\n");
		printf ("...\n\n");
		fgets (substfilename, LINE_LEN, stdin);
	}
        if (substfilename[0] == '-') { *polymorphfp = NULL; }
	else if ((*polymorphfp = fopen (substfilename, "r")) == NULL)
        {
                printf ("cannot open file %s \n", substfilename);
                exit (-1);
        }


	if (argc > 3) strcpy (outfilename, argv[3]);
	else
	{
		printf ("Enter name of outfile\n");
		fgets (outfilename, LINE_LEN, stdin);
	}

        strcpy (errorfilename, outfilename);
        strcat (errorfilename, ".error");
        if ((errorfp = fopen (errorfilename, "w")) == NULL) {
                printf ("couldn't open file %s\n", errorfilename);
                exit (-1);
        }

	if (argc > 4) *seq_identity = atoi (argv[4]);
        else {
                printf ("100% sequence identity that will be filtered out by default because argument not set\n");
                *seq_identity = 100; /* default of 100% */
                /*scanf ("%d", seq_identity); */
        }


} /* end of getargs */
