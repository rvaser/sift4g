/* (C) Copyright 2000, Fred Hutchinson Cancer Research Center */
/* Use, modification or distribution of these programs is subject to */
/* the terms of the non-commercial licensing agreement in license.h. */

#ifndef _PN_BLOCKS_C_
#define _PN_BLOCKS_C_
/* 02/07/01 in copy_block, added copyoing weight of sequences */

#include <assert.h>

#include "PN_blocks.h"
#include "PN_convert.h" /* so that information_per_residue can still operate on
			wide blocks */

#include "Array_math.h"

Block*
copy_block (Block* block)
{
        Block* newblock;
	Residue* residue_pointer;
	int i, pos;

        newblock = (Block*) malloc (sizeof (Block));

         /* copy header stuff */
        strcpy (newblock->id, block->id);
        strcpy (newblock->ac, block->ac);
        strcpy (newblock->de, block->de);
        strcpy (newblock->bl, block->bl);
        strcpy (newblock->number, block->number);
        strcpy (newblock->motif, block->motif);
        newblock->width = block->width;
        newblock->strength = block->strength;
        newblock->max_sequences = block->max_sequences;
        newblock->num_sequences = block->num_sequences;
        newblock->max_clusters = 1; /* allocate space for a single cluster*/
        newblock->num_clusters = 1;
        newblock->min_prev = block->min_prev;
        newblock->max_prev = block->max_prev;

	newblock->clusters= (Cluster *) calloc (newblock->max_clusters,
                                                sizeof (Cluster));
        newblock->clusters[0].num_sequences = newblock->num_sequences;

        /* sequences is an array of Sequence structures */
        newblock->sequences = (Sequence *) calloc (block->max_sequences,
                                                sizeof (Sequence));
        residue_pointer = (Residue *) calloc (block->width * block->max_sequences,
		 sizeof (Residue));

        newblock->residues = (Residue**) calloc (block->max_sequences, sizeof
	                   (Residue *));


        newblock->clusters[0].sequences = &(newblock->sequences[0]);

/* initialize residues */
        for (i=0; i < block->max_sequences; i++) {
                newblock->sequences[i].sequence = &(residue_pointer[i*block->width]);
                newblock->residues[i] = &(residue_pointer[i*block->width]);
        }

        for (i = 0; i < block->num_sequences; i++) {
                /* copy information stuff */
                strcpy (newblock->sequences[i].name,
                            block->sequences[i].name);
                strcpy (newblock->sequences[i].info,
                            block->sequences[i].info);
                newblock->sequences[i].position =
                            block->sequences[i].position;
                newblock->sequences[i].length =
                            block->sequences[i].length;
                newblock->sequences[i].max_length=
                            block->sequences[i].max_length;
                newblock->sequences[i].type =
                            block->sequences[i].type;
		newblock->sequences[i].weight =
			    block->sequences[i].weight;

                        /* copy residues */
                for (pos = 0; pos < newblock->width; pos++) {
                            newblock->sequences[i].sequence[pos] =
                                        block->sequences[i].sequence[pos];
/*                            printf ("%c", aa_btoa[newblock->sequences[i].sequence[pos]] ); */
                }
        } /* end of for, gone through all sequences in oldblock*/
        return (newblock);

} /* end of copy_block */

/* returns sequences that have a real amino acid at position
The crude way but tried to do pointer reassignment and couldn't
get it to work.
*/
Block*
subblock_of_seqs_with_aa_at_pos (Block* originalblock, int test_position)
{
	Block* newblock;
	HashTable seqnamehash;
	int seq;
	char* key;

	seqnamehash = InitializeTable (100);
	for (seq = 0; seq < originalblock->num_sequences; seq++) {
                  if (originalblock->residues[seq][test_position] >= 1 &&
		      originalblock->residues[seq][test_position] <=  22) {
			Insert (originalblock->sequences[seq].name, seqnamehash);
  		  } else {
			printf ("not including %s with %c at %d\n",
				originalblock->sequences[seq].name,
			aa_btoa[originalblock->residues[seq][test_position]],
			test_position);

		}
        }
	newblock = extract_seqs_from_old_block (seqnamehash, originalblock);
	assert (newblock->num_sequences == seqnamehash->no_of_elements);
	DestroyTable (seqnamehash);
	pb_weights (newblock);
	assert (newblock->residues != NULL);
	assert (newblock->clusters != NULL);
	assert (newblock->sequences[0].sequence != NULL);
	assert (newblock->sequences != NULL);
	assert (newblock != NULL);

	return (newblock);

}

void
convert_gap_to_X_block (Block* block)
{
	int i, pos;

	for (i = 0; i < block->num_sequences; i++) {
		for (pos = 0; pos < block->width; pos++) {
			if (block->sequences[i].sequence[pos] == aa_atob['-']) {
				block->sequences[i].sequence[pos] =aa_atob['X'];
			}
		}
	}
} /* end of convert_gap_to_X */

void
add_block2_to_block1 (Block* block1, Block* block2)
{
	void *tmp_ptr;
	int starting_index, new_block_index, pos;
	int i, j, s;

	starting_index = block1->num_sequences;
	/* going to start adding sequences at this index */

	block1->max_sequences += block2->num_sequences;
	block1->num_sequences += block2->num_sequences;
	/* Pauline! not in resize_blocks */
	block1->clusters[0].num_sequences += block2->num_sequences;

	tmp_ptr = (Sequence *) realloc (block1->sequences,
		   sizeof (Sequence) * block1->max_sequences);

	block1->sequences = (Sequence *) tmp_ptr;

	tmp_ptr = (Residue *) realloc (block1->sequences[0].sequence,
		sizeof (Residue) * block1->width * block1->max_sequences);

	block1->sequences[0].sequence = (Residue *) tmp_ptr;
        tmp_ptr =   (Residue **) realloc(block1->residues,
                                sizeof(Residue *) *
                                block1->max_sequences);
	  block1->residues = (Residue **) tmp_ptr;
  /* re-initialize sequences and residues, same data, just different place */
  for(i=0; i<block1->max_sequences; i++) {
    block1->sequences[i].length = block1->width;
    block1->sequences[i].sequence =
      &(block1->sequences[0].sequence[i * block1->width]);
    block1->residues[i] =
      &(block1->sequences[0].sequence[i * block1->width]);
  }
  s=0;
  for (j=0; j<block1->num_clusters; j++)
  {
      block1->clusters[j].sequences = &(block1->sequences[s]);
      for (i = 0; i < block1->clusters[j].num_sequences; i++)
		s++;
   }
  /* add in the new sequences */
/*  printf ("%d new blockindex %d block1 num sequences%d\n", new_block_index,
			block1->num_sequences); */
    for (i = 0, new_block_index=starting_index; new_block_index
					 < block1->num_sequences;
					 i++, new_block_index++) {
            /* copy information stuff */
	strcpy (block1->sequences[new_block_index].name,
                       block2->sequences[i].name);
            strcpy (block1->sequences[new_block_index].info,
                       block2->sequences[i].info);
            block1->sequences[new_block_index].position =
                       block2->sequences[i].position;
            block1->sequences[new_block_index].length =
                       block2->sequences[i].length;
            block1->sequences[new_block_index].max_length=
                       block2->sequences[i].max_length;
            block1->sequences[new_block_index].type =
                       block2->sequences[i].type;
             /* copy residues */
/*       printf ("width %d\n", block1->width); */
       for (pos = 0; pos < block1->width; pos++) {
                       block1->sequences[new_block_index].sequence[pos] =
                                        block2->sequences[i].sequence[pos];
/*            	printf ("%c", aa_btoa[block1->sequences[new_block_index].sequence[pos]]); */
	 }
      } /* end of for, gone through all sequences in oldblock*/

	pb_weights (block1);


} /*  end of add block2 to block1 */

Block*
extract_seqs_from_old_block (HashTable seqs_names, Block* block)
{
        char key[SMALL_BUFF_LENGTH];
	int min_string_length;
	int new_block_index;
	Residue* residue_pointer; int i;
	int pos;
	Residue res;

	Block* newblock;

        newblock = (Block*) malloc (sizeof (Block));

	 /* copy header stuff */
	strcpy (newblock->id, block->id);
        strcpy (newblock->ac, block->ac);
        strcpy (newblock->de, block->de);
        strcpy (newblock->bl, block->bl);
        strcpy (newblock->number, block->number);
        strcpy (newblock->motif, block->motif);
        newblock->width = block->width;
        newblock->percentile = block->percentile;
        newblock->strength = block->strength;
	newblock->max_sequences = block->max_sequences;
	newblock->num_sequences = seqs_names->no_of_elements; /*HASH*/
	newblock->max_clusters = 1; /* allocate space for a single cluster*/
	newblock->num_clusters = 1;
	newblock->min_prev = block->min_prev;
	newblock->max_prev = block->max_prev;

	newblock->clusters= (Cluster *) calloc (newblock->max_clusters,
						sizeof (Cluster));
	newblock->clusters[0].num_sequences = newblock->num_sequences;

	new_block_index = 0;
	/* sequences is an array of Sequence structures */
	newblock->sequences = (Sequence *) calloc (block->max_sequences,
						sizeof (Sequence));
	residue_pointer = (Residue *) calloc (block->width * block->max_sequences, sizeof (Residue));

	newblock->residues = (Residue**) calloc (block->max_sequences, sizeof (Residue *));

	newblock->clusters[0].sequences = &(newblock->sequences[0]);

/* initialize residues */
	for (i=0; i < block->max_sequences; i++) {
		newblock->sequences[i].sequence = &(residue_pointer[i*block->width]);
		newblock->residues[i] = &(residue_pointer[i*block->width]);
	}

	for (i = 0; i < block->num_sequences; i++) {
		if (Exists (block->sequences[i].name, seqs_names) ) {
			/* copy information stuff */

			strcpy (newblock->sequences[new_block_index].name,
				    block->sequences[i].name);
                        strcpy (newblock->sequences[new_block_index].info,
                                    block->sequences[i].info);
			newblock->sequences[new_block_index].position =
				    block->sequences[i].position;
                        newblock->sequences[new_block_index].length =
                                    block->sequences[i].length;
                        newblock->sequences[new_block_index].max_length=
                                    block->sequences[i].max_length;
                        newblock->sequences[new_block_index].type =
                                    block->sequences[i].type;

			/* copy residues */
			for (pos = 0; pos < newblock->width; pos++) {
				res = block->sequences[i].sequence[pos];
				newblock->sequences[new_block_index].sequence[pos] = res;
	/*			printf ("%c", aa_btoa[newblock->sequences[new_block_index].sequence[pos]] ); */
			}
			new_block_index++;
		} /* end of if Exists */
	} /* end of for, gone through all sequences in oldblock*/

	/* assert that all the sequence names were found in the block */
	/* changed assert from == to >= Pauline 01/19/04 */
	assert (new_block_index >= seqs_names->no_of_elements);

	pb_weights (newblock);
	return (newblock);

} /* end of extract_seqs_from_old_block */

int
seqindex (Block* originalblock, char seqname[])
{
        int seqindex;
        int i, pos;
        int found;

        found = FALSE;
        seqindex = 0;
        while (seqindex < originalblock->num_sequences && !found) {
                if (strstr (originalblock->sequences[seqindex].name, seqname) !=
 NULL){
                        found = TRUE;
                        seqindex--; /* to compensate for increase in i
                                        after if */
                }
                seqindex++;
        }
        if (!found) {
                fprintf (errorfp, "couldn't find %s\n", seqname);
		exit (-1);
	}
	return seqindex;
}

void
add_seq_to_block (Block* block, Block* originalblock, char seqname[])
{
	int seqindex;
	int i, pos;
	int found;
	int newseqindex;

	found = FALSE;
	seqindex = 0;
	while (seqindex < originalblock->num_sequences && !found) {
		if (strstr (originalblock->sequences[seqindex].name, seqname) != NULL){
			found = TRUE;
			seqindex--; /* to compensate for increase in i
					after if */
		}
		seqindex++;
	}
	if (!found) {
		fprintf (errorfp, "couldn't find %s\n", seqname);
		exit (-1);
	}
	newseqindex = block->num_sequences;

	/* insert sequence as last item in the array */
	/* copy information stuff */
       strcpy (block->sequences[newseqindex].name,
                            originalblock->sequences[seqindex].name);
       strcpy (block->sequences[newseqindex].info,
                            originalblock->sequences[seqindex].info);
       block->sequences[newseqindex].position =
                            originalblock->sequences[seqindex].position;
       block->sequences[newseqindex].length =
                            originalblock->sequences[seqindex].length;
       block->sequences[newseqindex].max_length=
                            originalblock->sequences[seqindex].max_length;
       block->sequences[newseqindex].type =
                            originalblock->sequences[seqindex].type;
                        /* copy residues */
       for (pos = 0; pos < block->width; pos++) {
       /*  printf ("%c\n", aa_btoa[originalblock->sequences[seqindex].sequence[pos]]); */
          block->sequences[newseqindex].sequence[pos]=
                              originalblock->sequences[seqindex].sequence[pos];
       }
	block->num_sequences++;
	if (block->num_clusters != 1) {
                fprintf (errorfp, "there is more than one cluster\n");
                fprintf (errorfp, "add seq subroutine will fail\n");
                exit (-1);
        }
	block->clusters[0].num_sequences++;
	pb_weights(block);

} /* end of add_seq_to_block*/

Block* new_block_with_1_seq (Sequence seq)
{
	Block* block;
	int pos;

	block = new_block (seq.length, 1);
        /* copy information stuff */
       block->percentile = 100;
	block->clusters[0].num_sequences = 1;

	strcpy (block->sequences[0].name,
                            seq.name);
       strcpy (block->sequences[0].info,
                            seq.info);
       block->sequences[0].position =
                            seq.position;
       block->sequences[0].length =
                        seq.length;
       block->sequences[0].max_length=
                            seq.max_length;
       block->sequences[0].type =
                         seq.type;
	/* it matches itself */

                        /* copy residues */
       for (pos = 0; pos < seq.length; pos++) {
          block->sequences[0].sequence[pos]=
                              seq.sequence[pos];

	}
	return block;

}

void
remove_last_seq_of_block (Block* block)
{
	block->num_sequences--;
	block->clusters[0].num_sequences--;
	pb_weights (block);
}

double
information_per_residue (Block* block)
{
	Matrix* matrix;
	double ln2, e, r, hmax, totalweight, dtemp, totalr, temp;
        int aa, pos, seq, num_residues;
        double max_aa_height, aa_height;
	int AA, width;

	AA = 20;

	totalr = 0.0;
	matrix = PN_block_to_matrix (block, 2);
	ln2 = log (2.0);
        hmax = log((double) AA);

	 for(pos=0; pos < matrix->width; pos++)
        {
        /* =====Calculating error correction for small sample size ===*/
       /*  count actual number of residues in each column.
        Including the 20 aa ('residues' matrix values 1-20)
        and B (Asp or Asn, value 21) and Z (Glu or Gln, value 22)
        Excluding gaps (-, value 0) unidentified aa (X, value 23) and
        stop codon (*, value 24). */

        num_residues = matrix->num_sequences ;
	/* for sample sizes > 50 */

        for(seq=0; seq < matrix->num_sequences; seq++)
            if (matrix->block->residues[seq][pos] < 1 ||
                matrix->block->residues[seq][pos] > 22) num_residues-- ;
	if (num_residues > 0) { /* to prevent counting columns that are
                                all X's or all -'s */

                               /* correction factor for small sample size */
        e = (AA-1) / ((double) 2 * ln2 * num_residues) ;

        /* =====end of calculating error correction ================*/
        assert (matrix->weights != NULL);

        r = hmax ; /* R = log 20 */
        for(aa=1, totalweight=0.; aa < AA+1; aa++) {
           totalweight +=  matrix->weights[aa][pos]  ;
        }
        /* totalweight is the sum of all aa weights at a certain pos */
        for(aa=1; aa < AA+1; aa++)
                    /* start loop at 1 and end at AA+1 because aa values
                      at the weight matrix start at 1, 0 is the gap value */
       {
         if (matrix->weights[aa][pos] > 0)
            {
            dtemp = (double) matrix->weights[aa][pos] / totalweight;
           /* in Tom Schneider's logo eq., dtemp corresponds to the frequency
                of an amino acid */
        /* dtemp is the proportion of that particular aa relative to all
         the aa's (total weight ) */
 temp = dtemp * log (dtemp) / ln2;
            r += dtemp * log(dtemp); /* sum over all residues dtemp*log(dtemp)*/
                                    /* is H(pos) */
            } /* end of if weights > 0 */
         } /* end of for for calculating r*/

      r /= ln2 ;       /* convert to bits , R = 2 - H(pos) in bits now */
      r -= e ;                     /* a correction for small sample sizes */
	totalr += r;

	} /* end if num_residues > 0 */
   } /* end of for-loop for all positions */
	width = matrix->width;
	free_matrix (matrix);
	return (totalr / width);
} /* end of information_per_residue  */

double*
calculate_info_for_each_pos (Block* block, int error_correction_option)
{
        Matrix* matrix;
        double ln2, e, r, hmax, totalweight, dtemp, totalr, temp;
        int aa, pos, seq, num_residues;
        double max_aa_height, aa_height;
        int AA, width;
	double* info_array;

        AA = 20;

	info_array = (double *) calloc (block->width, sizeof (double));
        matrix = PN_block_to_matrix (block, 2);
        ln2 = log (2.0);
        hmax = log((double) AA);

         for(pos=0; pos < matrix->width; pos++)
        {
        /* =====Calculating error correction for small sample size ===*/
       /*  count actual number of residues in each column.
        Including the 20 aa ('residues' matrix values 1-20)
        and B (Asp or Asn, value 21) and Z (Glu or Gln, value 22)
        Excluding gaps (-, value 0) unidentified aa (X, value 23) and
        stop codon (*, value 24). */
        num_residues = matrix->num_sequences ;
        for(seq=0; seq < matrix->num_sequences; seq++) {
            if (matrix->block->residues[seq][pos] < 1 ||
                matrix->block->residues[seq][pos] > 22) num_residues-- ;
 	}

        if (num_residues > 0) { /* to prevent counting columns that are
                                all X's or all -'s */

          /* correction factor for small sample size */
        e = (AA-1) / ((double) 2 * ln2 * num_residues) ;

        /* =====end of calculating error correction ================*/
        assert (matrix->weights != NULL);

        r = hmax ; /* R = log 20 */
        for(aa=1, totalweight=0.; aa < AA+1; aa++) {
           totalweight +=  matrix->weights[aa][pos]  ;
        }
        /* totalweight is the sum of all aa weights at a certain pos */
        for(aa=1; aa < AA+1; aa++)
                    /* start loop at 1 and end at AA+1 because aa values
                      at the weight matrix start at 1, 0 is the gap value */
       {
	if (matrix->weights[aa][pos] > 0)
            {
            dtemp = (double) matrix->weights[aa][pos] / totalweight;
           /* in Tom Schneider's logo eq., dtemp corresponds to the frequency
                of an amino acid */
        /* dtemp is the proportion of that particular aa relative to all
         the aa's (total weight ) */
	 temp = dtemp * log (dtemp) / ln2;
            r += dtemp * log(dtemp); /* sum over all residues dtemp*log(dtemp)*/
                                    /* is H(pos) */
            } /* end of if weights > 0 */
         } /* end of for for calculating r*/

      r /= ln2 ;       /* convert to bits , R = 2 - H(pos) in bits now */
	if (error_correction_option) {
	      r -= e ;                     /* a correction for small sample sizes */
	}
	info_array[pos] = r;

        } /* end if num_residues > 0 */
   } /* end of for-loop for all positions */
	free_matrix (matrix);
	return info_array;

} /* end of information_per_residue  */

/* June 25, 2002, changed index from i=1 to i=0 so first sequence
will be assigned 100% identity */
void
percentage_identity_with_seq0_seqs (Sequence* seqs[MAXSEQ], int nseqs)
{
	int i;
	int identical, pos;
	int alignment_length;


	for (i = 0; i < nseqs; i++) {
		assert (seqs[i]->length == seqs[0]->length);
		/* this is an alignment, lengths should be equal */
		alignment_length = seqs[0]->length;
		identical = 0;
		for (pos = 0; pos < seqs[0]->length; pos++) {
			if (seqs[i]->sequence[pos] == seqs[0]->sequence[pos]) {
				identical++;
			}
			/* for partial alignments */
			if (seqs[i]->sequence[pos] == aa_atob['X']) {
				alignment_length--;
			}
		}
/*	printf("%d identical %d aligned\n", identical, alignment_length); */
	seqs[i]->undefined = (int) round ((double) identical/
						   (double)alignment_length * 100);

		printf ("%s has %d identity\n", seqs[i]->name, seqs[i]->undefined);
	}
}


/* puts names of sequences from block into seqnamehash that are >=
percentage identical */

void
percentage_identity_with_seq0 (Block* block, HashTable seqnamehash, double percentage)
{
	int identical;
	int width;
	int seq, pos;
	double decimal;

	Insert (block->sequences[0].name, seqnamehash);
	for (seq = 1; seq < block->num_sequences; seq++) {
		identical = 0;
		for (pos = 0; pos < block->width; pos++) {
			if (block->sequences[0].sequence[pos] ==
			    block->sequences[seq].sequence[pos]) {
				identical++;
			}
		} /* end of for all positions */
		decimal = (double) (identical) / (double) block->width;
/*printf ("Percentage for %s is %.1f\n", block->sequences[seq].name,
					decimal); */
		if ( decimal >= percentage) {
/*printf ("**Inserting %s\n", block->sequences[seq].name); */
			Insert(block->sequences[seq].name, seqnamehash);
		}
	}

} /* end percentage_identity_with_seq0 */

void
calculate_percentage_identity (Block* block)
{
	int identical;
	int min_percentage, percentage_identical;
	int s1, s2, pos;

	min_percentage = 100;

	for (s1 = 0; s1 < block->num_sequences; s1++) {
		for (s2 = s1 + 1 ; s2 < block->num_sequences; s2++) {
			identical = 0;
			for (pos = 0; pos < block->width; pos++) {
				if (block->sequences[s1].sequence[pos] ==
				    block->sequences[s2].sequence[pos]) {
					identical++;
				}
			} /* end of for all positions */
			percentage_identical =  (int) ((double) identical/ (double) block->width * 100);
			if (percentage_identical < min_percentage) {
				min_percentage = percentage_identical;
			}
		} /* end for s2 */
	} /* end for s1 */
	block->percentile = min_percentage;

} /* end of calculate_percentage_identity */

Block*
remove_seq0_Xes_from_block (Block* block)
{
        Block* newblock;
        Residue* residue_pointer;
        int i, pos;
	int non_X_length, newblockpos;

        non_X_length = 0;
	for (i=0; i < block->width; i++) {
		if (block->sequences[0].sequence[i] != aa_atob['X']) {
			non_X_length++;
		}
	}
/* printf ("entered remove X %d\n", non_X_length); */
	newblock = (Block*) malloc (sizeof (Block));

         /* copy header stuff */
        strcpy (newblock->id, block->id);
        strcpy (newblock->ac, block->ac);
        strcpy (newblock->de, block->de);
        sprintf (newblock->bl, "; width=%d; seqs=%d; ", non_X_length, block->num_sequences);
	strcpy (newblock->number, block->number);
        strcpy (newblock->motif, block->motif);
        newblock->width = non_X_length;  /* NEW WIDTH */
        newblock->strength = block->strength;
        newblock->max_sequences = block->max_sequences;
        newblock->num_sequences = block->num_sequences;
        newblock->max_clusters = 1; /* allocate space for a single cluster*/
        newblock->num_clusters = 1;
        newblock->min_prev = block->min_prev;
        newblock->max_prev = block->max_prev;

        newblock->clusters= (Cluster *) calloc (newblock->max_clusters,
                                                sizeof (Cluster));
        newblock->clusters[0].num_sequences = newblock->num_sequences;

        /* sequences is an array of Sequence structures */
        newblock->sequences = (Sequence *) calloc (block->max_sequences,
                                                sizeof (Sequence));
        residue_pointer = (Residue *) calloc (non_X_length * block->max_sequences,
                 sizeof (Residue));

        newblock->residues = (Residue**) calloc (block->max_sequences, sizeof
                           (Residue *));


        newblock->clusters[0].sequences = &(newblock->sequences[0]);

/* initialize residues */
        for (i=0; i < block->max_sequences; i++) {
                newblock->sequences[i].sequence = &(residue_pointer[i*newblock->width]);
                newblock->residues[i] = &(residue_pointer[i*newblock->width]);
        }

        for (i = 0; i < block->num_sequences; i++) {
                /* copy information stuff */
                strcpy (newblock->sequences[i].name,
                            block->sequences[i].name);
                strcpy (newblock->sequences[i].info,
                            block->sequences[i].info);
                newblock->sequences[i].position =
                            block->sequences[i].position;
                newblock->sequences[i].length =
                            newblock->width ;
                newblock->sequences[i].max_length=
                            block->sequences[i].max_length;
                newblock->sequences[i].type =
                            block->sequences[i].type;

                        /* copy residues */
                for (pos = 0, newblockpos = 0; pos < block->width; pos++) {
				if (block->sequences[0].sequence[pos] != aa_atob['X']){
	                               newblock->sequences[i].sequence[newblockpos] =
                                        block->sequences[i].sequence[pos];
/*                            printf ("%c", aa_btoa[newblock->sequences[i].sequence[newblockpos]] ); */
                			newblockpos++;
				}
		}
/*      printf ("\n"); */
	} /* end of for, gone through all sequences in oldblock*/
	strcat (newblock->de, "updated");
	assert (newblock->width == newblockpos);
	print_block (newblock);
	return (newblock);

} /* end of remove_seq0_Xes_from_block */

void
copy_sequence_identity_from_sequences (Block* block,
				       Sequence*seqs[MAXSEQ],
				       int nseqs)
/* copies sequence identity (in sequence->undefined) to block */
{
	int i, j;

	for (i = 1; i < block->num_sequences; i++)
	{
		j = 0;
		while (j < nseqs && strcmp (seqs[j]->name, block->sequences[i].name) != 0) {
			j++;
		}
		if ( j == nseqs) {
			fprintf (errorfp, "Could not find %s in sequences\n", block->sequences[i].name);
			exit (-1);
		}
		block->sequences[i].undefined = seqs[j]->undefined;
/*		printf ("%s has %d identity \n", block->sequences[i].name,
					  block->sequences[i].undefined); */
	}

} /* end of copy_sequence_identity_from_sequences */

void
remove_seqs_percent_identical_to_query (Sequence* seqs[MAXSEQ], int* no_of_seqs, double percent_identical)
{

	int nseqs;
	int i;
	int entered;


	nseqs = *(no_of_seqs);
	entered = 0;

	printf ("*** The following sequences have been removed because they ");
	printf (" were found to be over %d%% identical with your protein query: ", (int) percent_identical);

	fprintf (stderr, "*** The following sequences have been removed because they ");
    fprintf (stderr, " were found to be over %d%% identical with your protein query: ", (int) percent_identical);

    percentage_identity_with_seq0_seqs (seqs, nseqs);
    for (i = 1; i < nseqs; i++)  {
        if (seqs[i]->undefined >= percent_identical ) { /* >= identical, remove */
            entered = 1;
			printf (" %s,", seqs[i]->name);
			fprintf (stderr, " %s,", seqs[i]->name);

            /* below was removed to reduce verbiage on webpage */
            /* printf ("*** Sequence %s was found to be %d%% ",
                seqs[i]->name,
                seqs[i]->undefined);
            printf ("identical with the query sequence and was ");
            printf ("removed to reduce redundancy.***<BR>\n");
            fprintf (stderr, "*** Sequence %s was found to be ", seqs[i]->name);
            fprintf (stderr,"%d%% identical with the query ", seqs[i]->undefined);
            fprintf (stderr, "sequence and was removed to reduce ");
            fprintf (stderr, "redundancy.***<BR>\n"); */

            free_sequence(seqs[i]);
            seqs[i] = seqs[nseqs-1];
            seqs[nseqs-1] = NULL;

            nseqs--;
            /* bug found by Chris Saunders 08/08/01
            **  culling loop mod:
            **    decriment the index so that the new seq now
            **  addressed as seqs[i] can be checked against
            **  against the percent_identical threshold
            */
			i--;
		}
    }
	if (entered == 1) {
		printf (".\n\n"); fprintf ( stderr, ".<BR><BR>\n");
	} else {
		printf (" NONE.\n\n"); fprintf (stderr, " NONE.<BR><BR>\n");
	}

	if (*no_of_seqs != nseqs) {
/*		fprintf (stderr, "*** <B>%d</B> sequences were used to estimate probabilities.  ***<BR>\n", nseqs); */
	}
	*no_of_seqs = nseqs;
}


/* option suppress_pbweights. in order to save time, if block is already weighted,
pass in TRUE and pbweights will not be recalculated */

double
calculate_median_information_of_block (Block* block, int error_correction,
int suppress_pbweights)
{
	double* info_array;
	double median;
	int i;

	if (!suppress_pbweights) {
		pb_weights (block);
	}
	info_array =  calculate_info_for_each_pos (block,
                 error_correction);
/*	for (i=0; i < block->width; i++) {
		printf ("info %.2f\n", info_array[i]);
	} */
	median = median_of_array (info_array, block->width);
	free (info_array);
	return median;

}

double
calculate_average_information_of_block (Block* block, int error_correction,
int suppress_pbweights)
{
        double* info_array;
        double mean;
        int i;

        if (!suppress_pbweights) {
                pb_weights (block);
        }
        info_array =  calculate_info_for_each_pos (block,
                 error_correction);
/*      for (i=0; i < block->width; i++) {
                printf ("info %.2f\n", info_array[i]);
        } */
        mean = mean_of_array (info_array, block->width);
        free (info_array);
        return mean;

}
void add_sequence_to_block (Block* block, Sequence* seq, int weight_block)
{
	int index, pos;

	if (block->num_sequences >= block->max_sequences) {
		printf ("add sequence to block resizeing block\n");
		resize_block_sequences (block);
	}

  /* add in the new sequence */

	index = block->num_sequences;

        /* insert sequence as last item in the array */
        /* copy information stuff */
       strcpy (block->sequences[index].name,
                            seq->name);
       strcpy (block->sequences[index].info,
                            seq->info);
       block->sequences[index].position =
                            seq->position;
       block->sequences[index].length =
                            seq->length;
       block->sequences[index].max_length=
                            seq->max_length;
       block->sequences[index].type =
                            seq->type;
                        /* copy residues */
       for (pos = 0; pos < block->width; pos++) {
       /*  printf ("%c\n", aa_btoa[originalblock->sequences[seqindex].sequence[pos]]); */
          block->sequences[index].sequence[pos]=
                              seq->sequence[pos];
       }
        block->num_sequences++;
        if (block->num_clusters != 1) {
                fprintf (errorfp, "there is more than one cluster\n");
                fprintf (errorfp, "add seq subroutine will fail\n");
                exit (-1);
        }
        block->clusters[0].num_sequences++;
        if (weight_block) {
		pb_weights(block);
	}
} /* add_seq_to_block*/

void
print_block_sequences (Block* block, FILE* outfp)
{

	int i;

	for (i = 0; i < block->num_sequences; i++) {
		output_sequence (block->sequences[i], outfp);
	}

} /* end print_block_sequences */

void
print_block_ids (Block* block, FILE* fp, int print_first_id)
{
	int i;

	if (print_first_id) {
		fprintf (fp, "%s\n", block->sequences[0].name);
	}

	for (i = 1; i < block->num_sequences; i++) {
		fprintf (fp, "%s\n", block->sequences[i].name);
	}

} /* end print_block_ids */

void copy_weights (Block* block, Block* oldblock)

{
	int i;

      for (i = 0; i < block->num_sequences; i++) {
                block->sequences[i].weight =
                            oldblock->sequences[i].weight;
	}
}

double
calculate_median_info_for_pos_in_block (Block* block, int pos)
{
	Block* block_with_seqs_at_pos;
	double * info_array;
	double median;

        block_with_seqs_at_pos = subblock_of_seqs_with_aa_at_pos
                                                        (block, pos);
        info_array =  calculate_info_for_each_pos
                                (block_with_seqs_at_pos, FALSE);
        median = median_of_array (info_array, block->width);

	free_block (block_with_seqs_at_pos);
	free (info_array);

	return median;

}


#endif
