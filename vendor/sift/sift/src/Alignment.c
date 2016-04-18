/* (C) Copyright 2000, Fred Hutchinson Cancer Research Center */
/* Use, modification or distribution of these programs is subject to */
/* the terms of the non-commercial licensing agreement in license.h. */

/* Alignment.c
contains subroutines in dealing with alignment conversions
*/
/* 12/20/00 Bob found an error in clump.  Query sequence (first sequence
in seqs[MAXSEQ] is not necessarily the first sequence in the cluster
when there aren't any other sequences clustered with it.  Fixed so that
query sequence is moved to be the first cluster */
/* modified psiblast_pairwise to read in evalue */

/* #include "blocksprogs.h" has already been defined */

#ifndef _ALIGNMENT_C_
#define _ALIGNMENT_C_

#include <assert.h>

#include "PN_convert.h"
#include "Alignment.h"

struct cluster_pair {
        double score;
        int cluster;
};
#define INDEX(n, col, row) (col*n - (col*(col+3))/2 - 1 + row)

/*========================================================================
 Try clustal format, which has multiple lines per sequence with the
 sequence name in the first 16 columns of each line.
 Writes each set of sequence segments out to a temporary file, then
 reads them back in & appends them to sequence array.
 Assumes each set of segments is separarated by at least one blank line.
 Assumes each line of a set contains the sequence name, whitespace, residues:

CLUSTAL W(1.60) multiple sequence alignment

JC2395          NVSDVNLNK---YIWRTAEKMK---ICDAKKFARQHKIPESKIDEIEHNSPQDAAE----
KPEL_DROME      MAIRLLPLPVRAQLCAHLDAL-----DVWQQLATAVKLYPDQVEQISSQKQRGRS-----
FASA_MOUSE      NASNLSLSK---YIPRIAEDMT---IQEAKKFARENNIKEGKIDEIMHDSIQDTAE----

JC2395          -------------------------QKIQLLQCWYQSHGKT--GACQALIQGLRKANRCD
KPEL_DROME      -------------------------ASNEFLNIWGGQYN----HTVQTLFALFKKLKLHN
FASA_MOUSE      -------------------------QKVQLLLCWYQSHGKS--DAYQDLIKGLKKAECRR

JC2395          IAEEIQAM
KPEL_DROME      AMRLIKDY
FASA_MOUSE      TLDKFQDM
============================================================================*/
int try_clustal(FILE* ifp, Sequence* seqs[MAXSEQ], char desc[SMALL_BUFF_LENGTH])
{
   FILE *tmp;
   Sequence *tmp_seq;
   int nseq, this_seq, i, my_pid, db_type, seq_type;
   char line[MAXLEN], tmp_name[30], *name, *residues;
   Residue c;

   nseq = 0; seqs[0] = NULL;

   /*  Look for the word clustal on the first non-blank line, if find it,
        skip to next non-blank line and begin  */
   while (!feof(ifp) && fgets(line, MAXLEN, ifp) != NULL &&
          strlen(line) < (int) 4 )
                ;
   if (strstr(line, "CLUSTAL") == NULL && strstr(line, "Clustal") == NULL)
   {
      rewind(ifp);
      return(nseq);
   }
   else
   {
      fprintf(stderr, "Trying CLUSTAL format.\n");
      fgets(line, MAXLEN, ifp);
   }

   /*  Read sequences until the next blank line, or until a line
       with blanks at the start
        then start over and add to existing sequences, etc.  */
   sprintf(tmp_name, "mablock.%s", desc);
   tmp=fopen(tmp_name, "w");
   while (!feof(ifp))
   {
      /*   Write a set of sequence segments to a temp file */
      /*  NOTE: Blank lines from netscape can be "^M\n\0", len=3 */
      while (!feof(ifp) && fgets(line, MAXLEN, ifp) != NULL &&
           strlen(line) > (int) 3 )
      {
         if (line[0] != ' ')        /* some lines of asterisks, etc */
         {
            name = strtok(line, " \t");
            residues = strtok(NULL, " \t\r\n");
/*printf("%s %s\n", name, residues); */
            if (tmp != NULL)
            {
               if (name != NULL && residues != NULL)
               {
                   fprintf(tmp, ">%s\n", name);
                   fprintf(tmp, "%s\n", residues);
               }
               else
               {   fprintf(stderr,"ERROR reading Clustal file\n%s", line);  }
            }
         }
      } /* end while */
      /*   Encountered a blank line or eof */
      fclose(tmp);
      if ( (tmp = fopen(tmp_name, "r")) != NULL)
      {
         db_type = type_dbs (tmp, DbInfo);   seq_type = AA_SEQ;
         if (seqs[0] == NULL)
         {
         /*   initialize sequences */
            nseq = 0;
            if (db_type >= 0)
            {
               while ( nseq < MAXSEQ &&
                   (seqs[nseq] = read_a_sequence(tmp, db_type, seq_type)) != NULL)
               {
                  nseq++;
               }
            }
         }
         else
         /*     add on to existing sequences */
         {
            this_seq = 0;
            if (db_type >= 0)
            {
               while ( this_seq < MAXSEQ &&
                   (tmp_seq = read_a_sequence(tmp, db_type, seq_type)) != NULL)
               {
                  if (seqs[this_seq] != NULL &&
                      strcmp(seqs[this_seq]->name, tmp_seq->name) == 0 )
                  {
                     if ((seqs[this_seq]->length + tmp_seq->length) >
                          seqs[this_seq]->max_length)
                              resize_sequence(seqs[this_seq]);
                     for (i=0; i<tmp_seq->length; i++)
                     {
                        c = tmp_seq->sequence[i];
			seqs[this_seq]->sequence[ seqs[this_seq]->length++ ] =
                          c;
                     }
                  }
                  this_seq++;
                  free_sequence(tmp_seq);
               }
            }
         }
         fclose(tmp);
      }
      tmp=fopen(tmp_name, "w");
   }

   sprintf(line, "\\rm %s", tmp_name);
   system(line);

   return(nseq);
}  /*  end of try_clustal   */
/*============================================================================
MSF format:  Comments until alignment begins after "//" line:

//

            1                                                   50
P09254-1    .....KRQED AGYDICVPYN LYLKR..... NEFIKIVLPI IRDWDLQHPS
A37470-1    .TFAPKRDED AGYDIAMPYT AVL....... APGENLHVRL PVAYAADAHA
Q00030-2    DYFAPKRDED AGYDISAQTN ATI....... EPDESYFVEL PIVFSSSNPA
P28892-1    .....KRVED AGYDISAPED ATI....... DPDESHFVDL PIVFANSNPA
P10234-1    .....KREED AGFDIVVRRP VTV.P..... ANGTTVVQPS LRMLHADAGP
============================================================================*/
int try_msf(FILE* ifp, Sequence *seqs[MAXSEQ], char desc[SMALL_BUFF_LENGTH])
{
   FILE *tmp;
   Sequence *tmp_seq;
   int nseq, this_seq, i, my_pid, db_type, seq_type;
   char line[MAXLEN], tmp_name[30], *name, *residues;

   nseq = 0; seqs[0] = NULL;
   fprintf(stderr, "Trying  MSF format.\n");

   /*  Look for "//"     */
   while (!feof(ifp) && fgets(line, MAXLEN, ifp) != NULL &&
           strncmp(line, "//", 2) != 0)
                ;
   /*  Didn't find "//"  */
   if (feof(ifp)) { rewind(ifp); }

   /*  Read sequences until the next blank line, or until a line
       with blanks at the start
        then start over and add to existing sequences, etc.  */
   sprintf(tmp_name, "mablock.%s", desc);
   tmp=fopen(tmp_name, "w");
   while (!feof(ifp))
   {
      /*   Write a set of sequence segments to a temp file */
      /*  NOTE: Blank lines from netscape can be "^M\n\0", len=3 */
      while (!feof(ifp) && fgets(line, MAXLEN, ifp) != NULL &&
           strlen(line) > (int) 3 )
      {
         if (line[0] != ' ')        /* some lines of asterisks, etc */
         {
            name = strtok(line, " \t");
            residues = strtok(NULL, "\t\r\n");
/*printf("%s %s\n", name, residues);*/
            if (tmp != NULL)
            {
               if (name != NULL && residues != NULL)
               {
                   fprintf(tmp, ">%s\n", name);
                   fprintf(tmp, "%s\n", residues);
               }
               else
               {   fprintf(stderr,"ERROR reading MSF file\n%s", line);  }
            }
         }
      } /* end while */
      /*   Encountered a blank line or eof */
      fclose(tmp);
      if ( (tmp = fopen(tmp_name, "r")) != NULL)
      {
         db_type = type_dbs(tmp, DbInfo);   seq_type = AA_SEQ;

         if (seqs[0] == NULL)
         {
         /*   initialize sequences */
            nseq = 0;
            if (db_type >= 0)
            {
               while ( nseq < MAXSEQ &&
                   (seqs[nseq] = read_a_sequence(tmp, db_type, seq_type)) != NULL)
               {
                  nseq++;
               }
            }
         }
         else
         /*     add on to existing sequences */
         {
            this_seq = 0;
            if (db_type >= 0)
            {
               while ( this_seq < MAXSEQ &&
                   (tmp_seq = read_a_sequence(tmp, db_type, seq_type)) != NULL)
               {
                  if (seqs[this_seq] != NULL &&
                      strcmp(seqs[this_seq]->name, tmp_seq->name) == 0 )
                  {
                     if ((seqs[this_seq]->length + tmp_seq->length) >
                          seqs[this_seq]->max_length)
                              resize_sequence(seqs[this_seq]);
                     for (i=0; i<tmp_seq->length; i++)
                     {
			seqs[this_seq]->sequence[ seqs[this_seq]->length++ ] =
                          tmp_seq->sequence[i];
	             }
                  }
                  this_seq++;
                  free_sequence(tmp_seq);
               }
            }
         }
         fclose(tmp);
      }
      tmp=fopen(tmp_name, "w");
   }

   sprintf(line, "\\rm %s", tmp_name);
   system(line);

   return(nseq);
}  /*  end of try_msf   */

void
change_periods_to_dashes (int nseqs, Sequence* seqs[MAXSEQ])
{
	int pos;
	int seq;
	int length;
	length = seqs[0]->length;
	for (pos = 0; pos < length; pos++) {
		for (seq = 0; seq < nseqs; seq++) {
		    if (aa_btoa[seqs[seq]->sequence[pos]] == '*') {
			seqs[seq]->sequence[pos] = aa_atob['-'] ;
		    }
		}
	}
} /* end of change_periods_to_dashes */


void
change_Xes_to_dashes (int nseqs, Sequence* seqs[MAXSEQ])
{
        int pos;
        int seq;
        int length;
        length = seqs[0]->length;
        for (pos = 0; pos < length; pos++) {
                for (seq = 0; seq < nseqs; seq++) {
                    if (aa_btoa[seqs[seq]->sequence[pos]] == 'X') {
                        seqs[seq]->sequence[pos] = aa_atob['-'] ;
                    }
                }
        }
} /*end of change_Xes_to_dashes */
/* ===================================================================
remove_gaps
Takes multiple alignment in seqs and removes columns which correspond
to gaps "-" in QUERY sequence.  QUERY is checked to be in 1rst entry
=====================================================================*/

void
remove_gaps (Sequence* newseqs[MAXSEQ], Sequence* seqs[MAXSEQ], int nseqs, int query_length)
{
	int pos, i;
	int newpos;
	Residue c;

	assert (strstr (seqs[0]->name, "QUERY") != 0);

	initialize_seqs (newseqs, nseqs, query_length);

	for (i = 0; i < nseqs; i++) {
		strcpy (newseqs[i]->name, seqs[i]->name);
		strcpy (newseqs[i]->info, " ");
	}
	newpos = 0;
	for (pos = 0; pos < seqs[0]->length; pos++) {
		if (seqs[0]->sequence[pos] != aa_atob['-']) {
			for (i = 0; i < nseqs; i++) {
				c = seqs[i]->sequence[pos];
				newseqs[i]->sequence[newpos] =c;
				newseqs[i]->length++;
			}
			assert (newpos < query_length);
			newpos++;
		}
	} /* end for finished looking through gapped alignment */
/*	assert (newpos == query_length); */
/* Sept. 10th, removed assert because alignment may be shorter than entire
query sequence length */
}

/* =======================================================================
initialize_seqs (Sequence* seqs[MAXSEQ], int nseq, int query_length)

initializes the array of seqs, allocating space for nseq # of sequences
with query_length # of residues

=========================================================================*/
void
initialize_seqs (Sequence* seqs[MAXSEQ], int nseq, int query_length)
{
	int i;

	for (i = 0; i < nseq; i++) {
		seqs[i] = (Sequence *) malloc (sizeof (Sequence));
		seqs[i]->position = 0;
		seqs[i]->sequence = (Residue *) calloc (query_length,
						sizeof (Residue));
		seqs[i]->length = 0;
		seqs[i]->type = AA_SEQ;
	}
}

/*==========================================================================
free_seqs
deallocate memory as initialized in initialize_seqs
========================================================================*/
void
free_seqs (Sequence* seqs[MAXSEQ], int nseq)
{
	int i;

	for ( i = 0; i < nseq; i++) {
		if (seqs[i]->sequence != NULL) {
			free (seqs[i]->sequence);
		}
		if (seqs[i] != NULL) {
			free(seqs[i]);
		}
		seqs[i] = NULL;
	}
}

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
int flat_master_slave (FILE* fp, Sequence *seqs[MAXSEQ], char desc[SMALL_BUFF_LENGTH])
{
	FILE* tmp;
	Sequence *tmp_seq;
	int nseq, this_seq, i, my_pid, db_type, seq_type;
	char line[MAXLEN], tmp_name[30], *name, *residues;
	char* startpos, *endpos; /* string pointers to this numbers, */
				/* i'm not using them, so not bothering to*/
				/* convert to int */
	char* prev_name, *prev_residues;
	char* ptr;
	int pos;
	int done;

	nseq = 0; seqs[0] = NULL;
	pos = -1;

	done = false;

	prev_name = name = prev_residues = residues = NULL;
   sprintf(tmp_name, "%s.tmp", desc);
   tmp=fopen(tmp_name, "w");
   while (!feof(fp) && !done)
   {
     /*   Write a set of sequence segments to a temp file */
      while (!feof(fp) && fgets(line, MAXLEN, fp) != NULL &&
           line[0] != '\n' && strstr (line, "Lambda") == NULL
	   && strstr (line, "Database:") == NULL )
	{
	if (line[0] != ' ')        /* some lines of asterisks, etc */
         {
   	    name = strtok(line, " \t");
            startpos = strtok (NULL, " \t");
	    if (strstr (startpos, "------") != NULL) /*there are no
						residues aligned for this
						particular sequence */
	    {
		residues = startpos;
	    	ptr = strpbrk (startpos, "\n");
	     	*ptr = '\0';
	    } else {
		 if ( (strcmp (name, "QUERY") == 0) && pos == -1) {
			pos = atoi (startpos);
	         }
		residues = strtok(NULL, " \t\r\n");
	    } /* end if strstr ------------ */
	    if (prev_name != NULL && strcmp (prev_name, name) == 0) {
			printf ("MERGING*************%s  %s\n", prev_name, name);
			merge_res (residues, prev_residues);
			free (prev_name); free (prev_residues);
			prev_name = NULL; prev_residues = NULL;
	    }
            if (tmp != NULL)
            {
	       if (name != NULL && residues != NULL &&
		   prev_name != NULL && (strcmp (prev_name, name) != 0) )
               {
                   fprintf(tmp, ">%s\n", prev_name);
                   fprintf(tmp, "%s\n", prev_residues);
		   free (prev_name); free (prev_residues);
		   prev_name = NULL; prev_residues= NULL;
	       }
               else
               {
		   fprintf(stderr,"Please contact us: ERROR reading psiblast alignment\n%s", line);
		}
            } /* end of if (tmp != NULL) */
         /* update */
	 if (prev_name == NULL && prev_residues == NULL) {
	 	prev_name = copy_string (name);
	 	prev_residues = copy_string (residues);
         } else {
		fprintf(errorfp, "ERROR! Why haven't %s and %s been freed?\n", prev_name,
		prev_residues);
		fprintf (errorfp, "Current %s %s\n", name, residues);
		exit (-1);
	 }
	} /* end of if (line[0] != ' ') */
      } /* end while */

	if (strstr (line, "Database:") != NULL) {
		done = true;
	}
	/* print out last sequences */
      fprintf (tmp, ">%s\n", prev_name);
      fprintf (tmp, "%s\n", prev_residues);
      free (prev_name); free (prev_residues);
      prev_name = NULL; prev_residues = NULL;

      /*   Encountered a blank line or eof */
      fclose(tmp);
      if ( (tmp = fopen(tmp_name, "r")) != NULL)
      {
         db_type = type_dbs(tmp, DbInfo);
         seq_type = AA_SEQ;
         if (seqs[0] == NULL)
         {
         /*   initialize sequences */
            nseq = 0;
            if (db_type >= 0)
            {
               while ( nseq < MAXSEQ &&
                   (seqs[nseq] = read_a_sequence(tmp, db_type, seq_type)) != NULL )
               {
                  nseq++;
               }
            }
         } /* end of if(seqs[0] == NULL) */
	/*	sequences have already been initialized*/
         else
         /*     add on to existing sequences */
         {
            this_seq = 0;
            if (db_type >= 0)
            {
               while ( this_seq < MAXSEQ &&
                   (tmp_seq = read_a_sequence(tmp, db_type, seq_type)) != NULL)
               {
                  if (seqs[this_seq] != NULL &&
                      strcmp(seqs[this_seq]->name, tmp_seq->name) == 0 )
                  {
                     if ((seqs[this_seq]->length + tmp_seq->length) >
                          seqs[this_seq]->max_length)
                              resize_sequence(seqs[this_seq]);
                     for (i=0; i<tmp_seq->length; i++)
                     {
                        seqs[this_seq]->sequence[ seqs[this_seq]->length++ ] =
                          tmp_seq->sequence[i];
                     }
                  }
                  this_seq++;
                  free_sequence(tmp_seq);
               }
            }
         }
         fclose(tmp);
      } /* end of line[0] != '\n' && line != "Lambda" */
      seqs[0]->position = pos;
	tmp = fopen (tmp_name, "w");
      if (strstr (line, "Lambda") != NULL) {
	sprintf (line, "\\rm %s", tmp_name);
	/*system (line); */
	return (nseq);
      }
   } /* end of while feof (fp) */

   sprintf(line, "\\rm %s", tmp_name);
/*   system(line); */

   return(nseq);
} /* end of flatmaster slave */

/* copy string : allocates space in dest and copies character by
character from src so that the contents are the same, but addresses
 are different*/
char*
copy_string (char* src)
{
	int length, i;
	char* dest;

	length = strlen (src);
	dest = (char*) calloc (length + 1, sizeof (char));
	for (i = 0; i < length; i++) {
		*(dest + i)  =  src[i];
	}
	*(dest+ length) = '\0';

	return dest;

} /* copy_string */

/* subroutine to merge local alignments from psiblast output
If there is an aa in one pos and a - in the corresponding pos in the
other sequence, choses the aa over the -.  There should not be an aa
in both seq1 and seq2 at the same position, this will return an error.
Change seq1.
*/
void
merge_res (char* seq1, char* seq2)
{
	int length1, length2, i;
	char c;

	length1 = strlen (seq1);
	length2 = strlen (seq2);
	assert (length1 = length2);

	for (i = 0; i < length1; i++) {
		c = *(seq2 + i);
		if (c != '-') {
			if (*(seq1 + i) != '-') {
				printf ("%d %c %c, ", i, *(seq1+i), *(seq2+i) );			}
			*(seq1 + i) = c;
		}
	}

} /* end of merge_res */

/*========================================================================
makes a block from a region of sequences. block starts from firstpos
(the position in array, not sequence) and the width is passed in.
if seq_weights_option == TRUE, use Block's weights, otherwise if == FALSE,
calculate
========================================================================*/
Block *make_block(int width, int firstpos, int nseq, Sequence* seqs[MAXSEQ],
		  int use_seq_weights_option)
{
   Block *block;
   int bpos, spos, s;

   block = new_block(width, nseq);
   block->min_prev = 0; block->max_prev = 0;
   block->strength = 0;
   strcpy (block->id, "UNK_ID");
   strcpy (block->ac, "UNK_AC");
   strcpy (block->de, "UNK_DE");
   strcpy (block->bl, "UNK_BL");
   strcpy (block->number, "UNK_BLNUM");
   strcpy (block->motif, "UNK_MOTIF");

   for (s=0; s < nseq; s++)
   {
      strcpy(block->sequences[s].name, seqs[s]->name);
      block->sequences[s].position = firstpos + 1; /* position in sequence, */
                                                /* not in array */
      block->sequences[s].weight = 100.0;
      bpos = 0;
      for (spos = firstpos; spos < firstpos + width; spos++)
          block->sequences[s].sequence[bpos++] = seqs[s]->sequence[spos];
   } /* end of sequence s */

   if (use_seq_weights_option == TRUE) {
	for (s = 0; s < nseq; s++) {
		block->sequences[s].weight = seqs[s]->weight;
	}
   } else if (use_seq_weights_option == FALSE) {
   	pb_weights(block);
   }
   return(block);
}  /*  end of make_block */

Block* make_block_from_positions (int width, int* positions_with_real_aa,
				  Sequence* seqs[MAXSEQ],
				  int nseqs)
{
   Block *block;
   int bpos, spos, s;
   int numseq;
   int blockseq;

   numseq = 0;
   for (s = 0; s < nseqs; s++) {
	if (seqs[s]->undefined) {
		numseq++;
	}
   }

   block = new_block(width, numseq);
   block->min_prev = 0; block->max_prev = 0;
   block->strength = 0;
   strcpy (block->id, "UNK_ID");
   strcpy (block->ac, "UNK_AC");
   strcpy (block->de, "UNK_DE");
   strcpy (block->bl, "UNK_BL");
   strcpy (block->number, "UNK_BLNUM");
   strcpy (block->motif, "UNK_MOTIF");

   printf("entering for loop\n");
   blockseq = 0;
   for (s=0; s < nseqs; s++)
   {
      if (seqs[s]->undefined) {
	strcpy(block->sequences[blockseq].name, seqs[s]->name);
      	block->sequences[blockseq].position = 0;
        block->sequences[blockseq].weight = 100.0;
        bpos = 0;
        for (spos = 0; spos < seqs[s]->length; spos++) {
          if (positions_with_real_aa[spos]) {
		block->sequences[blockseq].sequence[bpos++] =
					seqs[s]->sequence[spos];
	  } /* end of if */
         } /* end for sequence */
        printf ("bpos %d withd %d\n", bpos, width);
	assert (bpos == width);
         blockseq++;
       } /* end if sequence defined */
   } /* end of for loop  s */

   assert (blockseq == numseq );

   pb_weights(block);
   return(block);
}  /*  end of make_block */


/* makes a block out of a seqs array.  gaps kept */
Block*
block_from_seq_given_pos (Sequence* seqs[MAXSEQ],
			int nseq, int beg, int end)
{
	int width;
	Block* block;

	assert (strcmp (seqs[0]->name, "QUERY") == 0);
	width = width_with_gaps (seqs[0], beg, end);
	block = make_block (width, beg, nseq, seqs, FALSE);
	return block;

}

int
width_with_gaps (Sequence* seq, int beg,int end)
{
	int longer_width; /* width with gaps */
	int pos;
	int aa_width; /* width with just amino acids */

	aa_width = end - beg + 1;

	longer_width = 0;
	pos = beg;
	while (aa_width > 0) {
		if (seq->sequence[pos] != aa_atob['-'] &&
			seq->sequence[pos] != aa_atob['.']) {
			aa_width--;
		}
		pos++;
		longer_width++;
	}
	return longer_width;

} /* end of width_with_gaps */

/* 3/2/01 aded option return_error.  TRUE -- exit program, FALSE, just return */

/* read_psiblast_header_until_first */
void
read_psiblast_header_until_first (FILE* fp)
{
        char line[MAXLEN];

	while (!feof (fp) && fgets (line, LINE_LEN, fp) != NULL &&
                strstr (line, "Sequences producing significant alignments") == NULL)
        {
                /* keep reading in lines */
        }
	if (feof (fp)) {
		fclose (fp);
			fprintf (errorfp, "ERROR! PSI-BLAST found no hits. Program terminating.\n");
                exit (-1);
	}

	fgets (line, LINE_LEN, fp); /* read in newline */

} /* end of read_psiblast_header_until_first */

/* for trent 3/23/01 */
void
read_psiblast_header_until_first_alignment (FILE* fp)
{
	char line[LARGE_BUFF_LENGTH];

	read_psiblast_header_until_first (fp);

	while (!feof (fp) && fgets (line, LINE_LEN, fp) != NULL &&
                strstr (line, "Alignments") == NULL ) {
		/* keep reading in lines */
	}
	fgets (line, LINE_LEN, fp); /* read in newline */
	fgets (line,LARGE_BUFF_LENGTH, fp);
	/*should be start of first alignment */
} /* end of read_psiblast_header_until_first_alingment */


int
read_psiblast_header_until_first_no_error (FILE* fp,int return_error )
{
        char line[MAXLEN];

        while (!feof (fp) && fgets (line, LINE_LEN, fp) != NULL &&
                strstr (line, "Sequences producing significant alignments") == NULL)
        {
                /* keep reading in lines */
        }
        if (feof (fp)) {
                fclose (fp);
                if (return_error) {
                        fprintf (errorfp, "ERROR! PSI-BLAST found no hits. Program terminating.\n");
                exit (-1);
                } else {
                        return (-1);
                }
        }

        fgets (line, LINE_LEN, fp); /* read in newline */
        return 0;

} /* end of read_psiblast_header_until_first_no_error*/

/*======================================================================
read_psiblast_header_until_last:
Processes lines until reach the beginning of last alignment.
Returns true if alignment has converged before then, false otherwise

=====================================================================*/

int
read_psiblast_header_until_last(FILE* fp, int max_iterations)
{
	char line[MAXLEN];
	int iteration;
	fpos_t filepos;
printf ("entered read_psiblastuntillat\n");
   if (max_iterations != 1) {
	while (!feof (fp) && fgets (line, LINE_LEN, fp) != NULL &&
		strstr (line, "Results from round") == NULL)
	{
		/* keep reading in lines */
	}
printf ("finished reading psiblast file\n");
	if (feof (fp)) {
		fclose (fp);
		fprintf (errorfp, "ERROR! PSI-BLAST found no hits. Program terminating.\n");
		exit (-1);
	}
	sscanf (line, "Results from round %d", &iteration);
	if (iteration == max_iterations) {
		/* get to the alignment */
		while (!feof (fp) && fgets (line, LINE_LEN, fp) != NULL &&
		strstr (line, "Sequences not found previously") == NULL)
		{
			/* keep reading in lines */
		}
		fgets (line, LINE_LEN, fp); /* this should be a newline */
		fgetpos (fp, &filepos);
		fgets (line, LINE_LEN, fp);
		if (line[0] == '\n') {
			/* converged */
			fgets (line, LINE_LEN, fp);
			/* Sept 3 2010 Pauline, no longer says converged
			   commented out */
			/* assert (strstr (line, "CONVERGED!") != NULL); */
			return TRUE;
		} else if (line[0] == '>') {
		/* 10-11-00 when there are more sequences found but not
		printed out because exceeded MAXSEQ already */
			fsetpos (fp, &filepos);
			return FALSE;
		} else { /* didn't converge */
			while (!feof (fp) &&
				 fgets (line, LINE_LEN, fp) != NULL &&
				 line[0] != '\n')
			{
			/* keep reading in lines until at a newline */
			}
			/* next line should be QUERY 1 */
			return FALSE; /* didn't converge, return FALSE */
		} /* end of else didn't converge */
	} else { /* maybe at this iteration, sequence converged */
		if (iteration != 1) {
			while (!feof (fp) &&
				fgets (line, LINE_LEN, fp) != NULL &&
                		strstr (line, "Sequences not found previously")
								 == NULL)
	                {
       		                 /* keep reading in lines */
       		         }
			fgets (line, LINE_LEN, fp); /*this should be a newline*/
			fgets (line, LINE_LEN, fp);
			if (line[0] == '\n') {
				fgets (line, LINE_LEN, fp);
				/*assert (strstr (line, "CONVERGED!") != NULL); */
				/* this is not true for psiblast ,it's just an empty line Pauline 08-23-2010*/
	printf ("converged and returning %s\n", line);
			return TRUE; /* converged */
			} else {
				return read_psiblast_header_until_last (fp,
							max_iterations);
			}
		} else { /* iteration == 1, no chance of early convergence */
			return read_psiblast_header_until_last (fp, max_iterations);
		}
	} /* end of else early convergence */
  } else {   /* end of if max_iterations != 1 */
	read_psiblast_header_until_first (fp);
	printf ("at this point lin is %s\n", line);
	while (	!feof (fp) &&	fgets (line, LINE_LEN, fp) != NULL &&
                                 line[0] != '\n')
       {
        printf ("skipping line %s\n", line);
	       /* keep reading in lines until at a newline */
       }
	/* next line should be start of first alignment */
	printf ("line is %s\n", line);
	return TRUE;
   } /* end else max_iterations == 1 */

} /* end of read_psiblast_header_until_last */

int
get_length (Sequence* seq)
{
	int i ;
	int length;

	length = 0;

	for (i = 0; i < seq->length; i++) {
		if (aa_btoa [seq->sequence[i]] != '-') {
			length++;
		}
	}
	return length;
}

/*======================================================================
   Fix up the sequence names

Chris saunders modifications for NCBI
=======================================================================*/
void fix_names_practice (int nseq, Sequence* seqs[MAXSEQ])
{
  Boolean done;
  int i, s, len, mlen, startindex, index;
  char ctemp[80];
  char* ptr;
  char* strptr;

  for (s=0; s < nseq; s++)
  {
    /*   Get rid of leading "P1;"    */
    if (strncmp(seqs[s]->name, "P1;", 3) == 0)
    {
        if ((int) strlen(seqs[s]->name) > 15) seqs[s]->name[15] = '\0';
        strcpy(ctemp, seqs[s]->name);
        strcpy(seqs[s]->name, ctemp+3);
    }
   /* Sept. 10, 2000 get gi or gb */
    strptr = strstr (seqs[s]->name, "gi|");
    if (strptr != NULL) {
        /* this is a NCBI database -- start modifications here 3/16/2002 */
	}

    /*  Check for duplicate sequence names here */
    i = s-1; done = FALSE;
    while (!done && i >= 0)
    {
       if (strcmp(seqs[s]->name, seqs[i]->name) == 0)
       {
          printf("Non-unique sequence name: %s\n", seqs[s]->name);
/*
          printf("First 18 characters of sequence name must be unique.\n");
*/
        printf ("nonunique problem\n");
  strcpy(ctemp, seqs[s]->name);
          len = strlen(seqs[s]->name);
          mlen = 9 - (int) log10((double) s);   /* #digits in s*/
          ctemp[len - mlen] = '\0'; /* Pauline's modification, want to keep
                                        name the same length (for PHYLIP) */
          sprintf(seqs[s]->name, "%s%d", ctemp, s);
          seqs[s]->name[SNAMELEN] = '\0';
          printf("Modified name to %s\n", seqs[s]->name);
          done = TRUE;
      } /* end of non-unique name */
      i--;
    }
  }  /* end of for s */
}  /*  end of fix_names */




/*======================================================================
   Fix up the sequence names
=======================================================================*/
void fix_names(int nseq, Sequence* seqs[MAXSEQ])
{
  Boolean done;
  int i, s, len, mlen, startindex, index;
  char ctemp[80];
  char* ptr;
  char* strptr;

  for (s=0; s < nseq; s++)
  {
    /*   Get rid of leading "P1;"    */
    if (strncmp(seqs[s]->name, "P1;", 3) == 0)
    {
        if ((int) strlen(seqs[s]->name) > 15) seqs[s]->name[15] = '\0';
        strcpy(ctemp, seqs[s]->name);
        strcpy(seqs[s]->name, ctemp+3);
    }
   /* Sept. 10, 2000 get gi or gb */
    strptr = strstr (seqs[s]->name, "gi|");
    if (strptr != NULL) {
	strncpy (ctemp, strptr +3 , 8);
	ctemp[8] = '\0';
/*	printf ("newname from %s to %s\n", seqs[s]->name, ctemp); */
	strcpy (seqs[s]->name, "gi");
        strcat (seqs[s]->name, ctemp);
    } /* else {
	strptr = strstr (seqs[s]->name, "gb|");
	if (strptr != NULL) {
		strncpy (ctemp, strptr + 3 , 8);
		ctemp[8] = '\0';
		strcpy (seqs[s]->name, "gb");
		strcat (seqs[s]->name, ctemp);
	}
    }
*/

    strptr = strstr (seqs[s]->name, "gnl|"); /* gb and gnl sequences, take the last
					   15 characters */
    if (strptr == NULL) { strptr = strstr (seqs[s]->name, "gb|");}
    if (strptr == NULL) {strptr = strstr (seqs[s]->name, "emb|"); }
    if (strptr != NULL) {
	len = (int) strlen(seqs[s]->name);
/*	printf ("currentname : %s ", seqs[s]->name); */
	startindex = len - 18;
	if (startindex < 0) { startindex = 0; }
	strcpy (ctemp, seqs[s]->name + startindex);
	ctemp[SNAMELEN] = '\0';
	strcpy (seqs[s]->name, ctemp);
     /* printf (" last 18 characters %s name\n", seqs[s]->name); */
    }

/* Pauline added Aug 23,2010 for psiblast processing 2nd word*/
/*    strptr = strtok (seqs[s]->name, " \t\r\n");
    if (strptr != NULL) {
         strptr = strtok ((char*) NULL, "\t\n\0");
         strcpy (ctemp, &strptr[0]);
	ctemp[SNAMELEN]='\0';
	strcpy (seqs[s]->name, ctemp);

    }
*/

/* Pauline's code : this is because psiblast can not parse incomplete
names like PIP5_RAT|P and by shortening it, was looking for field P...
which wasn't there, causing error. here I just get rid  of the |P
sept.15, 2000 changed so that | -> _ , rather than \0*/
/*    strcpy (ctemp, seqs[s]->name); index = 0; ;
printf ("name is 3423 %s\n", seqs[s]->name);
    for (i = 0; i < (int) strlen (ctemp); i++) {
       printf ("to behere %d\n", (int) strlen (ctemp) );
	 if (ctemp[i] != '|') {
       printf ("huh %s\n", ctemp[i]);
         seqs[s]->name[index] = ctemp[i];
                index++;
        }
    }
*/
    if ((int) strlen(seqs[s]->name) > SNAMELEN)
        seqs[s]->name[SNAMELEN] = '\0'; /* ensure ony 18 characters long */

ptr = strstr (seqs[s]->name, "|");
if (ptr != NULL) {
	*ptr = '\0';
}
/* commented out 12/11/00 , names were too long
this part was used for genome sequence data */
/*    ptr = strstr (seqs[s]->name, "|");
    while (ptr != NULL) {
	*ptr = '_';
        ptr = strstr (seqs[s]->name, "|");
    }
*/
/* end Pauline's code */
    /*  Check for duplicate sequence names here */
    i = s-1; done = FALSE;
    while (!done && i >= 0)
    {
       if (strcmp(seqs[s]->name, seqs[i]->name) == 0)
       {
          printf("Non-unique sequence name: %s\n", seqs[s]->name);
/*
          printf("First 18 characters of sequence name must be unique.\n");
*/
        printf ("nonunique problem\n");
  strcpy(ctemp, seqs[s]->name);
          len = strlen(seqs[s]->name);
          mlen = 9 - (int) log10((double) s);   /* #digits in s*/
          ctemp[len - mlen] = '\0'; /* Pauline's modification, want to keep
					name the same length (for PHYLIP) */
          sprintf(seqs[s]->name, "%s%d", ctemp, s);
          seqs[s]->name[SNAMELEN] = '\0';
          printf("Modified name to %s\n", seqs[s]->name);
          done = TRUE;
      } /* end of non-unique name */
      i--;
    }
  }  /* end of for s */
}  /*  end of fix_names */

/* copies start pos through gapped_end pos to newseqs from seqs.
   allocates the space.
	called in psiblast_res_to_fasta_db.c */
/************************************************************/
void
copy_region (Sequence* newseqs[MAXSEQ], Sequence* seqs[MAXSEQ], int nseqs,
		int gapped_start_pos, int gapped_end_pos)
{
	int i, pos;
	Residue c;

	for (i = 0 ;i < nseqs; i ++) {
		newseqs[i] = (Sequence*) malloc (sizeof (Sequence));
		strcpy (newseqs[i]->name, seqs[i]->name);
		strcpy (newseqs[i]->info, seqs[i]->info);
		newseqs[i]->length = gapped_end_pos - gapped_start_pos + 1;
		newseqs[i]->max_length = seqs[i]->max_length;
		newseqs[i]->sequence = (Residue*) calloc
			(newseqs[i]->max_length, sizeof (Residue));
		for (pos = gapped_start_pos; pos <= gapped_end_pos; pos++) {
			c = seqs[i]->sequence[pos];
			newseqs[i]->sequence[pos - gapped_start_pos] =c;
		}
	}

} /* end copy_region */

/*******************************************************************/
/* copies regions from gapped_start_pos to gapped_end_pos and adds it to
	newseqs*/

void
add_region (Sequence* newseqs[MAXSEQ], Sequence* seqs[MAXSEQ], int nseqs,
                int gapped_start_pos, int gapped_end_pos)

{
	int i, pos;
	Residue c;

	assert (newseqs[0] != NULL);
	for (i = 0; i < nseqs; i++) {
		assert (strcmp(newseqs[i]->name, seqs[i]->name) == 0);
		for (pos = gapped_start_pos; pos <= gapped_end_pos; pos++) {
			c = seqs[i]->sequence[pos];
			newseqs[i]->sequence[newseqs[i]->length++] = c;
		}
	}

} /* end add_region */

void
extract_seqs (Sequence* newseqs[MAXSEQ], Sequence* oldseqs[MAXSEQ], int nseqs,
 char names[MAXSEQ][SMALL_BUFF_LENGTH],
		int new_no_of_seq)
{

	int i, new_index, j;
	int string_length;

	new_index = 0;

	for (i = 0; i < new_no_of_seq; i++) {
		newseqs[i] = (Sequence *) malloc (sizeof (Sequence));
	}

	for (i = 0; i < nseqs; i++) {
		for (j = 0; j < new_no_of_seq; j++) {
			string_length = strlen (oldseqs[i]->name);
			if (strlen (names[j]) < string_length) {
				string_length = strlen (names[j]);
			}
			if (strncmp (oldseqs[i]->name, names[j], string_length)
					 == 0 ) {
                                copy_sequence (newseqs[new_index],  oldseqs[i]);
				new_index++;
			}
		}
	}

	assert (new_index == new_no_of_seq);

} /* end of extract_seqs */

void
copy_sequence (Sequence* newseq, Sequence* oldseq)
{
         int i;
	Residue c;

	strcpy (newseq->name, oldseq->name);
         strcpy (newseq->info, oldseq->info);
		newseq->position = oldseq->position;
		newseq->length = oldseq->length;
		newseq->max_length = oldseq->max_length;

/*         printf ("type %d\n", oldseq->type); */
		newseq->type = oldseq->type;

/*         printf ("%.3f weight %d lenght\n", oldseq->weight, newseq->length); */
		newseq->weight = oldseq->weight;

	if (newseq->sequence == NULL) {
		newseq->sequence = (Residue *) calloc (newseq->length,
					sizeof (Residue));
	}
	 for (i = 0; i < oldseq->length; i++) {
		c = oldseq->sequence[i];
		newseq->sequence[i] = c;
	}
} /* end of copy_sequence */

int
compare_alignments (Sequence* align1[MAXSEQ], Sequence* align2[MAXSEQ],
		    int align1_nseq, 		int align2_nseq)
/* make pairwise comparisons between align1 and align2
   returns TRUE if alignment matches, FALSE otherwise */
{

/*	int i, pos;
	char name[SMALL_BUFF_LENGTH];

	for (i = 0; i < align1_nseq; i++) {
		for ( j = i + 1; j < align1_nseq; j++) {
			strcpy (name1, align1[i]->name);
			strcpy (name2, align1[j]->name);
			index1 =get_index_from_seqs(name1, align2, align2_nseq);
			index2 =get_index_from_seqs(name2, align2, align2_nseq);
			status = compare_pairwise_alignments (align1[i],
				 align1[j], align2[index1], align2[index2]);
			if (status == FALSE) {
				fprintf (errofp, "alignment not the same between ");
				fprintf (errorfp, "%s and %s in the two alignments \n",
					name1, name2);
				exit (-1);
			}
		}
	}
	return TRUE;
*/
    return FALSE;
} /* end of compare_alignments */

int
compare_pairwise_alignments (Sequence* aln1_seq1, Sequence* aln1_seq2,
				Sequence* aln2_seq1, Sequence* aln2_seq2)
/* this sucks if sequences are not of same length and don't start and end
at the same place */
{
/*	assert (strcmp (aln1_seq1->name, aln2_seq1->name) == 0);
	assert (strcmp (aln1_seq2->name, aln2_seq2->name) == 0);
*/
/*	remove_gaps_from_pairwise (gapless_aln1_seq1, gapless_aln1_seq2,
				           aln1_seq1,     aln1_seq2);

	remove_gaps_from_pairwise (gapless_aln2_seq1, gapless_aln2_seq2,
                                       aln2_seq1,     aln2_seq2);
*/
/*	min_length = gapless_aln1_seq1->length;
	if (gapless_aln1_seq1->length > gapless_aln2_seq2->length) {
		printf ("WARNING: longer alignment in alignment 1\n");
		min_length = gapless_aln2_seq2->length;
	} else if (gapless_aln1_seq1->length < gapless_aln2_seq2->length){
		printf ("WARNING: longer alignment in alignment 2\n");
		min_length = gapless_aln1_seq1->length;
	}

	for (i = 0; i < min_length; i++) {
		assert (gapless_aln1_seq1->sequence[i] ==
			gapless_aln2_seq1->sequence[i]);
		if (gapless_aln1_seq2->sequence[i] !=
			    gapless_aln2_seq2->sequence[i]) {
			printf ("not matching res %d\n",i);
			printf ("for sequences %s and %s\n", gapless_aln1_seq2->name, gapless_aln2_seq1->name);
			return FALSE;
		}
	}

	return TRUE;
*/
    return FALSE;
} /* end of compare_pairwise_alignments */

int
get_index_from_seqs (char name[], Sequence* seqs[MAXSEQ],
			int nseq)
{
	int i;

	for (i = 0; i < nseq; i++) {
		if (strcmp (name, seqs[i]->name) == 0) {
			return i;
		}
	}
	fprintf (errorfp, "EXITING -- couldn't find %s in get_index_from_seqs\n", name);
	exit (-1);
}

void
cleanup_short_length_seqs (Sequence* seqs[MAXSEQ], int *nseqs)
{
	int newnseqs, oldnseqs, i;
	int min_length;
	int length;

	oldnseqs = *nseqs;
	newnseqs = oldnseqs;
	min_length = (int) ( 0.8* (double) seqs[0]->length);

	for (i = 0; i < oldnseqs && i < newnseqs; i++) {
		length = get_length (seqs[i]);
		if (length < min_length ) {
			free_sequence (seqs[i]);
			seqs[i] = (Sequence *) malloc (sizeof (Sequence) );
			seqs[i]->sequence = NULL;
			if (i != newnseqs -1) { /* not the last one, copy the
						last seq to this location */
				copy_sequence (seqs[i], seqs[newnseqs-1]);
				free_sequence (seqs[newnseqs-1]);
				i--; /* have to check this sequence */
			}
			newnseqs--;
		}
	}
	*nseqs = newnseqs;
} /* end cleanup_zero_length_seqs */
int
get_gapped_pos (Sequence* seq, int aa_pos)
{
	int index, aa_count;

	index = 0;
	aa_count = 0;

	while (aa_count < aa_pos) {
		if (seq->sequence[index] != aa_atob['-']) {
			aa_count++;
		}
		index++;
	}
	return index;
} /* end of get_gapped_pos */

void
remove_gaps_in_sequence (Sequence* seq)
{
	Residue* new_sequence;
	int new_length;
	int new_i, i;

	new_length = 0;
	new_i = 0;
	new_sequence = (Residue*) calloc (seq->length, sizeof (Residue));

	for (i = 0; i < seq->length; i++) {
		if (seq->sequence[i] != aa_atob['-']) {
			new_sequence[new_i++] = seq->sequence[i];
			new_length++;
		}
	}
	/* reassign gapless sequence to seq */
	free (seq->sequence);
	seq->sequence = new_sequence;
	seq->length = new_length;
}

void
make_query_first (Sequence* seqs[MAXSEQ], int nseqs)
{
        int i;
        Sequence* tmp;

        i=0;

        while ( i < nseqs && ((strstr (seqs[i]->name , "QUERY") == NULL)
			  && (strstr (seqs[i]->name, "CONSENSUS0") == NULL)) ) {
                i++;
        }
        assert (i < nseqs); /* if i == nseqs, then QUERY isn't in list*/
        tmp = seqs[0];
        seqs[0] = seqs[i];
        seqs[i] = tmp;
}


int
psiblast_pairwise (FILE* fp, Sequence * seqs[MAXSEQ], int length, int option_carve_gaps)
{
	char* strptr;
	char name[SMALL_BUFF_LENGTH];
	int done;
	int nseq;
	fpos_t filepos;
	char line[LARGE_BUFF_LENGTH];
	char info[SMALL_BUFF_LENGTH];
	Sequence* subject;
	int identity;
	double evalue;
	int int_tmp;

	nseq = 0;
	done = false;
	while (!feof (fp) && !done) {
		/* starting at at a protein match */
		/* Pauline added Aug. 23,2010 in case there's extra whitespace lines */
		 while ( !feof (fp) &&   fgets (line, LINE_LEN, fp) != NULL &&
                                 line[0] == '\n')
       		{
		        printf ("skipping line %s\n", line);
               /* keep reading in lines until at a newline */
       		}

		/* end add */
	/*	fgets (line, LARGE_BUFF_LENGTH, fp); Pauline commented out because at a nonzero line*/
 printf ("here 3423 %s\n", line);
		if (strncmp (line , ">lcl", 4) == 0) {

        printf ("in 2343 here %s\n" ,line);
                	strptr = strtok (line, "| \t\r\n");
			strptr = strtok ((char*) NULL, "| \t\n\0");
			strcpy (name, &strptr[0]);
                }
		else {
printf ("in here this else");
			strptr = strtok (line, " \t\r\n");
/* Pauline added Aug 23,2010 for psiblast processing 2nd word
  and then commented out for blast2.2.24 Oct 8,2010 */

/*		strptr = strtok ((char*) NULL, " \t\n\0"); */
		strcpy (name, &strptr[1]);  /* don't need > */
		}
printf ("name here is 23422 %s\n", name);
	printf ("name is 233333 %s\n", name);
/*		strcpy (name, &strptr[1]); */ /* don't need > */
printf ("name %s\n",name);


		if ((strptr = strtok ((char*) NULL, "\t\n\0")) != NULL) {
			strcpy (info, strptr);
		} else {info[0] = '\0'; }
		name[SMALL_BUFF_LENGTH -1] = '\0';
		info[SMALL_BUFF_LENGTH -1 ] = '\0';

/* get evalue added 06-24-02 */
		while (strstr (line, "Expect") == NULL) {
			fgets (line, LARGE_BUFF_LENGTH, fp);
		}
		strptr = NULL;
		strptr = strstr (line, "Expect = ");
		if (strptr == NULL) {
			strptr = strstr (line, "Expect(");
			sscanf (strptr, "Expect(%d) = %lf", &int_tmp, &evalue);
		} else if (strptr != NULL) {
/*printf ("strptr is %s\n", strptr);; */
		sscanf (strptr, "Expect = %lf", &evalue);
		} else {
			fprintf (errorfp, "Unable to read e-value, parsing incorrect");
                exit (-1);
        }
/* printf ("evalue is %lf\n", evalue); */
/* end get evalue */

		while (strstr (line, "Identities") == NULL) {
			fgets (line, LARGE_BUFF_LENGTH, fp);
		}

		/* get identity score added 07-31-00 */
	/*	printf ("%s dafds\n", line);  */
		strptr  = strtok (line, " \t\r\n");
		while (strstr (strptr, "(") == NULL) {
			strptr = strtok ((char *) NULL, " \t\r\n");
		}
		sscanf (strptr, "(%d)", &identity);
	/*	printf ("got identity %d\n", identity); */

	        fgets (line, LARGE_BUFF_LENGTH, fp); /* read newline */
		/* Aug.16, 2000 check if line says Frame for DNA blast results*/
		if (strstr (line, "Frame") != NULL) {
			printf ("%s here\n", line);
			fgets (line, LARGE_BUFF_LENGTH, fp);
		}
		/* end of what was added Aug. 16, 2000 */
		subject = read_1st_4lines_of_pairwise (fp, name, length, option_carve_gaps);
		subject->undefined = identity; /* aug. 16 2000 */
		subject->undefined_dbl = evalue;
		strcpy (subject->info, info);
		while (read_4lines_of_pairwise (fp, subject, option_carve_gaps))
		{
		}
		fgetpos (fp, &filepos);
		fgets (line, LARGE_BUFF_LENGTH, fp);
/*	printf ("linedfaw %s\n", line);*/
		/* March 23, 2001, read in extra newlines until
		encounter line with something written on it .  This
		was added to read Trent's blast format */

		while (line[0] == '\n' || line[0] == '\r') {
			fgets (line,LARGE_BUFF_LENGTH, fp);
		}

		/* This line could either be "Database:" (done) or
		   "Score = " (another alignment with the same subject)
		or >SUBJECT_MATCH (next match */
		/* Aug.26, 2000, added "Searching" as done so can check
		results from round to round */
/*printf ("andthe line is %s\n", line);	  */
		/* Pauline added Lambda August 23,2010 */
		if ( (strstr (line, "Database:") != NULL) ||
		    (strstr (line, "Searching") != NULL)  ||
		     (strstr (line, "Lambda") != NULL)  ) {
			seqs[nseq] = subject;
			nseq++;
			done = true;
		} else if (strstr (line, "Score") != NULL) {
			printf ("%s has more thanone alignment\n", name);
			psiblast_pairwise_more_alignments
					(fp, subject, &done, option_carve_gaps);
			seqs[nseq] = subject;
			nseq++;
			/* more alignments */
		} else {
			assert (strstr (line, ">") != NULL);
			printf ("%s name\n", subject->name);
			seqs[nseq] = subject;
			nseq++;
			fsetpos (fp, &filepos); /* rewind position */
		}
	} /* end while feof */
	return nseq;
}


/* Check that both sequences do not have amino acids at the same position.
To be used to determine whether alignment returned by blast for the same
sequence are two different local alignments or overlapping alignments (as
is the case for repeats ) */

int
nonoverlapping (Sequence* seq1, Sequence* seq2)
{
	int i;

	assert (seq1->length == seq2->length);
	if (strcmp (seq1->name, seq2->name ) != 0) {
		fprintf (errorfp, "seeing whether %s and %s overlap, should be same name\n", seq1->name, seq2->name);
		exit (-1);
	}
	for (i = 0; i < seq1->length; i++) {
		/* check that both sequence do not have amino acids at the
			same position */
		if (seq1->sequence[i] != aa_atob['-'] &&
		    seq2->sequence[i] != aa_atob['-']) {
			return false;
		}
	}
	return true;
}

/* given that the sequences are nonoverlapping & the same length,
takes seq2 residues and copies it to sequence 1 gap residues.

To be used for joining two local alignments */
void
merge_sequence_res (Sequence* seq1, Sequence* seq2)
{
	int i;

	for (i = 0; i < seq1->max_length; i++) {
		if (seq2->sequence[i] != aa_atob['-']) {
			seq1->sequence[i] = seq2->sequence[i];
		}
	}
} /* end merge_sequence_res */

/* at another sequence alignment for the same protein.  at line Score */
void
psiblast_pairwise_more_alignments (FILE* fp, Sequence* subject, int* done,
				int option_carve_gaps)
{
	char line[LARGE_BUFF_LENGTH];
	Sequence* newsubject;
	fpos_t filepos;

	fgets (line, LARGE_BUFF_LENGTH, fp);
        fgets (line, LARGE_BUFF_LENGTH, fp);
	/* July 3, 2002 check if line says Frame for DNA blast results*/
        if (strstr (line, "Frame") != NULL) {
                   /*   printf ("%s here\n", line); */
                        fgets (line, LARGE_BUFF_LENGTH, fp);
        }


	newsubject = read_1st_4lines_of_pairwise (fp, subject->name,
				subject->max_length, option_carve_gaps);
        while (read_4lines_of_pairwise (fp, newsubject,option_carve_gaps)) {
        }

	if (option_carve_gaps == true && nonoverlapping (subject, newsubject)) {
/*		printf ("nonoverlapping: MERGING %s\n", subject->name);
		print_sequence (subject);
		print_sequence (newsubject); */
		merge_sequence_res (subject, newsubject);
        }

	fgetpos (fp, &filepos);
        fgets (line, LARGE_BUFF_LENGTH, fp);

	/* Pauline 11/22/02 for unknown # of blanklines before second
		alignment */
	while (line[0] == '\n') {
		fgets (line, LARGE_BUFF_LENGTH, fp);
	}

        /* This line could either be "Database:" (done) or
                   "Score = " (another alignment with the same subject)
                or >SUBJECT_MATCH (next match */

        if (strstr (line, "Database:") != NULL) {
              *done = true;
        } else if (strstr (line, "Score") != NULL) {
                psiblast_pairwise_more_alignments
                              (fp, subject, done, option_carve_gaps);
                        /* more alignments */
        } else {
                assert (strstr (line, ">") != NULL);
                /* new protein match */
	        fsetpos (fp, &filepos); /* rewind position */
        	return;
	}
}

int
read_4lines_of_pairwise (FILE* fp, Sequence* subject, int option_carve_gaps)
{
	char line[LARGE_BUFF_LENGTH];
	int query_start, seq_length, i;
	char query_seq[LARGE_BUFF_LENGTH], subject_seq[LARGE_BUFF_LENGTH];
	char* strptr;
	int gap_shift;
	int subject_start;

	assert (fp != NULL);

	fgets (line, LARGE_BUFF_LENGTH, fp);
	if (line[0] == '\n') {
		return false; /* end of alignment */
	}

	if (strstr (line, "Query") == NULL) {
		printf ("%s line should have Query:\n", line);
		exit (-1);
	}
        strptr = strtok (line, " \t\r\n"); /* this should be Query */
        strptr = strtok ((char *) NULL, " \t\r\n"); /* this should be the start
                                                        pos relative to query*/
        query_start = atoi (strptr);
	query_start--; /* location in ARRAY */
	strptr = strtok ((char *) NULL, " \t\r\n"); /* query seq */
        strcpy (query_seq, strptr);
        seq_length = strlen (query_seq);

        fgets (line, LARGE_BUFF_LENGTH, fp); /* identity & + line */
        fgets (line, LARGE_BUFF_LENGTH, fp); /* subject line */
        assert (strstr (line, "Sbjct") != NULL);
        strptr = strtok (line, " \t\r\n"); /* this should be Sbjct */
        strptr = strtok ((char *) NULL, " \t\r\n"); /* this should be the start
                                                        pos relative to sbject*/
        subject_start = atoi (strptr);
	subject_start--; /* location in array*/
	strptr = strtok ((char *) NULL, " \t\r\n"); /* subject seq */
        strcpy (subject_seq, strptr);
	if (option_carve_gaps == TRUE) {
/* assign sequence aa only if there is NOT a gap inthe query sequence */
		/* query_start is location in array */
		gap_shift= 0;
       		 for (i = 0; i < seq_length; i++) {
       		         if (query_seq[i] != '-') {
				subject->sequence[query_start + i - gap_shift] =
                                         aa_atob[subject_seq[i]];
               		 } else {
				gap_shift++;
			 }

        	}
	} else if (option_carve_gaps == FALSE) {
	/* allocate more space if subject length is longer than the query's
		protein length */
		if (subject->length + seq_length > subject->max_length) {
			/* print_sequence (subject); */
			resize_sequence(subject);
		}
		gap_shift = 0;
		for (i = 0; i < seq_length; i++) {
			if (subject_seq[i] != '-') {
				subject->sequence[subject_start+i-gap_shift
					- subject->position  ] =
						aa_atob[subject_seq[i]];
				subject->length++;
			} else {
				gap_shift++;
			}
		}
	}

        fgets (line, LARGE_BUFF_LENGTH, fp); /* newline */
        assert (line[0] = '\n');

        return true;

} /* end of read_4lines_of_pairwise*/

Sequence*
read_1st_4lines_of_pairwise (FILE* fp, char subject_name[SMALL_BUFF_LENGTH],
				int length, int option_carve_gaps)
{
	Sequence* subject;
	char line[LARGE_BUFF_LENGTH];
	char* strptr;
	int query_start, subject_start;
	char query_seq[LARGE_BUFF_LENGTH]; char subject_seq[LARGE_BUFF_LENGTH];
	int i, seq_length, gap_shift;
	int subject_end, query_end;

	CheckMem (subject = (Sequence*) malloc (sizeof (Sequence)));
	CheckMem (subject->sequence=(Residue*)calloc(length, sizeof (Residue)));
	subject->max_length = length;
	subject->length = length;
	subject->info[0] = '\0';
	subject->type = AA_SEQ;
	subject->position = 0;
	subject->weight = 0.0;
	subject->undefined = 0;
	subject->undefined_ptr = NULL;

	for (i = 0; i < length; i++) {
		subject->sequence[i] = aa_atob['-'];
	}
	strcpy (subject->name, subject_name);
	assert (fp != NULL);
	fgets (line, LARGE_BUFF_LENGTH, fp);
	if (strstr (line, "Query") == NULL) {
		printf ("%s line should have Query:\n", line);
		exit (-1);
	}
	strptr = strtok (line, " \t\r\n"); /* this should be Query */
	strptr = strtok ((char *) NULL, " \t\r\n"); /* this should be the start
							pos relative to query*/
	query_start = atoi (strptr);
	query_start--; /* location in ARRAY */
	strptr = strtok ((char *) NULL, " \t\r\n"); /* query seq */
	strcpy (query_seq, strptr);
	seq_length = strlen (query_seq);
	strptr = strtok ((char *) NULL, " \t\r\n");
	query_end = atoi (strptr);
	fgets (line, LARGE_BUFF_LENGTH, fp); /* identity & + line */
	fgets (line, LARGE_BUFF_LENGTH, fp); /* subject line */
/*	printf ("here: %s\n", line); */
	assert (strstr (line, "Sbjct") != NULL);
        strptr = strtok (line, " \t\r\n"); /* this should be Sbjct */
        strptr = strtok ((char *) NULL, " \t\r\n"); /* this should be the start
                                                        pos relative to sbject*/
	subject_start = atoi(strptr);
	subject_start--; /*location in array*/
	subject->position = subject_start;
	strptr = strtok ((char *) NULL, " \t\r\n"); /* subject seq */
	strcpy (subject_seq, strptr);
	strptr = strtok ((char *) NULL, " \t\r\n");
	subject_end = atoi (strptr);
        if (option_carve_gaps == TRUE) {
/* assign sequence aa only if there is NOT a gap inthe query sequence */
                /* query_start is location in array */
                gap_shift= 0;
                 for (i = 0; i < seq_length; i++) {
                         if (query_seq[i] != '-') {
                                subject->sequence[query_start + i - gap_shift] =
                                         aa_atob[subject_seq[i]];
                         } else {
                                gap_shift++;
                         }

                }
        } else if (option_carve_gaps == FALSE) {
                subject->length = 0;
		gap_shift = 0;
		for (i = 0; i < seq_length; i++) {
                        if (subject_seq[i] != '-') {
				subject->sequence[subject_start + i -gap_shift
						- subject->position  ]=
						 aa_atob[subject_seq[i]];
                		subject->length++;
			} else {
				gap_shift++;
			}
		}
        }


	fgets (line, LARGE_BUFF_LENGTH, fp); /* newline */
	assert (line[0] = '\n');
	return subject;


} /* end of read_1st_3lines_of_pairwise*/

void
convert_gap_to_X (Sequence* seq) {

	int i;

	for (i = 0 ; i < seq->length; i++) {
		if (seq->sequence[i] == aa_atob['-']) {
			seq->sequence[i] = aa_atob['X'];
		}
	}
}

void
output_sequence_without_gaps_or_Xes (Sequence* seq, FILE* osfp)
{

  int k;
  int length_printed;

  length_printed = 0;
  fprintf(osfp, ">%s  %s\n", seq->name, seq->info);

  if (seq->type == AA_SEQ) {
    for(k=0; k<seq->length; k++) {
      if ((aa_btoa[seq->sequence[k]] != '-') &&
			(aa_btoa[seq->sequence[k]] != 'X') ) {
		fprintf(osfp, "%c", aa_btoa[seq->sequence[k]]);
      		length_printed++;
		if ((length_printed+1)%60 == 0) {
        	fprintf(osfp, "\n");
      		}
    	}
    } /* end for */
  }
  else if (seq->type == NA_SEQ) {
    for(k=0; k<seq->length; k++) {
      fprintf(osfp, "%c", nt_btoa[seq->sequence[k]]);
      if ((k+1)%60 == 0) {
        fprintf(osfp, "\n");
      }
    }
  }
  else {
    for(k=0; k<seq->length; k++) {
      if ((k+1)%60 == 0) {
        fprintf(osfp, "\n");
      }
      fprintf(osfp, "%c", aa_btoa[seq->sequence[k]]);
    }
  }

  fprintf(osfp, "\n");
} /* end output_sequence_without_gapes_or_Xes */

void
convert_sequence_with_beg_and_end_gaps_to_X (Sequence* seq)
{
	int pos;
	Boolean done;

	done = FALSE; pos = 0;
	while (pos < seq->length && !done ) {
		if (aa_btoa[seq->sequence[pos]] == '-') {
			seq->sequence[pos] = aa_atob['X'];
 		} else {
			done = TRUE;
		}
		pos++;
	}
	done = FALSE; pos = seq->length - 1;
	while (pos > 0 && !done) {
		if (aa_btoa[seq->sequence[pos]] == '-') {
			seq->sequence[pos] = aa_atob['X'];
		} else {
			done = TRUE;
		}
		pos--;
	}
}

void
output_sequence_without_beg_and_end_X (Sequence* seq, FILE* osfp)
{

  int k;
  int length_printed;
  int beg, end;
  Boolean done;

  done = FALSE; k = 0;
   while (k < seq->length && !done) {
	if (aa_btoa[seq->sequence[k]] != 'X') {
		done = TRUE;
		k--;
	}
	k++;
  }
  beg = k;
  done = FALSE; k = seq->length - 1;
  while (k > 0 && !done) {
	if (aa_btoa[seq->sequence[k]] != 'X') {
		done = TRUE;
		k++; /* to compensate for following k--*/
	}
	k--;
  }
 end = k;

  length_printed = 0;
  fprintf(osfp, ">%s  %s\n", seq->name, seq->info);

  if (seq->type == AA_SEQ) {
    for(k=beg; k<= end; k++) {
         fprintf(osfp, "%c", aa_btoa[seq->sequence[k]]);
         length_printed++;
         if ((length_printed+1)%60 == 0) {
                fprintf(osfp, "\n");
         }
    } /* end for */
  }
  else if (seq->type == NA_SEQ) {
    for(k=0; k<seq->length; k++) {
      fprintf(osfp, "%c", nt_btoa[seq->sequence[k]]);
      if ((k+1)%60 == 0) {
        fprintf(osfp, "\n");
      }
    }
  }
  else {
    for(k=0; k<seq->length; k++) {
      if ((k+1)%60 == 0) {
        fprintf(osfp, "\n");
      }
      fprintf(osfp, "%c", aa_btoa[seq->sequence[k]]);
    }
  }

  fprintf(osfp, "\n");
} /* end output_sequence_without_gapes_or_Xes */

int
get_length_nogap_noX (Sequence* seq)
{
	int i;
	int length;

	length = 0;

	for (i = 0; i < seq->length; i++) {
		if (seq->sequence[i] != aa_atob['X'] &&
			seq->sequence[i] != aa_atob['-']) {
			length++;
		}
	}
	return length;

}

void output_sequence_clean (seq, outfp)
     Sequence *seq;
     FILE *outfp;
{
  int k;


  fprintf(outfp, ">%s  %s\n", seq->name, seq->info);

  if (seq->type == AA_SEQ) {
    for(k=0; k<seq->length; k++) {
      fprintf(outfp, "%c", aa_btoa[seq->sequence[k]]);
      if (  ((k+1)%60 == 0) && (k != seq->length - 1)) {
        fprintf(outfp, "\n");
      }
    }
  }
  else {
    for(k=0; k<seq->length; k++) {
      if ((k+1)%60 == 0) {
        fprintf(outfp, "\n");
      }
      fprintf(outfp, "%c", aa_btoa[seq->sequence[k]]);
    }
  }

  fprintf(outfp, "\n");
}

int
cluster (double clus, Sequence* seqs[MAXSEQ], int nseq, FILE* outfp,
	FILE* keyoutfp, int option)
{

   int iclus, npair, threshold, s1, s2, l1, l2, px, i, i1, i2, first, j;
   int check1, check2, temp_cluster_no, tempint;
   int nclus[MAXSEQ], icluster[MAXSEQ], minclus, oldclus, width;
   int match, no_of_clumps;
   struct cluster_pair *pairs;
  int alignment_width;
  int nseqs_cluster;
  Sequence* seqs_cluster[MAXSEQ]; Block* block; Matrix* pssm;
        Sequence* consensus_sequence;
        char comment[SMALL_BUFF_LENGTH];

   no_of_clumps = 0;
   npair = nseq*(nseq-1)/2;
   pairs = (struct cluster_pair *) malloc(npair * sizeof(struct cluster_pair));

   /*    Compute scores for all possible pairs of sequences            */
   for (s1=0; s1<nseq-1; s1++)                  /* col = 0, n-2     */
   {
      for (s2=s1+1; s2<nseq; s2++)      /* row = col+1, n-1 */
      {
         alignment_width = seqs[0]->length;
         match = 0;
         px = INDEX(nseq, s1, s2);
         pairs[px].score = 0;
         pairs[px].cluster = -1;
         for (i=0; i< seqs[0]->length; i++)
         {
                if ( seqs[s1]->sequence[i] == seqs[s2]->sequence[i]) {
                        if (seqs[s1]->sequence[i] == aa_atob['X'] ||
                              seqs[s1]->sequence[i]== aa_atob['-']) {
                                alignment_width--;
                        } else {
                                match++;
                        }
                } /* end of if letters match */
          } /* end of for width */
         pairs[px].score = ((double) match)/ ((double) alignment_width);
      }  /* end of s2 */
   }  /* end of s1 */

   /*-------Cluster if score exceeds threshold by scanning cols (s1) */
   for (s1=0; s1<nseq; s1++)
   {
      icluster[s1] = -1;                        /* clear out old values */
      nclus[s1] = 0;
   }
   iclus = 0;                                   /* cluster number */
   for (s1=0; s1<nseq-1; s1++)                  /* col = 0, n-2     */
      for (s2=s1+1; s2<nseq; s2++)      /* row = col+1, n-1 */
      {
         px = INDEX(nseq, s1, s2);
         if (pairs[px].score >= clus)      /*  cluster this pair */
         {
            if (icluster[s1] < 0)          /* s1 not yet clustered */
            {
               if (icluster[s2] < 0)       /* new cluster */
               {
                 icluster[s1] = iclus++;
                  icluster[s2] = icluster[s1];
                }
               else                             /* use s2's cluster  */
                  icluster[s1] =  icluster[s2];
            }
            /*  use s1's cluster if it has one and s2 doesn't */
            else if (icluster[s1] >= 0 && icluster[s2] < 0)
               icluster[s2] = icluster[s1];
            /* merge the two clusters into the lower number */
            else if (icluster[s1] >= 0 && icluster[s2] >= 0)
            {
               minclus = icluster[s1]; oldclus = icluster[s2];
               if (icluster[s2] < icluster[s1])
               {
                  minclus = icluster[s2]; oldclus = icluster[s1];
               }
               for (i1=0; i1<nseq; i1++)
                 if (icluster[i1] == oldclus)
                     icluster[i1] = minclus;
            }
         }  /* end of if pairs */
      }  /* end of s2 */

   /*---  Set ncluster, get rid of negative cluster numbers --*/
   for (s1=0; s1<nseq; s1++)
   {
      if (icluster[s1] < 0)
          icluster[s1] = iclus++;
   }
   for (s1=0; s1<nseq; s1++)
          nclus[icluster[s1]] += 1;

   i2 = 0;
/* printf ("no of clusters %d\n", iclus); */

   /* move QUERY to first cluster */
	if (icluster[0] != 0) { /* QUERY is not first in cluster */
printf ("entered if statment\n");
		check1 = 0; check2 = 0;
		temp_cluster_no = icluster[0];
		assert (temp_cluster_no != 0);
		tempint = nclus[temp_cluster_no];
				/* # of sequences grouped with QUERY */
		nclus[temp_cluster_no] = nclus[0];
		nclus[0] = tempint;
		for (s1 =0; s1 < nseq; s1++) {
			if (icluster[s1] == temp_cluster_no) {
				/* same cluster as seq 0*/
				icluster[s1] = 0; check1++;
			} else if (icluster[s1] == 0) {
				/* in cluster 0, moving it to where
				cluster number where query was */
			icluster[s1] = temp_cluster_no;
				check2++;
			}
		}
		assert (check1 == nclus[0]); /* same # of sequences with query*/
		assert (check2 == nclus[temp_cluster_no]) ;
			/* move all the sequences in cluster 0 to query */
	} /* end of if QUERY not infirst cluster */
   for (i=0; i< iclus; i++)
   {
      if (nclus[i] > 1 || (nclus[i] > 0 && option == 0))
	 /* print clumps that have 2 or more sequences or
	if option 0, print out clumps with 1 or more sequences */
      {
        i2++;
         nseqs_cluster = 0;
        for (s1=0; s1 < nseq; s1++)
         {
            if (icluster[s1] == i)
            {
		seqs_cluster[nseqs_cluster] = seqs[s1];
                nseqs_cluster++;
/*              output_sequence (seqs[s1], outfp); */
            }
         } /* end for nseqs */
        /* get consensus for cluster i */
        block = make_block (seqs[0]->length, 0, nseqs_cluster, seqs_cluster, FALSE);
        pb_weights (block);
        pssm = PN_block_to_matrix (block, 2);
        consensus_sequence = get_consensus (pssm);
/*              fprintf (outfp, "CONSENSUS for cluster %d \n", i); */
        /* August 28, 2000 added 1 to i below so that name starts as
	CONSENSUS1, CONSENSUS2, ... because QUERY will be represented
	by itself */
	sprintf (comment, "%d", i);
        strcat (consensus_sequence->name, comment);
        convert_gap_to_X (consensus_sequence);
	if (seq_not_all_Xes (consensus_sequence) ) {
	        output_sequence_without_beg_and_end_X (consensus_sequence, outfp);
		no_of_clumps++;
	}

/* now print out consensus key */
	fprintf (keyoutfp, "%s\n", consensus_sequence->name);
	for (j = 0; j < nseqs_cluster; j++) {
		fprintf (keyoutfp, "%s %s\n", seqs_cluster[j]->name, seqs_cluster[j]->info);
	}
	fprintf (keyoutfp, "\n");
        free_block (block);
        free_matrix (pssm);
        free_sequence (consensus_sequence);
     } /* end if nclus > 1 || nclus[i] > 0 && option == 0 */
   }


 printf ("Clumped %d sequences into %d clusters\n", nseq, i2);
 /* keep this print because get # of clustersin awk statement */
  free(pairs);
printf ("exiting clump\n");
 return no_of_clumps;
}  /* end of cluster */

int
seq_not_all_Xes (Sequence* seq)
{
	int pos;

	pos = 0;
	while ( seq->sequence[pos] == aa_atob['X']) {
		pos++;
	}
	if (pos < seq->length) {
		return TRUE;
	} else {
		/* gone through all of sequence, and it's all Xes*/
		return FALSE;
	}
} /* end of seq_not_all_Xes */

void sort_seq_by_identity (Sequence* seqs[MAXSEQ], int num_seqs)
{
/* assumes identity is in seq->undefined
 implementation crude -- wanted to code this up quickly */

	int j, P;
	Sequence* tmp_seq;

	for (P = 1; P < num_seqs; P++) {
		tmp_seq = seqs[P];
		for (j = P ; j > 0 && seqs[j-1]->undefined < tmp_seq->undefined; j--) {
			seqs[j] = seqs[j-1];
		}
		seqs[j] = tmp_seq;
	}
} /* end of sort_seq_by_identity */


#endif
