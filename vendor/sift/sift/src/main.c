#ifdef SIFT_MAIN_

#define EXTERN

#include <stdlib.h>
#include <stdio.h>

#include "info_on_seqs_pid.h"

#define MAXSEQ 400 /* Maximum number of sequences */
#define LINE_LEN 800

FILE* errorfp;
char errorfilename[LARGE_BUFF_LENGTH];

void getargs (int argc, char* argv[], FILE** seqfp, FILE** polymorphfp,
	char outfilename[LARGE_BUFF_LENGTH], int* seq_identity);

int main (int argc, char* argv[])
{
	FILE* seqfp; FILE* outfp; FILE* polymorphismfp;
	char outfilename[LARGE_BUFF_LENGTH], currentoutfile[LARGE_BUFF_LENGTH];
	char tempname[LARGE_BUFF_LENGTH];
	char desc[SMALL_BUFF_LENGTH]; int desc_length; char* strptr;
	Sequence *seqs[MAXSEQ];
	int nseqs, i, pos;
	int db_type;
	int seq_type;
	int seq_identity;
	int original_aa;

	ErrorLevelReport = 5;

	init_frq_qij();
	printf ("tell me i've entered\n");
	getargs (argc, argv, &seqfp, &polymorphismfp, outfilename, &seq_identity);

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
        int h;
        for (h = 0; h < seqs[0]->length; ++h) {
            printf("%c", seqs[0]->sequence[h] + 'A');
        }
        printf("\n");
    }
    else /* CLUSTAL or MSF? */
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
		{
			nseqs = try_msf(seqfp, seqs, desc);
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

	generate_predictions(seqs, nseqs, polymorphismfp, seq_identity, outfp);

	if (polymorphismfp != NULL) close(polymorphismfp);
	fclose(seqfp);
	fclose(outfp);
	fclose (errorfp);
	rm_file (errorfilename);
	exit (0);
} /* end main */

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
    if (substfilename[0] == '-')
	{
		*polymorphfp = NULL;
	}
	else if ((*polymorphfp = fopen (substfilename, "r")) == NULL)
    {
        printf ("cannot open file %s \n", substfilename);
        exit (-1);
    }


	if (argc > 3)
	{
		strcpy (outfilename, argv[3]);
	}
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

	if (argc > 4)
	{
		*seq_identity = atoi (argv[4]);
	}
    else
	{
        printf ("100% sequence identity that will be filtered out by default because argument not set\n");
        *seq_identity = 100; /* default of 100% */
        /*scanf ("%d", seq_identity); */
    }
} /* end of getargs */

#endif
