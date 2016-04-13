/* (C) Copyright 2000, Fred Hutchinson Cancer Research Center */
/* Use, modification or distribution of these programs is subject to */
/* the terms of the non-commercial licensing agreement in license.h. */

#ifndef MATRIX_INFO_C_
#define MATRIX_INFO_C_

#include <string.h>

#include "Matrix_Info.h"
#include "Array_math.h"

#define THRESHOLD 0.05
#define TRUE 1
#define FALSE 0


/* removed read_matrix_header() and read_matrix_20_aa_body() see july 2nd copy if
need this again (for Kevin Karplus pssms)

July 17, 2001 from Chris Saunders -- he was getting a core dump for info_on_seqs
but he was still able to get results.  cleaning up according to him, check this when I get back from D.C.  in free_Matrix_Info,

1) commented out free (matrix->logo_heights)
2) free (matrix->info)
3) free_matrix (matrix->block_matrix)

*/

Matrix_Info* initialize_matrix_info (Matrix* matrix)
{
	Matrix_Info* new_matrix_info;
	int block_width;
	int aa, pos;

	new_matrix_info = (Matrix_Info*) calloc (1, sizeof (Matrix_Info) );
	new_matrix_info->block_matrix = matrix;

	block_width = matrix->width;

	new_matrix_info->R = (double*) calloc(block_width, sizeof (double));
	new_matrix_info->avg_R_per_residue = (double*) calloc (block_width,
						sizeof (double));

	for (pos=0; pos < block_width; pos++) {
		new_matrix_info->R[pos] = 0;
		new_matrix_info->avg_R_per_residue[pos] = 0;
	}

	/* allocate space for info array */
	for (aa=0; aa < MATRIX_AA_WIDTH; aa++) {
		new_matrix_info->info[aa] = (int *) calloc (block_width,
								sizeof (int));
		new_matrix_info->logo_heights[aa] = (double *) calloc
						(block_width,sizeof (double));
	}

	for (aa=0; aa < MATRIX_AA_WIDTH; aa++) {
		for (pos = 0; pos < block_width; pos++) {
			new_matrix_info->info[aa][pos] = no_data;
			new_matrix_info->logo_heights[aa][pos] = 0;
		}
	}

	return new_matrix_info;
}

void free_Matrix_Info (Matrix_Info* matrix)
{
	int aa;

	for (aa = 0; aa < MATRIX_AA_WIDTH; aa++) {
		free (matrix->info[aa]);
		free (matrix->logo_heights[aa]);
	}


	free( matrix->R);
        free (matrix->avg_R_per_residue);

/*July 17, 20001 commented out recommended by Chris Sanders
on linux
 	free (matrix->logo_heights);
	free (matrix->info);

	free_matrix (matrix->block_matrix);
*/
	matrix->block_matrix = NULL;
	free (matrix);
	matrix = NULL;

}

/*======================================================================
calculate_information_R (Matrix_Info* matrix,int error_correction_option)

calculates the total information at position in R array and fills in
height of aa at position pos height = f (aa, pos) * R (pos) in
2-dim. array logo_heights
error_correction_option : whether or not to include error correction
=======================================================================*/

double calculate_information_R (Matrix_Info* info_matrix,
				int error_correction_option)
{
        double ln2, e, r, hmax, totalweight, dtemp;
        int aa, pos, seq, num_residues;
        double max_aa_height, aa_height;
	Matrix* matrix;
	int num_diff_residues;
	double total_R;

	total_R = 0.0;

	matrix = info_matrix->block_matrix;
        ln2 = log (2.0);
        hmax = log((double) AAs);
     for(pos=0; pos < matrix->width; pos++)
     {
        /* =====Calculating error correction for small sample size ===*/
       /*  count actual number of residues in each column.
        Including the 20 aa ('residues' matrix values 1-20)
        and B (Asp or Asn, value 21) and Z (Glu or Gln, value 22)
        Excluding gaps (-, value 0) unidentified aa (X, value 23) and
        stop codon (*, value 24). */
	num_residues = matrix->num_sequences ;

/* correction factor for PSSM assumes no X's.
  Basically, subtract from num_residues the X's and gaps
  but in PSSM, don't know these counts because don't
  have the block.
*/
  if (matrix->block != NULL) {

        for(seq=0; seq < matrix->num_sequences; seq++)
            if (matrix->block->residues[seq][pos] < 1 ||
                matrix->block->residues[seq][pos] > 22) num_residues-- ;
  }
                                /* correction factor for small sample size */
        e = (AAs-1) / ((double) 2 * ln2 * num_residues) ;

        /* =====end of calculating error correction ================*/
        assert (matrix->weights != NULL);

        r = hmax ; /* R = log 20 */
	 for(aa=1, totalweight=0.; aa < AAs+1; aa++) {
           totalweight +=  matrix->weights[aa][pos]  ;
        }
        /* totalweight is the sum of all aa weights at a certain pos */
        for(aa=1; aa < AAs+1; aa++)
                    /* start loop at 1 and end at AAs+1 because aa values
                      at the weight matrix start at 1, 0 is the gap value */
       {
         if (matrix->weights[aa][pos] > 0)
            {
            dtemp = (double) matrix->weights[aa][pos] / totalweight;
           /* in Tom Schneider's logo eq., dtemp corresponds to the frequency
                of an amino acid */
        /* dtemp is the proportion of that particular aa relative to all
          the aa's (total weight ) */
            r += dtemp * log(dtemp); /* sum over all residues dtemp*log(dtemp)*/
                                    /* is H(pos) */
            } /* end of if weights > 0 */
         } /* end of for for calculating r*/

      r /= ln2 ;       /* convert to bits , R = 2 - H(pos) in bits now */
/* if (r > 4) { printf ("pos %d exceeded over 4 bits\n", pos); } */
	if (error_correction_option) {
	      r -= e ;                     /* a correction for small sample sizes */
	}

	info_matrix->R[pos] = r;
	total_R += r;

	num_diff_residues = 0;
      for(aa=1; aa < AAs+1; aa++) {
         if (matrix->weights[aa][pos] > 0 ) {
		num_diff_residues++;
	 }
	 aa_height = r * ((double) matrix->weights[aa][pos])/ (double) 100;
         info_matrix->logo_heights[aa][pos] = aa_height;
	}
/* printf ("r %.2f %d num_diff_residues\n", r, num_diff_residues); */
	info_matrix->avg_R_per_residue[pos] = r/ (double) num_diff_residues;
/*	printf ("avg %.2f \n", info_matrix->avg_R_per_residue[pos]); */

   } /* end of for-loop for all positions */
   return total_R;

} /* end of calculate_information_R */

double
calculate_max_info_possible (int num_sequences)
{
        double r, e;
	double ln2;

	ln2 = log(2.0);

	r = log((double) AAs)/ ln2;
        e = (AAs-1) / ((double) 2 * ln2 * (double) num_sequences) ;
        r -= e ;                     /* a correction for small sample sizes */
	printf ("for %d sequences, max info %.2f\n", num_sequences, r);
	return r;

}


/* converts allowed, less_severe, more_severe, and less_severe to numerical
values so that a Pearson correlation coefficient can be calculated
*/

double
correlation_coefficient (Matrix* matrix_x, Matrix* matrix_y )
{
	double SumX, SumY, SumXsquared, SumYsquared, SumXY, x, y, r, num, denom;
	int n, aa, pos;

	SumX = SumY = SumXsquared = SumYsquared = SumXY = x = y = 0.0;
	n = 0;

	assert (matrix_x->width == matrix_y->width);

	printf ("POS\tX\tY\n");
	for (aa = 0; aa < MATRIX_AA_WIDTH; aa++) {
		for (pos = 0; pos < matrix_x->width; pos++) {
			x = matrix_x->weights[aa][pos];
			y = matrix_y->weights[aa][pos];
			if (x != NO_DATA && y != NO_DATA){
				printf ("%d\t%.3f\t%.3f\n", pos, x, y);
				n++;
				SumX += x;
				SumY += y;
				SumXsquared += (x*x);
				SumYsquared += (y*y);
				SumXY += (y*x);
			}
		}
	}

	num = SumXY - (SumX * SumY / n);
	denom = (SumXsquared - (SumX * SumX / n)) * (SumYsquared - (SumY * SumY / n));
	denom = sqrt (denom);
	r = num/denom;

	return r;

} /* end of correlation coefficient */

void
convert_experimental_info_to_values (Matrix_Info* info_matrix)
{
   int aa, pos;
   int phenotype;

        for (aa=0; aa < MATRIX_AA_WIDTH; aa++) {
                for (pos = 0; pos < info_matrix->block_matrix->width; pos++) {
                        phenotype = info_matrix->info[aa][pos];
			switch (phenotype) {
			   case allowed:
				info_matrix->block_matrix->weights[aa][pos] = 1.0; break;
			   case less_severe:
				info_matrix->block_matrix->weights[aa][pos] = 2.0/3.0; break;
			   case more_severe:
				info_matrix->block_matrix->weights[aa][pos] = 1.0/3.0; break;
			   case not_allowed:
				info_matrix->block_matrix->weights[aa][pos] = 0.0; break ;
			    case no_data:
				info_matrix->block_matrix->weights[aa][pos] = NO_DATA ; break ;
			    default:
				fprintf (errorfp, "Error!  what type of data is this? pos %d aa %c phenotype %d\n", pos, aa_btoa[aa], phenotype);
				exit (-1);
			} /* end of switch */
                }
        }

}


void read_experimental_info (Matrix_Info* matrix, FILE* infofp,
			int block_beg_pos_relative_to_seq,
			int block_end_pos_relative_to_seq,
			const fpos_t start_of_file)
{
	int aa, i, aa_pos, pos;
	char line[LARGE_BUFF_LENGTH];
	char* stringp;
	char *string_pos, result[10];
	char original_aa, substituted_aa;
	char word1[10], word2[10];

/*	printf ("entered read_experimental_info\n"); */

	fsetpos (infofp, &start_of_file);
	while (fgets (line, LARGE_BUFF_LENGTH, infofp) != NULL) {
		if (line[0] != '%') {
		word1[0] = '\0'; word2[0] = '\0', result[0]='\0';
		stringp = strtok (line, " \r\n\t");
	         strcpy (word1, stringp);
		if (( stringp = strtok ((char *) NULL, " \r\n\t" )) != NULL) {
                        strcpy (result, stringp);
                }
	/* parse original amino acid, position, and new amino acid */
	/* Example of what file should read: Y2P => tyrosine in original
	sequence at position 2 substituted for Proline */
		original_aa = word1[0];
		strcpy (word2, &word1[1]);
		pos = atoi (word2);
		stringp = strpbrk (word2, "ACDEFGHIKLMNPQRSTVWY");
		substituted_aa = *stringp;

	/* make sure POSITION IS IN RANGE OF PSSM */
/*printf ("pos %d original_aa %c substituted %c\n", pos, original_aa, substituted_aa);; */
	if (pos >= block_beg_pos_relative_to_seq &&
		pos <= block_end_pos_relative_to_seq) {

		if (original_aa == substituted_aa) {
			if (strcmp ("+" , result) != 0);
				printf ("pos %d original subst not allowed %s\n", pos, result);
/*				exit (-1); */
		}

		if (strstr (result, "+-") != NULL) {
			matrix->info[aa_atob[substituted_aa]][pos -
					block_beg_pos_relative_to_seq]
						= less_severe;
		} else if (strstr (result, "-+") != NULL) {
			matrix->info[aa_atob[substituted_aa]][pos -
				block_beg_pos_relative_to_seq] = more_severe;
		} else if (strstr (result, "-") != NULL) {
			matrix->info[aa_atob[substituted_aa]][pos -
					block_beg_pos_relative_to_seq]
							= not_allowed;
		} else if (strstr (result, "+" ) != NULL) {
			matrix->info[aa_atob[substituted_aa]][pos -
					block_beg_pos_relative_to_seq ]
							= allowed;
		}
	} /* end of check if position is in range */
      } /* end of if not a comment % line */
	} /* end of while reading file */
/*	printf ("exited reading experimantl info\n"); */
} /* end of read_experimental_info */

/**********************************************************************************/
void reclassify_Matrix_info (Matrix_Info* info_matrix, int severity_option)
{
/* severity_option == 1 => group +- phenotype with +, i.e. classify
	less_severe as allowed
   severity_option == 2 => group +- phenotype with -, i.e. classify
	less_severe as not_allowed
*/
 Matrix* matrix;
 int l; char c;

 matrix = info_matrix->block_matrix;

  for (l=0; l<matrix->width; l++) {
    for (c='A'; c <= 'Z'; c++) {
       if ((c != 'J') && (c != 'O') && (c != 'U') && (c != 'B') &&
           (c != 'X')  && (c != 'Z') ) {
           if (info_matrix->info[aa_atob[c]][l] ==
				less_severe && severity_option == 1){
                        info_matrix->info[aa_atob[c]][l] = allowed;
	   } else if (info_matrix->info[aa_atob[c]][l] ==
				less_severe &&
					 severity_option == 2) {
			info_matrix->info[aa_atob[c]][l] = not_allowed;
	   } else if (info_matrix->info[aa_atob[c]][l] == more_severe) {
			info_matrix->info[aa_atob[c]][l] = not_allowed;
	   }
	}
    } /* end of for alphabet */
  } /* end of for positons in matrix */
}

void output_Matrix_Info (Matrix_Info* info_matrix, FILE* omfp, int option, int severity_option,
			int style)
{

  Matrix* matrix;
  int error_fp, error_tn, correct, total;
  int notobserved_intolerant, observed_tolerant;
  int pos_observed_tolerant, pos_error_tn, pos_notobserved_intolerant,
	pos_error_fp, pos_total;
  char c; int l;
  int print_correct, print_incorrect, print_X, print_type;

	print_correct = 1; print_incorrect = 2; print_X =3;
  correct = 0; error_fp=0; error_tn = 0, total=0;
  notobserved_intolerant = 0; observed_tolerant = 0;

  matrix = info_matrix->block_matrix;
  fprintf (omfp, "severity option %d\n", severity_option);

  fprintf(omfp, "ID   %s\n", matrix->id);
  fprintf(omfp, "AC   %s\n", matrix->ac);
  fprintf(omfp, "DE   %s\n", matrix->de);
  fprintf(omfp, "MA   %s\n", matrix->ma);

  if (style == FLOAT_OUTPUT) {
	output_Matrix_Info_float (info_matrix, omfp, option, severity_option);
	return;
  }

    for (c='A'; c <= 'Z'; c++) {
      if ((c != 'J') && (c != 'O') && (c != 'U') && (c != 'B') && (c != 'X')
		&& (c != 'Z') )
	{
          fprintf(omfp, "  %c ", c);
        }
     } /* end of for loop,finished printing amino acids */
    if (style == FLOAT_OUTPUT) fprintf(omfp, "    *         -\n");
    else                       fprintf(omfp, "  R   -\n");

    	for (l=0; l<matrix->width; l++) {
         /* score accuracy for position */
	pos_observed_tolerant = 0; pos_notobserved_intolerant = 0;
	 pos_error_tn = 0; pos_error_fp = 0; pos_total = 0;
	 for (c='A'; c <= 'Z'; c++) {
           if ((c != 'J') && (c != 'O') && (c != 'U') && (c != 'B') &&
		(c != 'X')  && (c != 'Z') ) {
  	   switch (option) {
	     case 1: /* for PSSMs with negative scores.  predict scores <
			as intolerant.  for log odds and matrices.
			0 is neutral (therefore tolerant) */
			/* NO DATA */
		if (info_matrix->info[aa_atob[c]][l] == no_data ) {
			/* fprintf (omfp, "no data\n"); */
			print_type = print_X;
                } else if (info_matrix->info[aa_atob[c]][l] == allowed) {
			if ( matrix->weights[aa_atob[c]][l] < 0)
			{    /* ERROR of NOT OBSERVED, BUT TOLERANT true negative */
				error_tn++; print_type = print_incorrect;
				pos_error_tn++;
			} else if ( matrix->weights[aa_atob[c]][l] >= 0)
			{ /* predicted correctly */
				 observed_tolerant++; pos_observed_tolerant++;
				print_type = print_correct;
			}
	        	else { printf ("missing somethingdsaff\n"); }
		} else if (info_matrix->info[aa_atob[c]][l] == not_allowed) {
			if ( matrix->weights[aa_atob[c]][l] >= 0 )
			{ /* ERROR of OBSERVED, but INTOLERANT false positive */
				print_type = print_incorrect;
				pos_error_fp++; error_fp++;
			} else if ( matrix->weights[aa_atob[c]][l] < 0)
			{ /* PREDICTED CORRECTLY */
				notobserved_intolerant++;
				pos_notobserved_intolerant++;
				print_type = print_correct;
			} else {printf ("dfafaw\n"); }
		} else if (info_matrix->info[aa_atob[c]][l] == less_severe) {
			switch (severity_option) {
				case 1: /* include +- phenotype with + */
					if ( matrix->weights[aa_atob[c]][l] >= 0) {
						pos_observed_tolerant++;
		printf("less sever stuff\n");
						observed_tolerant++;
						print_type=print_correct;
					} else if (matrix->weights[aa_atob[c]][l] < 0){
						pos_error_tn++; error_tn++;
						 print_type = print_incorrect;
					printf ("incorrect less severe %d %d %c\n", l, (int) round (matrix->weights[aa_atob[c]][l]), c);
					 }
					break;
				case 2: /* +- phenotype is - */
					if ( matrix->weights[aa_atob[c]][l] >= 0) {
						pos_error_fp++;	error_fp++;
		printf ("incorrect less severe %d %d %c\n", l, matrix->weights[aa_atob[c]][l], c);
						print_type = print_incorrect;
					} else if (matrix->weights[aa_atob[c]][l] < 0){
					pos_notobserved_intolerant++;
			printf("less severe pos %d aa %c %d\n", l, c, round (matrix->weights[aa_atob[c]][l]));
					notobserved_intolerant++;
					print_type = print_correct;
					}
					break;
				default:
					fprintf (errorfp, "something wrong, severity option %d\n", severity_option);
					exit (-1);
			} /* end of switch (severity option) */

		} /* end of if == less severe */
	        break;
             case 2: /* for PSSMs with scores >= 0 i.e. counts
			scores == 0 predicted to be intolerant, > 0 tolerant*/
                        /* NO DATA */
                if (info_matrix->info[aa_atob[c]][l] == no_data ) {
                        print_type = print_X;
                } else if (info_matrix->info[aa_atob[c]][l] == allowed) {
                        if ( matrix->weights[aa_atob[c]][l] <= 0.0)
                        {    /* ERROR of NOT OBSERVED, BUT TOLERANT true negative */
				pos_error_tn++;
				error_tn++; print_type = print_incorrect;
                        } else if ( matrix->weights[aa_atob[c]][l] > 0.0)
                        { /* predicted correctly */
                                 pos_observed_tolerant++;
				observed_tolerant++; print_type = print_correct;
                        }
                } else if (info_matrix->info[aa_atob[c]][l] == not_allowed) {
                        if ( matrix->weights[aa_atob[c]][l] > 0.0 )
                        { /* ERROR of OBSERVED, but INTOLERANT false positive */
				print_type = print_incorrect;
                                pos_error_fp++; error_fp++;
                        } else if ( matrix->weights[aa_atob[c]][l] <= 0.0)
                        { /* PREDICTED CORRECTLY */
                                pos_notobserved_intolerant++;
				notobserved_intolerant++;
				print_type = print_correct;
                        }
                } else if (info_matrix->info[aa_atob[c]][l] == less_severe) {
                        switch (severity_option) {
                                case 1: /* include +- phenotype with + */
					if ( matrix->weights[aa_atob[c]][l] > 0.0) {
                                                pos_observed_tolerant++;
						observed_tolerant++;
						print_type=print_correct;
                                        } else if (matrix->weights[aa_atob[c]][l] <= 0.0){
                                                pos_error_tn++;
						error_tn++;
						 print_type = print_incorrect;
                                        } else {printf ("something is wrong\n");}
                                	break;
				case 2: /* +- phenotype is - */
                                        if ( matrix->weights[aa_atob[c]][l] > 0.0) {
                                               pos_error_fp++;
						error_fp++;
						 print_type = print_incorrect;
                                       	} else if (matrix->weights[aa_atob[c]][l] <= 0.0){
                                                pos_notobserved_intolerant++;
						notobserved_intolerant++;
						print_type = print_correct;
                                        }else {printf ("what?\n"); }
                                        break;
                        	default:
					fprintf(errorfp, "severity option error %d\n", severity_option);
					exit (-1);
				} /* end of switch (severity option) */

                }
                break;
	    default:
		fprintf (errorfp, "ERROR: Unknown style type %d\n", option);
		exit (-1);
	    } /* end of switch (option ) */
            if (print_type == print_correct) {
		if (style == INT_OUTPUT) {
			fprintf(omfp, "%3d ",
                                 (int) round(matrix->weights[aa_atob[c]][l]));
	    	} else if (style == FLOAT_OUTPUT) {
			fprintf(omfp, "%.1f ",
                                 matrix->weights[aa_atob[c]][l]);
		}
	    } else if (print_type == print_incorrect) {
		if (style == INT_OUTPUT) {
			fprintf(omfp, "%2d* ",
                                 (int) round(matrix->weights[aa_atob[c]][l]));
	        } else if (style == FLOAT_OUTPUT) {
			fprintf(omfp, "%.1f* ",
                                  matrix->weights[aa_atob[c]][l]);
		}
            } else if (print_type == print_X) {
	    	 fprintf(omfp, "  X ");
            }
	} /* end of if c != J, c!= O, ... */
         } /* end of for amino acids */
	pos_total = pos_observed_tolerant + pos_error_tn +
			pos_notobserved_intolerant + pos_error_fp;

         fprintf(omfp, "%.2f ", info_matrix->R[l]);
	fprintf (omfp, "%.2f", info_matrix->avg_R_per_residue[l]);
/*	fprintf (omfp, "%.2f ", (double) (pos_observed_tolerant + pos_notobserved_intolerant ) / (double) pos_total * 100); */
/*        fprintf (omfp, "%.2f ", (double) pos_notobserved_intolerant / (double) pos_total * 100); */
/*        fprintf (omfp, "%.2f ", (double) (pos_error_fp + pos_error_tn ) / (double) pos_total * 100); */
/*        fprintf (omfp, "%.2f ", (double) pos_error_tn / (double) pos_total * 100 ); */
/*         fprintf(omfp, "%3d\n", (int) round(matrix->weights[aa_atob['-']][l]));
*/
	fprintf (omfp, "%d ", pos_total);
	fprintf (omfp, "\n");
      } /* end of for entire matrix */

    total = notobserved_intolerant + observed_tolerant + error_fp + error_tn;
    fprintf(omfp, "//\n");
    fprintf (omfp, "correct: observed and tolerant    %5d\n", observed_tolerant);
    fprintf (omfp, "         not observed, intolerant %5d\n", notobserved_intolerant);
    fprintf (omfp, "incorrect:  observed, intolerant: %5d\n", error_fp);
    fprintf (omfp, "          not observed, tolerant: %5d\n", error_tn);
    fprintf (omfp, "\n");
    fprintf (omfp, "total counts : %d \n\n", total);
    fprintf (omfp, "In percentages -- \n");
    fprintf (omfp, "correct: observed and tolerant %3.1f\n", (double) observed_tolerant / (double) total * 100);
    fprintf (omfp, "      not observed, intolerant %3.1f\n", (double) notobserved_intolerant / (double) total * 100);
    fprintf (omfp, "incorrect: observed, intolerant: %3.1f\n", ((double) error_fp/ (double) total)*100);
    fprintf (omfp, "           not observed, tolerant: %3.1f\n", ((double)error_tn/(double) total)*100 );
    fprintf (omfp, "of tolerant residues, %3.1f predicted to be tolerant\n", (double)observed_tolerant/(double)
										(observed_tolerant + error_tn) * 100);
    fprintf (omfp, "of intolerant residues, %3.1f predicted to be intolerant\n", (double) notobserved_intolerant /
										(double) (notobserved_intolerant + error_fp) * 100);
   fprintf (omfp, "%3.1f of predicted intolerant residues are actually intolerant\n", (double) notobserved_intolerant/
										(double) (notobserved_intolerant + error_tn) * 100);

}  /*  end of output_matrix */

void output_Matrix_Info_float (Matrix_Info* info_matrix, FILE* omfp, int option, int severity_option)

{

  Matrix* matrix;
  int error_fp, error_tn, correct, total;
  int notobserved_intolerant, observed_tolerant;
  int pos_observed_tolerant, pos_error_tn, pos_notobserved_intolerant,
        pos_error_fp, pos_total;
  char c; int l;
  int print_correct, print_incorrect, print_X, print_type;

  print_correct = 1; print_incorrect = 2; print_X =3;
  correct = 0; error_fp=0; error_tn = 0, total=0;
  notobserved_intolerant = 0; observed_tolerant = 0;

  matrix = info_matrix->block_matrix;
  fprintf (omfp, "severity option %d\n", severity_option);

   for (c='A'; c <= 'Z'; c++) {
      if ((c != 'J') && (c != 'O') && (c != 'U') && (c != 'B') && (c != 'X')
                && (c != 'Z') )
        {
          fprintf(omfp, "  %c ", c);
        }
    } /* end of for loop,finished printing amino acids */
     fprintf(omfp, "    *         -\n");

  printf("pos\tpredicted tolerant correct\texper. tolerant\tpredicted intolerant corrent\texper. intolerant\n");

        for (l=0; l<matrix->width; l++) {
         /* score accuracy for position */
        pos_observed_tolerant = 0; pos_notobserved_intolerant = 0;
         pos_error_tn = 0; pos_error_fp = 0; pos_total = 0;
         for (c='A'; c <= 'Z'; c++) {
           if ((c != 'J') && (c != 'O') && (c != 'U') && (c != 'B') &&
                (c != 'X')  && (c != 'Z') ) {
           switch (option) {
             case 1: /* for PSSMs with negative scores.  predict scores <
                        as intolerant.  for log odds and matrices.
                        0 is neutral (therefore tolerant) */
                        /* NO DATA */
                if (info_matrix->info[aa_atob[c]][l] == no_data ) {
                        /* fprintf (omfp, "no data\n"); */
                        print_type = print_X;
                } else if (info_matrix->info[aa_atob[c]][l] == allowed) {
                        if ( matrix->weights[aa_atob[c]][l] < 0.01)
                        {    /* ERROR of NOT OBSERVED, BUT TOLERANT true negative */
                                error_tn++; print_type = print_incorrect;
                                pos_error_tn++;
                        } else if ( matrix->weights[aa_atob[c]][l] >= 0.01)
                        { /* predicted correctly */
                                 observed_tolerant++; pos_observed_tolerant++;
                                print_type = print_correct;
                        }
                        else { printf ("missing somethingdsaff\n"); }
                } else if (info_matrix->info[aa_atob[c]][l] == not_allowed) {
                        if ( matrix->weights[aa_atob[c]][l] >= 0.01 )
                        { /* ERROR of OBSERVED, but INTOLERANT false positive */
                                print_type = print_incorrect;
                                pos_error_fp++; error_fp++;
                        } else if ( matrix->weights[aa_atob[c]][l] < 0.01)
                        { /* PREDICTED CORRECTLY */
                                notobserved_intolerant++;
                                pos_notobserved_intolerant++;
                                print_type = print_correct;
                        } else {fprintf (errorfp, "dfafaw\n"); exit (-1); }
                } else if (info_matrix->info[aa_atob[c]][l] == less_severe) {
                        switch (severity_option) {
                                case 1: /* include +- phenotype with + */
                                        if ( matrix->weights[aa_atob[c]][l] >= 0.01) {
                                                pos_observed_tolerant++;
                                                observed_tolerant++;
                                                print_type=print_correct;
                                        } else if (matrix->weights[aa_atob[c]][l] < 0.01){
                                                pos_error_tn++; error_tn++;
                                                 print_type = print_incorrect;
                                         }
                                        break;
                                case 2: /* +- phenotype is - */
                                        if ( matrix->weights[aa_atob[c]][l] >= 0.01) {
                                                pos_error_fp++; error_fp++;
                                                print_type = print_incorrect;
                                        } else if (matrix->weights[aa_atob[c]][l] < 0.01){
                                        pos_notobserved_intolerant++;
                                        notobserved_intolerant++;
                                        print_type = print_correct;
                                        }
                                        break;
                                default:
                                        fprintf (errorfp, "something wrong, severity option %d\n", severity_option);
                                        exit (-1);
                        } /* end of switch (severity option) */

                } /* end of if == less severe */
                break;
             case 2: /* for PSSMs with scores >= 0 i.e. counts
                        scores == 0 predicted to be intolerant, > 0 tolerant*/
                        /* NO DATA */
                if (info_matrix->info[aa_atob[c]][l] == no_data ) {
                        print_type = print_X;
                } else if (info_matrix->info[aa_atob[c]][l] == allowed) {
                        if ( matrix->weights[aa_atob[c]][l] <= 0.01)
                        {    /* ERROR of NOT OBSERVED, BUT TOLERANT true negative */
                                pos_error_tn++;
                                error_tn++; print_type = print_incorrect;
                        } else if ( matrix->weights[aa_atob[c]][l] > 0.01)
                        { /* predicted correctly */
                                 pos_observed_tolerant++;
                                observed_tolerant++; print_type = print_correct;
                        }
                } else if (info_matrix->info[aa_atob[c]][l] == not_allowed) {
                        if ( matrix->weights[aa_atob[c]][l] > 0.01 )
                        { /* ERROR of OBSERVED, but INTOLERANT false positive */
                                print_type = print_incorrect;
                                pos_error_fp++; error_fp++;
                        } else if ( matrix->weights[aa_atob[c]][l] <= 0.01)
                        { /* PREDICTED CORRECTLY */
                                pos_notobserved_intolerant++;
                                notobserved_intolerant++;
                                print_type = print_correct;
                        }
                } else if (info_matrix->info[aa_atob[c]][l] == less_severe) {
                        switch (severity_option) {
                                case 1: /* include +- phenotype with + */
                                        if ( matrix->weights[aa_atob[c]][l] > 0.01) {
                                                pos_observed_tolerant++;
                                                observed_tolerant++;
                                                print_type=print_correct;
                                        } else if (matrix->weights[aa_atob[c]][l] <= 0.01){
                                                pos_error_tn++;
                                                error_tn++;
                                                 print_type = print_incorrect;
                                        } else {printf ("something is wrong\n");}
                                        break;
                                case 2: /* +- phenotype is - */
                                        if ( matrix->weights[aa_atob[c]][l] > 0.01) {
                                               pos_error_fp++;
                                                error_fp++;
                                                 print_type = print_incorrect;
                                        } else if (matrix->weights[aa_atob[c]][l] <= 0.01){
                                                pos_notobserved_intolerant++;
                                                notobserved_intolerant++;
                                                print_type = print_correct;
                                        }else {fprintf (errorfp, "what?\n"); exit (-1); }
                                        break;
                                default:
                                        fprintf(errorfp, "severity option error %d\n", severity_option);
                                        exit (-1);
                                } /* end of switch (severity option) */

                }
                break;
            default:
                fprintf (errorfp, "ERROR: Unknown style type %d\n", option);
                exit (-1);
            } /* end of switch (option ) */
            if (print_type == print_correct) {
                        fprintf(omfp, "%.2f ",
                                 matrix->weights[aa_atob[c]][l]);
            } else if (print_type == print_incorrect) {
                        fprintf(omfp, "%.2f* ",
                                  matrix->weights[aa_atob[c]][l]);
            } else if (print_type == print_X) {
                 fprintf(omfp, "  X ");
            }
        } /* end of if c != J, c!= O, ... */
         } /* end of for amino acids */
        pos_total = pos_observed_tolerant + pos_error_tn +
                        pos_notobserved_intolerant + pos_error_fp;
	if (matrix->block != NULL) {
		printf ("%d\t", l + matrix->block->sequences[0].position);
	}
	printf ("%d\t%d\t%d\t%d\n", pos_observed_tolerant, pos_error_tn + pos_observed_tolerant, pos_notobserved_intolerant, pos_notobserved_intolerant + pos_error_fp);

         fprintf(omfp, "%.2f ", info_matrix->R[l]);
        fprintf (omfp, "%.2f", info_matrix->avg_R_per_residue[l]);
        fprintf (omfp, "%d ", pos_total);
        fprintf (omfp, "\n");
      } /* end of for entire matrix */

    total = notobserved_intolerant + observed_tolerant + error_fp + error_tn;
    fprintf(omfp, "//\n");
    fprintf (omfp, "correct: observed and tolerant    %5d\n", observed_tolerant);
    fprintf (omfp, "         not observed, intolerant %5d\n", notobserved_intolerant);
    fprintf (omfp, "incorrect:  observed, intolerant: %5d\n", error_fp);
    fprintf (omfp, "          not observed, tolerant: %5d\n", error_tn);
    fprintf (omfp, "\n");
    fprintf (omfp, "total counts : %d \n\n", total);
    fprintf (omfp, "In percentages -- \n");
    fprintf (omfp, "correct: observed and tolerant %3.1f\n", (double) observed_tolerant / (double) total * 100);
    fprintf (omfp, "      not observed, intolerant %3.1f\n", (double) notobserved_intolerant / (double) total * 100);
    fprintf (omfp, "incorrect: observed, intolerant: %3.1f\n", ((double) error_fp/ (double) total)*100);
    fprintf (omfp, "           not observed, tolerant: %3.1f\n", ((double)error_tn/(double) total)*100);
   fprintf (omfp, "of tolerant residues, (%d/%d) ", observed_tolerant,
						observed_tolerant + error_tn);
   fprintf (omfp, " %3.1f predicted to be tolerant\n",
				(double)observed_tolerant/(double)
                                (observed_tolerant + error_tn) * 100);
    fprintf (omfp, "of intolerant residues, (%d/%d) ", notobserved_intolerant,
					notobserved_intolerant + error_fp);
    fprintf (omfp, "%3.1f predicted to be intolerant\n",
			(double) notobserved_intolerant /
                        (double) (notobserved_intolerant + error_fp) * 100);
   fprintf (omfp, "%3.1f (%d/%d) predicted intolerant residues are actually intolerant\n",
		(double) notobserved_intolerant/
                (double) (notobserved_intolerant + error_tn) * 100,
		notobserved_intolerant,
		notobserved_intolerant + error_tn);


}  /*  end of output_matrix */


void read_matrix_20_aa_body(mfp, matrix)
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
printf ("matrix width %d\n", matrix->width);
  /* skip the letters */
  fgets(Buffer, LARGE_BUFF_LENGTH, mfp);
printf ("skiipping letters %s\n", Buffer);
  /* read in the numbers */
  while (fgets(Buffer, LARGE_BUFF_LENGTH, mfp) &&
         !blank_line(Buffer) &&
         !((Buffer[0] == '/') && (Buffer[1] == '/'))) {
/*
    if (len > matrix->max_length) {
      resize_matrix(matrix);
    }
*/
    num_buf = get_token(Buffer);
    /*  NOTE:  Assumes matrix is in integer format !  */
    for (c='A'; c<='Z'; c++) {
      if ((c != 'J') && (c != 'O') && (c != 'U') && (c != 'B') && (c != 'X') &&
	  (c != 'Z')  ) {
        sscanf(num_buf, "%f", &dtemp);
        matrix->weights[aa_atob[c]][len] = (MatType) dtemp;
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

    len++;
  }
}  /* end of read_matrix_body */

void
normalize_matrix (Matrix* matrix)
{

	int pos, aa;
	double total;

	for (pos = 0; pos < matrix->width; pos++) {
		total = 0.0;
		for (aa = 0; aa < AAS; aa++) {
			total += matrix->weights[aa][pos];
		}
		if (total > 0.0) {
			for (aa = 0; aa < AAS; aa++) {
				matrix->weights[aa][pos] /= total;
			}
		}
	}
} /* end of normalize_matrix */

void
allow_min_and_above (Matrix* matrix)
{
	int min_aa, pos, aa;

	assert (matrix->block != NULL);

	for (pos = 0; pos < matrix->width; pos++) {
		min_aa = min_aa_in_column (matrix, pos);
		for (aa = 1; aa < AAS; aa++) {
			matrix->weights[aa][pos] -= min_aa;
		}
	}
}

int
min_aa_in_column (Matrix * matrix, int pos)
{
	double min;
	int min_aa, aa;

	min = 1000;
	for (aa = 0; aa < AAS; aa++) {
		if (matrix->weights[aa][pos] < min) {
			min_aa = aa;
			min = matrix->weights[aa][pos];
		}
	}
	return min_aa;
}

static void pssm_ratio_mutated_to_normal (Matrix* matrix, double threshold)
{
  int pos, normal_aa, aa;
  Block* block;
  double normal_aa_weight;

  assert (matrix->block != NULL);
  block = matrix->block;
  normalize_matrix (matrix);
  for (pos = 0; pos < matrix->width; pos++)
  {
    normal_aa = block->residues[0][pos];
    normal_aa_weight = matrix->weights[normal_aa][pos];
    /* calculate ratio of aa to normal aa observed in sequence */
    for (aa = 0; aa < AAS; aa++) {
	matrix->weights[aa][pos] /= normal_aa_weight;
        matrix->weights[aa][pos] -= threshold;
   }

  } /* end of for all positions */
}

void output_Matrix_Info_error_for_each_subst (
Matrix_Info* info_matrix, FILE* omfp, int severity_option)
{
	Matrix* matrix;
	int original_aa_index , new_aa_index;
	int total_array[20][20];
	int correct_array[20][20];
	double performance_array[20][20];
	int l, i, j;
	char c;
	matrix = info_matrix->block_matrix;

	for (i = 0; i < 20; i++) {
		for (j = 0; j < 20; j++) {
			correct_array[i][j] = 0;
			total_array[i][j] = 0;
		}
	}

	reclassify_Matrix_info (info_matrix, severity_option);


       for (l=0; l<matrix->width; l++) {
        printf ("%c original seq checkthi is a char\n", aa_btoa[matrix->block->sequences[0].sequence[l]]);

	original_aa_index = index_for_char( aa_btoa[matrix->block->sequences[0].sequence[l]]);

	for (c='A'; c <= 'Z'; c++) {
           if ((c != 'J') && (c != 'O') && (c != 'U') && (c != 'B') &&
                (c != 'X')  && (c != 'Z') ) {
		new_aa_index = index_for_char (c);
/*		if (info_matrix->info[aa_atob[c]][l] == allowed) { */

		if (info_matrix->info[aa_atob[c]][l] == not_allowed) {
			total_array[original_aa_index][new_aa_index]++;
			total_array[new_aa_index][original_aa_index]++;
/*			if ( matrix->weights[aa_atob[c]][l] >= THRESHOLD) { */
			if (matrix->weights[aa_atob[c]][l] < THRESHOLD) {
				correct_array[original_aa_index][new_aa_index]++;
				correct_array[new_aa_index][original_aa_index]++;
			}
		} /* end if allowd */
	  } /* end of if reall aa */
	} /* end of for loop through letters */
	} /* end of through entire matrix */
	for (i = 0; i < 20; i++) {
		for (j = 0; j < 20; j++) {
			if (total_array[i][j] > 0) {
				performance_array[i][j] = ((double) correct_array[i][j])/
						  ((double) total_array[i][j]);
			} else {
				performance_array[i][j] = -100;
			}
			fprintf (omfp, "%.2f  ", performance_array[i][j]);
		}
		fprintf (omfp, "\n");
	}

}

int
index_for_char (char c)
{
	switch (c) {
		case 'A': return 0;
		case 'C': return 1;
		case 'D': return 2;
		case 'E': return 3;
		case 'F': return 4;
		case 'G': return 5;
		case 'H': return 6;
		case 'I': return 7;
		case 'K': return 8;
		case 'L': return 9;
		case 'M': return 10;
		case 'N': return 11;
		case 'P': return 12;
		case 'Q': return 13;
		case 'R': return 14;
		case 'S': return 15;
		case 'T': return 16;
		case 'V': return 17;
		case 'W': return 18;
		case 'Y' : return 19;

		default: printf ("unknown character %c\n", c );
	} /* end switch */
}

void output_Matrix_Info_float_uncertaintybin (Matrix_Info* info_matrix,
						FILE* omfp,
						int option, int severity_option,
						int style,
						double upper_uncertainty_threshold,
						double lower_uncertainty_threshold)
{
 Matrix* matrix;
  int error_fp, error_tn, correct, total;
  int notobserved_intolerant, observed_tolerant;
  int pos_observed_tolerant, pos_error_tn, pos_notobserved_intolerant,
        pos_error_fp, pos_total;
  char c; int l;
  int print_correct, print_incorrect, print_X, print_type;
  int uncertain_aa, uncertain_tolerant, uncertain_intolerant;
  int pos_uncertain_aa, pos_uncertain_tolerant, pos_uncertain_intolerant;
  int total_correct, total_counted;

  reclassify_Matrix_info (info_matrix, severity_option);
  print_correct = 1; print_incorrect = 2; print_X =3;
  correct = 0; error_fp=0; error_tn = 0, total=0;
  notobserved_intolerant = 0; observed_tolerant = 0;
  uncertain_aa = 0; uncertain_tolerant = 0; uncertain_intolerant = 0;

  matrix = info_matrix->block_matrix;
  fprintf (omfp, "severity option %d\n", severity_option);

   for (c='A'; c <= 'Z'; c++) {
      if ((c != 'J') && (c != 'O') && (c != 'U') && (c != 'B') && (c != 'X')
                && (c != 'Z') )
        {
          fprintf(omfp, "  %c ", c);
        }
    } /* end of for loop,finished printing amino acids */
     fprintf(omfp, "    *         -\n");

  printf("pos\tpredicted tolerant correct\texper. tolerant\tpredicted intolerant correct\texper. intolerant\n");

        for (l=0; l<matrix->width; l++) {
         /* score accuracy for position */
        pos_observed_tolerant = 0; pos_notobserved_intolerant = 0;
         pos_error_tn = 0; pos_error_fp = 0; pos_total = 0;
         pos_uncertain_aa = 0; pos_uncertain_tolerant = 0; pos_uncertain_intolerant = 0;

	for (c='A'; c <= 'Z'; c++) {
           if ((c != 'J') && (c != 'O') && (c != 'U') && (c != 'B') &&
                (c != 'X')  && (c != 'Z') ) {
           switch (option) {
             case 1: /* for PSSMs with negative scores.  predict scores <
                        as intolerant.  for log odds and matrices.
                        0 is neutral (therefore tolerant) */
                        /* NO DATA */
                if (info_matrix->info[aa_atob[c]][l] == no_data ) {
                        /* fprintf (omfp, "no data\n"); */
                        print_type = print_X;
                } else if (matrix->weights[aa_atob[c]][l] < upper_uncertainty_threshold &&
			  matrix->weights[aa_atob[c]][l] > lower_uncertainty_threshold) {
			pos_uncertain_aa++; uncertain_aa++;
			if (info_matrix->info[aa_atob[c]][l] == allowed) {
				pos_uncertain_tolerant++; uncertain_tolerant++;
			} else if (info_matrix->info[aa_atob[c]][l] == not_allowed) {
				pos_uncertain_intolerant++; uncertain_intolerant++;
			} else { fprintf (errorfp, "spmething wrong inintolerant bin %c %d\n", c, l);
				 exit (-1);
			}
			/* UNCERTAINTY BIN */
		} else if (info_matrix->info[aa_atob[c]][l] == allowed) {
                        if ( matrix->weights[aa_atob[c]][l] < THRESHOLD)
                        {    /* ERROR of NOT OBSERVED, BUT TOLERANT true negative */
                                error_tn++; print_type = print_incorrect;
                                pos_error_tn++;
                        } else if ( matrix->weights[aa_atob[c]][l] >= THRESHOLD)
                       { /* predicted correctly */
                                 observed_tolerant++; pos_observed_tolerant++
;
                                print_type = print_correct;
                        }
                        else { printf ("missing somethingdsaff\n"); }
                } else if (info_matrix->info[aa_atob[c]][l] == not_allowed) {
                        if ( matrix->weights[aa_atob[c]][l] >= THRESHOLD )
                        { /* ERROR of OBSERVED, but INTOLERANT false positive */
                                print_type = print_incorrect;
                                pos_error_fp++; error_fp++;
                        } else if ( matrix->weights[aa_atob[c]][l] < THRESHOLD)
                        { /* PREDICTED CORRECTLY */
                                notobserved_intolerant++;
                                pos_notobserved_intolerant++;
                                print_type = print_correct;
                        } else {fprintf (errorfp, "dfafaw\n"); exit (-1); }
                } else { fprintf (errorfp, "something wrong at aa %c pos %d\n", c, l); exit (-1);}
                break;
             case 2: /* for PSSMs with scores >= 0 i.e. counts
                        scores == 0 predicted to be intolerant, > 0 tolerant */
                        /* NO DATA */
                if (info_matrix->info[aa_atob[c]][l] == no_data ) {
                        print_type = print_X;
                } else if (info_matrix->info[aa_atob[c]][l] == allowed) {
                        if ( matrix->weights[aa_atob[c]][l] <= 0.05)
                        {    /* ERROR of NOT OBSERVED, BUT TOLERANT true nega
tive */
                                pos_error_tn++;
                                error_tn++; print_type = print_incorrect;
                        } else if ( matrix->weights[aa_atob[c]][l] > 0.05)
                        { /* predicted correctly */
                                 pos_observed_tolerant++;
                                observed_tolerant++; print_type = print_correct;
                        }
                } else if (info_matrix->info[aa_atob[c]][l] == not_allowed) {
                        if ( matrix->weights[aa_atob[c]][l] > 0.05 )
                        { /* ERROR of OBSERVED, but INTOLERANT false positive
 */
                                print_type = print_incorrect;
                                pos_error_fp++; error_fp++;
                        } else if ( matrix->weights[aa_atob[c]][l] <= 0.05)
                        { /* PREDICTED CORRECTLY */
                                pos_notobserved_intolerant++;
                                notobserved_intolerant++;
                                print_type = print_correct;
                        }
                } else { fprintf (errorfp, "something very wrong at pos %d aa %c\n", l, c); exit (-1);
                }
                break;
            default:
                fprintf (errorfp, "ERROR: Unknown style type %d\n", option);
                exit (-1);
            } /* end of switch (option ) */
            if (print_type == print_correct) {
                        fprintf(omfp, "%.2f\t",
                                 matrix->weights[aa_atob[c]][l]);
            } else if (print_type == print_incorrect) {
                        fprintf(omfp, "%.2f*\t",
                                  matrix->weights[aa_atob[c]][l]);
            } else if (print_type == print_X) {
                 fprintf(omfp, "X\t");
            }
        } /* end of if c != J, c!= O, ... */
         } /* end of for amino acids */
        pos_total = pos_observed_tolerant + pos_error_tn +
                        pos_notobserved_intolerant + pos_error_fp;
        if (matrix->block != NULL) {
                printf ("%d\t", l + matrix->block->sequences[0].position);
        }
        printf ("%d\t%d\t%d\t%d\n", pos_observed_tolerant, pos_error_tn +
		pos_observed_tolerant, pos_notobserved_intolerant,
		pos_notobserved_intolerant +
		pos_error_fp);

         fprintf(omfp, "%.2f ", info_matrix->R[l]);
        fprintf (omfp, "%.2f", info_matrix->avg_R_per_residue[l]);
        fprintf (omfp, "%d ", pos_total);
        fprintf (omfp, "\n");
      } /* end of for entire matrix */

    total = notobserved_intolerant + observed_tolerant + error_fp + error_tn +
	    uncertain_tolerant + uncertain_intolerant;
    fprintf(omfp, "//\n");
    fprintf (omfp, "correct: observed and tolerant    %5d\n",
	observed_tolerant);
    fprintf (omfp, "         not observed, intolerant %5d\n",
					notobserved_intolerant);
    fprintf (omfp, "incorrect:  observed, intolerant: %5d\n", error_fp);
    fprintf (omfp, "          not observed, tolerant: %5d\n", error_tn);
    fprintf (omfp, "\n");
    fprintf (omfp, "total counts : %d \n\n", total);
    fprintf (omfp, "In percentages -- \n");
    fprintf (omfp, "correct: observed and tolerant %3.1f\n",
			(double) observed_tolerant / (double) total * 100);
    fprintf (omfp, "      not observed, intolerant %3.1f\n",
			(double) notobserved_intolerant / (double) total * 100);
    fprintf (omfp, "incorrect: observed, intolerant: %3.1f\n",
					((double) error_fp/ (double) total)*100);
    fprintf (omfp, "           not observed, tolerant: %3.1f\n",
					((double)error_tn/(double) total)*100);
   fprintf (omfp, "of tolerant residues, (%d/%d) ", observed_tolerant,
                                                observed_tolerant + error_tn);
   fprintf (omfp, " of the residues we chose to predict on ");
  fprintf (omfp, " %3.1f predicted to be tolerant\n\n",
                                (double)observed_tolerant/(double)
                                (observed_tolerant + error_tn) * 100);
    fprintf (omfp, "of intolerant residues, (%d/%d) ", notobserved_intolerant,
                                        notobserved_intolerant + error_fp);
    fprintf (omfp, "%3.1f predicted to be intolerant\n",
                        (double) notobserved_intolerant /
                        (double) (notobserved_intolerant + error_fp) * 100);
   fprintf (omfp, "%3.1f (%d/%d) predicted intolerant residues are actually intolerant\n",
                (double) notobserved_intolerant/
                (double) (notobserved_intolerant + error_tn) * 100,
                notobserved_intolerant,
                notobserved_intolerant + error_tn);
	total_correct = notobserved_intolerant + observed_tolerant;
	total_counted = notobserved_intolerant + error_fp +
			observed_tolerant + error_tn ;
   fprintf (omfp, "TOTAL PREDICTION ACCURACY %.3f (%d/%d)\n",
		 ((double) total_correct/ (double ) total_counted),
		total_correct, total_counted );

  fprintf (omfp, "UNCERTAINTY BIN*** (unclassified)\n");
  fprintf (omfp, "There were %d residues unclassified, %d of these are actually tolerated\n",
		uncertain_aa, uncertain_tolerant);
  fprintf (omfp, "while %d residues not tolerated \n", uncertain_intolerant);

}  /*  end of output_matrix */

void graph_scores_without_label (Matrix_Info* info_matrix)
{
        int l; char c;
        Matrix* matrix;

        matrix = info_matrix->block_matrix;

        printf ("POS\tTOLERATED\tLESS SEVERE\tMORE SEVERE\tINTOLERANT\n");
        for (l=0; l<matrix->width; l++) {
        for (c='A'; c <= 'Z'; c++) {
           if ((c != 'J') && (c != 'O') && (c != 'U') && (c != 'B') &&
                (c != 'X')  && (c != 'Z') ) {
                if (info_matrix->info[aa_atob[c]][l] != no_data ) {
                        if (info_matrix->info[aa_atob[c]][l] == allowed) {
                                printf ("%d\t%.2f\n", l, matrix->weights[aa_atob[c]][l]);
                        } else if (info_matrix->info[aa_atob[c]][l] ==
                                                        less_severe) {
                                printf ("%d\t\t%.2f\n", l, matrix->weights[aa_atob[c]][l]);
                        } else if (info_matrix->info[aa_atob[c]][l] ==
                                                        more_severe) {
                                printf ("%d\t\t\t%.2f\n", l,
                                                matrix->weights[aa_atob[c]][l]);
                        } else if (info_matrix->info[aa_atob[c]][l] ==
                                                not_allowed) {
                                printf ("%d\t\t\t\t%.2f\n", l, matrix->weights[aa_atob[c]][l]);
                        }
                }
            } /* end of if ! characters */
        } /* end of for characters */
        } /* end of for positions */
} /* end of graph_scores  */


/* prints median, mean, s.d. of info_matrix values based on phenotype */
/* 03/06/01 to be used for training SIFT v.2 */
void
print_stat_info_on_pssm (Matrix_Info* info_matrix)
{
	int pos, c;
	Matrix* matrix;
	/* for each position */
	double severe_at_pos[20]; int num_severe_at_pos;
	double less_severe_at_pos[20]; int num_less_severe_at_pos;
	double more_severe_at_pos[20]; int num_more_severe_at_pos;
	double tolerant_at_pos[20]; int num_tolerant_at_pos;

	/* over entire protein */
	/*double* severe;
	double* less_severe;
	double* more_severe;
	double* tolerant;
	int num_tolerant;
	int num_more_severe;
	int num_severe; int num_less_severe;
*/

	matrix =info_matrix->block_matrix;
/*	num_tolerant = 0;
	num_more_severe = 0;
	num_severe = 0; num_less_severe = 0;
	for (aa = 0 ; aa < 20; aa++) {
		severe[aa] = 0.0;
		less_severe[aa] = 0.0;
		more_severe[aa] = 0.0;
		tolerant[aa] = 0.0;
	}
*/

/*	severe = (double *) calloc (matrix->width, sizeof (double));
	less_severe = (double *) calloc (matrix->width, sizeof (double));
	more_severe = (double *) calloc (matrix->width, sizeof (double));
	tolerant = (double *) calloc (matrix->width, sizeof (double));
*/
	 for (pos = 0; pos < matrix->width; pos++) {
		num_tolerant_at_pos = 0;
		num_more_severe_at_pos = 0;
		num_severe_at_pos = 0;
		num_less_severe_at_pos = 0;

		for (c = 'A'; c<= 'Z'; c++) {
                   if (info_matrix->info[aa_atob[c]][pos] == allowed) {
                        tolerant_at_pos[num_tolerant_at_pos++] =
				matrix->weights[aa_atob[c]][pos];
		   } else if (info_matrix->info[aa_atob[c]][pos] ==
                                                        less_severe) {
                        less_severe_at_pos[num_less_severe_at_pos++] =
				matrix->weights[aa_atob[c]][pos];
                   } else if (info_matrix->info[aa_atob[c]][pos] ==
                                                        more_severe) {
                        more_severe_at_pos[num_more_severe_at_pos++] =
				matrix->weights[aa_atob[c]][pos];
                   } else if (info_matrix->info[aa_atob[c]][pos] ==
                                                not_allowed) {
			severe_at_pos[num_severe_at_pos++] =
				matrix->weights[aa_atob[c]][pos];
		  }
	        } /* end of for */

		/* finished looking at all amino acids at this position.
		now calculate median & mean */
		printf ("%d", pos+1);
		if (num_tolerant_at_pos != 0) {
			printf ("\t%.2f\t%.2f",
			mean_of_array (tolerant_at_pos, num_tolerant_at_pos),
			median_of_array(tolerant_at_pos, num_tolerant_at_pos)) ;
		} else {printf ("\tx\tx");}
		if (num_less_severe_at_pos != 0) {
		/*mean_of_array (less_severe_at_pos, num_less_severe_at_pos)*/;
		}
		if (num_more_severe_at_pos != 0) {
/*			mean_of_array (more_severe_at_pos, num_more_severe_at_pos); */
		}
		if (num_severe_at_pos != 0) {
			printf ("\t%.2f\t%.2f",
			mean_of_array(severe_at_pos, num_severe_at_pos),
			median_of_array (severe_at_pos, num_severe_at_pos));
		} else {printf ("\tx\tx"); }
		printf ("\n");
	} /* end of all pos */
}


void
graph_scores (Matrix_Info* info_matrix)
{
	int l; char c;
	Matrix* matrix;

	matrix = info_matrix->block_matrix;

	printf ("POS\tTOLERATED\tLESS SEVERE\tMORE SEVERE\tINTOLERANT\n");
        for (l=0; l<matrix->width; l++) {
        for (c='A'; c <= 'Z'; c++) {
           if ((c != 'J') && (c != 'O') && (c != 'U') && (c != 'B') &&
                (c != 'X')  && (c != 'Z') ) {
                if (info_matrix->info[aa_atob[c]][l] != no_data ) {
                        if (info_matrix->info[aa_atob[c]][l] == allowed) {
                                printf ("%d\t%.2f:T\n", l, matrix->weights[aa_atob[c]][l]);
                        } else if (info_matrix->info[aa_atob[c]][l] ==
							less_severe) {
                                printf ("%d\t\t%.2f:L\n", l, matrix->weights[aa_atob[c]][l]);
                        } else if (info_matrix->info[aa_atob[c]][l] ==
							more_severe) {
                                printf ("%d\t\t\t%.2f:S\n", l,
						matrix->weights[aa_atob[c]][l]);
			} else if (info_matrix->info[aa_atob[c]][l] ==
						not_allowed) {
                                printf ("%d\t\t\t\t%.2f:I\n", l, matrix->weights[aa_atob[c]][l]);
                        }
		}
	    } /* end of if ! characters */
	} /* end of for characters */
	} /* end of for positions */
} /* end of graph_scores  */

double
R_at_pos (Matrix* matrix, int pos)
{
	int seq, aa, num_residues;
	double ln2, e, r, dtemp, totalweight;

	num_residues = matrix->num_sequences ;
	ln2 = log(2.0);

	if (matrix->block != NULL) {
        	for(seq=0; seq < matrix->num_sequences; seq++)
            		if (matrix->block->residues[seq][pos] < 1 ||
                	matrix->block->residues[seq][pos] > 22) num_residues-- ;
  	}
                                /* correction factor for small sample size */
        e = (AAs-1) / ((double) 2 * ln2 * num_residues) ;

        /* =====end of calculating error correction ================*/
        assert (matrix->weights != NULL);

        r = log ((double) AAs);  /* R = log 20 */
        for(aa=1, totalweight=0.; aa < AAs+1; aa++) {
           totalweight +=  matrix->weights[aa][pos]  ;
        }
        /* totalweight is the sum of all aa weights at a certain pos */
        for(aa=1; aa < AAs+1; aa++)
                    /* start loop at 1 and end at AAs+1 because aa values
                      at the weight matrix start at 1, 0 is the gap value */
       {
         if (matrix->weights[aa][pos] > 0)
         {
            dtemp = (double) matrix->weights[aa][pos] / totalweight;
            r += dtemp * log(dtemp);
         } /* end of if weights > 0 */
        } /* end of for for calculating r*/

      r /= ln2 ;       /* convert to bits , R = 2 - H(pos) in bits now */

/* if (r > 4) { printf ("pos %d exceeded over 4 bits\n", pos); } */

      r -= e ;                     /* a correction for small sample sizes */

	return r;
}

void
free_AAnode (AAnode aanode)
{
	AAnode current;
	AAnode tmp;


	current = aanode;
	if (current != NULL) {

	while (current->next != NULL) {
		tmp = current->next;
		free (current);
		current = tmp;
	}
	free (current);
	} /* end check that current != NULL */
}

void
print_multiple_polymorphism_data (AAnode* aanodes, int length)
{
	int i;
	AAnode current;

	for (i=0; i < length; i++) {
		current = aanodes[i]->next; /* nothing in first node */
		while (current != NULL)
		{
			printf ("%d%c\n", i, aa_btoa[current->aa]);
			current = current->next;
		}
	}
}


AAnode*
read_multiple_polymorphism_data (FILE* infofp, Sequence* seq)
{
       int aa, i, aa_pos, pos;
        char line[LARGE_BUFF_LENGTH];
        char* stringp;
        char *string_pos;
        char original_aa, substituted_aa;
        char word1[10], word2[10];
        AAnode* polymorph_data;
        int warn_aasubst_not_matching_with_query ;
	AAnode newAAnode; AAnode tmp_AAnode;

        warn_aasubst_not_matching_with_query  = 0;
        polymorph_data = (AAnode*) calloc (seq->length, sizeof (AAnode));
        for (i = 0; i < seq->length; i++) {
                polymorph_data[i] = (AAnode) malloc (sizeof (struct aanode));
		polymorph_data[i]->aa = aa_atob['-'];
        	polymorph_data[i]->next = NULL;
	}
	if (infofp == NULL) {
		// fprintf (errorfp, "no substitutions read\n");
		return polymorph_data;
	} /* no file to read*/
	// printf ("before seg fault?\n");
	while (fgets (line, LARGE_BUFF_LENGTH, infofp) != NULL) {
                word1[0] = '\0'; word2[0] = '\0';
                stringp = strtok (line, " \r\n\t");
                if (stringp != NULL && stringp[0] != '#') {
                                /* not an empty line or a comment*/
                         strcpy (word1, stringp);
        /* parse original amino acid, position, and new amino acid */
        /* Example of what file should read: Y2P => tyrosine in original
        sequence at position 2 substituted for Proline */
                original_aa = word1[0];
                strcpy (word2, &word1[1]);
                pos = atoi (word2);
                pos--; /* location in array */
                stringp = strpbrk (word2, "ACDEFGHIKLMNPQRSTVWY");
                substituted_aa = *stringp;

                if ( pos >= 0  && original_aa != aa_btoa[seq->sequence[pos]] ) {
                        warn_aasubst_not_matching_with_query =1;
                        fprintf (errorfp, "<font color=red>ERROR!!  ");
                        fprintf (errorfp, "Substitution entered says that pos %d is %c",pos +1, original_aa);
                        fprintf (errorfp, " but the sequence file has %c at this position</font><BR><BR>\n", aa_btoa[seq->sequence[pos]]);
                }
		tmp_AAnode = polymorph_data[pos];
		while (tmp_AAnode->next != NULL) {
			tmp_AAnode = tmp_AAnode->next;
		}
		newAAnode = (AAnode) malloc (sizeof (struct aanode));
		newAAnode->next = NULL;
                newAAnode->aa = aa_atob[substituted_aa];
		tmp_AAnode->next = newAAnode;

                } /* end of if nothng on line */
        } /* finished reading file */

        if (warn_aasubst_not_matching_with_query == 1) {
                fprintf (errorfp, "<font color=red>*****Please check that you have entered your substitutions and the sequence correctly.*****</font><BR><BR>\n");
        }

        return polymorph_data;


} /* end read_multiple_polymprhism_data*/

Residue*
read_polymorphism_data (FILE* infofp, Sequence* seq)
{
        int aa, i, aa_pos, pos;
        char line[LARGE_BUFF_LENGTH];
        char* stringp;
        char *string_pos;
        char original_aa, substituted_aa;
        char word1[10], word2[10];
        Residue* polymorph_data;
	int warn_aasubst_not_matching_with_query ;

	warn_aasubst_not_matching_with_query  = 0;

        polymorph_data = (Residue*) calloc (seq->length, sizeof (Residue));
        for (i = 0; i < seq->length; i++) {
                polymorph_data[i] = aa_atob['-'];
        }
        if (infofp == NULL) { return polymorph_data;} /* no file to read*/
        while (fgets (line, LARGE_BUFF_LENGTH, infofp) != NULL) {
                word1[0] = '\0'; word2[0] = '\0';
                stringp = strtok (line, " \r\n\t");
                if (stringp != NULL && stringp[0] != '#') {
				/* not an empty line or a comment*/
                         strcpy (word1, stringp);
        /* parse original amino acid, position, and new amino acid */
        /* Example of what file should read: Y2P => tyrosine in original
        sequence at position 2 substituted for Proline */
                original_aa = word1[0];
                strcpy (word2, &word1[1]);
                pos = atoi (word2);
                pos--; /* location in array */
                stringp = strpbrk (word2, "ACDEFGHIKLMNPQRSTVWY*");
                substituted_aa = *stringp;

                if (pos >= 0 && original_aa != aa_btoa[seq->sequence[pos]] ) {
			warn_aasubst_not_matching_with_query =1;
                        fprintf (stderr, "<font color=red>ERROR!!  ");
                        fprintf (stderr, "Substitution entered says that pos %d is %c",pos +1, original_aa);
			fprintf (stderr, " but the sequence file has %c at this position</font><BR><BR>\n", aa_btoa[seq->sequence[pos]]);
                }
		if (substituted_aa != '*') {
                polymorph_data[pos] = aa_atob[substituted_aa];
                } /* don't record stop codons */
		} /* end of if nothng on line */
        } /* finished reading file */

	if (warn_aasubst_not_matching_with_query == 1) {
                fprintf (stderr, "<font color=red>*****Please check that you have entered your substitutions and the sequence correctly.*****</font><BR><BR>\n");
        }

        return polymorph_data;

} /* end of read_polymorphism_data */

void
remove_data_for_low_fraction_stored_pos (Matrix_Info* matrix_info, double* fraction_stored)
{
	int pos, aa;

	for (pos = 0;pos < matrix_info->block_matrix->width; pos++) {
		if (fraction_stored[pos] < 0.5) {
			for (aa =0; aa < MATRIX_AA_WIDTH; aa++) {
				matrix_info->info[aa][pos] = no_data;
			}
		}
	}
} /* end of remove data */
#endif
