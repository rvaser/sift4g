/* (C) Copyright 1993-2000, Fred Hutchinson Cancer Research Center */
/* Use, modification or distribution of these programs is subject to */
/* the terms of the non-commercial licensing agreement in license.h. */

/* blimps.c: main module of the BLIMPS blocks searcher.  */
/* Written by: Bill Alford */
/* Change log information is at the end of the file. */

#define EXTERN

/*	system headers not in global.h  */
#include <math.h>
#include <signal.h>

/*	blimps headers */
#include <global.h>
#include <version.h>
#include <residues.h>
#include <files.h>
#include <blocks.h>	/* includes sequences.h, output.h */
#include <matrix.h>	/* includes pattern.h */
#include <convert.h>
#include <frequency.h>
#include <options.h>

/*	headers in current directory */
#include "blimps.h"
#include "blimps-mem.h"
#include "config.h"
#include "scores.h"
#include "scoring.h"		/* for Alignments_Done and Scores_Done */
#include "lists.h"



/*
 * Exported variables and data structures
 */

/* variables set by the configuration file */

int StrandsToSearch;
int NumberToReport;
int SearchType;
int SequenceType;
Boolean RepeatsAllowed;
Boolean SavedScoresFlag;
int GeneticCodeInitializer;
int SiteSpecificScoringMatrixType;
int BlockToMatrixConversionMethod; /* default method is two */
int SequenceMatrixScoringMethod;   /* default method is zero */



/*
 * Local variables and data structures
 */


static Block *block;
static Matrix *matrix;
static Sequence *sequence;
static Sequence *trans1, *trans2, *trans3;
				/* translated DNA reading frames */
static Sequence *trans_1, *trans_2, *trans_3;
				/* translated DNA reading frames */

static unsigned char gcode[64], revgcode[64];
				/* genetic codes for translation */

static int records_searched;


/*
 * sequence_vs_blocks
 *   When SearchType is BLOCK, there is one sequence against a blocks
 *   database.
 *   Parameters:
 *     FILE *emfp: the file to output the matrix information to.
 *   Error codes:
 */

static void sequence_vs_blocks(emfp)
     FILE *emfp;
{
  char prev_number[20];

  /* make sure we are at the beginning of the sequence data */
  rewind_file(SEQUENCE_FILES);

  /* read in the first sequence (if there are others, won't be processed */
  sequence = read_a_sequence(get_file(SEQUENCE_FILES),
			     get_sequence_db_db_type(),
			     get_sequence_db_seq_type());

 if (sequence == NULL) {
    sprintf(ErrorBuffer,
	    "Unable to read the sequence for scoring against the blocks.");
    ErrorReport(SERIOUS_ERR_LVL);
    return;
  }
 if (sequence->length <= 0) {
    sprintf(ErrorBuffer,
	    "Query sequence has zero length, not scoring it.\n");
    ErrorReport(SERIOUS_ERR_LVL);
    return;
  }

  /*   Just translate the query sequence once  */
  /*   May take too much memory for long DNA sequences ... */
  trans1 = trans2 = trans3 = trans_1 = trans_2 = trans_3 = NULL;
  if (sequence->type == NA_SEQ)
  {
     if (StrandsToSearch > 0)
     {
        trans1 = translate_sequence(sequence, 1, gcode, revgcode);
        trans2 = translate_sequence(sequence, 2, gcode, revgcode);
        trans3 = translate_sequence(sequence, 3, gcode, revgcode);
     }
     if (StrandsToSearch < 0 || StrandsToSearch == 2)
     {
        trans_1 = translate_sequence(sequence, -1, gcode, revgcode);
        trans_2 = translate_sequence(sequence, -2, gcode, revgcode);
        trans_3 = translate_sequence(sequence, -3, gcode, revgcode);
     }
  }

  if (NumberToReport == 0) {
    NumberToReport = DEFAULT_BLOCK_SEARCH_NUMBER_TO_REPORT(sequence);
  }

  /* inform what sequence the search is being done with */
/*	Trying to speed things up ...
  sprintf(ErrorBuffer,
	  "Searching sequence %s against the blocks database.\n",
	  sequence->name);
  ErrorReport(INFO_ERR_LVL);
*/

    records_searched = 0;

    /* make sure we are at the begining of the block data */
    rewind_file(BLOCK_FILES);

    /* ------------------------loop through all the blocks */
    prev_number[0] = '\0';
    /*  Assuming the database is correctly formatted  */
    while ((block = read_a_block_faster(get_file(BLOCK_FILES))) != NULL)
    {
        /* inform which block scoring against */
/*  Trying to speed things up ...
        sprintf(ErrorBuffer, "Scoring vs block %s.\n", block->number);
	ErrorReport(INFO_ERR_LVL);
*/
/*		Check for new family
	if ((int) strlen(prev_number) < 7 ||
                  strncmp(block->number, prev_number, 7) != 0)
        {
		strcpy(prev_number, block->number);
        }
*/
        /* Compute a PSSM from the block */
	matrix = block_to_matrix(block, BlockToMatrixConversionMethod);
	if (emfp != NULL) { output_matrix(matrix, emfp); }


	/* score and enter the data into the list */
/*>>>>NOTE: just call default_scoring_method() instead of
	score_and_enter() since no other methods are implemented
 <<<<<< */

        if (sequence->type == NA_SEQ)
        {
           if (StrandsToSearch > 0)
           {
/*
	      score_and_enter(trans1, matrix, 1, RepeatsAllowed, SearchType);
	      score_and_enter(trans2, matrix, 2, RepeatsAllowed, SearchType);
	      score_and_enter(trans3, matrix, 3, RepeatsAllowed, SearchType);
*/
              default_scoring_method(trans1, matrix, 1, RepeatsAllowed, SearchType, FALSE);
              default_scoring_method(trans2, matrix, 2, RepeatsAllowed, SearchType, FALSE);
              default_scoring_method(trans3, matrix, 3, RepeatsAllowed, SearchType, FALSE);
           }
           if (StrandsToSearch < 0 || StrandsToSearch == 2)
           {
/*
	      score_and_enter(trans_1, matrix, -1, RepeatsAllowed, SearchType);
	      score_and_enter(trans_2, matrix, -2, RepeatsAllowed, SearchType);
	      score_and_enter(trans_3, matrix, -3, RepeatsAllowed, SearchType);
*/
              default_scoring_method(trans_1, matrix, -1, RepeatsAllowed, SearchType, FALSE);
              default_scoring_method(trans_2, matrix, -2, RepeatsAllowed, SearchType, FALSE);
              default_scoring_method(trans_3, matrix, -3, RepeatsAllowed, SearchType, FALSE);
           }
        }
        else  /* AA  */
        {
/*
           score_and_enter(sequence, matrix, 0, RepeatsAllowed, SearchType);
*/
           default_scoring_method(sequence, matrix, 0, RepeatsAllowed, SearchType, FALSE);
        }

	records_searched++;
	free_block(block);
	free_matrix(matrix);
    } /* end while there are blocks */

    if (sequence->type == NA_SEQ)
    {
       if (trans1 != NULL) free_sequence(trans1);
       if (trans2 != NULL) free_sequence(trans2);
       if (trans3 != NULL) free_sequence(trans3);
       if (trans_1 != NULL) free_sequence(trans_1);
       if (trans_2 != NULL) free_sequence(trans_2);
       if (trans_3 != NULL) free_sequence(trans_3);
    }

/*  free(sequence);  This is a global var, used later !!! */

}  /* end of sequence_vs_blocks */


/*
 * block_vs_sequences
 *   When SearchType is MATRIX, there is one block against a sequences
 *   database.
 *   Parameters:
 *     FILE *emfp: the file to output the matrix information to.
 *   Error codes:
 */

static void block_vs_sequences(emfp)
     FILE *emfp;
{
  /* make sure we are at the begining of the block data */
  rewind_file(BLOCK_FILES);

  /* loop through all the blocks (one) */
  block = read_a_block(get_file(BLOCK_FILES));

  if (block == NULL) {
    sprintf(ErrorBuffer,
	    "Unable to read block for scoring against the sequences.");
    ErrorReport(SERIOUS_ERR_LVL);
    sprintf(ErrorBuffer,
	    "Not scoring the block.\n");
    ErrorReport(SERIOUS_ERR_LVL);
    return;
  }

  if (NumberToReport == 0) {
    NumberToReport = DEFAULT_MATRIX_SEARCH_NUMBER_TO_REPORT(block);
  }

  /* inform what block the search is being done with */
  sprintf(ErrorBuffer,
	  "Searching block %s against the sequence database.\n",
	  block->number);
  ErrorReport(INFO_ERR_LVL);

  /* get the block matrix */
  matrix = block_to_matrix(block, BlockToMatrixConversionMethod);
  if (emfp != NULL) {
    output_matrix(matrix, emfp);
  }

  /* make sure we are at the beginning of the sequence data */
  rewind_file(SEQUENCE_FILES);

  records_searched = 0;

  /* loop through all the sequences */
  while ((sequence = read_a_sequence(get_file(SEQUENCE_FILES),
				     get_sequence_db_db_type(),
				     get_sequence_db_seq_type())) != NULL) {

    /* inform which sequence scoring against */
    sprintf(ErrorBuffer,
	    "Scoring vs sequence %s.\n",
	    sequence->name);
    ErrorReport(INFO_ERR_LVL);

    /* is the sequence a NA_SEQ?  if so, translate */
    if (sequence->type == NA_SEQ) {
      /*   Search forward direction  */
      if (StrandsToSearch > 0)
      {
         /* translate */
         trans1 = translate_sequence(sequence, 1, gcode, revgcode);
         /* score and enter the data into the list */
         score_and_enter(trans1, matrix, 1, RepeatsAllowed, SearchType);
         /* free the translated sequences */
         free_sequence(trans1);

         trans2 = translate_sequence(sequence, 2, gcode, revgcode);
         score_and_enter(trans2, matrix, 2, RepeatsAllowed, SearchType);
         free_sequence(trans2);

         trans3 = translate_sequence(sequence, 3, gcode, revgcode);
         score_and_enter(trans3, matrix, 3, RepeatsAllowed, SearchType);
         free_sequence(trans3);
      }
      /*   Search reverse direction  */
      if (StrandsToSearch == 2 || StrandsToSearch < 0)
      {
	 trans_1 = translate_sequence(sequence, -1, gcode, revgcode);
	 score_and_enter(trans_1, matrix, -1, RepeatsAllowed, SearchType);
	 free_sequence(trans_1);

	 trans_2 = translate_sequence(sequence, -2, gcode, revgcode);
	 score_and_enter(trans_2, matrix, -2, RepeatsAllowed, SearchType);
	 free_sequence(trans_2);

	 trans_3 = translate_sequence(sequence, -3, gcode, revgcode);
	 score_and_enter(trans_3, matrix, -3, RepeatsAllowed, SearchType);
	 free_sequence(trans_3);
      }
    }  /* end of NA_SEQ  */
    /* otherwise, just use the sequence, assuming it is an AA_SEQ */
    else {
      /* score and enter the data into the list */
      score_and_enter(sequence, matrix, 0, RepeatsAllowed, SearchType);
    }

    free_sequence(sequence);

    records_searched++;

  } /* end while there are sequences */

  free_block(block);
  free_matrix(matrix);
}  /* end of block_vs_sequences */


/*
 * blocks_vs_sequences
 *   When SearchType is UNKNOWN, there is probably a blocks database against
 *   a sequences database.  Score all pairwise matches of blocks and
 *   sequences.
 *   Parameters:
 *     FILE *emfp: the file to output the matrix information to.
 *   Error codes:
 */

static void blocks_vs_sequences(emfp)
     FILE *emfp;
{
  /* NOTE: If there are speed problems, might look at the SearchType */
  /*       variable in files.c to see which of sequences or blocks is */
  /*       considered the database.  The one that is the database is most */
  /*       likely to have the most elements, and should be in the outer */
  /*       loop to reduce "block->matrix" and "if NA_SEQ; NA_SEQ->AA_SEQ" */
  /*       overhead. */

  /* set the number to report */
  if (NumberToReport == 0) {
    rewind_file(BLOCK_FILES);
    block = read_a_block(get_file(BLOCK_FILES));
    rewind_file(SEQUENCE_FILES);
    sequence = read_a_sequence(get_file(SEQUENCE_FILES),
			       get_sequence_db_db_type(),
			       get_sequence_db_seq_type());
    if (block == NULL) {
      if (sequence == NULL) {
	/* there are no comparisons to do */
	NumberToReport = 1;
      }
      else {
	NumberToReport = DEFAULT_BLOCK_SEARCH_NUMBER_TO_REPORT(sequence);
	free_sequence(sequence);
      }
    }
    else {
      if (sequence == NULL) {
	/* there are no comparisons to do */
	NumberToReport = DEFAULT_MATRIX_SEARCH_NUMBER_TO_REPORT(block);
	free_block(block);
      }
      else {
	NumberToReport =
	  max((DEFAULT_BLOCK_SEARCH_NUMBER_TO_REPORT(sequence)),
	      (DEFAULT_MATRIX_SEARCH_NUMBER_TO_REPORT(block)));
	free_block(block);
	free_sequence(sequence);
      }
    }
    rewind_file(BLOCK_FILES);
    rewind_file(SEQUENCE_FILES);
  }

  /* make sure we are at the begining of the block data */
  rewind_file(BLOCK_FILES);

  /* loop through all the blocks */
  while ((block = read_a_block(get_file(BLOCK_FILES))) != NULL) {

    /* get the block matrix */
    matrix = block_to_matrix(block, BlockToMatrixConversionMethod);
    if (emfp != NULL) {
      output_matrix_s(matrix, emfp, FLOAT_OUTPUT);
    }

    /* make sure we are at the beginning of the sequence data */
    rewind_file(SEQUENCE_FILES);

    /* loop through all the sequences */
    while ((sequence = read_a_sequence(get_file(SEQUENCE_FILES),
				       get_sequence_db_db_type(),
				       get_sequence_db_seq_type())) != NULL) {

      /* inform which block and sequence are being scored */
      sprintf(ErrorBuffer,
	      "Scoring sequence %s vs block %s.\n",
	      sequence->name, block->number);
      ErrorReport(INFO_ERR_LVL);

      /* is the sequence a NA_SEQ?  if so, translate */
      if (sequence->type == NA_SEQ) {
	/* translate */
	trans1 = translate_sequence(sequence, 1, gcode, revgcode);
	trans2 = translate_sequence(sequence, 2, gcode, revgcode);
	trans3 = translate_sequence(sequence, 3, gcode, revgcode);
	/* score and enter the data into the list */
	score_and_enter(trans1, matrix, 1, RepeatsAllowed, SearchType);
	score_and_enter(trans2, matrix, 2, RepeatsAllowed, SearchType);
	score_and_enter(trans3, matrix, 3, RepeatsAllowed, SearchType);
	/* free the translated sequences */
	free_sequence(trans1);
	free_sequence(trans2);
	free_sequence(trans3);
	/* if we are supposed to search both strands, translate and score */
	if (StrandsToSearch == 2) {
	  /* translate */
	  trans_1 = translate_sequence(sequence, -1, gcode, revgcode);
	  trans_2 = translate_sequence(sequence, -2, gcode, revgcode);
	  trans_3 = translate_sequence(sequence, -3, gcode, revgcode);
	  /* score and enter the data into the list */
	  score_and_enter(trans_1, matrix, -1, RepeatsAllowed, SearchType);
	  score_and_enter(trans_2, matrix, -2, RepeatsAllowed, SearchType);
	  score_and_enter(trans_3, matrix, -3, RepeatsAllowed, SearchType);
	  /* free the translated sequences */
	  free_sequence(trans_1);
	  free_sequence(trans_2);
	  free_sequence(trans_3);
	}
      }
      /* otherwise, just use the sequence, assuming it is an AA_SEQ */
      else {
	/* score and enter the data into the list */
	score_and_enter(sequence, matrix, 0, RepeatsAllowed, SearchType);
      }

      free_sequence(sequence);

      records_searched++;

    } /* end while there are sequences */

    free_block(block);
    free_matrix(matrix);

  } /* end while there are blocks */
}  /* end of blocks_vs_sequences */



/*
 * sequence_vs_matrices
 *   When SearchType is BLOCK, there is one sequence against a matrix
 *   database.
 *   Parameters:
 *     FILE *emfp: the file to output the matrix information to.
 *   Error codes:
 */

static void sequence_vs_matrices(emfp)
     FILE *emfp;
{
  /* make sure we are at the beginning of the sequence data */
  rewind_file(SEQUENCE_FILES);

  /* loop through all the sequences (one) */
  sequence = read_a_sequence(get_file(SEQUENCE_FILES),
			     get_sequence_db_db_type(),
			     get_sequence_db_seq_type());

 if (sequence == NULL) {
    sprintf(ErrorBuffer,
	    "Unable to read the sequence for scoring against the matrices.");
    ErrorReport(SERIOUS_ERR_LVL);
    sprintf(ErrorBuffer,
	    "Not scoring the sequence.\n");
    ErrorReport(SERIOUS_ERR_LVL);
    return;
  }

  if (NumberToReport == 0) {
    NumberToReport = DEFAULT_BLOCK_SEARCH_NUMBER_TO_REPORT(sequence);
  }

  /* inform what sequence the search is being done with */
  sprintf(ErrorBuffer,
	  "Searching sequence %s against the matrices database.\n",
	  sequence->name);
  ErrorReport(INFO_ERR_LVL);

  records_searched = 0;

  /* is the sequence a NA_SEQ?  if so, translate and score */
  if (sequence->type == NA_SEQ) {

    /* translate */
    trans1 = translate_sequence(sequence, 1, gcode, revgcode);
    trans2 = translate_sequence(sequence, 2, gcode, revgcode);
    trans3 = translate_sequence(sequence, 3, gcode, revgcode);

    /* if we are supposed to search both strands, translate and score */
    if (StrandsToSearch == 2) {

      /* translate */
      trans_1 = translate_sequence(sequence, -1, gcode, revgcode);
      trans_2 = translate_sequence(sequence, -2, gcode, revgcode);
      trans_3 = translate_sequence(sequence, -3, gcode, revgcode);

    } /* end if StrandsToSearch == 2 */

    /* make sure we are at the begining of the matrix data */
    rewind_file(MATRIX_FILES);

    if (StrandsToSearch != 2) {
      /* loop through all the matrices */
      while ((matrix = read_a_matrix(get_file(MATRIX_FILES))) != NULL) {
	/* inform which matrix scoring against */
	sprintf(ErrorBuffer,
		"Scoring vs matrix %s.\n",
		matrix->number);
	ErrorReport(INFO_ERR_LVL);
	if (emfp != NULL) {
	  output_matrix(matrix, emfp);
	}
	/* score and enter the data into the list */
	score_and_enter(trans1, matrix, 1, RepeatsAllowed, SearchType);
	score_and_enter(trans2, matrix, 2, RepeatsAllowed, SearchType);
	score_and_enter(trans3, matrix, 3, RepeatsAllowed, SearchType);

	records_searched++;

	free_matrix(matrix);
      } /* end while there are matrices */
    }
    else { /* StrandsToSearch == 2 */
      /* loop through all the matrices */
      while ((matrix = read_a_matrix(get_file(MATRIX_FILES))) != NULL) {
	/* inform which matrix scoring against */
	sprintf(ErrorBuffer,
		"Scoring vs matrix %s.\n",
		matrix->number);
	ErrorReport(INFO_ERR_LVL);
	if (emfp != NULL) {
	  output_matrix(matrix, emfp);
	}
	/* score and enter the data into the list */
	score_and_enter(trans1, matrix, 1, RepeatsAllowed, SearchType);
	score_and_enter(trans2, matrix, 2, RepeatsAllowed, SearchType);
	score_and_enter(trans3, matrix, 3, RepeatsAllowed, SearchType);
	score_and_enter(trans_1, matrix, -1, RepeatsAllowed, SearchType);
	score_and_enter(trans_2, matrix, -2, RepeatsAllowed, SearchType);
	score_and_enter(trans_3, matrix, -3, RepeatsAllowed, SearchType);

	records_searched++;

	free_matrix(matrix);
      } /* end while there are matrices */
    }

  } /* end if NA_SEQ */
  /* otherwise, just use the sequence, assuming it is an AA_SEQ */
  else {
    /* make sure we are at the begining of the matrix data */
    rewind_file(MATRIX_FILES);

    /* loop through all the matrices */
    while ((matrix = read_a_matrix(get_file(MATRIX_FILES))) != NULL) {
      /* inform which matrix scoring against */
      sprintf(ErrorBuffer,
	      "Scoring vs matrix %s.\n",
	      matrix->number);
      ErrorReport(INFO_ERR_LVL);
      if (emfp != NULL) {
	output_matrix(matrix, emfp);
      }
      /* score and enter the data into the list */
      score_and_enter(sequence, matrix, 0, RepeatsAllowed, SearchType);

      records_searched++;

      free_matrix(matrix);
    } /* end while there are matrices */
  } /* end else assuming an AA_SEQ */

/*  free(sequence); This is a global var. & is used later !!*/

}  /* end of sequence_vs_matrices */


/*
 * matrix_vs_sequences
 *   When SearchType is MATRIX, there is one matrix against a sequences
 *   database.
 *   Parameters:
 *     FILE *emfp: the file to output the matrix information to.
 *   Error codes:
 */

static void matrix_vs_sequences(emfp)
     FILE *emfp;
{
  /* make sure we are at the begining of the matrix data */
  rewind_file(MATRIX_FILES);

  /* loop through all the matrices (one) */
  matrix = read_a_matrix(get_file(MATRIX_FILES));

  if (matrix == NULL) {
    sprintf(ErrorBuffer,
	    "Unable to read matrix for scoring against the sequences.");
    ErrorReport(SERIOUS_ERR_LVL);
    sprintf(ErrorBuffer,
	    "Not scoring the matrix.\n");
    ErrorReport(SERIOUS_ERR_LVL);
    return;
  }

  if (NumberToReport == 0) {
    NumberToReport = DEFAULT_MATRIX_SEARCH_NUMBER_TO_REPORT_M(matrix);
  }

  /* inform what matrix the search is being done with */
  sprintf(ErrorBuffer,
	  "Searching matrix %s against the sequence database.\n",
	  matrix->number);
  ErrorReport(INFO_ERR_LVL);

  /* output the matrix */
  if (emfp != NULL) {
    output_matrix(matrix, emfp);
  }

  /* make sure we are at the beginning of the sequence data */
  rewind_file(SEQUENCE_FILES);

  records_searched = 0;

  /* loop through all the sequences */
  while ((sequence = read_a_sequence(get_file(SEQUENCE_FILES),
				     get_sequence_db_db_type(),
				     get_sequence_db_seq_type())) != NULL) {

    /* inform which sequence scoring against */
    sprintf(ErrorBuffer,
	    "Scoring vs sequence %s.\n",
	    sequence->name);
    ErrorReport(INFO_ERR_LVL);

    /* is the sequence a NA_SEQ?  if so, translate */
    if (sequence->type == NA_SEQ) {
      /* translate */
      trans1 = translate_sequence(sequence, 1, gcode, revgcode);
      trans2 = translate_sequence(sequence, 2, gcode, revgcode);
      trans3 = translate_sequence(sequence, 3, gcode, revgcode);
      /* score and enter the data into the list */
      score_and_enter(trans1, matrix, 1, RepeatsAllowed, SearchType);
      score_and_enter(trans2, matrix, 2, RepeatsAllowed, SearchType);
      score_and_enter(trans3, matrix, 3, RepeatsAllowed, SearchType);
      /* free the translated sequences */
      free_sequence(trans1);
      free_sequence(trans2);
      free_sequence(trans3);
      /* if we are supposed to search both strands, translate and score */
      if (StrandsToSearch == 2) {
	/* translate */
	trans_1 = translate_sequence(sequence, -1, gcode, revgcode);
	trans_2 = translate_sequence(sequence, -2, gcode, revgcode);
	trans_3 = translate_sequence(sequence, -3, gcode, revgcode);
	/* score and enter the data into the list */
	score_and_enter(trans_1, matrix, -1, RepeatsAllowed, SearchType);
	score_and_enter(trans_2, matrix, -2, RepeatsAllowed, SearchType);
	score_and_enter(trans_3, matrix, -3, RepeatsAllowed, SearchType);
	/* free the translated sequences */
	free_sequence(trans_1);
	free_sequence(trans_2);
	free_sequence(trans_3);
      }
    }
    /* otherwise, just use the sequence, assuming it is an AA_SEQ */
    else {
      /* score and enter the data into the list */
      score_and_enter(sequence, matrix, 0, RepeatsAllowed, SearchType);
    }

    free_sequence(sequence);

    records_searched++;

  } /* end while there are sequences */

  free_matrix(matrix);
}


/*
 * matrices_vs_sequences
 *   When SearchType is UNKNOWN, there is probably a matrices database against
 *   a sequences database.  Score all pairwise matches of matrices and
 *   sequences.
 *   Parameters:
 *     FILE *emfp: the file to output the matrix information to.
 *   Error codes:
 */

static void matrices_vs_sequences(emfp)
     FILE *emfp;
{
  /* NOTE: If there are speed problems, might look at the SearchType */
  /*       variable in files.c to see which of sequences or blocks is */
  /*       considered the database.  The one that is the database is most */
  /*       likely to have the most elements, and should be in the outer */
  /*       loop to reduce "block->matrix" and "if NA_SEQ; NA_SEQ->AA_SEQ" */
  /*       overhead. */

  /* set the number to report */
  if (NumberToReport == 0) {
    rewind_file(MATRIX_FILES);
    matrix = read_a_matrix(get_file(MATRIX_FILES));
    rewind_file(SEQUENCE_FILES);
    sequence = read_a_sequence(get_file(SEQUENCE_FILES),
			       get_sequence_db_db_type(),
			       get_sequence_db_seq_type());
    if (matrix == NULL) {
      if (sequence == NULL) {
	/* there are no comparisons to do */
	NumberToReport = 1;
      }
      else {
	NumberToReport = DEFAULT_BLOCK_SEARCH_NUMBER_TO_REPORT(sequence);
	free_sequence(sequence);
      }
    }
    else {
      if (sequence == NULL) {
	/* there are no comparisons to do */
	NumberToReport = DEFAULT_MATRIX_SEARCH_NUMBER_TO_REPORT_M(matrix);
	free_matrix(matrix);
      }
      else {
	NumberToReport =
	  max((DEFAULT_BLOCK_SEARCH_NUMBER_TO_REPORT(sequence)),
	      (DEFAULT_MATRIX_SEARCH_NUMBER_TO_REPORT_M(matrix)));
	free_block(matrix);
	free_sequence(sequence);
      }
    }
    rewind_file(MATRIX_FILES);
    rewind_file(SEQUENCE_FILES);
  }

  /* make sure we are at the begining of the matrix data */
  rewind_file(MATRIX_FILES);

  /* loop through all the matrices */
  while ((matrix = read_a_matrix(get_file(MATRIX_FILES))) != NULL) {

    /* output the matrix */
    if (emfp != NULL) {
      output_matrix(matrix, emfp);
    }

    /* make sure we are at the beginning of the sequence data */
    rewind_file(SEQUENCE_FILES);

    /* loop through all the sequences */
    while ((sequence = read_a_sequence(get_file(SEQUENCE_FILES),
				       get_sequence_db_db_type(),
				       get_sequence_db_seq_type())) != NULL) {

      /* inform which matrix and sequence are being scored */
      sprintf(ErrorBuffer,
	      "Scoring sequence %s vs matrix %s.\n",
	      sequence->name, matrix->number);
      ErrorReport(INFO_ERR_LVL);

      /* is the sequence a NA_SEQ?  if so, translate */
      if (sequence->type == NA_SEQ) {
	/* translate */
	trans1 = translate_sequence(sequence, 1, gcode, revgcode);
	trans2 = translate_sequence(sequence, 2, gcode, revgcode);
	trans3 = translate_sequence(sequence, 3, gcode, revgcode);
	/* score and enter the data into the list */
	score_and_enter(trans1, matrix, 1, RepeatsAllowed, SearchType);
	score_and_enter(trans2, matrix, 2, RepeatsAllowed, SearchType);
	score_and_enter(trans3, matrix, 3, RepeatsAllowed, SearchType);
	/* free the translated sequences */
	free_sequence(trans1);
	free_sequence(trans2);
	free_sequence(trans3);
	/* if we are supposed to search both strands, translate and score */
	if (StrandsToSearch == 2) {
	  /* translate */
	  trans_1 = translate_sequence(sequence, -1, gcode, revgcode);
	  trans_2 = translate_sequence(sequence, -2, gcode, revgcode);
	  trans_3 = translate_sequence(sequence, -3, gcode, revgcode);
	  /* score and enter the data into the list */
	  score_and_enter(trans_1, matrix, -1, RepeatsAllowed, SearchType);
	  score_and_enter(trans_2, matrix, -2, RepeatsAllowed, SearchType);
	  score_and_enter(trans_3, matrix, -3, RepeatsAllowed, SearchType);
	  /* free the translated sequences */
	  free_sequence(trans_1);
	  free_sequence(trans_2);
	  free_sequence(trans_3);
	}
      }
      /* otherwise, just use the sequence, assuming it is an AA_SEQ */
      else {
	/* score and enter the data into the list */
	score_and_enter(sequence, matrix, 0, RepeatsAllowed, SearchType);
      }

      free_sequence(sequence);

      records_searched++;

    } /* end while there are sequences */

    free_matrix(matrix);

  } /* end while there are matrices */
}





static void print_stats(ofp)
     FILE *ofp;
{
  fprintf(ofp, "Records Searched:   %d\n", records_searched);
  fprintf(ofp, "\n");
  fprintf(ofp, "Scores Done:        %16.0f\n", Scores_Done);
  fprintf(ofp, "\n");
  fprintf(ofp, "Alignments Done:    %16.0f\n", Alignments_Done);
  fprintf(ofp, "\n");

  if (DoHistogram) {
    print_histogram(ofp);
  }
}



/*
 * main
 *   controls flow of program
 *   Parameters: argc, argv
 *   Error codes:
 */

int main(argc, argv)
     int argc;
     char *argv[];
{
  FILE *ofp;			/* the OutputFile file pointer */
  FILE *emfp;			/* the ExportMatrixFile file pointer */
  FILE *fqij; 			/* the Qij file */
  int qargc;			/* args from the OP line in .cs file */
  char **qargv;
  char *blimps_dir;		/* BLIMPS_DIR environment variable */
  int i;

  /* setup version info */
  version_strings("BLIMPS (BLocks IMProved Searcher)",
		  VERSION,
		  "",
		  VERSION_DATE,
		  "(C) Copyright 1993-2000, Fred Hutchinson Cancer Research Center\n\n");

  /* read comand line parameters */
  /* NOTE:!!! Must decide if parameters or config file takes precedence */
  /*          will only use config file to start with, so the only arg */
  /*          should be the config file */

  signal(SIGABRT, ABRT_signal_handler);
#ifdef MALLOC_DEBUG		/* sun specific?  At least not on DECs. */
  malloc_debug(1);		/* needed to raise the error in free(), also */
				/* gives better error messages */
#endif

  /* print the version info to stdout */
  print_version(stdout);
  ErrorLevelReport = PROGRAM_ERR_LVL;     /* Most errors are suppressed */

  /* check that the correct number of argments are given */
  if (argc != 2) {
    fprintf(stderr, "Usage: %s config_file\n", argv[0]);
    exit(1);
  }

  /* read the configuration file, assuming it is arg 1 for now */
  /* NOTE: No errors will be placed in the error file until read_config_file */
  read_config_file(argv[1]);

  /*  See if the BLIMPS_DIR environment variable has been set */
  blimps_dir = getenv("BLIMPS_DIR");

  /* initialize the memory routines */
  init_reclaim_space(blimps_reclaim_space);

  /* initialize the list data structures used for scores, blocks, */
  /* matrices, and sequences */
  initialize_lists();

  /* load the frequencies for converting blocks to matrices */
  /* see if the frequency file has not been specified */
  if (number_of_files(FREQUENCY_FILE) <= 0) {
    /* get a sequence file and find out the type */
    /* open the default file */
    if (blimps_dir != NULL) sprintf(Buffer, "%s/docs/", blimps_dir);
    else Buffer[0] = '\0';
    /* When searching a block vs a DNA database, use codon frequency instead
       of amino acid frequency because most of translated seq is not coding */
    if( (SearchType == SEARCH_TYPE_MATRIX) &&
        (get_sequence_db_seq_type() == NA_SEQ) )
    {
      strcat(Buffer, LOCAL_CODON_FREQUENCY_FILE);
    }
    else {
      strcat(Buffer, LOCAL_AMINO_FREQUENCY_FILE);
    }
  }
  else {
    /* open the listed frequency file name */
    strcpy(Buffer, get_current_file_name(FREQUENCY_FILE));
  }

  sprintf(ErrorBuffer, "Using frequencies from %s.\n", Buffer);
  ErrorReport(INFO_ERR_LVL);
  /*   Creates global array frequency[]  */
  if (! load_frequencies(Buffer)) {
    sprintf(ErrorBuffer, "Error loading frequencies.\n");
    ErrorReport(SERIOUS_ERR_LVL);
  }

  /*  Load the substitution probability matrix for conversion method 3 */
  if (BlockToMatrixConversionMethod >= 3)
  {
     Qij = NULL;
     if (blimps_dir != NULL) sprintf(Buffer, "%s/docs/", blimps_dir);
     else Buffer[0] = '\0';
     strcat(Buffer, LOCAL_QIJ_FILE);
     RTot = LOCAL_QIJ_RTOT;
    /*     Get the arguments from the "OP alts:" line if there is one */
    qargc = -1;
    if (! get_option_args("alts", & qargc, & qargv)) {
      sprintf(ErrorBuffer,
	      "Using default.qij\n");
      ErrorReport(WARNING_ERR_LVL);
    }
    else
    {
       if (qargc > 0)
       {
          RTot = atof(qargv[0]);
          if (qargc > 1)  strcpy(Buffer, qargv[1]);
       }
    }

    if ((fqij = fopen(Buffer, "r")) == NULL) {
      sprintf(ErrorBuffer,
	      "Cannot open 'alts' file \"%s\".  Exiting...\n", Buffer);
      ErrorReport(FATAL_ERR_LVL);
    }

    Qij = load_qij(fqij);
    fclose(fqij);
    sprintf(ErrorBuffer,
	  "Using pseudo-count parameters %f and %s.\n", RTot, Buffer);
    ErrorReport(INFO_ERR_LVL);
  }   /*   end of qij file stuff */


  /* initialize translation codes */
  if ((0>GeneticCodeInitializer) ||
      (GeneticCodeInitializer>=NUMBER_OF_GENETIC_CODES)) {
    sprintf(ErrorBuffer,
	    "Unknown genetic code %d.  Setting to Standard code (0).",
	    GeneticCodeInitializer);
    ErrorReport(WARNING_ERR_LVL);
    GeneticCodeInitializer = 0;
  }
  sprintf(ErrorBuffer,
	  "Using genetic code %d; %s.\n", GeneticCodeInitializer,
	  gcodes[GeneticCodeInitializer].name);
  ErrorReport(INFO_ERR_LVL);
  init_gcode(&gcodes[GeneticCodeInitializer], gcode, revgcode);

  /* open the output file for writing */
  if (OutputFile[0] != '\0') {
    /* open the output file for reading */
    ofp = fopen(OutputFile, "w");

    if (ofp == NULL) {
      sprintf(ErrorBuffer,
	      "Unable to open output file \"%s\" for writing, aborting.",
	      OutputFile);
      ErrorReport(FATAL_ERR_LVL);
    }
  }
  else {
    ofp = stdout;
  }


  /* open the export file for writing */
  if (ExportMatrixFile[0] != '\0') {
    /* open the output file for reading */
    emfp = fopen(ExportMatrixFile, "w");

    if (emfp == NULL) {
      sprintf(ErrorBuffer,
	      "Unable to open output file \"%s\" for writing.",
	      ExportMatrixFile);
      ErrorReport(SERIOUS_ERR_LVL);
      sprintf(ErrorBuffer,
	      "Not outputing the matrix output.\n");
      ErrorReport(SERIOUS_ERR_LVL);
      emfp = NULL;
    }
  }
  else {
    emfp = NULL;
  }


  /*
   * read through the blocks and sequences
   */
  switch (SearchType) {
  case SEARCH_TYPE_BLOCK :
    if (SiteSpecificScoringMatrixType == SSSM_BLOCK) {
      sequence_vs_blocks(emfp);
    }
    else if (SiteSpecificScoringMatrixType == SSSM_PRECOMP_MAT) {
      sequence_vs_matrices(emfp);
    }
    break;
  case SEARCH_TYPE_MATRIX :
    if (SiteSpecificScoringMatrixType == SSSM_BLOCK) {
      block_vs_sequences(emfp);
    }
    else if (SiteSpecificScoringMatrixType == SSSM_PRECOMP_MAT) {
      matrix_vs_sequences(emfp);
    }
    break;
  case SEARCH_TYPE_UNKNOWN :
    if (SiteSpecificScoringMatrixType == SSSM_BLOCK) {
      blocks_vs_sequences(emfp); /* note plural blocks */
    }
    else if (SiteSpecificScoringMatrixType == SSSM_PRECOMP_MAT) {
      matrices_vs_sequences(emfp); /* note plural matrices */
    }
    break;
  default:
    /* BIG ERROR, should never get here.  Memory corrupted */
    sprintf(ErrorBuffer,
	    "main(): SearchType variable corrupted, possible memory problems.");
    ErrorReport(PROGRAM_ERR_LVL);
    sprintf(ErrorBuffer,
	    "        Should never have reached this part of the program.");
    ErrorReport(PROGRAM_ERR_LVL);
    sprintf(ErrorBuffer,
	    "        The rest of the run may be invalid.");
    ErrorReport(PROGRAM_ERR_LVL);
    sprintf(ErrorBuffer,
	    "        Assuming SQ vs BL search.  Setting SearchType to UNKNOWN.\n");
    ErrorReport(PROGRAM_ERR_LVL);
    SearchType = SEARCH_TYPE_UNKNOWN;
    blocks_vs_sequences(); /* note plural blocks */
    break;
  }


  /* print out header information */

  /* print the version info */
  print_version(ofp);


  switch (SearchType) {
  case SEARCH_TYPE_BLOCK :
    if (sequence != NULL) {
      fprintf(ofp, "\n");
      fprintf(ofp, "Probe Sequence: %s %s\n",
             sequence->name, sequence->info);
      fprintf(ofp, "\n");
      if (sequence->type == NA_SEQ) {
	fprintf(ofp, "Probe Size: %d Base Pairs\n", sequence->length);
      }
      else {
	fprintf(ofp, "Probe Size: %d Amino Acids\n", sequence->length);
      }
      fprintf(ofp, "\n");
    }
    else {
      fprintf(ofp, "\n");
      fprintf(ofp, "Probe Sequence: NO SEQUENCE\n");
      fprintf(ofp, "\n");
      fprintf(ofp, "Probe Size: UNKNOWN\n");
      fprintf(ofp, "\n");
    }

    fprintf(ofp, "Probe File: %s\n", get_file_name(0, SEQUENCE_FILES));
    fprintf(ofp, "\n");
    if (SiteSpecificScoringMatrixType == SSSM_BLOCK) {
      fprintf(ofp, "Target File (s) : %s\n", get_file_name(0, BLOCK_FILES));
      for(i=1; i<number_of_files(BLOCK_FILES); i++) {
	fprintf(ofp, "                  %s\n", get_file_name(i, BLOCK_FILES));
      }
    }
    else if (SiteSpecificScoringMatrixType == SSSM_PRECOMP_MAT) {
      fprintf(ofp, "Target File (s) : %s\n", get_file_name(0, MATRIX_FILES));
      for(i=1; i<number_of_files(MATRIX_FILES); i++) {
	fprintf(ofp, "                  %s\n", get_file_name(i, MATRIX_FILES));
      }
    }

    fprintf(ofp, "\n");

    print_stats(ofp);

    fprintf(ofp, "AC#         Description");
    fprintf(ofp, "                                         Strength  Score ");
    fprintf(ofp, "RF    AA#\n");
    break;
  case SEARCH_TYPE_MATRIX :
    fprintf(ofp, "\n");
    fprintf(ofp, "\n");
    if (SiteSpecificScoringMatrixType == SSSM_BLOCK) {
      fprintf(ofp, "Block File: %s\n", get_file_name(0, BLOCK_FILES));
    }
    else if (SiteSpecificScoringMatrixType == SSSM_PRECOMP_MAT) {
      fprintf(ofp, "Matrix File: %s\n", get_file_name(0, MATRIX_FILES));
    }
    fprintf(ofp, "\n");
    fprintf(ofp, "Target File (s) : %s\n", get_file_name(0, SEQUENCE_FILES));
    for(i=1; i<number_of_files(SEQUENCE_FILES); i++) {
      fprintf(ofp, "                  %s\n", get_file_name(i, SEQUENCE_FILES));
    }
    fprintf(ofp, "\n");

    print_stats(ofp);

    fprintf(ofp, "AC#                    Description");
    fprintf(ofp, "                                                Score ");
    fprintf(ofp, "RF  AA# Length\n");
    break;
    break;
  case SEARCH_TYPE_UNKNOWN :
    if (SiteSpecificScoringMatrixType == SSSM_BLOCK) {
      fprintf(ofp, "Block File (s) :    %s\n", get_file_name(0, BLOCK_FILES));
      for(i=1; i<number_of_files(BLOCK_FILES); i++) {
	fprintf(ofp, "                   %s\n", get_file_name(i, BLOCK_FILES));
      }
    }
    else if (SiteSpecificScoringMatrixType == SSSM_PRECOMP_MAT) {
      fprintf(ofp, "Matrix File (s) :    %s\n",get_file_name(0, MATRIX_FILES));
      for(i=1; i<number_of_files(MATRIX_FILES); i++) {
	fprintf(ofp, "                   %s\n",get_file_name(i, MATRIX_FILES));
      }
    }
    fprintf(ofp, "\n");

    fprintf(ofp, "Sequence File (s) : %s\n", get_file_name(0, SEQUENCE_FILES));
    for(i=1; i<number_of_files(SEQUENCE_FILES); i++) {
      fprintf(ofp, "                   %s\n", get_file_name(i, SEQUENCE_FILES));
    }
    fprintf(ofp, "\n");

    print_stats(ofp);

    break;
  default:
    break;
  }
  /* output the data */
  output_scores(NumberToReport, ofp);

  /* close the output files */
  fclose(ofp);
  if (emfp != NULL) {
    fclose(emfp);
  }

  /* exit with a valid code, also closes any open files that were not closed */
  exit(0);

}  /*   end of main */



/* Change log information follows. */
/*
 Changes since version 3.5:
  7/22/02 Set ErrorLevelReport
 Changes since version 3.3.1:
 10/ 2/99 Added <string.h> & <math.h> for LINUX users.
 Changes since version 3.2.5:
  5/22/99 Only translate query sequence once in sequence_vs_blocks()
  2/23/99 Changed sequence_vs_blocks() to use read_a_block_faster() which
	  ignores clusters & assumes sequence weights are in the block.
  2/ 3/99 Clean up sequence_vs_blocks()
 Changes since version 3.2.4:
 12/12/98 Process StrandsToSearch == -1
          Used default.amino.frq instead of default.codon.frq for
          DNA seq vs blocks.
 Changes since version 3.1.4:
 5/13/98  Adjust block vs sequence search output for longer seq names.
 Changes since version 3.1:
 1/30/97  Changed use of default qij matrix error from SERIOUS to WARNING.
 1/20/97  Fixed bug in calls to get_file_name() when multiple files
 Changes since version 3.0.0:
 4/22/96  Added SequenceType global variable.
 4/11/96  Changed blocks_to_sequences() to call output_matrix_s(matrix,
          FLOAT_OUTPUT) instead of output_matrix(matrix). This is a bad
          kludge to print a matrix in floating point format by using an
	  "unknown" search type combination of BL with SQ.    JGH
 4/ 3/96  Changed to assume all conversion types >= 3 require loading a
	  Qij file.  JGH
*/
