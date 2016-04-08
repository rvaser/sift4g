/* (C) Copyright 1993-7, Fred Hutchinson Cancer Research Center */
/* Use, modification or distribution of these programs is subject to */
/* the terms of the non-commercial licensing agreement in license.h. */

/* blimps.h: */
/* Written by: Bill Alford */
/* Change log information is at the end of the file. */

#ifndef BLIMPS_H_
#define BLIMPS_H_

/*
 * Default score reporting settings
 */
/*
 * Algorithm for computing the default size of the output:
 *  Query = block:
 *     Save 200 results
 *  Query = sequence:
 *    Protein: save 400 results
 *    DNA:     roof(length in nts/1000) x 4, max = 4000, min = 400
 */
#define MIN_DEFAULT_MATRIX_SEARCH_NUM_TO_REPORT 200
#define MAX_DEFAULT_MATRIX_SEARCH_NUM_TO_REPORT 200
#define MIN_DEFAULT_BLOCK_SEARCH_NUM_TO_REPORT  400
#define MAX_DEFAULT_BLOCK_SEARCH_NUM_TO_REPORT 4000

#define DEFAULT_MATRIX_SEARCH_NUMBER_TO_REPORT(b)        \
  MAX_DEFAULT_MATRIX_SEARCH_NUM_TO_REPORT

#define DEFAULT_MATRIX_SEARCH_NUMBER_TO_REPORT_M(m)      \
  MAX_DEFAULT_MATRIX_SEARCH_NUM_TO_REPORT

#define DEFAULT_BLOCK_SEARCH_NUMBER_TO_REPORT(s)         \
  ((s->type == NA_SEQ) ?                                 \
   (min(                                                 \
	(max(                                            \
	     MIN_DEFAULT_BLOCK_SEARCH_NUM_TO_REPORT,     \
	     (int)(ceil((double)s->length/1000.0)) * 4)) \
	, MAX_DEFAULT_BLOCK_SEARCH_NUM_TO_REPORT)) :     \
   MIN_DEFAULT_BLOCK_SEARCH_NUM_TO_REPORT )



/* variables set when reading configuration file */

/*
 * search types allowed
 */
#define SEARCH_TYPE_UNSET  -1	/* search type not set yet */
#define SEARCH_TYPE_UNKNOWN 0	/* unknown search type given */
#define SEARCH_TYPE_MATRIX  1	/* sequences vs databases */
#define SEARCH_TYPE_BLOCK   2	/* blocks vs databases */
extern int SearchType;

/*
 *  SequenceType 
*/
extern int SequenceType;

/* 
 * site specific scoring matrix types
 */
#define SSSM_BLOCK       0
#define SSSM_PRECOMP_MAT 1
extern int SiteSpecificScoringMatrixType;

extern int StrandsToSearch;
extern int NumberToReport;

extern Boolean RepeatsAllowed;

extern int GeneticCodeInitializer;

extern int BlockToMatrixConversionMethod;
extern int SequenceMatrixScoringMethod;


#endif /*  BLIMPS_H_ */

/* Change log information follows. */
/* Changes since version 3.0.0:
  4/24/96 Addeded SequenceType
*/
