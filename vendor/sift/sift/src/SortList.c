/* (C) Copyright 2000, Fred Hutchinson Cancer Research Center */
/* Use, modification or distribution of these programs is subject to */
/* the terms of the non-commercial licensing agreement in license.h. */

#ifndef _SORTLIST_C_
#define _SORTLIST_C_

#include <assert.h>

#include "PN_convert.h"
#include "SortList.h"

extern FILE* errorfp;

int
lowest_scoring_aa (Matrix* matrix, int pos)
{
	struct working* col;
	int aa, min_aa = 1;
	double min;

	min = 10000;
	col = make_col();
	counts (matrix->block, col, pos);
	for (aa = 1 ; aa < AAS; aa++) {
		if (col->cnt[aa] > 0.0 && matrix->weights[aa][pos] < min) {
			min = matrix->weights[aa][pos];
			min_aa = aa;
		}
	}

	free (col);
	return min_aa;
} /* end of lowest_scoring_aa */

void print_list (struct aa_score list[20], Matrix* matrix, int pos)
{
	int i; struct working* col;
	FILE* outfp; char filename[LARGE_BUFF_LENGTH];

	sprintf (filename, "TEMP%d", pos);
	if ( (outfp = fopen (filename, "w")) == NULL) {
		fprintf (errorfp, "couldn't open %s\n", filename);
		exit(-1);
	}

	col = make_col();
	counts (matrix->block, col, pos);
	for (i = 0; i < 20; i++) {
		fprintf (outfp, "%c", aa_btoa[list[i].aa]);
		if (col->cnt[list[i].aa] > 0.0)  {
			fprintf (outfp, "*");
		}
		fprintf (outfp, "  %d  %.3f\n", i, list[i].score);
	}
	fclose (outfp);
	free(col);

} /* end of print_list */

void
copy_values (struct aa_score list[20], Matrix* matrix, int pos)
{
	char c;
	int i;

	i = 0;
	for (c = 'A'; c < 'Z'; c++) {
		if (c != 'B' && c != 'J' && c != 'O' && c != 'U' &&
							c != 'X' && c!= 'Z') {
			list[i].aa = aa_atob[c];
			list[i].score = matrix->weights[aa_atob[c]][pos];
			list[i].bool_variable = FALSE;
			i++;
		}
	}
	assert (i == 20);
} /* end of copy_values */

void one_tailed_CI_for_matrix (Matrix* matrix, const double confidence_interval)
{
	int pos;

	for (pos = 0; pos < matrix->width; pos++) {
		one_tailed_CI_for_pos (matrix, pos, confidence_interval);
	}
}

void one_tailed_CI_for_pos (Matrix* matrix, const int pos,
		    const double confidence_interval)
{
	struct aa_score list[20];
	double cum_prob;
	int i, aa;

	cum_prob = 0.0;
	assert (confidence_interval >= 0.0 && confidence_interval <= 1.0);
	copy_values (list, matrix, pos);
	sort_list (list);
	i = 0;
	while (cum_prob + list[i].score <= confidence_interval && i < 20) {
		cum_prob += list[i].score;
		list[i].bool_variable = TRUE;
		i++;
	}
	/* allow the last aa that just exceeds confidence interval */
	if (i != 20) {
		cum_prob += list[i].score;
		list[i].bool_variable = TRUE;
	}

	for (i = 0; i  < 20; i++) {
		if (list[i].bool_variable == FALSE) {
			matrix->weights[list[i].aa][pos] = 0.0;
		}
	}
} /* end of one_tailed_CI */

void
sort_list_in_PIMA_order (struct aa_score list[20])
{

	struct aa_score templist[20];
	char aa_order[20];
	int i, aa, index; double score;

	for (i = 0;  i < 20; i++) {
		score = list[i].score;
		aa = list[i].aa;
		index = PIMA_index(aa_btoa[aa]);
		templist[index].score = score;
		templist[index].aa = aa;
	}
	for (i = 0; i < 20; i++) {
		list[i].score = templist[i].score;
		list[i].aa = templist[i].aa;
	}
} /* end of sort list in PIMA order */

int
PIMA_index (char c)
{
	switch (c)
	{
		case 'D': return 0; break;
		case 'E': return 1;
		case 'K': return 2;
		case 'R': return 3;
		case 'H': return 4;
		case 'N': return 5;
		case 'Q': return 6;
		case 'S': return 7;
		case 'T': return 8;
		case 'I': return 9;
		case 'L': return 10;
		case 'V': return 11;
		case 'F': return 12;
		case 'W': return 13;
		case 'Y': return 14;
		case 'C': return 15;
		case 'M': return 16;
		case 'A': return 17;
		case 'G': return 18;
		case 'P': return 19;
	}

	return -1; // don't know which value to return if character is not in the list above
}

/* sort list of amino acid scores from highest to lowest */
void
sort_list (struct aa_score list[20])
{
	int j, P;
	double tmp;
	int tmp_aa;
	int tmp_bool_variable;

	for (P = 1; P < 20; P++) {
		tmp = list[P].score;
		tmp_aa = list[P].aa;
		tmp_bool_variable = list[P].bool_variable;
		for (j = P; j > 0 && list[j-1].score < tmp; j--) {
			list[j].score = list[j-1].score;
			list[j].aa = list[j-1].aa;
			list[j].bool_variable = list[j-1].bool_variable;
		}
		list[j].score = tmp;
		list[j].aa = tmp_aa;
		list[j].bool_variable = tmp_bool_variable;
	}
} /* end of sort_list */

#endif
