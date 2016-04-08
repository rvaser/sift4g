/* (C) Copyright 2000, Fred Hutchinson Cancer Research Center */
/* Use, modification or distribution of these programs is subject to */
/* the terms of the non-commercial licensing agreement in license.h. */

#ifndef _SORTLIST_H_
#define _SORTLIST_H_

#include "blimps/blocksprogs.h"

struct aa_score {
        int aa;
        int bool_variable;
	double score;
};

void print_list (struct aa_score list[20], Matrix* matrix, int pos);
void sort_list (struct aa_score list[20]);
int lowest_scoring_aa (Matrix* matrix, int pos);
void print_list (struct aa_score list[20], Matrix* matrix, int pos);
void copy_values (struct aa_score list[20], Matrix* matrix, int pos);

/* confidence interval */
void one_tailed_CI_for_matrix (Matrix* matrix, const double confidence_interval);
void one_tailed_CI_for_pos (Matrix* matrix, const int pos, const double confidence_interval);

/* PIMA stuff */
void sort_list_in_PIMA_order (struct aa_score list[20]);
int PIMA_index (char c);

#endif
