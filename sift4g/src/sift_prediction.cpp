/*!
 * @file sift_prediction.cpp
 *
 * @brief SIFT predictions source file
 */

#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <iostream>

#include "utils.hpp"
#include "sift_prediction.hpp"
#include "sift_prediction_scores.hpp"

constexpr uint32_t kMaxSequences = 400;

// declared extern in sift...
FILE* errorfp;


void siftPredictions(std::vector<std::vector<Chain*>>& alignment_strings,
    Chain** queries, int32_t queries_length, const std::string& subst_path,
    int32_t sequence_identity, const std::string& out_path) {

/* warning, this code changes the alignment_strings by adding query 
to the alignment  */
	fprintf(stderr, "** Generating SIFT predictions with sequence identity: %.2f%% **\n", (float) sequence_identity);
	
    std::string subst_extension = ".subst";
    std::string out_extension = ".SIFTprediction_new";

    char* error_file_name = createFileName("sift", out_path, ".err");
    errorfp = fopen(error_file_name, "w");

    for (int32_t i = 0; i < queries_length; ++i) {

	/* only keep first 399 hits, erase more distant ones */
	if (alignment_strings[i].size() > kMaxSequences -1 ) {
		for (int32_t j = kMaxSequences -1; j < int (alignment_strings[i].size()); j++) {
		alignment_strings[i].erase(alignment_strings[i].begin() + j);
		}
	}

	int query_length = chainGetLength (queries[i]); 
	remove_seqs_percent_identical_to_query(queries[i], alignment_strings[i], sequence_identity);

	/* add query sequence to the beginning of the alignment */
        alignment_strings[i].insert (alignment_strings[i].begin(), queries[i]);
	
	int total_seq = alignment_strings[i].size();

	std::vector<std::vector<double>> matrix (
		query_length, std::vector <double> (26));	

        std::vector<std::vector<double>> SIFTscores (
                query_length, std::vector <double> (26));

	std::vector <double> weights_1 (alignment_strings[i].size(), 1.0);
	std::vector<double> aas_stored (query_length);

	createMatrix (alignment_strings[i], queries[i], weights_1, matrix, aas_stored); 
	printMatrix (matrix, "tmp.txt");
	calcSIFTScores (alignment_strings[i], queries[i], matrix, SIFTscores); 	

	int num_seqs_in_alignment = alignment_strings[i].size();
	std::vector <double> seq_weights (num_seqs_in_alignment);
	std::vector <double> number_of_diff_aas (query_length);
	calcSeqWeights (alignment_strings[i], matrix, aas_stored, seq_weights,
	  number_of_diff_aas);

        query_log(i + 1, queries_length);
        char* subst_file_name = createFileName(chainGetName(queries[i]), subst_path, subst_extension);
	char* out_file_name = createFileName(chainGetName(queries[i]), out_path, out_extension);

	if ( exists(subst_file_name)) {
		std::list <std::string> subst_list;
		std::unordered_map <std::string, double> medianSeqInfoForPos;

		readSubstFile (subst_file_name, subst_list);
		hashPredictedPos (subst_list, medianSeqInfoForPos);
		addPosWithDelRef (queries[i],  SIFTscores, medianSeqInfoForPos);
		addMedianSeqInfo (alignment_strings[i], queries[i], matrix, medianSeqInfoForPos);
		printSubstFile (subst_list, medianSeqInfoForPos, SIFTscores, aas_stored, total_seq, queries[i], out_file_name); 
	} else {
		printMatrix (SIFTscores,  out_file_name); 	
	}

	} // end for loop through all queries
}
