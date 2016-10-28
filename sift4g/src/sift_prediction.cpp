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

// BLIMPS needs that...
#define EXTERN

extern "C" {
    #include "sift/info_on_seqs_pid.h"
    #include "sift/Alignment.h"
}

constexpr uint32_t kMaxSequences = 400;

// declared extern in sift...
FILE* errorfp;

Sequence* chainPtrToSequencePtr(Chain* chain_ptr, int32_t id) {

    Sequence* sequence_ptr = (Sequence*) malloc(sizeof(Sequence));

    std::string chain_name = std::to_string(id);
    strcpy(sequence_ptr->name, chain_name.c_str());
    sequence_ptr->position = 0;

    sequence_ptr->length = chainGetLength(chain_ptr);
    sequence_ptr->max_length = chainGetLength(chain_ptr);
    sequence_ptr->type = AA_SEQ;
    sequence_ptr->weight = 0.0;

    sequence_ptr->undefined = 0;
    sequence_ptr->undefined_dbl = 0.0;
    sequence_ptr->undefined_ptr = NULL;
    sequence_ptr->sequence = (Residue*) malloc(sequence_ptr->max_length * sizeof(Residue));

    for (int32_t i = 0; i < chainGetLength(chain_ptr); ++i) {
        int32_t stemp = aa_atob[(int) chainGetChar(chain_ptr, i)];
        assert(stemp > -1 && stemp < AAID_MAX + 1 && "error while transforming Chain to Sequence");
        sequence_ptr->sequence[i] = (Residue) stemp;
    }

    return sequence_ptr;
}

void siftPredictionsNoSIFT(std::vector<std::vector<Chain*>>& alignment_strings,
    Chain** queries, int32_t queries_length, const std::string& subst_path,
    int32_t sequence_identity, const std::string& out_path) {

/* warning, this code changes the alignment_strings by adding query 
to the alignment  */
    cout << "-------------------" << endl;
    cout << "GIS code starting now" << endl; 
    fprintf(stderr, "** Generating SIFT predictions with sequence identity: %.2f%% **\n", (float) sequence_identity);

    std::string subst_extension = ".subst";
    std::string out_extension = ".SIFTprediction_new";

    char* error_file_name = createFileName("sift", out_path, ".err");
    errorfp = fopen(error_file_name, "w");

    for (int32_t i = 0; i < queries_length; ++i) {

	/* add query sequence to the beginning of the alignment */
       /* alignment_strings[i].insert (alignment_strings[i].begin(), queries[i]);
*/
	/* only keep first 399 hits, erase more distant ones */
	if (alignment_strings[i].size() > kMaxSequences -1 ) {
		for (int32_t j = kMaxSequences -1; j < int (alignment_strings[i].size()); j++) {
		alignment_strings[i].erase(alignment_strings[i].begin() + j);
		}
	}

	cout << chainGetName(queries[i]) << endl;
	int query_length = chainGetLength (queries[i]); 

	cout << "new-SIFT before num of seq " << to_string (alignment_strings[i].size()) << endl;	
	remove_seqs_percent_identical_to_query(queries[i], alignment_strings[i], sequence_identity);

	/* add query sequence to the beginning of the alignment */
        alignment_strings[i].insert (alignment_strings[i].begin(), queries[i]);
//	printSeqNames (alignment_strings[i]);
	
	int total_seq = alignment_strings[i].size();
	cout << "new-SIFT num of seq " << total_seq << endl;

	std::vector<std::vector<double>> matrix (
		query_length, std::vector <double> (26));	

        std::vector<std::vector<double>> SIFTscores (
                query_length, std::vector <double> (26));

	std::vector <double> weights_1 (alignment_strings[i].size(), 1.0);
	std::vector<double> aas_stored (query_length);

	createMatrix (alignment_strings[i], queries[i], weights_1, matrix, aas_stored); 
	printMatrix (matrix, "tmp.txt");
	calcSIFTScores (alignment_strings[i], queries[i], matrix, SIFTscores); 	
/*	cout << "about to enter calcSeqWeights " << endl;
*/
	int num_seqs_in_alignment = alignment_strings[i].size();
	std::vector <double> seq_weights (num_seqs_in_alignment);
	std::vector <double> number_of_diff_aas (query_length);
	calcSeqWeights (alignment_strings[i], matrix, aas_stored, seq_weights,
	  number_of_diff_aas);

/* Shimin, at this point the alignment is identical to what is being passed in to old SIFT, and seq_weights have been calculated. put basic_matrix_construction, and print it out here. */

	/* printMatrix (matrix_shimin, "shimin.txt"); */

        query_log(i + 1, queries_length);
//	cout << "about to read substitution file " << endl;
        char* subst_file_name = createFileName(chainGetName(queries[i]), subst_path, subst_extension);
	char* out_file_name = createFileName(chainGetName(queries[i]), out_path, out_extension);

	if ( exists(subst_file_name)) {
		std::list <std::string> subst_list;
		std::unordered_map <std::string, double> medianSeqInfoForPos;

		readSubstFile (subst_file_name, subst_list);
		hashPredictedPos (subst_list, medianSeqInfoForPos);
		addPosWithDelRef (queries[i],  SIFTscores, medianSeqInfoForPos);
		addMedianSeqInfo (alignment_strings[i], queries[i], matrix, medianSeqInfoForPos);
//		cout << "hey there " << endl;
		printSubstFile (subst_list, medianSeqInfoForPos, SIFTscores, aas_stored, total_seq, queries[i], out_file_name); 
	} else {
		printMatrix (SIFTscores,  out_file_name); 	
	}

//	cout << "query " << to_string(i) << endl;	
	} // end for loop through all queries
//	cout << "done with all queries " << endl;
}

void siftPredictions(const std::vector<std::vector<Chain*>>& alignment_strings,
    Chain** queries, int32_t queries_length, const std::string& subst_path,
    int32_t sequence_identity, const std::string& out_path) {

    fprintf(stderr, "** Generating SIFT predictions with sequence identity: %.2f%% **\n", (float) sequence_identity);

    std::string subst_extension = ".subst";
    std::string out_extension = ".SIFTprediction";

    char* error_file_name = createFileName("sift", out_path, ".err");
    errorfp = fopen(error_file_name, "w");


    for (int32_t i = 0; i < queries_length; ++i) {

	query_log(i + 1, queries_length);
        char* subst_file_name = createFileName(chainGetName(queries[i]), subst_path, subst_extension);

        FILE* subst_fp = exists(subst_file_name) ? fopen(subst_file_name, "r") : nullptr;
        char* out_file_name = createFileName(chainGetName(queries[i]), out_path, out_extension);
        FILE* out_fp = fopen(out_file_name, "w");

	cout << "results in " << out_file_name <<endl ;
        int32_t num_sequences = (kMaxSequences - 1) < alignment_strings[i].size() ? (kMaxSequences - 1) : alignment_strings[i].size();

	/* adding the query to the alignment */
        Sequence** sequences = new Sequence*[num_sequences + 1];
        sequences[0] = chainPtrToSequencePtr(queries[i], 0);

        for (int32_t j = 0; j < num_sequences; ++j) {
            sequences[j + 1] = chainPtrToSequencePtr(alignment_strings[i][j], j + 1);
        }

        generate_predictions(sequences, num_sequences + 1, subst_fp, sequence_identity, out_fp);

        for (int32_t j = 0; j < num_sequences + 1; ++j) {
            if (sequences[j] != NULL) {
                free(sequences[j]->sequence);
                free(sequences[j]);
            }
        }
        delete[] sequences;

        fclose(out_fp);
        delete[] out_file_name;

        if (subst_fp != nullptr) fclose(subst_fp);
        delete[] subst_file_name;
    }

    fclose(errorfp);
    delete[] error_file_name;

    fprintf(stderr, "\n\n");
}
