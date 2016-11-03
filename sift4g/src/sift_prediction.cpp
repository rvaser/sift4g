/*!
 * @file sift_prediction.cpp
 *
 * @brief SIFT predictions source file
 */

#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <iostream>

#include <unistd.h>
#include <sys/wait.h>

#include "utils.hpp"
#include "sift_prediction.hpp"
#include "sift_prediction_scores.hpp"

constexpr uint32_t kMaxSequences = 400;
int qCount = 0;
// declared extern in sift...
FILE* errorfp;

void waitForChild(){
	while (true) {
    		int status;
    		pid_t done = wait(&status);
    		if (done == -1) {
        		if (errno == ECHILD) break; // no more child processes
    		} else {
        		if (!WIFEXITED(status) || WEXITSTATUS(status) != 0) {
            			cerr << "pid " << done << " failed" << endl;
            			//exit(1);
        		}
    		}
	}

}

void siftPredictions(std::vector<std::vector<Chain*>>& alignment_strings,
    Chain** queries, int32_t queries_length, const std::string& subst_path,
    int32_t sequence_identity, const std::string& out_path, int32_t num_threads) {

    /* warning, this code changes the alignment_strings by adding query 
    to the alignment  */
    fprintf(stderr, "** Generating SIFT predictions with sequence identity: %.2f%% **\n", (float) sequence_identity);

    char* error_file_name = createFileName("sift", out_path, ".err");
    errorfp = fopen(error_file_name, "w");

    /*FORK*/
    for (int32_t i = 0; i < queries_length; i++){

	pid_t pid;
	/*if exceed limit, (default 8), wait for all child to finish*/
	if(qCount > num_threads){
		waitForChild();
		qCount = 0;
	}

	query_log(i + 1, queries_length);
	qCount++; 
	pid = fork();
	
    	if (pid == 0){
        	// child process
        	processQuery(sequence_identity, alignment_strings[i], queries[i],queries_length, subst_path, out_path, i);
		
		exit(0);
    	}
    }
    /*wait for remaining child to finish up*/
    waitForChild();

    fclose(errorfp);
}

void processQuery(int32_t sequence_identity, std::vector<Chain*>& alignment_string, Chain* query,int32_t queries_length, const std::string& subst_path, const std::string& out_path, int32_t idx){

    	std::string subst_extension = ".subst";
    	std::string out_extension = ".SIFTprediction_new";

	/* only keep first 399 hits, erase more distant ones */
	if (alignment_string.size() > kMaxSequences -1 ) {
		for (int32_t j = kMaxSequences -1; j < int (alignment_string.size()); j++) {
		alignment_string.erase(alignment_string.begin() + j);
		}
	}

	int query_length = chainGetLength (query); 
	remove_seqs_percent_identical_to_query(query, alignment_string, sequence_identity);

	/* add query sequence to the beginning of the alignment */
        alignment_string.insert (alignment_string.begin(), query);
	
	int total_seq = alignment_string.size();

	std::vector<std::vector<double>> matrix (
		query_length, std::vector <double> (26));	

        std::vector<std::vector<double>> SIFTscores (
                query_length, std::vector <double> (26));

	std::vector <double> weights_1 (alignment_string.size(), 1.0);
	std::vector<double> aas_stored (query_length);

	createMatrix (alignment_string, query, weights_1, matrix, aas_stored); 

	calcSIFTScores (alignment_string, query, matrix, SIFTscores); 	

	int num_seqs_in_alignment = alignment_string.size();
	std::vector <double> seq_weights (num_seqs_in_alignment);
	std::vector <double> number_of_diff_aas (query_length);
	calcSeqWeights (alignment_string, matrix, aas_stored, seq_weights,
	  number_of_diff_aas);


        char* subst_file_name = createFileName(chainGetName(query), subst_path, subst_extension);
	char* out_file_name = createFileName(chainGetName(query), out_path, out_extension);

	if ( exists(subst_file_name)) {
		std::list <std::string> subst_list;
		std::unordered_map <std::string, double> medianSeqInfoForPos;

		readSubstFile (subst_file_name, subst_list);
		hashPredictedPos (subst_list, medianSeqInfoForPos);
		addPosWithDelRef (query,  SIFTscores, medianSeqInfoForPos);
		addMedianSeqInfo (alignment_string, query, matrix, medianSeqInfoForPos);
		printSubstFile (subst_list, medianSeqInfoForPos, SIFTscores, aas_stored, total_seq, query, out_file_name); 
	} else {
		printMatrix (SIFTscores,  out_file_name); 	
	}

}
