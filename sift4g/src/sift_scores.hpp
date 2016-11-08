/*!
 * @file sift_prediction.hpp
 *
 * @brief SIFT predictions header file
 *
 * @author: pauline-ng, angsm
 */

#pragma once

#include <stdint.h>
#include <vector>
#include <string>
#include <list>
#include <unordered_map>

#include "swsharp/swsharp.h"

void readSubstFile(char* substFilename, std::list<std::string>& substList);

void printSubstFile(const std::list<std::string>& substList, std::unordered_map<std::string, double>& medianSeqInfoForPos,
    const std::vector<std::vector<double>>& SIFTscores, const std::vector<double>& aas_stored,
    const int total_seq, Chain* query, const std::string outfile);

void createMatrix(const std::vector<Chain*>& alignment_string, Chain* query,
    const std::vector<double>& seq_weights, std::vector<std::vector<double>>& matrix,
    std::vector<double>& tot_pos_weight) ;

void remove_seqs_percent_identical_to_query(Chain *queries, std::vector<Chain*>& alignment_string, double seq_identity);

/* this has to be after createMatrix because it uses amino_acids_present */
void calcSeqWeights(const std::vector<Chain*>& alignment_string, std::vector<std::vector<double>>& matrix,
    std::vector<double>& amino_acids_present, std::vector<double>& seq_weights, std::vector<double>& number_of_diff_aas);

void seqs_without_X(const std::vector<Chain*>& alignment_string, int pos, std::vector<Chain*>& out_alignment);

void basic_matrix_construction(const std::vector<Chain*>& alignment_string,
    const std::vector<double>& seq_weights, std::vector<std::vector<double>>& matrix);

int aa_to_idx(char character);

bool valid_amino_acid(char aa);

double calculateMedianSeqInfo(std::vector<Chain*>& alignment_strings, std::unordered_map<int, int>& invalidSeq,
    std::vector<double>& seq_weights, std::vector<std::vector<double>>& matrix);

void hashPredictedPos(std::list<std::string>& substList, std::unordered_map<std::string, double>& medianSeqInfoForPos);

void addPosWithDelRef(Chain* query, std::vector<std::vector<double>>& SIFTscores,
    std::unordered_map<std::string, double>& medianSeqInfoForPos);

void addMedianSeqInfo(std::vector<Chain*>& alignment_string, Chain* query,
    std::vector<std::vector<double>>& matrix, std::unordered_map<std::string, double>& medianSeqInfoForPos);

void printSeqNames(std::vector<Chain*>& alignment_string);

void check_refaa_against_query(char ref_aa, int aa_pos, Chain* query, std::ofstream& outfp);
void printMatrix(std::vector<std::vector<double>>& matrix, std::string filename);

void calcDiri(std::vector<std::vector<double>>& weighted_matrix, std::vector<std::vector<double>>& diri_matrix);

void calcSIFTScores(std::vector<Chain*>& alignment_string, Chain* query,
    std::vector<std::vector <double>>& matrix, std::vector<std::vector <double>>& SIFTscores);

void add_diric_values(std::vector<double>& count_col, std::vector<double>& diric_col);

double add_logs(double logx, double logy);
