/*!
 * @file sift_prediction.cpp
 *
 * @brief SIFT predictions source file
 *
 * @author: pauline-ng, angsm
 */

#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <cstring>
#include <ostream>
#include <iostream>
#include <fstream>
#include <regex>
#include <iomanip>
#include <algorithm>

#include "utils.hpp"
#include "constants.hpp"
#include "sift_scores.hpp"

constexpr uint32_t kMaxSequences = 400;
constexpr double TOLERANCE_PROB_THRESHOLD = 0.05;
constexpr double ADEQUATE_SEQ_INFO =3.25;
// init_frq_qij
constexpr double LOCAL_QIJ_RTOT = 5.0;

void scale_matrix_to_max_aa(std::vector<std::vector <double>>& matrix, std::vector<int>& max_aa_array) {

    int query_length = matrix.size();
    int aas = matrix[0].size();

    for (int pos = 0; pos < query_length; pos++) {
        int max_aa = max_aa_array[pos];
        double max_aa_score = matrix[pos][max_aa];
        for (int aa_index = 0; aa_index < aas; aa_index++) {
            matrix[pos][aa_index] = matrix[pos][aa_index] / max_aa_score;
        }
    }
}

void find_max_aa_in_matrix(std::vector<std::vector<double>>& matrix, std::vector<int>& max_aa_index) {
    int query_length = matrix.size();
    int aas = matrix[0].size();

    for (int pos = 0; pos < query_length; pos++) {
        int max_aa = -1;
        double max_count = -1.0;
        for (int aa_index =0; aa_index < aas; aa_index++) {
            if (matrix[pos][aa_index] > max_count) {
                max_aa = aa_index;
                max_count = matrix[pos][aa_index];
            }
        }
        max_aa_index[pos] = max_aa;
    }
}

void calcEpsilon(std::vector<std::vector<double>>& weighted_matrix, std::vector<int>& max_aa_array,
    std::vector<double>& number_of_diff_aas, std::vector<double>& epsilon) {

    int query_length = weighted_matrix.size();
    /*vector <int> max_aa_array (query_length);
    find_max_aa_in_matrix (weighted_matrix, max_aa_array);*/

    for (int pos =0; pos < query_length; pos++) {
        if (number_of_diff_aas[pos] == 1) {
            epsilon[pos] = 0;
        } else {
            int max_aa = max_aa_array[pos];
            double sum = 0.0;
            double pos_tot = 0.0;
            for (char aa = 'A'; aa <= 'Z'; aa++) {
                if (valid_amino_acid(aa)) {
                    int aa_index = int(aa) - int('A');
                    int rank = rank_matrix[max_aa][aa_index];
                    sum += (double) rank * weighted_matrix[pos][aa_index];
                    pos_tot +=  weighted_matrix[pos][aa_index];
                }
            }
            sum = sum / pos_tot;
            epsilon[pos] = exp((double) sum);
        }
    }
}

void readSubstFile(char* substFilename, std::list<std::string>& substList) {

    std::ifstream infile(substFilename);
    std::string line;

    if (infile.good()) {
        while (std::getline(infile, line)) {
            substList.push_back(line);
        }
    }
    infile.close();
}

void addMedianSeqInfo(std::vector<Chain*>& alignment_string, Chain* query, std::vector<std::vector<double>>& matrix,
    std::unordered_map<std::string,double>& medianSeqInfoForPos) {

    int query_length = matrix.size();
    std::vector<double> number_of_diff_aas(query_length);

    for (auto it = medianSeqInfoForPos.begin(); it != medianSeqInfoForPos.end(); ++it) {
        int pos = stoi(it->first); // position in protein, start count at 1
        pos = pos -1;  // position in array
        if (it->second == -1) {
            /* get sequences that don't have X or invalid aa */
            std::unordered_map <int, int> invalidSeq;
            for (int i = 0; i < (int) alignment_string.size(); i++) {
                char c = chainGetChar(alignment_string[i], pos);
                if (!valid_amino_acid(c)) {
                    invalidSeq.insert({i,1}); //invalidSeq[i] = 1;
                }
            }
            std::vector<Chain*> alignment_with_noX_at_pos;
            seqs_without_X(alignment_string, pos, alignment_with_noX_at_pos );

            int num_seqs_in_alignment = alignment_with_noX_at_pos.size();
            if (num_seqs_in_alignment == 0) {
                medianSeqInfoForPos[it->first] = 0.0;
                continue;
            }

            /* have to recalculate sequence weights */
            std::vector<double> weights_1 (num_seqs_in_alignment, 1.0);
            /* raw counts of valid amino acids */
            std::vector<double> aas_stored_at_each_pos (query_length);
            std::vector<std::vector<double>> matrix_noX_raw(query_length, std::vector<double>(26, 0.0));

            createMatrix(alignment_with_noX_at_pos, query, weights_1, matrix_noX_raw, aas_stored_at_each_pos);
                std::vector <double> seq_weights (num_seqs_in_alignment);
            /* calcSeqWeights is the only function where all sequences
            and all positions have to be looked at at the same time.
            all other routines treat each position independently */
            calcSeqWeights (alignment_with_noX_at_pos,matrix_noX_raw, aas_stored_at_each_pos, seq_weights,
                number_of_diff_aas);
            std::vector<std::vector<double>> matrix_noX (query_length, std::vector<double>(26, 0.0));

            basic_matrix_construction(alignment_with_noX_at_pos, seq_weights, matrix_noX);
            // printMatrix (matrix, "tmp_basicmatrix.txt");
            double medianSeqInfo = calculateMedianSeqInfo(alignment_with_noX_at_pos, invalidSeq, seq_weights, matrix_noX);
            medianSeqInfoForPos[it->first] = medianSeqInfo;
        }
    }
}

double calculateMedianSeqInfo(std::vector<Chain*>& alignment_strings, std::unordered_map<int, int>& invalidSeq,
    std::vector<double>& seq_weights, std::vector<std::vector<double>>& matrix) {

    int query_len = chainGetLength(alignment_strings[0]);
    int amino_acid_num = 26;
    float median = kLog_2_20;
    double tmp = 0.0;

    double* amino_acid_nums = new double[amino_acid_num];
    for (int i = 0; i < amino_acid_num; ++i) {
        amino_acid_nums[i] = 0.0;
    }

    float* pos_freq = new float[query_len];
    for (int i = 0; i < query_len; ++i) {
        pos_freq[i] = 0.0;
    }

    // char c;
    // double valid_aa;
    for (int pos_index = 0; pos_index < query_len; pos_index++) {
        // valid_aa = 0.0;
        // re-initialize to 0
        double total_weight = 0.0;
        double r = 0.0;
        for (char aa = 'A'; aa <= 'Z'; aa++) {
            if (valid_amino_acid(aa)) {
                int idx = aa_to_idx(aa);
                amino_acid_nums[idx] = 0.0;
                total_weight += matrix[pos_index][idx];
            }
        }

        for (char aa= 'A';aa <= 'Z'; aa++) {
            int idx = aa_to_idx(aa);
            tmp = matrix[pos_index][idx]/total_weight;
            if (tmp > 0.0 && valid_amino_acid(aa)) {
                r += tmp * log(tmp);
            }
        }
        r = r / log(2.0);
        pos_freq[pos_index] = r + kLog_2_20;
    } /* end pos_index */
    median = getMedian(pos_freq, query_len);

    delete[] pos_freq;
    delete[] amino_acid_nums;

    return median;
}

bool checkSubsts(const std::list<std::string>& substList, Chain* query) {

    std::regex regexSubst("^([A-Z])([0-9]+)([A-Z])");  /*, std::regex_constants::basic); */
    std::smatch m;

    bool valid = true;
    uint32_t num_valid_lines = 0;
    for (std::list<std::string>::const_iterator it = substList.begin(); it != substList.end(); ++it) {
        if (regex_search(*it, m, regexSubst)) {
            ++num_valid_lines;
            char ref_aa = std::string(m[1])[0];
            int pos = stoi(std::string(m[2])) - 1;
            if (pos >= chainGetLength(query)) {
                fprintf(stderr, "* skipping prediction for [ %s ]: substitution list has a position out of bounds (line: %s, query length = %d) *\n",
                    chainGetName(query), (*it).c_str(), chainGetLength(query));
                valid = false;
                break;
            }
            if (chainGetChar(query, pos) != ref_aa) {
                fprintf(stderr, "* skipping prediction for [ %s ]: substitution list assumes wrong amino acid at position %d (line: %s, query amino acid = %c) *\n",
                    chainGetName(query), pos+1, (*it).c_str(), chainGetChar(query, pos));
                valid = false;
                break;
            }
        }
    }

    if (num_valid_lines == 0) {
        fprintf(stderr, "* skipping prediction for [ %s ]: substitution list contains zero valid lines *\n", chainGetName(query));
        valid = false;
    }

    return valid;
}

void hashPredictedPos(std::list<std::string>& substList, std::unordered_map<std::string, double>& medianSeqInfoForPos) {

    std::list<std::string>::const_iterator iterator;

    std::regex regexSubst ("^[A-Z]([0-9]+)[A-Z]");  /*, std::regex_constants::basic); */
    std::smatch m;

    for (iterator = substList.begin(); iterator != substList.end(); ++iterator)     {
        std::string substLine = *iterator;
        if (regex_search(substLine, m, regexSubst)) {
            std::string string_pos = m[1];
            medianSeqInfoForPos[string_pos] = -1;
        }
    }
}

void addPosWithDelRef(Chain* query, std::vector<std::vector<double>>& SIFTscores,
    std::unordered_map<std::string, double>& medianSeqInfoForPos) {

    int query_length = SIFTscores.size();

    for (int pos = 0; pos < query_length; pos++) {
        char ref_aa = chainGetChar(query, pos);
        int ref_aa_index = (int) ref_aa - (int) 'A';

        if (SIFTscores[pos][ref_aa_index] < TOLERANCE_PROB_THRESHOLD) {
            medianSeqInfoForPos[std::to_string(pos+1)] = -1; // pos + 1 because the substitution files are 1-based
        }
    }
}

void check_refaa_against_query(char ref_aa, int aa_pos, Chain* query, std::ofstream& outfp) {
    char aa = chainGetChar(query, aa_pos);
    if (aa != ref_aa) {
        outfp << "WARNING! Amino acid " << aa << " is at position " <<
            std::to_string(aa_pos + 1) << ", but your list of substitutions assumes it's a " << ref_aa << std::endl;
    }
}

std::string print_double(double num, int precision) {
    std::stringstream stream;
    stream << std::fixed << std::setprecision(precision) << num;
    return stream.str();
}

void printSubstFile(const std::list<std::string>& substList, std::unordered_map<std::string, double>& medianSeqInfoForPos,
    const std::vector<std::vector<double>>& SIFTscores, const std::vector<double>& aas_stored,
    const int total_seq, Chain* query, const std::string outfile) {

    std::list<std::string>::const_iterator iterator;
    std::regex regexSubst ("^([A-Z])([0-9]+)([A-Z])");  /*, std::regex_constants::basic); */
    std::smatch m;
    std::ofstream outfp;
    outfp.open(outfile, std::ios::out);
    int query_length = SIFTscores.size();

    for (int pos = 0; pos < query_length; pos++) {
        char ref_aa = chainGetChar(query, pos);
        int ref_aa_index = (int) ref_aa - (int) 'A';

        if (SIFTscores[pos][ref_aa_index] < TOLERANCE_PROB_THRESHOLD) {
            auto search = medianSeqInfoForPos.find(std::to_string(pos));
            double median = search->second;
            if (median < ADEQUATE_SEQ_INFO) {
                outfp << "WARNING! " << ref_aa << std::to_string(pos+1) << " not allowed! score: " <<
                    print_double( SIFTscores[pos][ref_aa_index],2) << " median: " <<
                    print_double(medianSeqInfoForPos[std::to_string(pos)],2) << " # of sequence: " <<
                    std::to_string((int) aas_stored[pos]) << std::endl;
            }
        } /* end if less than TOLERANCE_PROB_THRESHOLD */
    }

    for (iterator = substList.begin(); iterator != substList.end(); ++iterator) 	{
        std::string substLine = *iterator;
        std::stringstream ss(std::stringstream::in | std::stringstream::out);
        ss << substLine;
        std::string cleanSubst;
        ss >> cleanSubst;

        if (regex_search(substLine, m, regexSubst)) {
            char tmp[10];
            std::string tmp_string = m[1];
            std::strncpy(tmp, tmp_string.c_str(), tmp_string.size()+1);
            char ref_aa = tmp[0]; //.c_str()[0];
            std::string aa_pos_string = m[2];
            int aa_pos = stoi ( aa_pos_string); //.c_str());
            aa_pos = aa_pos -1; // position in array
            tmp_string = m[3];
            std::strncpy(tmp, tmp_string.c_str(), tmp_string.size()+1);
            char new_aa = tmp[0]; //.c_str()[0];

            int new_aa_index = (int) new_aa - (int) 'A';
            double score = SIFTscores[aa_pos][new_aa_index];

            check_refaa_against_query(ref_aa, aa_pos, query, outfp);
            outfp << cleanSubst << "\t";
            if (score >= TOLERANCE_PROB_THRESHOLD) {
                outfp << "TOLERATED\t" << print_double(score,2);
            } else {
                outfp << "DELETERIOUS\t" << print_double(score,2);
            }
            auto search = medianSeqInfoForPos.find(aa_pos_string);
            outfp << "\t" << print_double(search->second,2) << "\t" <<
                std::to_string((int) aas_stored[aa_pos]) << "\t" << std::to_string(total_seq) << std::endl;
        }
    }

    outfp.close();
}

bool valid_amino_acid(char aa) {
    if (aa == 'B' or aa == 'Z' or aa == 'J' or aa == 'O' or aa == 'U' or aa == 'X' or aa == '-' or aa == '*') {
        return false;
    } else {
        return true;
    }
}

void calcSIFTScores(std::vector<Chain*>& alignment_string, Chain* query, std::vector<std::vector<double>>& matrix,
    std::vector<std::vector<double>>& SIFTscores) {

    int query_length = matrix.size();

    std::vector<int> max_aa_array(query_length);
    std::vector<double> number_of_diff_aas(query_length);
    std::vector<double> epsilon(query_length);

    int num_seqs_in_alignment = alignment_string.size();
    std::vector<std::vector<double>> raw_count_matrix(query_length, std::vector<double>(26, 0.0));

    /* initialize an array with weight 1 */
    std::vector<double> weights_1(num_seqs_in_alignment, 1.0);

    /* raw counts of valid amino acids */
    std::vector<double> aas_stored_at_each_pos(query_length);
    createMatrix(alignment_string, query, weights_1, raw_count_matrix, aas_stored_at_each_pos);

    std::vector <double> seq_weights (num_seqs_in_alignment);

    /* calcSeqWeights is the only function where all sequences
    and all positions have to be looked at at the same time.
    all other routines treat each position independently */
    calcSeqWeights (alignment_string, matrix, aas_stored_at_each_pos, seq_weights,
        number_of_diff_aas);

    /* now construct matrix with weighted sequence values */
    std::vector<std::vector<double>> seq_weighted_matrix(query_length, std::vector<double>(26,0.0));
    std::vector<double> tot_weights_each_pos(query_length);

    createMatrix(alignment_string, query, seq_weights, seq_weighted_matrix, tot_weights_each_pos);

    find_max_aa_in_matrix(seq_weighted_matrix, max_aa_array);

    calcEpsilon(seq_weighted_matrix, max_aa_array, number_of_diff_aas, epsilon );

    /* pseudo_diri */
    std::vector<std::vector<double>> diric_matrix(query_length, std::vector<double>(26, 0.0));
    calcDiri(seq_weighted_matrix, diric_matrix);
    //printMatrix (diric_matrix, "diri_matrix.txt");

    for (int pos= 0; pos < query_length; pos++) {
        for (char aa = 'A'; aa <= 'Z'; aa++) {
            int aa_index = (int) aa - (int) 'A';
            SIFTscores[pos][aa_index] = seq_weighted_matrix[pos][aa_index] + epsilon[pos] * diric_matrix[pos][aa_index];
            SIFTscores[pos][aa_index] /= (tot_weights_each_pos[pos] + epsilon[pos]);
        }
    }
    /* have to find max aa again because it'schanged with newly added weights */
    find_max_aa_in_matrix(SIFTscores, max_aa_array);
    scale_matrix_to_max_aa(SIFTscores, max_aa_array);
    // printMatrix(SIFTscores, "siftscores_matrix.txt");
}

void calcDiri(std::vector<std::vector<double>>& count_matrix, std::vector<std::vector<double>>& diric_matrix) {

    int query_length = count_matrix.size();
    for (int pos = 0; pos < query_length; pos++) {
        add_diric_values(count_matrix[pos], diric_matrix[pos]);
    }
}

double add_logs (double logx, double logy) {
    if (logx > logy) {
        return (logx + log (1.0 + exp (logy -logx)));
    } else {
        return (logy + log (1.0 + exp (logx - logy)));
    }
}

void add_diric_values(std::vector<double>& count_col, std::vector<double>& diric_col) {

    double tmp;
    int diri_comp_num = (int) (sizeof (diri_altot)/sizeof (diri_altot[0]));
    std::vector<double> probn(diri_comp_num, 0.0);
    std::vector<double> probj(diri_comp_num,0.0);

    double pos_count_tot = 0.0;
    for (int j=0; j < (int) count_col.size(); j++) {
        pos_count_tot += count_col[j];
    }

    /*-----------   compute equation (3), Prob(n|j) ------------  */
    for (int j = 0; j < diri_comp_num; j++) {
        probn[j] = lgamma(pos_count_tot+1.0) + lgamma(diri_altot[j]);
        probn[j] -= lgamma(pos_count_tot + diri_altot[j]);

        for (char aa = 'A'; aa <= 'Z'; aa++) {
            int aa_index = (int) aa - (int) 'A';
            if (valid_amino_acid(aa)) {
                tmp = lgamma(count_col[aa_index] + diri_alpha[j][aa_index]);
                tmp -= lgamma(count_col[aa_index] + 1.0);
                tmp -= lgamma(diri_alpha[j][aa_index]);
                probn[j] += tmp;
            } /* end if greater than 0 */
        } /* end all amino acids */
    } /* end for all Diri components */

    /*------ compute sum qk * p(n|k) using logs & exponents ----------*/
    double denom = log(diri_q[0]) + probn[0];

    for (int j =1; j < diri_comp_num; j++) {
        double tmp = log(diri_q[j]) + probn[j];
        denom = add_logs(denom, tmp);
    }

    /*   compute equation (3), Prob(j|n)  */
    for (int j = 0; j < diri_comp_num; j++) {
        probj[j] = log(diri_q[j]) + probn[j] - denom;
    }

    double totreg = 0.0;
    for (char aa = 'A'; aa <= 'Z'; aa++) {
        if (valid_amino_acid(aa)) {
            int aa_index = (int) aa - (int) 'A';
            for (int j =0; j < (int) probj.size(); j++) {
                diric_col[aa_index] += exp (probj[j]) * diri_alpha[j][aa_index];
            } /* gone through all components */
            totreg += diric_col[aa_index];
        } /* end valid amino acid */
    }
    /* now normalize */
    for (char aa= 'A'; aa <= 'Z'; aa++) {
        int aa_index = (int) aa - (int) 'A';
        diric_col[aa_index] /= totreg;
    }
}

void calcSeqWeights(const std::vector<Chain*>& alignment_string, std::vector<std::vector<double>>& matrix,
    std::vector<double>& amino_acids_present, std::vector<double>& seq_weights, std::vector<double>& number_of_diff_aas) {

    int query_length = matrix.size();
    /* first calculate number of different aas at each position */
    /* std::vector <int> number_of_diff_aas (query_length); */

    /* initialize */
    for (int pos = 0; pos < query_length; pos++) {
        number_of_diff_aas[pos] = 0;
    }
    for (uint32_t seq_index = 0; seq_index < alignment_string.size(); seq_index++) {
        seq_weights[seq_index] = 0.0;
    }

    /* now tabulate # of unique amino acids at each position */
    for (int pos = 0; pos < query_length; pos++) {
        for (char aa = 'A'; aa <= 'Z'; aa++) {
            int aa_index = int (aa) - int ('A');
            if (valid_amino_acid(aa) &&  matrix[pos][aa_index] > 0.0) {
                number_of_diff_aas[pos] += 1.0f;
            }
        }
    }

    double tot = 0.0;
    /* now calculate position-based weights */
    for (uint32_t seq_index = 0; seq_index < alignment_string.size(); seq_index++) {
        for (int pos = 0; pos < query_length; pos++) {
            char aa = chainGetChar(alignment_string[seq_index], pos);
            int aa_index = (int) aa - (int) 'A';
            if (valid_amino_acid(aa) &&  matrix[pos][aa_index] > 0.0) {
                double tmp = number_of_diff_aas[pos] * matrix[pos][aa_index];
                seq_weights[seq_index] += 1.0/tmp;
            }
        }
        tot += seq_weights[seq_index];
    }

    /* normalize so weights sum up to the number of sequences */
    double new_tot_weight = 0.0;
    for (uint32_t seq_index = 0; seq_index<alignment_string.size(); seq_index++) {
        seq_weights[seq_index] = seq_weights[seq_index] / tot * alignment_string.size()  ;
        new_tot_weight += seq_weights[seq_index];
    }
}

void remove_seqs_percent_identical_to_query(Chain* queries, std::vector<Chain*>& alignment_string, double seq_identity) {

    double identity;
    double seqTotal;
    int lenOfQuery = chainGetLength(queries);
    int currPos = 0;

    /*iterate through elements in vector of chains*/
    while (currPos < int(alignment_string.size())){

        identity = 0;
        seqTotal = 0;
        int lenOfAlign = chainGetLength(alignment_string[currPos]);

        /*see if alignment_string is indeed made to match query seq length*/
        if (lenOfQuery == lenOfAlign){
            for (int m = 0; m < lenOfAlign; m++){
                /*compare each char of query to alignment_string char at same position in chain*/
                char qChar = chainGetChar(queries, m);
                char aChar = chainGetChar(alignment_string[currPos], m);

                if (aChar != 'X'){
                    if (valid_amino_acid(aChar) && valid_amino_acid(qChar)){
                        seqTotal++;
                        if (qChar == aChar){
                            identity++;
                        }
                    }
                }
            }
        } else {
            std::cout << "Length does not match!" << std::endl;
            exit(1);
        }

        double perc_similar = (identity / seqTotal) * 100;
        /*Read curr element, delete if beyond threshold, else move to next pos*/
        if (perc_similar >= seq_identity) {
            chainDelete(alignment_string[currPos]);
            alignment_string.erase(alignment_string.begin() + currPos);
        } else {
            currPos++;
        }
    }
}

void printSeqNames (std::vector<Chain*>& alignment_string) {
    for (int i=0; i < int(alignment_string.size()); i++){
        std::cout << "printName in here " << std::to_string(i) << std::endl;
        std::cout << "printName" << chainGetName(alignment_string[i]) << std::endl;
    }
}

/* add it with seq_weights 1 or pb weights
with 1 it's raw count, with pb_weights it's weighted */
void createMatrix(const std::vector<Chain*>& alignment_string, Chain* query,
    const std::vector<double>& seq_weights, std::vector<std::vector<double>>& matrix,
    std::vector<double>& tot_pos_weight) {

    int query_len = chainGetLength(query);
    for (uint32_t seq_index = 0; seq_index < alignment_string.size(); seq_index++) {
        for (int pos = 0; pos < query_len; pos++) {
            char aa = chainGetChar(alignment_string[seq_index], pos);
            if (valid_amino_acid(aa)) {
                int aa_index = (int) aa - (int) 'A';
                matrix[pos][aa_index] += seq_weights[seq_index] ; /*seq_weight[seq_index];*/
                tot_pos_weight[pos] += seq_weights[seq_index];
            }
        }
    }
}

void printMatrix(std::vector<std::vector<double>>& matrix, std::string filename) {

    int query_length = matrix.size();
    int aas = matrix[0].size();
    std::ofstream out_file;

    out_file.open(filename);

    /* print out header */
    for (int aa_index = 0; aa_index< aas; aa_index++) {
        out_file << "\t";
        out_file << (char) (aa_index + (int) 'A');
    }
    out_file << std::endl;

    for (int pos = 0; pos < query_length; pos++) {
        out_file << (pos + 1);
        for (int aa_index =0; aa_index < aas; aa_index++) {
            out_file << "\t" << matrix[pos][aa_index];
        }
        out_file << std::endl;
    }
    out_file.close();
}

void printMatrixOriginalFormat(std::vector<std::vector<double>>& matrix, std::string filename) {

    auto fp = fopen(filename.c_str(), "w");

    int query_length = matrix.size();
    int aas = matrix[0].size();

    // print out header
    fprintf(fp, "ID   UNK_ID; MATRIX\nAC   UNK_AC\nDE   UNK_DE\nMA   UNK_BL\n");
    fprintf(fp, " ");
    for (int aa_index = 0; aa_index < aas; aa_index++) {
        // ignore J O U
        if (aa_index != 9 && aa_index != 14 && aa_index != 20) {
            fprintf(fp, " %c  ", aa_index + 'A');
        }
    }
    fprintf(fp, " *   -\n");

    for (int pos = 0; pos < query_length; pos++) {
        for (int aa_index = 0; aa_index < aas; aa_index++) {
            if (aa_index != 9 && aa_index != 14 && aa_index != 20) {
                fprintf(fp, " %6.4f ", matrix[pos][aa_index]);
            }
        }
        fprintf(fp, " %6.4f  %6.4f\n", 0.0, 0.0);
    }
    fprintf(fp, "//\n");

    fclose(fp);
}

int aa_to_idx(char character){
    // A = 0 and etc
    return int(character) - int('A');
}

void basic_matrix_construction(const std::vector<Chain*>& alignment_string, const std::vector<double>& seq_weights,
    std::vector<std::vector<double>>& matrix) {

    // do partial calculation on aa that can be represented by B or Z as well
    double part_D = aa_frequency[aa_to_idx('D')] / ( aa_frequency[aa_to_idx('D')] + aa_frequency[aa_to_idx('N')] );
    double part_N = aa_frequency[aa_to_idx('N')] / ( aa_frequency[aa_to_idx('D')] + aa_frequency[aa_to_idx('N')] );
    double part_E = aa_frequency[aa_to_idx('E')] / ( aa_frequency[aa_to_idx('E')] + aa_frequency[aa_to_idx('Q')] );
    double part_Q = aa_frequency[aa_to_idx('Q')] / ( aa_frequency[aa_to_idx('E')] + aa_frequency[aa_to_idx('Q')] );

    int proteinLen = chainGetLength(alignment_string[0]);
    int aaLen =  matrix[0].size();

    // loop every position in chain length
    for (int pos = 0; pos < proteinLen; pos++){
        double total = 0.0;

        // read pos in each chain
        for (int seq = 0; seq < int(alignment_string.size()); seq++){
            char currChar = chainGetChar(alignment_string[seq],pos);

            // parition B between D and N
            if (currChar == 'B') {
                if (aa_frequency[aa_to_idx('D')] != 0.0) { /* try to protect div by zero */
                    double num = (part_D * seq_weights[seq]) / aa_frequency[aa_to_idx('D')];
                    matrix[pos][aa_to_idx('D')] += num;
                    total += num;
                }

                if (aa_frequency[aa_to_idx('N')] != 0.0) { /* try to protect div by zero */
                    double num = (part_N * seq_weights[seq]) / aa_frequency[aa_to_idx('N')];
                    matrix[pos][aa_to_idx('N')] += num;
                    total += num;
                }
                /* partition Z between E and Q */
            } else if (currChar == 'Z') {
                if (aa_frequency[aa_to_idx('E')] != 0.0) { /* try to protect div by zero */
                    double num = (part_E * seq_weights[seq]) / aa_frequency[aa_to_idx('E')];
                    matrix[pos][aa_to_idx('E')] += num;
                    total += num;
                }

                if (aa_frequency[aa_to_idx('Q')] != 0.0) { /* try to protect div by zero */
                    double num = (part_Q * seq_weights[seq]) / aa_frequency[aa_to_idx('Q')];
                    matrix[pos][aa_to_idx('Q')] += num;
                    total += num;
                }

            } else {
                if (aa_frequency[aa_to_idx(currChar)] != 0.0){
                    if (currChar != 'X' && currChar != '-' && currChar != '*'){
                        double num = seq_weights[seq] / aa_frequency[aa_to_idx(currChar)];
                        matrix[pos][aa_to_idx(currChar)] += num;
                        total += num;
                    }
                }
            }
        }

        // loop to calculate percentage, excludes index 26 (-) and index 27 (*)
        for (int n = 0; n < aaLen; n++){
            // A to Z excluding 'X' (alignment filler)
            if (n <= aa_to_idx('Z') || n != aa_to_idx('X')) {
                matrix[pos][n] = matrix[pos][n] * 100.0 / total; // X and * and -
            } else {
                matrix[pos][n] = aa_frequency[n];
            }
        }

        // lastly, tali B and Z rows
        matrix[pos][aa_to_idx('B')] = (matrix[pos][aa_to_idx('D')] * part_D) + (matrix[pos][aa_to_idx('N')] * part_N);
        matrix[pos][aa_to_idx('Z')] = (matrix[pos][aa_to_idx('E')] * part_E) + (matrix[pos][aa_to_idx('Q')] * part_Q);
    }
}

void seqs_without_X(const std::vector<Chain*>& alignment_string, int pos, std::vector<Chain*>& out_alignment) {

    char posChar;
    /*iterate vector of aligment and check pos*/
    for( int n = 0; n < int(alignment_string.size()); n++){
        posChar = chainGetChar(alignment_string[n], pos);

        /*check if char in pos is valid, and push it into new vector*/
        if (valid_amino_acid(posChar)){
            out_alignment.push_back(alignment_string[n]);
        }
    }
}
