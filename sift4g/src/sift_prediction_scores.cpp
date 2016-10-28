/*!
 * @file sift_prediction.cpp
 *
 * @brief SIFT predictions source file
 */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <ostream>
#include <iostream>
#include <fstream>
#include <vector>
#include <list>
#include <string>
#include <regex>
#include <unordered_map>
#include <iomanip>
#include "utils.hpp"
#include "sift_prediction_scores.hpp"
#include "select_alignments.hpp"
using namespace std;


constexpr uint32_t kMaxSequences = 400;
constexpr double TOLERANCE_PROB_THRESHOLD = 0.05; 
constexpr double  ADEQUATE_SEQ_INFO =3.25;
// init_frq_qij

constexpr double LOCAL_QIJ_RTOT = 5.0; 
 
void scale_matrix_to_max_aa (std::vector <std::vector <double>>& matrix, std::vector <int> max_aa_array) {
	int query_length = matrix.size();
        int aas = matrix[0].size();

	for (int pos = 0; pos < query_length; pos++) {
		int max_aa = max_aa_array[pos];
                double max_aa_score = matrix[pos][max_aa];
                for (int aa_index =0; aa_index < aas; aa_index++) {
			matrix[pos][aa_index] = matrix[pos][aa_index]/max_aa_score;
		}
	}
}

void find_max_aa_in_matrix (std::vector <std::vector <double>> matrix, std::vector <int>& max_aa_index) {
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

void calcEpsilon (std::vector <std::vector <double>>  weighted_matrix, std::vector <int> max_aa_array, std::vector <double>& number_of_diff_aas, std::vector <double>& epsilon  ) {

	int query_length = weighted_matrix.size();
	/*       vector <int> max_aa_array (query_length);
		 find_max_aa_in_matrix (weighted_matrix, max_aa_array);
	 */
	for (int pos =0; pos < query_length; pos++) {
		if (number_of_diff_aas[pos] == 1) {

			epsilon[pos] = 0;
		} else {
			int max_aa = max_aa_array[pos];
			double sum = 0.0;
			double pos_tot = 0.0;
			for (char aa = 'A'; aa <= 'Z'; aa++) {
				if (valid_amino_acid (aa)) {
					int aa_index = int (aa) - int ('A');
					int rank = rank_matrix[max_aa][aa_index];
					sum += (double) rank * weighted_matrix[pos][aa_index];
					pos_tot +=  weighted_matrix[pos][aa_index];
				}
			}
			sum = sum / pos_tot;
			epsilon[pos] = exp ((double) sum);
		}
		//cout << "Epsilon pos: " << to_string (pos) << " ep " << to_string (epsilon[pos]) << endl;
	}
}

void readSubstFile (char* substFilename, std::list <std::string> &substList) {

	std::ifstream infile (substFilename);
	std::string line;

	if (infile.good()) {
		while (std::getline (infile, line))
		{
			substList.push_back (line);
		}
	}
	infile.close();
}

void addMedianSeqInfo (std::vector<Chain*> alignment_string, Chain* query, std::vector<std::vector<double>> matrix, std::unordered_map <std::string, double>& medianSeqInfoForPos) {

        int query_length = matrix.size();
	vector <double> number_of_diff_aas (query_length);

	for (auto it = medianSeqInfoForPos.begin(); it != medianSeqInfoForPos.end(); ++it) {
		//cout << "pos " << it->first << endl;
		int pos = stoi( it->first); // position in protein, start count at 1
		pos = pos -1;  // position in array
		if (it->second == -1) { 

			/* get sequences that don't have X or invalid aa */
			std::unordered_map <int, int> invalidSeq;
			for (int i = 0; i < (int) alignment_string.size(); i++) {
				char c = chainGetChar (alignment_string[i], pos);
				if (!valid_amino_acid (c)) {
					invalidSeq.insert({i,1});//invalidSeq[i] = 1;
				}
			}  
			std::vector <Chain*> alignment_with_noX_at_pos;
			seqs_without_X  (alignment_string, pos, alignment_with_noX_at_pos );

			int num_seqs_in_alignment = alignment_with_noX_at_pos.size();
			//cout << "at pos " << to_string (pos ) << "now this number of seqs " << to_string (num_seqs_in_alignment) << endl;		
		        /* have to recalculate sequence weights */
			std::vector <double> weights_1 (num_seqs_in_alignment, 1.0);
		        /* raw counts of valid amino acids */
		        std::vector<double> aas_stored_at_each_pos (query_length);
		        std::vector<std::vector<double>>
					 matrix_noX_raw (query_length,
			                std::vector <double> (26, 0.0));

		        createMatrix (alignment_with_noX_at_pos, query, weights_1, matrix_noX_raw, aas_stored_at_each_pos);
        		std::vector <double> seq_weights (num_seqs_in_alignment);
        	/* calcSeqWeights is the only function where all sequences 
	        and all positions have to be looked at at the same time. 
       		all other routines treat each position independently */
       			calcSeqWeights (alignment_with_noX_at_pos,matrix_noX_raw, aas_stored_at_each_pos, seq_weights,
		          number_of_diff_aas);
			std::vector<std::vector<double>>
                                         matrix_noX (query_length,
                                        std::vector <double> (26, 0.0));
/*			for (int seq_index = 0; seq_index < num_seqs_in_alignment; seq_index++) {
				cout << "newSIFT seq weight " << to_string (seq_index) << " " << seq_weights[seq_index] << endl;
				for (int j = 0; j < query_length ; j++) {
            				cout << chainGetChar(alignment_with_noX_at_pos[seq_index], j);
        			}

			}
*/
			basic_matrix_construction (alignment_with_noX_at_pos, seq_weights, matrix_noX); 
//			printMatrix (matrix, "tmp_basicmatrix.txt");	
			double medianSeqInfo = calculateMedianSeqInfo (alignment_with_noX_at_pos, invalidSeq, seq_weights, matrix_noX);
			medianSeqInfoForPos[it->first] = medianSeqInfo;
			//cout << it->first << " with median " << std::to_string (medianSeqInfo) << endl;
		}
	}
	//cout << "done with median seq info " << endl;
}

double  calculateMedianSeqInfo (std::vector<Chain*>& alignment_strings, std::unordered_map <int, int> invalidSeq, std::vector <double> seq_weights, std::vector<std::vector<double>> matrix) {

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

	/*char c;
	double  valid_aa;
*/
	
	for (int pos_index = 0; pos_index < query_len; pos_index++) {
		/*valid_aa = 0.0; */
		/* re-initialize to 0 */
		double total_weight = 0.0;
		double r = 0.0;
		for (char aa = 'A'; aa <= 'Z'; aa++) {
			if (valid_amino_acid (aa)) {
				int idx = aa_to_idx (aa);
				amino_acid_nums[idx] = 0.0;
				total_weight += matrix[pos_index][idx];
			}
		}


		for (char aa= 'A';aa <= 'Z'; aa++) {
			int idx = aa_to_idx (aa);
			tmp = matrix[pos_index][idx]/total_weight;
			if (tmp > 0.0 && valid_amino_acid (aa)) {
				r += tmp * log (tmp);
			//cout << "newSIFT: freq " << to_string (aa) << " " << to_string (pos_index) << " weight " << to_string (matrix[pos_index][idx]) << " tmp " << to_string(tmp) << " total " << to_string(total_weight) << endl;
			}
		}
		r = r / log (2.0);
		
/*for (int seq_index = 0; seq_index < (int) alignment_strings.size(); seq_index++) {
			if (invalidSeq.find(seq_index) == invalidSeq.end()) {
				c = chainGetChar(alignment_strings[seq_index], pos_index);
				if (valid_amino_acid (c)) {
//					valid_aa++;
//					amino_acid_nums[(int) c - 'A']++;
					valid_aa += seq_weights[seq_index]; 
					amino_acid_nums[(int) c - 'A'] += seq_weights[seq_index];
//					cout << "adding seq weight " << to_string (seq_weights[seq_index]) << endl; 
				}
			}
		} */ /* end seq_index */
/*		for (int k = 0; k < amino_acid_num; ++k) {
			if (amino_acid_nums[k] != 0) {
				float f  = amino_acid_nums[k] / (float) valid_aa;
				cout << "newSIFT: freq " << to_string (k) << " " << to_string (pos_index) << " weight " << amino_acid_nums[k] << " freq " << f << " total " << valid_aa << endl; 
				pos_freq[pos_index] += f * log2f (f); 
			}
		}
*/
/*		pos_freq[pos_index] += kLog_2_20; */
		pos_freq[pos_index] = r + kLog_2_20;;
//	      cout << "newSIFT: pos " << std::to_string (pos_index)  << " string " << std::to_string (pos_freq[pos_index]) << endl; 
	} /* end pos_index */
	median = getMedian(pos_freq, query_len);


	delete[] pos_freq;
	delete[] amino_acid_nums;

	return median;

} 


void hashPredictedPos (std::list <std::string> &substList, std::unordered_map <string, double>& medianSeqInfoForPos) {

	std::list<std::string>::const_iterator iterator;

	std::regex regexSubst ("^[A-Z]([0-9]+)[A-Z]");  /*, std::regex_constants::basic); */ 
	std::smatch m;

	for (iterator = substList.begin(); iterator != substList.end(); ++iterator)     {
		string substLine = *iterator;
		if (regex_search (substLine, m, regexSubst)) {
			string string_pos = m[1];
						
			medianSeqInfoForPos[string_pos] = -1;
		}
	}

}

void addPosWithDelRef (Chain* query,  std::vector<std::vector<double>> SIFTscores, std::unordered_map <string, double>& medianSeqInfoForPos) {

	int query_length = SIFTscores.size();

        for (int pos = 0; pos < query_length; pos++) {
             char ref_aa = chainGetChar (query, pos);
             int ref_aa_index = (int) ref_aa - (int) 'A';

             if (SIFTscores[pos][ref_aa_index] < TOLERANCE_PROB_THRESHOLD) {
                     medianSeqInfoForPos[to_string(pos)] = -1;
             }
        }
	
}

void check_refaa_against_query (char ref_aa, int aa_pos, Chain* query, std::ofstream& outfp) {
	char aa = chainGetChar (query, aa_pos);
	if (aa != ref_aa) {
		outfp << "WARNING! Amino acid " << aa << " is at position " <<
			to_string (aa_pos + 1) << ", but your list of substitutions assumes it's a " << ref_aa << endl;
	}		

}

string print_double (double num, int precision) 
{
	stringstream stream;
	stream << fixed << std::setprecision (precision) << num;
	return stream.str();

}

void printSubstFile (const std::list <std::string> substList,std::unordered_map <string, double> medianSeqInfoForPos,const  std::vector<std::vector<double>> SIFTscores,const  std::vector <double> aas_stored,const  int total_seq,  Chain* query,const  std::string outfile) {
//	cout.flush() << "in printSubstFile2" << endl;
	std::list<std::string>::const_iterator iterator;
	std::regex regexSubst ("^([A-Z])([0-9]+)([A-Z])");  /*, std::regex_constants::basic); */
        std::smatch m;
	ofstream outfp;
	outfp.open (outfile, ios::out);
	int query_length = SIFTscores.size();

//	cout.flush() << "in printSubstFile " << endl;
	
	for (int pos = 0; pos < query_length; pos++) {
		 char ref_aa = chainGetChar (query, pos);
		 int ref_aa_index = (int) ref_aa - (int) 'A';
	//cout << "in ref del " << to_string (pos) <<  endl;
		 if (SIFTscores[pos][ref_aa_index] < TOLERANCE_PROB_THRESHOLD) {
			auto search = medianSeqInfoForPos.find (to_string(pos));
			double median = search->second; 
			if (median < ADEQUATE_SEQ_INFO) { 
		       outfp << "WARNING! " << ref_aa << to_string (pos+1) << " not allowed!  score: " << print_double( SIFTscores[pos][ref_aa_index],2) << " median: " << print_double (medianSeqInfoForPos[to_string(pos)],2) << " # of sequence: " << to_string ((int) aas_stored[pos]) << endl;
		
                        }
		} /* end if less than TOLERANCE_PROB_THRESHOLD */
	}

	for (iterator = substList.begin(); iterator != substList.end(); ++iterator) 	{
		 string substLine = *iterator;
//		cout << "reading substLine " << substLine << endl;
		 stringstream ss(stringstream::in|stringstream::out);
		ss << substLine;
		string cleanSubst;
		ss >> cleanSubst;	

		if (regex_search (substLine, m, regexSubst)) {
			char tmp[10];
			std::string tmp_string = m[1];
			strncpy (tmp, tmp_string.c_str(), tmp_string.size()+1);
                        char ref_aa = tmp[0]; //.c_str()[0];
			string aa_pos_string = m[2];
			int aa_pos = stoi ( aa_pos_string); //.c_str());
			aa_pos = aa_pos -1; // position in array
			tmp_string = m[3];
			strncpy (tmp, tmp_string.c_str(), tmp_string.size()+1);
			char new_aa = tmp[0]; //.c_str()[0];
			
			int new_aa_index = (int) new_aa - (int) 'A';
			double score = SIFTscores[aa_pos][new_aa_index]; 
	
			check_refaa_against_query (ref_aa, aa_pos, query, outfp);		
			outfp << cleanSubst << "\t"; 
			if (score >= TOLERANCE_PROB_THRESHOLD) {
				outfp << "TOLERATED\t" << print_double (score,2);   
			} else {
				outfp << "DELETERIOUS\t" << print_double (score,2);	
			}
			auto search = medianSeqInfoForPos.find (aa_pos_string);
			outfp << "\t" << print_double( search->second,2)
				 << "\t" <<
			to_string ((int) aas_stored[aa_pos]) << "\t" << to_string (total_seq)  << endl;
 
		
		}
//	std::cout << *iterator << std::endl;
	}	
	outfp.close();
}

bool valid_amino_acid (char aa) 
{
	if (aa == 'B' or aa == 'Z' or aa == 'J' or aa == 'O' or aa == 'U' or aa == 'X' or aa == '-' or aa == '*') {
		return false;
	} else {
		return true;
	}
}

void calcSIFTScores (std::vector <Chain *> alignment_string, Chain* query, std::vector <std::vector <double>>  matrix, std::vector <std::vector <double>> & SIFTscores) {

	int query_length = matrix.size();

	vector <int> max_aa_array (query_length);
	vector <double> number_of_diff_aas (query_length);
	vector <double> epsilon (query_length);

	int num_seqs_in_alignment = alignment_string.size();
        std::vector<std::vector<double>> raw_count_matrix (query_length, 
                std::vector <double> (26, 0.0));

	/* initialize an array with weight 1 */
	std::vector <double> weights_1 (num_seqs_in_alignment, 1.0);

	/* raw counts of valid amino acids */
        std::vector<double> aas_stored_at_each_pos (query_length);
        createMatrix (alignment_string, query, weights_1, raw_count_matrix, aas_stored_at_each_pos);
//        printMatrix (matrix, "tmp.txt");

        //cout << "about to enter calcSeqWeights " << endl;
        std::vector <double> seq_weights (num_seqs_in_alignment);

	/* calcSeqWeights is the only function where all sequences 
	and all positions have to be looked at at the same time. 
	all other routines treat each position independently */ 
       calcSeqWeights (alignment_string, matrix, aas_stored_at_each_pos, seq_weights,
          number_of_diff_aas);
	
	/* now construct matrix with weighted sequence values */
	std::vector<std::vector<double>> seq_weighted_matrix (query_length, 
	                std::vector <double> (26,0.0));
        std::vector <double> tot_weights_each_pos (query_length);

	createMatrix (alignment_string, query, seq_weights, seq_weighted_matrix, tot_weights_each_pos);
//	printMatrix (seq_weighted_matrix, "tmpweighted.txt");

//	cout << "find max aa in matrix " << endl;
	find_max_aa_in_matrix (seq_weighted_matrix, max_aa_array); 

//	cout << "calc Epsilon " << endl;
	calcEpsilon (seq_weighted_matrix, max_aa_array, number_of_diff_aas, epsilon ); 
	/* pseudo_diri */
//	cout << "scale matrix to max aa " << endl;
	std::vector<std::vector<double>> diric_matrix (query_length,
                std::vector <double> (26, 0.0));
	calcDiri (seq_weighted_matrix, diric_matrix);
//	printMatrix (diric_matrix, "diri_matrix.txt");
	
	for (int pos= 0; pos < query_length; pos++) {
		for (char aa = 'A'; aa <= 'Z'; aa++) {
                        int aa_index = (int) aa - (int) 'A';
			SIFTscores[pos][aa_index] = seq_weighted_matrix[pos][aa_index] + epsilon[pos] * diric_matrix[pos][aa_index];
//			cout << "newpos " << to_string (pos) << " aa " << aa << "colcont " << to_string (seq_weighted_matrix[pos][aa_index]) << " epsilon " << to_string (epsilon[pos]) << " colreg " << to_string (diric_matrix[pos][aa_index]) << endl;
			SIFTscores[pos][aa_index] /= (tot_weights_each_pos[pos] + epsilon[pos]);	
		}
//		cout << "tot_weights_each_pos " << to_string (tot_weights_each_pos[pos]) << " epsilon " << epsilon[pos] << endl;
	}	
	/* have to find max aa again because it'schanged with newly added weights */
	find_max_aa_in_matrix (SIFTscores, max_aa_array);
	scale_matrix_to_max_aa (SIFTscores, max_aa_array);
	printMatrix (SIFTscores, "siftscores_matrix.txt");	
//	cout << "done " << endl;
	
}


void calcDiri ( std::vector<std::vector <double>> & count_matrix, std::vector<std::vector<double>> & diric_matrix) {
	int query_length = count_matrix.size();
	for (int pos = 0; pos < query_length; pos++) {
//		cout << "pos for diri " << to_string (pos) << endl;
		add_diric_values (count_matrix[pos], diric_matrix[pos]);
	}
}

double add_logs (double logx, double logy)
{
	if (logx > logy) {
		return (logx + log (1.0 + exp (logy -logx)));
	}
	else {
		return (logy + log (1.0 + exp (logx - logy)));
	}
}

void add_diric_values (std::vector <double> count_col, std::vector <double>& diric_col) {
	double tmp;
	int diri_comp_num = (int) (sizeof (diri_altot)/sizeof (diri_altot[0]));
//	cout << "number of components " << to_string (diri_comp_num) << endl;
	std::vector <double>  probn (diri_comp_num, 0.0);
	std::vector <double> probj (diri_comp_num,0.0);

	double pos_count_tot = 0.0;
	for (int j=0; j < (int) count_col.size(); j++) {
		pos_count_tot += count_col[j];
	}

	/*-----------   compute equation (3), Prob(n|j) ------------  */
	for (int j = 0; j < diri_comp_num; j++) {
		probn[j] = lgamma (pos_count_tot+1.0) + lgamma (diri_altot[j]);
		probn[j] -= lgamma (pos_count_tot + diri_altot[j]);
	
		for (char aa = 'A'; aa <= 'Z'; aa++) {
			int aa_index = (int) aa - (int) 'A';
			if (valid_amino_acid (aa)) {
				tmp = lgamma (count_col[aa_index] + diri_alpha[j][aa_index]);
				tmp -= lgamma (count_col[aa_index] + 1.0);
				tmp -= lgamma (diri_alpha[j][aa_index]);
				probn[j] += tmp;
//				cout << " aa " << aa << "count " << count_col[aa_index] << " tmp " << to_string (tmp) << "comp " << to_string(j) << " probn " << to_string (probn[j]) << endl; 
			} /* end if greater than 0 */
		} /* end all amino acids */
	} /* end for all Diri components */ 	

	/*------ compute sum qk * p(n|k) using logs & exponents ----------*/
	double denom = log (diri_q[0]) + probn[0];

	for (int j =1; j < diri_comp_num; j++) {
		double tmp = log (diri_q[j]) + probn[j];
		denom = add_logs (denom, tmp);
	}
//	cout << "4gdiri denom " << denom << endl;

	/*   compute equation (3), Prob(j|n)  */
	for (int j = 0; j < diri_comp_num; j++) {
		probj[j] = log (diri_q[j]) + probn[j] - denom;
//		cout << "4gdiri j " << to_string (j) << " probn[j] " << to_string (probn[j]) << " probj[j] " << to_string (probj[j]) << endl;
	}

	double totreg =0.0;
	for (char aa = 'A'; aa <= 'Z'; aa++) {

		if (valid_amino_acid (aa)) {
			int aa_index = (int) aa - (int) 'A';
			for (int j =0; j < (int) probj.size(); j++) {
				diric_col[aa_index] += exp (probj[j]) * diri_alpha[j][aa_index];
			} /* gone through all components */
			totreg += diric_col[aa_index];
		} /* end valid amino acid */
	}
//	cout << "4gdiri total prob " << totreg << endl;	
	/* now normalize */
	for (char aa= 'A'; aa <= 'Z'; aa++) {
		int aa_index = (int) aa - (int) 'A';
		diric_col[aa_index] /= totreg;
	}
}


void calcSeqWeights (const std::vector <Chain *> alignment_string, std::vector<std::vector<double>> matrix, std::vector <double> amino_acids_present, std::vector<double> & seq_weights, std::vector <double>& number_of_diff_aas) {
	int query_length = matrix.size();
	/* first calculate number of different aas at each position */
/*	std::vector <int> number_of_diff_aas (query_length); */

	/* initialize */
	for (int pos = 0; pos < query_length; pos++) {
		number_of_diff_aas[pos] = 0;
	}
	for (uint32_t seq_index = 0; seq_index < alignment_string.size();
             seq_index++) {
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
//		cout << " pos " << to_string (pos) << " " << to_string (number_of_diff_aas[pos]) << endl;  
	}


	double tot = 0.0;	
	/* now calculate position-based weights */
	for (uint32_t seq_index = 0; seq_index < alignment_string.size();
	     seq_index++) {
		for (int pos = 0; pos < query_length; pos++) {
			char aa = chainGetChar (alignment_string[seq_index], pos);
			int aa_index = (int) aa - (int) 'A';
                        if (valid_amino_acid (aa) &&  matrix[pos][aa_index] > 0.0) {
				double tmp = number_of_diff_aas[pos] * matrix[pos][aa_index]; 
			/*	cout << "tmp " << to_string (tmp) << " diff " 
			<< to_string (number_of_diff_aas[pos]) << " ," <<
			to_string (matrix[aa_index][pos]) << endl; */
				seq_weights[seq_index] += 1.0/tmp;
			}
		}
//		cout << "seq weight " << to_string (seq_index) << " : " << to_string (seq_weights[seq_index]) << endl;
		tot += seq_weights[seq_index];
	} 

//	 cout << " total weight " << to_string (tot)  << endl;

	/* normalize so weights sum up to the number of sequences */	
	double new_tot_weight = 0.0;
	for (uint32_t seq_index = 0; seq_index<alignment_string.size(); seq_index++)
{
		seq_weights[seq_index]  = seq_weights[seq_index] / tot * alignment_string.size()  ;
//		cout << "seq weight " << to_string (seq_index) << " : " << to_string (seq_weights[seq_index]) << endl;
		new_tot_weight += seq_weights[seq_index];
	}
//	cout << "new total weight " << to_string (new_tot_weight)  << endl;

	/* print sequence weights for checking */
/*	for (uint32_t seq_index = 0; seq_index<alignment_string.size(); seq_index++) {
		cout << "seq index " << std::to_string (seq_index) << " weight " << std::to_string (seq_weights[seq_index]) << endl;
	}
*/
}

void remove_seqs_percent_identical_to_query(Chain *queries, std::vector<Chain *> &alignment_string, double seq_identity)

{

    double identity;

    double seqTotal;

//    cout << "remove_seqs iest entered here\n";

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
                }/*else{

                    cout << "Invalid amino acid found!" << endl;

                  exit(1); 

                }
		*/
            }

        }else{

            cout << "Length does not match!" << endl;

            exit(1);

        }



        double perc_similar = (identity / seqTotal) * 100;



        /*Read curr element, delete if beyond threshold, else move to next pos*/
	/// BUG, (int) round <value>
        if ( perc_similar >= seq_identity ){

            alignment_string.erase(alignment_string.begin() + currPos);

        }else{

            currPos++;

        }

    }

        /*iterates vector and prints out name*/

/*    for (std::vector<Chain *>::const_iterator i = alignment_string.begin(); i != alignment_string.end(); ++i){

        cout << "AFTER_NAME: " << chainGetName(*i) << endl;

    }
*/
}

void printSeqNames (std::vector <Chain *> alignment_string ) {
        for (int i=0; i < int(alignment_string.size()); i++){
		cout << "printName in here " << to_string (i) << endl;
                cout << "printName" << chainGetName (alignment_string[i]) << endl;
        }
}

/* add it with seq_weights 1 or pb weights
with 1 it's raw count, with pb_weights it's weighted */ 
void createMatrix (const std::vector <Chain *> alignment_string, 
			Chain* query, 
		    const std::vector <double> seq_weights, 
		    std::vector<std::vector<double>> & matrix, 
		   std::vector<double> & tot_pos_weight) {
	int query_len = chainGetLength(query);
	for (uint32_t seq_index = 0; 
	     seq_index < alignment_string.size();
             seq_index++) {
		for (int pos = 0; pos < query_len; pos++) {
			char aa = chainGetChar (alignment_string[seq_index], pos);
			if (valid_amino_acid (aa)) {
		/*		std::cout << "pos " << std::to_string(pos) << "sequence " << std::to_string (seq_index) << " aa" << std::to_string(aa) << std::endl; */
				int aa_index = (int) aa - (int) 'A'; 
				matrix[pos][aa_index] += seq_weights[seq_index] ; /*seq_weight[seq_index];*/ 
				tot_pos_weight[pos] += seq_weights[seq_index];
			}
		}
	}	
}

void printMatrix (std::vector<std::vector<double>> & matrix, std::string filename) {
	int query_length = matrix.size();
	int aas = matrix[0].size();
	std::ofstream out_file;

	out_file.open (filename);

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

int aa_to_idx (char character){

	// A = 0 and etc

	return int(character) - int('A');

}



void basic_matrix_construction (const std::vector <Chain*> alignment_string, const  std::vector <double> seq_weights, std::vector<std::vector<double>>& matrix){

	// do partial calculation on aa that can be represented by B or Z as well
	double part_D = aa_frequency[aa_to_idx('D')] / ( aa_frequency[aa_to_idx('D')] + aa_frequency[aa_to_idx('N')] );
	double part_N = aa_frequency[aa_to_idx('N')] / ( aa_frequency[aa_to_idx('D')] + aa_frequency[aa_to_idx('N')] );
	double part_E = aa_frequency[aa_to_idx('E')] / ( aa_frequency[aa_to_idx('E')] + aa_frequency[aa_to_idx('Q')] );
	double part_Q = aa_frequency[aa_to_idx('Q')] / ( aa_frequency[aa_to_idx('E')] + aa_frequency[aa_to_idx('Q')] );


	int proteinLen = chainGetLength(alignment_string[0]);	
	int aaLen =  matrix[0].size();
	/*int aaLen = 28; // MATRIX_AA_WIDTH in old code
	matrix.resize(aaLen);
*/
	// loop every position in chain length
	for (int pos = 0; pos < proteinLen; pos++){
		double total = 0.0;

                // expand by 1 and zero out col j
/*                for (int i = 0; i < aaLen; i++){
	                matrix[i].push_back(0.0);
                }
*/
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
			}else if (currChar == 'Z') {
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
			}else{
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
			if (n <= aa_to_idx('Z') || n != aa_to_idx('X')){
                        	matrix[pos][n] = matrix[pos][n] * 100.0 / total;
			// X and * and -
			}else{
                                matrix[pos][n] = aa_frequency[n];
			}
	        }

		// lastly, tali B and Z rows

		matrix[pos][aa_to_idx('B')] = (matrix[pos][aa_to_idx('D')] * part_D) + (matrix[pos][aa_to_idx('N')] * part_N);

                matrix[pos][aa_to_idx('Z')] = (matrix[pos][aa_to_idx('E')] * part_E) + (matrix[pos][aa_to_idx('Q')] * part_Q);

	}
}

void seqs_without_X  (const std:: vector <Chain*>  alignment_string, int pos, std::vector <Chain*>& out_alignment){


        char posChar;

        /*iterate vector of aligment and check pos*/

        for( int n = 0; n < int(alignment_string.size()); n++){
                /*cout << "NAME: " << chainGetName(alignment_string[n]) << endl;
                for (int m=0; m < chainGetLength(alignment_string[n]); m++){
                        cout << chainGetChar(alignment_string[n], m);
                }
                cout << endl;*/

                posChar = chainGetChar(alignment_string[n], pos);

                /*cout << posChar << " :: ";*/

                /*check if char in pos is valid, and push it into new vector*/
                if (valid_amino_acid(posChar)){
                        out_alignment.push_back(alignment_string[n]);
                        /*cout << "YAS" << endl;*/
                }/*else{

                        cout << "NO" << endl;

                }*/
        } 

}
