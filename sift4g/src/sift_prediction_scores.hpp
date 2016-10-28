/*!
 * @file sift_prediction.hpp
 *
 * @brief SIFT predictions header file
 */

#pragma once

#include <stdint.h>
#include <vector>
#include <string>
#include <list>
#include <unordered_map>
#include "swsharp/swsharp.h"
using namespace std;

void readSubstFile (char* substFilename, std::list <std::string> &substList);

//void printSubstFile (std::list <std::string> substList);
void printSubstFile (const std::list <std::string> substList,std::unordered_map <string, double> medianSeqInfoForPos,const std::vector<std::vector<double>> SIFTscores,const std::vector <double> aas_stored,const int total_seq,Chain* query,const std::string outfile);
 
void createMatrix (const std::vector <Chain *> alignment_string,Chain* query,const  std::vector <double> seq_weights, std::vector<std::vector<double>> & matrix, std::vector<double> & tot_pos_weight) ; 

void remove_seqs_percent_identical_to_query(Chain *queries, std::vector <Chain*>& alignment_string, double seq_identity);

/* this has to be after createMatrix because it uses amino_acids_present */
void calcSeqWeights (const std::vector <Chain *> alignment_string, std::vector<std::vector<double>> matrix, std::vector <double> amino_acids_present, std::vector<double> & seq_weights, std::vector <double>& number_of_diff_aas);

void seqs_without_X  (const std:: vector <Chain*>  alignment_string, int pos, std::vector <Chain*>& out_alignment);

void basic_matrix_construction (const std::vector <Chain*> alignment_string,const  std::vector <double> seq_weights, std::vector<std::vector<double>>& matrix);

int aa_to_idx (char character);

bool valid_amino_acid (char aa);

double  calculateMedianSeqInfo (std::vector<Chain*>& alignment_strings, std::unordered_map <int, int> invalidSeq, std::vector <double> seq_weights, std::vector<std::vector<double>> matrix);

void hashPredictedPos (std::list <std::string> &substList, std::unordered_map <std::string, double>& medianSeqInfoForPos);

void addPosWithDelRef (Chain* query,  std::vector<std::vector<double>> SIFTscores, std::unordered_map <string, double>& medianSeqInfoForPos);

void addMedianSeqInfo (std::vector<Chain*> alignment_string, Chain* query, std::vector<std::vector<double>> matrix, std::unordered_map <std::string, double>& medianSeqInfoForPos);

void printSeqNames (std::vector <Chain *> alignment_string ); 

void check_refaa_against_query (char ref_aa, int aa_pos, Chain* query, std::ofstream& outfp) ;
void printMatrix (std::vector<std::vector<double>> &  matrix, std::string filename);

void calcDiri (std::vector<std::vector<double>> &  weighted_matrix, std::vector<std::vector<double>> &  diri_matrix);

void calcSIFTScores (std::vector <Chain *> alignment_string, Chain* query, std::vector <std::vector <double>>  matrix, std::vector <std::vector <double>>& SIFTscores);

void add_diric_values (std::vector <double> count_col, std::vector <double>& diric_col);

double add_logs (double logx, double logy);

/* generate from default.rank, but re-orderd to alphabet index using
convert_rank_matrix.py . Matrix is assymetric. rank_matrix[row][col]*/
constexpr int rank_matrix[26][26] = 
   {
	{1,-1,6,17,10,19,3,16,12,-1,7,14,11,15,-1,9,8,13,2,4,-1,5,20,-1,18,-1},
{-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
{2,-1,1,19,20,10,12,16,6,-1,17,7,8,13,-1,14,15,18,5,4,-1,3,9,-1,11,-1},
{12,-1,17,1,2,18,9,8,15,-1,6,19,13,3,-1,10,5,11,4,7,-1,16,20,-1,14,-1},
{10,-1,20,3,1,18,14,6,19,-1,4,17,12,8,-1,11,2,5,7,9,-1,15,16,-1,13,-1},
{10,-1,12,19,18,1,16,8,6,-1,15,4,5,14,-1,20,17,13,11,9,-1,7,3,-1,2,-1},
{2,-1,14,5,10,17,1,9,20,-1,6,19,15,4,-1,11,8,12,3,7,-1,18,13,-1,16,-1},
{12,-1,18,9,5,10,14,1,20,-1,7,17,11,3,-1,15,4,6,8,13,-1,19,16,-1,2,-1},
{8,-1,7,16,17,5,20,19,1,-1,12,3,4,18,-1,13,14,15,10,6,-1,2,11,-1,9,-1},
{-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
{10,-1,19,8,4,20,13,9,17,-1,1,16,12,5,-1,11,3,2,6,7,-1,15,18,-1,14,-1},
{9,-1,8,19,16,5,20,15,3,-1,14,1,2,18,-1,17,11,12,13,7,-1,4,10,-1,6,-1},
{8,-1,12,20,16,5,19,15,3,-1,10,2,1,17,-1,18,6,11,14,7,-1,4,13,-1,9,-1},
{11,-1,15,2,8,17,9,4,18,-1,7,19,14,1,-1,12,5,10,3,6,-1,16,20,-1,13,-1},
{-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
{3,-1,16,8,6,19,11,12,15,-1,4,17,14,9,-1,1,7,10,2,5,-1,13,20,-1,18,-1},
{11,-1,19,8,2,20,14,5,18,-1,3,16,9,6,-1,12,1,4,7,10,-1,17,15,-1,13,-1},
{10,-1,20,11,4,18,15,5,19,-1,2,14,9,6,-1,13,3,1,7,8,-1,16,17,-1,12,-1},
{3,-1,12,8,6,18,9,13,17,-1,7,19,14,4,-1,11,5,10,1,2,-1,15,20,-1,16,-1},
{3,-1,11,12,10,19,16,18,9,-1,7,15,6,4,-1,13,8,14,2,1,-1,5,20,-1,17,-1},
{-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
{6,-1,7,20,14,8,19,18,2,-1,12,3,4,17,-1,13,11,15,10,5,-1,1,16,-1,9,-1},
{11,-1,7,20,16,3,10,8,12,-1,17,5,4,19,-1,18,6,13,14,9,-1,15,1,-1,2,-1},
{-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
{13,-1,17,20,15,2,19,4,8,-1,14,6,5,16,-1,18,9,12,11,10,-1,7,3,-1,1,-1},
{-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1}
   }; 


/*
    13-component Dirichlet
    Name = merge-opt.13comp
    Order of amino acids is A,B,C,D,...    
*/

const double diri_q[13] = { 0.10935200, 0.08624370, 0.09149870, 0.08177380, 0.07078570, 0.04727980, 0.06654990, 0.10472800, 0.09536710, 0.05250920, 0.03843090, 0.09601160, 0.05946920}; 

const double aa_frequency[26] = {0.07587,0,0.01666,0.05287,0.06377,0.04096,0.06847,0.02246,0.05817,0,0.05957,0.09427,0.02376,0.04456,0,0.04906,0.03976,0.05166,0.07127,0.05677,0,0.06587,0.01236,0,0.03186,0};

const double diri_altot[13] = { 1.30226000, 1.81388000, 1.70395000, 2.11182000, 4.96847000, 3.33809000, 2.19222000, 0.05904830, 0.08441720, 2.10287000, 1.10699000, 7.12351000, 2.18943000}; 

const double diri_alpha[13][26] = 
	{
	{0.388746,0,0.0367385,0.026518,0.0322076,0.0137612,0.206244,0.0059514,0.0158113,
0,0.0351076,0.0350778,0.0164858,0.0224268,0,0.0917171,0.0241792,0.0253054,0.187681,
0.0626166,0,0.0602319,0.00498456,0,0.0104666,0},
	{0.0517631,0,0.0133196,0.00631385,0.0132966,0.141159,0.0103533,0.0104886,0.264004
,0,0.0177537,0.807661,0.175596,0.0127646,0,0.0220077,0.0271782,0.0177749,0.0144998,
0.0299289,0,0.143995,0.016497,0,0.0175224,0},
	{0.0610492,0,0.0148718,0.315694,0.0908927,0.012609,0.222958,0.0672017,0.006044,0
,0.0877828,0.0135798,0.00485097,0.368744,0,0.0498296,0.0588252,0.044453,0.168811
,0.0652581,0,0.010088,0.00652446,0,0.0338823,0},
	{0.0796662,0,0.00761021,0.0136208,0.110994,0.0165387,0.0563167,0.0713015,0.0341013
,0,0.564383,0.0808963,0.0267238,0.077247,0,0.0627218,0.18985,0.475607,0.0651016
,0.0855045,0,0.0452026,0.0125425,0,0.0358879,0},
	{0.317696,0,0.110939,0.128582,0.170879,0.441344,0.206556,0.339981,0.246824,0,
0.181928,0.407255,0.138334,0.192446,0,0.127958,0.183911,0.217794,0.251084,0.245119,
0,0.314739,0.115968,0,0.629135,0},
	{0.0446044,0,0.0170118,0.0134169,0.014928,0.0818952,0.0279526,2.74064e-06,1.24306,
0,0.0235435,0.433818,0.0752952,0.0156684,0,0.00769209,0.00598883,0.0194969,
0.0237729,0.0605068,0,1.17986,0.00984592,0,0.039723,0},
	{0.142671,0,0.00172471,0.391867,0.796428,0.0127552,0.0507474,0.0337799,0.0169946
,0,0.143701,0.0315845,0.0130077,0.0560163,0,0.086296,0.189149,0.0377603,0.0701466
,0.0617127,0,0.0327986,0.00600271,0,0.0170779,0},
	{0.00281745,0,0.00552169,0.00417499,0.00218621,0.00287303,0.00978086,0.00226513, 
0.00229518,0,0.000288688,0.00312013,0.00106406,0.00244641,0,0.00279557,0.00165027
,0.00285057,0.00125284,0.00271746,0,0.00320849,0.00225774,0,0.00348149,0},
	{0.00322581,0,0.00013552,0.00597818,0.00438741,0.00193385,0.0158901,0.00243118,
0.00206323,0,0.00656549,0.00518845,0.000562418,0.00273016,0,0.0121928,0.00322314,
0.00579248,0.00329082,0.00288987,0,0.00206985,0.00227726,0,0.00158929,0},
	{0.14246,0,0.0306225,0.0639103,0.0396449,0.0154444,0.0394262,0.018336,0.0452637,
0,0.0606921,0.0412351,0.0225931,0.117688,0,0.0752369,0.0475692,0.038186,0.50948,
0.696177,0,0.0708571,0.0102097,0,0.0178433,0},
	{0.0177105,0,0.00908289,0.00274088,0.00444563,0.411122,0.00282536,1.4813e-06,
0.0384249,0,1.47392e-06,0.0931481,0.0186396,1.47581e-06,0,0.00676848,1.47629e-06,
1.47755e-06,0.012354,0.0111622,0,0.035578,0.0874503,0,0.35553,0},
	{0.668982,0,0.0472318,0.554692,0.890999,0.0440375,0.336612,0.107352,0.16552,0,
0.822984,0.267021,0.0936297,0.437926,0,0.242805,0.521498,0.464298,0.649512,0.507446
,0,0.252405,0.0151661,0,0.0333911,0},
	{0.304645,0,0.0610182,0.0139475,0.0452959,0.0312398,1.47952e-06,0.0113568,0.3088088
,0,0.0448228,0.209537,0.0639208,0.020799,0,0.0641272,0.0336194,0.0248684,
0.0434541,0.21585,0,0.677503,0.0052348,0,0.00937749,0}
	}; /* end 2-D array */




