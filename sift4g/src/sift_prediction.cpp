/*!
 * @file sift_prediction.cpp
 *
 * @brief SIFT predictions source file
 *
 * @author: rvaser, pauline-ng, angsm
 */

#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "utils.hpp"
#include "sift_scores.hpp"
#include "sift_prediction.hpp"

constexpr uint32_t kMaxSequences = 400;

class ThreadPredictionData {
public:
    ThreadPredictionData(std::vector<Chain*>& _alignment_strings, Chain* _query, const std::string& _subst_path,
        int32_t _sequence_identity, const std::string& _out_path)
            : alignment_strings(_alignment_strings), query(_query), sequence_identity(_sequence_identity),
            subst_path(_subst_path), out_path(_out_path) {
    }

    std::vector<Chain*>& alignment_strings;
    Chain* query;
    int32_t sequence_identity;
    std::string subst_path;
    std::string out_path;
};

void* threadSiftPredictions(void* params);

/*****************************************************************************
*****************************************************************************/

void siftPredictions(std::vector<std::vector<Chain*>>& alignment_strings,
    Chain** queries, int32_t queries_length, const std::string& subst_path,
    int32_t sequence_identity, const std::string& out_path) {

    fprintf(stderr, "** Generating SIFT predictions with sequence identity: %.2f%% **\n", (float) sequence_identity);

    std::vector<ThreadPoolTask*> thread_tasks(queries_length, nullptr);

    for (int32_t i = 0; i < queries_length; ++i) {

        auto thread_data = new ThreadPredictionData(alignment_strings[i], queries[i],
            subst_path, sequence_identity, out_path);

        thread_tasks[i] = threadPoolSubmit(threadSiftPredictions, (void*) thread_data);
    }

    for (int32_t i = 0; i < queries_length; ++i) {
        threadPoolTaskWait(thread_tasks[i]);
        threadPoolTaskDelete(thread_tasks[i]);
        queryLog(i + 1, queries_length);
    }

    fprintf(stderr, "\n\n");
}

/*****************************************************************************
*****************************************************************************/

void* threadSiftPredictions(void* params) {

    auto thread_data = (ThreadPredictionData*) params;

	std::string subst_extension = ".subst";
	std::string out_extension = ".SIFTprediction";

	// only keep first 399 hits, erase more distant ones
    if (thread_data->alignment_strings.size() > kMaxSequences - 1) {
        for (uint32_t j = kMaxSequences - 1; j < thread_data->alignment_strings.size(); ++j) {
            chainDelete(thread_data->alignment_strings[j]);
        }
        thread_data->alignment_strings.resize(kMaxSequences - 1);
    }

    int query_length = chainGetLength(thread_data->query);
    remove_seqs_percent_identical_to_query(thread_data->query,
        thread_data->alignment_strings, thread_data->sequence_identity);

    // add query sequence to the beginning of the alignment
    thread_data->alignment_strings.insert(thread_data->alignment_strings.begin(),
        chainCreateView(thread_data->query, 0, query_length - 1, 0));
    int total_seq = thread_data->alignment_strings.size();

    std::vector<std::vector<double>> matrix(query_length, std::vector<double>(26, 0.0));
    std::vector<std::vector<double>> SIFTscores(query_length, std::vector<double>(26, 0.0));

    std::vector<double> weights_1(thread_data->alignment_strings.size(), 1.0);
    std::vector<double> aas_stored(query_length, 0.0);

    createMatrix(thread_data->alignment_strings, thread_data->query, weights_1, matrix, aas_stored);

    calcSIFTScores(thread_data->alignment_strings, thread_data->query, matrix, SIFTscores);

    int num_seqs_in_alignment = thread_data->alignment_strings.size();
    std::vector <double> seq_weights (num_seqs_in_alignment);
    std::vector <double> number_of_diff_aas (query_length);
    calcSeqWeights(thread_data->alignment_strings, matrix, aas_stored, seq_weights, number_of_diff_aas);

    char* subst_file_name = createFileName(chainGetName(thread_data->query),
        thread_data->subst_path, subst_extension);
    char* out_file_name = createFileName(chainGetName(thread_data->query),
        thread_data->out_path, out_extension);

    if (isExtantPath(subst_file_name)) {
        std::list<std::string> subst_list;
        std::unordered_map<std::string, double> medianSeqInfoForPos;

        readSubstFile(subst_file_name, subst_list);
        if (checkSubsts(subst_list, thread_data->query) == true) {
            hashPredictedPos(subst_list, medianSeqInfoForPos);
            addPosWithDelRef(thread_data->query, SIFTscores, medianSeqInfoForPos);
            addMedianSeqInfo(thread_data->alignment_strings, thread_data->query,
                matrix, medianSeqInfoForPos);
            printSubstFile(subst_list, medianSeqInfoForPos, SIFTscores, aas_stored,
                total_seq, thread_data->query, out_file_name);
        }
    } else {
        // printMatrix(SIFTscores, out_file_name);
        printMatrixOriginalFormat(SIFTscores, out_file_name);
    }

    delete[] out_file_name;
    delete[] subst_file_name;

    delete thread_data;

    return nullptr;
}
