/*!
 * @file select_alignments.cpp
 *
 * @brief Alignment selection source file
 */

#include <math.h>
#include <assert.h>
#include <cstring>
#include <fstream>
#include <algorithm>

#include "utils.hpp"
#include "select_alignments.hpp"

constexpr double kLog_2_20 = 4.321928095;

class ThreadData {
public:
    ThreadData(std::vector<Chain*>& _dst, DbAlignment** _alignments, int _alignments_length,
        Chain* _query, float _threshold):
            dst(_dst), alignments(_alignments), alignments_length(_alignments_length),
            query(_query), threshold(_threshold) {
    }

    std::vector<Chain*>& dst;
    DbAlignment** alignments;
    int alignments_length;
    Chain* query;
    float threshold;
};
float getMedian(float* a, int len);

void aligmentStr(char** query_str, char** target_str, Alignment* alignment, const char gap_item);

void alignmentsExtract(std::vector<Chain*>& dst, Chain* query, DbAlignment** alignments,
    int alignments_length);

int alignmentsSelect(std::vector<Chain*>& alignment_strings, Chain* query, float threshold);

void* threadSelectAlignments(void* params);

/*****************************************************************************
*****************************************************************************/

void selectAlignments(std::vector<std::vector<Chain*>>& dst, DbAlignment*** alignments,
    int32_t* alignments_lengths, Chain** queries, int32_t queries_length,
    float threshold) {

    dst.resize(queries_length);

    fprintf(stderr, "** Selecting alignments with median threshold: %.2f **\n", threshold);

    std::vector<ThreadPoolTask*> thread_tasks(queries_length, nullptr);

    for (int32_t i = 0; i < queries_length; ++i) {

        auto thread_data = new ThreadData(dst[i], alignments[i], alignments_lengths[i],
            queries[i], threshold);

        thread_tasks[i] = threadPoolSubmit(threadSelectAlignments, (void*) thread_data);
    }

    for (int32_t i = 0; i < queries_length; ++i) {
        threadPoolTaskWait(thread_tasks[i]);
        threadPoolTaskDelete(thread_tasks[i]);
        query_log(i + 1, queries_length);
    }

    fprintf(stderr, "\n\n");
}

void outputSelectedAlignments(std::vector<std::vector<Chain*>>& alignment_strings,
    Chain** queries, int32_t queries_length, const std::string& out_path) {

    std::string out_extension = ".aligned.fasta";

    for (uint32_t i = 0; i < alignment_strings.size(); ++i) {

        std::ofstream out_file;
        char* out_file_name = createFileName(chainGetName(queries[i]), out_path, out_extension);

        int query_len = chainGetLength(queries[i]);

        out_file.open(out_file_name);
        out_file << ">QUERY" << std::endl;

        for (int j = 1; j < query_len + 1; ++j) {
            out_file << chainGetChar(queries[i], j - 1);
            if (j % 60 == 0) out_file << std::endl;
        }
        out_file << std::endl;

        for (uint32_t j = 0; j < alignment_strings[i].size(); ++j) {
            out_file << ">" << chainGetName(alignment_strings[i][j]) << std::endl;

            for (int k = 1; k < query_len + 1; ++k) {
                out_file << chainGetChar(alignment_strings[i][j], k - 1);
                if (k % 60 == 0) out_file << std::endl;
            }
            out_file << std::endl;
        }

        out_file.close();
        delete[] out_file_name;
    }
}

void deleteSelectedAlignments(std::vector<std::vector<Chain*>>& alignment_strings) {
    for (uint32_t i = 0; i < alignment_strings.size(); ++i) {
        for (uint32_t j = 0; j < alignment_strings[i].size(); ++j) {
            if (alignment_strings[i][j] != nullptr) {
                chainDelete(alignment_strings[i][j]);
                alignment_strings[i][j] = nullptr;
            }
        }
        std::vector<Chain*>().swap(alignment_strings[i]);
    }
}

/*****************************************************************************
*****************************************************************************/

float getMedian(float* a, int len) {

	std::sort(&a[0], &a[len - 1]);

	if (len % 2 == 0) {
        return (a[len / 2 - 1] + a[len / 2]) / 2.0;
    } else {
        return a[len / 2];
    }
}

void alignmentsExtract(std::vector<Chain*>& dst, Chain* query, DbAlignment** alignments,
    int alignments_length) {

	int query_len = chainGetLength(query);

	char* alignment_str = new char[query_len];
	int query_start, length;

	Chain* target = nullptr;

	char* query_str = nullptr;
	char* target_str = nullptr;
	Alignment* alignment = nullptr;

	const char gap_item = '-';

	for (int i = 0; i < alignments_length; ++i) {

		target = dbAlignmentGetTarget(alignments[i]);

		query_start = dbAlignmentGetQueryStart(alignments[i]);
		length = dbAlignmentGetPathLen(alignments[i]);

		alignment = dbAlignmentToAlignment(alignments[i]);
		aligmentStr(&query_str, &target_str, alignment, gap_item);
		alignmentDelete(alignment);

		int j = 0;
		for (; j < query_start; ++j) {
			alignment_str[j] = 'X';
		}

		for (int k = 0; k < length; ++k) {
			if (query_str[k] != gap_item) {
				if (target_str[k] != gap_item) {
					alignment_str[j++] = target_str[k];
				} else {
					alignment_str[j++] = 'X';
				}
			}
		}

		for (; j < query_len; ++j) {
			alignment_str[j] = 'X';
		}

        dst.push_back(chainCreate((char*) chainGetName(target), strlen(chainGetName(target)),
            alignment_str, query_len));

        delete[] query_str;
        delete[] target_str;
	}

    delete[] alignment_str;
}

int alignmentsSelect(std::vector<Chain*>& alignment_strings, Chain* query, float threshold) {

	int amino_acid_num = 26;
	float median = kLog_2_20;

    int* amino_acid_nums = new int[amino_acid_num];
	for (int i = 0; i < amino_acid_num; ++i) {
		amino_acid_nums[i] = 0;
	}

	int query_len = chainGetLength(query);

    float* pos_freq = new float[query_len];
	for (int i = 0; i < query_len; ++i) {
		pos_freq[i] = 0.0;
	}

	char c;
	int i, valid;
	for (i = 1; median > threshold && i <= (int) alignment_strings.size(); ++i) {

		for (int j = 0; j < query_len; ++j) {
			valid = 0;

			for (int k = 0; k < i; ++k) {
				c = chainGetChar(alignment_strings[k], j);
				if (c != 'X') {
					valid++;
					amino_acid_nums[(int) c - 'A']++;
				}
			}

			for (int k = 0; k < amino_acid_num; ++k) {
				if (amino_acid_nums[k] != 0) {
					pos_freq[j] += amino_acid_nums[k] / (float) valid *
						log2f(amino_acid_nums[k] / (float) valid);
				}
			}

			pos_freq[j] += kLog_2_20;

			for (int k = 0; k < amino_acid_num; ++k) {
				if (amino_acid_nums[k] != 0) {
					amino_acid_nums[k] = 0;
				}
			}
		}

		median = getMedian(pos_freq, query_len);

		for (int j = 0; j < query_len; ++j) {
			pos_freq[j] = 0.0;
		}
	}

    delete[] pos_freq;
    delete[] amino_acid_nums;

	return i - 1;
}

void aligmentStr(char** query_str, char** target_str, Alignment* alignment, const char gap_item) {
    /* code take from SW# (its a private function)*/
    Chain* query = alignmentGetQuery(alignment);
    int query_start = alignmentGetQueryStart(alignment);

    Chain* target = alignmentGetTarget(alignment);
    int target_start = alignmentGetTargetStart(alignment);

    int path_len = alignmentGetPathLen(alignment);

    int query_idx = query_start;
    int target_idx = target_start;

    *query_str = new char[path_len];
    *target_str = new char[path_len];

    int i;
    for (i = 0; i < path_len; ++i) {

        char query_chr;
        char target_chr;

        switch (alignmentGetMove(alignment, i)) {
        case MOVE_LEFT:

            query_chr = gap_item;
            target_chr = chainGetChar(target, target_idx);

            target_idx++;

            break;
        case MOVE_UP:

            query_chr = chainGetChar(query, query_idx);
            target_chr = gap_item;

            query_idx++;

            break;
        case MOVE_DIAG:

            query_chr = chainGetChar(query, query_idx);
            target_chr = chainGetChar(target, target_idx);

            query_idx++;
            target_idx++;

            break;
        default:
            // error
            return;
        }

        (*query_str)[i] = query_chr;
        (*target_str)[i] = target_chr;
    }
}

void* threadSelectAlignments(void* params) {

    auto thread_data = (ThreadData*) params;

    thread_data->dst.clear();

    alignmentsExtract(thread_data->dst, thread_data->query, thread_data->alignments,
        thread_data->alignments_length);

    uint32_t selected_alignments_length = alignmentsSelect(thread_data->dst,
        thread_data->query, thread_data->threshold);

    for (uint32_t i = selected_alignments_length; i < thread_data->dst.size(); ++i) {
        chainDelete(thread_data->dst[i]);
    }
    thread_data->dst.resize(selected_alignments_length);

    delete thread_data;

    return nullptr;
}
