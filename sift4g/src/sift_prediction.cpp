/*!
 * @file sift_prediction.cpp
 *
 * @brief SIFT predictions source file
 */

#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "utils.hpp"
#include "sift_prediction.hpp"

// BLIMPS needs that...
#define EXTERN

extern "C" {
    #include "sift/info_on_seqs_pid.h"
    #include "sift/Alignment.h"
}

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

        int32_t num_sequences = 400 < alignment_strings[i].size() ? 400 : alignment_strings[i].size();

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
