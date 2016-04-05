/*!
 * @file database_alignment.cpp
 *
 * @brief Database alignment source file
 */

#include <assert.h>

#include "database_alignment.hpp"

constexpr uint32_t database_chunk = 1000000000; /* ~1GB */
constexpr float log_step_percentage = 2.5;

void database_alignment_log(uint32_t part, float part_size, float percentage) {
    fprintf(stderr, "* aligning database part %u (size ~%.2f GB): %.2f/100.00%% *\r",
        part, part_size, percentage);
    fflush(stderr);
}

void valueFunction(double* values, int* scores, Chain* query, Chain** database,
    int databaseLen, int* cards, int cardsLen, void* param_ );

void createFilteredDatabase(std::vector<uint32_t>& used_indices, Chain*** filtered_database,
    std::vector<uint32_t>& indices, Chain** database, uint32_t database_length);

void alignDatabase(DbAlignment**** alignments, int** alignments_lengths, Chain*** _database,
    int32_t* _database_length, const std::string& database_path, Chain** queries,
    int32_t queries_length, std::vector<std::vector<uint32_t>>& indices,
    int32_t algorithm, EValueParams* evalue_params, double max_evalue,
    uint32_t max_alignments, Scorer* scorer, int32_t* cards,
    int32_t cards_length, const std::string& out_path,
    int32_t out_format) {

    fprintf(stderr, "** Aligning queries with candidate sequences **\n");

    Chain** database = nullptr;
    int database_length = 0;
    int database_start = 0;

    FILE* handle = nullptr;
    int serialized = 0;
    readFastaChainsPartInit(&database, &database_length, &handle, &serialized,
        database_path.c_str());

    uint32_t part = 1;
    float part_size = database_chunk / (float) 1000000000;
    uint32_t log_size = queries_length / (100. / log_step_percentage);

    while (true) {

        int status = 1;

        status &= readFastaChainsPart(&database, &database_length, handle,
            serialized, database_chunk);

        uint32_t log_counter = 0;
        float log_percentage = log_step_percentage;

        database_alignment_log(part, part_size, 0);

        /* have to malloc... */
        DbAlignment*** alignments_part = (DbAlignment***) malloc(queries_length * sizeof(DbAlignment**));
        int* alignments_part_lengths = (int*) malloc(queries_length * sizeof(int));

        std::vector<bool> used_sequences(database_length, false);

        for (int32_t i = 0; i < queries_length; ++i) {

            ++log_counter;
            if (log_counter % log_size == 0 && log_percentage < 100.) {
                database_alignment_log(part, part_size, log_percentage);
                log_percentage += log_step_percentage;
            }

            std::vector<uint32_t> used_indices;
            Chain** filtered_database = nullptr;
            createFilteredDatabase(used_indices, &filtered_database, indices[i],
                database, database_length);

            if (used_indices.empty()) {
                alignments_part[i] = nullptr;
                alignments_part_lengths[i] = 0;
                continue;
            }

            ChainDatabase* chain_database = chainDatabaseCreate(filtered_database, 0,
                used_indices.size(), cards, cards_length);

            alignDatabase(&alignments_part[i], &alignments_part_lengths[i], algorithm,
                queries[i], chain_database, scorer, max_alignments, valueFunction,
                (void*) evalue_params, max_evalue, nullptr, 0, cards,
                cards_length, nullptr);

            for (int32_t j = 0; j < alignments_part_lengths[i]; ++j) {
                used_sequences[used_indices[dbAlignmentGetTargetIdx(alignments_part[i][j])]] = true;
            }

            chainDatabaseDelete(chain_database);

            delete[] filtered_database;
        }

        if (*alignments == nullptr) {
            *alignments = alignments_part;
            *alignments_lengths = alignments_part_lengths;
        } else {
            dbAlignmentsMerge(*alignments, *alignments_lengths, alignments_part,
                alignments_part_lengths, queries_length, max_alignments);
            deleteShotgunDatabase(alignments_part, alignments_part_lengths, queries_length);
        }

        for (int32_t i = database_start; i < database_length; ++i) {
            if (!used_sequences[i] && database[i] != nullptr) {
                chainDelete(database[i]);
                database[i] = nullptr;
            }
        }

        database_alignment_log(part, part_size, 100);
        fprintf(stderr, "\n");
        ++part;

        if (status == 0) {
            break;
        }

        database_start = database_length;
    }

    fclose(handle);
    *_database = database;
    *_database_length = database_length;
}

void valueFunction(double* values, int* scores, Chain* query, Chain** database,
    int databaseLen, int* cards, int cardsLen, void* param_ ) {

    EValueParams* eValueParams = (EValueParams*) param_;
    eValues(values, scores, query, database, databaseLen, cards, cardsLen, eValueParams);
}

void createFilteredDatabase(std::vector<uint32_t>& used_indices, Chain*** filtered_database,
    std::vector<uint32_t>& indices, Chain** database, uint32_t database_length) {

    uint32_t database_end = database_length == 0 ? 0 : database_length - 1;

    uint32_t i = 0;
    for (; i < indices.size(); ++i) {
        if (indices[i] > database_end) {
            break;
        }
    }

    uint32_t used_indices_length = i;

    if (used_indices_length != 0) {
        used_indices.reserve(used_indices_length);
        *filtered_database = new Chain*[used_indices_length]();

        for (uint32_t j = 0; j < used_indices_length; ++j) {
            used_indices.emplace_back(indices[j]);
            (*filtered_database)[j] = database[indices[j]];
        }

        std::vector<uint32_t> tmp(indices.begin() + used_indices_length, indices.end());
        indices.swap(tmp);
    }
}
