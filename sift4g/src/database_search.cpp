/*!
 * @file database_search.cpp
 *
 * @brief Database search source file
 */

#include <algorithm>
#include <limits>

#include "hash.hpp"
#include "database_search.hpp"

#include "swsharp/swsharp.h"

class Candidate {

    Candidate(int _score, int _id) :
            score(_score), id(_id) {
    }

    int32_t score;
    int32_t id;
};

class ThreadData {

    ThreadData(Hash* _query_hash, Chain** _database, uint32_t _database_begin,
        uint32_t _database_end, std::vector<std::vector<Candidate>>* _candidates,
        std::vector<Mutex>* _candidates_mutexes, uint32_t _max_candidates):
            query_hash(_query_hash), database(_database), database_begin(_database_begin),
            database_end(_database_end), candidates(_candidates),
            candidates_mutexes(_candidates_mutexes),
            max_candidates(_max_candidates) {
    }

    Hash* query_hash;
    Chain** database;
    uint32_t database_begin;
    uint32_t database_end;
    std::vector<std::vector<Candidate>>* candidates;
    std::vector<Mutex>* candidates_mutexes;
    uint32_t max_candidates;
};

void* threadSearchDatabase(void* params);

int32_t longestIncreasingSubsequence(const std::vector<int32_t>& src);

uint64_t searchDatabase(std::vector<std::vector<uint32_t>>& dst,
    const std::string& database_path, uint32_t database_chunk,
    const std::string& query_path, uint32_t kmer_length,
    uint32_t max_candidates, uint32_t num_threads) {

    Chain** queries = nullptr;
    int queries_length = 0;
    readFastaChains(&queries, &queries_length, query_path.c_str());

    auto hash = createHash(queries, queries_length, 0, queries_length, kmer_length);

    Chain** database = nullptr;
    int database_length = 0;
    int database_start = 0;

    FILE* handle = nullptr;
    int serialized = 0;
    readFastaChainsPartInit(&database, &database_length, &handle, &serialized,
        database_path.c_str());

    uint64_t database_cells = 0;

    std::vector<std::vector<Candidate>> candidates(queries_length);
    std::vector<Mutex> candidates_mutexes(queries_length);
    for (int i = 0; i < queries_length; ++i) {
        mutexCreate(&candidates_mutexes[i]);
    }

    while (true) {

        int status = 1;

        status &= readFastaChainsPart(&database, &database_length, handle,
            serialized, database_chunk);

        std::vector<ThreadPoolTask*> thread_tasks(num_threads, nullptr);

        for (uint32_t i = 0; i < num_threads; ++i) {

        }

        for (uint32_t i = 0; i < num_threads; ++i) {
            //threadPoolTaskWait(thread_tasks[i]);
            //threadPoolTaskDelete(thread_tasks[i]);
        }

        for (int i = database_start; i < database_length; ++i) {
            database_cells += chainGetLength(database[i]);
            chainDelete(database[i]);
            database[i] = nullptr;
        }

        if (status == 0) {
            break;
        }

        database_start = database_length;
    }

    for (int i = 0; i < queries_length; ++i) {
        mutexDelete(&candidates_mutexes[i]);
    }

    fclose(handle);

    deleteFastaChains(database, database_length);
    deleteFastaChains(queries, queries_length);

    return database_cells;
}

void* threadSearchDatabase(void* params) {

   ThreadData* thread_data = static_cast<ThreadData*>(params);

   delete thread_data;

   return nullptr;
}

int32_t longestIncreasingSubsequence(const std::vector<int32_t>& src) {

    int32_t score = 0, l, u, temp;

    std::vector<int32_t> help(src.size() + 1, std::numeric_limits<int32_t>::max());
    help[0] = std::numeric_limits<int32_t>::min();

    for (uint32_t i = 0; i < src.size(); ++i) {
        l = 0;
        u = src.size();

        while (u > l) {
            temp = floor((l + u) / 2.0f);
            if (help[temp] < src[i]) {
                l = temp + 1;
            } else {
                u = temp;
            }
        }

        help[l] = src[i];
        score = score > l ? score : l;
    }

    return score;
}
