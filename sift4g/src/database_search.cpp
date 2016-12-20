/*!
 * @file database_search.cpp
 *
 * @brief Database search source file
 *
 * @author: rvaser
 */

#include <algorithm>
#include <limits>

#include "hash.hpp"
#include "utils.hpp"
#include "database_search.hpp"

constexpr uint32_t database_chunk = 250000000; /* ~250MB */
constexpr float log_step_percentage = 2.5;

class Candidate {
public:
    Candidate(float _score, int _id) :
            score(_score), id(_id) {
    }

    bool operator<(const Candidate& other) const {
        return this->score > other.score;
    }

    float score;
    int32_t id;
};

class ThreadSearchData {
public:
    ThreadSearchData(std::shared_ptr<Hash> _query_hash, uint32_t _queries_length,
        std::vector<float>& _min_scores, Chain** _database,
        uint32_t _database_begin, uint32_t _database_end,
        uint32_t _kmer_length, uint32_t _max_candidates,
        std::vector<std::vector<Candidate>>& _candidates,
        bool _log, uint32_t _part, float _part_size):
            query_hash(_query_hash), queries_length(_queries_length), min_scores(_min_scores),
            database(_database), database_begin(_database_begin), database_end(_database_end),
            kmer_length(_kmer_length), max_candidates(_max_candidates), candidates(_candidates),
            log(_log), part(_part), part_size(_part_size) {
    }

    std::shared_ptr<Hash> query_hash;
    uint32_t queries_length;
    std::vector<float>& min_scores;
    Chain** database;
    uint32_t database_begin;
    uint32_t database_end;
    uint32_t kmer_length;
    uint32_t max_candidates;
    std::vector<std::vector<Candidate>>& candidates;
    bool log;
    uint32_t part;
    float part_size;
};

void* threadSearchDatabase(void* params);

int32_t longestIncreasingSubsequence(const std::vector<int32_t>& src);

uint64_t searchDatabase(std::vector<std::vector<uint32_t>>& dst,
    const std::string& database_path, Chain** queries, int32_t queries_length,
    uint32_t kmer_length, uint32_t max_candidates, uint32_t num_threads) {

    fprintf(stderr, "** Searching database for candidate sequences **\n");

    std::shared_ptr<Hash> query_hash = createHash(queries, queries_length, 0,
        queries_length, kmer_length);

    Chain** database = nullptr;
    int database_length = 0;
    int database_start = 0;

    FILE* handle = nullptr;
    int serialized = 0;
    readFastaChainsPartInit(&database, &database_length, &handle, &serialized,
        database_path.c_str());

    uint64_t database_cells = 0;

    std::vector<float> min_scores(queries_length, 1000000.0);
    std::vector<std::vector<std::vector<Candidate>>> candidates(num_threads);

    uint32_t part = 1;
    float part_size = database_chunk / (float) 1000000000;

    while (true) {

        int status = 1;

        status &= readFastaChainsPart(&database, &database_length, handle,
            serialized, database_chunk);

        databaseLog(part, part_size, 0);

        uint32_t database_split_size = (database_length - database_start) / num_threads;
        std::vector<uint32_t> database_splits(num_threads + 1, database_start);
        for (uint32_t i = 1; i < num_threads; ++i) {
            database_splits[i] += i * database_split_size;
        }
        database_splits[num_threads] = database_length;

        std::vector<ThreadPoolTask*> thread_tasks(num_threads, nullptr);

        for (uint32_t i = 0; i < num_threads; ++i) {

            auto thread_data = new ThreadSearchData(query_hash, queries_length, min_scores,
                database, database_splits[i], database_splits[i + 1], kmer_length,
                max_candidates, candidates[i], i == num_threads - 1, part,
                part_size);

            thread_tasks[i] = threadPoolSubmit(threadSearchDatabase, (void*) thread_data);
        }

        for (uint32_t i = 0; i < num_threads; ++i) {
            threadPoolTaskWait(thread_tasks[i]);
            threadPoolTaskDelete(thread_tasks[i]);
        }

        for (int i = database_start; i < database_length; ++i) {
            database_cells += chainGetLength(database[i]);
            chainDelete(database[i]);
            database[i] = nullptr;
        }

        // merge candidates from all threads
        for (int32_t i = 0; i < queries_length; ++i) {
            for (uint32_t j = 1; j < num_threads; ++j) {
                if (candidates[j][i].empty()) {
                    continue;
                }
                candidates[0][i].insert(candidates[0][i].end(),
                    candidates[j][i].begin(), candidates[j][i].end());
                std::vector<Candidate>().swap(candidates[j][i]);
            }

            if (num_threads > 1) {
                std::sort(candidates[0][i].begin(), candidates[0][i].end());
                if (candidates[0][i].size() > max_candidates) {
                    std::vector<Candidate> tmp(candidates[0][i].begin(),
                        candidates[0][i].begin() + max_candidates);
                    candidates[0][i].swap(tmp);
                }
            }

            min_scores[i] = candidates[0][i].back().score;
        }

        databaseLog(part, part_size, 100);
        ++part;

        if (status == 0) {
            break;
        }

        database_start = database_length;
    }
    fprintf(stderr, "\n\n");

    fclose(handle);
    deleteFastaChains(database, database_length);

    dst.clear();
    dst.resize(queries_length);

    for (int32_t i = 0; i < queries_length; ++i) {
        dst[i].reserve(candidates[0][i].size());
        for (uint32_t j = 0; j < candidates[0][i].size(); ++j) {
            dst[i].emplace_back(candidates[0][i][j].id);
        }
        std::vector<Candidate>().swap(candidates[0][i]);
        std::sort(dst[i].begin(), dst[i].end());
    }

    return database_cells;
}

void* threadSearchDatabase(void* params) {

    auto thread_data = (ThreadSearchData*) params;

    thread_data->candidates.resize(thread_data->queries_length);

    std::vector<uint32_t> kmer_vector;
    std::vector<std::vector<int32_t>> hits(thread_data->queries_length);
    std::vector<float> min_scores(thread_data->min_scores);

    uint32_t log_counter = 0;
    uint32_t log_size = (thread_data->database_end - thread_data->database_begin) / (100. / log_step_percentage);
    float log_percentage = log_step_percentage;

    for (uint32_t i = thread_data->database_begin; i < thread_data->database_end; ++i) {

        if (thread_data->log && log_percentage < 100.0) {
            ++log_counter;
            if (log_size != 0 && log_counter % log_size == 0) {
                databaseLog(thread_data->part, thread_data->part_size, log_percentage);
                log_percentage += log_step_percentage;
            }
        }

        createKmerVector(kmer_vector, thread_data->database[i], thread_data->kmer_length);

        for (uint32_t j = 0; j < kmer_vector.size(); ++j) {
            if (j != 0 && kmer_vector[j] == kmer_vector[j - 1]) {
                continue;
            }

            Hash::Iterator begin, end;
            thread_data->query_hash->hits(begin, end, kmer_vector[j]);
            for (; begin != end; ++begin) {
                hits[begin->id].emplace_back(begin->position);
            }
        }

        for (uint32_t j = 0; j < thread_data->queries_length; ++j) {
            if (hits[j].empty()) {
                continue;
            }

            float similartiy_score = longestIncreasingSubsequence(hits[j]) /
                (float) chainGetLength(thread_data->database[i]);

            if (thread_data->candidates[j].size() < thread_data->max_candidates || similartiy_score > min_scores[j]) {
                thread_data->candidates[j].emplace_back(similartiy_score, i);
                min_scores[j] = std::min(min_scores[j], similartiy_score);
            }

            std::vector<int32_t>().swap(hits[j]);
        }
    }

    for (uint32_t i = 0; i < thread_data->queries_length; ++i) {
        std::sort(thread_data->candidates[i].begin(), thread_data->candidates[i].end());

        if (thread_data->candidates[i].size() > thread_data->max_candidates) {
            std::vector<Candidate> tmp(thread_data->candidates[i].begin(),
                thread_data->candidates[i].begin() + thread_data->max_candidates);
            thread_data->candidates[i].swap(tmp);
        }
    }

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
