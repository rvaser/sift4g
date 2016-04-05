/*!
 * @file database_search.hpp
 *
 * @brief Database search header file
 */

#pragma once

#include <stdint.h>
#include <vector>
#include <string>

uint64_t searchDatabase(std::vector<std::vector<uint32_t>>& dst,
    const std::string& database_path, const std::string& query_path,
    uint32_t kmer_length, uint32_t max_candidates, uint32_t num_threads);
