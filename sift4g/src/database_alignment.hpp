/*!
 * @file database_alignment.hpp
 *
 * @brief Database alignment header file
 *
 * @author: rvaser
 */

#pragma once

#include <stdint.h>
#include <vector>
#include <string>

#include "swsharp/evalue.h"
#include "swsharp/swsharp.h"

void alignDatabase(DbAlignment**** alignments, int** alignments_lengths, Chain*** database,
    int32_t* database_length, const std::string& database_path, Chain** queries,
    int32_t queries_length, std::vector<std::vector<uint32_t>>& indices,
    int32_t algorithm, EValueParams* evalue_params, double max_evalue,
    uint32_t max_alignments, Scorer* scorer, int32_t* cards,
    int32_t cards_length);
