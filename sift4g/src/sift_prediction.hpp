/*!
 * @file sift_prediction.hpp
 *
 * @brief SIFT predictions header file
 */

#pragma once

#include <stdint.h>
#include <vector>
#include <string>

#include "swsharp/swsharp.h"

void siftPredictions(std::vector<std::vector<Chain*>>& alignment_strings,
    Chain** queries, int32_t queries_length, const std::string& subst_path,
    int32_t sequence_identity, const std::string& out_path);
/* this function changes alignment_strings by adding query sequence to it */
