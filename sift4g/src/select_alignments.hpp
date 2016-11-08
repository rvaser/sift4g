/*!
 * @file select_alignments.hpp
 *
 * @brief Alignment selection header file
 *
 * @author: rvaser
 */

#pragma once

#include <stdint.h>
#include <vector>
#include <string>

#include "swsharp/swsharp.h"

void selectAlignments(std::vector<std::vector<Chain*>>& dst, DbAlignment*** alignments,
    int32_t* alignments_lengths, Chain** queries, int32_t queries_length,
    float threshold);

void outputSelectedAlignments(std::vector<std::vector<Chain*>>& alignment_strings,
    Chain** queries, int32_t queries_length, const std::string& out_path);

void deleteSelectedAlignments(std::vector<std::vector<Chain*>>& alignment_strings);
