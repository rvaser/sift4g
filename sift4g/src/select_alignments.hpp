/*!
 * @file select_alignments.hpp
 *
 * @brief Alignment selection header file
 */

#pragma once

#include <stdint.h>
#include <vector>
#include <string>

#include "swsharp/swsharp.h"

constexpr double kLog_2_20 = 4.321928095;

void selectAlignments(std::vector<std::vector<Chain*>>& dst, DbAlignment*** alignments,
    int32_t* alignments_lengths, Chain** queries, int32_t queries_length,
    float threshold);

void outputSelectedAlignments(std::vector<std::vector<Chain*>>& alignment_strings,
    Chain** queries, int32_t queries_length, const std::string& out_path);

float getMedian(float* a, int len);

void deleteSelectedAlignments(std::vector<std::vector<Chain*>>& alignment_strings);
