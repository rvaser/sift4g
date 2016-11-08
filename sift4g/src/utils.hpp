/*!
 * @file utils.hpp
 *
 * @brief Utils header file
 *
 * @author: rvaser
 */

#pragma once

#include <string>

bool isExtantPath(const char* path);

/* call delete[] after usage */
char* createFileName(const char* name, const std::string& path, const std::string& extension);

void queryLog(uint32_t part, uint32_t total);

void databaseLog(uint32_t part, float part_size, float percentage);
