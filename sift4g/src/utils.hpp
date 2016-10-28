/*!
 * @file utils.hpp
 *
 * @brief Utils header file
 */

#pragma once

#include <string>

bool exists(const char* path);

/* call delete[] after usage */
char* createFileName(const char* name, const std::string& path, const std::string& extension);

void query_log(uint32_t part, uint32_t total);

void database_log(uint32_t part, float part_size, float percentage);
