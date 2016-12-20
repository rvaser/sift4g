/*!
 * @file utils.hpp
 *
 * @brief Utils header file
 *
 * @author: rvaser
 */

#pragma once

#include <string>

#define ASSERT(expr, fmt, ...)\
    do {\
        if (!(expr)) {\
            fprintf(stderr, "[ERROR]: " fmt "\n", ##__VA_ARGS__);\
            exit(-1);\
        }\
    } while(0)

int isExtantPath(const char* path);

/* call delete[] after usage */
char* createFileName(const char* name, const std::string& path, const std::string& extension);

void queryLog(uint32_t part, uint32_t total);

void databaseLog(uint32_t part, float part_size, float percentage);
