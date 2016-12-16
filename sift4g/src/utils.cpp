/*!
 * @file utils.cpp
 *
 * @brief Utils source file
 *
 * @author: rvaser
 */

#include <sys/types.h>
#include <sys/stat.h>
#include <stdint.h>
#include <string.h>

#include "utils.hpp"

constexpr uint32_t kBufferSize = 4096;

int isExtantPath(const char* path) {
    struct stat buffer;
    if (stat(path, &buffer) == 0) {
        if (buffer.st_mode & S_IFDIR) {
            // directory
            return 0;
        } else if (buffer.st_mode & S_IFREG) {
            // file
            return 1;
        } else {
            // something
            return 2;
        }
    }
    return -1;
}

char* createFileName(const char* name, const std::string& path, const std::string& extension) {

    char* file_name = new char[kBufferSize];

    if (!path.empty()) {
        strcpy(file_name, path.c_str());
        strcat(file_name, "/");
        strcat(file_name, name);
    } else {
        strcpy(file_name, name);
    }

    strcat(file_name, extension.c_str());

    return file_name;
}

void queryLog(uint32_t part, uint32_t total) {
    fprintf(stderr, "* processing queries: %.2f/100.00%% *\r", 100 * part / (float) total);
    fflush(stderr);
}

void databaseLog(uint32_t part, float part_size, float percentage) {
    fprintf(stderr, "* processing database part %u (size ~%.2f GB): %.2f/100.00%% *\r",
        part, part_size, percentage);
    fflush(stderr);
}
