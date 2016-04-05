#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "database_search.hpp"

#include "swsharp/swsharp.h"

static struct option options[] = {
    {"query", required_argument, 0, 'i'},
    {"target", required_argument, 0, 'j'},
    {"kmer-length", required_argument, 0, 'k'},
    {"max-candidates", required_argument, 0, 'c'},
    {"threads", required_argument, 0, 't'},
    {"help", no_argument, 0, 'h'},
    {0, 0, 0, 0}
};

static void help();

int main(int argc, char* argv[]) {

    std::string query_path;
    std::string database_path;

    uint32_t kmer_length = 5;
    uint32_t max_candidates = 5000;

    uint32_t num_threads = 8;

    while (1) {

        char argument = getopt_long(argc, argv, "i:j:h", options, NULL);

        if (argument == -1) {
            break;
        }

        switch (argument) {
        case 'i':
            query_path = optarg;
            break;
        case 'j':
            database_path = optarg;
            break;
        case 'k':
            kmer_length = atoi(optarg);
            break;
        case 't':
            num_threads = atoi(optarg);
            break;
        case 'c':
            max_candidates = atoi(optarg);
            break;
        case 'h':
        default:
            help();
            return -1;
        }
    }

    assert(!query_path.empty() && "missing option -i (query file)");
    assert(!database_path.empty() && "missing option -j (database file)");

    assert(kmer_length > 2 && kmer_length < 6 && "kmer_length possible values = 3,4,5");

    threadPoolInitialize(num_threads);

    std::vector<std::vector<uint32_t>> indices;
    uint64_t cells = searchDatabase(indices, database_path, query_path, kmer_length,
        max_candidates, num_threads);

    threadPoolTerminate();

    return 0;
}

static void help() {
    printf(
    "usage: sift4g -i <query db file> -j <target db file> [arguments ...]\n"
    "\n"
    "arguments:\n"
    "    -i, --query <file>\n"
    "        (required)\n"
    "        input fasta database query file\n"
    "    -j, --target <file>\n"
    "        (required)\n"
    "        input fasta database target file\n"
    "    --kmer-length <int>\n"
    "        default: 5\n"
    "        length of kmers used for database search\n"
    "        possible values: 3, 4, 5\n"
    "    --max-candidates <int>\n"
    "        default: 5000\n"
    "        number of database sequences passed on to the Smith-Waterman part\n"
    "    --threads <int>\n"
    "        default: 8\n"
    "        number of threads used in thread pool\n"
    "    -h, -help\n"
    "        prints out the help\n");
}
