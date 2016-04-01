#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "hash.hpp"

#include "swsharp/swsharp.h"

static struct option options[] = {
    {"query", required_argument, 0, 'i'},
    {"target", required_argument, 0, 'j'},
    {"kmer-length", required_argument, 0, 'k'},
    {"help", no_argument, 0, 'h'},
    {0, 0, 0, 0}
};

static void help();

int main(int argc, char* argv[]) {

    std::string query_path;
    std::string database_path;

    int kmer_length = 5;

    while (1) {

        char argument = getopt_long(argc, argv, "i:j:k:h", options, NULL);

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
        case 'h':
        default:
            help();
            return -1;
        }
    }

    assert(!query_path.empty() && "missing option -i (query file)");
    assert(!database_path.empty() && "missing option -j (database file)");

    assert(kmer_length > 2 && kmer_length < 6 && "kmer_length possible values = 3,4,5");

    Chain** queries = nullptr;
    int queries_length = 0;
    readFastaChains(&queries, &queries_length, query_path.c_str());

    auto hash = createHash(queries, queries_length, 0, queries_length, kmer_length);

    deleteFastaChains(queries, queries_length);

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
    "    -h, -help\n"
    "        prints out the help\n");
}
