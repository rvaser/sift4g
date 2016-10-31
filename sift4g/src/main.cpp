#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "utils.hpp"
#include "database_search.hpp"
#include "database_alignment.hpp"
#include "select_alignments.hpp"
#include "sift_prediction.hpp"

#include "swsharp/evalue.h"
#include "swsharp/swsharp.h"

#define CHAR_INT_LEN(x) (sizeof(x) / sizeof(CharInt))
typedef struct CharInt {
    const char* format;
    const int code;
} CharInt;

static struct option options[] = {
    {"query", required_argument, 0, 'i'},
    {"target", required_argument, 0, 'j'},
    {"kmer-length", required_argument, 0, 'k'},
    {"max-candidates", required_argument, 0, 'C'},
    {"median-threshold", required_argument, 0, 'T'},
    {"cards", required_argument, 0, 'c'},
    {"gap-extend", required_argument, 0, 'e'},
    {"gap-open", required_argument, 0, 'g'},
    {"matrix", required_argument, 0, 'm'},
    {"out", required_argument, 0, 'o'},
    {"subst", required_argument, 0, 'S'},
    {"seq-id", required_argument, 0, 'I'},
    {"sub-results", no_argument, 0, 's'},
    {"outfmt", required_argument, 0, 'f'},
    {"evalue", required_argument, 0, 'E'},
    {"max-aligns", required_argument, 0, 'M'},
    {"algorithm", required_argument, 0, 'A'},
    {"threads", required_argument, 0, 't'},
    {"help", no_argument, 0, 'h'},
    {0, 0, 0, 0}
};

static CharInt outFormats[] = {
    { "bm0", SW_OUT_DB_BLASTM0 },
    { "bm8", SW_OUT_DB_BLASTM8 },
    { "bm9", SW_OUT_DB_BLASTM9 },
    { "light", SW_OUT_DB_LIGHT }
};

static CharInt algorithms[] = {
    { "SW", SW_ALIGN },
    { "NW", NW_ALIGN },
    { "HW", HW_ALIGN },
    { "OV", OV_ALIGN }
};

static void getCudaCards(int** cards, int* cardsLen, char* optarg);
static int getOutFormat(char* optarg);
static int getAlgorithm(char* optarg);
static void help();

int main(int argc, char* argv[]) {

    std::string query_path;
    std::string database_path;

    uint32_t kmer_length = 5;
    uint32_t max_candidates = 5000;

    int32_t gap_open = 10;
    int32_t gap_extend = 1;
    char* matrix = "BLOSUM_62";

    uint32_t max_alignments = 400;
    double max_evalue = 0.0001;

    int32_t* cards = nullptr;
    int32_t cards_length = 0;

    std::string out_path = "";
    bool sub_results = false;
    int32_t out_format = SW_OUT_DB_BLASTM9;

    int32_t algorithm = SW_ALIGN;

    float median_threshold = 2.75;
    std::string subst_path = "";
    int32_t sequence_identity = 100;

    uint32_t num_threads = 8;

    while (1) {

        char argument = getopt_long(argc, argv, "i:j:g:e:t:h", options, NULL);

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
        case 'C':
            max_candidates = atoi(optarg);
            break;
        case 'T':
            median_threshold = atof(optarg);
            break;
        case 'g':
            gap_open = atoi(optarg);
            break;
        case 'e':
            gap_extend = atoi(optarg);
            break;
        case 'c':
            getCudaCards(&cards, &cards_length, optarg);
            break;
        case 'o':
            out_path = optarg;
            break;
        case 'S':
            subst_path = optarg;
            break;
        case 's':
            sub_results = true;
            break;
        case 'I':
            sequence_identity = atoi(optarg);
            break;
        case 'f':
            out_format = getOutFormat(optarg);
            break;
        case 'M':
            max_alignments = atoi(optarg);
            break;
        case 'E':
            max_evalue = atof(optarg);
            break;
        case 'm':
            matrix = optarg;
            break;
        case 'A':
            algorithm = getAlgorithm(optarg);
            break;
        case 't':
            num_threads = atoi(optarg);
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
    assert(max_candidates > 0 && "invalid max candidates number");

    assert(max_evalue > 0 && "invalid evalue");
    assert(num_threads > 0 && "invalid thread number");

    if (!out_path.empty()) {
        assert(exists(out_path.c_str()) && "invalid out directory");
    }

    if (!subst_path.empty()) {
        assert(exists(subst_path.c_str()) && "invalid substitutions directory");
    }

    if (cards_length == -1) {
        cudaGetCards(&cards, &cards_length);
    }
    assert(cudaCheckCards(cards, cards_length) && "invalid cuda cards");

    threadPoolInitialize(num_threads);

    Chain** queries = nullptr;
    int32_t queries_length = 0;
    readFastaChains(&queries, &queries_length, query_path.c_str());

    std::vector<std::vector<uint32_t>> indices;
    uint64_t cells = searchDatabase(indices, database_path, queries, queries_length,
        kmer_length, max_candidates, num_threads);

    Scorer* scorer = nullptr;
    scorerCreateMatrix(&scorer, matrix, gap_open, gap_extend);

    EValueParams* evalue_params = createEValueParams(cells, scorer);

    DbAlignment*** alignments = nullptr;
    int* alignments_lenghts = nullptr;

    Chain** database = nullptr;
    int32_t database_length = 0;

    alignDatabase(&alignments, &alignments_lenghts, &database, &database_length,
        database_path, queries, queries_length, indices, algorithm, evalue_params,
        max_evalue, max_alignments, scorer, cards, cards_length);

    deleteEValueParams(evalue_params);
    scorerDelete(scorer);

    if (sub_results) {
        char* alignments_path = createFileName("alignments", out_path, ".txt");
        outputShotgunDatabase(alignments, alignments_lenghts, queries_length, alignments_path, out_format);
        delete[] alignments_path;
    }

    std::vector<std::vector<Chain*>> alignment_strings;
    selectAlignments(alignment_strings, alignments, alignments_lenghts, queries, queries_length, median_threshold);

    deleteShotgunDatabase(alignments, alignments_lenghts, queries_length);
    deleteFastaChains(database, database_length);

    if (sub_results) {
        outputSelectedAlignments(alignment_strings, queries, queries_length, out_path);
    }

    siftPredictions (alignment_strings, queries, 
				queries_length, subst_path,
			        sequence_identity, out_path);
    deleteSelectedAlignments(alignment_strings);
/*    The queries have been combined with alignments in siftPredictions
	deleteFastaChains(queries, queries_length);
*/
    threadPoolTerminate();

    free(cards);

    return 0;
}

static void getCudaCards(int** cards, int* cardsLen, char* optarg) {

    *cardsLen = strlen(optarg);
    *cards = (int*) malloc(*cardsLen * sizeof(int));

    for (int32_t i = 0; i < *cardsLen; ++i) {
        (*cards)[i] = optarg[i] - '0';
    }
}

static int getOutFormat(char* optarg) {

    for (uint32_t i = 0; i < CHAR_INT_LEN(outFormats); ++i) {
        if (strcmp(outFormats[i].format, optarg) == 0) {
            return outFormats[i].code;
        }
    }

    assert(0 && "unknown out format");
}

static int getAlgorithm(char* optarg) {

    for (uint32_t i = 0; i < CHAR_INT_LEN(algorithms); ++i) {
        if (strcmp(algorithms[i].format, optarg) == 0) {
            return algorithms[i].code;
        }
    }

    assert(0 && "unknown algorithm");
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
    "    -g, --gap-open <int>\n"
    "        default: 10\n"
    "        gap opening penalty, must be given as a positive integer \n"
    "    -e, --gap-extend <int>\n"
    "        default: 1\n"
    "        gap extension penalty, must be given as a positive integer and\n"
    "        must be less or equal to gap opening penalty\n"
    "    --matrix <string>\n"
    "        default: BLOSUM_62\n"
    "        similarity matrix, can be one of the following:\n"
    "            BLOSUM_45\n"
    "            BLOSUM_50\n"
    "            BLOSUM_62\n"
    "            BLOSUM_80\n"
    "            BLOSUM_90\n"
    "            BLOSUM_30\n"
    "            BLOSUM_70\n"
    "            BLOSUM_250\n"
    "    --evalue <float>\n"
    "        default: 0.0001\n"
    "        evalue threshold, alignments with higher evalue are filtered,\n"
    "        must be given as a positive float\n"
    "    --max-aligns <int>\n"
    "        default: 400\n"
    "        maximum number of alignments to be outputted\n"
    "    --algorithm <string>\n"
    "        default: SW\n"
    "        algorithm used for alignment, must be one of the following: \n"
    "            SW - Smith-Waterman local alignment\n"
    "            NW - Needleman-Wunsch global alignment\n"
    "            HW - semiglobal alignment\n"
    "            OV - overlap alignment\n"
    "    --cards <ints>\n"
    "        default: all available CUDA cards\n"
    "        list of cards should be given as an array of card indexes delimited with\n"
    "        nothing, for example usage of first two cards is given as --cards 01\n"
    "    --out <string>\n"
    "        default: current directory\n"
    "        output directory for SIFT predictions\n"
    "    --sub-results\n"
    "        prints sub results (alignment file and a file per query containing\n"
    "        its selected alignments forp rediction) to same directory defined\n"
    "        with --out\n"
    "    --outfmt <string>\n"
    "        default: bm9\n"
    "        out format for the alignment file, must be one of the following:\n"
    "            bm0      - blast m0 output format\n"
    "            bm8      - blast m8 tabular output format\n"
    "            bm9      - blast m9 commented tabular output format\n"
    "            light    - score-name tabbed output\n"
    "    --kmer-length <int>\n"
    "        default: 5\n"
    "        length of kmers used for database search\n"
    "        possible values: 3, 4, 5\n"
    "    --max-candidates <int>\n"
    "        default: 5000\n"
    "        number of database sequences passed on to the Smith-Waterman part\n"
    "    --median-threshold <float>\n"
    "        default: 2.75\n"
    "        represents alignment diversity, used to output only a set of alignments\n"
    "    --subst <string>\n"
    "        default: current directory\n"
    "        directory containing substitution files for each query (extension .subst)\n"
    "        files must have the same name as their corresponding query in FASTA file\n"
    "    --seq-id <int>\n"
    "        default: 100\n"
    "    -t, --threads <int>\n"
    "        default: 8\n"
    "        number of threads used in thread pool\n"
    "    -h, -help\n"
    "        prints out the help\n");
}
