/*!
 * @file hash.hpp
 *
 * @brief Hash class header file
 */

#pragma once

#include <stdint.h>
#include <memory>
#include <vector>

struct Chain;

void createKmerVector(std::vector<uint32_t>& dst, Chain* chain, uint32_t chain_length,
    uint32_t kmer_length);

class Hash;

std::unique_ptr<Hash> createHash(Chain** chains, uint32_t chains_length,
    uint32_t start, uint32_t length, uint32_t kmer_length);

class Hash {
public:

    ~Hash() {};

    using Iterator = std::vector<uint32_t>::iterator;
    void hits(Iterator& start, Iterator& end, uint32_t key);

    friend std::unique_ptr<Hash> createHash(Chain** chains, uint32_t chains_length,
        uint32_t start, uint32_t length, uint32_t kmer_length);

private:

    Hash(Chain** chains, uint32_t chains_length, uint32_t start, uint32_t length,
        uint32_t kmer_length);

    Hash(const Hash&) = delete;
    const Hash& operator=(const Hash&) = delete;

    std::vector<size_t> starts_;
    std::vector<uint32_t> hits_;
};
