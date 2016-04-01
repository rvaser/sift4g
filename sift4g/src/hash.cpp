/*!
 * @file hash.cpp
 *
 * @brief Hash class source file
 */

#include <assert.h>
#include <queue>

#include "hash.hpp"

#include "swsharp/swsharp.h"

constexpr uint32_t kProtMaxValue = 25;
constexpr uint32_t kProtBitLength = 5;
std::vector<uint32_t> kDelMasks = { 0, 0, 0, 0x7FFF, 0xFFFFF, 0x1FFFFFF };
std::vector<uint32_t> kNumDiffKmers = { 0, 0, 0, 26427, 845627, 27060027 };

void createKmerVector(std::vector<uint32_t>& dst, Chain* chain, uint32_t kmer_length) {

    uint32_t chain_length = chainGetLength(chain);
    if (chain_length < kmer_length) {
        return;
    }

    dst.clear();
    auto codes = chainGetCodes(chain);

    uint32_t kmer = 0;
    uint32_t del_mask = kDelMasks[kmer_length];

    for (uint32_t i = 0; i < kmer_length; ++i) {
        kmer = (kmer << kProtBitLength) | codes[i];
    }

    for (uint32_t i = kmer_length; i < chain_length; ++i) {
        kmer = ((kmer << kProtBitLength) | codes[i]) & del_mask;
        dst.emplace_back(kmer);
    }
}

std::unique_ptr<Hash> createHash(Chain** chains, uint32_t chains_length,
    uint32_t start, uint32_t length, uint32_t kmer_length) {

    assert(chains_length);
    assert(start < chains_length && start + length <= chains_length);
    assert(kmer_length && kmer_length < 6);

    return std::unique_ptr<Hash>(new Hash(chains, chains_length, start, length, kmer_length));
}

Hash::Hash(Chain** chains, uint32_t chains_length, uint32_t start, uint32_t length,
    uint32_t kmer_length)
        : starts_(kNumDiffKmers[kmer_length], 0) {

    std::vector<uint32_t> chain_kmers;
    for (uint32_t i = start; i < start + length; ++i) {

        createKmerVector(chain_kmers, chains[i], kmer_length);

        for (uint32_t j = 0; j < chain_kmers.size(); ++j) {
            ++starts_[chain_kmers[j] + 1];
        }
    }

    for (uint32_t i = 1; i < starts_.size() - 1; ++i) {
        starts_[i + 1] += starts_[i];
    }

    hits_.resize(starts_[starts_.size() - 1]);
    std::vector<uint32_t> tmp(starts_.begin(), starts_.end());

    for (uint32_t i = start; i < start + length; ++i) {

        createKmerVector(chain_kmers, chains[i], kmer_length);

        for (uint32_t j = 0; j < chain_kmers.size(); ++j) {
            hits_[tmp[chain_kmers[j]]++] = i - start;
        }
    }
}

void Hash::hits(Iterator& start, Iterator& end, uint32_t key) {
    start = hits_.begin() + starts_[key];
    end = hits_.begin() + starts_[key + 1];
}
