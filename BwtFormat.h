#ifndef BWT_FORMAT_H
#define BWT_FORMAT_H

#include <stdint.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>

#include <algorithm>
#include <stdexcept>
#include <string>
#include <vector>

namespace bwt_format {

enum Encoding {
    PLAIN = 0,
    RLE = 1,
    PACKED = 2
};

static const unsigned char MAGIC[8] = {'E', 'B', 'W', 'T', 'C', 'M', 'P', '1'};
static const unsigned char VERSION = 1;
static const size_t ADAPTIVE_BLOCK_BYTES = 1024 * 1024;

inline void fail(const std::string &message) {
    throw std::runtime_error(message);
}

inline void writeExact(FILE *file, const void *data, size_t size, const std::string &path) {
    if (size != 0 && fwrite(data, 1, size, file) != size)
        fail("cannot write compressed BWT: " + path);
}

inline void putU64LE(unsigned char *dst, uint64_t value) {
    for (unsigned int i = 0; i < 8; ++i) {
        dst[i] = static_cast<unsigned char>(value & 0xffU);
        value >>= 8;
    }
}

inline void putU32LE(unsigned char *dst, uint32_t value) {
    for (unsigned int i = 0; i < 4; ++i) {
        dst[i] = static_cast<unsigned char>(value & 0xffU);
        value >>= 8;
    }
}

inline Encoding parseEncoding(const char *value) {
    if (value == NULL || strcmp(value, "rle") == 0)
        return RLE;
    if (strcmp(value, "packed") == 0)
        return PACKED;
    if (strcmp(value, "plain") == 0)
        return PLAIN;
    fail("BCR_BWT_FORMAT must be one of: rle, packed, plain");
    return RLE;
}

inline const char *suffix(Encoding encoding) {
    if (encoding == RLE)
        return ".rl_bwt";
    if (encoding == PACKED)
        return ".pck_bwt";
    return ".ebwt";
}

class Writer {
public:
    Writer(const std::string &path, Encoding encoding, unsigned char terminator)
        : path_(path), file_(NULL), encoding_(encoding), terminator_(terminator),
          length_(0), checksum_(crc32(0L, Z_NULL, 0)), closed_(false),
          bit_buffer_(0), bit_count_(0) {
        file_ = fopen(path.c_str(), "wb");
        if (file_ == NULL)
            fail("cannot create BWT output: " + path);
        if (encoding_ != PLAIN)
            writeHeader(0, 0);
        if (encoding_ == RLE)
            adaptive_block_.reserve(ADAPTIVE_BLOCK_BYTES);
    }

    ~Writer() {
        if (file_ != NULL)
            fclose(file_);
    }

    void write(const unsigned char *data, size_t size) {
        if (closed_)
            fail("attempted to write a closed BWT output");
        for (size_t i = 0; i < size; ++i)
            validateSymbol(data[i]);
        updateMetadata(data, size);
        if (encoding_ == PLAIN) {
            writeExact(file_, data, size, path_);
            return;
        }
        if (encoding_ == RLE) {
            size_t offset = 0;
            while (offset < size) {
                const size_t available = ADAPTIVE_BLOCK_BYTES - adaptive_block_.size();
                const size_t take = std::min(available, size - offset);
                adaptive_block_.insert(
                    adaptive_block_.end(), data + offset, data + offset + take);
                offset += take;
                if (adaptive_block_.size() == ADAPTIVE_BLOCK_BYTES)
                    flushAdaptiveBlock();
            }
        } else {
            for (size_t i = 0; i < size; ++i)
                writePackedSymbol(data[i]);
        }
    }

    void close() {
        if (closed_)
            return;
        if (encoding_ == RLE)
            flushAdaptiveBlock();
        else if (encoding_ == PACKED)
            flushPacked();
        if (encoding_ != PLAIN) {
            if (fseek(file_, 0, SEEK_SET) != 0)
                fail("cannot finalize compressed BWT header: " + path_);
            writeHeader(length_, checksum_);
        }
        if (fclose(file_) != 0) {
            file_ = NULL;
            fail("cannot close BWT output: " + path_);
        }
        file_ = NULL;
        closed_ = true;
    }

    uint64_t length() const { return length_; }

private:
    std::string path_;
    FILE *file_;
    Encoding encoding_;
    unsigned char terminator_;
    uint64_t length_;
    uLong checksum_;
    bool closed_;
    uint64_t bit_buffer_;
    unsigned int bit_count_;
    std::vector<unsigned char> adaptive_block_;
    std::vector<unsigned char> rle_payload_;
    std::vector<unsigned char> packed_payload_;

    void validateSymbol(unsigned char symbol) const {
        if (symbol != terminator_ && symbol != 'A' && symbol != 'C' &&
            symbol != 'G' && symbol != 'T')
            fail("cannot encode symbol outside A,C,G,T and the configured terminator");
    }

    void updateMetadata(const unsigned char *data, size_t size) {
        size_t offset = 0;
        while (offset < size) {
            const size_t left = size - offset;
            const uInt chunk = left > UINT_MAX ? UINT_MAX : static_cast<uInt>(left);
            checksum_ = crc32(checksum_, data + offset, chunk);
            offset += chunk;
        }
        length_ += static_cast<uint64_t>(size);
    }

    void writeHeader(uint64_t length, uint64_t checksum) {
        unsigned char header[32];
        memset(header, 0, sizeof(header));
        memcpy(header, MAGIC, sizeof(MAGIC));
        header[8] = VERSION;
        header[9] = static_cast<unsigned char>(encoding_);
        header[10] = terminator_;
        putU64LE(header + 16, length);
        putU64LE(header + 24, checksum);
        writeExact(file_, header, sizeof(header), path_);
    }

    void appendVarint(std::vector<unsigned char> &output, uint64_t value) {
        do {
            unsigned char byte = static_cast<unsigned char>(value & 0x7fU);
            value >>= 7;
            if (value != 0)
                byte |= 0x80U;
            output.push_back(byte);
        } while (value != 0);
    }

    void encodeRleBlock() {
        rle_payload_.clear();
        rle_payload_.reserve(adaptive_block_.size() / 2 + 16);
        size_t start = 0;
        while (start < adaptive_block_.size()) {
            size_t end = start + 1;
            while (end < adaptive_block_.size() && adaptive_block_[end] == adaptive_block_[start])
                ++end;
            rle_payload_.push_back(adaptive_block_[start]);
            appendVarint(rle_payload_, static_cast<uint64_t>(end - start));
            start = end;
        }
    }

    void encodePackedBlock() {
        packed_payload_.assign((adaptive_block_.size() * 3 + 7) / 8, 0);
        uint64_t bits = 0;
        unsigned int bitCount = 0;
        size_t output = 0;
        for (size_t i = 0; i < adaptive_block_.size(); ++i) {
            bits |= static_cast<uint64_t>(packedCode(adaptive_block_[i])) << bitCount;
            bitCount += 3;
            while (bitCount >= 8) {
                packed_payload_[output++] = static_cast<unsigned char>(bits & 0xffU);
                bits >>= 8;
                bitCount -= 8;
            }
        }
        if (bitCount != 0)
            packed_payload_[output] = static_cast<unsigned char>(bits & 0xffU);
    }

    void flushAdaptiveBlock() {
        if (adaptive_block_.empty())
            return;
        encodeRleBlock();
        encodePackedBlock();
        const bool usePacked = packed_payload_.size() < rle_payload_.size();
        const std::vector<unsigned char> &payload = usePacked ? packed_payload_ : rle_payload_;
        if (adaptive_block_.size() > UINT32_MAX || payload.size() > 0x7fffffffU)
            fail("adaptive BWT block is too large");
        unsigned char blockHeader[8];
        putU32LE(blockHeader, static_cast<uint32_t>(adaptive_block_.size()));
        uint32_t storedMode = static_cast<uint32_t>(payload.size());
        if (usePacked)
            storedMode |= 0x80000000U;
        putU32LE(blockHeader + 4, storedMode);
        writeExact(file_, blockHeader, sizeof(blockHeader), path_);
        writeExact(file_, &payload[0], payload.size(), path_);
        adaptive_block_.clear();
    }

    unsigned int packedCode(unsigned char symbol) const {
        if (symbol == terminator_)
            return 0;
        if (symbol == 'A')
            return 1;
        if (symbol == 'C')
            return 2;
        if (symbol == 'G')
            return 3;
        return 4;
    }

    void writePackedSymbol(unsigned char symbol) {
        bit_buffer_ |= static_cast<uint64_t>(packedCode(symbol)) << bit_count_;
        bit_count_ += 3;
        while (bit_count_ >= 8) {
            const unsigned char byte = static_cast<unsigned char>(bit_buffer_ & 0xffU);
            writeExact(file_, &byte, 1, path_);
            bit_buffer_ >>= 8;
            bit_count_ -= 8;
        }
    }

    void flushPacked() {
        if (bit_count_ != 0) {
            const unsigned char byte = static_cast<unsigned char>(bit_buffer_ & 0xffU);
            writeExact(file_, &byte, 1, path_);
            bit_buffer_ = 0;
            bit_count_ = 0;
        }
    }
};

} // namespace bwt_format

#endif
