#ifndef BCR_COMPRESSED_IO_H
#define BCR_COMPRESSED_IO_H

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

#ifndef BCR_COMPRESSED_CYC
#define BCR_COMPRESSED_CYC 0
#endif

#ifndef BCR_COMPRESSED_PARTIALS
#define BCR_COMPRESSED_PARTIALS 0
#endif

#if BCR_COMPRESSED_CYC
typedef gzFile BCRCycFile;
#else
typedef FILE *BCRCycFile;
#endif

#if BCR_COMPRESSED_PARTIALS
struct BCRPackedPartial {
    FILE *file;
    bool writing;
    uint64_t decodedLength;
    uint64_t processed;
    std::vector<unsigned char> block;
    std::vector<unsigned char> packed;
    size_t blockPosition;
    size_t blockSize;
    std::string path;
};
typedef BCRPackedPartial *BCRPartialFile;
#else
typedef FILE *BCRPartialFile;
#endif

inline void bcrIoFail(const char *operation, const char *path) {
    std::cerr << operation << ": " << path << std::endl;
    exit(EXIT_FAILURE);
}

inline BCRCycFile bcrOpenCyc(const char *path, const char *mode) {
#if BCR_COMPRESSED_CYC
    const char *compressedMode = mode[0] == 'a' ? "ab1" : (mode[0] == 'w' ? "wb1" : "rb");
    BCRCycFile file = gzopen(path, compressedMode);
#else
    BCRCycFile file = fopen(path, mode);
#endif
    if (file == NULL)
        bcrIoFail("cannot open cyclic-input file", path);
    return file;
}

inline size_t bcrReadCyc(void *data, size_t size, BCRCycFile file) {
#if BCR_COMPRESSED_CYC
    size_t total = 0;
    unsigned char *out = static_cast<unsigned char *>(data);
    while (total < size) {
        const size_t left = size - total;
        const unsigned int chunk = left > (1U << 30) ? (1U << 30) : static_cast<unsigned int>(left);
        const int got = gzread(file, out + total, chunk);
        if (got <= 0)
            break;
        total += static_cast<size_t>(got);
    }
    return total;
#else
    return fread(data, 1, size, file);
#endif
}

inline size_t bcrWriteCyc(const void *data, size_t size, BCRCycFile file) {
#if BCR_COMPRESSED_CYC
    size_t total = 0;
    const unsigned char *in = static_cast<const unsigned char *>(data);
    while (total < size) {
        const size_t left = size - total;
        const unsigned int chunk = left > (1U << 30) ? (1U << 30) : static_cast<unsigned int>(left);
        const int wrote = gzwrite(file, in + total, chunk);
        if (wrote <= 0)
            break;
        total += static_cast<size_t>(wrote);
    }
    return total;
#else
    return fwrite(data, 1, size, file);
#endif
}

inline int bcrCloseCyc(BCRCycFile file) {
#if BCR_COMPRESSED_CYC
    return gzclose(file) == Z_OK ? 0 : EOF;
#else
    return fclose(file);
#endif
}

#if BCR_COMPRESSED_PARTIALS
static const unsigned char BCR_PARTIAL_MAGIC[8] = {'B', 'C', 'R', 'P', 'K', 'B', '0', '1'};
static const size_t BCR_PARTIAL_BLOCK_BYTES = 1024 * 1024;

inline void bcrStoreU32(unsigned char *dst, uint32_t value) {
    for (unsigned int i = 0; i < 4; ++i) {
        dst[i] = static_cast<unsigned char>(value & 0xffU);
        value >>= 8;
    }
}

inline uint32_t bcrLoadU32(const unsigned char *src) {
    uint32_t value = 0;
    for (int i = 3; i >= 0; --i)
        value = (value << 8) | src[i];
    return value;
}

inline void bcrStoreU64(unsigned char *dst, uint64_t value) {
    for (unsigned int i = 0; i < 8; ++i) {
        dst[i] = static_cast<unsigned char>(value & 0xffU);
        value >>= 8;
    }
}

inline uint64_t bcrLoadU64(const unsigned char *src) {
    uint64_t value = 0;
    for (int i = 7; i >= 0; --i)
        value = (value << 8) | src[i];
    return value;
}

struct BCRPackedTables {
    uint16_t encodePair[65536];
    uint32_t decodeFour[4096];
    bool decodeValid[4096];

    BCRPackedTables() {
        for (size_t i = 0; i < 65536; ++i)
            encodePair[i] = UINT16_MAX;
        const unsigned char symbols[6] = {'$', 'A', 'C', 'G', 'T', 'N'};
        for (unsigned int first = 0; first < 6; ++first)
            for (unsigned int second = 0; second < 6; ++second)
                encodePair[(static_cast<unsigned int>(symbols[first]) << 8) | symbols[second]] =
                    static_cast<uint16_t>(first | (second << 3));
        for (unsigned int packed = 0; packed < 4096; ++packed) {
            unsigned int bits = packed;
            uint32_t word = 0;
            bool valid = true;
            for (unsigned int i = 0; i < 4; ++i) {
                const unsigned int code = bits & 0x7U;
                bits >>= 3;
                if (code >= 6) {
                    valid = false;
                    break;
                }
                word |= static_cast<uint32_t>(symbols[code]) << (8 * i);
            }
            decodeFour[packed] = word;
            decodeValid[packed] = valid;
        }
    }
};

inline BCRPackedTables &bcrPackedTables() {
    static BCRPackedTables tables;
    return tables;
}

inline unsigned int bcrPackedCode(unsigned char symbol, const std::string &path) {
    if (symbol == '$') return 0;
    if (symbol == 'A') return 1;
    if (symbol == 'C') return 2;
    if (symbol == 'G') return 3;
    if (symbol == 'T') return 4;
    if (symbol == 'N') return 5;
    bcrIoFail("partial BWT contains a symbol unsupported by packed mode", path.c_str());
    return 0;
}

inline void bcrFlushPartialBlock(BCRPartialFile file) {
    if (file->blockSize == 0)
        return;
    const size_t packedSize = (file->blockSize * 3 + 7) / 8;
    file->packed.assign(packedSize, 0);
    BCRPackedTables &tables = bcrPackedTables();
    size_t input = 0;
    size_t output = 0;
    while (file->blockSize - input >= 8) {
        const uint16_t pair0 = tables.encodePair[
            (static_cast<unsigned int>(file->block[input]) << 8) | file->block[input + 1]];
        const uint16_t pair1 = tables.encodePair[
            (static_cast<unsigned int>(file->block[input + 2]) << 8) | file->block[input + 3]];
        const uint16_t pair2 = tables.encodePair[
            (static_cast<unsigned int>(file->block[input + 4]) << 8) | file->block[input + 5]];
        const uint16_t pair3 = tables.encodePair[
            (static_cast<unsigned int>(file->block[input + 6]) << 8) | file->block[input + 7]];
        if (pair0 == UINT16_MAX || pair1 == UINT16_MAX || pair2 == UINT16_MAX || pair3 == UINT16_MAX)
            bcrIoFail("partial BWT contains a symbol unsupported by packed mode", file->path.c_str());
        const uint32_t bits = static_cast<uint32_t>(pair0) |
            (static_cast<uint32_t>(pair1) << 6) |
            (static_cast<uint32_t>(pair2) << 12) |
            (static_cast<uint32_t>(pair3) << 18);
        file->packed[output] = static_cast<unsigned char>(bits & 0xffU);
        file->packed[output + 1] = static_cast<unsigned char>((bits >> 8) & 0xffU);
        file->packed[output + 2] = static_cast<unsigned char>((bits >> 16) & 0xffU);
        input += 8;
        output += 3;
    }
    uint32_t tail = 0;
    unsigned int tailBits = 0;
    while (input < file->blockSize) {
        tail |= static_cast<uint32_t>(bcrPackedCode(file->block[input++], file->path)) << tailBits;
        tailBits += 3;
    }
    while (tailBits != 0) {
        file->packed[output++] = static_cast<unsigned char>(tail & 0xffU);
        tail >>= 8;
        tailBits = tailBits > 8 ? tailBits - 8 : 0;
    }
    unsigned char header[4];
    bcrStoreU32(header, static_cast<uint32_t>(file->blockSize));
    if (fwrite(header, 1, sizeof(header), file->file) != sizeof(header) ||
        fwrite(&file->packed[0], 1, packedSize, file->file) != packedSize)
        bcrIoFail("cannot write packed partial-BWT block", file->path.c_str());
    file->blockSize = 0;
}

inline bool bcrLoadPartialBlock(BCRPartialFile file) {
    if (file->processed >= file->decodedLength)
        return false;
    unsigned char header[4];
    if (fread(header, 1, sizeof(header), file->file) != sizeof(header))
        bcrIoFail("truncated packed partial-BWT block header", file->path.c_str());
    const uint32_t decodedSize = bcrLoadU32(header);
    if (decodedSize == 0 || decodedSize > BCR_PARTIAL_BLOCK_BYTES ||
        decodedSize > file->decodedLength - file->processed)
        bcrIoFail("invalid packed partial-BWT block header", file->path.c_str());
    const size_t packedSize = (static_cast<size_t>(decodedSize) * 3 + 7) / 8;
    file->packed.resize(packedSize);
    if (fread(&file->packed[0], 1, packedSize, file->file) != packedSize)
        bcrIoFail("truncated packed partial-BWT block", file->path.c_str());
    file->block.resize(decodedSize);
    BCRPackedTables &tables = bcrPackedTables();
    size_t input = 0;
    size_t output = 0;
    while (decodedSize - output >= 8) {
        const uint32_t bits = static_cast<uint32_t>(file->packed[input]) |
            (static_cast<uint32_t>(file->packed[input + 1]) << 8) |
            (static_cast<uint32_t>(file->packed[input + 2]) << 16);
        const unsigned int low = bits & 0xfffU;
        const unsigned int high = (bits >> 12) & 0xfffU;
        if (!tables.decodeValid[low] || !tables.decodeValid[high])
            bcrIoFail("packed partial BWT contains an invalid code", file->path.c_str());
        const uint32_t word0 = tables.decodeFour[low];
        const uint32_t word1 = tables.decodeFour[high];
        file->block[output] = static_cast<unsigned char>(word0 & 0xffU);
        file->block[output + 1] = static_cast<unsigned char>((word0 >> 8) & 0xffU);
        file->block[output + 2] = static_cast<unsigned char>((word0 >> 16) & 0xffU);
        file->block[output + 3] = static_cast<unsigned char>((word0 >> 24) & 0xffU);
        file->block[output + 4] = static_cast<unsigned char>(word1 & 0xffU);
        file->block[output + 5] = static_cast<unsigned char>((word1 >> 8) & 0xffU);
        file->block[output + 6] = static_cast<unsigned char>((word1 >> 16) & 0xffU);
        file->block[output + 7] = static_cast<unsigned char>((word1 >> 24) & 0xffU);
        input += 3;
        output += 8;
    }
    uint32_t tail = 0;
    unsigned int tailBits = 0;
    unsigned int shift = 0;
    while (input < packedSize) {
        tail |= static_cast<uint32_t>(file->packed[input++]) << shift;
        shift += 8;
        tailBits += 8;
    }
    while (output < decodedSize) {
        const unsigned int code = tail & 0x7U;
        if (code >= 6)
            bcrIoFail("packed partial BWT contains an invalid code", file->path.c_str());
        static const unsigned char symbols[6] = {'$', 'A', 'C', 'G', 'T', 'N'};
        file->block[output++] = symbols[code];
        tail >>= 3;
        tailBits -= 3;
    }
    if (tail != 0)
        bcrIoFail("packed partial BWT has non-zero padding", file->path.c_str());
    file->blockPosition = 0;
    file->blockSize = decodedSize;
    return true;
}
#endif

inline BCRPartialFile bcrOpenPartial(const char *path, const char *mode) {
#if BCR_COMPRESSED_PARTIALS
    if (mode[0] == 'a')
        bcrIoFail("packed partial-BWT append mode is not supported", path);
    BCRPartialFile file = new BCRPackedPartial;
    file->writing = mode[0] == 'w';
    file->decodedLength = 0;
    file->processed = 0;
    file->blockPosition = 0;
    file->blockSize = 0;
    file->path = path;
    file->block.resize(BCR_PARTIAL_BLOCK_BYTES);
    file->file = fopen(path, file->writing ? "w+b" : "rb");
    if (file->file == NULL) {
        delete file;
        bcrIoFail("cannot open partial-BWT file", path);
    }
    if (file->writing) {
        unsigned char header[16] = {0};
        memcpy(header, BCR_PARTIAL_MAGIC, sizeof(BCR_PARTIAL_MAGIC));
        if (fwrite(header, 1, sizeof(header), file->file) != sizeof(header))
            bcrIoFail("cannot write packed partial-BWT header", path);
    } else {
        unsigned char header[16];
        if (fread(header, 1, sizeof(header), file->file) != sizeof(header) ||
            memcmp(header, BCR_PARTIAL_MAGIC, sizeof(BCR_PARTIAL_MAGIC)) != 0)
            bcrIoFail("invalid or truncated packed partial-BWT header", path);
        file->decodedLength = bcrLoadU64(header + 8);
    }
#else
    BCRPartialFile file = fopen(path, mode);
    if (file == NULL)
        bcrIoFail("cannot open partial-BWT file", path);
#endif
    return file;
}

inline size_t bcrReadPartial(void *data, size_t size, BCRPartialFile file) {
#if BCR_COMPRESSED_PARTIALS
    if (file->writing)
        bcrIoFail("cannot read a partial BWT opened for writing", file->path.c_str());
    unsigned char *output = static_cast<unsigned char *>(data);
    size_t total = 0;
    while (total < size && file->processed < file->decodedLength) {
        if (file->blockPosition == file->blockSize && !bcrLoadPartialBlock(file))
            break;
        const size_t available = file->blockSize - file->blockPosition;
        const size_t take = std::min(available, size - total);
        memcpy(output + total, &file->block[file->blockPosition], take);
        file->blockPosition += take;
        file->processed += static_cast<uint64_t>(take);
        total += take;
    }
    return total;
#else
    return fread(data, 1, size, file);
#endif
}

inline size_t bcrWritePartial(const void *data, size_t size, BCRPartialFile file) {
#if BCR_COMPRESSED_PARTIALS
    if (!file->writing)
        bcrIoFail("cannot write a partial BWT opened for reading", file->path.c_str());
    const unsigned char *input = static_cast<const unsigned char *>(data);
    size_t total = 0;
    while (total < size) {
        const size_t available = BCR_PARTIAL_BLOCK_BYTES - file->blockSize;
        const size_t take = std::min(available, size - total);
        memcpy(&file->block[file->blockSize], input + total, take);
        file->blockSize += take;
        total += take;
        if (file->blockSize == BCR_PARTIAL_BLOCK_BYTES)
            bcrFlushPartialBlock(file);
    }
    file->decodedLength += static_cast<uint64_t>(size);
    return size;
#else
    return fwrite(data, 1, size, file);
#endif
}

inline int bcrClosePartial(BCRPartialFile file) {
#if BCR_COMPRESSED_PARTIALS
    if (file->writing) {
        bcrFlushPartialBlock(file);
        unsigned char length[8];
        bcrStoreU64(length, file->decodedLength);
        if (fseek(file->file, 8, SEEK_SET) != 0 ||
            fwrite(length, 1, sizeof(length), file->file) != sizeof(length))
            bcrIoFail("cannot finalize packed partial-BWT header", file->path.c_str());
    }
    const int result = fclose(file->file);
    delete file;
    return result;
#else
    return fclose(file);
#endif
}

inline bool bcrPartialEof(BCRPartialFile file) {
#if BCR_COMPRESSED_PARTIALS
    return file->processed >= file->decodedLength;
#else
    return feof(file) != 0;
#endif
}

inline uint64_t bcrPartialDecodedSize(BCRPartialFile file) {
#if BCR_COMPRESSED_PARTIALS
    return file->decodedLength;
#else
    unsigned char buffer[65536];
    uint64_t total = 0;
    size_t got = 0;
    while ((got = bcrReadPartial(buffer, sizeof(buffer), file)) != 0)
        total += static_cast<uint64_t>(got);
    return total;
#endif
}

#endif
