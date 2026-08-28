#!/usr/bin/env python3
"""Streaming reader for plain and self-describing compressed eBWT files."""

from __future__ import annotations

from pathlib import Path
import struct
import zlib


MAGIC = b"EBWTCMP1"
VERSION = 1
RLE = 1
PACKED = 2
HEADER = struct.Struct("<8sBBBB4xQQ")
BLOCK_HEADER = struct.Struct("<II")
ADAPTIVE_BLOCK_BYTES = 1024 * 1024


class BwtReader:
    def __init__(self, path: Path, terminator: str):
        self.path = path
        self.terminator = ord(terminator)
        self.handle = path.open("rb")
        self.encoding = "plain"
        self.size = 0
        self.expected_checksum: int | None = None
        self.checksum = 0
        self.decoded = 0
        self.adaptive_block = b""
        self.adaptive_position = 0
        self.bit_buffer = 0
        self.bit_count = 0
        self.finished = False
        self._detect()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, traceback):
        self.handle.close()

    def _detect(self) -> None:
        probe = self.handle.read(len(MAGIC))
        if probe == MAGIC:
            rest = self.handle.read(HEADER.size - len(MAGIC))
            if len(rest) != HEADER.size - len(MAGIC):
                raise ValueError("truncated compressed eBWT header")
            magic, version, encoding, terminator, flags, size, checksum = HEADER.unpack(
                probe + rest
            )
            if magic != MAGIC or version != VERSION:
                raise ValueError("unsupported compressed eBWT version")
            if encoding not in (RLE, PACKED):
                raise ValueError("unsupported compressed eBWT encoding")
            if terminator != self.terminator:
                raise ValueError("compressed eBWT terminator does not match --terminator")
            if flags != 0:
                raise ValueError("compressed eBWT header uses unsupported flags")
            if checksum > 0xFFFFFFFF:
                raise ValueError("compressed eBWT header contains an invalid checksum")
            self.encoding = "rle" if encoding == RLE else "packed"
            self.size = size
            self.expected_checksum = checksum
            return

        if probe and MAGIC.startswith(probe):
            raise ValueError("truncated compressed eBWT header")
        self.handle.seek(0, 2)
        self.size = self.handle.tell()
        self.handle.seek(0)

    @staticmethod
    def _read_varint(payload: bytes, position: int) -> tuple[int, int]:
        value = 0
        for index in range(10):
            if position >= len(payload):
                raise ValueError("truncated RLE run length")
            byte = payload[position]
            position += 1
            if index == 9 and byte & 0xFE:
                raise ValueError("RLE run length overflows 64 bits")
            value |= (byte & 0x7F) << (7 * index)
            if not byte & 0x80:
                if value == 0:
                    raise ValueError("RLE run length must be positive")
                return value, position
        raise ValueError("invalid RLE run length")

    def _load_adaptive_block(self) -> None:
        raw_header = self.handle.read(BLOCK_HEADER.size)
        if len(raw_header) != BLOCK_HEADER.size:
            raise ValueError("truncated adaptive eBWT block header")
        decoded_size, stored_mode = BLOCK_HEADER.unpack(raw_header)
        packed = bool(stored_mode & 0x80000000)
        stored_size = stored_mode & 0x7FFFFFFF
        if not 0 < decoded_size <= ADAPTIVE_BLOCK_BYTES:
            raise ValueError("adaptive eBWT block has an invalid decoded size")
        if decoded_size > self.size - self.decoded:
            raise ValueError("adaptive eBWT block exceeds advertised decoded length")
        if stored_size == 0:
            raise ValueError("adaptive eBWT block has an empty payload")
        if packed:
            expected = (decoded_size * 3 + 7) // 8
            if stored_size != expected:
                raise ValueError("packed adaptive eBWT block has an invalid stored size")
        elif stored_size < 2 or stored_size > decoded_size * 2:
            raise ValueError("RLE adaptive eBWT block has an invalid stored size")
        payload = self.handle.read(stored_size)
        if len(payload) != stored_size:
            raise ValueError("truncated adaptive eBWT block payload")

        if packed:
            decoded = (self.terminator, ord("A"), ord("C"), ord("G"), ord("T"))
            output = bytearray()
            bits = 0
            bit_count = 0
            position = 0
            while len(output) < decoded_size:
                while bit_count < 3:
                    if position >= len(payload):
                        raise ValueError("truncated packed adaptive eBWT block")
                    bits |= payload[position] << bit_count
                    position += 1
                    bit_count += 8
                code = bits & 0x7
                bits >>= 3
                bit_count -= 3
                if code >= len(decoded):
                    raise ValueError("packed adaptive eBWT block contains an invalid symbol")
                output.append(decoded[code])
            if position != len(payload):
                raise ValueError("packed adaptive eBWT block contains trailing bytes")
            if bits:
                raise ValueError("packed adaptive eBWT block has non-zero padding bits")
            self.adaptive_block = bytes(output)
        else:
            output = bytearray()
            position = 0
            valid = {self.terminator, ord("A"), ord("C"), ord("G"), ord("T")}
            while len(output) < decoded_size:
                if position >= len(payload):
                    raise ValueError("truncated RLE adaptive eBWT block")
                symbol = payload[position]
                position += 1
                if symbol not in valid:
                    raise ValueError("RLE adaptive eBWT block contains an invalid symbol")
                run_length, position = self._read_varint(payload, position)
                if run_length > decoded_size - len(output):
                    raise ValueError("RLE run exceeds its adaptive eBWT block")
                output.extend(bytes((symbol,)) * run_length)
            if position != len(payload):
                raise ValueError("RLE adaptive eBWT block contains trailing bytes")
            self.adaptive_block = bytes(output)
        self.adaptive_position = 0

    def _read_rle(self, wanted: int) -> bytes:
        output = bytearray()
        while len(output) < wanted:
            if self.adaptive_position == len(self.adaptive_block):
                self._load_adaptive_block()
            available = len(self.adaptive_block) - self.adaptive_position
            take = min(wanted - len(output), available)
            output.extend(
                self.adaptive_block[self.adaptive_position : self.adaptive_position + take]
            )
            self.adaptive_position += take
        return bytes(output)

    def _read_packed(self, wanted: int) -> bytes:
        output = bytearray()
        decoded = (self.terminator, ord("A"), ord("C"), ord("G"), ord("T"))
        while len(output) < wanted:
            while self.bit_count < 3:
                raw = self.handle.read(1)
                if not raw:
                    raise ValueError("truncated packed eBWT payload")
                self.bit_buffer |= raw[0] << self.bit_count
                self.bit_count += 8
            code = self.bit_buffer & 0x7
            self.bit_buffer >>= 3
            self.bit_count -= 3
            if code >= len(decoded):
                raise ValueError("packed eBWT contains an invalid 3-bit symbol")
            output.append(decoded[code])
        return bytes(output)

    def read(self, wanted: int = 1024 * 1024) -> bytes:
        if wanted < 0:
            wanted = self.size - self.decoded
        wanted = min(wanted, self.size - self.decoded)
        if wanted == 0:
            return b""
        if self.encoding == "plain":
            chunk = self.handle.read(wanted)
            if len(chunk) != wanted:
                raise ValueError("plain eBWT changed size while being read")
        elif self.encoding == "rle":
            chunk = self._read_rle(wanted)
        else:
            chunk = self._read_packed(wanted)
        self.checksum = zlib.crc32(chunk, self.checksum)
        self.decoded += len(chunk)
        return chunk

    def finish(self) -> None:
        if self.finished:
            return
        if self.decoded != self.size:
            raise ValueError("eBWT decoder did not consume advertised length")
        if self.encoding == "rle" and self.adaptive_position != len(self.adaptive_block):
            raise ValueError("adaptive eBWT block exceeds advertised decoded length")
        if self.encoding == "packed" and self.bit_buffer:
            raise ValueError("packed eBWT has non-zero padding bits")
        if self.handle.read(1):
            raise ValueError("eBWT contains trailing payload bytes")
        if self.expected_checksum is not None and self.checksum != self.expected_checksum:
            raise ValueError("compressed eBWT checksum mismatch")
        self.finished = True

    def chunks(self, chunk_size: int = 1024 * 1024):
        while self.decoded < self.size:
            yield self.read(chunk_size)
        self.finish()


def read_all(path: Path, terminator: str, max_bytes: int) -> bytes:
    with BwtReader(path, terminator) as reader:
        if reader.size > max_bytes:
            raise ValueError(
                f"refusing in-memory inversion of {reader.size} decoded bytes; limit is {max_bytes}"
            )
        return b"".join(reader.chunks())
