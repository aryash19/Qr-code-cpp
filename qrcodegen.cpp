#include <algorithm>
#include <cassert>
#include <climits>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <sstream>
#include <utility>
#include "qrcodegen.hpp"

// The following macro and include are for generating PNG images.
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

using std::int8_t;
using std::uint8_t;
using std::size_t;
using std::vector;

namespace qrcodegen {

    /*---- Class QrSegment ----*/
    // QrSegment represents a part (or segment) of the QR code data.
    // This class encapsulates data (such as mode, character count, and bits)
    // along with functions to create segments from various data types.
    class QrSegment {
    public:
        // Nested class to hold encoding mode information.
        class Mode {
        public:
            // Constructor: Sets mode bits and how many bits are used to store character counts.
            Mode(int mode, int cc0, int cc1, int cc2) : modeBits(mode) {
                numBitsCharCount[0] = cc0;
                numBitsCharCount[1] = cc1;
                numBitsCharCount[2] = cc2;
            }
            // Return the mode bits.
            int getModeBits() const { return modeBits; }
            // Return the number of bits used for the character count field based on version.
            int numCharCountBits(int ver) const { return numBitsCharCount[(ver + 7) / 17]; }
        private:
            int modeBits;           // Mode bits that define the segment type (numeric, alphanumeric, etc.)
            int numBitsCharCount[3]; // Array storing bit lengths for different version ranges.
        };

        // Factory methods (static functions) to create segments from different types of data.
        // They hide the details of converting input text/data into the internal bit representation.
        static QrSegment makeBytes(const vector<uint8_t> &data);
        static QrSegment makeNumeric(const char *digits);
        static QrSegment makeAlphanumeric(const char *text);
        static vector<QrSegment> makeSegments(const char *text);
        static QrSegment makeEci(long assignVal);

        // Constructor overloads that use encapsulation: the mode, number of characters, and data bits are bundled.
        QrSegment(const Mode &md, int numCh, const std::vector<bool> &dt);
        QrSegment(const Mode &md, int numCh, std::vector<bool> &&dt);

        // Accessor functions allow safe retrieval of internal data.
        const Mode &getMode() const;
        int getNumChars() const;
        const vector<bool> &getData() const;

        // Computes total bits required for an array of segments for a given QR version.
        static int getTotalBits(const vector<QrSegment> &segs, int version);

        // Utility methods to check the type of string content.
        static bool isNumeric(const char *text);
        static bool isAlphanumeric(const char *text);
    private:
        const Mode *mode;       // Points to a Mode object, encapsulating the encoding type.
        int numChars;           // Stores the number of characters in the segment.
        vector<bool> data;      // The actual data bits representing the segment.

        // Character set for alphanumeric mode.
        static const char *ALPHANUMERIC_CHARSET;
    };

    // Implementation of QrSegment factory methods:
    QrSegment QrSegment::makeBytes(const vector<uint8_t> &data) {
        if (data.size() > static_cast<unsigned int>(INT_MAX))
            throw std::length_error("Data too long");
        BitBuffer bb;
        // Append each byte into the bit buffer as 8 bits.
        for (uint8_t b : data)
            bb.appendBits(b, 8);
        return QrSegment(Mode::BYTE, static_cast<int>(data.size()), std::move(bb));
    }

    QrSegment QrSegment::makeNumeric(const char *digits) {
        BitBuffer bb;
        int accumData = 0;
        int accumCount = 0;
        int charCount = 0;
        // Convert numeric characters into bits, grouping by 3 digits.
        for (; *digits != '\0'; digits++, charCount++) {
            char c = *digits;
            if (c < '0' || c > '9')
                throw std::domain_error("String contains non-numeric characters");
            accumData = accumData * 10 + (c - '0');
            accumCount++;
            if (accumCount == 3) {
                bb.appendBits(static_cast<uint32_t>(accumData), 10);
                accumData = 0;
                accumCount = 0;
            }
        }
        if (accumCount > 0)  // Handle remaining 1 or 2 digits.
            bb.appendBits(static_cast<uint32_t>(accumData), accumCount * 3 + 1);
        return QrSegment(Mode::NUMERIC, charCount, std::move(bb));
    }

    QrSegment QrSegment::makeAlphanumeric(const char *text) {
        BitBuffer bb;
        int accumData = 0;
        int accumCount = 0;
        int charCount = 0;
        // Convert alphanumeric characters into bits.
        for (; *text != '\0'; text++, charCount++) {
            const char *temp = std::strchr(ALPHANUMERIC_CHARSET, *text);
            if (temp == nullptr)
                throw std::domain_error("String contains unencodable characters in alphanumeric mode");
            accumData = accumData * 45 + static_cast<int>(temp - ALPHANUMERIC_CHARSET);
            accumCount++;
            if (accumCount == 2) {
                bb.appendBits(static_cast<uint32_t>(accumData), 11);
                accumData = 0;
                accumCount = 0;
            }
        }
        if (accumCount > 0)  // If one character remains, use 6 bits.
            bb.appendBits(static_cast<uint32_t>(accumData), 6);
        return QrSegment(Mode::ALPHANUMERIC, charCount, std::move(bb));
    }

    vector<QrSegment> QrSegment::makeSegments(const char *text) {
        // This method chooses the best encoding for the text.
        vector<QrSegment> result;
        if (*text == '\0');  // Empty text produces no segments.
        else if (isNumeric(text))
            result.push_back(makeNumeric(text));
        else if (isAlphanumeric(text))
            result.push_back(makeAlphanumeric(text));
        else {
            // For any other text, use byte mode.
            vector<uint8_t> bytes;
            for (; *text != '\0'; text++)
                bytes.push_back(static_cast<uint8_t>(*text));
            result.push_back(makeBytes(bytes));
        }
        return result;
    }

    QrSegment QrSegment::makeEci(long assignVal) {
        BitBuffer bb;
        if (assignVal < 0)
            throw std::domain_error("ECI assignment value out of range");
        else if (assignVal < (1 << 7))
            bb.appendBits(static_cast<uint32_t>(assignVal), 8);
        else if (assignVal < (1 << 14)) {
            bb.appendBits(2, 2);
            bb.appendBits(static_cast<uint32_t>(assignVal), 14);
        } else if (assignVal < 1000000L) {
            bb.appendBits(6, 3);
            bb.appendBits(static_cast<uint32_t>(assignVal), 21);
        } else
            throw std::domain_error("ECI assignment value out of range");
        return QrSegment(Mode::ECI, 0, std::move(bb));
    }

    // Constructors for QrSegment use encapsulation to initialize their data.
    QrSegment::QrSegment(const Mode &md, int numCh, const std::vector<bool> &dt)
        : mode(&md), numChars(numCh), data(dt) {
        if (numCh < 0)
            throw std::domain_error("Invalid value");
    }

    QrSegment::QrSegment(const Mode &md, int numCh, std::vector<bool> &&dt)
        : mode(&md), numChars(numCh), data(std::move(dt)) {
        if (numCh < 0)
            throw std::domain_error("Invalid value");
    }

    // Accessor methods providing read-only access to private members.
    const QrSegment::Mode &QrSegment::getMode() const { return *mode; }
    int QrSegment::getNumChars() const { return numChars; }
    const vector<bool> &QrSegment::getData() const { return data; }

    // Helper methods to check if a string contains only numeric or only allowed alphanumeric characters.
    bool QrSegment::isNumeric(const char *text) {
        for (; *text != '\0'; text++) {
            char c = *text;
            if (c < '0' || c > '9')
                return false;
        }
        return true;
    }

    bool QrSegment::isAlphanumeric(const char *text) {
        for (; *text != '\0'; text++) {
            if (std::strchr(ALPHANUMERIC_CHARSET, *text) == nullptr)
                return false;
        }
        return true;
    }

    // Computes the total number of bits needed to encode all segments for a specific QR version.
    int QrSegment::getTotalBits(const vector<QrSegment> &segs, int version) {
        int result = 0;
        for (const QrSegment &seg : segs) {
            int ccbits = seg.mode->numCharCountBits(version);
            if (seg.numChars >= (1L << ccbits))
                return -1;  // The segment's length doesn't fit within the allowed bit width.
            if (4 + ccbits > INT_MAX - result)
                return -1;  // Prevent integer overflow.
            result += 4 + ccbits;
            if (seg.data.size() > static_cast<unsigned int>(INT_MAX - result))
                return -1;
            result += static_cast<int>(seg.data.size());
        }
        return result;
    }

    // Define static constants for modes.
    const QrSegment::Mode QrSegment::Mode::NUMERIC     (0x1, 10, 12, 14);
    const QrSegment::Mode QrSegment::Mode::ALPHANUMERIC(0x2,  9, 11, 13);
    const QrSegment::Mode QrSegment::Mode::BYTE        (0x4,  8, 16, 16);
    const QrSegment::Mode QrSegment::Mode::KANJI       (0x8,  8, 10, 12);
    const QrSegment::Mode QrSegment::Mode::ECI         (0x7,  0,  0,  0);

    const char *QrSegment::ALPHANUMERIC_CHARSET = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ $%*+-./:";

    /*---- Class QrCode ----*/
    // QrCode is the main class that builds and represents a complete QR Code.
    // It encapsulates properties like version, size, error correction level, and the grid of modules.
    class QrCode {
    public:
        // Static methods that provide a simple interface to create QR Codes from text or binary data.
        static QrCode encodeText(const char *text, Ecc ecl);
        static QrCode encodeBinary(const vector<uint8_t> &data, Ecc ecl);
        static QrCode encodeSegments(const vector<QrSegment> &segs, Ecc ecl,
                                     int minVersion = MIN_VERSION, int maxVersion = MAX_VERSION,
                                     int mask = -1, bool boostEcl = true);

        // Getter methods to access QR Code properties.
        int getVersion() const;
        int getSize() const;
        Ecc getErrorCorrectionLevel() const;
        int getMask() const;
        bool getModule(int x, int y) const;

        // Method to output the QR Code as a PNG image.
        void toPng(const char *filename, int scale) const;
    private:
        int version;                    // QR Code version (affects size and data capacity)
        Ecc errorCorrectionLevel;       // Error correction level (robustness against errors)
        int size;                       // Dimension of the QR Code grid (derived from version)
        int mask;                       // Mask pattern number used to improve readability
        vector<vector<bool>> modules;   // 2D grid representing the QR Code (true = dark, false = light)
        vector<vector<bool>> isFunction;  // Marks cells that are reserved for function patterns

        // Private constructor to force the use of static encoding methods.
        QrCode(int ver, Ecc ecl, const vector<uint8_t> &dataCodewords, int msk);

        // Private methods that perform various steps of the QR Code generation.
        void drawFunctionPatterns();
        void drawFormatBits(int msk);
        void drawVersion();
        void drawFinderPattern(int x, int y);
        void drawAlignmentPattern(int x, int y);
        void setFunctionModule(int x, int y, bool isDark);
        bool module(int x, int y) const;
        vector<uint8_t> addEccAndInterleave(const vector<uint8_t> &data) const;
        void drawCodewords(const vector<uint8_t> &data);
        void applyMask(int msk);
        long getPenaltyScore() const;
        vector<int> getAlignmentPatternPositions() const;
        int getNumRawDataModules(int ver);
        int getNumDataCodewords(int ver, Ecc ecl);

        // Reed-Solomon error correction helper methods.
        vector<uint8_t> reedSolomonComputeDivisor(int degree);
        vector<uint8_t> reedSolomonComputeRemainder(const vector<uint8_t> &data,
                                                     const vector<uint8_t> &divisor);
        uint8_t reedSolomonMultiply(uint8_t x, uint8_t y);

        // Helper methods for penalty score calculations used in mask evaluation.
        int finderPenaltyCountPatterns(const std::array<int,7> &runHistory) const;
        int finderPenaltyTerminateAndCount(bool currentRunColor, int currentRunLength,
                                           std::array<int,7> &runHistory) const;
        void finderPenaltyAddHistory(int currentRunLength, std::array<int,7> &runHistory) const;
        static bool getBit(long x, int i);
    };

    // Implementation of static encoding methods.
    QrCode QrCode::encodeText(const char *text, Ecc ecl) {
        vector<QrSegment> segs = QrSegment::makeSegments(text);
        return encodeSegments(segs, ecl);
    }

    QrCode QrCode::encodeBinary(const vector<uint8_t> &data, Ecc ecl) {
        vector<QrSegment> segs{QrSegment::makeBytes(data)};
        return encodeSegments(segs, ecl);
    }

    QrCode QrCode::encodeSegments(const vector<QrSegment> &segs, Ecc ecl,
                                  int minVersion, int maxVersion, int mask, bool boostEcl) {
        // Check that version numbers and mask values are valid.
        if (!(MIN_VERSION <= minVersion && minVersion <= maxVersion && maxVersion <= MAX_VERSION) || mask < -1 || mask > 7)
            throw std::invalid_argument("Invalid value");

        // Determine the minimum version that can hold the data.
        int version, dataUsedBits;
        for (version = minVersion; ; version++) {
            int dataCapacityBits = getNumDataCodewords(version, ecl) * 8;  // available data bits
            dataUsedBits = QrSegment::getTotalBits(segs, version);
            if (dataUsedBits != -1 && dataUsedBits <= dataCapacityBits)
                break;  // Found a version that fits
            if (version >= maxVersion) {
                std::ostringstream sb;
                if (dataUsedBits == -1)
                    sb << "Segment too long";
                else {
                    sb << "Data length = " << dataUsedBits << " bits, ";
                    sb << "Max capacity = " << dataCapacityBits << " bits";
                }
                throw data_too_long(sb.str());
            }
        }
        assert(dataUsedBits != -1);

        // Optionally boost the error correction level if the data still fits.
        for (Ecc newEcl : {Ecc::MEDIUM, Ecc::QUARTILE, Ecc::HIGH}) {
            if (boostEcl && dataUsedBits <= getNumDataCodewords(version, newEcl) * 8)
                ecl = newEcl;
        }

        // Concatenate segment data into one BitBuffer.
        BitBuffer bb;
        for (const QrSegment &seg : segs) {
            bb.appendBits(static_cast<uint32_t>(seg.getMode().getModeBits()), 4);
            bb.appendBits(static_cast<uint32_t>(seg.getNumChars()), seg.getMode().numCharCountBits(version));
            bb.insert(bb.end(), seg.getData().begin(), seg.getData().end());
        }
        assert(bb.size() == static_cast<unsigned int>(dataUsedBits));

        // Append terminator and pad bits to form complete bytes.
        size_t dataCapacityBits = static_cast<size_t>(getNumDataCodewords(version, ecl)) * 8;
        assert(bb.size() <= dataCapacityBits);
        bb.appendBits(0, std::min(4, static_cast<int>(dataCapacityBits - bb.size())));
        bb.appendBits(0, (8 - static_cast<int>(bb.size() % 8)) % 8);
        assert(bb.size() % 8 == 0);

        // Add padding bytes until reaching the required capacity.
        for (uint8_t padByte = 0xEC; bb.size() < dataCapacityBits; padByte ^= 0xEC ^ 0x11)
            bb.appendBits(padByte, 8);

        // Pack bits into bytes (big endian).
        vector<uint8_t> dataCodewords(bb.size() / 8);
        for (size_t i = 0; i < bb.size(); i++)
            dataCodewords.at(i >> 3) |= (bb.at(i) ? 1 : 0) << (7 - (i & 7));

        // Create and return the QR Code object.
        return QrCode(version, ecl, dataCodewords, mask);
    }

    // Private constructor initializes the QR Code properties and builds the module grid.
    QrCode::QrCode(int ver, Ecc ecl, const vector<uint8_t> &dataCodewords, int msk)
        : version(ver), errorCorrectionLevel(ecl) {
        if (ver < MIN_VERSION || ver > MAX_VERSION)
            throw std::domain_error("Version value out of range");
        if (msk < -1 || msk > 7)
            throw std::domain_error("Mask value out of range");
        size = ver * 4 + 17;
        size_t sz = static_cast<size_t>(size);
        // Initialize the module grid (all light initially) and the function pattern markers.
        modules    = vector<vector<bool>>(sz, vector<bool>(sz));
        isFunction = vector<vector<bool>>(sz, vector<bool>(sz));

        // Draw fixed patterns (finder, timing, alignment, etc.).
        drawFunctionPatterns();
        const vector<uint8_t> allCodewords = addEccAndInterleave(dataCodewords);
        drawCodewords(allCodewords);

        // Determine and apply the best mask if not provided.
        if (msk == -1) {
            long minPenalty = LONG_MAX;
            for (int i = 0; i < 8; i++) {
                applyMask(i);
                drawFormatBits(i);
                long penalty = getPenaltyScore();
                if (penalty < minPenalty) {
                    msk = i;
                    minPenalty = penalty;
                }
                applyMask(i);  // Undo the mask (XOR operation is its own inverse)
            }
        }
        assert(0 <= msk && msk <= 7);
        mask = msk;
        applyMask(msk);        // Apply chosen mask.
        drawFormatBits(msk);   // Draw format information bits.

        // Free up memory from isFunction since it's no longer needed.
        isFunction.clear();
        isFunction.shrink_to_fit();
    }

    // Getter functions for various QR Code properties.
    int QrCode::getVersion() const { return version; }
    int QrCode::getSize() const { return size; }
    QrCode::Ecc QrCode::getErrorCorrectionLevel() const { return errorCorrectionLevel; }
    int QrCode::getMask() const { return mask; }

    // Returns whether the module (cell) at (x, y) is dark.
    bool QrCode::getModule(int x, int y) const {
        return 0 <= x && x < size && 0 <= y && y < size && module(x, y);
    }

    // Draws function patterns such as timing patterns, finder patterns, and alignment patterns.
    void QrCode::drawFunctionPatterns() {
        // Draw horizontal and vertical timing patterns.
        for (int i = 0; i < size; i++) {
            setFunctionModule(6, i, i % 2 == 0);
            setFunctionModule(i, 6, i % 2 == 0);
        }
        // Draw finder patterns at three corners.
        drawFinderPattern(3, 3);
        drawFinderPattern(size - 4, 3);
        drawFinderPattern(3, size - 4);

        // Draw alignment patterns based on version.
        const vector<int> alignPatPos = getAlignmentPatternPositions();
        size_t numAlign = alignPatPos.size();
        for (size_t i = 0; i < numAlign; i++) {
            for (size_t j = 0; j < numAlign; j++) {
                // Skip finder corners.
                if (!((i == 0 && j == 0) || (i == 0 && j == numAlign - 1) || (i == numAlign - 1 && j == 0)))
                    drawAlignmentPattern(alignPatPos.at(i), alignPatPos.at(j));
            }
        }

        // Draw placeholders for format and version bits.
        drawFormatBits(0);  // Dummy mask; will be updated later.
        drawVersion();
    }

    // Draws format bits which include error correction and mask information.
    void QrCode::drawFormatBits(int msk) {
        // Calculate combined format bits and then an error correction code for them.
        int data = getFormatBits(errorCorrectionLevel) << 3 | msk;
        int rem = data;
        for (int i = 0; i < 10; i++)
            rem = (rem << 1) ^ ((rem >> 9) * 0x537);
        int bits = (data << 10 | rem) ^ 0x5412;  // 15-bit result.
        assert(bits >> 15 == 0);

        // Draw the first copy of format bits near the top-left.
        for (int i = 0; i <= 5; i++)
            setFunctionModule(8, i, getBit(bits, i));
        setFunctionModule(8, 7, getBit(bits, 6));
        setFunctionModule(8, 8, getBit(bits, 7));
        setFunctionModule(7, 8, getBit(bits, 8));
        for (int i = 9; i < 15; i++)
            setFunctionModule(14 - i, 8, getBit(bits, i));

        // Draw the second copy in a mirrored position.
        for (int i = 0; i < 8; i++)
            setFunctionModule(size - 1 - i, 8, getBit(bits, i));
        for (int i = 8; i < 15; i++)
            setFunctionModule(8, size - 15 + i, getBit(bits, i));
        setFunctionModule(8, size - 8, true);  // Always dark.
    }

    // Draws version information (for version 7 or higher).
    void QrCode::drawVersion() {
        if (version < 7)
            return;
        int rem = version;
        for (int i = 0; i < 12; i++)
            rem = (rem << 1) ^ ((rem >> 11) * 0x1F25);
        long bits = static_cast<long>(version) << 12 | rem;
        assert(bits >> 18 == 0);

        // Draw two copies of the version information.
        for (int i = 0; i < 18; i++) {
            bool bit = getBit(bits, i);
            int a = size - 11 + i % 3;
            int b = i / 3;
            setFunctionModule(a, b, bit);
            setFunctionModule(b, a, bit);
        }
    }

    // Draws the finder pattern (the large square patterns in three corners).
    void QrCode::drawFinderPattern(int x, int y) {
        for (int dy = -4; dy <= 4; dy++) {
            for (int dx = -4; dx <= 4; dx++) {
                int dist = std::max(std::abs(dx), std::abs(dy));  // Use Chebyshev distance.
                int xx = x + dx, yy = y + dy;
                if (0 <= xx && xx < size && 0 <= yy && yy < size)
                    // Dark except in a ring around the center.
                    setFunctionModule(xx, yy, dist != 2 && dist != 4);
            }
        }
    }

    // Draws a smaller alignment pattern used in larger QR Codes.
    void QrCode::drawAlignmentPattern(int x, int y) {
        for (int dy = -2; dy <= 2; dy++) {
            for (int dx = -2; dx <= 2; dx++)
                setFunctionModule(x + dx, y + dy, std::max(std::abs(dx), std::abs(dy)) != 1);
        }
    }

    // Sets a module (cell) in the grid as part of the function patterns.
    void QrCode::setFunctionModule(int x, int y, bool isDark) {
        size_t ux = static_cast<size_t>(x);
        size_t uy = static_cast<size_t>(y);
        modules.at(uy).at(ux) = isDark;
        isFunction.at(uy).at(ux) = true;
    }

    // Returns the color (dark/light) of a module.
    bool QrCode::module(int x, int y) const {
        return modules.at(static_cast<size_t>(y)).at(static_cast<size_t>(x));
    }

    // Adds error correction codewords and interleaves the data as per QR Code standards.
    vector<uint8_t> QrCode::addEccAndInterleave(const vector<uint8_t> &data) const {
        if (data.size() != static_cast<unsigned int>(getNumDataCodewords(version, errorCorrectionLevel)))
            throw std::invalid_argument("Invalid argument");

        // Compute block parameters.
        int numBlocks = NUM_ERROR_CORRECTION_BLOCKS[static_cast<int>(errorCorrectionLevel)][version];
        int blockEccLen = ECC_CODEWORDS_PER_BLOCK[static_cast<int>(errorCorrectionLevel)][version];
        int rawCodewords = getNumRawDataModules(version) / 8;
        int numShortBlocks = numBlocks - rawCodewords % numBlocks;
        int shortBlockLen = rawCodewords / numBlocks;

        // Split data into blocks and compute ECC for each block.
        vector<vector<uint8_t>> blocks;
        const vector<uint8_t> rsDiv = reedSolomonComputeDivisor(blockEccLen);
        for (int i = 0, k = 0; i < numBlocks; i++) {
            vector<uint8_t> dat(data.cbegin() + k, data.cbegin() + (k + shortBlockLen - blockEccLen + (i < numShortBlocks ? 0 : 1)));
            k += static_cast<int>(dat.size());
            const vector<uint8_t> ecc = reedSolomonComputeRemainder(dat, rsDiv);
            if (i < numShortBlocks)
                dat.push_back(0);
            dat.insert(dat.end(), ecc.cbegin(), ecc.cend());
            blocks.push_back(std::move(dat));
        }

        // Interleave bytes from the blocks.
        vector<uint8_t> result;
        for (size_t i = 0; i < blocks.at(0).size(); i++) {
            for (size_t j = 0; j < blocks.size(); j++) {
                if (i != static_cast<unsigned int>(shortBlockLen - blockEccLen) || j >= static_cast<unsigned int>(numShortBlocks))
                    result.push_back(blocks.at(j).at(i));
            }
        }
        assert(result.size() == static_cast<unsigned int>(rawCodewords));
        return result;
    }

    // Draws the encoded codewords into the QR Code grid using a zigzag pattern.
    void QrCode::drawCodewords(const vector<uint8_t> &data) {
        if (data.size() != static_cast<unsigned int>(getNumRawDataModules(version) / 8))
            throw std::invalid_argument("Invalid argument");
        size_t i = 0;  // Bit index
        for (int right = size - 1; right >= 1; right -= 2) {
            if (right == 6)
                right = 5;
            for (int vert = 0; vert < size; vert++) {
                for (int j = 0; j < 2; j++) {
                    size_t x = static_cast<size_t>(right - j);
                    bool upward = ((right + 1) & 2) == 0;
                    size_t y = static_cast<size_t>(upward ? size - 1 - vert : vert);
                    if (!isFunction.at(y).at(x) && i < data.size() * 8) {
                        modules.at(y).at(x) = getBit(data.at(i >> 3), 7 - static_cast<int>(i & 7));
                        i++;
                    }
                }
            }
        }
        assert(i == data.size() * 8);
    }

    // Applies a mask pattern to the QR Code (to reduce problematic patterns).
    void QrCode::applyMask(int msk) {
        if (msk < 0 || msk > 7)
            throw std::domain_error("Mask value out of range");
        size_t sz = static_cast<size_t>(size);
        for (size_t y = 0; y < sz; y++) {
            for (size_t x = 0; x < sz; x++) {
                bool invert;
                switch (msk) {
                    case 0:  invert = (x + y) % 2 == 0; break;
                    case 1:  invert = y % 2 == 0; break;
                    case 2:  invert = x % 3 == 0; break;
                    case 3:  invert = (x + y) % 3 == 0; break;
                    case 4:  invert = (x / 3 + y / 2) % 2 == 0; break;
                    case 5:  invert = x * y % 2 + x * y % 3 == 0; break;
                    case 6:  invert = (x * y % 2 + x * y % 3) % 2 == 0; break;
                    case 7:  invert = ((x + y) % 2 + x * y % 3) % 2 == 0; break;
                    default: throw std::logic_error("Unreachable");
                }
                modules.at(y).at(x) = modules.at(y).at(x) ^ (invert & !isFunction.at(y).at(x));
            }
        }
    }

    // Computes the penalty score for the QR Code, used to determine the best mask.
    long QrCode::getPenaltyScore() const {
        long result = 0;
        // Check rows for adjacent modules with the same color.
        for (int y = 0; y < size; y++) {
            bool runColor = false;
            int runX = 0;
            std::array<int,7> runHistory = {};
            for (int x = 0; x < size; x++) {
                if (module(x, y) == runColor) {
                    runX++;
                    if (runX == 5)
                        result += PENALTY_N1;
                    else if (runX > 5)
                        result++;
                } else {
                    finderPenaltyAddHistory(runX, runHistory);
                    if (!runColor)
                        result += finderPenaltyCountPatterns(runHistory) * PENALTY_N3;
                    runColor = module(x, y);
                    runX = 1;
                }
            }
            result += finderPenaltyTerminateAndCount(runColor, runX, runHistory) * PENALTY_N3;
        }
        // Check columns similarly.
        for (int x = 0; x < size; x++) {
            bool runColor = false;
            int runY = 0;
            std::array<int,7> runHistory = {};
            for (int y = 0; y < size; y++) {
                if (module(x, y) == runColor) {
                    runY++;
                    if (runY == 5)
                        result += PENALTY_N1;
                    else if (runY > 5)
                        result++;
                } else {
                    finderPenaltyAddHistory(runY, runHistory);
                    if (!runColor)
                        result += finderPenaltyCountPatterns(runHistory) * PENALTY_N3;
                    runColor = module(x, y);
                    runY = 1;
                }
            }
            result += finderPenaltyTerminateAndCount(runColor, runY, runHistory) * PENALTY_N3;
        }
        // Check for 2x2 blocks of same color.
        for (int y = 0; y < size - 1; y++) {
            for (int x = 0; x < size - 1; x++) {
                bool color = module(x, y);
                if (color == module(x + 1, y) &&
                    color == module(x, y + 1) &&
                    color == module(x + 1, y + 1))
                    result += PENALTY_N2;
            }
        }
        // Balance dark and light modules.
        int dark = 0;
        for (const vector<bool> &row : modules) {
            for (bool color : row) {
                if (color)
                    dark++;
            }
        }
        int total = size * size;
        int k = static_cast<int>((std::abs(dark * 20L - total * 10L) + total - 1) / total) - 1;
        assert(0 <= k && k <= 9);
        result += k * PENALTY_N4;
        assert(0 <= result && result <= 2568888L);
        return result;
    }

    // Determines positions for alignment patterns.
    vector<int> QrCode::getAlignmentPatternPositions() const {
        if (version == 1)
            return vector<int>();
        else {
            int numAlign = version / 7 + 2;
            int step = (version * 8 + numAlign * 3 + 5) / (numAlign * 4 - 4) * 2;
            vector<int> result;
            for (int i = 0, pos = size - 7; i < numAlign - 1; i++, pos -= step)
                result.insert(result.begin(), pos);
            result.insert(result.begin(), 6);
            return result;
        }
    }

    // Calculates the total number of raw data modules.
    int QrCode::getNumRawDataModules(int ver) {
        if (ver < MIN_VERSION || ver > MAX_VERSION)
            throw std::domain_error("Version number out of range");
        int result = (16 * ver + 128) * ver + 64;
        if (ver >= 2) {
            int numAlign = ver / 7 + 2;
            result -= (25 * numAlign - 10) * numAlign - 55;
            if (ver >= 7)
                result -= 36;
        }
        assert(208 <= result && result <= 29648);
        return result;
    }

    // Computes how many data codewords can be stored.
    int QrCode::getNumDataCodewords(int ver, Ecc ecl) {
        return getNumRawDataModules(ver) / 8
            - ECC_CODEWORDS_PER_BLOCK[static_cast<int>(ecl)][ver]
            * NUM_ERROR_CORRECTION_BLOCKS[static_cast<int>(ecl)][ver];
    }

    // Computes the Reed-Solomon divisor polynomial.
    vector<uint8_t> QrCode::reedSolomonComputeDivisor(int degree) {
        if (degree < 1 || degree > 255)
            throw std::domain_error("Degree out of range");
        vector<uint8_t> result(static_cast<size_t>(degree));
        result.at(result.size() - 1) = 1;  // The polynomial starts with 1.
        uint8_t root = 1;
        for (int i = 0; i < degree; i++) {
            for (size_t j = 0; j < result.size(); j++) {
                result.at(j) = reedSolomonMultiply(result.at(j), root);
                if (j + 1 < result.size())
                    result.at(j) ^= result.at(j + 1);
            }
            root = reedSolomonMultiply(root, 0x02);
        }
        return result;
    }

    // Computes the remainder polynomial for Reed-Solomon error correction.
    vector<uint8_t> QrCode::reedSolomonComputeRemainder(const vector<uint8_t> &data, const vector<uint8_t> &divisor) {
        vector<uint8_t> result(divisor.size());
        for (uint8_t b : data) {
            uint8_t factor = b ^ result.at(0);
            result.erase(result.begin());
            result.push_back(0);
            for (size_t i = 0; i < result.size(); i++)
                result.at(i) ^= reedSolomonMultiply(divisor.at(i), factor);
        }
        return result;
    }

    // Multiplies two numbers in the finite field GF(2^8) used by Reed-Solomon.
    uint8_t QrCode::reedSolomonMultiply(uint8_t x, uint8_t y) {
        int z = 0;
        for (int i = 7; i >= 0; i--) {
            z = (z << 1) ^ ((z >> 7) * 0x11D);
            z ^= ((y >> i) & 1) * x;
        }
        assert(z >> 8 == 0);
        return static_cast<uint8_t>(z);
    }

    // Helper methods used in computing penalties for patterns.
    int QrCode::finderPenaltyCountPatterns(const std::array<int,7> &runHistory) const {
        int n = runHistory.at(1);
        assert(n <= size * 3);
        bool core = n > 0 && runHistory.at(2) == n && runHistory.at(3) == n * 3 && runHistory.at(4) == n && runHistory.at(5) == n;
        return (core && runHistory.at(0) >= n * 4 && runHistory.at(6) >= n ? 1 : 0)
             + (core && runHistory.at(6) >= n * 4 && runHistory.at(0) >= n ? 1 : 0);
    }

    int QrCode::finderPenaltyTerminateAndCount(bool currentRunColor, int currentRunLength, std::array<int,7> &runHistory) const {
        if (currentRunColor) {  // Terminate dark run.
            finderPenaltyAddHistory(currentRunLength, runHistory);
            currentRunLength = 0;
        }
        currentRunLength += size;  // Add light border.
        finderPenaltyAddHistory(currentRunLength, runHistory);
        return finderPenaltyCountPatterns(runHistory);
    }

    void QrCode::finderPenaltyAddHistory(int currentRunLength, std::array<int,7> &runHistory) const {
        if (runHistory.at(0) == 0)
            currentRunLength += size;
        std::copy_backward(runHistory.cbegin(), runHistory.cend() - 1, runHistory.end());
        runHistory.at(0) = currentRunLength;
    }

    // Helper to extract a bit from a number.
    bool QrCode::getBit(long x, int i) {
        return ((x >> i) & 1) != 0;
    }

    // Constants for penalty scores.
    const int QrCode::PENALTY_N1 = 3;
    const int QrCode::PENALTY_N2 = 3;
    const int QrCode::PENALTY_N3 = 40;
    const int QrCode::PENALTY_N4 = 10;

    // Tables for error correction parameters.
    const int8_t QrCode::ECC_CODEWORDS_PER_BLOCK[4][41] = {
        {-1,  7, 10, 15, 20, 26, 18, 20, 24, 30, 18, 20, 24, 26, 30, 22, 24, 28, 30, 28, 28, 28, 28, 30, 30, 26, 28, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30},
        {-1, 10, 16, 26, 18, 24, 16, 18, 22, 22, 26, 30, 22, 22, 24, 24, 28, 28, 26, 26, 26, 26, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28},
        {-1, 13, 22, 18, 26, 18, 24, 18, 22, 20, 24, 28, 26, 24, 20, 30, 24, 28, 28, 26, 30, 28, 30, 30, 30, 30, 28, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30},
        {-1, 17, 28, 22, 16, 22, 28, 26, 26, 24, 28, 24, 28, 22, 24, 24, 30, 28, 28, 26, 28, 30, 24, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30},
    };

    const int8_t QrCode::NUM_ERROR_CORRECTION_BLOCKS[4][41] = {
        {-1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 4,  4,  4,  4,  4,  6,  6,  6,  6,  7,  8,  8,  9,  9, 10, 12, 12, 12, 13, 14, 15, 16, 17, 18, 19, 19, 20, 21, 22, 24, 25},
        {-1, 1, 1, 1, 2, 2, 4, 4, 4, 5, 5,  5,  8,  9,  9, 10, 10, 11, 13, 14, 16, 17, 17, 18, 20, 21, 23, 25, 26, 28, 29, 31, 33, 35, 37, 38, 40, 43, 45, 47, 49},
        {-1, 1, 1, 2, 2, 4, 4, 6, 6, 8, 8,  8, 10, 12, 16, 12, 17, 16, 18, 21, 20, 23, 23, 25, 27, 29, 34, 34, 35, 38, 40, 43, 45, 48, 51, 53, 56, 59, 62, 65, 68},
        {-1, 1, 1, 2, 4, 4, 4, 5, 6, 8, 8, 11, 11, 16, 16, 18, 16, 19, 21, 25, 25, 25, 34, 30, 32, 35, 37, 40, 42, 45, 48, 51, 54, 57, 60, 63, 66, 70, 74, 77, 81},
    };

    // Custom exception for data too long.
    data_too_long::data_too_long(const std::string &msg) : std::length_error(msg) {}

    /*---- Class BitBuffer ----*/
    // BitBuffer extends std::vector<bool> to provide an easy way to append bits.
    // This is an example of inheritance in OOP.
    class BitBuffer : public std::vector<bool> {
    public:
        BitBuffer() : std::vector<bool>() {}
        // Appends a value as bits into the buffer.
        void appendBits(std::uint32_t val, int len) {
            if (len < 0 || len > 31 || val >> len != 0)
                throw std::domain_error("Value out of range");
            for (int i = len - 1; i >= 0; i--)  // Append each bit, most-significant bit first.
                this->push_back(((val >> i) & 1) != 0);
        }
    };

    // Method to output the QR Code as a PNG image.
    void QrCode::toPng(const char *filename, int scale) const {
        if (scale <= 0)
            throw std::invalid_argument("Scale must be positive");
        int border = 4;
        int dim = (size + border * 2) * scale;
        vector<uint8_t> image(dim * dim, 255);  // Start with a white image.

        // Draw dark modules (black pixels) based on the QR Code grid.
        for (int y = 0; y < size; y++) {
            for (int x = 0; x < size; x++) {
                if (getModule(x, y)) {
                    for (int dy = 0; dy < scale; dy++) {
                        for (int dx = 0; dx < scale; dx++) {
                            int index = ((border + y) * scale + dy) * dim + (border + x) * scale + dx;
                            image[index] = 0;
                        }
                    }
                }
            }
        }

        stbi_write_png(filename, dim, dim, 1, image.data(), dim);
    }

} // namespace qrcodegen
