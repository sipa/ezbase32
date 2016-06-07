#include <thread>
#include <mutex>

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <set>

#include <immintrin.h>

namespace {

const uint32_t tbl[60][5] = {
    /* length91 hd5 BCH codes */
    {0x053395c1, 0x09172b82, 0x115e51e4, 0x21d7ffc8, 0x2edf2570},
    {0x0549b4e1, 0x0817e9c2, 0x102f7784, 0x204ace28, 0x1281bd70},
    {0x07339121, 0x0c372242, 0x186e4024, 0x30ce9448, 0x2b8fa890},
    {0x0749ace1, 0x0d5dd9c2, 0x1975c784, 0x313b7ca8, 0x1bb88d50},
    {0x08af8ea1, 0x114b9862, 0x2013b0c4, 0x12276188, 0x244e6710},
    {0x08fac5a1, 0x11eed742, 0x20b6f464, 0x2f6d34c8, 0x33b1e990},
    {0x0a8ac161, 0x150796c2, 0x285f2924, 0x1aacd248, 0x354bb430},
    {0x0acc8aa1, 0x158792e2, 0x28df2264, 0x2ba0c368, 0x2e917570},
    {0x0cbc6e61, 0x1965b7a2, 0x317b6824, 0x175bbc48, 0x2d1af890},
    {0x0cdf7261, 0x19a9dd22, 0x31b48644, 0x3f990968, 0x23d592d0},
    {0x0eb9d061, 0x1d6dd4c2, 0x3915dd84, 0x0bfb4f08, 0x1438edb0},
    {0x0ed9d041, 0x1da78482, 0x39df0904, 0x233a9208, 0x16f1a130},
    {0x10bdb201, 0x2160e2e2, 0x2fb11f24, 0x321238a8, 0x0a247150},
    {0x11c71941, 0x205e3282, 0x3aa2e2a4, 0x0f4536e8, 0x1d5a6a70},
    {0x1243cc61, 0x24870cc2, 0x030e1984, 0x045eb308, 0x08afe2b0},
    {0x135237e1, 0x251468a2, 0x3f983d44, 0x0a9dfa88, 0x15269e70},
    {0x14677101, 0x28ce1602, 0x2b82aba4, 0x2ed550e8, 0x27aa55d0},
    {0x14712281, 0x28e24262, 0x27c46fa4, 0x3a383428, 0x026de850},
    {0x149b7321, 0x292dbca2, 0x3f30f944, 0x13112e88, 0x25525bf0},
    {0x1588eec1, 0x29f6e462, 0x0f1d74c4, 0x1cddd588, 0x39ac9710},
    {0x165eda01, 0x2ca99522, 0x0b47aa44, 0x141f51a8, 0x282a8750},
    {0x170aa0a1, 0x2c57c142, 0x12af1684, 0x254ca9a8, 0x02dbd350},
    {0x1864ede1, 0x30c97fc2, 0x3386dea4, 0x379d1c68, 0x3fbeb8d0},
    {0x1871be81, 0x30e378e2, 0x3fc64dc4, 0x237c2788, 0x1a1fcaf0},
    {0x18ebefe1, 0x31c5cf62, 0x2bdb0e64, 0x1ff49868, 0x3db930d0},
    {0x19e872e1, 0x307d8ea2, 0x16e69a24, 0x2dcd3328, 0x2e37e130},
    {0x1a4d5e21, 0x3484cfe2, 0x13096c64, 0x25dcacc8, 0x3277d990},
    {0x1b1abc81, 0x355ef902, 0x07d6ae04, 0x0cdd5ae8, 0x19a1ef30},
    {0x1c53cc41, 0x38a74482, 0x1f4e5504, 0x3df7f608, 0x169f36f0},
    {0x1d52af21, 0x38555ba2, 0x2eaa0ea4, 0x034398a8, 0x04773150},
    {0x1ede2a21, 0x3daed0e2, 0x331fb5c4, 0x2e7deb88, 0x16e9c7b0},
    {0x1fb49981, 0x3df93302, 0x2b76e324, 0x067d6768, 0x0ceeeff0},
    {0x20cdad01, 0x2f80da02, 0x32716ee4, 0x0ae20728, 0x15c408b0},
    {0x216e9561, 0x1259aac2, 0x24a7d0a4, 0x1b4f0548, 0x341a8a90},
    {0x22ae35e1, 0x0f4eebc2, 0x1cdfc724, 0x39ad9ee8, 0x3b19b970},
    {0x237d11a1, 0x3f34a342, 0x07b94124, 0x0cbcf648, 0x19679f30},
    {0x242c7b41, 0x024ae622, 0x0487dce4, 0x090f2dc8, 0x105cdb90},
    {0x243f63a1, 0x26659da2, 0x22cb3b44, 0x2b8df068, 0x3a70bcd0},
    {0x27bc2cc1, 0x3ad5d982, 0x03ab5f04, 0x04fbd568, 0x09eac6d0},
    {0x27cf34c1, 0x1379e982, 0x2414ef04, 0x162967e8, 0x2c45f630},
    {0x28811d41, 0x27023a82, 0x3bb47264, 0x02d80fa8, 0x05ad9830},
    {0x28f19d81, 0x0fe33b02, 0x1d3673e4, 0x389c5e28, 0x2f2f85b0},
    {0x2a2cb0a1, 0x2e47e142, 0x268f3684, 0x3700eaa8, 0x17d126f0},
    {0x2a3fa881, 0x066bd102, 0x0cc38604, 0x19870928, 0x319e1250},
    {0x2d0fa7a1, 0x23d1c8e2, 0x3e7365c4, 0x06e63f88, 0x0dcc78b0},
    {0x2db8e141, 0x0796fe82, 0x0ddd44e4, 0x195db5c8, 0x305ceb90},
    {0x2f1cbb41, 0x0ebdf3a2, 0x1d6fc664, 0x385bade8, 0x22a3dbd0},
    {0x2fc87981, 0x2a3d9f02, 0x2266b964, 0x32cd72c8, 0x13878ef0},
    {0x32b2d5c1, 0x3b651782, 0x2a3a2ae4, 0x0a63d028, 0x14c71c50},
    {0x33cb2541, 0x0afdca82, 0x15e0cfe4, 0x28b14528, 0x3f625650},
    {0x34ecc6e1, 0x3bcdace2, 0x271fd9c4, 0x1ebb9788, 0x3d63aa30},
    {0x3527dba1, 0x1fff5c22, 0x3c53d444, 0x0ea74488, 0x1d4e6510},
    {0x35744ba1, 0x36182ea2, 0x3227d8a4, 0x3a4f0d48, 0x2a899a90},
    {0x35ec0a41, 0x239a9022, 0x0f77a044, 0x1cbf4088, 0x396c9510},
    {0x3a9fc2c1, 0x0f21f622, 0x1d931fe4, 0x38f63868, 0x0bec70d0},
    {0x3b245f41, 0x2ab80762, 0x0b678b24, 0x143f13a8, 0x2869a2b0},
    {0x3b67cb41, 0x037f7de2, 0x055397c4, 0x091728e8, 0x119e51d0},
    {0x3b8f0e41, 0x1a759a62, 0x34eb3224, 0x07cde2a8, 0x0cf09fb0},
    {0x3ed25521, 0x0ba44642, 0x14f867e4, 0x29eda4a8, 0x2676c950},
    {0x3fbbb981, 0x3735f302, 0x263b76a4, 0x0664fde8, 0x0cc96fd0},
};

/* BCH codes over GF(2^5)
 */
uint32_t compute_checksum4(const uint8_t* data, int len, int gen) {
    uint32_t l = 0;

    while (len > 6) {
        uint8_t c = *(data++);
        len--;
        uint8_t t = (l >> 25) ^ c;
        l = (l & 0x1FFFFFFULL) << 5;
        l ^= -((t & 1) != 0) & tbl[gen][0];
        l ^= -((t & 2) != 0) & tbl[gen][1];
        l ^= -((t & 4) != 0) & tbl[gen][2];
        l ^= -((t & 8) != 0) & tbl[gen][3];
        l ^= -((t & 16) != 0) & tbl[gen][4];
    }

    uint32_t f = *(data++);
    f <<= 5;
    f |= *(data++);
    f <<= 5;
    f |= *(data++);
    f <<= 5;
    f |= *(data++);
    f <<= 5;
    f |= *(data++);
    f <<= 5;
    f |= *(data++);
    return l ^ f;
}

/*
inline uint32_t rstable(uint16_t x) {
    uint32_t t = ((x << 1) ^ x) << 1 ^ x;
    uint32_t h = x >> 4 | (t >> 7) << 10 | (t >> 9) << 20;
    return ((x & 0xF) << 6 | ((t & 0x7F) | (t & 0x1FF) << 8) << 13) ^ h ^ (h << 3);
}

// RS code over GF(2^10), max length 2046 (1023*2), HD 4
uint32_t compute_checksum3(const uint8_t* data, int len, int gen) {
    uint32_t l = 0;

    while (len > 6) {
        uint16_t c = *(data++);
        len--;
        if (len & 1) {
            c <<= 5;
            c |= *(data++);
            len--;
        }
        l = ((l & 0xFFFFF) << 10) ^ rstable((l >> 20) ^ c);
    }

    uint32_t f = *(data++);
    f <<= 5;
    f |= *(data++);
    f <<= 5;
    f |= *(data++);
    f <<= 5;
    f |= *(data++);
    f <<= 5;
    f |= *(data++);
    f <<= 5;
    f |= *(data++);
    return l ^ f;
}

// different CRC-6's over the individual bits of 5-bit characters, max length 57, HD 3
uint32_t compute_checksum2(const uint8_t* data, int len, int gen) {
    uint32_t ret = 0;
    while (len--) {
        ret = ((ret & 0x1FFFFFF) << 5) ^ (((ret >> 25) & 0x15) * 0x2000421) ^ (((ret >> 25) & 0xA) * 0x2100021) ^ *(data++);
    }
    return ret;
}

// A CRC-6 over the individual bits of 5-bit characters. max length 57, HD 3
uint32_t compute_checksum1(const uint8_t* data, int len, int gen) {
    uint32_t ret = 0;
    while (len--) {
        ret = ((ret & 0x1FFFFFF) << 5) ^ ((ret >> 25) * 0x2000421) ^ *(data++);
    }
    return ret;
}

// A CRC-30 over all the bits. max length infinity, HD 2
uint32_t compute_checksum0(const uint8_t* data, int len, int gen) {
    uint32_t ret = 0;
    while (len--) {
        uint8_t c = *(data++);
        for (int i = 0; i < 5; i++) {
            ret = ((ret & 0x1FFFFFFF) << 1) ^ ((ret >> 29) * 0x2030b9c7) ^ ((c >> i) & 1);
        }
    }
    return ret;
}
*/

/*
struct Checksum {
    const char* name;
    uint32_t (*fun)(const uint8_t*, int, int);
    int gen;
} static const checksums[] = {
//    {"CRC30", &compute_checksum0, 0},
//    {"CRC-6-sliced", &compute_checksum1, 0, 3},
//    {"2xCRC-6-sliced", &compute_checksum2, 0, 3},
//    {"RS_10", &compute_checksum3, 0},
};
*/

class Rander {
    unsigned char data[4096];
    int pos;

    void Step() {
        if ((pos & 0xFFF) == 0) {
            for (int i = 0; i < 4096; i += sizeof(unsigned long long)) {
                _rdrand64_step((unsigned long long*)(data + i));
            }
            pos = 0;
        }
    }

public:
    Rander() {
        memset(data, 0, sizeof(data));
        pos = 0;
    }

    uint8_t GetByte() {
        Step();
        uint8_t ret = data[pos];
        pos++;
        return ret;
    }

    uint8_t GetInt(uint8_t max, int bits) {
        uint8_t r;
        do {
            r = GetByte() & ((1 << bits) - 1);
        } while (r >= max);
        return r;
    }
};

#define LEN 91
#define LENBITS 7

#define CHECKSUMS (sizeof(tbl)/sizeof(tbl[0]))
// #define CHECKSUMS (sizeof(checksums)/sizeof(checksums[0]))
// #define CHECKSUMS 180

#define MINERR 4
#define MAXERR 8


struct CRCOutputs {
    uint32_t val[LEN][32][CHECKSUMS];

    CRCOutputs() {
        unsigned char data[LEN] = {0};
        for (int pos = 0; pos < LEN; pos++) {
            for (int v = 0; v < 32; v++) {
                data[pos] = v;
                for (unsigned int c = 0; c < CHECKSUMS; c++) {
//                  val[pos][v][c] = checksums[c].fun(data, LEN, checksums[c].gen);
                    val[pos][v][c] = compute_checksum4(data, LEN, c);
                }
            }
            data[pos] = 0;
        }
    }
};

const CRCOutputs outputs;

struct Results {
    uint64_t fails[CHECKSUMS];
    uint64_t count;

    Results() : fails{0}, count(0) {}

    Results& operator+=(const Results& x) {
        for (unsigned int i = 0; i < CHECKSUMS; i++) {
            fails[i] += x.fails[i];
        }
        count += x.count;
        return *this;
    }
};

void test(int errors, uint64_t loop, Results* ret, Rander& rng, const std::set<int>& which) {
    Results res = {};
    for (uint64_t i = 0; i < loop; i++) {
        uint32_t crc[CHECKSUMS] = {0};
        int errpos[MAXERR];
        for (int j = 0; j < errors; j++) {
            int ok;
            do {
                errpos[j] = rng.GetInt(LEN, LENBITS);
                ok = 1;
                for (int k = 0; k < j; k++) {
                    if (errpos[j] == errpos[k]) ok = 0;
                }
            } while (!ok);
            int mis = 1 + rng.GetInt(31, 5);
            for (int c : which) {
                crc[c] ^= outputs.val[errpos[j]][mis][c];
            }
        }
        for (int c : which) {
            res.fails[c] += (crc[c] == outputs.val[0][0][c]);
        }
        res.count++;
    }
    *ret = res;
}

static Results allresults[MAXERR-MINERR+1];
std::mutex cs_allresults;

void thread_crc() {
    Rander rng;
    std::set<int> which;
    for (unsigned int c = 0; c < CHECKSUMS; c++) {
        which.insert(c);
    }

    do {
        Results r[MAXERR-MINERR+1];
        for (int e = 0; e < MAXERR-MINERR+1; e++) {
            test(e + MINERR, 1 << 16, &r[e], rng, which);
        }
        {
            std::unique_lock<std::mutex> lock(cs_allresults);
            for (int e = 0; e < MAXERR-MINERR+1; e++) {
                allresults[e] += r[e];
            }
            for (auto it = which.begin(); it != which.end(); ) {
                bool keep = false;
                for (int e = 0; e < MAXERR-MINERR+1; e++) {
                    if (!allresults[e].fails[*it]) {
                        keep = true;
                        break;
                    }
                }
                if (!keep) {
                    which.erase(it++);
                } else {
                    ++it;
                }
            }
        }
    } while(true);
}

void thread_dump() {
    do {
        sleep(10);
        std::unique_lock<std::mutex> lock(cs_allresults);
        printf("#counts: ");
        for (int e = 0; e < MAXERR-MINERR+1; e++) {
            printf("%llu,", (unsigned long long)allresults[e].count);
        }
        printf("\n");
        for (unsigned int i = 0; i < CHECKSUMS; i++) {
            bool ok = true;
            for (int e = 0; e < MAXERR-MINERR+1; e++) {
                if (allresults[e].fails[i]) {
                    ok = false;
                    break;
                }
            }
            if (ok) {
                printf("\"BCH 0x%08x 0x%08x 0x%08x 0x%08x 0x%08x\"", tbl[i][0], tbl[i][1], tbl[i][2], tbl[i][3], tbl[i][4]);
                for (int e = 0; e < MAXERR-MINERR+1; e++) {
                    printf(",% 10g", ((double)allresults[e].fails[i]) / allresults[e].count * 1073741824.0);
                }
                printf("\n");
            }
        }
        printf("\n\n");
    } while(true);
}

}

int main(int argc, char** argv) {
    setbuf(stdout, NULL);
    int threads = argc > 1 ? strtol(argv[1], NULL, 10) : 1;
    for (int t = 0; t < threads; t++) {
        std::thread th(thread_crc);
        th.detach();
    }
    thread_dump();
    return 0;
}
