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

const uint32_t tbl[1][5] = {
  {0x3b6a57b2, 0x26508e6d, 0x1ea119fa, 0x3d4233dd, 0x2a1462b3}
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

static inline uint64_t rotl(const uint64_t x, int k) {
    return (x << k) | (x >> (64 - k));
}

class Rander {
    uint64_t s[2];
    uint32_t count;

    uint64_t next;
    int bits;

    void Produce(void) {
        const uint64_t s0 = s[0];
        uint64_t s1 = s[1];
        next = s0 + s1;
        bits = 64;

        s1 ^= s0;
        s[0] = rotl(s0, 55) ^ s1 ^ (s1 << 14); // a, b
        s[1] = rotl(s1, 36); // c
    }

    void Step() {
        static_assert(sizeof(unsigned long long) == sizeof(uint64_t), "Bad ULL length");
        if ((count & 0xFFFF) == 0) {
            s[0] = 0;
            s[1] = 0;
            _rdrand64_step((unsigned long long*)(s + 0));
            _rdrand64_step((unsigned long long*)(s + 1));
        }
        ++count;

        Produce();
    }

public:
    Rander() : count(0), bits(0) {}

    uint32_t GetBits(int bits_) {
        if (bits_ > bits) {
            Step();
        }

        uint32_t ret = next & ((1UL << bits_) - 1);
        next >>= bits_;
        bits -= bits_;
        return ret;
    }

    uint32_t GetInt(uint32_t range, int bits_) {
        do {
            uint32_t r = GetBits(bits_);
            if (r < range) return r;
        } while(true);
    }
};

#define LEN 71
#define LENBITS 9

#define MINSYM 5

#define CHECKSUMS (sizeof(tbl)/sizeof(tbl[0]))
// #define CHECKSUMS (sizeof(checksums)/sizeof(checksums[0]))
// #define CHECKSUMS 180

#define MINERR 5
#define MAXERR 10


struct CRCOutputs {
    uint32_t val[LEN][5][CHECKSUMS];

    CRCOutputs() {
        unsigned char data[LEN] = {0};
        for (int pos = 0; pos < LEN; pos++) {
            for (int v = 0; v < 5; v++) {
                data[pos] = (1 << v);
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
        char biterrs[LEN * 5] = {0};
        char symerrs[LEN] = {0};
        int syms = 0;
        for (int j = 0; j < errors && syms + errors - j >= MINSYM; j++) {
            while (true) {
                int bitpos = rng.GetInt(LEN * 5, LENBITS);
                if (biterrs[bitpos]) continue;
                int sympos = bitpos / 5;
                syms += 1 ^ symerrs[sympos];
                biterrs[bitpos] = 1;
                symerrs[sympos] = 1;
//                for (int c : which) {
                int c = 0;
                    crc[c] ^= outputs.val[sympos][bitpos % 5][c];
//                }
                break;
            }
        }
        if (syms < MINSYM) {
            continue;
        }
        for (int c : which) {
            res.fails[c] += (crc[c] == 0);
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
//            if (ok) {
                printf("\"BCH 0x%08x 0x%08x 0x%08x 0x%08x 0x%08x\"", tbl[i][0], tbl[i][1], tbl[i][2], tbl[i][3], tbl[i][4]);
                for (int e = 0; e < MAXERR-MINERR+1; e++) {
                    printf(",% 10g", ((double)allresults[e].fails[i]) / allresults[e].count * 1073741824.0);
                }
                printf("\n");
//            }
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
