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
  {0x53A0C81, 0x8F09902, 0x11E13204, 0x21526128, 0x12346650}
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

#define LEN 90
#define LENBITS 7

#define CHECKSUMS (sizeof(tbl)/sizeof(tbl[0]))
// #define CHECKSUMS (sizeof(checksums)/sizeof(checksums[0]))
// #define CHECKSUMS 180

#define MINERR 4
#define MAXERR 6


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
