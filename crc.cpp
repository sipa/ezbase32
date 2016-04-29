#include <thread>
#include <mutex>

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <immintrin.h>

namespace {

/* BCH code over GF(2^5)
 * gen 0: HD5 degree6 length93
 * gen 1: HD4 degree5 length1023 * degree1
 * gen 2: HD4 degree6 length1023
 */
uint32_t compute_checksum4(const uint8_t* data, int len, int gen) {
    uint32_t l = 0;
    static const uint32_t tbl[4][32] = {
        {0, 792092555, 373043126, 956985405, 745200617, 56394850, 978777183, 358656980, 315074418, 1039268089, 83593412, 734910287, 1051557019, 295379728, 680612653, 128388262, 630147681, 179901930, 867395031, 478493276, 166039944, 651415043, 534363710, 821027253, 927526163, 410956440, 561300133, 239246638, 455227130, 873752945, 220076364, 573064903},
        {0, 1020869833, 61667730, 1064682843, 77039385, 944683984, 121339531, 1004683842, 154077615, 905327974, 176418877, 911920372, 230198966, 828486271, 236222244, 852315117, 294687582, 759803799, 305949388, 786781701, 352785479, 702425230, 380366293, 712167708, 415210225, 610707000, 454255459, 667574186, 472391144, 552673569, 528802938, 593352883},
        {0, 624561345, 37132674, 655136067, 74265348, 559366085, 106670726, 593602119, 147382952, 770895465, 184449834, 801404907, 212193708, 698342765, 244533294, 732513519, 293852656, 884594993, 330648690, 915357875, 367978228, 819259957, 400064374, 853701559, 424359768, 1014053785, 461221594, 1044881947, 489030748, 941361309, 521182686, 975868191},
    };

    while (len > 6) {
        uint8_t c = *(data++);
        len--;
        l = ((l & 0x1FFFFFFULL) << 5) ^ tbl[gen][(l >> 25) ^ c];
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

inline uint32_t rstable(uint16_t x) {
    uint32_t t = ((x << 1) ^ x) << 1 ^ x;
    uint32_t h = x >> 4 | (t >> 7) << 10 | (t >> 9) << 20;
    return ((x & 0xF) << 6 | ((t & 0x7F) | (t & 0x1FF) << 8) << 13) ^ h ^ (h << 3);
}

/* RS code over GF(2^10), max length 2046 (1023*2), HD 4 */
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

/* 2 different CRC-6's over the individual bits of 5-bit characters, max length 57, HD 3 */
uint32_t compute_checksum2(const uint8_t* data, int len, int gen) {
    uint32_t ret = 0;
    while (len--) {
        ret = ((ret & 0x1FFFFFF) << 5) ^ (((ret >> 25) & 0x15) * 0x2000421) ^ (((ret >> 25) & 0xA) * 0x2100021) ^ *(data++);
    }
    return ret;
}

/* A CRC-6 over the individual bits of 5-bit characters. max length 57, HD 3 */
uint32_t compute_checksum1(const uint8_t* data, int len, int gen) {
    uint32_t ret = 0;
    while (len--) {
        ret = ((ret & 0x1FFFFFF) << 5) ^ ((ret >> 25) * 0x2000421) ^ *(data++);
    }
    return ret;
}

/* A CRC-30 over all the bits. max length infinity, HD 2 */
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

struct Checksum {
    const char* name;
    uint32_t (*fun)(const uint8_t*, int, int);
    int gen;
    int hd;
} static const checksums[] = {
    {"CRC-30", &compute_checksum0, 0, 2},
//    {"CRC-6-sliced", &compute_checksum1, 0, 3},
//    {"2xCRC-6-sliced", &compute_checksum2, 0, 3},
    {"RS-10bit", &compute_checksum3, 0, 4},
    {"BCH(m4d6l93)", &compute_checksum4, 0, 5},
    {"BCH(m3d5l1023+m1)", &compute_checksum4, 1, 4},
    {"BCH(m3d6l1023)", &compute_checksum4, 2, 4},
};

#define ROT(x,i) (((x) << (i)) | (x) >> (32-(i)))
#define QR(a,b,c,d) do { a += b; d ^= a; d = ROT(d,16); c += d; b ^= c; b = ROT(b,12); a += b; d ^= a; d = ROT(d,8); c += d; b ^= c; b = ROT(b,7); } while(false)

class Rander {
    uint32_t state[16];
    int pos;

    void Step() {
        uint32_t q[16];
        if ((pos & 0xFFFF) == 0) {
            // Add hardware randomness
            for (int i = 0; i < 16; i++) {
                _rdrand32_step(&q[i]);
            }
            for (int i = 0; i < 16; i++) {
                state[i] += q[i];
            }
        } else if ((pos & 0x1F) == 0) {
            // One double round of ChaCha20
            QR(state[0], state[4], state[8], state[12]);
            QR(state[1], state[5], state[9], state[13]);
            QR(state[2], state[6], state[10], state[14]);
            QR(state[3], state[7], state[11], state[15]);
            QR(state[0], state[5], state[10], state[15]);
            QR(state[1], state[6], state[11], state[12]);
            QR(state[2], state[7], state[8], state[13]);
            QR(state[3], state[4], state[9], state[14]);
        }
    }

public:
    Rander() {
        memset(state, 0, sizeof(state));
        pos = 0;
    }

    uint8_t GetByte() {
        Step();
        uint8_t ret = state[(pos >> 2) & 7] >> (pos & 3);
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

#define LEN 93
#define LENBITS 7

#define CHECKSUMS (sizeof(checksums)/sizeof(checksums[0]))

#define MAXERR 8

struct CRCOutputs {
    uint32_t val[LEN][32][CHECKSUMS];

    CRCOutputs() {
        unsigned char data[LEN] = {0};
        for (int pos = 0; pos < LEN; pos++) {
            for (int v = 0; v < 32; v++) {
                data[pos] = v;
                for (unsigned int c = 0; c < CHECKSUMS; c++) {
                    val[pos][v][c] = checksums[c].fun(data, LEN, checksums[c].gen);
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

void test(int errors, uint64_t loop, Results* ret, Rander& rng) {
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
            for (unsigned int c = 0; c < CHECKSUMS; c++) {
                crc[c] ^= outputs.val[errpos[j]][mis][c];
            }
        }
        for (unsigned int c = 0; c < CHECKSUMS; c++) {
            res.fails[c] += (crc[c] == outputs.val[0][0][c]);
        }
        res.count++;
    }
    *ret = res;
}

static Results allresults[MAXERR];
std::mutex cs_allresults;

void thread_crc() {
    Rander rng;
    do {
        Results r;
        int e = 1 + (rng() % MAXERR);
        test(e, 1 << 20, &r, rng);
        {
            std::unique_lock<std::mutex> lock(cs_allresults);
            allresults[e - 1] += r;
        }
    } while(true);
}

void thread_dump() {
    do {
        sleep(10);
        std::unique_lock<std::mutex> lock(cs_allresults);
        printf("#counts: ");
        for (int e = 0; e < MAXERR; e++) {
            printf("%llu,", (unsigned long long)allresults[e].count);
        }
        printf("\n");
        for (unsigned int i = 0; i < CHECKSUMS; i++) {
            printf("\"%s\"", checksums[i].name);
            for (int e = 0; e < MAXERR; e++) {
                printf(",%g", ((double)allresults[e].fails[i]) / allresults[e].count * 1000000000.0);
            }
            printf("\n");
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
