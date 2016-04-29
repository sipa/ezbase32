#include <thread>
#include <mutex>

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <immintrin.h>

namespace {

/* BCH codes over GF(2^5) M=2 N=93 DISTANCE=5 DEGREE=6
 */
uint32_t compute_checksum4(const uint8_t* data, int len, int gen) {
    uint32_t l = 0;
    static const uint32_t tbl[16][32] = {
        {0x0, 0x203a88dc, 0x1e629197, 0x3e58194b, 0x3cc52301, 0x1cffabdd, 0x22a7b296, 0x29d3a4a, 0x278a43e2, 0x7b0cb3e, 0x39e8d275, 0x19d25aa9, 0x1b4f60e3, 0x3b75e83f, 0x52df174, 0x251779a8, 0x13f3be24, 0x33c936f8, 0xd912fb3, 0x2daba76f, 0x2f369d25, 0xf0c15f9, 0x31540cb2, 0x116e846e, 0x3479fdc6, 0x1443751a, 0x2a1b6c51, 0xa21e48d, 0x8bcdec7, 0x2886561b, 0x16de4f50, 0x36e4c78c}, // F=(x^5 + x^3 + x^2 + x + 1) E=(e^2 + 19*e + 25) alpha=(18*e + 28) powers=3..6
        {0x0, 0x39dba3af, 0x2f50c2be, 0x168b6111, 0x2513cb3, 0x3b8a9f1c, 0x2d01fe0d, 0x14da5da2, 0x4a27949, 0x3d79dae6, 0x2bf2bbf7, 0x12291858, 0x6f345fa, 0x3f28e655, 0x29a38744, 0x107824eb, 0x9444e92, 0x309fed3d, 0x26148c2c, 0x1fcf2f83, 0xb157221, 0x32ced18e, 0x2445b09f, 0x1d9e1330, 0xde637db, 0x343d9474, 0x22b6f565, 0x1b6d56ca, 0xfb70b68, 0x366ca8c7, 0x20e7c9d6, 0x193c6a79}, // F=(x^5 + x^3 + x^2 + x + 1) E=(e^2 + 23*e + 10) alpha=(21*e + 21) powers=3..6
        {0x0, 0x38825597, 0x3b043f0b, 0x3866a9c, 0x3e587ab6, 0x6da2f21, 0x55c45bd, 0x3dde102a, 0x36a2e5e9, 0xe20b07e, 0xda6dae2, 0x35248f75, 0x8fa9f5f, 0x3078cac8, 0x33fea054, 0xb7cf5c3, 0x27455fd2, 0x1fc70a45, 0x1c4160d9, 0x24c3354e, 0x191d2564, 0x219f70f3, 0x22191a6f, 0x1a9b4ff8, 0x11e7ba3b, 0x2965efac, 0x2ae38530, 0x1261d0a7, 0x2fbfc08d, 0x173d951a, 0x14bbff86, 0x2c39aa11}, // F=(x^5 + x^2 + 1) E=(e^2 + 22*e + 22) alpha=(14*e + 27) powers=3..6
        {0x0, 0xea19f1f, 0x1d433abb, 0x13e2a5a4, 0x38d671f3, 0x3677eeec, 0x25954b48, 0x2b34d457, 0x3bac77c3, 0x350de8dc, 0x26ef4d78, 0x284ed267, 0x37a0630, 0xddb992f, 0x1e393c8b, 0x1098a394, 0x3f1aff26, 0x31bb6039, 0x2259c59d, 0x2cf85a82, 0x7cc8ed5, 0x96d11ca, 0x1a8fb46e, 0x142e2b71, 0x4b688e5, 0xa1717fa, 0x19f5b25e, 0x17542d41, 0x3c60f916, 0x32c16609, 0x2123c3ad, 0x2f825cb2}, // F=(x^5 + x^2 + 1) E=(e^2 + 30*e + 17) alpha=(7*e + 9) powers=3..6
        {0x0, 0x38825597, 0x3b043f0b, 0x3866a9c, 0x3e587ab6, 0x6da2f21, 0x55c45bd, 0x3dde102a, 0x36a2e5e9, 0xe20b07e, 0xda6dae2, 0x35248f75, 0x8fa9f5f, 0x3078cac8, 0x33fea054, 0xb7cf5c3, 0x27455fd2, 0x1fc70a45, 0x1c4160d9, 0x24c3354e, 0x191d2564, 0x219f70f3, 0x22191a6f, 0x1a9b4ff8, 0x11e7ba3b, 0x2965efac, 0x2ae38530, 0x1261d0a7, 0x2fbfc08d, 0x173d951a, 0x14bbff86, 0x2c39aa11}, // F=(x^5 + x^2 + 1) E=(e^2 + 30*e + 21) alpha=(27*e + 31) powers=3..6
        {0x0, 0x2b04cde3, 0x23b977c6, 0x8bdba25, 0x32df84ec, 0x19db490f, 0x1166f32a, 0x3a623ec9, 0x13a289d8, 0x38a6443b, 0x301bfe1e, 0x1b1f33fd, 0x217d0d34, 0xa79c0d7, 0x2c47af2, 0x29c0b711, 0x24f5138b, 0xff1de68, 0x74c644d, 0x2c48a9ae, 0x162a9767, 0x3d2e5a84, 0x3593e0a1, 0x1e972d42, 0x37579a53, 0x1c5357b0, 0x14eeed95, 0x3fea2076, 0x5881ebf, 0x2e8cd35c, 0x26316979, 0xd35a49a}, // F=(x^5 + x^4 + x^3 + x + 1) E=(e^2 + 2*e + 2) alpha=(23*e + 10) powers=3..6
        {0x0, 0x2ee347a2, 0x33c655a4, 0x1d251206, 0xafc7748, 0x241f30ea, 0x393a22ec, 0x17d9654e, 0x15e3b470, 0x3b00f3d2, 0x2625e1d4, 0x8c6a676, 0x1f1fc338, 0x31fc849a, 0x2cd9969c, 0x23ad13e, 0x28b768d7, 0x6542f75, 0x1b713d73, 0x35927ad1, 0x224b1f9f, 0xca8583d, 0x118d4a3b, 0x3f6e0d99, 0x3d54dca7, 0x13b79b05, 0xe928903, 0x2071cea1, 0x37a8abef, 0x194bec4d, 0x46efe4b, 0x2a8db9e9}, // F=(x^5 + x^4 + x^2 + x + 1) E=(e^2 + 17*e + 7) alpha=(20*e + 5) powers=3..6
        {0x0, 0x229bd67c, 0x17238df1, 0x35b85b8d, 0x2cd71bcb, 0xe4ccdb7, 0x3bf4963a, 0x196f4046, 0xbae32b6, 0x2935e4ca, 0x1c8dbf47, 0x3e16693b, 0x2779297d, 0x5e2ff01, 0x305aa48c, 0x12c172f0, 0x15d8e065, 0x37433619, 0x2fb6d94, 0x2060bbe8, 0x390ffbae, 0x1b942dd2, 0x2e2c765f, 0xcb7a023, 0x1e76d2d3, 0x3ced04af, 0x9555f22, 0x2bce895e, 0x32a1c918, 0x103a1f64, 0x258244e9, 0x7199295}, // F=(x^5 + x^3 + 1) E=(e^2 + 9*e + 13) alpha=(31*e + 25) powers=3..6
        {0x0, 0x330b47f5, 0x13bbe4b1, 0x20b0a344, 0x24daa559, 0x17d1e2ac, 0x376141e8, 0x46a061d, 0x3fa8ca89, 0xca38d7c, 0x2c132e38, 0x1f1869cd, 0x1b726fd0, 0x28792825, 0x8c98b61, 0x3bc2cc94, 0xafcfe72, 0x39f7b987, 0x19471ac3, 0x2a4c5d36, 0x2e265b2b, 0x1d2d1cde, 0x3d9dbf9a, 0xe96f86f, 0x355434fb, 0x65f730e, 0x26efd04a, 0x15e497bf, 0x118e91a2, 0x2285d657, 0x2357513, 0x313e32e6}, // F=(x^5 + x^4 + x^3 + x + 1) E=(e^2 + e + 27) alpha=(5*e + 11) powers=3..6
        {0x0, 0x1acda58a, 0x358fcb14, 0x2f426e9e, 0x3b9bb721, 0x215612ab, 0xe147c35, 0x14d9d9bf, 0x27b3eb62, 0x3d7e4ee8, 0x123c2076, 0x8f185fc, 0x1c285c43, 0x6e5f9c9, 0x29a79757, 0x336a32dd, 0x1ff777e4, 0x53ad26e, 0x2a78bcf0, 0x30b5197a, 0x246cc0c5, 0x3ea1654f, 0x11e30bd1, 0xb2eae5b, 0x38449c86, 0x2289390c, 0xdcb5792, 0x1706f218, 0x3df2ba7, 0x19128e2d, 0x3650e0b3, 0x2c9d4539}, // F=(x^5 + x^3 + 1) E=(e^2 + 15*e + 4) alpha=(9*e + 2) powers=3..6
        {0x0, 0x215413af, 0x3718203e, 0x164c3391, 0x1b9dc047, 0x3ac9d3e8, 0x2c85e079, 0xdd1f3d6, 0x3496ec8e, 0x15c2ff21, 0x38eccb0, 0x22dadf1f, 0x2f0b2cc9, 0xe5f3f66, 0x18130cf7, 0x39471f58, 0x1f2d351c, 0x3e7926b3, 0x28351522, 0x961068d, 0x4b0f55b, 0x25e4e6f4, 0x33a8d565, 0x12fcc6ca, 0x2bbbd992, 0xaefca3d, 0x1ca3f9ac, 0x3df7ea03, 0x302619d5, 0x11720a7a, 0x73e39eb, 0x266a2a44}, // F=(x^5 + x^4 + x^3 + x + 1) E=(e^2 + 2*e + 31) alpha=(20*e + 9) powers=3..6
        {0x0, 0x1d0bfcf3, 0x39d98ddb, 0x24d27128, 0xa7d9b8b, 0x17766778, 0x33a41650, 0x2eafeaa3, 0x14e5b0b6, 0x9ee4c45, 0x2d3c3d6d, 0x3037c19e, 0x1e982b3d, 0x393d7ce, 0x2741a6e6, 0x3a4a5a15, 0x29cb6151, 0x34c09da2, 0x1012ec8a, 0xd191079, 0x23b6fada, 0x3ebd0629, 0x1a6f7701, 0x7648bf2, 0x3d2ed1e7, 0x20252d14, 0x4f75c3c, 0x19fca0cf, 0x37534a6c, 0x2a58b69f, 0xe8ac7b7, 0x13813b44}, // F=(x^5 + x^4 + x^3 + x^2 + 1) E=(e^2 + 13*e + 15) alpha=(31*e + 6) powers=3..6
        {0x0, 0x31df19b6, 0x1a70b351, 0x2bafaae7, 0x34e1613f, 0x53e7889, 0x2e91d26e, 0x1f4ecbd8, 0x13c23643, 0x221d2ff5, 0x9b28512, 0x386d9ca4, 0x2723577c, 0x16fc4eca, 0x3d53e42d, 0xc8cfd9b, 0x24546b26, 0x158b7290, 0x3e24d877, 0xffbc1c1, 0x10b50a19, 0x216a13af, 0xac5b948, 0x3b1aa0fe, 0x37965d65, 0x64944d3, 0x2de6ee34, 0x1c39f782, 0x3773c5a, 0x32a825ec, 0x19078f0b, 0x28d896bd}, // F=(x^5 + x^4 + x^3 + x^2 + 1) E=(e^2 + 26*e + 29) alpha=(26*e + 13) powers=3..6
        {0x0, 0x2302f606, 0x16954d2c, 0x3597bb2a, 0x2d2a3e58, 0xe28c85e, 0x3bbf7374, 0x18bd8572, 0xad0f9b9, 0x29d20fbf, 0x1c45b495, 0x3f474293, 0x27fac7e1, 0x4f831e7, 0x316f8acd, 0x126d7ccb, 0x15a1575b, 0x36a3a15d, 0x3341a77, 0x2036ec71, 0x388b6903, 0x1b899f05, 0x2e1e242f, 0xd1cd229, 0x1f71aee2, 0x3c7358e4, 0x9e4e3ce, 0x2ae615c8, 0x325b90ba, 0x115966bc, 0x24cedd96, 0x7cc2b90}, // F=(x^5 + x^3 + 1) E=(e^2 + 21*e + 24) alpha=(7*e + 9) powers=3..6
        {0x0, 0x352f955d, 0x734aa8d, 0x321b3fd0, 0xd1953fa, 0x3836c6a7, 0xa2df977, 0x3f026c2a, 0x1959fd23, 0x2c76687e, 0x1e6d57ae, 0x2b42c2f3, 0x1440aed9, 0x216f3b84, 0x13740454, 0x265b9109, 0x31d8a646, 0x4f7331b, 0x36ec0ccb, 0x3c39996, 0x3cc1f5bc, 0x9ee60e1, 0x3bf55f31, 0xedaca6c, 0x28815b65, 0x1daece38, 0x2fb5f1e8, 0x1a9a64b5, 0x2598089f, 0x10b79dc2, 0x22aca212, 0x1783374f}, // F=(x^5 + x^4 + x^2 + x + 1) E=(e^2 + 31*e + 4) alpha=(13*e + 29) powers=3..6
        {0x0, 0x3b2c919a, 0x3e1ba311, 0x537328b, 0x3625c2a7, 0xd09533d, 0x83e61b6, 0x3312f02c, 0x264b15ee, 0x1d678474, 0x1850b6ff, 0x237c2765, 0x106ed749, 0x2b4246d3, 0x2e757458, 0x1559e5c2, 0x684abdc, 0x3da83a46, 0x389f08cd, 0x3b39957, 0x30a1697b, 0xb8df8e1, 0xebaca6a, 0x35965bf0, 0x20cfbe32, 0x1be32fa8, 0x1ed41d23, 0x25f88cb9, 0x16ea7c95, 0x2dc6ed0f, 0x28f1df84, 0x13dd4e1e}, // F=(x^5 + x^2 + 1) E=(e^2 + 14*e + 19) alpha=(26*e + 31) powers=3..6
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
} static const checksums[] = {
    {"CRC30", &compute_checksum0, 0},
//    {"CRC-6-sliced", &compute_checksum1, 0, 3},
//    {"2xCRC-6-sliced", &compute_checksum2, 0, 3},
    {"RS_10", &compute_checksum3, 0},
    {"BCH_0", &compute_checksum4, 0},
    {"BCH_1", &compute_checksum4, 1},
    {"BCH_2", &compute_checksum4, 2},
    {"BCH_3", &compute_checksum4, 3},
    {"BCH_4", &compute_checksum4, 4},
    {"BCH_5", &compute_checksum4, 5},
    {"BCH_6", &compute_checksum4, 6},
    {"BCH_7", &compute_checksum4, 7},
    {"BCH_8", &compute_checksum4, 8},
    {"BCH_9", &compute_checksum4, 9},
    {"BCH_A", &compute_checksum4, 9},
    {"BCH_B", &compute_checksum4, 10},
    {"BCH_C", &compute_checksum4, 11},
    {"BCH_D", &compute_checksum4, 12},
    {"BCH_E", &compute_checksum4, 13},
    {"BCH_F", &compute_checksum4, 14},
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
        Results r[MAXERR];
        for (int e = 0; e < MAXERR; e++) {
            test(e + 1, 1 << 16, &r[e], rng);
        }
        {
            std::unique_lock<std::mutex> lock(cs_allresults);
            for (int e = 0; e < MAXERR; e++) {
                allresults[e] += r[e];
            }
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
