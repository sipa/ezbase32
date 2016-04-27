#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <immintrin.h>

/* BCH code over GF(2^5)
 * gen 0: HD5 degree6 length93
 * gen 1: HD4 degree5 length1023 * degree1
 * gen 2: HD4 degree6 length1023
 */
static uint32_t compute_checksum4(const uint8_t* data, int len, int gen) {
    uint32_t l = 1;
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

static inline uint32_t rstable(uint16_t x) {
    uint32_t t = ((x << 1) ^ x) << 1 ^ x;
    uint32_t h = x >> 4 | (t >> 7) << 10 | (t >> 9) << 20;
    return ((x & 0xF) << 6 | ((t & 0x7F) | (t & 0x1FF) << 8) << 13) ^ h ^ (h << 3);
}

/* RS code over GF(2^10), max length 2046 (1023*2), HD 4 */
static uint32_t compute_checksum3(const uint8_t* data, int len) {
    uint32_t l = rstable(1);

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
/*static uint32_t compute_checksum2(const uint8_t* data, int len) {
    uint32_t ret = 3;
    while (len--) {
        ret = ((ret & 0x1FFFFFF) << 5) ^ (((ret >> 25) & 0x15) * 0x2000421) ^ (((ret >> 25) & 0xA) * 0x2100021) ^ *(data++);
    }
    return ret;
}*/

/* A CRC-6 over the individual bits of 5-bit characters. max length 57, HD 3 */
/*static uint32_t compute_checksum1(const uint8_t* data, int len) {
    uint32_t ret = 3;
    while (len--) {
        ret = ((ret & 0x1FFFFFF) << 5) ^ ((ret >> 25) * 0x2000421) ^ *(data++);
    }
    return ret;
}*/

uint32_t insecure_rand_Rz = 11;
uint32_t insecure_rand_Rw = 11;
static inline uint32_t insecure_rand(void)
{
    uint32_t rdr = 0;
    _rdrand32_step(&rdr);
    insecure_rand_Rz = 36969 * (insecure_rand_Rz & 65535) + (insecure_rand_Rz >> 16);
    insecure_rand_Rw = 18000 * (insecure_rand_Rw & 65535) + (insecure_rand_Rw >> 16);
    return (insecure_rand_Rw << 16) + (insecure_rand_Rw >> 16) + insecure_rand_Rz + rdr;
}

#define LEN 93

void test(int errors, int swaps, uint64_t loop) {
    uint8_t data[LEN];
    uint8_t olddata[LEN];
    for (int j = 0; j < LEN; j++) {
        data[j] = insecure_rand() % 32;
    }
    uint32_t crc3 = compute_checksum3(data, LEN);
    uint32_t crc40 = compute_checksum4(data, LEN, 0);
    uint32_t crc41 = compute_checksum4(data, LEN, 1);
    uint32_t crc42 = compute_checksum4(data, LEN, 2);
    uint64_t /*fail1 = 0, fail2 = 0,*/ fail3 = 0, fail40 = 0, fail41 = 0, fail42 = 0, fail43 = 0;
    for (uint64_t i = 0; i < loop; i++) {
        int errpos[8];
        memcpy(olddata, data, LEN);
        do {
            for (int j = 0; j < errors; j++) {
                int ok;
                do {
                    errpos[j] = insecure_rand() % LEN;
                    ok = 1;
                    for (int k = 0; k < j; k++) {
                        if (errpos[j] == errpos[k]) ok = 0;
                    }
                } while (!ok);
                int mis;
                do {
                    mis = insecure_rand() % 32;
                } while (mis == 0);
                data[errpos[j]] ^= mis;
            }
            for (int j = 0; j < swaps; j++) {
                int s;
                do {
                    s = insecure_rand() % (LEN - 1);
                } while (data[s] == data[s + 1]);
                uint8_t t = data[s];
                data[s] = data[s + 1];
                data[s + 1] = t;
            }
        } while (memcmp(olddata, data, LEN) == 0);
        uint32_t crc3n = compute_checksum3(data, LEN);
        uint32_t crc40n = compute_checksum4(data, LEN, 0);
        uint32_t crc41n = compute_checksum4(data, LEN, 1);
        uint32_t crc42n = compute_checksum4(data, LEN, 2);
        fail3 += (crc3n == crc3);
        fail40 += (crc40n == crc40);
        fail41 += (crc41n == crc41);
        fail42 += (crc42n == crc42);
        crc3 = crc3n;
        crc40 = crc40n;
        crc41 = crc41n;
        crc42 = crc42n;
    }
    fprintf(stderr, "Out of %llu (HD%i mistakes, %i swaps): fails={CRC3:%llu, CRC4:[%llu,%llu,%llu,%llu]}\n", (unsigned long long)loop, errors, swaps, (unsigned long long)fail3, (unsigned long long)fail40, (unsigned long long)fail41, (unsigned long long)fail42, (unsigned long long)fail43);
}

int main(int argc, char** argv) {
    if (argc > 1) { insecure_rand_Rz *= strtol(argv[1], NULL, 10); }
    uint64_t n = ((uint64_t)1) << 25;
    for (int i = 0; i < 40; i++) {
       test(3, 0, n);
       test(4, 0, n);
       test(5, 0, n);
       test(6, 0, n);
       n *= 2;
    }
}
