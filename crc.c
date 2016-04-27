#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static inline int rstable(uint16_t x) {
    uint32_t t = ((x << 1) ^ x) << 1 ^ x;
    uint32_t h = x >> 4 | (t >> 7) << 10 | (t >> 9) << 20;
    return ((x & 0xF) << 6 | ((t & 0x7F) | (t & 0x1FF) << 8) << 13) ^ h ^ (h << 3);
}

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

static uint32_t compute_checksum2(const uint8_t* data, int len) {
    uint32_t ret = 3;
    while (len--) {
        ret = ((ret & 0x1FFFFFF) << 5) ^ (((ret >> 25) & 0x15) * 0x2000421) ^ (((ret >> 25) & 0xA) * 0x2100021) ^ *(data++);
    }
    return ret;
}

uint32_t insecure_rand_Rz = 11;
uint32_t insecure_rand_Rw = 11;
static inline uint32_t insecure_rand(void)
{
    insecure_rand_Rz = 36969 * (insecure_rand_Rz & 65535) + (insecure_rand_Rz >> 16);
    insecure_rand_Rw = 18000 * (insecure_rand_Rw & 65535) + (insecure_rand_Rw >> 16);
    return (insecure_rand_Rw << 16) + (insecure_rand_Rw >> 16) + insecure_rand_Rz;
}

#define LEN 61

void test(int errors, int swaps, uint64_t loop) {
    uint8_t data[LEN];
    uint8_t olddata[LEN];
    for (int j = 0; j < LEN; j++) {
        data[j] = insecure_rand() % 32;
    }
    uint32_t crc2 = compute_checksum2(data, LEN);
    uint32_t crc3 = compute_checksum3(data, LEN);
    uint64_t fail1 = 0, fail2 = 0, fail3 = 0;
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
                    s = insecure_rand() % 39;
                } while (data[s] == data[s + 1]);
                uint8_t t = data[s];
                data[s] = data[s + 1];
                data[s + 1] = t;
            }
        } while (memcmp(olddata, data, LEN) == 0);
        uint32_t crc2n = compute_checksum2(data, LEN);
        uint32_t crc3n = compute_checksum3(data, LEN);
        fail2 += (crc2n == crc2);
        fail3 += (crc3n == crc3);
        crc2 = crc2n;
        crc3 = crc3n;
    }
    fprintf(stderr, "Out of %llu (HD%i mistakes, %i swaps): fails={CRC2:%llu, CRC3:%llu}\n", (unsigned long long)loop, errors, swaps, (unsigned long long)fail2, (unsigned long long)fail3);
}

int main(int argc, char** argv) {
    if (argc > 1) { insecure_rand_Rz *= strtol(argv[1], NULL, 10); }
    uint64_t n = 10000000;
    for (int i = 0; i < 40; i++) {
       test(3, 0, n);
       test(4, 0, n);
       test(5, 0, n);
       n *= 2;
    }
}
