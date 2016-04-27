#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <immintrin.h>

/* BCH code over GF(2^5), max length 93 (1023/11), HD 5 */
static uint32_t compute_checksum4(const uint8_t* data, int len, int gen) {
    uint32_t l = 1;
    static const uint32_t tbl[4][32] = {
        {0, 792092555, 373043126, 956985405, 745200617, 56394850, 978777183, 358656980, 315074418, 1039268089, 83593412, 734910287, 1051557019, 295379728, 680612653, 128388262, 630147681, 179901930, 867395031, 478493276, 166039944, 651415043, 534363710, 821027253, 927526163, 410956440, 561300133, 239246638, 455227130, 873752945, 220076364, 573064903},
        {0, 318034653, 636040511, 924706786, 59828827, 291780742, 645521252, 881659321, 91309107, 394630894, 546861324, 845980625, 115453544, 336952501, 587766615, 838617482, 145890371, 440848030, 760949116, 1068497825, 186840600, 433464517, 785105703, 1010766330, 230905968, 523733677, 673869135, 987676562, 240366123, 480731382, 733644564, 961435081},
        {0, 559158191, 924327998, 374092689, 463323207, 986305512, 746971257, 231863254, 882306190, 365100833, 59690160, 584769311, 789261513, 241123174, 403901687, 960962392, 523056412, 1048127155, 674567458, 157353613, 78706011, 635758324, 866702693, 318555850, 733731218, 183487037, 480508332, 1039657475, 807803349, 292686458, 121518571, 644491844},
        {0, 77639975, 155279950, 232903529, 282219580, 343074075, 429078130, 489916245, 564411485, 621055354, 686120467, 742780724, 829843553, 903272774, 943131183, 1016576776, 186153119, 263768504, 39394001, 117025782, 466234531, 527080836, 311119597, 371982282, 716981442, 773649893, 603793036, 660445099, 980291838, 1053729241, 858747568, 932168599},
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

#define LEN 40

void test(int errors, int swaps, uint64_t loop) {
    uint8_t data[LEN];
    uint8_t olddata[LEN];
    for (int j = 0; j < LEN; j++) {
        data[j] = insecure_rand() % 32;
    }
/*    uint32_t crc1 = compute_checksum1(data, LEN);
    uint32_t crc2 = compute_checksum2(data, LEN);*/
    uint32_t crc3 = compute_checksum3(data, LEN);
    uint32_t crc40 = compute_checksum4(data, LEN, 0);
    uint32_t crc41 = compute_checksum4(data, LEN, 1);
    uint32_t crc42 = compute_checksum4(data, LEN, 2);
    uint32_t crc43 = compute_checksum4(data, LEN, 3);
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
                    s = insecure_rand() % 39;
                } while (data[s] == data[s + 1]);
                uint8_t t = data[s];
                data[s] = data[s + 1];
                data[s + 1] = t;
            }
        } while (memcmp(olddata, data, LEN) == 0);
/*        uint32_t crc1n = compute_checksum1(data, LEN);
        uint32_t crc2n = compute_checksum2(data, LEN);*/
        uint32_t crc3n = compute_checksum3(data, LEN);
        uint32_t crc40n = compute_checksum4(data, LEN, 0);
        uint32_t crc41n = compute_checksum4(data, LEN, 1);
        uint32_t crc42n = compute_checksum4(data, LEN, 2);
        uint32_t crc43n = compute_checksum4(data, LEN, 3);
/*        fail1 += (crc1n == crc1);
        fail2 += (crc2n == crc2);*/
        fail3 += (crc3n == crc3);
        fail40 += (crc40n == crc40);
        fail41 += (crc41n == crc41);
        fail42 += (crc42n == crc42);
        fail43 += (crc43n == crc43);
/*        crc1 = crc1n;
        crc2 = crc2n;*/
        crc3 = crc3n;
        crc40 = crc40n;
        crc41 = crc41n;
        crc42 = crc42n;
        crc43 = crc43n;
    }
    fprintf(stderr, "Out of %llu (HD%i mistakes, %i swaps): fails={CRC3:%llu, CRC4:[%llu,%llu,%llu,%llu]}\n", (unsigned long long)loop, errors, swaps, (unsigned long long)fail3, (unsigned long long)fail40, (unsigned long long)fail41, (unsigned long long)fail42, (unsigned long long)fail43);
}

int main(int argc, char** argv) {
    if (argc > 1) { insecure_rand_Rz *= strtol(argv[1], NULL, 10); }
    uint64_t n = ((uint64_t)1) << 27;
    for (int i = 0; i < 40; i++) {
       test(4, 0, n);
       test(5, 0, n);
       test(6, 0, n);
       n *= 2;
    }
}
