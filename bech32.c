#include <stdio.h>
#include <stdint.h>
#include <string.h>

static uint32_t bech32_polymod_step(uint32_t pre, uint32_t encval, uint32_t decval) {
    int b = encval ^ (pre >> 25);
    return decval ^ ((pre & 0x1FFFFFF) << 5) ^
        (-((b >> 0) & 1) & 0x3b6a57b2UL) ^
        (-((b >> 1) & 1) & 0x26508e6dUL) ^
        (-((b >> 2) & 1) & 0x1ea119faUL) ^
        (-((b >> 3) & 1) & 0x3d4233ddUL) ^
        (-((b >> 4) & 1) & 0x2a1462b3UL);
}

void bech32_encode(char* output, const char* hrp, size_t hrp_len, const uint8_t* data, size_t data_len) {
    static const char* zbase32="ybndrfg8ejkmcpqxot1uwisza345h769";
    uint32_t chk = 0x3b6a57b2UL;
    size_t i;
    for (i = 0; i < hrp_len; ++i) {
        chk = bech32_polymod_step(chk, hrp[i] >> 5, 0);
    }
    chk = bech32_polymod_step(chk, 0, 0);
    for (i = 0; i < hrp_len; ++i) {
        chk = bech32_polymod_step(chk, *hrp & 0x1f, 0);
        *(output++) = *(hrp++);
    }
    *(output++) = '-';
    chk = bech32_polymod_step(chk, 0, 0);
    for (i = 0; i < data_len; ++i) {
        chk = bech32_polymod_step(chk, *data, 0);
        *(output++) = zbase32[*(data++)];
    }
    chk ^= 1;
    for (i = 0; i < 6; ++i) {
        *(output++) = zbase32[(chk >> ((5 - i) * 5)) & 0x1f];
    }
    *output = 0;
}

int bech32_decode(size_t *hrp_len, uint8_t* data, size_t* data_len, char* input, size_t input_len) {
    static const int zbase32_alpha[26] = {24,1,12,3,8,5,6,28,21,9,10,-1,11,2,16,13,14,4,22,17,19,-1,20,15,0,23};
    static const int zbase32_num[10] = {-1,18,-1,25,26,27,30,29,7,31};
    uint32_t chk = 1;
    size_t i;
    if (input_len < 8 || input_len > 89) {
        return -1;
    }
    *data_len = 0;
    while (*data_len < input_len && input[input_len - 1 - *data_len] != '-') {
        ++(*data_len);
    }
    *hrp_len = input_len - *data_len;
    if (*hrp_len < 1) {
        return -2;
    }
    for (i = 0; i < *hrp_len - 1; ++i) {
        if (input[i] < 32 || input[i] > 127) {
            return -3;
        }
        chk = bech32_polymod_step(chk, 0, input[i] >> 5);
    }
    chk = bech32_polymod_step(chk, 0, 0);
    for (i = 0; i < *hrp_len - 1; ++i) {
        chk = bech32_polymod_step(chk, 0, input[i] & 0x1f);
    }
    chk = bech32_polymod_step(chk, 0, 0);
    ++i;
    while (i < input_len) {
        int v = (input[i] >= '0' && input[i] <= '9') ? zbase32_num[input[i] - '0'] :
                (input[i] >= 'a' && input[i] <= 'z') ? zbase32_alpha[input[i] - 'a'] :
                -1;
        if (v == -1) {
            return -4;
        }
        chk = bech32_polymod_step(chk, 0, v);
        if (i + 6 < input_len) {
            data[i - *hrp_len - 1] = v;
        }
        ++i;
    }
    return chk == 1;
}

static int16_t exp10[1024], log10[1024];

void logtable10(void) {
    // Build table for GF(32)
    int8_t exp5[32], log5[32];
    int fmod = 41;
    log5[0] = -1;
    log5[1] = 0;
    exp5[0] = 1;
    exp5[31] = 1;
    int v = 1;
    for (int i = 1; i < 31; i++) {
        v = v << 1;
        if (v & 32) v ^= fmod;
        exp5[i] = v;
        log5[v] = i;
    }

    // Build table for GF(1024)
    log10[0] = -1;
    log10[1] = 0;
    exp10[0] = 1;
    exp10[1023] = 1;
    v = 1;
    for (int i = 1; i < 1023; i++) {
        int v0 = v & 31;
        int v1 = v >> 5;

        int v1n = (v1 ? exp5[(log5[v1] + log5[6]) % 31] : 0) ^ (v0 ? exp5[(log5[v0] + log5[9]) % 31] : 0);
        int v0n = (v1 ? exp5[(log5[v1] + log5[27]) % 31] : 0) ^ (v0 ? exp5[(log5[v0] + log5[15]) % 31] : 0);
        v = v1n << 5 | v0n;
        exp10[i] = v;
        log10[v] = i;
    }
}

static inline uint32_t syndrome(uint32_t fault) {
    uint32_t low = fault & 0x1f;
    return low ^ (low << 10) ^ (low << 20) ^
           (-((fault >> 5) & 1) & 0x31edd3c4UL) ^
           (-((fault >> 6) & 1) & 0x335f86a8UL) ^
           (-((fault >> 7) & 1) & 0x363b8870UL) ^
           (-((fault >> 8) & 1) & 0x3e6390c9UL) ^
           (-((fault >> 9) & 1) & 0x2ec72192UL) ^
           (-((fault >> 10) & 1) & 0x1046f79dUL) ^
           (-((fault >> 11) & 1) & 0x208d4e33UL) ^
           (-((fault >> 12) & 1) & 0x130ebd6fUL) ^
           (-((fault >> 13) & 1) & 0x2499fadeUL) ^
           (-((fault >> 14) & 1) & 0x1b27d4b5UL) ^
           (-((fault >> 15) & 1) & 0x04be1eb4UL) ^
           (-((fault >> 16) & 1) & 0x0968b861UL) ^
           (-((fault >> 17) & 1) & 0x1055f0c2UL) ^
           (-((fault >> 18) & 1) & 0x20ab4584UL) ^
           (-((fault >> 19) & 1) & 0x1342af08UL) ^
           (-((fault >> 20) & 1) & 0x24f1f318UL) ^
           (-((fault >> 21) & 1) & 0x1be34739UL) ^
           (-((fault >> 22) & 1) & 0x35562f7bUL) ^
           (-((fault >> 23) & 1) & 0x3a3c5bffUL) ^
           (-((fault >> 24) & 1) & 0x266c96f7UL) ^
           (-((fault >> 25) & 1) & 0x25c78b65UL) ^
           (-((fault >> 26) & 1) & 0x1b1f13eaUL) ^
           (-((fault >> 27) & 1) & 0x34baa2f4UL) ^
           (-((fault >> 28) & 1) & 0x3b61c0e1UL) ^
           (-((fault >> 29) & 1) & 0x265325c2UL);
}

static inline uint32_t mod1023(uint32_t x) { return (x & 0x3FF) + (x >> 10); }

int locate_errors2(uint32_t fault, int length, int maxnon1bit)
{
    if (fault == 0) {
        return 0;
    }

    uint32_t syn = syndrome(fault);
    int s0 = syn & 0x3FF;
    int s1 = (syn >> 10) & 0x3FF;
    int s2 = syn >> 20;

    int l_s0 = log10[s0], l_s1 = log10[s1], l_s2 = log10[s2];
    if (l_s0 != -1 && l_s1 != -1 && l_s2 != -1 && mod1023(mod1023(2 * l_s1 - l_s2 - l_s0 + 2047)) == 1) {
        int p1 = mod1023(l_s1 - l_s0 + 1024) - 1;
        if (p1 >= length) return -1;
        int l_e1 = mod1023(l_s0 + (1023 - 997) * p1 + 1) - 1;
        if (maxnon1bit == 0 && (l_e1 > 396 || (l_e1 % 99))) return -1;
        if (l_e1 % 33) return -1;
        return p1 + 1;
    }

    for (int p1 = 0; p1 < length; p1++) {
        int s2_s1p1 = s2 ^ (s1 == 0 ? 0 : exp10[mod1023(l_s1 + p1)]);
        if (s2_s1p1 == 0) continue;
        int s1_s0p1 = s1 ^ (s0 == 0 ? 0 : exp10[mod1023(l_s0 + p1)]);
        if (s1_s0p1 == 0) continue;
        int l_s1_s0p1 = log10[s1_s0p1];

        int p2 = mod1023(log10[s2_s1p1] - l_s1_s0p1 + 1023);
        if (p2 >= length || p1 == p2) continue;

        int s1_s0p2 = s1 ^ (s0 == 0 ? 0 : exp10[mod1023(l_s0 + p2)]);
        if (s1_s0p2 == 0) continue;

        int inv_p1_p2 = 1023 - log10[exp10[p1] ^ exp10[p2]];

        int l_e2 = mod1023(l_s1_s0p1 + inv_p1_p2 + (1023 - 997) * p2 + 1) - 1;
        if (l_e2 % 33) continue;
        int non1bit = l_e2 > 396 || (l_e2 % 99);
        int l_e1 = mod1023(log10[s1_s0p2] + inv_p1_p2 + (1023 - 997) * p1 + 1) - 1;
        if (l_e1 % 33) continue;
        non1bit += l_e1 > 396 || (l_e1 % 99);
        if (non1bit > maxnon1bit) return -1;

        if (p1 < p2) {
            return (p1 + 1) << 8 | (p2 + 1);
        } else {
            return (p2 + 1) << 8 | (p1 + 1);
        }
    }
    return -1;
}

#define min(a,b) ((a)<(b)?(a):(b))
#define max(a,b) ((a)>(b)?(a):(b))

int locate_errors3(uint32_t fault, int length)
{
    int result = locate_errors2(fault, length, 0);
    if (result != -1) {
        return result;
    }
    for (int non1bit = 0; non1bit < 4; ++non1bit) {
        int found = 0;
        for (int i_e1 = 0; i_e1 < (non1bit == 3 ? 31 : 5); i_e1++) {
            int e1 = 1 << i_e1;
            if (non1bit == 3) {
                e1 = 1 + i_e1;
            }
            uint32_t fault_e1 = e1;
            for (int p1 = 0; p1 < length; ++p1) {
                int res = locate_errors2(fault ^ fault_e1, length, min(2, non1bit));
//                printf("Try: %i@%i (fault=%x) => fault=%x => res=%x\n", e1, p1, fault_e1, fault ^ fault_e1, res);
                fault_e1 = bech32_polymod_step(fault_e1, 0, 0);
                if (res < 0x100) continue;

                int p2 = (res & 0xff) - 1;
                int p3 = (res >> 8) - 1;
                int ps[3];
                ps[0] = min(min(p1, p2), p3);
                ps[2] = max(max(p1, p2), p3);
                ps[1] = p1 + p2 + p3 - ps[0] - ps[2];
//                printf("* Res: %i %i %i\n", ps[0], ps[1], ps[2]);
                int comb = ((1 + ps[0]) << 16) | ((1 + ps[1]) << 8) | (1 + ps[2]);
                if (found == 0 || comb != result) {
                    found++;
                    result = comb;
                }
            }
        }
        if (found == 1) return result;
    }
    return -1;
}

int main(void) {
    logtable10();

    for (int i = 0; i < 5; ++i) {
        printf("log(%i) = %i\n", i, log10[1 << i]);
    }

/*
    Generate the fault -> syndrome table:

    for (int p = 0; p < 6; p++) {
        for (int i = 0; i < 5; i++) {
            int fault = 1 << i;
            for (int x = 0; x < p; x++) { fault = bech32_polymod_step(fault, 0, 0); }
            int s0 = 0, s1 = 0, s2 = 0;
            for (int i = 0; i < 6; i++) {
                int fi = (fault >> (5 * (5 - i))) & 31;
                s0 = mul10l(s0, 997) ^ fi;
                s1 = mul10l(s1, 998) ^ fi;
                s2 = mul10l(s2, 999) ^ fi;
            }
            printf("0x%x ", s2 << 20 | s1 << 10 | s0);
        }
    }
    printf("\n");
*/

    static uint32_t faults[89][31];

    for (int err = 1; err < 32; ++err) {
        faults[0][err - 1] = err;
        for (int pos = 1; pos < 89; ++pos) {
            faults[pos][err - 1] = bech32_polymod_step(faults[pos - 1][err - 1], 0, 0);
        }
    }


/*
    for (int pos1 = 0; pos1 < 89; ++pos1) {
        for (int err1 = 1; err1 < 32; ++err1) {
            uint32_t fault = faults[pos1][err1 - 1];
            int solve = locate_errors2(fault, 89, 0);
            if (solve != pos1 + 1) {
                printf("Fail: E%iP%i -> S%x\n", err1, pos1, solve);
            }
        }
    }

    for (int pos1 = 0; pos1 < 89; ++pos1) {
        for (int err1 = 1; err1 < 32; ++err1) {
            uint32_t fault1 = faults[pos1][err1 - 1];
            for (int pos2 = pos1 + 1; pos2 < 89; ++pos2) {
                for (int err2 = 1; err2 < 32; ++err2) {
                    uint32_t fault = fault1 ^ faults[pos2][err2 - 1];
                    int solve = locate_errors2(fault, 89, 0);
                    if (solve != (((pos1 + 1) << 8) | (pos2 + 1))) {
                        printf("Fail: E%i@P%i E%iP%i -> S%x\n", err1, pos1, err2, pos2, solve);
                    }
                }
            }
        }
    }
*/


    int len = 10;

    uint64_t iter[4][31*31*31] = {0};
    uint64_t fail[4][31*31*31] = {0};
    uint64_t unk[4][31*31*31] = {0};
    #pragma omp parallel for
    for (int es = 0; es < 31 * 31 * 31; es++) {
        int e1 = 1 + ((es) % 31);
        int e2 = 1 + ((es / 31) % 31);
        int e3 = 1 + ((es / 31 / 31) % 31);
        int non1bit = (!!(e1 & (e1-1))) + (!!(e2 & (e2-1))) + (!!(e3 & (e3-1)));
        for (int p1 = 0; p1 < len - 2; p1++) {
            for (int p2 = p1 + 1; p2 < len - 1; p2++) {
                for (int p3 = p2 + 1; p3 < len; p3++) {
                    int res = locate_errors3(faults[p1][e1 - 1] ^ faults[p2][e2 - 1] ^ faults[p3][e3 - 1], len);
                    if (res == -1) {
                        ++unk[non1bit][es];
                    } else {
                        int exp = ((1 + p1) << 16) | ((1 + p2) << 8) | (1 + p3);
                        fail[non1bit][es] += (res != exp);
                    }
                    ++iter[non1bit][es];
                }
            }
        }
    }

     for (int non1bit = 0; non1bit < 4; ++non1bit) {
         uint64_t iters = 0;
         uint64_t fails = 0;
         uint64_t unks = 0;
         for (int es = 0; es < 31*31*31; ++es) {
             iters += iter[non1bit][es];
             fails += fail[non1bit][es];
             unks += unk[non1bit][es];
         }
         printf("%i non-1-bit: Unk=%llu Fails=%llu Total=%llu\n", non1bit, (unsigned long long)unks, (unsigned long long)fails, (unsigned long long)iters);
     }
}
