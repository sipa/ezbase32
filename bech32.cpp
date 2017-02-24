#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <assert.h>

#include <map>
#include <algorithm>

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

static int16_t exptable[1024], logtable[1024];

void build_gftables(void) {
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
    logtable[0] = -1;
    logtable[1] = 0;
    exptable[0] = 1;
    exptable[1023] = 1;
    v = 1;
    for (int i = 1; i < 1023; i++) {
        int v0 = v & 31;
        int v1 = v >> 5;

        int v1n = (v1 ? exp5[(log5[v1] + log5[6]) % 31] : 0) ^ (v0 ? exp5[(log5[v0] + log5[9]) % 31] : 0);
        int v0n = (v1 ? exp5[(log5[v1] + log5[27]) % 31] : 0) ^ (v0 ? exp5[(log5[v0] + log5[15]) % 31] : 0);
        v = v1n << 5 | v0n;
        exptable[i] = v;
        logtable[v] = i;
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


int find_error_pos(uint32_t fault, int length)
{
    if (fault == 0) {
        return 0;
    }

    uint32_t syn = syndrome(fault);
    int s0 = syn & 0x3FF;
    int s1 = (syn >> 10) & 0x3FF;
    int s2 = syn >> 20;

    int l_s0 = logtable[s0], l_s1 = logtable[s1], l_s2 = logtable[s2];
    if (l_s0 != -1 && l_s1 != -1 && l_s2 != -1 && mod1023(mod1023(2 * l_s1 - l_s2 - l_s0 + 2047)) == 1) {
        int p1 = mod1023(l_s1 - l_s0 + 1024) - 1;
        if (p1 >= length) return -1;
        int e1 = exptable[mod1023(mod1023(l_s0 + (1023 - 997) * p1))];
        if (e1 >= 32) return -1;
        return p1 + 1;
    }

    for (int p1 = 0; p1 < length; p1++) {
        int s2_s1p1 = s2 ^ (s1 == 0 ? 0 : exptable[mod1023(l_s1 + p1)]);
        if (s2_s1p1 == 0) continue;
        int s1_s0p1 = s1 ^ (s0 == 0 ? 0 : exptable[mod1023(l_s0 + p1)]);
        if (s1_s0p1 == 0) continue;
        int l_s1_s0p1 = logtable[s1_s0p1];

        int p2 = mod1023(logtable[s2_s1p1] - l_s1_s0p1 + 1023);
        if (p2 >= length || p1 == p2) continue;

        int s1_s0p2 = s1 ^ (s0 == 0 ? 0 : exptable[mod1023(l_s0 + p2)]);
        if (s1_s0p2 == 0) continue;

        int inv_p1_p2 = 1023 - logtable[exptable[p1] ^ exptable[p2]];

        int e1 = exptable[mod1023(mod1023(logtable[s1_s0p2] + inv_p1_p2 + (1023 - 997)*p1))];
        if (e1 >= 32) continue;

        int e2 = exptable[mod1023(mod1023(l_s1_s0p1 + inv_p1_p2 + (1023 - 997)*p2))];
        if (e2 >= 32) continue;

        if (p1 < p2) {
            return (p1 + 1) << 8 | (p2 + 1);
        } else {
            return (p2 + 1) << 8 | (p1 + 1);
        }
    }
    return -1;
}

static int cmp3(const void* a, const void* b) {
    int* pa = (int*)a;
    int* pb = (int*)b;
    if (pa[0] < pb[0]) return -1;
    if (pa[0] > pb[0]) return 1;
    if (pa[1] < pb[1]) return -1;
    if (pa[1] > pb[1]) return 1;
    if (pa[2] < pb[2]) return -1;
    if (pa[2] > pb[2]) return 1;
    return 0;
}


int main(void) {
    build_gftables();

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

    for (int pos1 = 0; pos1 < 89; ++pos1) {
        for (int err1 = 1; err1 < 32; ++err1) {
            uint32_t fault = faults[pos1][err1 - 1];
            int solve = find_error_pos(fault, 89);
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
                    int solve = find_error_pos(fault, 89);
                    if (solve != (((pos1 + 1) << 8) | (pos2 + 1))) {
                        printf("Fail: E%i@P%i E%iP%i -> S%x\n", err1, pos1, err2, pos2, solve);
                    }
                }
            }
        }
    }

    int len = 10;
    uint64_t total[32] = {0};
    uint64_t fails[32] = {0};
    std::map<int, uint64_t> xfails[32];
    #pragma omp parallel for
    for (int err1 = 1; err1 < 32; ++err1) {
        int counts[89*31][3];
        for (int pos1 = 0; pos1 < len; ++pos1) {
            uint32_t fault1 = faults[pos1][err1 - 1];
            for (int pos2 = pos1 + 1; pos2 < len; ++pos2) {
                for (int err2 = 1; err2 < 32; ++err2) {
                    uint32_t fault2 = fault1 ^ faults[pos2][err2 - 1];
                    for (int pos3 = pos2 + 1; pos3 < len; ++pos3) {
                        for (int err3 = 1; err3 < 32; ++err3) {
                            uint32_t fault3 = fault2 ^ faults[pos3][err3 - 1];
                            int solve = find_error_pos(fault3, len);
                            ++total[err1];
                            if (solve != -1) {
                                fails[err1] += (solve != -1);
                            } else {
                                int ncounts = 0;
                                for (int pos4 = 0; pos4 < len; ++pos4) {
                                    for (int err4 = 1; err4 < 32; ++err4) {
                                        uint32_t fault4 = fault3 ^ faults[pos4][err4 - 1];
                                        int solvex = find_error_pos(fault4, len);
                                        if (solvex != -1) {
                                            unsigned int p1 = (solvex & 0xFF) - 1;
                                            unsigned int p2 = (solvex >> 8) - 1;
                                            assert(p1 < len && p2 < len);
                                            std::array<int, 3> pos = {pos4, p1, p2};
                                            std::sort(pos.begin(), pos.end());
                                            assert(ncounts < 89*31);
                                            counts[ncounts][0] = pos[0];
                                            counts[ncounts][1] = pos[1];
                                            counts[ncounts][2] = pos[2];
                                            ++ncounts;
                                        }
                                    }
                                }
                                qsort(counts, ncounts, sizeof(counts[0]), cmp3);
                                int diff = 0;
                                for (int xx = 0; xx < ncounts; ++xx) {
                                    if (xx == 0 || cmp3(counts[xx], counts[xx-1]) != 0) {
                                        ++diff;
                                    }
                                }
                                ++xfails[err1][diff];
                            }
                        }
                    }
                }
            }
        }
    }
    uint64_t o_total = 0;
    uint64_t o_fails = 0;
    std::map<int, int> o_xfails;
    for (int p = 0; p < 32; ++p) {
        o_total += total[p];
        o_fails += fails[p];
        for (auto const &e : xfails[p]) {
            o_xfails[e.first] += e.second;
        }
    }
    printf("%llu out of %llu HD3 errors are HD2 from another valid codeword\n", (unsigned long long)o_fails, (unsigned long long)o_total);
    for (auto const &e : o_xfails) {
        printf("%i HD3 interpretations occur %llu times\n", e.first, (unsigned long long)e.second);
    }
}
