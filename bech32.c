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

static uint32_t bech32_polymod_rstep(uint32_t pre, uint32_t encval, uint32_t decval) {
    int b = (encval ^ pre) & 0x1F;
    return (decval << 25) ^ (pre >> 5) ^
        (-((b >> 0) & 1) & 0x228bf1a8UL) ^
        (-((b >> 1) & 1) & 0x1703c750UL) ^
        (-((b >> 2) & 1) & 0x2c972fa9UL) ^
        (-((b >> 3) & 1) & 0xb2e5a72UL) ^
        (-((b >> 4) & 1) & 0x14d895edUL);
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

int main(void) {
    uint32_t test = bech32_polymod_rstep(bech32_polymod_step(0x123456,3,7),7,3);
    printf("res=%x\n", test);
    char out[16];
    uint8_t data[20] = {1,2,3};
    const char* hrp = "bc";
    size_t hrp_len = 2, data_len;
    bech32_encode(out, hrp, hrp_len, data, 3);
    printf("%s\n", out);
    printf("decode = %i\n", bech32_decode(&hrp_len, data, &data_len, out, strlen(out)));
}
