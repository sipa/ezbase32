#include <stdio.h>
#include <stdint.h>

static uint32_t bech32_polymod_step(uint32_t pre, uint32_t encval, uint32_t decval) {
    int b = encval ^ (pre >> 25);
    return decval ^ ((pre & 0x1FFFFFF) << 5) ^
        (-((b >> 0) & 1) & 0x53A0C81UL) ^
        (-((b >> 1) & 1) & 0x8F09902UL) ^
        (-((b >> 2) & 1) & 0x11E13204UL) ^
        (-((b >> 3) & 1) & 0x21526128UL) ^
        (-((b >> 4) & 1) & 0x12346650UL);
}

void bech32_encode(char* output, const char* hrp, size_t hrp_len, const uint8_t* data, size_t data_len) {
    static const char* zbase32="ybndrfg8ejkmcpqxot1uwisza345h769";
    uint32_t chk = 0x53A0C81UL;
    for (size_t i = 0; i < hrp_len; ++i) {
        chk = bech32_polymod_step(chk, hrp[i] >> 5, 0);
    }
    chk = bech32_polymod_step(chk, 0, 0);
    for (size_t i = 0; i < hrp_len; ++i) {
        chk = bech32_polymod_step(chk, *hrp & 0x1f, 0);
        *(output++) = *(hrp++);
    }
    *(output++) = '-';
    chk = bech32_polymod_step(chk, 0, 0);
    for (size_t i = 0; i < data_len; ++i) {
        chk = bech32_polymod_step(chk, *data, 0);
        *(output++) = ZBASE32[*(data++)];
    }
    chk ^= 1;
    for (size_t i = 0; i < 6; ++i) {
        *(output++) = ZBASE32[(chk >> ((5 - i) * 5)) & 0x1f];
    }
    *output = 0;
}

int main(void) {
    char out[16];
    uint8_t data[3] = {1,2,3};
    const char* hrp = "bc";
    bech32_encode(out, hrp, 2, data, sizeof(data));
    printf("%s\n", out);
}
