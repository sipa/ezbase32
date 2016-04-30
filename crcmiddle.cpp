#include <thread>
#include <mutex>
#include <algorithm>

#include <map>
#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <vector>

#include <immintrin.h>

namespace {

#define CHECKSUMBITS 30

/* BCH codes over GF(2^5)
 */
uint32_t compute_checksum4(const uint8_t* data, int len, const uint32_t* tbl) {
    uint32_t l = 0;

    while (len > 6) {
        uint8_t c = *(data++);
        len--;
        l = ((l & 0x1FFFFFFULL) << 5) ^ (tbl[(l >> 25) ^ c] & 0x3FFFFFFF);
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
    return ((l ^ f) * 1337) >> 2;
}

#define LEN 64
#define MAXERR 31

struct CRCOutputs {
    uint32_t val[LEN][MAXERR];

    CRCOutputs(const uint32_t *tbl) {
        unsigned char data[LEN] = {0};
        uint32_t none = compute_checksum4(data, LEN, tbl);
        for (int pos = 0; pos < LEN; pos++) {
            for (int v = 0; v < MAXERR; v++) {
                data[pos] = v + 1;
                val[pos][v] = compute_checksum4(data, LEN, tbl) ^ none;
            }
            data[pos] = 0;
        }
    }
};

/*
const uint32_t tbl[32] =  {
0, 0x61628a6c, 0xe6c59eb1, 0xab1789c8,
0x0e976928, 0x9ba2b753, 0x9c263a1a, 0xa0eeaa56,
0xfdb9a429, 0xc80a03eb, 0x7dab6118, 0x10919891,
0xc0dfb88a, 0xa7e81655, 0xe24098bf, 0x21c175df,
0x502437b2, 0x10e91bc4, 0xfcfdca8a, 0x52937cce,
0xc8825e02, 0x894fd3eb, 0xc08845cf, 0x70db58b2,
0xf01e4f7e, 0x2f34712b, 0xee67a508, 0xf2714345,
0x0ba51b75, 0x0b0e90cb, 0xb1b6c2fd, 0xb823972c,
};
*/

// const uint32_t tbl[32] = {0x0, 0x2f5651d3, 0x33dc7f91, 0x1c8a2e42, 0xad3a5f5, 0x2585f426, 0x390fda64, 0x16598bb7, 0x15a74bdd, 0x3af11a0e, 0x267b344c, 0x92d659f, 0x1f74ee28, 0x3022bffb, 0x2ca891b9, 0x3fec06a, 0x283e4d6d, 0x7681cbe, 0x1be232fc, 0x34b4632f, 0x22ede898, 0xdbbb94b, 0x11319709, 0x3e67c6da, 0x3d9906b0, 0x12cf5763, 0xe457921, 0x211328f2, 0x374aa345, 0x181cf296, 0x496dcd4, 0x2bc08d07}; // F=(x^5 + x^4 + x^2 + x + 1) E=(e^2 + 10*e + 14) alpha=(18*e + 15) powers=3..6

// const uint32_t tbl[32] = {0x0, 0x33174241, 0x3ade3d62, 0x9c97f23, 0x2babfac4, 0x18bcb885, 0x1175c7a6, 0x226285e7, 0xbb0cc68, 0x38a78e29, 0x316ef10a, 0x279b34b, 0x201b36ac, 0x130c74ed, 0x1ac50bce, 0x29d2498f, 0x159124d0, 0x26866691, 0x2f4f19b2, 0x1c585bf3, 0x3e3ade14, 0xd2d9c55, 0x4e4e376, 0x37f3a137, 0x1e21e8b8, 0x2d36aaf9, 0x24ffd5da, 0x17e8979b, 0x358a127c, 0x69d503d, 0xf542f1e, 0x3c436d5f}; // F=(x^5 + x^3 + x^2 + x + 1) E=(e^3 + 10*e^2 + 4*e + 3) alpha=(7*e^2 + 31*e + 3) powers=1..2

// BCH_0077
const uint32_t tbl[32] =  {0x0, 0x2f5651d3, 0x33dc7f91, 0x1c8a2e42, 0xad3a5f5, 0x2585f426, 0x390fda64, 0x16598bb7, 0x15a74bdd, 0x3af11a0e, 0x267b344c, 0x92d659f, 0x1f74ee28, 0x3022bffb, 0x2ca891b9, 0x3fec06a, 0x283e4d6d, 0x7681cbe, 0x1be232fc, 0x34b4632f, 0x22ede898, 0xdbbb94b, 0x11319709, 0x3e67c6da, 0x3d9906b0, 0x12cf5763, 0xe457921, 0x211328f2, 0x374aa345, 0x181cf296, 0x496dcd4, 0x2bc08d07}; // F=(x^5 + x^4 + x^2 + x + 1) E=(e^2 + 10*e + 14) alpha=(18*e + 15) powers=3..6

class IncMap {
    // Encoding: packed per-byte serialized:
    // * 00000000: end
    // * 00000001 - 01111111: value N
    // * 10000000 - 10111111: skip 1-64
    // * 11000000 - 11011111 (+ 1 byte): skip 65-8256
    // * 11100000 - 11101111 (+ 2 bytes): skip 8257-1056832
    // * 11110000 - 11110111 (+ 3 bytes): skip 1056833-135274560
    // * 11111000 (+ 4 bytes): 135274561-...

    std::vector<uint8_t> data;
    uint32_t lastkey;

    IncMap(const IncMap& x) = delete;
    IncMap& operator=(const IncMap& x) = delete;

    void WriteSkip(uint32_t skip) {
        if (skip == 0) {
            return;
        } else if (skip <= 64) {
            skip -= 1;
            data.push_back(0x80 | skip);
        } else if (skip <= 8256) {
            skip -= 65;
            data.push_back(0xC0 | (skip >> 8));
            data.push_back(skip & 0xFF);
        } else if (skip <= 1056832) {
            skip -= 8257;
            data.push_back(0xE0 | (skip >> 16));
            data.push_back((skip >> 8) & 0xFF);
            data.push_back(skip & 0xFF);
        } else if (skip <= 135274560) {
            skip -= 1056833;
            data.push_back(0xF0 | (skip >> 24));
            data.push_back((skip >> 16) & 0xFF);
            data.push_back((skip >> 8) & 0xFF);
            data.push_back(skip & 0xFF);
        } else {
            skip -= 135274561;
            data.push_back(0xF8);
            data.push_back((skip >> 24) & 0xFF);
            data.push_back((skip >> 16) & 0xFF);
            data.push_back((skip >> 8) & 0xFF);
            data.push_back(skip & 0xFF);
        }
    }

    void AppendOrdered(uint32_t key, uint8_t value) {
        WriteSkip(key - lastkey);
        lastkey = key + 1;
        data.push_back(value);
    }

public:

    IncMap() : data{0x00}, lastkey(0) {}
    IncMap(IncMap&& m) noexcept { data = std::move(m.data); lastkey = m.lastkey; }

    IncMap& operator=(IncMap&& m) noexcept { data = std::move(m.data); lastkey = m.lastkey; return *this; }

    class iterator {
        uint32_t key;
        std::vector<uint8_t>::const_iterator inner;

        void Next() {
            if ((*inner & 0x80) == 0) {
                return;
            } else if ((*inner & 0xC0) == 0x80) {
                key += 1 + (*inner & 0x3F);
                inner++;
            } else if ((*inner & 0xE0) == 0xC0) {
                key += 65 + (((uint32_t)(*inner & 0x1F)) << 8) + inner[1];
                inner += 2;
            } else if ((*inner & 0xF0) == 0xE0) {
                key += 8257 + (((uint32_t)(*inner & 0x0F)) << 16) + (((uint32_t)inner[1]) << 8) + inner[2];
                inner += 3;
            } else if ((*inner & 0xF8) == 0xF0) {
                key += 1056833 + (((uint32_t)(*inner & 0x0F)) << 24) + (((uint32_t)inner[1]) << 16) + (((uint32_t)inner[2]) << 8) + inner[3];
                inner += 4;
            } else {
                key += 135274560 + (((uint32_t)inner[1]) << 24) + (((uint32_t)inner[2]) << 16) + (((uint32_t)inner[3]) << 8) + inner[4];
                inner += 5;
            }
        }

        iterator(uint32_t key_, std::vector<uint8_t>::const_iterator inner_) : key(key_), inner(inner_) {
            Next();
        }

    public:
        uint32_t GetKey() {
            return key;
        }
        uint8_t GetValue() {
            return *inner;
        }
        bool Valid() {
            return *inner != 0x00;
        }
        void Increment() {
            inner++;
            key++;
            Next();
        }

        friend class IncMap;
    };

    iterator begin() const { return iterator((uint32_t)0, data.begin()); }

    IncMap(const IncMap& old, const std::vector<uint32_t>& sortednew) : lastkey(0) {
        iterator oldit = old.begin();
        std::vector<uint32_t>::const_iterator newit = sortednew.begin();
        uint32_t key = 0;
        data.reserve(old.data.size() + 3 * sortednew.size());
        while(true) {
            uint8_t value = 0;
            if (oldit.Valid() && oldit.GetKey() == key) {
                value += oldit.GetValue();
                oldit.Increment();
            }
            while (newit != sortednew.end() && *newit == key) {
                value++;
                newit++;
            }
            if (value > 0) {
                if (value > 0x7F) {
                    abort();
                }
//                printf("  * Inserting (%i,%i)\n", (int)key, (int)value);
                AppendOrdered(key, value);
            }
            if (newit == sortednew.end()) {
                if (!oldit.Valid()) {
                    break;
                } else {
                    key = oldit.GetKey();
                }
            } else {
                if (!oldit.Valid()) {
                    key = *newit;
                } else {
                    key = std::min(oldit.GetKey(), *newit);
                }
            }
        }
        data.push_back(0);
    }

protected:
    size_t Memusage() const {
        return data.size();
    }
};

class MutableIncMap : public IncMap {
    std::vector<uint32_t> unsorted;
    uint64_t total;

    void Merge() {
        std::sort(unsorted.begin(), unsorted.end());
        *(static_cast<IncMap*>(this)) = IncMap(*this, unsorted);
        unsorted.clear();
    }

public:
    MutableIncMap() : total(0) {}

    void Increment(uint32_t key) {
        unsorted.push_back(key);
        if (unsorted.size() * 4 > Memusage() + 65536) {
            Merge();
        }
        total++;
    }

    void Prepare() {
        if (unsorted.size() > 0) {
            Merge();
        }
    }

    uint64_t GetTotal() const {
        return total;
    }
};

const CRCOutputs outputs(tbl);

void Recurse(uint32_t accum, int errors, int minpos, int endpos, MutableIncMap& data) {
    if (errors == 0) {
        data.Increment(accum);
        if ((data.GetTotal() & 0xFFFFFF) == 0) {
            printf(".");
        }
        return;
    }
    for (int pos = minpos; pos < endpos - errors + 1; pos++) {
        for (int err = 1; err <= MAXERR; err++) {
            Recurse(accum ^ outputs.val[pos][err - 1], errors - 1, pos + 1, endpos, data);
        }
    }
}

void BuildEndMap(int errors, int endpos, MutableIncMap& data) {
    for (int err = 1; err <= MAXERR; err++) {
        Recurse(outputs.val[endpos][err - 1], errors - 1, 0, endpos, data);
    }
    data.Prepare();
}

void BuildBeginMap(int errors, int beginpos, MutableIncMap& data) {
    for (int err = 1; err <= MAXERR; err++) {
        Recurse(outputs.val[beginpos][err - 1], errors - 1, beginpos + 1, LEN, data);
    }
    data.Prepare();
}

struct Maps {
    std::vector<MutableIncMap> beginmaps;
    std::vector<MutableIncMap> endmaps;

    Maps(int errors) {
        uint64_t count = 0;
        printf("Building maps for %i errors...", errors);
        beginmaps.resize(LEN);
        endmaps.resize(LEN);
        for (int i = 0; i < LEN; i++) {
            BuildBeginMap(errors, i, beginmaps[i]);
            count += beginmaps[i].GetTotal();
            BuildEndMap(errors, i, endmaps[i]);
            count += endmaps[i].GetTotal();
        }
        printf("done (%llu combinations)\n", (unsigned long long)count);
    }

    uint64_t Failures() {
        uint64_t ret1 = 0, ret2 = 0;
        for (int i = 0; i < LEN; i++) {
            IncMap::iterator it1 = beginmaps[i].begin();
            if (it1.Valid() && it1.GetKey()==0) ret1 += it1.GetValue();
            IncMap::iterator it2 = endmaps[i].begin();
            if (it2.Valid() && it2.GetKey()==0) ret2 += it2.GetValue();
        }
        assert(ret1 == ret2);
        return ret1;
    }

    uint64_t FailuresCombined(const Maps& other) {
        uint64_t ret = 0;
        for (int i = 1; i < LEN - 1; i++) {
            const IncMap& endmap = endmaps[i];
            for (int j = i + 1; j < LEN; j++) {
                const IncMap& beginmap = beginmaps[j];
                IncMap::iterator eit = endmap.begin();
                IncMap::iterator bit = beginmap.begin();
                while (eit.Valid() && bit.Valid()) {
                    while (eit.Valid() && bit.Valid() && bit.GetKey() != eit.GetKey()) {
                        while (eit.Valid() && eit.GetKey() < bit.GetKey()) eit.Increment();
                        while (bit.Valid() && bit.GetKey() < eit.GetKey()) bit.Increment();
                    }
                    if (!bit.Valid() || !eit.Valid()) break;
                    assert(eit.GetKey() == bit.GetKey());
                    ret += ((uint64_t)bit.GetValue()) * eit.GetValue();
                    bit.Increment();
                    eit.Increment();
                }
            }
        }
        return ret;
    }
};

}

int main(void) {
    std::vector<Maps> list;
    setbuf(stdout, NULL);

/*
    std::map<uint32_t, uint16_t> test;
    MutableIncMap tes;
    for (int i = 0; i < 1000000; i++) {
        int r = random() % 100000;
        tes.Increment(r);
        test[r]++;
        printf("Adding %i\n", r);
        if (i % 2 == 0) {
            tes.Prepare();
            MutableIncMap::iterator it = tes.begin();
            while (it.Valid()) {
                printf("- Testing: (%i,%i), real: %i\n", (int)it.GetKey(), (int)it.GetValue(), (int)test[it.GetKey()]);
                assert(test[it.GetKey()] == it.GetValue());
                it.Increment();
            }
        }
    }
*/

    for (int i = 1; i <= 3; i++) {
        list.push_back(Maps(i));
        printf("Undetected HD%i = %llu\n", i, (unsigned long long)list[i - 1].Failures());
        for (int j = 1; j <= i; j++) {
            printf("Undetected HD%i = %llu\n", i + j, (unsigned long long)list[i - 1].FailuresCombined(list[j - 1]));
        }
    }
/*
        for (int j = 1; j <= i; j++) {
            uint64_t comb = GetCombinedError(list[j-1], j, list[i-1], i);
            printf("Undetected HD%i", i+j);
            for (int x = i+j-2; x >= 1; x -= 2) {
                printf(" + HD%i", x);
            }
            printf(" <= %llu\n", (unsigned long long)comb);
        }
*/
    return 0;
}
