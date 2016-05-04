#include <thread>
#include <mutex>
#include <algorithm>

#include <math.h>
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
uint32_t compute_bch(const uint8_t* data, int len, const uint32_t* tbl) {
    uint32_t l = 0;

    while (len > 6) {
        uint8_t c = *(data++);
        len--;
        l = ((l & 0x1FFFFFFULL) << 5) ^ tbl[(l >> 25) ^ c];
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

#define LEN 64
#define MAXERR 31

struct CRCOutputs {
    uint32_t val[LEN][MAXERR];

    CRCOutputs(const uint32_t *tbl) {
        unsigned char data[LEN] = {0};
        uint32_t none = compute_bch(data, LEN, tbl);
        for (int pos = 0; pos < LEN; pos++) {
            for (int v = 0; v < MAXERR; v++) {
                data[pos] = v + 1;
                val[pos][v] = compute_bch(data, LEN, tbl) ^ none;
            }
            data[pos] = 0;
        }
    }
};

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

uint64_t SimpleRecurse(uint32_t accum, int errors, int minpos, const CRCOutputs& outputs) {
    if (errors == 0) {
        return accum == 0;
    }
    uint64_t ret = 0;
    for (int pos = minpos; pos < LEN; pos++) {
        for (int err = 1; err <= MAXERR; err++) {
            ret += SimpleRecurse(accum ^ outputs.val[pos][err - 1], errors - 1, pos + 1, outputs);
        }
    }
    return ret;
}

void Recurse(uint32_t accum, int errors, int minpos, int endpos, MutableIncMap& data, const CRCOutputs& outputs) {
    if (errors == 0) {
        data.Increment(accum);
        if ((data.GetTotal() & 0xFFFFFF) == 0) {
            printf(".");
        }
        return;
    }
    for (int pos = minpos; pos < endpos - errors + 1; pos++) {
        for (int err = 1; err <= MAXERR; err++) {
            Recurse(accum ^ outputs.val[pos][err - 1], errors - 1, pos + 1, endpos, data, outputs);
        }
    }
}

void BuildEndMap(int errors, int endpos, MutableIncMap& data, const CRCOutputs& outputs) {
    for (int err = 1; err <= MAXERR; err++) {
        Recurse(outputs.val[endpos][err - 1], errors - 1, 0, endpos, data, outputs);
    }
    data.Prepare();
}

void BuildBeginMap(int errors, int beginpos, MutableIncMap& data, const CRCOutputs& outputs) {
    for (int err = 1; err <= MAXERR; err++) {
        Recurse(outputs.val[beginpos][err - 1], errors - 1, beginpos + 1, LEN, data, outputs);
    }
    data.Prepare();
}

struct Maps {
    std::vector<MutableIncMap> beginmaps;
    std::vector<MutableIncMap> endmaps;

    Maps(int errors, const CRCOutputs& outputs) {
        uint64_t count1 = 0, count2 = 0;
        printf("Building maps for %i errors...", errors);
        beginmaps.resize(LEN);
        endmaps.resize(LEN);
        for (int i = 0; i < LEN; i++) {
            BuildBeginMap(errors, i, beginmaps[i], outputs);
            count1 += beginmaps[i].GetTotal();
            BuildEndMap(errors, i, endmaps[i], outputs);
            count2 += endmaps[i].GetTotal();
        }
        assert(count1  == count2);
        printf("done (%llu combinations)\n", (unsigned long long)count1);
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
        uint64_t iter = 0;
        for (int i = 1; i < LEN - 1; i++) {
            const IncMap& endmap = endmaps[i];
            for (int j = i + 1; j < LEN; j++) {
                const IncMap& beginmap = other.beginmaps[j];
                IncMap::iterator eit = endmap.begin();
                IncMap::iterator bit = beginmap.begin();
                while (eit.Valid() && bit.Valid()) {
                    if (((++iter) % 0xFFFFFFF) == 0) printf(".");
                    if (eit.GetKey() < bit.GetKey()) {
                        eit.Increment();
                        continue;
                    }
                    if (bit.GetKey() < eit.GetKey()) {
                        bit.Increment();
                        continue;
                    }
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

double Combination(int k, int n) {
    double num = 1.0;
    double den = 1.0;
    if (n - k < k) k = n - k;
    for (int i = 1; i <= k; i++) {
        num *= (n - i + 1);
        den *= i;
    }
    return num / den;
}

}

#define COMPUTEDISTANCE 6

int main(int argc, char** argv) {
    if (argc != 33) {
        fprintf(stderr, "Usage: %s v0 v1 v2 v3... v31\n", argv[0]);
        return(1);
    }
    uint32_t tbl[32];
    for (int i = 0; i < 32; i++) {
        unsigned long long r = strtoul(argv[i + 1], NULL, 0);
        if (r >> CHECKSUMBITS) {
            fprintf(stderr, "Error: table entry %i is outside of range\n", i);
            return(1);
        }
        tbl[i] = r;
    }
    std::vector<Maps> list;
    uint64_t output[COMPUTEDISTANCE + 1] = {0};
    bool computed[COMPUTEDISTANCE + 1] = {false};

    setbuf(stdout, NULL);
    for (int i = 1; i <= (COMPUTEDISTANCE+1)/2; i++) {
        list.push_back(Maps(i, tbl));
        if (!computed[i]) {
            uint64_t directfail = (unsigned long long)list[i - 1].Failures();
            printf("Undetected HD%i... %llu\n", i, (unsigned long long)directfail);
            output[i] = directfail;
            computed[i] = true;
        }
        for (int j = 1; j <= i; j++) {
            if (j + i <= COMPUTEDISTANCE && !computed[j + i]) {
                printf("Undetected HD%i...", i + j);
                uint64_t compoundfail = list[i - 1].FailuresCombined(list[j - 1]);
                printf(" %llu\n", (unsigned long long)compoundfail);
                output[j + i] = compoundfail;
                computed[i] = true;
            }
        }
    }
    for (int i = 1; i <= COMPUTEDISTANCE; i++) {
        if (i > 1) printf(",");
        if (output[i]) {
            printf("%.16g", output[i] / (Combination(i, LEN) * pow(MAXERR, i)) * 1073741824);
        } else {
            printf("0");
        }
    }
    printf("\n");
    return 0;
}
