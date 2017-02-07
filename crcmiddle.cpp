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

#define FIELD 9

#define CHECKSYMBOLS 6

#define CHECKSUMBITS (CHECKSYMBOLS * 5)

#define BASEBITS (5*(CHECKSYMBOLS-1))

#define REQUIRE_ZEROES 1
//#define REQUIRE_ZEROES 0

/* BCH codes over GF(2^5)
 */
uint64_t compute_bch(const uint8_t* data, int len, const uint64_t* tbl) {
    uint64_t l = 0;

    while (len > CHECKSYMBOLS) {
        uint8_t c = *(data++);
        uint8_t e = (l >> BASEBITS) ^ c;
        l = (l & ((((uint64_t)1) << BASEBITS) - 1)) << 5;
        for (int i = 0; i < 5; i++) {
            l ^= (~(((uint64_t)((e >> i) & 1)) - 1)) & tbl[i];
        }
        len--;
    }

    uint64_t f = 0;
    for (int i = 0; i < CHECKSYMBOLS; i++) {
        f <<= 5;
        f |= *(data++);
    }
    return l ^ f;
}

#define LEN 256
#define MAXERR 31

struct CRCOutputs {
    uint64_t val[LEN][MAXERR];

    CRCOutputs(const uint64_t *tbl, int len) {
        unsigned char data[LEN] = {0};
        uint64_t none = compute_bch(data, len, tbl);
        for (int pos = 0; pos < len; pos++) {
            for (int v = 0; v < MAXERR; v++) {
                data[pos] = v + 1;
                val[pos][v] = compute_bch(data, len, tbl) ^ none;
            }
            data[pos] = 0;
        }
    }
};

class IncMap {
    std::vector<uint8_t> data;
    uint64_t lastkey;

    IncMap(const IncMap& x) = delete;
    IncMap& operator=(const IncMap& x) = delete;

    void WriteNum(uint64_t num) {
//        printf("[W %lu]", (unsigned long)num);
        while (num >= 128) {
            data.push_back(0x80 | (num & 0x7F));
            num >>= 7;
        }
        data.push_back(num);
    }

    static uint64_t ReadNum(std::vector<uint8_t>::const_iterator& it) {
        uint64_t ret = 0;
        int shift = 0;
        do {
            uint64_t r = *(it++);
            ret |= ((r & 0x7F) << shift);
            if (r & 0x80) {
                shift += 7;
            } else {
                break;
            }
        } while (true);
//        printf("[R %lu]", (unsigned long)ret);
        return ret;
    }

    void AppendOrdered(uint64_t key, uint64_t value) {
        assert(key != lastkey);
        WriteNum(key - lastkey);
        WriteNum(value - 1);
        lastkey = key;
    }

public:

    IncMap() : data{0x00}, lastkey(0) {}
    IncMap(IncMap&& m) noexcept { data = std::move(m.data); lastkey = m.lastkey; }

    IncMap& operator=(IncMap&& m) noexcept { data = std::move(m.data); lastkey = m.lastkey; return *this; }

    class iterator {
        bool done;
        uint64_t key;
        uint64_t value;
        std::vector<uint8_t>::const_iterator inner;

        void Next() {
            uint64_t skip = IncMap::ReadNum(inner);
            if (skip == 0) {
                done = true;
                return;
            }
            key += skip;
            value = IncMap::ReadNum(inner) + 1;
        }

        iterator(std::vector<uint8_t>::const_iterator inner_) : done(false), key((uint64_t)(-1)), inner(inner_) {
            Next();
        }

    public:
        uint64_t GetKey() {
            return key;
        }
        uint64_t GetValue() {
            return value;
        }
        bool Valid() {
            return !done;
        }
        void Increment() {
            if (!done) {
                Next();
            }
        }

        friend class IncMap;
    };

    iterator begin() const { return iterator(data.begin()); }

    IncMap(const IncMap& old, const std::vector<uint64_t>& sortednew) : lastkey((uint64_t)(-1)) {
        iterator oldit = old.begin();
        std::vector<uint64_t>::const_iterator newit = sortednew.begin();
        uint64_t key = 0;
        data.reserve(old.data.size() + 3 * sortednew.size());
        while(true) {
            uint64_t value = 0;
            if (oldit.Valid() && oldit.GetKey() == key) {
                value += oldit.GetValue();
                oldit.Increment();
            }
            while (newit != sortednew.end() && *newit == key) {
                value++;
                newit++;
            }
            if (value > 0) {
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
/*        auto it = begin();
        while (it.Valid()) {
            printf("[MAP %lu:%lu] ", (unsigned long)it.GetKey(), (unsigned long)it.GetValue());
            it.Increment();
        }*/
    }

protected:
    size_t Memusage() const {
        return data.size();
    }
};

class MutableIncMap : public IncMap {
    std::vector<uint64_t> unsorted;
    uint64_t total;

    void Merge() {
        std::sort(unsorted.begin(), unsorted.end());
        *(static_cast<IncMap*>(this)) = IncMap(*this, unsorted);
        unsorted.clear();
    }

public:
    MutableIncMap() : total(0) {}

    void Increment(uint64_t key) {
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

uint64_t SimpleRecurse(uint64_t accum, int errors, int minpos, const CRCOutputs& outputs) {
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

void Recurse(uint64_t accum, int errors, int minpos, int endpos, MutableIncMap& data, const CRCOutputs& outputs) {
    if (errors == 0) {
//        printf("[ADD %lu]", (unsigned long)accum);
        data.Increment(accum);
        if ((data.GetTotal() & 0xFFFFFF) == 0) {
            fprintf(stderr, ".");
        }
        return;
    }
    for (int pos = minpos; pos < endpos - errors + 1; pos++) {
        for (int err = 1; err <= MAXERR; err++) {
            Recurse(accum ^ outputs.val[pos][err - 1], errors - 1, pos + 1, endpos, data, outputs);
        }
    }
}

/* Always sets positions beginpos and endpos */
void BuildMap(int errors, int beginpos, int endpos, MutableIncMap& data, const CRCOutputs& outputs, bool reduce) {
//    [5~printf("[begin %i in %i..%i] ", errors, beginpos, endpos);
    if (errors == 1) {
        if (beginpos == endpos) {
            for (int err = 1; err <= (reduce ? 1 : MAXERR); err++) {
//                printf("[ADD %lu]", (unsigned long)outputs.val[endpos][err - 1]);
                data.Increment(outputs.val[endpos][err - 1]);
            }
        }
    } else {
        if (beginpos != endpos) {
            for (int err = 1; err <= MAXERR; err++) {
                for (int err2 = 1; err2 <= (reduce ? 1 : MAXERR); err2++) {
                    Recurse(outputs.val[beginpos][err - 1] ^ outputs.val[endpos][err2 - 1], errors - 2, beginpos + 1, endpos, data, outputs);
                }
            }
        }
    }
//    printf("[end %i in %i..%i] ", errors, beginpos, endpos);
    data.Prepare();
}

struct BeginMaps {
    std::vector<MutableIncMap> beginmaps;

    BeginMaps(int errors, int len, const CRCOutputs& outputs) {
        uint64_t count1 = 0;
#ifndef REQUIRE_ZEROES
        fprintf(stderr, "Building begin maps for %i errors...", errors);
#endif
        beginmaps.resize(len);
        for (int i = 0; i < len; i++) {
            BuildMap(errors, i, len - 1, beginmaps[i], outputs, true);
            count1 += beginmaps[i].GetTotal();
        }
#ifndef REQUIRE_ZEROES
        fprintf(stderr, "done (%llu combinations)\n", (unsigned long long)count1);
#endif
    }

    uint64_t Failures(int len) {
        uint64_t ret1 = 0;
        IncMap::iterator it1 = beginmaps[0].begin();
        if (it1.Valid() && it1.GetKey()==0) ret1 += it1.GetValue();
        return ret1 * MAXERR;
    }
};

struct EndMaps {
    int errors;
    int curlen;
    const CRCOutputs& outputs;
    std::vector<MutableIncMap> endmaps;

    EndMaps(int errors_, const CRCOutputs& outputs_) : errors(errors_), outputs(outputs_) {}

    void Extend(int len) {
        uint64_t count2 = 0;
#ifndef REQUIRE_ZEROES
        fprintf(stderr, "Building end maps for %i errors...", errors);
#endif
        endmaps.resize(len);
        while (curlen < len) {
            BuildMap(errors, 0, curlen, endmaps[curlen], outputs, false);
            count2 += endmaps[curlen].GetTotal();
            curlen++;
        }
#ifndef REQUIRE_ZEROES
        fprintf(stderr, "done (%llu combinations)\n", (unsigned long long)(count2));
#endif
    }

    uint64_t FailuresCombined(const BeginMaps& other, int len) {
        uint64_t ret = 0;
        uint64_t iter = 0;
        for (int i = 0; i < len - 1; i++) {
            const IncMap& endmap = endmaps[i];
            for (int j = i + 1; j < len; j++) {
                const IncMap& beginmap = other.beginmaps[j];
                IncMap::iterator eit = endmap.begin();
                IncMap::iterator bit = beginmap.begin();
                while (eit.Valid() && bit.Valid()) {
                    if (((++iter) & 0xFFFFFFF) == 0) fprintf(stderr, ".");
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
        return ret * MAXERR;
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

#define COMPUTEDISTANCE 5


int analyse(uint64_t code, int codelen, int maxtestlen, uint64_t num) {
    uint64_t tbl[5] = {code,0,0,0,0};
    for (int i = 1; i < 5; i++) {
        for (int j = 0; j < CHECKSYMBOLS; j++) {
            unsigned int prev = (tbl[i - 1] >> (5 * j)) & 0x1f;
            tbl[i] |= ((uint64_t)(((prev & 0xf) << 1) ^ (FIELD * (prev >> 4)))) << (5 * j);
        }
    }
    std::vector<EndMaps> endlist;
    CRCOutputs outputs(tbl, codelen);
    setbuf(stdout, NULL);
    for (int i = 1; i <= (COMPUTEDISTANCE)/2; i++) {
        endlist.push_back(EndMaps(i, outputs));
    }
    for (int testlen = 1; testlen <= maxtestlen; testlen++) {
        for (int i = 1; i <= (COMPUTEDISTANCE)/2; i++) {
            endlist[i - 1].Extend(testlen);
        }
        uint64_t output[COMPUTEDISTANCE + 1] = {0};
        bool computed[COMPUTEDISTANCE + 1] = {false};
        std::vector<BeginMaps> beginlist;
        for (int i = 1; i <= (COMPUTEDISTANCE+1)/2; i++) {
            beginlist.push_back(BeginMaps(i, testlen, outputs));
            if (!computed[i]) {
                uint64_t directfail = (unsigned long long)beginlist[i - 1].Failures(testlen);
#if REQUIRE_ZEROES
                if (directfail) return 0;
#else
                fprintf(stderr, "Undetected HD%i... %llu\n", i, (unsigned long long)directfail);
#endif
                output[i] = directfail;
                computed[i] = true;
            }
            for (int j = 1; j <= i; j++) {
                if (j + i <= COMPUTEDISTANCE && !computed[j + i]) {
                    uint64_t compoundfail = endlist[j - 1].FailuresCombined(beginlist[i - 1], testlen);
#if REQUIRE_ZEROES
                    if (compoundfail) return 0;
#else
                    fprintf(stderr, "Undetected HD%i...", i + j);
                    fprintf(stderr, " %llu\n", (unsigned long long)compoundfail);
#endif
                    output[j + i] = compoundfail;
                    computed[i] = true;
                }
            }
        }
#if REQUIRE_ZEROES
        if (testlen < maxtestlen) continue;
        printf("%llu ", (unsigned long long)num);
#endif
        printf("\"0x%lx 0x%lx 0x%lx 0x%lx 0x%lx\",%i,%i", (unsigned long)tbl[0], (unsigned long)tbl[1], (unsigned long)tbl[2], (unsigned long)tbl[3], (unsigned long)tbl[4], codelen, testlen);
        for (int i = 1; i <= COMPUTEDISTANCE; i++) {
            if (output[i]) {
                printf(",%.16g", output[i] / (Combination(i - 2, testlen - 2) * pow(MAXERR, i)) * 1073741824);
            } else {
                printf(",0");
            }
        }
        printf("\n");
    }
    return 1;
}

/*
int main(int argc, char** argv) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s codelen testlen gen\n", argv[0]);
        return(1);
    }
    int codelen = strtoul(argv[1], NULL, 0);
    int maxtestlen = strtoul(argv[2], NULL, 0);
    unsigned long long r = strtoul(argv[3], NULL, 0);
    if (r >> CHECKSUMBITS) {
        fprintf(stderr, "Error: generator %llx is outside of range\n", r);
        return(1);
    }
    analyse(r, codelen, maxtestlen, 0);
    return 0;
}
*/

int main(int argc, char** argv) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s codelen testlen\n", argv[0]);
        return(1);
    }
    int codelen = strtoul(argv[1], NULL, 0);
    int maxtestlen = strtoul(argv[2], NULL, 0);
    uint64_t num = 0;
    while (1) {
        char c[1024];
        if (!fgets(c, sizeof(c), stdin)) {
            break;
        }
        uint64_t r = strtoul(c, NULL, 0);
        assert((r >> CHECKSUMBITS) == 0);
        analyse(r, codelen, maxtestlen, num++);
    }
    return 0;
}
