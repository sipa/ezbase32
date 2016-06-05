#include <stdio.h>
#include <stdint.h>
#include <vector>
#include <set>
#include <string.h>
#include <immintrin.h>
#include <unistd.h>
#include <mutex>
#include <time.h>
#include <thread>
#include <math.h>
#include <tuple>

#define MAXLEN 80

static uint32_t bch(const uint8_t* data, size_t len, uint32_t x0, uint32_t x1, uint32_t x2, uint32_t x3, uint32_t x4) {
    uint32_t checksum = 0;
    while (len) {
        uint32_t b = (checksum >> 25) ^ *(data++);
        len--;
        checksum = (checksum & 0x1FFFFFF) << 5;
        checksum ^= -((b & 1) != 0) & x0;
        checksum ^= -((b & 2) != 0) & x1;
        checksum ^= -((b & 4) != 0) & x2;
        checksum ^= -((b & 8) != 0) & x3;
        checksum ^= -((b & 16) != 0) & x4;
    }
    return checksum;
}

static double timer(void) {
    struct timespec ret;
    clock_gettime(CLOCK_THREAD_CPUTIME_ID, &ret);
    return ret.tv_nsec * 0.000000001 + ret.tv_sec;
}

class Rander {
    unsigned char data[4096];
    int pos;

    void Step() {
        if ((pos & 0xFFF) == 0) {
            for (int i = 0; i < 4096; i += sizeof(unsigned long long)) {
                _rdrand64_step((unsigned long long*)(data + i));
            }
            pos = 0;
        }
    }

public:
    Rander() {
        memset(data, 0, sizeof(data));
        pos = 0;
    }

    uint8_t GetByte() {
        Step();
        uint8_t ret = data[pos];
        pos++;
        return ret;
    }

    uint8_t GetInt(uint8_t max, int bits) {
        uint8_t r;
        do {
            r = GetByte() & ((1 << bits) - 1);
        } while (r >= max);
        return r;
    }

    uint32_t GetBigInt(uint32_t max) {
        int bits = 8 * sizeof(unsigned long) - __builtin_clzl((unsigned long)max);
        uint32_t r;
        do {
            r = (((uint32_t)GetByte()) << 24 | ((uint32_t)GetByte()) << 16 | ((uint32_t)GetByte()) << 8 | ((uint32_t)GetByte())) & ((((uint32_t)1) << bits) - 1);
        } while (r >= max);
        return r;
    }
};

class BCHCode {
    uint32_t val[MAXLEN][31];

public:
    BCHCode(uint32_t x0, uint32_t x1, uint32_t x2, uint32_t x3, uint32_t x4) {
        uint8_t data[MAXLEN] = {0};
        uint32_t base = bch(data, MAXLEN, x0, x1, x2, x3, x4);
        for (int i = 0; i < MAXLEN; i++) {
            for (int j = 1; j <= 31; j++) {
                data[i] = j;
                val[i][j - 1] = bch(data, MAXLEN, x0, x1, x2, x3, x4) ^ base;
            }
            data[i] = 0;
        }
    }

    BCHCode(const BCHCode& code) {
        memcpy(val, code.val, sizeof(val));
    }

    inline uint32_t syndrome(size_t pos, int err) {
        return val[pos][err - 1];
    }
};

struct StateInfo {
    size_t len;
    size_t count;
    double rate;
    double progress;
};

class State {
    std::mutex mutex;
    const std::vector<BCHCode>* codes;

    size_t len;
    std::set<uint16_t> candidates;
    std::vector<uint16_t> vcandidates;
    double progress;

    double tim;
    double rate;
    double weight;

public:
    State(const std::vector<BCHCode>* codes_, size_t len_) : codes(codes_), len(len_ + 1), progress(0), tim(0), rate(0), weight(0) {}

    std::pair<size_t, std::vector<const BCHCode*>> GetCodes(size_t num, Rander& rander) {
        std::vector<const BCHCode*> ret;
        ret.reserve(num);

        std::unique_lock<std::mutex> lock(mutex);
        if (candidates.empty()) {
            for (size_t i = 0; i < codes->size(); i++) {
                candidates.insert(i);
            }
            len--;
            progress = 0;
        }

        if (vcandidates.size() != candidates.size()) {
            vcandidates.assign(candidates.begin(), candidates.end());
        }

        for (size_t i = 0; i < num; i++) {
            std::swap(vcandidates[i], vcandidates[i + rander.GetBigInt(codes->size() - i)]);
            ret.push_back(&((*codes)[vcandidates[i]]));
        }

        return std::make_pair(len, std::move(ret));
    }

    void Failed(const BCHCode* code, size_t lendone) {
        std::unique_lock<std::mutex> lock(mutex);
        if (lendone == len) {
            candidates.erase(code - &((*codes)[0]));
        }
    }

    StateInfo GetInfo() {
        std::unique_lock<std::mutex> lock(mutex);
        return StateInfo{len, candidates.size(), rate / weight, progress};
    }

    void Update(double timer, double count) {
        std::unique_lock<std::mutex> lock(mutex);

        double coef = pow(0.99, timer - tim);
        rate = coef * rate + (1.0 - coef) * count / (timer - tim);
        weight = weight * coef + (1.0 - coef);
        tim = timer;
        progress += count / candidates.size();
    }
};

void ThreadCheck(State* state, size_t num, size_t errs, size_t iterations) {
    Rander rander;
    do {
        auto x = state->GetCodes(num, rander);
        int bits = 8 * sizeof(unsigned int) - __builtin_clz((unsigned int)x.first);
        std::vector<BCHCode> codes;
        codes.reserve(x.second.size());
        for (auto p : x.second) {
            codes.emplace_back(*p);
        }

        std::vector<uint32_t> crc;
        std::vector<size_t> errpos;
        errpos.resize(errs);
        for (size_t i = 0; i < iterations; i++) {
            crc.assign(codes.size(), 0);
            for (size_t e = 0; e < errs; e++) {
                int ok;
                do {
                    errpos[e] = rander.GetInt(x.first, bits);
                    ok = 1;
                    for (size_t f = 0; f < e; f++) {
                        if (errpos[e] == errpos[f]) {
                            ok = 0;
                            break;
                        }
                    }
                } while (!ok);
                int mis = 1 + rander.GetInt(31, 5);
                for (size_t c = 0; c < codes.size(); c++) {
                    crc[c] ^= codes[c].syndrome(errpos[e], mis);
                }
            }
            for (size_t c = 0; c < codes.size(); c++) {
                if (crc[c] == 0) {
                    state->Failed(x.second[c], x.first);
                }
            }
        }
        state->Update(timer(), (double)iterations * (double)codes.size());
    } while(true);
}

void ThreadDump(State* state) {
    do {
        sleep(2);
        auto x = state->GetInfo();
        printf("%lu length-%i codes left (progress %g); %g TPF\n", (unsigned long)x.count, (int)x.len, x.progress / 1073741824.0, 1073741824.0 / x.rate);
    } while(true);
}

int main(void) {
    std::vector<BCHCode> codes;
    unsigned long l[5];
    while (scanf("0x%lx 0x%lx 0x%lx 0x%lx 0x%lx\n", &l[0], &l[1], &l[2], &l[3], &l[4]) == 5) {
        codes.emplace_back(l[0], l[1], l[2], l[3], l[4]);
    }
    printf("Got %lu codes\n", (unsigned long)codes.size());
    State state(&codes, 80);
    for (int i = 0; i < 4; i++) {
        std::thread th(ThreadCheck, &state, 2 << 5, 4, 2 << 16);
        th.detach();
    }
    ThreadDump(&state);
    return 0;
}
