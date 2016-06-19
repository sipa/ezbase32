#include <stdio.h>
#include <stdint.h>
#include <vector>
#include <set>
#include <string.h>
#include <assert.h>
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
    const std::string desc;

public:
    BCHCode(const std::string& desc_, uint32_t x0, uint32_t x1, uint32_t x2, uint32_t x3, uint32_t x4) : desc(desc_) {
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

    BCHCode(const BCHCode& code) : desc(code.desc) {
        memcpy(val, code.val, sizeof(val));
    }

    inline uint32_t syndrome(size_t pos, int err) {
        return val[pos][err - 1];
    }

    std::string GetDesc() const {
        return desc;
    }
};

struct BCHState {
    size_t codepos;
    uint64_t attempts;
    size_t len;
};

struct StateInfo {
    size_t len[4];
    size_t count[4];
    size_t firstiter;
    double rate;
    std::vector<std::string> names;
    std::vector<uint64_t> iterations;
};

struct BCHStatePtrComparator {
    bool operator()(const BCHState* a, const BCHState* b) {
        // Speed up the common case
        if (a == b) return false;

        // Codes that are candidates for higher lengths first
        if (a->len > b->len) return true;
        if (a->len < b->len) return false;

        // Code with fewer attempts at the same length first
        if (a->attempts < b->attempts) return true;
        if (a->attempts > b->attempts) return false;

        // Disambiguation
        if (a < b) return true;
        return false;
    }
};

class State {
    std::mutex mutex;
    const std::vector<BCHCode>* codes;
    std::vector<BCHState> states;
    std::set<BCHState*, BCHStatePtrComparator> stateptrs;

    double rate, weight, tim;

public:
    State(const std::vector<BCHCode>* codes_, size_t len_) : codes(codes_), rate(0), weight(0), tim(0) {
        states.resize(codes->size());
        for (size_t i = 0; i < states.size(); i++) {
            states[i].codepos = i;
            states[i].attempts = 0;
            states[i].len = len_;
            stateptrs.insert(&states[i]);
        }
    }

    const BCHCode& code(size_t pos) { return (*codes)[pos]; }

    std::pair<size_t, std::vector<size_t>> GetCodes(size_t num, size_t iterations) {
        std::vector<size_t> ret;
        ret.reserve(num);

        std::unique_lock<std::mutex> lock(mutex);

        auto it = stateptrs.begin();
        assert(it != stateptrs.end());
        size_t len = (*it)->len;
        ret.push_back((*it)->codepos);
        {
            auto oldit = it++;
            auto oldptr = *oldit;
            stateptrs.erase(oldit);
            oldptr->attempts += iterations;
            stateptrs.insert(oldptr);
        }

        while (it != stateptrs.end() && ret.size() < num) {
            if ((*it)->len < len) {
                break;
            }
            ret.push_back((*it)->codepos);
            {
                auto oldit = it++;
                auto oldptr = *oldit;
                stateptrs.erase(oldit);
                oldptr->attempts += iterations;
                stateptrs.insert(oldptr);
            }
        }


        return std::make_pair(len, std::move(ret));
    }

    void Failed(size_t codepos, size_t faillen) {
        std::unique_lock<std::mutex> lock(mutex);
        if (faillen <= states[codepos].len) {
            stateptrs.erase(&states[codepos]);
            states[codepos].attempts = 0;
            states[codepos].len = faillen - 1;
            stateptrs.insert(&states[codepos]);
        }
    }

    void Update(uint64_t count, double starttime) {
        std::unique_lock<std::mutex> lock(mutex);

        double t = timer();
        double duration = t - starttime;
        double coef = pow(0.99, t - tim);
        rate = coef * rate + count;
        weight = coef * weight + duration;
        tim = t;
    }

    StateInfo GetInfo() {
        std::unique_lock<std::mutex> lock(mutex);

        StateInfo ret = {};
        auto it = stateptrs.begin();
        int got = -1;
        size_t oldlen = 0;
        size_t count = 0;

        while (it != stateptrs.end() && got < 4) {
            if ((*it)->len != oldlen) {
                if (got >= 0) {
                    ret.len[got] = oldlen;
                    ret.count[got] = count;
                } else {
                    ret.firstiter = (*it)->attempts;
                }
                got++;
                count = 1;
                oldlen = (*it)->len;
            } else {
                count++;
            }
            if (got == 0) {
                ret.names.push_back(code((*it)->codepos).GetDesc());
                ret.iterations.push_back((*it)->attempts);
            }
            it++;
        }

        ret.rate = rate / weight;
        return ret;
    }
};

void ThreadCheck(State* state, size_t num, size_t errs, size_t iterations) {
    Rander rander;
    do {
        double start = timer();
        auto x = state->GetCodes(num, iterations);
        int bits = 8 * sizeof(unsigned int) - __builtin_clz((unsigned int)x.first);
        std::vector<BCHCode> codes;
        codes.reserve(x.second.size());
        for (auto p : x.second) {
            codes.emplace_back(state->code(p));
        }

        std::vector<uint32_t> crc;
        std::vector<size_t> errpos;
        errpos.resize(errs);
        for (size_t i = 0; i < iterations; i++) {
            size_t firstpos = 0;
            size_t lastpos = 0;
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
                if (e == 0) {
                    firstpos = errpos[e];
                    lastpos = errpos[e];
                } else {
                    firstpos = std::min(firstpos, errpos[e]);
                    lastpos = std::max(lastpos, errpos[e]);
                }
            }
            for (size_t c = 0; c < codes.size(); c++) {
                if (crc[c] == 0) {
                    state->Failed(x.second[c], lastpos - firstpos + 1);
                }
            }
        }
        state->Update(iterations * codes.size(), start);
    } while(true);
}

void ThreadDump(State* state) {
    do {
        sleep(10);
        auto x = state->GetInfo();
        printf("[%g Gi, %g Gi/s] {%i: %llu} {%i: %llu} {%i: %llu} {%i: %llu}\n", x.firstiter / 1073741824.0, x.rate / 1073741824.0, (int)x.len[0], (unsigned long long)x.count[0], (int)x.len[1], (unsigned long long)x.count[1], (int)x.len[2], (unsigned long long)x.count[2], (int)x.len[3], (unsigned long long)x.count[3]);
        FILE* file = fopen("crcelim.dump.tmp", "w");
        assert(file != NULL);
        for (size_t i = 0; i < x.names.size(); i++) {
            fprintf(file, "%s # %llu iterations\n", x.names[i].c_str(), (unsigned long long)x.iterations[i]);
        }
        fclose(file);
        rename("crcelim.dump.tmp", "crcelim.dump");
    } while(true);
}

int main(void) {
    setbuf(stdout, NULL);
    std::vector<BCHCode> codes;
    unsigned long l[5];
    char c[256];
    while (fgets(c, sizeof(c), stdin)) {
        if (sscanf(c, "0x%lx 0x%lx 0x%lx 0x%lx 0x%lx\n", &l[0], &l[1], &l[2], &l[3], &l[4]) == 5) {
            while (isspace(c[strlen(c) - 1])) {
                c[strlen(c) - 1] = 0;
            }
            codes.emplace_back(std::string(c), l[0], l[1], l[2], l[3], l[4]);
        }
    }
    printf("Got %lu codes\n", (unsigned long)codes.size());
    State state(&codes, MAXLEN);
    for (int i = 0; i < 8; i++) {
        std::thread th(ThreadCheck, &state, 2 << 5, 4, 2 << 16);
        th.detach();
    }
    ThreadDump(&state);
    return 0;
}
