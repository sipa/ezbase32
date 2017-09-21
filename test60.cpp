#include <stdint.h>
#include <stdio.h>
#include <vector>
#include <string.h>
#include <algorithm>
#include <tuple>
#include <set>
#include <map>
#include <thread>

#include "tinyformat.h"

#define DEGREE 12
#define LENGTH 64
#define ERR 6
#define TABLE_ERR 3
#define THREADS 48

static inline uint32_t rdrand() {
    uint32_t ret;
    __asm__ volatile(".byte 0x0f, 0xc7, 0xf0" : "=a"(ret) :: "cc");
    return ret;
}

template<class ForwardIt, class T, class Compare=std::less<>>
ForwardIt binary_find(ForwardIt first, ForwardIt last, const T& value, Compare comp={})
{
    first = std::lower_bound(first, last, value, comp);
    return first != last && !comp(value, *first) ? first : last;
}

constexpr uint64_t mul2(uint64_t x) {
    uint64_t high = 0x842108421084210 & x;
    uint64_t low = 0x7BDEF7BDEF7BDEF & x;
    return (low << 1) ^ (high >> 4) ^ (high >> 1);
}

uint64_t scramble(uint64_t x) {
    x = ((x << 5) | (x >> 55)) & 0xFFFFFFFFFFFFFFFULL;
    x = (x * 0xbff548e311f5fdbULL) & 0xFFFFFFFFFFFFFFFULL;
    return x;
}

constexpr uint64_t comb(uint64_t n, uint64_t k) {
    uint64_t ret = 1;
    for (uint64_t i = 1; i <= k; ++i) {
        ret *= n;
        ret *= 31;
        ret /= i;
        --n;
    }
    return ret;
}

static const char *charset = "0123456789ABCDEFGHIJKLMNOPQRSTUV";

uint64_t checksums[LENGTH][32];
uint64_t GEN1, GEN2, GEN4, GEN8, GEN16;

static uint64_t step(uint64_t pre) {
    uint8_t b = pre >> 55;
    return ((pre & 0x7FFFFFFFFFFFFF) << 5) ^
        (-(uint64_t)((b >> 0) & 1) & GEN1) ^
        (-(uint64_t)((b >> 1) & 1) & GEN2) ^
        (-(uint64_t)((b >> 2) & 1) & GEN4) ^
        (-(uint64_t)((b >> 3) & 1) & GEN8) ^
        (-(uint64_t)((b >> 4) & 1) & GEN16);
}

struct entry {
    uint64_t syndrome;
    signed char pos[3], err[3];
    signed char minpos, maxpos;

    entry() : syndrome(0), pos { -1, -1, -1 }, err { 0, 0, 0 }, minpos(-1), maxpos(-1) {}
};

bool operator<(const entry& a, const entry& b) { return a.syndrome < b.syndrome; }

static_assert(sizeof(entry) <= 16, "OO");

std::vector<double> table_mult;
std::vector<size_t> table_below;
std::vector<size_t> table_above;
std::vector<entry> table[TABLE_ERR + 1];

void build_table(uint64_t syndrome, int done, entry& ent, int minpos, int maxpos, int maxerr) {
    if (done == maxerr) {
        ent.syndrome = scramble(syndrome);
        table[maxerr].emplace_back(ent);
        return;
    }
    for (int p = minpos; p < maxpos; ++p) {
        for (int err = 1; err < 32; ++err) {
            if (done == 0) ent.minpos = p;
            if (done == maxpos - 1) ent.maxpos = p;
            ent.pos[done] = p;
            ent.err[done] = err;
            build_table(syndrome ^ checksums[p][err], done + 1, ent, p + 1, maxpos, maxerr);
        }
    }
}

uint64_t count_fixes(uint64_t syndrome, int err, int min_pos) {
    if (err <= TABLE_ERR) {
        entry ent;
        ent.syndrome = scramble(syndrome);
/*        if (ent.syndrome < table[err].front().syndrome) return 0;
        if (ent.syndrome > table[err].back().syndrome) return 0;
        size_t exp = (ent.syndrome - table[err].front().syndrome) * table_mult[err] + 0.5;
        size_t low = exp < table_below[err] ? 0 : exp - table_below[err];
        size_t high = exp + table_above[err] + 1 >= table[err].size() ? table[err].size() : exp + table_above[err] + 1;
        auto it = binary_find(table[err].cbegin() + low, table[err].cbegin() + high, ent);
        if (it != table[err].cbegin() + high) {
            return (it->minpos >= min_pos);
        }*/
        auto it = binary_find(table[err].cbegin(), table[err].cend(), ent);
        if (it != table[err].cend()) {
            return (it->minpos >= min_pos);
        }
        return 0;
    }
    uint64_t count = 0;
    for (int pos = min_pos; pos < LENGTH; ++pos) {
        for (int e = 1; e < 32; ++e) {
            count += count_fixes(syndrome ^ checksums[pos][e], err - 1, pos + 1);
        }
    }
    return count;
}

void thread_runner() {
    while(true) {
        uint64_t sum = 0;
        for (int i = 0; i < ERR; ++i) {
            sum ^= checksums[rdrand() % LENGTH][1 + (rdrand() % 31)];
        }
        std::string report = strprintf("[%x]", sum);
        for (int i = 0; i <= ERR; ++i) {
            report += strprintf(" %i:%lu", i, count_fixes(sum, i, 0));
        }
        report += "\n";
        printf("%s", report.c_str());
    }
}

int main(int argc, char** argv) {
    if (argc > 1) {
        if (strlen(argv[1]) == 12) {
            uint64_t gen = 0;
            for (int i = 0; i < 12; ++i) {
                const char *ptr = strchr(charset, toupper(argv[1][DEGREE - 1 - i]));
                if (ptr == nullptr) {
                    fprintf(stderr, "Unknown character '%c'\n", argv[1][DEGREE - 1 - i]);
                    return 1;
                }
                gen |= ((uint64_t)(ptr - charset)) << (i * 5);
            }
            GEN1 = gen;
            GEN2 = mul2(GEN1);
            GEN4 = mul2(GEN2);
            GEN8 = mul2(GEN4);
            GEN16 = mul2(GEN8);
        } else {
            fprintf(stderr, "Wrong length\n");
            return 1;
        }
    } else {
        fprintf(stderr, "Use: ./%s GEN%i\n", argv[0], DEGREE);
        return 1;
    }
    for (uint64_t err = 0; err <= 32; ++err) {
        uint64_t syndrome = err;
        for (size_t pos = 0; pos < LENGTH; ++pos) {
            checksums[pos][err] = syndrome;
            syndrome = step(syndrome);
        }
    }
    entry ent;
    for (int size = 0; size <= TABLE_ERR; ++size) {
        table[size].reserve(comb(LENGTH, size));
        build_table(0, 0, ent, 0, LENGTH, size);
        std::sort(table[size].begin(), table[size].end());

        double mult = size > 0 ? ((double)(table[size].size() - 1)) / (table[size].back().syndrome - table[size].front().syndrome) : 1;

        size_t below = 0, sbelow = 0;
        size_t above = 0, sabove = 0;
        for (size_t pos = 0; pos < table[size].size(); ++pos) {
            size_t exppos = ((table[size][pos].syndrome - table[size].front().syndrome) * mult + 0.5);
            if (pos > exppos) {
                sabove += pos - exppos;
                if (pos - exppos > above) above = pos - exppos;
            } else if (exppos > pos) {
                sbelow += exppos - pos;
                if (exppos - pos > below) below = exppos - pos;
            }
        }
        printf("%i errors\n", size);
        printf("* below=%lu above=%lu\n", (unsigned long)below, (unsigned long)above);
        printf("* abelow=%g aabove=%g\n", (double)sbelow / table[size].size(), (double)sabove / table[size].size());
        table_mult.push_back(mult);
        table_below.push_back(below);
        table_above.push_back(above);
    }

    std::vector<std::thread> threads;
    for (int i = 0; i < THREADS - 1; ++i) {
         threads.emplace_back(thread_runner);
    }
    thread_runner();
    return 0;
}

/*
static const char* charset = "qpzry9x8gf2tvdw0s3jn54khce6mua7l";

static const int8_t charset_rev[128] = {
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    15, -1, 10, 17, 21, 20, 26, 30,  7,  5, -1, -1, -1, -1, -1, -1,
    -1, 29, -1, 24, 13, 25,  9,  8, 23, -1, 18, 22, 31, 27, 19, -1,
     1,  0,  3, 16, 11, 28, 12, 14,  6,  4,  2, -1, -1, -1, -1, -1,
    -1, 29, -1, 24, 13, 25,  9,  8, 23, -1, 18, 22, 31, 27, 19, -1,
     1,  0,  3, 16, 11, 28, 12, 14,  6,  4,  2, -1, -1, -1, -1, -1
};


struct dpos {
    signed char pos1, pos2, pos3, pos4;
};

extern "C" {
struct table_entry {
    uint64_t syndrome;
    signed char pos1, pos2;
};
extern struct table_entry table[116059];
}

static bool compare_pair(const table_entry& p1, const table_entry& p2) {
    return p1.syndrome < p2.syndrome;
}

void buildtable() {
    std::set<uint64_t> entries;
    std::vector<table_entry> table_build;
    for (int e1 = 0; e1 < 2; ++e1) {
        uint64_t x1 = step(e1);
        for (signed char p1 = 1; e1 ? (p1 < 87) : (p1 == 1); ++p1) {
            for (int e2 = 0; e2 < 32; ++e2) {
                uint64_t x2 = e2;
                for (signed char p2 = 0; e2 ? (p2 < p1) : (p2 == 0); ++p2) {
                    if (entries.count(x1 ^ x2) == 0) {
                        table_build.emplace_back(table_entry{((x1 ^ x2) * 0x23cd3274b51129bULL) & 0xFFFFFFFFFFFFFFFULL, e1 ? p1 : (signed char)-1, e2 ? p2 : (signed char)-1});
                        entries.insert(x1 ^ x2);
                    }
                    x2 = step(x2);
                }
            }
            x1 = step(x1);
        }
    }
    std::sort(table_build.begin(), table_build.end(), compare_pair);
    uint64_t low = table_build.front().syndrome;
    uint64_t high = table_build.back().syndrome;
    std::map<size_t, size_t> off;
    size_t sum_off = 0;
    for (size_t pos = 0; pos < table_build.size(); ++pos) {
        size_t inter = (size_t)(((table_build[pos].syndrome - low) / (double)(high - low)) * (table_build.size() - 1) + 0.5);
        off[std::abs((long)(pos - inter))]++;
        sum_off += std::abs((long)(pos - inter));
    }
    printf("OFF: %g\n", ((double)sum_off) / table_build.size());
    for (const auto& pair : off) {
        printf("OFF %lu: %lu\n", (unsigned long)pair.first, (unsigned long)pair.second);
    }
    printf("#include <stdint.h>\n");
    printf("struct table_entry {\n");
    printf("    uint64_t syndrome;\n");
    printf("    signed char pos1, pos2;\n");
    printf("};\n");
    printf("const struct table_entry table[] = {\n");
    for (size_t i = 0; i < table_build.size(); ++i) {
        printf("    {0x% 14llxULL, %i, %i},\n", (unsigned long long)table_build[i].syndrome, table_build[i].pos1, table_build[i].pos2);
    }
    printf("};\n");
}

dpos locate(uint64_t x) {
    table_entry query;
    for (size_t loc = 0; loc < sizeof(table) / sizeof(table[0]); ++loc) {
        uint64_t entry = table[loc].syndrome;
        for (int i1 = 0; i1 < 31; ++i1) {
            query.syndrome = x ^ entry;
            for (int i2 = 0; i2 < 31; ++i2) {
                auto it = std::lower_bound(std::begin(table), std::end(table), query, compare_pair);
                if (it != std::end(table) && it->syndrome == query.syndrome) {
                    return dpos{table[loc].pos1, table[loc].pos2, it->pos1, it->pos2};
                }
                query.syndrome = gf_double(query.syndrome);
            }
            entry = gf_double(entry);
        }
    }

    return dpos{-2,-2,-2,-2};
}

int main(int argc, char** argv) {
    if (argc == 1) {
        uint64_t x = 1;
        for (int i = 0; i < 55; ++i) {
            uint8_t v = rdrand() & 0x1f;
            x = step(x) ^ v;
            printf("%c", charset[v]);
        }
        for (int i = 0; i < 12; ++i) {
            x = step(x);
        }
        x ^= 1;
        for (int i = 0; i < 12; ++i) {
            printf("%c", charset[(x >> ((11 - i) * 5)) & 0x1f]);
        }
        printf("\n");
        return 0;
    }

    buildtable();
    uint64_t x = 1;
    const char *s = argv[1];
    int len = strlen(s);
    printf("%s\n", s);
    while (*s) {
        x = step(x) ^ charset_rev[(int)*s];
        ++s;
    }
    x ^= 1;
    dpos pos = locate(x);
    for (int i = 0; i < len; ++i) {
        if (len - i - 1 == pos.pos1 || len - i - 1 == pos.pos2 || len - i - 1 == pos.pos3 || len - i - 1 == pos.pos4) {
            printf("^");
        } else {
            printf(" ");
        }
    }
    printf("\n");
    printf("Locs: %i, %i, %i, %i\n", (int)pos.pos1, (int)pos.pos2, (int)pos.pos3, (int)pos.pos4);
    return 0;
}
*/