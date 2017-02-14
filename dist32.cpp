#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <immintrin.h>

#include <string>
#include <vector>
#include <set>

static const char badpairs[][3] = {
    {'0', 'o', 9}, {'1', '7', 5}, {'1', 'i', 9}, {'1', 'l', 10},
    {'2', 'z', 2}, {'3', 'e', 1}, {'4', 'a', 1}, {'4', 'h', 1},
    {'4', '9', 2}, {'5', 's', 2}, {'6', 'b', 3}, {'8', 'b', 3},
    {'9', 'p', 1}, {'a', 'c', 2}, {'a', 'd', 2}, {'a', 'e', 2},
    {'a', 'o', 3}, {'a', 'h', 1}, {'b', 'd', 2}, {'b', 'h', 2},
    {'b', 'e', 1}, {'b', 'p', 2}, {'b', 'r', 2}, {'c', 'e', 2},
    {'c', 'o', 3}, {'c', 'g', 2}, {'d', 'o', 2}, {'e', 'o', 2},
    {'e', 'f', 2}, {'f', 'p', 1}, {'f', 't', 3}, {'g', 'q', 2},
    {'g', 'o', 1}, {'h', 'k', 1}, {'h', 'n', 4}, {'i', 'j', 5},
    {'i', 'l', 10}, {'k', 'x', 1}, {'m', 'n', 1}, {'n', 'r', 1},
    {'o', 'p', 1}, {'o', 'q', 4}, {'p', 'q', 1}, {'p', 'r', 2},
    {'u', 'v', 1}, {'v', 'w', 1}, {'v', 'x', 1}, {'v', 'y', 2},
    {'x', 'y', 1}
};

static inline uint64_t rotl(const uint64_t x, int k) {
    return (x << k) | (x >> (64 - k));
}

class Rander {
    uint64_t s[2];
    uint32_t count;

    uint64_t next;
    int bits;

    void Produce(void) {
        const uint64_t s0 = s[0];
        uint64_t s1 = s[1];
        next = s0 + s1;
        bits = 64;

        s1 ^= s0;
        s[0] = rotl(s0, 55) ^ s1 ^ (s1 << 14); // a, b
        s[1] = rotl(s1, 36); // c
    }

    void Step() {
        static_assert(sizeof(unsigned long long) == sizeof(uint64_t), "Bad ULL length");
        if ((count & 0xFFFFFUL) == 0) {
            s[0] = 0;
            s[1] = 0;
            _rdrand64_step((unsigned long long*)(s + 0));
            _rdrand64_step((unsigned long long*)(s + 1));
        }
        ++count;

        Produce();
    }

public:
    Rander() : count(0), bits(0) {}

    uint32_t GetBits(int bits_) {
        if (bits_ > bits) {
            Step();
        }

        uint32_t ret = next & ((1UL << bits_) - 1);
        next >>= bits_;
        bits -= bits_;
        return ret;
    }

    uint32_t GetInt(uint32_t range, int bits_) {
        do {
            uint32_t r = GetBits(bits_);
            if (r < range) return r;
        } while(true);
    }
};

static Rander rander;

int randtrans() {
    do {
        int x = rander.GetBits(10);
        if ((x >> 5) == 0) continue;
        if ((x >> 5) == 31) continue;
        if ((x & 31 & (~(x >> 5))) == 0) continue;
        return x;
    } while(true);
}

struct charset {
    char map1[128];
    char map2[32];

    charset(const char* str) {
        for (int i = 0; i < 128; i++) {
            map1[i] = -1;
        }
        for (int i = 0; i < 32; i++) {
            map2[i] = str[i];
            map1[str[i] + 0] = i;
        }
    }

    friend bool operator!=(const charset& o1, const charset& o2) {
        return memcmp(o1.map2, o2.map2, 32) != 0;
    }

    std::string print() const {
        std::string str;
        for (int i = 0; i < 32; i++) {
            str += (map2[i]);
        }
        return str;
    }

    static int transform(int trans, int val) {
        int mask = trans >> 5;
        int masked = val & mask;
        int need = trans & 31 & mask;
        if (masked == need) {
            val ^= trans & 31 & ~mask;
        }
        return val;
    }

    charset(const charset& other, int trans1) {
        for (int i = 0; i < 128; i++) {
            int p = other.map1[i];
            if (p != -1) {
                p = transform(trans1, p);
                map2[p] = i;
            }
            map1[i] = p;
        }
    }

    charset(const charset& other, int trans1, int trans2) {
        for (int i = 0; i < 128; i++) {
            int p = other.map1[i];
            if (p != -1) {
                p = transform(trans1, p);
                p = transform(trans2, p);
                map2[p] = i;
            }
            map1[i] = p;
        }
    }

    charset(const charset& other, int trans1, int trans2, int trans3) {
        for (int i = 0; i < 128; i++) {
            int p = other.map1[i];
            if (p != -1) {
                p = transform(trans1, p);
                p = transform(trans2, p);
                p = transform(trans3, p);
                map2[p] = i;
            }
            map1[i] = p;
        }
    }

    charset(const charset& other, int trans1, int trans2, int trans3, int trans4) {
        for (int i = 0; i < 128; i++) {
            int p = other.map1[i];
            if (p != -1) {
                p = transform(trans1, p);
                p = transform(trans2, p);
                p = transform(trans3, p);
                p = transform(trans4, p);
                map2[p] = i;
            }
            map1[i] = p;
        }
    }

    charset(const charset& other, int trans1, int trans2, int trans3, int trans4, int trans5) {
        for (int i = 0; i < 128; i++) {
            int p = other.map1[i];
            if (p != -1) {
                p = transform(trans1, p);
                p = transform(trans2, p);
                p = transform(trans3, p);
                p = transform(trans4, p);
                p = transform(trans5, p);
                map2[p] = i;
            }
            map1[i] = p;
        }
    }


    uint64_t score() const {
        uint64_t ret[5] = {0};
        for (int p = 0; p < 7; ++p) {
            for (int s = 0; s < 32; ++s) {
                int e = 1 << p;
                int d = map2[s];
                int md = map1[e ^ d];
                if (md != -1) {
                    int mw = __builtin_popcount(s ^ md);
                    ++ret[mw - 1];
                }
            }
        }

        uint64_t bad = 0;
        for (size_t p = 0 ; p < sizeof(badpairs)/sizeof(badpairs[0]); p++) {
            int v1 = map1[badpairs[p][0] + 0], v2 = map1[badpairs[p][1] + 0];
            if (v1 == -1 || v2 == -1 || __builtin_popcount(v1 ^ v2) == 1) {
                bad += badpairs[p][2];
            }
        }

        return ret[0] << 32 | bad << 40 | ret[1] << 24 | ret[2] << 16 | ret[3] << 8 | ret[4];
    }
};

int main(int argc, char** argv) {
    setbuf(stdout, NULL);

    int bestbad = 0;
    for (int d = 0; d < 2; d++) {
    for (int i1 = 0; i1 < 33; i1++) {
      char c1 = i1 < 10 ? '0' + i1 : 'a' + (i1 - 10);
      for (int i2 = i1 + 1; i2 < 34; i2++) {
        char c2 = i2 < 10 ? '0' + i2 : 'a' + (i2 - 10);
        for (int i3 = i2 + 1; i3 < 35; i3++) {
          char c3 = i3 < 10 ? '0' + i3 : 'a' + (i3 - 10);
          for (int i4 = i3 + 1; i4 < 36; i4++) {
            char c4 = i4 < 10 ? '0' + i4 : 'a' + (i4 - 10);
            int totalbad = 0;
            for (size_t i = 0; i < sizeof(badpairs) / sizeof(badpairs[0]); i++) {
              char b1 = badpairs[i][0], b2 = badpairs[i][1];
              if (c1 == b1 || c2 == b1 || c3 == b1 || c4 == b1 || c1 == b2 || c2 == b2 || c3 == b2 || c4 == b2) {
                totalbad += badpairs[i][2];
              }
            }
            if (d == 0 && totalbad > bestbad) {
              bestbad = totalbad;
            } else if (d == 1 && totalbad == bestbad) {
              printf("Best bad: %c %c %c %c (%x)\n", c1, c2, c3, c4, totalbad);
            }
          }
        }
      }
    }
    }

    int totalbad = 0;
    for (size_t i = 0; i < sizeof(badpairs) / sizeof(badpairs[0]); i++) {
        totalbad += badpairs[i][2];
    }
    printf("Total bad: %x\n", totalbad);

    charset start(argv[1]);
    uint64_t val = start.score();
    printf("Initial: %s (%llx)\n", start.print().c_str(), (unsigned long long)val);

    std::set<std::string> vars;
    for (int i = 0; i < 32; i++) {
        charset m(start, i);
        vars.insert(m.print());
    }
    for (const auto& x : vars) {
        printf("Var: %s\n", x.c_str());
    }

    bool cont = true;
    while (cont) {
        cont = false;
        for (int x = 0; x < 1048476; ++x) {
            charset m(start, x & 1023, x >> 10);
            uint64_t nval = m.score();
            if (nval > val) {
                start = m;
                val = nval;
                printf("Update %i: %s (%llx)\n", x, start.print().c_str(), , (unsigned long long)val);
                cont = true;
                break;
            }
        }
        if (cont) continue;
        cont = false;

        for (int c = 3; c <= 6 && !cont; c++) {
            for (int i = 0; i < 1048576; i++) {
                charset m(start, randtrans(), randtrans(), randtrans());
                for (int it = 3; it < c; it++) {
                    charset x(m, randtrans());
                    m = x;
                }
                uint64_t nval = m.score();
                if (nval > val) {
                    start = m;
                    val = nval;
                    printf("Update R%i: %s (%llx)\n", c, start.print().c_str(), (unsigned long long)val);
                    cont = true;
                    break;
                }
            }
        }
        if (cont) continue;
        cont = false;

/*
        charset m(start, randtrans(), randtrans(), randtrans());
        uint64_t nval = m.score();
        printf("Update X3: ");
        start.print();
        printf(" (%llx)\n", (unsigned long long)val);
        start = m;
        val = nval;
        cont = true;
*/
        while (!cont) {
            for (int i = 0; i < 65536; i++) {
                charset m(start, 0);
                while (true) {
                    charset m2(start, randtrans(), randtrans(), randtrans());
                    if (m2 != start) {
                        m = m2;
                        break;
                    }
                }
                uint64_t nval = m.score();
                if (nval >= val) {
                    start = m;
                    val = nval;
                    printf("Update E3: %s (%llx)\n", start.print().c_str(), (unsigned long long)val);
                    cont = true;
                    break;
                }
            }
            if (cont) break;
            cont = false;

            for (int i = 0; i < 65536 * 2; i++) {
                charset m(start, 0);
                while (true) {
                    charset m2(start, randtrans(), randtrans(), randtrans(), randtrans());
                    if (m2 != start) {
                        m = m2;
                        break;
                    }
                }
                uint64_t nval = m.score();
                if (nval >= val) {
                    start = m;
                    val = nval;
                    printf("Update E4: %s (%llx)\n", start.print().c_str(), (unsigned long long)val);
                    cont = true;
                    break;
                }
            }
            if (cont) break;
            cont = false;

            for (int i = 0; i < 65536 * 4; i++) {
                charset m(start, 0);
                while (true) {
                    charset m2(start, randtrans(), randtrans(), randtrans(), randtrans(), randtrans());
                    if (m2 != start) {
                        m = m2;
                        break;
                    }
                }
                uint64_t nval = m.score();
                if (nval >= val) {
                    start = m;
                    val = nval;
                    printf("Update E5: %s (%llx)\n", start.print().c_str(), (unsigned long long)val);
                    cont = true;
                    break;
                }
            }
            if (cont) break;
            cont = false;
        }
    }
    return 0;
}

