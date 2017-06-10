#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <algorithm>
#include <vector>
#include <array>
#include <tuple>
#include <math.h>
#include <thread>
#include <mutex>
#include <unistd.h>
#include <atomic>
#include <condition_variable>
#include <set>

#include "tinyformat.h"

#define RANDOMIZE_ON_ERROR 1
#define NUMCODES 1
#define POLYLEN 12
#define DEGREE 12
#define LENGTH 65
#define ERRORS 4
#define MAX_DEFICIENCY 2
#define THREADS 8

static inline uint64_t reduce2(uint64_t x) {
    uint64_t high = (x & 0xE0E0E0E0E0E0E0E0ULL) >> 5;
    uint64_t low = x & 0x1F1F1F1F1F1F1F1FULL;
    return low ^ high ^ (high << 3);
}

static inline uint64_t reduce3(uint64_t x) {
    uint64_t high1 = (x & 0x6060606060606060ULL) >> 5;
    uint64_t high2 = (x & 0x8080808080808080ULL) >> 7;
    uint64_t low = x & 0x1f1f1f1f1f1f1f1full;
    return low ^ high1 ^ ((high1 ^ high2) << 3) ^ high2 ^ (high2 << 2);
}

static uint64_t mul0(uint64_t x) { return 0; }
static uint64_t mul1(uint64_t x) { return x; }
static uint64_t mul2(uint64_t x) { return reduce2(x << 1); }
static uint64_t mul3(uint64_t x) { return reduce2(x ^ (x << 1)); }
static uint64_t mul4(uint64_t x) { return reduce2(x << 2); }
static uint64_t mul5(uint64_t x) { return reduce2(x ^ (x << 2)); }
static uint64_t mul6(uint64_t x) { return reduce2((x << 1) ^ (x << 2)); }
static uint64_t mul7(uint64_t x) { return reduce2(x ^ (x << 1) ^ (x << 2)); }
static uint64_t mul8(uint64_t x) { return reduce3(x << 3); }
static uint64_t mul9(uint64_t x) { return reduce3(x ^ (x << 3)); }
static uint64_t mul10(uint64_t x) { return reduce3((x << 1) ^ (x << 3)); }
static uint64_t mul11(uint64_t x) { return reduce3(x ^ (x << 1) ^ (x << 3)); }
static uint64_t mul12(uint64_t x) { return reduce3((x << 2) ^ (x << 3)); }
static uint64_t mul13(uint64_t x) { return reduce3(x ^ (x << 2) ^ (x << 3)); }
static uint64_t mul14(uint64_t x) { return reduce3((x << 1) ^ (x << 2) ^ (x << 3)); }
static uint64_t mul15(uint64_t x) { return reduce3(x ^ (x << 1) ^ (x << 2) ^ (x << 3)); }
static uint64_t mul16(uint64_t x) { return mul4(mul4(x)); }
static uint64_t mul17(uint64_t x) { uint64_t m4 = mul4(x); return reduce2(x ^ (m4 << 2)); }
static uint64_t mul18(uint64_t x) { uint64_t m4 = mul4(x); return reduce2((x << 1) ^ (m4 << 2)); }
static uint64_t mul19(uint64_t x) { uint64_t m4 = mul4(x); return reduce2(x ^ (x << 1) ^ (m4 << 2)); }
static uint64_t mul20(uint64_t x) { uint64_t m4 = mul4(x); return reduce2((x ^ m4) << 2); }
static uint64_t mul21(uint64_t x) { uint64_t m4 = mul4(x); return reduce2(x ^ ((x ^ m4) << 2)); }
static uint64_t mul22(uint64_t x) { uint64_t m4 = mul4(x); return reduce2((x << 1) ^ ((x ^ m4) << 2)); }
static uint64_t mul23(uint64_t x) { uint64_t m4 = mul4(x); return reduce2(x ^ (x << 1) ^ ((x ^ m4) << 2)); }
static uint64_t mul24(uint64_t x) { uint64_t m4 = mul4(x); return reduce2((m4 << 2) ^ (m4 << 1)); }
static uint64_t mul25(uint64_t x) { uint64_t m4 = mul4(x); return reduce2(x ^ (m4 << 2) ^ (m4 << 1)); }
static uint64_t mul26(uint64_t x) { uint64_t m4 = mul4(x); return reduce2((m4 << 2) ^ ((x ^ m4) << 1)); }
static uint64_t mul27(uint64_t x) { uint64_t m4 = mul4(x); return reduce2(x ^ (m4 << 2) ^ ((x ^ m4) << 1)); }
static uint64_t mul28(uint64_t x) { uint64_t m4 = mul4(x); return reduce2(m4 ^ (m4 << 1) ^ (m4 << 2)); }
static uint64_t mul29(uint64_t x) { uint64_t m4 = mul4(x); return reduce2(x ^ m4 ^ (m4 << 1) ^ (m4 << 2)); }
static uint64_t mul30(uint64_t x) { uint64_t m4 = mul4(x); return reduce2(m4 ^ ((x ^ m4) << 1) ^ (m4 << 2)); }
static uint64_t mul31(uint64_t x) { uint64_t m4 = mul4(x); return reduce2(x ^ m4 ^ ((x ^ m4) << 1) ^ (m4 << 2)); }

typedef uint64_t (*mulfun)(uint64_t);

static const mulfun mulfuns[32] = {
    mul0,mul1,mul2,mul3,mul4,mul5,mul6,mul7,
    mul8,mul9,mul10,mul11,mul12,mul13,mul14,mul15,
    mul16,mul17,mul18,mul19,mul20,mul21,mul22,mul23,
    mul24,mul25,mul26,mul27,mul28,mul29,mul30,mul31
};

class MulTable {
    uint8_t table[32][32];
    uint8_t divtable[32][32];

public:
    MulTable() {
        for (int i = 0; i < 32; ++i) {
            for (int j = 0; j < 32; ++j) {
                int p = mulfuns[i](j);
                table[i][j] = p;
                divtable[p][i] = j;
            }
        }
    }

    uint8_t mul(uint8_t a, uint8_t b) const {
        return table[a][b];
    }

    uint8_t div(uint8_t a, uint8_t b) const {
        return divtable[a][b];
    }

    const uint8_t* ptr(uint8_t a) const {
        return table[a];
    }
};

static const MulTable multable;

template<int N>
class Vector {
    uint8_t d[N];

public:
    Vector() {
        memset(d, 0, sizeof(d));
    }

    uint8_t& operator[](int a) { return d[a]; }
    const uint8_t& operator[](int a) const { return d[a]; }

    bool IsZero() const {
        for (int i = 0; i < N; ++i) {
            if (d[i] != 0) return false;
        }
        return true;
    }

    bool operator==(const Vector& a) const {
        for (int i = 0; i < N; ++i) {
            if (d[i] != a[i]) return false;
        }
        return true;
    }

    bool operator<(const Vector& a) const {
        for (int i = 0; i < N; ++i) {
            if (d[i] > a[i]) return false;
            if (d[i] < a[i]) return true;
        }
        return false;
    }

    Vector<N>& operator+=(const Vector<N>& a) {
        for (int i = 0; i < N; ++i) {
            d[i] ^= a[i];
        }
        return *this;
    }

    void SubMul(const Vector<N>& a, uint8_t v) {
        auto ptr = multable.ptr(v);
        for (int i = 0; i < N; ++i) {
            d[i] ^= ptr[a[i]];
        }
    }

    Vector<N>& operator*=(uint8_t a) {
        auto ptr = multable.ptr(a);
        for (int i = 0; i < N; ++i) {
            d[i] = ptr[d[i]];
        }
        return *this;
    }

    int Weight() const {
        int ret = 0;
        for (int i = 0; i < N; ++i) {
            ret += (d[i] != 0);
        }
        return ret;
    }

    void PolyMulXMod(const Vector<N>& mod) {
        auto ptr = multable.ptr(d[N - 1]);
        uint8_t over = 0;
        for (int i = 0; i < N; ++i) {
            uint8_t nover = d[i];
            d[i] = over ^ ptr[mod[i]];
            over = nover;
        }
    }

    template<int A>
    Vector<A> Low() const {
        Vector<A> ret;
        memcpy(&ret[0], d, A);
        return ret;
    }

    template<int A>
    Vector<A> High() const {
        Vector<A> ret;
        memcpy(&ret[0], d + N - A, A);
    }
};

template<int A>
uint8_t Multiply(const Vector<A>& a, const Vector<A>& b) {
    uint8_t ret = 0;
    for (int i = 0; i < A; ++i) {
        ret ^= multable.mul(a[i], b[i]);
    }
    return ret;
}

template<int A, int B>
Vector<A+B> Concat(const Vector<A>& a, const Vector<B>& b)
{
    Vector<A+B> ret;
    memcpy(&ret[0], &a[0], A);
    memcpy(&ret[A], &b[0], B);
    return ret;
}

template<int R, int C>
class Matrix {
    Vector<C> row[R];

public:
    const Vector<C>& operator[](int r) const { return row[r]; }
    Vector<C>& operator[](int r) { return row[r]; }

    void MakeIdentity() {
        static_assert(R == C, "Matrix must be square");
        for (int r = 0; r < R; ++r) {
            for (int c = 0; c < R; ++c) {
                row[r][c] = (r == c);
            }
        }
    }

    Matrix<C,R> Transpose() const {
        Matrix<C,R> ret;
        for (int r = 0; r < R; ++r) {
            for (int c = 0; c < C; ++c) {
                ret[c][r] = row[r][c];
            }
        }
        return ret;
    }

    void SwapRows(int r1, int r2) {
        std::swap(row[r1], row[r2]);
    }

    void MulRow(int r, uint8_t v) {
        row[r] *= v;
    }

    void SubMulRow(int r, int ro, uint8_t v) {
        row[r].SubMul(row[ro], v);
    }

    Matrix<R,C>& operator+=(const Matrix<R,C>& a) {
        for (int r = 0; r < R; ++r) {
            row[r] += a[r];
        }
        return *this;
    }

    Matrix<R,C>& operator*=(uint8_t v) {
        for (int r = 0; r < R; ++r) {
            row[r] *= v;
        }
        return *this;
    }


    int Invert(Matrix<R,R>& res) {
        static_assert(R == C, "Matrix must be square");
        res = *this;
        MakeIdentity();
        int r = 0;
        for (int c = 0; c < C; ++c) {
            if (res[r][c] == 0) {
                int r2 = r + 1;
                while (r2 < R) {
                    if (res[r2][c] != 0) {
                        res.SwapRows(r2, r);
                        SwapRows(r2, r);
                        break;
                    }
                    ++r2;
                }
                if (r2 < R) {
                    assert(res[r][c] != 0);
                } else {
                    assert(res[r][c] == 0);
                    continue;
                }
            }
            uint8_t i = multable.div(1, res[r][c]);
            res.MulRow(r, i);
            MulRow(r, i);
            assert(res[r][c] == 1);
            for (int r2 = 0; r2 < R; ++r2) {
                if (r2 != r) {
                    uint8_t i = res[r2][c];
                    res.SubMulRow(r2, r, i);
                    SubMulRow(r2, r, i);
                    assert(res[r2][c] == 0);
                }
            }
            ++r;
        }
        int ret = r;
        for (int c = 0; c < C; ++c) {
            int solvepos = -1;
            for (int r = 0; r < R; ++r) {
                if (res[r][c] != 0 && res[r][c] != 1) {
                    solvepos = -2;
                    break;
                }
                if (solvepos == -1 && res[r][c] == 1) {
                    solvepos = r;
                } else if (solvepos != -1 && res[r][c] != 0) {
                    solvepos = -2;
                    break;
                }
            }
            if (solvepos >= 0 && solvepos != c) {
                res.SwapRows(solvepos, c);
                SwapRows(solvepos, c);
            }
        }
        return ret;
    }
};

template<int RM, int CM>
Vector<RM> Multiply(const Matrix<RM,CM>& m, const Vector<CM>& v) {
    Vector<RM> ret;
    for (int t = 0; t < CM; ++t) {
        auto ptr = multable.ptr(v[t]);
        for (int r = 0; r < RM; ++r) {
            ret[r] ^= ptr[m[r][t]];
        }
    }
    return ret;
}

template<int N>
struct PartialSolution {
    Vector<N> constraints[MAX_DEFICIENCY];
    Matrix<N, N> solutions;
    Vector<N> freedom[MAX_DEFICIENCY];
    int deficiency;
};

template<int N>
PartialSolution<N> PartialSolve(const Matrix<N, N>& equations) {
    PartialSolution<N> ret;
    Matrix<N, N> inverse = equations;
    Matrix<N, N> residual;
    int rank = inverse.Invert(residual);
    assert(rank >= N - MAX_DEFICIENCY);
    int def_num = 0;
    for (int r = 0; r < N; ++r) {
        if (residual[r].IsZero()) {
            ret.constraints[def_num] = inverse[r];
            for (int c = 0; c < N; ++c) {
                ret.freedom[def_num][c] = residual[c][r] ^ (r == c);
            }
            ++def_num;
        } else {
            for (int c = 0; c < N; ++c) {
                ret.solutions[r][c] = inverse[r][c];
            }
        }
    }
    ret.deficiency = def_num;
    assert(def_num == N - rank);
    return ret;
}

template<int N>
uint64_t BaseSolution(Vector<N>& base_sol, const PartialSolution<N>& partial, const Vector<N>& knowns) {
    for (int def_num = 0; def_num < partial.deficiency; ++def_num) {
        uint8_t check = Multiply(partial.constraints[def_num], knowns);
        if (check) return 0;
    }
    base_sol = Multiply(partial.solutions, knowns);
    return ((uint64_t)1) << (5 * partial.deficiency);
}

template<int N>
Vector<N> ExtSolution(const Vector<N>& base_sol, const PartialSolution<N>& partial, uint64_t solnum) {
    Vector<N> ext_sol = base_sol;
    for (int def_num = 0; def_num < partial.deficiency; ++def_num) {
        ext_sol.SubMul(partial.freedom[def_num], (solnum >> (5 * def_num)) & 0x1F);
    }
    return ext_sol;
}

struct Result {
    Vector<DEGREE> fault;
    int min_pos;
    int max_pos;
    int num_err;
    int pos[ERRORS];

/*    Vector<ERRORS> err;
    std::array<int, ERRORS> pos;*/

    Result(const Vector<DEGREE>& fault_, int min_pos_, int max_pos_, int num_err_) : fault(fault_), min_pos(min_pos_), max_pos(max_pos_), num_err(num_err_) {}

    bool operator==(const Result& result) const {
        if (!(fault == result.fault)) return false;
        if (min_pos != result.min_pos) return false;
        if (max_pos != result.max_pos) return false;
        if (num_err != result.num_err) return false;
        return true;
    }

    bool operator<(const Result& result) const {
        if (fault < result.fault) return true;
        if (result.fault < fault) return false;
        if (min_pos < result.min_pos) return true;
        if (result.min_pos < min_pos) return false;
        if (max_pos < result.max_pos) return true;
        if (result.max_pos < max_pos) return false;
        if (num_err < result.num_err) return true;
        if (result.num_err < num_err) return false;
        return false;
    }
};

typedef std::vector<Result> result_type;
typedef std::vector<Vector<DEGREE>> basis_type;
typedef std::vector<std::pair<std::array<int, ERRORS>, PartialSolution<ERRORS>>> psol_type;

struct ErrCount {
    uint64_t count[2*ERRORS+1][LENGTH + 1];
    uint64_t total;

    ErrCount() : total(0) {
        memset(count, 0, sizeof(count));
    }

    void Inc(int errors, int length) {
        assert(length < LENGTH + 1);
        assert(errors < 2*ERRORS+1);
        for (int len = length; len < LENGTH + 1; ++len) {
            ++count[errors][len];
        }
    }

    void operator+=(const ErrCount& e) {
        for (int c = 0; c < 2*ERRORS+1; ++c) {
            for (int l = 0; l < LENGTH + 1; ++l) {
                count[c][l] += e.count[c][l];
            }
        }
        total += e.total;
    }
};

void RecursePositions(int idx, int min, std::array<int, ERRORS>& pos, psol_type& psol, const basis_type& basis) {
    if (idx == ERRORS) {
        Matrix<ERRORS,ERRORS> rr;
        for (int i = 0; i < ERRORS; ++i) {
            rr[i] = basis[pos[i]].Low<ERRORS>();
        }
        rr = rr.Transpose();
        psol.push_back(std::make_pair(pos, PartialSolve(rr)));
/*        if (psol.back().second.deficiency) {
            printf("Deficient %i: %i %i %i\n", psol.back().second.deficiency, pos[0], pos[1], pos[2]);
        }*/
        return;
    }
    for (pos[idx] = min; pos[idx] < LENGTH; ++pos[idx]) {
        RecursePositions(idx + 1, pos[idx] + 1, pos, psol, basis);
    }
}

long double Combination(int k, int n) {
    long double num = 1.0;
    long double den = 1.0;
    if (n - k < k) k = n - k;
    for (int i = 1; i <= k; i++) {
        num *= (n - i + 1);
        den *= i;
    }
    return num / den;
}

static int require_len = 0, require_err = 0;

static std::mutex results_mutex;
static ErrCount results;

struct ErrorLocations {
    int errors;
    int pos[ERRORS * 2];
    uint64_t progress;
};

static std::string code;

static long double total_comb() {
    long double ret = 0;
    for (int i = 1; i <= ERRORS; ++i) {
        ret += Combination(i, LENGTH) * powl(31, i - 1);
    }
    return ret;
}

bool RecurseShortFaults(int pos, bool allzerobefore, Vector<ERRORS>& fault, const psol_type& psol, const basis_type& basis, int part, uint64_t hash, std::atomic<bool>& abort, ErrorLocations& err) {
    if (pos == ERRORS) {
        if ((hash % ((uint64_t)THREADS)) != (uint64_t)part) return true;
        ErrCount errcount;
        result_type res;
        for (const auto& ps : psol) {
            if (abort.load(std::memory_order_relaxed)) return false;
            Vector<ERRORS> base_errors;
            uint64_t solcount = BaseSolution(base_errors, ps.second, fault);
            for (uint64_t sol = 0; sol < solcount; ++sol) {
                Vector<ERRORS> ext_errors = ExtSolution(base_errors, ps.second, sol);

                /* Filter out duplicates for fewer errors than max */
                bool consec = true;
                bool ok = true;
                for (int i = 0; i < ERRORS; ++i) {
                    consec = consec && (ps.first[i] == (i ? ps.first[i - 1] + 1 : 0));
                    if (ext_errors[i] == 0 && !consec) {
                        ok = false;
                        break;
                    }
                }
                if (!ok) continue;


                // Compute the full fault and verify it
                Vector<DEGREE> bigfault;
                for (int i = 0; i < ERRORS; ++i) {
                    bigfault.SubMul(basis[ps.first[i]], ext_errors[i]);
                }
                for (int i = 0; i < ERRORS; ++i) {
                    assert(bigfault[i] == fault[i]);
                }

                if (allzerobefore) {
                    for (int i = ERRORS; i < DEGREE; ++i) {
                        if (bigfault[i] != 0) {
                            ok = (bigfault[i] == 1);
                            break;
                        }
                    }
                    if (!ok) continue;
                }

                int num_error = 0;
                int min_pos = LENGTH;
                int max_pos = 0;
                for (int i = 0; i < ERRORS; ++i) {
                    if (ext_errors[i]) {
                        ++num_error;
                        min_pos = std::min(min_pos, ps.first[i]);
                        max_pos = std::max(max_pos, ps.first[i]);
                    }
                }

                if (num_error == 0) continue;
                ++errcount.total;
                res.emplace_back(bigfault, min_pos, max_pos, num_error);
                int nn = 0;
                for (int i = 0; i < ERRORS; ++i) {
                    if (ext_errors[i]) {
                        res.back().pos[nn] = ps.first[i];
                        ++nn;
                    }
                }
            }
        }
        std::sort(res.begin(), res.end());
        for (size_t pos = 0; pos < res.size(); ++pos) {
            if (abort.load(std::memory_order_relaxed)) return false;
            size_t pos1 = pos;
            auto &key = res[pos].fault;
            while (pos + 1 < res.size() && key == res[pos + 1].fault) {
                ++pos;
            }
            for (size_t posA = pos1; posA + 1 <= pos; ++posA) {
                for (size_t posB = posA + 1; posB <= pos; ++posB) {
                    if (res[posA].num_err <= res[posB].num_err && res[posA].num_err + 1 >= res[posB].num_err && res[posA].max_pos < res[posB].min_pos) {
                        int total_err = res[posA].num_err + res[posB].num_err;
                        const int length = res[posB].max_pos;
                        if (length + 1 - res[posA].min_pos <= require_len && total_err <= require_err) {
                            std::string line;
                            std::unique_lock<std::mutex> lock(results_mutex);
                            abort.store(true, std::memory_order_relaxed);
                            line += strprintf("%s: %i errors in a window of size %i: ", code.c_str(), total_err, length + 1 - res[posA].min_pos);
                            for (int nn = 0; nn < res[posA].num_err; ++nn) {
                                line += strprintf("%i ", res[posA].pos[nn]);
                                err.pos[err.errors++] = res[posA].pos[nn];
                            }
                            for (int nn = 0; nn < res[posB].num_err; ++nn) {
                                line += strprintf("%i ", res[posB].pos[nn]);
                                err.pos[err.errors++] = res[posB].pos[nn];
                            }
                            {
                                results += errcount;
                                line += strprintf(" # %Lg%% done\n", results.total / total_comb() * 100.0L);
                            }
                            printf("%s", line.c_str());
                            return false;
                        }
                        errcount.Inc(total_err, length + 1);
                    }
                }
            }
        }
        {
            std::unique_lock<std::mutex> lock(results_mutex);
            results += errcount;
        }
        return true;
    }
    int max = allzerobefore ? 2 : 32;
    int rrr = allzerobefore ? (random() & 1) : (random() & 0x1f);
    for (int x = 0; x < max; ++x) {
        fault[pos] = x ^ rrr;
        if (!RecurseShortFaults(pos + 1, allzerobefore && fault[pos] == 0, fault, psol, basis, part, hash * 9672876866715837601ULL + fault[pos], abort, err)) {
            return false;
        }
    }
    return true;
}

static const char* charset = "0123456789ABCDEFGHIJKLMNOPQRSTUV";

void run_thread(const psol_type* partials, const basis_type* basis, int part, std::atomic<bool>* aborter, ErrorLocations* locs) {
    Vector<ERRORS> faults;
    (void)RecurseShortFaults(0, true, faults, *partials, *basis, part, 0, *aborter, *locs);
}

using namespace std::chrono_literals;


void stat_thread() {
    while (true) {
        bool quit = false;
        std::this_thread::sleep_for(60000ms);
        std::unique_lock<std::mutex> lock(results_mutex);
        static const long double denom = 1.0L / total_comb();
        long double frac = results.total * denom;
        for (int l = 1; l <= LENGTH; ++l) {
            std::string line;
            line += strprintf("%s % 4i", code.c_str(), l);
            for (int e = 1; e <= ERRORS*2; ++e) {
                line += strprintf(" % 19.15f", (double)(results.count[e][l] / (Combination(e, l) * frac * powl(31.0, e - 1)) * powl(32.0, DEGREE)));
            }
            line += strprintf("  # %Lg%% done\n", frac * 100.0L);
            printf("%s", line.c_str());
        }
        if (quit) return;
    }
}

bool testalot(const psol_type* partials, const basis_type* basis, ErrorLocations* locs) {
    {
        std::unique_lock<std::mutex> lock(results_mutex);
        results = ErrCount();
    }
    std::atomic<bool> aborter{false};
    std::vector<std::thread> t;
    for (int part = 0; part < THREADS; ++part) {
        t.emplace_back(&run_thread, partials, basis, part, &aborter, locs);
    }
    for (int part = 0; part < THREADS; ++part) {
        t[part].join();
    }
    locs->progress = results.total;
    return !aborter.load();
}

std::string namecode(const Vector<DEGREE>& v) {
   std::string ret;
   ret.resize(DEGREE);
   for (int d = 0; d < DEGREE; ++d) {
       ret[d] = charset[v[d]];
   }
   return ret;
}


int main(int argc, char** argv) {
    srandom(time(NULL) * 7 + getpid() * 127);
    setbuf(stdout, NULL);
    Vector<POLYLEN> gen[NUMCODES];
    if (argc < 2 || strlen(argv[1]) != (POLYLEN + 1) * NUMCODES - 1) {
        fprintf(stderr, "Usage: %s GEN%i,...x%i\n", argv[0], POLYLEN, NUMCODES);
        return 1;
    }
    for (int c = 0; c < NUMCODES; ++c) {
        if (c && argv[1][(POLYLEN + 1) * c - 1] != ',') {
            fprintf(stderr, "Comma expected: '%c'\n", argv[1][(POLYLEN + 1) * c - 1]);
            return 1;
        }
        for (int i = 0; i < POLYLEN; ++i) {
            const char *ptr = strchr(charset, toupper(argv[1][(POLYLEN + 1) * c + POLYLEN - 1 - i]));
            if (ptr == nullptr) {
                fprintf(stderr, "Unknown character '%c'\n", argv[1][(POLYLEN + 1) * c + POLYLEN - 1 - i]);
                return 1;
            }
            gen[c][i] = ptr - charset;
        }
    }
    if (argc >= 3) { require_err = strtoul(argv[2], NULL, 0); }
    if (argc >= 4) { require_len = strtoul(argv[3], NULL, 0); }

    basis_type basis;
    basis.resize(LENGTH);
    Vector<POLYLEN> x;
    x[0] = 1;

    Matrix<DEGREE, DEGREE> rand;
    while(true) {
        Matrix<DEGREE, DEGREE> randi, res;
        for (int i = 0; i < DEGREE; ++i) {
            for (int j = 0; j < DEGREE; ++j) {
                rand[i][j] = random() & 0x1F;
            }
        }
        randi = rand;
        int rank = randi.Invert(res);
        if (rank == DEGREE) break;
    }

    std::vector<std::string> codes;
    codes.resize(LENGTH);
    code = std::string(argv[1]);

    for (int i = 0; i < LENGTH; ++i) {
        Vector<DEGREE> base = x.Low<DEGREE>();
        codes[i] = namecode(base);
        basis[i] = Multiply(rand, base);
        x.PolyMulXMod(gen[i % NUMCODES]);
    }

    std::thread th(&stat_thread);

    std::string best_code = code;
    std::vector<std::string> best_codes = codes;
    basis_type best_basis = basis;
    ErrorLocations best_locations;
    best_locations.progress = 0;

    do {
        srandom(random() * 13 + time(NULL) * 17 + getpid() * 19);
        psol_type partials;
        {
            std::array<int, ERRORS> pos;
            RecursePositions(0, 0, pos, partials, basis);
        }

        ErrorLocations locs;
        locs.errors = 0;
        testalot(&partials, &basis, &locs);
        if (locs.errors == 0) break;

        if (locs.progress <= best_locations.progress) {
            code = best_code;
            codes = best_codes;
            basis = best_basis;
            locs = best_locations;
        } else {
            best_code = code;
            best_codes = codes;
            best_basis = basis;
            best_locations = locs;
            std::string str = strprintf("%s: new best (%Lg%%)\n", best_code, best_locations.progress / total_comb() * 100.0L);
            printf("%s", str.c_str());
        }

        Vector<DEGREE> rv;
        for (int j = 0; j < DEGREE; ++j) {
            rv[j] = random() & 0x1F;
        }

        int pos;
        do {
            pos = locs.pos[random() % locs.errors];
        } while (pos < DEGREE);

        codes[pos] = namecode(rv);
        basis[pos] = Multiply(rand, rv);

        std::string rescode = codes[DEGREE];
        for (int l = DEGREE + 1; l < LENGTH; ++l) {
            rescode = rescode + "|" + codes[l];
        }
        code = rescode;

#if RANDOMIZE_ON_ERROR
    } while(true);
#else
    } while(false);
#endif
    return 0;
}
