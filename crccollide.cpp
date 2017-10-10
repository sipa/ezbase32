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

#define DEGREE 6
#define LENGTH 90
#define ERRORS 3
#define MAX_DEFICIENCY 2
#define THREADS 8

#define MIN_FACTOR_DEGREE 1

static inline uint32_t rdrand() {
    uint32_t ret;
    unsigned char ok;
    __asm__ volatile(".byte 0x0f, 0xc7, 0xf0; setc %1" : "=a"(ret), "=q"(ok) :: "cc");
    assert(ok);
    return ret;
}

static inline uint64_t reduce2(uint64_t x) {
    uint64_t high = (x & 0x180C06030180C060ULL) >> 5;
    uint64_t low = x & 0x7C3E1F0F87C3E1FULL;
    return low ^ high ^ (high << 3);
}

static inline uint64_t reduce4(uint64_t x) {
    uint64_t high1 = (x & 0x180C06030180C060ULL) >> 5;
    uint64_t high2 = (x & 0x6030180C06030180ULL) >> 7;
    uint64_t low = x & 0x7C3E1F0F87C3E1FULL;
    uint64_t tmp = high1 ^ high2;
    return low ^ tmp ^ (tmp << 3) ^ (high2 << 2);
}

static const int exptable[32] = {1,2,4,8,16,9,18,13,26,29,19,15,30,21,3,6,12,24,25,27,31,23,7,14,28,17,11,22,5,10,20,1};
static const int logtable[32] = {-1,0,1,14,2,28,15,22,3,5,29,26,16,7,23,11,4,25,6,10,30,13,27,21,17,18,8,19,24,9,12,20};

static uint64_t mul0(uint64_t x) { return 0; }
static uint64_t mul1(uint64_t x) { return x; }
static uint64_t mul2(uint64_t x) { return reduce2(x << 1); }
static uint64_t mul3(uint64_t x) { return reduce2(x ^ (x << 1)); }
static uint64_t mul4(uint64_t x) { return reduce2(x << 2); }
static uint64_t mul5(uint64_t x) { return reduce2(x ^ (x << 2)); }
static uint64_t mul6(uint64_t x) { return reduce2((x << 1) ^ (x << 2)); }
static uint64_t mul7(uint64_t x) { return reduce2(x ^ (x << 1) ^ (x << 2)); }
static uint64_t mul8(uint64_t x) { return reduce4(x << 3); }
static uint64_t mul9(uint64_t x) { return reduce4(x ^ (x << 3)); }
static uint64_t mul10(uint64_t x) { return reduce4((x << 1) ^ (x << 3)); }
static uint64_t mul11(uint64_t x) { return reduce4(x ^ (x << 1) ^ (x << 3)); }
static uint64_t mul12(uint64_t x) { return reduce4((x << 2) ^ (x << 3)); }
static uint64_t mul13(uint64_t x) { return reduce4(x ^ (x << 2) ^ (x << 3)); }
static uint64_t mul14(uint64_t x) { return reduce4((x << 1) ^ (x << 2) ^ (x << 3)); }
static uint64_t mul15(uint64_t x) { return reduce4(x ^ (x << 1) ^ (x << 2) ^ (x << 3)); }
static uint64_t mul16(uint64_t x) { return reduce4(x << 4); }
static uint64_t mul17(uint64_t x) { return reduce4((x << 4) ^ x); }
static uint64_t mul18(uint64_t x) { return reduce4((x << 4) ^ (x << 1)); }
static uint64_t mul19(uint64_t x) { return reduce4((x << 4) ^ (x << 1) ^ x); }
static uint64_t mul20(uint64_t x) { return reduce4((x << 4) ^ (x << 2)); }
static uint64_t mul21(uint64_t x) { return reduce4((x << 4) ^ (x << 2) ^ x); }
static uint64_t mul22(uint64_t x) { return reduce4((x << 4) ^ (x << 2) ^ (x << 1)); }
static uint64_t mul23(uint64_t x) { return reduce4((x << 4) ^ (x << 2) ^ (x << 1) ^ x); }
static uint64_t mul24(uint64_t x) { return reduce4((x << 4) ^ (x << 3)); }
static uint64_t mul25(uint64_t x) { return reduce4((x << 4) ^ (x << 3) ^ x); }
static uint64_t mul26(uint64_t x) { return reduce4((x << 4) ^ (x << 3) ^ (x << 1)); }
static uint64_t mul27(uint64_t x) { return reduce4((x << 4) ^ (x << 3) ^ (x << 1) ^ x); }
static uint64_t mul28(uint64_t x) { return reduce4((x << 4) ^ (x << 3) ^ (x << 2)); }
static uint64_t mul29(uint64_t x) { return reduce4((x << 4) ^ (x << 3) ^ (x << 2) ^ x); }
static uint64_t mul30(uint64_t x) { return reduce4((x << 4) ^ (x << 3) ^ (x << 2) ^ (x << 1)); }
static uint64_t mul31(uint64_t x) { return reduce4((x << 4) ^ (x << 3) ^ (x << 2) ^ (x << 1) ^ x); }


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
        for (unsigned int i = 0; i < 32; ++i) {
            for (uint64_t j = 0; j < 32; ++j) {
                uint64_t p = (i == 0 || j == 0) ? 0 : exptable[(logtable[i] + logtable[j]) % 31];
                assert(mulfuns[i](j) == p);
                assert(mulfuns[i](j << 9) == p << 9);
                assert(mulfuns[i](j << 18) == p << 18);
                assert(mulfuns[i](j << 27) == p << 27);
                assert(mulfuns[i](j << 36) == p << 36);
                assert(mulfuns[i](j << 45) == p << 45);
                assert(mulfuns[i](j << 54) == p << 54);
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
    static constexpr int L = (N + 6) / 7;
    uint64_t l[L];

public:
    Vector() {
        memset(l, 0, sizeof(l));
    }

    uint8_t operator[](int a) const { return (l[a / 7] >> (9 * (a % 7))) & 0x1F; }

    void Set(int a, uint8_t val) { l[a / 7] = (l[a / 7] & ~(((uint64_t)0x1f) << (9 * (a % 7)))) | (((uint64_t)val) << (9 * (a % 7))); }

    bool IsZero() const {
        for (int i = 0; i < L; ++i) {
            if (l[i] != 0) return false;
        }
        return true;
    }

    bool operator==(const Vector& a) const {
        for (int i = 0; i < L; ++i) {
            if (l[i] != a.l[i]) return false;
        }
        return true;
    }

    bool operator<(const Vector& a) const {
        for (int i = 0; i < L; ++i) {
            if (l[i] > a.l[i]) return false;
            if (l[i] < a.l[i]) return true;
        }
        return false;
    }

    Vector<N>& operator+=(const Vector<N>& a) {
        for (int i = 0; i < L; ++i) {
            l[i] ^= a.l[i];
        }
        return *this;
    }

    void SubMul(const Vector<N>& a, mulfun fn) {
        for (int i = 0; i < L; ++i) {
            l[i] ^= fn(a.l[i]);
        }
    }

    void SubMul(const Vector<N>& a, uint8_t v) {
        auto fn = mulfuns[v];
        for (int i = 0; i < L; ++i) {
            l[i] ^= fn(a.l[i]);
        }
    }

    Vector<N>& operator*=(uint8_t a) {
        auto fn = mulfuns[a];
        for (int i = 0; i < L; ++i) {
            l[i] = fn(l[i]);
        }
        return *this;
    }

    int Weight() const {
        int ret = 0;
        for (int i = 0; i < N; ++i) {
            ret += ((*this)[i] != 0);
        }
        return ret;
    }

    void PolyMulXMod(const Vector<N>& mod) {
        auto ptr = multable.ptr((*this)[N - 1]);
        uint8_t over = 0;
        for (int i = 0; i < N; ++i) {
            uint8_t nover = (*this)[i];
            Set(i, over ^ ptr[mod[i]]);
            over = nover;
        }
    }

    template<int A>
    Vector<A> Low() const {
        static_assert(A <= N, "ow");
        Vector<A> ret;
        for (int i = 0; i < A; ++i) {
            ret.Set(i, (*this)[i]);
        }
        return ret;
    }

    template<int A>
    Vector<A> High() const {
        static_assert(A <= N, "ow");
        Vector<A> ret;
        for (int i = 0; i < A; ++i) {
            ret.Set(i, (*this)[i + N - A]);
        }
        return ret;
    }
};

template<int A, int B>
uint8_t Divides(const Vector<A>& a, const Vector<B>& b) {
    Vector<A> m;
    m[0] = 1;
    for (int d = B - 1; d >= 0; --d) {
        m.PolyMulXMod(a);
        m[0] ^= b[d];
    }
    return m.IsZero();
}

template<int A, int B>
bool AnyFactor_(int pos, Vector<A>& a, const Vector<B>& b) {
    if (pos == A) {
        return Divides(a, b);
    }
    for (int i = 0; i < 32; ++i) {
        a[pos] = i;
        if (AnyFactor_(pos + 1, a, b)) return true;
    }
    return false;
}

template<int A, int B>
bool AnyFactor(const Vector<B>& b) {
    Vector<A> a;
    return AnyFactor_(0, a, b);
}

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
                row[r].Set(c, r == c);
            }
        }
    }

    Matrix<C,R> Transpose() const {
        Matrix<C,R> ret;
        for (int r = 0; r < R; ++r) {
            for (int c = 0; c < C; ++c) {
                ret[c].Set(r, row[r][c]);
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

template<int R, int C>
class Transform {
    Vector<R> column[C];

public:
    Transform() {}

    Transform(const Matrix<R, C>& mat) {
        for (int r = 0; r < R; ++r) {
            for (int c = 0; c < C; ++c) {
                column[c].Set(r, mat[r][c]);
            }
        }
    }

    Vector<R> Apply(const Vector<C>& in) const {
        Vector<R> out;
        for (int c = 0; c < C; ++c) {
            out.SubMul(column[c], in[c]);
        }
        return out;
    }
};

template<int RM, int CM>
Vector<RM> Multiply(const Matrix<RM,CM>& m, const Vector<CM>& v) {
    Vector<RM> ret;
    for (int t = 0; t < CM; ++t) {
        auto ptr = multable.ptr(v[t]);
        for (int r = 0; r < RM; ++r) {
            ret.Set(r, ret[r] ^ ptr[m[r][t]]);
        }
    }
    return ret;
}

template<int N>
struct PartialSolution {
    Vector<N> constraints[MAX_DEFICIENCY];
    Transform<N, N> solutions;
    Vector<N> freedom[MAX_DEFICIENCY];
    int deficiency;
};

template<int N>
PartialSolution<N> PartialSolve(const Matrix<N, N>& equations) {
    PartialSolution<N> ret;
    Matrix<N, N> inverse = equations;
    Matrix<N, N> residual;
    Matrix<N, N> solutions;
    int rank = inverse.Invert(residual);
    assert(rank >= N - MAX_DEFICIENCY);
    int def_num = 0;
    for (int r = 0; r < N; ++r) {
        if (residual[r].IsZero()) {
            ret.constraints[def_num] = inverse[r];
            for (int c = 0; c < N; ++c) {
                ret.freedom[def_num].Set(c, residual[c][r] ^ (r == c));
            }
            ++def_num;
        } else {
            for (int c = 0; c < N; ++c) {
                solutions[r].Set(c, inverse[r][c]);
            }
        }
    }
    ret.solutions = Transform<N, N>(solutions);
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
    base_sol = partial.solutions.Apply(knowns);
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

bool ComparePsol(const std::pair<std::array<int, ERRORS>, PartialSolution<ERRORS>>& a, const std::pair<std::array<int, ERRORS>, PartialSolution<ERRORS>>& b) {
    if (a.first == b.first) return false;
    if (a.second.deficiency > b.second.deficiency) return true;
    if (a.second.deficiency < b.second.deficiency) return false;
    for (int p = ERRORS - 1; p >= 0; --p) {
        if (a.first[p] < b.first[p]) return true;
        if (a.first[p] > b.first[p]) return false;
    }
    return false;
}

struct ErrCount {
    uint64_t count[2*ERRORS+1][LENGTH + 1];
    uint64_t total = 0;

    int errors = 0;
    int pos[ERRORS * 2];

    ErrCount() {
        memset(count, 0, sizeof(count));
    }

    void Inc(int errors, int length) {
        assert(length < LENGTH + 1);
        assert(errors < 2*ERRORS+1);
        for (int len = length; len < LENGTH + 1; ++len) {
            ++count[errors][len];
        }
    }

    uint64_t Update(const ErrCount& e) {
        for (int c = 0; c < 2*ERRORS+1; ++c) {
            for (int l = 0; l < LENGTH + 1; ++l) {
                count[c][l] += e.count[c][l];
            }
        }
        total += e.total;
        if (e.errors && !errors) {
            errors = e.errors;
            memcpy(pos, e.pos, sizeof(pos));
        }
        return total;
    }
};

struct LockedErrCount : public ErrCount {
#if THREADS > 1
    std::mutex lock;
#endif
    std::atomic<bool> cleanup{false};

    uint64_t Update(const ErrCount& e) {
#if THREADS > 1
        std::unique_lock<std::mutex> cs(lock);
#endif
        return ErrCount::Update(e);
    }
};

void RecursePositions(int idx, int errors, int min, int max, std::array<int, ERRORS>& pos, psol_type& psol, const basis_type& basis) {
    if (idx == errors) {
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
    for (pos[idx] = min; pos[idx] < max; ++pos[idx]) {
        RecursePositions(idx + 1, errors, pos[idx] + 1, max, pos, psol, basis);
    }
}

constexpr long double Combination(int k, int n) {
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

static std::string code;

static constexpr long double total_comb() {
    long double ret = 0;
    for (int i = 1; i <= ERRORS; ++i) {
        ret += Combination(i, LENGTH) * powl(31, i - 1);
    }
    return ret;
}

static constexpr long double denom = 1.0L / total_comb();

static void ExpandSolutions(result_type& res, const basis_type& basis, const psol_type& psol, const Vector<ERRORS>& fault, bool allzerobefore) {
    res.reserve(psol.size() * 2);
    for (const auto& ps : psol) {
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
}

static bool CompareResultPointer(const Result* a, const Result* b) {
    return (*a) < (*b);
}

bool RecurseShortFaults(int pos, bool allzerobefore, Vector<ERRORS>& fault, const psol_type& psol, const basis_type& basis, int part, uint64_t hash, LockedErrCount& err) {
    if (pos == ERRORS) {
        if ((hash % ((uint64_t)THREADS)) != (uint64_t)part) return true;
        if (err.cleanup.load(std::memory_order_relaxed)) {
            return false;
        }

        ErrCount local_err;
        result_type res;
        ExpandSolutions(res, basis, psol, fault, allzerobefore);
        local_err.total = res.size();

        std::vector<const Result*> sres[1024];
        for (int s = 0; s < 1024; ++s) {
            sres[s].reserve(res.size() / 300);
        }
        for (size_t i = 0; i < res.size(); ++i) {
            sres[res[i].fault[DEGREE - 1] + 32 * res[i].fault[DEGREE - 2]].push_back(&res[i]);
        }
        for (int s = 0; s < 1024; ++s) {
            std::sort(sres[s].begin(), sres[s].end(), CompareResultPointer);
            for (size_t pos = 0; pos < sres[s].size(); ++pos) {
                size_t pos1 = pos;
                auto &key = sres[s][pos]->fault;
                while (pos + 1 < sres[s].size() && key == sres[s][pos + 1]->fault) {
                    ++pos;
                }
                for (size_t posA = pos1; posA + 1 <= pos; ++posA) {
                    for (size_t posB = posA + 1; posB <= pos; ++posB) {
                        if (sres[s][posA]->num_err <= sres[s][posB]->num_err && sres[s][posA]->num_err + 1 >= sres[s][posB]->num_err && sres[s][posA]->max_pos < sres[s][posB]->min_pos) {
                            int total_err = sres[s][posA]->num_err + sres[s][posB]->num_err;
                            const int length = sres[s][posB]->max_pos;
                            if (length + 1 - sres[s][posA]->min_pos <= require_len && total_err <= require_err) {
                                err.cleanup.store(true, std::memory_order_relaxed);
                                for (int nn = 0; nn < sres[s][posA]->num_err; ++nn) {
                                    local_err.pos[local_err.errors++] = sres[s][posA]->pos[nn];
                                }
                                for (int nn = 0; nn < sres[s][posB]->num_err; ++nn) {
                                    local_err.pos[local_err.errors++] = sres[s][posB]->pos[nn];
                                }
                                err.Update(local_err);
                                return false;
                            }
                            local_err.Inc(total_err, length + 1);
                        }
                    }
                }
            }
        }
        err.Update(local_err);
        return true;
    }
    int max = allzerobefore ? 2 : 32;
    int rrr = allzerobefore ? 0 : (rdrand() & 0x1f);
    for (int x = 0; x < max; ++x) {
        fault.Set(pos, x ^ rrr);
        if (!RecurseShortFaults(pos + 1, allzerobefore && fault[pos] == 0, fault, psol, basis, part, hash * 9672876866715837617ULL + fault[pos], err)) {
            return false;
        }
    }
    return true;
}

static const char* charset = "0123456789ABCDEFGHIJKLMNOPQRSTUV";

void run_thread(const psol_type* partials, const basis_type* basis, int part, LockedErrCount* locs) {
    Vector<ERRORS> faults;
    (void)RecurseShortFaults(0, true, faults, *partials, *basis, part, 0, *locs);
}

using namespace std::chrono_literals;


void show_stats(const ErrCount& results) {
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
}

bool testalot(const basis_type* basis, ErrCount* res) {
    psol_type partials;
    std::array<int, ERRORS> pos;
    RecursePositions(0, ERRORS, 0, LENGTH, pos, partials, *basis);

    LockedErrCount ret;
#if THREADS > 1
    std::vector<std::thread> t;
    for (int part = 0; part < THREADS - 1; ++part) {
        t.emplace_back(&run_thread, &partials, basis, part, &ret);
    }
#endif
    run_thread(&partials, basis, THREADS - 1, &ret);
#if THREADS > 1
    for (int part = 0; part < THREADS - 1; ++part) {
        t[part].join();
    }
#endif
    *res = ret;
    if (res->errors) {
        std::string line;
        line += strprintf("%s: %i errors in a window of size %i: ", code.c_str(), res->errors, res->pos[res->errors - 1] + 1 - res->pos[0]);
        for (int e = 0; e < res->errors; ++e) {
            line += strprintf("%i ", res->pos[e]);
        }
        line += strprintf(" # %Lg%% done\n", res->total / total_comb() * 100.0L);
        printf("%s", line.c_str());
        return false;
    }
    return true;
}

std::string namecode(const Vector<DEGREE>& v) {
   std::string ret;
   ret.resize(DEGREE);
   for (int d = 0; d < DEGREE; ++d) {
       ret[d] = charset[v[DEGREE - 1 - d]];
   }
   return ret;
}

int main(int argc, char** argv) {
    setbuf(stdout, NULL);
    Vector<DEGREE> gen;
    if (argc < 2 || (strlen(argv[1]) != DEGREE && strlen(argv[1]) != 0)) {
        fprintf(stderr, "Usage: %s GEN%i\n", argv[0], DEGREE);
        return 1;
    }
    if (strlen(argv[1]) == 0) {
        do {
            for (int i = 0; i < DEGREE; ++i) {
                gen.Set(i, rdrand() & 0x1f);
            }
            code = namecode(gen);
            if (gen[0] == 0) { printf("%s: divisible by x\n", code.c_str()); continue; }
#if MIN_FACTOR_DEGREE > 1
            if (AnyFactor<1, DEGREE>(gen)) { printf("%s: degree 1 factor\n", code.c_str()); continue; }
#endif
#if MIN_FACTOR_DEGREE > 2
            if (AnyFactor<2, DEGREE>(gen)) { printf("%s: degree 2 factor\n", code.c_str()); continue; }
#endif
#if MIN_FACTOR_DEGREE > 3
            if (AnyFactor<3, DEGREE>(gen)) { printf("%s: degree 3 factor\n", code.c_str()); continue; }
#endif
#if MIN_FACTOR_DEGREE > 4
            if (AnyFactor<4, DEGREE>(gen)) { printf("%s: degree 4 factor\n", code.c_str()); continue; }
#endif
#if MIN_FACTOR_DEGREE > 5
            if (AnyFactor<5, DEGREE>(gen)) { printf("%s: degree 5 factor\n", code.c_str()); continue; }
#endif
#if MIN_FACTOR_DEGREE > 6
            if (AnyFactor<6, DEGREE>(gen)) { printf("%s: degree 6 factor\n", code.c_str()); continue; }
#endif
            break;
        } while(true);
    } else {
        for (int i = 0; i < DEGREE; ++i) {
            const char *ptr = strchr(charset, toupper(argv[1][DEGREE - 1 - i]));
            if (ptr == nullptr) {
                fprintf(stderr, "Unknown character '%c'\n", argv[1][DEGREE - 1 - i]);
                return 1;
            }
            gen.Set(i, ptr - charset);
        }
    }
    code = namecode(gen);
    printf("%s: starting\n", code.c_str());

    if (argc >= 3) { require_err = strtoul(argv[2], NULL, 0); }
    if (argc >= 4) { require_len = strtoul(argv[3], NULL, 0); }

    basis_type basis;
    basis.resize(LENGTH);
    Vector<DEGREE> x;
    x.Set(0, 1);

    Matrix<DEGREE, DEGREE> rand;
    while(true) {
        Matrix<DEGREE, DEGREE> randi, res;
        for (int i = 0; i < DEGREE; ++i) {
            for (int j = 0; j < DEGREE; ++j) {
                rand[i].Set(j, rdrand() & 0x1F);
            }
        }
        randi = rand;
        int rank = randi.Invert(res);
        if (rank == DEGREE) break;
    }

    assert(LENGTH >= DEGREE);

    for (int i = 0; i < LENGTH; ++i) {
        Vector<DEGREE> base = x.Low<DEGREE>();
        basis[i] = Multiply(rand, base);
        x.PolyMulXMod(gen);
    }

    ErrCount locs;
    locs.errors = 0;
    testalot(&basis, &locs);
    if (locs.errors == 0) {
        show_stats(locs);
    }

    return 0;
}
