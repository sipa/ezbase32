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
#include <x86intrin.h>

#include "tinyformat.h"

#define DEGREE 13
#define LENGTH 115
#define ERRORS 4
#define MAX_DEFICIENCY 2
#define THREADS 1

#define MIN_FACTOR_DEGREE 1

static inline uint32_t rdrand() {
    uint32_t ret;
    unsigned char ok;
    __asm__ volatile(".byte 0x0f, 0xc7, 0xf0; setc %1" : "=a"(ret), "=q"(ok) :: "cc");
    assert(ok);
    return ret;
}

static constexpr int exptable[32] = {1,2,4,8,16,9,18,13,26,29,19,15,30,21,3,6,12,24,25,27,31,23,7,14,28,17,11,22,5,10,20,1};
static constexpr int logtable[32] = {-1,0,1,14,2,28,15,22,3,5,29,26,16,7,23,11,4,25,6,10,30,13,27,21,17,18,8,19,24,9,12,20};

static constexpr uint8_t gf32_mul(uint8_t x, uint8_t y) { return (x == 0 || y == 0) ? 0 : exptable[(logtable[x] + logtable[y]) % 31]; }

struct MulTable {
    uint8_t table[32][32];
    uint8_t invtable[32];

    MulTable() {
        for (unsigned int i = 0; i < 32; ++i) {
            for (uint64_t j = 0; j < 32; ++j) {
                uint8_t p = gf32_mul(i, j);
                table[i][j] = p;
                if (p == 1) {
                    invtable[i] = j;
                }
            }
        }
        for (unsigned int i = 1; i < 32; ++i) {
            assert(invtable[invtable[i]] == i);
            assert(table[i][invtable[i]] == 1);
        }
    }
};

static const MulTable multable;

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


struct Char {
    static constexpr int Count = 1;
    uint8_t v;

    Char() : v(0) {}

    void Set(int pos, uint8_t val) {
        v = val;
    }

    uint8_t Get(int pos) const { return v; }
    bool IsZero() const { return v == 0; }
    bool IsZero(int pos) const { return v == 0; }
    bool IsOne(int pos) const { return v == 1; }
    int Cmp(Char x) const { if (v < x.v) return -1; if (v > x.v) return 1; return 0; }
    void Xor(Char x) { v ^= x.v; }

    void SubMul(Char x, const uint8_t val) {
        v ^= multable.table[val][x.v];
    }
};

template<typename T, int N>
struct Pack {
    static constexpr int Count = N;
    static_assert(std::is_unsigned<T>::value, "T must be unsigned");
    static_assert(std::numeric_limits<T>::max() >> (5 * N - 1), "T not large enough");

    T v;
    Pack() : v(0) {}
    void Set(int pos, uint8_t val) {
        v = (v & ~(((T)31) << (5 * pos))) | (((T)val) << (5 * pos));
    }

    uint8_t Get(int pos) const {
        return (v >> (5 * pos)) & 31;
    }

    bool IsZero() const { return v == 0; }
    bool IsZero(int pos) const { return ((v >> (5 * pos)) & 0x1f) == 0; }
    bool IsOne(int pos) const { return ((v >> (5 * pos)) & 0x1f) == 1; }
    int Cmp(Pack x) const { if (v < x.v) return -1; if (v > x.v) return 1; return 0; }
    void Xor(Pack x) { v ^= x.v; }

    void SubMul(Pack x, const uint8_t val) {
        for (int i = 0; i < Count; ++i) {
            v ^= (((T)(multable.table[val][(x.v >> (5 * i)) & 0x1f])) << (5 * i));
        }
    }
};

template<typename T, int N>
struct PackX {
    static constexpr int Count = N;
    static_assert(std::is_unsigned<T>::value, "T must be unsigned");
    static_assert(std::numeric_limits<T>::max() >> (5 * N - 1), "T not large enough");
    static_assert(N <= 7, "N is too large");

    T v;
    PackX() : v(0) {}
    void Set(int pos, uint8_t val) {
        v = (v & ~(((T)31) << (5 * pos))) | (((T)val) << (5 * pos));
    }

    uint8_t Get(int pos) const {
        return (v >> (5 * pos)) & 31;
    }

    bool IsZero() const { return v == 0; }
    bool IsZero(int pos) const { return ((v >> (5 * pos)) & 0x1f) == 0; }
    bool IsOne(int pos) const { return ((v >> (5 * pos)) & 0x1f) == 1; }
    int Cmp(PackX x) const { if (v < x.v) return -1; if (v > x.v) return 1; return 0; }
    void Xor(PackX x) { v ^= x.v; }

    void SubMul(PackX x, const uint8_t val) {
        constexpr uint64_t mask = 0x7C3E1F0F87C3E1FUL;
        __m128i xv = _mm_set_epi64x(0, _pdep_u64(x.v, mask));
        __m128i vv = _mm_set_epi64x(0, val);
        __m128i rv = _mm_clmulepi64_si128(xv, vv, 0);
        v ^= _pext_u64(reduce4(_mm_extract_epi64(rv, 0)), mask);
    }
};

template<typename I, int N>
struct BitsX {
    static constexpr int Count = N;
    I v;

    BitsX() : v(0) {
        static_assert(std::is_unsigned<I>::value, "I must be unsigned");
        static_assert(std::numeric_limits<I>::max() >> (9 * N - 1), "I must be large enough");
    }

    void Set(int pos, uint8_t val) {
        v = (v & ~(((I)31) << (9 * pos))) | (((I)val) << (9 * pos));
    }

    uint8_t Get(int pos) const {
        return (v >> (9 * pos)) & 0x1f;
    }

    bool IsZero() const { return v == 0; }
    bool IsZero(int pos) const { return ((v >> (9 * pos)) & 0x1f) == 0; }
    bool IsOne(int pos) const { return ((v >> (9 * pos)) & 0x1f) == 1; }

    int Cmp(BitsX x) const { if (v < x.v) return -1; if (v > x.v) return 1; return 0; }

    void Xor(BitsX x) { v ^= x.v; }

    void SubMul(BitsX x, uint8_t val) {
        __m128i xv = _mm_set_epi64x(0, x.v);
        __m128i vv = _mm_set_epi64x(0, val);
        __m128i rv = _mm_clmulepi64_si128(xv, vv, 0);
        v ^= reduce4(_mm_extract_epi64(rv, 0));
    }
};

template<typename A, typename B>
struct Cat {
    static constexpr int Count = A::Count + B::Count;
    A a;
    B b;
    void Set(int pos, uint8_t val) {
        if (pos < A::Count) {
            a.Set(pos, val);
        } else {
            b.Set(pos - A::Count, val);
        }
    }
    uint8_t Get(int pos) const {
        if (pos < A::Count) {
            return a.Get(pos);
        } else {
            return b.Get(pos - A::Count);
        }
    }
    bool IsZero() const { return a.IsZero() && b.IsZero(); }
    bool IsZero(int pos) const { if (pos < A::Count) return a.IsZero(pos); else return b.IsZero(pos - A::Count); }
    bool IsOne(int pos) const { if (pos < A::Count) return a.IsOne(pos); else return b.IsOne(pos - A::Count); }
    int Cmp(const Cat& x) const { int ret = a.Cmp(x.a); if (ret) return ret; return b.Cmp(x.b); }
    void Xor(const Cat& x) { a.Xor(x.a); b.Xor(x.b); }
    void SubMul(const Cat& x, uint8_t val) { a.SubMul(x.a, val); b.SubMul(x.b, val); }
};

template<typename A, int N>
struct Mult {
    static constexpr int Count = A::Count * N;
    A a[N];
    void Set(int pos, uint8_t val) {
        a[pos / A::Count].Set(pos % A::Count, val);
    }
    uint8_t Get(int pos) const {
        return a[pos / A::Count].Get(pos % A::Count);
    }
    bool IsZero() const {
        for (int i = 0; i < N; ++i) {
            if (!a[i].IsZero()) return false;
        }
        return true;
    }
    bool IsZero(int pos) const { return a[pos / A::Count].IsZero(pos % A::Count); }
    bool IsOne(int pos) const { return a[pos / A::Count].IsOne(pos % A::Count); }
    int Cmp(const Mult& x) const {
        for (int i = 0; i < N; ++i) {
            int ret = a[i].Cmp(x.a[i]);
            if (ret) return ret;
        }
        return false;
    }
    void Xor(const Mult& x) {
        for (int i = 0; i < N; ++i) {
            a[i].Xor(x.a[i]);
        }
    }
    void SubMul(const Mult& x, uint8_t val) {
        for (int i = 0; i < N; ++i) {
            a[i].SubMul(x.a[i], val);
        }
    }
};

#ifdef CONF
template<int N> struct Config { typedef CONF Elem; };
#else
template<int N> struct Config;
template<> struct Config<1> { typedef Char Elem; };
template<> struct Config<2> { typedef Cat<Char,Char> Elem; };
template<> struct Config<3> { typedef Mult<Char,3> Elem; };
template<> struct Config<4> { typedef Mult<Char,4> Elem; };
template<> struct Config<5> { typedef PackX<uint32_t, 5> Elem; };
template<> struct Config<6> { typedef PackX<uint32_t, 6> Elem; };
template<> struct Config<8> { typedef Mult<Char,8> Elem; };
template<> struct Config<9> { typedef Mult<Char,9> Elem; };
template<> struct Config<12> { typedef Pack<uint64_t, 12> Elem; };
template<> struct Config<13> { typedef Mult<Pack<uint16_t, 3>,5> Elem; };
#endif

template<int N>
class Vector {
    typedef typename Config<N>::Elem Elem;
    static constexpr int Count = Elem::Count;
    static_assert(N <= Count, "Internal data type insufficient for vector size");

    Elem d;

public:
    uint8_t operator[](int a) const { return d.Get(a); }
    void Set(int a, uint8_t val) { d.Set(a, val); }

    bool IsZero() const { return d.IsZero(); }
    bool IsZero(int pos) const { return d.IsZero(pos); }
    bool IsOne(int pos) const { return d.IsOne(pos); }
    bool operator==(const Vector& a) const { return d.Cmp(a.d) == 0; }
    bool operator<(const Vector& a) const { return d.Cmp(a.d) < 0; }
    Vector<N>& operator+=(const Vector& a) { d.Xor(a.d); return *this; }
    void SubMul(const Vector& a, uint8_t v) { d.SubMul(a.d, v); }
    Vector<N>& operator*=(uint8_t a) { Elem t; t.SubMul(d, a); d = t; return *this; }

    void PolyMulXMod(const Vector& mod) {
        auto ptr = multable.table[(*this)[N - 1]];
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
        ret ^= multable.table[a[i]][b[i]];
    }
    return ret;
}

template<int A, int B>
Vector<A+B> Concat(const Vector<A>& a, const Vector<B>& b)
{
    Vector<A+B> ret;
    for (int i = 0; i < A; ++i) {
        ret.Set(i, a[i]);
    }
    for (int i = 0; i < B; ++i) {
        ret.Set(i + A, b[i]);
    }
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
            if (res[r].IsZero(c)) {
                int r2 = r + 1;
                while (r2 < R) {
                    if (!res[r2].IsZero(c)) {
                        res.SwapRows(r2, r);
                        SwapRows(r2, r);
                        break;
                    }
                    ++r2;
                }
                if (r2 < R) {
                    assert(!res[r].IsZero(c));
                } else {
                    assert(res[r].IsZero(c));
                    continue;
                }
            }
            uint8_t i = multable.invtable[res[r][c]];
            assert(!res[r].IsZero(c));
            res.MulRow(r, i);
            MulRow(r, i);
            assert(res[r].IsOne(c));
            for (int r2 = 0; r2 < R; ++r2) {
                if (r2 != r) {
                    uint8_t i = res[r2][c];
                    res.SubMulRow(r2, r, i);
                    SubMulRow(r2, r, i);
                    assert(res[r2].IsZero(c));
                }
            }
            ++r;
        }
        int ret = r;
        for (int c = 0; c < C; ++c) {
            int solvepos = -1;
            for (int r = 0; r < R; ++r) {
                if (!res[r].IsZero(c) && !res[r].IsOne(c)) {
                    solvepos = -2;
                    break;
                }
                if (solvepos == -1 && res[r].IsOne(c)) {
                    solvepos = r;
                } else if (solvepos != -1 && !res[r].IsZero(c)) {
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

#ifndef BENCH
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
typedef std::vector<Vector<DEGREE-ERRORS>> extbasis_type;
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

static void ExpandSolutions(result_type& res, const extbasis_type& extbasis, const psol_type& psol, const Vector<ERRORS>& fault, bool allzerobefore) {
    res.reserve(psol.size() * 2);
    for (const auto& ps : psol) {
        Vector<ERRORS> base_errors;
        uint64_t solcount = BaseSolution(base_errors, ps.second, fault);
        for (uint64_t sol = 0; sol < solcount; ++sol) {
            Vector<ERRORS> ext_errors = ExtSolution(base_errors, ps.second, sol);

            // Filter out duplicates for fewer errors than max
            bool consec = true;
            bool ok = true;
            for (int i = 0; i < ERRORS; ++i) {
                consec = consec && (ps.first[i] == (i ? ps.first[i - 1] + 1 : 0));
                if (ext_errors.IsZero(i) && !consec) {
                    ok = false;
                    break;
                }
            }
            if (!ok) continue;

            // Compute the full fault and verify it
            Vector<DEGREE-ERRORS> extbigfault;
            for (int i = 0; i < ERRORS; ++i) {
                extbigfault.SubMul(extbasis[ps.first[i]], ext_errors[i]);
            }
            Vector<DEGREE> bigfault = Concat(fault, extbigfault);

            if (allzerobefore) {
                for (int i = ERRORS; i < DEGREE; ++i) {
                    if (!bigfault.IsZero(i)) {
                        ok = (bigfault.IsOne(i));
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

bool RecurseShortFaults(int pos, bool allzerobefore, Vector<ERRORS>& fault, const psol_type& psol, const extbasis_type& extbasis, int part, uint64_t hash, LockedErrCount& err) {
    if (pos == ERRORS) {
        if ((hash % ((uint64_t)THREADS)) != (uint64_t)part) return true;
        if (err.cleanup.load(std::memory_order_relaxed)) {
            return false;
        }

        ErrCount local_err;
        result_type res;
        ExpandSolutions(res, extbasis, psol, fault, allzerobefore);
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
        if (!RecurseShortFaults(pos + 1, allzerobefore && fault.IsZero(pos), fault, psol, extbasis, part, hash * 9672876866715837617ULL + fault[pos], err)) {
            return false;
        }
    }
    return true;
}

static const char* charset = "0123456789ABCDEFGHIJKLMNOPQRSTUV";

void run_thread(const psol_type* partials, const extbasis_type* extbasis, int part, LockedErrCount* locs) {
    Vector<ERRORS> faults;
    (void)RecurseShortFaults(0, true, faults, *partials, *extbasis, part, 0, *locs);
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

bool testalot(const basis_type* basis, const extbasis_type* extbasis, ErrCount* res) {
    psol_type partials;
    std::array<int, ERRORS> pos;
    RecursePositions(0, ERRORS, 0, LENGTH, pos, partials, *basis);

    LockedErrCount ret;
#if THREADS > 1
    std::vector<std::thread> t;
    for (int part = 0; part < THREADS - 1; ++part) {
        t.emplace_back(&run_thread, &partials, extbasis, part, &ret);
    }
#endif
    run_thread(&partials, extbasis, THREADS - 1, &ret);
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

template<typename I, typename J>
void ElemTest() {
    for (int ip = 0; ip < I::Count; ++ip) {
        for (int jp = 0; jp < J::Count; ++jp) {
            for (int x = 0; x < 32; ++x) {
                for (int y = 0; y < 32; ++y) {
                    for (int z = 0; z < 32; ++z) {
                        I i1, i2;
                        J j1, j2;
                        i1.Set(ip, x);
                        j1.Set(jp, x);
                        assert(i1.IsZero() == (x == 0));
                        assert(j1.IsZero() == (x == 0));
                        assert(i1.IsZero(ip) == (x == 0));
                        assert(j1.IsZero(jp) == (x == 0));
                        assert(i1.IsOne(ip) == (x == 1));
                        assert(j1.IsOne(jp) == (x == 1));
                        assert(j1.Get(jp) == i1.Get(ip));
                        i2.Set(ip, y);
                        j2.Set(jp, y);
                        assert(i2.IsZero() == (y == 0));
                        assert(j2.IsZero() == (y == 0));
                        assert(i2.IsZero(ip) == (y == 0));
                        assert(j2.IsZero(jp) == (y == 0));
                        assert(i2.IsOne(ip) == (y == 1));
                        assert(j2.IsOne(jp) == (y == 1));
                        assert(j2.Get(jp) == i2.Get(ip));
                        i1.Xor(i2);
                        j1.Xor(j2);
                        assert((i1.Cmp(i2) == 0) == (j1.Cmp(j2) == 0));
                        assert(i1.Get(ip) == j1.Get(jp));
                        assert(i1.Get(ip) == (x ^ y));
                        assert(i2.IsZero() == j2.IsZero());
                        assert(i2.IsZero(ip) == j2.IsZero(jp));
                        assert(i2.IsOne(ip) == j2.IsOne(jp));
                        i2.SubMul(i1, z);
                        j2.SubMul(j1, z);
                        assert(i1.Get(ip) == j1.Get(jp));
                        assert(i2.Get(ip) == j2.Get(jp));
                        assert((i1.Cmp(i2)) == (j1.Cmp(j2)));
                        assert(i2.IsZero() == j2.IsZero());
                        assert(i2.IsZero(ip) == j2.IsZero(jp));
                        assert(i2.IsOne(ip) == j2.IsOne(jp));
                    }
                }
            }
        }
    }
}


int main(int argc, char** argv) {
    setbuf(stdout, NULL);
    ElemTest<Char, Config<DEGREE>::Elem>();
    ElemTest<Char, Config<DEGREE-ERRORS>::Elem>();
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
    extbasis_type extbasis;
    basis.resize(LENGTH);
    extbasis.resize(LENGTH);
    Vector<DEGREE> x;
    x.Set(0, 1);

    Matrix<DEGREE, DEGREE> rand;
    Transform<DEGREE, DEGREE> rand_trans;
    while(true) {
        Matrix<DEGREE, DEGREE> randi, res;
        for (int i = 0; i < DEGREE; ++i) {
            for (int j = 0; j < DEGREE; ++j) {
                rand[i].Set(j, rdrand() & 0x1F);
            }
        }
        randi = rand;
        int rank = randi.Invert(res);
        if (rank == DEGREE) { rand_trans = Transform<DEGREE, DEGREE>(rand); break; }
    }

    assert(LENGTH >= DEGREE);

    for (int i = 0; i < LENGTH; ++i) {
        Vector<DEGREE> base = x.Low<DEGREE>();
        basis[i] = rand_trans.Apply(base);
        extbasis[i] = basis[i].High<DEGREE-ERRORS>();
        x.PolyMulXMod(gen);
    }

    ErrCount locs;
    locs.errors = 0;
    testalot(&basis, &extbasis, &locs);
    if (locs.errors == 0) {
        show_stats(locs);
    }

    return 0;
}
#else
int main() {
    printf("Size: %i\n", sizeof(Vector<NUM>));
    Matrix<NUM,NUM> mat;
    do {
        for (int i = 0; i < NUM; ++i) {
            for (int j = 0; j < NUM; ++j) {
                mat[i].Set(j, rdrand() & 0x1f);
            }
        }
        Matrix<NUM,NUM> mati = mat, res;
        if (mati.Invert(res) == NUM) break;
    } while(true);
    Transform<NUM,NUM> trans(mat);
    Vector<NUM> vec;
    for (int i = 0; i < NUM; ++i) {
        vec.Set(i, rdrand() & 0x1f);
    }
    for (int i = 0; i * NUM * NUM < 1000000000; ++i) {
        vec = trans.Apply(vec);
    }
    for (int i = 0; i < NUM; ++i) {
        printf("%i ", vec[i]);
    }
    printf("\n");
    return 0;
}

#endif
