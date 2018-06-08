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

static constexpr int DEGREE = 6;
static constexpr int LENGTH = 200;
static constexpr int ERRORS = 3;
static constexpr int MAX_DEFICIENCY = 2;
#define THREADS 8
static constexpr int MIN_FACTOR_DEGREE = 1;
#define EXTENSION 1
static constexpr int FIELD_MODULUS = 9;

static constexpr int BITS = 5 * EXTENSION;
static constexpr uint64_t FIELD_SIZE = 1 << BITS;

#if EXTENSION == 1
typedef uint8_t field_storage_type;
#elif EXTENSION <= 3
typedef uint16_t field_storage_type;
#else
typedef uint32_t field_storage_type;
#endif

static inline uint32_t rdrand() {
    uint32_t ret;
#ifdef NDEBUG
    unsigned char ok;
    __asm__ volatile(".byte 0x0f, 0xc7, 0xf0; setc %1" : "=a"(ret), "=q"(ok) :: "cc");
    assert(ok);
#else
    FILE* f = fopen("/dev/urandom", "r");
    fread(&ret, sizeof(ret), 1, f);
    fclose(f);
#endif
    return ret;
}

struct MulTable {
    field_storage_type table[FIELD_SIZE][FIELD_SIZE];
    field_storage_type invtable[FIELD_SIZE];

    MulTable() {
        int64_t exptable[FIELD_SIZE], logtable[FIELD_SIZE];
        int64_t logx = 0, expx = 1;
        while (logx != FIELD_SIZE) {
            exptable[logx] = expx;
            logtable[expx] = logx;
            ++logx;
            expx *= 2;
            if (expx & FIELD_SIZE) expx ^= (FIELD_SIZE ^ FIELD_MODULUS);
        }

        for (uint64_t i = 0; i < FIELD_SIZE; ++i) {
            for (uint64_t j = 0; j < FIELD_SIZE; ++j) {
                uint64_t p = (i == 0 || j == 0) ? 0 : exptable[(logtable[i] + logtable[j]) % (FIELD_SIZE - 1)];
                table[i][j] = p;
                if (p == 1) {
                    invtable[i] = j;
                }
            }
        }

        for (uint64_t i = 1; i < FIELD_SIZE; ++i) {
            assert(invtable[invtable[i]] == i);
            assert(table[i][invtable[i]] == 1);
        }
    }
};

static const MulTable multable;

class FieldElem {
    field_storage_type val;
public:
    FieldElem() : val(0) {}
    explicit FieldElem(field_storage_type val_) : val(val_) {}

    bool operator==(FieldElem x) const { return val == x.val; }
    bool operator!=(FieldElem x) const { return val != x.val; }
    bool operator<(FieldElem x) const { return val < x.val; }
    bool operator>(FieldElem x) const { return val > x.val; }

    FieldElem operator*(FieldElem x) const { return FieldElem(multable.table[x.val][val]); }
    FieldElem& operator*=(FieldElem x) { val = multable.table[x.val][val]; return *this; }
    FieldElem Inverse() const { return FieldElem(multable.invtable[val]); }
    FieldElem operator+(FieldElem x) const { return FieldElem(x.val ^ val); }
    FieldElem operator+=(FieldElem x) { val ^= x.val; return *this; }
    uint32_t Int() const { return val; }
    bool IsZero() const { return val == 0; }
    bool IsOne() const { return val == 1; }
};

template<int N>
class Vector {
    FieldElem val[N];
public:
    Vector<N>& operator*=(FieldElem y) {
        auto tbl = multable.table[y.Int()];
        for (int i = 0; i < N; ++i) {
            val[i] = FieldElem(tbl[val[i].Int()]);
        }
        return *this;
    }

    Vector<N>& operator+=(Vector<N> x) {
        for (int i = 0; i < N; ++i) {
            val[i] += x[i];
        }
        return *this;
    }

    void AddMul(const Vector<N>& x, FieldElem y) {
        auto tbl = multable.table[y.Int()];
        for (int i = 0; i < N; ++i) {
            val[i] += FieldElem(tbl[x[i].Int()]);
        }
    }

    bool IsZero() const {
        for (int i = 0; i < N; ++i) {
            if (!val[i].IsZero()) return false;
        }
        return true;
    }

    bool IsOne() const {
        for (int i = 0; i < N; ++i) {
            if (!val[i].IsOne()) return false;
        }
        return true;
    }

    FieldElem& operator[](int pos) { return val[pos]; }
    const FieldElem& operator[](int pos) const { return val[pos]; }

    bool operator==(const Vector& a) const {
        for (int i = 0; i < N; ++i) {
            if (val[i] != a.val[i]) return false;
        }
        return true;
    }

    bool operator!=(const Vector& a) const {
        for (int i = 0; i < N; ++i) {
            if (val[i] == a.val[i]) return false;
        }
        return true;
    }

    bool operator<(const Vector& a) const {
        for (int i = 0; i < N; ++i) {
            if (val[i] < a.val[i]) return true;
            if (val[i] > a.val[i]) return false;
        }
        return false;
    }

    bool operator>(const Vector& a) const {
        for (int i = 0; i < N; ++i) {
            if (val[i] > a.val[i]) return true;
            if (val[i] < a.val[i]) return false;
        }
        return false;
    }

    uint32_t Hash() const {
        uint32_t ret = 0;
        if (N >= 1) ret = 22937 * val[N - 1].Int();
        if (N >= 2) ret += 17167 * val[N - 2].Int();
        if (N >= 3) ret += 28411 * val[N - 3].Int();
        return ret & 0x7FFF;
    }
};

FieldElem RandFieldElem() {
    return FieldElem(rdrand() % FIELD_SIZE);
}

template<int N>
Vector<N> RandVector() {
    Vector<N> ret;
    for (int i = 0; i < N; ++i) {
        ret[i] = RandFieldElem();
    }
    return ret;
}

template<int N>
Vector<N> PolyMulXMod(const Vector<N>& x, const Vector<N>& mod) {
    Vector<N> ret;
    for (int i = 0; i < N; ++i) {
        ret[i] = (i ? x[i - 1] : FieldElem()) + mod[i] * x[N - 1];
    }
    return ret;
}

template<int A, int B>
Vector<A> Low(const Vector<B>& x) {
    static_assert(A <= B, "ow");
    Vector<A> ret;
    for (int i = 0; i < A; ++i) {
        ret[i] = x[i];
    }
    return ret;
}

template<int A, int B>
Vector<A> High(const Vector<B>& x) {
    static_assert(A <= B, "ow");
    Vector<A> ret;
    for (int i = 0; i < A; ++i) {
        ret[i] = x[i + B - A];
    }
    return ret;
}

template<int A, int B>
bool Divides(const Vector<A>& a, const Vector<B>& b) {
    Vector<A> m;
    m[0] = FieldElem(1);
    for (int d = B - 1; d >= 0; --d) {
        m = PolyMulXMod(m, a);
        m[0] += b[d];
    }
    return m.IsZero();
}

template<int A, int B>
bool AnyFactor_(int pos, Vector<A>& a, const Vector<B>& b) {
    if (pos == A) {
        return Divides(a, b);
    }
    for (uint64_t i = 0; i < FIELD_SIZE; ++i) {
        a[pos] = FieldElem(i);
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
FieldElem Multiply(const Vector<A>& a, const Vector<A>& b) {
    FieldElem ret;
    for (int i = 0; i < A; ++i) {
        ret += a[i] * b[i];
    }
    return ret;
}

template<int A, int B>
Vector<A+B> Concat(const Vector<A>& a, const Vector<B>& b)
{
    Vector<A+B> ret;
    for (int i = 0; i < A; ++i) {
        ret[i] = a[i];
    }
    for (int i = 0; i < B; ++i) {
        ret[i + A] = b[i];
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
                row[r][c] = FieldElem(r == c);
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

    void MulRow(int r, FieldElem v) {
        row[r] *= v;
    }

    void AddMulRow(int r, int ro, FieldElem v) {
        row[r].AddMul(row[ro], v);
    }

    Matrix<R,C>& operator+=(const Matrix<R,C>& a) {
        for (int r = 0; r < R; ++r) {
            row[r] += a[r];
        }
        return *this;
    }

    Matrix<R,C>& operator*=(FieldElem v) {
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
            if (res[r][c].IsZero()) {
                int r2 = r + 1;
                while (r2 < R) {
                    if (!res[r2][c].IsZero()) {
                        res.SwapRows(r2, r);
                        SwapRows(r2, r);
                        break;
                    }
                    ++r2;
                }
                if (r2 < R) {
                    assert(!res[r][c].IsZero());
                } else {
                    assert(res[r][c].IsZero());
                    continue;
                }
            }
            FieldElem inv = res[r][c].Inverse();
            assert(!res[r][c].IsZero());
            res.MulRow(r, inv);
            MulRow(r, inv);
            assert(res[r][c].IsOne());
            for (int r2 = 0; r2 < R; ++r2) {
                if (r2 != r) {
                    FieldElem i = res[r2][c];
                    res.AddMulRow(r2, r, i);
                    AddMulRow(r2, r, i);
                    assert(res[r2][c].IsZero());
                }
            }
            ++r;
        }
        int ret = r;
        for (int c = 0; c < C; ++c) {
            int solvepos = -1;
            for (int r = 0; r < R; ++r) {
                if (!res[r][c].IsZero() && !res[r][c].IsOne()) {
                    solvepos = -2;
                    break;
                }
                if (solvepos == -1 && res[r][c].IsOne()) {
                    solvepos = r;
                } else if (solvepos != -1 && !res[r][c].IsZero()) {
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
                column[c][r] = mat[r][c];
            }
        }
    }

    Vector<R> Apply(const Vector<C>& in) const {
        Vector<R> out;
        for (int c = 0; c < C; ++c) {
            out.AddMul(column[c], in[c]);
        }
        return out;
    }
};

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
                ret.freedom[def_num][c] = residual[c][r] + FieldElem(r == c);
            }
            ++def_num;
        } else {
            for (int c = 0; c < N; ++c) {
                solutions[r][c] = inverse[r][c];
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
        FieldElem check = Multiply(partial.constraints[def_num], knowns);
        if (!check.IsZero()) return 0;
    }
    base_sol = partial.solutions.Apply(knowns);
    return ((uint64_t)1) << (BITS * partial.deficiency);
}

template<int N>
Vector<N> ExtSolution(const Vector<N>& base_sol, const PartialSolution<N>& partial, uint64_t solnum) {
    Vector<N> ext_sol = base_sol;
    for (int def_num = 0; def_num < partial.deficiency; ++def_num) {
        ext_sol.AddMul(partial.freedom[def_num], FieldElem((solnum >> (BITS * def_num)) % FIELD_SIZE));
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
            rr[i] = Low<ERRORS>(basis[pos[i]]);
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
        ret += Combination(i, LENGTH) * powl(FIELD_SIZE - 1, i - 1);
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
                if (ext_errors[i].IsZero() && !consec) {
                    ok = false;
                    break;
                }
            }
            if (!ok) continue;

            // Compute the full fault and verify it
            Vector<DEGREE-ERRORS> extbigfault;
            for (int i = 0; i < ERRORS; ++i) {
                extbigfault.AddMul(extbasis[ps.first[i]], ext_errors[i]);
            }
            Vector<DEGREE> bigfault = Concat(fault, extbigfault);

            if (allzerobefore) {
                for (int i = ERRORS; i < DEGREE; ++i) {
                    if (!bigfault[i].IsZero()) {
                        ok = (bigfault[i].IsOne());
                        break;
                    }
                }
                if (!ok) continue;
            }

            int num_error = 0;
            int min_pos = LENGTH;
            int max_pos = 0;
            for (int i = 0; i < ERRORS; ++i) {
                if (!ext_errors[i].IsZero()) {
                    ++num_error;
                    min_pos = std::min(min_pos, ps.first[i]);
                    max_pos = std::max(max_pos, ps.first[i]);
                }
            }

            if (num_error == 0) continue;
            res.emplace_back(bigfault, min_pos, max_pos, num_error);
            int nn = 0;
            for (int i = 0; i < ERRORS; ++i) {
                if (!ext_errors[i].IsZero()) {
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
#if THREADS > 1
        if ((hash % ((uint64_t)THREADS)) != (uint64_t)part) return true;
#endif
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
            sres[res[i].fault.Hash() % 1024].push_back(&res[i]);
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
    int max = allzerobefore ? 2 : FIELD_SIZE;
    FieldElem rrr = allzerobefore ? FieldElem() : RandFieldElem();
    for (int x = 0; x < max; ++x) {
        fault[pos] = FieldElem(x) + rrr;
        if (!RecurseShortFaults(pos + 1, allzerobefore && fault[pos].IsZero(), fault, psol, extbasis, part, hash * 9672876866715837617ULL + fault[pos].Int(), err)) {
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
            line += strprintf(" % 19.15f", (double)(results.count[e][l] / (Combination(e, l) * frac * powl(FIELD_SIZE - 1, e - 1)) * powl(FIELD_SIZE, DEGREE)));
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
    std::vector<std::thread> t;
#if THREADS > 1
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

std::string FormatCode(const Vector<DEGREE>& v) {
    std::string ret;
    ret.resize(DEGREE * EXTENSION);
    for (int d = 0; d < DEGREE; ++d) {
        uint64_t i = v[DEGREE - 1 - d].Int();
        for (int e = 0; e < EXTENSION; ++e) {
            ret[d * EXTENSION + e] = charset[(i >> ((EXTENSION - e - 1) * 5)) & 0x1F];
        }
    }
    return ret;
}

Vector<DEGREE> ParseCode(const std::string& str) {
    Vector<DEGREE> ret;
    if (str.size() != DEGREE * EXTENSION) {
        throw std::runtime_error("Wrong code length");
    }
    for (int i = 0; i < DEGREE; ++i) {
        uint64_t v = 0;
        for (int e = 0; e < EXTENSION; ++e) {
            v = v * 32;
            const char *ptr = strchr(charset, toupper(str[e + (DEGREE - 1 - i) * EXTENSION]));
            if (ptr == nullptr) {
                throw std::runtime_error("Unknown character");
            }
            v += ptr - charset;
        }
        ret[i] = FieldElem(v);
    }
    return ret;
}

int main(int argc, char** argv) {
    setbuf(stdout, NULL);
    Vector<DEGREE> gen;
    if (argc < 2) {
        fprintf(stderr, "Usage: %s GEN%i\n", argv[0], DEGREE);
        return 1;
    }
    if (strlen(argv[1]) == 0) {
        do {
            gen = RandVector<DEGREE>();
            code = FormatCode(gen);
            if (gen[0].IsZero()) { printf("%s: divisible by x\n", code.c_str()); continue; }
            if (MIN_FACTOR_DEGREE > 1 && AnyFactor<1, DEGREE>(gen)) { printf("%s: degree 1 factor\n", code.c_str()); continue; }
            if (MIN_FACTOR_DEGREE > 2 && AnyFactor<2, DEGREE>(gen)) { printf("%s: degree 2 factor\n", code.c_str()); continue; }
            if (MIN_FACTOR_DEGREE > 3 && AnyFactor<3, DEGREE>(gen)) { printf("%s: degree 3 factor\n", code.c_str()); continue; }
            if (MIN_FACTOR_DEGREE > 4 && AnyFactor<4, DEGREE>(gen)) { printf("%s: degree 4 factor\n", code.c_str()); continue; }
            if (MIN_FACTOR_DEGREE > 5 && AnyFactor<5, DEGREE>(gen)) { printf("%s: degree 5 factor\n", code.c_str()); continue; }
            if (MIN_FACTOR_DEGREE > 6 && AnyFactor<6, DEGREE>(gen)) { printf("%s: degree 6 factor\n", code.c_str()); continue; }
            break;
        } while(true);
    } else {
        gen = ParseCode(std::string(argv[1]));
    }
    code = FormatCode(gen);
    printf("%s: starting\n", code.c_str());

    if (argc >= 3) { require_err = strtoul(argv[2], NULL, 0); }
    if (argc >= 4) { require_len = strtoul(argv[3], NULL, 0); }

    basis_type basis;
    extbasis_type extbasis;
    basis.resize(LENGTH);
    extbasis.resize(LENGTH);
    Vector<DEGREE> x;
    x[0] = FieldElem(1);

    Matrix<DEGREE, DEGREE> rand;
    Transform<DEGREE, DEGREE> rand_trans;
    while(true) {
        Matrix<DEGREE, DEGREE> randi, res;
        for (int i = 0; i < DEGREE; ++i) {
            for (int j = 0; j < DEGREE; ++j) {
                rand[i][j] = RandFieldElem();
            }
        }
        randi = rand;
        int rank = randi.Invert(res);
        if (rank == DEGREE) { rand_trans = Transform<DEGREE, DEGREE>(rand); break; }
    }

    assert(LENGTH >= DEGREE);

    for (int i = 0; i < LENGTH; ++i) {
        Vector<DEGREE> base = Low<DEGREE>(x);
        basis[i] = rand_trans.Apply(base);
        extbasis[i] = High<DEGREE-ERRORS>(basis[i]);
        x = PolyMulXMod(x, gen);
    }

    ErrCount locs;
    locs.errors = 0;
    testalot(&basis, &extbasis, &locs);
    if (locs.errors == 0) {
        show_stats(locs);
    }

    return 0;
}
