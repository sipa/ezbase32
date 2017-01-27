import random

def randlist(n, all):
    ret = []
    for i in range(n):
        ret.append(random.choice(all))
    return ret

def alllist(n, max):
    rets = []
    ret = [0 for x in range(n)]
    while True:
        rets.append(list(ret))
        p = 0;
        while p < n:
            if ret[p] + 1 < max:
                for q in range(p):
                    ret[q] = 0
                ret[p] += 1
                break
            p += 1
        if p == n:
            break
    return rets

def extexp(var, pow):
    if pow == 0:
        return 0
    return var**(pow - 1)

def polyfromarray(var, arr):
    ret = 0
    for a in arr:
        ret = ret * var + a
    return ret

used_F = {}

found = 0
used_E = {}

def attempt_exhaust(B,P,M,N,DISTANCE,DEGREE):
    Q = B**P
    format_len = int(math.ceil(log(B**(P*DEGREE)-1) /log(16)))
    assert(((Q**M)-1) % N == 0)
    assert(B.is_prime())

    # Find alphas
    alphalogs = []
    ZP = Integers(Q**M - 1)
    for i in range(Q**M - 1):
        if (ZP(i).additive_order() == N):
            alphalogs.append(i)

    # Find base field moduli
    BF.<bf> = GF(B)
    BP.<bp> = BF[]
    F_moduli = []
    for coefs in alllist(P, B):
        poly = polyfromarray(bp, [1] + coefs)
        if poly.is_primitive():
            F_moduli.append(poly)

    # Iterative over all base fields
    for F_modulus in F_moduli:
        F.<f> = GF(Q, repr='int', modulus=F_modulus)
        F_map = [extexp(f,i) for i in range(Q)]
        F_all = [x for x in F]
        F_from_int = {f.integer_representation() : f for f in F}
        FP.<x> = F[]

        # Find extension field modulus (any is fine)
        while True:
            E_modulus = polyfromarray(x, [F(1)] + randlist(M, F_all))
            E.<e> = F.extension(E_modulus)
            if E_modulus.is_primitive():
                break

        # Find primitive element in extension field (any is fine)
        while True:
            E_prim = polyfromarray(e, randlist(M, F_all))
            ok = True
            for di in (Q**M-1).divisors():
                if di != (Q**M-1) and (E_prim ** di) == 1:
                    ok = False
                    break
            if ok:
                break

        # Find c values
        cs = []
        alpha = E_prim ** alphalogs[0]
        mp = []
        num = DISTANCE - 1
        alphan = 1
        for i in range(1,N-num):
            alphan *= alpha
            mp.append(alphan.minpoly())
            if (i >= num):
                generator=lcm(mp[-num:])
                if (generator.degree()  == DEGREE):
                    cs.append(i-num+1)

        # Iterate over all alphas
        for alphalog in alphalogs:
            alpha = E_prim ** alphalog
            # Iterate over all c values
            for c in cs:
                c1 = alpha^c
                minpolys = [(c1*alpha^i).minpoly() for i in range(num)]
                generator = lcm(minpolys)
                table = []
                for p in range(P):
                    j = B**p
                    n = 0
                    for p in range(generator.degree()):
                        n = n * Q + (F_from_int[j] * generator.list()[generator.degree()-1-p]).integer_representation()
                    table.append(n)
                print "GEN {%s} F_mod=%r E_mod=%r E_primitive=%r alphalog=%r alpha=%r c=%r minpolys=%r" % (' '.join(['{0:#0{1}x}'.format(t, format_len + 2) for t in table]), polyfromarray(B, [int(cc) for cc in F_modulus.coefficients(sparse=False)]), E_modulus.coefficients(sparse=False), E_prim.list(), alphalog, alpha.list(), c, minpolys)
        break


def attempt(Q,M,N,DISTANCE,DEGREE,max):
    assert(((Q**M)-1) % N == 0)
    print "* ITERATION %i**%i / %i" % (Q, M, N)

    if Q in used_F:
        F = used_F[Q]
        f = F.gen()
    else:
        F.<f> = GF(Q, repr='int', modulus='random')
        used_F[Q] = F

    print "  * FIELD mod", F.modulus()

    F_all = [x for x in F]
    FP.<fp> = F[]
    F_from_int = {f.integer_representation() : f for f in F}

    if M == 1:
        E.<e> = F.extension(fp+1)
    else:
        if (Q, M) in used_E:
            print "* REUSING"
            E = used_E[(Q, M)]
            e = E.gen()
        else:
            while True:
                pol = polyfromarray(fp, [1] + randlist(M, F_all))
                if pol.is_primitive():
                    break
            E.<e> = F.extension(pol)
            used_E[(Q, M)] = E
    print "  * EXTFIELD mod", E.modulus()

    ok = False
    while not ok:
        alpha = polyfromarray(e, randlist(M, F_all)) ** ((Q**M-1) // N)
        alphan = 1
        n = 0
        if alpha^N != 1:
            continue
        for di in N.divisors():
            if di == N:
                ok = True
                break
            else:
                alphan *= (alpha ** (di - n))
                n = di
                if alphan == 1:
                    break

    print "  * ALPHA", alpha

    mp=[]
    ld={}
    print "  * LCM"
    num = DISTANCE-1 
    find = 0
    for i in range(1,N-num):
        mp.append((alpha^i).minpoly())
        if (i >= num):
            generator=lcm(mp[-num:])
            if (generator.degree() <= DEGREE):
                table=[]
                for j in range(Q):
                    n = 0
                    for p in range(generator.degree()):
                        n = n * Q + (F_from_int[j] * generator.list()[generator.degree()-1-p]).integer_representation()
                    table.append(n)
                print "      * TABLE: {%s}, // N=%i M=%i F=(%r) E=(%r) alpha=(%r) powers=%i..%i minpolys=%s gen=(%s)" % (' '.join(["0x%08x" % table[v] for v in [1,2,4,8,16]]), N, M, F.modulus(), E.modulus(), alpha, i-num+1, i, mp[-num:], generator)
                find += 1
                if (find == max):
                    return find
            else:
                 pass
#                print "      * POLY of degree %i" % generator.degree()

if False:
    Q=32
    Ns={}
    for M in range(1,10):
      for d in (Q**M-1).divisors():
        if d > 100 and d < 10000 and d not in Ns:
          Ns[d] = M
    for N in sorted(Ns.keys()):
      M = Ns[N]
      attempt(Q,M,N,7,11,1)
else:
    attempt_exhaust(2,5,2,341,7,11)
