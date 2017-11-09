import random

def swapgen(gen):
  coefs = gen.list()
  degree = gen.degree()
  var = gen.parent().gen()
  return sum(coefs[degree-i] * var^i for i in range(degree+1)) / coefs[0]

def expand(gen,F,Q):
  if Q == 32:
      R = [1,2,4,8,16]
  elif Q == 1024:
      R = [1,2,4,8,16,32,64,128,256,512]
  ret = []
  var = gen.parent().gen()
  tgen = swapgen(gen)
  for m in range(1,Q):
      for pow in R:
          ret.append(gen.subs(F.fetch_int(m)*var).map_coefficients(lambda c: c^pow).monic())
          ret.append(tgen.subs(F.fetch_int(m)*var).map_coefficients(lambda c: c^pow).monic())
  print (len(set(ret)))
  return ret

CHARSET = "0123456789ABCDEFGHIJKLMNOPQRSTUV"

def base32repr(gen,Q):
    gens = ""
    genlist = gen.list()
    degree = gen.degree()
    pow = 0
    while Q > 1:
        pow = pow + 1
        Q >>= 5
    for p in range(degree):
        v = genlist[degree-1-p].integer_representation()
        ch = ""
        for i in range(pow):
            ch = CHARSET[int(v) % 32] + ch
            v /= 32
        gens += ch
    return gens

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
    fil = open("deg%i_len%i_hd%i" % (DEGREE, N, DISTANCE), 'w')
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

    MP_cache = dict()

    # Iterative over all base fields
    for F_modulus in F_moduli:
        F.<f> = GF(Q, repr='int', modulus=F_modulus)
        F_map = [extexp(f,i) for i in range(Q)]
        F_all = [x for x in F]
        F_from_int = {f.integer_representation() : f for f in F}
        FP.<x> = F[]

        # Find extension field modulus (any is fine)
        E_modulus_list = [0 for i in range(M)]
        while True:
            if M == 1:
                E = F
                e = f
                E_modulus = x
                break
            E_modulus = polyfromarray(x, [F(1)] + [F_from_int[j] for j in E_modulus_list])
            if E_modulus.is_primitive():
                E.<e> = F.extension(E_modulus)
                break
            U = 0
            while True:
                E_modulus_list[M-U-1] = (E_modulus_list[M-U-1] + 1) % Q
                if (E_modulus_list[M-U-1] != 0):
                    break
                U = U + 1

        # Find primitive element in extension field (any is fine)
        E_prim = e
        print(E_prim)
        E_base = e ** (((Q**M)-1) / N)

        # Find c values
        cs = []
        alpha = E_base
        mp = []
        num = DISTANCE - 1
        alphan = 1
        for i in range(1,N-num):
            alphan *= alpha
            if M == 1:
                minpol = x - alphan
            else:
                minpol = alphan.minpoly()
            mp.append(minpol)
            if (i >= num):
                generator=lcm(mp[-num:])
                if (generator.degree()  == DEGREE):
                    cs.append(i-num+1)

        all_generators = set()
        all_generators_exp = dict()

        # Iterate over all alphas
        for alphalog in range(1, N):
            if gcd(alphalog, N) != 1:
                continue
            # Iterate over all c values
            for c in cs:
                minpolyset=set()
                minpolys=[]
                for i in range(num):
                    root = (c + i) * alphalog % N
                    if root in MP_cache:
                        minpoly = MP_cache[root]
                    else:
                        if M == 1:
                            minpoly = x - (E_base ** root)
                        else:
                            minpoly = (E_base ** root).minpoly()
                        MP_cache[root] = minpoly
                    minpolyset.add(minpoly)
                    minpolys.append(minpoly)
                generator = 1
                for minpoly in minpolyset:
                    generator *= minpoly
                gens = base32repr(generator,Q)
                if gens in all_generators:
                    continue
                dupof = None
                if gens in all_generators_exp:
                    dupof = all_generators_exp[gens]
                all_generators.add(gens)
                if not dupof:
                    for exp in expand(generator, F, Q):
                        all_generators_exp[base32repr(exp,Q)] = gens
                if not dupof:
                    fil.write("GEN=%s F_mod=%r E_mod=%r alphalog=%r c=%r minpolys=%r gen=(%r)%s\n" % (gens, polyfromarray(B, [int(cc) for cc in reversed(F_modulus.coefficients(sparse=False))]), E_modulus.coefficients(sparse=False), alphalog, c, minpolys, generator, " DUP="+dupof if dupof else ""))
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
                gens = base32repr(generator,Q)
                print "      * GEN=%s N=%i M=%i F=(%r) E=(%r) alpha=(%r) powers=%i..%i minpolys=%s gen=(%s)" % (gens, N, M, F.modulus(), E.modulus(), alpha, i-num+1, i, mp[-num:], generator)
                find += 1
                if (find == max):
                    return find
            else:
                 pass
#                print "      * POLY of degree %i" % generator.degree()

attempt_exhaust(2,10,1,341,7,6)
exit

if True:
    Q=1024
    Ns={}
    for M in range(1,4):
      for d in (Q**M-1).divisors():
        if d > 30 and d < 2000 and d not in Ns:
          Ns[d] = M
    for N in sorted(Ns.keys()):
      M = Ns[N]
      attempt(Q,M,N,7,6,1)
else:
    for (E,L) in [(2,1023),(2,341),(4,1025),(4,205),(4,165),(3,1057),(3,151)]:
        for (DIST,DEG) in [(7,12),(6,12),(5,12)]:
            attempt_exhaust(2,5,E,L,DIST,DEG)
