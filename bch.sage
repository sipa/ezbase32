import random

def randlist(n, all):
    ret = []
    for i in range(n):
        ret.append(random.choice(all))
    return ret

def polyfromarray(var, arr):
    ret = 0
    for a in arr:
        ret = ret * var + a
    return ret

used_F = {}

found = 0
used_E = {}

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
            if (generator.degree() == DEGREE):
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

Q=32
Ns={}

for M in range(1,6):
  for d in (Q**M-1).divisors():
    if d > 64 and d < 65537 and d not in Ns:
      Ns[d] = M

for N in sorted(Ns.keys()):
  M = Ns[N]
  attempt(Q,M,N,4,6,1)
