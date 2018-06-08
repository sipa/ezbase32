#!/usr/bin/python3

# * ITERATION 32**2 / 1023
#   * FIELD mod x^5 + x^2 + 1
#   * EXTFIELD mod e^2 + 26*e + 5
#   * ALPHA 2*e + 5
#   * LCM
#       * GEN=I5DFMIVBNVU4LISMBRK0QQA5F5S N=1023 M=2 F=(x^5 + x^2 + 1) E=(e^2 + 26*e + 5) alpha=(2*e + 5) powers=20..33 minpolys=[x^2 + 25*x + 25, x^2 + 21*x + 8, x^2 + 10*x + 14, x^2 + 26*x + 24, x^2 + 25*x + 18, x^2 + 19*x + 13, x^2 + 3*x + 19, x^2 + x + 23, x^2 + 26*x + 16, x^2 + 28*x + 28, x^2 + 12*x + 21, x^2 + 7*x + 1, x^2 + 17*x + 26, x + 26] gen=(x^27 + 18*x^26 + 5*x^25 + 13*x^24 + 15*x^23 + 22*x^22 + 18*x^21 + 31*x^20 + 11*x^19 + 23*x^18 + 31*x^17 + 30*x^16 + 4*x^15 + 21*x^14 + 18*x^13 + 28*x^12 + 22*x^11 + 11*x^10 + 27*x^9 + 20*x^8 + 26*x^6 + 26*x^5 + 10*x^4 + 5*x^3 + 15*x^2 + 5*x + 28)

FIELD=5
EXTFIELD=(26, 5)
CONST="I5DFMIVBNVU4LISMBRK0QQA5F5S"
ALPHA=2*32 + 5
GEN=(2, 5)
C=20
SYNDROMES=14

VALS=[ord(c) - ord('0') if c >= '0' and c <= '9' else ord(c) - ord('A') + 10 for c in CONST]

def red(x):
    while x >= 32:
        pos = 0
        y = x
        while y > 1:
            pos += 1
            y >>= 1
        x ^= (FIELD + 32) << (pos - 5)
    return x

def table():
    lt=[0] * 32
    et=[0] * 31
    e = 1
    lt[0] = -1
    for l in range(31):
        lt[e] = l
        et[l] = e
        e = red(e << 1)
    return (et, lt)

TABLE = table()

def mul(a,b):
    if (a == 0 or b == 0):
        return 0
    return TABLE[0][(TABLE[1][a] + TABLE[1][b]) % 31]

def exttable():
    lt = [0] * 1024
    et = [0] * 1023
    e = 1
    lt[0] = -1
    for l in range(1023):
        lt[e] = l
        et[l] = e
        e1 = e >> 5
        e0 = e & 31
        n2 = mul(e1, GEN[0])
        n1 = mul(e0, GEN[0]) ^ mul(e1, GEN[1]) ^ mul(n2, EXTFIELD[0])
        n0 = mul(e0, GEN[1]) ^ mul(n2, EXTFIELD[1])
        e = (n1 << 5) | n0
    return (et, lt)

EXTTABLE = exttable()

print(EXTTABLE)

def extmul(a,b):
    if (a == 0 or b == 0):
        return 0
    return EXTTABLE[0][(EXTTABLE[1][a] + EXTTABLE[1][b]) % 1023]

def extinv(a):
    return EXTTABLE[0][(1023 - EXTTABLE[1][a]) % 1023]

def toint(a, bits=5):
    r = 0
    for v in a:
        r = (r << bits) | v
    return r

GENS=[toint([red(x << i) for x in VALS]) for i in range(5)]

def polymod(values):
    chk = 1
    for value in values:
        top = chk >> 130
        chk = ((chk - (top << 130)) << 5) ^ value
        for i in range(5):
            chk ^= GENS[i] if ((top >> i) & 1) else 0
    return chk

def print_polymod():
    s = ""
    s += "def polymod(values):\n"
    s += "    generator = [0x%x, 0x%x, 0x%x, 0x%x, 0x%x]\n" % (GENS[0], GENS[1], GENS[2], GENS[3], GENS[4])
    s += "    chk = 1\n"
    s += "    for value in values:\n"
    s += "        top = chk >> 130\n"
    s += "        chk = (chk - top << 130) << 5 ^ value\n"
    s += "        for i in range(5):\n"
    s += "            chk ^= generator[i] if ((top >> i) & 1) else 0\n"
    s += "return chk\n"
    return s

def syndrome(values, a):
    chk = 1
    for v in values:
        chk = extmul(chk, a) ^ v
    return chk

def syndromes(values):
    a = 1
    for i in range(C):
        a = extmul(a, ALPHA)
    s = []
    for i in range(SYNDROMES):
        s.append(syndrome(values, a))
        a = extmul(a, ALPHA)
    return s

def print_syndromes():
    s = ""
    s += "def syndrome(p):\n"
    s += "    s = 0\n"
    for b in range(135):
        p = 1 << b
        values = [(p >> (26 - i) * 5) & 31 for i in range(27)]
        ss = syndromes(values)
        s += "    s ^= 0x%x if ((p >> %i) & 1) else 0\n" % (toint(ss, 10), b)
    s += "    return s\n"
    return s

def matrix(values, errors):
    s = syndromes(values)
    m = []
    for i in range(errors):
        v = []
        for j in range(errors + 1):
            v.append(s[i + j])
        m.append(v)
    return m

def extreduce(m):
    for rc in range(len(m)):
        found = None
        for r in range(rc, len(m)):
            if m[r][rc] != 0:
                found = r
                break
        if found is None:
            return None
        t = m[found]
        m[found] = m[rc]
        m[rc] = t
        i = extinv(m[rc][rc])
        for c in range(len(m[0])):
            m[rc][c] = extmul(i, m[rc][c])
        for r in range(len(m)):
            if r != rc:
                p = m[r][rc]
                for c in range(len(m[0])):
                    m[r][c] ^= extmul(p, m[rc][c])
    return m

def locate_errors(values, length):
    for v in range(SYNDROMES // 2, -1, -1):
        m = extreduce(matrix(values, v))
        if m is None:
            continue
        print(m)
        a = 1
        pos = []
        for p in range(length):
            d = 0
            for j in range(v):
                d = extmul(d, a) ^ m[j][v]
            d = extmul(d, a) ^ 1
            if d == 0:
                pos.append(p)
            a = extmul(a, extinv(ALPHA))
        return pos

CHARSET = "qpzry9x8gf2tvdw0s3jn54khce6mua7l"

DATA = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27]
CHKSUM = polymod(DATA + [0]*27) ^ 1 ^ 32 ^ 1024 ^ 32678 ^ 1048576 ^ 33554432 ^ 1073741824
XDATA = DATA + [(CHKSUM >> 5 * (26 - i)) & 31 for i in range(27)]

print(''.join(CHARSET[x] for x in XDATA))

print(XDATA)
print(polymod(XDATA))
print(locate_errors(XDATA, 1023))

print(print_polymod())
print(print_syndromes())
