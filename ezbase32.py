#!/usr/bin/env python3

import random

# GEN {0x0223776d 0x04464ffa 0x088c3efd 0x110cf8f3 0x209dd5cf} F_mod=37 E_mod=[29, 21, 1] E_primitive=[8, 30] alphalog=247 alpha=[6, 14] c=750 minpolys=[x^2 + 7*x + 4, x^2 + 3*x + 14, x^2 + 5*x + 21]

def bch(c, d):
    b = c >> 25
    c = d ^ (c & 0x1FFFFFF) << 5
    c ^= -((b >> 0) & 1) & 0x223776D
    c ^= -((b >> 1) & 1) & 0x4464FFA
    c ^= -((b >> 2) & 1) & 0x88C3EFD
    c ^= -((b >> 3) & 1) & 0x110CF8F3
    c ^= -((b >> 4) & 1) & 0x209DD5CF
    return c

def convertbits(data, frombits, tobits, pad=True):
    acc = 0
    bits = 0
    ret = []
    maxv = (1 << tobits) - 1
    for d in data:
        if d < 0 or (d >> frombits):
            return None
        acc = (acc << frombits) | d
        bits += frombits
        while (bits >= tobits):
            bits -= tobits
            ret.append((acc >> bits) & maxv)
    if (pad):
        if (bits):
            ret.append((acc << (tobits - bits)) & maxv)
    elif (acc << (tobits - bits)) & maxv:
        return None
    return ret

ZBASE32 = "bcdefghijkmnpqrstvwxyz0123456789"

INIT = 0x34056df6

def decode(s):
    r = []
    p = 0
    chk = 1
    for c in s:
        if c == '-':
            chk = bch(chk, p)
            r.append(-1)
        else:
            f = ZBASE32.find(c)
            if f == -1:
                return None
            chk = bch(chk, f)
            r.append(f)
    if c != INIT:
        return None
    
def baseencode(data):
    n = 0
    chk = [len(data)]
    for s in data:
        print("data: %r" % ','.join("%i" % x for x in s))
        n += len(s)
        chk.append(n % 32)
    print("meta: %r" % chk)
    for s in data:
        chk.extend(s)
    print("all: %r" % chk)
    checksum = bch(chk) ^ 0x0a8ce26f
    print("chk: %r" % chk)
    return ('-'.join([''.join([ZBASE32[x] for x in s]) for s in data])) + ''.join([ZBASE32[x] for x in convertbits([checksum], 30, 5)])

if __name__ == '__main__':
    random.seed()
    ff = [random.randrange(0,256) for x in xrange(32)]
    print("odata: %s" % ''.join("%02x" % x for x in ff))
    print("res: %r\n" % baseencode([[0,1],[0] + convertbits(ff, 8, 5)]))
