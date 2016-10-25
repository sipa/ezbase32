#!/usr/bin/env python3

# GEN {0x0223776d 0x04464ffa 0x088c3efd 0x110cf8f3 0x209dd5cf} F_mod=37 E_mod=[29, 21, 1] E_primitive=[8, 30] alphalog=247 alpha=[6, 14] c=750 minpolys=[x^2 + 7*x + 4, x^2 + 3*x + 14, x^2 + 5*x + 21]

def bch(checksum, data):
    for d in data:
        b = (checksum >> 25) ^ d
        checksum = (checksum & 0x1FFFFFF) << 5
        checksum ^= ((b & 1) != 0) * 35878765
        checksum ^= ((b & 2) != 0) * 71716858
        checksum ^= ((b & 4) != 0) * 143408893
        checksum ^= ((b & 8) != 0) * 286062835
        checksum ^= ((b & 16) != 0) * 547214799
    return checksum

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

ZBASE32 = "ybndrfg8ejkmcpqxot1uwisza345h769"

def encode(prefix, context, databytes):
    if len(prefix) > 30 or len(databytes) > 36:
        return None
    u5 = convertbits(databytes, 8, 5)
    checksum = bch(bch(bch(context, [1 + len(prefix)]), bytearray(prefix, "utf8")), u5)
    return prefix + ''.join([ZBASE32[x] for x in u5 + convertbits([checksum], 30, 5)])

def decode(context, ezbase32):
    prefixlen = None
    for p in range(31):
        if all(x in ZBASE32 for x in ezbase32[p:]):
            prefixlen = p
            break
    if prefixlen is None:
        return (None, None)
    datalen = (len(ezbase32) - prefixlen - 6) * 5 // 8
    if datalen < 0 or datalen > 36 or (datalen * 8 + 4) // 5 != len(ezbase32) - 6 - prefixlen:
        return (None, None)
    u5checksum = [ZBASE32.find(x) for x in ezbase32[prefixlen:]]
    checksum = bch(bch(bch(context, [1 + prefixlen]), bytearray(ezbase32[:prefixlen], "utf8")), u5checksum[:-6])
    [checksum2] = convertbits(u5checksum[-6:], 5, 30, False)
    if checksum != checksum2:
        return (None, None)
    return (ezbase32[:prefixlen], bytearray(convertbits(u5checksum[:-6], 5, 8, False)))

if __name__ == '__main__':
    ff = [100,200,300,400,500,600,700,800,900,1000]
    data = "hello"
    for x in range(100000):
        s = bytearray("hello" + str(x), "utf8")
        print("s: %s" % s)
        e = encode("btc:", 0, s)
        print("e: %r" % e)
        (prefix, d) = decode(0, e)
        assert prefix == "btc:"
        assert d == s
