#!/usr/bin/env python3

# {, // N=1057 M=3 F=(x^5 + x^3 + 1) E=(e^3 + 11*e^2 + 11*e + 4) alpha=(11*e^2 + 11*e + 11) powers=340..342 minpolys=[x^3 + 24*x^2 + 22*x + 1, x^3 + 14*x^2 + 19*x + 1, x^3 + 14*x^2 + 19*x + 1] gen=(x^6 + 22*x^5 + 24*x^4 + 13*x^3 + 4*x^2 + 5*x + 1)

def bch(checksum, data):
    for d in data:
        b = (checksum >> 25) ^ d
        checksum = (checksum & 0x1FFFFFF) << 5
        checksum ^= -((b & 1) != 0) & 763793569 # 0x2d8690a1
        checksum ^= -((b & 2) != 0) & 194847042 # 0x0b9d2142
        checksum ^= -((b & 4) != 0) & 364823172 # 0x15bec284
        checksum ^= -((b & 8) != 0) & 704226344 # 0x29f9a428
        checksum ^= -((b & 16) != 0) & 58181712 # 0x0377c850
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

B5EXP=[1,2,4,8,16,9,18,13,26,29,19,15,30,21,3,6,12,24,25,27,31,23,7,14,28,17,11,22,5,10,20,1]
B5LOG=[-1,0,1,14,2,28,15,22,3,5,29,26,16,7,23,11,4,25,6,10,30,13,27,21,17,18,8,19,24,9,12,20]



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
