#!/usr/bin/env python3

def rstable(x):
    t = (x << 1 ^ x) << 1 ^ x
    h = x >> 4 | (t >> 7) << 10 | (t >> 9) << 20
    l = (x & 0xF) << 6 | (t & 0x7F | (t & 0x1FF) << 8) << 13
    return l ^ h ^ h << 3

RSTABLE = [rstable(x) for x in range(1024)]

def rsupdate(crc, data):
    return RSTABLE[crc >> 20 ^ data] ^ (crc & 0xFFFFF) << 10

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

def auxchecksum(aux):
    crc = rsupdate(0, len(aux) + 1)
    for c in aux:
        crc = rsupdate(crc, c)
    return crc

def encode(auxstr, databytes):
    if len(auxstr) > 100 or len(databytes) > 512:
        return None
    crc = auxchecksum(bytearray(auxstr, "utf8"))
    u5 = convertbits(databytes, 8, 5)
    for u10 in [1 + len(databytes)] + convertbits(u5, 5, 10):
        crc = rsupdate(crc, u10)
    return ''.join([ZBASE32[x] for x in u5 + convertbits([crc], 30, 5)])

def decode(auxstr, ezbase32):
    datalen = (len(ezbase32) - 6) * 5 // 8
    if len(auxstr) > 100 or datalen < 0 or datalen > 512 or (datalen * 8 + 34) // 5 != len(ezbase32):
        return None
    crc = auxchecksum(bytearray(auxstr, "utf8"))
    u5 = [ZBASE32.find(x) for x in ezbase32]
    databytes = convertbits(u5[:-6], 5, 8, False)
    if databytes is None:
        return None
    for u10 in [1 + len(databytes)] + convertbits(u5[:-6], 5, 10):
        crc = rsupdate(crc, u10)
    [crc2] = convertbits(u5[-6:], 5, 30)
    if crc != crc2 or databytes is None:
        return None
    return bytearray(databytes)

if __name__ == '__main__':
    ff = [100,200,300,400,500,600,700,800,900,1000]
    for x in range(10):
        crc = auxchecksum(ff)
        ff = ff + convertbits([crc], 30, 10)
    print("ff=%r\n" % ff)
    data = "hello"
    s = [0,20,55,55,55,55,55,55,55,55,55,55,55,55,55,55,55,55,55,55,55,55]
    print("%r\n" % encode("btc:", s))
    s = [0,32,55,55,55,55,55,55,55,55,55,55,55,55,55,55,55,55,55,55,55,55,55,55,55,55,55,55,55,55,55,55,55,55]
    print("%r\n" % encode("btc:", s))
    for x in range(100000):
        s = bytearray("hello" + str(x), "utf8")
        e = encode("btc:", s)
        d = decode("btc:", e)
        assert s == d
