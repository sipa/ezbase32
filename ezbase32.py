#!/usr/bin/env python3

def rstable(x):
    t = ((x << 1) ^ x) << 1 ^ x
    h = x >> 4 | (t >> 7) << 10 | (t >> 9) << 20
    return ((x & 0xF) << 6 | ((t & 0x7F) | (t & 0x1FF) << 8) << 13) ^ h ^ (h << 3)

def rsupdate(crc, data):
    return ((crc & 0xFFFFF) << 10) ^ rstable((crc >> 20) ^ data)

def convert8to5(data):
    acc = 0
    bits = 0
    ret = []
    for byte in data:
        if byte < 0 or byte > 0xFF:
            raise Exception("Data bytes are outside [0..255]")
        acc = (acc << 8) | byte
        bits += 8
        while (bits >= 5):
            bits -= 5
            ret.append((acc >> bits) & 0x1F)
    if bits > 0:
        ret.append((acc << (5 - bits)) & 0x1F)
    return ret

def convert5to8(data):
    acc = 0
    bits = 0
    ret = []
    for u5 in data:
        if u5 < 0 or u5 > 31:
            raise Exception("Base32 data is outside [0..31]")
        acc = (acc << 5) | u5
        bits += 5
        if (bits >= 8):
            bits -= 8
            ret.append((acc >> bits) & 0xFF)
    if (acc << (8 - bits)) & 0xFF:
        return None
    return ret

ZBASE32 = "ybndrfg8ejkmcpqxot1uwisza345h769"

def tozbase32(vals):
    return ''.join([ZBASE32[x] for x in vals])

def fromzbase32(s):
    ret = []
    for c in s:
        p = ZBASE32.find(c)
        if p < 0:
            return None
        ret.append(p)
    return ret

def headerchecksum(header):
    crc = rsupdate(0, len(header) + 1)
    for c in header:
        crc = rsupdate(crc, c)
    return crc

def datachecksum(headercrc, data5):
    crc = rsupdate(headercrc, len(data5) + 1)
    pos = 0
    if (len(data5) & 1):
        crc = rsupdate(crc, data5[0])
        pos = 1
    while pos < len(data5):
        crc = rsupdate(crc, data5[pos] << 5 | data5[pos + 1])
        pos += 2
    return crc

def checksumto5(crc):
    return [(crc >> (5 * (5 - i))) & 0x1F for i in range(6)]

def checksumfrom5(data):
    return data[0] << 25 | data[1] << 20 | data[2] << 15 | data[3] << 10 | data[4] << 5 | data[5]

def encode(headerstr, databytes):
    if len(headerstr) > 100 or len(databytes) > 512:
        raise None
    crc = headerchecksum(headerstr.encode())
    u5 = convert8to5(databytes)
    if u5 is None:
        return None
    crc = datachecksum(crc, u5)
    return headerstr + tozbase32(u5 + checksumto5(crc))

def decode(ezbase32, headerlen):
    datalen = (len(ezbase32) - 6 - headerlen) * 5 // 8
    if headerlen < 0 or headerlen > 100 or datalen < 0 or datalen > 512 or (datalen * 8 + 4) // 5 + 6 + headerlen != len(ezbase32):
        return None
    crc = headerchecksum(ezbase32[0:headerlen].encode())
    u5 = fromzbase32(ezbase32[headerlen:])
    if u5 is None:
        return None
    data = convert5to8(u5[:-6])
    crc = datachecksum(crc, u5[:-6])
    crc2 = checksumfrom5(u5[-6:])
    if crc != crc2 or data is None:
        return None
    return ezbase32[0:headerlen], bytes(data)

if __name__ == '__main__':
    data = "hello"
    for x in range(100):
        print(tozbase32(convert8to5((data + str(x)).encode())))
        e = encode("btc:", (data + str(x)).encode())
        print(e)
        h,d = decode(e, 4)
        print(h,d.decode())

