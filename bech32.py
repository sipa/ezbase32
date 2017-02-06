#!/usr/bin/python3

import hashlib

def bech32_polymod(values):
  GEN = [0x3b6a57b2, 0x26508e6d, 0x1ea119fa, 0x3d4233dd, 0x2a1462b3]
  chk = 1
  for v in values:
    b = (chk >> 25)
    chk = (chk & 0x1ffffff) << 5 ^ v
    for i in range(5):
      chk ^= GEN[i] if ((b >> i) & 1) else 0
  return chk

def bech32_hrp_expand(s):
  return [ord(x) >> 5 for x in s] + [0] + [ord(x) & 31 for x in s] + [0]

def bech32_verify_checksum(hrp, data):
  return bech32_polymod(bech32_hrp_expand(hrp) + data) == 1

def bech32_create_checksum(hrp, data):
  values = bech32_hrp_expand(hrp) + data
  polymod = bech32_polymod(values + [0,0,0,0,0,0]) ^ 1
  return [(polymod >> 5 * (5 - i)) & 31 for i in range(6)]

ZBASE32="ybndrfg8ejkmcpqxot1uwisza345h769"

def bech32_encode(hrp, data):
  combined = data + bech32_create_checksum(hrp, data)
  return hrp + '-' + ''.join([ZBASE32[d] for d in combined])

def bech32_decode(s):
  if any (ord(x) < 31 or ord(x) > 127 for x in s):
    return (None, None)
  pos = s.rfind('-')
  if pos < 1 or pos + 7 > len(s) or len(s) > 89:
    return (None, None)
  if not all(x in ZBASE32 for x in s[pos+1:]):
    return (None, None)
  hrp = s[:pos]
  data = [ZBASE32.find(x) for x in s[pos+1:]]
  if not bech32_verify_checksum(hrp, data):
    return (None, None)
  return (hrp, data[:-6])

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

def segwit_addr_encode(testnet, witver, witprog):
  hrp = "bctest" if testnet else "bc"
  assert (witver >= 0 and witver <= 16)
  return bech32_encode(hrp, [witver] + convertbits(witprog, 8, 5))

def segwit_addr_decode(testnet, addr):
  if any (ord(x) < 31 or ord(x) > 127 for x in addr):
    return (None, None)
  hrp, data = bech32_decode(addr.lower())
  hrpexp = "bctest" if testnet else "bc"
  if hrp != hrpexp:
    return None
  decoded = convertbits(data[1:], 5, 8, False)
  if decoded is None or len(decoded) < 2 or len(decoded) > 40:
    return None
  if (data[0] > 16):
    return None
  if (data[0] == 0 and len(decoded) != 20 and len(decoded) != 32):
    return None
  return (data[0], decoded)

KEY="0279BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798".decode("hex")

def hash160(d):
  sha256 = hashlib.new("sha256")
  ripemd160 = hashlib.new("ripemd160")
  sha256.update(d)
  ripemd160.update(sha256.digest())
  return ripemd160.digest()

def hash256(d):
  sha256 = hashlib.new("sha256")
  sha256.update(d)
  return sha256.digest()

KEYID = [ord(x) for x in hash160(KEY)]
SCRIPTID = [ord(x) for x in hash256("21".decode("hex") + KEY + "ac".decode("hex"))]

ENC1 = segwit_addr_encode(False, 0, KEYID)
ENC2 = segwit_addr_encode(True, 0, KEYID)
ENC3 = segwit_addr_encode(False, 0, SCRIPTID)
ENC4 = segwit_addr_encode(True, 0, SCRIPTID)
DEC1 = segwit_addr_decode(False, ENC1)
DEC2 = segwit_addr_decode(True, ENC2)
DEC3 = segwit_addr_decode(False, ENC3)
DEC4 = segwit_addr_decode(True, ENC4)

print("ENC1: %r" % ENC1)
print("ENC2: %r" % ENC2)
print("ENC3: %r" % ENC3)
print("ENC4: %r" % ENC4)
print("DEC1: %r" % (DEC1,))
print("DEC2: %r" % (DEC2,))
print("DEC3: %r" % (DEC3,))
print("DEC4: %r" % (DEC4,))

print(bech32_encode("bc", [1,2,3]));
