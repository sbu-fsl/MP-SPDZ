#!/usr/bin/env python3

import os, sys
from copy import copy
# add MP-SPDZ dir to path so we can import from Compiler. It is assumed this file lives in MP-SPDZ/Programs/Source. 
sys.path.insert(0, os.path.dirname(sys.argv[0]) + '/../..') 
from Compiler.library import print_ln, vectorize, if_e, else_
from Compiler.types import cgf2n, sgf2n, Array, Matrix, VectorArray
from Compiler.compilerLib import Compiler

def apply_field_embedding_bd(in_bytes: list[cgf2n | sgf2n]) -> list[cgf2n | sgf2n]:
    '''
    Applies the field embedding f: GF(2^8) -> GF(2^40) given by x = y^5+1.
    Assumes irreducible polynomials x^8 + x^4 + x^3 + x + 1 (per FIPS 197) for GF(2^8) and 
    y^40 + y^20 + y^15 + y^10 + 1 for GF(2^40) (this is the irreducible polynomial chosen when an 
    MPC protocol is executed with --lg2 40)

    :param in_bytes: list[cgf2n | sgf2n]. Assumed to hold an element of GF(2^8) in bit decomposed form. 
    :returns: list[cgf2n | sgf2n]. Same type as in_bytes. Image of in_bytes in GF(2^40) under embedding f.
    '''
    out_bytes = [type(in_bytes[0])(0) for _ in range(8)] # will hold coefficients of 1, y^5, y^10,...,y^35 determined by embedding
    # embedding f can be computed as as:
    # f( \sum_{i=0}^7 a_i x^i ) 
    # = \sum_{i=0}^7 a_i (y^5+1)^i 
    # = \sum_{i=0}^7 a_i ( \sum_{k=0}^i \binom{i}{k} y^5k)
    out_bytes[0] = sum(in_bytes[0:8]) # a_0 + ... + a_7
    out_bytes[1] = in_bytes[1] + in_bytes[3] + in_bytes[5] + in_bytes[7] # a_1 + a_3 + a_5 + a_7
    out_bytes[2] = in_bytes[2] + in_bytes[3] + in_bytes[6] + in_bytes[7] # a_2 + a_3 + a_6 + a_7
    out_bytes[3] = in_bytes[3] + in_bytes[7] # a_3 + a_7
    out_bytes[4] = in_bytes[4] + in_bytes[5] + in_bytes[6] + in_bytes[7] # a_4 + a_5 + a_6 + a_7
    out_bytes[5] = in_bytes[5] + in_bytes[7] # a_5 + a_7
    out_bytes[6] = in_bytes[6] + in_bytes[7] # a_6 + a_7
    out_bytes[7] = in_bytes[7] # a_7
    return out_bytes

def apply_inverse_field_embedding_bd(in_bytes: list[cgf2n | sgf2n]) -> list[cgf2n | sgf2n]:
    '''
    Apply the left inverse f^{-1} of the field embedding f: GF(2^8) -> GF(2^40) given by x = y^5 + 1.
    Assumes irreducible polynomials x^8 + x^4 + x^3 + x + 1 (per FIPS 197) for GF(2^8) and 
    y^40 + y^20 + y^15 + y^10 + 1 for GF(2^40) (this is the irreducible polynomial chosen when an 
    MPC protocol is executed with --lg2 40)

    :param in_bytes: list[cgf2n | sgf2n]. Assumed to hold the image of some GF(2^8) element under the embedding f in bit decomposed form. 
    :returns: list[cgf2n | sgf2n]. Same type as in_bytes. Holds f^{-1}(y) 
    '''
    out_bytes = [type(in_bytes[0])(0) for _ in range(8)] # will hold coefficients of x^0, x^1,...,x^7 determined by inverse embedding
    # essentially undo the XORs from apply_field_embedding_bd
    out_bytes[7] = in_bytes[7]
    out_bytes[6] = in_bytes[6] + out_bytes[7]
    out_bytes[5] = in_bytes[5] + out_bytes[7]
    out_bytes[4] = in_bytes[4] + out_bytes[5] + out_bytes[6] + out_bytes[7]
    out_bytes[3] = in_bytes[3] + out_bytes[7]
    out_bytes[2] = in_bytes[2] + out_bytes[3] + out_bytes[6] + out_bytes[7]
    out_bytes[1] = in_bytes[1] +  out_bytes[3] + out_bytes[5] + out_bytes[7]
    out_bytes[0] = in_bytes[0] + sum(out_bytes[1:8])
    return out_bytes

def apply_field_embedding(x: cgf2n | sgf2n) -> cgf2n | sgf2n:
    '''
    Applies the field embedding f: GF(2^8) -> GF(2^40) given by x = y^5+1.
    Assumes irreducible polynomials x^8 + x^4 + x^3 + x + 1 (per FIPS 197) for GF(2^8) and 
    y^40 + y^20 + y^15 + y^10 + 1 for GF(2^40) (this is the irreducible polynomial chosen when an 
    MPC protocol is executed with --lg2 40)

    :param x: cgf2n | sgf2n. Assumed to hold an element of GF(2^8) in its lower 8 bits. Assumed cgf2n is 40 bits long
    :returns: cgf2n | sgf2n. Same type as x. Image of x in GF(2^40) under embedding f.
    '''
    in_bytes = x.bit_decompose(8) # select lower 8 bits of x into list[cgf2n] length 8. LSB first. 
    out_bytes = apply_field_embedding_bd(in_bytes)

    # now that we have the coefficients in out_bytes, need to multiply them by their respective y^{5k} and sum into a single cgf2n/sgf2n
    return type(x)(sum(out_bytes[idx] * (cgf2n(2) ** (5*idx)) for idx in range(8)))

def apply_inverse_field_embedding(y: cgf2n | sgf2n) -> cgf2n | sgf2n:
    '''
    Apply the left inverse f^{-1} of the field embedding f: GF(2^8) -> GF(2^40) given by x = y^5 + 1.
    Assumes irreducible polynomials x^8 + x^4 + x^3 + x + 1 (per FIPS 197) for GF(2^8) and 
    y^40 + y^20 + y^15 + y^10 + 1 for GF(2^40) (this is the irreducible polynomial chosen when an 
    MPC protocol is executed with --lg2 40)

    :param y: cgf2n | sgf2n. Assumed to hold the image of some GF(2^8) element under the embedding f. 
    :returns: cgf2n | sgf2n. Same type as y. Holds f^{-1}(y) in its lower 8 bits. 
    '''
    # select bits of y corresponding to coefficients of y^0, y^5, y^10,...,y^35. LSB first
    in_bytes = y.bit_decompose(bit_length=40, step=5)
    out_bytes = apply_inverse_field_embedding_bd(in_bytes)
   
    # now that we have the coefficients in out_bytes, need to multiply them by their respective x^k
    return type(y)(sum(out_bytes[idx] * (cgf2n(2) ** idx) for idx in range(8)))

class EmbeddedInverter():
    '''
    Used to invert GF(2^8) elements embedded in GF(2^40) via x = y^5+1 embedding. 
    Recall that the multiplicative group of GF(2^n) is cyclic, so for any 
    z in GF(2^8), z^254 = z^-1. We can compute z^254 much faster than naive approach 
    with exponentiation by squaring. 

    (NOTE: Is this really faster / worth doing vs. just using built in inversion method on GF(2^40) element? Probably?)

    Specifically, z^254 = z^2 * z^4 * z^8 * z^16 * z^32 * z^64 * z^128.
    Knowing this, it would be handy if we could precompute these powers for any z. 
    It suffices to precompute y^{2^i * j} in GF(2^40) for i=0 to 7, j=0 to 39 (aka the embedded_powers table).

    This is because we can leverage the following fact about arithmetic in GF(2^n):
        For any a_0,...,a_m in GF(2^n), (a_0 + ... + a_m)^{2^i} = (a_0^{2^i} + ... + a_m^{2^i})

    To be extra concrete, if we have a GF(2^40) element z = c_0 + c_1y + ... + c_{39}y^{39},
    then to compute z^{2^i}, we can take every non-zero c_j*y^j term in z and raise it to 2^i: c_j * y^{2^i * j},
    which just amounts to a lookup in the embedded_powers table. 
    '''

    def __init__(self, size=1):
        '''
        Only need to instantiate this class once in order to set up the precomputation table. 

        :param size: number of elements we want to invert at once. 
        '''
        # Here is what the rows of _embedded_powers correspond to if they were written as polynomials instead of bytes:
        # [y^0, y^1, y^2, y^3, ..., y^39] // no modular reduction needed because we never cross y^40
        # [y^0, y^2, y^4, y^6,..., y^38, y^40 = y^20+y^15+y^10+1 = 0x108401,...],
        # [y^0, y^4,...]
        # [y^0, y^8,...]
        # [y^0, y^16,...]
        # [y^0, y^32,...]
        # [y^0, y^64 = y^39+y^34+y^19+y^14+y^4 = 0x8400084010,...]
        # [y^0, y^128 = 0x108,...]
        _embedded_powers = [
            [0x1,0x2,0x4,0x8,0x10,0x20,0x40,0x80,0x100,0x200,0x400,0x800,0x1000,0x2000,0x4000,0x8000,0x10000,0x20000,0x40000,0x80000,0x100000,0x200000,0x400000,0x800000,0x1000000,0x2000000,0x4000000,0x8000000,0x10000000,0x20000000,0x40000000,0x80000000,0x100000000,0x200000000,0x400000000,0x800000000,0x1000000000,0x2000000000,0x4000000000,0x8000000000],
            [0x1,0x4,0x10,0x40,0x100,0x400,0x1000,0x4000,0x10000,0x40000,0x100000,0x400000,0x1000000,0x4000000,0x10000000,0x40000000,0x100000000,0x400000000,0x1000000000,0x4000000000,0x108401,0x421004,0x1084010,0x4210040,0x10840100,0x42100400,0x108401000,0x421004000,0x1084010000,0x4210040000,0x840008401,0x2100021004,0x8400084010,0x1000000842,0x4000002108,0x100021,0x400084,0x1000210,0x4000840,0x10002100],
            [0x1,0x10,0x100,0x1000,0x10000,0x100000,0x1000000,0x10000000,0x100000000,0x1000000000,0x108401,0x1084010,0x10840100,0x108401000,0x1084010000,0x840008401,0x8400084010,0x4000002108,0x400084,0x4000840,0x40008400,0x400084000,0x4000840000,0x8021004,0x80210040,0x802100400,0x8021004000,0x210802008,0x2108020080,0x1080010002,0x800008421,0x8000084210,0x108,0x1080,0x10800,0x108000,0x1080000,0x10800000,0x108000000,0x1080000000],
            [0x1,0x100,0x10000,0x1000000,0x100000000,0x108401,0x10840100,0x1084010000,0x8400084010,0x400084,0x40008400,0x4000840000,0x80210040,0x8021004000,0x2108020080,0x800008421,0x108,0x10800,0x1080000,0x108000000,0x800108401,0x10002108,0x1000210800,0x20004010,0x2000401000,0x42008020,0x4200802000,0x84200842,0x8420084200,0x2000421084,0x40000420,0x4000042000,0x10040,0x1004000,0x100400000,0x40108401,0x4010840100,0x1080200040,0x8021080010,0x2100421080],
            [0x1,0x10000,0x100000000,0x10840100,0x8400084010,0x40008400,0x80210040,0x2108020080,0x108,0x1080000,0x800108401,0x1000210800,0x2000401000,0x4200802000,0x8420084200,0x40000420,0x10040,0x100400000,0x4010840100,0x8021080010,0x40108421,0x1080000040,0x100421080,0x4200040100,0x1084200,0x842108401,0x1004210042,0x2008400004,0x4210000008,0x401080210,0x840108001,0x1000000840,0x100001000,0x840100,0x8401000000,0x800000001,0x84210800,0x2100001084,0x210802100,0x8001004210],
            [0x1,0x100000000,0x8400084010,0x80210040,0x108,0x800108401,0x2000401000,0x8420084200,0x10040,0x4010840100,0x40108421,0x100421080,0x1084200,0x1004210042,0x4210000008,0x840108001,0x100001000,0x8401000000,0x84210800,0x210802100,0x800000401,0x2100420080,0x8000004000,0x4010002,0x4000800100,0x842000420,0x8421084,0x421080210,0x80010042,0x10802108,0x800000020,0x1084,0x8401084010,0x1004200040,0x4000840108,0x100020,0x2108401000,0x8400080210,0x84210802,0x10802100],
            [0x1,0x8400084010,0x108,0x2000401000,0x10040,0x40108421,0x1084200,0x4210000008,0x100001000,0x84210800,0x800000401,0x8000004000,0x4000800100,0x8421084,0x80010042,0x800000020,0x8401084010,0x4000840108,0x2108401000,0x84210802,0x20,0x8000004210,0x2100,0x8401004,0x200800,0x802108420,0x21084000,0x4200842108,0x2000020000,0x1084210000,0x100421,0x1004010,0x10840008,0x108421080,0x1000200840,0x108001,0x8020004210,0x10040108,0x2108401004,0x1084210040],
            [0x1,0x108,0x10040,0x1084200,0x100001000,0x800000401,0x4000800100,0x80010042,0x8401084010,0x2108401000,0x20,0x2100,0x200800,0x21084000,0x2000020000,0x100421,0x10840008,0x1000200840,0x8020004210,0x2108401004,0x400,0x42000,0x4010000,0x421080000,0x21004,0x2008420,0x210800100,0x4200002,0x401000210,0x2108401084,0x8000,0x840000,0x80200000,0x8421000000,0x420080,0x40108400,0x4210002000,0x84000040,0x8020004200,0x2108400084]
        ]

        embedded_powers = VectorArray(8 * 40, cgf2n, size)
        for i,_list in enumerate(_embedded_powers):
            for j,x in enumerate(_list):
                embedded_powers[40 * i + j] = cgf2n(x, size=size)
        self.embedded_powers = embedded_powers
        self.size = size
    
    def repeated_squaring(self, bd_val: list[cgf2n | sgf2n], exponent: int) -> cgf2n | sgf2n:
        '''
        Compute bd_val^{2^exponent} using lookups into self.embedded_powers

        :param bd_val: list[cgf2n | sgf2n]. A bit-decomposed GF(2^40) value. 
        :param exponent: int. Constrained to 0 <= exponent <= 7, by self.embedded_powers lookup table. 
        '''
        return sum(self.embedded_powers[exponent * 40 + idx] * bd_val[idx] for idx in range(len(bd_val)))

    def invert(self, val: cgf2n | sgf2n) -> cgf2n | sgf2n:
        '''
        Compute val^254 via exponentiation by squaring. 

        :param val: cgf2n/sgf2n, assumed to be a GF(2^8) value embedded in GF(2^40)
        :returns: cgf2n/sgf2n, same type as val
        '''

        bd_val = val.bit_decompose(bit_length=40)
        powers_of_val = [0] * 129
        for idx in range(1,8):
            powers_of_val[2 ** idx] = self.repeated_squaring(bd_val, idx)
        return (powers_of_val[2] * powers_of_val[4] * powers_of_val[8] * powers_of_val[16] * powers_of_val[32] * powers_of_val[64] * powers_of_val[128])
    
if __name__ == "__main__":
    usage = "usage: %prog [options] [args]"
    compiler = Compiler(usage=usage)
    compiler.parse_args()

    @compiler.register_function("test_embedding")
    def test_embedding():
        # test case 1: (x^3+x^2+x+1 = 0xf) <=> (0x8000 = y^15)
        a = cgf2n(0xf)
        b = apply_field_embedding(a)
        b_preimage = apply_inverse_field_embedding(b)
        error_pattern = a + b_preimage
        @if_e(error_pattern)
        def _():
            print_ln("❌ TEST 1 FAILED:\nb_preimage=%s\nexpected=%s", b_preimage, a)
        @else_
        def _():
            print_ln("✅ TEST 1 PASSED")

        # test case 2: (x^6+x = 0x42) <=> (0x40100420 = y^30+y^20+y^10+y^5) 
        a = cgf2n(0x42)
        b = apply_field_embedding(a)
        b_preimage = apply_inverse_field_embedding(b)
        error_pattern = a + b_preimage
        @if_e(error_pattern)
        def _():
            print_ln("❌ TEST 2 FAILED:\nb_preimage=%s\nexpected=%s", b_preimage, a)
        @else_
        def _():
            print_ln("✅ TEST 2 PASSED")

        # test case 3: 0xf * 0xc7 = 0x1 in GF(2^8), so multiplication of their embeddings should also equal 0x1 in GF(2^40)
        a = cgf2n(0xf)
        b = cgf2n(0xc7)
        a_embedded = apply_field_embedding(a)
        b_embedded = apply_field_embedding(b)
        c_embedded = a_embedded * b_embedded
        error_pattern = c_embedded + cgf2n(1)
        @if_e(error_pattern)
        def _():
            print_ln("❌ TEST 3 FAILED:\nc_embedded=%s\nexpected=%s", c_embedded, cgf2n(1))
        @else_
        def _():
            print_ln("✅ TEST 3 PASSED")

    compiler.compile_func()