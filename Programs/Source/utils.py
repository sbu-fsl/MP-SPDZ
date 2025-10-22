#!/usr/bin/env python3

import os, sys
# add MP-SPDZ dir to path so we can import from Compiler
sys.path.insert(0, os.path.dirname(sys.argv[0]) + '/../..') 
from Compiler.library import print_ln, if_e, else_
from Compiler.types import sint, cint, Array, sgf2n, cgf2n, regint, _number, _secret
from Compiler.compilerLib import Compiler # only used for testing
from Compiler.oram import OptimalORAM, AbstractORAM

# we assume these modules reside in Programs/Source/ 
from embeddings import apply_field_embedding, apply_inverse_field_embedding

def pad_byte[T: cgf2n | sgf2n](byte: T, offset: int) -> T:
    '''
    Given a byte = b_0,...,b_7 with offset number of meaningful bits, 
    we pad byte with a 1 followed by zeroes. For example, if offset = 4,
    then we want b_0, b_1, b_2, b_3, 1, 0, 0, 0.
    '''
    t = type(byte)
    byte = byte.bit_decompose(8)
    byte.reverse()
    byte[offset:] = [t(1)] + ([t(0)] * len(byte[offset+1:]))
    byte.reverse()
    return t.bit_compose(byte)

def str_to_hex(x):
        ''' Convert a string into a list of hex values. Obviously the string should represent valid hex to begin with. '''
        return [int(x[i : i + 2], 16) for i in range(0, len(x), 2)]

def get_random_sgf2n(bit_length: int, size=1) -> sgf2n:
    return sgf2n.bit_compose([sgf2n.get_random_bit(size=size) for _ in range(bit_length)])

def poly_eval[S,T: _number](coeffs: list[S], x: T) -> S|T:
    '''
    Use Horner's method to evaluate the polynomial defined by coeffs at the point x.
    Assumes coeffs[0] holds constant term.
    '''
    if len(coeffs) == 1:
        return coeffs[0]
    return coeffs[0] + x * poly_eval(coeffs[1:], x)

def interpolate_zero[T](xs: list[T], ys: list[T], size=1) -> T:
    '''
    Lagrange interpolate the point at x=0 from the points given by zip(xs,ys)
    '''
    assert(len(xs) == len(ys))
    deg = len(xs)
    t = type(xs[0])
    res = t(0, size=size)
    for i in range(deg):
        prod = t(1, size=size) 
        for j in range(deg):
            if j != i:
                prod *= xs[j].field_div((xs[j] - xs[i]))
        res += ys[i] * prod 
    return res


def find_nonzero_secret_idx(arr: AbstractORAM) -> _secret:
    '''Return secret index of last nonzero element in arr.'''
    t: _secret = arr.index_type
    num_entries = arr.size
    res = t(0)
    for i in range(num_entries):
        b = (arr[t(i)] != t(0))
        res = b.cond_swap(res, t(i))[0] # if b is True, res = t(i)
    return res

def dot_product_random_preimage(r: list[sgf2n], y: sgf2n) -> list[sgf2n]:
    '''
    Given a vector r and a scalar y, randomly sample a dot product preimage x such that <x,r>=y.
    NOTE: this implementation uses ORAM.
    NOTE: only support sgf2n for now.
    NOTE: this assumes values of r and y are GF(2^8) values embedded in GF(2^40)
    NOTE: returned list is embedded
    '''
    # init oram instances
    size = len(r)
    r_oram = OptimalORAM(size, sgf2n)
    x_oram = OptimalORAM(size, sgf2n)
    for i in range(size):
        r_oram[i] = r[i]
        x_oram[i] = apply_field_embedding(get_random_sgf2n(8))
    # solve for x
    j = find_nonzero_secret_idx(r_oram)
    x_j, r_j = x_oram[j], r_oram[j]
    x_oram[j] = (y - (sum(x_oram[i] * r_oram[i] for i in range(size)) - (x_j * r_j))).field_div(r_j)
    return [x_oram[i] for i in range(size)]

if __name__ == "__main__":
    usage = "usage: %prog [options] [args]"
    compiler = Compiler(usage=usage)

    @compiler.register_function("test_utils")
    def test_utils():
        print_ln("TEST find_nonzero_secret_idx")
        v = [sgf2n(1), sgf2n(0), sgf2n(2), sgf2n(0)]
        v_oram = OptimalORAM(4, sgf2n)
        for i in range(4):
            v_oram[i] = v
        idx = find_nonzero_secret_idx(v_oram).reveal()
        error_pattern = cgf2n(3) - idx
        @if_e(error_pattern != cgf2n(0))
        def _():
            print_ln("❌ TEST 1 FAILED\nidx=%s\nexpected idx=%s", idx, cgf2n(3))
        @else_
        def _():
            print_ln("✅ TEST 1 PASSED")


        print_ln("TEST dot_product_random_preimage")
        r = [sgf2n(1), sgf2n(0), sgf2n(2), sgf2n(0)]
        r = [apply_field_embedding(r_i) for r_i in r]
        y = sgf2n(5)
        y = apply_field_embedding(y)
        x = dot_product_random_preimage(r, y)
        y_res = sum(x_i * r_i for (x_i, r_i) in zip(x,r)).reveal()
        error_pattern = y.reveal() - y_res
        @if_e(error_pattern != cgf2n(0))
        def _():
            print_ln("❌ TEST 1 FAILED\ny=%s\nexpected y=%s", y_res, y.reveal())
        @else_
        def _():
            print_ln("✅ TEST 1 PASSED")



    compiler.compile_func()