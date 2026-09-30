#!/usr/bin/env python3

import os, sys
from operator import mul
# add MP-SPDZ dir to path so we can import from Compiler
sys.path.insert(0, os.path.dirname(sys.argv[0]) + '/../..') 
from Compiler.library import print_ln, if_e, else_
from Compiler.types import sint, cint, Array, sgf2n, cgf2n, regint, _number, _secret
from Compiler.compilerLib import Compiler # only used for testing
from Compiler.oram import OptimalORAM, AbstractORAM
from Compiler.util import tree_reduce

# we assume these modules reside in Programs/Source/ 
from embeddings import apply_field_embedding, apply_inverse_field_embedding

def int_to_sgf2n_bytes(x: int, num_bytes: int) -> list[sgf2n]:
    '''
    Convert compile-time int x into a list[sgf2n] of length num_bytes, where each sgf2n holds a byte of x in its lower 8 bits.
    '''
    x_bin_lst = [int(bit) for bit in list(bin(x)[2:][::-1])]
    return [sgf2n.bit_compose(x_bin_lst[j*8:(j+1)*8]) for j in range(num_bytes)]

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
    ''' 
    Convert a string into a list of hex values. 
    It is assumed the string represents valid hex. 
    '''
    return [int(x[i : i + 2], 16) for i in range(0, len(x), 2)]

def get_random_sgf2n_bits(bit_length: int, size: int = 1) -> sgf2n:
    '''
    Sample an sgf2n element where lower **bit_length**-many bits are uniformly
    random and all upper bits are zero.
    '''
    return sgf2n.bit_compose([sgf2n.get_random_bit(size=size) for _ in range(bit_length)])

def get_random_sgf2n(size: int = 1) -> sgf2n:
    '''
    Sample a uniform secret element of the active GF(2^n) field.
    '''
    # The first factor of an MP-SPDZ multiplication triple is uniform over the
    # compiled field. The other two factors are intentionally consumed and
    # discarded.
    return sgf2n.get_random_triple(size=size)[0]

def poly_eval[S,T: _number](coeffs: list[S], x: T) -> S|T:
    '''
    Use Horner's method to evaluate the polynomial defined by coeffs at the point x.
    Assumes coeffs[0] holds constant term.
    '''
    if len(coeffs) == 1:
        return coeffs[0]
    return coeffs[0] + x * poly_eval(coeffs[1:], x)

def poly_eval_balanced[S,T: _number](coeffs: list[S], x: T) -> S|T:
    '''
    Evaluate a non-empty polynomial with logarithmic multiplication depth.

    Coefficients are split into the even and odd powers recursively, using
    p(x) = p_even(x^2) + x * p_odd(x^2). ``poly_eval()`` intentionally remains
    the Horner evaluator for callers that do not need this larger circuit.
    '''
    assert len(coeffs) > 0
    powers = [x]
    for _ in range((len(coeffs) - 1).bit_length() - 1):
        powers.append(powers[-1] * powers[-1])

    def evaluate(values, depth):
        if len(values) == 1:
            return values[0]
        return (evaluate(values[::2], depth + 1)
                + powers[depth] * evaluate(values[1::2], depth + 1))

    return evaluate(coeffs, 0)

def interpolate_zero[X,Y: _number](xs: list[X], ys: list[Y], size=1) -> Y:
    '''
    Lagrange interpolate the point at x=0 from the points given by zip(xs,ys)
    '''
    assert(len(xs) == len(ys))
    deg = len(xs)
    t = type(ys[0])
    res = t(0, size=size)
    for i in range(deg):
        prod = t(1, size=size) 
        for j in range(deg):
            if j != i:
                prod *= xs[j].field_div((xs[j] - xs[i]))
        res += ys[i] * prod 
    return res

def interpolation_weights(
    xs: list,
    size: int=1
) -> list[sgf2n] | list[cgf2n]:
    '''
    Lagrange weights at zero. 
    '''
    n = len(xs)
    t = type(xs[0])
    weights = []
    for i in range(n):
        factors = []
        for j in range(n):
            if j != i:
                factors.append(xs[j].field_div(xs[j] - xs[i]))
        weights.append(tree_reduce(mul, factors))
    return weights



def obliv_interpolation_weights(
    xs: list,
    valid_xs: list,
    size: int=1
) -> list[sgf2n] | list[cgf2n]:
    '''
    Lagrange weights at zero, excluding invalid coordinates.

    Reuse these weights for different values with the same coordinates and
    validity bitmap. An invalid coordinate has weight zero.

    Inputs must be nonempty and equally sized. Coordinates must be distinct
    field elements, and validity entries must be public or secret field bits
    (zero or one). Public validity bits and coordinates produce public weights.
    '''
    assert(len(xs) == len(valid_xs))
    n = len(xs)
    t = type(valid_xs[0])
    one = t(1, size=size)
    weights = []
    for i in range(n):
        factors = []
        for j in range(n):
            if j != i:
                rational_term = xs[j].field_div(xs[j] - xs[i])
                factors.append(valid_xs[j].if_else(rational_term, one))
        prod = tree_reduce(mul, factors) if factors else one
        weights.append(prod * valid_xs[i])
    return weights


def obliv_interpolate_zero(
    xs: list,
    valid_xs: list,
    ys: list,
    size: int=1
) -> sgf2n:
    '''Interpolate at zero, excluding coordinates with validity bit zero.'''
    assert(len(xs) == len(valid_xs) and len(xs) == len(ys))
    weights = obliv_interpolation_weights(xs, valid_xs, size=size)
    return sum((y * weight for y, weight in zip(ys, weights)),
               type(ys[0])(0, size=size))


def find_nonzero_secret_idx(arr: AbstractORAM) -> _secret:
    '''Return secret index of last nonzero element in arr.'''
    t: _secret = arr.index_type
    num_entries = arr.size
    res = t(0)
    for i in range(num_entries):
        b = (arr[t(i)] != t(0))
        res = b.cond_swap(res, t(i))[0] # if b is True, res = t(i)
    return res

def dot_product_random_preimage(
        r: list[sgf2n], 
        y: sgf2n, 
) -> list[sgf2n]:
    '''
    Given a vector r and a scalar y, randomly sample a dot product preimage x
    such that <x,r>=y.
    
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
        x_oram[i] = apply_field_embedding(get_random_sgf2n_bits(8))
    # solve for x
    j = find_nonzero_secret_idx(r_oram)
    x_j, r_j = x_oram[j], r_oram[j]
    x_oram[j] = (y - (sum(x_oram[i] * r_oram[i] for i in range(size)) - (x_j * r_j))).field_div(r_j)
    return [x_oram[i] for i in range(size)]

def mac(key: tuple[sgf2n, sgf2n], msg: list[sgf2n]) -> sgf2n:
    '''
    One-time MAC from polynomial hashing with logarithmic multiplication depth.

    :param key: a tuple (a,b) 
    :param msg: a list of sgf2n
    '''
    a, b = key
    return poly_eval_balanced([a] + msg, b)

def mac_verify(key: tuple[sgf2n, sgf2n], msg:list[sgf2n], tag: sgf2n) -> sgf2n:
    '''
    Canonical MAC verification. Outputs sgf2n(1) if tag verifies, else sgf2n(0).
    '''
    return tag.equal(mac(key, msg))

if __name__ == "__main__":
    usage = "usage: %prog [options] [args]"
    compiler = Compiler(usage=usage)

    @compiler.register_function("test_utils")
    def test_utils():

        print_ln("TEST mac")
        key = (sgf2n(1), sgf2n(2))
        msg = [sgf2n(3), sgf2n(6)]
        tag = mac(key, msg)
        b = mac_verify(key, msg, tag).reveal()
        @if_e(b == cgf2n(0))
        def _():
            print_ln("❌ TEST mac FAILED")
        @else_
        def _():
            print_ln("✅ TEST mac PASSED")
        

        # print_ln("TEST find_nonzero_secret_idx")
        # v = [sgf2n(1), sgf2n(0), sgf2n(2), sgf2n(0)]
        # v_oram = OptimalORAM(4, sgf2n)
        # for i in range(4):
        #     v_oram[i] = v
        # idx = find_nonzero_secret_idx(v_oram).reveal()
        # error_pattern = cgf2n(3) - idx
        # @if_e(error_pattern != cgf2n(0))
        # def _():
        #     print_ln("❌ TEST 1 FAILED\nidx=%s\nexpected idx=%s", idx, cgf2n(3))
        # @else_
        # def _():
        #     print_ln("✅ TEST 1 PASSED")


        # print_ln("TEST dot_product_random_preimage")
        # r = [sgf2n(1), sgf2n(0), sgf2n(2), sgf2n(0)]
        # r = [apply_field_embedding(r_i) for r_i in r]
        # y = sgf2n(5)
        # y = apply_field_embedding(y)
        # x = dot_product_random_preimage(r, y)
        # y_res = sum(x_i * r_i for (x_i, r_i) in zip(x,r)).reveal()
        # error_pattern = y.reveal() - y_res
        # @if_e(error_pattern != cgf2n(0))
        # def _():
        #     print_ln("❌ TEST 1 FAILED\ny=%s\nexpected y=%s", y_res, y.reveal())
        # @else_
        # def _():
        #     print_ln("✅ TEST 1 PASSED")



    compiler.compile_func()
