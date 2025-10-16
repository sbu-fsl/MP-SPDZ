#!/usr/bin/env python3

import os, sys
# add MP-SPDZ dir to path so we can import from Compiler
sys.path.insert(0, os.path.dirname(sys.argv[0]) + '/../..') 
from Compiler.library import print_ln, if_e, else_
from Compiler.types import sint, cint, Array, sgf2n, cgf2n, regint, _number
from Compiler.compilerLib import Compiler # only used for testing

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

if __name__ == "__main__":
    usage = "usage: %prog [options] [args]"
    compiler = Compiler(usage=usage)

    @compiler.register_function("test_utils")
    def test_utils():
        pass

    compiler.compile_func()