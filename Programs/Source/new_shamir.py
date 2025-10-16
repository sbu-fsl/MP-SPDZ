#!/usr/bin/env python3

import os, sys
# add MP-SPDZ dir to path so we can import from Compiler
sys.path.insert(0, os.path.dirname(sys.argv[0]) + '/../..') 
from Compiler.library import print_ln, if_e, else_
from Compiler.types import sint, cint, Array, sgf2n, cgf2n, regint, _number
from Compiler.compilerLib import Compiler # only used for testing

from utils import get_random_sgf2n, poly_eval, interpolate_zero

def shamir_share[T: (S, C), S: (sint, sgf2n), C: (cint, cgf2n)](
    msg: T, 
    threshold: int, 
    num_parties: int, 
    eval_points:list[C]=None, 
    rand:list[T]=None, 
    size:int=1,
) -> tuple[list[C], list[T]]:
    '''Perform textbook Shamir's secret sharing.'''
    assert threshold <= num_parties
    t = type(msg)
    ct = t if not hasattr(t, "clear_type") else t.clear_type

    # setup eval_points
    if eval_points is None:
        # by default, eval_points are 1,...,num_parties interpreted as clear type of msg type
        eval_points = [ct(i, size=size) for i in range(1, num_parties+1)] 
    
    # setup poly_coeffs
    poly_coeffs = []
    if rand:
        assert(len(rand) == threshold)
        poly_coeffs = rand
    else:
        if t == sgf2n:
            # TODO: how can we reliably get field bit length at compile time? Seems difficult since we set field at runtime...
            poly_coeffs = [get_random_sgf2n(128, size=size) for _ in range(threshold)] 
        elif t == sint:
            poly_coeffs = [sint.get_random(size=size) for _ in range(threshold)]
        elif t == cint:
            poly_coeffs = [cint(regint.get_random(128, size=size)) for _ in range(threshold)]
        else:
            raise TypeError(f"type {t} not yet supported")
    poly_coeffs[0] = msg
    
    # compute share values
    vals = [poly_eval(poly_coeffs, eval_points[i]) for i in range(num_parties)]
    return eval_points, vals

def shamir_reconstruct[T: (S, C), S: (sint, sgf2n), C: (cint, cgf2n)](
    vals: T,
    eval_points: list[C]=None,
    size=1,
) -> T:
    '''Shamir secret reconstruction.'''
    t = type(vals[0])
    ct = t if not hasattr(t, "clear_type") else t.clear_type
    # setup eval_points
    if eval_points is None:
        # by default, eval_points are 1,...,num_parties interpreted as clear type of msg type
        eval_points = [ct(i, size=size) for i in range(1, len(vals)+1)] 
    secret = interpolate_zero(eval_points, vals, size=size)
    return secret


if __name__ == "__main__":
    usage = "usage: %prog [options] [args]"
    compiler = Compiler(usage=usage)

    @compiler.register_function("test_shamir")
    def test_shamir():
        print_ln("SHAMIR TESTS")

        msg = sgf2n(1)
        threshold = 2
        num_parties = 3
        _,y = shamir_share(msg, threshold=threshold, num_parties=num_parties)
        secret: sgf2n = shamir_reconstruct(y)
        error_pattern = (secret - msg).reveal()
        @if_e(error_pattern != cgf2n(0))
        def _():
            print_ln("❌ TEST 1 FAILED\nsecret=%s\nexpected secret=%s", secret.reveal(), msg.reveal())
        @else_
        def _():
            print_ln("✅ TEST 1 PASSED")

       


    compiler.compile_func()
