#!/usr/bin/env python3

import os, sys
# add MP-SPDZ dir to path so we can import from Compiler
sys.path.insert(0, os.path.dirname(sys.argv[0]) + '/../..') 
from Compiler.library import print_ln, if_e, else_
from Compiler.types import sint, cint, Array, sgf2n, cgf2n, regint, _number
from Compiler.compilerLib import Compiler # only used for testing

# we assume these modules reside in Programs/Source/
from utils import get_random_sgf2n, poly_eval, interpolate_zero
from embeddings import apply_field_embedding, apply_inverse_field_embedding

def shamir_share[K,X](
    msg: K, 
    threshold: int, 
    num_parties: int, 
    eval_points:list[X]=None, 
    rand:list[K]=None, 
    size:int=1,
) -> tuple[list[X],list[K]]:
    '''
    Perform textbook Shamir's secret sharing. 
    
    :param eval_points: Defaults to 1,...,num_parties in clear type corresponding to msg type. 
    :param rand: Representing secret polynomial coefficients. Generated internally by default.
    '''
    assert threshold <= num_parties
    kt = type(msg)
    ct = kt if not hasattr(kt, "clear_type") else kt.clear_type

    # setup eval_points
    if eval_points is None:
        # by default, eval_points are 1,...,num_parties in clear type
        # corresponding to msg type.
        eval_points = [ct(i, size=size) for i in range(1, num_parties+1)] 
    
    # setup poly_coeffs
    poly_coeffs = []
    if rand:
        assert(len(rand) == threshold)
        poly_coeffs = rand
    else:
        if kt == sgf2n:
            poly_coeffs = [get_random_sgf2n(128, size=size) for _ in range(threshold)] 
        elif kt == sint:
            poly_coeffs = [sint.get_random(size=size) for _ in range(threshold)]
        elif kt == cint:
            poly_coeffs = [cint(regint.get_random(128, size=size)) for _ in range(threshold)]
        else:
            raise TypeError(f"type {kt} not yet supported")
    poly_coeffs[0] = msg
    
    # compute share values
    vals = [poly_eval(poly_coeffs, eval_points[i]) for i in range(num_parties)]
    return eval_points, vals

def shamir_reconstruct[K,X](
    vals: list[K],
    eval_points: list[X]=None,
    size=1,
) -> K:
    '''Shamir secret reconstruction.'''
    kt = type(vals[0])
    ct = kt if not hasattr(kt, "clear_type") else kt.clear_type
    # setup eval_points
    if eval_points is None:
        # by default, eval_points are 1,...,num_parties in clear type
        # corresponding to kt
        eval_points = [ct(i, size=size) for i in range(1, len(vals)+1)] 
    secret = interpolate_zero(eval_points, vals, size=size)
    return secret

if __name__ == "__main__":
    usage = "usage: %prog [options] [args]"
    compiler = Compiler(usage=usage)

    @compiler.register_function("test_shamir")
    def test_shamir():
        print_ln("SHAMIR TESTS")

        # Test 1: sgf2n, defaults
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

        # Test 2: sgf2n, vector
        msg = sgf2n([1,2,3,4,5])
        threshold = 2
        num_parties = 3
        size = 5
        _,y = shamir_share(msg, threshold=threshold, num_parties=num_parties, size=size)
        secret: sgf2n = shamir_reconstruct(y)
        error_pattern = (secret - msg).reveal()
        @if_e(error_pattern != cgf2n(0))
        def _():
            print_ln("❌ TEST 2 FAILED\nsecret=%s\nexpected secret=%s", secret.reveal(), msg.reveal())
        @else_
        def _():
            print_ln("✅ TEST 2 PASSED")

        # Test 3: sgf2n, embedding
        msg = sgf2n(64)
        msg_emb = apply_field_embedding(msg)
        threshold = 2
        num_parties = 3
        eval_points = [apply_field_embedding(cgf2n(i)) for i in range(1,num_parties+1)]
        rand = [apply_field_embedding(get_random_sgf2n(8)) for _ in range(threshold)]
        _,y = shamir_share(msg_emb, threshold=threshold, num_parties=num_parties, eval_points=eval_points, rand=rand)
        secret: sgf2n = apply_inverse_field_embedding(shamir_reconstruct(y))
        error_pattern = (secret - msg).reveal()
        @if_e(error_pattern != cgf2n(0))
        def _():
            print_ln("❌ TEST 3 FAILED\nsecret=%s\nexpected secret=%s", secret.reveal(), msg.reveal())
        @else_
        def _():
            print_ln("✅ TEST 3 PASSED")

        # Test 4: sgf2n, vector, embedding
        msg = sgf2n([1,2,3,4,5])
        msg_emb = apply_field_embedding(msg)
        size = 5
        threshold = 2
        num_parties = 3
        eval_points = [apply_field_embedding(cgf2n(i)) for i in range(1,num_parties+1)]
        rand = [apply_field_embedding(get_random_sgf2n(8, size=size)) for _ in range(threshold)]
        _,y = shamir_share(
            msg_emb, 
            threshold=threshold, 
            num_parties=num_parties, 
            eval_points=eval_points, 
            rand=rand, 
            size=size
        )
        secret: sgf2n = apply_inverse_field_embedding(shamir_reconstruct(y))
        error_pattern = (secret - msg).reveal()
        @if_e(error_pattern != cgf2n(0))
        def _():
            print_ln("❌ TEST 4 FAILED\nsecret=%s\nexpected secret=%s", secret.reveal(), msg.reveal())
        @else_
        def _():
            print_ln("✅ TEST 4 PASSED")


    


    compiler.compile_func()
