#!/usr/bin/env python3

import os, sys, math
from itertools import product, chain
# add MP-SPDZ dir to path so we can import from Compiler
sys.path.insert(0, os.path.dirname(sys.argv[0]) + '/../..') 
from Compiler.library import print_ln, if_e, else_
from Compiler.types import sgf2n, cgf2n, get_sgf2nuint
from Compiler.compilerLib import Compiler # only used for testing

# we assume these modules reside in Programs/Source/
from lrss import lr_share
from shamir import shamir_share, shamir_reconstruct, obliv_shamir_reconstruct
from utils import get_random_sgf2n, poly_eval, mac, mac_verify

def robust_lr_share(
        msg: sgf2n,
        threshold: int,
        num_parties: int,
        mu: int, 
        secpar: int=40, 
        size: int=1,
) -> list[tuple]:
    '''
    Robust variant of lr_share obtained with pairwise MAC trick of Rabin and Ben-Or.
    '''
    n = num_parties
    lr_shares = lr_share(msg, threshold, n, mu, secpar, size)
    pairs = list(product(range(n),repeat=2))
    keys = {
        (j,i): (get_random_sgf2n(128, size=size), get_random_sgf2n(128, size=size)) 
        for j,i in pairs
    }
    tags = {}
    for i in range(n):
        src, ct, ss, mst = lr_shares[i]
        for j in range(n):
            tags[(j,i)] = mac(keys[(j,i)], [*src, ct, *ss, *mst])
    robust_shares = [
        (
         lr_shares[i], 
         [keys[(i,j)] for j in range(n)], 
         [tags[(j,i)] for j in range(n)]
        ) 
        for i in range(n)
    ]
    return robust_shares

def robust_lr_rec(
        shares: list[tuple],
        size: int=1
) -> sgf2n:
    '''
    Robust variant of lr_share obtained with pairwise MAC trick of Rabin and
    Ben-Or. **shares** must have length equal to the number of parties. 
    '''
    n = len(shares) 
    lr_shares = [s[0] for s in shares]
    pairs = list(product(range(n),repeat=2))
    keys = {(i,j): shares[i][1][j] for i,j in pairs}
    tags = {(i,j): shares[j][2][i] for i,j in pairs}
    candidates = [None for _ in range(n)]
    for i in range(n):
        '''
        We only want to reconstruct using lr_shares[i] if the i-th share is
        valid. Ideally, we would just filter lr_shares into a list containing
        only the valid shares, and pass this list (along with corresponding
        filtered coords) into lr_rec. Unfortunately, share validity is a runtime
        property determined by the above MAC checks. This means we cannot
        produce such filtered lists as the list size (i.e., number of valid
        shares) is only known at runtime. Note this issue precludes the use of
        runtime data structures (e.g., Array) since their size must also be a
        compile-time value. 

        Our solution is to introduce a bitmap valid_coords such that
        valid_coords[i] == sgf2n(1) if the i-th share is valid, else sgf2n(0).
        We define an "oblivious" version of Shamir reconstruction taking this
        bitmap as input, and proceed as in lr_rec, using oblivious
        reconstruction instead of the normal Shamir reconstruction.
        '''
        valid_coords = []
        for j in range(n):
            src, ct, ss, mst = lr_shares[j]
            b = mac_verify(
                    keys[(i,j)], 
                    [*src, ct, *ss, *mst], 
                    tags[(i,j)]
                )
            valid_coords.append(b)
        # begin adapted lr_rec code (too lazy to generalize lr_rec)
        [sources, ct, seed_shares, mask_shares_transposed] = list(map(tuple, zip(*lr_shares)))
        mask_shares = list(map(list, zip(*mask_shares_transposed)))
        masks = [obliv_shamir_reconstruct(s, valid_coords, size=size) for s in mask_shares]
        seed_shares_transposed = list(map(list, zip(*seed_shares)))
        seed = [obliv_shamir_reconstruct(s, valid_coords, size=size) for s in seed_shares_transposed]
        ext_outputs = [sum(seed[j] * source[j] for j in range(len(seed))) for source in sources]
        intermediate_shares = [ct[i] + ext_outputs[i] + masks[i] for i in range(len(ct))]
        candidates[i] = obliv_shamir_reconstruct(intermediate_shares, valid_coords, size=size)
        # end adapted lr_rec code
    '''a little hack to get most frequent item in candidates'''
    sgf2nuint_logn = get_sgf2nuint(math.ceil(math.log2(n)+1))
    frqs = [
        sum(sgf2nuint_logn(candidates[i].equal(candidates[j]), size=size) for j in range(n)) 
        for i in range(n)
    ] # frqs: list[sgf2nuint_logn]
    indicators: list[sgf2n] = [frq >= math.ceil(n/2) for frq in frqs] # honest majority assumption
    res = sgf2n(0, size=size)
    for c, b in zip(candidates, indicators):
        res = b.cond_swap(res, c)[0]
    return res



if __name__ == "__main__":
    usage = "usage: %prog [options] [args]"
    compiler = Compiler(usage=usage)

    @compiler.register_function("test_robust_lrss")
    def test_robust_lrss():
        print_ln("ROBUST LRSS TESTS")

        print_ln("-----TEST 1: BASIC-----")
        msg = sgf2n(2)
        shares = robust_lr_share(
            msg=msg,
            threshold=2,
            num_parties=3,
            mu=1,
            secpar=40,
        )
        rec_msg = robust_lr_rec(shares)
        error_pattern = (rec_msg - msg).reveal()
        @if_e(error_pattern != cgf2n(0))
        def _():
            print_ln("❌ TEST 1 FAILED\nreconstructed message=%s\nexpected message=%s", rec_msg.reveal(), msg.reveal())
        @else_
        def _():
            print_ln("✅ TEST 1 PASSED")

        print_ln("-----TEST 2: Vectorized-----")
        msg = sgf2n(list(range(100)))
        size = 100
        shares = robust_lr_share(
            msg=msg,
            threshold=2,
            num_parties=3,
            mu=1,
            secpar=40,
            size=size
        )
        rec_msg = robust_lr_rec(shares, size=size)
        error_pattern = (rec_msg - msg).reveal()
        @if_e(error_pattern != cgf2n(0))
        def _():
            print_ln("❌ TEST 2 FAILED\nreconstructed message=%s\nexpected message=%s", rec_msg.reveal(), msg.reveal())
        @else_
        def _():
            print_ln("✅ TEST 2 PASSED")
    
    compiler.compile_func()

