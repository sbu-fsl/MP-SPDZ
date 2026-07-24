#!/usr/bin/env python3

import os, sys, math
from itertools import product, chain
# add MP-SPDZ dir to path so we can import from Compiler
sys.path.insert(0, os.path.dirname(sys.argv[0]) + '/../..') 
from Compiler.library import print_ln, if_e, else_
from Compiler.types import sint, cint, Array, sgf2n, cgf2n, regint, _number
from Compiler.compilerLib import Compiler # only used for testing
from Compiler.oram import OptimalORAM, AbstractORAM

# we assume these modules reside in Programs/Source/
from lrss import lr_share
from shamir import shamir_share, shamir_reconstruct
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
    # off_diagonal = {(j,i) in {0,...,n-1} x {0,...,n-1} : j != i}
    off_diagonal = filter(lambda pair : pair[0] != pair[1], product(range(n),repeat=2))
    keys = {(j,i): (get_random_sgf2n(128, size=size), get_random_sgf2n(128, size=size)) for j,i in off_diagonal}
    tags = {(j,i): mac(keys[(j,i)], list(chain.from_iterable(lr_shares[i]))) for j,i in off_diagonal}
    robust_shares = [
        (
         lr_share[i], 
         [keys[(i,j)] for j in range(n) if j != i], 
         [tags[(j,i)] for j in range(n) if j != i]
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

    **WARNING**: uses ORAM
    '''
    n = len(shares) 
    candidate_secrets = OptimalORAM(n, sgf2n)
    lr_shares = [s[0] for s in shares]
    # off_diagonal = {(j,i) in {0,...,n-1} x {0,...,n-1} : j != i}
    off_diagonal = filter(lambda pair : pair[0] != pair[1], product(range(n),repeat=2))
    keys = {(i,j): shares[i][1][j] for i,j in off_diagonal}
    tags = {(i,j): shares[j][2][i] for i,j in off_diagonal}
    for i in range(n):
        is_valid_share = {
            j: mac_verify(keys[(i,j)], lr_shares[j], tags[(i,j)]) 
            for j in range(n) if j != i
        } 
        is_valid_share[i] = sgf2n(1) # party i trusts his own share
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
        '''
        valid_coords = list(is_valid_share.values()) 




if __name__ == "__main__":
    usage = "usage: %prog [options] [args]"
    compiler = Compiler(usage=usage)

    @compiler.register_function("test_robust_lrss")
    def test_robust_lrss():
        print_ln("ROBUST LRSS TESTS")
    
    compiler.compile_func()

