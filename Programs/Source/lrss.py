#!/usr/bin/env python3

import os, sys, math
# add MP-SPDZ dir to path so we can import from Compiler
sys.path.insert(0, os.path.dirname(sys.argv[0]) + '/../..') 
from Compiler.library import print_ln, if_e, else_
from Compiler.types import sint, cint, Array, sgf2n, cgf2n, regint, _number
from Compiler.compilerLib import Compiler # only used for testing

from shamir import shamir_share, shamir_reconstruct
from utils import dot_product_random_preimage, apply_field_embedding, apply_inverse_field_embedding, get_random_sgf2n

def get_source_length(
    num_parties: int,
    mu: int,
    secpar: int, 
    k: int=128, 
) -> int:
    '''
    Compile-time helper function for determining the source length for the
    inner product extractor used in the SV LRSS scheme. See 3.2 in LRPSS paper
    for derivation.
    '''
    eps_ext = (2 ** (- secpar)) / (2 * num_parties)
    lower_bound = k + mu + 2 * math.log2(1 / eps_ext) - 2
    source_length = math.ceil(lower_bound / k) # eta / k
    return source_length



def lr_share(
        msg: sgf2n,
        threshold: int,
        num_parties: int,
        mu: int, 
        secpar: int=40, 
        field_bit_length: int=128,
        size: int=1,
    ) -> list[list]:
    '''
    LRShare algorithm of Srinivasan and Vasudevan strong local leakage-resilient
    secret sharing scheme. https://eprint.iacr.org/2018/1154.
    
    We hardcode Shamir's secret sharing with default eval points
    1,...,num_parties for the underlying threshold scheme, and instantiate the
    extractor with the inner product extractor. The source/seed length is
    determined by a desired statistical security parameter **secpar**, and a
    desired leakage budget **mu**.  
    
    :param msg: message to be shared.
    :param threshold: secret sharing threshold
    :param num_parties: number of shareholders
    :param mu: number of bits of local leakage allowed on a share
    :param secpar: 2**(- secpar) is desired statistical distance
    between views on share sets of any two messages (see paper). 
    :param field_bit_length: bit-length of field msg lives in. Only 128 for now. 
    :param size: the usual MP-SPDZ parallelization argument 
    '''
    source_length = get_source_length(num_parties, mu, secpar)
    seed_length = source_length

    # Shamir share msg to get intermediate shares
    eval_points = [cgf2n(i, size=size) for i in range(1,num_parties+1)]
    intermediate_shares = shamir_share(msg, threshold, num_parties, size=size)[1]

    # uniformly sample one extractor seed, num_parties sources, and num_parties masks
    seed = [get_random_sgf2n(field_bit_length, size=size) for _ in range(seed_length)]
    sources = [[get_random_sgf2n(field_bit_length, size=size) for _ in range(source_length)] for _ in range(num_parties)]
    masks = [get_random_sgf2n(field_bit_length, size=size) for _ in range(num_parties)]

    # double mask intermediate shares
    # crucially assumes characteristic-two field where addition is XOR
    ext_outputs = [sum(seed[j] * source[j] for j in range(seed_length)) for source in sources]
    ct = [intermediate_shares[i] + ext_outputs[i] + masks[i] for i in range(num_parties)]
    
    # share masks and seed.
    mask_shares = [shamir_share(masks[i], threshold, num_parties, size=size)[1] for i in range(num_parties)]
    mask_shares_transposed = list(map(list, zip(*mask_shares)))
    seed_shares_transposed = [shamir_share(seed[j], threshold, num_parties, size=size)[1] for j in range(seed_length)]
    seed_shares = list(map(list, zip(*seed_shares_transposed)))

    # output final shares
    shares = [[sources[i], ct[i], seed_shares[i], mask_shares_transposed[i]] for i in range(num_parties)]
    return shares

def lr_rec(
        shares: list[list],
        coords: list[cgf2n]=None,
        size: int=1
    ):
    '''
    LRRec algorithm of Srinivasan and Vasudevan strong local leakage-resilient
    secret sharing scheme. https://eprint.iacr.org/2018/1154. 

    :param shares: each list in shares represents an LRSS share of the form
    output by lr_share. 
    :param coords: coords[i] indicates which party shares[i] comes from. If
    None, assumes coords = [cgf2n(i) for i in range(1, len(shares)+1)]
    :param size: the usual parallelization parameter
    '''
    if not coords: coords = [cgf2n(i) for i in range(1, len(shares)+1)]
    [sources, ct, seed_shares, mask_shares_transposed] = list(map(list, zip(*shares)))

    # reconstruct masks and seed
    mask_shares = list(map(list, zip(*mask_shares_transposed)))
    masks = [shamir_reconstruct(s, coords, size=size) for s in mask_shares]
    seed_shares_transposed = list(map(list, zip(*seed_shares)))
    seed = [shamir_reconstruct(s, coords, size=size) for s in seed_shares_transposed]

    # unmask intermediate shares
    ext_outputs = [sum(seed[j] * source[j] for j in range(len(seed))) for source in sources]
    intermediate_shares = [ct[i] + ext_outputs[i] + masks[i] for i in range(len(ct))]
    msg = shamir_reconstruct(intermediate_shares, coords, size=size)
    return msg
    

if __name__ == "__main__":
    usage = "usage: %prog [options] [args]"
    compiler = Compiler(usage=usage)

    @compiler.register_function("test_lrss")
    def test_lrss():
        print_ln("LRSS TESTS")

        print_ln("-----TEST 1: Basic-----")
        msg = sgf2n(2)
        shares = lr_share(
            msg=msg,
            threshold=2,
            num_parties=3,
            mu=1,
            secpar=40,
        )
        # for i, share in enumerate(shares):
        #     [source, ct, seed_shares, mask_shares_transposed] = share
        #     source = [s.reveal() for s in source]
        #     ct = ct.reveal()
        #     seed_shares = [s.reveal() for s in seed_shares]
        #     mask_shares_transposed = [s.reveal() for s in mask_shares_transposed]
        #     print_ln("shares[%s] = %s, %s, %s, %s\n", i, source, ct, seed_shares, mask_shares_transposed)
        rec_msg = lr_rec(shares)
        error_pattern = (rec_msg - msg).reveal()
        @if_e(error_pattern != cgf2n(0))
        def _():
            print_ln("❌ TEST 1 FAILED\nreconstructed message=%s\nexpected message=%s", rec_msg.reveal(), msg.reveal())
        @else_
        def _():
            print_ln("✅ TEST 1 PASSED")

        print_ln("-----TEST 2: vectorized-----")
        msg = sgf2n(list(range(100)))
        size = 100
        shares = lr_share(
            msg=msg,
            threshold=2,
            num_parties=3,
            mu=1,
            secpar=40,
            size=100
        )
        rec_msg = lr_rec(shares)
        error_pattern = (rec_msg - msg).reveal()
        @if_e(error_pattern != cgf2n(0))
        def _():
            print_ln("❌ TEST 2 FAILED\nreconstructed message=%s\nexpected message=%s", rec_msg.reveal(), msg.reveal())
        @else_
        def _():
            print_ln("✅ TEST 2 PASSED")
    
    compiler.compile_func()

