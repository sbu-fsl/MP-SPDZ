#!/usr/bin/env python3

import os, sys, math
# add MP-SPDZ dir to path so we can import from Compiler
sys.path.insert(0, os.path.dirname(sys.argv[0]) + '/../..') 
from Compiler.library import print_ln, if_e, else_
from Compiler.types import sint, cint, Array, sgf2n, cgf2n, regint, _number
from Compiler.compilerLib import Compiler # only used for testing

from new_shamir import shamir_share, shamir_reconstruct
from utils import dot_product_random_preimage, apply_field_embedding, apply_inverse_field_embedding, get_random_sgf2n

def get_source_length(
    leakage_budget: int, 
    leakage_secpar: int, 
    field_bit_length: int, 
    num_parties: int
) -> int:
    '''
    Helper function for determining (at compile-time) the source length for the
    inner product extractor used in the CKOS LRSS
    '''
    # see CKOS Theorem 1 for ext_error
    # see leftover hash lemma for min_entropy
    # for source_length, observe source needs (min_entropy + leakage_budget) total
    # min entropy before leakage. Add 1 because of how we sample source. Divide
    # by field_bit_length to get source_length in terms of embedded field elements
    ext_error = (2 ** (- leakage_secpar)) / (6 * num_parties)
    min_entropy = field_bit_length + 2 * math.log2(1 / ext_error) - 2
    source_length = math.ceil( (min_entropy + leakage_budget + 1) / field_bit_length)
    return source_length



def lr_share[T: sgf2n](
        msg: T,
        threshold: int,
        num_parties: int,
        leakage_budget: int, 
        leakage_secpar: int=40, 
        field_bit_length: int=8,
        size: int=1,
        embedded: bool=True
    ) -> tuple[list, list, list]:
    '''
    LRShare algorithm of CKOS22 leakage-resilient secret sharing scheme.
    https://eprint.iacr.org/2022/216
    We hardcode Shamir's secret sharing for MShare and SdShare, and instantiate
    the strong linear extractor with the inner product extractor.
    The consequence is that the resulting LRSS scheme enjoys perfect privacy and perfect local uniformity, 
    and we can easily set the leakage error with a single parameter leakage_secpar.
    
    :param msg: message to be shared. If msg is embedded, set embedded flag.
    :param threshold: secret sharing threshold
    :param num_parties: number of shares
    :param leakage_budget: number of bits of local leakage allowed on a share
    :param leakage_secpar: 2^{- leakage_secpar} is desired statistical distance
    between leakage distributions on share sets of any two messages (see paper).
    :param field_bit_length: number of bits needed to represent an element of
    the field in which the message lives (check embedding flag)
    :param size: the usual MP-SPDZ parallelization parameter
    :param embedded: is msg an embedded field element. For now we only support
    embedded elements so this is always set to true. 
    '''
    # handles for message type and corresponding clear type
    t = type(msg)
    ct = t if not hasattr(t, "clear_type") else t.clear_type

    # determine appropriate seed length (= source length)
    source_length = get_source_length(leakage_budget, leakage_secpar, field_bit_length, num_parties)
    seed_length = source_length

    # Shamir share msg to get intermediate shares
    eval_points = [apply_field_embedding(ct(i, size=size)) for i in range(1,num_parties+1)]
    rand = [apply_field_embedding(get_random_sgf2n(field_bit_length, size=size)) for _ in range(threshold)]
    _, intermediate_shares = shamir_share(msg, threshold, num_parties, rand=rand, size=size)
    
    # sample uniform extractor seed
    seed = [apply_field_embedding(get_random_sgf2n(field_bit_length, size=size)) for _ in range(seed_length)]
    # (2,n) secret share the extractor seed
    rands = [[apply_field_embedding(get_random_sgf2n(field_bit_length, size=size)) for _ in range(2)] for _ in range(seed_length)]
    seed_shares = [shamir_share(seed[i], 2, num_parties, eval_points=eval_points, rand=rands[i], size=size)[1] for i in range(seed_length)]
    seed_shares_by_party = {party: [] for party in range(num_parties)}
    for i in range(seed_length):
        seed_shares = shamir_share(
            msg=seed[i],
            threshold=2,
            num_parties=num_parties,
            eval_points=eval_points,
            rand=rands[i],
            size=size
        )[1]
        for party, sh in enumerate(seed_shares):
            seed_shares_by_party[party].append(sh)

    # sample sources. *** uses ORAM ***
    sources = [dot_product_random_preimage(seed, intermediate_shares[i]) for i in range(num_parties)]



