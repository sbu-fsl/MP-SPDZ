#!/usr/bin/env python3

import os
import sys
from itertools import product
from operator import mul
# add MP-SPDZ dir to path so we can import from Compiler
sys.path.insert(0, os.path.dirname(sys.argv[0]) + '/../..') 
from Compiler.library import print_ln, if_e, else_
from Compiler.types import sgf2n, cgf2n, get_sgf2nuint
from Compiler.compilerLib import Compiler # only used for testing
from Compiler.util import tree_reduce

# we assume these modules reside in Programs/Source/
from lrss import lr_share
from utils import get_random_sgf2n, mac, mac_verify, obliv_interpolation_weights


def _parse_mac_message(share):
    source, ciphertext, seed_share, mask_share = share
    return [*source, ciphertext, *seed_share, *mask_share]

def _has_strict_majority(frequency, num_parties):
    return frequency >= num_parties // 2 + 1

def _frequency_bit_length(num_parties):
    # get_sgf2nuint() currently derives from MP-SPDZ's signed bit-integer type,
    # so the count needs an explicit sign bit.
    return num_parties.bit_length() + 1

def _frequency_sum(bits, size=1):
    """Count secret bits with a balanced tree of bounded integer adders.

    Generic sgf2n integer addition uses full-field bitwise AND for carries.
    Here all operands are known bits, so AND is simply field multiplication.
    Keep the sign bit required by MP-SPDZ's comparison implementation.
    """
    width = _frequency_bit_length(len(bits))

    def add(left, right):
        carry = 0
        result = []
        for a, b in zip(left, right):
            different = a + b
            result.append(different + carry)
            carry = a * b + different * carry
        return result

    # Explicit secret constants keep 1 + 1 in GF(2), not Python integers.
    words = [[sgf2n(bit, size=size)] + [0] * (width - 1) for bit in bits]
    return get_sgf2nuint(width).compose(tree_reduce(add, words))

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

    Party i receives (lr_shares[i], verification_keys, authentication_tags).
    Its keys authenticate each owner's share; its tags are indexed by verifier.
    '''
    n = num_parties
    lr_shares = lr_share(msg, threshold, n, mu, secpar, size)
    # keys[verifier, owner] authenticates the owner's share to the verifier.
    keys = {
        (verifier, owner): (
            get_random_sgf2n(size=size),
            get_random_sgf2n(size=size),
        )
        for verifier, owner in product(range(n), repeat=2)
    }
    tags = {}
    for owner, share in enumerate(lr_shares):
        message = _parse_mac_message(share)
        for verifier in range(n):
            tags[verifier, owner] = mac(keys[verifier, owner], message)
    return [
        (
            lr_shares[party],
            [keys[party, owner] for owner in range(n)],
            [tags[verifier, party] for verifier in range(n)],
        )
        for party in range(n)
    ]

def robust_lr_rec(
        shares: list[tuple],
        size: int=1
) -> sgf2n:
    '''
    Reconstruct from a nonempty list of robust shares in original party order.

    Include one share per party, using the layout returned by robust_lr_share.
    Each verifier produces a candidate using the shares its MAC checks accept.
    Return the strict-majority candidate, or zero if no strict majority exists.
    '''
    n = len(shares) 
    lr_shares = [s[0] for s in shares]
    pairs = list(product(range(n), repeat=2))
    keys = {(verifier, owner): shares[verifier][1][owner]
            for verifier, owner in pairs}
    tags = {(verifier, owner): shares[owner][2][verifier]
            for verifier, owner in pairs}
    coords = [cgf2n(j, size=size) for j in range(1, n + 1)]
    mac_messages = [_parse_mac_message(share) for share in lr_shares]

    # mac_verifies[owner][verifier] is a secret bit indicating verifier's acceptance 
    # of owner's share
    mac_verifies = [
        [
            mac_verify(keys[(verifier, owner)], mac_messages[owner], tags[(verifier, owner)]) 
            for verifier in range(n)
        ] 
        for owner in range(n)
    ]
    # valid_coords[i] == sgf2n(1) iff a majority of entries in mac_verifies[i] equal sgf2n(1)
    valid_coords = [
        _has_strict_majority(_frequency_sum(mac_verifies[i], size=size), n) for i in range(n)
    ]

    # reconstruct from lr_shares using only shares projected onto valid_coords
    weights = obliv_interpolation_weights(coords, valid_coords, size=size)
    sources, ct, seed_shares, mask_shares_transposed = zip(*lr_shares)
    dot = sgf2n.dot_product
    seed = [dot(column, weights) for column in zip(*seed_shares)]
    masks = [dot(column, weights) for column in zip(*mask_shares_transposed)]
    extractor_outputs = [dot(seed, src) for src in sources]
    intermediate_shares = [ct[i] + extractor_outputs[i] + masks[i] for i in range(len(ct))]
    msg = dot(intermediate_shares, weights)
    return msg


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
            threshold=1,
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
            threshold=1,
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

        print_ln("-----TEST 3: CORRUPTED SHARE AND TAGS-----")
        msg = sgf2n(42)
        shares = robust_lr_share(
            msg=msg,
            threshold=1,
            num_parties=3,
            mu=1,
            secpar=40,
        )
        lr_share_to_corrupt, keys, tags = shares[0]
        source, ciphertext, seed_shares, mask_shares = lr_share_to_corrupt
        shares[0] = (
            (source, ciphertext + sgf2n(1), seed_shares, mask_shares),
            keys,
            [tag + sgf2n(1) for tag in tags],
        )
        rec_msg = robust_lr_rec(shares)
        error_pattern = (rec_msg - msg).reveal()
        @if_e(error_pattern != cgf2n(0))
        def _():
            print_ln("❌ TEST 3 FAILED\nreconstructed message=%s\nexpected message=%s", rec_msg.reveal(), msg.reveal())
        @else_
        def _():
            print_ln("✅ TEST 3 PASSED")

        print_ln("-----TEST 4: EVEN PARTY COUNT, SIMD-----")
        size = 4
        msg = sgf2n([3, 5, 8, 13])
        shares = robust_lr_share(
            msg=msg,
            threshold=1,
            num_parties=4,
            mu=1,
            secpar=40,
            size=size,
        )
        rec_msg = robust_lr_rec(shares, size=size)
        error_pattern = (rec_msg - msg).reveal()
        @if_e(error_pattern != cgf2n(0))
        def _():
            print_ln("❌ TEST 4 FAILED\nreconstructed message=%s\nexpected message=%s", rec_msg.reveal(), msg.reveal())
        @else_
        def _():
            print_ln("✅ TEST 4 PASSED")
    
    compiler.compile_func()
