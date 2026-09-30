#!/usr/bin/env python3

import os
import sys
from itertools import product
# add MP-SPDZ dir to path so we can import from Compiler
sys.path.insert(0, os.path.dirname(sys.argv[0]) + '/../..') 
from Compiler.library import print_ln, if_e, else_, start_timer, stop_timer
from Compiler.types import sgf2n, cgf2n, MemValue, MPCThread
from Compiler.compilerLib import Compiler # only used for testing

# we assume these modules reside in Programs/Source/
from lrss import lr_share
from utils import get_random_sgf2n, mac, obliv_interpolation_weights


def _parse_mac_message(share):
    source, ciphertext, seed_share, mask_share = share
    return [*source, ciphertext, *seed_share, *mask_share]

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
    Reconstruct using the owners' shares accepted by a strict majority of verifiers.

    Each owner's worker reveals masked MAC differences and computes a public
    strict-majority validity bit.
    '''

    n = len(shares) 
    lr_shares = [s[0] for s in shares]
    pairs = list(product(range(n), repeat=2))
    # MemValue shares scalar or SIMD values between MP-SPDZ threads.
    keys = {(verifier, owner): tuple(MemValue(k) for k in shares[verifier][1][owner])
            for verifier, owner in pairs}
    tags = {(verifier, owner): MemValue(shares[owner][2][verifier])
            for verifier, owner in pairs}
    coords = [cgf2n(j, size=size) for j in range(1, n + 1)]
    mac_messages = [[MemValue(x) for x in _parse_mac_message(share)] for share in lr_shares]

    # Public masked MAC differences: zero means the verifier accepts the owner.
    start_timer(1)
    mac_verifies = [[MemValue(cgf2n(0, size=size)) for owner in range(n)]
                    for verifier in range(n)]
    valid_coords = [MemValue(cgf2n(0, size=size)) for owner in range(n)]

    def is_valid_coord(owner):
        message = [x.read() for x in mac_messages[owner]]
        for verifier in range(n):
            key = tuple(k.read() for k in keys[verifier, owner])
            tag = tags[verifier, owner].read()
            r = get_random_sgf2n(size=size)
            masked_difference = (mac(key, message) - tag) * r
            mac_verifies[verifier][owner].write(masked_difference.reveal())
        column = [mac_verifies[verifier][owner].read() for verifier in range(n)]
        # cgf2n equality checks the whole field and returns a public regint bit.
        frequency = sum(entry == 0 for entry in column)
        valid_coords[owner].write(cgf2n(frequency >= n // 2 + 1, size=size))

    # Each tape writes its owner's column and strict-majority validity bit.
    threads = [MPCThread(is_valid_coord, 'is_valid_coord', args=[owner], single_thread=True)
               for owner in range(n)]
    for thread in threads:
        thread.start()
    for thread in threads:
        thread.join()
    valid_coords = [bit.read() for bit in valid_coords]
    stop_timer(1)

    # reconstruct from lr_shares using only shares projected onto valid_coords
    weights = obliv_interpolation_weights(coords, valid_coords, size=size)
    sources, ct, seed_shares, mask_shares_transposed = zip(*lr_shares)
    mixed_dot = lambda x, y : sum(xi * yi for xi, yi in zip(x,y))
    sdot = sgf2n.dot_product
    seed = [mixed_dot(column, weights) for column in zip(*seed_shares)]
    masks = [mixed_dot(column, weights) for column in zip(*mask_shares_transposed)]
    extractor_outputs = [sdot(seed, src) for src in sources]
    intermediate_shares = [ct[i] + extractor_outputs[i] + masks[i] for i in range(len(ct))]
    msg = mixed_dot(intermediate_shares, weights)

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

        print_ln("-----TEST 5: STRICT MAJORITY BOUNDARY, SIMD-----")
        # In each lane, reject a different corrupted owner with exactly two
        # accepting verifiers. Every other owner has exactly three acceptances.
        for owner in range(4):
            share, keys, tags = shares[owner]
            source, ciphertext, seed_shares, mask_shares = share
            corrupted = [int(lane == owner) for lane in range(size)]
            share = (source, ciphertext + sgf2n(corrupted), seed_shares, mask_shares)
            message = _parse_mac_message(share)
            tags = []
            for verifier in range(4):
                mismatch = cgf2n([
                    int(verifier >= (2 if lane == owner else 3))
                    for lane in range(size)
                ])
                tags.append(mac(shares[verifier][1][owner], message) + mismatch)
            shares[owner] = (share, keys, tags)
        rec_msg = robust_lr_rec(shares, size=size)
        error_pattern = (rec_msg - msg).reveal()
        @if_e(error_pattern != cgf2n(0))
        def _():
            print_ln("❌ TEST 5 FAILED\nreconstructed message=%s\nexpected message=%s", rec_msg.reveal(), msg.reveal())
        @else_
        def _():
            print_ln("✅ TEST 5 PASSED")

        print_ln("-----TEST 6: PUBLIC ZERO CHECK, UPPER FIELD BITS-----")
        # Zero comparison must inspect the upper half of a 128-bit field too.
        column = [cgf2n([0, 1 << 64, 1 << 127, 1]), cgf2n([0, 0, 0, 0])]
        frequency = sum(entry == 0 for entry in column)
        valid = cgf2n(frequency >= 2, size=4)
        @if_e(valid != cgf2n([1, 0, 0, 0]))
        def _():
            print_ln("❌ TEST 6 FAILED: validity=%s", valid)
        @else_
        def _():
            print_ln("✅ TEST 6 PASSED")
    
    compiler.compile_func()
