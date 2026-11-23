#!/usr/bin/env python3

import os, sys
# add MP-SPDZ dir to path so we can import from Compiler
sys.path.insert(0, os.path.dirname(sys.argv[0]) + '/../..')

from Compiler.library import listen_for_clients, accept_client_connection, public_input, start_timer, stop_timer
from Compiler.compilerLib import Compiler

from lrss import get_source_length
from share_io import read_blocks, write_shares_to_socket
from robust_lrss import robust_lr_rec, robust_lr_share, get_mac_repetitions

usage = "usage: %prog [options] [args]"
compiler = Compiler(usage=usage)
compiler.parser.add_option("--threshold", dest="t")
compiler.parser.add_option("--num-parties", dest="n")
compiler.parser.add_option("--mu", dest="mu")
compiler.parser.add_option("--secpar", dest="secpar")
compiler.parser.add_option("--size", dest="size") # parallelization parameter
compiler.parse_args()
if not compiler.options.t:
    compiler.parser.error("--threshold required")
if not compiler.options.n:
    compiler.parser.error("--num-parties required")
if not compiler.options.mu:
    compiler.parser.error("--mu required")
if not compiler.options.secpar:
    compiler.parser.error("--secpar required")
if not compiler.options.size:
    compiler.parser.error("--size required")

@compiler.register_function('robust_lrpss')
def robust_lrpss():
    '''
    Refresh robust LRSS shares by reconstructing and resharing.
    Assumes 128-bit field. 
    '''
    opt = compiler.options
    args = ("t", "n", "mu", "secpar", "size")
    t, n, mu, secpar, size = tuple(map(lambda x : int(getattr(opt, x)), args))
    if min(t, n, secpar, size) <= 0 or mu < 0:
        raise ValueError(
            "threshold, num-parties, secpar, and size must be positive; "
            "mu must be non-negative"
        )
    if t >= n:
        raise ValueError("threshold must be less than num-parties")

    PORT_BASE = public_input()
    listen_for_clients(PORT_BASE)
    socket = accept_client_connection(PORT_BASE)

    source_length = get_source_length(n, mu, secpar)

    start_timer(11)
    old_shares = []
    base_share_length = 2 * source_length + n + 1
    c = get_mac_repetitions(t, base_share_length, mu, secpar)
    key_length = base_share_length + 1
    keys_end = base_share_length + n * c * key_length
    tags_end = keys_end + n * c
    key_count = n * n * c
    new_key_blocks = [0] * (2 * key_count)
    for i in range(n):
        blocks = read_blocks(i, tags_end + len(new_key_blocks), size)
        lrss_blocks = blocks[:base_share_length]
        mac_key_blocks = blocks[base_share_length:keys_end]
        tag_blocks = blocks[keys_end:tags_end]
        # Field addition is XOR. One independent honest contribution makes the
        # sum uniform, with every summand kept secret inside MPC.
        new_key_blocks = [total + value for total, value in zip(new_key_blocks, blocks[tags_end:])]
        old_shares.append(
            (
                (
                    lrss_blocks[:source_length],
                    lrss_blocks[source_length],
                    lrss_blocks[source_length + 1:2 * source_length + 1],
                    lrss_blocks[2 * source_length + 1:],
                ),
                [
                    [tuple(mac_key_blocks[j:j + key_length])
                     for j in range(owner * c * key_length, (owner + 1) * c * key_length, key_length)]
                    for owner in range(n)
                ],
                [tag_blocks[verifier * c:(verifier + 1) * c]
                 for verifier in range(n)],
            )
        )
    stop_timer(11)
    start_timer(12)
    secret = robust_lr_rec(old_shares, size=size, expanded_keys=True)
    stop_timer(12)
    start_timer(13)
    key_blocks = iter(new_key_blocks)
    new_keys = {(verifier, owner): [(next(key_blocks), next(key_blocks)) for _ in range(c)]
                for verifier in range(n) for owner in range(n)}
    new_shares = robust_lr_share(secret, t, n, mu, secpar, size=size, mac_keys=new_keys)
    stop_timer(13)

    start_timer(14)
    write_shares_to_socket(socket, new_shares)
    stop_timer(14)

if __name__ == "__main__":
    compiler.compile_func()
