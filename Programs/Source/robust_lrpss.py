#!/usr/bin/env python3

import os, sys
# add MP-SPDZ dir to path so we can import from Compiler
sys.path.insert(0, os.path.dirname(sys.argv[0]) + '/../..')

from Compiler.library import listen_for_clients, accept_client_connection, public_input
from Compiler.compilerLib import Compiler

from lrss import get_source_length
from share_io import read_blocks, write_shares_to_socket
from robust_lrss import robust_lr_rec, robust_lr_share

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

    # refresh
    old_shares = []
    base_share_length = 2 * source_length + n + 1
    for i in range(n):
        blocks = read_blocks(i, base_share_length + 3 * n, size)
        lrss_blocks = blocks[:base_share_length]
        mac_key_blocks = blocks[base_share_length:base_share_length + 2 * n]
        old_shares.append(
            (
                (
                    lrss_blocks[:source_length],
                    lrss_blocks[source_length],
                    lrss_blocks[source_length + 1:2 * source_length + 1],
                    lrss_blocks[2 * source_length + 1:],
                ),
                [
                    tuple(mac_key_blocks[j:j + 2])
                    for j in range(0, 2 * n, 2)
                ],
                blocks[base_share_length + 2 * n:],
            )
        )
    secret = robust_lr_rec(old_shares, size=size)
    new_shares = robust_lr_share(secret, t, n, mu, secpar, size=size)

    write_shares_to_socket(socket, new_shares)

if __name__ == "__main__":
    compiler.compile_func()
