#!/usr/bin/env python3

import os, sys
# add MP-SPDZ dir to path so we can import from Compiler
sys.path.insert(0, os.path.dirname(sys.argv[0]) + '/../..')

from Compiler.library import listen_for_clients, accept_client_connection, if_, public_input
from Compiler.types import cgf2n
from Compiler.compilerLib import Compiler

from lrss import get_source_length
from lrss_io import read_robust_share, reveal_and_encode_robust_share
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
    if t > n:
        raise ValueError("threshold cannot exceed num-parties")

    PORT_BASE = public_input()
    listen_for_clients(PORT_BASE)
    socket = accept_client_connection(PORT_BASE)

    # See lrss.py for notes on parsing input
    # Robust input parsing is shared with key generation through lrss_io.
    # Robust shares extend the LRSS layout with n MAC key pairs and n tags.
    # Every component of party i's share is read from player i.
    source_length = get_source_length(n, mu, secpar)

    # refresh
    old_shares = [
        read_robust_share(i, n, source_length, size)
        for i in range(n)
    ]
    secret = robust_lr_rec(old_shares, size=size)
    new_shares = robust_lr_share(secret, t, n, mu, secpar, size=size)

    # write back shares to appropriate parties
    # using term "block" for field element bc we are working in 128-bit
    # field (like AES)
    for i in range(n):
        # Flattening and socket encoding are shared with key generation in
        # lrss_io so both programs emit the same robust-share layout.
        write_back_vals = reveal_and_encode_robust_share(new_shares[i], i)
        @if_(i == socket)
        def _():
            cgf2n.write_to_socket(socket, write_back_vals)

if __name__ == "__main__":
    compiler.compile_func()
