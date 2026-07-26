#!/usr/bin/env python3

import os, sys
# add MP-SPDZ dir to path so we can import from Compiler
sys.path.insert(0, os.path.dirname(sys.argv[0]) + '/../..')

from Compiler.library import listen_for_clients, accept_client_connection, if_, public_input
from Compiler.types import cint, sgf2n
from Compiler.compilerLib import Compiler

from lrss import get_source_length
from robust_lrss import robust_lr_rec, robust_lr_share
from utils import get_random_sgf2n

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

@compiler.register_function('lrpss')
def lrpss():
    '''
    Refresh SV shares by reconstructing and resharing. 
    Assumes 128-bit field. 
    '''
    opt = compiler.options
    args = ("t", "n", "mu", "secpar", "size")
    t, n, mu, secpar, size = tuple(map(lambda x : int(getattr(opt, x)), args))

    PORT_BASE = public_input()
    listen_for_clients(PORT_BASE)
    socket = accept_client_connection(PORT_BASE)

    # See lrss.py for notes on parsing input
    source_length = seed_length = get_source_length(n, mu, secpar)
    share_length = source_length + 1 + seed_length + n
    def parse_input(i):
        src = [sgf2n.get_input_from(i, size=size) for _ in range(source_length)]
        ct = sgf2n.get_input_from(i, size=size)
        ss = [sgf2n.get_input_from(i, size=size) for _ in range(seed_length)]
        mst = [sgf2n.get_input_from(i, size=size) for _ in range(n)]
        lr_share = (src, ct, ss, mst)
        keys = [
            (sgf2n.get_input_from(i, size=size), sgf2n.get_input_from(i, size=size))
            for i in range(n)
        ]
        tags = [sgf2n.get_input_from(i, size=size) for i in range(n)]
        return (lr_share, keys, tags)

    # refresh
    old_shares = [parse_input(i) for i in range(n)]
    secret = robust_lr_rec(old_shares, size=size)
    new_shares = robust_lr_share(secret, t, n, mu, secpar, size=size)

    # write back shares to appropriate parties
    # using term "block" for field element bc we are working in 128-bit
    # field (like AES)
    for i in range(n):
        lr_share, keys, tags = new_shares[i]
        src, ct, ss, mst = lr_share 
        flat_share = [*src, ct, *ss, *mst, *keys, *tags]
        flat_share_personal = [block.reveal_to(i) for block in flat_share]
        @if_(i == socket)
        def _():
            write_back_vals = [cint(block._v) for block in flat_share_personal]
            cint.write_to_socket(socket, flat_share_personal)

if __name__ == "__main__":
    compiler.compile_func()
