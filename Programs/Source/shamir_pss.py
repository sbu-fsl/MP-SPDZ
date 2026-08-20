#!/usr/bin/env python3

import os, sys
# add MP-SPDZ dir to path so we can import from Compiler
sys.path.insert(0, os.path.dirname(sys.argv[0]) + '/../..')

from Compiler.library import listen_for_clients, accept_client_connection, public_input
from Compiler.types import sgf2n
from Compiler.compilerLib import Compiler

from share_io import read_blocks, write_shares_to_socket
from shamir import shamir_share

usage = "usage: %prog [options] [args]"
compiler = Compiler(usage=usage)
compiler.parser.add_option("--threshold", dest="t")
compiler.parser.add_option("--num-parties", dest="n")
compiler.parser.add_option("--size", dest="size", default="1")
compiler.parse_args()
if not compiler.options.t:
    compiler.parser.error("--threshold required")
if not compiler.options.n:
    compiler.parser.error("--num-parties required")

@compiler.register_function('shamir_pss')
def shamir_pss():
    '''
    Refresh Shamir shares using the "shares of zero" trick.
    Assumes a 128-bit field.
    '''
    opt = compiler.options
    args = ("t", "n", "size")
    t, n, size = (
        int(getattr(opt, name))
        for name in args
    )
    if min(t, n, size) <= 0:
        raise ValueError("threshold, num-parties, and size must be positive")
    if t >= n:
        raise ValueError("threshold must be less than num-parties")

    PORT_BASE = public_input()
    listen_for_clients(PORT_BASE)
    socket = accept_client_connection(PORT_BASE)

    input_shares = [read_blocks(i, 1, size)[0] for i in range(n)]
    _, zero_shares = shamir_share(sgf2n(0, size=size), t, n, size=size)
    refreshed_shares = [
        share + zero_share
        for share, zero_share in zip(input_shares, zero_shares)
    ]

    write_shares_to_socket(socket, refreshed_shares)

if __name__ == "__main__":
    compiler.compile_func()
