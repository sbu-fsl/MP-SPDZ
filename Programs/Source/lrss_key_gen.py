#!/usr/bin/env python3

import os, sys
# add MP-SPDZ dir to path so we can import from Compiler
sys.path.insert(0, os.path.dirname(sys.argv[0]) + '/../..')

from Compiler.library import listen_for_clients, accept_client_connection, if_, public_input
from Compiler.types import cint, sgf2n
from Compiler.compilerLib import Compiler

from lrss import lr_share, lr_rec, get_source_length
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
    Generate lrss shares of a uniform 128-bit key (i.e., one field element).  
    '''
    opt = compiler.options
    args = (t, n, mu, secpar, size)
    t, n, mu, secpar, size = tuple(map(lambda x : int(getattr(opt, x)), args))

    PORT_BASE = public_input()
    listen_for_clients(PORT_BASE)
    socket = accept_client_connection(PORT_BASE)

    secret = get_random_sgf2n(128, size=size)
    shares = lr_share(secret, t, n, mu, secpar, size=size)

    for party in range(n):
        @if_(party == socket)
        def _():
            # using "block" to remind each item in share is 128-bits (like AES)
            # TODO: check github issue to see if we can use cgf2n now. 
            share_values = [cint(block._v) for block in shares[party]]
            cint.write_to_socket(socket, share_values)

if __name__ == "__main__":
    compiler.compile_func()
