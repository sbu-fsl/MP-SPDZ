#!/usr/bin/env python3

import os, sys
# add MP-SPDZ dir to path so we can import from Compiler
sys.path.insert(0, os.path.dirname(sys.argv[0]) + '/../..')

from Compiler.library import listen_for_clients, accept_client_connection, public_input
from Compiler.compilerLib import Compiler

from share_io import write_shares_to_socket
from shamir import shamir_share
from utils import get_random_sgf2n

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

@compiler.register_function('shamir_key_gen')
def shamir_key_gen():
    '''Generate and distribute Shamir shares of ``size`` random 128-bit keys.'''
    t = int(compiler.options.t)
    n = int(compiler.options.n)
    size = int(compiler.options.size)
    if min(t, n, size) <= 0:
        raise ValueError("threshold, num-parties, and size must be positive")
    if t >= n:
        raise ValueError("threshold must be less than num-parties")

    PORT_BASE = public_input()
    listen_for_clients(PORT_BASE)
    socket = accept_client_connection(PORT_BASE)

    secret = get_random_sgf2n(128, size=size)
    _, shares = shamir_share(secret, t, n, size=size)
    write_shares_to_socket(socket, shares)

if __name__ == "__main__":
    compiler.compile_func()
