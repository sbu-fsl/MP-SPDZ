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
    Refresh SV shares by reconstructing and resharing. 
    Assumes 128-bit field. 
    '''
    opt = compiler.options
    args = (t, n, mu, secpar, size)
    t, n, mu, secpar, size = tuple(map(lambda x : int(getattr(opt, x)), args))

    PORT_BASE = public_input()
    listen_for_clients(PORT_BASE)
    socket = accept_client_connection(PORT_BASE)

    '''
    Parsing input shares is only mildly involved. We want to parse each
    share as (source, ct, seed_share, mask_shares_transposed) where source,
    seed_share, and mask_shares_transposed are all list[sgf2n]. We always have
    `len(mask_shares_transposed) == n` and `len(source) == len(seed_share)`, but
    value of source length is determined by `get_source_length(n,mu,secpar)`. 
    
    The only tool we have to read secret inputs is
    `sgf2n.get_input_from(i,size=size)`, which reads `size`-many field elements
    at a time from whatever input method is specified at runtime. See
    https://mp-spdz.readthedocs.io/en/latest/io.html#private-inputs-from-computing-parties.
    However, we want to reserve the size argument for parallelization, i.e.,
    size=1000 means we are running lrpss on a batch of 1000 inputs at once.
    Thus, for each party i we want to call sgf2n.get_input_from(i, size=size)
    `share_length`-many times, where `share_length = source_length + 1 +
    seed_length + n`.
    '''
    source_length = seed_length = get_source_length(n, mu, secpar)
    share_length = source_length + 1 + seed_length + n
    def parse_input(i):
        source = [sgf2n.get_input_from(i, size=size) for _ in range(source_length)]
        ct = sgf2n.get_input_from(i, size=size)
        seed_share = [sgf2n.get_input_from(i, size=size) for _ in range(seed_length)]
        mask_shares_transposed = [sgf2n.get_input_from(i, size=size) for _ in range(n)]
        return (source, ct, seed_share, mask_shares_transposed)

    old_shares = [parse_input(i) for i in range(n)]
    secret = lr_rec(old_shares, size=size)
    new_shares = lr_share(secret, t, n, mu, secpar, size=size)

    for party in range(n):
        @if_(party == socket)
        def _():
            # using "block" to remind each item in share is 128-bits (like AES)
            # TODO: check github issue to see if we can use cgf2n now. 
            new_share_values = [cint(block._v) for block in new_shares[party]]
            cint.write_to_socket(socket, new_share_values)

if __name__ == "__main__":
    compiler.compile_func()
