#!/usr/bin/env python3

import os, sys
# add MP-SPDZ dir to path so we can import from Compiler
sys.path.insert(0, os.path.dirname(sys.argv[0]) + '/../..')

from Compiler.library import listen_for_clients, accept_client_connection, if_, public_input
from Compiler.types import cint
from Compiler.compilerLib import Compiler

from lrss import lr_share, lr_rec, get_source_length
from lrss_io import read_share, reveal_and_encode_share

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

@compiler.register_function('lrpss')
def lrpss():
    '''
    Refresh SV shares by reconstructing and resharing. 
    Assumes 128-bit field. 
    '''
    opt = compiler.options
    args = ("t", "n", "mu", "secpar", "size")
    t, n, mu, secpar, size = (
        int(getattr(opt, name))
        for name in args
    )
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

    '''
    We want to parse each share as 
    (source, ct, seed_share, mask_shares_transposed) 
    where source, seed_share, and mask_shares_transposed
    are all list[sgf2n], and ct is a single sgf2n. 
    We always have `len(mask_shares_transposed) == n` and
    `len(source) == len(seed_share)`, but value of source length is determined
    by `get_source_length(n,mu,secpar)`. 
    
    The only tool we have to read secret sgf2n inputs is
    `sgf2n.get_input_from(i,size=size)`, which reads `size`-many field elements
    at a time from whatever input method is specified at runtime. See
    https://mp-spdz.readthedocs.io/en/latest/io.html#private-inputs-from-computing-parties.
    However, we want to reserve the size argument for parallelization, i.e.,
    size=1000 means we are running lrpss on a batch of 1000 inputs at once.
    Thus, for each party i we want to call sgf2n.get_input_from(i, size=size)
    `share_length`-many times, where `share_length = source_length + 1 +
    seed_length + n`.
    '''

    # The former inline parse_input() helper now lives in lrss_io.read_share().
    source_length = get_source_length(n, mu, secpar)

    # refresh is short and sweet
    old_shares = [
        read_share(i, n, source_length, size)
        for i in range(n)
    ]
    secret = lr_rec(old_shares, size=size)
    new_shares = lr_share(secret, t, n, mu, secpar, size=size)

    # write back shares to appropriate parties
    # using term "block" for field element bc we are working in 128-bit
    # field (like AES)
    # The former inline flatten/reveal logic now lives in
    # lrss_io.reveal_and_encode_share(), together with the encoding used for
    # socket output.
    for i in range(n):
        output = reveal_and_encode_share(new_shares[i], i, size)
        @if_(i == socket)
        def _():
            cint.write_to_socket(socket, output)

if __name__ == "__main__":
    compiler.compile_func()
