#!/usr/bin/env python3

import os, sys
# add MP-SPDZ dir to path so we can import from Compiler
sys.path.insert(0, os.path.dirname(sys.argv[0]) + '/../..')

from Compiler.library import listen_for_clients, accept_client_connection, if_, public_input
from Compiler.types import cint, sgf2n
from Compiler.compilerLib import Compiler

from shamir import shamir_share
from embeddings import apply_field_embedding, apply_inverse_field_embedding
from utils import get_random_sgf2n

usage = "usage: %prog [options] [args]"
compiler = Compiler(usage=usage)
compiler.parser.add_option("--key-length", dest="key_len")
compiler.parser.add_option("--threshold", dest="t")
compiler.parser.add_option("--num-parties", dest="n")
compiler.parse_args()
if not compiler.options.key_len:
    compiler.parser.error("--key-length required")
if not compiler.options.t:
    compiler.parser.error("--threshold required")
if not compiler.options.n:
    compiler.parser.error("--num-parties required")

@compiler.register_function('pss')
def pss():
    '''
    We now refresh the existing Shamir shares by adding a fresh random sharing of zero, instead of reconstructing the secret.
    '''
    key_len, t, n = int(compiler.options.key_len), int(compiler.options.t), int(compiler.options.n)
    num_bytes = key_len // 8

    PORT_BASE = public_input()
    listen_for_clients(PORT_BASE)
    socket = accept_client_connection(PORT_BASE)

    input_shares = [sgf2n.get_input_from(i, size=num_bytes) for i in range(n)]
    
    input_shares_embedded = [apply_field_embedding(share) for share in input_shares]
    eval_points_embedded = [apply_field_embedding(sgf2n(i)) for i in range(1, n + 1)]

    # Sample one sharing of zero and add it to the existing shares.
    zero_embedded = apply_field_embedding(sgf2n(0, size=num_bytes))
    randomness_embedded = [
        apply_field_embedding(get_random_sgf2n(8, size=num_bytes))
        for _ in range(t)
    ]
    _, zero_shares_embedded = shamir_share(
        msg=zero_embedded,
        threshold=t,
        num_parties=n,
        eval_points=eval_points_embedded,
        rand=randomness_embedded,
        size=num_bytes,
    )

    #  Addition of zero shares refreshes the shares while preserving the secret.
    refreshed_shares_by_party = [
        apply_inverse_field_embedding(
            input_shares_embedded[party] + zero_shares_embedded[party]
        ).reveal_to(party)
        for party in range(n)
    ]

    for party in range(n):
        @if_(party == socket)
        def _():
            refreshed_share = refreshed_shares_by_party[party]
            
            byte_values = [cint(refreshed_share[i]._v) for i in range(num_bytes)]
            cint.write_to_socket(socket, byte_values)

if __name__ == "__main__":
    compiler.compile_func()
