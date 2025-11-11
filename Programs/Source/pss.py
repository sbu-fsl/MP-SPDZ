#!/usr/bin/env python3

import os, sys
# add MP-SPDZ dir to path so we can import from Compiler
sys.path.insert(0, os.path.dirname(sys.argv[0]) + '/../..') 

from Compiler.library import print_ln, listen_for_clients, accept_client_connection, if_, public_input
from Compiler.types import cint, Array, sgf2n
from Compiler.compilerLib import Compiler

# we assume these modules reside in Programs/Source/ 
from shamir import shamir_share, shamir_reconstruct
from embeddings import apply_field_embedding, apply_inverse_field_embedding

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
    Perform the proactive secret sharing functionality on the set of input shares. That is,
    reconstruct the secret from the input shares and output fresh secret shares. Each client
    actually inputs a list of secret shares whose concatenation corresponds to a single share 
    of a key. PSS is performed independently "column-wise". In other words, visualize a matrix
    of shares where each row is a client's input list. Each column is a set of secret shares that 
    gets refreshed. 
    '''
    # constants
    key_len, t, n = int(compiler.options.key_len), int(compiler.options.t), int(compiler.options.n)
    num_bytes = key_len // 8

    # set up external client connections
    PORT_BASE = public_input()
    listen_for_clients(PORT_BASE)
    socket = accept_client_connection(PORT_BASE)

    input_shares = [sgf2n.get_input_from(i, size=num_bytes) for i in range(n)] # read from Player-Data/Input-P<player>-<thread> in HEX FORMAT
    input_shares_embedded = [[apply_field_embedding(y) for y in vectorized_share] for vectorized_share in input_shares] # need to "un-vectorize" sgf2n before applying embedding... for now
    
    # eval points need to be embedded since they participate in arithmetic with embedded key elements. 
    eval_points_embedded = Array(n, sgf2n).assign([apply_field_embedding(sgf2n(i)) for i in range(1,n+1)])

    # reconstruct secret
    secret_embedded = [shamir_reconstruct(list(ys), eval_points=eval_points_embedded) for ys in zip(*input_shares_embedded)]

    # reshare secret one byte at a time, and group new shares by party. 
    new_shares_by_party = {party: [] for party in range(n)}
    for byte_idx in range(num_bytes):
        # need to make sure random field elements used in shamir_share are also embedded field elements. 
        # have to do this inside for loop to ensure we don't reuse randomness across shamir_share() calls
        randomness_embedded = [apply_field_embedding(sgf2n.bit_compose([sgf2n.get_random_bit() for _ in range(8)])) for i in range(t)]
        randomness_embedded = Array(t,sgf2n).assign(randomness_embedded)  # convert to Array since shamir_share expects Array (for now)
        byte_shares = shamir_share(
            msg=secret_embedded[byte_idx], 
            threshold=t, 
            num_parties=n, 
            eval_points=eval_points_embedded,
            rand=randomness_embedded)[1]
        for party, new_share in enumerate(byte_shares):
            new_share = apply_inverse_field_embedding(new_share)
            new_shares_by_party[party].append(new_share.reveal_to(party)) 

    # write shares back to corresponding parties
    for party in range(n):
        @if_(party == socket)
        def _():
            byte_values = [cint(value._v) for value in new_shares_by_party[party]]
            cint.write_to_socket(socket, byte_values)
    
if __name__ == "__main__":
    compiler.compile_func()
