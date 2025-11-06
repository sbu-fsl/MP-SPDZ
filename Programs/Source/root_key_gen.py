#!/usr/bin/env python3

import os, sys
# add MP-SPDZ dir to path so we can import from Compiler
sys.path.insert(0, os.path.dirname(sys.argv[0]) + '/../..') 

from Compiler.library import print_ln, listen_for_clients, accept_client_connection, for_range, get_player_id, if_, public_input
from Compiler.types import sint, cint, regint, Array, Matrix, ClientMessageType, sgf2n
from Compiler.compilerLib import Compiler

# we assume these modules reside in Programs/Source/ 
from shamir import shamir_share
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

@compiler.register_function('root_key_gen')
def root_key_gen():
    '''
    Generate a root key and secret-share it according to the key generation parameters 
    that have been passed into the compiler via command-line args. Each secret share
    is sent back to its corresponding party. 
    '''
    # get key_gen_params from command-line args
    key_len, t, n = int(compiler.options.key_len), int(compiler.options.t), int(compiler.options.n)

    # set up external client connections
    PORT_BASE = public_input()
    listen_for_clients(PORT_BASE)
    socket = accept_client_connection(PORT_BASE)

    # generate uniformly random key and embed each byte of the key
    key = [sgf2n.bit_compose([sgf2n.get_random_bit() for _ in range(8)]) for _ in range(key_len // 8)]
    key_embedded = [ apply_field_embedding(byte) for byte in key]

    # eval points need to be embedded since they participate in arithmetic with embedded key elements. 
    eval_points_embedded = [apply_field_embedding(sgf2n(i)) for i in range(1,n+1)]
    eval_points_embedded = Array(n,sgf2n).assign(eval_points_embedded) # convert to Array since shamir_share expects Array (for now)

    # secret share each byte, then group shares by party 
    shares_by_party = {party: [] for party in range(n)}
    for byte_idx in range(key_len // 8):
        # need to make sure random field elements used in shamir_share are also embedded field elements. 
        # have to do this inside for loop to ensure we don't reuse randomness across shamir_share() calls
        randomness_embedded = [apply_field_embedding(sgf2n.bit_compose([sgf2n.get_random_bit() for _ in range(8)])) for i in range(t)]
        randomness_embedded = Array(t,sgf2n).assign(randomness_embedded)  # convert to Array since shamir_share expects Array (for now)
        byte_shares = shamir_share(
            msg=key_embedded[byte_idx], 
            threshold=t, 
            num_parties=n, 
            eval_points=eval_points_embedded,
            rand=randomness_embedded
        )[1] # only want polynomial evaluations, not evaluation points
        for party, share in enumerate(byte_shares):
            share = apply_inverse_field_embedding(share)
            shares_by_party[party].append(share.reveal_to(party)) 
    
    # write shares back to corresponding parties
    for party in range(n):
        @if_(party == socket)
        def _():
            byte_values = [cint(value) for value in shares_by_party[party]]
            cint.write_to_socket(socket, byte_values)

if __name__ == "__main__":
    compiler.compile_func()