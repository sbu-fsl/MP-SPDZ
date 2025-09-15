#!/usr/bin/env python3

import os, sys
# add MP-SPDZ dir to path so we can import from Compiler
sys.path.insert(0, os.path.dirname(sys.argv[0]) + '/../..') 

from Compiler.library import print_ln, listen_for_clients, accept_client_connection, for_range, get_player_id, if_, public_input
from Compiler.types import sint, cint, regint, Array, Matrix, ClientMessageType, sgf2n
from Compiler.compilerLib import Compiler

# we assume these modules reside in Programs/Source/ 
from shamir import shamir_share
from aes import apply_field_embedding, apply_inverse_field_embedding

MAX_KEY_BITS = 128
MAX_KEY_BYTES = MAX_KEY_BITS // 8
MAX_N = 3
MAX_T = 2

usage = "usage: %prog [options] [args]"
compiler = Compiler(usage=usage)

@compiler.register_function('root_key_gen')
def root_key_gen():
    '''
    Generate a root key and secret-share it according to the key generation parameters 
    that have been passed into the compiler via command-line args. Each secret share
    is sent back to its corresponding party. 
    '''
    # ---- Runtime public parameters ----
    key_len_bits = public_input()   # public int
    t = public_input()   # public int
    n = public_input()   # public int
    key_len_bytes = key_len_bits // 8

    print("Read from public_input")

    # set up external client connections
    PORT_BASE = 15000
    listen_for_clients(PORT_BASE)
    socket = accept_client_connection(PORT_BASE)

    # generate uniformly random key and embed each byte of the key
    key = [sgf2n.bit_compose([sgf2n.get_random_bit() for _ in range(8)]) for _ in range(MAX_KEY_BYTES)]
    key_embedded = [ apply_field_embedding(byte) for byte in key]

    # eval points need to be embedded since they participate in arithmetic with embedded key elements. 
    eval_points_embedded = [apply_field_embedding(sgf2n(i)) for i in range(1,MAX_N+1)]
    eval_points_embedded = Array(MAX_N,sgf2n).assign(eval_points_embedded) # convert to Array since shamir_share expects Array (for now)

    # secret share each byte, then group shares by party 
    shares_by_party = {party: [] for party in range(MAX_N)}
    for byte_idx in range(MAX_KEY_BYTES):
        # need to make sure random field elements used in shamir_share are also embedded field elements. 
        # have to do this inside for loop to ensure we don't reuse randomness across shamir_share() calls
        randomness_embedded = [apply_field_embedding(sgf2n.bit_compose([sgf2n.get_random_bit() for _ in range(8)])) for i in range(MAX_T)]
        randomness_embedded = Array(MAX_T,sgf2n).assign(randomness_embedded)  # convert to Array since shamir_share expects Array (for now)
        byte_shares = shamir_share(
            msg=key_embedded[byte_idx], 
            threshold=MAX_T, 
            num_parties=MAX_N, 
            eval_points=eval_points_embedded,
            rand=randomness_embedded
        )[1] # only want polynomial evaluations, not evaluation points
        for party, share in enumerate(byte_shares):
            share = apply_inverse_field_embedding(share)
            shares_by_party[party].append(share.reveal_to(party)) 
    
    # write shares back to corresponding parties
    for party in range(MAX_N):
        @if_((party < n) & (party == socket))
        def _send_to_socket():
            for i in range(MAX_KEY_BYTES):
                @if_(i < key_len_bytes)
                def _maybe_send():
                    cint.write_to_socket(socket, cint(shares_by_party[party][i]._v))

if __name__ == "__main__":
    compiler.compile_func()