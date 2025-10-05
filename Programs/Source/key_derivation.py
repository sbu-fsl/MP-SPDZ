#!/usr/bin/env python3

import os, sys
# add MP-SPDZ dir to path so we can import from Compiler
sys.path.insert(0, os.path.dirname(sys.argv[0]) + '/../..') 

from Compiler.library import print_ln, listen_for_clients, accept_client_connection, if_, public_input
from Compiler.types import cint, Array, sgf2n, cgf2n
from Compiler.compilerLib import Compiler

# we assume these modules reside in Programs/Source/ 
from shamir import shamir_share, shamir_reconstruct
from embeddings import apply_field_embedding, apply_inverse_field_embedding
from aes_ctr import aes_ctr_encrypt


usage = "usage: %prog [options] [args]"
compiler = Compiler(usage=usage)
compiler.parser.add_option("--root-key-length", dest="root_key_len")
compiler.parser.add_option("--child-key-length", dest="child_key_len")
compiler.parser.add_option("--read-nonce-from-file", action="store_true", dest="nonce")
# TODO: for now assuming same access structure for root key and child key
compiler.parser.add_option("--threshold", dest="t")
compiler.parser.add_option("--num-parties", dest="n")
compiler.parse_args()
if not compiler.options.root_key_len:
    compiler.parser.error("--root-key-length required")
if not compiler.options.child_key_len:
    compiler.parser.error("--child-key-length required")
if not compiler.options.t:
    compiler.parser.error("--threshold required")
if not compiler.options.n:
    compiler.parser.error("--num-parties required")

@compiler.register_function('derive')
def derive():
    '''
    Derive a child key from a secret-shared root key and a child key ID using AES-CTR. 
    Outputs shares of the child key to each party using the same access structure as the root key. 
    '''
    # constants
    opts = compiler.options
    [root_key_len, child_key_len, t, n] = [int(x) for x in [opts.root_key_len, opts.child_key_len, opts.t, opts.n]] 
    num_bytes_root_key = root_key_len // 8
    num_bytes_child_key = child_key_len // 8

    # set up external client connections
    PORT_BASE = 15000
    listen_for_clients(PORT_BASE)
    socket = accept_client_connection(PORT_BASE)

    # runtime inputs
    # nonces should be written to Programs/Public-Input/<progname> in integer format (public_input returns cint)
    nonce = [cgf2n(public_input()) for _ in range(12)] if opts.nonce else None
    # TODO: make sure child_key_id is 16-byte (i.e. 128-bit) aligned.
    child_key_id = [cgf2n(public_input()) for _ in range(num_bytes_child_key)] 

    ### Step 1: Reconstruct root key
    input_shares = [sgf2n.get_input_from(i, size=num_bytes_root_key) for i in range(n)] # read from Player-Data/Input-P<player>-<thread> in HEX FORMAT
    input_shares_embedded = [[apply_field_embedding(y) for y in vectorized_share] for vectorized_share in input_shares] # need to "un-vectorize" sgf2n before applying embedding... for now
    # eval points need to be embedded since they participate in arithmetic with embedded key elements. 
    eval_points_embedded = Array(n, sgf2n).assign([apply_field_embedding(sgf2n(i)) for i in range(1,n+1)])
    # reconstruct root key
    root_key = [apply_inverse_field_embedding(shamir_reconstruct(list(ys), eval_points=eval_points_embedded)) for ys in zip(*input_shares_embedded)]

    ### Step 2: use root key and aes encryption to derive child key (aes encryption plays the role of PRF)
    nonce, child_key = aes_ctr_encrypt(root_key, child_key_id, nonce=nonce)
    # dbg_key = [x.reveal() for x in child_key]
    # print_ln("dbg_key=%s", dbg_key)
    child_key = [apply_field_embedding(x) for x in child_key] # embed child key for next step

    ### Step 3: Secret share child_key (client will have to manually get shares from the nodes)
    # reshare secret one byte at a time, and group new shares by party. 
    child_key_shares_by_party = {party: [] for party in range(n)}
    for byte_idx in range(num_bytes_child_key):
        # need to make sure random field elements used in shamir_share are also embedded field elements. 
        # have to do this inside for loop to ensure we don't reuse randomness across shamir_share() calls
        randomness_embedded = [apply_field_embedding(sgf2n.bit_compose([sgf2n.get_random_bit() for _ in range(8)])) for i in range(t)]
        randomness_embedded = Array(t,sgf2n).assign(randomness_embedded)  # convert to Array since shamir_share expects Array (for now)
        byte_shares = shamir_share(
            msg=child_key[byte_idx], 
            threshold=t, 
            num_parties=n, 
            eval_points=eval_points_embedded,
            rand=randomness_embedded)[1]
        for party, byte_share in enumerate(byte_shares):
            byte_share = apply_inverse_field_embedding(byte_share)
            # child_key_shares_by_party[party].append(byte_share.reveal())
            child_key_shares_by_party[party].append(byte_share.reveal_to(party)) 
    
    # for party in range(n):
    #     print_ln("child_key_shares_by_party=%s", child_key_shares_by_party[party])

    ## Step 4: Write nonce and child key shares back to corresponding parties
    for party in range(n):
        @if_(party == socket)
        def _():
            for nonce_byte in nonce:
                cint.write_to_socket(socket, cint(nonce_byte))
            child_key_share = child_key_shares_by_party[party]
            for byte_share in child_key_share:
                cint.write_to_socket(socket, cint(byte_share._v)) 


if __name__ == "__main__":
    compiler.compile_func()
