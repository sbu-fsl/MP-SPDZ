#!/usr/bin/env python3

import os, sys
# add MP-SPDZ dir to path so we can import from Compiler
sys.path.insert(0, os.path.dirname(sys.argv[0]) + '/../..') 

from Compiler.library import print_ln, listen_for_clients, accept_client_connection, for_range, get_player_id, if_
from Compiler.types import sint, cint, regint, Array, Matrix, ClientMessageType, sgf2n
from Compiler.compilerLib import Compiler

# we assume these modules reside in Programs/Source/ 
from shamir import shamir_share, shamir_reconstruct
from aes import apply_field_embedding, apply_inverse_field_embedding

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
    # get key_gen_params from command-line args
    key_len, t, n = int(compiler.options.key_len), int(compiler.options.t), int(compiler.options.n)

    # set up external client connections
    PORT_BASE = 15000
    listen_for_clients(PORT_BASE)
    socket = accept_client_connection(PORT_BASE)

    # get a list of shares from each client, think it has to be through Player-Data/ since 
    # only method in client interface for sending private data to parties is client.send_private_inputs(values:list) ,
    # but this assumes client is connected to all parties. 
    input_shares = [sgf2n.get_input_from(i) for i in range(n)]


    
    # secret share each byte, then group shares by party 
    shares_by_party = {party: [] for party in range(n)}
    for byte_idx in range(key_len // 8):
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
            for share in shares_by_party[party]:
                cint.write_to_socket(socket, cint(share._v))

if __name__ == "__main__":
    compiler.compile_func()
