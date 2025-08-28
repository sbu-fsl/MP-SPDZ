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
# compiler.parser.add_option("--key-length", dest="key_len")
# compiler.parser.add_option("--threshold", dest="t")
# compiler.parser.add_option("--num-parties", dest="n")
# compiler.parse_args()
# if not compiler.options.key_len:
#     compiler.parser.error("--key-length required")
# if not compiler.options.t:
#     compiler.parser.error("--threshold required")
# if not compiler.options.n:
#     compiler.parser.error("--num-parties required")

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

    # key_len, t, n = int(compiler.options.key_len), int(compiler.options.t), int(compiler.options.n)
    key_len, t, n = 128, 2, 3
    # TODO: in order to use key_len, we want to pass size= key_len // 8 into the size parameter for get_input_from below.
    # Unfortunately, dependencies like field embeddings and shamir do not support vectorized types like this yet. 
    # Can we coerce vector into list? Maybe we can do Array.get_vector() first and then put into list? 

    # set up external client connections
    # PORT_BASE = 15000
    # listen_for_clients(PORT_BASE)
    # socket = accept_client_connection(PORT_BASE)

    input_shares = [sgf2n.get_input_from(i, size=1) for i in range(n)] # read from Player-Data/Input-P<player>-<thread>
    input_shares_embedded = [apply_field_embedding(x) for x in input_shares]

    # reconstruct secret
    secret_embedded = shamir_reconstruct(input_shares_embedded)

    # reshare secret
    new_shares_embedded = shamir_share(secret_embedded, t, n)[1]
    new_shares = [apply_inverse_field_embedding(x) for x in new_shares_embedded]
    new_shares_personal = [share.reveal_to(i+1) for i, share in enumerate(new_shares)]

    # write shares back to corresponding parties
    # for party in range(n):
    #     @if_(party == socket)
    #     def _():
    #         new_share = new_shares_personal[party]
    #         cint.write_to_socket(socket, cint(new_share._v)) 
    
if __name__ == "__main__":
    compiler.compile_func()
