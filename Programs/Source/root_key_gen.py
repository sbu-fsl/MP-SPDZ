#!/usr/bin/env python3

import os, sys
# add MP-SPDZ dir to path so we can import from Compiler
sys.path.insert(0, os.path.dirname(sys.argv[0]) + '/../..') 

from Compiler.library import print_ln, listen_for_clients, accept_client_connection, for_range, get_player_id, if_
from Compiler.types import sint, cint, regint, Array, ClientMessageType, sgf2n
from Compiler.compilerLib import Compiler

# we assume these modules reside in Programs/Source/ 
from shamir import shamir_share
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
    PORT_BASE = 15000
    listen_for_clients(PORT_BASE)
    socket = accept_client_connection(PORT_BASE)

    # generate and secret share key
    key = [
            apply_field_embedding(
                sgf2n.bit_compose([sgf2n.get_random_bit() for _ in range(8)])
            ) 
            for _ in range(key_len//8)
          ]

    # key = sint.get_random_int(bits=key_len)

    # key = sint(val=1) # debugging
    eval_list = [
            apply_field_embedding(sgf2n(i)) for i in range(n)]

    key_eval_points = Array(n,sgf2n).assign(eval_list)

    rand_list = [
        apply_field_embedding(sgf2n.bit_compose([sgf2n.get_random_bit() for _ in range(8)])) for i in range(t)]
    
    key_rand = Array(t,sgf2n).assign(rand_list)

    all_poly_evals = []
    for i in range(key_len//8):
        poly_evals = shamir_share(msg=key[i], threshold=t, num_parties=n, eval_points=key_eval_points, rand=key_rand)[1]
        for j in range(len(poly_evals)):
            poly_evals[j].update(apply_inverse_field_embedding(poly_evals[j]))
        all_poly_evals.append(poly_evals)
    # write Shamir shares of key back to client. 
    for i in range(n):
        @if_(regint(i) == socket) # think this is equiv to @if_(i == regint(get_player_id()._v)) 
        def _():
            for j in range(key_len//8):
                poly_eval_personal = all_poly_evals[j][i].reveal_to(i)
                cint.write_to_socket(socket, cint(poly_eval_personal._v))

    # print_ln("KEY=%s\n", key.reveal()) # debugging

if __name__ == "__main__":
    compiler.compile_func()
