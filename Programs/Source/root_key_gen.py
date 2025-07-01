#!/usr/bin/env python3

import os, sys
# add MP-SPDZ dir to path so we can import from Compiler
sys.path.insert(0, os.path.dirname(sys.argv[0]) + '/../..') 

from Compiler.library import print_ln, listen_for_clients, accept_client_connection, for_range, get_player_id, if_
from Compiler.types import sint, cint, regint, Array, ClientMessageType
from Compiler.compilerLib import Compiler

# we assume these modules reside in Programs/Source/ 
from shamir import shamir_share

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
    key = sint.get_random_int(bits=key_len)
    # key = sint(val=1) # debugging
    eval_points, poly_evals = shamir_share(msg=key, threshold=t, num_parties=n)
    eval_points = eval_points.reveal() # nothing secret about these as all. Should consider changing impl of shamir to make them public. 

    # write Shamir shares of key back to client. 
    @for_range(regint(n))
    def _(i):
        poly_eval_personal = poly_evals[i].reveal_to(i)
        @if_(i == socket) # think this is equiv to @if_(i == regint(get_player_id()._v)) 
        def _():
            cint.write_to_socket(socket, [eval_points[i], poly_eval_personal._v])

    # print_ln("KEY=%s\n", key.reveal()) # debugging

if __name__ == "__main__":
    compiler.compile_func()
