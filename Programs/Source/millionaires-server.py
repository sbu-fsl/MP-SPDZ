'''
NOTE: I wrote this example trying to stick to the "one-to-one" approach of personal-client-example, 
where each client is connected to exactly one MPC party that it trusts. This is why
we use personal.read_int_from_socket and cint.write_to_socket. 
After giving it some thought though, I am not sure it matters too much whether we 
do the one-to-one approach vs. the all-to-all approach of bankers_bonus. The reason it likely does
not matter is because the private input is eventually going to be secret shared across all parties either way. 
The all-to-all approach does the secret sharing right away when client calls send_private_input and server
calls sint.receive_from_client, while the one-to-one approach procrastinates the secret sharing until 
the server casts the personal values to sint. 
'''
import os, sys
# add MP-SPDZ dir to path so we can import from Compiler
sys.path.insert(0, os.path.dirname(sys.argv[0]) + '/../..') 

from Compiler.compilerLib import Compiler
from Compiler.library import listen_for_clients, accept_client_connection, for_range
from Compiler.types import sint, personal, Array, cint


usage = "usage: %prog" # not taking any arguments for now
compiler = Compiler(usage=usage)


@compiler.register_function('millionaires')
def main():
    # listen+accept client connections
    PORT_BASE = 15000
    listen_for_clients(PORT_BASE)
    socket = accept_client_connection(PORT_BASE)

    # read input from clients
    a = personal.read_int_from_socket(0, socket, 1) # player_id, socket, size=1
    b = personal.read_int_from_socket(1, socket, 1)
    c = personal.read_int_from_socket(2, socket, 1)

    client_values = Array(3, sint)
    client_values.assign([sint(a), sint(b), sint(c)])

    max_value = Array(1, sint)
    max_value[0] = client_values[0]
    winner_id = Array(1, sint)
    winner_id[0] = sint(0)

    @for_range(2)
    def _(i):
        is_new_max = max_value[0] < client_values[i+1]
        max_value[0] = is_new_max.if_else(client_values[i+1], max_value[0])
        winner_id[0] = is_new_max.if_else(sint(i+1), winner_id[0])
    res = winner_id[0].reveal() # res type is cint

    # write back to clients. 
    cint.write_to_socket(socket, res)
    cint.write_to_socket(socket, res)
    cint.write_to_socket(socket, res)

if __name__ == "__main__":
    compiler.compile_func()