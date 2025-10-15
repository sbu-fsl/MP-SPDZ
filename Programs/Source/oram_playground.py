#!/usr/bin/env python3

import os, sys
# add MP-SPDZ dir to path so we can import from Compiler
sys.path.insert(0, os.path.dirname(sys.argv[0]) + '/../..') 

from Compiler.library import print_ln, listen_for_clients, accept_client_connection, for_range, get_player_id, if_
from Compiler.types import sint, cint, regint, Array, Matrix, ClientMessageType, sgf2n, _secret
from Compiler.compilerLib import Compiler
from Compiler.oram import OptimalORAM, AbstractORAM

from embeddings import apply_field_embedding, apply_inverse_field_embedding

usage = "usage: %prog [options] [args]"
compiler = Compiler(usage=usage)
compiler.parse_args()

@compiler.register_function("test_oram")
def test_oram():
    num_entries = 3
    vals_per_entry = 1 # TODO: how to get this working with vals_per_entry > 1?
    oram = OptimalORAM(num_entries, sgf2n, vals_per_entry)
    idx_0 = sgf2n(0)
    idx_1 = sgf2n(1)
    idx_2 = sgf2n(2)
    a = sgf2n(2)
    b = sgf2n(4)
    c = sgf2n(8)
    write = sgf2n(True)
    # TODO: do we have to use access w/ secret write flag, or can we just use read/write? I think the former. 
    oram.access(idx_0, a, write)
    oram.access(idx_1, b, write)
    oram.access(idx_2, c, write)
    oram_dbg = [oram[i].reveal() for i in [sgf2n(0), sgf2n(1), sgf2n(2)]]
    print_ln("oram_dbg=%s", oram_dbg)

    # so it looks like __getitem__ and __setitem__ can be used. Does not hide the type of operation, but does hide the index. 
    # r = OptimalORAM(5, sgf2n)
    # r[0] = sgf2n(0)
    # r[1] = sgf2n(0)
    # r[2] = sgf2n(2)
    # r[3] = sgf2n(3)
    # r[4] = sgf2n(5)


    def find_nonzero_secret_idx(arr: AbstractORAM) -> _secret:
        t: _secret = arr.index_type
        num_entries = arr.size
        res = t(0)
        for i in range(num_entries):
            b = (arr[t(i)] != t(0))
            res = b.cond_swap(res, t(i))[0]
        return res

    def get_random_sgf2n(bit_length: int) -> sgf2n:
        return sgf2n.bit_compose([sgf2n.get_random_bit() for _ in range(bit_length)])
    
    x = OptimalORAM(5, sgf2n)
    r = OptimalORAM(5, sgf2n)
    y = sgf2n(1)
    for i in range(5):
        x[i] = apply_field_embedding(get_random_sgf2n(8))
        r[i] = apply_field_embedding(get_random_sgf2n(8))
    j = find_nonzero_secret_idx(r)
    x_j, r_j = r[j], x[j]
    x[j] = (y - (sum(x[i] * r[i] for i in range(5)) - (x_j * r_j))).field_div(r_j)
    x = [apply_inverse_field_embedding(x[i]).reveal() for i in range(5)]
    print_ln("x=%s", x)


if __name__ == "__main__":
    compiler.compile_func()