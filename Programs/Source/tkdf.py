#!/usr/bin/env python3

import os, sys
# add MP-SPDZ dir to path so we can import from Compiler
sys.path.insert(0, os.path.dirname(sys.argv[0]) + '/../..') 

from Compiler.library import print_ln, listen_for_clients, accept_client_connection, for_range, get_player_id, if_, public_input
from Compiler.types import regint, sgf2n, cgf2n, cint
from Compiler.compilerLib import Compiler

# we assume these modules reside in Programs/Source/ 
from new_shamir import shamir_reconstruct, shamir_share
from embeddings import apply_field_embedding, apply_inverse_field_embedding
from kdf_ctr import kdf_ctr, BLOCK_SIZE
from utils import get_random_sgf2n

usage = "usage: %prog [options] [args]"
compiler = Compiler(usage=usage)
compiler.parser.add_option("--root-key-length", dest="root_key_len")
compiler.parser.add_option("--child-key-length", dest="child_key_len")
compiler.parser.add_option("--context-length", dest="ctx_len")
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

@compiler.register_function('tkdf')
def tkdf():
    # constants
    opts = compiler.options
    [root_key_len, child_key_len, ctx_len, t, n] = [int(x) for x in [opts.root_key_len, opts.child_key_len, opts.ctx_len, opts.t, opts.n]] 
    assert(root_key_len % 8 == 0)
    assert(child_key_len % 8 == 0)
    assert(ctx_len % 8 == 0)
    num_bytes_root_key = root_key_len // 8
    num_bytes_child_key = child_key_len // 8
    ctx_len = ctx_len // 8

    # set up external client connections
    PORT_BASE = public_input()
    listen_for_clients(PORT_BASE)
    socket = accept_client_connection(PORT_BASE)

    # context should be written to Programs/Public-Input/<progname> in decimal integer format (public_input returns cint)
    # each integer should be in range [0,255]. There should be exactly ctx_len integers. 
    # we cast to sgf2n not because context is secret info, but because kdf_ctr expects a list[sgf2n]
    context = [sgf2n(public_input()) for _ in range(ctx_len)]

    # reconstruct root key
    input_shares = [sgf2n.get_input_from(i, size=num_bytes_root_key) for i in range(n)] # read from Player-Data/Input-P<player>-<thread> in HEX FORMAT
    input_shares_emb = [apply_field_embedding(share) for share in input_shares]
    eval_points_emb = [apply_field_embedding(sgf2n(i)) for i in range(1,n+1)]
    root_key_emb: sgf2n = shamir_reconstruct(
        input_shares_emb, 
        eval_points=eval_points_emb, 
        size=num_bytes_root_key
    )
    root_key: sgf2n = apply_inverse_field_embedding(root_key_emb)
    root_key: list[sgf2n] = [byte for byte in root_key] # "unvectorize" root_key into list[sgf2n] for use with kdf_ctr
    
    # derive new key using root_key as kdk and label+context from public_input
    # hardcode label as sgf2n(0x01) for now, as we only have one purpose for child keys (use as a DEK)
    child_key: list[sgf2n] = kdf_ctr(root_key, h=BLOCK_SIZE, r=1, L=num_bytes_child_key, label=[sgf2n(0x01)], context=context)

    # secret share and reveal each share to appropriate party
    child_key_shares_by_party = {party: None for party in range(n)}
    child_key_emb = apply_field_embedding(sgf2n(child_key, size=num_bytes_child_key)) # vectorize and embed before secret sharing
    rand_emb = [apply_field_embedding(get_random_sgf2n(8, size=num_bytes_child_key)) for _ in range(t)]
    _, child_key_shares_emb = shamir_share(
        msg=child_key_emb, 
        threshold=t, 
        num_parties=n, 
        eval_points=eval_points_emb, 
        rand=rand_emb,
        size=num_bytes_child_key
    )
    for party, child_key_share_emb in enumerate(child_key_shares_emb):
        child_key_share = apply_inverse_field_embedding(child_key_share_emb)
        child_key_shares_by_party[party] = child_key_share.reveal_to(party)
    
    # write back child key shares to corresponding parties
    for party in range(n):
        @if_(party == socket)
        def _():
            child_key_share = child_key_shares_by_party[party]
            # cint.write_to_socket(socket, cint(child_key_share._v))
            for byte_share in child_key_share:
                cint.write_to_socket(socket, cint(byte_share._v))


if __name__ == "__main__":
    compiler.compile_func()