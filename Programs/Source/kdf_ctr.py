#!/usr/bin/env python3

import os, sys
# add MP-SPDZ dir to path so we can import from Compiler
sys.path.insert(0, os.path.dirname(sys.argv[0]) + '/../..') 

from Compiler.library import print_ln, listen_for_clients, accept_client_connection, if_, public_input
from Compiler.types import cint, Array, sgf2n, cgf2n
from Compiler.compilerLib import Compiler
from math import ceil

# we assume these modules reside in Programs/Source/ 
from embeddings import apply_field_embedding, apply_inverse_field_embedding
from cmac import aes_cmac, BLOCK_SIZE
from utils import int_to_sgf2n_bytes, str_to_hex

def kdf_ctr(kdk: list[sgf2n], h: int, r: int,  L: int, label: list[sgf2n], context: list[sgf2n]) -> list[sgf2n]:
    '''
    KDF in counter mode as described in NIST SP 800-108r1-upd1. Note the specification calls for h,r,L to 
    represent a number of bits, but since everything we do is byte-aligned, we have h,r,L represent lengths in 
    bytes. Length of (context, label, L) as a single list[sgf2n] must not exceed 

    :param kdk: key derivation key as unembedded list[sgf2n] 
    :param h: compile-time int representing length of output of a single PRF invocation in bytes.
    :param r: compile-time int representing length of a binary counter in bytes. 
    :param L: compile-time int representing length of the derived key in bytes
    :param label: unembedded list[sgf2n] encoding the purpose for the derived keying material. Caller is responsible for the encoding method. 
    :param context: unembedded list[sgf2n] encoding information related to derive keying material (e.g., entity identifiers, session IDs)
    '''
    n = ceil(L / h) 
    assert(n < (2**(r*8)) - 1) # make sure n doesn't exceed max possible value of counter. 

    sep = sgf2n(0x00) # separator between label and context
    L_bytes = int_to_sgf2n_bytes(L, ceil(L.bit_length() / 8))
    res = []
    for i in range(n):
        ctr = int_to_sgf2n_bytes(i, r)
        cmac_input = ctr + label + [sep] + context + L_bytes
        res += aes_cmac(kdk, cmac_input, h)
    # get L leftmost bytes of result
    return res[:L]

if __name__ == "__main__":
    usage = "usage: %prog [options] [args]"
    compiler = Compiler(usage=usage)

    @compiler.register_function("test_kdf")
    def test_kdf():
        print_ln("TEST kdf_ctr")
        # test 1: easiest case
        kdk = [sgf2n(byte) for byte in str_to_hex("2B7E151628AED2A6ABF7158809CF4F3C")]
        h = BLOCK_SIZE
        r = 1
        L = BLOCK_SIZE
        label = [sgf2n(0x01)] # e.g., enum specifying DEK=0x01
        context = [sgf2n(byte) for byte in str_to_hex("897F5978B0F011F0A0D9D6F009A4B5AF")] # uuid.uuid1().hex.upper()
        dek = kdf_ctr(kdk, h, r, L, label, context)
        print_ln("dek=%s", [x.reveal() for x in dek])
        # TODO: need test cases
    compiler.compile_func()
