
from aes import AESCipher, WORDS_PER_BLOCK, BYTES_PER_WORD
from Compiler.library import print_ln, if_e, else_
from Compiler.types import cgf2n, sgf2n, regint
from Compiler.compilerLib import Compiler
from math import ceil
from copy import copy

# assume these modules live in Programs/Source/
from utils import str_to_hex

BLOCK_SIZE = WORDS_PER_BLOCK * BYTES_PER_WORD



def aes_cmac(key: list[sgf2n], m: list[sgf2n | cgf2n], tlen: int) -> list[sgf2n | cgf2n]: 
    '''
    CMAC(K,M,Tlen) as as described in NIST SP 800-38B (with AES cipher).
    For ease of implementation, we enforce that the message and tag are byte-aligned. 
    Additionally, we enforce that tlen <= BLOCK_SIZE.

    :param key: MAC key represented as unembedded list[sgf2n]
    :param m: the message to be authenticated
    :param tlen: compile-time int representing the length of the outputted tag in bytes.
    '''
    assert(tlen <= BLOCK_SIZE)
    
    # define subkey generation function using aes as the cipher. 
    aes = AESCipher(key)
    def subkey() -> tuple[list[sgf2n], list[sgf2n]]:
        '''
        SUBK(k) function as described in NIST SP 800-38B for use in CMAC.

        NOTE: As implemented here subkey does not take in a key parameter, 
        since we construct the aes instance outside of the function scope. 
        This is only to save on multiple constructions of aes. Ideally there's 
        a better way without having to construct aes more than once...
        '''
        # degree 128 irreducible poly defined by 0^{120} || 10000111 as in CMAC standard
        R = [sgf2n(0) for _ in range(15)] + [sgf2n(0x87)] 
        zero_block = [sgf2n(0) for _ in range(BLOCK_SIZE)] # no need to embed 0
        L = aes.cipher(zero_block)

        # compute k_1 as (L<<1) if msb(L) == 0, else (L<<1) XOR R
        L = [byte.bit_decompose(8) for byte in L]
        L_msb = L[0][-1]
        L_shifted_lsb = [L[i+1][-1] for i in range(BLOCK_SIZE - 1)] + [sgf2n(0)] # lsb of L_shifted[i] is msb of L[i+1]
        L_shifted = [[lsb] + byte[:-1] for lsb, byte in zip(L_shifted_lsb,L)] # rotate each decomposed byte one bit to the RIGHT (since LSB first in each word).
        L_shifted = [sgf2n.bit_compose(byte) for byte in L_shifted] # L_shifted := L << 1
        k_1 = [byte + (L_msb * r) for byte, r in zip(L_shifted, R)]

        # compute k_2 as (k_1<<1) if msb(k_1) == 0, else (k_1<<1) XOR R
        k_1_dec = [byte.bit_decompose(bit_length=8) for byte in k_1]
        k_1_msb = k_1_dec[0][-1]
        k_1_shifted_lsb = [k_1_dec[i+1][-1] for i in range(BLOCK_SIZE - 1)] + [sgf2n(0)]
        k_1_shifted_dec = [[k_1_shifted_lsb[i]] + byte[:-1] for i, byte in enumerate(k_1_dec)]
        k_1_shifted = [sgf2n.bit_compose(byte) for byte in k_1_shifted_dec]
        k_2 = [byte + (k_1_msb * r) for byte, r in zip(k_1_shifted, R)]

        return k_1, k_2

    # generate subkeys
    k_1, k_2 = subkey() 

    # set last block of m; pad if necessary
    n = ceil(len(m) / (BLOCK_SIZE)) if len(m) != 0 else 1 # number of blocks in m
    print(f"len(m)={len(m)}, n={n}")
    m = copy(m) # avoid mutating argument
    last_block = m[(n-1)*BLOCK_SIZE:]
    if len(last_block) == BLOCK_SIZE:
        last_block = [k_1[i] + last_block[i] for i in range(BLOCK_SIZE)]
    else: # need to pad!
        padding = [sgf2n(0)] * (BLOCK_SIZE - len(last_block))
        first_padding_byte = padding[0].bit_decompose(8)
        first_padding_byte[-1] = sgf2n(1)
        padding[0] = sgf2n.bit_compose(first_padding_byte)
        last_block = last_block + padding
        last_block = [k_2[i] + last_block[i] for i in range(BLOCK_SIZE)]
    m[(n-1)*BLOCK_SIZE:] = last_block
    assert(len(m) == n * BLOCK_SIZE)

    # cipher block chaining
    c = [sgf2n(0)] * BLOCK_SIZE
    for i in range(n):
        block = m[i*BLOCK_SIZE : (i+1)*BLOCK_SIZE]
        c = aes.cipher([c[i] + block[i] for i in range(BLOCK_SIZE)])

    # extract tag from tlen most significant bits of c
    tag = c[:tlen]
    return tag



if __name__ == "__main__":
    usage = "usage: %prog [options] [args]"
    compiler = Compiler(usage=usage)

    @compiler.register_function("test_cmac")
    def test_cmac():
        # see https://csrc.nist.gov/CSRC/media/Projects/Cryptographic-Standards-and-Guidelines/documents/examples/AES_CMAC.pdf
        # for cmac correctness tests

        print_ln("TEST aes_cmac")
        # test 1: empty message
        # expected k_2 = F7DDAC30 6AE266CC F90BC11E E46D513B
        key = [sgf2n(byte) for byte in str_to_hex("2B7E151628AED2A6ABF7158809CF4F3C")]
        msg = []
        tag = [t.reveal() for t in aes_cmac(key, msg, BLOCK_SIZE)]
        expected_tag = [cgf2n(byte) for byte in str_to_hex("BB1D6929E95937287FA37D129B756746")]
        error_pattern = [x - y for x,y in zip(tag, expected_tag)]
        @if_e(sum(error_pattern) != cgf2n(0))
        def _():
            print_ln("❌ TEST 1 FAILED\ntag=%s\nexpected tag=%s", tag, expected_tag)
        @else_
        def _():
            print_ln("✅ TEST 1 PASSED")
        
        # test 2: m = 6BC1BEE2 2E409F96 E93D7E11 7393172A
        key = [sgf2n(byte) for byte in str_to_hex("2B7E151628AED2A6ABF7158809CF4F3C")]
        msg = [sgf2n(byte) for byte in str_to_hex("6BC1BEE22E409F96E93D7E117393172A")]
        tag = [t.reveal() for t in aes_cmac(key, msg, BLOCK_SIZE)]
        expected_tag = [cgf2n(byte) for byte in str_to_hex("070A16B46B4D4144F79BDD9DD04A287C")]
        error_pattern = [x - y for x,y in zip(tag, expected_tag)]
        @if_e(sum(error_pattern) != cgf2n(0))
        def _():
            print_ln("❌ TEST 2 FAILED\ntag=%s\nexpected tag=%s", tag, expected_tag)
        @else_
        def _():
            print_ln("✅ TEST 2 PASSED")

    compiler.compile_func()
