#!/usr/bin/env python3

import os, sys
from copy import copy
# add MP-SPDZ dir to path so we can import from Compiler. It is assumed this file lives in MP-SPDZ/Programs/Source. 
sys.path.insert(0, os.path.dirname(sys.argv[0]) + '/../..') 
from Compiler.library import print_ln, vectorize
from Compiler.types import cgf2n, sgf2n, Array, Matrix, VectorArray
from Compiler.compilerLib import Compiler

from embeddings import *

# CONSTANTS
BYTES_PER_WORD = 4
WORDS_PER_BLOCK = 4 # Nb = 4 words = 16 bytes = 128 bits


class SBox():
    '''
    Implementation of the S-Box, and its inverse, as in FIPS 197. 
    '''

    def __init__(self):
        self.EI = EmbeddedInverter()
        self.matrix = [
            [1,0,0,0,1,1,1,1],
            [1,1,0,0,0,1,1,1],
            [1,1,1,0,0,0,1,1],
            [1,1,1,1,0,0,0,1],
            [1,1,1,1,1,0,0,0],
            [0,1,1,1,1,1,0,0],
            [0,0,1,1,1,1,1,0],
            [0,0,0,1,1,1,1,1]
        ]
        self.matrix_inv = [
            [0,0,1,0,0,1,0,1],
            [1,0,0,1,0,0,1,0],
            [0,1,0,0,1,0,0,1],
            [1,0,1,0,0,1,0,0],
            [0,1,0,1,0,0,1,0],
            [0,0,1,0,1,0,0,1],
            [1,0,0,1,0,1,0,0],
            [0,1,0,0,1,0,1,0]
        ]
        self.affine_constant = [1,1,0,0,0,1,1,0]

    def apply(self, byte: cgf2n | sgf2n) -> cgf2n | sgf2n:
        '''
        Applies the S-Box to an embedded byte

        :param byte: cgf2n | sgf2n, assumed to be an embedded byte. 
        '''
        byte_inv_bd = self.EI.invert(byte).bit_decompose(bit_length=40, step=5)
        b_tilde = apply_inverse_field_embedding_bd(byte_inv_bd) # named b_tilde to match FIPS 197 notation
        affine_transform = [ sum(b_tilde[idx] * row[idx] for idx in range(8)) + c for row, c in zip(self.matrix, self.affine_constant) ]

        # affine_transform holds our result in unembedded, bit_decomposed form (list of length 8)
        # so we basically just need the last two lines of apply_field_embedding
        out_bytes = apply_field_embedding_bd(affine_transform)
        return type(byte)(sum(out_bytes[idx] * (cgf2n(2) ** (5*idx)) for idx in range(8)))


class AESCipher():
    '''
    Implementation of the AES cipher, per FIPS 197.
    '''
    def __init__(self, key: list[sgf2n], sbox=None):
        '''
        Create an instance of AESCipher. This involves setting up basic parameters, an SBox, and 
        performing key expansion. 

        :param key: list[sgf2n]. An unembedded key of length 16, 24, or 32 bytes (128, 192, or 256 bits).
        :param sbox: SBox, optional. If not provided, a default SBox will be used. 

        TODO: add nparallel as a parameter, and pass it as an arg every time we instantiate a runtime data type (e.g., cgf2n, sgf2n, VectorArray)
        '''
        assert(len(key) in [16, 24, 32]) # key must be 16, 24, or 32 bytes long
        num_rounds_from_key_length = {16: 10, 24: 12, 32: 14}
        self.num_rounds = num_rounds_from_key_length[len(key)] # set num_rounds based on key length

        self.key_length = len(key) // BYTES_PER_WORD # Nk = length of key in words
        self.sbox = sbox if sbox else SBox()
        key = [apply_field_embedding(x) for x in key] # embed the key, so it is in GF(2^40)
        self.key_schedule = self.key_expansion(key) # 4*(Nr+1) words = 16*(Nr+1) bytes

    def key_expansion(self, key: list[sgf2n]) -> list[sgf2n]:
        '''
        Perform FIPS 197 key expansion on an embedded key.

        :param key: list[sgf2n]. Assumed to be an embedded key of length self.key_length words
        :return: list[sgf2n]. A key schedule consisting of (self.num_rounds + 1) embedded round keys.
        '''
        # initialize round constants
        rcon_raw = [0x00, 0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, 0x80, 0x1b, 0x36] # 0th byte should never be used, it's just that rcon is 1-indexed in FIPS 197
        # TODO: why does VectorArray keep cropping up? Can't we just do rcon = [cgf2n(x) for x in rcon_raw]?
        rcon = VectorArray(len(rcon_raw), cgf2n, 1)
        for idx in range(len(rcon_raw)):
            rcon[idx] = cgf2n(rcon_raw[idx])
        
        # number of round keys is (Nr + 1) blocks
        key_schedule = [sgf2n(0)] * ((WORDS_PER_BLOCK * BYTES_PER_WORD) * (self.num_rounds + 1))
        temp = [sgf2n(0)] * BYTES_PER_WORD
        
        # first Nk words of key schedule are just the key
        for i in range(self.key_length):
            for j in range(BYTES_PER_WORD):
                key_schedule[i * BYTES_PER_WORD + j] = key[i * BYTES_PER_WORD + j]
        
        for i in range(self.key_length, WORDS_PER_BLOCK * (self.num_rounds + 1)): 
            # init temp to (i-1)-th word of key schedule
            for j in range(BYTES_PER_WORD): # TODO: can we use slicing to clean this up?
                temp[j] = key_schedule[(i-1) * BYTES_PER_WORD + j]

            # print_ln("i =%s, temp init = %s", i, [apply_inverse_field_embedding(x.reveal()) for x in temp]) # temp is initialized just fine

            # do stuff to temp based on Nk
            if i % self.key_length == 0:
                temp = self.sub_word(self.rot_word(temp)) 
                temp[0] += apply_field_embedding(rcon[i // self.key_length])
            elif self.key_length > 6 and i % self.key_length == 4:
                temp = self.sub_word(temp)
            # next word of key_schedule depends on key_schedule Nk words back, and whatever temp is.
            for j in range(BYTES_PER_WORD):
                key_schedule[i * BYTES_PER_WORD + j] = key_schedule[(i - self.key_length) * BYTES_PER_WORD + j] + temp[j]
        return key_schedule
    
    def cipher(self, input: list[sgf2n]) -> list[sgf2n | cgf2n]:
        '''
        Apply the AES block cipher to a 128-bit input.

        :param input: Assumed to be a list[sgf2n] of length 16, where each element holds an unembedded byte in its lower 8 bits.
        :return: Resulting cipher output of length 16, where each element holds an unembedded byte in its lower 8 bits.
        '''
        state = [apply_field_embedding(x) for x in input] # embed input and copy to state vector
        round_key = self.key_schedule[0 : (WORDS_PER_BLOCK * BYTES_PER_WORD)] # each round key is 4 words of key schedule
        self.add_round_key(state, round_key)
        for round in range(1, self.num_rounds):
            self.sub_bytes(state)
            self.shift_rows(state)
            self.mix_columns(state)
            round_key = self.key_schedule[round * WORDS_PER_BLOCK * BYTES_PER_WORD : ((round+1) * WORDS_PER_BLOCK * BYTES_PER_WORD)]
            self.add_round_key(state, round_key)
        self.sub_bytes(state)
        self.shift_rows(state)
        round_key = self.key_schedule[self.num_rounds * WORDS_PER_BLOCK * BYTES_PER_WORD : ]
        self.add_round_key(state, round_key)
        return [apply_inverse_field_embedding(x) for x in state]
    
    def cipher_inverse(self, ciphertext):
        '''
        '''
        pass

    def sub_word(self, word: list[sgf2n]) -> list[sgf2n]:
        '''
        Applies the S-Box to each byte in word.

        :param word: list[sgf2n]. Assumed to hold 1 word = 4 bytes.
        '''
        return [self.sbox.apply(word[i]) for i in range(BYTES_PER_WORD)]
        pass

    def rot_word(self, word: list[sgf2n]) -> list[sgf2n]:
        '''
        Rotates word one byte to the left: [a_0, a_1, a_2, a_3] -> [a_1, a_2, a_3, a_0]

        :param word: list[sgf2n]. Assumed to hold 1 word = 4 bytes
        '''
        return [word[1], word[2], word[3], word[0]]
    
    def rot_word_left(self, word: list[sgf2n], offset: int) -> list[sgf2n]:
        '''
        Rotates word offset bytes to the left. 
        For example: [a_0, a_1, a_2, a_3] -> [a_1, a_2, a_3, a_0] for offset 1. 

        :param word: list[sgf2n]. Assumed to hold 1 word = 4 bytes
        :param offset: int. We assume 0 <= offset < 4.
        '''
        # return [word[offset % 4], word[(offset + 1) % 4], word[(offset + 2) % 4], word[(offset + 3) % 4]]
        return word[offset:] + word[:offset]

    def sub_bytes(self, state: list[sgf2n]):
        '''
        Apply S-Box to every element (an embedded sgf2n byte) of state in-place

        :param state: list[sgf2n]. We assume elements are embedded. Modified in-place. 
        '''
        for i in range(len(state)):
            state[i] = self.sbox.apply(state[i])

    def shift_rows(self, state: list[sgf2n]):
        '''
        Shifts (i.e. rotated) row of the state. The zeroth row is unchanged. The first row is rotated 
        one to the left. The second row is rotated two to the left. The third is rotated three to the 
        left. The state is modified in-place.

        :param state: list[sgf2n]. We assume elements are embedded, and 4x4 state matrix is stored in column-major order, as specified by FIPS 197. Modified in-place.
        '''
        for i in range(WORDS_PER_BLOCK): # WORDS_PER_BLOCK = 4 = length of each row, i.e. number of columns.
            state[i::WORDS_PER_BLOCK] = self.rot_word_left(state[i::WORDS_PER_BLOCK],i)

    def mix_columns(self, state: list[sgf2n]):
        '''
        Multiply state (viewed as a matrix of embedded values) on the left with 
        an embedded version of the following fixed MDS matrix:
        [
            [0x02, 0x03, 0x01, 0x01],
            [0x01, 0x02, 0x03, 0x01],
            [0x01, 0x01, 0x02, 0x03],
            [0x03, 0x01, 0x01, 0x02]
        ]

        :param state: list[sgf2n]. Elements assumed to be embedded. Modified in-place. 

        NOTE: In this implementation, we avoid matrix multiplication altogether by exploiting the 
        specific values of the fixed MDS matrix and properties of binary fields. 
        '''

        def mix_column(column: list[sgf2n]) -> list[sgf2n]:
            '''
            Helper function for computing a single column of mix_columns matrix multiplication.
            '''
            temp = copy(column)
            doubles = [apply_field_embedding(cgf2n(2)) * t for t in temp]
            column[0] = doubles[0] + (temp[1] + doubles[1]) + temp[2] + temp[3]
            column[1] = temp[0] + doubles[1] + (temp[2] + doubles[2]) + temp[3]
            column[2] = temp[0] + temp[1] + doubles[2] + (temp[3] + doubles[3])
            column[3] = (temp[0] + doubles[0]) + temp[1] + temp[2] + doubles[3]
        
        for i in range(WORDS_PER_BLOCK):
            column = []
            for j in range(BYTES_PER_WORD):
                column.append(state[i*BYTES_PER_WORD+j])
            mix_column(column)
            for j in range(BYTES_PER_WORD):
                state[i*BYTES_PER_WORD+j] = column[j]
    
    def add_round_key(self, state: list[sgf2n], round_key: list[sgf2n]):
        '''
        XOR the state with round_key. Modifies state in-place. 
        '''
        for i in range(len(state)):
            state[i] = state[i] + round_key[i]



if __name__ == "__main__":
    usage = "usage: %prog [options] [args]"
    compiler = Compiler(usage=usage)
    compiler.parse_args()

    @compiler.register_function("test_inverse")
    def test_inverse():
        EI = EmbeddedInverter()
        b = cgf2n(0x8000) # embedding of 0xf
        b_inv = EI.invert(b)
        b_inv_alt = cgf2n(1).field_div(b)
        a_inv = apply_inverse_field_embedding(b_inv)
        print_ln("b=%s, b_inv=%s, b_inv_alt=%s, a_inv=%s", b, b_inv, b_inv_alt, a_inv) # b_inv = 0x802008401 = b_inv_alt. a_inv = 0xc7

        a = cgf2n(0x8d)
        b = apply_field_embedding(a)
        b_inv = EI.invert(b)
        b_inv_alt = cgf2n(1).field_div(b)
        a_inv = apply_inverse_field_embedding(b_inv)
        print_ln("b=%s, b_inv=%s, b_inv_alt=%s, a_inv=%s", b, b_inv, b_inv_alt, a_inv) # a_inv = 0x2

    compiler.compile_func()

    @compiler.register_function("test_sbox")
    def test_sbox():
        sbox = SBox()
        a = apply_field_embedding(cgf2n(0x61))
        b = sbox.apply(a)
        b_unembed = apply_inverse_field_embedding(b)
        print_ln("a=%s, b=%s, b_unembed=%s, expected_b_unembed=%s", a, b, b_unembed, cgf2n(0xef))

    compiler.compile_func()

    # @compiler.register_function("test_mix_columns")
    # def test_mix_columns():
    #     # TODO: consider making mix_columns a static / class method. Don't want/need to compile all of AES for this test. 
    #     key_raw = "2b7e151628aed2a6abf7158809cf4f3c"
    #     key = [apply_field_embedding(sgf2n(x)) for x in str_to_hex(key_raw)]
    #     aes = AESCipher(key)

    #     aes.mix_columns(state)

    # compiler.compile_func()

    def str_to_hex(x):
        ''' Convert a string into a list of hex values. Obviously the string should represent valid hex to begin with. '''
        return [int(x[i : i + 2], 16) for i in range(0, len(x), 2)]


    @compiler.register_function("test_key_expansion")
    def test_key_expansion():
        # FIPS 197 Appendix A.1 example
        key_raw = "2b7e151628aed2a6abf7158809cf4f3c"
        key = [sgf2n(x) for x in str_to_hex(key_raw)]
        aes = AESCipher(key)
        key_schedule = [apply_inverse_field_embedding(x.reveal()) for x in aes.key_schedule]
        expected_key_schedule_raw = "2b7e151628aed2a6abf7158809cf4f3c" + "a0fafe17" + "88542cb1" + "23a33939" + "2a6c7605" + "f2c295f2" + "7a96b943" + "5935807a" + "7359f67f" + "3d80477d" + "4716fe3e" + "1e237e44" + "6d7a883b" + "ef44a541" + "a8525b7f" + "b671253b" + "db0bad00" + "d4d1c6f8" + "7c839d87" + "caf2b8bc" + "11f915bc" + "6d88a37a" + "110b3efd" + "dbf98641" + "ca0093fd" + "4e54f70e" + "5f5fc9f3" + "84a64fb2" + "4ea6dc4f" + "ead27321" + "b58dbad2" + "312bf560" + "7f8d292f" + "ac7766f3" + "19fadc21" + "28d12941" + "575c006e" + "d014f9a8" + "c9ee2589" + "e13f0cc8" + "b6630ca6"
        expected_key_schedule = [cgf2n(x) for x in str_to_hex(expected_key_schedule_raw)]
        error_pattern = [x + y for x,y in zip(key_schedule, expected_key_schedule)]
        print_ln("expanded key = %s\nexpected key = %s\nerror pattern = %s", key_schedule, expected_key_schedule, error_pattern)

    compiler.compile_func()

    @compiler.register_function("test_cipher")
    def test_cipher():
        # FIPS 197 Appendix B example
        key_raw = "2b7e151628aed2a6abf7158809cf4f3c"
        msg_raw = "3243f6a8885a308d313198a2e0370734"
        expected_ct_raw = "3925841d02dc09fbdc118597196a0b32"
        key = [sgf2n(x) for x in str_to_hex(key_raw)]
        msg = [sgf2n(x) for x in str_to_hex(msg_raw)]
        expected_ct = [cgf2n(x) for x in str_to_hex(expected_ct_raw)]
        aes = AESCipher(key)
        ct = [x.reveal() for x in aes.cipher(msg)]
        error_pattern = [x + y for x,y in zip(expected_ct, ct)]
        print_ln("EX1: ciphertext = %s\nexpected ciphertext = %s\nerror_pattern = %s", ct, expected_ct, error_pattern)

        # Original MP-SPDZ aes.mpc example
        msg_raw = "6bc1bee22e409f96e93d7e117393172a"
        expected_ct_raw = "3ad77bb40d7a3660a89ecaf32466ef97"
        msg = [sgf2n(x) for x in str_to_hex(msg_raw)]
        expected_ct = [cgf2n(x) for x in str_to_hex(expected_ct_raw)]
        ct = [x.reveal() for x in aes.cipher(msg)]
        error_pattern = [x + y for x,y in zip(expected_ct, ct)]
        print_ln("EX2: ciphertext = %s\nexpected ciphertext = %s\nerror_pattern = %s", ct, expected_ct, error_pattern)

        # TODO: why doesn't this work? Get list index out of range in mix_columns: column.append(state[i*4+j])
        # would be nice if we could get this to work - should give shorter compile times. 

        # key_raw = "2b7e151628aed2a6abf7158809cf4f3c"
        # key = [apply_field_embedding(sgf2n(x)) for x in str_to_hex(key_raw)]
        # aes = AESCipher(key)
        # messages_hex: list[list[int]] = [str_to_hex(x) for x in ["3243f6a8885a308d313198a2e0370734", "6bc1bee22e409f96e93d7e117393172a"]]
        # expected_ciphertexts_raw: list[list[int]] = [str_to_hex(x) for x in ["3925841d02dc09fbdc118597196a0b32", "3ad77bb40d7a3660a89ecaf32466ef97"]]
        # msgs = [apply_field_embedding(sgf2n(list(x))) for x in zip(messages_hex)]
        # expected_cts = [cgf2n(list(x)) for x in zip(expected_ciphertexts_raw)]
        # cts = [apply_inverse_field_embedding(x.reveal()) for x in aes.cipher(msgs)]
        # error_patterns = [x + y for x,y in zip(expected_cts, cts)]
        # print_ln("ciphertexts = %s\nexpected ciphertexts = %s\nerror_patterns = %s", cts, expected_cts, error_patterns)

    compiler.compile_func()





