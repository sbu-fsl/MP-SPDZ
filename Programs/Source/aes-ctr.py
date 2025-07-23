from aes import AESCipher, apply_field_embedding, apply_inverse_field_embedding, BLOCK_SIZE, BYTES_PER_WORD
from Compiler.library import print_ln
from Compiler.types import cgf2n, sgf2n, regint
from Compiler.compilerLib import Compiler
from Compiler.GC.types import sbits # requires binary circuits. Have to compile with -X or -Y option for daBit / edaBit

MAX_BLOCKS = 2 ** 32  # maximum number of plaintext blocks we can handle in AES-CTR mode

class AES_CTR():
    def __init__(self, key: list[sgf2n], num_blocks: int):
        '''
        Initialize nonce and counters for AES-CTR mode for num_blocks of plaintext.
        1 block = 16 bytes = 128 bits.

        :param key: 16-byte key for AES (unembedded)
        :param num_blocks: number of blocks of plaintext to encrypt/decrypt
        '''
        num_rounds_from_key_length = {16: 10, 24: 12, 32: 14}
        assert(len(key) in num_rounds_from_key_length)  # key must be 16, 24, or 32 bytes long
        self.aes = AESCipher(num_rounds_from_key_length[len(key)], key)

        assert(num_blocks <= MAX_BLOCKS)
        self.num_blocks = num_blocks
        # nonce should be 96 bits = 12 bytes, and counter should be 32 bits = 4 bytes.
        self.nonce = [cgf2n(regint.get_random(bit_length=8)) for _ in range(12)]
        # TODO: for parallelism, maybe we can wrap counters in a single cgf2n?
        self.counters = [self.nonce + [cgf2n(x) for x in i.to_bytes(length=4)] for i in range(num_blocks)]

    def encrypt(self, plaintext: list[sgf2n]) -> tuple[list[sgf2n], list[sgf2n]]:
        '''
        Encrypt plaintext using AES-CTR mode.

        :param plaintext: plaintext to encrypt. Length of plaintext should be self.num_blocks * BLOCK_SIZE * BYTES_PER_WORD, since we don't yet support padding.
        :return: list of sgf2n values representing the ciphertext
        '''
        assert(len(plaintext) == self.num_blocks * BLOCK_SIZE * BYTES_PER_WORD) # Too lazy to deal with padding.
        # chunk plaintext into blocks and copy to ciphertext
        plaintext = [plaintext[i : i + BLOCK_SIZE * BYTES_PER_WORD] for i in range(self.num_blocks)]
        ciphertext = plaintext
        for i in range(self.num_blocks):
            cipher_output = self.aes.cipher(self.counters[i])
            for j in range(BLOCK_SIZE * BYTES_PER_WORD):
                # is this technically incorrect bc not indexing in column major order?
                ciphertext[i][j] = plaintext[i][j] + cipher_output[j]
        # return nonce and flattened ciphertext
        return self.nonce, [c for block in ciphertext for c in block]
    
    def decrypt(self, ciphertext: list[sgf2n], nonce) -> list[sgf2n]:
        '''
        '''
        pass


if __name__ == "__main__":
    usage = "usage: %prog [options] [args]"
    compiler = Compiler(usage=usage)
    compiler.parse_args()

    def str_to_hex(x):
        ''' Convert a string into a list of hex values. Obviously the string should represent valid hex to begin with. '''
        return [int(x[i : i + 2], 16) for i in range(0, len(x), 2)]

    @compiler.register_function("test_aes_ctr")
    def test_aes_ctr():
        key_raw = "2b7e151628aed2a6abf7158809cf4f3c"
        key = [sgf2n(x) for x in str_to_hex(key_raw)]
        aes_ctr = AES_CTR(key, num_blocks=4)
        msg_raw = "68656C6C6F2074686572652C207468697320697320612074657374206D65737361676520666F7220656E6372797074696F6E20707572706F7365732E2E2E2E2E"
        msg = [sgf2n(x) for x in str_to_hex(msg_raw)]
        nonce, ct = aes_ctr.encrypt(msg)
        ct = [c.reveal() for c in ct]
        print_ln("nonce=%s\n ct=%s", nonce, ct)


        
        

        


    compiler.compile_func()