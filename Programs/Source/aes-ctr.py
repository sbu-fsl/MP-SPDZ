from aes import AESCipher, apply_field_embedding, apply_inverse_field_embedding, BLOCK_SIZE, BYTES_PER_WORD
from Compiler.library import print_ln
from Compiler.types import cgf2n, sgf2n, regint
from Compiler.compilerLib import Compiler

MAX_BLOCKS = 2 ** 32  # maximum number of plaintext blocks we can handle in AES-CTR mode


def aes_ctr_encrypt(key: list[sgf2n], plaintext: list[sgf2n], nonce: list[cgf2n] = None) -> tuple[list[sgf2n], list[sgf2n]]:
    '''
    Encrypt plaintext using AES-CTR mode.

    :param plaintext: plaintext to encrypt. Length of plaintext should be multiple of BLOCK_SIZE * BYTES_PER_WORD, since we don't yet support padding.
    :return: A tuple with the nonce as the first coordinate, and the ciphertext as the second coordinate. 
    '''
    # validate plaintext length and set up nonce + counters based on length of plaintext
    assert(len(plaintext) % (BLOCK_SIZE * BYTES_PER_WORD) == 0) # Too lazy to deal with padding.
    num_blocks = int(len(plaintext) / (BLOCK_SIZE * BYTES_PER_WORD))
    assert(num_blocks <= MAX_BLOCKS)
    if(nonce is None):
        nonce = [cgf2n(regint.get_random(bit_length=8)) for _ in range(12)]
    counters = [[cgf2n(x) for x in i.to_bytes(length=4)] for i in range(num_blocks)]

    # AES block cipher setup (performs key expansion)
    aes = AESCipher(key)

    # encrypt
    plaintext = [plaintext[i * (BLOCK_SIZE*BYTES_PER_WORD) : (i+1) * (BLOCK_SIZE*BYTES_PER_WORD)] for i in range(num_blocks)]
    ciphertext = plaintext
    for i in range(num_blocks):
        cipher_output = aes.cipher(nonce + counters[i])
        # print_ln("[Encrypt] cipher_output in block %s = %s", i, [x.reveal() for x in cipher_output])
        # print_ln("[Encrypt] plaintext in block %s = %s", i, [x.reveal() for x in plaintext[i]])
        # print_ln("[Encrypt] ciphertext before XOR in block %s = %s", i, [x.reveal() for x in ciphertext[i]])
        for j in range(BLOCK_SIZE * BYTES_PER_WORD):
            ciphertext[i][j] = plaintext[i][j] + cipher_output[j]
        # print_ln("[Encrypt] ciphertext after XOR in block %s = %s", i, [x.reveal() for x in ciphertext[i]])
    return nonce, [c for block in ciphertext for c in block]
    
def aes_ctr_decrypt(key: list[sgf2n], ciphertext: list[sgf2n], nonce: list[cgf2n]) -> list[sgf2n]:
    '''
    Decrypt ciphertext with nonce and key using AES-CTR mode. 
    '''
    # validate ciphertext length and nonce length and set up counters based on ciphertext length
    assert(len(ciphertext) % (BLOCK_SIZE * BYTES_PER_WORD) == 0) # Too lazy to deal with padding.
    assert(len(nonce) == 12) 
    num_blocks = int(len(ciphertext) / (BLOCK_SIZE * BYTES_PER_WORD))
    assert(num_blocks <= MAX_BLOCKS)
    counters = [[cgf2n(x) for x in i.to_bytes(length=4)] for i in range(num_blocks)]

    # AES block cipher setup (performs key expansion)
    num_rounds_from_key_length = {16: 10, 24: 12, 32: 14}
    assert(len(key) in num_rounds_from_key_length)  # key must be 16, 24, or 32 bytes long
    aes = AESCipher(key)

    # decrypt
    ciphertext = [ciphertext[i * (BLOCK_SIZE*BYTES_PER_WORD) : (i+1) * (BLOCK_SIZE*BYTES_PER_WORD)] for i in range(num_blocks)]
    plaintext = ciphertext
    for i in range(num_blocks):
        cipher_output = aes.cipher(nonce + counters[i])
        # print_ln("[Decrypt] cipher_output in block %s = %s", i, [x.reveal() for x in cipher_output])
        # print_ln("[Decrypt] ciphertext in block %s = %s", i, [x.reveal() for x in ciphertext[i]])
        # print_ln("[Decrypt] plaintext before XOR in block %s = %s", i, [x.reveal() for x in plaintext[i]])
        for j in range(BLOCK_SIZE * BYTES_PER_WORD):
            plaintext[i][j] = ciphertext[i][j] + cipher_output[j]
        # print_ln("[Decrypt] plaintext after XOR in block %s = %s", i, [x.reveal() for x in plaintext[i]])
    return [b for block in plaintext for b in block]


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
        msg_raw = "68656C6C6F2074686572652C207468697320697320612074657374206D65737361676520666F7220656E6372797074696F6E20707572706F7365732E2E2E2E2E"
        msg = [sgf2n(x) for x in str_to_hex(msg_raw)]
        # test_nonce = [cgf2n(x) for x in [0xb2, 0xf, 0x14, 0xbd, 0x91, 0x25, 0xd8, 0x48, 0xa7, 0xa7, 0x30, 0x1a]]
        nonce, ct = aes_ctr_encrypt(key, msg)
        pt = aes_ctr_decrypt(key, ct, nonce)
        print_ln("msg=%s", [x.reveal() for x in msg])
        print_ln("final plaintext = %s", [x.reveal() for x in pt])
        print_ln("error pattern = %s", [x.reveal() + y.reveal() for x,y in zip(msg, pt)])

    compiler.compile_func()