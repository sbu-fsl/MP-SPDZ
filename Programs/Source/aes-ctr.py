from aes import AESCipher, apply_field_embedding, apply_inverse_field_embedding, BLOCK_SIZE, BYTES_PER_WORD
from Compiler.library import print_ln
from Compiler.types import cgf2n, sgf2n
from Compiler.compilerLib import Compiler
from Compiler.GC.types import sbits # requires binary circuits. Have to compile with -X or -Y option for daBit / edaBit

# need to set up a random counter 
# sgf2n has no random method, so maybe we can use _secret.get_random_bit() and sgf2n.bit_compose? 
# have to test it out.

class AES_CTR():
    def __init__(self, aes: AESCipher, key: bytes):
        self.aes = aes
        self.key = key
        self.counter = [sgf2n.get]

if __name__ == "__main__":
    usage = "usage: %prog [options] [args]"
    compiler = Compiler(usage=usage)
    compiler.parse_args()

    @compiler.register_function("main")
    def main():
        # sb128 = sbits.get_type(128) # don't know if we can use this bc don't know whether we can convert to sgf2n somehow
        # sgf2n.bit_type = sgf2n... what does this mean/imply? 

        # _secret.bit_compose(bits) says it works whenever bits is an iterable of any type convertible to sint
        #interestingly, if we check the code for _secret.bit_compose
        a = sgf2n.get_random_bit()
        print_ln("a=%s", a.reveal())

    compiler.compile_func()