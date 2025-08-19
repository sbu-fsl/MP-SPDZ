# doing this as a python script for educational purposes. 

# Chris: had to add following two lines so that Compiler module can be found
# Chris: compile as: python3 Programs/Source/arithmetic.mpc
# Chris: run as: Scripts/mascot.sh testmpc
import os, sys
sys.path.insert(0, os.path.dirname(sys.argv[0]) + '/../..')

from Compiler.library import print_ln
from Compiler.types import sint, Matrix
from Compiler.compilerLib import Compiler


usage = "usage: %prog" # not taking any arguments for now
compiler = Compiler(usage=usage)


@compiler.register_function('arithmetic')
def main():
    a = sint(1)
    b = sint(2)
    print_ln("a/b=%s", a.field_div(b).reveal())

    # print_ln("prime=%s", program.prime) # None
    # # https://github.com/data61/MP-SPDZ/issues/271 looks like default prime is 128 bits. Where is it though? Seems like determined at runtime?
    # print_ln("security=%s", program.security)

    # example input
    # echo 65 > Player-Data/Input-P0-0
    # echo 72 > Player-Data/Input-P1-0
    # echo 97 > Player-Data/Input-P2-0

    # prints at compile time!
    # import numpy as np
    # a = np.asarray([1,2])
    # print(a)

    # create a matrix: create list. instantiate Matrix. assign list to Matrix.
    v = [
        [1, 1, 1, 1],
        [1, 2, 4, 8],
        [1, 3, 9, 27]
    ]
    V = sint.Matrix(3,4) # we choose sint because we later concatenate a secret column to V. 
    V.assign(v)

    # Example: create a random (secret) vector and compute V*c
    c = sint.Array(4)
    c.randomize() # randomize modifies c in place and returns None
    res = V.dot(c)

    # use Player-Data to construct vector s, create augmented matrix [V|s]
    s0 = sint.get_input_from(0)
    s1 = sint.get_input_from(1)
    s2 = sint.get_input_from(2)
    s = Matrix.create_from([s0, s1, s2])
    V_aug_s = V.concat_columns(s)
    # print_ln("%s", V_aug_s.reveal())

    # it seems like numpy doesn't really work with run-time data types
    # a1 = np.asarray([s0], dtype=object)
    # a2 = np.asarray([s1], dtype=object)
    # a3 = a1 + a2
    # a4 = a3.tolist()
    # a5 = Matrix.create_from(a4) # Compiler.exceptions.CompilerError: cannot create Array of <class 'list'>
    # print_ln("a5=%s", a5.reveal())


if __name__ == "__main__":
    compiler.compile_func()