#!/usr/bin/env python3

import os, sys
# add MP-SPDZ dir to path so we can import from Compiler
sys.path.insert(0, os.path.dirname(sys.argv[0]) + '/../..') 
from Compiler.library import print_ln
from Compiler.types import cgf2n
from Compiler.compilerLib import Compiler


if __name__ == "__main__":
    compiler = Compiler()

    @compiler.register_function("gf2n")
    def test_gf2n():
        '''
        Can set the "bit length" n for GF(2^n) at runtime with 
        the -lg2 argument (default 128). 
        The statistical security parameter -S (default 40) must be at most n.  For example: 
        $> Scripts/mascot.sh -lg2 4 -S 4 gf2n
        '''
        a = cgf2n(1) # 0x01
        print_ln("a=%s", a)
        b = cgf2n(2) # 0x02
        print_ln("b=%s", b)
        c = cgf2n(3) # 0x03
        print_ln("c=%s", c)
        d = cgf2n(15) # 0x0f
        print_ln("d=%s", d)
        e = cgf2n(16) # 0x10
        print_ln("e=%s", e)
        f = cgf2n(17) # 0x11
        print_ln("f=%s", f)
        g = cgf2n(18) # 0x12
        print_ln("g=%s", g)
        h = cgf2n(48) # 0x30
        print_ln("h=%s", h)
        i = cgf2n(0x12)
        print_ln("i=%s", i)

        '''
        For -lg2 4 -S 4: (field on 16 elements)
        a=1
        b=2
        c=3
        d=f
        e=0
        f=1
        g=2
        h=0
        i=2
        
        For -lg 5 -S 5: (field on 32 elements)
        a=1
        b=2
        c=3
        d=f
        e=10
        f=11
        g=12
        h=10
        i=12
        '''

        cc = c.bit_decompose(4)
        print_ln("cc=%s", cc) # [1, 1, 0, 0]. just bit_decompose() has [1,1,0...0] length 40.

    compiler.compile_func()
