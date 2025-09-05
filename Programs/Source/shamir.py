#!/usr/bin/env python3

import os, sys
# add MP-SPDZ dir to path so we can import from Compiler
sys.path.insert(0, os.path.dirname(sys.argv[0]) + '/../..') 
from Compiler.library import print_ln
from Compiler.types import sint, cint, Array, sgf2n, cgf2n
from Compiler.compilerLib import Compiler # only used for testing

# we assume these modules reside in Programs/Source/ 
from linalg import LUSolver, create_vandermonde_matrix
from aes import apply_field_embedding, apply_inverse_field_embedding

def shamir_share(msg: sint|sgf2n, threshold: int, num_parties: int, eval_points:Array=None, rand:Array=None) -> tuple[Array, Array]:
    '''
    Perform textbook Shamir's secret sharing

    :param msg: Secret message to be secret shared, numerical type interpreted according to computation domain (e.g., prime field for arithmetic circuit domain).
    :param threshold: Reconstruction / privacy threshold. Must be less than num_parties.
    :param num_parties: Number of shareholders.
    :param eval_points: Optional public Array of explicit evaluation points.
    - If given, the length of the list must equal :param num_parties
    - If eval_pts=None, we default to eval_pts=[1,...,num_parties]. Note these integers will be interpreted according to the computation domain.
    :param rand: Optional secret Array of random coefficients to use. Length must equal threshold
    
    :returns: A 2-tuple, where the first item holds evaluation points, and the second iten holds corresponding polynomial evaluations. 
    In other words, if we have (eval_points, poly_evals) = shamir_share(...), then (eval_points[i], poly_evals[i]) is the i-th Shamir 
    share of the form (x, p(x))
    '''
    assert threshold <= num_parties
    msg_type = type(msg)
    if eval_points is None:
        eval_points = Array(num_parties, msg_type).assign([i for i in range(1,num_parties+1)])
    
    V = create_vandermonde_matrix(num_parties, threshold - 1, msg_type, eval_points)
    poly_coeffs = Array(threshold, msg_type)
    if rand:
        poly_coeffs.assign(rand)
    else:
        if msg_type == sgf2n:
            for i in range(len(poly_coeffs)):
                poly_coeffs[i] = sgf2n.bit_compose([sgf2n.get_random_bit() for _ in range(128)])
        else:
            poly_coeffs.randomize()
    poly_coeffs[0] = msg
    poly_evals = V.dot(poly_coeffs)
    return eval_points, poly_evals

def shamir_reconstruct(poly_evals: Array | list[sint | sgf2n], eval_points:Array=None) -> sint | sgf2n:
    '''
    Attempts best-effort reconstruction of Shamir secret shares from :param shares
     
    :param poly_evals: Polynomial evaluations corresponding to :param eval_points. 
    :param eval_points: Optional public Array of explicit evaluation points. Defaults to 1,2,...,len(poly_evals)
    In other words, (eval_points[i], poly_evals[i]) is a single Shamir secret share of the form (x, p(x))
    :returns: A reconstructed secret, interpreted as a field element in the computation domain
    '''
    # workaround to make this work regardless of whether poly_evals is list or Array
    msg_type = type(poly_evals[0])
    n = 0
    if type(poly_evals) == list:
        n = len(poly_evals)
    else:
        n = poly_evals.length
    poly_evals = Array(n, msg_type).assign(poly_evals)
    if eval_points is None:
        eval_points = Array(n, msg_type).assign([i for i in range(1,n+1)])
    V = create_vandermonde_matrix(n, n - 1, msg_type)
    solver = LUSolver(V)
    poly_coeffs = solver.solve(poly_evals)
    return poly_coeffs[0]

if __name__ == "__main__":
    usage = "usage: %prog [options] [args]"
    compiler = Compiler(usage=usage)

    @compiler.register_function("test_shamir")
    def test_shamir():
        print_ln("SHAMIR TESTS")

        print_ln("----Test 1----")
        msg = sint(246432)
        rand = sint.Array(2).assign([4, 5])
        print_ln("msg=%s", msg.reveal())
        eval_points, poly_evals = shamir_share(msg, 2, 3, rand=rand)
        print_ln("eval_points=%s \npoly_evals=%s", eval_points.reveal(), poly_evals.reveal())
        reconstructed_msg = shamir_reconstruct(eval_points, poly_evals)
        print_ln("reconstructed_msg=%s", reconstructed_msg.reveal())

        print_ln("----Test 2 (sgf2n with a random array)----")
        msg = sgf2n(700)
        rand = Array(2, sgf2n).assign([88, 57])
        print_ln("msg=%s", msg.reveal())
        eval_points, poly_evals = shamir_share(msg, 2, 3, rand=rand)
        print_ln("eval_points=%s \npoly_evals=%s", eval_points.reveal(), poly_evals.reveal())
        reconstructed_msg = shamir_reconstruct(eval_points, poly_evals)
        print_ln("reconstructed_msg=%s", sgf2n(reconstructed_msg.reveal()).reveal())


        print_ln("----Test 3 (sgf2n)----")
        msg = sgf2n(1283890184043)
        print_ln("msg=%s", msg.reveal())
        eval_points, poly_evals = shamir_share(msg, 2, 3)
        print_ln("eval_points=%s \npoly_evals=%s", eval_points.reveal(), poly_evals.reveal())
        reconstructed_msg = shamir_reconstruct(eval_points, poly_evals)
        print_ln("reconstructed_msg=%s", sgf2n(reconstructed_msg.reveal()).reveal())

        print_ln("----Test 4 (buggy case discovered from pss.py)----")
        msg = apply_field_embedding(sgf2n(114))
        print_ln("msg_unembedded=%s, msg_embedded=%s", cgf2n(114), msg.reveal())
        eval_points_embedded = Array(3, sgf2n).assign([apply_field_embedding(sgf2n(i)) for i in range(1,4)])
        randomness_embedded = [apply_field_embedding(sgf2n.bit_compose([sgf2n.get_random_bit() for _ in range(8)])) for i in range(2)]
        randomness_embedded = Array(2,sgf2n).assign(randomness_embedded) 
        eval_points, poly_evals = shamir_share(msg, 2, 3, eval_points=eval_points_embedded, rand=randomness_embedded)
        print_ln("eval_points=%s \npoly_evals_embedded=%s, poly_evals_unembedded=%s,", eval_points.reveal(), poly_evals.reveal(), [apply_inverse_field_embedding(x).reveal() for x in poly_evals])
        reconstructed_msg = shamir_reconstruct(poly_evals, eval_points=eval_points_embedded)
        print_ln("reconstructed_msg_embedded=%s, reconstructed_msg_unembedded=%s", reconstructed_msg.reveal(), apply_inverse_field_embedding(reconstructed_msg).reveal())

        print_ln("----Test 5 (a simple case)----")
        msg = apply_field_embedding(sgf2n(1))
        print_ln("msg_unembedded=%s, msg_embedded=%s", cgf2n(1), msg.reveal())
        eval_points_embedded = Array(3, sgf2n).assign([apply_field_embedding(sgf2n(i)) for i in range(1,4)])
        randomness_embedded = [apply_field_embedding(sgf2n.bit_compose([sgf2n.get_random_bit() for _ in range(8)])) for i in range(2)]
        randomness_embedded = Array(2,sgf2n).assign(randomness_embedded) 
        eval_points, poly_evals = shamir_share(msg, 2, 3, eval_points=eval_points_embedded, rand=randomness_embedded)
        print_ln("eval_points=%s \npoly_evals_embedded=%s, poly_evals_unembedded=%s,", [x.reveal() for x in eval_points_embedded], poly_evals.reveal(), [apply_inverse_field_embedding(x).reveal() for x in poly_evals])
        reconstructed_msg = shamir_reconstruct(poly_evals, eval_points=eval_points_embedded)
        print_ln("reconstructed_msg_embedded=%s, reconstructed_msg_unembedded=%s", reconstructed_msg.reveal(), apply_inverse_field_embedding(reconstructed_msg).reveal())

        print_ln("----Test 6 (another bug from pss.py)----")
        eval_points_embedded = Array(3, sgf2n).assign([apply_field_embedding(sgf2n(i)) for i in range(1,4)])
        poly_evals_embedded = [apply_field_embedding(x) for x in [sgf2n(198), sgf2n(64), sgf2n(203)]]
        print_ln("poly_evals_embedded=%s", [x.reveal() for x in poly_evals_embedded])
        reconstructed_msg = shamir_reconstruct(poly_evals_embedded, eval_points=eval_points_embedded)
        print_ln("reconstructed_msg_embedded=%s, reconstructed_msg_unembedded=%s", reconstructed_msg.reveal(), apply_inverse_field_embedding(reconstructed_msg).reveal())

    compiler.compile_func()  