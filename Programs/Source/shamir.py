#!/usr/bin/env python3

import os, sys
# add MP-SPDZ dir to path so we can import from Compiler
sys.path.insert(0, os.path.dirname(sys.argv[0]) + '/../..') 
from Compiler.library import print_ln
from Compiler.types import sint, cint, Array
from Compiler.compilerLib import Compiler # only used for testing

# we assume these modules reside in Programs/Source/ 
from linalg import LUSolver, create_vandermonde_matrix

def shamir_share(msg, threshold, num_parties, eval_points=None, rand=None):
    '''
    Perform textbook Shamir's secret sharing

    :param msg: Secret message to be secret shared, numerical type interpreted according to computation domain (e.g., prime field for arithmetic circuit domain).
    :param threshold: Reconstruction / privacy threshold. Must be less than num_parties.
    :param num_parties: Number of shareholders.
    :param eval_pts: Optional public Array of explicit evaluation points.
    - If given, the length of the list must equal :param num_parties
    - If eval_pts=None, we default to eval_pts=[1,...,num_parties]. Note these integers will be interpreted according to the computation domain.
    :param rand: Optional secret Array of random coefficients to use. Length must equal threshold
    
    :returns: A 2-tuple, where the first item holds evaluation points, and the second iten holds corresponding polynomial evaluations. 
    In other words, if we have (eval_points, poly_evals) = shamir_share(...), then (eval_points[i], poly_evals[i]) is the i-th Shamir 
    share of the form (x, p(x))
    '''
    assert threshold <= num_parties
    if eval_points is None:
        eval_points = Array(num_parties, sint).assign([i for i in range(1,num_parties+1)]) # TODO: do we need sint if we return eval_points as part of tuple?
    
    V = create_vandermonde_matrix(num_parties, threshold - 1, sint, eval_points)
    poly_coeffs = Array(threshold, sint)
    if rand:
        poly_coeffs.assign(rand)
    else:
        poly_coeffs.randomize()
    poly_coeffs[0] = msg
    poly_evals = V.dot(poly_coeffs)
    # shares = Matrix(2, num_parties, sint)
    # shares[0] = eval_points
    # shares[1] = poly_evals
    return eval_points, poly_evals

def shamir_reconstruct(eval_points, poly_evals):
    '''
    Attempts best-effort reconstruction of Shamir secret shares from :param shares
    
    :param eval_points: Evaluation points used during the sharing phase. 
    :param poly_evals: Polynomial evaluations corresponding to :param eval_points. 
    In other words, (eval_points[i], poly_evals[i]) is a single Shamir secret share of the form (x, p(x))
    :returns: A reconstructed secret, interpreted as a field element in the computation domain
    '''
    V = create_vandermonde_matrix(eval_points.length, eval_points.length - 1, cint)
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
        msg = sint(12)
        rand = sint.Array(2).assign([4, 5])
        print_ln("msg=%s", msg.reveal())
        eval_points, poly_evals = shamir_share(msg, 2, 3, rand=rand)
        print_ln("eval_points=%s \npoly_evals=%s", eval_points.reveal(), poly_evals.reveal())
        reconstructed_msg = shamir_reconstruct(eval_points, poly_evals)
        print_ln("reconstructed_msg=%s", reconstructed_msg.reveal())

    compiler.compile_func()  