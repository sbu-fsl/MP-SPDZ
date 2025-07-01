#!/usr/bin/env python3

import os, sys, math
# add MP-SPDZ dir to path so we can import from Compiler
sys.path.insert(0, os.path.dirname(sys.argv[0]) + '/../..') 
from Compiler.library import print_ln, for_range
from Compiler.types import sint, cint, Matrix, Array
from Compiler.compilerLib import Compiler # only used for testing

# we assume these modules reside in Programs/Source/ 
from linalg import LUSolver
from shamir import shamir_share, shamir_reconstruct


def lr_share(msg, threshold, num_parties, leakage_budget, leakage_error, prime_modulus_bit_length):
    '''
    LRShare algorithm of CKOS22 leakage-resilient secret sharing scheme.
    https://eprint.iacr.org/2022/216
    We hardcode Shamir's secret sharing for MShare and SdShare.
    The consequence is that the resulting LRSS scheme enjoys perfect privacy and perfect local uniformity, 
    and we can easily set the leakage error with a single parameter :param leakage_error.

    :param msg: Secret message to be secret shared, numerical type interpreted according to computation domain (e.g., prime field for arithmetic circuit domain).
    :param threshold: Reconstruction / privacy threshold. Must be less than num_parties.
    :param num_parties: Number of shareholders.
    :param leakage_budget: Amount of leakage (in bits) from each share the scheme should support (see mu in CKOS22)
    :param leakage_error: For any leakage function f in the local leakage family specified by :param leakage_budget, 
        and for any two messages m, m', :param leakage_error is the maximum possible statistical distance between f(Share(m)) and f(Share(m')) (see epsilon_lr in CKOS22).
    :param prime_modulus_bit_length: Number of bits needed to represent the size of the prime field being used. 
    
    :returns: Leakage-resilient secret shares as a 3-tuple. The first item in the tuple is an Array of evaluation points. These are the evaluation points
    used in every invocation of Shamir's secret sharing within the LRSS scheme. 
    The second item is a Matrix of extractor "sources" (w_i in CKOS22). Each of the :param num_parties rows corresponds to a 
    different w_i. 
    The third item is a Matrix of Shamir shares of the seed. Each of the :param num_parties rows corresponds to a different seed share. 
    The i-th LRSS share is the concatenation of the i-th coordinate of each item in the returned tuple. That is, if we have 
    (eval_points, sources, seed_shares) = lr_share(...), then (eval_points[i], sources[i], seed_shares[i]) is the i-th LRSS share. 
    
    '''
    
    #first Shamir secret share msg 
    eval_points, msg_shares = shamir_share(msg, threshold, num_parties)
    
    # The next step is supposed to be generating a random seed for use with InvExt. But how long is this seed supposed to be?:
    # Because our extractor output should be a single field element, i.e. l = 1 * log |F|,
    # then by rules of matrix multiplication our seed length d = n = source length. 
    # (In other words, our extractor is just a dot product where one of the vectors is a random seed)
    # To determine source length, we have to determine a bunch of other parameters, starting with the extractor error.
    # NOTE: the following params need to be computed at compile-time (i.e., with normal python types) because sizes of Container types depend on these. 
    ext_error = leakage_error / (6 * num_parties) # see theorem 1 CKOS22
    # By leftover hash lemma, we must have min_entropy >= l + 2log(1/ext_error) - 2 = log|F| + 2log(1/ext_error) - 2
    min_entropy = prime_modulus_bit_length + 2 * math.log2(1 / ext_error) - 2 
    # In CKOS22, source length in bits is equal to min_entropy + leakage_budget + 1. Why is this?:
    # First notice our source must retain at least 'min_entropy' bits of min-entropy in order for our extractor to be secure by leftover hash lemma. 
    # Second, notice the source is output as part of the final secret share, so it will be subject to leakage, and it must retain the 
    # required min-entropy even AFTER the leakage, so we are at least forced to increase the length by leakage_budget bits.
    # Finally, notice when we actually sample the source below by solving a linear system with one equation in source_length unknowns,
    # there are (source_length - 1) free variables that we set to random field elements. 
    # This is because exactly one variable is determined by the system, so this variable does not contribute any entropy to the source.
    # As a result, we add 1 to the source length in order to get a source with min-entropy (min_entropy + leakage_budget). 
    # This way, no matter how the adversary attains leakage_budget bits of information from the source, its min-entropy can only decrease
    # by at most leakage_budget bits.
    # Also, for us, source_length must correspond to a number of field elements, so divide by log|F| and round up.
    source_length = math.ceil( (min_entropy + leakage_budget) / prime_modulus_bit_length )
    seed_length = source_length

    # generate random secret seed
    seed = Matrix(1, seed_length, sint) # needs to be Matrix if using LUSolver
    seed.randomize()
    # print_ln("seed=%s", seed.reveal())

    # run InvExt on each msg share (just the p(x) part - 2nd row) to obtain source w_i. Collect {w_i} into sources
    # NOTE: because we are using inner product extractor, every share must have a preimage, so no need to check for failures in InvExt.
    # NOTE: if above is true and always exists a solution, then might be able to get same leakage_error with bigger ext_error => potentially smaller seed
    sources = Matrix(num_parties, source_length, sint)
    seed_solver = LUSolver(seed)
    @for_range(num_parties)
    def _(i):
        sources[i] = seed_solver.solve( Array(1, sint).assign_all(msg_shares[i]), free_vars='rand') 
    
    # Shamir secret share the seed. Use same eval points as msg
    seed_shares_T = Matrix(seed_length, num_parties, sint) # _T to indicate we will transpose this later
    @for_range(seed_length)
    def _(i):
        seed_shares_T[i] = shamir_share(seed[0][i], 2, num_parties, eval_points)[1]
    seed_shares = seed_shares_T.transpose() # num_parties x seed_length. i-th row of seed_shares is now the i-th "share" of the seed that goes into i-th LRSS share. Corresponds to i-th eval point.
    
    return eval_points, sources, seed_shares # i-th LRSS share is collection of i-th coord of each of these

def lr_reconstruct(eval_points, sources, seed_shares):
    '''
    Attempts best-effort reconstruction of LRSS shares. 

    :param eval_points: Array of evaluation points used in every invocation of Shamir's secret sharing within the LRSS scheme.
    :param sources: Matrix of extractor "sources" (w_i in CKOS22). Each row corresponds to a different w_i. 
    :param seed_shares: Matrix of Shamir shares of the seed. Each row corresponds to a different seed share.

    The reconstruction function assumes that the i-th coordinate of each parameter belongs to the i-th LRSS share. 
    That is, (eval_points[i], sources[i], seed_shares[i]) is the i-th LRSS share. 

    :returns A single field element as an sint.
    '''
    
    # Do all inputs have same length?
    num_shares = eval_points.length
    assert sources.sizes[0] == num_shares
    assert seed_shares.sizes[0] == num_shares

    # Attempt best-effort seed reconstruction (we will reconstruct something, just no guarantee it's the seed that was used.)
    seed_shares_T = seed_shares.transpose() # seed_length x num_shares. For convenience in what we do next. 
    seed_length = seed_shares_T.sizes[0]
    seed = Array(seed_length, sint)
    @for_range(seed_length)
    def _(i):
        seed[i] = shamir_reconstruct(eval_points, seed_shares_T[i])

    # Run extractor on each source using seed (just dot product of seed with each source)
    # TODO: would it be simpler to just do shamir_poly_evals = sources.dot(seed) all in one go?
    shamir_poly_evals = Array(num_shares, sint)
    @for_range(num_shares)
    def _(i):
        shamir_poly_evals[i] = seed.dot(sources[i])[0]

    # Attempt best-effort reconstruction of message.
    msg = shamir_reconstruct(eval_points, shamir_poly_evals)
    return msg

if __name__ == "__main__":
    usage = "usage: %prog [options] [args]"
    compiler = Compiler(usage=usage)
    
    @compiler.register_function("test_lrss")
    def test_lrss():
        print_ln("TEST CKOS22 LRSS")
        PRIME_MODULUS = compiler.prog.prime
        PRIME_MODULUS_BIT_LENGTH = math.ceil(math.log2(PRIME_MODULUS))


        print_ln("----Test 1----") # basic easy test
        msg = cint(1)
        eval_points, sources, seed_shares = lr_share(
            msg=msg,
            threshold=3, 
            num_parties=3, 
            leakage_budget=1, 
            leakage_error=0.25, 
            prime_modulus_bit_length=PRIME_MODULUS_BIT_LENGTH)
        print_ln("eval_points=%s", eval_points.reveal())
        print_ln("sources=%s", sources.reveal())
        print_ln("seed_shares=%s", seed_shares.reveal())
        reconstructed_msg = lr_reconstruct(eval_points, sources, seed_shares)
        print_ln("reconstructed_msg=%s", reconstructed_msg.reveal())

        print_ln("----Test 2----") # higher leakage budget and much stricter leakage error
        msg = cint(2025)
        eval_points, sources, seed_shares = lr_share(
            msg=msg, 
            threshold=4, 
            num_parties=5, 
            leakage_budget=32, 
            leakage_error=1 / (2**128), 
            prime_modulus_bit_length=PRIME_MODULUS_BIT_LENGTH)
        print_ln("eval_points=%s", eval_points.reveal())
        print_ln("sources=%s", sources.reveal())
        print_ln("seed_shares=%s", seed_shares.reveal())
        reconstructed_msg = lr_reconstruct(eval_points, sources, seed_shares)
        print_ln("reconstructed_msg=%s", reconstructed_msg.reveal())

    compiler.compile_func()