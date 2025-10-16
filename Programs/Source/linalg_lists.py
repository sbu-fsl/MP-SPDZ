#!/usr/bin/env python3

import os, sys
# add MP-SPDZ dir to path so we can import from Compiler
sys.path.insert(0, os.path.dirname(sys.argv[0]) + '/../..') 
from Compiler.library import print_ln, for_range, while_do, break_loop, if_, if_e, else_, print_ln_if
from Compiler.types import sint, cint, Matrix, Array, sgf2n, cgf2n, regint, _secret
from Compiler.compilerLib import Compiler # only used for testing

def create_vandermonde_matrix(num_rows: int, num_cols: int, value_type: cint | sint | cgf2n | sgf2n, eval_points: list = None) -> list[list]:
    '''
    Creates a Vandermonde matrix from the given parameters. Recall that a Vandermonde 
    matrix V allows for evaluation of a polynomial p at multiple points via matrix multiplication.
    We choose to support left multiplication: given a Vandermonde matrix V over evaluation points 
    x_1,...,x_{num_rows}},  and a column matrix P of coefficients c_1,...,c_{num_cols} corresponding to a degree 
    d = num_cols-1 polynomial p, the matrix multiplication V*P yields a column matrix with the points
    p(x_1),...,p(x_{num_rows}).

    :param num_rows: Number of rows; corresponds to the number of evaluation points.
    :type num_rows: int
    :param num_cols: Number of columns; corresponds to degree of the polynomial being evaluated: num_cols = degree + 1.
    :type num_cols: int
    :param value_type: Runtime MP-SPDZ data type of entries
    :type value_type: cint, sint, cgf2n, sgf2n
    :param eval_pts: Optional list of explicit evaluation points. The type of the points must be the clear type corresponding to value_type (e.g., if value_type == sgf2n, then eval_pts must be list[cgf2n]). If given, the length of the list must equal num_rows. If eval_pts=None, we default to eval_pts=[1,...,num_rows] (where integers are actually the clear type versions of value_type). 
    :type eval_pts: list[cint] if value_type == sint, list[cgf2n] if value_type == sgf2n, or None

    :return: A Vandermonde matrix with num_rows rows and num_cols columns.
    :rtype: list[list[value_type]]
    '''
    assert(value_type in (sint, cint, sgf2n, cgf2n))
    clear_from_value = {sint: cint, sgf2n: cgf2n, cint: cint, cgf2n: cgf2n}
    clear_type = clear_from_value[value_type]

    if eval_points is None:
        eval_points = [clear_type(i) for i in range(1, num_rows + 1)]
    else:
        assert all(type(x) == clear_type for x in eval_points)
        assert len(eval_points) == num_rows
    
    V = [[value_type(1) for _ in range(num_cols)] for _ in range(num_rows)]
    for row in range(num_rows):
        for col in range(1, num_cols):
            V[row][col] = V[row][col-1] * eval_points[row]

    return V


######## Testing ########
if __name__ == "__main__":
    usage = "usage: %prog [options] [args]"
    compiler = Compiler(usage=usage)

    compiler.parser.add_option(
        "--vandermonde", 
        dest="vandermonde", 
        action="store_false",
        default=True,
        help="Disable Vandermonde tests"
    )
    compiler.parse_args()


    @compiler.register_function('test_linalg_lists')
    def test_linalg():

        if compiler.options.vandermonde:
            print_ln("VANDERMONDE TESTS")
            print_ln("---- Test 1 (3x3, sint) ----")
            V = create_vandermonde_matrix(3, 3, sint)
            V = [[x.reveal() for x in row] for row in V]
            expected_V = [
                [cint(1), cint(1), cint(1)],
                [cint(1), cint(2), cint(4)],
                [cint(1), cint(3), cint(9)]
            ]
            error_pattern = [[x-y for x,y in zip(row1, row2)] for row1,row2 in zip(V, expected_V)]
            @if_e(sum([x for row in error_pattern for x in row]))
            def _():
                print_ln("FAILED\nV=%s\nexpected_V=%s", V, expected_V)
            @else_
            def _():
                print_ln("PASSED")

            print_ln("---- Test 2 (3x5, sint) ----")
            V = create_vandermonde_matrix(3, 5, sint)
            V = [[x.reveal() for x in row] for row in V]
            expected_V = [
                [cint(1), cint(1), cint(1), cint(1), cint(1)],
                [cint(1), cint(2), cint(4), cint(8), cint(16)],
                [cint(1), cint(3), cint(9), cint(27), cint(81)]
            ]
            error_pattern = [[x-y for x,y in zip(row1, row2)] for row1,row2 in zip(V, expected_V)]
            @if_e(sum([x for row in error_pattern for x in row]))
            def _():
                print_ln("FAILED\nV=%s\nexpected_V=%s", V, expected_V)
            @else_
            def _():
                print_ln("PASSED")


            print_ln("---- Test 3 (5x3, sint) ----")
            V = create_vandermonde_matrix(5, 3, sint)
            V = [[x.reveal() for x in row] for row in V]
            expected_V = [
                [cint(1), cint(1), cint(1)],
                [cint(1), cint(2), cint(4)],
                [cint(1), cint(3), cint(9)],
                [cint(1), cint(4), cint(16)],
                [cint(1), cint(5), cint(25)]
            ]
            error_pattern = [[x-y for x,y in zip(row1, row2)] for row1,row2 in zip(V, expected_V)]
            @if_e(sum([x for row in error_pattern for x in row]))
            def _():
                print_ln("FAILED\nV=%s\nexpected_V=%s", V, expected_V)
            @else_
            def _():
                print_ln("PASSED")


            print_ln("---- Test 4 (3x3, cint, eval_points) ----")
            # test 3x3 cint matrix with explicit evaluation points
            eval_points = [cint(3), cint(5), cint(7)]
            V = create_vandermonde_matrix(3, 3, cint, eval_points=eval_points)
            expected_V = [
                [cint(1), cint(3), cint(9)],
                [cint(1), cint(5), cint(25)],
                [cint(1), cint(7), cint(49)]
            ]
            error_pattern = [[x-y for x,y in zip(row1, row2)] for row1,row2 in zip(V, expected_V)]
            @if_e(sum([x for row in error_pattern for x in row]))
            def _():
                print_ln("FAILED\nV=%s\nexpected_V=%s", V, expected_V)
            @else_
            def _():
                print_ln("PASSED")


            print_ln("---- Test 6 (3x3, sgf2n) ----")
            V = create_vandermonde_matrix(3,3,sgf2n)
            V = [[x.reveal() for x in row] for row in V]
            expected_V = [
                [cgf2n(1), cgf2n(1), cgf2n(1)],
                [cgf2n(1), cgf2n(2), cgf2n(4)],
                [cgf2n(1), cgf2n(3), cgf2n(5)] 
            ]
            error_pattern = [[x-y for x,y in zip(row1, row2)] for row1,row2 in zip(V, expected_V)]
            @if_e(sum([x for row in error_pattern for x in row]))
            def _():
                print_ln("FAILED\nV=%s\nexpected_V=%s", V, expected_V)
            @else_
            def _():
                print_ln("PASSED")

    compiler.compile_func()