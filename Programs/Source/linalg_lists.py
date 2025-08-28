#!/usr/bin/env python3

import os, sys
# add MP-SPDZ dir to path so we can import from Compiler
sys.path.insert(0, os.path.dirname(sys.argv[0]) + '/../..') 
from Compiler.library import print_ln, for_range, while_do, break_loop, if_, if_e, else_, print_ln_if
from Compiler.types import sint, cint, Matrix, Array, sgf2n, cgf2n, regint
from Compiler.compilerLib import Compiler # only used for testing

class LUSolver:
    '''
    Perform LU decomposition and solve linear systems.

    :param M: arbitrary Matrix to decompose

    WARNING: the LU decomposition procedure (which is performed at instantiation, or with static method lu) may reveal up to 
    a bit of information for each element on or below the diagonal of M.
    '''

    def __init__(self, M):
        self.M = M # store unmodified original M
        M_copy = Matrix.create_from(M)
        self.P, self.L, self.U = LUSolver.lu(M_copy) # modifies M_copy in place. 

    @staticmethod
    def lu(M):
        '''
        Compute the LUP decomposition of M: PM = LU. See https://en.wikipedia.org/wiki/LU_decomposition. Modifies M in place.
        Return P,L,U

        WARNING: the LU decomposition procedure (which is performed at instantiation, or with static method lu) may reveal up to 
        a bit of information for each element on or below the diagonal of M.
        '''

        def swap_rows(M, i, j):
            '''Swap rows i and j of matrix M in place'''
            # need deeper copy than just tmp = M[i]
            tmp = M[i].same_shape()
            tmp.assign(M[i])
            M[i] = M[j]
            M[j] = tmp

        def create_identity_matrix(n, value_type):
            I = Matrix(n,n,value_type)
            I.assign_all(0)
            @for_range(n)
            def _(i):
                I[i][i] = 1
            return I
        
        # LU procedure is basically identical to Gaussian elimination, except information about the Gaussian elimination process gets "stored" in P and L. 
        (num_rows, num_cols) = (M.sizes[0], M.sizes[1])
        P = create_identity_matrix(num_rows, M.value_type)
        L = create_identity_matrix(num_rows, M.value_type)

        h = cint(0) # pivot row
        k = cint(0) # pivot col
        @while_do(lambda: (h < num_rows) & (k < num_cols))
        def _():
            # print_ln("h=%s, k=%s, M=%s", h, k, M.reveal())
            # choose a pivot by finding non-zero element in column k, and swapping its row with row h.
            @for_range(h, num_rows)
            def _(i):
                @if_((M[i][k] != 0).reveal()) # WARNING: revealing information about M
                def _():
                    swap_rows(M, h, i)
                    swap_rows(P, h, i)
                    break_loop()
            # if no pivot exists, pass to next column
            @if_e((M[h][k] == 0).reveal()) # WARNING: revealing information about M
            def _():
                k.update(k+1)
            # for all rows below the pivot, zero out elements in column below pivot using row ops
            @else_
            def _():
                @for_range(h+1, num_rows)
                def _(i):
                    scale_factor = M[i][k].field_div(M[h][k])
                    L[i][k] = scale_factor
                    @for_range(k, num_cols)
                    def _(j):
                        M[i][j] = M[i][j] - (M[h][j] * scale_factor)
                # increase pivot row and column
                h.update(h+1)
                k.update(k+1)
            return 1 # keep looping
        # print_ln("P=%s\nL=%s\nU=%s", P.reveal(), L.reveal(), M.reveal())
        return (P,L,M)

    def solve(self, b, free_vars=None):
        '''
        Solve the system Mx = b for x, where M = self.M is the matrix passed into the constructor. 

        :param b: Array of appropriate length
        :param free_vars: If system has more than one solution, free_vars is used to set the values of the free variables. 
            - If free_vars=None (default), then assign 0 to every free variable
            - Else if free_vars="rand", then assign a random value to every free variable
            - Else, assign free_vars[i] to the i-th free variable
            
        :returns: Array of appropriate length. If no solution exists, this Array is filled with False.

        WARNING: May reveal up to a bit of information for each element of b
        '''
        assert b.length == self.M.sizes[0]
        y = self._solve_lower(b) # must be the case that this always has a solution
        x = self._solve_upper(y, free_vars)
        return x

    def _solve_lower(self, b):
        '''
        Solve the system Ly = Pb for y with forward substitution, where L = self.L and P = self.P are obtained from 
        LU decomposition of self.M. In particular, this means L is a square lower triangular matrix 
        with ones on the diagonal. This means the system will always have a unique solution. 

        :param b: Array of appropriate length
        :returns: Array of appropriate length
        '''
        # initialize vars
        P = self.P
        L = self.L
        b = P.dot(b) # apply permutation to b
        num_rows = L.sizes[0] # L is square
        y = Array(num_rows, b.value_type).assign_all(0) # solution vector
        # begin forward substitution
        @for_range(num_rows)
        def _(i):
            # compute y[i] = ( b[i] - \sum_{j < i} L[i][j] * y[j] ) / L[i][i]
            y[i] = b[i] 
            @for_range(i) # for all j < i ...
            def _(j):
                y[i] = y[i] - ( L[i][j] * y[j] )
            y[i] = y[i].field_div(L[i][i])
        print("Lower: %s", y)
        return y

    def _solve_upper(self, y, free_vars):
        '''
        Solve the system Ux = y for x with backsubstitution, where U = self.U is the row-echelon form matrix obtained from LU decomposition of self.M.

        In the case that the system has no solution, return an Array filled with False.

        In the case that the system has more than one solution, use free_vars to assign values to the free variables.
            - If free_vars=None, then assign 0 to every free variable
            - Else if free_vars="rand", then assign a random value to every free variable
            - Else, assign free_vars[i] to the i-th free variable
        '''
        U = self.U
        num_rows, num_cols = U.sizes[0], U.sizes[1]
        # initialize solution vector x
        if free_vars is None:
            x = Array(num_cols, y.value_type).assign_all(0)
        elif free_vars == "rand":
            if y.value_type == cgf2n:
                # Bit Compose will work for any field up to size 128
                x = Array(num_cols, cgf2n)
                print_ln("Before Bit compose: %s", x.reveal())
                for i in range(len(x)):
                    x[i] = (cgf2n.bit_compose([cgf2n(regint.get_random(bit_length=64)) for _ in range(2)]))
                print_ln("After Bit compose: %s", x.reveal())
            elif y.value_type == sgf2n:
                # Bit Compose will work for any field up to size 128
                x = Array(num_cols, sgf2n)
                print_ln("Before Bit compose: %s", x.reveal())
                for i in range(0,len(x)):
                    x[i] = sgf2n.bit_compose([sgf2n.get_random_bit() for _ in range(128)])
                print_ln("After Bit compose: %s", x.reveal())
            else:
                x = Array(num_cols, y.value_type)
                print_ln("Array: %s", x.reveal())
                x.randomize()
                print_ln("Randomized Array: %s", x.reveal())
        else:
            x = Array(num_cols, y.value_type).assign(free_vars)

        # the while loop below sets the following three variables
        last_non_zero_row_idx = cint(num_rows - 1) # index of the last non-zero row in U
        last_pivot_idx = cint(-1) # if we find a non-zero row, we will find a pivot in that row. 
        is_unsolvable = sint(0) # boolean flag for whether solution exists. 
        @while_do(lambda: last_non_zero_row_idx >= 0)
        def _():
            # test U[last_non_zero_row_idx] to see if it is all zeroes
            @for_range(last_non_zero_row_idx, num_cols) # we assume U is row echelon form (upper triangular), so only need to check from diagonal and to the right
            def _(j):
                @if_((U[last_non_zero_row_idx][j] != 0).reveal()) # WARNING: leaks info about U
                def _():
                    last_pivot_idx.update(j)
                    break_loop()
            # if we found a non-zero element, break
            @if_(last_pivot_idx != -1)
            def _():
                break_loop()
            # U[last_non_zero_row_idx] is all zeros. If y[last_non_zero_row_idx] != 0, then unsolvable system
            @if_((y[last_non_zero_row_idx] != 0).reveal()) # WARNING: reveals a bit of information about y
            def _():
                is_unsolvable.update(1)
                break_loop()
            # U[last_non_zero_row_idx] is all zeros, but system is still solvable. Go up one row.
            last_non_zero_row_idx.update(last_non_zero_row_idx - 1)
        # print_ln("last_non_zero_row_idx=%s, last_pivot_idx=%s, is_unsolvable=%s", last_non_zero_row_idx, last_pivot_idx, is_unsolvable.reveal())
        # check to see if no solution exists
        @if_e(is_unsolvable.reveal())
        def _():
            x.assign_all(False) # TODO: could we instead just return False?
        @else_
        def _():
            # It is possible U is all zeros, so any vector is vacuously a solution. We return free_vars.
            @if_e(last_non_zero_row_idx == -1)
            def _():
                pass
            @else_
            # at this point, we must have a valid pivot in U[last_non_zero_row_idx][last_pivot_idx], and a solution is guaranteed to exist. 
            # notice that valid pivots only exist on the diagonal, so U[last_non_zero_row_idx][last_pivot_idx] = U[last_pivot_idx][last_pivot_idx] and we can forget about last_non_zero_row_idx
            # begin backsubstitution.
            def _():
                @for_range(last_pivot_idx, -1, -1)
                def _(i):
                    # compute x[i] = ( y[i] - \sum_{j > i}{ U[i][j] * x[j] } ) / U[i][i]
                    x[i] = y[i]
                    @for_range(i + 1, num_cols)
                    def _(j):
                        x[i] = x[i] - ( U[i][j] * x[j]) 
                    x[i] = x[i].field_div(U[i][i])
        return x

def create_vandermonde_matrix(num_eval_points: int, degree: int, value_type, eval_points=None):
    '''
    Creates a Vandermonde matrix from the given parameters. Recall that a Vandermonde 
    matrix V allows for evaluation of a polynomial p at multiple points via matrix multiplication.
    We choose to support left multiplication: given a Vandermonde matrix V over evaluation points 
    x_1,...,x_n,  and a column matrix P of coefficients c_1,...,c_{d+1} corresponding to a degree 
    d polynomial p, the matrix multiplication V*P yields a size n column matrix with the points
    p(x_1),...,p(x_n).

    :param num_eval_points: Number of evaluation points the Vandermonde matrix supports.
    :param degree: Polynomial degree the Vandermonde matrix supports. 
    :param value_type: Run-time data type of entries (e.g., sint, cint)
    :param eval_pts: Optional Array of explicit evaluation points. Value type must match :param value_type
    - If given, the length of the list must equal :param num_eval_points
    - If eval_pts=None, we default to eval_pts=[1,...,num_eval_points]. Note these integers will be interpreted according to the computation domain.

    :returns: A Vandermonde matrix with :param num_eval_points rows and :param degree + 1 columns.
    '''
    if eval_points is None:
        eval_points = Array(num_eval_points, value_type).assign([i for i in range(1,num_eval_points+1)])
    else:
        assert eval_points.value_type == value_type 
        assert eval_points.length == num_eval_points
    
    V = Matrix(num_eval_points, degree + 1, value_type)
    @for_range(num_eval_points)
    def _(i):
        V[i][0] = value_type(1)
        @for_range(1,degree+1)
        def _(j):
            V[i][j] = V[i][j-1] * eval_points[i]
    # print_ln("V=%s", V.reveal())
    return V

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
        "--lu", 
        dest="lu", 
        action="store_false",
        default=True,
        help="Disable LU decomposition tests"
    )
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
        if compiler.options.lu:
            print_ln("LU TESTS")
            m = [
                [0, 3, 6],
                [1, 4, 7],
                [2, 5, 8]
            ]
            M = sint.Matrix(3,3)
            M.assign(m)
            # print_ln("M=%s\n", M.reveal())

            print_ln("----Test 1----")
            # find a preimage of b under M when we know a preimage (x) exists. 
            # The preimage that is returned may not be x, and that is okay. This just means M is not a bijection.
            x = sint.Array(3)
            x.assign([1, 2, 3])
            b = M.dot(x)
            print_ln("b=%s", b.reveal())
            solver = LUSolver(M)
            preimage = solver.solve(b)
            print_ln("preimage of b=%s under M=%s is %s", b.reveal(), M.reveal(), preimage.reveal())
            print_ln("check against b: M * preimage = %s\n", M.dot(preimage).reveal())

            print_ln("----Test 2----")
            # given a vector, find its preimage under M, if it exists
            b = Array(3, sint).assign([4, 5, 6])
            preimage = solver.solve(b)
            print_ln("preimage of b=%s under M=%s is %s", b.reveal(), M.reveal(), preimage.reveal())
            print_ln("check against b: M * preimage = %s\n", M.dot(preimage).reveal())

            print_ln("----Test 3----")
            # find a random element of the kernel (preimage of 0 under M)
            b = Array(3, sint).assign([0,0,0])
            preimage = solver.solve(b, free_vars="rand")
            print_ln("preimage of b=%s under M=%s is %s", b.reveal(), M.reveal(), preimage.reveal())
            print_ln("check against b: M * preimage = %s\n", M.dot(preimage).reveal())

            print_ln("----Test 4----")
            # buggy case from Shamir test 1
            m = [[1, 1, 1], [1, 2, 4], [1, 3, 9]]
            M = cint.Matrix(3,3)
            M.assign(m)
            solver = LUSolver(M)
            b = Array(3, sint).assign([17, 22, 27])
            preimage = solver.solve(b)
            print_ln("preimage of b=%s under M=%s is %s", b.reveal(), M.reveal(), preimage.reveal())
            print_ln("check against b: M * preimage = %s\n", M.dot(preimage).reveal())

            print_ln("----Test 5----")
            # solving for random dot product - relevant to CKOS22 LRSS impl
            m = [1, 2]
            M = sint.Matrix(1, 2)
            M[0].assign(m)
            solver = LUSolver(M)
            b = Array(1, sint).assign([6])
            preimage = solver.solve(b, free_vars='rand')
            print_ln("preimage of b=%s under M=%s is %s", b.reveal(), M.reveal(), preimage.reveal())
            print_ln("check against b: M * preimage = %s\n", M.dot(preimage).reveal())
        

            print_ln("----Test 6 (Matrix with cgf2n)----")
            m = [1, 2]
            M = cgf2n.Matrix(1,2)
            M[0].assign(m)
            print_ln("%s", M.reveal())

            print_ln("----Test 7 (Matrix with sgf2n)----")
            i = sgf2n(400).get_random_bit()
            l = sgf2n(2).get_random_bit()
            m = [1, 2]
            M = sgf2n.Matrix(1,2)
            M[0].assign(m)
            print_ln("%s, %s, %s", M.reveal(), i.reveal(), l.reveal())


            print_ln("----Test 8 (Randomize with cgf2n)----")
            m = [[1, 1, 1], [1, 2, 4], [1, 3, 9]]
            M = cgf2n.Matrix(3,3)
            M.assign(m)
            solver = LUSolver(M)
            print_ln("%s",M)
            b = Array(3, cgf2n).assign([17, 22, 27])
            preimage = solver.solve(b, free_vars='rand')
            print_ln("%s and %s and %s",solver,b,preimage)
            print_ln("preimage of b=%s under M=%s is %s", b.reveal(), M.reveal(), preimage.reveal())
            print_ln("check against b: M * preimage = %s\n", M.dot(preimage).reveal())



            print_ln("----Test 9 (Randomize with sgf2n)----")
            m = [[31, 17, 771], [123, 839, 401], [165, 332, 912]]
            M = sgf2n.Matrix(3,3)
            M.assign(m)
            solver = LUSolver(M)
            b = Array(3, sgf2n).assign([17, 22, 27])
            preimage = solver.solve(b, free_vars='rand')
            print_ln("preimage of b=%s under M=%s is %s", b.reveal(), M.reveal(), preimage.reveal())
            print_ln("check against b: M * preimage = %s\n", M.dot(preimage).reveal())

            print_ln("----Test 10----")
            # given a vector, find its preimage under M, if it exists
            b = Array(3, sgf2n).assign([4, 5, 6])
            preimage = solver.solve(b)
            print_ln("preimage of b=%s under M=%s is %s", b.reveal(), M.reveal(), preimage.reveal())
            print_ln("check against b: M * preimage = %s\n", M.dot(preimage).reveal())
            # if M.dot(preimage).reveal() == b.reveal():
            #     print_ln("Test Passed!")
            # else:
            #     print_ln("Test Failed")

            print_ln("----Test 8 (Array Test)----")
            # solving for random dot product - relevant to CKOS22 LRSS impl
            m = [27, 49]
            M = sint.Matrix(1, 2)
            M[0].assign(m)
            solver = LUSolver(M)
            b = Array(1, sint).assign([6])
            preimage = solver.solve(b, free_vars='rand')
            print_ln("preimage of b=%s under M=%s is %s", b.reveal(), M.reveal(), preimage.reveal())
            print_ln("check against b: M * preimage = %s\n", M.dot(preimage).reveal())


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