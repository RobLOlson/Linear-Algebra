from decimal import Decimal, getcontext
from copy import deepcopy

from vector import Vector
from plane import Plane
from math import log10

import pdb

getcontext().prec = 30

class Paramaterization(object):
    """"""
    def __init__(self, basepoint, vector_space):
        self.basepoint = Vector(basepoint)
        self.vector_space = vector_space

class VectorSpace(object):
    """For a given possibly redundant basis,
    generates a spanned VectorSpace."""

    def __init__(self, basis_vectors, delta=Decimal('.0000001')):
        self.basis = [Vector(e, delta) for e in basis_vectors]
        self.originals = [Vector(e, delta) for e in basis_vectors]
        self.delta = delta

    @property
    def dimension(self):
        """The dimension of a VECTOR SPACE is the minimal number of
        basis vectors IN that space, **NOT the dimension of the vectors**."""
        system = deepcopy(self)
        system.orthonormalize()
        return len(system.basis)

    @property
    def prec(self):
        return round(log10(1/self.delta))-1

    def __setitem__(self, index, val):
        self.basis[index] = Vector(val, self.delta)

    def __getitem__(self, index):
        return self.basis[index]

    def __repr__(self):
        return f"{self.basis}"

    def orthonormalize(self):
        system = deepcopy(self.basis)

        for ti, target in enumerate(system):
            previous_targets = deepcopy(self[:ti])
            for prev in previous_targets:
                if abs(prev):
                    unit_prev = prev * (1/abs(prev))
                    amount = (target*unit_prev)
                    target = target - amount * unit_prev
            if target:
                unit_target = target * (1/abs(target))
            self[ti] = unit_target

        return self

# End of class VectorSpace(...) ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

class LinearSystem(object):
    """public attributes:
    dimension
    planes[]

    TODO
    paramaterization
    """

    ALL_PLANES_MUST_BE_IN_SAME_DIM_MSG = 'All planes in the system should live in the same dimension'
    NO_SOLUTIONS_MSG = 'No solutions'
    INF_SOLUTIONS_MSG = 'Infinitely many solutions'

    def __init__(self, planes, delta=Decimal('0.0001')):
        self.delta = Decimal(delta)
        self.permutations = []
        try:
            d = planes[0].dimension
            for p in planes:
                assert p.dimension == d

            self.planes = planes
            self.dimension = d

        except AssertionError:
            raise Exception(self.ALL_PLANES_MUST_BE_IN_SAME_DIM_MSG)

    @classmethod
    def from_column_space(cls, vector_space, target_vector, delta=Decimal('.0001')):
        """return My_Class(params...)"""

        planes = []
        target_vector = Vector(target_vector)
        assert vector_space[0].dimension == target_vector.dimension

        for pln_i in range(len(vector_space[0].dimension)):
            planes.append(Plane([pln[pln_i] for pln in vector_space], target_vector[pln_i]))

        return LinearSystem(planes, delta)

    @classmethod
    def diagonal(cls, diagonals=Vector([1, 1, 1])):
        """returns a square LinearSystem() with only diagonal entries non-zero"""

        slices = []
        for i, val in enumerate(diagonals):
            normal = Vector([0]*len(diagonals))
            normal[i] = val
            slices.append(Plane(normal_vector=normal, constant_term=0), delta=self.delta)

        return LinearSystem(slices)

    @property
    def prec(self):
        return round(log10(1/self.delta))-1

    def swap_rows(self, row1, row2):
        self.planes[row1], self.planes[row2] = self.planes[row2], self.planes[row1]

    def scalar_multiply(self, coefficient, row):
        self.planes[row] = Decimal(coefficient) * self.planes[row]

    def linear_combine(self, amount, row_to_add, destination):
        a = destination
        b = row_to_add
        self.planes[a] += Decimal(amount)*self.planes[b]

    def pivots_by_column(self):
        """For each column, if it has a pivot, add the pivot's row
        number to the list (-1 if no pivot in column)."""

        system = deepcopy(self)
        system = system.rref()
        row_pivots = system.pivots_by_row()
        pivots = [-1]*system.planes[0].dimension

        # for each column
        for n in range(system.planes[0].dimension):
            # if that column is a pivot column
            if n in row_pivots:
                # add the pivots row index to the list
                pivots[n] = row_pivots.index(n)
            else:
                pivots[n] = -1

        return pivots

    def pivots_by_row(self):
        """After gaussian elimination, find the 1st non-zero entry
        in each row and put it in a list.  -1 for a 0-row."""

        system = deepcopy(self)
        system = system.rref()
        pivots = system.indices_of_first_nonzero_terms_in_each_row()
        for a, b in reversed(system.permutations):
            pivots[a], pivots[b] = pivots[b], pivots[a]

        return pivots

    def null_space(self, orthonormal=False):
        """Returns a VectorSpace containing all the vectors for which
        the system produces 0."""

        system = deepcopy(self)
        system = system.rref()
        pivots = system.pivots_by_column()

        basis_vectors = []

        # with the matrix in RREF form
        # for each column (pivots = system.pivots_by_column)
        for coef_i, pln_i in enumerate(pivots):
            # if free column
            if pln_i == -1:
                # there is a null vector corresponding to that free column
                null_vector = [0]*system.planes[0].dimension

                # The null vector has elements that come from identity entries
                # and from the entries of the free column associated with it.
                # The pattern of identity vs. free entries in the null vector
                # is given by the inverse of that pattern in the columns of rref().

                # I.E., if the nth column is free, then the nth entry of the null
                # vector is an identity element; and if the nth column is an identity
                # column, then the nth entry will be (the * inverse of) a free entry.

                # We initialize the null vector with 0's on the assumption that
                # there are no pivot (identity) columns.  Then we only have to
                # fix entries corresponding to pivot columns, plus put a 1 in
                # entry corresponding to the free column we started this with.

                # for each column in the matrix
                for i, e in enumerate(pivots):
                    # if it's the i'th pivot column
                    if e > -1:
                        # change the i'th in the null vector to the i'th
                        # entry of the free column * -1
                        null_vector[i] = -1*system[e][coef_i]

                # the null vector corresponding to the i'th free column has
                # a 1 in its i'th entry
                null_vector[coef_i] = Decimal(1)

                basis_vectors.append(null_vector)

        if orthonormal:
            return VectorSpace(basis_vectors).orthonormalize()
        return VectorSpace(basis_vectors)
    nullspace = null_space

# End of null_space(...) ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    def column_space(self, orthonormal=False):
        """Returs a VectorSpace containing all possible linear
        combinations of the system's column vectors."""

        system = deepcopy(self)
        system = system.rref()
        pivots = system.pivots_by_row()

        basis_vectors = []

        # for a matrix in RREF
        # for each row
        for pln_i, coef_i in enumerate(pivots):
            piv_column = []
            # if it's a pivot row
            if coef_i > -1:
                # go thru that pivot's column
                for j in range(len(system.planes)):
                    piv_column.append(self[j][coef_i])

                basis_vectors.append(piv_column)

        if orthonormal:
            return VectorSpace(basis_vectors).orthonormalize()
        return VectorSpace(basis_vectors)
    columnspace = column_space

# End of column_space(...) ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    def row_space(self, orthonormal=False):
        """Returns a VectorSpace containing all possible linear
        combinations of the system's row (normal) vectors."""
        system = deepcopy(self)
        system = system.rref()
        pivots = system.indices_of_first_nonzero_terms_in_each_row()
        for a, b in reversed(system.permutations):
            pivots[a], pivots[b] = pivots[b], pivots[a]

        basis_vectors = []

        # for each row
        for pln_i, coef_i in enumerate(pivots):
            # if it's a pivot row
            if coef_i > -1:
                # then it's in the row space
                basis_vectors.append(self[pln_i].normal_vector)

        if orthonormal:
            return VectorSpace(basis_vectors).orthonormalize()
        return VectorSpace(basis_vectors)
    rowspace = row_space

# End of row_space(...) ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    def solution(self):
        system = deepcopy(self)
        system = system.compute_rref()

        r = system.pivot_num()
        m = len(system.planes)
        n = system.dimension

        solution = Vector.of_length(n)
        one_solution = True
        many_solutions = False
        no_solution = False

        non_zeroes = system.indices_of_first_nonzero_terms_in_each_row()

        # for a matrix in RREF
        # for each row (remembering that pivot rows are bubbled to the top)
        for pln_i, coef_i in enumerate(non_zeroes):
            # if it's a pivot row
            if pln_i < r:
                # then the right-hand side is part of the solution
                solution[coef_i] = system[pln_i].constant_term

            # if a row has all 0's
            if coef_i == -1:
                # then solution only possible if rhs is also 0
                if abs(system[pln_i].constant_term) > self.delta:
                    no_solution = True
                    break

        if no_solution:
            return False

        # non-empty nullspace
        if n > r:
            return Paramaterization(solution, system.null_space())

        return Vector(solution)

# end of LinearSystem.solution() ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    def make_zeroes_above(self, pln_i, coef_i):
        system = deepcopy(self)
        coefs = system.coefs_above(pln_i, coef_i)
        pivot = system[pln_i][coef_i]
        if pivot == 0:
            raise Exception("Pivot must be non-zero to eliminate above.")
        top = len(coefs)
        for i, target in enumerate(coefs):
            if abs(target) > .0001:
                ratio = Decimal(-target / pivot)
                system.linear_combine(ratio, pln_i, i)

        return system

    def compute_rref(self):
        system = deepcopy(self)
        system = system.triangular()

        for i in range(len(system.planes)):
            piv_columns = system.indices_of_first_nonzero_terms_in_each_row()
            if piv_columns[i] > -1:
                system = system.make_zeroes_above(i, piv_columns[i])
                if system[i][piv_columns[i]] != 1:
                    system[i] = system[i] * (1/system[i][piv_columns[i]])
                    # system.scalar_multiply(1/piv_columns[i], i)
        return system

    rref = compute_rref

    def pivot_num(self):
        pivots = self.pivots_by_row()
        return len([e for e in pivots if e > -1])

    def show_pivots(self):
        system = deepcopy(self)
        system = system.rref()
        # pivots = system.pivots()
        pivots = system.indices_of_first_nonzero_terms_in_each_row()
        for a, b in reversed(system.permutations):
            pivots[a], pivots[b] = pivots[b], pivots[a]

        pivot_elements = []
        for pln_i, coef_i in enumerate(pivots):
            if coef_i > -1:
                pivot_elements.append((pln_i, coef_i))
        return self.matrix_form(highlights=pivot_elements)

    def pivots(self):
        """Returns a list of (pivot_value, pivot_index) tuples.
        Pivots are ordered from top to bottom.
        Rows with no pivot return (None, None)."""

        pivots = [(None, None) for i in range(len(self.planes))]

        #these two loops scan through the matrix left->right, top->bottom
        #when they find a non-zero matrix entry, it becomes a pivot
        #entries in that row can no longer be pivots
        for pln_i, pln in enumerate(self.planes):
            for coef_i, coef in enumerate(pln.normal_vector):
                if coef != Decimal(0) and coef_i not in [e[1] for e in pivots]:
                    pivots[pln_i] = (coef, coef_i)

                    break
        return pivots

    def coefs_below(self, pln_i, coef_i):
        coefs = []
        pln_i += 1

        while pln_i < len(self.planes):
            coefs.append(self[pln_i][coef_i])
            pln_i += 1

        return coefs

    def coefs_above(self, pln_i, coef_i):
        coefs = []
        counter = 0
        while counter < pln_i:
            coefs.append(self[counter][coef_i])
            counter+=1

        return coefs

    def matrix_form(self, highlights=[]):
        largest_coef = 1
        largest_const = 1
        for pln in self.planes:
            if len(f'{self.simple_num(Decimal(round(pln.constant_term, self.prec)))}') > largest_const:
                largest_const = len(f'{self.simple_num(Decimal(round(pln.constant_term, self.prec)))}')

            for coef in pln.normal_vector:
                if len(f'{self.simple_num(Decimal(round(coef, self.prec)))}') > largest_coef:
                    largest_coef = len(f'{self.simple_num(Decimal(round(coef, self.prec)))}')
        result = ""
        # if largest_coef > 6:
        #     largest_coef = 6
        extra = 0
        s0= ""

        if self.dimension > len(self.planes):
            extra = self.dimension - len(self.planes)
            for e in range(extra):
                s0 += "\n"+("  "+" "*largest_coef+"  ")*self.dimension+f" [x{e+1}]"

        for i, pln in enumerate(self.planes):
            s = ""
            for j, coef in enumerate(pln.normal_vector):
                if not j:
                    s+= "\n"
                s += "["
                if (i, j) in highlights:
                    s += "*"
                else:
                    s += " "
                s += "{val:^{width}}".format(width=largest_coef, val=self.simple_num(Decimal(round(coef, self.prec))))
                if (i, j) in highlights:
                    s += "*"
                else:
                    s += " "
                s += "]"
            if i < self.dimension:
                s2 = " [x{}]".format(i+1+extra)
            else:
                s2 = "     "
            s3 = " || ["+" {:^"+str(self.simple_num(Decimal(round(largest_const, self.prec))))+"} "+"]"
            s3 = s3.format(self.planes[i].constant_term)

            result += s.format(width=largest_coef, *pln.normal_vector)+s2+s3
        result = s0 + result
        return result

# end of LinearSystem.matrix_form(...) ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    def index_of_first_non_zero_below(self, pln_i, coef_i):
        below = self.coefs_below(pln_i, coef_i)
        for i, e in enumerate(below):
            if e:
                return i

    def triangular(self):
        """Implement the 1st stage of gaussian elimination.
        Returns the matrix after it has been upper-triangularized."""

        system = deepcopy(self)

        # scan matrix entries left->right ; top->bottom
        # pln_i = m
        # coef_i = n
        pln_i = 0
        skipped_columns = 0
        while pln_i < len(system.planes):
            coef_i = pln_i + skipped_columns

            # find a non-zero column ; skip empty columns
            finding_column = True
            while coef_i < system.dimension and finding_column:
                #if current column not empty
                if(any(system.coefs_below(pln_i-1, coef_i))):
                    finding_column = False
                    piv_i = pln_i + system.index_of_first_non_zero_below(pln_i-1, coef_i)
                    if piv_i != pln_i:
                        system.swap_rows(piv_i, pln_i)
                        system.permutations.append((piv_i, pln_i))
                        piv_i = pln_i
                    piv_val = system.planes[pln_i][coef_i]
                    below = system.coefs_below(pln_i, coef_i)
                    for below_i, below_co in enumerate(below):
                        system.linear_combine(-below_co/piv_val,
                                              piv_i,
                                              pln_i + below_i + 1)
                #current column was empty, so check next
                else:
                    coef_i += 1
                    skipped_columns += 1
            pln_i += 1

        return system

# end of LinearSystem.triangular(...) ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    def indices_of_first_nonzero_terms_in_each_row(self):
        num_equations = len(self)
        num_variables = self.dimension

        indices = [-1] * num_equations

        for i, p in enumerate(self.planes):
            try:
                indices[i] = p.first_nonzero_index(p.normal_vector)
            except Exception as e:
                if str(e) == Plane.NO_NONZERO_ELTS_FOUND_MSG:
                    continue
                else:
                    raise e

        return indices

    @staticmethod
    def simple_num(num):
        return num.to_integral() if num == num.to_integral() else num.normalize()

    def __len__(self):
        return len(self.planes)

    def __delitem__(self, target):
        del self.planes[target]

    def __add__(self, target):
        new_system = self.planes
        if target.dimension == self.dimension:
            if type(target) is Plane:
                new_system.append(target)
            elif type(target) is LinearSystem:
                new_system += target.planes

            return LinearSystem(new_system)
        else:

            raise Exception("Equations must have same dimension.")

    def __getitem__(self, i):
        return self.planes[i]


    def __setitem__(self, i, x):
        try:
            assert x.dimension == self.dimension
            self.planes[i] = x

        except AssertionError:
            raise Exception(self.ALL_PLANES_MUST_BE_IN_SAME_DIM_MSG)


    def __str__(self):
        return self.show_pivots()
        ret = 'Linear System:\n'
        temp = ['Equation {}: {}'.format(i+1,p) for i,p in enumerate(self.planes)]
        ret += '\n'.join(temp)
        return ret

    __repr__ = __str__
# end of class LinearSystem(...)

class MyDecimal(Decimal):
    def is_near_zero(self, eps=1e-10):
        return abs(self) < eps

p0 = Plane(normal_vector=Vector(['1','1','1']), constant_term='1')
p1 = Plane(normal_vector=Vector(['0','1','0']), constant_term='2')
p2 = Plane(normal_vector=Vector(['1','1','-1']), constant_term='3')
p3 = Plane(normal_vector=Vector(['1','0','-2']), constant_term='2')

p00 = Plane([0,0,0], 0)

s = LinearSystem([p0,p1,p2,p3])
m = LinearSystem([deepcopy(p00), deepcopy(p00), deepcopy(p00)])

if __name__ == "__main__":
    s.swap_rows(0,1)
    if not (s[0] == p1 and s[1] == p0 and s[2] == p2 and s[3] == p3):
        print('test case 1 failed')

    s.swap_rows(1,3)
    if not (s[0] == p1 and s[1] == p3 and s[2] == p2 and s[3] == p0):
        print('test case 2 failed')

    s.swap_rows(3,1)
    if not (s[0] == p1 and s[1] == p0 and s[2] == p2 and s[3] == p3):
        print('test case 3 failed')

    s.scalar_multiply(1,0)
    if not (s[0] == p1 and s[1] == p0 and s[2] == p2 and s[3] == p3):
        print('test case 4 failed')

    s.scalar_multiply(-1,2)
    if not (s[0] == p1 and
            s[1] == p0 and
            s[2] == Plane(normal_vector=Vector(['-1','-1','1']), constant_term='-3') and
            s[3] == p3):
        print('test case 5 failed')
    s.scalar_multiply(10,1)
    if not (s[0] == p1 and
            s[1] == Plane(normal_vector=Vector(['10','10','10']), constant_term='10') and
            s[2] == Plane(normal_vector=Vector(['-1','-1','1']), constant_term='-3') and
            s[3] == p3):
        print('test case 6 failed')

    s.linear_combine(0,0,1)
    if not (s[0] == p1 and
            s[1] == Plane(normal_vector=Vector(['10','10','10']), constant_term='10') and
            s[2] == Plane(normal_vector=Vector(['-1','-1','1']), constant_term='-3') and
            s[3] == p3):
        print('test case 7 failed')

    s.linear_combine(1,0,1)
    if not (s[0] == p1 and
            s[1] == Plane(normal_vector=Vector(['10','11','10']), constant_term='12') and
            s[2] == Plane(normal_vector=Vector(['-1','-1','1']), constant_term='-3') and
            s[3] == p3):
        print('test case 8 failed')

    s.linear_combine(-1,1,0)
    if not (s[0] == Plane(normal_vector=Vector(['-10','-10','-10']), constant_term='-10') and
            s[1] == Plane(normal_vector=Vector(['10','11','10']), constant_term='12') and
            s[2] == Plane(normal_vector=Vector(['-1','-1','1']), constant_term='-3') and
            s[3] == p3):
        print('test case 9 failed')


    p1 = Plane(normal_vector=Vector(['1','1','1']), constant_term='1')
    p2 = Plane(normal_vector=Vector(['0','1','1']), constant_term='2')
    s1 = LinearSystem([p1,p2])
    t1 = s1.triangular()
    if not (t1[0] == p1 and
            t1[1] == p2):
        print('test case 2.1 failed')

    p1 = Plane(normal_vector=Vector(['1','1','1']), constant_term='1')
    p2 = Plane(normal_vector=Vector(['1','1','1']), constant_term='2')
    s2 = LinearSystem([p1,p2])
    t2 = s2.triangular()
    if not (t2[0] == p1 and
            t2[1] == Plane(constant_term='1')):
        print('test case 2.2 failed')

    p1 = Plane(normal_vector=Vector(['1','1','1']), constant_term='1')
    p2 = Plane(normal_vector=Vector(['0','1','0']), constant_term='2')
    p3 = Plane(normal_vector=Vector(['1','1','-1']), constant_term='3')
    p4 = Plane(normal_vector=Vector(['1','0','-2']), constant_term='2')
    s3 = LinearSystem([p1,p2,p3,p4])
    t3 = s3.triangular()
    if not (t3[0] == p1 and
            t3[1] == p2 and
            t3[2] == Plane(normal_vector=Vector(['0','0','-2']), constant_term='2') and
            t3[3] == Plane()):
        print('test case 2.3 failed')

    p1 = Plane(normal_vector=Vector(['0','1','1']), constant_term='1')
    p2 = Plane(normal_vector=Vector(['1','-1','1']), constant_term='2')
    p3 = Plane(normal_vector=Vector(['1','2','-5']), constant_term='3')
    s4 = LinearSystem([p1,p2,p3])
    t4 = s4.triangular()
    if not (t4[0] == Plane(normal_vector=Vector(['1','-1','1']), constant_term='2') and
            t4[1] == Plane(normal_vector=Vector(['0','1','1']), constant_term='1') and
            t4[2] == Plane(normal_vector=Vector(['0','0','-9']), constant_term='-2')):
        print('test case 2.4 failed')

    p1 = Plane(normal_vector=Vector(['1','1','1']), constant_term='1')
    p2 = Plane(normal_vector=Vector(['0','1','1']), constant_term='2')
    s5 = LinearSystem([p1,p2])
    r5 = s5.compute_rref()
    if not (r5[0] == Plane(normal_vector=Vector(['1','0','0']), constant_term='-1') and
            r5[1] == p2):
        print('test case 2.1 failed')

    p1 = Plane(normal_vector=Vector(['1','1','1']), constant_term='1')
    p2 = Plane(normal_vector=Vector(['1','1','1']), constant_term='2')
    s6 = LinearSystem([p1,p2])
    r6 = s6.compute_rref()
    if not (r6[0] == p1 and
            r6[1] == Plane(constant_term='1')):
        print('test case 2.2 failed')

    p1 = Plane(normal_vector=Vector(['1','1','1']), constant_term='1')
    p2 = Plane(normal_vector=Vector(['0','1','0']), constant_term='2')
    p3 = Plane(normal_vector=Vector(['1','1','-1']), constant_term='3')
    p4 = Plane(normal_vector=Vector(['1','0','-2']), constant_term='2')
    s7 = LinearSystem([p1,p2,p3,p4])
    r7 = s7.compute_rref()
    if not (r7[0] == Plane(normal_vector=Vector(['1','0','0']), constant_term='0') and
            r7[1] == p2 and
            r7[2] == Plane(normal_vector=Vector(['0','0','-2']), constant_term='2') and
            r7[3] == Plane()):
        print('test case 2.3 failed')

    p1 = Plane(normal_vector=Vector(['0','1','1']), constant_term='1')
    p2 = Plane(normal_vector=Vector(['1','-1','1']), constant_term='2')
    p3 = Plane(normal_vector=Vector(['1','2','-5']), constant_term='3')
    s8 = LinearSystem([p1, p2, p3])
    r8 = s8.compute_rref()
    if not (r8[0] == Plane(normal_vector=Vector(['1','0','0']), constant_term=Decimal('23')/Decimal('9')) and
            r8[1] == Plane(normal_vector=Vector(['0','1','0']), constant_term=Decimal('7')/Decimal('9')) and
            r8[2] == Plane(normal_vector=Vector(['0','0','1']), constant_term=Decimal('2')/Decimal('9'))):
        print('test case 2.4 failed')

    p1 = Plane([5.862, 1.178, -10.366], -8.15)
    p2 = Plane([-2.931, -.589, 5.183], -4.075)
    s9 = LinearSystem([p1, p2])
    sol9 = s9.solution()
    if sol9:
        print('test case 3.1 failed')

    p1 = Plane([8.631, 5.112, -1.816], -5.113)
    p2 = Plane([4.315, 11.132, -5.27], -6.775)
    p3 = Plane([-2.158, 3.01, -1.727], -.831)
    s10 = LinearSystem([p1, p2, p3])
    sol10 = s10.solution()

    if sol10:
        print('test case 3.2 failed')

    p1 = Plane([5.262, 2.739, -9.878], -3.441)
    p2 = Plane([5.111, 6.358, 7.638], -2.152)
    p3 = Plane([2.016, -9.924, -1.367], -9.278)
    p4 = Plane([2.167, -13.543, -18.883], -10.567)
    s11 = LinearSystem([p1, p2, p3, p4])
    sol11 = s11.solution()
    if sol11 != Vector([-1.1772, 0.7071, -0.0827]):
        print('test case 3.3 failed')

    p1 = Plane([.786, .786, .588], -.714)
    p2 = Plane([-.138, -.138, .244], .319)
    s12 = LinearSystem([p1, p2])
    sol12 = s12.solution()

    p1 = Plane([8.631, 5.112, -1.816], -5.113)
    p2 = Plane([4.315, 11.132, -5.27], -6.775)
    p3 = Plane([-2.158, 3.01, -1.727], -0.831)
    s13 = LinearSystem([p1, p2, p3])

    p1 = Plane([0.935, 1.76, -9.365], -9.955)
    p2 = Plane([0.187, 0.352, -1.873], -1.991)
    p3 = Plane([.374, .704, -3.746], -3.982)
    p4 = Plane([-.561, -1.056, 5.619], 5.973)
    s14 = LinearSystem([p1, p2, p3, p4])

    vs = VectorSpace([[1, 1, 1], [1, 2, 3]])
    vs1 = VectorSpace([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    vs2 = VectorSpace([[1, 10, 100], [2, 4, 6], [3, 9, 27]])
    vs3 = VectorSpace([[1, 1, 1, 1], [1, 2, 3, 4], [1, 4, 9, 16], [1, 8, 27, 64]])
