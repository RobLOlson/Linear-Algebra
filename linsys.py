from decimal import Decimal, getcontext
from copy import deepcopy

from vector import Vector
from plane import Plane
from math import log10

import pdb

getcontext().prec = 30

class LinearSystem(object):
    """public attributes:
    dimension
    planes[]

    TODO
    extract solution from rref
    """

    ALL_PLANES_MUST_BE_IN_SAME_DIM_MSG = 'All planes in the system should live in the same dimension'
    NO_SOLUTIONS_MSG = 'No solutions'
    INF_SOLUTIONS_MSG = 'Infinitely many solutions'

    def __init__(self, planes, delta='0.001'):
        self.delta = Decimal(delta)
        try:
            d = planes[0].dimension
            for p in planes:
                assert p.dimension == d

            self.planes = planes
            self.dimension = d

        except AssertionError:
            raise Exception(self.ALL_PLANES_MUST_BE_IN_SAME_DIM_MSG)

    @property
    def prec(self):
        return round(log10(1/self.delta))

    def swap_rows(self, row1, row2):
        self.planes[row1], self.planes[row2] = self.planes[row2], self.planes[row1]

    def scalar_multiply(self, coefficient, row):
        self.planes[row] = Decimal(coefficient) * self.planes[row]

    def linear_combine(self, amount, row_to_add, destination):
        a = destination
        b = row_to_add
        self.planes[a] += Decimal(amount)*self.planes[b]

    def solve_system(self):
        system = deepcopy(self)
        system = system.compute_rref()


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

    @staticmethod
    def pivot_num(self):
        pivots = self.pivots()
        empty_rows = [e for e in pivots if e[0] == None]

        return len(pivots) - len(empty_rows)

    def show_pivots(self):
        system = deepcopy(self)
        pivots = system.pivots()
        pivot_elements = []
        for pln_i, pln in enumerate(pivots):
            if pln[1] != None:
                pivot_elements.append((pln_i, pln[1]))
        return system.matrix_form(highlights=pivot_elements)

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
        # pivots = [e for e in pivots if e]
        print("pivots = {}".format(pivots))
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
                s2 = " [x{}]".format(i+1)
            else:
                s2 = "     "
            s3 = " || ["+" {:^"+str(self.simple_num(Decimal(round(largest_const, self.prec))))+"} "+"]"
            s3 = s3.format(self.planes[i].constant_term)

            result += s.format(width=largest_coef, *pln.normal_vector)+s2+s3
        return result
            # print(s.format(width=largest_coef,*pln.normal_vector)+s2+s3)

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
    # breakpoint()
    r7 = s7.compute_rref()
    # breakpoint()
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

