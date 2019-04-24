from decimal import Decimal, getcontext

from vector import Vector

import unittest, random

getcontext().prec = 5

class Line(object):
    """Defined by a normal vector and a constant term.
    If ax + by = c, then

    the normal vector is [a, b] and the constant term is c.

    TO IMPLEMENT:

    *Determine if two lines are parallel
    *Determine if two lines are equal
    *Compute the intersection of two lines"""
    NO_NONZERO_ELTS_FOUND_MSG = 'No nonzero elements found'


    def __init__(self, normal_vector=None, constant_term=None):
        """
        public facing instance attributes
        ---------------------------------
        self.normal_vector
        self.basepoint"""

        self.dimension = 2

# if no normal_vector is given, make it a 0 vector
        if not normal_vector:
            all_zeros = ['0']*self.dimension
            normal_vector = Vector(all_zeros)
        self.normal_vector = Vector(normal_vector)

        if not constant_term:
            constant_term = Decimal('0')
        self.constant_term = Decimal(constant_term)

# basepoint is any point on the line, or in this
# case, a point at which the line intersects
# either the x or y axis
        self.set_basepoint()
    #end of __init__

    def set_basepoint(self, vec=None):
        if vec:
            self.basepoint = Vector(vec)

# If we know a point on the line (basepoint) and
# we know the normal coefficients, we can calculate
# the constant term by simply plugging the paratmeters
# Ax + By = c (where x and y are given by the basep)
            x = self.normal_vector[0]*self.basepoint[0]
            y = self.normal_vector[1]*self.basepoint[1]

            self.constant_term = x+y

        try:
            n = self.normal_vector
            c = self.constant_term
            basepoint_coords = ['0']*self.dimension

            initial_index = Line.first_nonzero_index(n)
            initial_coefficient = n[initial_index]
            basepoint_coords[initial_index] = c/Decimal(initial_coefficient)
            self.basepoint = Vector(basepoint_coords)

        except Exception as e:
            if str(e) == Line.NO_NONZERO_ELTS_FOUND_MSG:
                self.basepoint = None
            else:
                raise e
    #end of set_basepoint(self, vec=None)

    def __getitem__(self, index):
        return self.normal_vector[index]

    def __setitem__(self, index, val):
        temp = list(self.normal_vector)
        temp[index] = val
        self.normal_vector = Vector(temp)

    def __str__(self):

        num_decimal_places = 3

        def write_coefficient(coefficient, is_initial_term=False):
            coefficient = round(coefficient, num_decimal_places)
            if coefficient % 1 == 0:
                coefficient = int(coefficient)

            output = ''

            if coefficient < 0:
                output += '-'
            if coefficient > 0 and not is_initial_term:
                output += '+'

            if not is_initial_term:
                output += ' '

            if abs(coefficient) != 1:
                output += '{}'.format(abs(coefficient))

            return output

        n = self.normal_vector

        try:
            initial_index = Line.first_nonzero_index(n)
            terms = [write_coefficient(n[i], is_initial_term=(i==initial_index)) + 'x_{}'.format(i+1)
                     for i in range(self.dimension) if round(n[i], num_decimal_places) != 0]
            output = ' '.join(terms)

        except Exception as e:
            if str(e) == self.NO_NONZERO_ELTS_FOUND_MSG:
                output = '0'
            else:
                raise e

        constant = round(self.constant_term, num_decimal_places)
        if constant % 1 == 0:
            constant = int(constant)
        output += ' = {}'.format(constant)

        return output
    # end of __str__

    __repr__ = __str__

    def __eq__(self, target):
        if self.normal_vector.is_zero():
            if target.normal_vector.is_zero():
                return abs(self.constant_term - target.constant_term) < 0.001
            else:
                return False
        elif target.normal_vector.is_zero():
            return False

        return self.same_line_Q(target)

    def direct_vector(self):
        return Vector([-self[1], self[0]])

    def parallel_Q(self, geometry):
        if type(geometry) is Line:
            l = geometry
            return self.normal_vector.parallel_Q(l.normal_vector)

        else:
            v = Vector(geometry)
            return self.normal_vector.parallel_Q(v)

    def same_line_Q(self, l):
        if type(l) is Line:
            if self.parallel_Q(l):
                return (abs(self[0]/l[0] - self.constant_term/l.constant_term) < .001)
            else:
                return False
        else:
            raise Exception("Cannot compare lines to other objects.")

    def random_point(self, interval=1):
        p = self.direct_vector()*Decimal(random.random()*interval)
        p += self.basepoint

        return p

    def intersection(self, geometry):

        if type(geometry) is Line:
            l = geometry

            print("TEST\n")
            print(self.parallel_Q(l))
            print(self.same_line_Q(l))
            if self.parallel_Q(l):
                if self.same_line_Q(l):
                    return self
                else:
                    return False

            A = Decimal(self.normal_vector[0])
            B = Decimal(self.normal_vector[1])
            C = Decimal(l.normal_vector[0])
            D = Decimal(l.normal_vector[1])
            k1 = Decimal(self.constant_term)
            k2 = Decimal(l.constant_term)

            x = (D*k1 - B*k2) / (A*D - B*C)
            y = (-C*k1 + A*k2) / (A*D - B*C)

            return Vector([x, y])

        #else geometry NOT a line
        geometry = Vector(geometry)

        test_point = geometry
        #move test_point to the origin
        test_point -= self.basepoint
        return test_point.orthogonal_Q(self.normal_vector)

    @staticmethod
    def first_nonzero_index(iterable):
        for k, item in enumerate(iterable):
            if not MyDecimal(item).is_near_zero():
                return k
        raise Exception(Line.NO_NONZERO_ELTS_FOUND_MSG)


class MyDecimal(Decimal):
    def is_near_zero(self, eps=1e-10):
        return abs(self) < eps


class TestLineMethods(unittest.TestCase):
    def setUp(self):
        self.A = Line([4.046, 2.836], 1.21)
        self.B = Line([10.115, 7.09], 3.025)
        self.C = Line([7.204, 3.182], 8.68)
        self.D = Line([8.172, 4.114], 9.883)
        self.E = Line([1.182, 5.562], 6.744)
        self.F = Line([1.773, 8.343], 9.525)

    def test_parallel_Q(self):
        self.assertEqual(True,
                         self.A.parallel_Q(self.B))

        self.assertEqual(False,
                         self.C.parallel_Q(self.D))

        self.assertEqual(True,
                         self.E.parallel_Q(self.F))

    def test_same_line_Q(self):
        self.assertTrue(self.A.same_line_Q(self.B))
        self.assertFalse(self.C.same_line_Q(self.D))
        self.assertFalse(self.E.same_line_Q(self.F))

    def test_random_point(self):
        self.assertTrue(self.A.intersection(self.A.random_point()))

    def test_intersection(self):
        self.assertEqual(self.A,
                         self.A.intersection(self.B))

        self.assertEqual(Vector([1.173, 0.073]),
                         self.C.intersection(self.D))

        self.assertEqual(False,
                         self.E.intersection(self.F))


A = Line([4.046, 2.836], 1.21)
B = Line([10.115, 7.09], 3.025)
C = Line([7.204, 3.182], 8.68)
D = Line([8.172, 4.114], 9.883)
E = Line([1.182, 5.562], 6.744)
F = Line([1.773, 8.343], 9.525)


a = Line([1, 1], 1)
b = Line([2, 2], 1)
c = Line(normal_vector=[4, 5], constant_term=6)
d = Line([1, 2], 3)

if __name__ == '__main__':
    unittest.main()
