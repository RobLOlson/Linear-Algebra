from decimal import Decimal, getcontext

from vector import Vector

import unittest, random

getcontext().prec = 30


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

        #if no normal_vector is given, make it a 0 vector
        if not normal_vector:
            all_zeros = ['0']*self.dimension
            normal_vector = Vector(all_zeros)
        self.normal_vector = Vector(normal_vector)

        if not constant_term:
            constant_term = Decimal('0')
        self.constant_term = Decimal(constant_term)

        #basepoint is any point on the line, or in this
        #case, a point at which the line intersects
        #either the x or y axis
        self.set_basepoint()


    def set_basepoint(self, vec=None):
        if vec:
            self.basepoint = Vector(vec)

            #the basepoint vector must have a component that
            #is perpendicular to the line; the length of that
            #component must be the constant term (when the line
            #passes through 0, the constant term is 0)
            #zero = self.basepoint.projected_onto(self.normal_vector)
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

    def parallel_Q(self, l):
        return self.normal_vector.parallel_Q(l.normal_vector)

    def same_line_Q(self, l):
        return self.random_point().parallel_Q(l.random_point())

    def random_point(self, interval=1):
        direct_vector = Vector([self.normal_vector[1],
                         -self.normal_vector[0]])

        direct_vector *= random.random()*interval
        direct_vector += self.basepoint

        return direct_vector

    def intersection(self, geometry):


        if type(geometry) is Line:
            l = geometry
            A = Decimal(self.normal_vector[0])
            B = Decimal(self.normal_vector[1])
            C = Decimal(l.normal_vector[0])
            D = Decimal(l.normal_vector[1])
            k1 = Decimal(self.constant_term)
            k2 = Decimal(l.constant_term)

            x = (D*k1 - B*k2) / (A*D - B*C)
            y = (-C*k1 + A*k2) / (A*D - B*C)

            return Vector([x, y])

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
            self.C.parallel_Q(self.D))

        self.assertEqual(False,
            self.C.parallel_Q(self.D))

        self.assertEqual(True,
            self.C.parallel_Q(self.D))

    def test_same_line_Q(self):
        self.assertTrue(self.A.same_line_Q(self.B))
        self.assertFalse(self.C.same_line_Q(self.D))
        self.assertFalse(self.E.same_line_Q(self.F))

    def test_random_point(self):
        self.assertTrue(self.intersection(self.random_point()))

    def test_intersection(self):
        self.assertEqual("infinity",
            self.A.intersection(self.B))

        self.assertEqual([1.173, 0.073],
            self.C.intersection(self.D))

        self.assertEqual(False,
            self.C.intersection(self.D))

A = Line([1,1], 1)
B = Line([2,2], 1)
C = Line(normal_vector=[4, 5], constant_term=6)

if __name__ == '__main__':
    unittest.main()
