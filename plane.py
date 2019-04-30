from decimal import Decimal, getcontext

from vector import Vector

import unittest, random

getcontext().prec = 5


class Plane(object):

    NO_NONZERO_ELTS_FOUND_MSG = 'No nonzero elements found'

    def __init__(self, normal_vector=None, constant_term=None):
        """public attributes:
        dimension
        normal_vector
        constant_term

        TO DO
        watch more lecture
        """
        self.dimension = 3

        # if no normal vector is given make it the 0-plane
        if not normal_vector:
            all_zeros = ['0']*self.dimension
            normal_vector = Vector(all_zeros)
        self.normal_vector = Vector(normal_vector)

        if not constant_term:
            constant_term = Decimal('0')
        self.constant_term = Decimal(constant_term)

        # basepoint is any point on the plane; as implemented
        # that means a point on the plane that passes thru
        # one of the coordinate axes
        self.set_basepoint()

    def set_basepoint(self):
        """set the basepoint"""

        try:
            n = self.normal_vector
            c = self.constant_term
            basepoint_coords = ['0']*self.dimension
            # set all of {x,y,z} to 0
            initial_index = Plane.first_nonzero_index(n)
            # if any({x,y,z}) has non-zero coefficient, divide
            # the equation (e.g., Ax + 0y + 0z = c) by that
            # coefficient.
            initial_coefficient = Decimal(n[initial_index])
            # e.g., x then becomes c / A
            basepoint_coords[initial_index] = Decimal(c/initial_coefficient)
            # finally, e.g., the basepoint is [x,y,z]==[c/A, 0, 0]
            self.basepoint = Vector(basepoint_coords)

        except Exception as e:
            if str(e) == Plane.NO_NONZERO_ELTS_FOUND_MSG:
                self.basepoint = None
            else:
                raise e

    def __str__(self):
        """Represent the Plane in the form of a 3 term equation."""

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
            initial_index = Plane.first_nonzero_index(n)
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
    # end of __str__(self)

    __repr__ = __str__

    def __eq__(self, target):
        return self.equal_Q(target)

    def __getitem__(self, index):
        return self.normal_vector[index]

    def __setitem__(self, index, value):
        temp = list(self.normal_vector)
        temp[index] = Decimal(str(value))
        self.normal_vector = Vector(temp)

    def __contains__(self):
        pass

    def __mul__(self, val):
        return Plane(Vector(self.normal_vector)*Decimal(val),
                     self.constant_term*Decimal(val))

    __rmul__ = __mul__

    def __sub__(self, val):
        return self+-1*val

    __rsub__ = __sub__

    def __add__(self, val):
        if self.dimension != val.dimension:
            raise Exception("Can only add Planes() of same .dimension.")
        try:
            new_normal = self.normal_vector + val.normal_vector
            new_constant = self.constant_term + val.constant_term

            return Plane(new_normal, new_constant)
        except TypeError:
            raise TypeError("Can only add Planes() to other Planes().")

    __radd__ = __add__

    def parallel_Q(self, geometry):
        if type(geometry) is Plane:
            pln = geometry
            return self.normal_vector.parallel_Q(pln.normal_vector)

        if type(geometry) is Line:
            l = geometry
            return self.normal_vector.parallel_Q(l.normal_vector)

    def equal_Q(self, pln):
        if self.parallel_Q(pln):
            if self[0]*pln.constant_term == pln[0]*self.constant_term:
                return True
            else:
                return False
        else:
            return False

    def random_point(self, interval):
        pass

    def random_basis(self):
        pass

    @staticmethod
    def first_nonzero_index(iterable):
        for k, item in enumerate(iterable):
            if not MyDecimal(item).is_near_zero():
                return k
        raise Exception(Plane.NO_NONZERO_ELTS_FOUND_MSG)
#end Plane class

class MyDecimal(Decimal):
    def is_near_zero(self, eps=1e-10):
        return abs(self) < eps

class TestPlane(unittest.TestCase):
    def setUp(self):
        self.A = Plane([-0.412, 3.806, 0.728], -3.46)
        self.B = Plane([1.03, -9.515, -1.82], 8.65)
        self.C = Plane([2.611, 5.528, 0.283], 4.6)
        self.D = Plane([7.715, 8.306, 5.342], 3.76)
        self.E = Plane([-7.926, 8.625, -7.212], -7.952)
        self.F = Plane([-2.642, 2.875, -2.404], -2.443)

    def test_parallel_Q(self):
        self.assertTrue(self.A.parallel_Q(self.B))
        self.assertFalse(self.C.parallel_Q(self.D))
        self.assertTrue(self.E.parallel_Q(self.F))

    def test_equal_Q(self):
        self.assertTrue(self.A.equal_Q(self.B))
        self.assertFalse(self.C.equal_Q(self.D))
        self.assertFalse(self.E.equal_Q(self.F))

    def test_contains(self):
        pass
# end UnitTest class


A = Plane([-0.412, 3.806, 0.728], -3.46)
B = Plane([1.03, -9.515, -1.82], 8.65)
C = Plane([2.611, 5.528, 0.283], 4.6)
D = Plane([7.715, 8.306, 5.342], 3.76)
E = Plane([-7.926, 8.625, -7.212], -7.952)
F = Plane([-2.642, 2.875, -2.505], -2.443)

if __name__ == "__main__":
    unittest.main()
