from math import sqrt, acos, pi, log10
from decimal import Decimal, setcontext, Context

import unittest

setcontext(Context(prec=10))
Dec = Decimal

class Vector(object):
    def __init__(self, coordinates, delta=Dec(.001)):
        self.delta = delta
        try:
            if not coordinates:
                raise ValueError
            self.coordinates = tuple([Dec(str(e)) for e in coordinates])
            self.dimension = len(coordinates)

        except ValueError:
            raise ValueError('The coordinates must be nonempty')

        except TypeError:
            raise TypeError('The coordinates must be an iterable')

    @property
    def prec(self):
        return round(log10(1/self.delta))

    def __len__(self):
        return len(self.coordinates)

    def __eq__(self, v):
        try:
            zipped = zip(self.coordinates, v.coordinates)
            diffs = [round(e[0] - e[1], ndigits=self.prec) for e in zipped]
            for diff in diffs:
                if abs(diff) > self.delta:
                    return False

            return True
            # return all([abs(Dec(round(e[0] - e[1], ndigits=self.prec))) < self.delta
            #             for e in zip(self.coordinates, v.coordinates)])
        except IndexError:
            raise Exception("Cannot compare vectors of different dimension")

    def __add__(self, v):
        zipped = list(zip(self.coordinates, v.coordinates))
        sums = [round(e[0] + e[1], ndigits=self.prec) for e in zipped]
        return Vector([Dec(round(e[0]+e[1], ndigits=self.prec)) for
                       e in zip(self.coordinates, v.coordinates)])

    def __sub__(self, v):
        return Vector([Dec(round(e[0]-e[1], ndigits=self.prec)) for
                       e in zip(self.coordinates, v.coordinates)])

    def __mul__(self, q):
        if type(q) is Vector:
            return Dec(
                           round(
                                 sum([e[0]*e[1] for
                                     e in
                                     zip(self.coordinates, q.coordinates)]
                                     ),
                                 ndigits=self.prec
                                 )
                           )
        else:
            prods = [round(Dec(q)*e, ndigits=self.prec) for e in self.coordinates]
            return Vector(prods)

    __rmul__ = __mul__

    def __getitem__(self, index):
        return self.coordinates[index]

    def __setitem__(self, index, value):
        temp = list(self.coordinates)
        temp[index] = round(Dec(value), ndigits=self.prec)
        self.coordinates = tuple(temp)

    def __abs__(self):
        return Dec(round(sqrt(sum([elem * elem for elem in self.coordinates])), self.prec))

    def __repr__(self):
        return "V"+str([float(e) for e in self.coordinates])

    __str__ = __repr__

    def normalize(self):
        try:
            return Vector([elem / Dec(abs(self)) for elem in self.coordinates])

        except ZeroDivisionError:
            raise Exception("Cannot normalize 0 vector.")

    def angle(self, v=[0]):
        if type(v) is not Vector:
            v = Vector(v)

        #if the 0-vector is given as an argument
        #instead returnn the angle of self with respect
        #to the x-axis
        if abs(v) < self.delta:
            v = [0] * len(v)
            v[0] = 1
            v = Vector(v)
        cosine = round(self*v/(abs(self)*abs(v)), ndigits=self.prec)
        return Dec(acos(cosine))

    def projected_part(self, v):
        return (v*self)/(Dec(abs(v)*abs(v)))*v

    def orthogonal_part(self, v):
        return self - self.projected_part(v)

    def decompose(self, v):
        return (self.projected_part(v), self.orthogonal_part(v))

    def cross(self, v):
        """Gives cross-product of self and v."""

        if len(self) > 3 | len(v) > 3:
            raise Exception("only 3-vectors have defined cross product")

        a = self
        b = v

        cz = a[0]*b[1] - a[1]*b[0]
        cy = a[2]*b[0] - a[0]*b[2]
        cx = a[1]*b[2] - a[2]*b[1]

        return Vector([cx, cy, cz])

    def parallel_Q(self, target):
        """Returns true if the target is parallel with self."""

        # the 0 vector is parallel to all vectors
        if abs(self) < self.delta or abs(target) < self.delta:
            return True

        # if target is a ratio multiple of self then it must be parallel
        i = Vector.find_non_zero_indeces(self, target)
        if i >= 0:
            ratio = self[i] / target[i]
            return self == target*ratio
        else:
            return False

        # if the target's angle is 0 or 180 it must be parallel
        # if (self.angle(target) < self.delta or
        #     abs(self.angle(target) - Dec(3.14159265)) < self.delta):
        #     return True
        # else:
        #     return False

    def orthogonal_Q(self, v2):
        """Returns true if v2 is orthogonal to self."""

        # if the dot product is close to 0 they must be orthogonal
        if abs(self*v2) < self.delta:
            return True
        else:
            return False

    def is_zero(self):
        return abs(self) < self.delta

    @staticmethod
    def find_non_zero_indeces(a, b):
        for i, l in enumerate(zip(a,b)):
            if l[0] and l[1]:
                return i

        return None


class TestVectorMethods(unittest.TestCase):
    def setUp(self):
        self.A = Vector([1, 1])
        self.B = Vector([2, 2])
        self.C = Vector([1, -1])
        self.D = Vector([1, 2])
        self.E = Vector([-10, -10])
        self.F = Vector([1.234, 5.678])
        self.G = Vector([9.876, 5.432])
        self.ux = Vector([1, 0])
        self.uy = Vector([0, 1])
        self.test_vectors = [self.A,
                             self.B,
                             self.C,
                             self.D,
                             self.E,
                             self.F,
                             self.G,
                             self.ux,
                             self.uy]

    def test_angle(self):
        self.assertEqual(self.A.angle(self.B),
                         0)

        self.assertEqual(self.A.angle(self.C),
                         pi/2)

        self.assertEqual(round(self.A.angle(self.D), 3),
                         round(Dec('0.32176'), 3))

        for i in self.test_vectors:
            self.assertEqual(i.angle(i), 0)

    def test_projected_part(self):
        self.assertEqual(self.A.projected_part(self.C),
                         Vector([0, 0]))

        self.assertEqual(self.A.projected_part(self.ux),
                         Vector([1, 0]))

        self.assertEqual(self.A.projected_part(self.E),
                         self.A)

        self.assertEqual(self.F.projected_part(self.G),
                         Vector([3.345, 1.840]))

    def test_orthogonal_part(self):
        self.assertEqual(self.A.orthogonal_part(self.B),
                         Vector([0, 0]))

        self.assertEqual(self.A.orthogonal_part(self.uy),
                         Vector([1, 0]))

        self.assertEqual(self.F.orthogonal_part(self.G),
                         Vector([-2.111, 3.838]))

    def test_equality(self):
        self.assertFalse(self.A == self.B)
        self.assertTrue(self.A == self.A)
        self.assertTrue(self.A*-10 == self.E)

    def test_cross(self):
        pass

    def test_parallel_Q(self):

        #c = -2.5
        self.assertTrue(Vector([-0.412, 3.806, 0.728]).parallel_Q(
                        Vector([1.03, -9.515, -1.82])))

    def test_orthogonal_Q(self):
        pass


A = Vector([1, 1])
B = Vector([2, 2])
C = Vector([1, -1])
D = Vector([1, 2])
E = Vector([-10, -10])
F = Vector([1.234, 5.678])
G = Vector([9.876, 5.432])
ux = Vector([1, 0])
uy = Vector([0, 1])

if __name__ == "__main__":
    unittest.main()
