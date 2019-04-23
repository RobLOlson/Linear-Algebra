from math import sqrt, acos
from decimal import Decimal


class Vector(object):
    def __init__(self, coordinates):
        try:
            if not coordinates:
                raise ValueError
            self.coordinates = tuple([Decimal(e) for e in coordinates])
            self.dimension = len(coordinates)

        except ValueError:
            raise ValueError('The coordinates must be nonempty')

        except TypeError:
            raise TypeError('The coordinates must be an iterable')

    def __str__(self):
        return 'Vector: {}'.format(self.coordinates)

    def __len__(self):
        return len(self.coordinates)

    def __eq__(self, v):
        return self.coordinates == v.coordinates

    def __add__(self, v):
        return Vector([e[0]+e[1] for
                       e in zip(self.coordinates, v.coordinates)])

    def __sub__(self, v):
        return Vector([e[0]-e[1] for
                       e in zip(self.coordinates, v.coordinates)])

    def __mul__(self, c):
        if type(c) is Vector:
            return Decimal(sum([e[0]*e[1] for
                                e in
                                zip(self.coordinates, c.coordinates)]))
        else:
            return Vector([Decimal(c)*elem for elem in self.coordinates])

    __rmul__ = __mul__

    def __getitem__(self, index):
        return self.coordinates[index]

    def __setitem__(self, index, value):
        temp = list(self.coordinates)
        temp[index] = value
        self.coordinates = tuple(temp)

    def __abs__(self):
        return sqrt(sum([elem * elem for elem in self.coordinates]))

    def __repr__(self):
        return str([float(e) for e in self.coordinates])

    def normalize(self):
        try:
            return Vector([elem / Decimal(abs(self)) for elem in self.coordinates])

        except ZeroDivisionError:
            raise Exception("Cannot normalize 0 vector.")

    def angle(self, v=[0]):
        if type(v) is not Vector:
            v = Vector(v)

        if abs(v) < 0.0001:
            v = [0] * len(v)
            v[0] = 1
            v = Vector(v)

        cosine = self*v/(Decimal(abs(self))*Decimal(abs(v)))
        return acos(cosine)

    def projected_onto(self, v):
        return (v*self)/(Decimal(abs(v)*abs(v)))*v

    def orthogonal_to(self, v):
        return self - self.projected_onto(v)

    def decompose(self, v):
        return (self.projected_onto(v), self.orthogonal_to(v))

    def cross(self, v):
        if len(self) > 3 | len(v) > 3:
            raise Exception("only 3-vectors have defined cross product")

        a = self
        b = v

        #u1 = Vector([1,0,0])
        #u2 = Vector([0,1,0])
        #u3 = Vector([0,0,1])

        #ax = len(self.projected_onto(u1))
        #ay = len(self.projected_onto(u2))
        #az = len(self.projected_onto(u3))

        #bx = len(v.projected_onto(u1))
        #by = len(v.projected_onto(u2))
        #bz = len(v.projected_onto(u3))

        cz = a[0]*b[1] - a[1]*b[0]
        cy = a[2]*b[0] - a[0]*b[2]
        cx = a[1]*b[2] - a[2]*b[1]

        return Vector([cx, cy, cz])

    def parallel_Q(self, v2):
        """Returns true if v2 is parallel with self."""

        if abs(self) < 0.001 or abs(v2) < 0.0001:
            return True
        if self.angle(v2) < 0.001 or abs(self.angle(v2) - 3.14159265) < 0.001:
            return True
        else:
            return False

    def orthogonal_Q(self, v2):
        """Returns true if v2 is orthogonal to self."""

        if abs(self*v2) < 0.001:
            return True
        else:
            return False

A = Vector([1, 2])
B = Vector([3, 4])
A.parallel_Q
