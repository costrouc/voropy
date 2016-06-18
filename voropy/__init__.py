"""
This is a library to calculate Voronoi cells and access their information.

Basic Process
~~~~~~~~~~~~~

  - Create a :class:`Container` object, using information about your system.
      - a  :class:`Container` is a `list` of :class:`Cell` objects
  - Access the :class:`Cell` methods to get information about them

Example
~~~~~~~

    >>> from tess import Container
    >>> c = Container([[1,1,1], [2,2,2]], limits=(3,3,3), periodic=False)
    >>> [round(v.volume(), 3) for v in c]
    [13.5, 13.5]
"""

from ._voro import (
    Container as _Container,
    ContainerPoly as _ContainerPoly,
    ContainerPeriodic as _ContainerPeriodic,
    ContainerPeriodicPoly as _ContainerPeriodicPoly,
    Cell
)


import numpy as np
from math import acos, sqrt, cos, sin


class Container:
    @staticmethod
    def _lattice_to_lowtri(lattice):
        """ Converts lattice to lower triangular form with lattice[i][i] > 0

        This is a requirement of non-orthogonal containers in voro++
        """
        def angle(v1, v2):
            return acos(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2)))

        a, b, c = map(np.linalg.norm, lattice)
        alpha, beta, gamma = (
            angle(lattice[1], lattice[2]),
            angle(lattice[0], lattice[2]),
            angle(lattice[0], lattice[1])
        )
        temp = (cos(alpha) - cos(gamma)*cos(beta)) / sin(gamma)
        return np.array([
            [a, 0.0, 0.0],
            [b * cos(gamma), b * sin(gamma), 0.0],
            [c * cos(beta), c * temp, c * sqrt(1.0 - (cos(beta)**2 - temp**2))]
        ])

    @property
    def orthogonal(self):
        """ Determines whether lattice is orthogonal

        """
        return all([(abs(e) < 1e-10 or i==j) for i, row in enumerate(self.lattice) for j, e in enumerate(row)])

    def _put_points(self, points, radii):
        if radii and len(points) != len(radii):
            raise ValueError('Number of radii must match number of points')

        if radii:
            for i, (x,y,z), r in zip(enumerate(points), radii):
                self._container.put(i, x, y, z, r)
        else:
            for i, (x,y,z) in enumerate(points):
                self._container.put(i, x, y, z)

    def __init__(self, points, origin=0.0, lattice=1.0, periodic=False, radii=None, blocks=None):
        """Get the voronoi cells for a given set of points."""
        periodic = np.array(periodic * np.array([1, 1, 1]), dtype=np.bool)
        lattice = lattice * np.eye(3)

        self.origin = origin * np.array([1.0, 1.0, 1.0])
        self.lattice = self._lattice_to_lowtri(lattice)

        if not all(periodic) and not self.orthogonal:
            raise ValueError("Box must be all periodic for non-orthogonal containers")

        if blocks is None:
            volume = 9 * self.lattice[0][0] * self.lattice[1][1] * self.lattice[1][1]
            Nthird = pow(len(points)/volume, 1.0/3.0)
            blocks = (
                round(Nthird * self.lattice[0][0]),
                round(Nthird * self.lattice[1][1]),
                round(Nthird * self.lattice[2][2])
            )
        self.blocks = np.array(blocks)

        if radii is None and self.orthogonal:
            self._container = _Container(
                self.origin[0], self.lattice[0][0],
                self.origin[1], self.lattice[1][1],
                self.origin[2], self.lattice[2][2],
                *blocks, *periodic, len(points)
            )
        elif radii is None and not self.orthogonal:
            self._container = _ContainerPeriodic(
                self.lattice[0][0],
                self.lattice[1][0], self.lattice[1][1],
                self.lattice[2][0], self.lattice[2][1], self.lattice[2][2],
                *blocks, len(points)
            )
        elif self.orthogonal:   #radii not None implied
            self._container = _ContainerPeriodicPoly(
                self.lattice[0][0],
                self.lattice[1][0], self.lattice[1][1],
                self.lattice[2][0], self.lattice[2][1], self.lattice[2][2],
                *blocks, len(points)
            )
        else:
            self._container = _Container(
                self.origin[0], self.lattice[0][0],
                self.origin[1], self.lattice[1][1],
                self.origin[2], self.lattice[2][2],
                *blocks, *periodic, len(points)
            )

        for point in points:
            frac_point = np.dot(point, np.linalg.inv(lattice))
            for i in range(3):
                if periodic[i]:
                    frac_point[i] = frac_point[i] % 1
                elif not (0 <= frac_point[i] <= 1):
                    raise ValueError('Point {} not within box in non-periodic axis'.format(point))
            point = np.dot(frac_point, self.lattice)

        self._put_points(points, radii)
