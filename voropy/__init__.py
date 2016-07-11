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

    def __init__(self, points, origin=0.0, lattice=1.0, periodic=False, radii=None, blocks=None, frac_points=False):
        """Get the voronoi cells for a given set of points.

        Args:
            points (list or np.array): [N x 3] points to construct voronoi from
            origin (np.array): origin of container. A float will be cast to 3x1
                array.
            lattice (np.array): [3 x 3] lattice vectors for container based at origin. A
                float or list will be cast to diagonal of 3x3 array.
            periodic (list): [3 x 1] periodicity of each edge of the lattice. For orthogonal
                containers all faces can be either True or False. For non-orthogonal
                containers all faces must be Periodic.
            radii (np.array): [N x 1] radii of each point. Dimension of points
                must be same as dimension of radii.
            blocks (list): [3 x 1] blocking to use for unit cell. Usually leave as None
                since default constructor should work for most systems. Blocking works
                as expected for orthogonal systems. For non-orthogonal systems blocking
                is best described by looking at how it is initialized in __init__.
            frac_points (bool): Whether the points specified are in factional coordinates.
                Often times fractional coordinates are easier to work with.
        """
        self.origin = origin * np.array([1.0, 1.0, 1.0])

        self.lattice = np.array(lattice)
        if lattice.shape != (3, 3):
            self.lattice = lattice * np.eye(3)

        self._lattice = self._lattice_to_lowtri(self.lattice)

        self.periodic = np.array(periodic * np.array([1, 1, 1]), dtype=np.bool)

        if not all(periodic) and not self.orthogonal:
            raise ValueError("Box must be all periodic for non-orthogonal containers")

        if blocks is None:
            volume = 9 * self.lattice[0][0] * self.lattice[1][1] * self.lattice[1][1]
            Nthird = pow(len(points)/volume, 1.0/3.0)
            blocks = (
                round(Nthird * self._lattice[0][0]),
                round(Nthird * self._lattice[1][1]),
                round(Nthird * self._lattice[2][2])
            )
        self.blocks = np.array(blocks)

        if radii is None and self.orthogonal:
            self._container = _Container(
                self.origin[0], self._lattice[0][0],
                self.origin[1], self._lattice[1][1],
                self.origin[2], self._lattice[2][2],
                *blocks, *periodic, len(points)
            )
        elif self.orthogonal: #radii implied
            self._container = _ContainerPoly(
                self.origin[0], self._lattice[0][0],
                self.origin[1], self._lattice[1][1],
                self.origin[2], self._lattice[2][2],
                *blocks, *periodic, len(points)
            )
        elif radii is None:
            self._container = _ContainerPeriodic(
                self._lattice[0][0],
                self._lattice[1][0], self._lattice[1][1],
                self._lattice[2][0], self._lattice[2][1], self._lattice[2][2],
                *blocks, len(points)
            )
        else:
            self._container = _ContainerPeriodicPoly(
                self._lattice[0][0],
                self._lattice[1][0], self._lattice[1][1],
                self._lattice[2][0], self._lattice[2][1], self._lattice[2][2],
                *blocks, len(points)
            )

        for point in points:
            frac_point = point if frac_points else np.dot(point, np.linalg.inv(lattice))

            for i in range(3):
                if periodic[i]:
                    frac_point[i] = frac_point[i] % 1
                elif not (0 <= frac_point[i] <= 1):
                    raise ValueError('Point {} not within box in non-periodic axis'.format(point))
            point = np.dot(frac_point, self.lattice)

        self._put_points(points, radii)
