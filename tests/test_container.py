import pytest
from voropy import Container
import numpy as np


def test_container_init_defaults():
    points = [[0.5, 0.5, 0.5]]
    c = Container(points)
    assert np.allclose(c.origin, [0.0, 0.0, 0.0])
    assert np.allclose(c.lattice, np.eye(3))
    assert all(c.periodic == [False, False, False])


def test_container_init_lattice_float():
    points = [[0.5, 0.5, 0.5]]
    lattice = 5.0
    c = Container(points, lattice=lattice)
    assert np.allclose(c.origin, [0.0, 0.0, 0.0])
    assert np.allclose(c.lattice, lattice * np.eye(3))
    assert all(c.periodic == [False, False, False])
