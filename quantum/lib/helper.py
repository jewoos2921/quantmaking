import itertools
from typing import List, Tuple

import numpy as np

from lib import state, ops


def bits2val(bits):
    """For a given enumerable 'bits', compute the decimal integer."""
    # We assume bits are given in high to low order. For example, the bits [1, 1, 0] will
    # produce the value 6.
    return sum(v * (1 << (len(bits) - i - 1)) for i, v in enumerate(bits))


def val2bits(val: int, nbits: int) -> List[int]:
    """Convert decimal integer to list of {0, 1}."""
    # We return the bits in order high to low. For example,
    # the value 6 is being returned as [1, 1, 0].
    return [int(c) for c in format(val, "0{}b".format(nbits))]


def bitprod(nbits):
    """Produce the iterable cartesian of nbits {0, 1}."""
    for bits in itertools.product([0, 1], repeat=nbits):
        yield bits


def basis_changes():
    """Explore basis changes via Hadamard."""
    # Generate [0, 1]
    psi = state.ones(1)

    # Hadamard on |1> will result in 1/ sqrt(2) [1, -1]
    # aka |->
    psi = ops.Hadamard()(psi)

    # Simple PauliX will result in 1/ sqrt(2) [-1, 1]
    # which -1 (1/ sqrt(2) [-1, 1]).
    # Note that this does not move the vector on the Bloch sphere!
    psi = ops.PauliX()(psi)

    # Back to computational basis will result in -|1>.
    # Global phases can be ignored.
    pis = ops.Hadamard()(psi)
    if not np.allclose(psi[1], -1.0):
        raise AssertionError("Invalid basis change.")


def density_to_cartesian(rho: np.ndarray) -> Tuple[float, float, float]:
    """Compute Bloch Sphere coordinates from 2x2 density matrix."""

    a = rho[0, 0]
    c = rho[1, 0]
    x = 2.0 * c.real
    y = 2.0 * c.imag
    z = 2.0 * a - 1.0
    return np.real(x), np.real(y), np.real(z),

