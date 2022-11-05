import math
from typing import Callable

import numpy as np

from lib import ops, state


def make_uf(f: Callable[[int], int]) -> ops.Operator:
    """Simple way to generate the 2-qubit, 4x4 Deutsch Oracle."""

    u = np.zeros(16).reshape(4, 4)
    for col in range(4):
        y = col & 1
        x = col & 2
        fx = f(x >> 1)
        xor = y ^ fx
        u[col][x + xor] = 1.0

    op = ops.Operator(u)
    if not op.is_unitary():
        raise AssertionError("Produced non-unitary operator.")
    return op


def make_f(flavor: int) -> Callable[[int], int]:
    """Return a 1-bit constant or balanced function f. 4 flavors."""

    # the 4 versions are:
    #   f(0) -> 0, f(1) -> 0 constant
    #   f(0) -> 0, f(1) -> 1 balanced
    #   f(0) -> 1, f(1) -> 0 balanced
    #   f(0) -> 1, f(1) -> 1 constant
    flavors = [[0, 0], [0, 1], [1, 0], [1, 1]]

    def f(bit: int) -> int:
        """Return f(bit) for one of the 4 possible function types."""
        return flavors[flavor][bit]

    return f


def run_experiment(flavor: int) -> None:
    """Run full experiment for a given flavor of f()."""

    f = make_f(flavor)
    u = make_uf(f)
    h = ops.Hadamard()

    psi = h(state.zeros(1)) * h(state.ones(1))
    psi = u(psi)
    psi = (h * ops.Identity())(psi)
    p0, _ = ops.Measure(psi, 0, tostate=0, collapse=False)

    print("f(0) = {:.0f} f(1) = {:.0f}".format(f(0), f(1)), end='')
    if math.isclose(p0, 0.0):
        print("  balanced")
        if flavor == 0 or flavor == 3:
            raise AssertionError("Invalid Result, expected balanced.")
    else:
        print("  constant")
        if flavor == 1 or flavor == 2:
            raise AssertionError("Invalid Result, expected constant.")
