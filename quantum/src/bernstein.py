from typing import Tuple

import numpy as np

from lib import ops, state, helper


def make_c(nbits: int) -> Tuple[bool]:
    """Make a random constant c from {0, 1}, the c we try to find."""

    constant_c = [0] * nbits
    for idx in range(nbits - 1):
        constant_c[idx] = int(np.random.random() < 0.5)
    return tuple(constant_c)


def make_u(nbits: int, constant_c: Tuple[bool]) -> ops.Operator:
    """Make general Bernstein oracle."""

    op = ops.Identity(nbits)
    for idx in range(nbits - 1):
        if constant_c[idx]:
            op = ops.Identity(idx) * ops.Cnot(idx, nbits - 1) @ op

    if not op.is_unitary():
        raise AssertionError("Constructed non-unitary operator.")
    return op


def run_experiments(nbits: int) -> None:
    """Run full experiment for a given number of bits."""

    c = make_c(nbits - 1)
    u = make_u(nbits, c)

    psi = state.zeros(nbits - 1) * state.ones(1)
    psi = ops.Hadamard(nbits)(psi)
    psi = u(psi)
    psi = ops.Hadamard(nbits)(psi)
    check_result(nbits, c, psi)


def check_result(nbits: int, c: Tuple[bool], psi: state.State) -> None:
    """Check expected vs. achieved results."""

    print(f"Expected: {c}")

    # The state with the 'flipped' bits will have probability 1.0
    # It will be found on the very first try.
    for bits in helper.bitprod(nbits):
        if psi.prob(*bits) > 0.1:
            print("Found   : {}, with prob: {:.1f}".format(bits[:-1], psi.prob(*bits)))
            if bits[:-1] != c:
                raise AssertionError("Invalid result")


# Alternative way to achieve the same result, using the Deutsch Oracle UF.
def make_oracle_f(c: Tuple[bool]) -> ops.Operator:
    """Return a function computing the dot product mod 2 of bits, c."""

    const_c = c

    def f(bit_string: Tuple[int]) -> int:
        val = 0
        for idx in range(len(bit_string)):
            val += const_c[idx] * bit_string[idx]
        return val % 2

    return f


def run_oracle_experiment(nbits: int) -> None:
    """Run full experiment for a given number of bits."""

    c = make_c(nbits - 1)
    f = make_oracle_f(c)
    u = ops.OracleUF(nbits, f)

    psi = state.zeros(nbits - 1) * state.ones(1)
    psi = ops.Hadamard(nbits)(psi)
    psi = u(psi)
    psi = ops.Hadamard(nbits)(psi)

    check_result(nbits, c, psi)
   