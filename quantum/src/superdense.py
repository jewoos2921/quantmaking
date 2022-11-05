import math

from lib import state, ops


def alice_manipulate(psi: state.State,
                     bit0: int, bit1: int) -> state.State:
    """Alice encodes 2 classical bits in her 1 qubit."""

    # Note: this logic applies the Z-gate and X-gate to qubit0
    ret = ops.Identity(2)(psi)
    if bit0:
        ret = ops.PauliX()(ret)
    if bit1:
        ret = ops.PauliZ()(ret)
    return ret


def bob_measures(psi: state.State,
                 expect0: int, expect1: int) -> None:
    """Bob measures both bits (in computational basis)."""

    # Change Hadamard basis back to computational basis.
    psi = ops.Cnot(0, 1)(psi)
    psi = ops.Hadamard()(psi)
    p0, _ = ops.Measure(psi, 0, tostate=expect1)
    p1, _ = ops.Measure(psi, 1, tostate=expect0)

    if not math.isclose(p0, 1.0, abs_tol=1e-6) or not math.isclose(p1, 1.0, abs_tol=1e-6):
        raise AssertionError(f"Invalid result p0 {p0} p1 {p1}")

    print(f'Expected/matched: {expect0}{expect1}.')
