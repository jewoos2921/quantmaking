from __future__ import annotations

import cmath
import math
from typing import Optional, Union, List

import numpy as np
from lib import tensor, state


class Operator(tensor.Tensor):
    """Operator are represented by square, unitary matrices."""

    def __repr__(self) -> str:
        s = 'Operator('
        s += super(Operator, self).__str__().replace("\n", "\n", " " * len(s))
        s += ')'
        return s

    def __str__(self) -> str:
        s = f'Operator for {self.nbits}-qubit state space.'
        s += ' Tensor:\n'
        s += super(Operator, self).__str__()
        return s

    def adjoint(self) -> Operator:
        return Operator(np.conj(self.transpose()))

    def dump(self, description: Optional[str] = None,
             zeros: bool = False) -> None:
        res = ''
        if description:
            res += f'{description} ({self.nbits}-qubits operator)\n'
        for row in range(self.shape[0]):
            for col in range(self.shape[1]):
                val = self[row, col]
                res += f'{val.real:+.1f}{val.imag:+.1f}j '
            res += '\n'
        if not zeros:
            res = res.replace("+0.0j", '     ')
            res = res.replace("+0.0", ' -   ')
            res = res.replace("-0.0", ' -   ')
            res = res.replace("+", ' ')
        print(res)

    def __call__(self, arg: Union[state.State, Operator], idx: int = 0) -> state.State:
        return self.apply(arg, idx)

    def apply(self, arg: Union[state.State, Operator], idx: int = 0) -> state.State:
        """Apply operator to a state or operator."""

        if isinstance(arg, Operator):
            arg_bits = arg.nbits
            if idx > 0:
                arg = Identity().kpow(idx) * arg

            if self.nbits > arg.nbits:
                arg = arg * Identity().kpow(self.nbits - idx - arg_bits)

            if self.nbits != arg.nbits:
                raise AssertionError("Operator(0) with mis-matched dimensions.")
            # Note: We reverse the order in this matmul. So:
            # X(Y) == Y @ X
            # This is to mirror that for a circuit like this:
            # --- X --- Y --- psi
            #
            # Incrementally updating states we would write:
            # psi = X(psi)
            # psi = Y(psi)
            #
            # But in a combined operator matrix, Y comes first:
            # (YX) (psi)
            #
            # The function call should mirror this semantic, since parameter are typically evaluated first
            # (and this mirrors the left to right in the circuit notation):
            # X(Y) = YX
            #
            return arg @ self

        if not isinstance(arg, state.State):
            raise AssertionError("Invalid parameter, expected State.")

        op = self
        if idx > 0:
            op = Identity().kpow(idx) * op

        if arg.nbits - idx - self.nbits > 0:
            op = op * Identity().kpow(arg.nbits - idx - self.nbits)

        # Note the reversed order compared to above.
        return state.State(np.matmul(op, arg))


def Identity(d: int = 1) -> Operator:
    return Operator(np.array([[1.0, 0.0], [0.0, 1.0]])).kpow(d)


def PauliX(d: int = 1) -> Operator:
    return Operator(np.array([[0.0j, 1.0], [1.0, 0.0j]])).kpow(d)


def PauliY(d: int = 1) -> Operator:
    return Operator(np.array([[0.0, -1.0j],
                              [1.0j, 0.0]])).kpow(d)


def PauliZ(d: int = 1) -> Operator:
    return Operator(np.array([[1.0, 0.0],
                              [0.0, -1.0]])).kpow(d)


# Cache Pauli matrices for performance reasons.
_PAULI_X = PauliX()
_PAULI_Y = PauliY()
_PAULI_Z = PauliZ()


def Rotation(v: np.ndarray, theta: float) -> np.ndarray:
    """Produce the single-qubit rotation operator."""
    v = np.asarray(v)
    if v.shape != (3,) or not math.isclose(v @ v, 1) or not np.all(np.isreal(v)):
        raise ValueError("Rotation vector must be 3D real unit vector.")

    return np.cos(theta / 2) * Identity() - 1j * np.sin(theta / 2) * (v[0] * _PAULI_X +
                                                                      v[1] * _PAULI_Y +
                                                                      v[2] * _PAULI_Z)


def RotationX(theta: float) -> Operator:
    return Rotation([1., 0., 0.], theta)


def RotationY(theta: float) -> Operator:
    return Rotation([0., 1., 0.], theta)


def RotationZ(theta: float) -> Operator:
    return Rotation([0., 0., 1.], theta)


# Phase gate, also called S or Z90. Rotate by 90 deg around z-axis
def Phase(d: int = 1) -> Operator:
    return Operator(np.array([[1.0, 0.0], [0.0, 1.0j]])).kpow(d)


# Phase gate is also called S-gate
def Sgate(d: int = 1) -> Operator:
    return Phase(d)


# Rk is one of the rotation gates used in Quantum Fourier Transform
def Rk(k: int, d: int = 1) -> Operator:
    return Operator(np.array([
        (1.0, 0.0),
        (0.0, cmath.exp(2.0 * cmath.pi * 1j / 2 ** k))
    ])).kpow(d)


# U1(lamda) gate
def U1(lam: float, d: int = 1) -> Operator:
    return Operator(np.array([
        (1.0, 0.0),
        (0.0, cmath.exp(1j * lam))
    ])).kpow(d)


# V-gate, which is sqrt(X)
# VV* = I
# V^2 = X.
def Vgate(d: int = 1) -> Operator:
    return Operator(0.5 *
                    np.array([
                        (1 + 1j, 1 - 1j),
                        (1 - 1j, 1 + 1j,)])).kpow(d)


# T-gate is equivalent to a 45' phase around the z-axis
def Tgate(d: int = 1) -> Operator:
    """T-gate is sqrt(S-gate)."""
    return Operator(
        np.array([[1.0, 0.0],
                  [0.0, cmath.exp(cmath.pi * 1j / 4)]])
    ).kpow(d)


def Yroot(d: int = 1) -> Operator:
    """Root of Y-gate"""
    return Operator(0.5 *
                    np.array([
                        (1 + 1j, -1 - 1j),
                        (1 + 1j, 1 + 1j)])).kpow(d)


def Projector(psi: state.State) -> Operator:
    """Construct projection Operator for basis state."""
    return Operator(psi.density())


def Hadamard(d: int = 1) -> Operator:
    return Operator(1 / np.sqrt(2) *
                    np.array([[1.0, 1.0],
                              [1.0, -1.0]])).kpow(d)


# Note on indices for controlled operators:
#
# The important aspects are direction and difference, not absolute values.
# In that regard, these are equivalent:
#     ControlledU(0, 3, U) == ControlledU(1, 4, U)
#     ControlledU(2, 0, U) == ControlledU(4, 2, U)
# We could have used -3 and 3, but felt this representation was more intuitive.
#
# Operator matrices are stored with all intermittent qubits (as Identities).
# When applying an operator, the starting qubit index can be specified
def ControlledU(idx0: int, idx1: int, u: Operator) -> Operator:
    """Control qubit at idx1 via controlling qubit at idx0."""
    if idx0 == idx1:
        raise ValueError("Control and controlled qubit must not be equal")

    p0 = Projector(state.zeros(1))
    p1 = Projector(state.ones(1))
    # Space between qubits
    ifill = Identity(abs(idx1 - idx0) - 1)
    # 'width' of U in terms of Identity matrices
    ufill = Identity().kpow(u.nbits)
    if idx1 > idx0:
        if idx1 - idx0 > 1:
            op = p0 * ifill * ufill + p1 * ifill * u
        else:
            op = p0 * ufill + p1 * u

    else:
        if idx0 - idx1 > 1:
            op = ufill * ifill * p0 + u * ifill * p1
        else:
            op = ufill * p0 + u * p1
    return op


def Cnot(idx0: int = 0, idx1: int = 1) -> Operator:
    """Controlled-Not between idx0 and idx1, controlled by |1>. """
    return ControlledU(idx0, idx1, PauliX())


def Cnot0(idx0: int = 0, idx1: int = 1) -> Operator:
    """Controlled-Not between idx0 and idx1, controlled by |0>."""
    if idx1 > idx0:
        x2 = PauliX() * Identity(idx1 - idx0)
    else:
        x2 = Identity(idx0 - idx1) * PauliX()
    return x2 @ ControlledU(idx0, idx1, PauliX()) @ x2


# Make Toffoli gate out of 2 controlled Cnot's.
#   idx1 and idx2 define the 'inner' cnot
#   idx0 defines the 'outer' cnot.
#
# For a Toffoli gate to control with qubit 5
# a Cnot from 4 and 1:
#   Toffolii(5, 4, 1)
#
def Toffoli(idx0: int, idx1: int, idx2: int) -> Operator:
    """Make a Toffoli gate."""
    cnot = Cnot(idx1, idx2)
    toffoli = ControlledU(idx0, idx1, cnot)
    return toffoli


def Swap(idx0: int = 0, idx1: int = 1) -> Operator:
    """Swap qubits at idx0 and idx1 via combination of Cnot gates."""
    return Cnot(idx1, idx0) @ Cnot(idx0, idx1) @ Cnot(idx1, idx0)


def TraceOutSingle(rho: Operator, index: int) -> Operator:
    """Trace out single qubit from density matrix."""

    nbits = int(math.log2(rho.shape[0]))
    if index > nbits:
        raise AssertionError(
            "Error in TraceOutSingle, invalid index (>nbits).")

    eye = Identity()
    zero = Operator(np.array([1.0, 0.0]))
    one = Operator(np.array([0.0, 1.0]))

    p0 = p1 = tensor.Tensor(1.0)
    for idx in range(nbits):
        if idx == index:
            p0 = p0 * zero
            p1 = p1 * one
        else:
            p0 = p0 * eye
            p1 = p1 * eye

    rho0 = p0 @ rho
    rho0 = rho0 @ p0.transpose()
    rho1 = p1 @ rho
    rho1 = rho1 @ p1.transpose()
    rho_reduced = rho0 + rho1
    return rho_reduced


def TraceOut(rho: Operator, index_set: List[int]) -> Operator:
    """Trace out multiple qubits from density matrix."""

    for idx, val in enumerate(index_set):
        nbits = int(math.log2(rho.shape[0]))
        if val > nbits:
            raise AssertionError(
                "Error in TraceOut, invalid index (>nbits).")
        rho = TraceOutSingle(rho, val)

        # Tracing out a bit of means that rho is now 1 bit smaller, the
        # indices right to the traced out qubit need to shift left by 1.
        # Example, to trace out bits 2, 4
        # Before:
        #   qubit 0 1 2 3 4 5
        #         a b c d e f
        # Trace out 2:
        #   qubit 0 1 <- 3 4 5
        #   qubit 0 1  2 3 4
        #         a b d  e f
        # Trace out 4 (is now 3)
        #   qubit 0 1 2 <- 4
        #   qubit 0 1 2 3
        #         a b d f
        for i in range(idx + 1, len(index_set)):
            index_set[i] = index_set[i] - 1
    return rho


def Measure(psi: state.State, idx: int,
            tostate: int = 0, collapse: bool = True) -> (float, state.State):
    """Measure a qubit via a projector on the density matrix."""

    # Compute probability of qubit(idx) to be in state 0 / 1.
    rho = psi.density()
    op = Projector(state.zero) if tostate == 0 else Projector(state.one)

    # Construct full matrix to apply to density matrix:
    if idx > 0:
        op = Identity().kpow(idx) * op

    if idx < psi.nbits - 1:
        op = op * Identity().kpow(psi.nbits - idx - 1)

    # Probability is the trace.
    prob0 = np.trace(np.matmul(op, rho))

    # Collapse state (don't forget to normalize if norm != 0)
    if collapse:
        mvmul = np.dot(op, psi)
        divisor = np.real(np.linalg.norm(mvmul))
        if divisor > 1e-10:
            normed = mvmul / divisor
        else:
            raise AssertionError("Measure() collapsed to 0.0 probability state.")
        return np.real(prob0), state.State(normed)

    # Return original state, enable chaining.
    return np.real(prob0), psi


def Pauli(d: int = 1) -> Operator:
    return Identity(d), PauliX(d), PauliY(d), PauliZ(d)
