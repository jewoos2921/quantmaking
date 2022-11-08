import random
from typing import Callable
import sys
import numpy as np

from lib import state, ops, tensor

sys.path.append("xgates.cpp")

# import xgates
try:
    import libxgates as xgates

    apply1 = xgates.apply1
    applyc = xgates.applyc
except:
    print("""
  **************************************************************
  WARNING: Could not find 'libxgates.so'.
  Please build it and point PYTHONPATH to it.
  Execution is being re-directed to a Python implementation,
  performance may suffer greatly.
  **************************************************************
  """)


    def apply1(psi, gate, nbits, qubit, bitwidth=0):
        """Apply a single-qubit gate via explicit indexing."""

        qubit = nbits - qubit - 1
        two_q = 2 ** qubit
        for g in range(0, 2 ** nbits, 2 ** (qubit + 1)):
            for i in range(g, g + two_q):
                t1 = gate[0] * psi[i] + gate[1] * psi[i + two_q]
                t2 = gate[2] * psi[i] + gate[3] * psi[i + two_q]
                psi[i] = t1
                psi[i + two_q] = t2
        return psi


    def applyc(psi, gate, nbits, control, target, bitwidth=64):
        """Apply a controlled 2-qubit gate via explicit indexing."""

        qubit = nbits - target - 1
        two_q = 2 ** qubit
        control = nbits - control - 1
        for g in range(0, 2 ** nbits, 2 ** (qubit + 1)):
            for i in range(g, g + two_q):
                idx = g * 2 ** nbits + i
                if idx & (1 << control):
                    t1 = gate[0] * psi[i] + gate[1] * psi[i + two_q]
                    t2 = gate[2] * psi[i] + gate[3] * psi[i + two_q]
                    psi[i] = t1
                    psi[i + two_q] = t2
        return psi


def id(gate: ops.Operator) -> ops.Operator:
    return gate


def adjoint(gate: ops.Operator) -> ops.Operator:
    return gate.adjoint()


class qc:
    """Wrapper class to maintain state + operators."""

    def __init__(self, name=None):
        self.name = name
        self.psi = 1.0
        state.reset()

    def reg(self, size: int, it, *, name: str = None):
        ret = state.Reg(size, it, self.global_reg)
        self.global_reg = self.global_reg + size
        self.psi = self.psi * ret.psi()
        return ret

    def qubit(self,
              alpha: np.complexfloating = None,
              beta: np.complexfloating = None) -> None:
        self.psi = self.psi * state.qubit(alpha, beta)
        self.global_reg = self.global_reg + 1

    def zeros(self, n: int) -> None:
        self.psi = self.psi * state.zeros(n)
        self.global_reg = self.global_reg + n

    def ones(self, n: int) -> None:
        self.psi = self.psi * state.ones(n)
        self.global_reg = self.global_reg + n

    def bitstring(self, *bits) -> None:
        self.psi = self.psi * state.bitstring(*bits)
        self.global_reg = self.global_reg + len(bits)

    def arange(self, n: int) -> None:
        self.zeros(n)
        for i in range(0, 2 ** n):
            self.psi[i] = float(i)
        self.global_reg = self.global_reg + n

    def rand(self, n: int) -> None:
        self.psi = self.psi * state.rand(n)
        self.global_reg = self.global_reg + n

    @property
    def nbits(self) -> int:
        return self.psi.nbits

    def apply1(self, gate: ops.Operator, idx: int,
               name: str = None, *, val: float = None):
        if isinstance(idx, state.Reg):
            for reg in range(idx.nbits):
                xgates.apply1(self.psi, gate.reshape(4), self.psi.nbits,
                              idx[reg], tensor.tensor_width)
            return
        xgates.apply1(self.psi, gate.reshape(4), self.psi.nbits, idx,
                      tensor.tensor_width)

    def applyc(self, gate: ops.Operator, ctl: int, idx: int,
               name: str = None, *, val: float = None):
        if isinstance(idx, state.Reg):
            raise AssertionError("controlled register not supported")
        xgates.applyc(self.psi, gate.reshape(4), self.psi.nbits, ctl,
                      idx, tensor.tensor_width)

    def cv(self, idx0: int, idx1: int) -> None:
        self.applyc(ops.Vgate(), idx0, idx1, "cv")

    def cv_adj(self, idx0: int, idx1: int) -> None:
        self.applyc(ops.Vgate().adjoint(), idx0, idx1, "cv_adj")

    def cx(self, idx0: int, idx1: int) -> None:
        self.applyc(ops.PauliX(), idx0, idx1, "cx")

    def cy(self, idx0: int, idx1: int) -> None:
        self.applyc(ops.PauliY(), idx0, idx1, "cy")

    def cz(self, idx0: int, idx1: int) -> None:
        self.applyc(ops.PauliZ(), idx0, idx1, "cz")

    def cu1(self, idx0: int, idx1: int, value) -> None:
        self.applyc(ops.U1(value), idx0, idx1, "cu1", val=value)

    def crk(self, idx0: int, idx1: int, value) -> None:
        self.applyc(ops.Rk(value), idx0, idx1, "crk", val=value)

    def ccx(self, idx0: int, idx1: int, idx2: int) -> None:
        """Sleator-Weinfurter Construction."""

        self.cv(idx0, idx2)
        self.cx(idx0, idx1)
        self.cv_adj(idx1, idx2)
        self.cx(idx0, idx1)
        self.cv(idx1, idx2)

    def toffoli(self, idx0: int, idx1: int, idx2: int) -> None:
        self.ccx(idx0, idx1, idx2)

    def h(self, idx: int) -> None:
        self.apply1(ops.Hadamard(), idx, "h")

    def t(self, idx: int) -> None:
        self.apply1(ops.Tgate(), idx, "t")

    def u1(self, idx: int, val) -> None:
        self.apply1(ops.U1(val), idx, "u1", val=val)

    def v(self, idx: int) -> None:
        self.apply1(ops.Vgate(), idx, "v")

    def x(self, idx: int) -> None:
        self.apply1(ops.PauliX(), idx, "x")

    def y(self, idx: int) -> None:
        self.apply1(ops.PauliY(), idx, "y")

    def z(self, idx: int) -> None:
        self.apply1(ops.PauliZ(), idx, "z")

    def s(self, idx: int) -> None:
        self.apply1(ops.Sgate(), idx, "s")

    # def s(self, idx: int, trans: Callable = id) -> None:
    #     self.apply1(trans(ops.Sgate()), idx, 's')

    def sdag(self, idx: int) -> None:
        self.apply1(ops.Sgate().adjoint(), idx, "sdag")

    def yroot(self, idx: int) -> None:
        self.apply1(ops.Yroot(), idx, "yroot")

    def rx(self, idx: int, theta: float) -> None:
        self.apply1(ops.RotationX(theta), idx, 'rx')

    def ry(self, idx: int, theta: float) -> None:
        self.apply1(ops.RotationY(theta), idx, 'ry')

    def rz(self, idx: int, theta: float) -> None:
        self.apply1(ops.RotationZ(theta), idx, 'rz')

    # This doesn't go through the apply functions. Don't use.
    # def unitary(self, op, idx: int) -> None:
    #     self.psi = ops.Operator(op)(self.psi, idx, 'u')

    def measure_bit(self, idx: int, tostate: int = 0,
                    collapse: bool = True) -> (float, state.State):
        return ops.Measure(self.psi, idx, tostate, collapse)

    def sample_state(self, prob_state0: float):
        if prob_state0 < random.random():
            return 1  # corresponds to |1>
        return 0  # corresponds to |0>

    def swap(self, idx0: int, idx1: int) -> None:
        self.cx(idx1, idx0)
        self.cx(idx0, idx1)
        self.cx(idx1, idx0)

    def cswap(self, ct1: int, idx0: int, idx1: int) -> None:
        self.ccx(ct1, idx1, idx0)
        self.ccx(ct1, idx0, idx1)
        self.ccx(ct1, idx1, idx0)

    def mult_control(self, ct1, idx1, aux, gate, desc: str):
        """Multi-Controlled gate, using aux as ancilla"""

        # This is a simpler version that requires n-1 ancilla, instead of n-2.
        # The benefit is that the gate can be used as a single-controlled gate,
        # which means we don't need to take the root (no need to include scipy).
        # This construction also makes the Controlled-By-0 gates a little bit easier,
        # those controllers are being passed as single-element lists, eg.:
        #       ct1 = [1, 2, [3], [4], 5]
        #
        # This can be optimized (later) to turn into a space-optimized n-2 version.
        #
        # We also generalize to the case where ct1 is empty or only has 1 control qubit.
        # This is very flexible and practically any gate could be expressed this way.
        # This would make bulk control of whole gate sequences straightforward, but
        # changes the trivial IR we're working with here.
        # Something to keep in mind
        with self.scope(self.ir, f"multi({ct1}, {idx1}) # {desc}"):
            if len(ct1) == 0:
                self.apply1(gate, idx1, desc)
                return
            if len(ct1) == 1:
                self.applyc(gate, ct1[0], idx1, desc)
                return

            # compute the predicate
            self.ccx(ct1[0], ct1[1], aux[0])
            aux_idx = 0
            for i in range(2, len(ct1)):
                self.ccx(ct1[i], aux[aux_idx], aux[aux_idx + 1])
                aux_idx = aux_idx + 1

                # Use predicate to single-control qubit at idx1.
                self.applyc(gate, aux[aux_idx], idx1, desc)

                # Uncompute predicate.
                aux_idx = aux_idx - 1
                for i in range(len(ct1) - 1, 1, -1):
                    self.ccx(ct1[i], aux[aux_idx], aux[aux_idx + 1])
                    aux_idx = aux_idx - 1
                self.ccx(ct1[0], ct1[1], aux[0])
