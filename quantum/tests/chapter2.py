import math
import random
import unittest

import numpy as np
import scipy.stats

from lib import state, helper, ops, bell


class MytestForChapter2(unittest.TestCase):
    def test_bit_length(self):
        psi = state.qubit(alpha=1.0)  # corresponds to |0>
        phi = state.qubit(beta=1.0)  # corresponds to |1>
        comb = psi * phi
        print(comb)

    def test_bit_equation(self):
        p1 = state.qubit(alpha=random.random())
        x1 = state.qubit(alpha=random.random())
        psi = p1 * x1

        # inner product of full state
        self.assertTrue(np.allclose(np.inner(psi.conj(), psi), 1.0))
        # inner product of the constituents multiplied
        self.assertTrue(np.allclose(np.inner(p1.conj(), p1) * (np.inner(x1.conj(), x1)), 1.0))

    def test_pis_bit_string(self):
        psi = state.bitstring(1, 0, 1, 0)
        psi.ampl(1, 0, 1, 1)
        psi.prob(1, 0, 1, 1)
        for bits in helper.bitprod(4):
            print(psi.prob(*bits))

        # umat =scipy.stats.unitary_group.rvs(2 **nbits)

    def test_ops(self):
        psi = state.bitstring(0, 0, 0)
        op = ops.Identity() * ops.PauliX() * ops.Identity()
        psi = op(psi)
        psi.dump()

        psi = state.bitstring(0, 0, 0)
        opx = ops.Identity() * ops.PauliX() * ops.Identity()
        opy = ops.Identity() * ops.Identity() * ops.PauliY()
        big_op = opx(opy)
        psi = big_op(psi)

        psi = state.bitstring(0, 0, 0)
        opx = ops.Identity() * ops.PauliX() * ops.Identity()
        psi = opx(psi)
        opy = ops.Identity() * ops.Identity() * ops.PauliY()
        psi = opy(psi)
        psi = state.bitstring(0, 0, 0)
        opxy = ops.Identity() * ops.PauliX() * ops.PauliY()
        psi = opxy(psi)

        psi = state.bitstring(0, 0, 0)
        psi = ops.PauliX()(psi, 1)

    def test_rotation(self):
        rz = ops.RotationZ(math.pi)
        rz.dump("RotationZ pi/2")
        rs = ops.Sgate()
        rs.dump('S-gate')

        psi = state.qubit(random.random())
        psi.dump("Random state")
        ops.Sgate()(psi).dump("After applying S-gate")
        ops.RotationZ(math.pi)(psi).dump("After applying RotationZ")

    def test_rk_ul(self):
        for i in range(10):
            u1 = ops.U1(2 * math.pi / (2 ** i))
            rk = ops.Rk(i)
            self.assertTrue(u1.is_close(rk))

    def test_rk(self):
        rk0 = ops.Rk(0)
        self.assertTrue(rk0.is_close(ops.Identity()))

        rk1 = ops.Rk(1)
        self.assertTrue(rk1.is_close(ops.PauliZ()))

        rk2 = ops.Rk(2)
        self.assertTrue(rk2.is_close(ops.Sgate()))

        rk3 = ops.Rk(3)
        self.assertTrue(rk3.is_close(ops.Tgate()))

    def test_t_gate(self):
        """Test that T^2 == S."""
        t = ops.Tgate()
        self.assertTrue(t(t).is_close(ops.Phase()))
        self.assertTrue(t(t).is_close(ops.Sgate()))

    def test_v_gate(self):
        """Test that V^2 == X."""
        v = ops.Vgate()
        self.assertTrue(v(v).is_close(ops.PauliX()))

    def test_yroot_gate(self):
        """Test that Yroot^2 == Y."""

        yr = ops.Yroot()
        self.assertTrue(yr(yr).is_close(ops.PauliY()))

        computed_yroot = scipy.linalg.sqrtm(ops.PauliY())
        self.assertTrue(ops.Yroot().is_close(computed_yroot))

    def test_braket(self):
        psi = state.zeros(1)
        psi = ops.PauliX()(psi)
        psi = ops.PauliY()(psi)
        psi = ops.PauliZ()(psi)

        psi = state.zeros(1)
        op = ops.PauliX(ops.PauliY(ops.PauliZ()))
        psi = op(psi)

        psi = state.zeros(1)
        psi = (ops.PauliZ() @ (ops.PauliZ() @ ops.PauliX()))(psi)

    def test_bloch(self):
        psi = state.zeros(1)
        x, y, z = helper.density_to_cartesian(psi.density())
        self.assertEqual(x, 0.0)
        self.assertEqual(y, 0.0)
        self.assertEqual(z, 1.0)

        psi = ops.PauliX()(psi)
        x, y, z = helper.density_to_cartesian(psi.density())
        self.assertEqual(x, 0.0)
        self.assertEqual(y, 0.0)
        self.assertEqual(z, -1.0)

        psi = ops.Hadamard()(psi)
        x, y, z = helper.density_to_cartesian(psi.density())
        self.assertTrue(math.isclose(x, -1.0, abs_tol=1e-6))
        self.assertTrue(math.isclose(y, 0.0, abs_tol=1e-6))
        self.assertTrue(math.isclose(z, 0.0, abs_tol=1e-6))

    def test_entangler(self):
        psi = state.zeros(2)
        op = ops.Hadamard() * ops.Identity()
        psi = op(psi)
        print(psi)
        psi = ops.Cnot(0, 1)(psi)
        print(psi)

    def test_trace(self):
        q0 = state.qubit(alpha=0.5)
        q1 = state.qubit(alpha=0.8660254)
        psi = q0 * q1
        reduced = ops.TraceOut(psi.density(), [1])
        self.assertTrue(math.isclose(np.real(np.trace(reduced)), 1.0))
        self.assertTrue(math.isclose(np.real(reduced[0, 0]),
                                     0.25, abs_tol=1e-6))
        self.assertTrue(math.isclose(np.real(reduced[1, 1]),
                                     0.75, abs_tol=1e-6))

        reduced = ops.TraceOut(psi.density(), [0])
        self.assertTrue(math.isclose(np.real(np.trace(reduced)), 1.0))
        self.assertTrue(math.isclose(np.real(reduced[0, 0]),
                                     0.75, abs_tol=1e-6))
        self.assertTrue(math.isclose(np.real(reduced[1, 1]),
                                     0.25, abs_tol=1e-6))
        psi = bell.bell_state(0, 0)
        reduced = ops.TraceOut(psi.density(), [1])
        self.assertTrue(math.isclose(np.real(np.trace(reduced)),
                                     1.0, abs_tol=1e-6))
        self.assertTrue(math.isclose(np.real(reduced[0, 0]),
                                     0.5, abs_tol=1e-6))
        self.assertTrue(math.isclose(np.real(reduced[1, 1]),
                                     0.5, abs_tol=1e-6))

    def test_measure(self):
        psi = state.bitstring(1, 0, 1, 0)
        psi.dump()
        p0, _ = ops.Measure(psi, 1)
        print(p0)
        # p0, _ = ops.Measure(psi, 1, tostate=1)
        # print(p0)
        psi = bell.bell_state(0, 0)
        psi.dump()

        psi = bell.bell_state(0, 0)
        p0, _ = ops.Measure(psi, 0, 0, collapse=False)
        print("Probability: ", p0)
        psi.dump()

        psi = bell.bell_state(0, 0)
        p0, psi = ops.Measure(psi, 0, 0, collapse=True)
        print("Probability: ", p0)
        psi.dump()


if __name__ == '__main__':
    unittest.main()
