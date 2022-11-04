import math
import unittest

from lib import ops, state


class MyTestCase(unittest.TestCase):
    def test_gate_equivalences(self):
        psi = ops.Hadamard()(state.zero)
        psi.dump()

    def test_t_gate(self):
        """Test that T^2 == S."""
        s = ops.Tgate() @ ops.Tgate()
        self.assertTrue(s.is_close(ops.Phase()))

    def test_s_gate(self):
        """Test that S^2 == Z."""
        x = ops.Sgate() @ ops.Sgate()
        self.assertTrue(x.is_close(ops.PauliZ()))

    def test_v_gate(self):
        """Test that V^2 == X."""
        v = ops.Vgate() @ ops.Vgate()
        self.assertTrue(v.is_close(ops.PauliX()))

    def test_had_cnot_had(self):
        h2 = ops.Hadamard(2)
        cnot = ops.Cnot(0, 1)
        op = h2(cnot(h2))
        self.assertTrue(op.is_close(ops.Cnot(1, 0)))
        (ops.Hadamard(2) @ ops.Cnot(0, 1) @ ops.Hadamard(2)).dump()
        ops.Cnot(1, 0).dump()

    def test_controlled_z(self):
        z0 = ops.ControlledU(0, 1, ops.PauliZ())
        z1 = ops.ControlledU(1, 0, ops.PauliZ())
        self.assertTrue(z0.is_close(z1))

    def test_xyx(self):
        x = ops.PauliX()
        y = ops.PauliY()
        print(y)
        op = x(y(x))
        print(op)
        self.assertTrue(op.is_close(-1.0 * y))

    def test_equalities(self):
        # Generate the Pauli and Hadamard matrices.
        _, x, y, z = ops.Pauli()
        h = ops.Hadamard()

        # Check equalities.
        op = h(x(h))
        self.assertTrue(op.is_close(z))

        op = h(y(h))
        self.assertTrue(op.is_close(-1.0 * y))

        op = h(z(h))
        self.assertTrue(op.is_close(x))

        op = x(z)
        self.assertTrue(op.is_close(1.0j * y))

    def test_global_phase(self):
        h = ops.Hadamard()
        op = h(ops.Tgate()(h))

        # If equal up to a global phase, all values should be eqaul.
        phase = op / ops.RotationX(math.pi / 4)
        self.assertTrue(math.isclose(phase[0, 0].real,
                                     phase[0, 1].real,
                                     abs_tol=1e-6))
        self.assertTrue(math.isclose(phase[0, 0].imag,
                                     phase[0, 1].imag,
                                     abs_tol=1e-6))
        self.assertTrue(math.isclose(phase[0, 0].real,
                                     phase[1, 0].real,
                                     abs_tol=1e-6))
        self.assertTrue(math.isclose(phase[0, 0].imag,
                                     phase[1, 0].imag,
                                     abs_tol=1e-6))
        self.assertTrue(math.isclose(phase[0, 0].real,
                                     phase[1, 1].real,
                                     abs_tol=1e-6))
        self.assertTrue(math.isclose(phase[0, 0].imag,
                                     phase[1, 1].imag,
                                     abs_tol=1e-6))

    def test_v_vdag_v(self):
        # Make Toffoli out of V = sqrt(X)
        #
        v = ops.Vgate()  # Could be any unitary, in principle!
        ident = ops.Identity()
        cnot = ops.Cnot(0, 1)

        o0 = ident * ops.ControlledU(1, 2, v)
        c2 = cnot * ident
        o2 = (ident * ops.ControlledU(1, 2, v.adjoint()))
        o4 = ops.ControlledU(0, 2, v)
        final = o4 @ c2 @ o2 @ c2 @ o0

        v2 = v @ v
        cv1 = ops.ControlledU(1, 2, v2)
        cv0 = ops.ControlledU(0, 1, cv1)
        self.assertTrue(final.is_close(cv0))

    def test_control_equalities(self):
        """Exercise 4.13 Nielson, Chuana."""

        i, x, y, z = ops.Pauli()
        x1 = x * i
        x2 = i * x
        y1 = y * i
        y2 = i * y
        z1 = z * i
        z2 = i * z
        c = ops.Cnot(0, 1)
        theta = 25.0 * math.pi / 180.0
        rx2 = i * ops.RotationX(theta)
        rz1 = ops.RotationZ(theta) * i

        self.assertTrue(c(x1(c)).is_close(x1(x2)))
        self.assertTrue((c @ x1 @ c).is_close(x1 @ x2))
        self.assertTrue((c @ y1 @ c).is_close(y1 @ x2))
        self.assertTrue((c @ z1 @ c).is_close(z1))
        self.assertTrue((c @ x2 @ c).is_close(x2))
        self.assertTrue((c @ y2 @ c).is_close(z1 @ y2))
        self.assertTrue((c @ z2 @ c).is_close(z1 @ z2))
        self.assertTrue((rz1 @ c).is_close(c @ rz1))
        self.assertTrue((rx2 @ c).is_close(c @ rx2))


if __name__ == '__main__':
    unittest.main()
