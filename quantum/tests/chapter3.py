import math
import unittest

from lib import ops, state, bell
from src import arith_classic, superdense, deutsch, deutsch_jozsa, simon, simon_general


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

    def test_classic(self):
        arith_classic.add_classic()

    def test_swap(self):
        print("Swap test. 0.5 means different, 1.0 means similar")
        arith_classic.run_experiment(1., 0.0, 0.5)
        arith_classic.run_experiment(0., 1.0, 0.5)
        arith_classic.run_experiment(1., 1.0, 1.)
        arith_classic.run_experiment(0., 0.0, 1.0)
        arith_classic.run_experiment(.1, 0.9, 0.65)

    def test_teleportation(self):
        # Step 1: Alice and Bob share an entangled pair, and separate.
        psi = bell.bell_state(0, 0)
        # Step 2: Alice wants to teleport a qubit |x> to Bob,
        #           which is in the state:
        #           |x> = a|0> + b|1> (with a^2 + b^2 == 1)
        a = 0.6
        b = math.sqrt(1.0 - a * a)
        x = state.qubit(a, b)
        print("Quantum Teleportation")
        print("Start with EPR Pair a={:.2f}, b={:.2f}".format(a, b))
        # Produce combined state 'alice'
        alice = x * psi
        # Alice lets qubit 0 (|x>) interact with qubit 1, which is her part of
        # the entangled state with Bob.
        alice = ops.Cnot(0, 1)(alice)
        # Now she applies a Hadamard to qubit 0. Bob still owns qubit 2.
        alice = ops.Hadamard()(alice, idx=0)

        # Alice measures and communicates the result |00>, |01>, ... to Bob.
        arith_classic.alice_measures(alice, a, b, 0, 0)
        arith_classic.alice_measures(alice, a, b, 0, 1)
        arith_classic.alice_measures(alice, a, b, 1, 0)
        arith_classic.alice_measures(alice, a, b, 1, 1)

    def test_superdense(self):
        # Step 1: Alice and Bob share an entangled pair, and separate.
        psi = bell.bell_state(0, 0)
        # Alice manipulates her qubit and sends her 1 qubit back to Bob,
        # who measures. In the Hadamard basis he would get b00, b01, etc.
        # but we're measuring in the computational basis by reverse
        # applying Hadamard and Cnot.

        for bit0 in range(2):
            for bit1 in range(2):
                psi_alice = superdense.alice_manipulate(psi, bit0, bit1)
                superdense.bob_measures(psi_alice, bit0, bit1)

    def test_deutsch(self):
        for i in range(4):
            f = deutsch.make_f(i)
            u = deutsch.make_uf(f)
            print(f"Flavor {i:02b} : {u}")

        print("\n")
        deutsch.run_experiment(0)
        deutsch.run_experiment(1)
        deutsch.run_experiment(2)
        deutsch.run_experiment(3)

    def test_deutsch_jozsa(self):
        for qubits in range(2, 8):
            result = deutsch_jozsa.run_experiment(qubits, deutsch_jozsa.exp_constant)
            print("Found: {} ({} qubits) (expected: {})".format(result, qubits, deutsch_jozsa.exp_constant))

            if result != deutsch_jozsa.exp_constant:
                raise AssertionError("Error, expected {}".format(deutsch_jozsa.exp_constant))

            result = deutsch_jozsa.run_experiment(qubits, deutsch_jozsa.exp_balanced)
            print("Found: {} ({} qubits) (expected: {})".format(result, qubits, deutsch_jozsa.exp_balanced))

            if result != deutsch_jozsa.exp_balanced:
                raise AssertionError("Error, expected {}".format(deutsch_jozsa.exp_balanced))

    def test_simon(self):
        simon.run_experiment()
        print("\n")
        # Note: Running with 5 (10 total) qubits takes about 5 minutes.
        simon_general.run_experiment(5)


if __name__ == '__main__':
    unittest.main()
