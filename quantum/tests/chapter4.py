import unittest

from lib import state
from lib import circuit


class MyTestCase(unittest.TestCase):
    def test_something(self):
        data = state.Reg(4, (1, 0, 1, 1), 0)  # 0b1011
        ancilla = state.Reg(3, 7, 4)  # 0b111

    def test_qc(self):
        qc = circuit.qc("test")
        qc.arange(4)
        print(qc.psi)

    def test_qc2(self):
        # Let's try this for qubits 0 to 3
        for idx in range(4):
            qc = circuit.qc("test")

            # Populate vector with values 0 to 15.
            qc.arange(4)

            # Apply X-gate to qubit at index 'idx'.
            qc.x(idx)
            print('Applied X to qubit {}: {}'.format(idx, qc.psi))


if __name__ == '__main__':
    unittest.main()
