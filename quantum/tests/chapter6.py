import unittest

from lib import state, ops, helper


class MyTestCase(unittest.TestCase):
    def test_phasekit(self):
        psi = state.bitstring(0, 0, 1)
        psi = ops.Hadamard(2)(psi)
        psi = ops.ControlledU(0, 2, ops.Sgate())(psi)
        psi = ops.ControlledU(1, 2, ops.Tgate())(psi, 1)
        psi.dump()

    def test_binary_fraction(self):
        val = helper.bit2sfrac((0,))
        print(val)

        val = helper.bit2sfrac((1,))
        print(val)

        val = helper.bit2sfrac((0, 1))
        print(val)

        val = helper.bit2sfrac((1, 0))
        print(val)

        val = helper.bit2sfrac((1, 1))
        print(val)

    def test_qft(self):
        psi = state.bitstring(0, 0)
        psi = ops.Hadamard()(psi)
        psi = ops.ControlledU(0, 1, ops.Sgate())(psi)
        psi = ops.Hadamard()(psi, 1)
        psi.dump()

    def test_sphere(self):
        psi = state.bitstring(1, 1)
        psi = ops.Qft(2)(psi)

        rho0 = ops.TraceOut(psi.density(), [1])
        rho1 = ops.TraceOut(psi.density(), [0])

        x0, y0, z0 = helper.density_to_cartesian(rho0)
        x1, y1, z1 = helper.density_to_cartesian(rho1)

        print("x0: {:.1f} y0: {:.1f} z0: {:.1f}".format(x0, y0, z0))
        print("x1: {:.1f} y1: {:.1f} z1: {:.1f}".format(x1, y1, z1))

if __name__ == '__main__':
    unittest.main()
