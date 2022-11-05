import random
import timeit

from lib import ops, state


def apply_single_gate(gate, qubit, psi):
    """Apply a single qubit gate via explicit indexing."""

    qubit = psi.nbits - qubit - 1
    two_q = 2 ** qubit
    for g in range(0, 2 ** psi.nbits, 2 ** (qubit + 1)):
        for i in range(g, g + two_q):
            t1 = gate[0, 0] * psi[i] + gate[0, 1] * psi[i + two_q]
            t2 = gate[1, 0] * psi[i] + gate[1, 1] * psi[i + two_q]
            psi[i] = t1
            psi[i + two_q] = t2
    return psi


def single_gate_complexity() -> None:
    """Compare times for full matmul vs single-gate."""

    nbits = 12
    qubit = random.randint(0, nbits - 1)
    gate = ops.PauliX()

    def with_matmul():
        psi = state.zeros(nbits)
        op = ops.Identity(qubit) * gate * ops.Identity(nbits - qubit - 1)
        psi = op(psi)

    def apply_single():
        psi = state.zeros(nbits)
        psi = apply_single_gate(gate, qubit, psi)

    print("Time with full matmul :{:.3f} secs".format(timeit.timeit(with_matmul, number=1)))
    print("Time with single gate :{:.3f} secs".format(timeit.timeit(apply_single, number=1)))
