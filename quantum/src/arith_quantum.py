from lib import helper
from lib import circuit


def arith_quantum(n: int, init_a: int, init_b: int,
                  factor: float = 1.0, dumpit: bool = False) -> None:
    """Run a quantum add experiment."""
    qc = circuit.qc("qadd")
    a = qc.reg(n + 1, helper.val2bits(init_a, n)[::-1], name='a')
    b = qc.reg(n + 1, helper.val2bits(init_b, n)[::-1], name='b')

    for i in range(n + 1):
        qft(qc, a, n - i)

    for i in range(n + 1):
        evolve(qc, a, b, n - i, factor)

    for i in range(n + 1):
        inverse_qft(qc, a, i)
