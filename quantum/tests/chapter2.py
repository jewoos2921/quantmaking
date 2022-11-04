import random
import unittest

import numpy as np

from lib import state


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


if __name__ == '__main__':
    unittest.main()

# %%
