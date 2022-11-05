import cmath
import math
import random
from typing import Optional, List

import numpy as np

from lib import tensor, ops
from lib import helper


class State(tensor.Tensor):
    """Produce a given state for a single qubit."""

    def __repr__(self) -> str:
        s = 'State('
        s += super(State, self).__str__().replace("\n", "\n", " " * len(s))
        s += ')'
        return s

    def __str__(self) -> str:
        s = f'{self.nbits}-qubit state.'
        s += ' Tensor:\n'
        s += super(State, self).__str__()
        return s

    def ampl(self, *bits) -> np.complexfloating:
        """Return amplitude for state indexed by 'bits'."""
        idx = helper.bits2val(bits)  # in helper.py
        return self[idx]

    def prob(self, *bits) -> float:
        """Return probability for state indexed by 'bits'."""
        amplitude = self.ampl(*bits)
        return np.real(amplitude.conj() * amplitude)

    def maxprob(self) -> (List[float], float):
        """Find state with the highest probability."""
        maxbits = []
        maxprobs = 0.0
        for bits in helper.bitprod(self.nbits):
            cur_prob = self.prob(*bits)
            if cur_prob > maxprobs:
                maxbits, maxprobs = bits, cur_prob
        return maxbits, maxprobs

    def normalize(self) -> None:
        """Re-normalize the state. Sum of norms == 1.0."""
        dprod = np.conj(self) @ self
        if dprod.is_close(0.0):
            raise AssertionError("Normalizing to zero-probability state.")
        self /= np.sqrt(np.real(dprod))

    def phase(self, *bits) -> float:
        """Compute phase of a state from the complex amplitude."""
        amplitude = self.ampl(*bits)
        return math.degrees(cmath.phase(amplitude))

    def dump(self, desc: Optional[str] = None,
             prob_only: bool = True) -> None:
        dump_state(self, desc, prob_only)

    def density(self) -> tensor.Tensor:
        return tensor.Tensor(np.outer(self, self.conj()))

    def apply1(self, gate, index: int) -> None:
        """Apply single-qubit gate to this state."""

        # To maintain qubit ordering in this infrastructure
        # index needs to be reversed.
        #
        index = self.nbits - index - 1
        pow_2_index = 1 << index
        g00 = gate[0, 0]
        g01 = gate[0, 1]
        g10 = gate[1, 0]
        g11 = gate[1, 1]
        for g in range(0, 1 << self.nbits, 1 << (index + 1)):
            for i in range(g, g + pow_2_index):
                t1 = g00 * self[i] + g01 * self[i + pow_2_index]
                t2 = g10 * self[i] + g11 * self[i + pow_2_index]
                self[i] = t1
                self[i + pow_2_index] = t2

    def applyc(self, gate, ctl: int, target: int) -> None:
        """Apply a controlled 2-qubit gate via explicit indexing."""

        # To maintain qubit ordering in this infrastructure
        # index needs to be reversed.
        #
        qbit = self.nbits - target - 1
        pow_2_index = 2 ** qbit
        ctl = self.nbits - ctl - 1
        g00 = gate[0, 0]
        g01 = gate[0, 1]
        g10 = gate[1, 0]
        g11 = gate[1, 1]
        for g in range(0, 1 << self.nbits, 1 << (qbit + 1)):
            idx_base = g * (1 << self.nbits)
            for i in range(g, g + pow_2_index):
                idx = idx_base + i
                if idx & (1 << ctl):
                    t1 = g00 * self[i] + g01 * self[i + pow_2_index]
                    t2 = g10 * self[i] + g11 * self[i + pow_2_index]
                    self[i] = t1
                    self[i + pow_2_index] = t2


def state_to_strings(bits) -> str:
    """Convert state to string like |010>."""
    s = ''.join(str(i) for i in bits)
    return '|{:s}> (|{:d}>)'.format(s, int(s, 2))


def dump_state(psi, desc: Optional[str] = None,
               prob_only: bool = True) -> None:
    """Dump probabilities for a state. as well as local qubit state."""
    if desc:
        print("/", end='')
        for i in range(psi.nbits - 1, -1, -1):
            print(i % 10, end='')
        print(f"> \'{desc}\'")

    state_list: List[str] = []
    for bits in helper.bitprod(psi.nbits):
        if prob_only and (psi.prob(*bits) < 10e-6):
            continue
        state_list.append('{:s}: ampl: {:+.2f} prob: {:.2f} Phase: {:5.1f}'
                          .format(state_to_strings(bits),
                                  psi.ampl(*bits),
                                  psi.prob(*bits),
                                  psi.phase(*bits)))
    state_list.sort()
    print(*state_list, sep='\n')


def qubit(alpha: Optional[np.complexfloating] = None,
          beta: Optional[np.complexfloating] = None) -> State:
    """Produce a given state for a single qubit."""
    if alpha is None and beta is None:
        raise ValueError("alpha, or beta, or both are required")

    if beta is None:
        beta = math.sqrt(1.0 - np.conj(alpha) * alpha)

    if alpha is None:
        alpha = math.sqrt(1.0 - np.conj(beta) * beta)

    if not math.isclose(np.conj(alpha) * alpha +
                        np.conj(beta) * beta, 1.0):
        raise ValueError("Qubit probabilities do not add to 1.")

    qb = np.zeros(2, dtype=tensor.tensor_type())
    qb[0] = alpha
    qb[1] = beta
    return State(qb)


"""
The functions zeros() and ones() produce the all-zero or all-one computational basis vector 
for 'd' qubits, i.e., 
    |000...0> or 
    |111...1> 
    
The result of this tensor product is always [1, 0, 0, ..., 0] ^ T or [0, 0, 0, ..., 1] ^ T 
"""


def zeros_or_ones(d: int = 1, idx: int = 0) -> State:
    """Produce the all-0/1 basis vector for 'd' qubits."""
    if d < 1:
        raise ValueError("Rank must be at least 1.")
    shape = 2 ** d
    t = np.zeros(shape, dtype=tensor.tensor_type())
    t[idx] = 1
    return State(t)


def zeros(d: int = 1) -> State:
    """Produce state with 'd' |0>, eg., |0000>."""
    return zeros_or_ones(d, 0)


def ones(d: int = 1) -> State:
    """Produce state with 'd' |1>, eg., |1111>."""
    return zeros_or_ones(d, 2 ** d - 1)


def bitstring(*bits) -> State:
    """Produce a state from a given bit sequence, eg., |0101>."""

    d = len(bits)
    if d == 0:
        raise ValueError("Rank must be at least 1.")

    t = np.zeros(1 << d, dtype=tensor.tensor_type())
    t[helper.bits2val(bits)] = 1
    return State(t)


def rand(n: int) -> State:
    """Produce random combination of |0> and |1>."""
    bits = [0] * n
    for i in range(n):
        bits[i] = random.randint(0, 1)
    return bitstring(*bits)


# These two are used so commonly, make them constants.
zero = zeros(1)
one = ones(1)


class Reg:
    def __init__(self, size: int, it=0, global_reg: int = None):
        self.size = size
        self.global_idx = list(range(global_reg, global_reg + size))

        self.val = [0] * size

        if it:
            if isinstance(it, int):
                it = format(it, "0{}b".format(size))
            if isinstance(it, (str, tuple, list)):
                for idx, val in enumerate(it):
                    if val == '1' or val == 1:
                        self.val[idx] = 1

    def __str__(self) -> str:
        s = '/'
        for _, val in enumerate(self.val):
            s += f'{val}'
        return s + ">"

    def __getitem__(self, idx: int) -> int:
        return self.global_idx[idx]

    def __setitem__(self, idx: int, value: int) -> None:
        self.val[idx] = value

    @property
    def nbits(self) -> int:
        return self.size

    def psi(self) -> State:
        return bitstring(*self.val)
