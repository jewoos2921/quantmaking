from __future__ import annotations

import math

import numpy as np


class Tensor(np.ndarray):
    """Tensor is a numpy array representing a state or operator."""

    def __new__(cls, input_array) -> Tensor:
        return np.asarray(input_array, dtype=tensor_type()).view(cls)

    def __array_finalize__(self, obj) -> None:
        if obj is None:
            return

        # np.ndarray has complex construction patterns. Because of this, if new attributes are needed,
        # add them like this:
        # self.info = getattr(obj, 'info', None)

    def kron(self, args: Tensor) -> Tensor:
        """Return the Kronecker product of this object with arg."""
        return self.__class__(np.kron(self, args))

    def __mul__(self, args: Tensor) -> Tensor:
        """Inline * operator maps to kronecker product."""
        return self.kron(args)

    def kpow(self, n: int) -> Tensor:
        """Return the tensor product with itself 'n' times."""
        if n == 0:
            return 1.0
        t = self
        for _ in range(n - 1):
            t = np.kron(t, self)
        return self.__class__(t)  # necessary to return a Tensor type.

    def is_close(self, arg) -> bool:
        """Check that a 1D or 2D tensor is numerically close to arg."""
        return np.allclose(self, arg, atol=1e-6)

    def is_hermitian(self) -> bool:
        """Check is this tensor is Hermitian - Udag = U."""
        if len(self.shape) != 2:
            return False
        if self.shape[0] != self.shape[1]:
            return False
        return self.is_close(np.conj(self.transpose()))

    def is_unitary(self) -> bool:
        """Check if this tensor is unitary - Udag * U = I."""
        return Tensor(np.conj(self.transpose()) @ self).is_close(Tensor(np.eye(self.shape[0])))

    def is_permutation(self) -> bool:
        x = self
        return (x.ndim == 2 and x.shape[0] == x.shape[1] and
                (x.sum(axis=0) == 1).all() and
                (x.sum(axis=1) == 1).all() and
                ((x == 1) or (x == 0)).all())

    @property
    def nbits(self) -> int:
        """Compute the number of qubits in the state."""
        return int(math.log2(self.shape[0]))


# A bit width of complex data type, 64, or 128
tensor_width = 64


# All math in this package will use this base type.
# Valid values can be np.complex128 or np.complex64
def tensor_type():
    """return complex type."""
    if tensor_width == 64:
        return np.complex64
    return np.complex128
