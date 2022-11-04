from lib import state, ops


def bell_state(a: int, b: int) -> state.State:
    """Make one of the four bell states with a, b from {0,1}."""

    if a not in [0, 1] or b not in [0, 1]:
        raise ValueError("Bell state arguments are bits and must be 0 or 1.")
    psi = state.bitstring(a, b)
    psi = ops.Hadamard()(psi)
    return ops.Cnot()(psi)


def ghz_state(nbits: int) -> state.State:
    """Make a maximally entangled nbits state (GHZ State)."""
    # Simple construction via:
    #
    # |0> --- H --- o ---------    --- H --- o ----o----
    # |0> ----------X --- o --- or ----------X --- | ---
    # |0> ----------------X ---    ----------------X ---
    #
    psi = state.zeros(nbits)
    psi = ops.Hadamard()(psi)
    for offset in range(nbits - 1):
        psi = ops.Cnot(0, 1)(psi, offset)
    return psi

