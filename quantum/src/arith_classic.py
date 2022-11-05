import math

import numpy as np

from lib import state, ops


def fulladder_matrix(psi: state.State):
    """Non-quantum-exploiting, classic full adder."""

    psi = ops.Cnot(0, 3)(psi, 0)
    psi = ops.Cnot(1, 3)(psi, 1)
    psi = ops.ControlledU(0, 1, ops.Cnot(1, 4))(psi, 0)
    psi = ops.ControlledU(0, 2, ops.Cnot(2, 4))(psi, 0)
    psi = ops.ControlledU(1, 2, ops.Cnot(2, 4))(psi, 1)
    psi = ops.Cnot(2, 3)(psi, 2)
    return psi


def experiment_matrix(a: int, b: int, cin: int,
                      expected_sum: int, expected_count: int):
    """Run a simple classic experiment, check results."""
    psi = state.bitstring(a, b, cin, 0, 0)
    psi = fulladder_matrix(psi)

    bsum, _ = ops.Measure(psi, 3, tostate=1, collapse=False)
    bout, _ = ops.Measure(psi, 4, tostate=1, collapse=False)
    print(f'a: {a} b:{b} cin: {cin} sum: {bsum} cout: {bout}')
    if bsum != expected_sum or bout != expected_count:
        raise AssertionError("invalid results")


def add_classic():
    """Full eval of the full adder."""

    for exp_function in [experiment_matrix]:
        exp_function(0, 0, 0, 0, 0)
        exp_function(0, 1, 0, 1, 0)
        exp_function(1, 0, 0, 1, 0)
        exp_function(1, 1, 0, 0, 1)


def run_experiment(a1: np.complexfloating, a2: np.complexfloating,
                   target: float) -> None:
    """Construct swap test circuit and measure"""

    # |0> --- H --- o --- H --- Measure
    #               |
    # a1 ---------- x ---------
    #               |
    # a2 ---------- x ---------
    psi = state.bitstring(0) * state.qubit(a1) * state.qubit(a2)
    psi = ops.Hadamard()(psi, 0)
    psi = ops.ControlledU(0, 1, ops.Swap(1, 2))(psi)
    psi = ops.Hadamard()(psi, 0)

    # Measure qubit 0 once.
    p0, _ = ops.Measure(psi, 0)

    if abs(p0 - target) > 0.05:
        raise AssertionError("Probability {:.2f} off more than 5 pct from target {:.2f}".format(p0, target))
    print("Similarity of a1: {:.2f}, a2: {:.2f} ==> \%: {:.2f}".format(a1, a2, 100.0 * p0))


def alice_measures(alice: state.State, expect0: np.complexfloating, expect1: np.complexfloating,
                   qubit0: np.complexfloating, qubit1: np.complexfloating):
    """Force measure and get teleported qubit."""

    # Alices measure her state and gets a collapsed |qubit0 qubit1>
    # She let's Bob know which one of the 4 combinations she obtained.
    # We force measurement here, collapsing to a state with the
    # first two qubits collapsed. Bob's qubit is still unmeasured.
    _, alice0 = ops.Measure(alice, 0, tostate=qubit0)
    _, alice1 = ops.Measure(alice0, 1, tostate=qubit1)

    # Depending on what was measured and communicated, Bob has to do
    # one of these things to his qubit2:
    if qubit0 == 0 and qubit1 == 0:
        pass
    if qubit0 == 0 and qubit1 == 1:
        alice1 = ops.PauliX()(alice1, idx=2)
    if qubit0 == 1 and qubit1 == 0:
        alice1 = ops.PauliZ()(alice1, idx=2)
    if qubit0 == 1 and qubit1 == 1:
        alice1 = ops.PauliX()(ops.PauliZ()(alice1, idx=2), idx=2)

    # Now Bob measures his qubit (2) (without collapse, so we can 'measure' it twice.
    # This is not necessary, but good to double-check the maths).
    p0, _ = ops.Measure(alice1, 2, tostate=0, collapse=False)
    p1, _ = ops.Measure(alice1, 2, tostate=1, collapse=False)

    # Alice should now have 'teleported' the qubit in state 'x'.
    # We sqrt() the probability, we want to show (original) amplitudes.
    bob_a = math.sqrt(p0.real)
    bob_b = math.sqrt(p1.real)

    print("Teleported (|{:d}{:d}>)   a={:.2f}, b={:.2f}".format(
        qubit0, qubit1, bob_a, bob_b))

    if not math.isclose(expect0, bob_a, abs_tol=1e-6) or not math.isclose(expect1, bob_b, abs_tol=1e-6):
        raise AssertionError("Invalid result.")
