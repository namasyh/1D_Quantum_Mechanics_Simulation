"""Statevector-friendly QFT split-step circuits for the educational notebooks.

The position register uses Qiskit's little-endian integer convention: basis
state ``|j>`` represents ``x[j]``.  The Fourier labels use
``2*pi*fftfreq(N, dx)``; consequently their high half represents negative
wave numbers.  The unused work qubit is retained only to make its reset state
explicit during validation--the diagonal construction needs no arithmetic
ancilla.
"""
from __future__ import annotations

import numpy as np
from qiskit import QuantumCircuit, QuantumRegister
from qiskit.circuit.library import DiagonalGate, QFTGate
from qiskit.quantum_info import Statevector


def periodic_grid(num_qubits: int, length: float = 8.0) -> np.ndarray:
    """Uniform periodic grid on [-length/2, length/2), with N=2**num_qubits."""
    if num_qubits < 1 or length <= 0:
        raise ValueError("num_qubits and length must be positive")
    return np.linspace(-length / 2, length / 2, 2**num_qubits, endpoint=False)


def wavenumbers(x: np.ndarray) -> np.ndarray:
    """FFT-ordered wave numbers, including the negative-frequency half."""
    dx = x[1] - x[0]
    return 2 * np.pi * np.fft.fftfreq(len(x), d=dx)


def _diagonal(circuit: QuantumCircuit, values: np.ndarray, position: QuantumRegister) -> None:
    circuit.append(DiagonalGate(np.asarray(values, dtype=complex)), position)


def qft_split_step_circuit(
    psi: np.ndarray,
    x: np.ndarray,
    potential: np.ndarray,
    dt: float,
    *,
    mass: float = 1.0,
    hbar: float = 1.0,
    include_work_qubit: bool = True,
    steps: int = 1,
) -> tuple[QuantumCircuit, QuantumRegister, QuantumRegister | None]:
    """One or more second-order Strang steps, ``V(dt/2) F^-1 T(dt) F V(dt/2)``.

    QFTGate's sign convention is the conjugate of NumPy's ``fft``.  Because
    ``k**2`` is even, the kinetic diagonal is identical under that reversal.
    """
    psi, x, potential = map(np.asarray, (psi, x, potential))
    n = int(np.log2(len(psi)))
    if 2**n != len(psi) or x.shape != psi.shape or potential.shape != psi.shape:
        raise ValueError("psi, x and potential require the same power-of-two length")
    if mass <= 0 or hbar <= 0:
        raise ValueError("mass and hbar must be positive")
    if not isinstance(steps, (int, np.integer)) or steps < 1:
        raise ValueError("steps must be a positive integer")
    if not np.isfinite(psi).all() or not np.isfinite(x).all() or not np.isfinite(potential).all():
        raise ValueError("inputs must be finite")
    if np.linalg.norm(psi) == 0:
        raise ValueError("initial state must have nonzero norm")
    if len(x) < 2 or x[1] <= x[0] or not np.allclose(np.diff(x), x[1] - x[0]):
        raise ValueError("x must be uniformly increasing")
    position = QuantumRegister(n, "x")
    work = QuantumRegister(1, "work") if include_work_qubit else None
    circuit = QuantumCircuit(position, work) if work else QuantumCircuit(position)
    circuit.initialize(psi / np.linalg.norm(psi), position)
    potential_phase = np.exp(-0.5j * dt * potential / hbar)
    kinetic_phase = np.exp(-0.5j * dt * hbar * wavenumbers(x) ** 2 / mass)
    # Repeated coherent unitary evolution, with only one state initialization.
    for _ in range(steps):
        _diagonal(circuit, potential_phase, position)
        circuit.append(QFTGate(n), position)
        _diagonal(circuit, kinetic_phase, position)
        circuit.append(QFTGate(n).inverse(), position)
        _diagonal(circuit, potential_phase, position)
    return circuit, position, work


def statevector_position_and_work(
    circuit: QuantumCircuit, position: QuantumRegister, work: QuantumRegister | None
) -> tuple[np.ndarray, float]:
    """Return position amplitudes and probability that an optional work qubit is 1."""
    state = np.asarray(Statevector.from_instruction(circuit).data)
    if work is None:
        return state, 0.0
    size = 2 ** len(position)
    return state[:size], float(np.sum(np.abs(state[size:]) ** 2))


def phase_aligned_errors(actual: np.ndarray, expected: np.ndarray) -> tuple[float, float]:
    """Maximum amplitude and probability errors after removing global phase."""
    phase = np.vdot(expected, actual)
    aligned = actual * np.exp(-1j * np.angle(phase))
    return (float(np.max(np.abs(aligned - expected))),
            float(np.max(np.abs(np.abs(actual) ** 2 - np.abs(expected) ** 2))))
