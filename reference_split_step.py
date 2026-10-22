"""Independent NumPy baseline for periodic 1D Schrödinger evolution.

Not a verification of the historical Qiskit circuits. Grid amplitudes have
sum(abs(psi)**2) = 1. Units default to m = hbar = 1.
"""
import numpy as np


def normalize_state(amplitudes):
    state = np.asarray(amplitudes, dtype=complex)
    norm = np.linalg.norm(state)
    if not np.isfinite(norm) or norm == 0:
        raise ValueError("State must have a finite, nonzero norm")
    return state / norm


def gaussian_packet(x, x0, sigma, k0):
    """Gaussian with probability-density std. dev. sigma, phase exp(i*k0*x)."""
    if sigma <= 0:
        raise ValueError("sigma must be positive")
    x = np.asarray(x, dtype=float)
    return normalize_state(np.exp(-((x - x0)**2) / (4 * sigma**2) + 1j*k0*x))


def split_step(psi, x, potential, dt, steps, *, mass=1.0, hbar=1.0):
    """Second-order Strang split on a uniform periodic grid, endpoint excluded."""
    state = np.asarray(psi, dtype=complex).copy()
    x = np.asarray(x, dtype=float)
    v = np.asarray(potential, dtype=float)
    if state.ndim != 1 or len(state) < 2 or x.shape != state.shape or v.shape != state.shape:
        raise ValueError("psi, x, potential must be equally sized 1D arrays")
    if not np.isfinite(state).all() or not np.isfinite(x).all() or not np.isfinite(v).all():
        raise ValueError("Inputs must be finite")
    dx = x[1] - x[0]
    if dx <= 0 or not np.allclose(np.diff(x), dx, rtol=1e-10, atol=1e-12):
        raise ValueError("x must be uniformly increasing")
    if mass <= 0 or hbar <= 0:
        raise ValueError("mass and hbar must be positive")
    if not isinstance(steps, (int, np.integer)) or steps < 0:
        raise ValueError("steps must be a nonnegative integer")
    k = 2*np.pi*np.fft.fftfreq(len(state), d=dx)
    kinetic = np.exp(-1j*dt*hbar*k**2/(2*mass))
    potential_half = np.exp(-0.5j*dt*v/hbar)
    for _ in range(steps):
        state = potential_half * state
        state = np.fft.ifft(kinetic*np.fft.fft(state))
        state = potential_half * state
    return state
