"""Qiskit Statevector checks against the independent NumPy reference."""
import unittest
import numpy as np

from qiskit_split_step import (periodic_grid, qft_split_step_circuit,
                               statevector_position_and_work, phase_aligned_errors)
from reference_split_step import gaussian_packet, split_step


class QiskitSplitStepTests(unittest.TestCase):
    def _one_step(self, potential, *, n=4, dt=0.03, k0=1.5):
        x = periodic_grid(n)
        psi = gaussian_packet(x, -1.5, 0.55, k0)
        circuit, position, work = qft_split_step_circuit(psi, x, potential(x), dt)
        actual, work_one_probability = statevector_position_and_work(circuit, position, work)
        expected = split_step(psi, x, potential(x), dt, 1)
        return actual, expected, work_one_probability

    def test_free_one_step_small_registers(self):
        for n in (3, 4, 5):
            actual, expected, work = self._one_step(lambda x: np.zeros_like(x), n=n)
            amplitude, probability = phase_aligned_errors(actual, expected)
            self.assertLess(amplitude, 2e-14)
            self.assertLess(probability, 2e-14)
            self.assertLess(work, 1e-15)

    def test_constant_potential_phase(self):
        x = periodic_grid(4)
        psi = gaussian_packet(x, -1, 0.6, 1.0)
        dt, value = 0.04, 2.5
        circuit, position, work = qft_split_step_circuit(psi, x, np.full(len(x), value), dt)
        actual, leakage = statevector_position_and_work(circuit, position, work)
        free = split_step(psi, x, np.zeros_like(x), dt, 1)
        np.testing.assert_allclose(actual, free * np.exp(-1j * value * dt), atol=2e-14)
        self.assertLess(leakage, 1e-15)

    def test_step_and_barrier_masks_match_reference(self):
        step = lambda x: np.where(x >= 0, 1.5, 0.0)
        barrier = lambda x: np.where((x >= -0.5) & (x < 0.5), 3.0, 0.0)
        for potential in (step, barrier):
            actual, expected, work = self._one_step(potential)
            amplitude, probability = phase_aligned_errors(actual, expected)
            self.assertLess(amplitude, 2e-14)
            self.assertLess(probability, 2e-14)
            self.assertLess(work, 1e-15)

    def test_norm_preserved(self):
        actual, _, _ = self._one_step(lambda x: np.where(x >= 0, 1.5, 0.0))
        self.assertAlmostEqual(float(np.sum(np.abs(actual) ** 2)), 1.0, places=14)

    def test_multistep_qft_matches_numpy(self):
        """Eight coherent Strang steps for free, step and barrier potentials."""
        for n in (3, 5):
            x = periodic_grid(n)
            psi = gaussian_packet(x, -0.5, 0.6, 2.0)
            for v in (
                np.zeros_like(x),
                np.where(x >= 0, 1.5, 0.0),
                np.where((x >= -0.5) & (x < 0.5), 3.0, 0.0),
            ):
                qc, position, work = qft_split_step_circuit(psi, x, v, 0.04, steps=8)
                actual, leakage = statevector_position_and_work(qc, position, work)
                expected = split_step(psi, x, v, 0.04, 8)
                amp_error, probability_error = phase_aligned_errors(actual, expected)
                self.assertLess(amp_error, 5e-12)
                self.assertLess(probability_error, 5e-12)
                self.assertAlmostEqual(float(sum(abs(actual)**2)), 1.0, places=12)
                self.assertLess(leakage, 1e-15)

    def test_optional_work_qubit_and_invalid_step_count(self):
        x = periodic_grid(3)
        psi = gaussian_packet(x, -0.5, 0.6, 1.0)
        v = np.where(x >= 0, 1.5, 0.0)
        qc, position, work = qft_split_step_circuit(
            psi, x, v, 0.04, include_work_qubit=False, steps=3
        )
        self.assertIsNone(work)
        actual, leakage = statevector_position_and_work(qc, position, work)
        expected = split_step(psi, x, v, 0.04, 3)
        self.assertLess(phase_aligned_errors(actual, expected)[0], 5e-12)
        self.assertEqual(leakage, 0.0)
        with self.assertRaises(ValueError):
            qft_split_step_circuit(psi, x, v, 0.04, steps=0)
