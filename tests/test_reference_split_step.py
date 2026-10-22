"""Test only the NumPy baseline, not Qiskit circuit correctness."""
import unittest
import numpy as np
from reference_split_step import gaussian_packet, split_step


class ReferenceEvolutionTests(unittest.TestCase):
    def setUp(self):
        self.n = 64
        self.length = 8.0
        self.x = np.linspace(-4, 4, self.n, endpoint=False)

    def test_packet_normalization(self):
        psi = gaussian_packet(self.x, -1, 0.4, 3)
        self.assertAlmostEqual(float(np.sum(abs(psi)**2)), 1, places=13)

    def test_free_plane_wave_phase(self):
        k = 2*np.pi*3/self.length
        psi = np.exp(1j*k*self.x)/np.sqrt(self.n)
        dt = 0.23
        actual = split_step(psi, self.x, np.zeros(self.n), dt, 1)
        np.testing.assert_allclose(actual, psi*np.exp(-1j*k**2*dt/2), atol=1e-13)

    def test_constant_potential_phase(self):
        psi = gaussian_packet(self.x, 0, 0.7, 0)
        dt = 0.07
        actual = split_step(psi, self.x, np.full(self.n, 2.5), dt, 5)
        free = split_step(psi, self.x, np.zeros(self.n), dt, 5)
        np.testing.assert_allclose(actual, free*np.exp(-1j*2.5*dt*5), atol=1e-12)

    def test_barrier_norm_and_reversal(self):
        psi = gaussian_packet(self.x, -2, 0.5, 4)
        barrier = np.where((self.x >= 0) & (self.x < 1), 3.0, 0.0)
        evolved = split_step(psi, self.x, barrier, 0.02, 30)
        self.assertAlmostEqual(float(np.sum(abs(evolved)**2)), 1, places=12)
        reversed_psi = split_step(evolved, self.x, barrier, -0.02, 30)
        np.testing.assert_allclose(reversed_psi, psi, atol=1e-12)


if __name__ == "__main__":
    unittest.main()
