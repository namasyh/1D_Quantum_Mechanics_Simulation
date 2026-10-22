# Preliminary physics audit (2026-10-22)

Static inspection of `Step_potential.ipynb` and `quantum_tunnelling.ipynb`. The very large free-particle notebook was not fully retrieved, inspected, or executed. This is **not** Qiskit statevector validation. Original notebooks and figures remain untouched.

1. **Aliasing:** Both inspected notebooks use 64 points from -0.5 to 4 inclusive and `exp(300j*x)`. Thus dx = 4.5/63 ≈ 0.07143, Nyquist k = π/dx ≈ 43.98, and k=300 aliases to approximately 36.11. Fix grid and physical parameters together.
2. **Gaussian width:** The actual amplitude envelope is `exp[-0.5*(x-x0)^2/0.016]`. Its amplitude-envelope Gaussian scale is sqrt(0.016) ≈ 0.1265, corresponding to probability-density standard deviation sqrt(0.008) ≈ 0.0894 on a broad grid—not the claimed sigma=0.4. The discrete normalization is mathematically reasonable.
3. **Step potential:** The step code uses `QuantumCircuit(aux_reg, q_reg, clas_reg)` and `qc.p(V*phi,4)`. Circuit qubit 4 is position bit 3, so the phase affects indices 8–15, 24–31, 40–47, 56–63. The figure instead marks a single step at index 16. Positive phase also conflicts with `exp(-i V dt / hbar)` for positive V and dt.
4. **Finite barrier:** `qc.cx(4,0); qc.cp(V*phi,0,5); qc.cx(4,0)` phases indices 24–31 and 56–63, rather than just the depicted interval around 23–32. The notebook's global `V=1.7` is also superseded by local `V=5`.
5. **Kinetic operator:** Positive phase gates in `Time_Evolution`, a global `phi`, QFT/inverse-QFT ordering, X gates around circuit index `n`, signed momentum interpretation, and ancilla uncomputation require explicit statevector checks. Histogram appearances cannot verify these.
6. **Splitting/boundaries:** The loop structurally applies a first-order product of potential and kinetic steps; the *independent NumPy reference* here uses second-order Strang splitting. Fourier grids impose periodic boundaries and can cause wavepacket wraparound.

## Bounded local validation tasks

- Set explicit `m`, `hbar`, `dt`, initial center/width/momentum, domain, and a periodic grid with `endpoint=False`.
- Build numerical step/barrier potential arrays first; make circuit phase masks and plot overlays match exactly.
- Compare **unmeasured Qiskit statevectors** (3–6 position qubits, including ancilla) against the reference for free, constant, step and barrier potentials. Check amplitude and probability error, global phase, register indexing and ancilla reset.
- Reconcile positive/negative phase signs, QFT direction, signed momentum, numerical convergence, boundary wraparound and energy relative to barrier height.
- Only then regenerate representative figures and update scientific claims. Probability beyond a barrier alone does not establish tunneling.

The added unit tests check only the standalone NumPy reference and do not show Qiskit agreement, hardware results, or quantum advantage.
