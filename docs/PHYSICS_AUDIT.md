# Preliminary physics audit (2026-10-22)

Initial static inspection of all three historical notebooks, followed by bounded local Qiskit Statevector validation of their corrected replacements. No hardware runs occurred.

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

## Correction and validation outcome

`qiskit_split_step.py` fixes the selected convention: `m=hbar=1`, `N=2^n`, `x in [-L/2,L/2)`, FFT-ordered signed wave numbers, and a Strang step with negative potential and kinetic phases. Qiskit's position-register bit order is explicitly treated as the little-endian integer index `j`. The QFT kinetic diagonal is even in signed `k`, so QFT sign reversal does not alter `k^2`. The former arithmetic-ancilla kinetic network was replaced by an explicit diagonal, avoiding ambiguous phase decomposition; an optional untouched work qubit has `P(|1>)=0` in every Statevector check.

Local tests compare Qiskit Statevector with the independent NumPy reference for free registers of 3, 4, and 5 qubits, constant potential, a step `x >= 0`, and a barrier `-0.5 <= x < 0.5`. All eight tests pass. At 5 qubits and `dt=0.03`, maximum phase-aligned amplitude error is below `7.1e-16`; probability error is below `3.4e-16`; norm remains one to floating-point precision. These are operator-equivalence checks for the discretized model, not a continuum error bound.

The reference solver was retained: it already implements endpoint-excluded periodic FFT Strang evolution with the same units and phase convention. A barrier refinement study at 256 points showed decreasing fine-step differences, but does not establish full convergence. Finite periodic domains wrap; plots and interpretation must remain before return flux. For the illustrated barrier packet, substantial above-barrier spectral weight remains, so right-side probability is numerical transmission containing above-barrier components, not a clean tunnelling measurement.


## Follow-up: multi-step statevector comparisons

The Qiskit helper accepts `steps`, implementing a single initialization followed by repeated unitary Strang steps. Each corrected notebook now includes a small-register plot comparing **multi-step Qiskit** and **multi-step NumPy** probabilities. The large illustrative wavepacket trajectories continue to use the NumPy reference and must not be described as long Qiskit-circuit output.

Two additional test methods cover 8-step statevector comparisons for free, step and barrier potentials and optional work-qubit behavior. In the existing `qiskit-qc` environment (Qiskit 2.6.0), all 10 tests passed; notebook multi-step comparisons reached maximum amplitude errors below `5.1e-15`, probability errors below `2.6e-15`, and zero work-qubit leakage.

Three original February 2026 notebooks are retained unchanged in `archive/Notebooks/`, including large historical outputs, for provenance rather than execution.
