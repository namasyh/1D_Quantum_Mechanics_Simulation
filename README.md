# 1D Quantum Mechanics Simulation

A self-directed educational exploration of 1D, single-particle Schrödinger evolution using Qiskit circuits. The notebooks study free Gaussian wavepackets, a step potential, and a finite barrier. **QFT/split-operator methods are established approaches; this is not a new simulation algorithm.**

> **Physics status:** Corrected notebook circuits now pass local Qiskit Statevector checks against the independent NumPy split-step reference. This validates the chosen finite, periodic-grid discretization—not continuum accuracy, hardware execution, or quantum advantage. See [the audit](docs/PHYSICS_AUDIT.md).

## Model and numerical method

The time-independent Hamiltonian is `H = T + V`, with `T = -ħ²∂²/(2m∂x²)`. All corrected examples set `m=ħ=1`, use a periodic grid `x_j∈[-L/2,L/2)` with `N=2^n`, and normalize amplitudes with `Σ_j|ψ_j|²=1`. `sigma` is the probability-density standard deviation and the travelling factor is `exp(i k0 x)`. Fourier labels are `2π fftfreq(N, dx)`, so the upper half encodes negative wave numbers. Evolution is the symmetric **second-order Strang split** `exp(-iVdt/2ħ) F⁻¹ exp(-iTdt/ħ) F exp(-iVdt/2ħ)`. Potential and kinetic diagonal phases both carry the physical negative sign for positive time.

## Corrected notebooks and historical originals

| Case | Notebook | Status |
|---|---|---|
| Free packet | [time_evolution_free_particle.ipynb](Notebooks/time_evolution_free_particle.ipynb) | One-step validation and new multi-step comparison example |
| Step potential | [Step_potential.ipynb](Notebooks/Step_potential.ipynb) | Exact `x >= 0` phase mask and plotted mask |
| Finite barrier | [quantum_tunnelling.ipynb](Notebooks/quantum_tunnelling.ipynb) | Exact `-0.5 <= x < 0.5` mask; no claim of confirmed sub-barrier tunneling |

The [archive/Notebooks](archive/Notebooks/) directory preserves byte-for-byte copies of all three historical notebooks from the original February 2026 main branch (including their large, **unvalidated** outputs). The corrected `Notebooks/` versions are the supported examples. See [archive/README.md](archive/README.md).

Preserved original free-particle illustrations (historical, not proof of physical accuracy):

| Initial stationary | Later stationary | Initial travelling | Later travelling |
|:---:|:---:|:---:|:---:|
| ![](images/image.png) | ![](images/image-1.png) | ![](images/image-2.png) | ![](images/image-3.png) |

## Validation and reproduction

NumPy alone suffices for the baseline reference, but the complete test suite requires Qiskit. Dependencies target Qiskit 2.6 and do not require Aer, IBM credentials or real hardware.

```bash
python -m venv .venv
# Activate: source .venv/bin/activate (macOS/Linux) or .venv\Scripts\Activate.ps1 (PowerShell)
python -m pip install -r requirements.txt
python -m unittest discover -s tests -v
```

Tested locally in the existing `qiskit-qc` Conda environment: Python 3.14.6, NumPy 2.5.1, Matplotlib 3.11.1, Qiskit 2.6.0, and ipykernel 7.2.0. `python -m pip check` reported no broken requirements. The independent [reference_split_step.py](reference_split_step.py) required no physics change: it already uses the selected endpoint-excluded FFT grid, negative phase, and Strang split. [qiskit_split_step.py](qiskit_split_step.py) implements the same one-step operator using QFT and explicit diagonal gates; it contains no IBM credentials, sampling, or arithmetic ancilla.

The local `qiskit-qc` validation pass reported **10 passing tests** with `python -m unittest discover -s tests -v`: analytic free plane wave, constant potential, normalization/reversal, free QFT steps for 3–5 qubits, one- and eight-step free/step/barrier QFT comparisons, and optional-work-qubit behavior. In the local 5-qubit one-step comparisons (`dt=0.03`), maximum global-phase-aligned amplitude errors were `6.3e-16`–`7.1e-16`; maximum probability errors were `2.5e-16`–`3.3e-16`; norm was 1 to floating-point precision; work-qubit `|1>` probability was 0. Eight-step notebook comparisons produced maximum amplitude errors below `5.1e-15` and probability errors below `2.6e-15`. These values compare two implementations of the same finite-grid operator, not the continuum solution.

For a 256-point barrier check (`L=16`, `x0=-4`, `sigma=0.5`, `k0=2.5`, `V0=6`, width 1, `t=2`), reducing `dt` from 0.02 to 0.01 changed maximum amplitude by `5.86e-3` and probability by `1.57e-3`; 0.01 to 0.005 changed them by `1.37e-3` and `2.10e-4`. This is expected product-formula refinement, not a proof of convergence in space or time. The packet has mean discrete kinetic energy 3.625 but about 0.200 of its discrete momentum probability lies above `V0`; observed right-side probability therefore cannot be identified as purely sub-barrier tunnelling. Periodic boundaries also wrap outgoing probability: a free 64-point `L=8` packet centred at `-2.5` with `k0=3` has its modal index return to `x=-0.625` by `t=2.5` after reaching the right side.

**Follow-up checks:** Two additional test methods and matching comparison cells in all three notebooks exercise multiple coherent Qiskit Statevector steps (including free, step and barrier masks) versus the NumPy reference. The 10-test suite and all corrected notebook code cells passed in `qiskit-qc`; there are no GitHub CI checks yet. The large wavepacket plots still use NumPy; the separately labeled small-register comparison plots use multi-step Qiskit evolution.

For an existing Miniconda environment rather than a new virtual environment:

```bash
conda activate qiskit-qc
python -m pip install -r requirements.txt
python -m unittest discover -s tests -v
python -m ipykernel install --user --name qiskit-qc --display-name "Python (qiskit-qc)"
```

**Limitations:** The QFT method represents a periodic, spectral discretization. Avoid interpreting late-time edge population after wraparound. Diagonal gates are locally statevector-simulated here; no compilation cost, noise, shots, or real-hardware behavior has been validated. Transmission beyond the barrier is not by itself proof of quantum tunnelling, especially for a packet with above-barrier spectral weight. No quantum-advantage, novel-algorithm, publication, or hardware-execution claims are made.

## Attribution

- Giuliano Benenti and Giuliano Strini, [*Quantum simulation of the single-particle Schrödinger equation*](https://arxiv.org/abs/0709.1704) (2008): established quantum simulation methods.
- Mohamed Abouelela, [*Quantum Simulation of the Schrödinger Equation using IBM's Quantum Computers*](https://fount.aucegypt.edu/capstone/20/) (2020): related Qiskit educational exploration.
- Lohit Potnuru, *Simulate Quantum Particles with Python and Qiskit* (2022): educational influence; citation retained from the original notebooks.

This project represents independent study/adaptation of pre-existing methods to illustrative 1D potentials. Further convergence analysis and physical scattering/tunneling observables remain future work; statevector equivalence for a finite discretization is not a continuum validation.
