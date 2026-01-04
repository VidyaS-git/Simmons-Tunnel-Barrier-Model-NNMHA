# Simmons-Tunnel-Barrier-Model-NNMHA
NGSpice-based implementation of a Simmons/Pickett tunnel-barrier memristor model with Python scripts for IV, potentiation/depression, and dynamic route map (DRM) characterization. Designed for reproducible education-oriented experiments.

# Simmons / Pickett Tunnel-Barrier Memristor (NGSpice + Python)

This repository contains an educational, reproducible workflow to simulate and characterize a memristor model in **NGSpice**, with **Python** used only for execution and post-processing.

The scripts generate the three standard characterization experiments required in the voluntary modeling project:
1. **I–V curve** under a triangular stimulus  
2. **Potentiation / Depression** via write–read pulse trains  
3. **Dynamic Route Map (DRM)** showing state dynamics

> Note: A PySpice-based implementation exists, but on Windows it may fail due to `ngspice.dll` loading issues (Python ABI / VC runtime). Therefore, the default workflow here uses **NGSpice batch mode** for maximum robustness and reproducibility.

---

## Repository structure

- `mem_simmons_pickett.subckt`  
  NGSpice subcircuit implementing the memristor behavior.
- `simulate_SimmonsMemristor.py`  
  Runs NGSpice in batch mode, exports waveforms to CSV, and generates plots.
- `iv_run.cir` *(optional / generated)*  
  Example netlist for the triangular I–V experiment.
- `iv.csv`, `pulses.csv` *(generated)*  
  Exported simulation results.
- `results/` *(generated)*  
  Plots (`fig_1.png`, `fig_2.png`, …) and any saved output.

---

## Requirements

- **NGSpice** (recommended v38+; tested with v45.x)
- **Python 3.9+**
- Python packages:
  - `numpy`
  - `matplotlib`

Install Python packages:
```bash
pip install numpy matplotlib
