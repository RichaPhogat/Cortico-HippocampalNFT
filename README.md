# Code for: A unifying framework for healthy and pathological cortico-hippocampal interactions

This repository contains the MATLAB code used to generate the results and figures for the manuscript:

**“A unifying framework for healthy and pathological cortico-hippocampal interactions”**

The code implements neural field models of cortico-hippocampal interactions as proposed in the manuscript. 
---

## Repository Structure

The repository is organized by figure and experiment for clarity and reproducibility:

- **`UncoupledCortexHippo/`**  
  Simulations of uncoupled cortico-thalamic and hippocampo-septal systems. These simulations give rise to the cortical alpha and hippocampal theta rhythm dynamics.

- **`Figure1/`**  
  Scripts used to generate all panels of Figure 1 (cortico-thalamic system dynamics).

- **`Figure2/`**  
  Scripts used to generate all panels of Figure 2 (hippocampo-septal system dynamics).

- **`Figure3/`**  
  Scripts used to generate all panels of Figure 3 (coupling the cortico-hippocampal dynamics).

- **`SeizureP1/`, `SeizureP2/`, `SeizureP3/`**  
  Patient-specific seizure modeling and parameter fitting scripts.

Each folder contains a dedicated `README.md` describing:
- the purpose of the folder,
- the main entry-point scripts,
- required inputs,
- generated outputs.

---

## Requirements

- MATLAB (tested on R2022b and later)

---

## Running the Code

1. Clone the repository:
   ```bash
   git clone https://github.com/RichaPhogat/Cortico-HippocampalNFT.git
