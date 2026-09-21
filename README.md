# Multi-Satellite Multi-Connectivity Beam-Hopping Scheduling for NTN D2C Communications

This repository contains the simulation code for multi-satellite beam scheduling in NTN direct-to-cell networks.

## Related Works

- **[1]** *Multi-Satellite Multi-Connectivity Beam-Hopping Scheduling for NTN D2C Communications*. IEEE Transactions on Wireless Communications (Major Revision).
- **[2]** *Beam Scheduling for Multi-Connectivity NTN LEO D2C Satellite Networks*. IEEE VTC 2026-Spring.
- **[3]** Z. Lin, Z. Ni, L. Kuang, C. Jiang, and Z. Huang, “Multi-satellite beam hopping based on load balancing and interference avoidance for NGSO satellite communication systems,” *IEEE Transactions on Communications*, vol. 71, no. 1, pp. 282–295, 2023.
- **[4]** M. Zhao, N. Ye, Q. Ouyang, Y. Jin, Y. Jin, and L. Zhao, “Multi-satellite cooperative communication: Exploiting time asynchrony in non-orthogonal transmissions,” *IEEE Transactions on Vehicular Technology*, vol. 72, no. 5, pp. 6868–6873, 2023.
- **[5]** X. Zhang, S. Sun, M. Tao, Q. Huang, and X. Tang, “Multi-satellite cooperative networks: Joint hybrid beamforming and user scheduling design,” *IEEE Transactions on Wireless Communications*, vol. 23, no. 7, pp. 7938–7952, 2024.

## Repository Structure

```text
Project_root/
├── data/                         # External and generated scenario data
├── figures/                      # Generated figures
├── matlab/
│   ├── main/                     # Simulation scripts
│   ├── module/                   # Scheduling algorithms
│   └── util/                     # Utility functions
├── python/
│   ├── scenario_generation/      # Scenario generation
│   └── plotting_scripts/         # Result visualization
└── results/                      # Generated simulation results
```

## Environment

- Python 3.10
- MATLAB R2024a
- Gurobi Optimizer 12.0.0 with a valid license

Install the Python dependencies:

```bash
pip install numpy pandas matplotlib shapely pyproj cartopy h3 folium geopandas fiona astropy "poliastro==0.17.0"
```

Configure Gurobi for MATLAB and replace `YOUR_GUROBI_MATLAB_PATH` in the MATLAB modules with your local Gurobi MATLAB path.

## Workflow

### 1. Generate Scenarios

Download [Kontur Population for China](https://geodata-eu-central-1-kontur-public.s3.amazonaws.com/kontur_datasets/kontur_population_CN_20220630.gpkg.gz), extract it, and place `kontur_population_CN_20220630.gpkg` in `data/`.

```bash
cd python/scenario_generation
python generate_cell.py
python relativeDensity.py
```

In `visibility_generate.py`, set `constellation_config_mode` to `1` and `2`, running the script once for each configuration:

```bash
python visibility_generate.py
```

### 2. Run MATLAB Simulations

In MATLAB:

```matlab
cd matlab/main
addpath('../module', '../util')
fig1to6
fig7to8
fig9to10
fig11to12
fig13to14
fig15
fig16
fig17
fig18
```

Results are saved in `results/`.

### 3. Generate Figures

```bash
cd python/plotting_scripts
python figureFile.py
```

Figures are saved in `figures/`.

## Usage

This code is provided solely for academic evaluation and reproducibility.
