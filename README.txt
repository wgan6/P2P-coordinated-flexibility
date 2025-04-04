# Coordinated Flexibility via P2P – Impact on LVDN Network Flows

This MATLAB project simulates the impact of peer-to-peer (P2P) coordinated flexibility on low-voltage distribution networks under UK future energy scenarios (UK FES).

## Contents

- `CoordinatedFlex_viaP2P_ImpactOnNetworks.m`: Main simulation and optimization script  
- `Case_Parameters_from_UK_FES.m`: Scenario parameters (PV, load, EV, DR)  
- `gendist.m`: Helper function for stochastic distributions  
- `*.mat`: Network data, PV/load profiles

## Usage

1. Open the main script in MATLAB.
2. Ensure YALMIP and MOSEK (or another compatible solver) are installed.
3. Run the simulation. Results will be saved as `.mat` files.

## Notes

- Statistically network cases: `SSNG_Keydata_network_1/2/3`
- PV/load profiles: Typical UK residential datasets








