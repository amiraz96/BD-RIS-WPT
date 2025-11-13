# Beyond-Diagonal RIS Wireless Power Transfer (BD-RIS-WPT)

This repository collects MATLAB scripts used to simulate and optimize reconfigurable intelligent surface (RIS) architectures for wireless power transfer experiments associated with the preprint [arXiv:2502.19176](https://arxiv.org/abs/2502.19176). The code explores iterative and convex-optimization-based beamforming strategies for beyond-diagonal RIS (BD-RIS) structures across Rayleigh and Rician fading channels.

## Requirements

- **MATLAB** (tested with recent releases). The scripts rely on base MATLAB functionality and complex arithmetic.
- **CVX** (with an SDP-capable backend such as MOSEK) for semidefinite programs solved inside the optimization scripts. Activate the solver by uncommenting the `cvx_solver mosek;` lines when MOSEK is available.

## Repository structure

| File | Description |
| --- | --- |
| `MONTE_CARLO_ITWF.m` | Monte Carlo driver for the iterative transmit waveform (ITWF) algorithm. Alternates between RIS impedance updates and waveform re-optimization under frequency-selective Rayleigh fading, subject to power-transfer metrics with polynomial diode coefficients `K2` and `K4`. |
| `Monte_CARLO_OPT_DiagRIS.m` | Evaluates a diagonal-RIS benchmark with successive convex approximation and SDP relaxation under Rician channels. Includes joint updates of waveform amplitudes and RIS phases using CVX and Takagi-based randomization. |
| `Monte_CARLO_OPT_SDRBDRIS.m` | Implements a semidefinite relaxation with Gaussian randomization for beyond-diagonal RIS designs, iteratively improving waveform and impedance variables. |
| `Monte_CARLO_OPT_SDP_BDRIS.m` | Similar to the SDR script but enforces additional nuclear-norm constraints to promote low-rank structure during the SDP iterations. |
| `Rayleigh_channel.m`, `Rayleigh_channel_flat.m`, `Rayleigh_channel_time.m` | Channel generators producing frequency-flat or frequency-selective Rayleigh fading responses (the time-varying version derives taps from coherence bandwidth and delay-spread parameters). |
| `Rician_channel.m` | Builds Rician fading matrices with user-defined K-factor and deterministic line-of-sight steering based on planar-array geometry. |
| `aprox_fixed_w.m`, `aprox_fixed_w_new.m` | Evaluate fourth-order diode-based harvested power approximations and their gradients with respect to transmit weights, for complex and magnitude-only variants respectively. |
| `takagi_phase_randomization_approximation.m` | Searches random phase rotations of a symmetric matrix’s Takagi factors to maximize harvested power metrics. |
| `gaussian_randomization_cholesky.m` | Generates candidate rank-one factorizations from an SDP solution by applying Cholesky-based Gaussian randomization. |
| `Rect_K.m` | Utility returning rectifier coefficients used in diode-based power models. |

## Running the simulations

1. Add the repository folder (and any required helper utilities) to your MATLAB path.
2. Ensure CVX is installed and initialized (`cvx_setup`). Uncomment solver selection lines if you have licensed solvers such as MOSEK.
3. Adjust scenario parameters inside the desired Monte Carlo script (number of RIS elements `N`, subcarriers `N_f`, fading parameters, etc.).
4. Execute the script from MATLAB. Output `.mat` files store optimized waveforms, impedance matrices, and harvested power statistics for each random channel realization.

The channel-generation functions and optimization routines can also be reused independently to prototype new RIS architectures or evaluate alternative waveform optimization strategies.
