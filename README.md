# Multi-DOF Flexible Joint Robot: Stiffness & Inertia Matrix Identification

Frequency-domain identification of joint stiffness $k_i$ and the **signed** inverse inertia matrix $M_l^{-1}$ for multi-DOF flexible joint robots (FJR) — without requiring a prior CAD-based mass matrix.

## Key Contribution

The conventional diagonal-FRF analysis (Vieta's formulas) recovers $k_i$ and $|M_l^{-1}_{ij}|$, but cannot determine the **signs** of off-diagonal elements. We resolve this via:

1. **Eigenvalue constraint** — reduces $2^3 = 8$ sign candidates to 4 (Klein four-group symmetry)
2. **Off-diagonal FRF fitting** — uniquely selects the correct sign by comparing model FRFs against measurement (>500× NRMSE separation)

This bypasses the need for off-diagonal anti-resonance extraction (which fails at practical frequency resolution) and SVD-based mode shape methods (which give mode shapes ≠ D' eigenvectors).

## Repository Structure

```
├── README.md
├── matlab/
│   ├── For_Share_Stiffness_Mmatrix_ID_MCL_plant_3DOF_urdf.m          # Original
│   └── For_Share_Stiffness_Mmatrix_ID_MCL_plant_3DOF_urdf_fixed.m    # Fixed (Methods 1-4)
├── python/
│   └── frf_identify_3dof.py             # Python equivalent (pinocchio + numpy)
├── paper/
│   ├── src/main.tex                     # LaTeX source
│   ├── figures/                         # All paper figures (PDF)
│   └── FJR_offdiag_fitting_paper.pdf    # Compiled paper
├── docs/
│   ├── SKILL_algorithms.md              # Tom Oomen algorithm reference (12 algorithms)
│   ├── SKILL.md                         # Skill metadata
│   ├── paper_list.md                    # Oomen publication tracking
│   └── cmif_methodology.md             # CMIF methodology document
└── urdf/
    └── SEA_URDF/                        # 7-DOF SEA arm URDF + meshes
```

## Methods Implemented

| Method | Source | Signs? | Off-diag Error |
|--------|--------|--------|----------------|
| 1. Roots (symbolic) | Diagonal FRF polynomial roots | No (magnitude only) | 0.04–0.31% |
| 2. CMIF + Vieta | CMIF peak + parabolic interpolation | No (magnitude only) | 0.03–0.07% |
| 3. Matrix Pencil | SVD mode shape → D' reconstruction | Yes (but inaccurate) | 16–132% |
| 4. Iterative D' | Refinement of Method 3 | Yes (same as 3) | 16–132% |
| **5. Off-diag FRF fitting** | **Proposed** | **Yes (accurate)** | **0.03–0.07%** |

## Quick Start

### Python (recommended)

```bash
pip install pin numpy scipy matplotlib

# Set pinocchio path (if installed via pip install pin)
export PYTHONPATH="/path/to/cmeel.prefix/lib/python3.12/site-packages:$PYTHONPATH"

cd python/
python frf_identify_3dof.py
```

### MATLAB

```matlab
cd matlab/
For_Share_Stiffness_Mmatrix_ID_MCL_plant_3DOF_urdf_fixed
```

Requires: Robotics System Toolbox, Symbolic Math Toolbox.

### Paper Compilation

```bash
cd paper/src/
pdflatex main.tex && pdflatex main.tex
```

## Robot Parameters (SEA Arm)

| Joint | $J_m$ (×10⁻⁴) | $N$ | $k$ (Nm/rad) |
|-------|----------------|-----|---------------|
| J2 | 0.6793 | 100 | 2500 |
| J4 | 0.2639 | 50 | 1800 |
| J6 | 0.3057 | 50 | 800 |

Pose 1: $q = [0, \pi/2, -\pi/2]^T$ for joints J2, J4, J6.

## Algorithm Summary

### Step 1–4: Diagonal FRF Analysis (established)

**Step 1** — Motor inertia from high-frequency asymptote: $|G_{ii}(j\omega)| \to 1/(J_{m,i}\omega)$.

**Step 2** — Common resonance frequencies $\omega_{R,r}$ via CMIF peak detection with parabolic refinement:

$$\hat{\omega}_R = \omega_k + \frac{\sigma_{k-1} - \sigma_{k+1}}{2(\sigma_{k-1} - 2\sigma_k + \sigma_{k+1})} \Delta\omega$$

**Step 3** — Channel-specific anti-resonance $\omega_{AR,r,i}$ from magnitude dips of each diagonal FRF $G_{ii}$.

**Step 4** — Vieta's formulas on the polynomial difference $\delta_i(\lambda) = P_{N_i'}(\lambda) - P_{D'}(\lambda)$:

$$\boxed{k_i = J_{\text{eff},i} \left( \sum_r \omega_{R,r}^2 - \sum_r \omega_{AR,r,i}^2 \right)}$$

$$W_j = \frac{A_k + A_l - A_j}{2}, \quad |\Gamma_{jk}|^2 = W_j W_k - B_i \quad (\{j,k\} \neq i)$$

At this point $k_i$, $W_i$, $|\Gamma_{ij}|$, $(M_l^{-1})_{ii}$ are all determined. The remaining unknowns are the $\binom{n}{2}$ signs $\sigma_{ij} = \text{sign}(\Gamma_{ij})$.

### Step 5: Eigenvalue Constraint (8 → 4 candidates)

For each sign combination $\sigma = (\sigma_{12}, \sigma_{13}, \sigma_{23}) \in \{+1,-1\}^3$, construct:

$$D'_\sigma = \begin{pmatrix} W_1 & \sigma_{12}|\Gamma_{12}| & \sigma_{13}|\Gamma_{13}| \\ \sigma_{12}|\Gamma_{12}| & W_2 & \sigma_{23}|\Gamma_{23}| \\ \sigma_{13}|\Gamma_{13}| & \sigma_{23}|\Gamma_{23}| & W_3 \end{pmatrix}$$

Filter: keep only $\sigma$ satisfying $\text{eig}(D'_\sigma) = \{\omega_{R,r}^2\}$.

**Klein four-group reduction**: The similarity transform $D' \to S D' S$ with $S = \text{diag}(s_1, s_2, s_3)$, $s_i \in \{+1,-1\}$, preserves eigenvalues but maps $\Gamma_{ij} \to s_i s_j \Gamma_{ij}$. This partitions the 8 combinations into 2 orbits of size 4; exactly 4 pass the eigenvalue test.

**General $n$-DOF**: $2^{n(n-1)/2}$ total → $2^{n-1}$ candidates after filtering.

### Step 6: Off-Diagonal FRF Fitting (4 → 1 unique selection)

**Key insight**: All 4 surviving candidates produce identical diagonal FRFs $G_{ii}$, but yield **completely different off-diagonal FRFs** $G_{ij}$ because the cofactor $C_{ij}(s^2)$ depends on the signed off-diagonal elements of $D'$.

For each surviving $\sigma$:

**6a. Reconstruct** the signed inverse inertia matrix:

$$(M_l^{-1})_{ii} = \frac{W_i - \omega_{m,i}^2}{k_i}, \qquad (M_l^{-1})_{ij} = \frac{\sigma_{ij}|\Gamma_{ij}|}{\sqrt{k_i k_j}}$$

Then invert: $M_{l,\sigma} = [(M_l^{-1})_\sigma]^{-1}$.

**6b. Compute model FRF** via the $2n \times 2n$ impedance matrix:

$$Z_{CL}(s) = \begin{pmatrix} J_m s^2 + K_N & -K_N \\ -K_N & M_{l,\sigma} s^2 + K_N \end{pmatrix}, \quad G_\sigma(j\omega) = j\omega \cdot [Z_{CL}^{-1}(j\omega)]_{1:n,\,1:n}$$

**6c. Evaluate NRMSE** over all off-diagonal channels:

$$J_\sigma = \frac{\sqrt{\sum_k \sum_{i \neq j} |G_{\sigma,ij}(j\omega_k) - G_{ij}^{\text{meas}}(j\omega_k)|^2}}{\sqrt{\sum_k \sum_{i \neq j} |G_{ij}^{\text{meas}}(j\omega_k)|^2}} \times 100\%$$

**6d. Select**: $\sigma^* = \arg\min_\sigma J_\sigma$.

Result: correct candidate NRMSE = 0.18%, incorrect ≥ 99.4% (>500× separation).

## Citation

```
@unpublished{mcl2026fjr,
  title={Identification of Joint Stiffness and Signed Inertia Matrix 
         for Multi-DOF Flexible Joint Robots via Off-Diagonal FRF Fitting},
  author={Motion Control Lab, DGIST},
  year={2026}
}
```

## License

MIT License (code). See individual files for details.
