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

$$k_i = J_{\text{eff},i} \left( \sum_r \omega_{R,r}^2 - \sum_r \omega_{AR,r,i}^2 \right)$$

### Step 5: Eigenvalue Constraint

For 3-DOF: 8 sign combinations → 4 candidates (Klein four-group $\mathbb{Z}_2 \times \mathbb{Z}_2$).

### Step 6: Off-Diagonal FRF Fitting (proposed)

$$\sigma^* = \arg\min_\sigma \sum_k \sum_{i \neq j} |G_{\sigma,ij}(j\omega_k) - G_{ij}^{\text{meas}}(j\omega_k)|^2$$

Result: correct candidate NRMSE = 0.18%, incorrect ≥ 99.4%.

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
