#!/usr/bin/env python3
"""
3-DOF FJR Stiffness & Inertia Matrix Identification
====================================================
MATLAB → Python 변환.  pinocchio로 URDF mass matrix 계산,
sympy로 전달함수 심볼릭 유도, numpy/scipy로 FRF & 식별.

방법 1: roots 기반 (심볼릭 다항식의 근)
방법 2: CMIF 기반 (SVD peak + parabolic interpolation)
방법 3: Matrix Pencil (D' reconstruct from CMIF mode shapes)
방법 4: Iterative D' Refinement

사용법:
    export PYTHONPATH="/usr/local/lib/python3.12/dist-packages/cmeel.prefix/lib/python3.12/site-packages:$PYTHONPATH"
    python3 frf_identify_3dof.py
"""

import os
import sys
import numpy as np
from numpy.linalg import inv, det, eig, svd, norm
import matplotlib.pyplot as plt
from scipy.signal import TransferFunction, freqs

# --------------- pinocchio import (cmeel prefix) ---------------
_cmeel_prefix = "/usr/local/lib/python3.12/dist-packages/cmeel.prefix/lib/python3.12/site-packages"
if _cmeel_prefix not in sys.path:
    sys.path.insert(0, _cmeel_prefix)
import pinocchio as pin

# ==============================================================
#  Robot Parameters  (SEA arm — getRobotParams_MCL('SEA'))
# ==============================================================
# 7-DOF SEA arm parameters (0-indexed in Python)
# MATLAB: params.Jm = [0.7, 0.67929, 0.66, 0.26386, 0.46, 0.30574, 0.411]*1e-4
# MATLAB: params.Nm = [100 100 100 50 50 50 50]
# MATLAB: params.K  = [4500, 2500, 5000, 1800, 2400, 800, 3000]
Jm_param = np.array([0.7, 0.67929, 0.66, 0.26386, 0.46, 0.30574, 0.411]) * 1e-4
Nm_param = np.array([100, 100, 100, 50, 50, 50, 50])
K_param  = np.array([4500, 2500, 5000, 1800, 2400, 800, 3000], dtype=float)

# MATLAB: selected_joints = [2,4,6] (1-indexed) → Python 0-indexed = [1,3,5]
SELECTED_JOINTS = [1, 3, 5]   # J2, J4, J6
DOF = len(SELECTED_JOINTS)

# ==============================================================
#  URDF load & mass matrix
# ==============================================================
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
URDF_PATH = os.path.join(SCRIPT_DIR, "..", "urdf", "SEA_URDF", "urdf", "SEA arm_woEE.urdf")

print(f"Loading URDF: {URDF_PATH}")
model = pin.buildModelFromUrdf(URDF_PATH)
data  = model.createData()
print(f"  nq={model.nq}, nv={model.nv}, joints={[model.names[i] for i in range(model.njoints)]}")

Jm_local = Jm_param[SELECTED_JOINTS]
Nm_local = Nm_param[SELECTED_JOINTS]
Jr_local = Jm_local * Nm_local**2   # reflected inertia
K_true   = K_param[SELECTED_JOINTS]

# ==============================================================
#  Transfer function construction (numerical, not symbolic)
# ==============================================================
def build_impedance_matrix_s_coeffs(M_l, Jm_vec, K_vec, N_vec):
    """
    Z_CL(s) = M_total * s^2 + K_total  →  returns (A2, A0) such that Z_CL = A2*s^2 + A0.
    M_total = [M_l, 0; 0, Jm],  K_total with gear ratio.
    """
    n = len(Jm_vec)
    Jm = np.diag(Jm_vec)
    K  = np.diag(K_vec)
    Ninv = np.diag(1.0 / N_vec)

    M_total = np.block([[M_l,             np.zeros((n,n))],
                        [np.zeros((n,n)), Jm             ]])

    K_total = np.block([[ K,              -K @ Ninv           ],
                        [-Ninv @ K,        Ninv @ Ninv @ K    ]])
    return M_total, K_total


def compute_frf_matrix(M_l, Jm_vec, K_vec, N_vec, omega_vec):
    """
    motor torque → motor velocity FRF matrix at each frequency.
    G_{ij}(jω) = jω * cofactor(Z, motor_i_row, motor_j_col) / det(Z)
    Returns: G  shape (DOF, DOF, N_freq) complex
    """
    n = len(Jm_vec)
    N_freq = len(omega_vec)
    G = np.zeros((n, n, N_freq), dtype=complex)

    M_total, K_total = build_impedance_matrix_s_coeffs(M_l, Jm_vec, K_vec, N_vec)

    for k, w in enumerate(omega_vec):
        s = 1j * w
        Z = M_total * s**2 + K_total   # 6x6

        det_Z = det(Z)
        if abs(det_Z) < 1e-300:
            continue

        for i_out in range(n):
            for i_in in range(n):
                row_idx = n + i_out   # motor side output
                col_idx = n + i_in    # motor side input
                # cofactor = (-1)^(r+c) * det(minor)
                minor = np.delete(np.delete(Z, row_idx, axis=0), col_idx, axis=1)
                sign_ij = (-1) ** (row_idx + col_idx)
                G_pos = sign_ij * det(minor) / det_Z
                G[i_out, i_in, k] = s * G_pos   # position → velocity: multiply by s

    return G


def roots_from_poly(coeffs):
    """Return roots of polynomial (numpy convention: highest degree first)."""
    r = np.roots(coeffs)
    return r


def extract_resonances_from_tf(G_diag_frf, omega_vec, n_dof):
    """
    From diagonal FRF (complex vector), extract resonance and anti-resonance
    by finding poles and zeros of the underlying rational TF.
    For simulation data, we have the exact TF, but here we use peak/dip picking
    from the FRF magnitude as a more general approach.
    """
    pass  # not needed for the direct numerical approach below

# ==============================================================
#  Peak detection utilities
# ==============================================================
def parabolic_peak(sig, idx, x_vec):
    """Parabolic interpolation for sub-grid peak refinement."""
    if idx <= 0 or idx >= len(sig) - 1:
        return x_vec[idx]
    sm1, s0, sp1 = sig[idx-1], sig[idx], sig[idx+1]
    denom = 2.0 * (sm1 - 2*s0 + sp1)
    if abs(denom) < 1e-15:
        return x_vec[idx]
    delta = (sm1 - sp1) / denom
    dx = x_vec[1] - x_vec[0]
    return x_vec[idx] + delta * dx


def find_peaks(sig, min_prominence=3.0, min_distance_idx=10):
    """Simple peak finder (prominence-based)."""
    peaks = []
    for i in range(1, len(sig)-1):
        if sig[i] > sig[i-1] and sig[i] > sig[i+1]:
            # check prominence
            left_min = np.min(sig[max(0,i-100):i])
            right_min = np.min(sig[i+1:min(len(sig),i+101)])
            prom = sig[i] - max(left_min, right_min)
            if prom >= min_prominence:
                peaks.append(i)
    # enforce minimum distance
    if len(peaks) > 1:
        filtered = [peaks[0]]
        for p in peaks[1:]:
            if p - filtered[-1] >= min_distance_idx:
                filtered.append(p)
        peaks = filtered
    return peaks


def recover_Minv(A_vals, B_vals, omega_m_sq, k_vec):
    """Recover M^{-1} diagonal and off-diagonal from polynomial coefficient comparison."""
    W_est = np.zeros(3)
    W_est[0] = (A_vals[1] + A_vals[2] - A_vals[0]) / 2
    W_est[1] = (A_vals[0] + A_vals[2] - A_vals[1]) / 2
    W_est[2] = (A_vals[0] + A_vals[1] - A_vals[2]) / 2
    diag_est = (W_est - omega_m_sq) / k_vec
    m23 = np.sqrt(max(0, (W_est[1]*W_est[2] - B_vals[0]) / (k_vec[1]*k_vec[2])))
    m13 = np.sqrt(max(0, (W_est[0]*W_est[2] - B_vals[1]) / (k_vec[0]*k_vec[2])))
    m12 = np.sqrt(max(0, (W_est[0]*W_est[1] - B_vals[2]) / (k_vec[0]*k_vec[1])))
    return diag_est, m12, m13, m23


# ==============================================================
#  Frequency grid
# ==============================================================
ROUND_DIGITS = 3
delta_f   = 10**(-ROUND_DIGITS)
f_vec     = np.arange(0.1, 100 + delta_f/2, delta_f)
omega_vec = 2 * np.pi * f_vec
N_freq    = len(omega_vec)
print(f"Frequency grid: {f_vec[0]:.4f} – {f_vec[-1]:.4f} Hz, Δf={delta_f} Hz, N={N_freq}")

# ==============================================================
#  Pose definitions
# ==============================================================
q_active_log = np.array([
    [0,      0,       0,       0      ],   # J2
    [np.pi/2, 0.0443, -0.4154, 1.5329],   # J4
    [-np.pi/2, 0.5236, 1.3909, -1.2658],  # J6
])
N_poses = q_active_log.shape[1]

# ==============================================================
#  Main computation: per-pose FRF generation
# ==============================================================
FRF_log      = []
M_log        = []
omega_R_log  = []
omega_AR_log = []

for p in range(N_poses):
    print(f"\n--- Pose {p+1}/{N_poses} ---")
    q_active = q_active_log[:, p]
    q_full = np.zeros(7)
    q_full[SELECTED_JOINTS] = q_active

    # Mass matrix via pinocchio
    pin.computeAllTerms(model, data, q_full, np.zeros(7))
    M_full = np.array(data.M)
    M_l_num = M_full[np.ix_(SELECTED_JOINTS, SELECTED_JOINTS)]
    M_log.append(M_l_num)
    print(f"  M_l =\n{M_l_num}")

    # FRF computation
    print("  Computing FRF matrix...", end="", flush=True)
    G = compute_frf_matrix(M_l_num, Jm_local, K_true, Nm_local, omega_vec)
    FRF_log.append(G)
    print(" done.")

    # Extract resonances / anti-resonances from diagonal FRF (roots of det)
    # We use the direct analytical approach: build Z_CL, find det polynomial roots
    M_total, K_total = build_impedance_matrix_s_coeffs(M_l_num, Jm_local, K_true, Nm_local)
    # det(Z_CL) is polynomial in s^2 → roots give resonance
    # For simplicity, use peak picking on FRF magnitude

    omega_AR_p = np.zeros((DOF, DOF))
    omega_R_p  = None

    for i_diag in range(DOF):
        mag = np.abs(G[i_diag, i_diag, :])
        mag_db = 20 * np.log10(mag + 1e-30)

        # Anti-resonances (dips)
        dip_idx = find_peaks(-mag_db, min_prominence=3, min_distance_idx=int(0.5/delta_f))
        ar_freqs = []
        for di in dip_idx:
            ar_freqs.append(parabolic_peak(-mag_db, di, omega_vec))
        ar_freqs = sorted(ar_freqs)[:DOF]
        for r in range(min(len(ar_freqs), DOF)):
            omega_AR_p[i_diag, r] = ar_freqs[r]

        # Resonances (peaks) — use first diagonal channel
        if i_diag == 0:
            pk_idx = find_peaks(mag_db, min_prominence=3, min_distance_idx=int(0.5/delta_f))
            res_freqs = []
            for pi_ in pk_idx:
                res_freqs.append(parabolic_peak(mag_db, pi_, omega_vec))
            omega_R_p = np.array(sorted(res_freqs)[:DOF])

        print(f"  G_{i_diag+1}{i_diag+1} AR: [{', '.join(f'{w/(2*np.pi):.4f}' for w in omega_AR_p[i_diag,:])}] Hz")

    print(f"  공진: [{', '.join(f'{w/(2*np.pi):.4f}' for w in omega_R_p)}] Hz")
    omega_R_log.append(omega_R_p)
    omega_AR_log.append(omega_AR_p)


# ==============================================================
#  Identification (pose_id = 0, first pose)
# ==============================================================
pose_id = 0
G_frf    = FRF_log[pose_id]
omega_AR = omega_AR_log[pose_id]
omega_R  = omega_R_log[pose_id]
M_l_id   = M_log[pose_id]
M_inv_gt = inv(M_l_id)
k_loc    = K_true.copy()
omega_m_sq = k_loc / Jr_local

print("\n" + "="*60)
print(f"  Identification at Pose {pose_id+1}")
print("="*60)
print(f"[Ground Truth M^{{-1}}]\n{M_inv_gt}")
print(f"[Ground Truth K] = {K_true}")

# ---------------------------------------------------------------
#  Method 1: roots-based (using peak-picked frequencies)
# ---------------------------------------------------------------
print("\n===== [방법 1] roots 기반 =====")

K_est_roots = np.zeros(DOF)
for i in range(DOF):
    K_est_roots[i] = Jr_local[i] * (np.sum(omega_R**2) - np.sum(omega_AR[i,:]**2))
    err = abs(K_est_roots[i] - K_true[i]) / K_true[i] * 100
    print(f"  Joint {SELECTED_JOINTS[i]+1}: True={K_true[i]:.2f}, Est={K_est_roots[i]:.2f}, Err={err:.2f}%")

# Polynomial coefficient comparison for M^{-1}
omega_ref = np.exp(np.mean(np.log(omega_R)))  # geometric mean
A_vals = np.zeros(DOF)
B_vals = np.zeros(DOF)
for i in range(DOF):
    lam_R_n  = (omega_R / omega_ref)**2
    lam_AR_n = (omega_AR[i,:] / omega_ref)**2
    # characteristic polynomial P(λ) = Π(λ - λ_k)
    P_R  = np.poly(lam_R_n)    # highest degree first
    P_AR = np.poly(lam_AR_n)
    delta_c = P_AR - P_R
    A_vals[i] = -delta_c[2] / delta_c[1] * omega_ref**2
    B_vals[i] =  delta_c[3] / delta_c[1] * omega_ref**4

Minv_diag_r, m12_r, m13_r, m23_r = recover_Minv(A_vals, B_vals, omega_m_sq, k_loc)
print("[대각항 (M^{-1})_ii]")
for i in range(DOF):
    err = abs(Minv_diag_r[i] - M_inv_gt[i,i]) / abs(M_inv_gt[i,i]) * 100
    print(f"  ({i+1},{i+1}): True={M_inv_gt[i,i]:+.6f}  Est={Minv_diag_r[i]:+.6f}  Err={err:.2f}%")
print("[비대각항 |(M^{-1})_ij|]")
for (ii,jj,est) in [(0,1,m12_r),(0,2,m13_r),(1,2,m23_r)]:
    tv = abs(M_inv_gt[ii,jj])
    err = abs(est - tv) / tv * 100
    print(f"  |({ii+1},{jj+1})|: True={tv:.6f}  Est={est:.6f}  Err={err:.2f}%")


# ---------------------------------------------------------------
#  Method 2: CMIF-based
# ---------------------------------------------------------------
print("\n===== [방법 2] CMIF 기반 =====")

sigma1 = np.zeros(N_freq)
for k in range(N_freq):
    sv = svd(G_frf[:,:,k], compute_uv=False)
    sigma1[k] = sv[0]

sigma1_db = 20 * np.log10(sigma1 + 1e-30)
pk_idx_cmif = find_peaks(sigma1_db, min_prominence=3, min_distance_idx=int(0.5/delta_f))
omega_R_cmif = np.array([parabolic_peak(sigma1_db, pi_, omega_vec) for pi_ in pk_idx_cmif[:DOF]])
omega_R_cmif = np.sort(omega_R_cmif)

omega_AR_cmif = np.zeros((DOF, DOF))
for i in range(DOF):
    mag_db = 20 * np.log10(np.abs(G_frf[i,i,:]) + 1e-30)
    dip_idx = find_peaks(-mag_db, min_prominence=3, min_distance_idx=int(0.5/delta_f))
    ar_f = sorted([parabolic_peak(-mag_db, di, omega_vec) for di in dip_idx])[:DOF]
    for r in range(len(ar_f)):
        omega_AR_cmif[i, r] = ar_f[r]

print(f"공진 (CMIF): [{', '.join(f'{w/(2*np.pi):.4f}' for w in omega_R_cmif)}] Hz")
for i in range(DOF):
    print(f"G_{i+1}{i+1} 반공진: [{', '.join(f'{w/(2*np.pi):.4f}' for w in omega_AR_cmif[i,:])}] Hz")

K_est_cmif = np.zeros(DOF)
for i in range(DOF):
    K_est_cmif[i] = Jr_local[i] * (np.sum(omega_R_cmif**2) - np.sum(omega_AR_cmif[i,:]**2))
    err = abs(K_est_cmif[i] - K_true[i]) / K_true[i] * 100
    print(f"  Joint {SELECTED_JOINTS[i]+1}: True={K_true[i]:.2f}, Est={K_est_cmif[i]:.2f}, Err={err:.2f}%")

# M^{-1} via polynomial coefficients
P_R_cmif = np.poly(omega_R_cmif**2)
A_cmif = np.zeros(DOF)
B_cmif = np.zeros(DOF)
for i in range(DOF):
    P_AR = np.poly(omega_AR_cmif[i,:]**2)
    dc = P_AR - P_R_cmif
    A_cmif[i] = -dc[2] / dc[1]
    B_cmif[i] =  dc[3] / dc[1]

Minv_diag_c, m12_c, m13_c, m23_c = recover_Minv(A_cmif, B_cmif, omega_m_sq, k_loc)
print("[대각항]")
for i in range(DOF):
    err = abs(Minv_diag_c[i] - M_inv_gt[i,i]) / abs(M_inv_gt[i,i]) * 100
    print(f"  ({i+1},{i+1}): True={M_inv_gt[i,i]:+.6f}  Est={Minv_diag_c[i]:+.6f}  Err={err:.2f}%")


# ---------------------------------------------------------------
#  Method 3: Matrix Pencil
# ---------------------------------------------------------------
print("\n===== [방법 3] Matrix Pencil 기반 =====")

V_modes = np.zeros((DOF, DOF))
for r in range(DOF):
    k_res = np.argmin(np.abs(omega_vec - omega_R_cmif[r]))

    G_mat = G_frf[:, :, k_res]
    scale = Jm_local * (1j * omega_vec[k_res])
    G_tilde = np.diag(scale) @ G_mat

    _, _, Vh = svd(G_tilde)
    v_complex = Vh[0, :].conj()   # numpy svd returns V^H, first row = first right singular vector

    # Fix 1: extract real mode shape
    max_idx = np.argmax(np.abs(v_complex))
    phase_corr = np.angle(v_complex[max_idx])
    v_real = np.real(v_complex * np.exp(-1j * phase_corr))
    V_modes[:, r] = v_real / norm(v_real)

print(f"[Mode shapes (columns)]\n{V_modes}")

Lambda_R = np.diag(omega_R_cmif**2)
D_prime = np.real(V_modes @ Lambda_R @ inv(V_modes))

asym_err = norm(D_prime - D_prime.T, 'fro') / norm(D_prime, 'fro') * 100
print(f"[D' reconstructed]  (비대칭도: {asym_err:.4f}%)")
print(D_prime)
D_prime = (D_prime + D_prime.T) / 2   # symmetrize

# Fix 2: eig(N_i') → anti-resonances → K
print("\n--- eig(N_i') 기반 반공진 ---")
e_mat = np.eye(DOF)
omega_AR_mp = np.zeros((DOF, DOF))
for i in range(DOF):
    N_i_prime = D_prime - omega_m_sq[i] * np.outer(e_mat[:,i], e_mat[:,i])
    eigvals = np.sort(np.real(eig(N_i_prime)[0]))
    if np.any(eigvals < -1e-6):
        print(f"  [경고] G_{i+1}{i+1}: 음수 고유값 {np.min(eigvals):.4e}")
    eigvals = np.maximum(eigvals, 0)
    omega_AR_mp[i, :] = np.sqrt(eigvals)
    print(f"  G_{i+1}{i+1} AR(MP): [{', '.join(f'{w/(2*np.pi):.4f}' for w in omega_AR_mp[i,:])}] Hz")

K_est_mp = np.zeros(DOF)
for i in range(DOF):
    K_est_mp[i] = Jr_local[i] * (np.sum(omega_R_cmif**2) - np.sum(omega_AR_mp[i,:]**2))
    err = abs(K_est_mp[i] - K_true[i]) / K_true[i] * 100
    print(f"  K_{SELECTED_JOINTS[i]+1}: True={K_true[i]:.2f}, Est={K_est_mp[i]:.2f}, Err={err:.2f}%")

# Fix 3: M^{-1} using estimated K
W_mp = np.diag(D_prime)
Gamma_mp = D_prime - np.diag(W_mp)
Minv_diag_mp = (W_mp - omega_m_sq) / K_est_mp
Minv_mp = np.diag(Minv_diag_mp)
for i in range(DOF):
    for j in range(i+1, DOF):
        val = Gamma_mp[i,j] / np.sqrt(K_est_mp[i] * K_est_mp[j])
        Minv_mp[i,j] = val
        Minv_mp[j,i] = val

print(f"\n[M^{{-1}} Matrix Pencil (부호 포함)]\n{Minv_mp}")
print(f"[Ground Truth]\n{M_inv_gt}")

for i in range(DOF):
    err = abs(Minv_mp[i,i] - M_inv_gt[i,i]) / abs(M_inv_gt[i,i]) * 100
    print(f"  ({i+1},{i+1}): True={M_inv_gt[i,i]:+.6f}  Est={Minv_mp[i,i]:+.6f}  Err={err:.2f}%")
for (ii,jj) in [(0,1),(0,2),(1,2)]:
    err = abs(Minv_mp[ii,jj] - M_inv_gt[ii,jj]) / abs(M_inv_gt[ii,jj]) * 100
    print(f"  ({ii+1},{jj+1}): True={M_inv_gt[ii,jj]:+.6f}  Est={Minv_mp[ii,jj]:+.6f}  Err={err:.2f}%")


# ---------------------------------------------------------------
#  Method 4: Iterative D' Refinement
# ---------------------------------------------------------------
print("\n===== [방법 4] Iterative D' Refinement =====")

ITER_MAX = 50
ITER_TOL = 1e-12

D_iter = D_prime.copy()
K_hist    = np.zeros((ITER_MAX, DOF))
Minv_hist = [None] * ITER_MAX
conv_hist = np.zeros(ITER_MAX)
n_iter_done = 0

for it in range(ITER_MAX):
    D_old = D_iter.copy()

    # Step 1: eig(N_i') → anti-resonances
    omega_AR_iter = np.zeros((DOF, DOF))
    for i in range(DOF):
        N_i = D_iter - omega_m_sq[i] * np.outer(e_mat[:,i], e_mat[:,i])
        ev = np.sort(np.real(eig(N_i)[0]))
        ev = np.maximum(ev, 0)
        omega_AR_iter[i, :] = np.sqrt(ev)

    # Step 2: Vieta → K
    K_iter = np.zeros(DOF)
    for i in range(DOF):
        K_iter[i] = Jr_local[i] * (np.sum(omega_R_cmif**2) - np.sum(omega_AR_iter[i,:]**2))
    K_iter = np.maximum(K_iter, 1.0)   # clamp to positive

    # Step 3: D' elements → M^{-1}
    W_iter = np.diag(D_iter)
    Gamma_iter = D_iter - np.diag(W_iter)
    Minv_diag_iter = (W_iter - omega_m_sq) / K_iter
    Minv_iter = np.diag(Minv_diag_iter)
    for i in range(DOF):
        for j in range(i+1, DOF):
            val = Gamma_iter[i,j] / np.sqrt(K_iter[i] * K_iter[j])
            Minv_iter[i,j] = val
            Minv_iter[j,i] = val

    # Step 4: K, M^{-1} → D' reconstruction
    D_new = np.zeros((DOF, DOF))
    for i in range(DOF):
        D_new[i,i] = omega_m_sq[i] + K_iter[i] * Minv_iter[i,i]
        for j in range(i+1, DOF):
            D_new[i,j] = np.sqrt(K_iter[i] * K_iter[j]) * Minv_iter[i,j]
            D_new[j,i] = D_new[i,j]
    D_iter = (D_new + D_new.T) / 2

    # Step 5: convergence
    rel_change = norm(D_iter - D_old, 'fro') / (norm(D_old, 'fro') + 1e-30)
    K_hist[it, :] = K_iter
    Minv_hist[it] = Minv_iter.copy()
    conv_hist[it] = rel_change
    n_iter_done = it + 1

    if rel_change < ITER_TOL:
        print(f"  수렴: iter {it+1}, ‖ΔD'‖/‖D'‖ = {rel_change:.2e}")
        break

if rel_change >= ITER_TOL:
    print(f"  [경고] {ITER_MAX}회 후 미수렴. 최종 = {rel_change:.2e}")

K_est_iter    = K_hist[n_iter_done-1, :]
Minv_iter_fin = Minv_hist[n_iter_done-1]

# Final anti-resonances
omega_AR_iter_fin = np.zeros((DOF, DOF))
for i in range(DOF):
    N_i = D_iter - omega_m_sq[i] * np.outer(e_mat[:,i], e_mat[:,i])
    ev = np.sort(np.real(eig(N_i)[0]))
    ev = np.maximum(ev, 0)
    omega_AR_iter_fin[i, :] = np.sqrt(ev)

print(f"\n--- Stiffness (Iterative, {n_iter_done} iter) ---")
for i in range(DOF):
    err = abs(K_est_iter[i] - K_true[i]) / K_true[i] * 100
    print(f"  Joint {SELECTED_JOINTS[i]+1}: True={K_true[i]:.2f}, Est={K_est_iter[i]:.2f}, Err={err:.2f}%")

print(f"\n[M^{{-1}} Iterative]\n{Minv_iter_fin}")
print(f"[Ground Truth]\n{M_inv_gt}")

for i in range(DOF):
    err = abs(Minv_iter_fin[i,i] - M_inv_gt[i,i]) / abs(M_inv_gt[i,i]) * 100
    print(f"  ({i+1},{i+1}): True={M_inv_gt[i,i]:+.6f}  Est={Minv_iter_fin[i,i]:+.6f}  Err={err:.2f}%")
for (ii,jj) in [(0,1),(0,2),(1,2)]:
    err = abs(Minv_iter_fin[ii,jj] - M_inv_gt[ii,jj]) / abs(M_inv_gt[ii,jj]) * 100
    print(f"  ({ii+1},{jj+1}): True={M_inv_gt[ii,jj]:+.6f}  Est={Minv_iter_fin[ii,jj]:+.6f}  Err={err:.2f}%")


# ==============================================================
#  Final comparison table
# ==============================================================
print("\n" + "="*70)
print("  방법 1 vs 2 vs 3 vs 4 비교")
print("="*70)
header = f"{'':12s} {'roots':>10s} {'CMIF':>10s} {'MP':>10s} {'Iter':>10s} {'True':>10s}"
print(header)
for i in range(DOF):
    print(f"K_{SELECTED_JOINTS[i]+1:d}         {K_est_roots[i]:10.3f} {K_est_cmif[i]:10.3f} "
          f"{K_est_mp[i]:10.3f} {K_est_iter[i]:10.3f} {K_true[i]:10.3f}")

print("\n(M^{-1})_{ii}:")
for i in range(DOF):
    print(f"  ({i+1},{i+1})  roots={Minv_diag_r[i]:+.6f}  CMIF={Minv_diag_c[i]:+.6f}  "
          f"MP={Minv_mp[i,i]:+.6f}  Iter={Minv_iter_fin[i,i]:+.6f}  True={M_inv_gt[i,i]:+.6f}")

print("(M^{-1})_{ij} (부호 포함, MP & Iter):")
for (ii,jj) in [(0,1),(0,2),(1,2)]:
    print(f"  ({ii+1},{jj+1})  MP={Minv_mp[ii,jj]:+.6f}  Iter={Minv_iter_fin[ii,jj]:+.6f}  True={M_inv_gt[ii,jj]:+.6f}")


# ==============================================================
#  Plots
# ==============================================================
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# (1) CMIF
ax = axes[0, 0]
ax.semilogx(f_vec, sigma1_db, 'b-', linewidth=1.5)
for w in omega_R_cmif:
    ax.axvline(w/(2*np.pi), color='r', linestyle='--', alpha=0.7)
ax.set_xlabel('Frequency (Hz)')
ax.set_ylabel('σ₁ (dB)')
ax.set_title(f'CMIF — Pose {pose_id+1}')
ax.set_xlim([0.1, 100])
ax.grid(True)

# (2) Diagonal Bode
ax = axes[0, 1]
colors = ['b', 'r', 'g']
for i in range(DOF):
    mag_db = 20*np.log10(np.abs(G_frf[i,i,:]) + 1e-30)
    ax.semilogx(f_vec, mag_db, color=colors[i], linewidth=1, label=f'G_{i+1}{i+1}')
ax.set_xlabel('Frequency (Hz)')
ax.set_ylabel('Magnitude (dB)')
ax.set_title('Diagonal FRFs')
ax.legend()
ax.set_xlim([0.1, 100])
ax.grid(True)

# (3) Convergence
ax = axes[1, 0]
ax.semilogy(np.arange(1, n_iter_done+1), conv_hist[:n_iter_done], 'b-o', markersize=3)
ax.axhline(ITER_TOL, color='r', linestyle='--', label=f'tol={ITER_TOL:.0e}')
ax.set_xlabel('Iteration')
ax.set_ylabel('‖ΔD\'‖ / ‖D\'‖')
ax.set_title('Iterative D\' Convergence')
ax.legend()
ax.grid(True)

# (4) K error history
ax = axes[1, 1]
for i in range(DOF):
    K_err = np.abs(K_hist[:n_iter_done, i] - K_true[i]) / K_true[i] * 100
    ax.plot(np.arange(1, n_iter_done+1), K_err, '-o', markersize=3,
            label=f'K_{SELECTED_JOINTS[i]+1}')
ax.set_xlabel('Iteration')
ax.set_ylabel('K Error (%)')
ax.set_title('Stiffness Estimation Error')
ax.legend()
ax.grid(True)

plt.tight_layout()
plt.savefig(os.path.join(SCRIPT_DIR, 'identification_results.png'), dpi=150)
print(f"\n[Plot saved: identification_results.png]")
plt.show()
