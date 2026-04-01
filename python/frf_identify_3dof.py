#!/usr/bin/env python3
"""
3-DOF FJR Stiffness & Inertia Matrix Identification (v2)
=========================================================
5 methods + multi-pose validation.

Methods:
  1-2. CMIF + Vieta (magnitude only)
  3. Matrix Pencil (SVD mode shape, signed but inaccurate)
  4. Iterative D' Refinement
  5. Off-diagonal FRF fitting (proposed, signed + accurate)

Usage:
    pip install pin numpy scipy matplotlib
    python frf_identify_3dof.py
"""
import os, sys, time
import numpy as np
from numpy.linalg import inv, det, eig, svd, norm
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

try:
    import pinocchio as pin
    if not hasattr(pin, 'buildModelFromUrdf'):
        raise ImportError
except ImportError:
    _cmeel = "/usr/local/lib/python3.12/dist-packages/cmeel.prefix/lib/python3.12/site-packages"
    if os.path.isdir(_cmeel):
        sys.path.insert(0, _cmeel)
    import pinocchio as pin

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# ================================================================
#  Robot Parameters (SEA arm)
# ================================================================
Jm_param = np.array([0.7, 0.67929, 0.66, 0.26386, 0.46, 0.30574, 0.411]) * 1e-4
Nm_param = np.array([100, 100, 100, 50, 50, 50, 50], dtype=float)
K_param  = np.array([4500, 2500, 5000, 1800, 2400, 800, 3000], dtype=float)
SEL = [1, 3, 5]; DOF = len(SEL)
Jm = Jm_param[SEL]; Nm = Nm_param[SEL]; K_true = K_param[SEL]; Jr = Jm*Nm**2

URDF_PATH = os.path.join(SCRIPT_DIR, '..', 'urdf', 'SEA_URDF', 'urdf', 'SEA arm_woEE.urdf')
model = pin.buildModelFromUrdf(URDF_PATH)
data  = model.createData()

Q_POSES = np.array([[0, 0, 0, 0],
                     [np.pi/2, 0.0443, -0.4154, 1.5329],
                     [-np.pi/2, 0.5236, 1.3909, -1.2658]])
N_POSES = Q_POSES.shape[1]

# ================================================================
def compute_Ml(q3):
    q = np.zeros(7); q[SEL] = q3
    pin.computeAllTerms(model, data, q, np.zeros(7))
    return np.array(data.M, dtype=float)[np.ix_(SEL, SEL)]

def build_G(Ml, K_v, omega):
    s = 1j*omega; n = DOF
    Mt = np.block([[Ml, np.zeros((n,n))],[np.zeros((n,n)), np.diag(Jm)]])
    Ni = np.diag(1.0/Nm); Kd = np.diag(K_v)
    Kt = np.block([[Kd,-Kd@Ni],[-Ni@Kd,Ni@Ni@Kd]])
    Z = Mt*s**2+Kt; dZ = complex(det(Z))
    G = np.zeros((n,n),dtype=complex)
    for io in range(n):
        for ii in range(n):
            minor = np.delete(np.delete(Z,n+io,0),n+ii,1)
            G[io,ii] = s*((-1)**(io+ii))*complex(det(minor))/dZ
    return G

def parabolic_peak(sig, idx, x):
    if idx<=0 or idx>=len(sig)-1: return x[idx]
    a,b,c = sig[idx-1], sig[idx], sig[idx+1]
    d = 2.0*(a-2*b+c)
    if abs(d)<1e-15: return x[idx]
    return x[idx] + (a-c)/d*(x[1]-x[0])

def find_peaks(sig, prom=3.0, dist=10):
    pk = []
    for i in range(1,len(sig)-1):
        if sig[i]>sig[i-1] and sig[i]>sig[i+1]:
            lm = np.min(sig[max(0,i-200):i])
            rm = np.min(sig[i+1:min(len(sig),i+201)])
            if sig[i]-max(lm,rm) >= prom: pk.append(i)
    if len(pk)>1:
        f=[pk[0]]
        for p in pk[1:]:
            if p-f[-1]>=dist: f.append(p)
        pk=f
    return pk

# ================================================================
#  Core identification for a single pose
# ================================================================
def identify(M_l):
    Minv_gt = inv(M_l)
    omega_m_sq = K_true/Jr

    # FRF (high-res)
    df = 0.001; f_vec = np.arange(0.5, 50+df/2, df); Nf = len(f_vec)
    omega_vec = 2*np.pi*f_vec
    G_frf = np.zeros((DOF,DOF,Nf),dtype=complex)
    for k in range(Nf):
        G_frf[:,:,k] = build_G(M_l, K_true, omega_vec[k])

    # CMIF → resonance
    sig1 = np.zeros(Nf)
    for k in range(Nf):
        sig1[k] = svd(G_frf[:,:,k],compute_uv=False)[0]
    sig1_db = 20*np.log10(sig1+1e-30)
    pk = find_peaks(sig1_db, prom=3, dist=int(0.5/df))
    omega_R = np.sort([parabolic_peak(sig1_db,p,omega_vec) for p in pk[:DOF]])

    # Diagonal anti-resonance
    omega_AR = np.zeros((DOF,DOF))
    for i in range(DOF):
        mdb = 20*np.log10(np.abs(G_frf[i,i,:])+1e-30)
        dips = find_peaks(-mdb, prom=3, dist=int(0.5/df))
        ar = sorted([parabolic_peak(-mdb,d,omega_vec) for d in dips])[:DOF]
        for r in range(len(ar)): omega_AR[i,r] = ar[r]

    # ---- Method 2: Vieta ----
    P_R = np.poly(omega_R**2)
    A = np.zeros(3); B = np.zeros(3)
    for i in range(3):
        dc = np.poly(omega_AR[i,:]**2) - P_R
        A[i]=-dc[2]/dc[1]; B[i]=dc[3]/dc[1]
    W = np.array([(A[1]+A[2]-A[0])/2,(A[0]+A[2]-A[1])/2,(A[0]+A[1]-A[2])/2])
    Gabs = np.array([np.sqrt(max(0,W[0]*W[1]-B[2])),
                     np.sqrt(max(0,W[0]*W[2]-B[1])),
                     np.sqrt(max(0,W[1]*W[2]-B[0]))])
    K_est = np.array([float(Jr[i])*(float(np.sum(omega_R**2))-float(np.sum(omega_AR[i,:]**2))) for i in range(3)])
    Minv_diag = (W-omega_m_sq)/K_est

    # ---- Method 3: Matrix Pencil ----
    V_m = np.zeros((DOF,DOF))
    for r in range(DOF):
        kr = np.argmin(np.abs(omega_vec-omega_R[r]))
        Gt = np.diag(Jm*1j*omega_vec[kr])@G_frf[:,:,kr]
        _,_,Vh = svd(Gt); vc = Vh[0,:].conj()
        mi = np.argmax(np.abs(vc))
        vr = np.real(vc*np.exp(-1j*np.angle(vc[mi])))
        V_m[:,r] = vr/norm(vr)
    D_mp = np.real(V_m@np.diag(omega_R**2)@inv(V_m))
    D_mp = (D_mp+D_mp.T)/2
    e_mat = np.eye(DOF)
    oAR_mp = np.zeros((DOF,DOF))
    for i in range(DOF):
        Ni = D_mp - omega_m_sq[i]*np.outer(e_mat[:,i],e_mat[:,i])
        ev = np.sort(np.real(eig(Ni)[0]))
        oAR_mp[i,:] = np.sqrt(np.maximum(ev,0))
    K_mp = np.array([float(Jr[i])*(float(np.sum(omega_R**2))-float(np.sum(oAR_mp[i,:]**2))) for i in range(3)])
    K_mp = np.maximum(K_mp,1.0)
    Wmp = np.diag(D_mp); Gmp = D_mp-np.diag(Wmp)
    Minv_mp = np.diag((Wmp-omega_m_sq)/K_mp)
    for i in range(DOF):
        for j in range(i+1,DOF):
            v = Gmp[i,j]/np.sqrt(K_mp[i]*K_mp[j])
            Minv_mp[i,j]=v; Minv_mp[j,i]=v

    # ---- Method 5: Off-diagonal FRF fitting ----
    target_eig = np.sort(omega_R**2)
    cands = []
    for s12 in [+1,-1]:
        for s13 in [+1,-1]:
            for s23 in [+1,-1]:
                Dc = np.diag(W).copy()
                Dc[0,1]=s12*Gabs[0]; Dc[1,0]=Dc[0,1]
                Dc[0,2]=s13*Gabs[1]; Dc[2,0]=Dc[0,2]
                Dc[1,2]=s23*Gabs[2]; Dc[2,1]=Dc[1,2]
                te = np.sort(np.real(eig(Dc)[0]))
                if norm(te-target_eig)/norm(target_eig)<1e-3:
                    cands.append((s12,s13,s23,Dc))

    df_fit = 0.5
    f_fit = np.arange(2.0, 45.0, df_fit)
    Nfit = len(f_fit)
    G_fit = np.zeros((DOF,DOF,Nfit),dtype=complex)
    for k in range(Nfit):
        G_fit[:,:,k] = build_G(M_l, K_true, 2*np.pi*f_fit[k])

    nrmse = np.zeros(len(cands))
    Minv_cands = []
    for ci,(s12,s13,s23,Dc) in enumerate(cands):
        Mc = np.diag([(W[i]-omega_m_sq[i])/K_est[i] for i in range(3)])
        for i in range(3):
            for j in range(i+1,3):
                v = float(Dc[i,j])/np.sqrt(float(K_est[i])*float(K_est[j]))
                Mc[i,j]=v; Mc[j,i]=v
        Ml_c = inv(Mc)
        Minv_cands.append(Mc.copy())
        en=0.0; ed=0.0
        for k in range(Nfit):
            Gm = build_G(Ml_c, K_est, 2*np.pi*f_fit[k])
            for ii in range(3):
                for jj in range(3):
                    if ii!=jj:
                        d = Gm[ii,jj]-G_fit[ii,jj,k]
                        en+=float((d*d.conjugate()).real)
                        ed+=float((G_fit[ii,jj,k]*G_fit[ii,jj,k].conjugate()).real)
        nrmse[ci] = (en/ed)**0.5*100.0

    bi = np.argmin(nrmse)
    best_signs = (cands[bi][0],cands[bi][1],cands[bi][2])
    Minv_m5 = Minv_cands[bi]

    return dict(
        M_l=M_l, M_inv_gt=Minv_gt, omega_R=omega_R, omega_AR=omega_AR,
        K_est=K_est, K_mp=K_mp, Minv_diag=Minv_diag, Gabs=Gabs, W=W,
        Minv_mp=Minv_mp, Minv_m5=Minv_m5, signs_m5=best_signs,
        cands=cands, nrmse=nrmse,
    )

# ================================================================
def print_results(r, pid):
    gt = r['M_inv_gt']
    gt_signs = (int(np.sign(gt[0,1])),int(np.sign(gt[0,2])),int(np.sign(gt[1,2])))
    print(f"\n  Resonances (Hz): {(r['omega_R']/(2*np.pi)).round(4)}")
    print(f"  K: {r['K_est'].round(2)} (true: {K_true})")
    print(f"  Signs: M5={r['signs_m5']} GT={gt_signs} {'✓' if r['signs_m5']==gt_signs else '✗'}")
    print(f"  NRMSE: ", end="")
    for ci,c in enumerate(r['cands']):
        tag = "←" if ci==np.argmin(r['nrmse']) else " "
        print(f"({c[0]:+d},{c[1]:+d},{c[2]:+d}):{r['nrmse'][ci]:7.2f}%{tag} ", end="")
    print()
    print(f"  M⁻¹ off-diag errors:")
    for (i,j) in [(0,1),(0,2),(1,2)]:
        e5 = abs(r['Minv_m5'][i,j]-gt[i,j])/abs(gt[i,j])*100
        emp = abs(r['Minv_mp'][i,j]-gt[i,j])/abs(gt[i,j])*100
        print(f"    ({i+1},{j+1}): M5={e5:.4f}%  MP={emp:.1f}%")

# ================================================================
def plot_all(results):
    N = len(results)
    fig, axes = plt.subplots(2, 2, figsize=(12, 9))

    # (1) NRMSE
    ax = axes[0,0]; cols = ['#e41a1c','#377eb8','#4daf4a','#984ea3']
    x = np.arange(N); w = 0.18
    for ci in range(4):
        vals = [r['nrmse'][ci] for r in results]
        lab = '(%+d,%+d,%+d)' % results[0]['cands'][ci][:3]
        ax.bar(x+ci*w-1.5*w, vals, w, color=cols[ci], label=lab, edgecolor='k', lw=0.3)
    ax.set_xticks(x); ax.set_xticklabels([f'Pose {i+1}' for i in range(N)])
    ax.set_yscale('log'); ax.set_ylabel('NRMSE (%)'); ax.set_title('Method 5: Sign Candidate NRMSE')
    ax.legend(fontsize=6, ncol=2); ax.grid(True,alpha=0.3,axis='y')

    # (2) Off-diag M^{-1} error
    ax = axes[0,1]; pairs = [(0,1),(0,2),(1,2)]
    e5=[]; emp=[]; xl=[]
    for r in results:
        for (i,j) in pairs:
            gt = r['M_inv_gt'][i,j]
            e5.append(abs(r['Minv_m5'][i,j]-gt)/abs(gt)*100)
            emp.append(abs(r['Minv_mp'][i,j]-gt)/abs(gt)*100)
            xl.append(f'P{r["pose_id"]+1}({i+1},{j+1})' if 'pose_id' in r else f'({i+1},{j+1})')
    x2 = np.arange(len(e5))
    ax.bar(x2-0.2, emp, 0.35, label='Matrix Pencil', color='#ff7f0e', edgecolor='k', lw=0.2)
    ax.bar(x2+0.2, e5, 0.35, label='Method 5', color='#2ca02c', edgecolor='k', lw=0.2)
    ax.set_xticks(x2); ax.set_xticklabels(xl, fontsize=5, rotation=45)
    ax.set_yscale('log'); ax.set_ylabel('Error (%)'); ax.set_ylim([1e-3,500])
    ax.set_title('Off-Diagonal $M^{-1}_{ij}$ Error'); ax.legend(fontsize=7); ax.grid(True,alpha=0.3,axis='y')

    # (3) K error
    ax = axes[1,0]; kc = ['#1f77b4','#ff7f0e','#2ca02c']
    for pi,r in enumerate(results):
        for ji in range(DOF):
            err = abs(r['K_est'][ji]-K_true[ji])/K_true[ji]*100
            ax.bar(pi*DOF+ji, err, color=kc[ji], edgecolor='k', lw=0.2)
    ax.set_xticks([pi*DOF+1 for pi in range(N)])
    ax.set_xticklabels([f'Pose {i+1}' for i in range(N)])
    ax.set_ylabel('K Error (%)'); ax.set_title('Stiffness ID Error'); ax.grid(True,alpha=0.3,axis='y')

    # (4) M_l off-diag signs
    ax = axes[1,1]
    for pi,r in enumerate(results):
        Ml = r['M_l']
        ax.bar(pi*3+0, Ml[0,1], 0.8, color='#1f77b4', edgecolor='k', lw=0.2)
        ax.bar(pi*3+1, Ml[0,2], 0.8, color='#ff7f0e', edgecolor='k', lw=0.2)
        ax.bar(pi*3+2, Ml[1,2], 0.8, color='#2ca02c', edgecolor='k', lw=0.2)
    ax.axhline(0, color='k', lw=0.5)
    ax.set_xticks([pi*3+1 for pi in range(N)])
    ax.set_xticklabels([f'Pose {i+1}' for i in range(N)])
    ax.set_ylabel('$M_{l,ij}$'); ax.set_title('Link Mass Matrix Off-Diagonal')
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color='#1f77b4',label='$M_{12}$'),
                       Patch(color='#ff7f0e',label='$M_{13}$'),
                       Patch(color='#2ca02c',label='$M_{23}$')],fontsize=7)
    ax.grid(True,alpha=0.3,axis='y')

    plt.tight_layout()
    out = os.path.join(SCRIPT_DIR,'multipose_results.png')
    plt.savefig(out, dpi=150); plt.close()
    print(f"\n[Plot saved: {out}]")

# ================================================================
if __name__ == '__main__':
    print("="*60)
    print("  3-DOF FJR Identification — Methods 1-5 + Multi-Pose")
    print("="*60)
    results = []
    for p in range(N_POSES):
        qa = Q_POSES[:,p]
        Ml = compute_Ml(qa)
        print(f"\n{'='*60}\n  Pose {p+1}/{N_POSES}: q = {np.rad2deg(qa).round(1)}°\n{'='*60}")
        t0 = time.time()
        r = identify(Ml)
        r['pose_id'] = p; r['q'] = qa.copy()
        print(f"  [{time.time()-t0:.1f}s]")
        print_results(r, p)
        results.append(r)

    # Summary
    print("\n\n" + "="*70)
    print("  MULTI-POSE SUMMARY")
    print("="*70)
    hdr = f"{'Pose':>5s} {'Signs_M5':>12s} {'OK?':>4s} {'Best%':>9s} {'2nd%':>9s} {'Ratio':>7s} {'(1,2)%':>8s} {'(1,3)%':>8s} {'(2,3)%':>8s}"
    print(hdr)
    for r in results:
        gt_s = (int(np.sign(r['M_inv_gt'][0,1])),int(np.sign(r['M_inv_gt'][0,2])),int(np.sign(r['M_inv_gt'][1,2])))
        sn = np.sort(r['nrmse'])
        ratio = sn[1]/sn[0] if sn[0]>0 else float('inf')
        ok = '✓' if r['signs_m5']==gt_s else '✗'
        errs = [abs(r['Minv_m5'][i,j]-r['M_inv_gt'][i,j])/abs(r['M_inv_gt'][i,j])*100 for (i,j) in [(0,1),(0,2),(1,2)]]
        ss = str(r['signs_m5'])
        print(f"  {r['pose_id']+1:>3d} {ss:>12s} {ok:>4s} {sn[0]:>8.3f}% {sn[1]:>8.1f}% {ratio:>6.0f}x {errs[0]:>7.4f} {errs[1]:>7.4f} {errs[2]:>7.4f}")

    plot_all(results)
