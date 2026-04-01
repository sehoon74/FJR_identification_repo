# Multi-DOF FJR Inertia Matrix Identification: Off-Diagonal FRF Fitting Method

## 개요

다자유도 유연관절 로봇(FJR)의 관절 강성 $k_i$와 **부호 포함** 역관성행렬 $M_l^{-1}$을 
주파수 응답 함수(FRF) 측정으로부터 식별하는 방법론.

**핵심 기여**: 대각 FRF의 Vieta 분석으로 크기($|Γ_{ij}|$)를 구하고, 
비대각 FRF 전체 형상 피팅으로 부호를 유일하게 결정 — 
비대각 반공진 직접 측정 불필요.

---

## 문제 정의

### 시스템 모델

$n$-DOF FJR의 motor torque → motor velocity 전달함수:

$$G_{ij}(s) = \frac{s \cdot C_{ij}(s^2)}{\prod_l J_{m,l} s \cdot \det(D(s))}$$

분모 행렬 $D(s) \in \mathbb{R}^{n \times n}$:

$$D(s) = s^2 I + D', \quad D' = \text{diag}(W) + \Gamma$$

파라미터 정의:
- $W_i = \omega_{m,i}^2 + k_i (M_l^{-1})_{ii}$ — D' 대각 원소
- $\Gamma_{ij} = \sqrt{k_i k_j} \cdot (M_l^{-1})_{ij}$ — D' 비대각 원소 (대칭)
- $\omega_{m,i}^2 = k_i / (J_{m,i} N_i^2)$ — motor-side 고유진동수²

### 식별 목표

FRF 측정 $G_{ij}(j\omega)$ ($i,j = 1,...,n$)로부터:
1. 관절 강성 $k_i$
2. **부호 포함** 역관성행렬 $M_l^{-1}$ (대칭, $n(n+1)/2$ 원소)

---

## 알고리즘 (6 Steps)

### Step 1: Motor Inertia

고주파 점근선에서 $|G_{ii}(j\omega)| \to 1/(J_{m,i}\omega)$.

### Step 2: Common Resonance (CMIF)

매 주파수에서 SVD: $G(j\omega_k) = U \Sigma V^H$

CMIF = $\sigma_1(\omega)$ (최대 특이값). 피크 = 공진주파수 $\omega_{R,r}$.

포물선 보간으로 격자 해상도 이상의 정밀도:

$$\hat{\omega}_R = \omega_k + \frac{\sigma_{k-1} - \sigma_{k+1}}{2(\sigma_{k-1} - 2\sigma_k + \sigma_{k+1})} \Delta\omega$$

### Step 3: Channel-Specific Anti-Resonance

각 대각 FRF $G_{ii}$의 크기 dip에서 반공진 $\omega_{AR,r,i}$ 추출 (동일 보간 적용).

### Step 4: Vieta's Formulas

**강성 식별** (핵심 수식):

$$\boxed{k_i = J_{\text{eff},i} \left(\sum_r \omega_{R,r}^2 - \sum_r \omega_{AR,r,i}^2\right)}$$

여기서 $J_{\text{eff},i} = J_{m,i} N_i^2$.

**원리**: $P_{D'}(\lambda) = \det(D' - \lambda I)$와 $P_{N_i'}(\lambda) = \det(N_i' - \lambda I)$의 
다항식 차 $\delta_i(\lambda) = P_{N_i'} - P_{D'}$의 계수비를 이용.

$\delta_i(\lambda) = \alpha_i \lambda^{n-1} + \beta_i \lambda^{n-2} + \gamma_i \lambda^{n-3} + \cdots$

$A_i = -\beta_i/\alpha_i$로부터 **개별 $W_j$** 복원:

$$W_j = \frac{A_k + A_l - A_j}{2}, \quad \{j,k,l\} = \{1,2,3\}$$

$B_i = \gamma_i/\alpha_i$로부터 **비대각 크기** 복원:

$$|\Gamma_{jk}|^2 = W_j W_k - B_i, \quad \{j,k\} \neq i$$

이 시점에서 결정된 것: $k_i$, $W_i$, $|\Gamma_{ij}|$, $(M_l^{-1})_{ii}$.
**미결정**: $\text{sign}(\Gamma_{ij})$ — $\binom{n}{2}$개의 부호.

### Step 5: Eigenvalue Constraint (8 → 4 후보)

$D'$의 고유값 = $\omega_{R,r}^2$ (알고 있음).

각 부호 조합 $\sigma = (\sigma_{12}, \sigma_{13}, \sigma_{23}) \in \{+1,-1\}^3$에 대해
$D'_\sigma$를 구성하고 고유값이 $\omega_R^2$과 일치하는지 검사.

**Proposition (Klein four-group 대칭)**: $n=3$에서 정확히 4개 조합이 통과.

증명: 유사 변환 $D' \to SDS$ ($S = \text{diag}(\pm 1, \pm 1, \pm 1)$)는 고유값을 보존하면서
$\Gamma_{ij} \to s_i s_j \Gamma_{ij}$. Kernel은 크기 4 → 8개 조합이 2개 orbit으로 분할.

**일반화**: $n$-DOF에서 $2^{n(n-1)/2}$ 조합 → $2^{n-1}$ 후보.

| $n$ | 전체 부호 | 고유값 필터 후 | FRF sweep 수 |
|-----|----------|-------------|-------------|
| 3 | 8 | 4 | 4 |
| 4 | 64 | 8 | 8 |
| 5 | 1024 | 16 | 16 |
| 6 | 32768 | 32 | 32 |
| 7 | >10⁶ | 64 | 64 |

### Step 6: Off-Diagonal FRF Fitting (4 → 1 유일 결정)

**핵심 관찰**: 4개 후보는 동일한 대각 FRF를 생성하지만, 
**비대각 FRF $G_{ij}(j\omega)$의 전체 형상이 완전히 다름**.

각 후보 $\sigma$에 대해:
1. $D'_\sigma \to M_l^{-1}{}_\sigma \to M_{l,\sigma}$ 재구성
2. $M_{l,\sigma}$로 모델 FRF $G_{\sigma,ij}(j\omega_k)$ 계산
3. 측정값과의 NRMSE:

$$J_\sigma = \frac{\sqrt{\sum_k \sum_{i \neq j} |G_{\sigma,ij}(j\omega_k) - G_{ij}^{\text{meas}}(j\omega_k)|^2}}{\sqrt{\sum_k \sum_{i \neq j} |G_{ij}^{\text{meas}}(j\omega_k)|^2}} \times 100\%$$

4. $\sigma^* = \arg\min_\sigma J_\sigma$ 선택.

**계산량**: $O(2^{n-1} \cdot (2n)^3 \cdot N_f)$. $n=3$, $N_f=86$ ($\Delta f=0.5$ Hz): ~1초.

---

## 왜 기존 방법이 실패하는가

### Matrix Pencil의 한계

공진에서 $\tilde{G} = \text{diag}(J_m j\omega) \cdot G$의 SVD leading vector를 
D'의 고유벡터로 사용하려는 시도.

**실패 원인**: $\tilde{G}$의 SVD vector는 $\text{adj}(Z_{CL})$의 motor-motor subblock의 
rank-1 방향에 대응하는데, 이것은 $D'$의 고유벡터와 **수학적으로 동치가 아님**.
$6 \times 6$ 임피던스 행렬의 cofactor 구조가 $3 \times 3$ D' 고유값 문제와 다르기 때문.

시뮬레이션: $\sigma_1/\sigma_2 > 40000$ (극도의 rank-1 지배)에서도 MAC = 0.867.
→ 비대각항 오차 131%.

### Iterative D' Refinement의 한계

$D' \to (K, M^{-1}) \to D'$ 재구성 루프는 **대수적 항등사상**.
1회 반복에서 $\|ΔD'\|/\|D'\| = 10^{-17}$로 수렴 — 초기 오차를 수정하는 
새로운 정보원이 없음.

### 비대각 반공진 직접 측정의 한계

비대각 FRF $G_{ij}$의 반공진은 closely spaced 공진 사이에 위치할 수 있음.

| $\Delta f$ | 반공진 dip 깊이 | 검출 가능? |
|-----------|--------------|----------|
| 0.001 Hz | -34.9 dB | ✓ |
| 0.01 Hz | -11.8 dB | △ |
| 0.1 Hz | +1.4 dB | ✗ |
| 0.5 Hz | 소실 | ✗ |

---

## 시뮬레이션 검증 결과

### 시스템

7-DOF SEA arm의 J2, J4, J6 서브시스템.

| Joint | $J_m$ (×10⁻⁴) | $N$ | $k$ (Nm/rad) |
|-------|---------------|-----|-------------|
| J2 | 0.6793 | 100 | 2500 |
| J4 | 0.2639 | 50 | 1800 |
| J6 | 0.3057 | 50 | 800 |

### 부호 판별 결과 (4 poses)

| Pose | $q$ (deg) | 부호 | 정답? | Best NRMSE | 2nd NRMSE | 분리비 |
|------|----------|-----|------|-----------|----------|-------|
| 1 | (0, 90, -90) | (-,-,+) | ✓ | 0.22% | 99.4% | 449× |
| 2 | (0, 3, 30) | (-,+,-) | ✓ | 0.82% | 105.8% | 129× |
| 3 | (0, -24, 80) | (-,+,-) | ✓ | 2.27% | 101.4% | 45× |
| 4 | (0, 88, -73) | (-,-,-) | ✓ | 2.73% | 70.7% | 26× |

**전 포즈에서 부호 판별 100% 성공**. 최소 분리비 26×.

### $M_l^{-1}$ 비대각항 오차 비교

| Pose | 원소 | Method 5 | Matrix Pencil |
|------|-----|---------|--------------|
| 1 | (1,2) | **0.09%** | 131.7% |
| 1 | (1,3) | **0.03%** | 45.6% |
| 1 | (2,3) | **0.09%** | 16.6% |
| 2 | (1,2) | **0.03%** | 16745% |
| 2 | (2,3) | **0.01%** | 171% |
| 4 | (1,2) | **0.28%** | 131.3% |

### 노이즈 강건성 (Pose 1)

| SNR | 정답 NRMSE | 오답 최소 | 분리비 |
|-----|----------|---------|-------|
| Clean | 0.18% | 99.4% | 553× |
| 60 dB | 0.21% | 99.4% | 473× |
| 40 dB | 1.00% | 99.5% | 100× |
| 20 dB | 10.2% | 99.6% | 10× |

---

## 구현 파일

| 파일 | 설명 |
|-----|-----|
| `python/frf_identify_3dof.py` | Python 전체 구현 (Methods 1-5, multi-pose) |
| `matlab/..._fixed.m` | MATLAB 구현 (Methods 1-4, Matrix Pencil 버그 수정) |
| `paper/src/main.tex` | LaTeX 논문 소스 |
| `paper/FJR_offdiag_fitting_paper.pdf` | 컴파일된 논문 (4페이지) |
| `urdf/SEA_URDF/` | 로봇 URDF + 메쉬 |

### Python 실행

```bash
pip install pin numpy scipy matplotlib
cd python/
python frf_identify_3dof.py
```

출력: 4-pose 전체 결과 + `multipose_results.png` 플롯.

---

## 연결 참조

- `SKILL_algorithms.md`: Tom Oomen 알고리즘 12개 상세 레퍼런스
- `cmif_methodology.md`: CMIF 기반 공통 공진 추출 방법론
- `paper_list.md`: Oomen 논문 목록 및 분석 상태

---

## 핵심 수식 요약

1. **강성**: $k_i = J_{\text{eff},i}(\Sigma \omega_R^2 - \Sigma \omega_{AR,i}^2)$
2. **대각 M⁻¹**: $(M_l^{-1})_{ii} = (W_i - \omega_{m,i}^2) / k_i$
3. **비대각 크기**: $|\Gamma_{ij}|^2 = W_i W_j - B_k$ ($k \neq i,j$)
4. **비대각 부호**: $\sigma^* = \arg\min_\sigma \sum_{i \neq j} |G_{\sigma,ij} - G_{ij}^{\text{meas}}|^2$
5. **역관성**: $(M_l^{-1})_{ij} = \sigma_{ij} |\Gamma_{ij}| / \sqrt{k_i k_j}$
