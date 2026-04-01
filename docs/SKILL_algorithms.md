---
name: tom-oomen-sysid-algorithms
description: >
  Tom Oomen (TU Eindhoven) 교수의 system identification 알고리즘의 수식적·알고리즘적 세부사항을 담은 스킬.
  주파수 영역 FRF 식별, MIMO 실험 설계, coprime factor identification, robust control 연계,
  LPV FRF 식별, nD-LPM 등 핵심 알고리즘의 수학적 구조와 구현 패턴을 포함.
  트리거: "Oomen 알고리즘", "FRF 식별 수식", "coprime factor identification",
  "control-relevant identification", "ETFE 확장", "LPV FRF", "nD-LPM",
  "optimal experiment design MIMO", "rank-constrained SDP", "element-wise constraints",
  "dual-Youla uncertainty", "MIMO FRF identification algorithm" 등.
  5편의 논문 전문 분석에 기반한 알고리즘 수준의 스킬.
---

# Tom Oomen System Identification — 알고리즘 상세 스킬

## 개요

이 스킬은 Tom Oomen 교수의 5편의 핵심 논문을 전문 분석하여 추출한 
**수식적·알고리즘적** 세부사항을 정리한 것입니다.

분석 논문:
1. **Connecting SysID and Robust Control** (TCST 2014, Award paper)
2. **FRF ID for MV motion control: OED with element-wise constraints** (Mechatronics 2020)  
3. **FRF ID of periodically scheduled LPV systems** (MSSP 2021)
4. **Accurate FRF ID of LPV systems: nD-LPM** (TCST 2017)
5. **LTI Model ID for Mechatronic Systems Based on MIMO FRF Data** (TMECH 2023)

---

## 알고리즘 1: Closed-Loop FRF 식별 — ETFE 기반

### 문제 설정

Closed-loop에서 plant P_o를 식별. 제어기 C^exp가 시스템을 안정화.

### 실험 구조

```
r1 (excitation) → C^exp → u → P_o → y
                    ↑                  |
                    └──── feedback ────┘
```

### 주파수 영역 데이터 모델

j번째 실험에서:

[Y^<j>; U^<j>] = [P_o; I] * (I + C^exp * P_o)^{-1} * R_1^<j>

여기서 R_1^<j> ∈ C^{n_u×1}은 excitation의 DFT.

### Closed-loop transfer function T(P,C) 추정

2개 실험 (j=1,2)으로:

T̃(P_o, C^exp) = [Y^<1> Y^<2>; U^<1> U^<2>] * [R_1^<1> R_1^<2>]^{-1} * [C^exp I]

유효 조건: det([R_1^<1> R_1^<2>]) ≠ 0 (→ 입력 공간을 span해야 함)

### Open-loop 복원

T(P,C) = [T11 T12; T21 T22]이면:
P̃_o = T̃_12 * T̃_22^{-1}

### 입력 신호 설계

Random-phase multisine:
r_1(t) = Q * Σ_k a_k sin(ω_k t + φ_k)

- ω_k ∈ Ω^id (식별 주파수 그리드)
- Q = I (순차 x, y 방향 가진)
- 주기적 신호 → windowing 불필요, bias error 제거

---

## 알고리즘 2: Robust-Control-Relevant Identification

### 핵심 아이디어

식별 기준을 제어 목적에 직접 연결:

**Control-relevant identification criterion:**

P̂ = argmin_P ‖W(T(P_o, C^exp) - T(P, C^exp))V‖_∞

이것은 "좋은 모델"이 아닌 "좋은 제어기를 설계할 수 있는 모델"을 찾는 것.

### Coprime Factor Identification으로의 변환

위 H∞ 문제를 coprime factor 도메인으로 변환:

min_{N̂,D̂} ‖W([N_o; D_o] - [N̂; D̂])Ñ_e‖_∞  s.t. N̂, D̂ ∈ RH_∞

여기서:
- P̂ = N̂ D̂^{-1} (right coprime factorization)
- [N_o; D_o]는 true system의 coprime factor (비모수적으로 추정)
- Ñ_e는 co-inner → H∞ norm에 영향 없이 제거 가능

**핵심 결과**: 4-block control-relevant 문제가 2-block coprime factor 문제로 축소됨.

### 비모수적 coprime factor 추정

[Ñ_o; D̃_o] = T̃(P_o, C^exp) * V * Ñ_e*   (for ω_i ∈ Ω^id)

### 파라메트릭 모델 구조

P̂(θ) = N̂(θ) D̂(θ)^{-1} = B(θ) A(θ)^{-1}

- B, A ∈ R^{2×2}[z]: 다항식 행렬 (matrix fraction description, MFD)
- Full polynomial form → 채널 간 공통 동역학 포착 → 최소 McMillan degree

**wafer stage 결과**: McMillan degree 8의 모델로 rigid body 4 states + 2 resonances 4 states

### Dual-Youla-Kučera Uncertainty Structure

Model set:
P^DY = {P | P = (N̂ + D_c Δ_u)(D̂ - N_c Δ_u)^{-1}, Δ_u ∈ Δ_u}

**핵심 성질**: worst-case performance가 Δ_u의 affine 함수:
J_WC(P^DY, C^exp) = sup_{Δ_u} ‖M̂_22 + M̂_21 Δ_u M̂_12‖_∞

→ 항상 유계 (additive/multiplicative uncertainty와 달리)

### Robust-Control-Relevant Model Set P^RCR

특정 coprime factorization 선택 시:

J_WC(P^RCR, C^exp) ≤ J(P̂, C^exp) + γ

여기서 γ = sup σ̄(Δ_u(ω_i))는 validation 기반으로 결정.

**핵심**: γ를 최소화하는 것이 nominal model의 control-relevant identification과 정확히 일치.

### Robust Controller Synthesis

Unstructured perturbation → μ-synthesis의 D-step이 **보수성 없이** 해결 가능 (μ-simple).
이것이 structured uncertainty를 사용하는 기존 방법 대비 핵심 장점.

실험 결과: PID 40Hz bandwidth → robust controller 69Hz bandwidth, 오차 표준편차 ~50% 감소.

---

## 알고리즘 3: MIMO Optimal Experiment Design (Element-wise Constraints)

### 문제 정의

Closed-loop에서 ETFE의 covariance 최소화:

ETFE: Ĝ(k) = Y(k) U^{-1}(k)

Covariance:
C_{Ĝ}(k) = (S(k) Σ_e W^[e](k) W^[e]H(k) S^H(k))^{-1} ⊗ C_Y(k)

여기서 S(k) = (I + K(k)G(k))^{-1} (input sensitivity)

### Cost function (A-optimality)

J(W) = Σ_k Tr(M(k) C_{Ĝ}(k))

### Element-wise power constraints

각 actuator/sensor 채널별 독립적 power 제약:

P^[e]_{ξ_i} = Σ_k |Ξ^[e]_i(k)|^2 ≤ c_{ξ_i}, ∀i,e

### 명시적 OED 문제 (NLP)

minimize Σ_k γ(k) Tr(S(k) Σ_e W^[e](k) W^[e]H(k) S^H(k))^{-1}
s.t. Σ_k ⟨H_{ξ_i}(k), W^[e](k) W^[e]H(k)⟩ ≤ c_{ξ_i}, ∀i,e

여기서 H_{ξ_i}(k) = G^H_{ξ_i}(k) G_{ξ_i}(k)

**이 문제는 non-convex, NP-hard.**

### 해법 1: Rank-Constrained SDP → Semi-Definite Relaxation

Lifting: Φ^[e]_w(k) = W^[e](k) W^[e]H(k), rank(Φ) = 1

Rank constraint 제거 → convex SDR:

minimize Σ_k γ(k) Tr(Z(k))
s.t. [Z(k), S^{-H}(k); S^{-1}(k), Σ_e Φ^[e]_w(k)] ≽ 0, ∀k
    Σ_k ⟨H_{ξ_i}(k), Φ^[e]_w(k)⟩ ≤ c_{ξ_i}, ∀i,e

SDR의 optimal cost f†_SDR은 원래 문제의 lower bound.

### 해법 2: SSDR (Sequential Semi-Definite Relaxation) Algorithm

Rank constraint를 eigenvalue 기반 affine constraint로 대체:
⟨W, X⟩ ≤ ε, W = V_2 V_2^H (smallest n_u-1 eigenvectors)

반복적으로:
1. Convex SDP 풀기 (SSDR subproblem)
2. Direction matrix W 업데이트
3. Slack variable ε 단조 감소
4. **수렴 보장**: local minimizer로 수렴 (Theorem 7)

### 해법 3: Relaxation and Randomization (RR) Algorithm

1. SDR을 한 번만 풀기 → Φ†_w(k) = Δ†(k) Δ†H(k) (Cholesky)
2. Random Haar unitary R(k) 생성
3. W(k) = Δ†(k) R(k) √T (feasibility 보장을 위한 scaling)
4. 여러 번 반복, 최소 cost 선택

**Performance bound** (Theorem 9):
f*_NLP ≤ (1 + √(max σ²{P^[e]_{ξ_i}} · (n_u n_ξ - 1))) · f†_SDR

**대규모 N에서의 기대값 bound** (Theorem 11):
E{f*_NLP} ≤ (1 + √(2 · (n_u-1)/(n_u+1) · ln(n_u n_ξ))) · f†_SDR

### 실험 결과 (7×8 wafer stage)

| Method | #simultaneous inputs | Improvement factor | Computation time |
|--------|---------------------|--------------------|-----------------|
| Preliminary SIMO | 1 | 1 | - |
| Optimized SIMO | 1 | 2.5 | 1s |
| SSDR | 4 | 6.3 | 774s |
| **RR** | **8** | **14.2** | **125s** |

---

## 알고리즘 4: LPV FRF Identification (Global Approach)

### 대상 시스템 클래스

DT SISO LPV input-output model:

a(ρ(t), q^{-1}) y(t) = b(ρ(t), q^{-1}) u(t)

여기서:
a(ρ,q^{-1}) = a_0(q^{-1}) + Σ_{i=1}^{n_φ} φ_i(ρ(t)) a_i(q^{-1})
b(ρ,q^{-1}) = b_0(q^{-1}) + Σ_{i=1}^{n_ψ} ψ_i(ρ(t)) b_i(q^{-1})

### 핵심 가정

- ρ(t)는 N_ρ-periodic scheduling signal
- u(t)는 N_u-periodic 입력
- → LPV 시스템이 LPTV (Linear Periodic Time-Varying)로 동작

### 주파수 영역 표현 (Harmonic Transfer Function)

DFT 적용 후:

A₀Y(k) + Σ_i (Φ_i(k) ⊛ A_i Y(k)) = B₀U(k) + Σ_i (Ψ_i(k) ⊛ B_i U(k))

여기서 ⊛은 circular convolution. Lifted 형태:

A·Y = B·U,  G = A^{-1}B ∈ C^{N×N}

G는 HTF (Harmonic Transfer Function)의 DFT grid 위의 값.

### Identifiability 조건 (Theorem 1)

G(θ)가 globally identifiable ⟺

1. rank([P₀, Φ₁, ..., Φ_{n_φ}]) = n_φ + 1
2. rank([P₀, Ψ₁, ..., Ψ_{n_ψ}]) = n_ψ + 1  
3. [θ]_s = 1 for some s ∈ [1,N]

여기서 P₀ = [1, 0, ..., 0]^T, Φ_i는 scheduling 신호의 DFT.

**물리적 의미**: scheduling 신호가 충분히 informative해야 함 (주파수 성분이 독립).

### Parameter Estimation: WNLS

E(θ) = Y - G(θ)U (output error)

θ̂ = argmin_θ Σ_i E^{⟨i⟩}(θ)^H E^{⟨i⟩}(θ)

여기서 E̲ = W·E, W는 weighting matrix.

**Maximum Likelihood**: W W^T = C_V^{-1}이면 MLE.
→ consistent, asymptotically efficient, Cramér-Rao bound 달성.

### Jacobian (Theorem 2) — Kronecker product 제거

J^{⟨i⟩}(θ) = -W A^{-1}(θ) [T_J^{⟨i⟩}(θ)]

T_J = [-T_U(Z^{⟨i⟩}), T_W(U^{⟨i⟩})]

여기서 Z^{⟨i⟩} = A^{-1}B U^{⟨i⟩} (simulated output)

**핵심**: Kronecker product가 나타나지 않아 메모리/계산 효율적.

### 최적화 전략

1. **IV (Instrumental Variable) method**: linearized cost의 반복 해
   - 장점: local minima 통과 가능
   - 단점: global 수렴 보장 없음
2. **Levenberg-Marquardt**: IV 결과를 초기값으로 GN 기반 정제
   - damping으로 global/local 수렴 균형

### 실험 결과

- Flexible beam + LPV controller (parameter-varying stiffness)
- Scheduling: ρ(t) = sin(8πt/N) + sin(10πt/N) → 불안정 영역을 주기적으로 방문
- **BFR (Best Fit Ratio)**: identification 95.86%, validation 92.12%
- 기존 ETFE (frozen): 95.24% (LTI에서만), LPV FRF는 시변 예측까지 가능

---

## 알고리즘 5: nD-LPM (n-Dimensional Local Polynomial Method)

### 핵심 아이디어

LPM이 주파수 방향의 smoothness를 활용하듯, **scheduling 방향의 smoothness도 동시 활용**.

### Standard LPM (1D, 주파수 방향만)

**Step 1**: Transient 추정 (non-excited frequencies)

Y(αP+r) = T_y(Ω_{αP}) + Σ_{x=1}^{Q_t} t_x r^x + Ṽ(αP+r)

→ least squares로 T̂_y 추정

**Step 2**: Plant 추정 (excited frequencies)

Y_c((α+w)P) = (G(Ω_{αP}) + Σ_{x=1}^{Q_g} g_x w^x) U((α+w)P) + Ṽ

→ 다시 least squares로 Ĝ 추정

### nD-LPM 확장 (주파수 + scheduling)

nD Taylor expansion:

G(Ω_{(α+w)P}, θ_β + z) = Σ ... Σ (w^{i1} Π z_x^{i_{x+1}}) / (Π i_x!) 
                          × ∂^{Σi_x} G / ∂Ω^{i1} Π ∂θ_x^{i_{x+1}}

**2D 특수 경우** (1 scheduling parameter):

G(Ω_{(α+w)P}, θ₁+z) = G(Ω_{αP}, θ₁) + Σ_{i1=1}^{Q_z} Σ_{i2=1}^{Q_g-i1} g_{i1,i2} w^{i2} z^{i1}

### Closed-loop 확장

Z(k,θ) = [Y(k,θ); U(k,θ)]

G_wz: w → z = [G; I](I + CG)^{-1}

Open-loop 복원: G = G_wy G_wu^{-1}

### Covariance 분석

Noise covariance:
^n C_z = (μ/q) R_t R_t^H

Total covariance (noise + nonlinearities):
^t C_z = (1/q) R_g R_g^H

Plant covariance 매핑:
Ĉ_{G_wz} = S̄^H S ⊗ Ĉ_z

Open-loop 매핑:
Ĉ_G = (G_wu ⊗ [I, -G]) Ĉ_{G_wz} (...)^H

### Bias-Variance Trade-off

- nD-LPM: 더 많은 데이터 → **variance 감소**
- scheduling 방향 interpolation → **bias 증가 가능**
- polynomial order Q_z로 조절: Q_z 증가 시 bias 감소, 
  variance 감소가 멈추면 trade-off 균형점

### 실험 결과 (Medical X-ray system)

- 2×2 MIMO + 3 acceleration outputs
- 13 frozen positions, θ₂ ∈ [-90°, 90°]
- **결과**: 2D-LPM의 variance가 LPM 대비 ~20dB 감소 (동일 데이터)

---

## 알고리즘 6: MIMO LTI Model Identification from FRF Data

### 모델 구조

G_{k_o k_i}(s) = [Σ_{k=1}^{n_c} (β_k s + α_k)/(s²+2ζ_k ω_k s+ω_k²) 
                + Σ_{k=1}^{n_r} γ_k/(s+p_k) + δ₀ + δ₁s + ...] · e^{-T_d s}

**핵심 제약**: 모든 i/o 채널에서 pole 위치 공통 (physical system의 characteristic equation).

### 6-Step Algorithm

**Step 1: Delay 추정 및 제거**
- τ_d 범위에서 grid search
- H_i = [exp(-jω_i T_d)]^{-1} H'_i

**Step 2: Lightly-damped modes (ω_k, ζ_k) 추정**
- CMIF (Complex Mode Indicator Function)로 공진 주파수 탐지
- 각 모드 근방에서 local LS fitting:

  G(ω) ≈ (jωβ + α)/(ω_k² - ω² + j2ζ_k ω_k ω) + (r + jq)

- 변수 변환: u = ω_k², v = 2ζ_k ω_k
- **MIMO 확장**: θ = [θ_uv (공통); ξ^{koki} (채널별)]
  → regressor를 block-diagonal로 구성, LS로 공통 u, v 추정

**Step 3: Participation factors (α, β) 추정**
- 각 모드에서 residual dynamics를 (jω)^{-2}ā + (jω)^{-1}b̄ + c̄ + (jω)d̄로 보정
- 보조 FRF: G₂(ω) = 1/(ω_k² - ω² + j2ζ_k ω_k ω)
- LS로 α, β 추정 (ā, b̄, c̄, d̄는 auxiliary → 버림)

**Step 4: Remainder dynamics — Rational Fraction Polynomial**
- 식별된 공진 모드 기여분 제거:
  H_rem(ω) = H(ω) - Σ_k Ĉ_k(ω)
- 나머지를 RFP로 fitting: T(s) = (b₀s^m + ...)/( s^n + a₁s^{n-1} + ...)
- **MIMO**: 분모 공통, 분자 채널별 독립
  → block 구조의 LS 문제

**Step 5: 전체 통합 및 participation factor 재최적화**
- 모든 pole (complex + real)과 direct terms를 통합
- 전 주파수 대역에서 participation factors 동시 LS fitting
- Regressor: Φ = [Φ_{αβ} | Φ_γ | Φ^{koki}_δ]

**Step 6: Pole 위치 비선형 최적화**
- 탐색 변수: x = [ω₁, ζ₁, ..., ω_{nc}, ζ_{nc}, p₁, ..., p_{nr}]
- 각 후보 pole set에 대해 Step 5로 최적 participation factors 계산
- Objective: J = rms(E), E = 전 채널 concatenated error
- fmincon (SQP) + multiple starting conditions

### 실험 결과 (T-type gantry)

- tfest 대비 RMS fitting error **~100배 개선** (0.001 vs 0.22 mm/V)
- tfest: 92 states (pole 비공통) vs 제안 방법: 46 states (pole 공통)

---

## Oomen 연구의 알고리즘적 Meta-Patterns

### Pattern A: Convex Relaxation → Non-convex 문제 해결

모든 논문에서 반복:
1. 원래 문제는 non-convex (rank constraint, H∞ norm, etc.)
2. Convex relaxation으로 lower bound + 초기해
3. 반복적 refinement (SSDR, LM, global search)

### Pattern B: MIMO에서 "공통 구조" 강제

- Pole commonality (알고리즘 6)
- 공통 coprime factor denominator (알고리즘 2)
- Element-wise constraints에서 공통 spectrum (알고리즘 3)

### Pattern C: 주파수 영역의 국소 구조 활용

- LPM/LRM: 주파수 방향 Taylor expansion
- nD-LPM: scheduling 방향까지 확장
- Coprime factors: 주파수별 uncertainty quantification

### Pattern D: Least Squares의 창조적 활용

거의 모든 알고리즘의 inner loop이 LS:
- Step 2, 3, 5의 participation factor estimation
- LPM의 transient/plant estimation
- LPV FRF의 WNLS

Non-convex 최적화는 pole 위치에만 집중, 나머지는 LS로 해결.

---

## FJR/로봇 제어 연구에의 구체적 적용 가이드

### 1. FJR 공진/반공진 식별

**알고리즘 6 (MIMO FRF ID) 직접 적용**:
- 모터 토크 → 링크 위치/모터 위치의 2×1 또는 2×2 FRF 측정
- CMIF로 공진/반공진 주파수 탐지
- Step 2에서 ω_k, ζ_k 공통 추정 → 커플링된 모드 식별
- Step 4의 RFP로 rigid body + damped dynamics 포착

### 2. 강성 추정을 위한 structured identification

**알고리즘 4 (LPV FRF)의 아이디어 활용**:
- Payload/configuration 변화 → scheduling parameter ρ
- a(ρ, q^{-1})의 구조에 강성 변화 모델 내장
- Global experiment로 전체 configuration 범위 식별

### 3. DOB 설계를 위한 nominal model + uncertainty

**알고리즘 2 (Robust-Control-Relevant ID) 적용**:
- Coprime factor ID로 DOB의 nominal model 추정
- Dual-Youla uncertainty로 Q-filter 설계 가이드
- γ bound가 DOB의 robustness margin 직접 제공

### 4. 다축 로봇의 실험 설계

**알고리즘 3 (OED) 적용**:
- 각 관절 토크 제한 → element-wise constraints
- RR algorithm으로 최적 multivariable excitation 설계
- 실험 시간 대폭 절감 (SIMO 대비 ~14배 FRF 품질 향상)

---

## 참고: 핵심 수식 용어 사전

| 기호 | 의미 |
|------|------|
| S(k) | Input sensitivity (I + KG)^{-1} |
| Φ_w^[e](k) | Excitation spectrum W^[e] W^[e]H |
| {N̂, D̂} | Right coprime factorization of P̂ |
| Ñ_e | Co-inner left coprime factor of [C^exp V₂, V₁] |
| HTF | Harmonic Transfer Function (LPTV의 주파수 표현) |
| CMIF | Complex Mode Indicator Function (SVD of FRF) |
| RFP | Rational Fraction Polynomial |
| BFR | Best Fit Ratio (simulation accuracy metric) |
| A-optimality | Tr(C_{Ĝ}) 최소화 |
