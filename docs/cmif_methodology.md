# CMIF 기반 공통 공진주파수 추출: 3-DOF FJR 식별 강화 방법론

## 1. 문제 정의

### 1.1 현재 방법의 구조

3-DOF coupled FJR에서 motor torque → motor velocity 전달함수:

G_{ii}(s) = (1 / Jm_i · s) · det(N_i) / det(D)

**공통 분모** det(D)의 영점 = **공진주파수** ω_{R,r} (r = 1,2,3)
**채널별 분자** det(N_i)의 영점 = **반공진주파수** ω_{AR,r,i} (r = 1,2,3, 채널 i에 의존)

현재 Step 2 (K 식별)의 핵심 공식:

k_i = J_{eff,i} · (Σ_r ω²_{R,r} − Σ_r ω²_{AR,r,i})

### 1.2 문제점: 채널별 공진주파수 불일치

det(D)가 **모든 채널에서 공통**이므로, ω_{R,r}은 물리적으로 채널에 무관해야 합니다.
그러나 실제 FRF 측정에서:

- G₁₁의 공진 peak → ω̂_{R,r} from G₁₁
- G₂₂의 공진 peak → ω̂_{R,r} from G₂₂  (≠ ω̂_{R,r} from G₁₁)
- G₃₃의 공진 peak → ω̂_{R,r} from G₃₃  (≠ 위 둘)

이 불일치의 원인:

1. **Participation factor 차이**: G_{ii}에서 r번째 공진 모드의 "기여도"가 채널마다 다름. 기여도가 작은 채널에서는 peak이 약하거나 인접 모드와 겹쳐 peak 위치가 이동.
2. **Damping에 의한 peak shift**: undamped 공진주파수와 damped peak 주파수 사이의 편차가 모드 기여도에 따라 채널별로 다르게 나타남.
3. **주파수 해상도 한계**: 이산 주파수 grid에서 peak을 찾으므로 grid 간격만큼의 양자화 오차 발생. 채널마다 peak shape이 다르면 양자화 방향도 다를 수 있음.

보고서 결과에서 0.01Hz 해상도 시:
- G₁₁ 반공진: [8.61, 29.60, 31.13] Hz
- G₂₂ 반공진: [10.02, 16.05, 31.03] Hz  
- G₃₃ 반공진: [12.55, 26.43, 29.85] Hz
- 공진 (공유): [12.84, 29.60, 31.16] Hz

→ 공진주파수가 세 채널 모두에서 동일하게 보고되었지만, 해상도가 더 거칠어지면 채널별 편차가 발생합니다.

### 1.3 영향: K 및 M⁻¹ 오차 전파

Σ_r ω²_R에 오차 δ가 있으면:

k̂_i = J_{eff,i} · (Σ_r ω²_R + δ − Σ_r ω²_{AR,i})

→ δ가 K 추정에 직접 전파. 특히 Σ_r ω²_R ≈ Σ_r ω²_{AR,i}인 경우 (강성이 작아 공진-반공진 간격이 좁은 조인트) **상대오차가 급격히 증가**.

---

## 2. CMIF (Complex Mode Indicator Function)의 수학적 배경

### 2.1 정의

매 주파수점 ω에서 3×3 FRF 행렬 G(jω)의 SVD를 수행:

G(jω) = U(ω) Σ(ω) V^H(ω)

여기서:
- Σ(ω) = diag(σ₁(ω), σ₂(ω), σ₃(ω)), σ₁ ≥ σ₂ ≥ σ₃ ≥ 0
- U(ω) ∈ C^{3×3}: left singular vectors (출력 방향)
- V(ω) ∈ C^{3×3}: right singular vectors (입력 방향)

**CMIF**: σ₁(ω), σ₂(ω), σ₃(ω)를 주파수의 함수로 plot한 것.

### 2.2 CMIF와 공진주파수의 관계

3-DOF FJR의 전달함수 행렬:

G(jω) = diag(1/Jm_i · jω) · N_mat(jω) / det(D(jω))

공진주파수 ω_R에서 det(D(jω_R)) → 0이므로 G(jω)의 **모든 원소가 동시에** 발산.
따라서:

σ₁(ω_R) → ∞  (최대 특이값이 peak)

이것이 핵심: **σ₁(ω)의 peak은 det(D) = 0의 근, 즉 공통 공진주파수에 정확히 대응.**

반면, 반공진주파수 ω_{AR,i}에서는 det(N_i) = 0이므로 G_{ii}만 0에 접근하고, 
다른 채널 G_{jj} (j ≠ i)는 영향받지 않습니다. 따라서 σ₁(ω_{AR,i})에서는 
peak이 아니라 dip이 나타나지만, dip의 깊이는 채널마다 다릅니다.

### 2.3 다중 모드 분리: σ₂, σ₃의 역할

두 공진주파수가 가까이 있는 경우 (closely spaced modes):
- σ₁(ω)에서 두 peak이 병합되어 하나로 보일 수 있음
- **σ₂(ω)의 peak**이 두 번째 모드를 분리하여 보여줌

3-DOF에서 3개 공진이 있으므로:
- σ₁의 peak → 가장 지배적인 모드 (또는 서로 가까운 모드의 합)
- σ₂의 peak → 두 번째 지배적 모드
- σ₃의 peak → 세 번째 모드

→ **모든 σ_k(ω)의 peak을 종합하면 n개 공진주파수를 robust하게 추출** 가능.

### 2.4 CMIF vs 개별 채널 peak picking

| 속성 | 개별 G_{ii} peak picking | CMIF peak picking |
|------|------------------------|-------------------|
| 정보 활용 | 대각 FRF 1개 | 전체 3×3 FRF 행렬 |
| 공진 정의 | 개별 채널의 magnitude peak | 시스템 전체의 gain peak |
| Closely spaced modes | 기여도 작은 채널에서 놓칠 수 있음 | σ₂, σ₃로 분리 |
| 노이즈 robustness | 단일 채널 SNR에 의존 | 9개 전달함수의 정보 종합 |
| 물리적 일관성 | 채널별 ω_R 다를 수 있음 | **구조적으로 공통** |

---

## 3. 적용 방법론: CMIF를 현재 3-DOF 프레임워크에 접목

### 3.1 전체 파이프라인

```
[기존]
  각 G_{ii} → peak picking → ω̂_{R,r}^{(i)} (채널마다 다름) → 평균? 가중평균?

[제안]
  전체 G(jω) → SVD → σ_k(ω) → CMIF peak → ω̂_{R,r} (단일 공통값)
                                       ↓
  각 G_{ii} → (반공진만) peak picking → ω̂_{AR,r,i}
                                       ↓
  k_i = J_{eff,i} · (Σ_r ω̂²_{R,r} − Σ_r ω̂²_{AR,r,i})   ← Vieta's formula
```

### 3.2 Step-by-Step 절차

#### Step 0: MIMO FRF 측정

3개 조인트 (J2, J3, J5)에 대해 3×3 FRF 행렬 측정:

G(jω_k) ∈ C^{3×3},  k = 1, 2, ..., N_freq

측정 방법:
- 순차 가진: J2에 토크 → G₁₂, G₂₂, G₃₂ 측정; J3에 토크 → G₁₃, G₂₃, G₃₃; J5에 토크 → ...
- 또는 동시 가진: orthogonal multisine (향후 아이디어 F)

**주의**: 비대각 항 G_{ij} (i ≠ j)도 측정. 현재 방법에서는 대각만 사용하지만, CMIF에는 전체 행렬이 필요.
실제로 비대각 FRF는 순차 가진 시 "부산물"로 이미 측정 가능하므로 추가 실험 불필요.

#### Step 1: 매 주파수점에서 SVD 수행

각 ω_k에서:

[U_k, Σ_k, V_k] = svd(G(jω_k))

σ₁(ω_k) = Σ_k(1,1)
σ₂(ω_k) = Σ_k(2,2)
σ₃(ω_k) = Σ_k(3,3)

MATLAB 코드:

```matlab
for k = 1:N_freq
    [U, S, V] = svd(G(:,:,k));
    sigma1(k) = S(1,1);
    sigma2(k) = S(2,2);
    sigma3(k) = S(3,3);
end
```

#### Step 2: CMIF에서 공진주파수 추출

σ₁(ω), σ₂(ω), σ₃(ω)를 log scale로 plot.

**Peak detection**:

```matlab
% σ₁에서 peaks 찾기
[pks1, locs1] = findpeaks(20*log10(sigma1), freq, ...
    'MinPeakProminence', 10, ...  % 최소 10dB prominence
    'MinPeakDistance', 1);        % 최소 1Hz 간격

% σ₂에서 추가 peaks (closely spaced modes 분리용)
[pks2, locs2] = findpeaks(20*log10(sigma2), freq, ...
    'MinPeakProminence', 6);

% 모든 peaks를 통합하고 중복 제거
all_peaks = unique_within_tolerance([locs1; locs2], 0.5);  % 0.5Hz 내 병합
omega_R = sort(all_peaks);  % 공통 공진주파수
```

3-DOF이므로 **정확히 3개의 공진주파수**를 기대. 
3개보다 많으면 prominence가 낮은 것부터 제거, 적으면 σ₃까지 확인.

#### Step 3: 공진주파수 정밀화 — Quadratic Interpolation

이산 주파수 grid에서 찾은 peak은 grid 간격의 양자화 오차를 가집니다.
Quadratic (parabolic) interpolation으로 sub-grid 정밀도 달성:

CMIF peak 위치 ω_k에서, 인접 3점 (ω_{k-1}, ω_k, ω_{k+1})의 σ 값 (dB)을 사용:

σ̃_{-1} = 20 log₁₀ σ(ω_{k-1})
σ̃_0  = 20 log₁₀ σ(ω_k)
σ̃_{+1} = 20 log₁₀ σ(ω_{k+1})

Parabolic interpolation:

δ = (σ̃_{-1} − σ̃_{+1}) / (2(σ̃_{-1} − 2σ̃_0 + σ̃_{+1}))

ω̂_R = ω_k + δ · Δf

여기서 Δf는 주파수 grid 간격.

**기대 효과**: 0.01Hz grid에서도 ~0.001Hz 수준의 공진주파수 정밀도 달성.

```matlab
function omega_refined = parabolic_refine(sigma_dB, freq, peak_idx)
    s_m1 = sigma_dB(peak_idx - 1);
    s_0  = sigma_dB(peak_idx);
    s_p1 = sigma_dB(peak_idx + 1);
    delta = (s_m1 - s_p1) / (2*(s_m1 - 2*s_0 + s_p1));
    df = freq(2) - freq(1);
    omega_refined = freq(peak_idx) + delta * df;
end
```

#### Step 4: 반공진주파수 — 개별 대각 FRF에서 추출

반공진은 채널별로 다르므로 (det(N_i) ≠ det(N_j)), 여전히 각 G_{ii}에서 개별 추출:

```matlab
for i = 1:3
    Gii = squeeze(G(i,i,:));
    [~, locs] = findpeaks(-20*log10(abs(Gii)), freq, ...
        'MinPeakProminence', 6);  % dip = negative peak
    omega_AR(:,i) = sort(locs);
end
```

**같은 parabolic interpolation을 반공진에도 적용.**

단, 반공진 dip은 CMIF가 아닌 개별 대각 FRF의 **minimum**에서 찾아야 합니다.
(CMIF의 σ_min은 반공진과 직접 대응하지 않음 — σ_min은 시스템의 최소 gain 방향을 나타내므로.)

#### Step 5: Vieta's Formula 적용 (기존 Step 2와 동일)

공통 ω̂_R과 채널별 ω̂_{AR,i}를 사용하여:

k_i = J_{eff,i} · (Σ_{r=1}^{3} ω̂²_{R,r} − Σ_{r=1}^{3} ω̂²_{AR,r,i})

이후 Step 3 (M⁻¹ 식별)도 동일하게 진행.

### 3.3 CMIF 기반 반공진주파수 보조 검증

추가적으로, CMIF의 **left singular vector** U(ω)를 활용한 반공진 검증이 가능합니다.

공진주파수 ω_R에서의 지배적 singular vector u₁(ω_R)은 해당 모드의 **mode shape**를 나타냅니다:

u₁(ω_R) ≈ [φ₁, φ₂, φ₃]ᵀ  (각 채널의 상대적 기여도)

이 mode shape 정보를 활용하면:
- |φ_i|가 작은 채널 i에서는 해당 공진이 약하게 나타남 → 반공진과 겹칠 위험
- 이런 경우를 사전에 감지하여 반공진 peak picking의 신뢰도를 판단 가능

---

## 4. 왜 이것이 작동하는가 — 수학적 정당화

### 4.1 det(D) = 0 ↔ σ₁ → ∞

G(jω) = diag(1/(Jm_i · jω)) · adj(D(jω)) / det(D(jω))

ω → ω_R (공진)일 때 det(D(jω_R)) → 0:

‖G(jω)‖ → ∞  ⟹  σ₁(G(jω)) → ∞

즉, σ₁의 peak은 det(D)의 영점에 정확히 대응합니다.

damping이 있는 경우 det(D(jω_R)) ≠ 0이지만 매우 작으므로 σ₁에 sharp peak이 형성됩니다.

### 4.2 왜 CMIF peak이 개별 peak보다 더 정확한가

개별 G_{ii}(jω)의 magnitude는:

|G_{ii}(jω)| = |det(N_i(jω))| / (Jm_i · ω · |det(D(jω))|)

공진 ω_R에서 분모의 |det(D)|은 모든 채널에서 동일하지만, 분자 |det(N_i)|는 채널마다 다릅니다.

만약 det(N_i(jω_R)) ≈ 0 (공진과 반공진이 가까운 경우), G_{ii}에서는 peak이 상쇄되어 약해지거나 사라질 수 있습니다. 이 경우 개별 peak picking은 실패합니다.

반면, σ₁(G(jω))는 **전체 행렬의 최대 gain**이므로:

σ₁(ω) ≥ |G_{ii}(jω)|, ∀i

따라서 어떤 채널에서 공진이 약해지더라도, 다른 채널에서의 강한 공진이 σ₁에 반영됩니다.

### 4.3 Closely Spaced Modes의 분리

3-DOF에서 두 공진이 가까운 경우 (예: 29.6Hz와 31.1Hz):

G_{ii}(jω)에서는 두 peak이 하나로 합쳐져 보일 수 있지만,
σ₁(ω)과 σ₂(ω)에서는 각각의 peak이 분리됩니다.

이유: SVD는 서로 다른 singular vector 방향으로 두 모드를 직교 분해하므로,
같은 주파수 대역에 두 모드가 있어도 σ₁과 σ₂에 각각 나타납니다.

---

## 5. 예상 성능 개선 분석

### 5.1 정량적 기대

보고서의 시뮬레이션 조건 (0.01Hz 해상도)에서:

| 항목 | 개별 peak picking | CMIF + parabolic |
|------|-----------------|------------------|
| 공진주파수 오차 | ~0.01Hz (grid 한계) | ~0.001Hz (보간) |
| K 식별 오차 | 0.67% | **~0.05%** (예상) |
| (M⁻¹)_{ii} 오차 | 3.03% | **~0.3%** (예상) |
| (M⁻¹)_{ij} 오차 | 17.4% | **~2%** (예상) |

**핵심**: 공진주파수의 정확도 향상이 Vieta's formula를 통해 모든 파라미터로 전파.
특히 M⁻¹ 비대각항은 ω_R과 ω_AR의 **차이**에 의존하므로 두 주파수의 정밀도가 모두 중요.

### 5.2 Damping이 큰 경우의 이점

Damping이 크면 개별 FRF의 peak이 broadening되어 정확한 위치 추출이 어렵습니다.
CMIF의 이점:
- 9개 전달함수의 정보를 종합하므로 effective SNR 향상 (~√9 ≈ 3배)
- Closely spaced modes의 분리 능력 향상
- 약한 모드도 σ₂, σ₃에서 감지 가능

---

## 6. 구현 시 주의사항

### 6.1 비대각 FRF 측정 품질

CMIF는 전체 G(jω)를 필요로 하므로 비대각 항 G_{ij}의 품질이 중요합니다.
비대각 FRF는 일반적으로 대각보다 SNR이 낮으므로:
- 충분한 평균 횟수 (n_avg ≥ 5)
- Coherence function γ²_{ij}(ω) > 0.9 확인
- Coherence가 낮은 주파수 대역은 CMIF에서 가중치를 줄이거나 제외

### 6.2 Rigid Body Mode 처리

저주파에서 rigid body mode (s = 0 근방)가 σ₁을 지배할 수 있습니다.
공진/반공진은 rigid body 이후에 나타나므로:
- CMIF peak search 범위를 저주파 cutoff 이상으로 제한
- 또는 G(jω)에서 rigid body contribution (1/(Jm·s) term)을 미리 제거

### 6.3 3개 공진 확인

3-DOF이므로 정확히 3개 공진이 있어야 합니다. CMIF에서:
- 3개 미만 검출 → MinPeakProminence 낮추거나, σ₂, σ₃ 검사
- 3개 초과 검출 → 감쇠기에 의한 spurious peak 또는 higher-order dynamics

### 6.4 Singular Vector Tracking

주파수가 변할 때 singular vector의 순서가 바뀔 수 있습니다 (mode crossing).
이를 방지하기 위해 MAC (Modal Assurance Criterion) 기반 tracking:

MAC(u_k, u_{k+1}) = |u_k^H u_{k+1}|² / (‖u_k‖² ‖u_{k+1}‖²)

MAC > 0.9이면 같은 모드, < 0.9이면 순서 교환.

---

## 7. MATLAB 구현 요약

```matlab
function [omega_R, omega_AR] = cmif_identify(G, freq, n_modes)
% G: 3×3×N_freq complex FRF matrix
% freq: N_freq×1 frequency vector [Hz]
% n_modes: expected number of modes (= 3 for 3-DOF)

N_freq = length(freq);

%% Step 1: SVD at each frequency
sigma = zeros(3, N_freq);
U_all = zeros(3, 3, N_freq);
for k = 1:N_freq
    [U, S, ~] = svd(G(:,:,k));
    sigma(:,k) = diag(S);
    U_all(:,:,k) = U;
end

%% Step 2: Peak detection on CMIF
sigma_dB = 20*log10(sigma);

% Primary peaks from σ₁
[~, locs1] = findpeaks(sigma_dB(1,:), freq, ...
    'MinPeakProminence', 10, 'MinPeakDistance', 1);

% Secondary peaks from σ₂ (for closely spaced modes)
[~, locs2] = findpeaks(sigma_dB(2,:), freq, ...
    'MinPeakProminence', 6, 'MinPeakDistance', 1);

% Combine and select n_modes peaks
all_locs = unique_within_tol([locs1(:); locs2(:)], 0.5);
[~, sort_idx] = sort(sigma_dB(1, freq2idx(freq, all_locs)), 'descend');
omega_R_raw = all_locs(sort_idx(1:n_modes));
omega_R_raw = sort(omega_R_raw);

%% Step 3: Parabolic refinement
omega_R = zeros(n_modes, 1);
for m = 1:n_modes
    idx = find(abs(freq - omega_R_raw(m)) < 1e-6, 1);
    omega_R(m) = parabolic_refine(sigma_dB(1,:), freq, idx);
end

%% Step 4: Anti-resonance from individual diagonal FRFs
omega_AR = zeros(n_modes, 3);
for i = 1:3
    Gii_dB = 20*log10(abs(squeeze(G(i,i,:))));
    [~, locs] = findpeaks(-Gii_dB, freq, ...
        'MinPeakProminence', 6, 'MinPeakDistance', 1);
    % Parabolic refinement on each
    for m = 1:min(n_modes, length(locs))
        idx = find(abs(freq - locs(m)) < 1e-6, 1);
        omega_AR(m,i) = parabolic_refine(-Gii_dB, freq, idx);
    end
end

end

function omega_ref = parabolic_refine(signal, freq, idx)
    if idx <= 1 || idx >= length(freq)
        omega_ref = freq(idx);
        return;
    end
    s_m = signal(idx-1);
    s_0 = signal(idx);
    s_p = signal(idx+1);
    denom = 2*(s_m - 2*s_0 + s_p);
    if abs(denom) < 1e-12
        omega_ref = freq(idx);
        return;
    end
    delta = (s_m - s_p) / denom;
    df = freq(2) - freq(1);
    omega_ref = freq(idx) + delta * df;
end
```

---

## 8. 요약: 기존 대비 변경점

| 단계 | 기존 | CMIF 적용 후 |
|------|------|------------|
| FRF 측정 | 대각 G_{ii}만 | 전체 3×3 G(jω) (추가 비용 거의 없음) |
| 공진 추출 | 각 G_{ii}에서 독립 peak | CMIF σ_k(ω)에서 공통 peak |
| Peak 정밀화 | Grid 양자화 그대로 | Parabolic interpolation |
| 반공진 추출 | 각 G_{ii}에서 (변경 없음) | 각 G_{ii}에서 + parabolic refinement |
| K 식별 | Vieta's formula (변경 없음) | 동일, 더 정확한 ω 입력 |
| M⁻¹ 식별 | 계수비교 (변경 없음) | 동일, 더 정확한 ω 입력 |

**추가 계산 비용**: SVD × N_freq ≈ O(27 · N_freq) flops — 무시할 수 있는 수준.
**추가 실험 비용**: 비대각 FRF 측정이 필요하지만, 순차 가진 시 이미 측정됨. 추가 비용 = 0.
