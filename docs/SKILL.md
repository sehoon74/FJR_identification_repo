---
name: tom-oomen-sysid
description: >
  Tom Oomen (TU Eindhoven) 교수의 system identification 알고리즘 설계 철학과 기법을 체계화한 스킬.
  주파수 영역 식별, 제어 지향 모델링, robust control 연계, 비인과적 학습 제어에 대한 연구 패턴을 포함.
  로봇/메카트로닉스 시스템의 식별 및 제어 설계 시 Oomen 스타일의 접근법을 적용하고 싶을 때 사용.
  트리거: "Oomen 스타일", "주파수 영역 식별", "control-oriented identification",
  "FRF 기반 식별", "robust control을 위한 식별", "ILC를 위한 모델링",
  "비모수적 식별", "local rational model", "kernel-based identification",
  "multirate identification", "structured identification" 등의 키워드가 등장할 때 사용.
---

# Tom Oomen System Identification 스킬

## 개요

Tom Oomen 교수(TU Eindhoven, Control Systems Technology)의 system identification 연구에서
반복적으로 나타나는 알고리즘 설계 철학, 핵심 기법, 그리고 실험-설계-제어를 잇는 프레임워크를 정리한 스킬.

이 스킬은 논문 전문 분석이 아닌 **publication 목록 기반 패턴 추출**임을 명시합니다.
더 깊은 분석이 필요하면 `references/paper_list.md`의 논문 PDF를 직접 읽고 보완할 것.

---

## 선별된 핵심 논문 20편 (System Identification 중심, 2018~2026)

아래는 Oomen 교수의 journal 논문 중 system identification과 직접 관련된 최근 20편입니다.

| # | Year | Title (약칭) | Journal | 핵심 키워드 |
|---|------|-------------|---------|------------|
| 1 | 2026 | Structured identification of multivariable modal systems | MSSP | MIMO modal, structured parametric |
| 2 | 2026 | Closed-loop data-enabled predictive control ≡ subspace predictive control | Automatica | data-driven, subspace, closed-loop |
| 3 | 2025 | Recursive identification of structured systems: IV approach | EJC | recursive IV, structured models |
| 4 | 2025 | Locating nonlinearities: frequency-domain dynamic network | MSSP | nonlinearity detection, network |
| 5 | 2025 | Identification of additive CT systems (open/closed-loop) | Automatica | continuous-time, additive structure |
| 6 | 2025 | Local rational modeling beyond Nyquist: wafer stage | TCST | local rational, beyond-Nyquist |
| 7 | 2025 | Performance analysis of multirate systems: direct FD approach | MSSP | multirate, frequency-domain |
| 8 | 2025 | Lifted FD identification of closed-loop multirate: HDD | Mechatronics | lifted, closed-loop, multirate |
| 9 | 2025 | System identification beyond Nyquist: kernel-regularized | CEP | kernel, beyond-Nyquist, regularization |
| 10 | 2025 | Identification of additive multivariable CT systems | L-CSS | MIMO, additive, continuous-time |
| 11 | 2025 | Non-parametric system norm estimation of MIMO systems | CEP | H∞ norm, non-parametric, MIMO |
| 12 | 2025 | GraFIT: FRF identification for gravitational wave detectors | Rev. Sci. Instrum. | FRF toolbox, fast identification |
| 13 | 2024 | Sampling in parametric/nonparametric sysid: aliasing & consistency | L-CSS | sampling theory, aliasing |
| 14 | 2024 | Kernel-based identification using Lebesgue-sampled data | Automatica | kernel, event-triggered sampling |
| 15 | 2024 | Statistical analysis of block coordinate descent for CT sysid | L-CSS | continuous-time, statistical |
| 16 | 2024 | Reset-free data-driven gain estimation: reversed-circulant | Automatica | data-driven, H∞ norm, no reset |
| 17 | 2023 | LTI model identification from MIMO FRF data | TMECH | MIMO, FRF, mechatronic |
| 18 | 2023 | Wavelet-based FRF identification from incomplete data | TIM | wavelet, missing data, FRF |
| 19 | 2023 | Beyond Nyquist in FRF identification: slow-sampled systems | L-CSS | dual-rate, beyond-Nyquist |
| 20 | 2018 | Non-parametric identification of MV systems: local rational | MSSP | local rational, MIMO, benchmark |

---

## 핵심 설계 철학 (Design Philosophy)

### 1. "Identification FOR Control" — 제어 지향 식별

Oomen의 가장 근본적인 원칙. 식별의 목적은 좋은 모델 자체가 아니라 **좋은 제어기를 설계할 수 있는 모델**을 얻는 것.

- **Uncertainty model 동시 추정**: 모델과 함께 model uncertainty (additive/multiplicative)를 
  주파수 영역에서 직접 추정 → robust control (H∞, μ-synthesis)에 바로 연결
- **Control-relevant validation**: 식별된 모델의 품질을 prediction error가 아닌
  **달성 가능한 제어 성능** 관점에서 평가
- **9-block control problem과의 연계**: inferential control (측정하지 않는 성능 변수 제어)을
  위한 식별 프레임워크

> **적용 지침**: 식별 목적을 명시하라. "이 모델로 어떤 제어기를 설계할 것인가?"를 
> 식별 실험 설계 단계에서부터 고려할 것.

### 2. 주파수 영역 우선 (Frequency-Domain First)

시간 영역 데이터를 바로 파라메트릭 모델에 피팅하는 대신, 
**비모수적 FRF 추정 → 구조 분석 → 파라메트릭 피팅** 의 2단계 접근을 선호.

- **이유**: FRF는 모델 구조 가정 없이 시스템의 동특성을 직접 보여줌
- **Local Rational Model (LRM)**: FRF 추정에서 leakage/transient 영향을 줄이기 위한
  국소 다항식/유리함수 모델링 — Oomen의 핵심 도구
- **주파수 영역에서의 uncertainty quantification**: coherence 함수 외에
  LRM의 variance 추정으로 모델 불확실성을 주파수별로 정량화

> **적용 지침**: 시간 영역 데이터를 바로 `tfest`에 넣지 말 것.
> 먼저 FRF를 추정하고 시각적으로 확인한 뒤, 모델 차수와 구조를 결정하라.

### 3. 구조적 사전 지식의 적극적 활용 (Structured Identification)

"Black-box" 식별을 기본으로 하되, 알려진 물리적 구조를 가능한 한 활용:

- **Additive structure**: MIMO 시스템을 독립적인 서브시스템의 합으로 분해
  (e.g., rigid body + flexible modes)
- **Modal structure**: 기계 시스템의 modal decomposition을 식별 단계에서 강제
- **Stability by construction**: 파라미터화 자체에 안정성 제약을 내장
  (unconstrained optimization으로도 항상 안정한 모델을 얻음)
- **Continuous-time identification**: 물리적 의미가 있는 연속시간 모델을
  이산시간 데이터로부터 직접 추정 (sampling 효과를 명시적으로 처리)

> **적용 지침**: 식별 대상 시스템의 물리적 구조를 먼저 파악하라.
> Rigid body mode 수, 예상되는 공진 주파수 범위, MIMO coupling 구조 등을
> 모델 파라미터화에 반영하면 추정 정확도가 크게 향상된다.

### 4. Sampling과 Multirate의 명시적 처리

표준 식별 이론이 가정하는 "uniform sampling at Nyquist 이상"을 넘어서는 상황을 체계적으로 다룸:

- **Beyond-Nyquist identification**: 샘플링 주파수보다 높은 주파수의 동특성을
  multirate 또는 kernel-regularized 방법으로 식별
- **Lifted system approach**: multirate 시스템을 single-rate lifted 시스템으로 변환하여
  표준 식별 이론을 적용
- **Event-triggered (Lebesgue) sampling**: 비균일 샘플링 데이터에서의 kernel-based 식별
- **Aliasing 분석**: 샘플링으로 인한 aliasing이 파라메트릭/비파라메트릭 추정기에
  미치는 영향의 엄밀한 분석

> **적용 지침**: 로봇 시스템에서 센서/액추에이터의 샘플링 레이트가 다르거나,
> 성능 변수가 저속으로 측정되는 경우, multirate/lifted 프레임워크를 고려하라.

### 5. 데이터 효율성과 실험 설계 (Experiment Design)

산업 현장에서 데이터 수집 시간은 비용 → 최소한의 실험으로 최대한의 정보를 추출:

- **Optimal experiment design**: Fisher information 기반으로 입력 신호의
  주파수 스펙트럼을 최적화
- **Element-wise constraints**: MIMO 시스템에서 각 입출력 채널별로
  다른 신호 크기 제약 하에서의 최적 설계
- **Multisine excitation**: 주기적 멀티사인 입력의 설계와
  유한 샘플에서의 통계적 성질 분석

> **적용 지침**: 로봇 시스템의 식별 실험 시, 관절별 토크 한계와 안전 범위를
> 고려한 constrained optimal experiment design을 적용하라.

---

## 핵심 알고리즘/기법 카탈로그

### A. Local Rational Model (LRM) — 비모수적 FRF 추정

```
목적: 주파수 응답 함수(FRF)를 transient/leakage 영향 없이 추정
입력: 시간 영역 입출력 데이터 (주기적 또는 비주기적)
출력: FRF 추정값 + 주파수별 variance 추정

핵심 아이디어:
- 각 주파수점 ω_k 근방에서 국소 유리함수로 시스템 + transient를 동시 모델링
- 전역 파라메트릭 모델 없이도 정확한 FRF + uncertainty 제공
- MIMO 확장, LPV 확장 가능 (nD-LPM)

Oomen 논문: #20 (MSSP 2018), #6 (TCST 2025), #12 (Rev.Sci.Instrum 2025)
```

### B. Structured Parametric Identification

```
목적: 물리적 구조를 반영한 파라메트릭 모델 추정
방법 카테고리:
  - Additive decomposition: 서브시스템별 독립 추정 (연속시간)
  - Modal decomposition: 공진/반공진 모드별 구조화
  - Instrumental Variable (IV): bias-free 추정, 재귀적 확장
  - Stability-by-construction: 안정성이 보장되는 파라미터화

Oomen 논문: #1 (MSSP 2026), #3 (EJC 2025), #5 (Automatica 2025), #10 (L-CSS 2025)
```

### C. Beyond-Nyquist / Multirate Identification

```
목적: 샘플링 레이트보다 높은 주파수의 동특성 식별
방법:
  - Kernel-regularized estimation: TC/DC 커널로 impulse response 사전분포 설정
  - Lifted frequency-domain: multirate → single-rate 변환
  - Local rational model 확장: inter-sample behavior 포착

Oomen 논문: #6, #7, #8, #9 (2025), #13, #19 (2023-2024)
```

### D. Data-Driven System Norm Estimation

```
목적: 식별된 모델 없이 입출력 데이터로부터 H∞ norm 직접 추정
방법:
  - Power iteration on reversed-circulant matrices
  - Reset-free (시스템 초기화 없이 반복 실험)
  - Non-parametric MIMO 확장

Oomen 논문: #11 (CEP 2025), #16 (Automatica 2024)
```

### E. Kernel-Based Identification

```
목적: 정규화를 통한 bias-variance tradeoff 최적화
방법:
  - TC (Tuned-Correlated), DC (Diagonal-Correlated) 커널
  - Lebesgue-sampled (event-triggered) 데이터에의 확장
  - Bayesian interpretation: 사전분포로서의 커널 선택

Oomen 논문: #9 (CEP 2025), #14 (Automatica 2024)
```

---

## Oomen 연구의 반복 패턴 (Meta-Pattern)

### Pattern 1: 이론 → 산업 실증의 완결성

거의 모든 논문이 다음 구조를 따름:
1. 수학적으로 엄밀한 문제 정의
2. 알고리즘 제안 + 통계적 성질 증명 (consistency, efficiency)
3. **실제 산업 시스템에서의 실험 검증** (wafer stage, printer, HDD, 중력파 검출기 등)

### Pattern 2: 주파수 영역 ↔ 시간 영역의 이중성 활용

- 분석과 설계는 주파수 영역에서
- 실험과 데이터 수집은 시간 영역에서
- 두 도메인을 명시적으로 연결하는 프레임워크를 항상 제공

### Pattern 3: 비모수적 → 파라메트릭의 점진적 복잡화

1. 먼저 FRF (비모수적)로 시스템의 "얼굴"을 확인
2. 구조적 사전 지식을 활용한 파라메트릭 모델 구축
3. 제어 목적에 필요한 만큼만 복잡도를 증가

### Pattern 4: 식별-제어-학습의 삼각 연결

```
        Identification
           /        \
          /          \
    Control  ←→  Learning (ILC/RC)
```

식별된 모델이 ILC/RC 설계에 사용되고, ILC/RC의 반복 데이터가 다시 모델을 개선하는
**closed-loop learning** 프레임워크가 다수의 논문에서 반복됨.

### Pattern 5: 제약 하 최적성 (Constrained Optimality)

- 무제약 최적 해 → 실제 제약(안전, 에너지, 시간) 하에서의 수정
- 항상 "왜 이 제약이 필요한가"와 "제약이 성능에 미치는 영향"을 정량적으로 분석

---

## 교수님 연구와의 연결점 (로봇 제어 관점)

### FJR (Flexible Joint Robot)에의 적용 가능성

| Oomen 기법 | FJR 적용 시나리오 |
|-----------|-----------------|
| Local Rational Model | 조인트 FRF 측정에서 공진/반공진 모드 정확 추정 |
| Structured (modal) ID | MIMO 연성 진동 모드 분리 및 강성 추정 |
| Beyond-Nyquist ID | 고속 모터 드라이브 vs 저속 링크 센서의 multirate 상황 |
| Additive CT ID | Rigid body + flexibility의 분리 식별 |
| IV-based recursive ID | 온라인 강성 변화 추정 (payload 변화 시) |

### DOB 프레임워크와의 시너지

Oomen의 식별 결과물(FRF + uncertainty)은 DOB 설계에 직접 활용 가능:
- FRF에서 추정된 nominal plant → DOB의 Q-filter 설계
- Uncertainty bound → DOB의 robust stability margin 분석
- Multirate 식별 → 멀티레이트 DOB 구현의 기반

---

## 사용법

이 스킬을 활용하여:

1. **식별 실험 설계**: Oomen 스타일의 주파수 영역 중심 실험 설계를 적용할 때
2. **FRF 분석**: 측정된 FRF 데이터에서 모델 구조를 결정할 때
3. **제어 지향 모델링**: robust control이나 ILC/feedforward를 위한 모델을 만들 때
4. **논문 작성**: Oomen 그룹의 연구를 참조하거나 비교할 때

### 보완이 필요한 부분

이 스킬은 **논문 제목, 초록, 알려진 연구 패턴 기반**으로 작성되었습니다.
다음 사항을 보완하면 정확도가 크게 향상됩니다:

- [ ] 핵심 논문 5편의 PDF를 직접 읽고 알고리즘 세부사항 보강
- [ ] MATLAB/Python 코드 예제 추가 (Oomen 그룹의 공개 코드 참조)
- [ ] 각 기법의 수학적 정의 (전달함수 형태, 최적화 목적함수 등) 추가
- [ ] 실제 FJR 데이터에 적용한 case study 추가

---

## 참고 자료

- Tom Oomen's publication page: https://toomen.eu/publications.html
- Google Scholar: https://scholar.google.com/citations?user=BCMH_VYAAAAJ
- Overview paper (2018): "Advanced motion control for precision mechatronics", IEEJ JIA
- Inaugural lecture (2024): "Learning in control"
- 공개 코드/데이터: https://gitlab.tue.nl/oomenpublic/
