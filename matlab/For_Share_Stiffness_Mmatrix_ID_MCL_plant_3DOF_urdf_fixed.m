%% =======================================================================
%  IROS 이후 이론 디벨롭. URDF로 전달함수 도출 후 K 및 M행렬 식별 테스트
%  분모의 대각항의 합 = 공진주파수^2의 합의 성질을 이용
% ========================================================================
close all; clear; clc;
sympref('PolynomialDisplayStyle', 'descend');

%% 공통 설정
MANIPULATOR     = 'SEA';
MAX_DOF         = 7;
selected_joints = [2,4,6];
DOF             = length(selected_joints);
is_gear_ratio   = 1;
round_digits    = 3;          % 반올림 자릿수 (FRF 주파수 해상도 모사)

% % [단일 채널 확인용 - 필요시 주석 해제]
% TF_input  = 1;   % local index
% TF_output = 1;   % local index
% input_joint_global  = selected_joints(TF_input);
% output_joint_global = selected_joints(TF_output);
% fprintf('Subsystem Input  (Local %d) is Global Joint %d\n', TF_input,  input_joint_global);
% fprintf('Subsystem Output (Local %d) is Global Joint %d\n', TF_output, output_joint_global);

%% 로봇 파라미터
params   = getRobotParams_MCL(MANIPULATOR);
Jm_param = params.Jm;   % 1x7
Nm_param = params.Nm;   % 1x7
K_param  = params.K;    % 1x7
K_true   = K_param(selected_joints).';  % 3x1

%% -----------------------------------------------------------------------
%  심볼릭 Z_CL 및 전체 G_ij (3×3) 한 번만 구성
% -----------------------------------------------------------------------
syms s
Jm_sym_full  = sym('Jm',  [1, MAX_DOF]);
M_l_sym_full = sym('M',   [MAX_DOF, MAX_DOF]);
k_sym_full   = sym('k',   [1, MAX_DOF]);
N_sym_full   = sym('N',   [1, MAX_DOF]);

M_l   = M_l_sym_full(selected_joints, selected_joints);
J_m   = diag(Jm_sym_full(selected_joints));
K_sym = diag(k_sym_full(selected_joints));
N_sym = diag(N_sym_full(selected_joints));

M_total = [M_l,        zeros(DOF);
           zeros(DOF), J_m       ];
if is_gear_ratio
    K_total = [ K_sym,            -K_sym*inv(N_sym);
               -inv(N_sym)*K_sym,  inv(N_sym)^2*K_sym];
else
    K_total = [ K_sym, -K_sym;
               -K_sym,  K_sym];
end
Z_CL = M_total*s^2 + K_total;   % (6x6)

% 전체 G_ij 심볼릭 (motor velocity / motor torque)
G_sym = cell(DOF, DOF);
for i_out = 1:DOF
    for i_in = 1:DOF
        row_idx = DOF + i_out;   % motor side output
        col_idx = DOF + i_in;    % motor side input
        M_minor            = Z_CL;
        M_minor(row_idx,:) = [];
        M_minor(:,col_idx) = [];
        sign_ij            = (-1)^(row_idx + col_idx);
        G_pos_ij           = sign_ij * det(M_minor) / det(Z_CL);
        G_sym{i_out, i_in}(s) = s * G_pos_ij;
    end
end
fprintf('Symbolic 3x3 TF matrix constructed.\n');

% % [단일 채널 심볼릭 - 필요시 주석 해제]
% TF_input = 1; TF_output = 1;
% G_vel_sym_single(s) = G_sym{TF_output, TF_input}(s);
% fprintf('Single-channel TF (local in=%d, out=%d) ready.\n', TF_input, TF_output);

%% -----------------------------------------------------------------------
%  URDF 로봇 로드
% -----------------------------------------------------------------------
baseDir = fileparts(mfilename('fullpath'));
urdfPath = fullfile(baseDir, '..', 'urdf', 'SEA_URDF', 'urdf', 'SEA arm_woEE.urdf');
robot = importrobot(urdfPath);
robot.DataFormat = 'column';

Jm_local = Jm_param(selected_joints);
Nm_local = Nm_param(selected_joints);
Jr_local = Jm_local .* (Nm_local.^2);   % reflected inertia

%% -----------------------------------------------------------------------
%  주파수 격자 설정
% -----------------------------------------------------------------------
delta_f      = 10^(-round_digits);
f_vec        = 0.1 : delta_f : 100;
omega_vec    = 2*pi*f_vec;

%% -----------------------------------------------------------------------
%  포즈 루프: M(q) → 전체 G_ij 수치화 → FRF 저장 → Bode (대각)
% -----------------------------------------------------------------------
q_active_log = [0     0       0       0     ;
                pi/2  0.0443  -0.4154 1.5329;
               -pi/2  0.5236  1.3909 -1.2658];
N_poses = size(q_active_log,2);

M_log        = cell(N_poses, 1);
sub_log      = cell(N_poses, 1);
FRF_log      = cell(N_poses, 1);   % {p}{i_out, i_in} = complex FRF vector
omega_AR_log = cell(N_poses, 1);   % roots 기반 (대각만)
omega_R_log  = cell(N_poses, 1);   % roots 기반

fig_robot = figure('Name', 'Robot Poses', 'Position', [100,100,1200,400]);

for p = 1:N_poses
    %% 포즈 구성
    q_active = q_active_log(:, p);
    q_full   = zeros(7,1);
    q_full(selected_joints) = q_active;

    figure(fig_robot);
    subplot(1, N_poses, p);
    show(robot, q_full);
    title(sprintf('Pose %d', p));
    view([45, 30]);

    %% M(q) 계산
    M_full   = massMatrix(robot, q_full);
    M_l_num  = M_full(selected_joints, selected_joints);
    M_log{p} = M_l_num;

    %% sub_struct 구성
    sub = struct();
    for i_local = 1:DOF
        i_global = selected_joints(i_local);
        sub.(string(Jm_sym_full(i_global))) = Jm_param(i_global);
        sub.(string(k_sym_full(i_global)))  = K_param(i_global);
        sub.(string(N_sym_full(i_global)))  = Nm_param(i_global);
    end
    for i_local = 1:DOF
        i_global = selected_joints(i_local);
        for j_local = 1:DOF
            j_global = selected_joints(j_local);
            sub.(string(M_l_sym_full(i_global, j_global))) = M_l_num(i_local, j_local);
        end
    end
    sub_log{p} = sub;

    %% 전체 G_ij 수치화 + FRF 저장
    FRF_p = cell(DOF, DOF);
    omega_AR_p = zeros(DOF, DOF);
    omega_R_p  = [];

    fig_bode = figure('Name', sprintf('Bode - Pose %d', p));

    for i_out = 1:DOF
        for i_in = 1:DOF
            G_real(s)      = subs(G_sym{i_out, i_in}(s), sub);
            [num_s, den_s] = numden(G_real);
            [C_d, ~]       = coeffs(formula(den_s), s);
            lc             = C_d(1);
            num_n          = sym2poly(vpa(formula(num_s / lc), 32));
            den_n          = sym2poly(vpa(formula(den_s / lc), 32));

            % FRF 저장 (CMIF용)
            sys       = tf(num_n, den_n);
            [mag, phase_d] = bode(sys, omega_vec);
            mag_sq    = squeeze(mag);
            phase_rad = squeeze(phase_d) * pi/180;
            FRF_p{i_out, i_in} = mag_sq .* exp(1j * phase_rad);

            % 대각 채널: roots 기반 ω_AR / ω_R 추출
            if i_out == i_in
                i_diag = i_out;
                
                % 반공진
                z_d      = roots(num_n);
                ar_cands = sort(abs(imag(z_d)));
                ar_cands = ar_cands(ar_cands > 1e-3);
                ar_cands = uniquetol(ar_cands, 1e-6);
                ar_cands = round(ar_cands / (2*pi), round_digits) * (2*pi);
                omega_AR_p(i_diag, :) = ar_cands(:)';

                % 공진
                p_d     = roots(den_n);
                r_cands = sort(abs(imag(p_d)));
                r_cands = r_cands(r_cands > 1e-3);
                r_cands = uniquetol(r_cands, 1e-6);
                r_cands = round(r_cands / (2*pi), round_digits) * (2*pi);
                omega_R_p = uniquetol(r_cands, 1e-6);

                fprintf('Pose %d  G_%d%d  반공진: [', p, i_diag, i_diag);
                fprintf('%.4f ', omega_AR_p(i_diag,:) / (2*pi));
                fprintf('] Hz\n');

                %% Bode 플롯 (대각만)
                mag_db    = 20*log10(mag_sq);
                phase_deg = wrapTo180(squeeze(phase_d));

                subplot(DOF, 2, (i_diag-1)*2 + 1);
                semilogx(f_vec, mag_db, 'LineWidth', 1.5); hold on; grid on;
                ylabel(sprintf('G_{%d%d}  Mag (dB)', i_diag, i_diag));
                xlabel('Frequency (Hz)'); xlim([0.1, 100]);

                subplot(DOF, 2, (i_diag-1)*2 + 2);
                semilogx(f_vec, phase_deg, 'LineWidth', 1.5); hold on; grid on;
                ylabel('Phase (deg)'); xlabel('Frequency (Hz)'); xlim([0.1, 100]);
            end
        end
    end

    fprintf('Pose %d  공진: [', p);
    fprintf('%.4f ', omega_R_p / (2*pi));
    fprintf('] Hz\n\n');

    FRF_log{p}      = FRF_p;
    omega_AR_log{p} = omega_AR_p;
    omega_R_log{p}  = omega_R_p;
end

%% -----------------------------------------------------------------------
%  식별 블록 (루프 밖)
%  공통: p_id 포즈 선택
% -----------------------------------------------------------------------
pose_id         = 1;
omega_AR_all = omega_AR_log{pose_id};
omega_R_id   = omega_R_log{pose_id};
M_l_num_id   = M_log{pose_id};
M_inv_gt     = inv(M_l_num_id);
k_loc        = K_param(selected_joints);

%% -----------------------------------------------------------------------
%  [방법 1] roots() 직접 추출 기반 식별
% -----------------------------------------------------------------------
fprintf('===== [방법 1] roots 기반 =====\n');

% K 식별
fprintf('--- Stiffness Identification ---\n');
K_est_roots = zeros(1, DOF);
for i = 1:DOF
    K_est_roots(i) = Jr_local(i) * (sum(omega_R_id.^2) - sum(omega_AR_all(i,:).^2));
    err_pct = abs(K_est_roots(i) - K_true(i)) / K_true(i) * 100;
    fprintf('Joint %d: True = %.2f,  Est = %.2f,  Err = %.2f%%\n', ...
        selected_joints(i), K_true(i), K_est_roots(i), err_pct);
end

% Step 3: 다항식 차분
syms lambda

%     % No 정규화
%     P_W   = prod(lambda - omega_R_id.^2);
%     P_W_c = sym2poly(expand(P_W));
%     omega_m_sq = k_loc ./ (Jm_param(selected_joints) .* Nm_param(selected_joints).^2);
%     
%     A_vals = zeros(1, DOF);
%     B_vals = zeros(1, DOF);
%     for i = 1:DOF
%         P_Ni    = prod(lambda - omega_AR_all(i,:).^2);
%         P_Ni_c  = sym2poly(expand(P_Ni));
%     
%         delta_c = P_Ni_c - P_W_c;
%         A_vals(i) = -delta_c(3) / delta_c(2);
%         B_vals(i) =  delta_c(4) / delta_c(2);
%     end
    
    % 정규화
    omega_ref = geomean(omega_R_id);          % 기하평균을 정규화 기준으로
    lambda_n  = lambda / omega_ref^2;
    
    P_W_n   = prod(lambda_n - (omega_R_id  / omega_ref).^2);
    P_W_n_c = sym2poly(expand(P_W_n));
    
    omega_m_sq = k_loc ./ (Jm_param(selected_joints) .* Nm_param(selected_joints).^2);
    
    A_vals = zeros(1, DOF);
    B_vals = zeros(1, DOF);
    for i = 1:DOF
        P_Ni_n   = prod(lambda_n - (omega_AR_all(i,:) / omega_ref).^2);
        P_Ni_n_c = sym2poly(expand(P_Ni_n));
    
        delta_c   = P_Ni_n_c - P_W_n_c;
    
        % 역정규화: A×ω_ref², B×ω_ref⁴
        A_vals(i) = -delta_c(3) / delta_c(2) * omega_ref^2;
        B_vals(i) =  delta_c(4) / delta_c(2) * omega_ref^4;
    end

% Step 4: M^{-1} 복원
[Minv_diag_est_r, Minv_12_r, Minv_13_r, Minv_23_r] = ...
    recoverMinv(A_vals, B_vals, omega_m_sq, k_loc, DOF);
printMinvResult(Minv_diag_est_r, Minv_12_r, Minv_13_r, Minv_23_r, M_inv_gt, DOF);

%% -----------------------------------------------------------------------
%  [방법 2] CMIF 기반 식별
% -----------------------------------------------------------------------
fprintf('\n===== [방법 2] CMIF 기반 =====\n');

FRF_p = FRF_log{pose_id};

% 3×3 FRF 행렬 → 매 주파수점 SVD → σ₁(ω)
N_freq  = length(omega_vec);
sigma1  = zeros(N_freq, 1);
for k = 1:N_freq
    G_mat = zeros(DOF, DOF);
    for i_out = 1:DOF
        for i_in = 1:DOF
            G_mat(i_out, i_in) = FRF_p{i_out, i_in}(k);
        end                                                        
    end
    sv          = svd(G_mat);
    sigma1(k)   = sv(1);
end

% σ₁ peak → 공진 (parabolic interpolation)
sigma1_db   = 20*log10(sigma1);
[~, pk_idx] = findpeaks(sigma1_db, 'MinPeakProminence', 3);
omega_R_cmif = zeros(size(pk_idx));
for r = 1:length(pk_idx)
    omega_R_cmif(r) = parabolicPeak(sigma1_db, pk_idx(r), omega_vec);
end

% 각 G_ii notch → 반공진 (parabolic interpolation on -mag)
omega_AR_cmif = zeros(DOF, DOF);
for i_diag = 1:DOF
    mag_diag    = abs(FRF_p{i_diag, i_diag});
    mag_diag_db = 20*log10(mag_diag);
    [~, nt_idx] = findpeaks(-mag_diag_db, 'MinPeakProminence', 3);
    for r = 1:length(nt_idx)
        omega_AR_cmif(i_diag, r) = parabolicPeak(-mag_diag_db, nt_idx(r), omega_vec);
    end
end

fprintf('공진 (CMIF): ['); fprintf('%.4f ', omega_R_cmif/(2*pi)); fprintf('] Hz\n');
for i = 1:DOF
    fprintf('G_%d%d 반공진: [', i, i);
    fprintf('%.4f ', omega_AR_cmif(i,:)/(2*pi));
    fprintf('] Hz\n');
end

% K 식별
fprintf('--- Stiffness Identification ---\n');
K_est_cmif = zeros(1, DOF);
for i = 1:DOF
    K_est_cmif(i) = Jr_local(i) * (sum(omega_R_cmif.^2) - sum(omega_AR_cmif(i,:).^2));
    err_pct = abs(K_est_cmif(i) - K_true(i)) / K_true(i) * 100;
    fprintf('Joint %d: True = %.2f,  Est = %.2f,  Err = %.2f%%\n', ...
        selected_joints(i), K_true(i), K_est_cmif(i), err_pct);
end

% Step 3 & 4
P_W_c_cmif = sym2poly(expand(prod(lambda - omega_R_cmif.^2)));
A_cmif = zeros(1, DOF);
B_cmif = zeros(1, DOF);
for i = 1:DOF
    P_Ni_c  = sym2poly(expand(prod(lambda - omega_AR_cmif(i,:).^2)));
    delta_c = P_Ni_c - P_W_c_cmif;
    A_cmif(i) = -delta_c(3) / delta_c(2);
    B_cmif(i) =  delta_c(4) / delta_c(2);
end
[Minv_diag_est_c, Minv_12_c, Minv_13_c, Minv_23_c] = ...
    recoverMinv(A_cmif, B_cmif, omega_m_sq, k_loc, DOF);
printMinvResult(Minv_diag_est_c, Minv_12_c, Minv_13_c, Minv_23_c, M_inv_gt, DOF);

%% CMIF σ₁ 플롯
figure('Name', sprintf('CMIF - Pose %d', pose_id));
semilogx(f_vec, sigma1_db, 'LineWidth', 1.5); grid on;
xlabel('Frequency (Hz)'); ylabel('\sigma_1 (dB)');
title(sprintf('CMIF — Pose %d', pose_id));
xline(omega_R_cmif/(2*pi), 'r--', 'LineWidth', 1);
xlim([0.1, 100]);

%% -----------------------------------------------------------------------
%  비교 출력
% -----------------------------------------------------------------------
fprintf('\n===== 방법 1 vs 방법 2 비교 =====\n');
fprintf('%-12s %10s %10s %10s\n', '', 'roots', 'CMIF', 'True');
for i = 1:DOF
    fprintf('K_%d        %10.3f %10.3f %10.3f\n', ...
        selected_joints(i), K_est_roots(i), K_est_cmif(i), K_true(i));
end
fprintf('\n(M^{-1})_{ii}:\n');
for i = 1:DOF
    fprintf('  (%d,%d)  roots=%+.6f  CMIF=%+.6f  True=%+.6f\n', ...
        i, i, Minv_diag_est_r(i), Minv_diag_est_c(i), M_inv_gt(i,i));
end

%% -----------------------------------------------------------------------
%  [방법 3] Matrix Pencil 기반 식별
%  수정사항:
%    Fix 1: 복소 singular vector → global phase 제거 → 실수 mode shape
%    Fix 2: eig(N_i') 기반 반공진 → Vieta's formula로 K 식별
%    Fix 3: 식별된 K_est_mp를 M^{-1} 계산에 사용 (ground truth 미사용)
% -----------------------------------------------------------------------
fprintf('\n===== [방법 3] Matrix Pencil 기반 =====\n');

% --- 스케일 보정된 G_tilde 에서 공진점 right singular vector 추출 ---
V_modes = zeros(DOF, DOF);   % 열: 각 공진의 실수 mode shape
for r = 1:DOF
    % 공진 주파수에 가장 가까운 격자점
    [~, k_res] = min(abs(omega_vec - omega_R_cmif(r)));

    % G_tilde = diag(Jm_i * jω) · G  (Jm·s 스케일 제거 → adj(D)/det(D) 근사)
    G_mat_raw = zeros(DOF, DOF);
    for i_out = 1:DOF
        for i_in = 1:DOF
            G_mat_raw(i_out, i_in) = FRF_p{i_out, i_in}(k_res);
        end
    end
    scale = Jm_local(:) .* (1j * omega_vec(k_res));   % (3x1)
    G_tilde = diag(scale) * G_mat_raw;

    % SVD → 최대 singular value의 right singular vector
    [~, ~, Vsvd] = svd(G_tilde);
    v_complex = Vsvd(:, 1);

    % [Fix 1] 복소 singular vector → 실수 mode shape 변환
    % global phase를 제거: 최대 크기 원소의 위상을 기준으로 회전
    [~, max_idx] = max(abs(v_complex));
    phase_corr = angle(v_complex(max_idx));
    v_real = real(v_complex * exp(-1j * phase_corr));
    V_modes(:, r) = v_real / norm(v_real);
end

fprintf('[추출된 mode shapes (열 = 각 공진 모드)]\n');
disp(V_modes);

% --- D' 재구성 ---
Lambda_R  = diag(omega_R_cmif.^2);
D_prime   = real(V_modes * Lambda_R / V_modes);   % V * Λ * V^{-1}, 실수화

fprintf('[재구성된 D prime]\n');
disp(D_prime);

% D' 대칭성 검증 (물리적으로 대칭이어야 함: Γ_{ij} = Γ_{ji})
asym_err = norm(D_prime - D_prime', 'fro') / norm(D_prime, 'fro') * 100;
fprintf('  D'' 비대칭도: %.4f%%\n', asym_err);
if asym_err > 1
    fprintf('  [경고] D''의 비대칭도가 큽니다. Mode shape 추정 확인 필요.\n');
end

% 대칭화 (물리적 제약 강제)
D_prime = (D_prime + D_prime') / 2;

% --- D' 원소에서 W_i, Gamma_ij 추출 ---
W_mp     = diag(D_prime)';          % (1x3)
Gamma_mp = D_prime - diag(diag(D_prime));   % 비대각 부분

% --- [Fix 2] eig(N_i') → 반공진 → Vieta's formula로 K 식별 ---
fprintf('\n--- eig(N_i'') 기반 반공진 추출 ---\n');
e_mat = eye(DOF);
omega_AR_mp = zeros(DOF, DOF);   % (i, r): 채널 i의 r번째 반공진 [rad/s]
for i = 1:DOF
    N_i_prime = D_prime - omega_m_sq(i) * (e_mat(:,i) * e_mat(:,i)');
    ar_sq     = sort(real(eig(N_i_prime)));

    % 음수 고유값 진단
    if any(ar_sq < -1e-6)
        fprintf('  [경고] G_%d%d: 음수 고유값 발견 (%.4e). D'' 부정확 가능.\n', ...
            i, i, min(ar_sq));
    end
    ar_sq = max(ar_sq, 0);    % 수치 안전
    omega_AR_mp(i, :) = sqrt(ar_sq);   % rad/s

    fprintf('  G_%d%d 반공진 (MP): [', i, i);
    fprintf('%.4f ', omega_AR_mp(i,:) / (2*pi));
    fprintf('] Hz   (CMIF: [');
    fprintf('%.4f ', omega_AR_cmif(i,:) / (2*pi));
    fprintf('] Hz)\n');
end

% K 식별: Vieta's formula with Matrix Pencil 반공진
fprintf('\n--- Stiffness Identification (Matrix Pencil) ---\n');
K_est_mp = zeros(1, DOF);
for i = 1:DOF
    K_est_mp(i) = Jr_local(i) * (sum(omega_R_cmif.^2) - sum(omega_AR_mp(i,:).^2));
    err_pct = abs(K_est_mp(i) - K_true(i)) / K_true(i) * 100;
    fprintf('  Joint %d: True = %.2f,  Est = %.2f,  Err = %.2f%%\n', ...
        selected_joints(i), K_true(i), K_est_mp(i), err_pct);
end

% --- [Fix 3] 식별된 K로 M^{-1} 식별 (ground truth 미사용) ---
fprintf('\n--- M^{-1} Identification (Matrix Pencil, 식별된 K 사용) ---\n');

% 대각항: (M^{-1})_{ii} = (W_i - omega_m_sq_i) / k_i
Minv_diag_mp = (W_mp - omega_m_sq) ./ K_est_mp;

% 비대각항: Γ_{ij} = sqrt(k_i * k_j) * (M^{-1})_{ij}
%          → (M^{-1})_{ij} = Γ_{ij} / sqrt(k_i * k_j)
%          부호도 Γ_{ij}로부터 직접 결정됨!
Minv_12_mp = Gamma_mp(1,2) / sqrt(K_est_mp(1) * K_est_mp(2));
Minv_13_mp = Gamma_mp(1,3) / sqrt(K_est_mp(1) * K_est_mp(3));
Minv_23_mp = Gamma_mp(2,3) / sqrt(K_est_mp(2) * K_est_mp(3));

% --- 결과 출력 ---
fprintf('\n[대각항 (M^{-1})_{ii}]\n');
for i = 1:DOF
    true_v = M_inv_gt(i,i);
    est_v  = Minv_diag_mp(i);
    err    = abs(est_v - true_v) / abs(true_v) * 100;
    fprintf('  (%d,%d): True=%+.6f  Est=%+.6f  Err=%.2f%%\n', ...
        i, i, true_v, est_v, err);
end

fprintf('[비대각항 (M^{-1})_{ij} — 부호 포함!]\n');
off_mp = {1,2,Minv_12_mp; 1,3,Minv_13_mp; 2,3,Minv_23_mp};
for r = 1:3
    ii = off_mp{r,1}; jj = off_mp{r,2}; est = off_mp{r,3};
    true_v = M_inv_gt(ii,jj);   % 부호 포함 비교
    err = abs(est - true_v) / abs(true_v) * 100;
    fprintf('  (%d,%d): True=%+.6f  Est=%+.6f  Err=%.2f%%\n', ...
        ii, jj, true_v, est, err);
end

% 전체 M^{-1} 행렬 복원 (부호 포함)
Minv_est_mp = diag(Minv_diag_mp) + ...
    [0,          Minv_12_mp, Minv_13_mp;
     Minv_12_mp, 0,          Minv_23_mp;
     Minv_13_mp, Minv_23_mp, 0          ];
fprintf('\n[복원된 M^{-1} (Matrix Pencil, 부호 포함)]\n');
disp(Minv_est_mp);
fprintf('[Ground Truth M^{-1}]\n');
disp(M_inv_gt);

% --- 3가지 방법 비교 ---
fprintf('\n===== 방법 1 vs 방법 2 vs 방법 3 비교 =====\n');
fprintf('%-12s %10s %10s %10s %10s\n', '', 'roots', 'CMIF', 'MatPencil', 'True');
for i = 1:DOF
    fprintf('K_%d        %10.3f %10.3f %10.3f %10.3f\n', ...
        selected_joints(i), K_est_roots(i), K_est_cmif(i), K_est_mp(i), K_true(i));
end
fprintf('\n(M^{-1})_{ii}:\n');
for i = 1:DOF
    fprintf('  (%d,%d)  roots=%+.6f  CMIF=%+.6f  MP=%+.6f  True=%+.6f\n', ...
        i, i, Minv_diag_est_r(i), Minv_diag_est_c(i), Minv_diag_mp(i), M_inv_gt(i,i));
end
fprintf('(M^{-1})_{ij} (부호 포함, MP만):\n');
for r = 1:3
    ii = off_mp{r,1}; jj = off_mp{r,2}; est = off_mp{r,3};
    true_v = M_inv_gt(ii,jj);
    fprintf('  (%d,%d)  MP=%+.6f  True=%+.6f\n', ii, jj, est, true_v);
end

%% -----------------------------------------------------------------------
%  [방법 4] Iterative D' Refinement
%  
%  원리:
%    방법 3의 D'는 CMIF singular vector 기반 → 근사적.
%    식별된 K, M^{-1}로 D'를 물리 파라미터로부터 직접 재구성하면
%    더 정확한 D'를 얻을 수 있음.
%    이 D'에서 다시 eig(N_i') → ω_AR → K, M^{-1} 를 반복.
%
%  수렴 조건:
%    ‖D'_new − D'_old‖ / ‖D'_old‖ < tol
%
%  알고리즘:
%    Init: D'⁽⁰⁾ = 방법 3의 SVD 기반 D'
%    Loop k = 1, 2, ...:
%      1. eig(N_i'⁽ᵏ⁻¹⁾) → ω_AR⁽ᵏ⁾
%      2. Vieta → K⁽ᵏ⁾
%      3. D' 원소 → M⁻¹⁽ᵏ⁾ (K⁽ᵏ⁾ 사용)
%      4. K⁽ᵏ⁾, M⁻¹⁽ᵏ⁾ → D'⁽ᵏ⁾ 재구성
%      5. 수렴 체크
% -----------------------------------------------------------------------
fprintf('\n===== [방법 4] Iterative D'' Refinement =====\n');

iter_max = 50;
iter_tol = 1e-12;

% --- 초기화: 방법 3의 D'를 시작점으로 사용 ---
D_iter = D_prime;   % 방법 3에서 얻은 D' (대칭화된 것)

% 반복 이력 저장
K_hist     = zeros(iter_max, DOF);
Minv_hist  = cell(iter_max, 1);
conv_hist  = zeros(iter_max, 1);
D_hist     = cell(iter_max, 1);

for iter = 1:iter_max
    D_old = D_iter;
    
    % === Step 1: eig(N_i') → 반공진 ===
    omega_AR_iter = zeros(DOF, DOF);   % (i, r) [rad/s]
    for i = 1:DOF
        N_i_iter = D_iter - omega_m_sq(i) * (e_mat(:,i) * e_mat(:,i)');
        ar_sq = sort(real(eig(N_i_iter)));
        ar_sq = max(ar_sq, 0);
        omega_AR_iter(i, :) = sqrt(ar_sq);
    end
    
    % === Step 2: Vieta's formula → K ===
    K_iter = zeros(1, DOF);
    for i = 1:DOF
        K_iter(i) = Jr_local(i) * (sum(omega_R_cmif.^2) - sum(omega_AR_iter(i,:).^2));
    end
    
    % K가 비물리적(음수)이면 보정
    if any(K_iter <= 0)
        fprintf('  [경고] iter %d: 음수 K 발견. 이전 값 유지.\n', iter);
        K_iter = max(K_iter, 1);   % 최소 양수로 clamp
    end
    
    % === Step 3: D' 원소 → M^{-1} ===
    W_iter = diag(D_iter)';
    Gamma_iter = D_iter - diag(diag(D_iter));
    
    Minv_diag_iter = (W_iter - omega_m_sq) ./ K_iter;
    
    Minv_offdiag_iter = zeros(DOF, DOF);
    for i = 1:DOF
        for j = i+1:DOF
            Minv_offdiag_iter(i,j) = Gamma_iter(i,j) / sqrt(K_iter(i) * K_iter(j));
            Minv_offdiag_iter(j,i) = Minv_offdiag_iter(i,j);
        end
    end
    Minv_iter = diag(Minv_diag_iter) + Minv_offdiag_iter;
    
    % === Step 4: K, M^{-1} → D' 재구성 ===
    % D'(i,i) = W_i = omega_m_sq_i + K_i * (M^{-1})_{ii}
    % D'(i,j) = Γ_{ij} = sqrt(K_i * K_j) * (M^{-1})_{ij}
    D_new = zeros(DOF, DOF);
    for i = 1:DOF
        D_new(i,i) = omega_m_sq(i) + K_iter(i) * Minv_iter(i,i);
        for j = i+1:DOF
            D_new(i,j) = sqrt(K_iter(i) * K_iter(j)) * Minv_iter(i,j);
            D_new(j,i) = D_new(i,j);
        end
    end
    
    % 대칭화 (수치 안전)
    D_iter = (D_new + D_new') / 2;
    
    % === Step 5: 수렴 체크 ===
    rel_change = norm(D_iter - D_old, 'fro') / norm(D_old, 'fro');
    
    % 이력 저장
    K_hist(iter, :)  = K_iter;
    Minv_hist{iter}  = Minv_iter;
    conv_hist(iter)   = rel_change;
    D_hist{iter}      = D_iter;
    
    if rel_change < iter_tol
        fprintf('  수렴: iter %d, ‖ΔD''‖/‖D''‖ = %.2e\n', iter, rel_change);
        break;
    end
end
n_iter_done = min(iter, iter_max);

if rel_change >= iter_tol
    fprintf('  [경고] %d회 반복 후 미수렴. 최종 ‖ΔD''‖/‖D''‖ = %.2e\n', ...
        iter_max, rel_change);
end

% --- 최종 결과 추출 ---
K_est_iter    = K_hist(n_iter_done, :);
Minv_iter_fin = Minv_hist{n_iter_done};
D_iter_fin    = D_hist{n_iter_done};

% 최종 반공진 (검증용)
omega_AR_iter_fin = zeros(DOF, DOF);
for i = 1:DOF
    N_i_fin = D_iter_fin - omega_m_sq(i) * (e_mat(:,i) * e_mat(:,i)');
    ar_sq = sort(real(eig(N_i_fin)));
    ar_sq = max(ar_sq, 0);
    omega_AR_iter_fin(i, :) = sqrt(ar_sq);
end

% --- 결과 출력 ---
fprintf('\n--- Stiffness (Iterative, %d iter) ---\n', n_iter_done);
for i = 1:DOF
    err_pct = abs(K_est_iter(i) - K_true(i)) / K_true(i) * 100;
    fprintf('  Joint %d: True = %.2f,  Est = %.2f,  Err = %.2f%%\n', ...
        selected_joints(i), K_true(i), K_est_iter(i), err_pct);
end

fprintf('\n[대각항 (M^{-1})_{ii}]\n');
for i = 1:DOF
    true_v = M_inv_gt(i,i);
    est_v  = Minv_iter_fin(i,i);
    err    = abs(est_v - true_v) / abs(true_v) * 100;
    fprintf('  (%d,%d): True=%+.6f  Est=%+.6f  Err=%.2f%%\n', ...
        i, i, true_v, est_v, err);
end

fprintf('[비대각항 (M^{-1})_{ij} — 부호 포함]\n');
off_pairs = [1,2; 1,3; 2,3];
for r = 1:3
    ii = off_pairs(r,1); jj = off_pairs(r,2);
    true_v = M_inv_gt(ii,jj);
    est_v  = Minv_iter_fin(ii,jj);
    err    = abs(est_v - true_v) / abs(true_v) * 100;
    fprintf('  (%d,%d): True=%+.6f  Est=%+.6f  Err=%.2f%%\n', ...
        ii, jj, true_v, est_v, err);
end

fprintf('\n[복원된 M^{-1} (Iterative)]\n');
disp(Minv_iter_fin);
fprintf('[Ground Truth M^{-1}]\n');
disp(M_inv_gt);

fprintf('[최종 D'' (Iterative)]\n');
disp(D_iter_fin);

% 반공진 비교
fprintf('[반공진 비교 (Hz)]\n');
fprintf('%-8s %12s %12s %12s %12s\n', '', 'roots', 'CMIF', 'MP(init)', 'Iterative');
for i = 1:DOF
    for r = 1:DOF
        fprintf('  AR_%d_%d  %12.4f %12.4f %12.4f %12.4f\n', i, r, ...
            omega_AR_all(i,r)/(2*pi), ...
            omega_AR_cmif(i,r)/(2*pi), ...
            omega_AR_mp(i,r)/(2*pi), ...
            omega_AR_iter_fin(i,r)/(2*pi));
    end
end

% --- 수렴 이력 플롯 ---
figure('Name', 'Iterative D'' Refinement — Convergence');

subplot(2,1,1);
semilogy(1:n_iter_done, conv_hist(1:n_iter_done), 'b-o', 'LineWidth', 1.5, 'MarkerSize', 4);
grid on; xlabel('Iteration'); ylabel('‖ΔD''‖ / ‖D''‖');
title('D'' Relative Change');
yline(iter_tol, 'r--', sprintf('tol = %.0e', iter_tol));

subplot(2,1,2);
hold on; grid on;
colors = lines(DOF);
for i = 1:DOF
    K_err_hist = abs(K_hist(1:n_iter_done, i) - K_true(i)) / K_true(i) * 100;
    plot(1:n_iter_done, K_err_hist, '-o', 'Color', colors(i,:), ...
        'LineWidth', 1.5, 'MarkerSize', 4, ...
        'DisplayName', sprintf('K_%d', selected_joints(i)));
end
xlabel('Iteration'); ylabel('K Error (%)');
title('Stiffness Estimation Error');
legend('Location', 'best');

% --- 4가지 방법 최종 비교 ---
fprintf('\n===== 방법 1 vs 방법 2 vs 방법 3 vs 방법 4 비교 =====\n');
fprintf('%-12s %10s %10s %10s %10s %10s\n', '', 'roots', 'CMIF', 'MP', 'Iterative', 'True');
for i = 1:DOF
    fprintf('K_%d        %10.3f %10.3f %10.3f %10.3f %10.3f\n', ...
        selected_joints(i), K_est_roots(i), K_est_cmif(i), K_est_mp(i), K_est_iter(i), K_true(i));
end
fprintf('\n(M^{-1})_{ii}:\n');
for i = 1:DOF
    fprintf('  (%d,%d)  roots=%+.6f  CMIF=%+.6f  MP=%+.6f  Iter=%+.6f  True=%+.6f\n', ...
        i, i, Minv_diag_est_r(i), Minv_diag_est_c(i), ...
        Minv_diag_mp(i), Minv_iter_fin(i,i), M_inv_gt(i,i));
end
fprintf('(M^{-1})_{ij} (부호 포함):\n');
for r = 1:3
    ii = off_pairs(r,1); jj = off_pairs(r,2);
    fprintf('  (%d,%d)  MP=%+.6f  Iter=%+.6f  True=%+.6f\n', ...
        ii, jj, off_mp{r,3}, Minv_iter_fin(ii,jj), M_inv_gt(ii,jj));
end

fig_robot.HandleVisibility = 'off'; 
setFigurePositions(4,500,500)
%% =======================================================================
%  로컬 함수
% =======================================================================

function omega_peak = parabolicPeak(sig, idx, omega_vec)
% 포물선 보간으로 peak 위치 정밀화
    if idx == 1 || idx == length(sig)
        omega_peak = omega_vec(idx);
        return;
    end
    sm1 = sig(idx-1);  s0 = sig(idx);  sp1 = sig(idx+1);
    delta = (sm1 - sp1) / (2*(sm1 - 2*s0 + sp1));
    domega = omega_vec(2) - omega_vec(1);
    omega_peak = omega_vec(idx) + delta * domega;
end

function [diag_est, m12, m13, m23] = recoverMinv(A_vals, B_vals, omega_m_sq, k_loc, DOF)
    W_est    = zeros(1, DOF);
    W_est(1) = (A_vals(2) + A_vals(3) - A_vals(1)) / 2;
    W_est(2) = (A_vals(1) + A_vals(3) - A_vals(2)) / 2;
    W_est(3) = (A_vals(1) + A_vals(2) - A_vals(3)) / 2;
    diag_est = (W_est - omega_m_sq) ./ k_loc;
    m23 = sqrt(max(0, (W_est(2)*W_est(3) - B_vals(1)) / (k_loc(2)*k_loc(3))));
    m13 = sqrt(max(0, (W_est(1)*W_est(3) - B_vals(2)) / (k_loc(1)*k_loc(3))));
    m12 = sqrt(max(0, (W_est(1)*W_est(2) - B_vals(3)) / (k_loc(1)*k_loc(2))));
end

function printMinvResult(diag_est, m12, m13, m23, M_inv_gt, DOF)
    fprintf('[대각항]\n');
    for i = 1:DOF
        err = abs(diag_est(i) - M_inv_gt(i,i)) / abs(M_inv_gt(i,i)) * 100;
        fprintf('  (%d,%d): True=%+.6f  Est=%+.6f  Err=%.2f%%\n', ...
            i, i, M_inv_gt(i,i), diag_est(i), err);
    end
    fprintf('[비대각항]\n');
    pairs = {1,2,m12; 1,3,m13; 2,3,m23};
    for r = 1:3
        ii=pairs{r,1}; jj=pairs{r,2}; est=pairs{r,3};
        true_v = abs(M_inv_gt(ii,jj));
        err = abs(est - true_v) / true_v * 100;
        fprintf('  |(%d,%d)|: True=%.6f  Est=%.6f  Err=%.2f%%\n', ...
            ii, jj, true_v, est, err);
    end

    fprintf('\n[복원된 M^{-1} 행렬 (부호 미결정)]\n');
    Minv_est_mat = diag(diag_est) + ...
        [0,    m12,  m13;
         m12,  0,    m23;
         m13,  m23,  0   ];
    disp(Minv_est_mat);
    
    fprintf('[Ground Truth M^{-1}]\n');
    disp(M_inv_gt);
end

function setFigurePositions(cols, width, height)
    % Set Positions for All Figures
    % cols : 숫자, 열의 수 / width : 너비 / height : 높이
    
    % Get all figure handles
    fig_handles = findall(groot, 'Type', 'figure', 'HandleVisibility', 'on');

    % Retrieve the figure numbers (internal identifiers)
    fig_numbers = arrayfun(@(h) h.Number, fig_handles);
    
    % Sort the handles by their figure numbers
    [~, idx] = sort(fig_numbers);  % Sort by Figure number
    fig_handles = fig_handles(idx);
    
    % Number of figures
    num_figures = length(fig_handles);
    
    % Positioning parameters
%     width = 1000;
%     height = 500;
    h_margin = 10;  % Increased horizontal margin
    v_margin = 100;   % Vertical margin remains the same

    % Compute the number of rows based on the number of figures and columns
    rows = ceil(num_figures / cols);

    % Compute positions for figures
    positions = zeros(num_figures, 4);

    for i = 1:num_figures
        row = floor((i-1) / cols);
        col = mod(i-1, cols);
        positions(i, :) = [col*(width + h_margin), (rows-row-1)*(height + v_margin) + 50, width, height];
    end

    % Apply positions to figures
    for i = 1:num_figures
        set(fig_handles(i), 'Position', positions(i, :));
    end
end

function params = getRobotParams_MCL(MANIPULATOR)
    % Fixed parameters
    if strcmp(MANIPULATOR, 'RA')
        params.Nm = [100 100 50 50 50 50 50];
        params.Jm = [0.7, 1.3, 0.66, 0.31, 0.46, 0.339, 0.411]*1e-4;
        params.K  = [4500, 8500, 5000, 4000, 2400, 3700, 3000];

        % Jl matrix
        % RA
        params.Jl.Config1 ...
               = [  3.1422    1.4226    0.2588
                    1.4226    0.7503    0.1541
                    0.2588    0.1541    0.0530]; 
        
        params.Jl.Config2 ...
               = [  2.7482    1.2256    0.2280
                    1.2256    0.7503    0.1541
                    0.2280    0.1541    0.0530]; %cond = 327
                    
        params.Jl.Config3 ...
               = [  1.7974    0.7502    0.1539
                    0.7502    0.7503    0.1541
                    0.1539    0.1541    0.0530]; %cond = 233.7
        
        params.Jl.Config4 ...
               = [  1.3859    0.4432   -0.0518
                    0.4432    0.5478    0.0528
                   -0.0518    0.0528    0.0530]; %cond = 45
    
        params.Jl.Config5 ...
               = [  1.8052    0.6532    0.1578
                    0.6532    0.5485    0.0532
                    0.1578    0.0532    0.0530]; %cond = 45

    else
        params.Nm = [100 100 100 50 50 50 50];
        params.Jm = [0.7, 0.67929, 0.66, 0.26386, 0.46, 0.30574, 0.411]*1e-4; 
        params.K  = [4500, 2500, 5000, 1800, 2400, 800, 3000];

        params.Jl.Config1 ...
               = [  2.48304746205680	1.12445830906160	0.202024356828311
                    1.12445830906160	0.592581357954115	0.120016294424166
                    0.202024356828311	0.120016294424166	0.0402896408754157]; 
        
        params.Jl.Config2 ...
               = [  1.09783323788665	0.352283291803278	-0.0415596731532293
                    0.352283291803278	0.433445547607615	0.0404483892509157
                   -0.0415596731532293	0.0404483892509157	0.0402896408754157];
                    
        params.Jl.Config3 ...
               = [  1.42131175508408	0.593590455575245	0.120179585445488
                    0.593590455575245	0.592581357954115	0.120016294424166
                    0.120179585445488	0.120016294424166	0.0402896408754157]; 
        
        params.Jl.Config4 ...
               = [  2.17290823715324	0.969388696609827	0.178120215650590
                    0.969388696609827	0.592581357954115	0.120016294424166
                    0.178120215650590	0.120016294424166	0.0402896408754157];
    
        params.Jl.Config5 ...
               = [  1.42523049400122	0.515664423109568	0.122138954904061
                    0.515664423109568	0.432810554105615	0.0401308924999157
                    0.122138954904061	0.0401308924999157	0.0402896408754157];

        params.Jl.Config6 ...
               = [  2.43546604575354	1.09006558024047	0.178637551108934
                    1.09006558024047	0.571377316615099	0.109414273754658
                    0.178637551108934	0.109414273754658	0.0402896408754157]; % [0, 0.0443, 0.5236]

        params.Jl.Config7 ...
               = [  2.02718909009732	0.792313187999455	-0.0208795948093955
                    0.792313187999455	0.384149487789295	0.0158003593417559
                    -0.0208795948093955	0.0158003593417559	0.0402896408754157]; % 0, 0.1515, 1.885

        params.Jl.Config8 ...
               = [  2.20346323496412	0.919361025221459	0.100832686106499
                    0.919361025221459	0.461971017366503	0.0547111241303598
                    0.100832686106499	0.0547111241303598	0.0402896408754157]; % 0, -0.4154, 1.3909

        params.Jl.Config10 ...
               = [  2.17378784927603	0.946364750690084	0.178560604847175
                    0.946364750690084	0.545653853991843	0.0965525424430301
                    0.178560604847175	0.0965525424430301	0.0402896408754157]; % [0, 0.7854,-0.7854]
    end
end