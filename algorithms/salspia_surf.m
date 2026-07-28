function [P, info] = salspia_surf(Au, Av, Q, P0, tol, maxiter, params)
% SALSPIA_SURF  Paper-faithful SaLSPIA for tensor-product surface fitting.
%   [P, info] = salspia_surf(Au, Av, Q, P0, tol, maxiter, params)
%
%   论文公式 (式3-12) 的曲面版本: 内积/范数在所有 d 个坐标上聚合,
%   谱上界 C_H = C_H^u * C_H^v.  trial weight 用二次方程正根.
%
%   Input:
%     Au,Av  - u/v 方向配置矩阵
%     Q      - cell{d}, 每个 (mu+1)x(mv+1)
%     P0     - cell{d}, 每个 (nu+1)x(nv+1)
%     params - .c(=sigma,1e-4), .M(10), .delta(1e-8),
%              .eps_saf(1e-30)

if nargin < 5 || isempty(tol),     tol = 1e-6;     end
if nargin < 6 || isempty(maxiter), maxiter = 10000; end
if nargin < 7, params = struct(); end

sigma = gfd(params, 'c',     1e-4);
M     = gfd(params, 'M',     10);
delta = gfd(params, 'delta', 1e-8);
eps_saf = gfd(params, 'eps_saf', 1e-30);

d = numel(Q);
CH = max(sum(Au,1)) * max(sum(Av,1));   % C_H^u * C_H^v
alpha0 = 1.0 / CH;

%% 初始化:  D_c = Au^T (Q_c - Au P_c Av^T) Av
[D, g0_norm2, e_val] = compute_DG(Au, Av, P0, Q);
P = P0;

err_hist  = nan(maxiter+2, 1);
Ek_hist   = nan(maxiter+2, 1);
halv_hist = zeros(maxiter+2, 1);
err_hist(1) = e_val;
Ek_hist(1)  = 1.0;
n_updates = 0;

%% 初始 safeguard: rho_0 足够小时，P0 已满足停止条件
if g0_norm2 <= eps_saf
    info.iter_count  = 0;
    info.err_history = err_hist(1);
    info.Ek_history  = Ek_hist(1);
    info.cpu_time    = 0;
    info.halvings    = halv_hist(1);
    info.CH          = CH;
    info.eps_saf     = eps_saf;
    return;
end

%% 第一步固定 alpha0
D_prev     = D;
alpha_prev = alpha0;
for c = 1:d, P{c} = P{c} + alpha0 * D{c}; end
[D, g2, e_val] = compute_DG(Au, Av, P, Q);
n_updates = n_updates + 1;
err_hist(n_updates+1) = e_val;
Ek_hist(n_updates+1)  = g2 / g0_norm2;

Phi_queue = [0.5*compute_err(Au,Av,P0,Q); 0.5*e_val];
Phi_queue = Phi_queue(max(1, end-M+1):end);

tic;
while n_updates < maxiter
    rho_km1 = sum_fro2(D_prev);
    rho_k   = sum_fro2(D);

    % Algorithm 1: while rho_k > eps_saf and E_k > tol
    if rho_k <= eps_saf || Ek_hist(n_updates+1) < tol
        break;
    end

    zeta_k  = 0;
    for c = 1:d, zeta_k = zeta_k + sum(sum(D_prev{c} .* D{c})); end

    b_k = (rho_km1 - zeta_k) / alpha_prev;
    c_k = (rho_km1 - 2*zeta_k + rho_k) / (alpha_prev^2);

    if b_k <= eps_saf || c_k <= eps_saf
        alpha_bar = alpha0;
    else
        a1 = rho_km1 / b_k;
        a2 = b_k / c_k;
        mu_k   = b_k / (CH * rho_km1);
        mu_hat = min(max(mu_k, 0), 1);
        m_k    = max(1 - mu_hat, delta);
        if (1 - mu_hat) <= delta
            alpha_bar = a2;
        else
            % 稳定计算式(2.10)的正根，避免普通二次公式的消去误差
            Aq = (1 - m_k) * c_k;
            Bq = (2*m_k - 1) * b_k;
            Dq = Bq^2 + 4*Aq*m_k*rho_km1;
            sqrt_Dq = sqrt(max(Dq, 0));

            if Aq <= eps_saf
                alpha_tilde = a1;
            elseif Bq >= 0
                alpha_tilde = 2*m_k*rho_km1 / (Bq + sqrt_Dq);
            else
                alpha_tilde = (sqrt_Dq - Bq) / (2*Aq);
            end
            alpha_bar = alpha_tilde;
        end
    end

    %% GLL 线搜索
    Fk     = max(Phi_queue);
    rho_ls = sum_fro2(D);
    alpha  = alpha_bar;
    n_halv = 0;
    while true
        et = 0;
        for c = 1:d
            Rt = Q{c} - Au * (P{c} + alpha*D{c}) * Av';
            et = et + norm(Rt,'fro')^2;
        end
        if 0.5*et <= Fk - sigma*alpha*rho_ls
            break;
        end
        alpha  = alpha / 2;
        n_halv = n_halv + 1;
        if n_halv > 100, break; end
    end
    halv_hist(n_updates+1) = n_halv;

    alpha_prev = alpha;
    D_prev     = D;
    for c = 1:d, P{c} = P{c} + alpha * D{c}; end
    [D, g2, e_val] = compute_DG(Au, Av, P, Q);

    Phi = 0.5*e_val;
    if numel(Phi_queue) < M
        Phi_queue = [Phi_queue; Phi];
    else
        Phi_queue = [Phi_queue(2:end); Phi];
    end

    n_updates = n_updates + 1;
    err_hist(n_updates+1) = e_val;
    Ek_hist(n_updates+1)  = g2 / g0_norm2;

end
cpu = toc;

info.iter_count  = n_updates;
info.err_history = err_hist(1:n_updates+1);
info.Ek_history  = Ek_hist(1:n_updates+1);
info.cpu_time    = cpu;
info.halvings    = halv_hist(1:max(1,n_updates));
info.CH          = CH;
info.eps_saf     = eps_saf;
end

%% ---- helpers ----
function [D, g2, e] = compute_DG(Au, Av, P, Q)
    d = numel(Q); D = cell(d,1); g2 = 0; e = 0;
    for c = 1:d
        Rc = Q{c} - Au * P{c} * Av';
        D{c} = Au' * Rc * Av;
        g2 = g2 + norm(D{c},'fro')^2;
        e  = e  + norm(Rc,'fro')^2;
    end
end

function e = compute_err(Au, Av, P, Q)
    d = numel(Q); e = 0;
    for c = 1:d
        Rc = Q{c} - Au * P{c} * Av';
        e = e + norm(Rc,'fro')^2;
    end
end

function s = sum_fro2(C)
    s = 0;
    for c = 1:numel(C), s = s + norm(C{c},'fro')^2; end
end

function v = gfd(s, f, d)
    if isfield(s,f), v = s.(f); else, v = d; end
end
