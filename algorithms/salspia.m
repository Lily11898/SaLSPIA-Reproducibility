function [P, info] = salspia(A, Q, P0, tol, maxiter, params)
% SALSPIA  Spectrally-adaptive LSPIA (SaLSPIA), paper-faithful version.
%   [P, info] = salspia(A, Q, P0, tol, maxiter, params)
%
%   严格按论文公式 (1)-(12) 实现:
%     - 两个局部权重 alpha^(1)=rho/b (long BB), alpha^(2)=b/c (short BB)  [式6]
%     - mu_k = b_k/(C_H*rho_{k-1}),  m_k = max(1-mu_hat, delta)           [式7-8]
%     - trial weight = 二次方程 (1-m)c a^2+(2m-1)b a - m rho=0 的正根      [式9-11]
%     - GLL 非单调线搜索接受/折半                                          [式12]
%
%   Input:
%     A,Q,P0,tol,maxiter  同 LSPIA
%     params.c (=sigma, default 1e-4), params.M (default 10),
%     params.delta (default 1e-8), params.eps_saf (default 1e-30)
%   Output:
%     info: iter_count, err_history, Ek_history, cpu_time, halvings, CH

if nargin < 4 || isempty(tol),     tol = 1e-6;     end
if nargin < 5 || isempty(maxiter), maxiter = 10000; end
if nargin < 6, params = struct(); end

sigma = gfd(params, 'c',     1e-4);   % GLL sufficient-decrease (论文 sigma)
M     = gfd(params, 'M',     10);     % GLL memory length
delta = gfd(params, 'delta', 1e-8);   % 截断阈值 delta (论文式8)
eps_saf = gfd(params, 'eps_saf', 1e-30); % safeguard tolerance (Algorithm 1)

%% CH = max_i sum_j B_i(t_j)   (式1),  alpha0 = 1/CH
CH = max(sum(A, 1));
alpha0 = 1.0 / CH;

%% 初始化:  r = Q - A P,  d = A^T r   (论文 D_i)
P = P0;
r = Q - A * P;
d = A' * r;
rho0 = norm(d, 'fro')^2;
g0_norm2 = rho0;                       % ||A^T(Q-AP0)||^2  用于 Ek

err_hist  = nan(maxiter+2, 1);
Ek_hist   = nan(maxiter+2, 1);
halv_hist = zeros(maxiter+2, 1);
err_hist(1) = norm(r, 'fro')^2;
Ek_hist(1)  = 1.0;
n_updates = 0;

%% 初始 safeguard: rho_0 足够小时，P0 已满足停止条件
if rho0 <= eps_saf
    info.iter_count  = 0;
    info.err_history = err_hist(1);
    info.Ek_history  = Ek_hist(1);
    info.cpu_time    = 0;
    info.halvings    = halv_hist(1);
    info.CH          = CH;
    info.eps_saf     = eps_saf;
    return;
end

%% 第一步固定用 alpha0:  P^(1)=P^(0)+alpha0 d^(0)   (论文式)
d_prev     = d;
alpha_prev = alpha0;
P = P + alpha0 * d;
r = Q - A * P;
d = A' * r;
n_updates = n_updates + 1;
err_hist(n_updates+1) = norm(r,'fro')^2;
Ek_hist(n_updates+1)  = norm(d,'fro')^2 / g0_norm2;

%% GLL 拟合误差队列 Phi = 0.5||r||^2
Phi_queue = [0.5*norm(Q - A*P0,'fro')^2; 0.5*norm(r,'fro')^2];
Phi_queue = Phi_queue(max(1, end-M+1):end);

tic;
while n_updates < maxiter
    %% rho_{k-1}, rho_k, zeta_k   (式3,4)
    rho_km1 = norm(d_prev,'fro')^2;
    rho_k   = norm(d,'fro')^2;

    % Algorithm 1: while rho_k > eps_saf and E_k > tol
    if rho_k <= eps_saf || Ek_hist(n_updates+1) < tol
        break;
    end

    zeta_k  = sum(sum(d_prev .* d));

    %% b_k, c_k   (式5)
    b_k = (rho_km1 - zeta_k) / alpha_prev;
    c_k = (rho_km1 - 2*zeta_k + rho_k) / (alpha_prev^2);

    %% trial weight  (式6-11)
    if b_k <= eps_saf || c_k <= eps_saf
        alpha_bar = alpha0;                 % 保护分支
    else
        a1 = rho_km1 / b_k;                 % alpha^(1) long BB
        a2 = b_k / c_k;                     % alpha^(2) short BB

        mu_k   = b_k / (CH * rho_km1);      % (式7)
        mu_hat = min(max(mu_k, 0), 1);
        m_k    = max(1 - mu_hat, delta);    % (式8)

        if (1 - mu_hat) <= delta
            alpha_bar = a2;                 % (式11)
        else
            % 稳定计算式(2.10)的正根，避免普通二次公式的消去误差
            Aq = (1 - m_k) * c_k;
            Bq = (2*m_k - 1) * b_k;
            Dq = Bq^2 + 4*Aq*m_k*rho_km1;
            sqrt_Dq = sqrt(max(Dq, 0));

            if Aq <= eps_saf
                alpha_tilde = a1;           % 极限值 rho_{k-1}/b_k
            elseif Bq >= 0
                alpha_tilde = 2*m_k*rho_km1 / (Bq + sqrt_Dq);
            else
                alpha_tilde = (sqrt_Dq - Bq) / (2*Aq);
            end
            alpha_bar = alpha_tilde;
        end
    end

    %% GLL 非单调线搜索  (式12):  0.5||Q-A(P+alpha d)||^2 <= Fk - sigma*alpha*rho_k
    Fk      = max(Phi_queue);
    rho_ls  = norm(d,'fro')^2;
    alpha   = alpha_bar;
    n_halv  = 0;
    while true
        r_try = Q - A * (P + alpha * d);
        if 0.5*norm(r_try,'fro')^2 <= Fk - sigma*alpha*rho_ls
            break;
        end
        alpha  = alpha / 2;
        n_halv = n_halv + 1;
        if n_halv > 100, break; end
    end
    halv_hist(n_updates+1) = n_halv;

    %% 更新:  P^{k+1} = P^k + alpha d^k
    alpha_prev = alpha;
    d_prev     = d;
    P = P + alpha * d;
    r = Q - A * P;
    d = A' * r;

    Phi = 0.5*norm(r,'fro')^2;
    if numel(Phi_queue) < M
        Phi_queue = [Phi_queue; Phi];
    else
        Phi_queue = [Phi_queue(2:end); Phi];
    end

    n_updates = n_updates + 1;
    err_hist(n_updates+1) = norm(r,'fro')^2;
    Ek_hist(n_updates+1)  = norm(d,'fro')^2 / g0_norm2;

end
cpu = toc;

%% 输出
info.iter_count  = n_updates;
info.err_history = err_hist(1:n_updates+1);
info.Ek_history  = Ek_hist(1:n_updates+1);
info.cpu_time    = cpu;
info.halvings    = halv_hist(1:max(1,n_updates));
info.CH          = CH;
info.eps_saf     = eps_saf;
end

function v = gfd(s, f, d)
    if isfield(s,f), v = s.(f); else, v = d; end
end
