function [P, info] = ablation_lspia(A, Q, P0, tol, maxiter, params)
% ABLATION_LSPIA  SaLSPIA ablation variants for curve fitting.
%
%   params.mode:
%       'bb1'         BB1 trial weight + GLL halving
%       'bb2'         BB2 trial weight + GLL halving
%       'interp-only' SaLSPIA interpolated trial weight, no line search
%       'salspia'     SaLSPIA interpolated trial weight + GLL halving
%
%   Residual and gradient are stored as r = A*P - Q and g = A'*r.
%   Hence the update is P_{k+1} = P_k - alpha_k*g_k.

if nargin < 4 || isempty(tol),     tol = 1e-6; end
if nargin < 5 || isempty(maxiter), maxiter = 10000; end
if nargin < 6 || isempty(params),  params = struct(); end

mode = lower(strrep(gfd(params, 'mode', 'salspia'), '_', '-'));
sigma       = gfd(params, 'c', 1e-4);
M           = gfd(params, 'M', 10);
delta       = gfd(params, 'delta', gfd(params, 'eps_m', 1e-8));
eps_saf     = gfd(params, 'eps_saf', 1e-30);
p_thr       = gfd(params, 'p', 3);
e_tol       = gfd(params, 'e_tol', 1e-15);
max_halving = gfd(params, 'max_halving', 100);
alpha_scale = gfd(params, 'alpha_scale', 1.0);

use_line_search = ~strcmp(mode, 'interp-only');

CH = max(sum(A, 1));
if CH <= 0
    error('C_H must be positive.');
end

P = P0;
r = A * P - Q;
g = A' * r;
f_val = 0.5 * norm(r, 'fro')^2;
e_val = 2 * f_val;
g_norm2 = norm(g, 'fro')^2;
g0_norm2 = g_norm2;

F_queue = f_val;
nc = 0;
n_updates = 0;
status = 'maxiter';

err_hist   = nan(maxiter + 1, 1);
Ek_hist    = nan(maxiter + 1, 1);
alpha_hist = nan(maxiter, 1);
trial_hist = nan(maxiter, 1);
bb1_hist   = nan(maxiter, 1);
bb2_hist   = nan(maxiter, 1);
mu_hist    = nan(maxiter, 1);
m_hist     = nan(maxiter, 1);
halv_hist  = zeros(maxiter, 1);

err_hist(1) = e_val;
Ek_hist(1) = initial_Ek(g0_norm2);

P_prev = P;
g_prev = g;
e_prev = e_val;

tic;
for k = 0:maxiter-1
    if k >= 1
        if ~isfinite(e_val) || ~isfinite(g_norm2)
            status = 'failed';
            break;
        end

        if abs(e_val - e_prev) < e_tol
            nc = nc + 1;
        else
            nc = 0;
        end

        Ek_cur = safe_ratio(g_norm2, g0_norm2);
        if g0_norm2 > 0 && Ek_cur < tol
            status = 'converged';
            break;
        end
        if nc >= p_thr
            status = 'stagnated';
            break;
        end
    end

    if k == 0
        alpha_hat = 1.0 / CH;
        step = empty_step();
    else
        s = P - P_prev;
        y = g - g_prev;
        [alpha_hat, step] = choose_trial_step(s, y, CH, delta, eps_saf, mode);
    end

    if ~isfinite(alpha_hat) || alpha_hat <= 0
        alpha_hat = 1.0 / CH;
        step.used_fallback = true;
    end
    alpha_hat = alpha_scale * alpha_hat;

    Ag = A * g;
    alpha = alpha_hat;
    n_halv = 0;

    if use_line_search
        Fk = max(F_queue);
        while true
            r_trial = r - alpha * Ag;
            f_trial = 0.5 * norm(r_trial, 'fro')^2;
            if isfinite(f_trial) && f_trial <= Fk - sigma * alpha * g_norm2
                break;
            end

            alpha = alpha / 2;
            n_halv = n_halv + 1;
            if n_halv >= max_halving
                status = 'line-search-limit';
                break;
            end
        end
    end

    P_prev = P;
    g_prev = g;
    e_prev = e_val;

    P = P - alpha * g;
    r = r - alpha * Ag;
    g = A' * r;

    f_val = 0.5 * norm(r, 'fro')^2;
    e_val = 2 * f_val;
    g_norm2 = norm(g, 'fro')^2;

    n_updates = n_updates + 1;
    err_hist(n_updates + 1) = e_val;
    Ek_hist(n_updates + 1) = safe_ratio(g_norm2, g0_norm2);
    alpha_hist(n_updates) = alpha;
    trial_hist(n_updates) = alpha_hat;
    bb1_hist(n_updates) = step.alpha_bb1;
    bb2_hist(n_updates) = step.alpha_bb2;
    mu_hist(n_updates) = step.mu_hat;
    m_hist(n_updates) = step.m;
    halv_hist(n_updates) = n_halv;

    if length(F_queue) < M
        F_queue = [F_queue; f_val];
    else
        F_queue = [F_queue(2:end); f_val];
    end

    if strcmp(status, 'line-search-limit')
        break;
    end
end
cpu = toc;

if strcmp(status, 'maxiter') && n_updates < maxiter
    status = 'stopped';
end

info.mode          = mode;
info.iter_count    = n_updates;
info.err_history   = err_hist(1:n_updates + 1);
info.Ek_history    = Ek_hist(1:n_updates + 1);
info.alpha_history = alpha_hist(1:n_updates);
info.trial_alpha   = trial_hist(1:n_updates);
info.alpha_BB1     = bb1_hist(1:n_updates);
info.alpha_BB2     = bb2_hist(1:n_updates);
info.mu_history    = mu_hist(1:n_updates);
info.m_history     = m_hist(1:n_updates);
info.halvings      = halv_hist(1:n_updates);
info.total_halving = sum(halv_hist(1:n_updates));
info.cpu_time      = cpu;
info.CH            = CH;
info.status        = status;
info.line_search   = use_line_search;
info.alpha_scale   = alpha_scale;
end

function [alpha_hat, step] = choose_trial_step(s, y, CH, delta, eps_saf, mode)
step = empty_step();

a = norm(s, 'fro')^2;
b = sum(s(:) .* y(:));
c = norm(y, 'fro')^2;

if a <= eps_saf || b <= eps_saf || c <= eps_saf
    alpha_hat = 1.0 / CH;
    step.used_fallback = true;
    return;
end

alpha_bb1 = a / b;
alpha_bb2 = b / c;
mu = b / (CH * a);
mu_hat = min(max(mu, 0), 1);
m = max(1.0 - mu_hat, delta);

switch mode
    case 'bb1'
        alpha_hat = alpha_bb1;
    case 'bb2'
        alpha_hat = alpha_bb2;
    case {'interp-only', 'salspia', 'full'}
        if 1.0 - mu_hat <= delta
            alpha_hat = alpha_bb2;
        else
            Aq = (1.0 - m) * c;
            Bq = (2.0 * m - 1.0) * b;
            Dq = Bq^2 + 4.0 * Aq * m * a;
            sqrt_Dq = sqrt(max(Dq, 0));
            if Aq <= eps_saf
                alpha_hat = alpha_bb1;
            elseif Bq >= 0
                alpha_hat = 2.0 * m * a / (Bq + sqrt_Dq);
            else
                alpha_hat = (sqrt_Dq - Bq) / (2.0 * Aq);
            end
        end
    otherwise
        error('Unknown ablation mode: %s', mode);
end

step.alpha_bb1 = alpha_bb1;
step.alpha_bb2 = alpha_bb2;
step.mu = mu;
step.mu_hat = mu_hat;
step.m = m;
end

function step = empty_step()
step.alpha_bb1 = NaN;
step.alpha_bb2 = NaN;
step.mu = NaN;
step.mu_hat = NaN;
step.m = NaN;
step.used_fallback = false;
end

function y = safe_ratio(a, b)
if b > 0
    y = a / b;
else
    y = 0;
end
end

function y = initial_Ek(g0_norm2)
if g0_norm2 > 0
    y = 1.0;
else
    y = 0.0;
end
end

function v = gfd(s, f, d)
if isfield(s, f)
    v = s.(f);
else
    v = d;
end
end
