function [P, info] = alspia(A, Q, P0, tol, maxiter)
% ALSPIA  Asynchronous LSPIA with Chebyshev semi-iterative step sizes.
%   [P, info] = alspia(A, Q, P0, tol, maxiter)
%
%   Reference: Wu & Liu, "Asynchronous progressive iterative approximation
%   method for least squares fitting", CAGD 2024.
%
%   Update rule:  P^{k+1} = P^k + omega_k * A^T * (Q - A * P^k)
%   Step sizes from Chebyshev roots (Theorem 3.1, Eq. 5):
%     omega_{ell} = 2 / ((v+u) + (v-u)*cos((2*ell+1)/(2*K)*pi))
%   Deterministic strategy: ell_k = mod(k-1, K), K = 10^6.
%
%   Input:
%     A       - (m+1)-by-(n+1) collocation matrix
%     Q       - (m+1)-by-d data point matrix
%     P0      - (n+1)-by-d initial control points
%     tol     - tolerance for relative gradient norm E_k (default 1e-6)
%     maxiter - maximum iterations (default 10000)
%
%   Output:
%     P       - final control points
%     info    - struct with iter_count, err_history, Ek_history, cpu_time

if nargin < 4, tol = 1e-6;   end
if nargin < 5, maxiter = 10000; end

%% Precomputation
eigvals = eig(A' * A);
eigvals_pos = eigvals(eigvals > 1e-14);
lam_max = max(eigvals_pos);

rank_tol = 1e-10 * lam_max;
is_full_rank = (min(eigvals) > rank_tol);

%% Precompute Chebyshev step sizes
K = min(1e6, maxiter);

if is_full_rank
    % Theorem 3.1, Eq. (5)
    u = min(eigvals_pos);
    v = lam_max;
    ell = (0:K-1)';
    omega_all = 2.0 ./ ((v + u) + (v - u) * cos((2*ell + 1) / (2*K) * pi));
else
    % Theorem 3.2, Eq. (6)
    v = lam_max;
    r_Kp1 = cos((2*K + 1) / (2*(K + 1)) * pi);
    ell = (0:K-1)';
    omega_all = (1 - r_Kp1) ./ ...
        (v * (cos((2*ell + 1) / (2*(K+1)) * pi) - r_Kp1));
end

%% Initialize
P = P0;
r = A * P - Q;
g = A' * r;
g0_norm2 = norm(g, 'fro')^2;

err_hist = zeros(maxiter+1, 1);
Ek_hist  = zeros(maxiter+1, 1);
err_hist(1) = norm(r, 'fro')^2;
Ek_hist(1)  = 1.0;

%% Iteration
tic;
for k = 1:maxiter
    wk = omega_all(mod(k-1, K) + 1);   % deterministic Chebyshev index

    P = P - wk * g;      % update
    r = A * P - Q;        % residual
    g = A' * r;           % gradient

    err_hist(k+1) = norm(r, 'fro')^2;
    Ek_hist(k+1)  = norm(g, 'fro')^2 / g0_norm2;

    if Ek_hist(k+1) < tol
        err_hist = err_hist(1:k+1);
        Ek_hist  = Ek_hist(1:k+1);
        break;
    end
end
cpu = toc;

%% Output
info.iter_count   = min(k, maxiter);
info.err_history  = err_hist;
info.Ek_history   = Ek_hist;
info.cpu_time     = cpu;
info.is_full_rank = is_full_rank;
info.K_chebyshev  = K;
end
