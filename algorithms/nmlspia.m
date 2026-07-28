function [P, info] = nmlspia(A, Q, P0, tol, maxiter)
% NMLSPIA  LSPIA accelerated by Nesterov's momentum.
%   [P, info] = nmlspia(A, Q, P0, tol, maxiter)
%
%   Reference: Liu, Wu, Li & Hu, "Two novel iterative approaches for
%   improved LSPIA convergence", CAGD 2024. NmLSPIA is Eq. (6)/(10);
%   optimal step sizes are given by Theorem 2. Among the two variants of
%   that paper, NmLSPIA is chosen because Corollary 1 proves its asymptotic
%   rate is faster than PmLSPIA.
%
%   Matrix form (Eq. 10), with mbar = P^k - P^{k-1} the previous momentum:
%     Plook   = P^k + eta * mbar
%     P^{k+1} = P^k + zeta * A^T( Q - A*Plook ) + eta * mbar
%
%   Optimal parameters (Theorem 2), requiring the extreme eigenvalues of
%   A^T*A:
%     zeta = 1 / lambda_max
%     eta  = (sqrt(lambda_max) - sqrt(lambda_min)) ...
%            / (sqrt(lambda_max) + sqrt(lambda_min))
%
%   Input / output conventions are identical to lspia.m and alspia.m so
%   that this function is a drop-in comparison method.
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
%     info    - struct with fields:
%                 iter_count, err_history, Ek_history, cpu_time,
%                 zeta, eta, lam_min, lam_max

if nargin < 4, tol = 1e-6;   end
if nargin < 5, maxiter = 10000; end

%% Precomputation (extreme eigenvalues of A^T*A)
eigvals = eig(A' * A);
eigvals = eigvals(eigvals > 1e-14);
lam_min = min(eigvals);
lam_max = max(eigvals);

% Theorem 2 optimal parameters.
% NOTE: in Liu (2024) the eigenvalues are sorted in non-increasing order,
% so lambda_0 = lam_max. Hence zeta = 1/lam_max (NOT 1/lam_min).
zeta = 1.0 / lam_max;
eta  = (sqrt(lam_max) - sqrt(lam_min)) / (sqrt(lam_max) + sqrt(lam_min));

%% Initialize
P      = P0;            % P^{k}
P_prev = P0;            % P^{k-1};  P^{1} = P^{0} so first momentum is zero

r = A * P - Q;          % residual at P^0
g = A' * r;             % gradient at P^0
g0_norm2 = norm(g, 'fro')^2;

err_hist = zeros(maxiter+1, 1);
Ek_hist  = zeros(maxiter+1, 1);
err_hist(1) = norm(r, 'fro')^2;
Ek_hist(1)  = 1.0;

%% Iteration
tic;
for k = 1:maxiter
    mbar  = P - P_prev;                 % previous momentum  P^k - P^{k-1}
    Plook = P + eta * mbar;             % Nesterov look-ahead point

    g_look = A' * (A * Plook - Q);      % gradient at the look-ahead point
    P_new  = P - zeta * g_look + eta * mbar;   % Eq. (10) update

    % advance
    P_prev = P;
    P      = P_new;

    % termination test uses the TRUE gradient at P^{k}, matching lspia.m
    r = A * P - Q;
    g = A' * r;

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
info.iter_count  = min(k, maxiter);
info.err_history = err_hist;
info.Ek_history  = Ek_hist;
info.cpu_time    = cpu;
info.zeta        = zeta;
info.eta         = eta;
info.lam_min     = lam_min;
info.lam_max     = lam_max;
end