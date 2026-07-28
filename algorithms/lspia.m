function [P, info] = lspia(A, Q, P0, tol, maxiter)
% LSPIA  Least Squares Progressive Iterative Approximation.
%   [P, info] = lspia(A, Q, P0, tol, maxiter)
%
%   Reference: Deng & Lin, "Progressive and iterative approximation for
%   least squares B-spline curve and surface fitting", CAD 2014.
%
%   Update rule:  P^{k+1} = P^k + omega * A^T * (Q - A * P^k)
%   Optimal step: omega = 2 / (lambda_max + lambda_min) of A^T*A
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
%                 iter_count, err_history, Ek_history, cpu_time

if nargin < 4, tol = 1e-6;   end
if nargin < 5, maxiter = 10000; end

%% Precomputation
eigvals = eig(A' * A);
eigvals = eigvals(eigvals > 1e-14);
lam_min = min(eigvals);
lam_max = max(eigvals);
omega   = 2.0 / (lam_max + lam_min);

%% Initialize
P = P0;
r = A * P - Q;          % residual
g = A' * r;             % gradient
g0_norm2 = norm(g, 'fro')^2;

err_hist = zeros(maxiter+1, 1);
Ek_hist  = zeros(maxiter+1, 1);
err_hist(1) = norm(r, 'fro')^2;
Ek_hist(1)  = 1.0;

%% Iteration
tic;
for k = 1:maxiter
    P = P - omega * g;   % update
    r = A * P - Q;       % residual
    g = A' * r;          % gradient 

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
info.omega       = omega;
info.lam_min     = lam_min;
info.lam_max     = lam_max;
end
