function [P, info] = mlspia(A, Q, P0, tol, maxiter)
% MLSPIA  Progressive Iterative Approximation with Memory for Least Squares.
%   [P, info] = mlspia(A, Q, P0, tol, maxiter)
%
%   Reference: Huang & Wang, "On a progressive and iterative approximation
%   method with memory for least square fitting", CAGD 2020.
%
%   Iteration (Eq. 8 in paper), at step k:
%     (a) Lambda^{k+1} = (1-omega)*Lambda^k  -  gamma*nu * A*(A^T*Lambda^k)
%                        + omega*(Q - A*P^k)
%     (b) P^{k+1}      = P^k  +  nu * A^T * Lambda^k
%
%   Optimal weights (Theorem 6, Eq. 27):
%     omega* = gamma* = 4*s1*sr / (s1+sr)^2,   nu* = 1/(s1*sr)
%   where s1, sr are the largest/smallest positive singular values of A.
%
%   Initial Lambda^0  (initial control points II in paper, Section 6.2):
%     Lambda^0 = omega * (Q - A*P^0)
%
%   Stopping: relative gradient norm
%     E_k = ||A^T(A*P^k - Q)||_F^2 / ||A^T(A*P^0 - Q)||_F^2 < tol
%
%   Input:
%     A       - (m+1) x (n+1) collocation matrix
%     Q       - (m+1) x d  data points
%     P0      - (n+1) x d  initial control points
%     tol     - tolerance for E_k  (default 1e-6)
%     maxiter - max iterations      (default 10000)
%
%   Output:
%     P       - converged control points
%     info    - struct: iter_count, err_history, Ek_history, cpu_time,
%               omega, gamma_w, nu, lam_min, lam_max, s1, sr
if nargin < 4, tol     = 1e-6;   end
if nargin < 5, maxiter = 10000;  end
%% ---- Optimal weights ---------------------------------------------------
eigvals     = eig(A' * A);
eigvals_pos = sort(eigvals(eigvals > 1e-14));   % ascending
s1 = sqrt(eigvals_pos(end));   % largest  singular value of A
sr = sqrt(eigvals_pos(1));     % smallest positive singular value of A
omega = 4 * s1 * sr / (s1 + sr)^2;   % = gamma*
gamma = omega;
nu    = 1.0 / (s1 * sr);
lam_max = eigvals_pos(end);
lam_min = eigvals_pos(1);
%% ---- Initialize --------------------------------------------------------
P      = P0;
r      = A * P - Q;           % r^0 = A*P^0 - Q,  size (m+1) x d
Lambda = omega * (-r);        % Lambda^0 = omega*(Q - A*P^0)
g      = A' * r;              % gradient at P^0
g0_norm2 = norm(g, 'fro')^2;
err_hist = zeros(maxiter+1, 1);
Ek_hist  = zeros(maxiter+1, 1);
err_hist(1) = norm(r, 'fro')^2;
Ek_hist(1)  = 1.0;
%% ---- Main iteration ----------------------------------------------------
% At the start of iteration k we hold:  P = P^k,  r = r^k,  Lambda = Lambda^k
tic;
for k = 1:maxiter
    %% Shared computation reused by both updates
    AtL  = A' * Lambda;       % (n+1) x d  :  A^T * Lambda^k
    AAtL = A  * AtL;          % (m+1) x d  :  A*A^T * Lambda^k
    %% (a) Update Lambda (uses current P^k via r = A*P^k - Q)
    Lambda = (1 - omega) * Lambda  -  gamma * nu * AAtL  -  omega * r;
    %% (b) Update P
    P = P + nu * AtL;
    %% Update residual and gradient:
    %%   r^{k+1} = A*P^{k+1} - Q = (A*P^k - Q) + nu*A*(A^T*Lambda^k)
    r = r + nu * AAtL;
    g = A' * r;
    %% Record
    err_hist(k+1) = norm(r, 'fro')^2;
    Ek_hist(k+1)  = norm(g, 'fro')^2 / g0_norm2;
    %% Convergence check
if Ek_hist(k+1) < tol
        err_hist = err_hist(1:k+1);
        Ek_hist  = Ek_hist(1:k+1);
break;
end
end
cpu = toc;
%% ---- Output ------------------------------------------------------------
info.iter_count  = min(k, maxiter);
info.err_history = err_hist;
info.Ek_history  = Ek_hist;
info.cpu_time    = cpu;
info.omega       = omega;
info.gamma_w     = gamma;
info.nu          = nu;
info.lam_min     = lam_min;
info.lam_max     = lam_max;
info.s1          = s1;
info.sr          = sr;
end