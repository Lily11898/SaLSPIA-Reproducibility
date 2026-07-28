function [P_cell, info] = mlspia_surf(Au, Av, Q_cell, P0_cell, tol, maxiter)
% MLSPIA_SURF  MLSPIA extended to tensor-product surface fitting.
%   [P_cell, info] = mlspia_surf(Au, Av, Q_cell, P0_cell, tol, maxiter)
%
%   Reference: Huang & Wang (2020), Section 5 (Theorem 7).
%
%   For each coordinate c, the fitted surface satisfies:
%     min_{Pc}  (1/2) * ||Au * Pc * Av^T - Qc||_F^2
%
%   Iteration (Theorem 7, Eq. 37 weights):
%     For each coordinate c:
%       R^{k}_{c}      = Au * P^k_c * Av^T - Q_c        (residual)
%       Lambda^{k+1}_c = (1-omega)*Lambda^k_c
%                        - gamma*nu * Au*(Au^T*Lambda^k_c*Av)*Av^T
%                        + omega*(-R^k_c)
%       P^{k+1}_c      = P^k_c + nu * Au^T * Lambda^k_c * Av
%
%   Optimal weights (Theorem 7, Eq. 37):
%     omega* = gamma* = 4*S1*Sr / (S1+Sr)^2
%     nu*    = 1 / (S1*Sr)
%   where S1 = s1(Au)*mu1(Av),  Sr = sr(Au)*mus(Av)
%   (products of extreme singular values of Au and Av).
%
%   Stopping: relative gradient norm,
%     E_k = sum_c ||Au^T*R^k_c*Av||_F^2  /  sum_c ||Au^T*R^0_c*Av||_F^2
%
%   Input:
%     Au, Av    - collocation matrices (m1+1)x(n1+1) and (m2+1)x(n2+1)
%     Q_cell    - cell{d} of (m1+1)x(m2+1) data coordinate matrices
%     P0_cell   - cell{d} of (n1+1)x(n2+1) initial control point matrices
%     tol       - tolerance for E_k  (default 1e-6)
%     maxiter   - max iterations      (default 10000)
%
%   Output:
%     P_cell    - cell{d} of converged control point matrices
%     info      - struct: iter_count, Ek_history, err_history, cpu_time,
%                 omega, gamma_w, nu

if nargin < 5, tol     = 1e-6;   end
if nargin < 6, maxiter = 10000;  end

d = numel(Q_cell);

%% ---- Optimal weights (Theorem 7) --------------------------------------
sv_Au = svd(Au);  sv_Au = sv_Au(sv_Au > 1e-14);
sv_Av = svd(Av);  sv_Av = sv_Av(sv_Av > 1e-14);

s1  = max(sv_Au);   sr  = min(sv_Au);   % singular values of Au
mu1 = max(sv_Av);   mus = min(sv_Av);   % singular values of Av

S1 = s1 * mu1;   % largest  effective singular value
Sr = sr * mus;   % smallest effective singular value

omega = 4 * S1 * Sr / (S1 + Sr)^2;
gamma = omega;
nu    = 1.0 / (S1 * Sr);

%% ---- Initialize --------------------------------------------------------
P_cell      = P0_cell;
R_cell      = cell(d, 1);   % residuals
Lambda_cell = cell(d, 1);   % Lambda matrices

g0_norm2 = 0;
for c = 1:d
    R_cell{c}      = Au * P_cell{c} * Av' - Q_cell{c};
    Lambda_cell{c} = omega * (-R_cell{c});
    Gc             = Au' * R_cell{c} * Av;     % gradient at P^0
    g0_norm2       = g0_norm2 + norm(Gc, 'fro')^2;
end

err_hist = zeros(maxiter+1, 1);
Ek_hist  = zeros(maxiter+1, 1);

% Initial error
e0 = 0;
for c = 1:d, e0 = e0 + norm(R_cell{c}, 'fro')^2; end
err_hist(1) = e0;
Ek_hist(1)  = 1.0;

%% ---- Main iteration ----------------------------------------------------
tic;
for k = 1:maxiter

    g_norm2 = 0;
    e_val   = 0;

    for c = 1:d
        %% Shared:  Au^T*Lambda^k_c  and  Au*(Au^T*Lambda^k_c)*Av^T
        AtL  = Au' * Lambda_cell{c} * Av;    % (n1+1)x(n2+1): Au^T*Lam*Av
        AAtL = Au  * AtL * Av';              % (m1+1)x(m2+1): Au*(Au^T*Lam*Av)*Av^T

        %% Update Lambda (uses current R = A*P^k - Q)
        Lambda_cell{c} = (1 - omega) * Lambda_cell{c} ...
                         - gamma * nu * AAtL ...
                         - omega * R_cell{c};

        %% Update P
        P_cell{c} = P_cell{c} + nu * AtL;

        %% Update residual efficiently
        R_cell{c} = R_cell{c} + nu * AAtL;

        %% Accumulate gradient norm and fitting error
        Gc       = Au' * R_cell{c} * Av;
        g_norm2  = g_norm2 + norm(Gc, 'fro')^2;
        e_val    = e_val   + norm(R_cell{c}, 'fro')^2;
    end

    err_hist(k+1) = e_val;
    Ek_hist(k+1)  = g_norm2 / g0_norm2;

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
info.S1          = S1;
info.Sr          = Sr;
end