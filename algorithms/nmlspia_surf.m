function [P, info] = nmlspia_surf(Au, Av, Q, P0, tol, maxiter)
% NMLSPIA_SURF  NmLSPIA for tensor-product B-spline surface fitting.
%   [P, info] = nmlspia_surf(Au, Av, Q, P0, tol, maxiter)
%
%   This is the tensor-product surface counterpart of nmlspia.m. For each
%   coordinate c, the least-squares problem is
%       min 1/2 || Au * P_c * Av' - Q_c ||_F^2.
%   The Nesterov momentum update is applied through tensor-product matrix
%   products without forming kron(Av, Au).
%
%   Optimal parameters follow Liu, Wu, Li & Hu (2024), Theorem 2, using the
%   extreme eigenvalues of Hs = kron(Av'*Av, Au'*Au):
%       zeta = 1 / lambda_max(Hs),
%       eta  = (sqrt(lambda_max) - sqrt(lambda_min)) /
%              (sqrt(lambda_max) + sqrt(lambda_min)).
%
%   Input conventions match lspia_surf.m, alspia_surf.m, and mlspia_surf.m.

if nargin < 5 || isempty(tol),     tol = 1e-6;    end
if nargin < 6 || isempty(maxiter), maxiter = 10000; end

d = numel(Q);

%% Extreme eigenvalues of Hs = kron(Av'*Av, Au'*Au)
Hu = Au' * Au;
Hv = Av' * Av;
eig_u = eig(Hu); eig_u = eig_u(eig_u > 1e-14);
eig_v = eig(Hv); eig_v = eig_v(eig_v > 1e-14);

all_eig = kron(eig_u, eig_v);
lam_min = min(all_eig);
lam_max = max(all_eig);

zeta = 1.0 / lam_max;
eta  = (sqrt(lam_max) - sqrt(lam_min)) / ...
       (sqrt(lam_max) + sqrt(lam_min));

%% Initialize
P = P0;
P_prev = P0;

[g0_norm2, err0] = compute_grad_error(Au, Av, P, Q);

err_hist = zeros(maxiter + 1, 1);
Ek_hist  = zeros(maxiter + 1, 1);
err_hist(1) = err0;
Ek_hist(1)  = 1.0;

%% Iteration
tic;
for k = 1:maxiter
    P_new = cell(d, 1);

    for c = 1:d
        mbar = P{c} - P_prev{c};
        Plook = P{c} + eta * mbar;
        Rlook = Au * Plook * Av' - Q{c};
        Glook = Au' * Rlook * Av;
        P_new{c} = P{c} - zeta * Glook + eta * mbar;
    end

    P_prev = P;
    P = P_new;

    [g_norm2, err_val] = compute_grad_error(Au, Av, P, Q);
    err_hist(k + 1) = err_val;
    Ek_hist(k + 1) = g_norm2 / g0_norm2;

    if Ek_hist(k + 1) < tol
        err_hist = err_hist(1:k + 1);
        Ek_hist  = Ek_hist(1:k + 1);
        break;
    end
end
cpu = toc;

info.iter_count  = min(k, maxiter);
info.err_history = err_hist;
info.Ek_history  = Ek_hist;
info.cpu_time    = cpu;
info.zeta        = zeta;
info.eta         = eta;
info.lam_min     = lam_min;
info.lam_max     = lam_max;
end

function [g_norm2, err_val] = compute_grad_error(Au, Av, P, Q)
    d = numel(Q);
    g_norm2 = 0;
    err_val = 0;
    for c = 1:d
        R = Au * P{c} * Av' - Q{c};
        G = Au' * R * Av;
        g_norm2 = g_norm2 + norm(G, 'fro')^2;
        err_val = err_val + norm(R, 'fro')^2;
    end
end
