function [P, info] = lspia_surf(Au, Av, Q, P0, tol, maxiter)
% LSPIA_SURF  LSPIA for tensor-product B-spline surface fitting.
%   [P, info] = lspia_surf(Au, Av, Q, P0, tol, maxiter)
%
%   Surface residual:   R_c = Au * P_c * Av' - Q_c
%   Surface gradient:   G_c = Au' * R_c * Av
%   Surface Hessian:    Hs = kron(Av'*Av, Au'*Au)
%   Step size:          omega = 2 / (lam_max(Hs) + lam_min(Hs))
%
%   Input:
%     Au  - (mu+1)-by-(nu+1) collocation matrix in u-direction
%     Av  - (mv+1)-by-(nv+1) collocation matrix in v-direction
%     Q   - cell array of d matrices, each (mu+1)-by-(mv+1)
%     P0  - cell array of d matrices, each (nu+1)-by-(nv+1) initial ctrl pts
%     tol, maxiter - convergence parameters
%
%   Output:
%     P    - cell array of d matrices (final control points)
%     info - struct with iter_count, err_history, Ek_history, cpu_time

if nargin < 5, tol = 1e-6;   end
if nargin < 6, maxiter = 10000; end

d = length(Q);  % number of coordinates

%% Eigenvalues of Hs = kron(Av'*Av, Au'*Au)
Hu = Au' * Au;
Hv = Av' * Av;
eig_u = eig(Hu);  eig_u = eig_u(eig_u > 1e-14);
eig_v = eig(Hv);  eig_v = eig_v(eig_v > 1e-14);
% Eigenvalues of Kronecker product are all products
all_eig = kron(eig_u, eig_v);
lam_min = min(all_eig);
lam_max = max(all_eig);
omega   = 2.0 / (lam_max + lam_min);

%% Initialize
P = P0;
[g0_norm2, ~] = compute_grad_norm(Au, Av, P, Q);
err0 = compute_fitting_error(Au, Av, P, Q);

err_hist = zeros(maxiter+1, 1);
Ek_hist  = zeros(maxiter+1, 1);
err_hist(1) = err0;
Ek_hist(1)  = 1.0;

%% Iteration
tic;
for k = 1:maxiter
    % Gradient and update for each coordinate
    for c = 1:d
        Rc = Au * P{c} * Av' - Q{c};
        Gc = Au' * Rc * Av;
        P{c} = P{c} - omega * Gc;
    end

    % Errors
    [gk_norm2, ~] = compute_grad_norm(Au, Av, P, Q);
    err_hist(k+1) = compute_fitting_error(Au, Av, P, Q);
    Ek_hist(k+1)  = gk_norm2 / g0_norm2;

    if Ek_hist(k+1) < tol
        err_hist = err_hist(1:k+1);
        Ek_hist  = Ek_hist(1:k+1);
        break;
    end
end
cpu = toc;

info.iter_count  = min(k, maxiter);
info.err_history = err_hist;
info.Ek_history  = Ek_hist;
info.cpu_time    = cpu;
end

%% Helper functions
function [g_norm2, G] = compute_grad_norm(Au, Av, P, Q)
    d = length(Q);
    g_norm2 = 0;
    G = cell(d,1);
    for c = 1:d
        Rc = Au * P{c} * Av' - Q{c};
        G{c} = Au' * Rc * Av;
        g_norm2 = g_norm2 + norm(G{c}, 'fro')^2;
    end
end

function e = compute_fitting_error(Au, Av, P, Q)
    d = length(Q);
    e = 0;
    for c = 1:d
        Rc = Au * P{c} * Av' - Q{c};
        e = e + norm(Rc, 'fro')^2;
    end
end
