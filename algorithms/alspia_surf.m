function [P, info] = alspia_surf(Au, Av, Q, P0, tol, maxiter)
% ALSPIA_SURF  ALSPIA for tensor-product B-spline surface fitting.
%   [P, info] = alspia_surf(Au, Av, Q, P0, tol, maxiter)
%
%   Uses Chebyshev semi-iterative step sizes (Wu & Liu 2024).
%
%   Input:
%     Au, Av - collocation matrices in u and v directions
%     Q      - cell array of d matrices (data points per coordinate)
%     P0     - cell array of d matrices (initial control points)
%     tol, maxiter - convergence parameters

if nargin < 5, tol = 1e-6;   end
if nargin < 6, maxiter = 10000; end

d = length(Q);

%% Eigenvalues of Hs = kron(Av'*Av, Au'*Au)
Hu = Au' * Au;
Hv = Av' * Av;
eig_u = eig(Hu);  eig_u = eig_u(eig_u > 1e-14);
eig_v = eig(Hv);  eig_v = eig_v(eig_v > 1e-14);
all_eig = kron(eig_u, eig_v);

lam_min = min(all_eig);
lam_max = max(all_eig);

is_full_rank = (lam_min > 1e-10 * lam_max);

%% Precompute Chebyshev step sizes
K = min(1e6, maxiter);
u_val = lam_min;
v_val = lam_max;

if is_full_rank
    omega_all = zeros(K, 1);
    for ell = 0:K-1
        omega_all(ell+1) = 2.0 / ((v_val + u_val) + ...
            (v_val - u_val) * cos((2*ell + 1) / (2*K) * pi));
    end
else
    r_Kp1 = cos((2*K + 1) / (2*(K+1)) * pi);
    omega_all = zeros(K, 1);
    for ell = 0:K-1
        omega_all(ell+1) = (1 - r_Kp1) / ...
            (v_val * (cos((2*ell+1)/(2*(K+1))*pi) - r_Kp1));
    end
end

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
    ell_k = mod(k-1, K);
    wk    = omega_all(ell_k + 1);

    for c = 1:d
        Rc = Au * P{c} * Av' - Q{c};
        Gc = Au' * Rc * Av;
        P{c} = P{c} - wk * Gc;
    end

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

function [g_norm2, G] = compute_grad_norm(Au, Av, P, Q)
    d = length(Q); g_norm2 = 0; G = cell(d,1);
    for c = 1:d
        Rc = Au * P{c} * Av' - Q{c};
        G{c} = Au' * Rc * Av;
        g_norm2 = g_norm2 + norm(G{c}, 'fro')^2;
    end
end

function e = compute_fitting_error(Au, Av, P, Q)
    d = length(Q); e = 0;
    for c = 1:d
        Rc = Au * P{c} * Av' - Q{c};
        e = e + norm(Rc, 'fro')^2;
    end
end
