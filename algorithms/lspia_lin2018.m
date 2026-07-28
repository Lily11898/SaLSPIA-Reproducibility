function [P, info] = lspia_lin2018(A, Q, P0, tol, maxiter)
% LSPIA_LIN2018  Reduced active-set implementation of Lin et al. (2018a)
% for singular least-squares fitting.
%
%   It applies the diagonal-weighted LSPIA only to active control-point
%   groups, i.e. columns i with sum_j A(j,i) > tol_col.
%   Inactive control points are kept fixed at their initial values.
%
%   Iteration on active set:
%       P_a^{k+1} = P_a^k + Lambda_a * A_a^T * (Q - A_a P_a^k)
%   where Lambda_a(ii) = 1 / sum_j A_a(j,i)

if nargin < 4, tol = 1e-6; end
if nargin < 5, maxiter = 10000; end

%% Column-group activity detection
colsum = sum(A, 1)';                  % (n+1)-by-1
tol_col = 1e-12 * max(1, max(colsum)); % relative threshold
active = (colsum > tol_col);
inactive = ~active;

n_all = size(A, 2);
n_act = nnz(active);

if n_act == 0
    error('lspia_lin2018:noActiveCols', ...
        'No active control-point groups remain after missing-data removal.');
end

A_act = A(:, active);
P_act = P0(active, :);
lambda_diag = 1 ./ colsum(active);

%% Initialize on active set
P = P0;                         % full output
r = A_act * P_act - Q;
g = A_act' * r;
g0_norm2 = norm(g, 'fro')^2;

err_hist = zeros(maxiter+1, 1);
Ek_hist  = zeros(maxiter+1, 1);
err_hist(1) = norm(r, 'fro')^2;

if g0_norm2 <= eps
    Ek_hist(1) = 0.0;
    info.iter_count   = 0;
    info.err_history  = err_hist(1);
    info.Ek_history   = Ek_hist(1);
    info.cpu_time     = 0.0;
    info.lambda_diag  = lambda_diag;
    info.active_mask  = active;
    info.inactive_mask = inactive;
    info.n_active     = n_act;
    info.n_inactive   = n_all - n_act;
    return;
else
    Ek_hist(1) = 1.0;
end

%% Iteration on active variables only
tic;
for k = 1:maxiter
    % P_act = P_act - Lambda * g, row-wise
    P_act = P_act - bsxfun(@times, lambda_diag, g);

    r = A_act * P_act - Q;
    g = A_act' * r;

    err_hist(k+1) = norm(r, 'fro')^2;
    Ek_hist(k+1)  = norm(g, 'fro')^2 / g0_norm2;

    if Ek_hist(k+1) < tol
        err_hist = err_hist(1:k+1);
        Ek_hist  = Ek_hist(1:k+1);
        break;
    end
end
cpu = toc;

%% Write back to full variable
P(active, :) = P_act;

%% Output
info.iter_count    = min(k, maxiter);
info.err_history   = err_hist;
info.Ek_history    = Ek_hist;
info.cpu_time      = cpu;
info.lambda_diag   = lambda_diag;
info.active_mask   = active;
info.inactive_mask = inactive;
info.n_active      = n_act;
info.n_inactive    = n_all - n_act;
end