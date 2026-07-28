%% ========================================================================
%  run_table2_parameter_sensitivity_surface.m
%  -----------------------------------------------------------------------
%  Parameter sensitivity of SaLSPIA on a TENSOR-PRODUCT SURFACE example,
%  reproducing the Peaks Ex. 4.5 row in the manuscript's Table 2.
%
%  The Peaks surface (Example 4.5, full-rank) is used, sweeping the three
%  user-defined parameters of SaLSPIA, exactly as in the curve study:
%
%        sigma  (GLL sufficient-decrease, params.c)
%        M      (GLL memory length,       params.M)
%        delta  (truncation threshold,    params.delta)
%
%  For each value we report:  E_inf (final relative gradient norm),
%  IT (iteration count), Halvings (total GLL backtracking steps).
%
%  The method called is salspia_surf.m  (paper-faithful SaLSPIA, eqs 3-12).
%  Defaults: sigma = 1e-4, M = 10, delta = 1e-8  (same as the paper).
%  ========================================================================
clear; clc;

script_dir = fileparts(mfilename('fullpath'));
repo_dir = fileparts(script_dir);
addpath(fullfile(repo_dir, 'algorithms'));
addpath(fullfile(repo_dir, 'utils'));

out_dir = fullfile(repo_dir, 'results', 'parameter_sensitivity_surface');
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

p_deg   = 3;
tol     = 1e-6;
maxiter = 10000;

% ----- default parameters (held fixed except for the one being swept) ---
def.c     = 1e-4;     % sigma
def.M     = 10;       % memory length
def.delta = 1e-8;     % truncation threshold
def.eps_saf = 1e-30;  % safeguard tolerance (Algorithm 1)

% ----- sweep ranges (mirror the curve table) ---------------------------
sigma_values = [1/2, 1e-1, 1e-2, 1e-4, 1e-6];
M_values     = [1, 3, 5, 10, 20];
delta_values = [1e-10, 1e-8, 1e-6, 1e-4, 1e-2];

%% ----- build the Peaks surface system (Example 4.5) --------------------
ex_id = 5;  m1u = 121;  m1v = 121;  n1u = 41;  n1v = 41;

[Qx, Qy, Qz, u_raw, v_raw] = gen_peaks(m1u, m1v);
u_par = (u_raw - u_raw(1)) / (u_raw(end) - u_raw(1));
v_par = (v_raw - v_raw(1)) / (v_raw(end) - v_raw(1));

knots_u = make_clamped_knots(n1u, p_deg);
knots_v = make_clamped_knots(n1v, p_deg);
Au = build_collocation(knots_u, p_deg, u_par);
Av = build_collocation(knots_v, p_deg, v_par);

Q_cell  = {Qx, Qy, Qz};
P0_cell = cell(3,1);
for c = 1:3, P0_cell{c} = select_initial_surf(Q_cell{c}, n1u, n1v); end

fprintf('Peaks surface (Ex 4.5): (m1,m2,n1,n2) = (%d,%d,%d,%d)\n\n', ...
        m1u, m1v, n1u, n1v);

%% ----- sweep sigma ------------------------------------------------------
fprintf('===== Sensitivity to sigma (M=%d, delta=%g fixed) =====\n', def.M, def.delta);
sig = run_sweep(Au, Av, Q_cell, P0_cell, tol, maxiter, def, 'c', sigma_values);

%% ----- sweep M ----------------------------------------------------------
fprintf('\n===== Sensitivity to M (sigma=%g, delta=%g fixed) =====\n', def.c, def.delta);
Mres = run_sweep(Au, Av, Q_cell, P0_cell, tol, maxiter, def, 'M', M_values);

%% ----- sweep delta ------------------------------------------------------
fprintf('\n===== Sensitivity to delta (sigma=%g, M=%d fixed) =====\n', def.c, def.M);
del = run_sweep(Au, Av, Q_cell, P0_cell, tol, maxiter, def, 'delta', delta_values);

save(fullfile(out_dir, 'peaks_parameter_sensitivity_results.mat'), ...
    'sig', 'Mres', 'del', 'sigma_values', 'M_values', 'delta_values', 'def');

%% ----- LaTeX rows for the M table (to slot into Table 2) ---------------
fprintf('\n\n%% ---- LaTeX row for the M-sensitivity table (Peaks surface) ----\n');
fprintf('Peaks Ex.~4.5\n');
fprintf('& $E_\\infty$');
for i = 1:numel(M_values), fprintf(' & %s', scifmt(Mres(i).Einf)); end
fprintf(' \\\\\n& IT');
for i = 1:numel(M_values), fprintf(' & %d', Mres(i).IT); end
fprintf(' \\\\\n& Halvings');
for i = 1:numel(M_values), fprintf(' & %d', Mres(i).halv); end
fprintf(' \\\\\n');
fprintf('\nSaved outputs to:\n  %s\n', out_dir);

%% ========================================================================
function res = run_sweep(Au, Av, Q, P0, tol, maxiter, def, field, vals)
    fprintf('%-12s %-12s %-8s %-10s\n', field, 'E_inf', 'IT', 'Halvings');
    fprintf('%s\n', repmat('-', 1, 44));
    n = numel(vals);
    res = struct('val',cell(n,1),'Einf',[],'IT',[],'halv',[]);
    for i = 1:n
        params = def;
        params.(field) = vals(i);
        [~, info] = salspia_surf(Au, Av, Q, P0, tol, maxiter, params);
        Einf = info.Ek_history(end);
        IT   = info.iter_count;
        halv = sum(info.halvings);
        res(i).val  = vals(i);
        res(i).Einf = Einf;
        res(i).IT   = IT;
        res(i).halv = halv;
        fprintf('%-12g %-12.2e %-8d %-10d\n', vals(i), Einf, IT, halv);
    end
end

function [Qx,Qy,Qz,u_par,v_par] = gen_peaks(m1u, m1v)
    th1 = linspace(-3,3,m1u); th2 = linspace(-4,4,m1v);
    [T1,T2] = meshgrid(th1,th2); T1 = T1'; T2 = T2';
    F = 3*(1-T1).^2.*exp(-T1.^2-(T2+1).^2) ...
      -10*(T1/5-T1.^3-T2.^5).*exp(-T1.^2-T2.^2) ...
      -1/3*exp(-(T1+1).^2-T2.^2);
    Qx = T1; Qy = T2; Qz = F;
    u_par = th1(:); v_par = th2(:);
end

function knots = make_clamped_knots(n1, p)
    n = n1-1; n_int = n-p;
    knots = zeros(1, n1+p+1);
    knots(1:p+1) = 0; knots(end-p:end) = 1;
    if n_int > 0
        internal = linspace(0,1,n_int+2);
        knots(p+2:end-p-1) = internal(2:end-1);
    end
end

function P0 = select_initial_surf(Qc, n1u, n1v)
    [mu1,mv1] = size(Qc);
    idx_u = round(linspace(1, mu1, n1u));
    idx_v = round(linspace(1, mv1, n1v));
    P0    = Qc(idx_u, idx_v);
end

function s = scifmt(x)
    if x == 0, s = '$0$'; return; end
    e = floor(log10(abs(x)));
    m = x / 10^e;
    s = sprintf('$%.2f\\times10^{%d}$', m, e);
end
