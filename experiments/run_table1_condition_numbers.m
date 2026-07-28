%% ========================================================================
%  run_table1_condition_numbers.m
%  Compute ranks and effective condition numbers of the collocation matrices
%  used in the manuscript experiments.
%
%  For full-column-rank matrices, kappa_eff(A) is the usual 2-norm condition
%  number kappa_2(A). For rank-deficient matrices, the usual condition number
%  is infinite, so we report
%
%      kappa_eff(A) = sigma_max(A) / sigma_r(A),
%
%  where r = rank(A) and sigma_r(A) is the smallest nonzero singular value.
%
%  Outputs:
%    results/condition_numbers/condition_numbers_summary.csv
%    results/condition_numbers/condition_numbers_table.tex
%  ========================================================================
clear; clc; close all;

script_dir = fileparts(mfilename('fullpath'));
repo_dir   = fileparts(script_dir);
addpath(fullfile(repo_dir, 'algorithms'));
addpath(fullfile(repo_dir, 'utils'));

p_deg = 3;
rows = {};

out_dir = fullfile(repo_dir, 'results', 'condition_numbers');
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

%% ===== Full-rank curve examples: Examples 4.1--4.3 =====================
curve_cfg = {
    'Ex.~4.1', 'full curve', 15001, 3001, 1;
    'Ex.~4.2', 'full curve', 10001, 2001, 2;
    'Ex.~4.3', 'full curve',  1000,  288, 3;
};

for i = 1:size(curve_cfg, 1)
    ex_label = curve_cfg{i, 1};
    setting  = curve_cfg{i, 2};
    m1       = curve_cfg{i, 3};
    n1       = curve_cfg{i, 4};
    ex_id    = curve_cfg{i, 5};

    Q = generate_curve_data(ex_id, m1);
    [t, knots] = make_params_knots(Q, n1 - 1, p_deg);
    A = build_collocation(knots, p_deg, t);
    rows(end+1, :) = append_diag(ex_label, setting, sprintf('%d', m1), ...
        sprintf('%d', n1), size(A, 1), size(A, 2), matrix_diagnostics(A)); %#ok<SAGROW>
end

%% ===== Full-rank tensor-product surface examples: Examples 4.4--4.6 ====
surface_cfg = {
    'Ex.~4.4', 'full surface',  81,  81,  31,  31;
    'Ex.~4.5', 'full surface', 121, 121,  41,  41;
    'Ex.~4.6', 'full surface', 860, 600, 121, 121;
};

for i = 1:size(surface_cfg, 1)
    ex_label = surface_cfg{i, 1};
    setting  = surface_cfg{i, 2};
    m1u      = surface_cfg{i, 3};
    m1v      = surface_cfg{i, 4};
    n1u      = surface_cfg{i, 5};
    n1v      = surface_cfg{i, 6};

    u_par = linspace(0, 1, m1u);
    v_par = linspace(0, 1, m1v);
    Au = build_collocation(make_clamped_knots(n1u, p_deg), p_deg, u_par);
    Av = build_collocation(make_clamped_knots(n1v, p_deg), p_deg, v_par);

    % For tensor-product data, A_s = Av \otimes Au. Singular values of
    % A_s are pairwise products of those of Au and Av, so the condition
    % number is the product of the one-dimensional condition numbers.
    d = tensor_product_diagnostics(Au, Av);
    rows(end+1, :) = append_diag(ex_label, setting, sprintf('%d x %d', m1u, m1v), ...
        sprintf('%d x %d', n1u, n1v), m1u*m1v, n1u*n1v, d); %#ok<SAGROW>
end

%% ===== Rank-deficient missing-data curve examples ======================
[A, data_size, ctrl_size] = build_missing_blob_curve(p_deg);
rows(end+1, :) = append_diag('Ex.~4.1', 'missing curve', sprintf('%d', data_size), ...
    sprintf('%d', ctrl_size), size(A, 1), size(A, 2), matrix_diagnostics(A)); %#ok<SAGROW>

[A, data_size, ctrl_size] = build_missing_helix_curve(p_deg);
rows(end+1, :) = append_diag('Ex.~4.2', 'missing curve', sprintf('%d', data_size), ...
    sprintf('%d', ctrl_size), size(A, 1), size(A, 2), matrix_diagnostics(A)); %#ok<SAGROW>

[A, data_size, ctrl_size] = build_missing_g_loop_curve(p_deg);
rows(end+1, :) = append_diag('Ex.~4.7', 'missing curve', sprintf('%d', data_size), ...
    sprintf('%d', ctrl_size), size(A, 1), size(A, 2), matrix_diagnostics(A)); %#ok<SAGROW>

%% ===== Rank-deficient missing-data surface example: Example 4.8 ========
[A, data_size, ctrl_size] = build_missing_franke_surface(p_deg);
rows(end+1, :) = append_diag('Ex.~4.8', 'missing surface', sprintf('%d', data_size), ...
    sprintf('%d', ctrl_size), size(A, 1), size(A, 2), matrix_diagnostics(A)); %#ok<SAGROW>

%% ===== Save =============================================================
T = cell2table(rows, 'VariableNames', ...
    {'Example', 'Setting', 'DataSize', 'ControlSize', 'Rows', 'Columns', ...
     'Rank', 'KappaEff', 'ClassicalKappa'});

disp(T);
writetable(T, fullfile(out_dir, 'condition_numbers_summary.csv'));
write_latex_table(T, fullfile(out_dir, 'condition_numbers_table.tex'));

fprintf('\nSaved condition-number diagnostics to:\n  %s\n', out_dir);

%% ========================================================================
%  Local example builders
%% ========================================================================
function Q = generate_curve_data(ex_id, m1)
switch ex_id
    case 1
        theta = linspace(0, 2*pi, m1)';
        r = 2 + 4*cos(2*theta + pi/4) + cos(3*theta + pi/4);
        Q = [r .* cos(theta), r .* sin(theta)];
    case 2
        theta = linspace(0, 2*pi, m1)';
        Q = [10*cos(theta*pi/3), 10*sin(theta*pi/3), theta*pi/3];
    case 3
        Q = load(salspia_data('cur_data deer'));
        if m1 ~= size(Q, 1)
            idx = round(linspace(1, size(Q, 1), m1));
            Q = Q(idx, :);
        end
    otherwise
        error('Unknown curve example id: %d', ex_id);
end
end

function [A, data_size, ctrl_size] = build_missing_blob_curve(p_deg)
n = 100;
m_full = 330;
theta = linspace(0, 2*pi, m_full + 1)';
r = 2 + 4*cos(2*theta + pi/4) + cos(3*theta + pi/4);
Q_full = [r .* cos(theta), r .* sin(theta)];
[t_full, knots] = make_params_knots(Q_full, n, p_deg);
hole_centers = [0.20, 0.42, 0.62, 0.82];
hole_width = 0.02;
[~, Q_obs, ~, t_obs] = remove_by_param_intervals(Q_full, t_full, hole_centers, hole_width);
A = build_collocation(knots, p_deg, t_obs);
data_size = size(Q_obs, 1);
ctrl_size = n + 1;
end

function [A, data_size, ctrl_size] = build_missing_helix_curve(p_deg)
n = 200;
m_full = 520;
theta = linspace(0, 5*6, m_full + 1)';
Q_full = [10*cos(theta*pi/3), 10*sin(theta*pi/3), theta*pi/3];
[t_full, knots] = make_params_knots(Q_full, n, p_deg);
hole_centers = [0.25, 0.50, 0.75];
hole_width = 0.02;
[~, Q_obs, ~, t_obs] = remove_by_param_intervals(Q_full, t_full, hole_centers, hole_width);
A = build_collocation(knots, p_deg, t_obs);
data_size = size(Q_obs, 1);
ctrl_size = n + 1;
end

function [A, data_size, ctrl_size] = build_missing_g_loop_curve(p_deg)
n = 100;
Q_full = load(salspia_data('s_loop_curve_data.txt'));
[t_full, knots] = make_params_knots(Q_full, n, p_deg);
hole_centers = [0.15, 0.38, 0.62, 0.99];
hole_width = 0.03;
[~, Q_obs, ~, t_obs] = remove_by_param_intervals(Q_full, t_full, hole_centers, hole_width);
A = build_collocation(knots, p_deg, t_obs);
data_size = size(Q_obs, 1);
ctrl_size = n + 1;
end

function [A_obs, data_size, ctrl_size] = build_missing_franke_surface(p_deg)
% This follows the Franke setup in run_table7_missing_surface.m.
m1u = 101;
m1v = 101;
n1u = 41;
n1v = 41;
u_par = linspace(0, 1, m1u)';
v_par = linspace(0, 1, m1v)';

hole_specs = [ ...
    0.45 0.25 0.11 0.12;
    0.72 0.62 0.11 0.11;
    0.25 0.55 0.12 0.11;
    0.60 0.85 0.11 0.11 ];

Au = build_collocation(make_clamped_knots(n1u, p_deg), p_deg, u_par);
Av = build_collocation(make_clamped_knots(n1v, p_deg), p_deg, v_par);
[U, V] = ndgrid(u_par, v_par);
hole_mask = false(size(U));
for h = 1:size(hole_specs, 1)
    c = hole_specs(h, :);
    hole_mask = hole_mask | (((U - c(1)) / c(3)).^2 + ((V - c(2)) / c(4)).^2 <= 1);
end
obs_mask = ~hole_mask;
A_obs = build_observed_surface_collocation(Au, Av, obs_mask);
data_size = nnz(obs_mask);
ctrl_size = n1u * n1v;
end

function [keep, Q_keep, Q_miss, t_keep] = remove_by_param_intervals( ...
    Q_full, t_full, hole_centers, hole_half_width)
keep = true(size(t_full));
for h = 1:length(hole_centers)
    lo = hole_centers(h) - hole_half_width;
    hi = hole_centers(h) + hole_half_width;
    keep((t_full >= lo) & (t_full <= hi)) = false;
end
Q_keep = Q_full(keep, :);
Q_miss = Q_full(~keep, :);
t_keep = t_full(keep);
end

function knots = make_clamped_knots(n1, p)
n = n1 - 1;
n_int = n - p;
knots = zeros(1, n1 + p + 1);
knots(1:p+1) = 0;
knots(end-p:end) = 1;
if n_int > 0
    internal = linspace(0, 1, n_int + 2);
    knots(p+2:end-p-1) = internal(2:end-1);
end
end

function [A_obs, obs_iu, obs_iv] = build_observed_surface_collocation(Au, Av, obs_mask)
[obs_iu, obs_iv] = find(obs_mask);
n_obs = numel(obs_iu);
n_ctrl = size(Au, 2) * size(Av, 2);
A_obs = zeros(n_obs, n_ctrl);
for k = 1:n_obs
    A_obs(k, :) = kron(Av(obs_iv(k), :), Au(obs_iu(k), :));
end
end

%% ========================================================================
%  Diagnostics
%% ========================================================================
function d = matrix_diagnostics(A)
% Compute singular values directly from A.
%
% Do NOT obtain them from eig(A'*A): forming the Gram matrix squares the
% condition number and can turn exact zero singular values into small
% positive values through roundoff.  That gave a spurious rank for the
% rank-deficient Franke example.  Direct SVD also matches MATLAB rank(A),
% whose default tolerance is max(size(A))*eps(sigma_max).
s = svd(full(A), 'econ');

if isempty(s) || s(1) == 0
    d.rank = 0;
    d.kappa_eff = Inf;
    d.classical_kappa = Inf;
    return;
end

% Singular values returned by SVD are in descending order.
smax = s(1);
tol = max(size(A)) * eps(smax);
d.rank = sum(s > tol);

if d.rank == 0
    d.kappa_eff = Inf;
    d.classical_kappa = Inf;
    return;
end

% sigma_r is the smallest singular value above MATLAB's rank tolerance.
sigma_r = s(d.rank);
d.kappa_eff = smax / sigma_r;

if d.rank < size(A, 2)
    d.classical_kappa = Inf;
else
    d.classical_kappa = d.kappa_eff;
end
end

function d = tensor_product_diagnostics(Au, Av)
du = matrix_diagnostics(Au);
dv = matrix_diagnostics(Av);
d.rank = du.rank * dv.rank;
d.kappa_eff = du.kappa_eff * dv.kappa_eff;
if isinf(du.classical_kappa) || isinf(dv.classical_kappa)
    d.classical_kappa = Inf;
else
    d.classical_kappa = d.kappa_eff;
end
end

function row = append_diag(ex_label, setting, data_size, ctrl_size, n_rows, n_cols, d)
row = {ex_label, setting, data_size, ctrl_size, n_rows, n_cols, ...
    d.rank, d.kappa_eff, d.classical_kappa};
fprintf('%-8s %-16s data=%-12s ctrl=%-10s rank=%5d/%-5d kappa_eff=%.3e\n', ...
    ex_label, setting, data_size, ctrl_size, d.rank, n_cols, d.kappa_eff);
end

function write_latex_table(T, fname)
fid = fopen(fname, 'w');
if fid == -1
    warning('Cannot write LaTeX table: %s', fname);
    return;
end

fprintf(fid, '\\begin{table}[H]\n');
fprintf(fid, '\\centering\n');
fprintf(fid, '\\caption{Ranks and effective condition numbers of the collocation matrices. For rank-deficient cases, the usual condition number is infinite and $\\kappa_{\\rm eff}=\\sigma_{\\max}/\\sigma_r$ is reported over the nonzero singular spectrum.}\n');
fprintf(fid, '\\label{tab:condition_numbers}\n');
fprintf(fid, '\\begin{tabular}{llcccc}\n');
fprintf(fid, '\\hline\n');
fprintf(fid, 'Example & Setting & Data size & Control size & Rank & $\\kappa_{\\rm eff}$ \\\\\n');
fprintf(fid, '\\hline\n');
for i = 1:height(T)
    fprintf(fid, '%s & %s & %s & %s & %d/%d & %.2e \\\\\n', ...
        T.Example{i}, T.Setting{i}, T.DataSize{i}, T.ControlSize{i}, ...
        T.Rank(i), T.Columns(i), T.KappaEff(i));
end
fprintf(fid, '\\hline\n');
fprintf(fid, '\\end{tabular}\n');
fprintf(fid, '\\end{table}\n');
fclose(fid);
end
