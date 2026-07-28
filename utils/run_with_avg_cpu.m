function [P, info, cpu_samples] = run_with_avg_cpu(run_fun, n_repeats)
%RUN_WITH_AVG_CPU  Run a deterministic solver repeatedly and average CPU time.
%   [P, INFO] = RUN_WITH_AVG_CPU(RUN_FUN, N_REPEATS) runs RUN_FUN
%   N_REPEATS times. RUN_FUN must return [P, INFO], where INFO contains
%   the field cpu_time. The first run supplies P and the convergence history.
%   INFO.cpu_time and INFO.cpu_time_mean contain the five-run mean,
%   INFO.cpu_time_std contains the sample standard deviation, and
%   INFO.cpu_time_samples stores the individual timings.

if nargin < 2 || isempty(n_repeats)
    n_repeats = 5;
end

[P, info] = run_fun();
cpu_samples = zeros(n_repeats, 1);
cpu_samples(1) = info.cpu_time;

for ir = 2:n_repeats
    [~, info_ir] = run_fun();
    cpu_samples(ir) = info_ir.cpu_time;
end

info.cpu_time_mean = mean(cpu_samples);
info.cpu_time = info.cpu_time_mean; % backward-compatible alias
info.cpu_time_std = std(cpu_samples, 0);
info.cpu_time_samples = cpu_samples;
info.cpu_n_repeats = n_repeats;
end
