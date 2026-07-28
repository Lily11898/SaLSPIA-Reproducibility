function [t, knots] = make_params_knots(Q, n, p)
% MAKE_PARAMS_KNOTS  Compute chord-length parameters and knot vector.
%   [t, knots] = make_params_knots(Q, n, p)
%
%   Following Piegl & Tiller "The NURBS Book", Eqs. (9.5) and (9.69).
%
%   Input:
%     Q - (m+1)-by-d matrix of data points
%     n - n+1 is the number of control points (0-indexed max index)
%     p - B-spline degree
%
%   Output:
%     t     - parameter values of length m+1 in [0,1]
%     knots - clamped knot vector of length (n+1)+p+1

[m1, ~] = size(Q);   % m+1
m = m1 - 1;

%% Chord-length parameterization (Eq. 9.5)
chords = sqrt(sum(diff(Q).^2, 2));
total  = sum(chords);

t = zeros(m1, 1);
t(1) = 0;
for j = 2:m1
    t(j) = t(j-1) + chords(j-1) / total;
end
t(end) = 1.0;

%% Knot vector (Eq. 9.69 for least squares)
% Clamped: first p+1 = 0, last p+1 = 1, with n-p internal knots
n_internal = n - p;   % number of internal knots (n+1 control pts, degree p)
knots = zeros(1, (n+1) + p + 1);

% Clamped ends
knots(1:p+1)       = 0;
knots(end-p:end)    = 1;

if n_internal > 0
    d = (m + 1) / (n - p + 1);
    for j = 1:n_internal
        i     = floor(j * d);       % 1-based index into t
        alpha = j * d - i;
        % t(i) corresponds to t_{i-1} (0-indexed), t(i+1) to t_i (0-indexed)
        knots(p + 1 + j) = (1 - alpha) * t(i) + alpha * t(i+1);
    end
end
end
