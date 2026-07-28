function A = build_collocation(knots, p, t)
% BUILD_COLLOCATION  Construct the B-spline collocation matrix.
%   A = build_collocation(knots, p, t)
%
%   Input:
%     knots - knot vector of length (n+1) + p + 1
%     p     - B-spline degree
%     t     - parameter values, column or row vector of length m+1
%
%   Output:
%     A     - (m+1)-by-(n+1) collocation matrix with A(j,i) = N_i^p(t(j))

t = t(:)';  % ensure row vector
m1 = length(t);          % m+1
n1 = length(knots) - p - 1;  % n+1 (number of basis functions)

A = sparse(m1, n1);

for j = 1:m1
    span = findspan(n1, p, t(j), knots);
    N    = basisfuns(span, t(j), p, knots);
    cols = (span - p):span;
    A(j, cols) = N;
end

A = full(A);
end
