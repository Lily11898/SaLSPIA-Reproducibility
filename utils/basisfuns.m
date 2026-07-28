function N = basisfuns(span, t, p, knots)
% BASISFUNS  Evaluate all non-zero B-spline basis functions at parameter t.
%   N = basisfuns(span, t, p, knots)
%
%   Based on Piegl & Tiller "The NURBS Book", Algorithm A2.2.
%   All indices are 1-based (MATLAB convention).
%
%   Input:
%     span  - knot span index (from findspan)
%     t     - parameter value
%     p     - B-spline degree
%     knots - knot vector
%
%   Output:
%     N     - row vector of p+1 non-zero basis function values
%             corresponding to basis functions (span-p) to span (1-based)

N     = zeros(1, p+1);
left  = zeros(1, p+1);
right = zeros(1, p+1);
N(1)  = 1.0;

for j = 1:p
    left(j)  = t - knots(span + 1 - j);
    right(j) = knots(span + j) - t;
    saved = 0.0;
    for r = 0:j-1
        temp   = N(r+1) / (right(r+1) + left(j - r));
        N(r+1) = saved + right(r+1) * temp;
        saved  = left(j - r) * temp;
    end
    N(j+1) = saved;
end
end
