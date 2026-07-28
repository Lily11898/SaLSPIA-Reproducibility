function span = findspan(n, p, t, knots)
% FINDSPAN  Find the knot span index for parameter t.
%   span = findspan(n, p, t, knots)
%
%   Based on Piegl & Tiller "The NURBS Book", Algorithm A2.1.
%   All indices are 1-based (MATLAB convention).
%
%   Input:
%     n     - number of basis functions (= number of control points)
%     p     - B-spline degree
%     t     - parameter value
%     knots - knot vector of length n+p+1
%
%   Output:
%     span  - knot span index (1-based) such that knots(span) <= t < knots(span+1)

% Handle right endpoint
if t >= knots(n+1)
    span = n;
    return;
end

% Handle left endpoint
if t <= knots(p+1)
    span = p + 1;
    return;
end

% Binary search
low  = p + 1;
high = n + 1;
mid  = floor((low + high) / 2);

while t < knots(mid) || t >= knots(mid+1)
    if t < knots(mid)
        high = mid;
    else
        low = mid;
    end
    mid = floor((low + high) / 2);
end

span = mid;
end
