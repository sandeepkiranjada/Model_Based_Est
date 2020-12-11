function [ H ] = householder( a,b )
%HOUSEHOLDER Get H, in Ha = b
%   Detailed explanation goes here

eps = 2*(norm(a)-norm(b))/(norm(a)+norm(b));

if abs(eps)>1e-10
    error('Householder() called with vectors of unequal magnitudes')
end

w = (a-b);

beta = 0.5*(w'*w);

H = eye(length(a)) - (w*w')/beta;

end

