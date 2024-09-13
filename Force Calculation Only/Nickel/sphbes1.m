% This function computes the spherical bessel function of the first kind.
% Source:

% Daniel Vangheluwe (2024). spherical Besselfunctions 
% (https://www.mathworks.com/matlabcentral/fileexchange/8488-spherical-besselfunctions), 
% MATLAB Central File Exchange. Retrieved August 9, 2024.

function js = sphbes1(nu, x)
% returns the spherical Bessel functions jnu(x)
% x is a vector or it may be a matrix if nu is a scalar
% if nu is a row and x a column vector, the output js is a matrix

[nnu lnu] = size(nu);
[nx lx] = size(x);
xm = repmat(x, 1, lnu);
js = sqrt(pi ./(2* xm)) .* besselj(nu + 0.5, x);

end


