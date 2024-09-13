% For the commented version, see the folder Force Calculation Only/Nickel

function js = sphbes2(nu, x)

[nnu lnu] = size(nu);
[nx lx] = size(x);
xm = repmat(x, 1, lnu);
js = sqrt(pi ./(2* xm)) .* bessely(nu + 0.5, x);

end


