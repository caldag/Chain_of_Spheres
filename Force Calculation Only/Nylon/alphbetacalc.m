% This function computes the scattering coefficient as described in Eq. (7)
% in:

% Hasegawa, T. (1979). Acoustic radiation force on a sphere in a 
% quasistationary wave field—theory. The Journal of the Acoustical Society 
% of America, 65(1), 32-40.

% The code returns the real (alphn) and imaginary (betan) parts of the
% scattering coefficient.

function [alphn, betan] = alphbetacalc(Fn,jn,jn1,djn,nn,nn1,dnn,x,n)

Gn=(Fn-n).*jn+x.*jn1;
Hn=(Fn-n).*nn+x.*nn1;

alphn=-Gn.^2/(Gn.^2+Hn.^2);
betan=-Gn.*Hn./(Gn.^2+Hn.^2);
end