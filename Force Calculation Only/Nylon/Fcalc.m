% This function computes F_n, needed to evaluate the radiation force
% function. F_n is used to evaluate the scattering coefficient c_n.

% The formulation is taken from:

% Hasegawa, Takahi, and Katuya Yosioka. "Acoustic‐radiation force on a solid
% elastic sphere." The Journal of the Acoustical Society of America 
% 46.5B (1969): 1139-1143.

function Fn=Fcalc(jn,jn1,djn,ddjn,x,x1,x2,rhop,rholiq,n)

jnx1=subs(jn,x,x1);
djnx1=subs(djn,x,x1);
ddjnx1=subs(ddjn,x,x1);
jnx2=subs(jn,x,x2);
djnx2=subs(djn,x,x2);
ddjnx2=subs(ddjn,x,x2); % See hasegawa79_rad_force_fct.m for definitions
% Here the code only substitutes non-dimensional sphere radius scaled with 
% compression (x1) or shear wave-based (x2) wavenumbers.

jn1x1=subs(jn1,x,x1);
jn1x2=subs(jn1,x,x2);

An=(n.*jnx1-x1.*jn1x1)./...
    ((n-1).*jnx1-x1.*jn1x1);

Bn=(2.*n.*(n+1).*jnx2)./...
    ((2.*n.^2-x2.^2-2).*jnx2+2.*x2.*jn1x2);

Dn=((x2.^2./2-n.*(n-1)).*jnx1-2.*x1.*jn1x1)./...
    ((n-1).*jnx1-x1.*jn1x1);

En=(2.*n.*(n+1).*((1-n).*jnx2+x2.*jn1x2))./...
    ((2.*n.^2-x2.^2-2).*jnx2+2.*x2.*jn1x2);

Fn=((x2.^2).*rholiq.*(An-Bn)) ./...
    (2.*rhop.*(Dn-En));
end