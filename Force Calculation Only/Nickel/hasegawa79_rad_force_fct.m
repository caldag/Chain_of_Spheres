% This function calculates the acoustic radiation force function on a spherical
% particle under standing wave. The formulation is taken from:
%
% Hasegawa, T. (1979) Acoustic radiation force on a sphere in a
% quasistationary wave field-theory.
%
% The formulation here includes the correction pointed out by Mitri (2005)
% and Glynne-Jones (2013).

function Yst=hasegawa79_rad_force_fct(xeval)

% Nickel properties
rhop=8900; % Particle density, kg/m3
c1=5639; % Compression wave velocity of particle, m/s
c2=2970; % Shear wave velocity of particle, m/s

% Water properties
rholiq=1000; % Liquid density
cliq=1480; % Liquid speed of sound

x=sym('x'); % x=k*a, non-dimensional sphere radius

% x1 and x2 are also non-dimensional sphere radii but they use different
% scales:
x1=x*cliq/c1; % Uses compression wave velocity
x2=x*cliq/c2; % Uses shear wave velocity

n=sym('n'); % Bessel order
n1=sym('n1'); % Bessel order, n+1

%% Bessel functions:

% Of first kind:
jn=sphbes1(n,x); % n
jn1=sphbes1(n1,x); % n+1
% First and second derivatives, pre-computed symbolically:
djn=(2.^(1/2).*pi.^(1/2).*besselj(n - 1/2, x).*(1./x).^(1./2))./2 - (2.^(1./2).*pi.^(1./2).*besselj(n + 1/2, x).*(n + 1).*(1./x).^(3./2))./2;
ddjn=(2.^(1./2).*pi.^(1/2).*besselj(n + 1/2, x).*(1./x).^(1./2).*(n.^2 + 3.*n - x.^2 + 2))./(2.*x.^2) - 2.^(1./2).*pi.^(1/2).*besselj(n - 1./2, x).*(1./x).^(3./2);

% Of second kind:
nn=sphbes2(n,x); % n
nn1=sphbes2(n1,x); % n+1
% First derivative needed only (pre-computed):
dnn=(2.^(1./2).*pi.^(1./2).*bessely(n - 1./2, x).*(1./x).^(1./2))./2 - (2.^(1./2).*pi.^(1./2).*bessely(n + 1./2, x).*(n + 1).*(1./x).^(3./2))./2;

%% Evaluation of the key parameters:

% F_n is used in computing the real and imaginary components of the
% scattering coefficient, c_n=\alpha_n+1i*\beta_n

Fn=Fcalc(jn,jn1,djn,ddjn,x,x1,x2,rhop,rholiq,n);

% alphn: \alpha_n, real part of scattering coefficient
% betan: \beta_n, imaginary part of the scattering coefficient

[alphn,betan]=alphbetacalc(Fn,jn,jn1,djn,nn,nn1,dnn,x,n);

% This is normally an infinite summation but we can truncate it down to 20
% with great accuracy. The loop substitutes values in symbolic expressions
% to obtain numerical values.
for ii=0:20
   alphn_ev(ii+1)=subs(alphn,{n n1},{ii ii+1});
   betan_ev(ii+1)=subs(betan,{n n1},{ii ii+1});
end

%% Hasegawa (1979), radiation force function under standing wave

Yphi=0;

for ii=1:20 % Eq. (21) in Caldag & Yesilyurt (2020) without the division
   Yphi=Yphi+(ii)*(-1)^(ii)*(betan_ev(ii)*(1+2*alphn_ev(ii+1))-betan_ev(ii+1)*(1+2*alphn_ev(ii))); 
end

Yphi=Yphi.*8/x^2; % Division operation added in completes Eq. (21)

%% Substitution: Substitute any x=k*a value to evaluate force function

aa=subs(Yphi,x,xeval);
Yst=-double(aa);

end