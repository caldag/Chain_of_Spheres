%This script calculates the radiation force function Yphi based on Hasegawa's
%1977 study.

%The code also includes Maidanik (1957) formulation (commented out) that
%should give the same result.

function Y_phi=hasegawa77_tr_calc(ka,MASSsc,LENGTHsc,TIMEsc)

%Nickel
rhop=8900*LENGTHsc^2/(MASSsc)*TIMEsc;; % Particle density
c1=5639/(LENGTHsc*TIMEsc); % Compression wave velocity of particle
c2=2970/(LENGTHsc*TIMEsc); % Shear wave velocity of particle

rholiq=1000*LENGTHsc^2/(MASSsc)*TIMEsc; % Liquid density, water
cliq=1480/(LENGTHsc*TIMEsc); % Liquid speed of sound, water

% rholiq=1264*LENGTHsc^2/(MASSsc)*TIMEsc; % Liquid density, glycerine
% cliq=1962/(LENGTHsc*TIMEsc); % Liquid speed of sound, glycerine

x=sym('x'); % x=k*a
x1=x*cliq/c1; 
x2=x*cliq/c2;
n=sym('n'); %Bessel order
n1=sym('n1'); %Bessel order, n+1

%% Bessel functions
jn=sphbes1(n,x);
jn1=sphbes1(n1,x);
djn=(2.^(1/2).*pi.^(1/2).*besselj(n - 1/2, x).*(1./x).^(1./2))./2 - (2.^(1./2).*pi.^(1./2).*besselj(n + 1/2, x).*(n + 1).*(1./x).^(3./2))./2;
ddjn=(2.^(1./2).*pi.^(1/2).*besselj(n + 1/2, x).*(1./x).^(1./2).*(n.^2 + 3.*n - x.^2 + 2))./(2.*x.^2) - 2.^(1./2).*pi.^(1/2).*besselj(n - 1./2, x).*(1./x).^(3./2);

nn=sphbes2(n,x);
nn1=sphbes2(n1,x);
dnn=(2.^(1./2).*pi.^(1./2).*bessely(n - 1./2, x).*(1./x).^(1./2))./2 - (2.^(1./2).*pi.^(1./2).*bessely(n + 1./2, x).*(n + 1).*(1./x).^(3./2))./2;

%% Key variables
Fn=Fcalc(jn,jn1,djn,ddjn,x,x1,x2,rhop,rholiq,n);

[alphn,betan]=alphbetacalc(Fn,jn,jn1,djn,nn,nn1,dnn,x,n);

for ii=0:20
   alphn_ev(ii+1)=subs(alphn,{n n1},{ii ii+1});
   betan_ev(ii+1)=subs(betan,{n n1},{ii ii+1});
end


Yphi=0;

%% Haesgawa (1977), ARF under travelling wave

for ii=1:20
   Yphi=Yphi+(ii)*(alphn_ev(ii)+alphn_ev(ii+1)+2*alphn_ev(ii)*alphn_ev(ii+1)+2*betan_ev(ii)*betan_ev(ii+1)); 
end

Yphi=-Yphi.*4/x^2;

%% Maidanik (1957), ARF under travelling wave. Comment out the Hasegawa part above before running.

% for ii=1:20
%    Yphi=Yphi+(2*(ii-1)+1)*alphn_ev(ii)+2*(ii)*(alphn_ev(ii)*alphn_ev(ii+1)+betan_ev(ii)*betan_ev(ii+1)); 
% end
% 
% Yphi=-Yphi.*4/x^2;


%% Substitution: Substitute any x=k*a value to evaluate ARF
aa=subs(Yphi,x,ka);
Y_phi=double(aa);

end