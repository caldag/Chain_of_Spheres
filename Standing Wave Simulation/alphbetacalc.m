% For the commented version, see the folder Force Calculation Only/Nickel/

function [alphn, betan] = alphbetacalc(Fn,jn,jn1,djn,nn,nn1,dnn,x,n)

Gn=(Fn-n).*jn+x.*jn1;
Hn=(Fn-n).*nn+x.*nn1;

alphn=-Gn.^2/(Gn.^2+Hn.^2);
betan=-Gn.*Hn./(Gn.^2+Hn.^2);
end