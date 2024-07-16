function xprime = LotkaVolterra(t,x,p);
%Computes the rates for a system of 2 competing herbivores
%(competitive Lotka-Volterra equations)

r(1)=p(1); r(2)=p(2); K(1)=p(3); K(2)=p(4); alpha=p(5); beta=p(6);
xprime=[r(1)*x(1)*(1-( x(1) + alpha*x(2) )/K(1) )  ; r(2)*x(2)*(1-( x(2) + beta*x(1) )/K(2) )];


end
