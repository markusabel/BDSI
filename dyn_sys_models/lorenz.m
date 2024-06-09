function xprime = lorenz(t,x,p);
%Computes the coordinates in the Lorenz equations.

sig=p(1); bet=p(2); rh=p(3);
xprime=[-sig*x(1) + sig*x(2); rh*x(1) - x(2) - x(1)*x(3); -bet*x(3) + x(1)*x(2)];

end
