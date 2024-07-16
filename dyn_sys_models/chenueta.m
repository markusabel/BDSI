function xprime = chenueta(t,x,p);
%Computes the derivatives in the Chen-Ueta system, Lynch p 285.

sig=p(1); beta=p(2); rho=p(3);
xprime=[-sig*x(1) + sig*x(2); (rho-sig)*x(1) + rho*x(2) - x(1)*x(3); -beta*x(3) + x(1)*x(2)];


end

%appears that xdot, zdot same as Lorenz; only distinction is in ydot
