function xprime = galerkin3D(t,x,p);
%Computes the derivatives in the 3D Galerkin system given by M Schlegel.
%(Noack 2003)

sigma=p(1); mu=p(2);  rho=p(3);
xprime=[sigma*x(1)-x(2)-x(3)*x(1);mu*x(2)+x(1)-x(3)*x(2);-rho*x(3)+x(1)^2+x(2)^2];


end
