function xprime = Vance(t,x,p);
%Computes the rates for a system of 2 prey and 1 predator
%(mod Lotka-Volterra equations)
%reported by Vance 1978, with lumped K terms
%reported as chaotic by Gilpin 1979


r(1)=p(1); r(2)=p(2); r(3)=p(3); 
alpha(1,1)=p(4); alpha(1,2)=p(5); alpha(1,3)=p(6); 
alpha(2,1)=p(7); alpha(2,2)=p(8); alpha(2,3)=p(9); 
alpha(3,1)=p(10); alpha(3,2)=p(11); alpha(3,3)=p(12);

%xprime=[r(1)*x(1) - x(1)*( alpha(1,1)*x(1) + alpha(1,2)*x(2) + alpha(1,3)*x(3) ); ... 
%        r(2)*x(2) - x(2)*( alpha(2,1)*x(1) + alpha(2,2)*x(2) + alpha(2,3)*x(3) ); ...
%        r(3)*x(3) - x(3)*( alpha(3,1)*x(1) + alpha(3,2)*x(2) + alpha(3,3)*x(3) )];

xprime=[r(1)*x(1) - alpha(1,1)*x(1)*x(1) - alpha(1,2)*x(1)*x(2) - alpha(1,3)*x(1)*x(3); ... 
        r(2)*x(2) - alpha(2,1)*x(2)*x(1) - alpha(2,2)*x(2)*x(2) - alpha(2,3)*x(2)*x(3); ...
        r(3)*x(3) - alpha(3,1)*x(3)*x(1) - alpha(3,2)*x(3)*x(2) - alpha(3,3)*x(3)*x(3) ];


end
