function xprime = FeIIFeIII(t,x,p);
%Computes the Fe(II)/Fe(III) dyn system of Ord etal 2012 eq. 6.2
%4 variables P, FeII, FeIII, C
%Old non-stoich form 

k=p;  %params = rate consts k1, k2, k3, k4 
xprime=[-k(1)*x(1); k(1)*x(1)-k(2)*x(2)-k(3)*x(2)*x(3)^2; ...
    k(2)*x(2)+k(3)*x(2)*x(3)^2-k(4)*x(3); k(4)*x(3)];

end
