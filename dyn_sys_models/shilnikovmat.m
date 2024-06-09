function xprime = shilnikovmat(t,x,p);
%Computes the coordinates in the equations of Shil'nikov, Shil'nikov and
%Turaev, International Journal of Bifurcation and Chaos, Vol. 3, No.5 (1993) 1123-1139

alf=p(1); lam=p(2); B=p(3);
m=size(x,1); 
won=ones(m,1);
xprime=[x(:,2), ... 
    x(:,1).*(won(:)-x(:,3))-B*(x(:,1).*x(:,1).*x(:,1))-lam*x(:,2), ...
    -alf*(x(:,3)-(x(:,1).*x(:,1))) ];
end
