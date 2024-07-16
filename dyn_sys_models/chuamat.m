function xprime = chuamat(t,x,p);
%Computes the derivatives in the Chua equations; Lynch p287

a=p(1); b=p(2); c=p(3); d=p(4);
g=c*x(:,1) + (1/2)*(d-c)*( abs(x(:,1)+1) - abs(x(:,1)-1) );
xprime=[a*(x(:,2)-x(:,1)-g), x(:,1)-x(:,2)+x(:,3), -b*x(:,2)];


end
