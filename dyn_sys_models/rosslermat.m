function xprime = rosslermat(t,x,p);
%Computes the derivatives for the Rossler equations; Lynch p280; ; Abel matrix form
a=p(1); b=p(2); c=p(3);
xprime=[-x(:,2)-x(:,3), x(:,1)+a*x(:,2), b+x(:,1).*x(:,3)-c*x(:,3)];


end
