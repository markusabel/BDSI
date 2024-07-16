function xprime = predatorprey(t,x,p);
%Computes the rates for a system of predator prey model
%with x(1)=N, x(2)=P

r=p(1); K=p(2); kappa=p(3); lambda=p(4);
xprime=[r*x(1)*(1-x(1)/K)-kappa*x(2); lambda*x(1)];


end
