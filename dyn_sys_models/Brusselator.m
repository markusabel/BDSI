function xprime = Brusselator(t,x,p);
%Computes the Brusselator chemical chaotic system 
%2 variables X, Y; 2 parameters A, B,

A=p(1); B=p(2); %C=p(3); D:=p(4);  %don't need params D, E, in model
xprime=[A+x(1)^2*x(2)-B*x(1)-x(1)  ; B*x(1)-x(1)^2*x(2)];

end
