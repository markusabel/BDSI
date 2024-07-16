function xprime = photosyn3spec(t,x,p);
%Computes the photosyn dyn system of Zupanonic + Sarah Hall
%4 variables X, Y, Z
%New stoich form 
%1st order coefficients

%format long e
numrn=8;
numspec=3;
%p=[1 2 3 4 1/5 2/5 3/5 4/5]                %test values
%x=[10,100,200]      %test values

k=p;

X=[x(1) x(2) x(3)];             %dummy var

%kinetic power matrix (reactants only)
KT=[1 1 0; 1 1 0; 1 1 0; 1 1 0; 0 1 1; 0 1 1; 1 0 1; 1 0 1;];               

%rate vector substituted
for i = 1:numrn
   Q(i,1)=k(i)*prod(X.^KT(i,:));
end
%makes column vector of rate terms, one for each rate eq

%stoichiometric matrix
ST=[-1 1 0; 1 -1 0; 1 -1 0; -1 1 0; 0 -1 1; 0 1 -1; 1 0 -1; -1 0 1;]; 
S=ST.';

xprime=S*Q;

end
