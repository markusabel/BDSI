function noisex = invErlangmat(u,s,lambda);
%inversion of standard 1D Erlang distribution
%integration in Maple / uses Matlab inverse function

%u is matrix over [0,1]
%valid for matrix input u, matrix output noisex

%uses inverse of CDF, which is monotonic and varies between 0 and 1
%inversion concept explained in 
%https://www.comsol.com/blogs/sampling-random-numbers-from-probability-distribution-functions

%consts
%s = k = parameter of the incomplete gamma distribution
%lambda = parameter of the Erlang distribution

%pdf is
%f  = piecewise(x < 0, 0, 0 <= x, lambda^k*x^(k - 1)*exp(-lambda*x)/(k -1)!); 

%CDF is
%F = gamma(k, lambda*x)/(k - 1)!;   %standard Wiki notation
% = (GAMMA(k) - GAMMA(k, lambda*x))/(k - 1)!;  %in Maple notation

%inverting the CDF to get x=invF(u) is tricky, since the incomplete gamma
%function is defined differently in Maple and Matlab; also with terms
%reversed!
% 
% I obtain
%x = gammaincinv(y*(k - 1)!/Gamma(k), k)/lambda
%seems to not need cutoff in y

noisex = gammaincinv(u*factorial(s - 1)/gamma(s), s)/lambda;




