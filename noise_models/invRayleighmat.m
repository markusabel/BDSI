function noisex = invRayleighmat(u,sig);
%inversion of standard 1D Rayleigh distribution
%integration in Maple

%u is matrix over [0,1]
%valid for matrix input u, matrix output noisex

%uses inverse of CDF, which is monotonic and varies between 0 and 1
%inversion concept explained in 
%https://www.comsol.com/blogs/sampling-random-numbers-from-probability-distribution-functions

%consts
%sig= sigma is a shift parameter

%pdf is
%f  = piecewise(x < 0, 0, 0 <= x, x*exp(-x^2/(2*sigma^2))/sigma^2)

%CDF is
%F = piecewise(x <= 0, 0, 0 < x, -exp(-x^2/(2*sigma^2)) + 1)

%inverting the CDF to get x=invF(u) gives
%x = x = sqrt(-2*ln(-y + 1))*sigma
%seems to handle it without needing cutoff in y

noisex = sqrt(-2*log(-u + 1))*sig


