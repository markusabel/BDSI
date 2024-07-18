function noisex = invGaussianmat(u,mu,var);
%inversion of standard 1D Gaussian distribution
%integration in Maple; uses inverse error function of MATLAB

%u is matrix over [0,1]
%valid for matrix input u, matrix output noisex

%uses inverse of CDF, which is monotonic and varies between 0 and 1
%inversion concept explained in 
%https://www.comsol.com/blogs/sampling-random-numbers-from-probability-distribution-functions

%consts
sig=sqrt(var);

%pdf is
%f  = sqrt(2)*exp(-(x - mu)^2/(2*var))/(2*sqrt(Pi*var))
%   =exp(-(x - mu)^2/(2*sigma^2))/sqrt(2*Pi*sigma^2)
%with var = sigma^2

%CDF is
%F = 1/2 - erf(sqrt(2)*(-x + mu)/(2*sigma))/2

%inverting the CDF to get x=invF(u) gives
%x = -erfinv(-2*u + 1)*sqrt(2)*sigma + mu

noisex=-erfinv(-2*u + 1)*sqrt(2)*sig + mu; 


