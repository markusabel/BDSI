function noisex = invLaplacemat(u,mu,var);
%inversion of standard 1D Laplace distribution
%integration in Maple

%u is matrix over [0,1]
%valid for matrix input u, matrix output noisex

%uses inverse of CDF, which is monotonic and varies between 0 and 1
%inversion concept explained in 
%https://www.comsol.com/blogs/sampling-random-numbers-from-probability-distribution-functions

%consts
sig=sqrt(var);

%pdf is
%f= exp(-sqrt(2)*abs(mu - x)/sqrt(var)) / (sqrt(2*var));

%CDF is
%F = piecewise(x <= mu, exp(-sqrt(2)*(mu - x)/sigma)/2, mu < x, -exp(sqrt(2)*(mu - x)/sigma)/2 + 1)

%inverting the CDF to get x=invF(u) gives
%piecewise(y <= 1/2, ((sqrt(2)*mu + ln(2*y)*sigma)*sqrt(2))/2, 1/2 < y, ((sqrt(2)*mu - ln(2 - 2*y)*sigma)*sqrt(2))/2)

%loops
[mm nn]=size(u); % find row and col number
%noisex=zeros(size(u));   %in main code

for i=1:mm
    for j=1:nn
        if u(i,j)<= 1/2
            noisex(i,j)=(((sqrt(2)*mu + log(2*u(i,j))*sig)*sqrt(2))/2);
        elseif u(i,j)> 1/2
            noisex(i,j)=(((sqrt(2)*mu - log(2-2*u(i,j)*sig)*sqrt(2))/2));
        end
    end
end

