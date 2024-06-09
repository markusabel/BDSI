function [fh,vf,ve,Residcol,Regcol,HR,Delta,niter]=VBA6(g,H,ae0,be0,af0,bf0,ntol);
% Author: Ali Mohammad-Djafari
% Date: April 2015
% Follows theoretical algorithm of A-MD and M Dumitru, in 2015 paper
% Corrected internal loop + additional exports + annotations added by
%   Robert Niven, 2021-2014
% Also corrected the ve, vf terms in the algorithm
% Also uses logdet function of Lin 2008
%
% Notation
% g=H * f + e    equiv to    y = Theta * Xi + epsilon
%
% m = no of data points in time series
% n = dimensions of indep parameters 
% o = dimensions of dep parameters (for dyn sys = n); must iterate over j in {1...o}
% c = variables in alphabet
%
% Formulation
% g = Y_j or dX_j = dependent variable(s), jth column from m x o (for dyn sys, jth column from m x n)
% H = Theta = alphabet applied to data, m x c 
% f = Xi_j = param matrix, jth column from c x o
% e = epsilon_j = error, jth column from m x o
%
% Calc quantities
% fh= \hat{Xi}_j = predicted param matrix, jth column from c x o
% ve = diag(V_epsilon_j) = diag of variance matrix of e, m x 1, valid for jth column of Y
% vf = diag(V_Xi_j) = diag of variance matrix of f, c x 1, valid for jth column of Y
% Residcol = (y_j - Theta*\hat{Xi}_j)' * iVe * (y_j - Theta*\hat{Xi}_j) = Gaussian norm for likelihood
% Regcol = \hat{Xi}_j ' *iVf * \hat{Xi}_j = Gaussian norm for prior
% HR = (Delta)^{-1} = inverse of covariance matrix of posterior, valid for jth column of Y
% Delta = covariance matrix of posterior, valid for jth column of Y
%
% VBA: approximate p(f,ve,vf|g)=p(f,ve,vf|g)=p(g|f)p(f|vf)p(ve)p(vf)/p(g)
% by q(f,ve,vf)=q_1(f) q_2(ve) q_3(vf)
% q_1(f)=N(f|fh,Delta), Delta=inv(H'*iVe*H+iVf); fh=Delta*(H'*iVe*g);
% q_2(ve)IG(ve|ae,ve),  ae=ae0+.5; be=be0+.5*(dg.*dg+hth.*sigma);
% q_3(vf)IG(vf|af,vf),  af=af0+.5; bf=bf0+.5*(fh.*fh+sigma);
%
 % Initialization
 [C,N]=size(H);%HtH=H'*H; %hth=diag(HtH);

 %initial values are scalars
 ve0=be0/(ae0-1);  %corrected
 vf0=bf0/(af0-1);  %corrected
 lambda=ve0/vf0;

 Delta=inv(H'*H+lambda*eye(N,N));
 %sigma=diag(Delta);  %not needed here
 fh=Delta*(H'*g);

 %Internal iterations
 niter=0;
 HR=eye(C,C);
 logdetHR=logdet(HR);     %logdet function, Lin 2008, uses LU factorization
 calctol=1E10;   %initial value
 
 while calctol>ntol
 %for iter=1:niter
 %'internal iteration'
  niter=niter+1    %these are true iteration numbers

  gh=H*fh;dg=g-gh;

  %pause;

  %then be, bf turned into vectors; ae, af kept as scalars
  ae=ae0+0.5;
  be=be0+0.5*(dg.*dg+diag(H*Delta*H'));
  %be=be0+.5*(dg.*dg+trace(H*Delta*H'));  %a/c to paper
  %ve=be./(ae-1);
  ve=be./(ae);     %corrected July 2023
  iVe=diag(1./ve);

  af=af0+0.5;
  bf=bf0+0.5*(fh.*fh+diag(Delta));
  %bf=bf0+.5*(fh.*fh+trace(Delta));  %a/c to paper
  %vf=bf./(af-1);
  vf=bf./(af);     %corrected July 2023
  iVf=diag(1./vf);

  Delta=inv(H'*iVe*H+iVf);
  fh=Delta*(H'*iVe*g);

  %f1h=abs(fh);fhmax=max(f1h);tsh=.1*fhmax;tih=find(abs(fh)<tsh);fh(tih)=0;

  logdetoldHR=logdetHR;
  logdetHR=logdet(HR);
  calctol= abs((logdetHR -logdetoldHR)/ logdetoldHR);
 end
 
 %Delta
 HR=inv(Delta);

 %Gaussian norms
 Kern  = g - H*fh;
 Residcol = Kern'*iVe*Kern;
 Regcol   = fh'*iVf*fh;
 
 
