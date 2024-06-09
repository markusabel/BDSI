function [fh,vf,ve,Residcol,Regcol,HR,Delta,niter]=JMAP6(g,H,ae0,be0,af0,bf0,ntol);
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
% p(g|f)=N(f|H f,ve I), p(ve)=IG(ve|ae0,be0)
% p(f|vf)=N(f|0,vf I),  p(vf)=IG(ve|ae0,be0)
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
% JMAP: maximizes p(f,ve,vf|g)=p(g|f)p(f|vf)p(ve)p(vf)/p(g)
% with respect to f, ve and vf 
% Must be applied columnwise to multivariate data
% 
 % Initialization
 [C,N]=size(H); 
 HtH=H'*H;   %need to calc inside fn
 %Htg=H'*g;  %not used

 %initial values are scalars
 ve0=be0/(ae0-1);  %corrected
 vf0=bf0/(af0-1);  %corrected
 lambda=ve0/vf0;
 fh=(HtH+lambda*eye(N,N))\(H'*g);

 %Internal iterations
 niter=0;
 HR=eye(C,C);
 logdetHR=logdet(HR);     %logdet function, Lin 2008, uses LU factorization
 calctol=1E10;   %initial value

 while calctol>ntol
 %for iter=1:niter
  %'internal iteration'
  niter=niter+1    %these are true iteration numbers
  
  gh=H*fh;  dg=g-gh;

  %then here be, bf turned into vectors; ae, af kept as scalars
  ae=ae0+1.5;  %corrected
  be=be0+0.5*(dg.*dg);  %same size as dg (= m x 1)
  ve=be./ae;    %not a matrix; same size as dg (= m x 1)
  iVe=diag(1./ve);  %matrix m x m

  af=af0+1.5;   %corrected
  bf=bf0+0.5*(fh.*fh);  %same size as fh (= c x 1)
  vf=bf./af;    %not a matrix; same size as fh (= c x 1)
  iVf=diag(1./vf);   %matrix c x c

  %pause;

  %inverse Delta = HR
  HR=(H'*iVe*H+iVf);
  fh=HR\(H'*iVe*g);    %matrix left division A\B = (A inverse) * B
  %fh=lsqr(HR,Htg,[],50,[],[],fh);
  
  %f1h=abs(fh);fhmax=max(f1h);tsh=.1*fhmax;tih=find(abs(fh)<tsh);fh(tih)=0;

  logdetoldHR=logdetHR;
  logdetHR=logdet(HR);
  calctol= abs((logdetHR -logdetoldHR)/ logdetoldHR);
 end
 
 Delta=inv(HR);
 %lambda=ve/vf;
 
 %Gaussian norms
 Kern  = g - H*fh;
 Residcol = Kern'*iVe*Kern;    %Gaussian norm for likelihood
 Regcol   = fh'*iVf*fh;        %Gaussian norm for prior

 


