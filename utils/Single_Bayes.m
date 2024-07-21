%% Bayesian method - single calculation
% Written by Mohammad-Djafari & Dumitru 2015
% Modifications by Robert Niven, 2023-2024

%% Single iteration
format shortEng

%now setup for aev, bev, afv, bfv as o-vectors, so can pass output from
%optBayes as well as a single Bayes analysis

if iterateBayes==0
        %i.e. single iteration
        %keep ae0; be0; af0; bf0; but preserve dimension
        aev=ae0*ones(1,o); 
        bev=be0*ones(1,o); 
        afv=af0*ones(1,o); 
        bfv=bf0*ones(1,o);
elseif iterateBayes>=10
        %use optima
        aev=aeopt; 
        bev=beopt; 
        afv=afopt;
        bfv=bfopt;
end

%parameters

    lambda(1:o)=NaN          %allows for output dimension o

    %inverse gamma criteria; over dimensions j
    %careful: aev, bev etc are o-vectors here
    psie=bev./(aev-1);  
    modee=bev./(aev+1);  
    omegae=bev.*bev/((aev-1).*(aev-1).*(aev-2)); 
    
    psif=bfv./(afv-1);  
    modef=bfv./(afv+1);  
    omegaf=bfv.*bfv./((afv-1).*(afv-1).*(afv-2));  

for j=1:o  %loop over output dimensions o
     %internal loop variables fh, vf, ve, Sigmaf, sigmaf 
     'Dimension *****' 
     j
     g=dx(:,j);

     if Bayesmeth==1
         if jmap=='JMAP6' 
             [fh,vf,ve,Residcol,Regcol,HR,Delta,niter]=...
                     JMAP6(g,H,aev(j),bev(j),afv(j),bfv(j),ntol);
         else
         %   do nothing
         end;
     elseif Bayesmeth==2
         if vba=='VBA6'
            [fh,vf,ve,Residcol,Regcol,HR,Delta,niter]=...
                    VBA6(g,H,aev(j),bev(j),afv(j),bfv(j),ntol);
         else
         %   do nothing
         end;
     end;

     numintit(j) = niter;

     %calc of variables ******

         %construct for jth column and kth iteration
         Xi(:,j)=fh;
         if knownXi=='Y'       %calc data - known Xi
            diffXi(:,j)=true_Xi(:,j) - Xi(:,j);    
         elseif knownXi=='N'    %unknown Xi
            diffXi(:,j) = NaN;          %protect against mismatch for posterior calc; will mess with graphs
            %diffXi=Xi(:,i)-Xi(:,i);    %to show errors; need to implement
         end;
     
         %'Calc ve, mean(ve), std(ve)' 
         %ve
         mean_ve(j)=mean(ve);
         std_ve(j)=std(ve);

         %'Calc vf, mean(vf), std(vf)' 
         %vf
         mean_vf(j)=mean(vf);
         std_vf(j)=std(vf);
     
         lambdamat(:,j)=mean_ve(j)*ones(cc,1)./vf;  %c-vector over alphabet, jth column
         lambda(j)=min(lambdamat(:,j));  %single value for each dimension

         %2norm quantities
         Resid2(j)= (norm(dx(:,j) - H*Xi(:,j)))^2;      %2norm for likelihood
         Reg2(j)= (norm(Xi(:,j)))^2;                    %2norm for prior
         Jtilde2(j)=Resid2(j) + lambda(j) * Reg2(j);   %2norm objective function incl lambda, each j
         %problematic; 2-norm not representative unless Gaussian isotropic
     
         %Gen Gaussian norm quantities - from inside JMAP5 or VBA5
         ResidGnorm(j)=Residcol;           %Gaussian norm for likelihood
         RegGnorm(j)=Regcol;               %Gaussian norm for prior
         %JGnorm(j)=Residcol+Regcol ;       %sum of norms, not much meaning
     
         %Gen Gaussian norm quantities - again
         %go direct to det using log simplification, assuming diag matrices
         %(assumption is correct - these are constructed as diagonal)
         logve=log(ve);
         logvf=log(vf);
         logdetve(j)=sum(logve);
         logdetvf(j)=sum(logvf);
         Resid(j)=-1/2*Residcol - m/2*log(2*pi) - 1/2*logdetve(j);  %log likelihood exact
         Reg(j)=-1/2*Regcol - cc/2*log(2*pi) - 1/2*logdetvf(j);  %log prior exact
         %J(j)=Resid(j) + Reg(j);  %sum of log likelihood + log prior, does not include log evidence

         %Posterior covariance calc
         %both HR, Delta are c x c, fairly small
         HR                           %inverse of posterior covariance matrix, jth column, discard each version
         Delta                        %posterior covariance matrix = Delta_i, jth column, discard each version
         Vf(:,j)=diag(Delta)          %diag of posterior covariance matrix, jth column, r
         Vfsqrt(:,j)=sqrt(Vf(:,j))          %std devs from posterior covariance matrix, jth column, retain columns       

         %need some ranking of non-diag of covariance matrix

         %log Posterior direct calc
         logdetDelta(j)=-logdet(HR)                    %logdet function, Lin 2008, uses LU factorization
                                                         %applied to HR = inverse Delta, then converted, to reduce ill-cond problems
         PostGnorm(j)=diffXi(:,j)'*HR*diffXi(:,j)     %Gaussian norm for posterior 
         logPost(j)=-1/2*PostGnorm(j) - cc/2*log(2*pi) ...
             - 1/2*logdetDelta(j);  %log posterior exact

         PostGnormcheck(j)=PostGnorm(j)-fval(j)     %check on Gaussian norm for posterior; only for optBayes

         %Evidence quantities
         estY(:,j)=  zeros(m,1);   % H*mu_prior(:,j);    %evidence mean, jth column, retain columns; assumes zero prior mean 
         Omega=H*diag(vf)*H'+diag(ve);                   %evidence cov matrix, jth column, discard each version   
         invOmega = inv(Omega);                          %inverse of evidence cov matrix, jth column, discard each version
         logdetOmega(j)=-logdet(invOmega)              %logdet function, Lin 2008, uses LU factorization
                                                         %applied to inverse Omega, then conv, to reduce ill-cond problems
         diffY(:,j)=dx(:,j) - estY(:,j);                 %for dynamical system, dx is the same as Y
         EvidGnorm(j)=diffY(:,j)'*invOmega*diffY(:,j)       %Gaussian norm for evidence      
         logEvid(j)=-1/2*EvidGnorm(j) - m/2*log(2*pi) ...
             - 1/2*logdetOmega(j);  %log evidence exact
  
         %alternative evidence calc from other terms
         EvidGnormalt(j)= ResidGnorm(j) + RegGnorm(j) - PostGnorm(j);   %Gnorms
         logEvidalt(j)= Resid(j) +  Reg(j) - logPost(j); %logs of full Gaussians

         %Bayesian check mult and sum
         Gnormcheck(j) = ResidGnorm(j) + RegGnorm(j) - PostGnorm(j) ...
             - EvidGnorm(j)   %Gnorms, will not be normalised
         logBayescheck(j) = Resid(j) + Reg(j) - logPost(j) - logEvid(j)     %logs of full Gaussians; should be 0
 
         %other metrics
         aic(j)= -2*Resid(j)+2*cc;                       %Bishop eq 1.73 mod with -2* factor
         bic(j)= -2*Resid(j)+cc*log(m);                  %Bishop eq 4.139 mod with -2* factor

         %corrected AIC; note Mangan etal 2017 uses incorrect correction in code
         aicc_err= 2*(cc+1)*(cc+2)/(m-cc-2);     %error in correction
         aicc(j)= aic(j) + aicc_err;             %correction

         %from Mangan etal 2017, with BIC corrected
         %note SINDy_BIC is incorrect in code (extra factor of 2 on cc*log(m) term)
         %also Chen & Li 2021
         logL2(j)=(-m/2*log(Resid2(j)/m));
         aic2norm(j)= -2*logL2(j)   +2*cc;
         bic2norm(j)= -2*logL2(j)   +cc*log(m);

         aicc2norm(j)= aic2norm(j) + aicc_err;    %corrected

         %Zhang & Lin "error bar"
         CV=Vf(:,j)./(Xi(:,j).*Xi(:,j));
         EB(j)=sum(CV);

     'end of Dimension *****' 

end

    
%Outputs: export of params and graphs

    if Bayesmeth==1
        outputs_JMAP

    elseif Bayesmeth==2
        outputs_VBA
        
    end



    





