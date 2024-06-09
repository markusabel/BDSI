%outputs from JMAP; trying to standardise
%Written by Robert Niven, 2023-2024

        numintit_JMAP=numintit
        
        %all iterations
        aev_JMAP=aev
        bev_JMAP=bev
        afv_JMAP=afv    
        bfv_JMAP=bfv    
        
        psie_JMAP=psie
        modee_JMAP=modee
        omegae_JMAP=omegae
        psif_JMAP=psif
        modef_JMAP=modef
        omegaf_JMAP=omegaf

        %crit0_JMAP=crit0

        mean_ve_JMAP=mean_ve
        std_ve_JMAP=std_ve

        mean_vf_JMAP=mean_vf
        std_vf_JMAP=std_vf

        true_Xi
        Xi_JMAP=Xi
        diffXi_JMAP=diffXi

        lambdamat_JMAP=lambdamat
        lambda_JMAP=lambda

        Resid2_JMAP=Resid2    %2norms to compare to SINDy
        Reg2_JMAP=Reg2     
        Jtilde2_JMAP=Jtilde2
        %problematic; cannot use 2-norm unless Gaussian isotropic

        ResidGnorm_JMAP=ResidGnorm      %Gaussian norms
        RegGnorm_JMAP=RegGnorm   
        %JGnorm_JMAP=JGnorm   

        logdetve_JMAP=logdetve
        logdetvf_JMAP=logdetvf

        Resid_JMAP=Resid        %log full Gaussians
        Reg_JMAP=Reg   
        %J_JMAP=J 

        Vf_JMAP=Vf                  %For Lorenz, is 9 x 3 variance matrix on parameters in Xi
        Vfsqrt_JMAP=Vfsqrt

        logdetDelta_JMAP = logdetDelta    %posterior terms
        PostGnorm_JMAP = PostGnorm                            
        logPost_JMAP = logPost            %log full Gaussian

        logdetOmega_JMAP=logdetOmega          %evidence terms
        EvidGnorm_JMAP=EvidGnorm                                  
        logEvid_JMAP=logEvid               %log full Gaussian                   

        EvidGnormalt_JMAP = EvidGnormalt
        logEvidalt_JMAP = logEvidalt

        Gnormcheck_JMAP=Gnormcheck
        logBayescheck_JMAP=logBayescheck

        aic_JMAP=aic              %metrics
        bic_JMAP=bic
        aicc_err_JMAP=aicc_err            
        aicc_JMAP=aicc             

        aic2norm_JMAP=aic2norm
        bic2norm_JMAP=bic2norm
        aicc2norm_JMAP=aicc2norm

        EB_JMAP=EB

        %Title_text='JMAP'; 
        Fighead_text='Bayes iterations: JMAP';
        Comp_name='Comp_JMAP';
        


