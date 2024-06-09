%outputs from VBA; trying to standardise
%Written by Robert Niven, 2023-2024

        numintit_VBA=numintit
        
        %all iterations
        aev_VBA=aev
        bev_VBA=bev
        afv_VBA=afv    
        bfv_VBA=bfv    
        
        psie_VBA=psie
        modee_VBA=modee
        omegae_VBA=omegae
        psif_VBA=psif
        modef_VBA=modef
        omegaf_VBA=omegaf

        %crit0_VBA=crit0

        mean_ve_VBA=mean_ve
        std_ve_VBA=std_ve

        mean_vf_VBA=mean_vf
        std_vf_VBA=std_vf

        true_Xi
        Xi_VBA=Xi
        diffXi_VBA=diffXi

        lambdamat_VBA=lambdamat
        lambda_VBA=lambda

        Resid2_VBA=Resid2    %2norms to compare to SINDy
        Reg2_VBA=Reg2     
        Jtilde2_VBA=Jtilde2
        %problematic; cannot use 2-norm unless Gaussian isotropic

        ResidGnorm_VBA=ResidGnorm      %Gaussian norms
        RegGnorm_VBA=RegGnorm   
        %JGnorm_VBA=JGnorm   

        logdetve_VBA=logdetve
        logdetvf_VBA=logdetvf

        Resid_VBA=Resid        %log full Gaussians
        Reg_VBA=Reg   
        %J_VBA=J 

        Vf_VBA=Vf                  %For Lorenz, is 9 x 3 variance matrix on parameters in Xi
        Vfsqrt_VBA=Vfsqrt

        logdetDelta_VBA = logdetDelta    %posterior terms
        PostGnorm_VBA = PostGnorm                            
        logPost_VBA = logPost            %log full Gaussian

        logdetOmega_VBA=logdetOmega          %evidence terms
        EvidGnorm_VBA=EvidGnorm                                  
        logEvid_VBA=logEvid               %log full Gaussian                   

        EvidGnormalt_VBA = EvidGnormalt
        logEvidalt_VBA = logEvidalt

        Gnormcheck_VBA=Gnormcheck
        logBayescheck_VBA=logBayescheck

        aic_VBA=aic              %metrics
        bic_VBA=bic
        aicc_err_VBA=aicc_err            
        aicc_VBA=aicc             

        aic2norm_VBA=aic2norm
        bic2norm_VBA=bic2norm
        aicc2norm_VBA=aicc2norm

        EB_VBA=EB

        %Title_text='VBA'; 
        Fighead_text='Bayes iterations: VBA';
        Comp_name='Comp_VBA';
        


