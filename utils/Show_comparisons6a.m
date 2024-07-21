%% Show_comparisons6a
%% Show and compare the coefficients and generated time series
%% and also the generated derivatives
% Modified R Niven June 2018 to Jul 2023
% With sqrt variance terms (Dec 2018)
% Also with \xi symbols, July 2019

% This is for 3D systems; need other codes for 2D, 4D, 7D or 8D

% can only specify true_Xi, if known to problem
% otherwise true_Xi and diffXi are dummies
% this code may not handle unknown Xi - must check
% clearvars diffXi;  %Clear variables from memory

%if knownXi=='Y';     %calc data - known param values in Xi
%    diffXi=true_Xi-Xi     %known function
%elseif knownXi=='N'
%    diffXi=Xi-Xi          %to show errors; unknown function
%end;
%NOW in Iterate_SINDy, Iterate_Bayes, and just after single SINDy

fighead3=strcat(fighead, ": ", Title_text);

%set number of derivative param subplots
numsub = max(3,o);

%plot of figure to contain parameter values ***
%applied if there is an iterative method ***************
if bayesian==0
    if (SINDy==1) && (iterateSINDy==1) 
        format shortEng
        str8c=['SINDy method'];
        str10=['Iteration k= ', int2str(k), ...
            ', Consts: lambda=', num2str(lambda(k),' %.6g ')]; 
        figure
        annotation('textbox',...
        [0.05 0.05 0.85 0.35],...
        'String',{str0, str0a, str8c, str10},...
        'FontSize',14,...
        'FontName','Helvetica');  
        fignum=fignum+1; 
        saveas(gcf,fullfile(pwd, dirname,strcat('Fig',int2str(fignum))),'epsc')
    elseif (LASSO==1) && (iterateLASSO==1)        
        format shortEng
        str8c=['LASSO method'];
        str10=['Iteration k= ', int2str(k), ...
            ', Consts: kappa=', num2str(kappa(k),' %.6g ')]; 
        figure
        annotation('textbox',...
        [0.05 0.05 0.85 0.35],...
        'String',{str0, str0a, str8c, str10},...
        'FontSize',14,...
        'FontName','Helvetica');  
        fignum=fignum+1; 
        saveas(gcf,fullfile(pwd, dirname,strcat('Fig',int2str(fignum))),'epsc')
     elseif (Ridge==1) && (iterateRidge==1)        
        format shortEng
        str8c=['Ridge regression method'];
        str10=['Iteration k= ', int2str(k), ...
            ', Consts: theta=', num2str(theta(k),' %.6g ')]; 
        figure
        annotation('textbox',...
        [0.05 0.05 0.85 0.35],...
        'String',{str0, str0a, str8c, str10},...
        'FontSize',14,...
        'FontName','Helvetica');  
        fignum=fignum+1; 
        saveas(gcf,fullfile(pwd, dirname,strcat('Fig',int2str(fignum))),'epsc')
   end
elseif bayesian==1
  if iterateBayes>=1
    format shortEng
    
    if Bayesmeth==1
       str8a=['Bayesian method: ', jmap];
    elseif Bayesmeth==2
       str8a=['Bayesian method: ', vba];
    end

    if iterateBayes==1
       str8b=['Iterate Bayes on error'];
    elseif iterateBayes==2
       str8b=['Iterate Bayes on model'];
    elseif iterateBayes==3
       str8b=['Iterate Bayes on error and model (equal)'];     
    elseif iterateBayes==4
       str8b=['Iterate Bayes on error (vary b, fixed a)'];     
    elseif iterateBayes==5
       str8b=['Iterate Bayes on error (vary a, fixed b)'];     
    elseif iterateBayes==10
       str8b=['Optimise Bayes on error (over a_e, b_e)'];     
    end


    if iterateBayes==0 
        str9=['Single Bayes:', ...
            ', Consts for j=1 to: ', int2str(o), ...
            ', a_{e0}=', num2str(ae0,' %.6g '), ...
            ', b_{e0}=', num2str(be0,' %.6g '), ...
            ', a_{f0}=', num2str(af0,' %.6g '), ...
            ', b_{f0}=', num2str(bf0,' %.6g ')];    
    elseif iterateBayes>=1 & iterateBayes<10
        str9=['Iterate Bayes: iteration k= ', int2str(k), ...
            ', Consts: a_{ev}=', num2str(aev(k,:),' %.6g '), ...
            ', b_{ev}=', num2str(bev(k,:),' %.6g '), ...
            ', a_{fv}=', num2str(afv(k,:),' %.6g '), ...
            ', b_{fv}=', num2str(bfv(k,:),' %.6g '), ...
            ', active b/(a-1)=', num2str(crit0(k),' %.6g ')];
    elseif iterateBayes==10
        str9=['Optimised Bayes:',  ...
            ', Consts for j=1 to: ', int2str(o), ...
            ', a_{eopt}=', num2str(aev,' %.6g '), ...
            ', b_{eopt}=', num2str(bev,' %.6g '), ...
            ', a_{f0}=', num2str(afv,' %.6g '), ...
            ', b_{f0}=', num2str(bfv,' %.6g ')];    
    elseif iterateBayes==11
        str9=['Optimised Bayes:',  ...
            ', Consts for j=1 to: ', int2str(o), ...
            ', a_{eopt}=', num2str(aev,' %.6g '), ...
            ', b_{eopt}=', num2str(bev,' %.6g '), ...
            ', a_{fopt}=', num2str(afv,' %.6g '), ...
            ', b_{fopt}=', num2str(bfv,' %.6g ')];    
    end
   
    figure
    annotation('textbox',...
    [0.05 0.05 0.85 0.35],...
    'String',{str0, str0a, str8a, str8b, str9},...
    'FontSize',14,...
    'FontName','Helvetica'); 
    fignum=fignum+1; 
    saveas(gcf,fullfile(pwd, dirname,strcat('Fig',int2str(fignum))),'epsc')
  end
end

% Then plot the coefficients xi_ij ***********************
% for labels see setlabels.m
% for colours see setlabels.m

%limits
for w=1:o
    if(bayesian==1)     
        ximax(w)=max([max(true_Xi(:,w)) max(Xi(:,w)+Vfsqrt(:,w))]);
        ximin(w)=min([min(true_Xi(:,w)) min(Xi(:,w)-Vfsqrt(:,w))]);
    else %not Bayesian
        ximax(w)=max([max(true_Xi(:,w)) max(Xi(:,w))]);
        ximin(w)=min([min(true_Xi(:,w)) min(Xi(:,w))]);
    end
end

% true Xi values, now for any dimension o; fixed y scale
figure 
for v=1:o-1 
    subplot(numsub,1,v),...
        plot(true_Xi(:,v), plotsym{v}, 'Color', plotcol{v}, 'LineWidth',1),...
        xlim([1 cc]),...
        %ylim([ximin(v) ximax(v)]),...
        ylabel(coefflab(v),'Interpreter','latex','FontSize',18),...
        ax = gca; ax.YAxis.FontSize = 16;
        ax.XTick = seqcols; ax.XTickLabel = [];  
end
subplot(numsub,1,o),...
    plot(true_Xi(:,o), plotsym{o}, 'Color', plotcol{o}, 'LineWidth',1),...
    xlabel('Dictionary','FontSize',18),...
    xlim([1 cc]),...
    %ylim([ximin(o) ximax(o)]),...
    ylabel(coefflab(o),'Interpreter','latex','FontSize',18),...
    axc = gca; axc.XAxis.FontSize = xlabsize; axc.YAxis.FontSize = 16;
    axc.XTick = seqcols; axc.XTickLabel = dictlist;  
    axc.XTickLabelRotation = xlabrot; 
figtitle(fighead3, 'fontweight','bold','FontSize',18)  
fignum=fignum+1; 
saveas(gcf,fullfile(pwd, dirname,strcat('Fig',int2str(fignum))),'epsc')

%pause

if(bayesian==1)     
    %Xi values; fixed y scale
    figure  
    for v=1:o-1  
        subplot(numsub,1,v),...
            errorbar(Xi(:,v),Vfsqrt(:,v), plotsym{v}, 'Color', plotcol{v}, 'LineWidth',1),...
            xlim([1 cc]),...
            %ylim([ximin(v) ximax(v)]), ...
            ylabel(coeffhatlab(v),'Interpreter','latex','FontSize',18),...
            ax = gca; ax.YAxis.FontSize = 16;
            ax.XTick = seqcols; ax.XTickLabel = [];  
    end
    subplot(numsub,1,o),...
        errorbar(Xi(:,o),Vfsqrt(:,o), plotsym{o}, 'Color', plotcol{o}, 'LineWidth',1),...
        xlabel('Dictionary','FontSize',18),...
        xlim([1 cc]),...
        %ylim([ximin(o) ximax(o)]), ...
        ylabel(coeffhatlab(o),'Interpreter','latex','FontSize',18),...
        axa = gca; axa.XAxis.FontSize = xlabsize; axa.YAxis.FontSize = 16;
        axa.XTick = seqcols; axa.XTickLabel = dictlist;  
        axa.XTickLabelRotation = xlabrot;  
    figtitle(fighead3, 'fontweight','bold','FontSize',18)  
    fignum=fignum+1; 
    saveas(gcf,fullfile(pwd, dirname,strcat('Fig',int2str(fignum))),'epsc')

    %pause;

    %diffXi values; floating y scale
    if knownXi=='Y'     
        figure 
        for v=1:o-1  
            subplot(numsub,1,v),...
                errorbar(diffXi(:,v),Vfsqrt(:,v), plotsym{v}, 'Color', plotcol{v}, 'LineWidth',1),...
                xlim([1 cc]), ... 
                ylabel(coeffdifflab(v),'Interpreter','latex','FontSize',18),...
                ax = gca; ax.YAxis.FontSize = 16;
                ax.XTick = seqcols; ax.XTickLabel = [];  
        end
        subplot(numsub,1,o),...
            errorbar(diffXi(:,o),Vfsqrt(:,o), plotsym{o}, 'Color', plotcol{o}, 'LineWidth',1),...
            xlabel('Dictionary','FontSize',18),...
            xlim([1 cc]), ...
            ylabel(coeffdifflab(o),'Interpreter','latex','FontSize',18),...
            axb = gca; axb.XAxis.FontSize = xlabsize; axb.YAxis.FontSize = 16; 
            axb.XTickLabelRotation = xlabrot; 
            axb.XTick = seqcols; axb.XTickLabel = dictlist;  
        figtitle(fighead3, 'fontweight','bold','FontSize',18)  
        fignum=fignum+1; 
        saveas(gcf,fullfile(pwd, dirname,strcat('Fig',int2str(fignum))),'epsc')
    end
 
    %centred diffXi values; floating y scale
    figure 
    for v=1:o-1  
        subplot(numsub,1,v),...
            errorbar(Xi(:,v)-Xi(:,v),Vfsqrt(:,v), plotsym{v}, 'Color', plotcol{v}, 'LineWidth',1),...
            xlim([1 cc]),...
            ylabel(vanishlab(v),'Interpreter','latex','FontSize',18),...
            ax = gca; ax.YAxis.FontSize = 16;
            ax.XTick = seqcols; ax.XTickLabel = [];  
    end
    subplot(numsub,1,o),...
        errorbar(Xi(:,o)-Xi(:,o),Vfsqrt(:,o), plotsym{o}, 'Color', plotcol{o}, 'LineWidth',1),...
        xlabel('Dictionary','FontSize',18), xlim([1 cc]),...
        ylabel(vanishlab(o),'Interpreter','latex','FontSize',18),...
        axb = gca; axb.XAxis.FontSize = xlabsize; axb.YAxis.FontSize = 16; 
        axb.XTickLabelRotation = xlabrot; 
        axb.XTick = seqcols; axb.XTickLabel = dictlist;  
    figtitle(fighead3, 'fontweight','bold','FontSize',18)  
    fignum=fignum+1; 
    saveas(gcf,fullfile(pwd, dirname,strcat('Fig',int2str(fignum))),'epsc')
   
        
else  %not bayesian
    %Xi values; fixed y scale
    figure 
    for v=1:o-1  
        subplot(numsub,1,v),...
            plot(Xi(:,v), plotsym{v}, 'Color', plotcol{v}, 'LineWidth',1),...
            xlim([1 cc]), ...
            %ylim([ximin(v) ximax(v)]), ...
            ylabel(coeffhatlab(v),'Interpreter','latex','FontSize',18),...
            ax = gca; ax.YAxis.FontSize = 16;
            ax.XTick = seqcols; ax.XTickLabel = []; 
    end
    subplot(numsub,1,o),...
        plot(Xi(:,o), plotsym{o}, 'Color', plotcol{o}, 'LineWidth',1),...
        xlabel('Dictionary','FontSize',18),...
        xlim([1 cc]), ...
        %ylim([ximin(v) ximax(v)]), ...
        ylabel(coeffhatlab(o),'Interpreter','latex','FontSize',18),...
        axc = gca; axc.XAxis.FontSize = xlabsize; axc.YAxis.FontSize = 16;
        axc.XTick = seqcols; axc.XTickLabel = dictlist;  
        axc.XTickLabelRotation = xlabrot; 
    figtitle(fighead3, 'fontweight','bold','FontSize',18)  
    fignum=fignum+1; 
    saveas(gcf,fullfile(pwd, dirname,strcat('Fig',int2str(fignum))),'epsc')

    %diffXi values; floating y scale   
    if knownXi=='Y'
        figure 
        for v=1:o-1  
            subplot(numsub,1,v),plot(diffXi(:,v), plotsym{v}, 'Color', plotcol{v}, 'LineWidth',1),...
                xlim([1 cc]),...
                ylabel(coeffdifflab(v),'Interpreter','latex','FontSize',18),...
                ax = gca; ax.YAxis.FontSize = 16;
                ax.XTick = seqcols; ax.XTickLabel = []; 
        end
        subplot(numsub,1,o),plot(diffXi(:,o), plotsym{o}, 'Color', plotcol{o}, 'LineWidth',1),...
            xlabel('Dictionary','FontSize',18),xlim([1 cc]),...
            ylabel(coeffdifflab(o),'Interpreter','latex','FontSize',18),...
            axd = gca; axd.XAxis.FontSize = xlabsize; axd.YAxis.FontSize = 16;
            axd.XTick = seqcols; axd.XTickLabel = dictlist; 
            axd.XTickLabelRotation = xlabrot; 
        figtitle(fighead3, 'fontweight','bold','FontSize',18)  
        fignum=fignum+1; 
        saveas(gcf,fullfile(pwd, dirname,strcat('Fig',int2str(fignum))),'epsc')
    end
          
end

%pause


% compare Xi with the original used to simulate data ******************
col_true_Xi=arrayfun(@colorlog,true_Xi);
col_Xi=arrayfun(@colorlog,Xi);

%first fig will plot as blank plot for Xi not known
figure 
subplot(1,2,1),imagesc(col_true_Xi),...
    title('True \Xi','FontSize',16),...
    xlabel('Variable','FontSize',18), ...
    xlim([0.5 o+0.5]), ...
    ylabel('Dictionary','FontSize',18), ...
    ylim([0.5 cc+0.5]),...
    colorbar,...
    ax = gca; ax.XAxis.FontSize = 16; ax.YAxis.FontSize = xlabsize;
    ax.XTick = seqY; ax.XTickLabel = varlist;  
    ax.YTick = seqcols; ax.YTickLabel = dictlist;  
    c=colorbar; c.FontSize = 14;
subplot(1,2,2),imagesc(col_Xi),...
    title('Est. \Xi','FontSize',16),...
    xlabel('Variable','FontSize',18),...
    xlim([0.5 o+0.5]), ...
    ylabel('Dictionary','FontSize',18),...
    ylim([0.5 cc+0.5]),...
    colorbar,...
    ax = gca; ax.XAxis.FontSize = 16; ax.YAxis.FontSize = xlabsize;
    ax.XTick = seqY; ax.XTickLabel = varlist;  
    ax.YTick = seqcols; ax.YTickLabel = dictlist; 
    c=colorbar; c.FontSize = 14;
fignum=fignum+1; 
saveas(gcf,fullfile(pwd, dirname,strcat('Fig',int2str(fignum))),'epsc')

%pause

%Computation of dynamical system from predicted coefficients ***************

%t_h = zeros(m,1);      %protects against breakage in calcs
%x_h = zeros(m,n); 

if isequal(checktype,'sp')     
   
    p

    %assignment of predicted coefficients to original function call
    %uses parameter block p_h
    %this is really clunky but would take effort to automate
    %%better to use checktype 'mat' to replace by generic fn call to 
    %%polynomial function, using entire Xi matrix

    if isequal(eqname,'lorenz') 
        sigma_h=Xi(3-1+inone,1);
        beta_h =-Xi(4-1+inone,3);
        rho_h  =Xi(2-1+inone,2);
        p_h=[sigma_h,beta_h,rho_h]
    elseif isequal(eqname,'rossler')
        if (inone==1)
            b_h=Xi(1-1+inone,3)
        else
            b_h=NaN
        end
        a_h=Xi(3-1+inone,2)
        c_h=-Xi(4-1+inone,3); 
        p_h=[a_h,b_h,c_h]  
    elseif isequal(eqname,'chenueta')
        ss_h=Xi(3-1+inone,1); 
        bb_h=-Xi(4-1+inone,3); 
        rr_h=Xi(3-1+inone,2);
        p_h=[ss_h,bb_h,rr_h]   
    elseif isequal(eqname,'lorenzchenueta') 
        %solution from Maple
        a_h = -(Xi(2-1+inone,1)*rho-sigma*Xi(3-1+inone,2)+sigma*Xi(2-1+inone,2)...
            +Xi(2-1+inone,1))/(Xi(2-1+inone,1)+Xi(3-1+inone,2)-Xi(2-1+inone,2)...
            +rho+sigma+1);
        b_h = (beta*Xi(2-1+inone,1)+beta*Xi(3-1+inone,2)-Xi(4-1+inone,3)*rho...
            -Xi(4-1+inone,3)*sigma-beta*Xi(2-1+inone,2)-Xi(4-1+inone,3))/(Xi(2-1+inone,1)+Xi(3-1+inone,2)-Xi(2-1+inone,2)+rho+sigma+1); 
        c_h = -(-Xi(3-1+inone,2)*rho-sigma*Xi(3-1+inone,2)+Xi(2-1+inone,1)...
            -Xi(2-1+inone,2))/(Xi(2-1+inone,1)+Xi(3-1+inone,2)-Xi(2-1+inone,2)+rho+sigma+1);
        phi_h = -(Xi(2-1+inone,1)+Xi(3-1+inone,2)-Xi(2-1+inone,2))/(rho+sigma+1);
        p_h=[sigma,beta,rho,a_h,b_h,c_h,phi_h]   
        %rather clunky, insufficient Xi parameters to solve explicitly
        %better to incorporate whole matrix Xi
    elseif isequal(eqname,'Vance') 
        alpha(1, 1) = -Xi(4 + inone, 1); 
        alpha(1, 2) = -Xi(5 + inone, 1); 
        alpha(1, 3) = - Xi(6 + inone, 1); 
        r(1) = Xi(1 + inone, 1); 
        alpha(2, 1) = -Xi(5 + inone, 2); 
        alpha(2, 2) = -Xi(7 + inone, 2); 
        alpha(2, 3) = -Xi(8 + inone, 2); 
        r(2) = Xi(2 + inone, 2); 
        alpha(3, 1) = -Xi(6 + inone, 3); 
        alpha(3, 2) = -Xi(8 + inone, 3); 
        alpha(3, 3) = -Xi(9 + inone, 3); 
        r(3) = Xi(3 + inone, 3);
        p_h=[r(1), r(2), r(3), alpha(1,1), alpha(1,2), alpha(1,3), ...
            alpha(2,1), alpha(2,2), alpha(2,3), ...
            alpha(3,1), alpha(3,2), alpha(3,3)];
    elseif isequal(eqname,'shilnikov')
        lambda = -Xi(2+inone, 2); 
        alpha = Xi(3+inone, 3); 
        B = - Xi(10+inone, 2); 
        p_h=[alpha, lambda, B]
    else
        'checktype sp not implemented for this system' 
    end

    if (mods==1)
        clear t_h x_h dx_h;
        %odeopt = odeset('RelTol',1e-8,'AbsTol',1e-10,'Events', @events_diverge); %Niven options
        %
        %choice of pde solver and solution of pde
        [t_h,x_h]=ode45(eqsys,tspan,x0,odeopt,p_h);       %nonstiff, 1 step, 1st try
        %[t,x]=ode23(eqsys,tspan,x0,odeopt,p);       %nonstiff, 1 step, mild stiffness
        %[t,x]=ode113(eqsys,tspan,x0,odeopt,p);       %nonstiff, multistep
        %[t,x]=ode15s(eqsys,tspan,x0,odeopt,p);       %stiff, multistep
        %[t,x]=ode23s(eqsys,tspan,x0,odeopt,p);       %stiff, 1 step
        %[t,x]=ode23t(eqsys,tspan,x0,odeopt,p);       %mod stiff, trap rule
        %[t,x]=ode23tb(eqsys,tspan,x0,odeopt,p);       %mod stiff, RK rule
      
        %protects against breakage - use length of x_h instead of x
        m_h = length(x_h); 
        % compute Derivative
        if isequal(eqname,'lorenz')
            dx_h = lorenzmat(t_h,x_h,p_h);  %Abel speedup
        elseif isequal(eqname,'rossler')
            dx_h = rosslermat(t_h,x_h,p_h);  
        elseif isequal(eqname,'chenueta')
            dx_h = chenuetamat(t_h,x_h,p_h);  
        elseif isequal(eqname,'Vance')
            dx_h = Vancemat(t_h,x_h,p_h);  
        else
           for i=1:m_h   
                dx_h(i,:) = eqsys(t_h(i),x_h(i,:),p_h);  %Niven call
           end
        end
        xcut = x(1:m_h,:);
        dxcut = dx(1:m_h,:);
    end;

elseif isequal(checktype,'mat')
    %compute directly from chosen alphabet and inferred parameters
    %%using generic fn call to polynomial function, using entire Xi matrix
    % already sped up!
    clear t_h x_h dx_h;
    
    %odeopt = odeset('RelTol',1e-15,'AbsTol',1e-5*ones(1,n)); 
    [t_h,x_h]=ode45(@genericpoly,tspan,x0,odeopt,Xi,n,inone,polyorder,...
        laurent,absshift,sineorder,intime,sintime);       
    %nonstiff, 1 step, 1st try
    %follows revised genericpoly commands, Nov 2021
    %[t_h,x_h]=ode23(@genericpoly,tspan,x0,odeopt,Xi,n,polyorder,inone,...
        %sineorder,laurent,intime,sintime);       %nonstiff, 1 step, mild stiffness
    %[t,x]=ode113(eqsys,tspan,x0,odeopt,p);       %nonstiff, multistep
    %[t,x]=ode15s(eqsys,tspan,x0,odeopt,p);       %stiff, multistep
    %[t,x]=ode23s(eqsys,tspan,x0,odeopt,p);       %stiff, 1 step
    %[t,x]=ode23t(eqsys,tspan,x0,odeopt,p);       %mod stiff, trap rule
    %[t,x]=ode23tb(eqsys,tspan,x0,odeopt,p);       %mod stiff, RK rule
        
    %protects against breakage - use length of x_h instead of x
    m_h = length(x_h); 
    for i=1:m_h
        dx_h(i,:) = genericpoly(t_h(i),x_h(i,:),Xi,n,inone,polyorder,...
        laurent,absshift,sineorder,intime,sintime);
    end
    xcut = x(1:m_h,:);
    dxcut = dx(1:m_h,:);

end

% variable plots, any n  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% labels from setlabels.m
% hardwired x axis to all time domain
% NOT hardwired y axis limits, third plot

% These are plots of recalc x and dx, to compare with original x and dx

if altplot == 1
    for v=1:o
        xmax=max(max(altscale.*x, [],'all'), max(altscale.*x_h, [],'all'));
        xmin=min(min(altscale.*x, [],'all'), min(altscale.*x_h, [],'all'));
    end  
    %fixed y scale for x, x_h plots
    figure 
    subplot(3,1,1),plot(t,altscale.*x),ylabel('Variable','FontSize',18),...
        legend(altxlist,'Interpreter','latex','FontSize',14),...
        title('Noisy data: coordinates','FontSize',16),xticklabels([])
        ax = gca; ax.XAxis.FontSize = 16; ax.YAxis.FontSize = 16;
        ax.XLim = [0 tfinal];  ax.YLim =([xmin,xmax]); 
    subplot(3,1,2),plot(t_h,altscale.*x_h),...
        ylabel('Variable','FontSize',18),...
        legend(altxestlist,'Interpreter','latex','FontSize',14),...
        title('Est. data: coordinates','FontSize',16),xticklabels([])
        ax = gca; ax.XAxis.FontSize = 16; ax.YAxis.FontSize = 16;
        ax.XLim = [0 tfinal]; ax.YLim =([xmin,xmax]);
    subplot(3,1,3),plot(t_h,altscale.*(xcut-x_h)),...
        xlabel('t','FontSize',18),ylabel('Variable','FontSize',18),...
        legend(altdiffxestlist,'Interpreter','latex','FontSize',14),...
        title('Est. data: coordinate error','FontSize',16)
        ax = gca; ax.XAxis.FontSize = 16; ax.YAxis.FontSize = 16; 
        ax.XLim = [0 tfinal]; 
        %ax.YLim = [-1000 1000]; 
    fignum=fignum+1; 
    saveas(gcf,fullfile(pwd, dirname,strcat('Fig',int2str(fignum))),'epsc')
else
    for v=1:o
        xmax=max(max(x,[],'all'), max(x_h,[],'all'));
        xmin=min(min(x,[],'all'), min(x_h,[],'all'));
    end
    %fixed y scale for x, x_h plots
    figure 
    subplot(3,1,1),plot(t,x),ylabel('Variable','FontSize',18),...
        legend(xlist,'Interpreter','latex','FontSize',14),...
        title('Noisy data: coordinates','FontSize',16),xticklabels([])
        ax = gca; ax.XAxis.FontSize = 16; ax.YAxis.FontSize = 16;
        ax.XLim = [0 tfinal]; ax.YLim =([xmin,xmax]);
    subplot(3,1,2),plot(t_h,x_h),ylabel('Variable','FontSize',18),...
        legend(xestlist,'Interpreter','latex','FontSize',14),...
        title('Est. data: coordinates','FontSize',16),xticklabels([])
        ax = gca; ax.XAxis.FontSize = 16; ax.YAxis.FontSize = 16;
        ax.XLim = [0 tfinal]; ax.YLim =([xmin,xmax]);
    subplot(3,1,3),plot(t_h,xcut-x_h),xlabel('t','FontSize',18),ylabel('Variable','FontSize',18),...
        legend(diffxestlist,'Interpreter','latex','FontSize',14),...
        title('Est. data: coordinate error','FontSize',16)
        ax = gca; ax.XAxis.FontSize = 16; ax.YAxis.FontSize = 16; 
        ax.XLim = [0 tfinal]; 
        %ax.YLim = [-1000 1000]; 
    fignum=fignum+1; 
    saveas(gcf,fullfile(pwd, dirname,strcat('Fig',int2str(fignum))),'epsc')
end


%pause

%variable phase plots, n = 2 or 3 (or more)
if n == 2
    figure 
    subplot(1,2,1),plot(x(:,1),x(:,2)),...
        xlabel(xlist(1),'Interpreter','latex','FontSize',18),...
        ylabel(xlist(2),'Interpreter','latex','FontSize',18),...
        title('Noisy data: coord. phase space','FontSize',16), view([-140 20]), ...
        box on, grid on
        ax = gca; ax.XAxis.FontSize = 16; ax.YAxis.FontSize = 16; ax.ZAxis.FontSize = 16;
    subplot(1,2,2),plot(x_h(:,1),x_h(:,2)),...
        xlabel(xestlist(1),'Interpreter','latex','FontSize',18),...
        ylabel(xestlist(2),'Interpreter','latex','FontSize',18),...
        title('Est. data: coord. phase space','FontSize',16), view([-140 20]), ...
        box on, grid on
        ax = gca; ax.XAxis.FontSize = 16; ax.YAxis.FontSize = 16; ax.ZAxis.FontSize = 16;
    fignum=fignum+1; 
    set(gcf, 'renderer', 'painters');
    saveas(gcf,fullfile(pwd, dirname,strcat('Fig',int2str(fignum))),'epsc')
elseif n >= 3
    %3D phase plots
    %for >3D, cycles through 3D plots in order of parameters
    for u=1:n-2   
        figure 
        subplot(1,2,1),plot3(x(:,u),x(:,u+1),x(:,u+2)),...
            xlabel(xlist(u),'Interpreter','latex','FontSize',18),...
            ylabel(xlist(u+1),'Interpreter','latex','FontSize',18),...
            zlabel(xlist(u+2),'Interpreter','latex','FontSize',18),...
            title('Noisy data: coord. phase space','FontSize',16), view([-140 20]), ...
            box on, grid on
            ax = gca; ax.XAxis.FontSize = 16; ax.YAxis.FontSize = 16; ax.ZAxis.FontSize = 16;
        subplot(1,2,2),plot3(x_h(:,u),x_h(:,u+1),x_h(:,u+2)),...
            xlabel(xestlist(u),'Interpreter','latex','FontSize',18),...
            ylabel(xestlist(u+1),'Interpreter','latex','FontSize',18),...
            zlabel(xestlist(u+2),'Interpreter','latex','FontSize',18),...
            title('Est. data: coord. phase space','FontSize',16), view([-140 20]), ...
            box on, grid on
            ax = gca; ax.XAxis.FontSize = 16; ax.YAxis.FontSize = 16; ax.ZAxis.FontSize = 16;
        fignum=fignum+1; 
        set(gcf, 'renderer', 'painters');
        saveas(gcf,fullfile(pwd, dirname,strcat('Fig',int2str(fignum))),'epsc')
        %
        %extra plot at normal scale for figs; can delete
        figure
            plot3(x_h(:,u),x_h(:,u+1),x_h(:,u+2)),...
            xlabel(xestlist(u),'Interpreter','latex','FontSize',18),...
            ylabel(xestlist(u+1),'Interpreter','latex','FontSize',18),...
            zlabel(xestlist(u+2),'Interpreter','latex','FontSize',18),...
            title('Est. data: coord. phase space','FontSize',16), view([-140 20]), ...
            box on, grid on
            ax = gca; ax.XAxis.FontSize = 16; ax.YAxis.FontSize = 16; ax.ZAxis.FontSize = 16;
        fignum=fignum+1; 
        set(gcf, 'renderer', 'painters');
        saveas(gcf,fullfile(pwd, dirname,strcat('Fig',int2str(fignum))),'epsc')
    end
end

% deriv plots, any o  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% labels from setlabels.m
% hardwired x axis to all time domain
% NOT hardwired y axis limits, third plot
for v=1:o
    dxmax=max(max(dx,[],'all'), max(dx_h,[],'all'));
    dxmin=min(min(dx,[],'all'), min(dx_h,[],'all'));
end
%fixed y scale for dx, dx_h plots
figure 
subplot(3,1,1),plot(t,dx),ylabel('Derivative','FontSize',18),...
    legend(dxlist,'Interpreter','latex','FontSize',14),...
    title('Noisy data: derivatives','FontSize',16),xticklabels([])
    ax = gca; ax.XAxis.FontSize = 16; ax.YAxis.FontSize = 16; 
    ax.XLim = [0 tfinal]; ax.YLim =([dxmin,dxmax]);
subplot(3,1,2),plot(t_h,dx_h),ylabel('Derivative','FontSize',18),...
    legend(dxestlist,'Interpreter','latex','FontSize',14),...
    title('Est. data: derivatives','FontSize',16),xticklabels([])
    ax = gca; ax.XAxis.FontSize = 16; ax.YAxis.FontSize = 16; 
    ax.XLim = [0 tfinal]; ax.YLim =([dxmin,dxmax]);
subplot(3,1,3),plot(t_h,dxcut-dx_h),xlabel('t','FontSize',18),ylabel('Deriv. error','FontSize',18),...
    legend(diffdxlist,'Interpreter','latex','FontSize',14),...
    title('Est. data: derivative error','FontSize',16)
    ax = gca; ax.XAxis.FontSize = 16; ax.YAxis.FontSize = 16; 
    ax.XLim = [0 tfinal]; 
fignum=fignum+1; 
saveas(gcf,fullfile(pwd, dirname,strcat('Fig',int2str(fignum))),'epsc')

%deriv phase plots, o = 2 or 3 (or more)
if o == 2
    figure 
    subplot(1,2,1),plot(dx(:,1),dx(:,2)),...
        xlabel(dxlist(1),'Interpreter','latex','FontSize',18),...
        ylabel(dxlist(2),'Interpreter','latex','FontSize',18),...
        title('Noisy data: deriv. phase space','FontSize',16), view([-140 20]), ...
        box on, grid on
        ax = gca; ax.XAxis.FontSize = 16; ax.YAxis.FontSize = 16; ax.ZAxis.FontSize = 16;
    subplot(1,2,2),plot(dx_h(:,1),dx_h(:,2)),...
        xlabel(dxestlist(1),'Interpreter','latex','FontSize',18),...
        ylabel(dxestlist(2),'Interpreter','latex','FontSize',18),...
        title('Est. data: deriv. phase space','FontSize',16), view([-140 20]), ...
        box on, grid on
        ax = gca; ax.XAxis.FontSize = 16; ax.YAxis.FontSize = 16; ax.ZAxis.FontSize = 16;
    fignum=fignum+1; 
    set(gcf, 'renderer', 'painters');
    saveas(gcf,fullfile(pwd, dirname,strcat('Fig',int2str(fignum))),'epsc') 
elseif o >= 3
    %3D phase plots
    %for >3D, cycles through 3D plots in order of parameters
    for u=1:o-2   
        figure 
        subplot(1,2,1),plot3(dx(:,u),dx(:,u+1),dx(:,u+2)),...
            xlabel(dxlist(u),'Interpreter','latex','FontSize',18),...
            ylabel(dxlist(u+1),'Interpreter','latex','FontSize',18),...
            zlabel(dxlist(u+2),'Interpreter','latex','FontSize',18),...
            title('Noisy data: deriv. phase space','FontSize',16), view([-140 20]), ...
            box on, grid on
            ax = gca; ax.XAxis.FontSize = 16; ax.YAxis.FontSize = 16; ax.ZAxis.FontSize = 16;
        subplot(1,2,2),plot3(dx_h(:,u),dx_h(:,u+1),dx_h(:,u+2)),...
            xlabel(dxestlist(u),'Interpreter','latex','FontSize',18),...
            ylabel(dxestlist(u+1),'Interpreter','latex','FontSize',18),...
            zlabel(dxestlist(u+2),'Interpreter','latex','FontSize',18),...
            title('Est. data: deriv. phase space','FontSize',16), view([-140 20]), ...
            box on, grid on
            ax = gca; ax.XAxis.FontSize = 16; ax.YAxis.FontSize = 16; ax.ZAxis.FontSize = 16;
        fignum=fignum+1; 
        set(gcf, 'renderer', 'painters');
        saveas(gcf,fullfile(pwd, dirname,strcat('Fig',int2str(fignum))),'epsc')
    end
end




