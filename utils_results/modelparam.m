% Copyright 2018-2024, All Rights Reserved
% Code by Robert Niven
% Looks up equation system and specifies parameters
% Need to go through and connect params with Xi

%Xi indices below include a row representing 1s
%if inone=0, this will be stripped out by the formula

    sbr=[10 8/3 28]; %params sigma, beta, rho of std Lorenz; 
    %also in shilnikov model
    sigma = sbr(1);  
    beta = sbr(2);  
    rho = sbr(3);
    %

    altplot = 0; %to override for 4D systems
    %altplot = 1; %4D systems

%note tspan now in main code
%label specs are in setlabels.m
   
if n==2 
    if isequal(eqname,'LotkaVolterra')
        x0=[10 10];   %ICs
        %tspan=[0,500];  %time axis
        %p=[0.1 0.1 100 150 1 -2];  %params - one dies 
        %p=[0.5 0.1 100 500 -3 -0.3];  %params - stable
        %p=[0.5 0.1 100 2000 -3 -0.3];  %params - stable
        %p=[5 0.1 200 500 0.1 0.2];  %params - stable with switchover
        %p=[5 0.1 2000 500 2 -0.3];  %params - stable with switchover
        %p=[0.5 0.1 2000 500 1 0.4];  %params - rise then fall
        p=[0.2 0.1 2000 2000 2 -2];  %params - one rises then dies
        %varlab=['x, y, xdot, ydot'];
        %xlab2={'x'; 'y'}; 
        %xlab2deriv={'xdot'; 'ydot'};
        %xlab3=xlab2;
        %fig2on=1;
        %xscale = [1; 1; 1];
    elseif isequal(eqname,'predatorprey')
        x0=[1 10];   %ICs
        %tspan=[0,50];  %time axis
        %p=[0.5 1000 2 1];  %params - non cyclic
        %p=[0 1000 1 1];  %params - perfect cyclic phase diag
        %p=[0.1 100 10 1];  %params - exploding phase diag
        %p=[-0.1 1000 1 1];  %params - converging phase diag
        p=[0.1 1000 1.1 0.1];  %params - converging phase diag
        %varlab=['x, y, xdot, ydot'];
        %xlab2={'x'; 'y'}; 
        %xlab2deriv={'xdot'; 'ydot'};
        %xlab3=xlab2;
        %fig2on=1;
        %xscale = [1; 1; 1];
    elseif isequal(eqname,'Brusselator')
        x0=[1 1];   %ICs
        %tspan=[0,100];  %time axis
        %p=[1.0 3.0];  %params 
        p=[1.0 2.0];  %alternative params 
        eps=eps/100; %less error added to this model (small ranges)
        %varlab=['x, y, xdot, ydot'];
        %xlab2={'x'; 'y'}; 
        %xlab2deriv={'xdot'; 'ydot'};
        %xlab3=%xlab2;
        %xscale = [1; 1; 1];
    end 
elseif n==3
    if isequal(eqname,'lorenz') | isequal(eqname,'lorenzmat')
        x0=[-2 0 27];   %ICs
        %x0=[-8; 8; 27];  % Initial condition
        %tspan=[0,100];  
        %tspan=[t0:tstep:tfinal]; %time axis
        %tspan now in main code
        p=sbr;  %params sigma, beta, rho of std Lorenz
        %%varlab=[' x, y, z, xdot, ydot, zdot'];
        %%xlab2={'x'; 'y'; 'z'};
        %xnlab2={'noisyx'; 'noisyy'; 'noisyz'};
        %%xlab2deriv={'xdot'; 'ydot'; 'zdot'};
        %%xscale = [200; 1; 1];
        %    
        true_Xi(2-1+inone,1)=-sigma; 
        true_Xi(2-1+inone,2)=rho; 
        true_Xi(3-1+inone,1)=sigma; 
        true_Xi(3-1+inone,2)=-1; 
        true_Xi(4-1+inone,3)=-beta;
        true_Xi(6-1+inone,3)=1; 
        true_Xi(7-1+inone,2)=-1;
    elseif isequal(eqname,'rossler') | isequal(eqname,'rosslermat')
        x0=[0 0 0];   %ICs
        %tspan=[0,1000];  %time axis
        p=[0.2 0.2 6.3];  %params
        a=p(1); b=p(2); c=p(3); 
        %varlab=[' x, y, z, xdot, ydot, zdot'];
        %xlab2={'x'; 'y'; 'z'};
        %xlab2deriv={'xdot'; 'ydot'; 'zdot'};
        %    
        if (inone==1)
            true_Xi(1-1+inone,3)=b; 
        end
        true_Xi(3-1+inone,2)=a; 
        true_Xi(4-1+inone,3)=-c; 
        true_Xi(2-1+inone,2)=1; 
        true_Xi(3-1+inone,1)=-1; 
        true_Xi(4-1+inone,1)=-1; 
        true_Xi(7-1+inone,3)=1;
        %    
    elseif isequal(eqname,'chenueta') | isequal(eqname,'chenuetamat')
        %x0=[-2 0 5];   %ICs
        x0=[20 20 30];   %ICs, Chenueta B
        %tspan=[0,100];  %time axis
        %p=sbr;  %params same as Lorenz; does not work?
        %p=[35 5 22.05];  %params in sbr order, chen & ueta A, 1st attractor
        p=[35 3 28];  %Chenueta B, params in sbr order, chen & ueta, Lynch p285
        %p=[35 1 25.264];  %params in sbr order, chen & ueta C, 2nd attractor
        ss=p(1); bb=p(2); rr=p(3);
        %varlab=[' x, y, z, xdot, ydot, zdot'];
        %xlab2={'x'; 'y'; 'z'};
        %xlab2deriv={'xdot'; 'ydot'; 'zdot'};
        %
        true_Xi(2-1+inone,1)=-ss; %y 
        true_Xi(2-1+inone,2)=rr-ss; %y
        true_Xi(3-1+inone,1)=ss;  %y
        true_Xi(3-1+inone,2)=rr;  %y
        true_Xi(4-1+inone,3)=-bb; %y
        true_Xi(6-1+inone,3)=+1; %y
        true_Xi(7-1+inone,2)=-1;  %y    
    elseif isequal(eqname,'Vance' ) | isequal(eqname,'Vancemat')
        x0=[100 100 100];   %educ guess
        %tspan=[0,100];  
        %tspan=[t0:tstep:tfinal]; %time axis
        %tspan now in main code
        p=[1 1 -1 0.001 0.001 0.01 0.0015 0.001 0.001 -0.005 -0.0005 0];  
        %params r(1)=p(1); r(2)=p(2); r(3)=p(3); 
        %alpha(1,1)=p(4); alpha(1,2)=p(5); alpha(1,3)=p(6); 
        %alpha(2,1)=p(7); alpha(2,2)=p(8); alpha(2,3)=p(9); 
        %alpha(3,1)=p(10); alpha(3,2)=p(11); alpha(3,3)=p(12);
        %eps=eps*10; %more error added to this model 
        %%varlab=[' x, y, z, xdot, ydot, zdot'];
        %%xlab2={'x'; 'y'; 'z'};
        %xnlab2={'noisyx'; 'noisyy'; 'noisyz'};
        %%xlab2deriv={'xdot'; 'ydot'; 'zdot'};
        %viewp=view([-140 20]);
        %%xscale = [200; 1; 1];
        %   from Maple
        true_Xi(4 + inone, 1) = -p(4);
        true_Xi(5 + inone, 1) = -p(5);
        true_Xi(6 + inone, 1) = -p(6);
        true_Xi(1 + inone, 1) = p(1);
        true_Xi(5 + inone, 2) = -p(7);
        true_Xi(7 + inone, 2) = -p(8);
        true_Xi(8 + inone, 2) = -p(9);
        true_Xi(2 + inone, 2) = p(2);
        true_Xi(6 + inone, 3) = -p(10);
        true_Xi(8 + inone, 3) = -p(11);
        true_Xi(9 + inone, 3) = -p(12);
        true_Xi(3 + inone, 3) = p(3);
    elseif isequal(eqname,'chua') | isequal(eqname,'chuamat')
        x0=[-1.6 0 1.6];   %ICs
        %tspan=[0,100];  %time axis
        p=[15 25.58 -5/7 -8/7];  %params - biattractor case
        %eps=eps/20; %less error added to this model (small ranges)
        %varlab=[' x, y, z, xdot, ydot, zdot'];
        %xlab2={'x'; 'y'; 'z'};
        %xlab2deriv={'xdot'; 'ydot'; 'zdot'};
        aq=p(1); bq=p(2); cq=p(3); dq=p(4);
        true_Xi(1+ inone, 1) = -aq*(cq + 1); 
        true_Xi(1+ inone, 2) = 1; 
        true_Xi(2+ inone, 1) = aq; 
        true_Xi(2+ inone, 2) = -1; 
        true_Xi(2+ inone, 3) = -bq; 
        true_Xi(3+ inone, 2) = 1; 
        if and(absshift==1, polyorder == 1)
            true_Xi(4+ inone, 1) = -aq*(cq - dq)/2; 
            true_Xi(5+ inone, 1) = aq*(cq - dq)/2;
        end
        %
    elseif isequal(eqname,'shilnikov') | isequal(eqname,'shilnikovmat')
        x0=[1 1 1];   %ICs
        %tspan=[0,1000];  %time axis
        sig=sbr(1); beta=sbr(2); rho=sbr(3);
        alp= beta/sqrt(sig*(rho-1)); lam=(1+sig)/sqrt(sig*(rho-1)); B=beta/(2*sig-beta); 
        p=[alp, lam, B]
        %params of shilnikov, passed to function
        %must avoid clash with lambda used in SINDy
        %eps=eps/10; %less error added to this model (small ranges)
        %p=[0.1623    0.6694    0.1538];
        %p=[0.16 0.66 0.15];
        %varlab=[' x, y, z, xdot, ydot, zdot'];
        %xlab2={'x'; 'y'; 'z'};
        %xnlab2={'noisyx'; 'noisyy'; 'noisyz'};
        %xlab2deriv={'xdot'; 'ydot'; 'zdot'};
        %%xscale = [200; 1; 1];  
        %
        true_Xi(1+inone, 2) = 1; 
        true_Xi(2+inone, 1) = 1; 
        true_Xi(2+inone, 2) = -lam; 
        true_Xi(3+inone, 3) = -alp; 
        true_Xi(4+inone, 3) = alp; 
        true_Xi(6+inone, 2) = -1; 
        if (polyorder >= 3)
            true_Xi(10+inone, 2) = -B; 
        end    
    elseif isequal(eqname,'galerkin3D')
        x0=[0.01 0.01 0.01];   %ICs
        %tspan=[0,100];  %time axis
        p=[1/10 1/10 1];  %params 
        eps=eps/1000; %less error added to this model (small ranges)
        %varlab=[' x, y, z, xdot, ydot, zdot'];
        %xlab2={'x'; 'y'; 'z'};
        %xlab2deriv={'xdot'; 'ydot'; 'zdot'};
    elseif or(isequal(eqname,'photosyn3spec'), isequal(eqname,'photosyn3speclin'))
        x0=[1 1 1];   %ICs
        %tspan=[0,0.1];  %time axis
        %p=[k(L) k(Lmin) k(d) k(dmin) k(plus1) k(min1) k(plus2) k(min2)];
        %p=[86.76 0.0016 1.18e8 0.00046 1.35E10 4.9E-10 238.3 4.74];  
        p=[100,0.01,1000,0.01,1000,0.00001,100,1];
        %varlab=['X, Y, Z, dotX, dotY, dotZ'];
        %xlab2={'X'; 'Y'; 'Z'}; 
        %xlab2deriv={'dotX'; 'dotY'; 'dotZ'}; 
        %%fig2on = 1;
        %%xscale = [1/200; 1; 1; 1/200];
        %%xlab3={'P/200'; 'Fe^{2+}'; 'Fe^{3+}'; 'C/200'};
    end
elseif n==4
    %generic labels (must have something)
    %xvarlist = {'w(t)', 'x(t)','y(t)','z(t)'};
    %dxvarlist = {'w^{\prime}(t)','x^{\prime}(t)','y^{\prime}(t)','z^{\prime}(t)'};
    %xvarlist={'w'; 'x'; 'y'; 'z'}; 
    %dxvarlist = {'w^{\prime}'; 'x^{\prime}'; 'y^{\prime}'; 'z^{\prime}'};
    %
    if isequal(eqname,'FeIIFeIII')
        x0=[500 0 0 0];   %ICs
        eps=eps/10; %less error added to this model (small ranges)
        %tspan=[0,1000];  %time axis
        %p=[0.01 0.01 0.0001 0.01];  %params
        p=[0.002 0.08 1 1];  %params - Hobbs settings
        %p=[0.002 0.08 0.4 1];  %e.g. params 
        %p=[0.002 0.08 0.5 1];  %e.g. params - extinguish osc over t=1:1000 for k(1)>0.02; k(4)<0.65; k(4)>3.5
        %varlab=['P, Fe^{2+}, Fe^{3+}, C, dotP, dotFe^{2+}, dotFe^{3+}, dotC'];
        %xvarlist={'P'; 'Fe^{2+}'; 'Fe^{3+}'; 'C'}; 
        %dxvarlist={'dotP'; 'dotFe^{2+}'; 'dotFe^{3+}'; 'dotC'};
        altplot = 1;
        %altscale = [1/200; 1; 1; 1/200];  %check altxlist, altxnoisylist in setlabels.m
        altscale = [1/200 1 1 1/200];  %check altxlist, altxnoisylist in setlabels.m
        %xlab3={'P/200'; 'Fe^{2+}'; 'Fe^{3+}'; 'C/200'};
        %truelab for x, x^2, x^3 terms
        true_Xi(1+inone, 1) = -p(1);
        true_Xi(1+inone, 2) = p(1);
        true_Xi(2+inone, 2) = -p(2);
        true_Xi(2+inone, 3) = p(2);
        true_Xi(3+inone, 3) = -p(4);
        true_Xi(3+inone, 4) = p(4);
        if polyorder == 3
          true_Xi(28+inone, 2) = -p(3);
          true_Xi(28+inone, 3) = p(3);
        end
    end
else
    p=[];
end
    

