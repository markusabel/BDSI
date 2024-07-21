%% Bayesian Dynamical System Identification (BDSI), version 9d
%% Extended version, with SINDy, JMAP, VBA, Lasso, ridge regression

% In this code we implement SINDy, two Bayesian algorithms (JMAP and VBA),
% and two traditional methods (LASSO and ridge regression),
% for the identification of a dynamical system from time series data.
% The Bayesian methods allow the quantification of model uncertainties.
%
% The program is structured to call a separate routine which specifies the
% dynamical system. A number of 2D, 3D, 4D dynamical systems are
% implemented.
%
% The Bayesian algorithms for a single set of hyperparameters are implemented
% for 2D, 3D, 4D systems; however the Bayesian iterations over hyperparameters
% are hardwired for 3D systems.
%
% This code was written by Robert Niven and Ali Mohammad-Djafari 2018-2024,
% with contributions from Laurent Cordier, Markus Abel and Markus Quade.
% It incorporates:
% - a previous ode solver written by Robert Niven
% - parts of the regularization code and SINDy routine of Brunton etal 2016 PNAS
% - the JMAP and VBA routines of Ali Mohammad-Djafari and Mircea Dumitru 2015

clear all; close all; clc;
addpath('utils')
changeplot

mods=1;  %1 for this code
addpath('dyn_sys_models')  %dyn system model library
addpath('noise_models')    %noise model library

%% Representations of data

%Dyn systems: contain [x1, x2,...,xn] and their derivs [x1dot, x2dot, ...xndot],
%all as functions of time t

%% Choice of system

dynsys=1;        %is a dynamical system model
RRmodel=0;       %not a rainfall-runoff model (NOT IMPLEMENTED HERE)
hydrolmodel=0;   %not a hydrological model (NOT IMPLEMENTED HERE)

knownXi='Y';     %calc data - known param values in Xi; can use checktype'sp' or 'mat'
%knownXi='N';     %calc data - Xi not known; omits or changes some plots

%%2D dyn systems  NOT FULLY IMPLEMENTED
%n = 2; %dim of system
%eqsys=@LotkaVolterra    %(competitive 2 herbivores) %not chaotic
%eqsys=@predatorprey
%eqsys=@Brusselator;

%%3D dyn systems
n = 3; %dim of system
eqsys=@lorenz
%eqsys=@chenueta             %sim to Lorenz; variants A, B, C; nice mix with lorenz
%eqsys=@chua                 %nice flat double attractor; contains shifted abs values
%eqsys=@rossler              %contains const term
%eqsys=@shilnikov            %contains cubic term, must setup Xi carefully
%eqsys=@Vance
%eqsys=@galerkin3D           %periodic
%eqsys=@photosyn3spec;

%%4D dyn systems   NOT FULLY IMPLEMENTED
%n = 4; %dim of system
%eqsys=@FeIIFeIII;

%% Function handling ****************************************
eqname = func2str(eqsys); %To construct character vector from function handle.
disp(['Equation name: ' eqname])% to display the eqname

%% Import Data
%not implemented here

%% Selection of options and parameter settings
%for all dynamical systems here
n;       %dimensions of input data (X)
o = n;   %dimensions of output data (here Xdot); %(if Y, must specify);

timeseries = 1;    %handle as a time series
%timeseries = 0;    %not a time series

posseries = 0;     %not a position series
%posseries = 1;      %handle as a 1D position series
%posseries = 2;      %handle as a 2D position series   %%NOT YET IMPLEMENTED
%posseries = 3;      %handle as a 3D position series   %%NOT YET IMPLEMENTED

%polyorder = 1;
polyorder = 2;
%polyorder = 3;    %IMPLEMENTED
%polyorder = 4;    %IMPLEMENTED

inone=0;        %do not have basis column of 1s
%inone=1;       %include basis column of 1s

sineorder = 0;    % no trig functions
%sineorder = 1;    % with trig functions, to sin(x), cos(x)
%sineorder = 2;    % with trig functions, only up to sin(2x), cos(2x)

laurent = 0;    % no reciprocal functions
%laurent = 1;    % with reciprocal functions
%laurent = 2;    % with 2nd order reciprocal functions

absshift=0;   %shifted absolute values (linear)
%absshift=1;   %shifted absolute values (linear), for Chua

intime=0;       %do not include time in alphabet
%intime=1;      %include time in alphabet

sintime=0;       %do not include sin/cos time in alphabet
%sintime=1;      %include sin/cos time in alphabet

%linear or log transformed variables
vartype='lin';
%vartype='log';    %NOT IMPLEMENTED

derivtype='call';
%here use a function call to calculate the derivatives precisely
%other choices not set up

%checktype='sp';  %uses original dyn eq with coeffs from Xi matrix, for check calcs (seems v slow)
checktype='mat'; %uses full Xi matrix for check calcs

special='nospec';          %no special condition

if timeseries == 1
    %tstep=0.05;   %
    tstep=0.01   %Lorenz, Chenueta B
    %tstep=0.02  %Vance, Shilnikov, others, 0.05 FeIIFeIII
    %tfinal=1.0  %Lorenz, special
    %tfinal=10  %Lorenz fast
    %tfinal=50  %Shilnikov fast
    tfinal=100  %Lorenz main, Chenueta,
    %tfinal=200   %Chua
    %tfinal=250   %Rossler
    %tfinal=500 %Vance, Shilnikov
    %tfinal=1000 % Shilnikov, FeIIFeIII
    t0=tstep;
    tspan=[t0:tstep:tfinal]; %time axis
    %tsamp=tstep;
    m=length(tspan); %length of time series
end

%*** noise model and parameters ***********************************************

%noisetype=0;   %no added noise
%noisetype=1;   %added Gaussian noise
noisetype=2;   %added Laplace noise
%noisetype=3;    %added Rayleigh noise
%noisetype=4;    %added Erlang noise

%noisy=0; %don't incorporate noisy x into derivative
noisy=1;  %incorporate noisy x into derivative

%for Gaussian, Laplace, Rayleigh, Erlang models
%eps is the multiplier of magnitude
%eps = 0 %no noise
%eps = 0.02; % noise, Shilnikov, Chua
%eps = 0.1; % noise,
eps = 0.2; % noise, Lorenz
%eps = 0.3; % noise, Chenueta B
%eps = 0.5; % noise,
%eps = 1.0; % noise, Lorenz special
%eps = 2.0; % noise, Lorenz special
%eps = 5.0; % noise, worse case

%for Gaussian, Laplace models: mean value
mu=0; %zero mean noise
%mu=1; %shifted mean noise
%special='shiftnoise1';      %shifted noise

%for Gaussian, Laplace, Rayleigh noise: std deviation or shift
sigg=1;

%for Erlang noise: parameters
ss=1;
lamb=1;

%*** SINDy parameters ***********************************************
%runSINDy=0;     %do not run SINDy
runSINDy=1;     %run SINDy

%same lambda0 across all columns
lambda0=10    %for iterations, used if iterating
%if not iterating:
%lambda0=0.01    %default
%lambda0=0.001   %Shilnikov
%lambda0=0.0001  %Vance
%subjective choice based on known values of coeffs

%iterateSINDy=0;      %no iteration
iterateSINDy=1;     %iterate to convergence in lambda criterion; IMPLEMENTED

if iterateSINDy==0
    lambdalim=NaN;    %dummy
elseif iterateSINDy==1
    lambdalim=1.0e-8;  %iterate to crit0 criterion, 1.0e-5 most; 1.0e-8 FeII
    lambdafact=(10)^(1/2);    %mult or division factor
end


%*** Bayesian routines and parameters ********************************
%JMAP, VBA functions
jmapfn=@JMAP6;        %regular algorithm, uses det tolerance
vbafn=@VBA6;          %regular algorithm, uses det tolerance
jmap=func2str(jmapfn);
vba=func2str(vbafn);

ntol=0.01;  %tolerance on inverse Delta calc, internal JMAP6, VBA6
%ntol=1.0e-6;  logntol=log10(ntol); %high precision

%iterateBayes=0;    %no iteration
iterateBayes=1;     %values 1=iterate over error (default)
                    %       2=iterate over model
                    %       3=iterate over both model and error (both a and b, equal values)
                    %       4=iterate over error (vary be, fixed ae)
                    %       5=iterate over error (vary ae, fixed be)

ae0=10; be0=0.1;          %large, allows many iterations;
                          %used also as high fixed values for iterateBayes=2

af0=10; bf0=0.1;          %broad f prior; use for iterateBayes=1 (error iteration)
                          %also use as start value for iterateBayes=2 (model iteration)
                          %can use for iterateBayes=3 (model and error iteration)

%crit0, only for iterateBayes
    %crit0_tol=1e-12;    %very good for run testing
    crit0_tol=1e-30;   %iterate to crit0 criterion, most dyn sys

%param_fact, only for IterateBayes
param_fact=10;    %mult or division factor

%*** LASSO parameters *****************

runLASSO=0;     %do not run LASSO
%runLASSO=1;     %run LASSO

kappa0=10    %for iterations
%if not iterating
%kappa0=0.01    %default, used
%subjective choice based on known values of coeffs

%iterateLASSO=0;      %no iteration
iterateLASSO=1;     %iterate to convergence in kappa criterion; IMPLEMENTED

if iterateLASSO==0
    kappalim=NaN;    %dummy
elseif iterateLASSO==1
    kappalim=1.0e-8;  %iterate to crit0 criterion, 1.0e-5 most; 1.0e-8 FeII
    kappafact=(10)^(1/2);    %mult or division factor
end


%*** Ridge regression parameters *****************

runRidge=0;     %do not run Ridge
%runRidge=1;     %run Ridge

theta0=10    %for iterations
%if not iterating
%theta0=0.01    %default
%subjective choice based on known values of coeffs

scaleRidge=0;  %normalises data and rescales; returns constant term in output
%scaleRidge=1;  %DEFAULT, normalises data, does not rescale

%iterateRidge=0;      %no iteration
iterateRidge=1;     %iterate to convergence in theta criterion; IMPLEMENTED

if iterateRidge==0
    thetalim=NaN;    %dummy
elseif iterateRidge==1
    thetalim=1.0e-17;  %iterate to crit0 criterion, need v low for ridge
    thetafact=(10);    %mult or division factor
end

%%%%%%%%%%%%%%

if runLASSO==1
    LAS='_LASSO'
else
    LAS=''
end

%options encoded as output directory *****************************
dirname = ['ex_', eqname, '_dim', num2str(n,1), ...
    '_ep', num2str(eps,3),'_mu', num2str(mu,3),...
    '_deriv', derivtype, '_inone', int2str(inone), ...
    '_poly', int2str(polyorder), ...
    '_absshift', int2str(absshift), ...
    '_noisetype', int2str(noisetype), ...
      '_', checktype, '_', special,  ...
    '_iterSINDy', num2str(iterateSINDy,1),...
    '_iterBayes', num2str(iterateBayes,2),...
    '_t', num2str(tfinal,4), '_tstep', num2str(tstep,4), LAS]

mkdir([dirname])      %directory for exported eps files


%% Output of figure headings and parameter values
% This outputs a special text figure with all options listed

fignum=1;

if strcmp(eqname,'chenueta')==1
    Eqname='Chen-Ueta';
else
    Eqname = eqname;
    Eqname(1)=upper(eqname(1));
end;

fighead=strcat(Eqname, ' system (', num2str(n,1), 'D)')
fighead2=[];

str0=[fighead, '; polynomial order ', int2str(polyorder)];
str0a=[fighead2];

if noisetype ==0
    str1=['No added noise'];
elseif noisetype ==1
    str1=['Added Gaussian noise: \mu=', num2str(mu,3),' \sigma=', ...
        num2str(sigg,3), ' (scaled by \epsilon=', num2str(eps,3), ')'];
elseif noisetype ==2
    str1=['Added Laplace noise: \mu=', num2str(mu,3), ' \sigma=', ...
        num2str(sigg,3), ' (scaled by \epsilon=', num2str(eps,3), ')'];
elseif noisetype ==3
    str1=['Added Rayleigh noise: \sigma=', num2str(sigg,3), ' (scaled by \epsilon=', num2str(eps,3), ')'];
elseif noisetype ==4
    str1=['Added Erlang noise: s=', num2str(ss,3),', \lambda=', ...
        num2str(lamb,3), ' (scaled by \epsilon=', num2str(eps,3), ')'];
end;

str2=['Derivatives calculated by ' derivtype];
if(inone==0)
    str3=['Not incl. 1s; '];
elseif(inone==1)
    str3=['Incl. 1s; '];
end
if(sineorder==0)
    str4=[str3, 'not incl. sin/cos functions'];
elseif(sineorder>=1)
    str4=[str3, 'incl. sin/cos to order ' int2str(sineorder)];
end
if(laurent==0)
    str5=[str4, '; not incl. reciprocal functions'];
elseif(laurent>=1)
    str5=[str4, '; incl. reciprocals to order ' int2str(laurent)];
end

if(absshift==0)
    str6a=['Not incl. shifted abs values; '];
elseif(absshift==1)
    str6a=['Incl. shifted abs values; '];
end

if(intime==0)
    str6b=['not incl. time; '];
elseif(intime==1)
    str6b=['incl. time; '];
end

if(sintime==0)
    str6=[str6a, str6b, 'not incl. sin/cos time; '];
elseif(sintime==1)
    str6=[str6a, str6b, 'incl. sin/cos time; '];
end
if isequal(checktype,'sp')
    str7=['Inferred time series eval. using sparse extracted coefficients'];
elseif isequal(checktype,'mat')
    str7=['Inferred time series eval. using full matrix of coefficients'];
end

str8=['time =', num2str(tfinal,5), ', tstep =', num2str(tstep,5)];

if runSINDy == 1
    str9=['SINDy: consts lambda0=', num2str(lambda0,6), ...
    ', lambdalim=' num2str(lambdalim,6), ...
       ', iterate=',num2str(iterateSINDy,1)];
else
    str9=['SINDy not run'];
end

str10=['Bayes: consts ae0=', num2str(ae0,4), ...
    ', be0=', num2str(be0,4), ...
    ', af0=', num2str(af0,4), ...
    ', bf0=', num2str(bf0,4), ...
    ', crit0tol=' num2str(crit0_tol,6), ...
    ', iterate=',num2str(iterateBayes,2)];

if runLASSO == 1
    str11=['LASSO: consts kappa0=', num2str(kappa0,6), ...
    ', kappalim=' num2str(kappalim,6), ...
       ', iterate=',num2str(iterateLASSO,1)];
else
    str11=['LASSO not run'];
end

if runRidge == 1
    str12=['Ridge: consts theta0=', num2str(theta0,6), ...
    ', thetalim=' num2str(thetalim,6), ...
       ', iterate=',num2str(iterateRidge,1),...
       ', scale param = ', num2str(scaleRidge,1)];
else
    str12=['Ridge regression not run'];
end

figure
    annotation('textbox',...
    [0.05 0.05 0.9 0.55],...
    'String',{str0, str0a, str1, str2, str5, str6, str7, str8, ...
        str9, str10, str11, str12},...
    'FontSize',14,...
    'FontName','Helvetica');
    saveas(gcf,fullfile(pwd, dirname,strcat('Fig',int2str(fignum))),'epsc')

pause


%% Setup of alphabet and output of X(t) and dotX(t) from dynamical system

%number of alphabet functions cc = Xi_rows
%uses figurate numbers; accounts for sin/cos, abshift and laurent
%alphabet based on n is correct
cc = nchoosek(polyorder+n, n) - 1 + inone + (2*n*sineorder) ...
    + nchoosek(laurent+n, n)-1 + 2*n*absshift + intime + 2*sintime

%number output format
%format shortEng;
format longEng;

true_Xi=zeros(cc,o);   %protects against missing Xi data; o allows for dimension of Y data
seqcols=[1:1:cc];
colslist=string(seqcols);
seqY=[1:1:o];

%set up model coefficients
modelparam; %looks up model, specifies parameters x0, p, true Xi variables;
%can override here if do correctly
x0
p

%pause

%calculation of data - integration of dyn sys ********************
if(mods==1)  %Niven library
    odeopt = odeset('RelTol',1e-12,'AbsTol',1e-5*ones(1,n));
    %
    %choice of pde solver and solution of pde
    [t,x]=ode45(eqsys,tspan,x0,odeopt,p);       %nonstiff, 1 step, 1st try
    %[t,x]=ode23(eqsys,tspan,x0,odeopt,p);       %nonstiff, 1 step, mild stiffness
    %[t,x]=ode113(eqsys,tspan,x0,odeopt,p);       %nonstiff, multistep
    %[t,x]=ode15s(eqsys,tspan,x0,odeopt,p);       %stiff, multistep
    %[t,x]=ode23s(eqsys,tspan,x0,odeopt,p);       %stiff, 1 step
    %[t,x]=ode23t(eqsys,tspan,x0,odeopt,p);       %mod stiff, trap rule
    %[t,x]=ode23tb(eqsys,tspan,x0,odeopt,p);       %mod stiff, RK rule
    %comment from Karsten - use stiff solver e.g. Rosenbrock, for 2 time scales
    %
    % compute Derivative
    %numx=length(x); %this is just m
    if isequal(eqname,'lorenz')
        dx = lorenzmat(t,x,p);  %Abel speedup
    elseif isequal(eqname,'rossler')
        dx = rosslermat(t,x,p);
    elseif isequal(eqname,'chenueta')
        dx = chenuetamat(t,x,p);
    elseif isequal(eqname,'Vance')
        dx = Vancemat(t,x,p);
    elseif isequal(eqname,'shilnikov')
        dx = shilnikovmat(t,x,p);
    elseif isequal(eqname,'chua')
        dx = chuamat(t,x,p);
    else
        for i=1:m
           dx(i,:) = eqsys(t(i),x(i,:),p);  %Niven call
        end
    end
end

%labels, colours, symbols of plots ********************
setlabels;

%plot generated data ********************
paramhead ='Raw data: coordinates';
paramsphead = 'Raw data: coordinate phase space';
fnhead = 'Raw data: derivatives';
fnsphead = 'Raw data: derivative phase space';
withnoise=0;
X = x; dX = dx;
Xlist = xlist; dXlist = dxlist;
if altplot == 1
    altXlist = altxlist;
end;
plotdynsys; %must come after setlabels

%pause

%% Addition of noise

% add noise to the state variables
% uses array based on size(x)
% now extended to include shift mu, and other noise functions

rng(10); %controls seed
u=rand(size(x));   %random number matrix, values in [0,1]
noisex = zeros(size(u));
yv = 0:0.01:1;

if noisetype ==0
    %nothing
    invcdf=zeros(size(yv));

elseif noisetype ==1     %added Gaussian noise
    invcdf=invGaussianmat(yv,mu,sigg);
    %noisex = randn(size(x));  Gaussian with 0 mean, 1 variance
    noisex = invGaussianmat(u,mu,sigg);

elseif noisetype ==2    %added Laplace noise
    invcdf=invLaplacemat(yv,mu,sigg);
    noisex = invLaplacemat(u,mu,sigg);

elseif noisetype ==3    %added Rayleigh noise
    invcdf=invRayleighmat(yv,sigg);
    noisex = invRayleighmat(u,sigg);

elseif noisetype ==4    %added Erlang noise
    invcdf=invErlangmat(yv,ss,lamb);
    noisex = invErlangmat(u,ss,lamb);
end;

figure
    plot(yv,invcdf);
    xlabel('Cumulative probability','FontSize',18),...
    ylabel('x','FontSize',18),...
    %xlim(crit0lim), ...
    title(str1,'FontSize',18)
    ax = gca; ax.XAxis.FontSize = 16; ax.YAxis.FontSize = 16;
    fignum=fignum+1;
    saveas(gcf,fullfile(pwd, dirname,strcat('Fig',int2str(fignum))),'epsc')

%pause;

xnoisy = x + eps*noisex;

figure
    histogram(noisex,50,'Normalization','probability');
    xlabel('$x_{noise}$','Interpreter','latex','FontSize',18),...
    ylabel('frequency','FontSize',18),...
    %xlim(crit0lim), ...
    title(str1,'FontSize',18)
    ax = gca; ax.XAxis.FontSize = 16; ax.YAxis.FontSize = 16;
    fignum=fignum+1;
    saveas(gcf,fullfile(pwd, dirname,strcat('Fig',int2str(fignum))),'epsc')

pause;

%incorporate noise into derivative by function call
if (noisy==1)
    %x = xnoisy;
    if(mods==1)  %Niven library
        if isequal(eqname,'lorenz')
            dxnoisy = lorenzmat(t,xnoisy,p);  %Abel speedup routines
        elseif isequal(eqname,'rossler')
            dxnoisy = rosslermat(t,xnoisy,p);
        elseif isequal(eqname,'chenueta')
            dxnoisy = chenuetamat(t,xnoisy,p);
        elseif isequal(eqname,'Vance')
            dxnoisy = Vancemat(t,xnoisy,p);
        elseif isequal(eqname,'shilnikov')
            dxnoisy = shilnikovmat(t,xnoisy,p);
        elseif isequal(eqname,'chua')
            dxnoisy = chuamat(t,xnoisy,p);
        else
            for i=1:m
                dxnoisy(i,:) = eqsys(t(i),xnoisy(i,:),p);  %Niven call
            end
       end
    end
end


%plot of generated data with added noise ********************
paramhead ='Noisy data: coordinates';
paramsphead = 'Noisy data: coordinate phase space';
fnhead = 'Noisy data: derivatives';
fnsphead = 'Noisy data: derivative phase space';
withnoise=1;
X = xnoisy; dX = dxnoisy; diffX=xnoisy-x; diffdX=dxnoisy-dx;
Xlist = xnoisylist; dXlist = dxnoisylist;
if altplot == 1
    altXlist = altxnoisylist;
    altdiffXlist = altdiffxlist;
end;
diffXlist = diffxlist; diffdXlist = diffdxlist;

plotdynsys;

%assign noisy data to data
x = xnoisy;
dx = dxnoisy;
clear xnoisy dxnoisy;


%pause

%% Build the Dictionary
% pool Data  (i.e., build library of nonlinear time series)
H = poolData3(t,x,n,inone,polyorder,laurent,absshift,sineorder,intime,sintime);

%use of n = no of input variables is correct
%cc = size(H,2)   #already defined

%pause

HtH=H'*H;

%set up column colouring for matrices
colM=arrayfun(@colorlog,H);

%rotation of labels
%keep as x description for consistency across subroutines
if cc >16
    xlabrot=90;
    xlabsize=10;
elseif cc >10
    xlabrot=45;
    xlabsize=16;
else
    xlabrot=0;
    xlabsize=16;
end;

%plots of H, HTH
figure,imagesc(colM),title({fighead; 'Matrix H'},'FontSize',18),...
    xlabel('Dictionary','FontSize',18),xlim([0.5 cc+0.5]),...
    ylabel('Time step','FontSize',18),...
    colorbar,...
    %ylabel(colorbar,'shifted log scale','FontSize',16)
    %caxis([-maxH maxH]);
    ax = gca; ax.XAxis.FontSize = xlabsize; ax.YAxis.FontSize = 16;
    ax.XTick = seqcols; ax.XTickLabel = dictlist;
    ax.XTickLabelRotation = xlabrot;
    c=colorbar; c.FontSize = 14;
    fignum=fignum+1;
    saveas(gcf,fullfile(pwd, dirname,strcat('Fig',int2str(fignum))),'epsc')


colM=arrayfun(@colorlog,HtH);

figure,imagesc(colM),title({fighead; 'Matrix H^tH'},'FontSize',18),...
    xlabel('Dictionary','FontSize',18),xlim([0.5 cc+0.5]),...
    ylabel('Dictionary','FontSize',18),ylim([0.5 cc+0.5]),...
    %xticklabels(dictlist),yticklabels(dictlist), colorbar,...
    %ylabel(colorbar,'shifted log scale','FontSize',16)
    %caxis([-maxHtH maxHtH]);
    ax = gca; ax.XAxis.FontSize = xlabsize; ax.YAxis.FontSize = xlabsize;
    ax.XTick = seqcols; ax.XTickLabel = dictlist;
    ax.XTickLabelRotation = xlabrot;
    ax.YTick = seqcols; ax.YTickLabel = dictlist;
    c=colorbar; c.FontSize = 14;
    fignum=fignum+1;
    saveas(gcf,fullfile(pwd, dirname,strcat('Fig',int2str(fignum))),'epsc')

%pause




%% SINDy Method: Brunton etal 2016 PNAS modified

if runSINDy==0
    %fignum=196;

elseif runSINDy==1

    bayesian=0; SINDy=1; LASSO=0; Ridge=0;

    %format shortEng;
    format longEng;

    if iterateSINDy==0
        numit_SINDy=1
        %same lambda0 across all columns
        lambda_SINDy=lambda0
        true_Xi
        Xi = sparsifyDynamics(H,dx,lambda0,o);   %allow for separate Y dimension
        Xi_SINDy=Xi

        if knownXi=='Y';     %calc data - known param values in Xi
            diffXi=true_Xi-Xi     %known function
        elseif knownXi=='N'
            diffXi=Xi-Xi          %to show errors; unknown function
        end;

        %2norm quantities
        for i=1:o   %loop over Y dimensions o
            Resid2_SINDy(i)= (norm(dx(:,i) - H*Xi(:,i)))^2;
            Reg2_SINDy(i)= (norm(Xi(:,i)))^2;
            Jtilde2_SINDy(i)=Resid2_SINDy(i) + lambda0*Reg2_SINDy(i);  %Tik obj function
        end
        Resid2_SINDy
        Reg2_SINDy
        Jtilde2_SINDy

        % Show and compare the coefficients and generated time series
        Title_text='SINDy';
        if dynsys==1
             Show_comparisons6a;
        elseif hydrolmodel==1 | RRmodel==1
             Show_comparisons8;
        end;

    elseif iterateSINDy==1

        Iterate_SINDy

    end
end

%fignum=196;

%%pause

%% JMAP and VBA methods: Mohammad-Djafari & Dumitru 2015 modified
bayesian=1; SINDy=0; LASSO=0;  Ridge=0;

for Bayesmeth=1:2   %1=JMAP, 2=VBA
%for Bayesmeth=1:1   %JMAP only
%for Bayesmeth=2:2   %VBA only

    if Bayesmeth==1
        Title_text='JMAP';
    elseif Bayesmeth==2
        Title_text='VBA';
    end

    %format shortEng;
    format longEng;

    Xi=zeros(cc,o);       %initialise; allow for separate Y dimension
    Vf=zeros(size(Xi));   %initialise
    fval=zeros(1,o);      %initialise to protect nonopt pathway
    numintit=zeros(1,o);  %to avoid some traps in loops
    numextit=0;           %to avoid some traps in loops

    if iterateBayes==0
        'Isolated Bayes'
        Single_Bayes
        if dynsys==1
             Show_comparisons6a;
        elseif hydrolmodel==1 | RRmodel==1
             Show_comparisons8;
        end;

    elseif iterateBayes>=1 & iterateBayes<10
        'Iterate Bayes'
        %need to pick out the actual solution that wish to use
        Iterate_Bayes    %includes Show_comparisons for each iteration

    end

end


fignum=578;

%% LASSO Method: from Ali code 2015, modified

if runLASSO==0
    %fignum=196;

elseif runLASSO==1

    bayesian=0; SINDy=0; LASSO=1;  Ridge=0;

    %format shortEng;
    format longEng;

    if iterateLASSO==0
        numit_LASSO=1
        %same kappa0 across all columns
        kappa_LASSO=kappa0
        true_Xi
        Xi=zeros(cc,n);
        for i=1:o
            g=dx(:,i);
            [B,Lasso_FitInfo]=lasso(H,g,'Lambda',kappa0);
            Xi(:,i)=B;
        end
        Lasso_FitInfo
        Xi_LASSO=Xi

        if knownXi=='Y';     %calc data - known param values in Xi
            diffXi=true_Xi-Xi     %known function
        elseif knownXi=='N'
            diffXi=Xi-Xi          %to show errors; unknown function
        end;

        %2norm and Frobenius norm quantities
        for i=1:o   %loop over Y dimensions o
            Resid2_LASSO(i)= (norm(dx(:,i) - H*Xi(:,i)))^2;
            Reg2_LASSO(i)= (norm(Xi(:,i)))^2;
            Jtilde2_LASSO(i)=Resid2_LASSO(i) + kappa0*Reg2_LASSO(i);  %Tik obj function
            %
            Reg1_LASSO(i)=norm(Xi(:,i),1);
            J_LASSO(i)=Resid2_LASSO(i) + 2*m*kappa0*Reg1_LASSO(i);
        end

        Resid2_LASSO
        Reg2_LASSO
        Jtilde2_LASSO

        Reg1_LASSO
        J_LASSO

        % Show and compare the coefficients and generated time series
        Title_text='LASSO';
        if dynsys==1
             Show_comparisons6a;
        elseif hydrolmodel==1 | RRmodel==1
             Show_comparisons8;
        end;

    elseif iterateLASSO==1

        Iterate_LASSO

    end
end

%%pause

%% Ridge Regression Method

if runRidge==0
    %fignum=196;

elseif runRidge==1

    bayesian=0; SINDy=0; LASSO=0; Ridge=1;

    %format shortEng;
    format longEng;

    if iterateRidge==0
        numit_Ridge=1
        %same theta0 across all columns
        theta_Ridge=theta0
        true_Xi
        Xi=zeros(cc,n);
        for i=1:o
            g=dx(:,i);
            Br=ridge(g,H,theta0,scaleRidge); %note is opposite notation to lasso
            if scaleRidge==0        %normalises data and rescales; returns constant term in output
                Xi(:,i)=Br(2:cc+1);
            elseif scaleRidge==1    %DEFAULT, normalises data, does not rescale
                Xi(:,i)=Br;
            end
        end
        Xi_Ridge=Xi

        if knownXi=='Y';     %calc data - known param values in Xi
            diffXi=true_Xi-Xi     %known function
        elseif knownXi=='N'
            diffXi=Xi-Xi          %to show errors; unknown function
        end;

        %2norm quantities
        for i=1:o   %loop over Y dimensions o
            Resid2_Ridge(i)= (norm(dx(:,i) - H*Xi(:,i)))^2;
            Reg2_Ridge(i)= (norm(Xi(:,i)))^2;
            Jtilde2_Ridge(i)=Resid2_Ridge(i) + theta0*Reg2_Ridge(i);  %Tik obj function
            %this really does match the problem!
        end

        Resid2_Ridge
        Reg2_Ridge
        Jtilde2_Ridge

        % Show and compare the coefficients and generated time series
        Title_text='Ridge';
        if dynsys==1
             Show_comparisons6a;
        elseif hydrolmodel==1 | RRmodel==1
             Show_comparisons8;
        end;

    elseif iterateRidge==1

        Iterate_Ridge

    end
end
