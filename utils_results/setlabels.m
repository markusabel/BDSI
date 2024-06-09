%% Set up labels + colours + symbols
% this version is hardwired

% Copyright 2019-2024, All Rights Reserved
% Code by Robert Niven

%if dynsys==1
    %xlab = 't (h)';         %NOT IMPL
    %if vartype == 'lin'
        %ylab1 = '$Q$'; 
        %ylablist1 = {'$Q$'}; 
        %ylab2 = '$\hat{Q}$';
        %ylab3 = '$Q-\hat{Q}$';
    %end;
%end

%colours, symbols for data plots; up to n=12 and o=12
plotcolord = [0 0 1; 1 0 0; 0 0.8 0; 0 0.6 1; 0 1 0; 1 0 1; ...
    0.5 0.5 0.5; 0 0 0; 0.4660 0.6740 0.1880; 0.5 0 0; 0 0 0.5; 1 1 0];
%plotcolord= [0, 0.8, 0; 0.5, 0, 0; 0, 0, 0.5; 0 0.7 1; 0.5, 0.5, 0.5; ...
%      0.4660 0.6740 0.1880];
set(0,'DefaultAxesColorOrder',plotcolord); 
%this overrides default for entire session
%see http://matlab.izmiran.ru/help/techdoc/creating_plots/axes_p19.html

%colours, symbols for xi_ij plots
darkGreen = [0, 0.8, 0];
darkRed = [0.5, 0, 0];
darkBlue = [0, 0, 0.5];
midBlue = [0 0.7 1];
grey = [0.5, 0.5, 0.5];
khaki = [0.4660 0.6740 0.1880];

plotsym = {'^','o','x','s','d','v','*','>','<','+','.','p'};
plotcol = {'blue', 'red', darkGreen, midBlue, ...
    'green', 'magenta', grey, 'black', ...
    khaki, darkRed, darkBlue, 'yellow'};

%labels for Xi plots
coefflab={'${\xi}_{\ell1}$','${\xi}_{\ell2}$','${\xi}_{\ell3}$','${\xi}_{\ell4}$',...
    '${\xi}_{\ell5}$','${\xi}_{\ell6}$','${\xi}_{\ell7}$','${\xi}_{\ell8}$',...
    '${\xi}_{\ell9}$','${\xi}_{\ell10}$','${\xi}_{\ell11}$','${\xi}_{\ell12}$'};

coeffhatlab={'$\hat{\xi}_{\ell1}$', '$\hat{\xi}_{\ell2}$', '$\hat{\xi}_{\ell3}$', ...
    '$\hat{\xi}_{\ell4}$','$\hat{\xi}_{\ell5}$','$\hat{\xi}_{\ell6}$',...
    '$\hat{\xi}_{\ell7}$','$\hat{\xi}_{\ell8}$','$\hat{\xi}_{\ell9}$',...
    '$\hat{\xi}_{\ell10}$','$\hat{\xi}_{\ell11}$','$\hat{\xi}_{\ell12}$'};

coeffdifflab={'${\xi}_{\ell1}-\hat{\xi}_{\ell1}$','${\xi}_{\ell2}-\hat{\xi}_{\ell2}$',...
    '${\xi}_{\ell3}-\hat{\xi}_{\ell3}$','${\xi}_{\ell4}-\hat{\xi}_{\ell4}$',...
    '${\xi}_{\ell5}-\hat{\xi}_{\ell5}$','${\xi}_{\ell6}-\hat{\xi}_{\ell6}$',...
    '${\xi}_{\ell7}-\hat{\xi}_{\ell7}$','${\xi}_{\ell8}-\hat{\xi}_{\ell8}$',...
    '${\xi}_{\ell9}-\hat{\xi}_{\ell9}$','${\xi}_{\ell10}-\hat{\xi}_{\ell10}$',...
    '${\xi}_{\ell11}-\hat{\xi}_{\ell11}$','${\xi}_{\ell12}-\hat{\xi}_{\ell12}$'};

vanishlab = {'$\hat{\xi}_{\ell1}-\hat{\xi}_{\ell1}$','$\hat{\xi}_{\ell2}-\hat{\xi}_{\ell2}$',...
    '$\hat{\xi}_{\ell3}-\hat{\xi}_{\ell3}$','$\hat{\xi}_{\ell4}-\hat{\xi}_{\ell4}$',...
    '$\hat{\xi}_{\ell5}-\hat{\xi}_{\ell5}$','$\hat{\xi}_{\ell6}-\hat{\xi}_{\ell6}$',...
    '$\hat{\xi}_{\ell7}-\hat{\xi}_{\ell7}$','$\hat{\xi}_{\ell8}-\hat{\xi}_{\ell8}$',...
    '$\hat{\xi}_{\ell9}-\hat{\xi}_{\ell9}$','$\hat{\xi}_{\ell10}-\hat{\xi}_{\ell10}$',...
    '$\hat{\xi}_{\ell11}-\hat{\xi}_{\ell11}$','$\hat{\xi}_{\ell12}-\hat{\xi}_{\ell12}$'};


%for variable plots, deriv plots and alphabet construction
if n==2
    %varlist={'x''','y'''};
    varlist = {'x^{\prime}', 'y^{\prime}'};
    xlist = {'$x(t)$','$y(t)$'};
    dxlist = {'$x^{\prime}(t)$','$y^{\prime}(t)$'};
    xnoisylist = {'$x_{noisy}(t)$','$y_{noisy}(t)$'};
    dxnoisylist = {'$x_{noisy}^{\prime}(t)$','$y_{noisy}^{\prime}(t)$'};
    diffxlist = {'$x_{noisy}(t)-x(t)$','$y_{noisy}(t)-y(t)$'};
    diffdxlist = {'$x_{noisy}^{\prime}(t)-x^{\prime}(t)$', ...
        '$y_{noisy}^{\prime}(t)-y^{\prime}(t)$'};
    xestlist = {'$\hat{x}(t)$','$\hat{y}(t)$'};
    dxestlist = {'$\hat{x}^{\prime}(t)$','$\hat{y}^{\prime}(t)$'}
    diffxestlist = {'$x(t)-\hat{x}(t)$','$y(t)-\hat{y}(t)$'};
    diffdxestlist = {'$x^{\prime}(t)-\hat{x}^{\prime}(t)$',...
       '$y^{\prime}(t)-\hat{y}^{\prime}(t)$'}
    %
    lablist1={'x','y'};
    lablist2={'x^2','xy', 'y^2'};
    lablist3={'x^3','x^2y','xy^2','y^3'};
    lablist4={'x^4'}; 
    lablist5={'x^5'}; 
    laurlist1={'1/x','1/y'};
    laurlist2={'1/x^2','1/xy','1/y^2'};
    laurlist3={'1/x^3'};
    sinlist1={'sin x','cos x','sin y','cos y'};
    sinlist2={'sin 2x','cos 2x','sin 2y','cos 2y'};
    absshiftlist1={'|x - 1|', '|x + 1|', '|y - 1|', '|y + 1|'};
elseif n==3
    %varlist={'x''','y''','z'''};
    varlist = {'x^{\prime}', 'y^{\prime}', 'z^{\prime}'};
    xlist = {'$x(t)$','$y(t)$','$z(t)$'};
    dxlist = {'$x^{\prime}(t)$','$y^{\prime}(t)$','$z^{\prime}(t)$'};
    xnoisylist = {'$x_{noisy}(t)$','$y_{noisy}(t)$','$z_{noisy}(t)$'};
    dxnoisylist = {'$x_{noisy}^{\prime}(t)$','$y_{noisy}^{\prime}(t)$',...
        '$z_{noisy}^{\prime}(t)$'};
    diffxlist = {'$x_{noisy}(t)-x(t)$', ...
        '$y_{noisy}(t)-y(t)$', '$z_{noisy}(t)-z(t)$'};
    diffdxlist = {'$x_{noisy}^{\prime}(t)-x^{\prime}(t)$', ...
        '$y_{noisy}^{\prime}(t)-y^{\prime}(t)$', ...
        '$z_{noisy}^{\prime}(t)-z^{\prime}(t)$'};
    xestlist = {'$\hat{x}(t)$','$\hat{y}(t)$','$\hat{z}(t)$'};
    dxestlist = {'$\hat{x}^{\prime}(t)$','$\hat{y}^{\prime}(t)$',...
       '$\hat{z}^{\prime}(t)$'}
    diffxestlist = {'$x(t)-\hat{x}(t)$','$y(t)-\hat{y}(t)$','$z(t)-\hat{z}(t)$'};
    diffdxestlist = {'$x^{\prime}(t)-\hat{x}^{\prime}(t)$',...
       '$y^{\prime}(t)-\hat{y}^{\prime}(t)$',...
       '$z^{\prime}(t)-\hat{y}^{\prime}(t)$'}
    %
    lablist1={'x','y','z'};
    lablist2={'x^2','xy','xz', 'y^2', 'yz', 'z^2'};
    lablist3={'x^3','x^2y','x^2z','xy^2','xyz','xz^2','y^3','y^2z','yz^2','z^3'};
    lablist4={'x^4'}; 
    lablist5={'x^5'}; 
    laurlist1={'1/x','1/y','1/z'};
    laurlist2={'1/x^2','1/xy','1/xz', '1/y^2', '1/yz', '1/z^2'};
    laurlist3={'1/x^3'};
    sinlist1={'sin x','cos x','sin y','cos y', 'sin z', 'cos z'};
    sinlist2={'sin 2x','cos 2x','sin 2y','cos 2y', 'sin 2z', 'cos 2z'};
    absshiftlist1={'|x - 1|', '|x + 1|', '|y - 1|', '|y + 1|', '|z - 1|', '|z + 1|'};
elseif n==4
    %varlist={'w''','x''','y''','z'''};
    varlist = {'w^{\prime}', 'x^{\prime}', 'y^{\prime}', 'z^{\prime}'};
    xlist = {'$w(t)$','$x(t)$','$y(t)$','$z(t)$'};
    altxlist = {'$w(t)/200$','$x(t)$','$y(t)$','$z(t)/200$'};  %hardwired
    dxlist = {'$w^{\prime}(t)$','$x^{\prime}(t)$','$y^{\prime}(t)$','$z^{\prime}(t)$'};
    xnoisylist = {'$w_{noisy}(t)$','$x_{noisy}(t)$','$y_{noisy}(t)$','$z_{noisy}(t)$'};
    altxnoisylist = {'$w_{noisy}(t)/200$','$x_{noisy}(t)$','$y_{noisy}(t)$', ...
        '$z_{noisy}(t)/200$'}; %hardwired
    dxnoisylist = {'$w_{noisy}^{\prime}(t)$','$x_{noisy}^{\prime}(t)$',...
        '$y_{noisy}^{\prime}(t)$','$z_{noisy}^{\prime}(t)$'};
    diffxlist = {'$w_{noisy}(t)-w(t)$','$x_{noisy}(t)-x(t)$', ...
        '$y_{noisy}(t)-y(t)$', '$z_{noisy}(t)-z(t)$'};
    altdiffxlist = {'$[w_{noisy}(t)-w(t)]/200$','$x_{noisy}(t)-x(t)$', ...
        '$y_{noisy}(t)-y(t)$', '$[z_{noisy}(t)-z(t)]/200$'};  %hardwired
    diffdxlist = {'$w_{noisy}^{\prime}(t)-x^{\prime}(t)$', ...
        '$x_{noisy}^{\prime}(t)-x^{\prime}(t)$', ...
        '$y_{noisy}^{\prime}(t)-y^{\prime}(t)$', ...
        '$z_{noisy}^{\prime}(t)-z^{\prime}(t)$'};
    xestlist = {'$\hat{w}(t)$','$\hat{x}(t)$','$\hat{y}(t)$','$\hat{z}(t)$'};
    altxestlist = {'$\hat{w}(t)/200$','$\hat{x}(t)$','$\hat{y}(t)$',...
        '$\hat{z}(t)/200$'}; %hardwired
    dxestlist = {'$\hat{w}^{\prime}(t)$','$\hat{x}^{\prime}(t)$',...
        '$\hat{y}^{\prime}(t)$','$\hat{z}^{\prime}(t)$'}
    diffxestlist = {'$w(t)-\hat{w}(t)$','$x(t)-\hat{x}(t)$',...
        '$y(t)-\hat{y}(t)$','$z(t)-\hat{z}(t)$'};
    altdiffxestlist = {'[$w(t)-\hat{w}(t)]/200$','$x(t)-\hat{x}(t)$',...
        '$y(t)-\hat{y}(t)$','$[z(t)-\hat{z}(t)]/200$'}; %hardwired
    diffdxestlist = {'$w^{\prime}(t)-\hat{w}^{\prime}(t)$',...
       '$x^{\prime}(t)-\hat{x}^{\prime}(t)$',...
       '$y^{\prime}(t)-\hat{y}^{\prime}(t)$',...
       '$z^{\prime}(t)-\hat{y}^{\prime}(t)$'}
    %     
    lablist1={'w','x','y','z'};
    lablist2={'w^2','wx','wy','wz','x^2','xy','xz', 'y^2', 'yz', 'z^2'};
    lablist3={'w^3','w^2x','w^2y','w^2z','wx^2','wxy','wxz','wy^2','wyz','wz^2',...       
              'x^3','x^2y','x^2z','xy^2','xyz','xz^2','y^3','y^2z','yz^2','z^3'};  
    lablist4={'w^4'}; 
    lablist5={'w^5'}; 
    laurlist1={'1/w','1/x','1/y','1/z'};
    laurlist2={'1/w^2','1/wx','1/wy','1/wz','1/x^2','1/xy','1/xz', '1/y^2', '1/yz', '1/z^2'};
    laurlist3={'1/w^3'};
    laur1poly2list={'w^2/x', 'w^2/y', 'w^2/z', 'wx/y', 'wx/z', 'wy/x', ...
        'wy/z', 'wz/x', 'wz/y', 'x^2/w', 'x^2/y', 'x^2/z', 'xy/w', 'xy/z',...
        'xz/w', 'xz/y', 'y^2/w', 'y^2/x', 'y^2/z', 'yz/w', 'yz/x', 'z^2/w', ...
        'z^2/x', 'z^2/y'};  %eval in Maple
    sinlist1={'sin w', 'cos w', 'sin x','cos x','sin y','cos y', 'sin z', 'cos z'};
    sinlist2={'sin 2w', 'cos 2w', 'sin 2x','cos 2x','sin 2y','cos 2y', 'sin 2z', 'cos 2z'};
    absshiftlist1={'|w - 1|', '|w + 1|', '|x - 1|', '|x + 1|', '|y - 1|', ...
        '|y + 1|', '|z - 1|', '|z + 1|'};
end

if (polyorder>=0)
    if (inone==0)
        dictlist=[];
    elseif (inone==1)
        dictlist=['1'];
    end
end

if (polyorder>=1)
    dictlist=[dictlist, lablist1];
end
if (polyorder>=2)
    dictlist=[dictlist,lablist2];
end    
if (polyorder>=3)
    dictlist=[dictlist,lablist3];
end  
if (polyorder>=4)
    dictlist=[dictlist,lablist4];
end
if (polyorder>=5)
    dictlist=[dictlist,lablist5];
end
%breaks if >5;

%laurent
if (laurent>=1)
        dictlist=[dictlist,laurlist1];
end
if (laurent>=2)
        dictlist=[dictlist,laurlist2];
end
if (laurent>=3)
        dictlist=[dictlist,laurlist3];
end
%breaks if >3

%shifted absolute values
if (absshift>=1)
        dictlist=[dictlist,absshiftlist1];
end
%breaks if >1; need to automate

%trig
if (sineorder>=1)
        dictlist=[dictlist,sinlist1];
end
if (sineorder>=2)
        dictlist=[dictlist,sinlist2];
end
%breaks if >2; need to automate

%time
if (intime==1)
        dictlist=[dictlist, 't'];
end

%sin cos time
if (sintime==1)
        dictlist=[dictlist, 'sin t', 'cos t'];
end

dictlist





