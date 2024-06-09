% Copyright 2019-2024, All Rights Reserved
% Code by Robert Niven
% to automate plotting of dynamical system and derivatives

%must come after setlabels.m, for default colours

%variable plots, any n
figure,plot(t,X),...
    xlabel('t','FontSize',18),ylabel('Variable','FontSize',18),...
    legend(Xlist,'Interpreter','latex','FontSize',16),...
    title({fighead; paramhead},'FontSize',18)
    ax = gca; ax.XAxis.FontSize = 16; ax.YAxis.FontSize = 16;
    fignum=fignum+1; 
    saveas(gcf,fullfile(pwd, dirname,strcat('Fig',int2str(fignum))),'epsc')
if withnoise==1
    figure,plot(t,diffX),xlabel('t','FontSize',18),...
        ylabel('Variable Error','FontSize',18),...
        legend(diffXlist,'Interpreter','latex','FontSize',16),...
        title({fighead; paramhead},'FontSize',18)
        ax = gca; ax.XAxis.FontSize = 16; ax.YAxis.FontSize = 16;
        fignum=fignum+1; 
        set(gcf, 'renderer', 'painters');
        saveas(gcf,fullfile(pwd, dirname,strcat('Fig',int2str(fignum))),'epsc')
end;

%pause

%variable plots, any n, alt scale
if altplot == 1
    figure,plot(t,altscale.*X),xlabel('t','FontSize',18),...
        ylabel('Variable','FontSize',18),...
        legend(altXlist,'Interpreter','latex','FontSize',16),...
        title({fighead; paramhead},'FontSize',18)
        ax = gca; ax.XAxis.FontSize = 16; ax.YAxis.FontSize = 16;
        fignum=fignum+1; 
        set(gcf, 'renderer', 'painters');
        saveas(gcf,fullfile(pwd, dirname,strcat('Fig',int2str(fignum))),'epsc')
    if withnoise==1
        figure,plot(t,diffX),xlabel('t','FontSize',18),...
            ylabel('Variable Error','FontSize',18),...
            legend(diffXlist,'Interpreter','latex','FontSize',16),...
            title({fighead; paramhead},'FontSize',18)
            ax = gca; ax.XAxis.FontSize = 16; ax.YAxis.FontSize = 16;
            fignum=fignum+1; 
            set(gcf, 'renderer', 'painters');
           saveas(gcf,fullfile(pwd, dirname,strcat('Fig',int2str(fignum))),'epsc')
        figure,plot(t,altscale.*diffX),xlabel('t','FontSize',18),...
            ylabel('Variable Error','FontSize',18),...
            legend(altdiffXlist,'Interpreter','latex','FontSize',16),...
            title({fighead; paramhead},'FontSize',18)
            ax = gca; ax.XAxis.FontSize = 16; ax.YAxis.FontSize = 16;
            fignum=fignum+1; 
            set(gcf, 'renderer', 'painters');
           saveas(gcf,fullfile(pwd, dirname,strcat('Fig',int2str(fignum))),'epsc')
    end;
end;

%variable derivative or Y plots, any o
figure,plot(t,dX),xlabel('t','FontSize',18),ylabel('Derivative','FontSize',18),...
    legend(dXlist,'Interpreter','latex','FontSize',16),...
    title({fighead; fnhead},'FontSize',18)
    ax = gca; ax.XAxis.FontSize = 16; ax.YAxis.FontSize = 16;
    fignum=fignum+1; 
    set(gcf, 'renderer', 'painters');
    saveas(gcf,fullfile(pwd, dirname,strcat('Fig',int2str(fignum))),'epsc')
if withnoise==1
    figure,plot(t,diffdX),xlabel('t','FontSize',18),...
        ylabel('Derivative Error','FontSize',18),...
        legend(diffdXlist,'Interpreter','latex','FontSize',16),...
        title({fighead; fnhead},'FontSize',18)
        ax = gca; ax.XAxis.FontSize = 16; ax.YAxis.FontSize = 16;
        fignum=fignum+1; 
        set(gcf, 'renderer', 'painters');
        saveas(gcf,fullfile(pwd, dirname,strcat('Fig',int2str(fignum))),'epsc')
end;

%variable phase plots, n = 2 or 3 (or more)
if n == 2
    figure,plot(X(:,1),X(:,2)),...
        xlabel(Xlist(1),'Interpreter','latex','FontSize',18),...
        ylabel(Xlist(2),'Interpreter','latex','FontSize',18),...
        title({fighead; paramsphead},'FontSize',18), ...
        ax = gca; ax.XAxis.FontSize = 16; ax.YAxis.FontSize = 16; 
        fignum=fignum+1; 
        set(gcf, 'renderer', 'painters');
        saveas(gcf,fullfile(pwd, dirname,strcat('Fig',int2str(fignum))),'epsc')
elseif n >= 3
   %3D phase plots
   %for >3D, cycles through 3D plots in order of parameters
   for u=1:n-2   
     figure,plot3(X(:,u),X(:,u+1),X(:,u+2)),...
        xlabel(Xlist(u),'Interpreter','latex','FontSize',18),...
        ylabel(Xlist(u+1),'Interpreter','latex','FontSize',18),...
        zlabel(Xlist(u+2),'Interpreter','latex','FontSize',18),...
        %view([-140 20]), ...  %default
        view([-150 10]), ...  %Vance
        title({fighead; paramsphead},'FontSize',18), ...
        box on, grid on
        ax = gca; ax.XAxis.FontSize = 16; ax.YAxis.FontSize = 16; ax.ZAxis.FontSize = 16;
        fignum=fignum+1; 
        set(gcf, 'renderer', 'painters');
        saveas(gcf,fullfile(pwd, dirname,strcat('Fig',int2str(fignum))),'epsc')
   end
end

%derivative phase plots, o = 2 or 3 (or more)
if o == 2
    figure,plot(dX(:,1),dX(:,2)), ...
        xlabel(dXlist(1),'Interpreter','latex','FontSize',18),...
        ylabel(dXlist(2),'Interpreter','latex','FontSize',18),...
        title({fighead; fnsphead},'FontSize',18), ...
        ax = gca; ax.XAxis.FontSize = 16; ax.YAxis.FontSize = 16; ax.ZAxis.FontSize = 16;
        fignum=fignum+1; 
        set(gcf, 'renderer', 'painters');
        saveas(gcf,fullfile(pwd, dirname,strcat('Fig',int2str(fignum))),'epsc')
elseif o >= 3
   %3D phase plots
   %for >3D, cycles through 3D plots in order of parameters
   for u=1:o-2   
     figure,plot3(dX(:,u),dX(:,u+1),dX(:,u+2)), ...
        xlabel(dXlist(u),'Interpreter','latex','FontSize',18),...
        ylabel(dXlist(u+1),'Interpreter','latex','FontSize',18),...
        zlabel(dXlist(u+2),'Interpreter','latex','FontSize',18),...
        view([-140 20]), ...  %default, incl Vance
        title({fighead; fnsphead},'FontSize',18), ...
        box on, grid on
        ax = gca; ax.XAxis.FontSize = 16; ax.YAxis.FontSize = 16; ax.ZAxis.FontSize = 16;
        fignum=fignum+1; 
        set(gcf, 'renderer', 'painters');
        saveas(gcf,fullfile(pwd, dirname,strcat('Fig',int2str(fignum))),'epsc')
   end
end


