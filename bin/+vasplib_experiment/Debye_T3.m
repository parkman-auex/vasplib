function [FittingData,Hc_brandnew,fig,ax]=Debye_T3(Hc,T,Data_density,titlename,T_range)
% 
import vasplib_experiment.*
% -----init
[fig,ax] = creat_figure();
hold(ax,'all');
titlename1  = strcat(titlename,'- Dsbye T3 fit ');
title(ax,titlename1);

nHc = length(Hc);
nT = length(T);
% [minT,locT1] = min(T);
 [Tmax,~] = max(T);
Hc_max = max(Hc);
% --check 
if nT ~= nHc
    error('length of M and T different');
end

if nargin<5
     fig2 = figure();
      ax2 = gca();
     loglog(ax2,T,Hc);
    xlabel(ax2,'Temperature (K)');
    ylabel(ax2,'C_p/J (mol^-1K^-1)');
    fprintf('No recommend T_range, the loglog plot tells you should select a range which is nearly a straight line\n');
    fprintf('just give T range:\n');
        Tmin = input('The T min is \n:');
        Tmax = input('The T max is \n:');
        T_range = [Tmin,Tmax];
    delete(fig2);
    delete(ax2);
end

% re shape data
Findlabel = find(T>=T_range(1)&T<=T_range(2));

Hc_new = Hc(Findlabel );
T_new = T(Findlabel );
Deliter = nT/Data_density;
maker_idx = [round(1:Deliter:nT),nT];
%% fit third

p=fittype('a *x^3 ');

f=fit(T_new,Hc_new,p,'Lower',[0]);

FittingData.alpha = f.a;

TheName = "Hc = "+string(FittingData.alpha)+" T^3 ";

legend(ax);

Xmin = 0; Xmax = Tmax + 10;
Ymin = 0; Ymax = Hc_max;
set(ax,'xlim',[Xmin,Xmax]);%y×ø±êÖá·¶Î§
set(ax,'ylim',[Ymin,Ymax]);%y×ø±êÖá·¶Î§

CheckT = linspace(0,Tmax,100);
CheckHc=CheckT.^3 .* FittingData.alpha ;

%% plot last
%yyaxis(ax,'left');
plot(ax,T,Hc ,'DisplayName','Hc-T',...
    'MarkerIndices',maker_idx,...
    'MarkerFaceColor',[0 0.447058826684952 0.74117648601532],...
    'MarkerEdgeColor',[0 0.447058826684952 0.74117648601532],...
    'MarkerSize',6,...
    'Marker','o',...
    'LineStyle','none');
% give xylabel
xlabel(ax,'Temperature (K)');
ylabel(ax,'C_p/J (mol^-1K^-1)');
plot(ax,CheckT,CheckHc ,'DisplayName',TheName,'LineWidth',3,'Color','r');
%% caculate 
for i =1:nT
    Hc_brandnew(i) = Hc(i)-FittingData.alpha*T(i)^3;
end

end


