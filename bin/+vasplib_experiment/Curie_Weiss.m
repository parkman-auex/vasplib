function [FittingData,fig2,ax2]=Curie_Weiss(M,H,T,num_sample,Range,Data_density,titlename)
% 
import vasplib_experiment.*
import vasplib_tool.*
% -----init
[fig,ax] = creat_figure();
[fig2,ax2] = creat_figure();
titlename2  = strcat(titlename,'- Curie Weiss fit with \chi_0');
titlename1  = strcat(titlename,'- Curie Weiss fit ignore \chi_0');
title(ax,titlename1);
title(ax2,titlename2);
X = (M./H)/num_sample;
Xm1 = num_sample*H./M;
nX = length(Xm1);
nT = length(T);
% [minT,locT1] = min(T);
 [Tmax,~] = max(T);

% --check 
if nT ~= nX
    error('length of M and T different');
end

if nargin<5
    Range = Range_choose();
end

% re shape data
Findlabel = find(T>=Range(1)&T<=Range(2));
Xm1_new = Xm1(Findlabel );
X_new = X(Findlabel );
Xm1_max= max(Xm1);
X_max = max(X);
T_new = T(Findlabel );

Deliter = nT/Data_density;
maker_idx = [round(1:Deliter:nT),nT];



%yyaxis(ax,'right');
% plot(ax,T,X ,'DisplayName','\chi^{-1}-T',...
%     'MarkerFaceColor',[0 0.447058826684952 0.74117648601532],...
%     'MarkerEdgeColor',[0 0.447058826684952 0.74117648601532],...
%     'MarkerSize',3,...
%     'Marker','o',...
%     'LineStyle','none');
% % give title

% 
%% fit second
F = polyfit(Xm1_new,T_new,1);
Tc = F(2);
C =F(1);
TheName = "T = "+string(C)+" / \chi +"+string(Tc);
TheName = strrep(TheName,'+-','- ');
% coeffnames(fun) % 可以校验欲拟合的系数


p1=fittype('c/(x-b)');
f1=fit(T_new,Xm1_new,p1,'StartPoint',[Tc,C]);
% disp(f);
CheckXm1 = 0:100*Data_density:Xm1_max ; 
CheckT=polyval(F,CheckXm1 ); 
plot(ax,CheckT,CheckXm1 ,'DisplayName',TheName,'LineWidth',3);
%
Xmin = Tc- 10;Xmax = Tmax + 10;
Ymin = 0; Ymax = Xm1_max;
set(ax,'xlim',[Xmin,Xmax]);%y坐标轴范围
set(ax,'ylim',[Ymin,Ymax]);%y坐标轴范围
%line(ax,[0 0],[0 Ymax],'Color','k');
legend(ax);
%ax.XAxisLocation = 'origin';
%ax.YAxisLocation = 'origin';
FittingData(1).Tc = F(2);
FittingData(1).C = F(1);
FittingData(1).fittype = p1;
FittingData(1).fitresult = f1;
FittingData(1).chi0 = 0;
%% fit third

p2=fittype('a + c/(x-b) ');

f2=fit(T_new,X_new,p2,'StartPoint',[0,Tc,C],'Lower',[0,-1000,0]);

FittingData(2).Tc = f2.b;
FittingData(2).C = f2.c;
FittingData(2).fittype = p2;
FittingData(2).fitresult = f2;
FittingData(2).chi0 = f2.a;

TheName = "T = "+string(FittingData(2).C)+" / (\chi-\chi_0) +"+string(FittingData(2).Tc);
TheName = strrep(TheName,'+-','- ');
% CheckT  = T;
% CheckX = feval(f2,CheckT);
% plot(ax2,CheckT,CheckX ,'DisplayName',TheName,'LineWidth',3);
legend(ax2);

Xmin = Tc- 10;Xmax = Tmax + 10;
Ymin = 0; Ymax = Xm1_max;
set(ax2,'xlim',[Xmin,Xmax]);%y坐标轴范围
set(ax2,'ylim',[Ymin,Ymax]);%y坐标轴范围

Xm2 = 1 ./(X-FittingData(2).chi0);
Xm2_max = 0:100*Data_density:Xm2_max ; 
CheckXm2 = 0:100*Data_density:Xm2_max ; 
CheckT= (Xm2- FittingData(2).Tc/FittingData(2).C).* FittingData(2).C;
plot(ax2,CheckT,CheckX ,'DisplayName',TheName,'LineWidth',3);
%% plot last
%yyaxis(ax,'left');
plot(ax,T,Xm1 ,'DisplayName','\chi^{-1}-T',...
    'MarkerIndices',maker_idx,...
    'MarkerFaceColor',[0 0.447058826684952 0.74117648601532],...
    'MarkerEdgeColor',[0 0.447058826684952 0.74117648601532],...
    'MarkerSize',6,...
    'Marker','o',...
    'LineStyle','none');
% give xylabel
xlabel(ax,'Temperature (K)');
ylabel(ax,'\chi^{-1} mol');

plot(ax2,T,Xm2 ,'DisplayName','\chi-T',...
    'MarkerIndices',maker_idx,...
    'MarkerFaceColor',[0 0.447058826684952 0.74117648601532],...
    'MarkerEdgeColor',[0 0.447058826684952 0.74117648601532],...
    'MarkerSize',6,...
    'Marker','o',...
    'LineStyle','none');
% give xylabel
xlabel(ax2,'Temperature (K)');
ylabel(ax2,'{\chi-\chi_0}^{-1} mol');

%% caculate 
FittingData(1).Mag = sqrt(8*FittingData(1).C);
FittingData(2).Mag = sqrt(8*FittingData(2).C);
annotation(fig,'textbox',...
    [0.592071428571429 0.190476190476191 0.227571428571429 0.152380952380952],...
        'FontSize',24,...
    'String','Mag =  '+  string(FittingData(1).Mag) +' \mu_B / f.u',...
    'FitBoxToText','off',...
    'EdgeColor','none');
annotation(fig2,'textbox',...
    [0.592071428571429 0.190476190476191 0.227571428571429 0.152380952380952],...
        'FontSize',24,...
    'String',['Mag =  ' + string(FittingData(2).Mag)+" \mu_B / f.u",...
    "\chi_0 = "+string(FittingData(2).chi0)],...
    'FitBoxToText','off',...
    'EdgeColor','none');

end


function Range = Range_choose(M,H,T)
    %init
    nM = length(M);
    nT = length(T);
    [minT,locT1] = min(T);
    [maxT,locT2] = max(T);
    Tcontrol = 30;
    % Recommand data
    RecoT1 = maxT - Tcontrol;
    RecoT2 = maxT - Tcontrol*2;
    %[~,loc_RecoT1 ] = find(T = );
    %[~,loc_RecoT1 ] = find(T = );
    % print
    fprint('We have a %d length data, \n with %d th minT :%f K; %d th maxT :%f K; \n',...
        nT,locT1,minT,locT2,maxT);   
    fprint('We recommand the range [%f,%f]\n ', RecoT1, RecoT2);
    inputchar = input('however, we can still do a range search fot Curie_Weiss fitting, y for continue, n for jump\n');
    if strcmp(inputchar,'y')
        disp('waiting');
    elseif strcmp(inputchar,'n')
        print('just give T range:\n');
        Tmin = input('The T min is \n:');
        Tmax = input('The T max is \n:');
        Range = [Tmin,Tmax];
    end
    
end