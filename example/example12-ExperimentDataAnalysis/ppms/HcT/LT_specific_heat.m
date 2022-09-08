function [FittingData,Hc_brandnew,fig,ax]=LT_specific_heat(Hc,T,Data_density,titlename,C_sch_input,T_range_extend,Low_value,High_value,Start_point)
%
import vasplib_experiment.*
% -------------- check --------------
[fig,ax] = creat_figure();
hold(ax,'all');
titlename1  = strcat(titlename,'Schottky -  fit ');
title(ax,titlename1);
nHc = length(Hc);
nT = length(T);
% [minT,locT1] = min(T);
% [Tmax,~] = max(T);
Hc_max = max(Hc);
%  -------------- C_sch --------------
Sample_num = C_sch_input.n ;       % mol 
%R = C_sch_input.R  ;               % J/(K * mol)                    
DegenRate = C_sch_input.DegenRate ;
R = 8.31446261815324  ;          % J/(mol*K)
kB = 8.617333262145e-2;            % meV/K
Theta_begin = C_sch_input.Delta/kB;      % K  
Csh_constant = Sample_num*DegenRate*R;
fprintf('The init Theta_begin is %f K\n',Theta_begin);
% 1 eV= k T   ; 1 meV =  k' T ;
% k = 8.6173324(78)e-5    eV/K
%  -------------- check --------------
if nT ~= nHc
    error('length of M and T different');
end
Deliter = nT/Data_density;
maker_idx = [round(1:Deliter:nT),nT];
% ----------------
if nargin<6
    fig2 = figure();
    ax2 = gca();
    loglog(ax2,T,Hc);
    xlabel(ax2,'Temperature (K)');
    ylabel(ax2,'C_p/J (mol^-1K^-1)');
    fprintf('No recommend T_range, the loglog plot tells you should select a range which is nearly a straight line\n');
    % 
    fprintf('just give T range 1 (low T):\n');
    Tmin = input('The T min is \n:');
    Tmax = input('The T max is \n:');
    T_range_extend(1,:) = [Tmin,Tmax];
    %
     fprintf('just give T range 2 (Mid T):\n');
    Tmin = input('The T min is \n:');
    Tmax = input('The T max is \n:');
    T_range_extend(2,:) = [Tmin,Tmax];
    %
    fprintf('just give T range 3 (High T):\n');
    Tmin = input('The T min is \n:');
    Tmax = input('The T max is \n:');
    T_range_extend(3,:) = [Tmin,Tmax];
    delete(fig2);
    delete(ax2);
end


% ------- reshape data -------
% low
Findlabel = find(T>=T_range_extend(1,1)&T<=T_range_extend(1,2));
Hc_low = Hc(Findlabel );
T_low = T(Findlabel );
% mid
Findlabel = find(T>=T_range_extend(1,1)&T<=T_range_extend(1,2));
Hc_mid = Hc(Findlabel );
T_mid = T(Findlabel );
% high
Findlabel = find(T>=T_range_extend(1,1)&T<=T_range_extend(1,2));
Hc_high = Hc(Findlabel );
T_high = T(Findlabel );
%% fit low - mid - high
if nargin < 9
    Start_point = [1e-4 1e-5 1e-8 Theta_begin];
end
% use lsqcurvefit

xdata = T_low';
ydata = Hc_low';

fun = @(fitdata,xdata)fitdata(1)*xdata+fitdata(2)*xdata.^3+fitdata(3)*xdata.^5+...
                      Csh_constant*(fitdata(4)*xdata.^-1).^2.*...
                      exp(fitdata(4)*xdata.^-1).*(1+DegenRate*exp(fitdata(4)*xdata.^-1)).^-2;
                  
x0 = Start_point;
lb = Low_value;
ub = High_value;
fitdata = lsqcurvefit(fun,x0,xdata,ydata,lb,ub);


% fittype_string = '0 * gamma * x + B3 * x^3 + B5 * x^5 + '...
%                  +string(Csh_constant)+'* (Theta/x)^2 * exp(Theta/x) /(1 + '+string(DegenRate)+'*exp(Theta/x))^2';
% p = fittype(fittype_string);
% 
% f=fit(T_low,Hc_low,p,'Lower',Low_value,'Upper',High_value,'Start',Start_point);

Tmax = max(T_low);
% ---------------------------
% FittingData.gamma = f.gamma;
% FittingData.B3 = f.B3;
% FittingData.B5 = f.B5;
% FittingData.Theta = f.Theta;
FittingData.gamma = fitdata(1);
FittingData.B3 =fitdata(2);
FittingData.B5 =fitdata(3);
FittingData.Theta = fitdata(4);

gamma = FittingData.gamma;
B3 = FittingData.B3;
B5 = FittingData.B5;
Theta = FittingData.Theta;
fprintf('The fit Theta is %f K\n',Theta);
% ---------------------------
TheName = "Hc = "+string(gamma)+" T +"+...
                  string(B3)+" T^3 +"+...
                  string(B5)+" T^5 +"+...
                  "Csh";

legend(ax);

Xmin = 0; Xmax = Tmax + 10;
Ymin = 0; Ymax = Hc_max;
set(ax,'xlim',[Xmin,Xmax]);%y坐标轴范围
set(ax,'ylim',[Ymin,Ymax]);%y坐标轴范围

CheckT  = linspace(0,Tmax,100);
Csh     = Csh_constant*(Theta*CheckT.^-1).^2 .* exp(Theta*CheckT.^-1).*(1 + DegenRate*exp(Theta*CheckT.^-1)).^-2;
CheckHc = CheckT*gamma + CheckT.^3 *B3 + CheckT.^5 *B5 +Csh;

%% plot last
%yyaxis(ax,'left');
plot(ax,T,Hc ,'DisplayName','Hc-T',...
    'MarkerFaceColor',[0 0.447058826684952 0.74117648601532],...
    'MarkerEdgeColor',[0 0.447058826684952 0.74117648601532],...
    'MarkerSize',6,...
    'Marker','o',...
    'LineStyle','none');
% give xylabel
xlabel(ax,'Temperature (K)');
ylabel(ax,'C_p/J (mol^-1K^-1)');

plot(ax,CheckT,CheckHc ,'DisplayName',TheName,'LineWidth',3,'Color','r');
plot(ax,CheckT,Csh ,'DisplayName','Csh','LineWidth',3,'Color','y');

set(ax,'xlim',[0,max(T_low)]);%y坐标轴范围
set(ax,'ylim',[0,max(Hc_low)]);
%% caculate
for i =1:nT
    Csh_one     = Csh_constant*(Theta/T(i))^2 * exp(Theta/T(i)) /(1 + DegenRate*exp(Theta/T(i)))^2;
    CheckHc_one = T(i)*gamma + T(i)^3 *B3 + T(i)^5 *B5 +Csh_one;
    Hc_brandnew(i) = CheckHc_one  ;
end

end


