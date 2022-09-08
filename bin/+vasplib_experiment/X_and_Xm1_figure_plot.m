function [fig,ax] = X_and_Xm1_figure_plot(M1,H1,T1,M2,H2,T2,num_sample,Trange,Data_density,titlename)
%---
import vasplib_tool.*;
import vasplib_experiment.*;
%----------------
if nargin < 8
    disp('whole range');
    Trange = [max([min(T1),min(T2)]),min([max(T1),max(T2)])];
end
if nargin < 9 
    Data_density = 100;
end
% init data
Tmax = Trange(2);
X_ZFC = (H1./M1)/num_sample;
X_FC = (H2./M2)/num_sample;
Xm1_ZFC = num_sample*M1./H1;
Xm1_FC = num_sample*M2./H2;

Xm1_ZFC_max= max(Xm1_ZFC);
Xm1_FC_max= max(Xm1_ZFC);
Xm1_max = max([Xm1_FC_max,Xm1_ZFC_max]);

X_ZFC_max= max(X_ZFC);
X_FC_max= max(X_ZFC);
X_max = max([X_FC_max,X_ZFC_max]);

% refersh
Findlabel1 = find(T1>=Trange(1)&T1<=Trange(2));
Findlabel2 = find(T2>=Trange(1)&T2<=Trange(2));
Xm1_ZFC_new = Xm1_ZFC(Findlabel1 );
Xm1_FC_new = Xm1_FC(Findlabel2 );
X_ZFC_new = X_ZFC(Findlabel1 );
X_FC_new = X_FC(Findlabel2 );
T1_new = T1(Findlabel1);
T2_new = T2(Findlabel2 );
maker_idx1 = round(linspace(1, length(T1_new),Data_density));
maker_idx2 = round(linspace(1, length(T2_new),Data_density));
%
Xmin = 0; Xmax = Tmax + 10;
Ymin = 0; Ymax = Xm1_max;

% create figure
[fig,ax] = creat_figure();
yyaxis(ax,'right');
plot(ax,T1_new,Xm1_ZFC_new,'DisplayName','\chi^{-1} ZFC',...
    'MarkerIndices',maker_idx1,...
    'MarkerFaceColor',[0.850980401039124 0.325490206480026 0.0980392172932625],...
    'MarkerEdgeColor',[0.850980401039124 0.325490206480026 0.0980392172932625],...
    'MarkerSize',10,...
    'Marker','o',...
    'LineWidth',0.1,...
    'Color',[0.9 0.9 0.9]);
plot(ax,T2_new,Xm1_FC_new,'DisplayName','\chi^{-1} FC',...
    'MarkerIndices',maker_idx2,...
    'MarkerFaceColor',[1 0.843137264251709 0],...
    'MarkerEdgeColor',[1 0.843137264251709 0],...
    'MarkerSize',10,...
    'Marker','^',...
    'LineWidth',0.1,...
    'Color',[0.9 0.9 0.9]);
xlabel(ax,'Temperature (K)');
ylabel(ax,'\chi^{-1} mol');
set(ax,'xlim',[Xmin,Xmax]);%y×ø±êÖá·¶Î§
set(ax,'ylim',[Ymin,Ymax]);%y×ø±êÖá·¶Î§

% 
Ymin = 0; Ymax = X_max;
yyaxis(ax,'left');
set(ax,'ylim',[Ymin,Ymax]);
plot(ax,T1_new,X_ZFC_new,'DisplayName','\chi ZFC',...
    'MarkerIndices',maker_idx1,...
    'MarkerFaceColor',[0 0.439215689897537 0.737254917621613],...
    'MarkerEdgeColor',[0 0.439215689897537 0.737254917621613],...
    'MarkerSize',10,...
    'Marker','square',...
    'LineWidth',0.1,...
    'Color',[0.9 0.9 0.9]);
plot(ax,T2_new,X_FC_new,'DisplayName','\chi FC','LineWidth',0.1,...
    'MarkerIndices',maker_idx2,...
    'MarkerFaceColor',[0 0.498039215803146 0],...
    'MarkerEdgeColor',[0 0.498039215803146 0],...
    'MarkerSize',10,...
    'Marker','hexagram',...
    'LineWidth',0.1,...
    'Color',[0.9 0.9 0.9]);
ylabel(ax,'\chi mol^{-1}');

% title
title(ax,strcat(titlename,': \chi^{-1} -T & \chi -T'));

%

%line(ax,[0 0],[0 Ymax],'Color','k');
legend(ax);

end