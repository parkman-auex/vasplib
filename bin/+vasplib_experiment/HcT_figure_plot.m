function [fig,ax] = HcT_figure_plot(Hc_cell,T_cell,H_list,Data_density,T_cut,cmap)
if nargin <4
    Data_density = 100;
end
if nargin < 5
    disp('whole range');
    for i = 1:length(T_cell)
        T_cut(i,:) = [min(T_cell{i}),max(T_cell{i})];
    end
end
if nargin < 6
    cmap = hsv;
end
%maker_idx = 1:30:length(x);
%plot(x,y,'-s',)

%---
import vasplib_tool.*;
import vasplib_experiment.*;
%-----check
nHc_cell = length(Hc_cell);
nT_cell =  length(T_cell);
nH_list=length(H_list);

ncmap=length(cmap);
Dcmap = round(ncmap/nH_list);
for i = 1:nH_list
    cmap_new(i,:) =cmap((Dcmap-1)*i+i,:);
end
cmap = cmap_new;
Maker_list = ["o";'^';'square';'pentagram';'hexagram';'diamond';'>';'<'];

if nH_list ~= nHc_cell 
    error('num error,nTcell ~= nH_cell ');
end
if nT_cell ~= nHc_cell 
    error('num error,nHccell ~= nH_cell ');
end


% create figure
[fig,ax] = creat_figure();
xlabel(ax,'T (K)');
ylabel(ax,'C_p/J (mol^-1K^-1)');
 
for i= 1: nH_list
    % init data
    Hc = Hc_cell{i};
    T = T_cell{i};
    color =cmap(i,:);

    % refersh
    Findlabel = find(T>=T_cut(i,1) & T<=T_cut(i,2));
    Hc_new = Hc(Findlabel );
    T_new = T(Findlabel);
    nT_new = length(T_new );
    maker_idx = round(linspace(1,nT_new,Data_density));
    %alpha_values = linspace(0.8,0.8,length(H_new))';
    Name = "H = "+string(H_list(i))+" Ost";
    % ----
%     M_new(end) = NaN;
    %
    plot(ax,T_new,Hc_new,'DisplayName',Name, ...
    'MarkerIndices',maker_idx,...
    'MarkerFaceColor',color,...
    'MarkerEdgeColor',color,...
    'MarkerSize',4,...
    'Marker',Maker_list(mod(i,length(Maker_list)),:),...
    'LineWidth',4,...
    'Color',color);
     clear('Findlabel');
     clear('Hc_new ');
     clear('T_new');
%     patch('EdgeColor',color,...
%     'FaceVertexAlphaData',alpha_values,'AlphaDataMapping','none',...
%     'EdgeAlpha','interp')
end


axis(ax,'auto');

% title
title(ax,'M-H');

%

%line(ax,[0 0],[0 Ymax],'Color','k');
legend(ax);

end


