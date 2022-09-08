function [fig,ax] = MH_figure_plot(M_cell,H_cell,T_list,Data_density,H_cut,cmap)
if nargin < 5
    disp('whole range');
    for i = 1:length(H_cell)
        H_cut(i,:) = [min(H_cell{i}),max(H_cell{i})];
    end
end
if nargin < 6
    cmap = jet;
end
%maker_idx = 1:30:length(x);
%plot(x,y,'-s',)

%---
import vasplib_tool.*;
import vasplib_experiment.*;
%-----check
 nM_cell = length(M_cell);
nH_cell =  length(H_cell);
nT_list=length(T_list);

ncmap=length(cmap);
Dcmap = round(ncmap/nT_list);
for i = 1:nT_list
    cmap_new(i,:) =cmap((Dcmap-1)*i+i,:);
end
cmap = cmap_new;
Maker_list = ["o";'^';'square';'pentagram';'hexagram';'diamond';'>';'<'];

if nT_list ~= nH_cell 
    error('num error,nTcell ~= nH_cell ');
end
if nM_cell ~= nH_cell 
    error('num error,nMcell ~= nH_cell ');
end


% create figure
[fig,ax] = creat_figure();
xlabel(ax,'H (Oe)');
ylabel(ax,'M (emu.)');
 
for i= 1: nT_list
    % init data
    M = M_cell{i};
    H = H_cell{i};
    color =cmap(i,:);

    % refersh
    Findlabel = find(H>=H_cut(i,1) & H<=H_cut(i,2));
    M_new = M(Findlabel );
    H_new = H(Findlabel);
    nH_new = length(H_new );
    maker_idx = round(linspace(1,nH_new,Data_density));
    %alpha_values = linspace(0.8,0.8,length(H_new))';
    Name = "T = "+string(T_list(i))+"K";
    % ----
%     M_new(end) = NaN;
    %
    plot(ax,H_new,M_new,'DisplayName',Name, ...
    'MarkerIndices',maker_idx,...
    'MarkerFaceColor',color,...
    'MarkerEdgeColor',color,...
    'MarkerSize',10,...
    'Marker',Maker_list(mod(i,length(Maker_list)),:),...
    'LineWidth',0.1,...
    'Color',[0.9 0.9 0.9]);
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


