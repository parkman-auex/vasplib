%% 
function [fig,ax]  = POSCAR_plot(Rm,sites,Atom_name,Atom_num,options)
arguments
    Rm =[];
    sites=[];
    Atom_name =[];
    Atom_num =[];
    options.fig = figure('WindowState','maximized');
    options.ax = [];
    options.boundary = [0,1;0,1;0,1];
    options.box = [-2,-2,-1;3,3,2];
    options.title = 'stucture';
    options.vectorL = [0,0,0];
    options.TwoD = false;
    options.scale = 1;
    options.atomscale = 0.3;
    options.fast = false;
    options.view = [30,60];
    options.tight = true;
end
% 
scale =options.scale;
%
if nargin <1
    [Rm,sites,Atom_name,Atom_num,~]=POSCAR_readin('POSCAR');
end
% 
elements_table=element_table();
% 创建 figure
fig = options.fig;
% 创建 axes
if isempty(options.ax)
    ax = axes(fig);
else
    ax = options.ax;
end

% 
hold(ax,'on');
% axis(ax,'equal');
axis image on;
grid on;
%% plot box
Rm = double(Rm)*scale;
if options.TwoD
   Rm(3,:) = [0,0,0.1];
end
XYZ = diag(Rm);
if ~options.tight
    scatter3(ax,XYZ(1)*options.box (1,1),...
        XYZ(2)*options.box (1,2),...
        XYZ(3)*options.box (1,3));
    scatter3(ax,XYZ(1)*options.box (2,1),...
        XYZ(2)*options.box (2,2),...
        XYZ(3)*options.box (2,3));
    plotAxis(Rm,'ax',ax,'W',0.5,'H',1,'box',options.box );
end
%% plot Lattice 
plotBox(Rm,'ax',ax,'vectorL',options.vectorL);

if isa(sites,'struct')
    position   = [[sites.rc1].',[sites.rc2].',[sites.rc3].'];
else
    position = sites;
end
if isa(Atom_name,'char') || isstring(Atom_name)
    tempseq=0;
    for i=1:length(Atom_num)
        for j=1:Atom_num(i)
            tempseq    = tempseq+1;
            sites_name = Atom_name(i);
            AtomList(tempseq,1) = find( elements_table.element_name==sites_name);
        end
    end
else
    AtomList = Atom_name;
end

%% plot sites
if ~options.fast
    for i=1:length(AtomList)
        %disp(temprows);
        markercolor = table2array(elements_table(AtomList(i),{'rgb1','rgb2','rgb3'}));
        radius =  table2array(elements_table(AtomList(i),{'element_radius'}))*options.atomscale;
        Atoms = plotAtoms(ax,Rm,'position',position(i,:),'radius',radius,'color', markercolor);
        Atom_handle=handle(Atoms);
        Atom_handle.DisplayName = string(i);
    end
else
    markercolor = ones(length(AtomList),3);
    markersize = ones(length(AtomList),1)*50;
    for i=1:length(AtomList)
        markercolor(i,:) = table2array(elements_table(AtomList(i),{'rgb1','rgb2','rgb3'}));
        markersize(i,:) = table2array(elements_table(AtomList(i),{'element_radius'}))*options.atomscale*50;
    end
    real_position=position*Rm;
    scatter3(ax,real_position(:,1), real_position(:,2), real_position(:,3),markersize, markercolor,'filled');
end 
titlename = options.title;
% for i =1:length(Atom_name)
%     titlename=titlename+Atom_name(i)+Atom_num(i);
% end

title(char(titlename));
%% lighting
l = light(ax);
l.Color = [1 1 1];
l.Style = 'local';
l.Position = XYZ;
view(options.view(1),options.view(2));
axis equal;
end

%% function
function [a,b,c,alpha,beta,gamma]=Rm2abc(Rm)
    a1 =Rm(1,:);a2 =Rm(2,:);a3 =Rm(3,:);
    a=norm(a1);b=norm(a2);c=norm(a3);
end
function plotBox(Rm,options)
    arguments
        Rm;
        options.vectorL = [0,0,0];
        options.ax = gca();
        options.LineSpec ='-';
        options.Color ='k';
        options.LineWidth = 1.5;
    end
    ax = options.ax;
    hold(ax,'on');
    %fig = figure();
    import vasplib_tool_outer.*;
    %box on;
    a1 =Rm(1,:);a2 =Rm(2,:);a3 =Rm(3,:);
    vectorList = double(options.vectorL);
    % 8site
    vertex = [0, 0, 0;...
        [a1];...
        [a2];...
        [a1+a2];...
        [a3];...
        [a1+a3];...
        [a2+a3];...
        [a1+a2+a3]];
    %vertex = vertex ;
    plotVertex(ax,vertex,'LineSpec',options.LineSpec,'Color', options.Color, 'LineWidth', options.LineWidth);
    for i =1:size(vectorList,1)
        if ~isequal( vectorList(i,:),[0,0,0])
            vertex_tmp = vertex + vectorList(i,:)*Rm;
            plotVertex(ax,vertex_tmp,'LineSpec','--','Color', options.Color, 'LineWidth', options.LineWidth/5);
        end
    end
    % 
%     Rm_=-Rm;
%     [x_min, x_max, y_min, y_max, z_min, z_max] = deal(min(Rm_(:,1)), max(Rm(:,1)), ...
%                                                       min(Rm_(:,2)), max(Rm(:,2)), ...
%                                                       min(Rm_(:,3)), max(Rm(:,3)));
%     x_len = (x_max - x_min)/2;
%     y_len = (y_max - y_min)/2;
%     z_len = (z_max - z_min)/2;
%    
%     axis([x_min-x_len x_max+x_len y_min-y_len y_max+y_len z_min-z_len z_max+z_len]);
     
end
function plotVertex(ax,vertex,options)
    arguments
        ax handle;
        vertex ;
        options.LineSpec ='-';
        options.Color ='k';
        options.LineWidth = 1.5;
    end
    % 12
    plotLine(ax,vertex(1,:), vertex(2,:),'LineSpec',options.LineSpec,'Color', options.Color, 'LineWidth', options.LineWidth);
    plotLine(ax,vertex(1,:), vertex(3,:),'LineSpec',options.LineSpec,'Color', options.Color, 'LineWidth', options.LineWidth);
    plotLine(ax,vertex(2,:), vertex(4,:),'LineSpec',options.LineSpec,'Color', options.Color, 'LineWidth', options.LineWidth);
    plotLine(ax,vertex(3,:), vertex(4,:),'LineSpec',options.LineSpec,'Color', options.Color, 'LineWidth', options.LineWidth);
    plotLine(ax,vertex(5,:), vertex(6,:),'LineSpec',options.LineSpec,'Color', options.Color, 'LineWidth', options.LineWidth);
    plotLine(ax,vertex(5,:), vertex(7,:),'LineSpec',options.LineSpec,'Color', options.Color, 'LineWidth', options.LineWidth);
    plotLine(ax,vertex(6,:), vertex(8,:),'LineSpec',options.LineSpec,'Color', options.Color, 'LineWidth', options.LineWidth);
    plotLine(ax,vertex(7,:), vertex(8,:),'LineSpec',options.LineSpec,'Color', options.Color, 'LineWidth', options.LineWidth);
    plotLine(ax,vertex(1,:), vertex(5,:),'LineSpec',options.LineSpec,'Color', options.Color, 'LineWidth', options.LineWidth);
    plotLine(ax,vertex(2,:), vertex(6,:),'LineSpec',options.LineSpec,'Color', options.Color, 'LineWidth', options.LineWidth);
    plotLine(ax,vertex(3,:), vertex(7,:),'LineSpec',options.LineSpec,'Color', options.Color, 'LineWidth', options.LineWidth);
    plotLine(ax,vertex(4,:), vertex(8,:),'LineSpec',options.LineSpec,'Color', options.Color, 'LineWidth', options.LineWidth);
end
function plotMoreBox(Rm,vectorL,options)
    arguments
        Rm;
        vectorL ;
        options.ax = gca();
    end
    import vasplib_tool_outer.*;
    ax = options.ax;
    hold(ax,'on');
    %fig = figure();
    for i =1:size(vectorL,1)
    end

end
function Axis = plotAxis(Rm,options)
arguments
    Rm
    options.ax = gca();
    options.W = 1;
    options.H = 2;
    options.box = [-2,-2,-1;3,3,2];
end
import vasplib_tool_outer.*;
W = options.W;
H = options.H;
ax = options.ax;
hold(ax,'on');
%fig = figure();

%box on;
XYZ = diag(Rm).' ;
XYZ(1) = XYZ(1)*options.box(1,1);
XYZ(2) = XYZ(2)*options.box(1,2);
XYZ(3) = XYZ(3)*options.box(1,3);
a1 =Rm(1,:);a2 =Rm(2,:);a3 =Rm(3,:);
if length(W) == 3
    Arrow = Plot_arrow(XYZ,XYZ+a1*0.5...
        ,'ax',ax...
        ,'LineWidth',1 ...
        ,'W',W(1) ...
        ,'H',H(1) ...
        ,'DisplayName','a_1'...
        );
    Axis(1,:) = Arrow.';
    Arrow = Plot_arrow(XYZ,XYZ+a2*0.5...
        ,'ax',ax...
        ,'LineWidth',1 ...
        ,'W',W(2) ...
        ,'H',H(2) ...
        ,'DisplayName','a_2'...
        );
    Axis(2,:) = Arrow.';
        Arrow = Plot_arrow(XYZ,XYZ+a3*0.5...
        ,'ax',ax...
        ,'LineWidth',1 ...
        ,'W',W(3) ...
        ,'H',H(3) ...
        ,'DisplayName','a_3'...
        );
    Axis(3,:) = Arrow.';
else
    Arrow = Plot_arrow(XYZ,XYZ+a1*0.5...
        ,'ax',ax...
        ,'LineWidth',1 ...
        ,'W',W(1) ...
        ,'H',H(1) ...
        ,'DisplayName','a_1'...
        );
    Axis(1,:) = Arrow.';
    Arrow = Plot_arrow(XYZ,XYZ+a2*0.5...
        ,'ax',ax...
        ,'LineWidth',1 ...
        ,'W',W(1) ...
        ,'H',H(1) ...
        ,'DisplayName','a_2'...
        );
    Axis(2,:) = Arrow.';
        Arrow = Plot_arrow(XYZ,XYZ+a3*0.5...
        ,'ax',ax...
        ,'LineWidth',1 ...
        ,'W',W(1) ...
        ,'H',H(1) ...
        ,'DisplayName','a_3'...
        );
    Axis(3,:) = Arrow.';
end
end
function Atoms = plotAtoms(ax,Rm,options)
arguments
    ax handle;
    Rm ;
    options.position = [0,0,0];
    options.radius = 1;
    options.color = 'g';
    options.n = 100;
end
    % a1 =Rm(1,:);a2 =Rm(2,:);a3 =Rm(3,:);
    real_position=options.position*Rm;
%     cur_point = real_position;
    for i  =1 :size(real_position,1)
        Atoms(i) = plotBall(ax,...
            'position',real_position,...
            'radius',options.radius,...
            'color',options.color,...
            'n',options.n...
            );
    end
%     plot3(ax,cur_point(:,1), cur_point(:,2), cur_point(:,3), 'ok', 'linewidth', 1.5, 'markersize', markersize, 'markerfacecolor', markercolor)

end
function Ball = plotBall(ax,options)
arguments
    ax handle
    options.position = [0,0,0];
    options.radius = 1;
    options.color = 'g';
    options.n = 100;
end
[X,Y,Z]  = sphere(options.n);
r = options.radius;
X2 = X * r;
Y2 = Y * r;
Z2 = Z * r;
X1 = options.position(:,1);
Y1 = options.position(:,2);
Z1 = options.position(:,3);
% C = zeros();
Ball = surf(ax,X2+X1,Y2+Y1,Z2+Z1);
Ball.FaceColor = options.color;
Ball.EdgeColor = 'none';
end
function plotLine(ax,x1, x2,options)
arguments
    ax handle;
    x1 ;
    x2 ;
    options.LineSpec ='-'; 
    options.Color ='k'; 
    options.LineWidth = 1.5;
end
%   '-o','Color','b','MarkerSize',10,'MarkerFaceColor','#D9FFFF'
    plot3(ax,[x1(1) x2(1)], [x1(2), x2(2)], [x1(3), x2(3)],  options.LineSpec,'Color', options.Color, 'LineWidth', options.LineWidth);
end
%%
function elements_table=element_table()
%% Initialize variables.
filename = 'elements.dat';
delimiter = ' ';

%% Format for each line of text:
%   column1: text (%s)
%	column2: double (%f)
%   column3: double (%f)
%	column4: double (%f)
%   column5: double (%f)
%	column6: double (%f)
%   column7: double (%f)
%	column8: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%s%f%f%f%f%f%f%f%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'TextType', 'string',  'ReturnOnError', false);

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% Create output variable
elements_table = table(dataArray{1:end-1}, 'VariableNames', {'element_name','element_value','element_radius','var1','var2','rgb1','rgb2','rgb3'});


%% Clear temporary variables
clearvars filename delimiter formatSpec fileID dataArray ans;
%%
end