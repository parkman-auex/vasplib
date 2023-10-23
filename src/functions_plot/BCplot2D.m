function varargout = BCplot2D(BCCAR,Grid,Rm,options)
arguments
    BCCAR=[];
    Grid=[];
    Rm=[];
    options.ColorCut double= 1;
    options.ColorCutMinus double= -1;
    options.ax =  handle([]);
    options.BZ logical= false;
    options.BZmode {mustBeMember(options.BZmode,{'3D','2D'})} = '3D';
    options.BZlabel logical = true; 
    options.shading logical= false;
    options.scalefactor double= 100;
    options.GkOrigin = [0 0 0];
    options.method {mustBeMember(options.method,{'linear','nearest','natural','v4'})} = 'linear';
    options.plotmode {mustBeMember(options.plotmode,{'surf','contour'})} = 'surf'
end
%
if isempty(options.ax)
    Fig =  Figs(1,1);
    ax = Fig.axes(1);
else
    if ishandle(options.ax)
        ax = options.ax;
    else
    end
end
if options.ColorCutMinus == -1
    options.ColorCutMinus = options.ColorCut;
end
%BCCAR = normalize(BCCAR,'scale');
BCCAR = real(BCCAR)*options.scalefactor;
maxBCplus = options.ColorCut*max((BCCAR),[],'all');
maxBCMinus = options.ColorCutMinus*min((BCCAR),[],'all');
maxBC = max(abs(maxBCplus),abs(maxBCMinus));
BCCAR(BCCAR>maxBCplus) = maxBCplus;
BCCAR(BCCAR<maxBCMinus) = maxBCMinus;
%             BCCAR = normalize(BCCAR);
if ~isvector(BCCAR)
    if size(Grid,3) == 1
        sizemesh = size(BCCAR);
        Grid1 = reshape(Grid(:,1),sizemesh);
        Grid2 = reshape(Grid(:,2),sizemesh);
        Grid3 = reshape(Grid(:,3),sizemesh);
    else
        Grid1 = Grid(:,:,1);Grid2 = Grid(:,:,2);Grid3 = Grid(:,:,3);
    end
    %
    h = surf(ax,Grid1,Grid2,Grid3,BCCAR,'EdgeColor','none');
else
    % use griddate to generate surf
    X = Grid(:,1); Y = Grid(:,2); Z = Grid(:,3);
    % Meshing the scatter solution
    % The key pointis that getting the mesh by X and Y.  Firstly, unique X and Y
    Grid1 = unique(X);
    Grid2 = unique(Y);
    Grid3 = unique(Z);
    if length(Grid3) == 1
        % secondly, obtain the mesh by Xcoord and Ycoord
        [meshX, meshY, meshZ] = meshgrid(Grid1, Grid2, Grid3);
        % Finally, 利用griddata来插值，从xyz生成栅格数据
        %最后一个为插值方法，包linear cubic natural nearest和v4等方法
        %v4方法耗时长但较为准确
        meshU = griddata(X,Y,BCCAR,meshX,meshY,options.method);
    else
        % secondly, obtain the mesh by Xcoord and Ycoord
        [meshX,meshY,meshZ] = meshgrid(Grid1,Grid2,Grid3);
        % Finally, 利用griddata来插值，从xyz生成栅格数据
        meshU = griddata(X,Y,Z,BCCAR,meshX,meshY,meshZ,options.method);
    end
    switch options.plotmode
        case 'surf'
            h = surf(ax,meshX,meshY,meshZ,meshU,'EdgeColor','none');
        case 'contour'
            h = contourf(ax,meshX,meshY,meshU);
    end
end
% colormap(ax,redbluecmap);
%colormap(ax,ColorMap.Matplotlib('coolwarm'));
colormap(ax,ColorMap.redblue);
% ref
Gk = (2*pi*eye(3)/Rm).';
Gknorm = norm(Gk(1,:))/10;
offsetX = ones(2,2)*Gknorm;
offsetY = ones(2,2)*Gknorm;
offsetY2 = ones(2,2)*2*Gknorm;
minx = min(Grid1,[],'all');
miny = min(Grid2,[],'all');
minz = min(Grid3,[],'all');
maxx = max(Grid1,[],'all');
maxy = max(Grid2,[],'all');
maxz = max(Grid3,[],'all');
h2 = surf(ax,minx-offsetX,maxy+offsetY ,zeros(2),-ones(2,2)*maxBC,'EdgeColor','none');
h3 = surf(ax,minx-offsetX,maxy+offsetY2,zeros(2),ones(2,2)*maxBC,'EdgeColor','none');
try
    colorbar(ax,...
        'Ticks',[-roundn(maxBC,-2),0,roundn(maxBC,-2)])
catch

end
%view(0,90);
axis(ax,'equal');
if options.shading
    shading(ax,'interp');
end
if options.BZ
    ax = BZplot(Rm,'color','none','ax',ax,'blackwhite',true,'mode',options.BZmode,'label',options.BZlabel,'OriginPoint',options.GkOrigin) ;
    view(ax,2);
end

%-------- return --------
if nargout  == 3
    varargout{1} = ax.Parent;
    varargout{2} = ax;
    varargout{3} = h;
end
if nargout  == 2
    varargout{1} = ax.Parent;
    varargout{2} = ax;
end
if nargout  == 1
    varargout{1} = ax;
end
end
