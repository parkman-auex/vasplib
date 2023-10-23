function varargout = heatplot3D(DOSCAR,klist1,klist2,Elist,options)
arguments
    DOSCAR
    klist1 {mustBeVector}
    klist2 {mustBeVector}
    Elist {mustBeVector}
    options.xslice = [];
    options.yslice = [];
    options.zslice = [];
    options.style {mustBeMember(options.style,{'Users','Tetris','Tetris_2'})}= 'Users';%  Tetris
    options.Alpha = 1;
    options.EdgeColor = 'none';
    options.method = 'nearest';
    options.ax =  handle([]);
    options.view = [60,60];
    options.log = true;
    options.square = false;
    options.loglog = false;
    options.xlabel = 'k_x';
    options.ylabel = 'k_z';
    options.zlabel = 'E(eV)';
    options.Accuracy = -1;
    options.fast = true;
end
%
%--------  init  --------

%
V = DOSCAR;
%
if options.Accuracy > 0
    V(V<options.Accuracy) = nan;
end
if options.log
    if options.loglog
        V = log(DOSCAR);
        V (V < 0 ) =nan;
        V = log(V);
    else
        V = log(DOSCAR);
    end
else
    %V = (DOSCAR);
end
%
if options.square
    V = sign(V).*(V.^2);
end
%
% style
if strcmp(options.style ,'Users') ...
        && isempty(options.xslice) ...
        && isempty(options.yslice)...
        && isempty(options.zslice)
    style = 'Tetris';
else
    style = options.style;
end
switch  char(style)
    case 'Users'
        xslice = options.xslice;
        yslice = options.yslice;
        zslice = options.zslice;
    case 'Tetris'
        Nx = length(klist1);
        Ny = length(klist2);
        Nz = length(Elist);
        V(1:round((Nx-3)/2),round((Ny+1)/2):Ny,:) = nan;
        V(1:round((Nx-3)/2),1:round((Ny-1)/2),round((Nz+1)/2):Nz) = nan;
        V(round((Nx-1)/2):Nx,round((Ny+1)/2):Ny,round((Nz+1)/2):Nz) = nan;
        if options.fast
            xslice = klist1([1,round((Nx-1)/2),Nx]);
            yslice = klist2([1,round((Ny-1)/2),Ny]);
            zslice = Elist([1,round((Nz-1)/2),Nz]);%Elist;%[]
        else
            xslice = klist1;
            yslice = klist2;
            zslice = Elist;%Elist;%[]
        end
    case 'Tetris_2'
        Nx = length(klist1);
        Ny = length(klist2);
        Nz = length(Elist);
        V(1:round((Nx-1)/2),round((Ny+1)/2):Ny,:) = nan;
        V(1:round((Nx-1)/2),1:round((Ny-1)/2),round((Nz+1)/2):Nz) = nan;
        V(round((Nx+1)/2):Nx,round((Ny+1)/2):Ny,round((Nz+1)/2):Nz) = nan;
        if options.fast
            xslice = klist1([1,round((Nx-1)/2),Nx]);
            yslice = klist2([1,round((Ny-1)/2),Ny]);
            zslice = Elist([1,round((Nz-1)/2),Nz]);%Elist;%[]
        else
            xslice = klist1;
            yslice = klist2;
            zslice = Elist;%Elist;%[]
        end
    otherwise
        xslice = options.xslice;
        yslice = options.yslice;
        zslice = options.zslice;
end
%
if  isempty(options.ax)
    Fig =  Figs(1,1);
    ax = Fig.axes(1);
else
    if ishandle(options.ax)
        ax = options.ax;
    else
    end
end
[X,Y,Z] = meshgrid(klist1,klist2,Elist);

h = slice(ax,X,Y,Z,V,xslice,yslice,zslice,options.method);
set(h,'FaceColor','interp','EdgeColor',options.EdgeColor,'FaceAlpha',options.Alpha );
switch style
    case 'Users'
        view(options.view(1),options.view(2));
    case 'Tetris'
        view(50,27);
    otherwise
        view(50,27);
end
xlabel(ax,options.xlabel);ylabel(ax,options.ylabel);zlabel(ax,options.zlabel);

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
