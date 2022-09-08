function varargout = quiverplot(R,A,options)
%--------  init  --------
arguments
    R
    A
    options.color = 'b';
    options.displayname = 'vector field';
    options.ax = handle([]);
    options.scale = 0.5;
    options.ShowArrowHead = 'on';
    options.Rm = POSCAR_readin;
    options.BZ logical= false;
    options.BZmode {mustBeMember(options.BZmode,{'3D','2D'})} = '2D';
    options.BZlabel logical = true; 
end
import vasplib_plot.*;
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
scale = options.scale;
displayname = options.displayname;
color = options.color;


X = R(:,1);
Y = R(:,2);
U = A(:,1);
V = A(:,2);
quiver(ax,X,Y,U,V,scale,'Displayname',displayname,'Color',color,'lineWidth',1,'ShowArrowHead',options.ShowArrowHead);


if options.BZ
    ax =vasplib_plot.BZplot(options.Rm,'color','none','ax',ax,'blackwhite',true,'mode',options.BZmode,'label',options.BZlabel) ;
end
view(0,90);
axis(ax,'equal');
%-------- return --------
if nargout  == 2
    varargout{1} = ax.Parent;
    varargout{2} = ax;
end
if nargout  == 1
    varargout{1} = ax;
end

end