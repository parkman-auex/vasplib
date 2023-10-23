function varargout = quiverplot3(R,A,color,displayname,ax,scale)
%--------  init  --------
import vasplib_plot.*
if nargin <6
    scale = 0.5;
end
if nargin < 5
    Fig =  Figs(1,1);
    ax = Fig.axes(1);
end
if nargin < 4
displayname = 'vector field';
end
if nargin < 3
color = 'b';
end
X = R(:,1);
Y = R(:,2);
Z = R(:,3);
U = A(:,1);
V = A(:,2);
W = A(:,3);
quiver3(ax,X,Y,Z,U,V,W,scale,'Displayname',displayname,'Color',color,'lineWidth',1);
view(45,30);
    %-------- return --------
    if nargout  == 2
        varargout{1} = ax.Parent;
        varargout{2} = ax;
    end
    if nargout  == 1
        varargout{1} = ax;
    end

end