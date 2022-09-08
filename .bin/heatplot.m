function [fig,ax] = heatplot(DOSCAR,Y,X,kpoints_l,kpoints_name,cmap,mode,fig,ax)
import vasplib_tool.*;
if nargin <7
    mode = 'surf';
end
if nargin < 8
    [fig,ax] = creat_figure();
end
if nargin <6
    cmap = turbo;
end
if nargin <3
    POSCAR_read;
    [~,klist_l,~,kpoints_l,kpoints_name]=kpathgen3D(Rm);
    X =klist_l;
end
fontname = 'Helvetica';
switch mode
    case 'surf'
        if isvector(X)&& isvector(Y)
            [X_grid,Y_grid] = meshgrid(X,Y);
            Xcut = [X(1),X(end)];
            Ecut = [Y(1),Y(end)];
        else
            X_grid = X;
            Y_grid = Y;
            Xcut = [X(1,1),X(end,1)];
            Ecut = [Y(1,1),Y(1,end)];
        end
        colormap(ax,cmap);
        surf(ax,X_grid,Y_grid,DOSCAR,'FaceColor','interp','EdgeColor','none');
        colormap(ax,cmap);
        try
        [fig,ax] = set_reference(kpoints_l,kpoints_name,Xcut,Ecut,...
            'fig',fig,...
            'ax',ax,...
            'fontname',fontname ...
            );
        catch
        end
    case 'arc'
        if isvector(X)&& isvector(Y)
            [X_grid,Y_grid] = meshgrid(X,Y);
            %Xcut = [X(1),X(end)];
            %Ecut = [Y(1),Y(end)];
        else
            X_grid = X;
            Y_grid = Y;
            %Xcut = [X(1,1),X(end,1)];
            %Ecut = [Y(1,1),Y(1,end)];
        end
        colormap(ax,cmap);
        surf(ax,X_grid,Y_grid,DOSCAR,'FaceColor','interp','EdgeColor','none');
        %Xcut = [X(1),X(end)];
        %Ecut = [Y(1),Y(end)];
        xlabel(ax,kpoints_l);
        ylabel(ax,kpoints_name);
        colormap(ax,cmap);
        %[fig,ax] = set_reference(fig,ax,'band',Xcut,Ecut,kpoints_l,kpoints_name);
    case 'parm'
        [X_grid,Y_grid] = meshgrid(X,Y);
        colormap(ax,jet);
        %contourf(ax,X_grid,Y_grid,DOSCAR,'LineStyle','none','LevelStep',10,'Fill','on');
        surf(ax,X_grid,Y_grid,DOSCAR,'FaceColor','interp','EdgeColor','none');
        %surf(ax,X_grid,Y_grid,DOSCAR);
        Xcut = [X(1),X(end)];
        Ecut = [Y(1),Y(end)];
        colormap(ax,cmap);
        [fig,ax] = set_reference([],[],Xcut,Ecut,...
            'fig',fig,...
            'ax',ax,...
            'fontname',fontname ...
            );
        
end
axis(ax,'tight')

end

