function varargout = heatplot(DOSCAR,Y,X,kpoints_l,kpoints_name,options)
arguments
    DOSCAR = [];
    Y  = [];
    X  = [];
    kpoints_l=[];
    kpoints_name = [];
    options.mode {mustBeMember(options.mode,{'surf','arc','parm'})}= 'surf';%  Tetris
    options.KPOINTS = 'KPOINTS';
    options.POSCAR = 'POSCAR';
    options.Alpha = 1;
    options.fontname = 'Helvetica';
    options.cmap = parula;
    options.ax =  handle([]);
    options.view = [0,90];
    options.log = true;
    options.enforceRef = true;
    options.square = false;
    options.loglog = false;
    options.xlabel = 'k_x';
    options.ylabel = 'k_z';
    options.zlabel = 'E(eV)';
    options.Accuracy = -1;
    options.fast = true;
    options.fin_dir = 3;
    options.cartisian = false;
end
import vasplib_plot.*;
%
V = DOSCAR;
cmap = options.cmap;
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
if isempty(options.ax)
    Fig =  Figs(1,1);
    ax = Fig.axes(1);
else
    if ishandle(options.ax)
        ax = options.ax;
    else
        
    end
end
%
if options.enforceRef
    vasplibobj = vasplib;
    vasplibobj.Rm = POSCAR_readin(options.POSCAR);
    vasplibobj = kpathgen3D(vasplibobj,options.KPOINTS);
    [~,kpoints_l,kpoints_name] = kpath_information(vasplibobj);
end
switch options.mode
    case 'surf'
        if isempty(X)
            vasplibobj = vasplib;
            vasplibobj.Rm = POSCAR_readin(options.POSCAR);
            vasplibobj = kpathgen3D(vasplibobj,options.KPOINTS);
            [klist_l,kpoints_l,kpoints_name] = kpath_information(vasplibobj);
            X =klist_l;
        end
        if length(Y) == 3
            Y = linspace(Y(1),Y(2),Y(3));
        end
        if isvector(X) && isvector(Y)
            [X_grid,Y_grid] = meshgrid(X,Y);
        else
            X_grid = X;
            Y_grid = Y;
        end
        colormap(ax,cmap);
        h = surf(ax,X_grid,Y_grid,V,'FaceColor','interp','EdgeColor','none');
        Xcut = [X(1),X(end)];
        Ecut = [Y(1),Y(end)];
        colormap(ax,cmap);
        ax = set_reference(kpoints_l,kpoints_name,Xcut,Ecut,...
            'ax',ax,...
            'fontname',options.fontname ...
            );
    case 'arc'
        if length(X) == 2
            vasplibobj = vasplib;
            vasplibobj.Rm = POSCAR_readin(options.POSCAR);
            kfermiarc = Y;
            kmesh = X;
            [klist,klist1,klist2]=kmesh3D(vasplibobj,kmesh,kfermiarc,'fermiarc');
            if options.cartisian
                %options.output = 'refined';
                klist = klist*vasplibobj.Gk;
                switch options.fin_dir
                    case 1
                        klist1 =klist(:,2);
                        klist2 =klist(:,3);
                    case 2
                        klist1 =klist(:,1);
                        klist2 =klist(:,3);
                    case 3
                        klist1 =klist(:,1);
                        klist2 =klist(:,2);
                end
                klist1 = reshape(klist1,kmesh);
                klist2 = reshape(klist2,kmesh);
            else
                switch options.fin_dir
                    case 1
                        klist1 =klist1(:,2);
                        klist2 =klist2(:,3);
                    case 2
                        klist1 =klist1(:,1);
                        klist2 =klist2(:,3);
                    case 3
                        klist1 =klist1(:,1);
                        klist2 =klist2(:,2);
                end
            end
            X = klist2;
            Y = klist1;
        end
        if isvector(X) && isvector(Y)
            [X_grid,Y_grid] = meshgrid(X,Y);
            Xcut = [X(1),X(end)];
            Ecut = [Y(1),Y(end)];
        elseif ~isvector(X) && ~isvector(Y)
            X_grid = X; Y_grid = Y;
            Xcut = [X(1,1),X(1,end)];
            Ecut = [Y(1,end),Y(end,end)];
        else
            warning('wrong input\n');
            return;
        end
        colormap(ax,cmap);
        h = surf(ax,X_grid,Y_grid,V,'FaceColor','interp','EdgeColor','none');
        xlim(ax,Xcut);
        ylim(ax,Ecut);
        xlabel(ax,kpoints_l);
        ylabel(ax,kpoints_name);
        colormap(ax,cmap);
        %axis(ax,'square');
        %[fig,ax] = set_reference(fig,ax,'band',Xcut,Ecut,kpoints_l,kpoints_name);
    case 'parm'
        [X_grid,Y_grid] = meshgrid(X,Y);
        colormap(ax,cmap);
        %contourf(ax,X_grid,Y_grid,DOSCAR,'LineStyle','none','LevelStep',10,'Fill','on');
        h = surf(ax,X_grid,Y_grid,DOSCAR,'FaceColor','interp','EdgeColor','none');
        %surf(ax,X_grid,Y_grid,DOSCAR);
        Xcut = [X(1),X(end)];
        Ecut = [Y(1),Y(end)];
        colormap(ax,cmap);
        ax = set_reference([],[],Xcut,Ecut,...
            'ax',ax,...
            'fontname',options.fontname ...
            );

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
