function varargout = bandplot3d(EIGENCAR_3D,klistX,klistY,options)
arguments
    EIGENCAR_3D = [];
    klistX double = [];
    klistY double = [];
    options.refined = true;
    options.Ecut = [-3,3];
    options.WEIGHTCAR = [];
    options.ax =  handle([]);
    options.fontname = 'Helvetica';
    options.KPOINTS = 'KPOINTS';
    options.POSCAR = 'POSCAR';
    options.cmap =  turbo;
    options.density = 1;
    options.title = '';
    options.klist_l = [];
    options.kpoints_l = [];
    options.kpoints_name = [];
    options.xlabel='k_1';
    options.ylabel='k_2';
    options.zlabel='E(eV)';
    options.LineSpec = '-';
    options.LineWidth = 1;
    options.MarkerSize = 3;
    options.EdgeColor = 'none';
    options.FaceAlpha= 0.8;
    options.view = [0,90];
end
if strcmp(options.title,'')
    [~,titlestring]=fileparts(pwd);
    titlestring = strrep(titlestring,'_','-');
else
    titlestring  = options.title;
end
if isempty(options.ax)
    figs = create_figure();
    ax = figs.axes(1);
else
    ax = options.ax;
end
if isempty(options.WEIGHTCAR)
    WEIGHTCAR = abs(EIGENCAR_3D);
else
    WEIGHTCAR = options.WEIGHTCAR;
end
hold(ax,'on');
if options.refined
    for i = 1:size(EIGENCAR_3D,3)
        h = surf(ax,klistX,klistY,EIGENCAR_3D(:,:,i),WEIGHTCAR(:,:,i));
        set(h,'EdgeColor',options.EdgeColor,'FaceAlpha',options.FaceAlpha);
    end
end
%
zlim(ax,options.Ecut);
xlim(ax,[min(min(klistX)),max(max(klistX))]);
ylim(ax,[min(min(klistY)),max(max(klistY))]);
xlabel(ax,options.xlabel);
ylabel(ax,options.ylabel);
zlabel(ax,options.zlabel);
colormap(ax,options.cmap);
view(ax,options.view(1),options.view(2));
%
%-------- return --------
if nargout  == 2
    varargout{1} = ax.Parent;
    varargout{2} = ax;
end
if nargout  == 1
    varargout{1} = ax;
end

end
