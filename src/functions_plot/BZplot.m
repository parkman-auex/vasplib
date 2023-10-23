function varargout = BZplot(Rm,options)
arguments
    Rm = POSCAR_readin;
    options.ax = handle([]);
    options.color ='none';
    options.blackwhite = false;
    options.alpha = 0.2;
    options.mode = '3D';
    options.Gk  = false;
    options.light;
    options.OriginPoint = [0 0 0];
    options.Rotation = eye(3);
    options.label = true;
    options.KPOINTS = '';
    options.title = '';
end
import vasplib_plot.*
import vasplib_tool_outer.*;
X0 = options.OriginPoint(1);
Y0 = options.OriginPoint(2);
Z0 = options.OriginPoint(3);
%             if isequal(options.OriginPoint,[0,0,0] )
%                 showk3 = true;
%             else
%                 showk3 = false;
%             end
if isempty(options.ax)
    Fig =  Figs(1,1);
    ax = Fig.axes(1);
else
    ax = options.ax;
end
if isa(Rm,'vasplib')
    Rm = Rm.Rm;
end
if isequal(size(Rm),[2,2])
    Rm = [Rm,[0;0]];
end
color  = options.color ;
alpha  = options.alpha ;
switch options.mode
    case '3D'
        if options.Gk
            Gk = Rm;
        else
            Gk = (eye(3)*2*pi/Rm).';
        end
        if size(Gk,1) <3 || size(Gk,2) <3
            Gk(3,3) = 1;
        end
        %Points = [0,0,0];
        vector = [];
        for i = linspace(-2,2,5)
            for j = linspace(-2,2,5)
                for k = linspace(-2,2,5)
                    vector = [vector ;i j k];
                end
            end
        end
        % rotation first
        Gk = options.Rotation*Gk;
        vectorR = vector*Gk ;
        [~,Line_000]=ismember([0,0,0],vectorR,'rows');
        %vectorR = unique(vectorR);
        [v,c] = voronoin(vectorR );
        V = v(c{Line_000},:);
        %V = V ;
        ax = plotPolyhedron(V,color,alpha,'ax',ax,'OriginPoint', options.OriginPoint);
        %
        scatter3(ax,X0,Y0,Z0,'filled','LineWidth',300,'MarkerEdgeColor','none','MarkerFaceColor','k');
        if strcmp(options.KPOINTS,'')
            text(ax,X0,Y0,Z0,"\Gamma",'FontSize',24);
        else
            [kpoints,~,~,kpoints_name] = KPOINTS_read(options.KPOINTS);
            [kpoints,label_unique] = unique(kpoints,'rows');
            kpoints_name = kpoints_name(label_unique);
            kpoints_name = strrep(kpoints_name,'GAMMA','Γ');
            kpoints_name = strrep(kpoints_name,'Gamma','Γ');
            kpoints_name = strrep(kpoints_name,'G','Γ');
            HighK = kpoints*Gk + options.OriginPoint;
            for i = 1:length(HighK)
                scatter3(ax,HighK(i,1),HighK(i,2),HighK(i,3),'filled','LineWidth',300);
                text(ax,HighK(i,1),HighK(i,2),HighK(i,3),kpoints_name(i),'FontSize',24);
            end
        end
        %
        xlabel('k_x');
        ylabel('k_y');
        zlabel('k_z');
        if options.label
            %
            if options.blackwhite
                quiver3(ax,X0,Y0,Z0,Gk(1,1),Gk(1,2),Gk(1,3),'Color','k','DisplayName','k_1','AutoScale','off');
                quiver3(ax,X0,Y0,Z0,Gk(2,1),Gk(2,2),Gk(2,3),'Color','k','DisplayName','k_2','AutoScale','off');
                quiver3(ax,X0,Y0,Z0,Gk(3,1),Gk(3,2),Gk(3,3),'Color','k','DisplayName','k_3','AutoScale','off');
            else
                quiver3(ax,X0,Y0,Z0,Gk(1,1),Gk(1,2),Gk(1,3),'Color','r','DisplayName','k_1','AutoScale','off');
                quiver3(ax,X0,Y0,Z0,Gk(2,1),Gk(2,2),Gk(2,3),'Color','g','DisplayName','k_2','AutoScale','off');
                quiver3(ax,X0,Y0,Z0,Gk(3,1),Gk(3,2),Gk(3,3),'Color','b','DisplayName','k_3','AutoScale','off');
            end
            %
            % add shift
            Gk = Gk + options.OriginPoint;
            text(ax,Gk(1,1),Gk(1,2),Gk(1,3)," k_1",'FontSize',24);
            text(ax,Gk(2,1),Gk(2,2),Gk(2,3)," k_2",'FontSize',24);
            text(ax,Gk(3,1),Gk(3,2),Gk(3,3)," k_3",'FontSize',24);
        end
    case '2D'
        if size(Rm,1) ==3
            Rm = Rm(1:2,:);
        elseif size(Rm,1) == 2 && size(Rm,2) == 3
        end
        if options.Gk
            Gk = Rm;
        else
            Gk =  (eye(3)*2*pi/Rm).';
        end
        % attention here
        Gk = Gk * options.Rotation;
        % change to local Gk
        % x equals Gk(1,:)
        % z equals Gk(1,:) cross Gk(2,:)
        % y equals z cross x
        Gk_x = Gk(1,:)/norm(Gk(1,:));
        Gk_z = cross(Gk_x,Gk(2,:));
        Gk_z = Gk_z/norm(Gk_z);
        Gk_y = cross(Gk_z,Gk_x);
        Gk_2d_Pmat = [Gk_x;Gk_y];
        %
        Gk_2d = Gk*Gk_2d_Pmat.';
        vector = [];
        for i = linspace(-3,3,7)
            for j = linspace(-3,3,7)
                vector = [vector ;i j];
            end
        end
        vectorR = vector*Gk_2d;
        [~,Line_000]=ismember([0,0],vectorR,'rows');
        [v,c] = voronoin(vectorR );
        V = v(c{Line_000},:);
        %
        % V = V*Gk_2d_Pmat;
        % V = V + options.baseOrigin;
        ax = plotPolyhedron(V,color,alpha,'ax',ax,...
            'OriginPoint', options.OriginPoint,'Orientation',Gk_2d_Pmat);
        scatter3(ax,X0,Y0,Z0,'filled','LineWidth',300);
        text(ax,X0,Y0,Z0,"\Gamma",'FontSize',24);
        xlabel('k_x');
        ylabel('k_y');
        zlabel('k_z');
        if options.label
            quiver3(ax,X0,Y0,Z0,Gk(1,1),Gk(1,2),Gk(1,3),'Color','r','DisplayName','k_1','AutoScale','off');
            quiver3(ax,X0,Y0,Z0,Gk(2,1),Gk(2,2),Gk(2,3),'Color','g','DisplayName','k_2','AutoScale','off');
            %
            Gk = Gk + options.OriginPoint;
            text(ax,Gk(1,1),Gk(1,2),Gk(1,3)," k_1",'FontSize',24);
            text(ax,Gk(2,1),Gk(2,2),Gk(2,3)," k_2",'FontSize',24);
            view(ax,0,90);
        end
    case '1D'

        %annotation(lineType,x,y)
        %  arrow3([0,0,0],Gk(1,:));
        %  arrow3([0,0,0],Gk(2,:));
        %  arrow3([0,0,0],Gk(3,:));
end

    %-------- return --------
    if nargout  == 2
        varargout{1} = ax.Parent;
        varargout{2} = ax;
    end
    if nargout  == 1
        varargout{1} = ax;
    end
end