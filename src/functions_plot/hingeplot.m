function  ax = hingeplot(orb_list,WaveFunc,mode,options)
    arguments
        orb_list = [];
        WaveFunc  = zeros(size(orb_list,1),1);
        mode {mustBeMember(mode,{'hex','square'})}= 'square';
        options.ax =  handle([]);
        options.Rm = [];
        options.h = 3;
        options.lnum = 100;
        options.hnum = 2;
        options.POSCAR = 'POSCAR';
        options.WaveMin = 1e-3;
        options.rmultiplealpha = 1.5;
        options.WaveColor = 'r';
        options.WaveSize = 1;
        options.OrbColor = 'k';
        options.OrbSize = 1;
        options.twoD = false;
        options.delta_r = 0; 
        options.delta_r_factor = 10; 
        options.delta_surfr_factor = 1;
        options.baseMuply = 10;
        options.sign = true;
        options.view = [24,24]; 
        options.cmap = parula;
        options.efsmall = 1e-3;
    end
    
    if isempty(options.Rm)
        try
            [Rm,~,~,~,~]=POSCAR_read(options.POSCAR);
        catch
            filename = input('please give a POSCAR file :');
            [Rm,~,~,~,~]=POSCAR_read(filename);
        end
    else
        Rm = options.Rm;
    end
    if isempty(options.ax)
        Fig =  Figs(1,1);
        ax = Fig.axes(1);
    else
        ax = options.ax;
    end
    [Nwave,Nlist] = size(WaveFunc);
    Norb = length(orb_list);
    if  Norb ~= Nwave
        error('Orbital list length is not equal to WaveFunc');
    end
    %
    orb_list_2D = orb_list(:,1:3);
    orb_list_2D(:,3) = 0;
    orb_r_list = orb_list_2D*Rm;
    % Base
    Base_X = orb_r_list(:,1);
    Base_Y = orb_r_list(:,2);
    Base_Z = orb_r_list(:,3);
    WFplot_list = options.baseMuply*sum(WaveFunc.*conj(WaveFunc),2);
    if options.delta_r == 0
        delta_r = (norm(Rm(1,:))/options.lnum)/2*options.delta_r_factor;
    else
        delta_r = options.delta_r*options.delta_r_factor;
    end
    delta_r2 =delta_r^2;
    efsmall =1/options.lnum + options.efsmall;
    switch mode
        case 'square'
            % para
            minx = min(orb_list(:,1));
            miny = min(orb_list(:,2));
            maxx = max(orb_list(:,1));
            maxy = max(orb_list(:,2));
            %Center_plus = [minx miny 0.5;maxx maxy 0.5;];
            %Center_minus = [minx maxy 0.5;maxx miny 0.5];
            % Go by clockwise
            FourCorner1 = [minx miny 0.0];
            FourCorner2 = [maxx maxy 0.0];
            FourCorner3 = [minx maxy/2 0.0];
            FourCorner4 = [maxx/2 maxy 0.0];
            Center = (FourCorner2 + FourCorner1)/2*Rm;
            r = norm((FourCorner2 - FourCorner1)*Rm)/2;
            %
            [X_seed,Y_seed,Z_seed] = cylinder(r,8);
            X_seed = X_seed(:,[2,4,6,8,2]);
            Y_seed = Y_seed(:,[2,4,6,8,2]);
            Z_seed = Z_seed(:,[2,4,6,8,2]);
            Z_seed = Z_seed*options.h;
            %
            Bulk_X = [];
            Bulk_Z = [];
            Bulk_Y = [];
            for i = 1:size(X_seed,2)-1
                Bulk_X = [Bulk_X,linspace(X_seed(1,i),X_seed(2,i+1),options.lnum)];
                Bulk_Y = [Bulk_Y,linspace(Y_seed(1,i),Y_seed(2,i+1),options.lnum)];
            end
            Bulk_X = Bulk_X + Center(1);
            Bulk_Y = Bulk_Y + Center(2);
            Bulk_Z = repmat(linspace(Z_seed(1,1),Z_seed(2,1),options.hnum).',[1 size(Bulk_X,2)]);
            Bulk_Z = Bulk_Z*options.h + Center(3);
            Corlor3D_surf = zeros(size(Bulk_X));
            % information 3D
            for i =1:numel(Corlor3D_surf)
                Corlor3D_surfx = Bulk_X(i);
                Corlor3D_surfy = Bulk_Y(i);
                if sign((Corlor3D_surfx - Center(1))*(Corlor3D_surfy - Center(2)))>=0
                    Tmpsign = 1;
                elseif  options.sign
                    Tmpsign = -1;
                else
                    Tmpsign = 1;
                end
                tmpLabel = find((Base_X-Corlor3D_surfx).^2 + (Base_Y-Corlor3D_surfy).^2  <= delta_r2*options.delta_surfr_factor);
                if isempty(tmpLabel)
                else
                    WFtemp = mean(WFplot_list(tmpLabel))*options.delta_surfr_factor;
                    if WFtemp >options.WaveMin
                        Corlor3D_surf(i) = Tmpsign*mean(WFplot_list(tmpLabel));
                    end
                end
            end
            % 3D expand
            Bulk_Y = repmat(Bulk_Y,[ size(Bulk_Z,1) 1]);
            Bulk_X = repmat(Bulk_X,[ size(Bulk_Z,1) 1]);
            Corlor3D_surf = repmat(Corlor3D_surf,[ size(Bulk_Z,1) 1]);
            % 2D
            Bottom_X = linspace(FourCorner1(1),FourCorner2(1),2*options.lnum);
            Bottom_Y = linspace(FourCorner1(2),FourCorner2(2),2*options.lnum);
            [Mesh_Bottom_X,Mesh_Bottom_Y] = meshgrid(Bottom_X,Bottom_Y);
            Mesh_Bottom_X_r = Mesh_Bottom_X*Rm(1,1)+Mesh_Bottom_Y*Rm(2,1);
            Mesh_Bottom_Y_r = Mesh_Bottom_X*Rm(1,2)+Mesh_Bottom_Y*Rm(2,2);
            Mesh_Bottom_Z_r = zeros(size(Mesh_Bottom_X_r));
            Mesh_Top_Z_r = ones(size(Mesh_Bottom_X_r))*Bulk_Z(end,1);
            % Color2D
            Corlor2D_Bottom = zeros(size(Mesh_Top_Z_r));
            % information 2D
            for i =1:numel(Corlor2D_Bottom)
                if ~isnan(Corlor2D_Bottom(i))
                    Corlor2D_Bottomx = Mesh_Bottom_X_r(i);
                    Corlor2D_Bottomy = Mesh_Bottom_Y_r(i);
                    if sign((Corlor2D_Bottomx - Center(1))*(Corlor2D_Bottomy - Center(2)))>=0
                        Tmpsign = 1;
                    elseif options.sign
                        Tmpsign = -1;
                    else
                        Tmpsign = 1;
                    end 
                    tmpLabel = find((Base_X-Corlor2D_Bottomx).^2 + (Base_Y-Corlor2D_Bottomy).^2  <= delta_r2);
                    if isempty(tmpLabel)
                    else
                        WFtemp = mean(WFplot_list(tmpLabel));
                        if WFtemp >options.WaveMin
                            Corlor2D_Bottom(i) = Tmpsign * mean(WFplot_list(tmpLabel));
                        end
                    end
                end
            end
            Corlor2D_Top = Corlor2D_Bottom;
            % plot 2D
            BottomSurf = surf(ax,Mesh_Bottom_X_r,Mesh_Bottom_Y_r,Mesh_Bottom_Z_r,Corlor2D_Bottom,'EdgeColor','none');
            if options.twoD
                axis(ax,'equal');
                shading(ax,'interp');
                view(ax,0,90);
                colormap(ax,options.cmap);
                axis(ax,'off');
                return;
            end
            % plot 3D
            EdgeSurf = surf(ax,Bulk_X,Bulk_Y,Bulk_Z,Corlor3D_surf,'EdgeColor','none');
            % plot 2D
            BottomSurf = surf(ax,Mesh_Bottom_X_r,Mesh_Bottom_Y_r,Mesh_Bottom_Z_r,Corlor2D_Bottom,'EdgeColor','none');
            TopSurf = surf(ax,Mesh_Bottom_X_r,Mesh_Bottom_Y_r,Mesh_Top_Z_r,Corlor2D_Top,'EdgeColor','none');

            axis(ax,'equal');
            shading(ax,'interp');
            view(ax,options.view(1),options.view(2));
            colormap(ax,options.cmap);
            axis(ax,'off');
        case 'hex'
% para
            minx = min(orb_list(:,1));
            miny = min(orb_list(:,2));
            maxx = max(orb_list(:,1));
            maxy = max(orb_list(:,2));
            %Center_plus = [minx miny 0.5;maxx maxy 0.5;];
            %Center_minus = [minx maxy 0.5;maxx miny 0.5];
            % Go by clockwise
            SixCorner1 = [minx miny 0.0];
            SixCorner2 = [minx maxy/2 0.0];
            SixCorner3 = [maxx/2 maxy 0.0];
            SixCorner4 = [maxx maxy 0.0];
            SixCorner5 = [maxx maxy/2 0.0];
            SixCorner6 = [maxx/2 miny 0.0];
            Center = (SixCorner4 + SixCorner1)/2*Rm;
            r = norm((SixCorner4 - SixCorner1)*Rm)/2;

%             Set_2D = [SixCorner1;SixCorner2;SixCorner3;...
%                 SixCorner4;SixCorner5;SixCorner6]*Rm;
%             AS_2D = alphaShape(Set_2D,r*options.rmultiplealpha);
%             plot(AS_2D);
            %pgon1 = nsidedpoly(6,'Center',[5 0],'SideLength',3);
            % bulk
            [X_seed,Y_seed,Z_seed] = cylinder(r,6);
            Z_seed = Z_seed*options.h;
            %
            Bulk_X = [];
            Bulk_Z = [];
            Bulk_Y = [];
            for i = 1:size(X_seed,2)-1
                Bulk_X = [Bulk_X,linspace(X_seed(1,i),X_seed(2,i+1),options.lnum)];
                Bulk_Y = [Bulk_Y,linspace(Y_seed(1,i),Y_seed(2,i+1),options.lnum)];
            end
            Bulk_X = Bulk_X + Center(1);
            Bulk_Y = Bulk_Y + Center(2);
            Bulk_Z = repmat(linspace(Z_seed(1,1),Z_seed(2,1),options.hnum).',[1 size(Bulk_X,2)]);
            Bulk_Z = Bulk_Z*options.h + Center(3);
            Corlor3D_surf = zeros(size(Bulk_X));
            % information 3D
            for i =1:numel(Corlor3D_surf)
                Corlor3D_surfx = Bulk_X(i);
                Corlor3D_surfy = Bulk_Y(i);
                tmpLabel = find((Base_X-Corlor3D_surfx).^2 + (Base_Y-Corlor3D_surfy).^2  <= delta_r2*options.delta_surfr_factor);
                if isempty(tmpLabel)
                else
                    WFtemp = mean(WFplot_list(tmpLabel))*options.delta_surfr_factor;
                    if WFtemp >options.WaveMin
                        Corlor3D_surf(i) = mean(WFplot_list(tmpLabel));
                    end
                end
            end
            % 3D expand
            Bulk_Y = repmat(Bulk_Y,[ size(Bulk_Z,1) 1]);
            Bulk_X = repmat(Bulk_X,[ size(Bulk_Z,1) 1]);
            Corlor3D_surf = repmat(Corlor3D_surf,[ size(Bulk_Z,1) 1]);
            % 2D
            Bottom_X = linspace(SixCorner1(1),SixCorner4(1),2*options.lnum);
            Bottom_Y = linspace(SixCorner1(2),SixCorner4(2),2*options.lnum);
            [Mesh_Bottom_X,Mesh_Bottom_Y] = meshgrid(Bottom_X,Bottom_Y);
            Mesh_Bottom_X_r = Mesh_Bottom_X*Rm(1,1)+Mesh_Bottom_Y*Rm(2,1);
            Mesh_Bottom_Y_r = Mesh_Bottom_X*Rm(1,2)+Mesh_Bottom_Y*Rm(2,2);
            Mesh_Bottom_Z_r = zeros(size(Mesh_Bottom_X_r));
            Mesh_Top_Z_r = ones(size(Mesh_Bottom_X_r))*Bulk_Z(end,1);
            % Color2D
            Corlor2D_Bottom = zeros(size(Mesh_Top_Z_r));
            % Two-point Form (y-y2)/ (y1-y2) = (x-x2)/ (x1-x2)
            % y- (y1-y2) / (x1-x2) *x > y2-x2* (y1-y2) / (x1-x2)
            a = -(SixCorner2(2)-SixCorner3(2))/(SixCorner2(1)-SixCorner3(1));
            b = a*SixCorner3(1)+SixCorner3(2);
            exceed_label1 = a*Mesh_Bottom_X + Mesh_Bottom_Y> b+efsmall;
            %  -y2*(x1-x2)/ (y1-y2) + x2 < x+ -(y)*(x1-x2)/ (y1-y2);
            a2 = -(SixCorner6(1)-SixCorner5(1))/(SixCorner6(2)-SixCorner5(2));
            b2 = a*SixCorner5(2)+SixCorner5(1);
            exceed_label2 = Mesh_Bottom_X +a2* Mesh_Bottom_Y> b2+efsmall;
            Corlor2D_Bottom(exceed_label1) = nan;
            Corlor2D_Bottom(exceed_label2) = nan;
            % information 2D
            for i =1:numel(Corlor2D_Bottom)
                if ~isnan(Corlor2D_Bottom(i))
                    Corlor2D_Bottomx = Mesh_Bottom_X_r(i);
                    Corlor2D_Bottomy = Mesh_Bottom_Y_r(i);
                    tmpLabel = find((Base_X-Corlor2D_Bottomx).^2 + (Base_Y-Corlor2D_Bottomy).^2  <= delta_r2);
                    if isempty(tmpLabel)
                    else
                        WFtemp = mean(WFplot_list(tmpLabel));
                        if WFtemp >options.WaveMin
                        Corlor2D_Bottom(i) = mean(WFplot_list(tmpLabel));
                        end
                    end
                end
            end
            Corlor2D_Top = Corlor2D_Bottom;
            % plot 2D
            BottomSurf = surf(ax,Mesh_Bottom_X_r,Mesh_Bottom_Y_r,Mesh_Bottom_Z_r,Corlor2D_Bottom,'EdgeColor','none');
            if options.twoD
                axis(ax,'equal');
                shading(ax,'interp');
                view(ax,0,90);
                colormap(ax,options.cmap);
                axis(ax,'off');
                return;
            end
            % plot 3D
            EdgeSurf = surf(ax,Bulk_X,Bulk_Y,Bulk_Z,Corlor3D_surf,'EdgeColor','none');
            % plot 2D
            BottomSurf = surf(ax,Mesh_Bottom_X_r,Mesh_Bottom_Y_r,Mesh_Bottom_Z_r,Corlor2D_Bottom,'EdgeColor','none');
            TopSurf = surf(ax,Mesh_Bottom_X_r,Mesh_Bottom_Y_r,Mesh_Top_Z_r,Corlor2D_Top,'EdgeColor','none');

            axis(ax,'equal');
            shading(ax,'interp');
            view(ax,options.view(1),options.view(2));
            colormap(ax,options.cmap);
            axis(ax,'off');
    end
end