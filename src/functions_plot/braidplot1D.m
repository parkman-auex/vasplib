function varargout=braidplot1D(EIGENCAR,Ecut,klist_l,kpoints_l,kpoints_name,options)
arguments
    EIGENCAR = EIGENVAL_read();
    Ecut double= [-3,3;-3,3];
    klist_l double= [];
    kpoints_l double= [];
    kpoints_name string= [];
    options.ax =  handle([]);
    options.FontName = 'Helvetica';
    options.FontSize = 24;
    options.KPOINTS = 'KPOINTS';
    options.POSCAR = 'POSCAR';
    options.Color =  @parula;
    options.title = '';
    options.klist_l = [];
    options.kpoints_l = [];
    options.kpoints_name = [];
    options.xlabel='';
    options.ylabel='E (eV)';
    options.LineSpec = '-';
    options.LineWidth = 3;
    options.MarkerSize = 3;
    options.MarkerEdgeColor = 'none';
    options.MarkerFaceColor = 'none';
    options.set_reference = true;
    options.set_plane = true;
    options.Units = 'pixels';
    options.Position = [];
    options.flip = true;
    options.Mirror = -1;
end

%--------  init  --------
    import vasplib_tool.*
    nOutputs = nargout;
    varargout = cell(1,nOutputs);
%--------  narg  -------- 
    if isempty(klist_l) 
       Rm=POSCAR_read(options.POSCAR);
       [~,klist_l,~,kpoints_l,kpoints_name]=kpathgen3D(Rm,options.KPOINTS);
    end
    Nbands=size(EIGENCAR,1);
    if isstring(options.Color)
        if isrow
            colorL = options.Color.';
        end
    elseif isnumeric(options.Color)
        x = linspace(0,1,size(options.Color,1));
        xq = linspace(0,1,Nbands);
        colorL = [...
            interp1(x,options.Color(:,1),xq).',...
            interp1(x,options.Color(:,2),xq).',...
            interp1(x,options.Color(:,3),xq).',...
            ];
    elseif isa(options.Color,'function_handle')
        colorL = options.Color(Nbands);
    else
        for i = 1:Nbands
            colorL(i,:) = [rand,rand,rand];
        end
    end
    if options.flip
        colorL = flip(colorL,1);
    end
    if isempty(options.MarkerEdgeColor)
        MarkerEdgeColor = colorL;
    else
        MarkerEdgeColor = options.MarkerEdgeColor;
    end
    if isempty(options.MarkerFaceColor)
        MarkerFaceColor = colorL;
    else
        MarkerFaceColor = options.MarkerFaceColor;
    end
    if strcmp(options.title,'')
        [~,titlestring]=fileparts(pwd);
        titlestring = strrep(titlestring,'_','-');
    else
        titlestring  = options.title;
    end
    if isempty(options.ax)
        Fig =  Figs(1,1);
        ax = Fig.axes(1);
    else
        if isvalid(options.ax)
            ax = options.ax;
        else
            Fig =  Figs(1,1);
            ax = Fig.axes(1);
        end
    end
%--------  juge  --------
    if isempty(options.klist_l)
        klist=klist_l;
    else
        klist=options.klist_l;
    end
    if isempty(options.kpoints_l)
%         kpoints_l=kpoints_l;
    else
        kpoints_l=options.kpoints_l;
    end
    if isempty(options.kpoints_name)
%         kpoints_name=kpoints_name;
    else
        kpoints_name=options.kpoints_name;
    end
    xmin=klist(1);
    xmax=klist(end);
    Xcut = [xmin,xmax];
    if options.set_reference
        %--------reference -------
        %ax.BoxStyle = 'full';
        %grid(ax,'on');
        Ycut = Ecut(1,:);
        Zcut = Ecut(2,:);
        %view(ax,-16,24);
        %set(ax,'Projection','perspective');
        axis(ax,'square');
        ylabel(ax,'Re(E(z))');
        zlabel(ax,'Im(E(z))');
        xlabel(ax,'kpath');
        % z plane
        xticks(ax,kpoints_l);
        xticklabels(ax,kpoints_name);
        xlim(ax,Xcut);
        ylim(ax,Ycut);
        zlim(ax,Zcut);
        
        if options.set_plane
            %count = 0;
            nkl = length(kpoints_l);
            for i = 1:length(kpoints_l)
                V{i} = [
                    kpoints_l(i),Ycut(1),Zcut(1);
                    kpoints_l(i),Ycut(1),Zcut(2);
                    kpoints_l(i),Ycut(2),Zcut(2);
                    kpoints_l(i),Ycut(2),Zcut(1);
                    ];
                patch(ax,'Faces',[1,2,3,4],'Vertices',V{i},...
                    'EdgeColor','none',...
                    'FaceColor',[102 102 102]./255,'FaceAlpha',0.1,...
                    'DisplayName',['T_',num2str(i)]);
            end  

        end
    end
    
    %--------color -------
    %    version = ver('matlab');
    %    version_num = version.Version;
    % %    if version_num > 9.6
    %         colororder(ax,color);
    %         %--------plotwhole--------
    %         for Ei=1:Nbands
    %             plot(ax,klist,EIGENCAR(Ei,:),'LineWidth',1.0,'DisplayName',num2str(Ei));
    %             hold(ax,'on');
    %             %hold on;
    %         end
    %    else
    RealEIGENCAR = real(EIGENCAR);
    ImagEIGENCAR = imag(EIGENCAR);
    nKinEIGENCAR = size(EIGENCAR,2);
    for Ei=1:Nbands
        plot3(ax,klist(1:nKinEIGENCAR),RealEIGENCAR(Ei,:),options.Mirror*ImagEIGENCAR(Ei,:),options.LineSpec,...
            'LineWidth',options.LineWidth,...
            'Color',colorL(Ei,:),...
            'MarkerSize',options.MarkerSize,...
            'MarkerEdgeColor',MarkerEdgeColor,...
            'MarkerFaceColor',MarkerFaceColor,...
            'DisplayName',['E_',num2str(Ei)]);
        %hold(ax,'on');
        %hold on;
    end
    if options.set_plane
        HSV = rgb2hsv(colorL);
        HSV(:,2) = HSV(:,2);
        HSV(:,3) = HSV(:,3)-0.2;
        ModifycolorL = hsv2rgb(HSV);
        for i = 1:length(kpoints_l)
            SelectL = find(kpoints_l(i) == klist);
            if ~isempty(SelectL)
                try
                    scatter3(ax,ones(1,Nbands)*kpoints_l(i),RealEIGENCAR(:,SelectL(1)),ImagEIGENCAR(:,SelectL(1)),...
                        ones(1,Nbands)*options.LineWidth*20,ModifycolorL,"filled",'DisplayName',['T_',num2str(i)]);
                catch
                end
            end
            %L = legend(ax);
        %L.String([1:nkl,nkl+Nbands-1:2*nkl-2+Nbands]) = '';
    end
    %--------  fbug  --------
    if isa(titlestring,'cell')
        titlestring= cell2mat(titlestring);
    end
    %--------  title  -------
    title(ax,titlestring);
%    end

    if length(varargout)>2
    %--------  save  --------
        plot_data_coll.EIGENCAR = EIGENCAR;
        plot_data_coll.Ecut = Ecut;
        plot_data_coll.titlestring = titlestring;
        plot_data_coll.color = colorL;
        plot_data_coll.klist_l = klist_l;
        plot_data_coll.kpoints_l = kpoints_l;
        plot_data_coll.kpoints_name = kpoints_name;
        plot_data_coll.fontname = options.fontname;
        %disp(class(titlestring));
    %     mkdir('plot_results');
        [~,ax] = save_figure(ax.Parent.Parent,titlestring+".eps",ax);
        varargout{3} = plot_data_coll;
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

    

