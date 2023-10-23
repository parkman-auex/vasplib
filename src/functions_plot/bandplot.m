function varargout=bandplot(EIGENCAR,Ecut,klist_l,kpoints_l,kpoints_name,options)
arguments
    EIGENCAR = EIGENVAL_read();
    Ecut double= [-3,3];
    klist_l double= [];
    kpoints_l double= [];
    kpoints_name string= [];
    options.ax =  handle([]);
    options.FontName = 'Helvetica';
    options.FontSize = 24;
    options.KPOINTS = 'KPOINTS';
    options.POSCAR = 'POSCAR';
    options.Color =  @jet;
    options.title = '';
    options.klist_l = [];
    options.kpoints_l = [];
    options.kpoints_name = [];
    options.xlabel='';
    options.ylabel='E (eV)';
    options.LineSpec = '-';
    options.LineWidth = 1;
    options.MarkerSize = 3;
    options.MarkerEdgeColor = [];
    options.MarkerFaceColor = [];
    options.set_reference = true;
    options.Units = 'pixels';
    options.Position = [];
end

%--------  init  --------
    import vasplib_tool.*
    nOutputs = nargout;
    varargout = cell(1,nOutputs);
%--------  narg  -------- 
    if isempty(klist_l) 
       Rm=POSCAR_readin(options.POSCAR);
       [~,klist_l,~,kpoints_l,kpoints_name]=kpathgen3D(Rm,options.KPOINTS);
    end
    if iscell(EIGENCAR)
        Nbands=size(EIGENCAR{1},1);
        color = options.Color ;
    else
        Nbands=size(EIGENCAR,1);
        if ischar(options.Color) || isnumeric(options.Color)
            color = options.Color;
        else
            color = [rand,rand,rand];
        end
        
    end
    if isempty(options.MarkerEdgeColor)
        MarkerEdgeColor = color;
    else
        MarkerEdgeColor = options.MarkerEdgeColor;
    end
    if isempty(options.MarkerFaceColor)
        MarkerFaceColor = color;
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
        [ax] = vasplib_plot.set_reference(kpoints_l,kpoints_name,Xcut,Ecut,...
            'ax',ax,...
            'fontname',options.FontName,...
            'FontSize',options.FontSize,...
            'xlabel',options.xlabel,...
            'ylabel',options.ylabel ...
            );
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
    if iscell(EIGENCAR)
        Npara = length(EIGENCAR);
        if isa(color,'function_handle')
            Colormap = color(Npara);
        else
            Colormap = color;
        end
        for j = 1:Npara
            switch j
                case {1,Npara,floor((Npara+1)/2),ceil((Npara+1)/4),floor((Npara+1)*3/4)}
                    LineWidth = 3.0;
                otherwise
                    LineWidth = 1;
            end
            Nbands=size(EIGENCAR{j},1);
            for Ei=1:Nbands
                
                plot(ax,klist,EIGENCAR{j}(Ei,:),options.LineSpec,'LineWidth',LineWidth,'Color',Colormap(j,:),'DisplayName',num2str(j)+"_"+num2str(Ei));
                hold(ax,'on');
                %hold on;
            end
        end
    else
        if isa(color,'function_handle')
            color = color(1);
        end
        for Ei=1:Nbands
            plot(ax,klist,EIGENCAR(Ei,:),options.LineSpec,...
                'LineWidth',options.LineWidth,...
                'Color',color,...
                'MarkerSize',options.MarkerSize,...
                'MarkerEdgeColor',MarkerEdgeColor,...
                'MarkerFaceColor',MarkerFaceColor,...
                'DisplayName',num2str(Ei));
            %hold(ax,'on');
            %hold on;
        end
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
        plot_data_coll.color = color;
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

    

