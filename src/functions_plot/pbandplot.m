%% pbandplot
%
% plotpband
%
%
% * Label: plot
%
%% Description of the Function:
%%
%% Usage:
%
% * [fig,ax]=pbandplot(WEIGHTCAR_struct,EIGENCAR,Ecut,titlestring,cmap,klist_l,kpoints_l,kpoints_name,fontname,fig,ax)
% * [fig,ax]=pbandplot(WEIGHTCAR_struct,EIGENCAR,Ecut,titlestring,cmap,klist_l,kpoints_l,kpoints_name,fontname)
% * [fig,ax]=pbandplot(WEIGHTCAR_struct,EIGENCAR,Ecut,titlestring,cmap,klist_l,kpoints_l,kpoints_name)
% * [fig,ax]=pbandplot(WEIGHTCAR_struct,EIGENCAR,Ecut,titlestring,cmap)
% * [fig,ax]=pbandplot(WEIGHTCAR_struct,EIGENCAR,Ecut,titlestring)
% * [fig,ax]=pbandplot(WEIGHTCAR_struct,EIGENCAR,Ecut)
% * [fig,ax]=pbandplot(WEIGHTCAR_struct,EIGENCAR)
% * [fig,ax]=pbandplot(WEIGHTCAR_struct,EIGENCAR)
% * [fig,ax]=pbandplot(WEIGHTCAR_struct,EIGENCAR)
% * [fig,ax]=pbandplot(WEIGHTCAR_struct)
% * [fig,ax]=pbandplot()
%% Input:
%
% # input1:
% # input2:
% # input3:
%
%% Output:
%
% # fig:
% # ax:
%
%% Note:
%
%  Take advantage of the scope of application of the function.
%
%% Change log
%
% * Document Date: 2020/12/04
% * Creation Date: 2020/12/04
% * Last updated : 2023/09/04
%
%% Source code :
%
function varargout = pbandplot(WEIGHTCAR_struct,EIGENCAR,klist_l,kpoints_l,kpoints_name,options,optionsplot)
arguments
    WEIGHTCAR_struct = [];
    EIGENCAR = EIGENVAL_read();
    klist_l double= [];
    kpoints_l double= [];
    kpoints_name string= [];
    options.ColorCut double= 1;
    options.ColorCutMinus double= -1;
    options.Ecut = [-3,3];
    options.fontname = 'Helvetica';
    options.KPOINTS = 'KPOINTS';
    options.POSCAR = 'POSCAR';
    options.cmap =  @jet;
    options.title = '';
    options.klist_l = [];
    options.kpoints_l = [];
    options.kpoints_name = [];
    options.xlabel='';
    options.ylabel='E(eV)';
    options.LineSpec = '-';
    options.LineWidth = 1;
    options.MarkerEdgeColor = [];
    options.MarkerFaceColor = [];
    options.silent = false;   % If you do not want to be asked to choose the projections, turn on it.
    optionsplot.density = 1;
    optionsplot.WEIGHTCAR_factor = 1;
    optionsplot.ax = handle([]);
    optionsplot.filled = true;   
end
%--------  init  --------
import vasplib_tool.*
%--------  narg  --------
if isempty(klist_l)
    Rm=POSCAR_read(options.POSCAR);
    [~,klist_l,~,kpoints_l,kpoints_name]=kpathgen3D(Rm,options.KPOINTS);
end
if strcmp(options.title,'')
    [~,titlestring]=fileparts(pwd);
    titlestring = strrep(titlestring,'_','-');
else
    titlestring  = options.title;
end
if options.ColorCutMinus == -1
    options.ColorCutMinus = options.ColorCut;
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
%--------  juge 2   --------
if isempty(WEIGHTCAR_struct)
    % WEIGHTCAR_in_a_file;
    [WEIGHTCAR_struct_cell,Name_list,~] = WEIGHTCAR_gen_dir('PBAND');
    [~,nWEIGHTCAR_struct_cell] = size(WEIGHTCAR_struct_cell);
    WEIGHTCAR_TOT = struct();
    % test 2021a new feature
    % propertyCell = namedargs2cell(options);
    %
    for i = 1:nWEIGHTCAR_struct_cell        
        element_name = split(Name_list(i,:),'_');
        titlestring_new = titlestring+"-"+element_name(1);
        if ~options.silent
            disp("Please choose which projections of "+element_name(1)+" will be plotted:");
            disp("(array format like 1:4, [1,4,5]... is supported)");
        end

        options.title = titlestring_new;
        propertyCell = namedargs2cell(options);
        WEIGHTCAR_struct = WEIGHTCAR_struct_cell{i};
        % rm the last
        WEIGHTCAR_struct(end) = [];
        %             nWEIGHTCAR = length(WEIGHTCAR_struct);
        %             cmap = cmap_function(nWEIGHTCAR);
        ax_temp=vasplib.pbandplot(WEIGHTCAR_struct,EIGENCAR,klist_l,kpoints_l,kpoints_name,propertyCell{:});
        %--------  fbug  --------
        fig(i) = ax_temp.Parent;
        ax(i) = ax_temp;
        % if i ==1
        %     WEIGHTCAR_TOT = WEIGHTCAR_struct;
        % else
        %     WEIGHTCAR_TOT = [WEIGHTCAR_TOT ,WEIGHTCAR_struct];
        % end
    end
    %         nWEIGHTCAR = length(WEIGHTCAR_TOT);
    %         cmap = cmap_function(nWEIGHTCAR);
    %--------  fbug  --------
    % options.title = 'TOT';
    % propertyCell = namedargs2cell(options);

    % fig(i+1) = ax_temp.Parent;ax(i+1) = ax_temp;
    %-------- return --------
    if nargout  == 2
        varargout{1} = fig;
        varargout{2} = ax;
    end
    if nargout  == 1
        varargout{1} = ax;
    end
    return;
end
%--------  norm  --------
if isa(WEIGHTCAR_struct,'double')
    %warning('This method has not been revised!')
    nWEIGHTCAR = length(WEIGHTCAR_struct);
    pbandmode = 'patch'         ;
    WEIGHTCAR = WEIGHTCAR_struct;
    maxBCplus = options.ColorCut*max((WEIGHTCAR),[],'all');
    maxBCMinus = options.ColorCutMinus*min((WEIGHTCAR),[],'all');
    maxBC = max(abs(maxBCplus),abs(maxBCMinus));
    WEIGHTCAR(WEIGHTCAR>maxBCplus) = maxBCplus;
    WEIGHTCAR(WEIGHTCAR<maxBCMinus) = maxBCMinus;
    if isa(options.cmap,'function_handle')
        try
            cmap = options.cmap(64,[-maxBC,maxBC]);
        catch
            cmap = options.cmap(64);
        end
    else
        cmap = options.cmap;
    end

    %cmap = cmap_function();
elseif isa(WEIGHTCAR_struct,'cell')
    warning('This method has not been revised!')

    cmap = options.cmap;
    if length(WEIGHTCAR_struct) > 1
        pbandmode = 'bubble_rough' ;
        WEIGHTCAR_cell = WEIGHTCAR_struct;
        nWEIGHTCAR = length(WEIGHTCAR_cell);
        Name_list(nWEIGHTCAR ,:) = "";
        cmap = options.cmap(nWEIGHTCAR);
    else
        pbandmode = 'bubble_only' ;
        WEIGHTCAR_mat = WEIGHTCAR_struct{1};
        nWEIGHTCAR = length(WEIGHTCAR_mat);
        Name_list(nWEIGHTCAR ,:) = "";
    end
    for i = 1:nWEIGHTCAR
        Name_list(i,:) = strcat('A',string(i));
    end
elseif isa(WEIGHTCAR_struct,'struct')
    % it means we can rebuild code for highly user-like plot
    pbandmode = 'bubble_refine';
    nWEIGHTCAR = length(WEIGHTCAR_struct);
    Name_list(nWEIGHTCAR ,:) = "";
    WEIGHTCAR_cell{nWEIGHTCAR} = WEIGHTCAR_struct(nWEIGHTCAR).WEIGHTCAR;
    for i = 1 : nWEIGHTCAR
        WEIGHTCAR_cell{i} =  WEIGHTCAR_struct(i).WEIGHTCAR;
        Name_list(i ,:) = WEIGHTCAR_struct(i).displayname;
        % if color
        % if alpha
    end
    %--------  bug  --------
    if isa(options.cmap,'function_handle')
        cmap = options.cmap(nWEIGHTCAR);
    else
        cmap = options.cmap;
    end
    nWEIGHTCAR = length(WEIGHTCAR_cell);
end
%--------  init  --------
Nbands=size(EIGENCAR,1);
xmin=klist_l(1);
xmax=klist_l(length(klist_l));
Xcut = [xmin,xmax];
%--------projections----                    

if nWEIGHTCAR ~= 9
   % warning("It seems like that the PBAND info does not come from vaspkit!")
end
if ~options.silent && (strcmp(pbandmode,"bubble_only") || strcmp(pbandmode,"bubble_rough") ||strcmp(pbandmode,"bubble_refine"))
    prompt = "0) all(or just press Enter)\n1) s    2) py    3) pz    4) px\n"...
        +"5) dxy  6) dyz   7) dz2   8) dxz   9) dx2-y2:\n";
    Selected_projs = input(prompt);
    if isempty(Selected_projs)
        Selected_projs = 1:nWEIGHTCAR;
    end
else
    Selected_projs = 1:nWEIGHTCAR;
end
%--------  init ax ------
if isempty(optionsplot.ax)
    Fig =  Figs(1,1);
    ax = Fig.axes(1);
else
    ax = optionsplot.ax;
end
%--------plotweight--------
switch pbandmode
    case 'patch'
        %             for Ei=1:Nbands
        %                 plot(ax,klist,EIGENCAR(Ei,:),'LineWidth',1.0,'Color',[0.1 0.1 0.1],'DisplayName',num2str(Ei));
        %                 %hold on;
        %             end
        ax = pband_plot_one_patch(klist,EIGENCAR,WEIGHTCAR,'cmap',cmap,'ax',ax,'LineWidth',options.LineWidth);
    case 'bubble_only'
        optionsplot.ax =ax;
        optionsplotcell = namedargs2cell(optionsplot);
        ax = pband_plot_one_bubble(klist,EIGENCAR,WEIGHTCAR_mat,cmap,Name_list,optionsplotcell{:});
    case 'bubble_rough'
        for Ei=1:Nbands
            plot(ax,klist,EIGENCAR(Ei,:),'LineWidth',1.0,'Color',[0.1 0.1 0.1],'DisplayName',num2str(Ei));
            %hold on;
        end
        ax = pband_plot_set(klist,EIGENCAR,WEIGHTCAR_cell,Name_list,cmap,options.density,ax);
    case 'bubble_refine'
        %--------baseline------
        for Ei=1:Nbands
            plot(ax,klist,EIGENCAR(Ei,:),...
                'LineWidth',1.0,'Color',[0.1 0.1 0.1],'DisplayName',num2str(Ei),...
                'HandleVisibility','off');
            %hold on;
        end
        optionsplot.ax =ax;
        optionsplotcell = namedargs2cell(optionsplot);
        ax = pband_plot_set(klist,EIGENCAR,WEIGHTCAR_cell,Name_list,Selected_projs,...
            'cmap',cmap,optionsplotcell{:});% waiting
        %---------legend-------
        %legend_all = ["s", "p_y", "p_z", "p_x", "d_{xy}", "d_{yz}", "d_{z^2}", "d_{xz}", "d_{x^2-y^2}"];
        %legend(ax,legend_all(Selected_projs))
end
%--------reference -------
ax = set_reference(kpoints_l,kpoints_name,Xcut,options.Ecut,...
    'ax',ax,...
    'fontname',options.fontname ...
    );
%[fig,ax] = set_reference(fig,ax,'band',Xcut,Ecut,kpoints_l,kpoints_name,fontname);
%--------  fbug  --------
if iscell(titlestring )
    titlestring = cell2mat(titlestring);
end
%--------  title  -------
title(ax,titlestring);
if nargout>2
    %--------  save  --------
    plot_data_coll.EIGENCAR = EIGENCAR;
    plot_data_coll.Ecut = Ecut;
    plot_data_coll.titlestring = titlestring;
    plot_data_coll.color = color;
    plot_data_coll.klist_l = klist_l;
    plot_data_coll.kpoints_l = kpoints_l;
    plot_data_coll.kpoints_name = kpoints_name;
    plot_data_coll.fontname = options.fontname;
    plot_data_coll.WEIGHTCAR_struct = WEIGHTCAR_struct;
    varargout{3} = plot_data_coll;
    varargout{1} = ax.Parent;
    varargout{2} = ax;
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