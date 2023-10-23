classdef MaterialAccom
    properties
        ID ; % the unique ID
        SG ;
        Source = 'unkown';
        checkindate    ;
        isInsulator    ;
        isRealChern    ;
        SYM_operation;
        Rm ;
        position ;
    end
    properties
        Atom_num_total ;
        isInversion    ;
        isC2z          ;
        isC2x          ;
        isC2y          ;
        SYM_num        ;
        neckname_3D    ;
        neckname_2D    ;
        neckname       ;
    end
    properties
        BAND_gap1 ;
        BAND_gap2 ;
        HOMO      ;
        TRIM      ;
        RealChern ;
        EIGENCAR  ;
        klist_l   ;
        kpoints_l ;
        kpoints_name ;
        symmetry_dataset;
    end
    
    properties(Hidden = true )
        Atom_num ; % list
        Atom_name; % list
        infinit_small = 1e-3;
        tolerence     = 2e-2;
        sites ;
    end
    properties(Access = private )
        POSCAR_file;
        POSCAR_syminfor_file;

    end
    methods
        function MA = MaterialAccom(POSCAR_file,POSCAR_syminfor_file,Source)
            if nargin < 3
                Source = 'unkown';
            end
            if nargin < 2
                POSCAR_syminfor_file = "SymInfor" + POSCAR_file;
            end
            MA.ID          = str2double(strrep(POSCAR_file,'POSCAR.',''));
            MA.POSCAR_file = POSCAR_file;
            MA.POSCAR_syminfor_file = POSCAR_syminfor_file;
            MA.Source = Source ;
            MA.checkindate = date;
            [MA.SYM_operation,MA.SG] = MaterialAccom.sym_mat_gen(MA.POSCAR_syminfor_file);
            % read POSCAR
            [MA.Rm,MA.sites,MA.Atom_name,MA.Atom_num,~]=POSCAR_readin(MA.POSCAR_file,'vasp');
            MA.position = MA.position_gen(MA.sites);
        end
        
    end
    methods % Dynamic
        function neckname_3D = get.neckname_3D(MA)
            neckname_3D = MA.SG+"-"+MA.SYM_num+"-"+MA.Atom_num_total+"-"+MA.ID;
        end
        function neckname_2D = get.neckname_2D(MA)
            neckname_2D = MA.SG+"-"+MA.SYM_num+"-"+MA.Atom_num_total+"-"+MA.isInversion+"-"+MA.isC2z+"-"+MA.ID;
        end
        function neckname = get.neckname(MA)
            neckname = "";
            for i = 1:length(MA.Atom_name)
                neckname = neckname+MA.Atom_name(i)+num2str(MA.Atom_num(i));
            end
        end
        function Atom_num_total = get.Atom_num_total(MA)
            Atom_num_total = sum(MA.Atom_num);
        end
        function isRealChern  = get.isRealChern(MA)
            if isempty(MA.RealChern)
                isRealChern = 0;
            elseif sum(MA.RealChern) > 0
                isRealChern = 1;
            else
                isRealChern = 0;
            end
        end
        function isInsulator = get.isInsulator(MA)
            if isempty(MA.BAND_gap1) 
                isInsulator = nan;
                return
            end
            if MA.BAND_gap1 <0
                isInsulator = 0;
                return
            end
            if isempty(MA.BAND_gap2)
                isInsulator = 1;
                return
            end
            if MA.BAND_gap2 >0
                isInsulator = 1;
                return
            elseif MA.BAND_gap2 <= 0 
                isInsulator = 2;
                return
            end
        end
        function isInversion = get.isInversion(MA)
            %
            target_rotation = [...
                -1,0,0;...
                0,-1,0;...
                0,0,-1];
            target_translation =[0,0,0];
            %
            isInversion = MA.isSymOper(target_rotation,target_translation);
        end
        function isC2z = get.isC2z(MA)
            %
            target_rotation = [...
                -1,0,0;...
                0,-1,0;...
                0,0,1];
            target_translation =[0,0,0];
            %
            isC2z = MA.isSymOper(target_rotation,target_translation);
        end
        function isC2x = get.isC2x(MA)
            %
            target_rotation = [...
                1,0,0;...
                0,-1,0;...
                0,0,-1];
            target_translation =[0,0,0];
            %
            isC2x = MA.isSymOper(target_rotation,target_translation);
        end
        function isC2y = get.isC2y(MA)
            %
            target_rotation = [...
                -1,0,0;...
                0,1,0;...
                0,0,-1];
            target_translation =[0,0,0];
            %
            isC2y = MA.isSymOper(target_rotation,target_translation);
        end
        function SYM_num = get.SYM_num(MA)
            SYM_num = size(MA.SYM_operation,3);
        end
        function isSymOper_label = isSymOper(MA,target_rotation,target_translation)
            isSymOper_label = 0;
            sym_mat = MA.SYM_operation;
            for i =1:size( sym_mat,3)
                if all(sym_mat(1:3,:,i) == target_rotation)
                    if norm(sym_mat(4,:,i)-target_translation)<MA.infinit_small
                        isSymOper_label = 1;
                        return;
                    end
                    
                end
            end
            
        end
    end
    
    methods % overload
        function label = eq(A,B)
            if isa(A,'MaterialAccom') && isa(B,'MaterialAccom')
                label = 1;
                if A.SG ~= B.SG
                    label = 0;
                    return
                end
                if A.Atom_num_total ~= B.Atom_num_total
                    label = 0;
                    return
                end
                tolerence_tmp = (A.tolerence+B.tolerence)/2;
                for i = 1:3
                    for j = 1:3
                        if abs(A.Rm(i,j) - B.Rm(i,j)) > tolerence_tmp
                            label = 0;
                            return
                        end
                    end
                end
                for i =1:A.Atom_num_total
                    if norm(A.position-B.position) > tolerence_tmp
                        label = 0;
                        return
                    end
                end
                if ~isequal(A.Atom_name,B.Atom_name)
                    label = 0;
                end
                if ~isequal(A.Atom_num,B.Atom_num)
                    label = 0;
                end
            else
                error('! not MaterialAccom')
            end
        end
    end
    methods % tools
        function neckname = neckname_ppt(MA)
            import mlreportgen.ppt.*;
            neckname = Paragraph("");
            for i = 1:length(MA.Atom_name)
                subs_num = Text(num2str(MA.Atom_num(i)));
                subs_num.Subscript = true;
                neckname.append(MA.Atom_name(i));
                neckname.append(subs_num);
            end
        end
        function MA = POSCAR_gen(MA,filename)
            if nargin < 2
                filename = "POSCAR_"+MA.neckname_3D;
            end
            POSCAR_gen( MA.Rm,MA.sites,MA.Atom_name,MA.Atom_num,filename);
        end
        function [fig,ax] = POSCAR_plot(MA,filename_w_suffix,options)
            arguments
                MA
                filename_w_suffix
                options.fig = figure();
                options.ax = gca();
            end
            if nargin >1
                if length(MA.Atom_name) == 1
                    MA.Atom_name = [MA.Atom_name,"H"];
                    MA.Atom_num = [MA.Atom_num,0];
                end
                [fig,ax]  = POSCAR_plot(MA.Rm,MA.sites,[MA.Atom_name],MA.Atom_num,...
                    'atomscale',0.3,...
                    'title' ,'',...
                    'fig',options.fig,...
                    'ax',options.ax,...
                    'box',[-0.5,-0.5,-0.2;1.5,1.5,1.2],...
                    'view',[0,90],...
                    'fast',false);...
                    grid(ax,'off');
                axis(ax,'off');
                %                 exportgraphics(ax,filename_w_suffix(1),'BackgroundColor','none');
                exportgraphics(ax,filename_w_suffix(1));
                cla(ax);
                [fig,ax]  = POSCAR_plot(MA.Rm,MA.sites,[MA.Atom_name],MA.Atom_num,...
                    'atomscale',0.3,...
                    'title' ,'',...
                    'fig',fig,...
                    'ax',options.ax,...
                    'box',[-0.5,-0.5,-0.2;1.5,1.5,1.2],...
                    'view',[90,0],...
                    'fast',false);...
                    grid(ax,'off');
                axis(ax,'off');
                %                 exportgraphics(ax,filename_w_suffix(2),'BackgroundColor','none');
                exportgraphics(ax,filename_w_suffix(2));
                cla(ax);
                axis(ax,'on');
                return;
            end
            [fig,ax]  = POSCAR_plot(MA.Rm,MA.sites,MA.Atom_name,MA.Atom_num,...
                'atomscale',0.3,...
                'title' ,'',...
                'fig',options.fig,...
                'ax',options.ax,...
                'box',[-0.5,-0.5,-0.2;1.5,1.5,1.2],...
                'view',[30,90]);...
                grid(ax,'off');
                axis(ax,'off');
                
        end
        function R_struct = Rm2abc(MA)
            R_struct.a = norm(MA.Rm(1,:));
            R_struct.b = norm(MA.Rm(2,:));
            R_struct.c = norm(MA.Rm(3,:));
            R_struct.alpha = acos(dot(MA.Rm(2,:),MA.Rm(3,:))/(norm(R_struct.b)*norm(R_struct.c)))/pi*180;
            R_struct.beta  = acos(dot(MA.Rm(3,:),MA.Rm(1,:))/(norm(R_struct.c)*norm(R_struct.a)))/pi*180;
            R_struct.gamma = acos(dot(MA.Rm(1,:),MA.Rm(2,:))/(norm(R_struct.a)*norm(R_struct.b)))/pi*180;
        end
        function [fig,ax] = bandplot(MA,options)
            arguments
                MA;
                options.Ecut = [-6,6];
                options.title = ' ';
                options.color = [rand rand rand];
                options.savefile  = '';
                options.fig =  handle([]);
                options.ax =  handle([]);
            end
            if isempty(options.fig) && isempty(options.ax)
                [fig,ax] = vasplib_tool.create_figure();
            else
                fig = options.fig;
                ax = options.ax;
            end
            
            if isempty(MA.klist_l)
                cla(ax);
                if ~strcmp(options.savefile,'')
                    exportgraphics(ax,options.savefile);
                    cla(ax);
                end
                return
                
            end
            Kpoints_name = strrep( MA.kpoints_name,'GAMMA','G');
            [fig,ax] = bandplot(...
                MA.EIGENCAR,...
                options.Ecut,...
                MA.klist_l,...
                MA.kpoints_l,...
                Kpoints_name,...
                'Color',options.color,...
                'fig',fig,...
                'ax',ax,...
                'title',options.title ...
                );
            if ~strcmp(options.savefile,'')
                exportgraphics(ax,options.savefile);
                cla(ax);
            end
            
            
        end
    end 
    methods(Static)
        function [sym_mat,SG] = sym_mat_gen(filename,startRow, endRow)
            %% --------------------------------------------------------
            delimiter = {',',':'};
            if nargin<=2
                startRow = 3;
                endRow = 3;
            end
            formatSpec = '%*s%f%*s%*s%[^\n\r]';
            fileID = fopen(filename,'r');
            dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines', startRow(1)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
            for block=2:length(startRow)
                frewind(fileID);
                dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines', startRow(block)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
                dataArray{1} = [dataArray{1};dataArrayBlock{1}];
            end
            fclose(fileID);
            SG = [dataArray{1:end-1}];
            %% --------------------------------------------------------
            delimiter = {',','[',']','#'};
            if nargin <=3
                startRow = 6;
                endRow = inf;
            end
            formatSpec = '%s%s%s%s%[^\n\r]';
            fileID = fopen(filename,'r');
            
            dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines', startRow(1)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
            for block=2:length(startRow)
                frewind(fileID);
                dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines', startRow(block)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
                for col=1:length(dataArray)
                    dataArray{col} = [dataArray{col};dataArrayBlock{col}];
                end
            end
            fclose(fileID);
            raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
            for col=1:length(dataArray)-1
                raw(1:length(dataArray{col}),col) = mat2cell(dataArray{col}, ones(length(dataArray{col}), 1));
            end
            numericData = NaN(size(dataArray{1},1),size(dataArray,2));
            for col=[2,3,4]
                rawData = dataArray{col};
                for row=1:size(rawData, 1)
                    regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
                    try
                        result = regexp(rawData(row), regexstr, 'names');
                        numbers = result.numbers;
                        invalidThousandsSeparator = false;
                        if numbers.contains(',')
                            thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
                            if isempty(regexp(numbers, thousandsRegExp, 'once'))
                                numbers = NaN;
                                invalidThousandsSeparator = true;
                            end
                        end
                        if ~invalidThousandsSeparator
                            numbers = textscan(char(strrep(numbers, ',', '')), '%f');
                            numericData(row, col) = numbers{1};
                            raw{row, col} = numbers{1};
                        end
                    catch
                        raw{row, col} = rawData{row};
                    end
                end
            end
            rawNumericColumns = raw(:, [2,3,4]);
            rawStringColumns = string(raw(:, 1));
            R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),rawNumericColumns);
            rawNumericColumns(R) = {NaN};
            sym_operation = table;
            sym_operation.Label_name = rawStringColumns(:, 1);
            sym_operation.rc1 = cell2mat(rawNumericColumns(:, 1));
            sym_operation.rc2 = cell2mat(rawNumericColumns(:, 2));
            sym_operation.rc3 = cell2mat(rawNumericColumns(:, 3));
            % sym_operation= readtable(filename, opts);
            %
            cut_label = find([sym_operation.Label_name]'=='atom_mapping:');
            %
            rc1 = [sym_operation.rc1];
            rc2 = [sym_operation.rc2];
            rc3 = [sym_operation.rc3];
            sym_operation_list = [rc1,rc2,rc3];
            sym_operation_list(cut_label:end,:) = [];
            operation_num = size(sym_operation_list,1)/5;
            sym_mat = zeros(4,3,operation_num);
            for i = 1:operation_num
                cut1 = (i-1)*5 +2;
                cut2 = (i-1)*5 +5;
                sym_mat(:,:,i) = sym_operation_list(cut1:cut2,:);
            end        
        end            
        function position = position_gen(sites)
            position_init = [[sites.rc1].',[sites.rc2].',[sites.rc3].'];
            position = sortrows(position_init,[3,1,2]) ;
        end

    end
    
end

