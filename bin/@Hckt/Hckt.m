classdef Hckt < matlab.mixin.CustomDisplay
    %UNTITLED 此处提供此类的摘要
    % suppose GND is the last
    properties
        %dim = 2;
        title = "Created by Hckt "+date;
        Nports ; % which is the homecell Nnodes except GND
        MeshNnodes ; % which is the netlist Nnodes after genmesh
        vectorAll = [0,0]; % the first line is always homecell
    end
    properties % with seq
        HnodeL  logical= [];
        ScktL Subckt = Subckt([]);
        vectorL = [1];
        PortInCell = {[]}; % The subckt must IN - OUT
        PortOutCell = {[]};
        DescriptionL = "";
    end
    % netlist mesh
    properties
        NetlistMesh;
        fin_dir = [];
        mesh = [];
        Port2VectorDist ;
        Port2PortDist ;
        Port2PortDistForVectorise ;
        Port2PortDist_BinBin ;
        AddtionInformation ='';
    end
    % other
    properties
        magnitude = 'p';
        Lib = [".inc 'Hckt.lib' ";];
        IC = ["";];
        Ipulse = ["";];
        Vac =["";];
        Options = [...
            ".option post=2 probe";...
            "*.option parhier=global";...
            ".option accurate=1";...
            "*.option INTERP";]
    end
    % Vectorise
    properties
        VectorizeScktL Subckt = Subckt([]) ;% attention the number of this ScktL depends on fin_dir;
    end
    properties
        pinned = false;
    end
    properties(Dependent)
        NRPTS;
        HomeCell Subckt  ;
        dim;
        Components;
        innerNetlist;
        NetlistStr;
        Netlist;
        pinnedNetlist;
        OPT;
    end
    properties
        pinnedNetlistStr;
    end
    % data
    properties(Transient =true)

    end
    %% Define which properties show
    methods (Access = protected)
        function propgrp = getPropertyGroups(~)
            proplist = {'title','dim','Nports','Components','vectorAll','Netlist'};
            propgrp = matlab.mixin.util.PropertyGroup(proplist);
        end
    end
    %% Construction
    methods
        function Hcktobj = Hckt(ScktL,vectorAll,propArgs,options)
            arguments
                ScktL = 1;
                vectorAll = [0,0];
                propArgs.?Hckt;%  https://ww2.mathworks.cn/help/matlab/matlab_prog/function-argument-validation-1.html#mw_1b62b6d6-a445-4c55-a9b9-9c70becfdbe6
                options.test = true;
            end
            if ScktL ~=1
                propArgs.HomeCell = ScktL(1);
                propArgs.HoppingScktL = ScktL(2:end);
            else
            end
            if ~isequal(vectorAll,[0,0])
                propArgs.vectorAll = vectorAll;
            else
            end
            fieldnamesL = fieldnames(propArgs);
            if ~isempty(fieldnamesL)
                for i =1:numel(fieldnamesL)
                    Hcktobj.(fieldnamesL{i}) = propArgs.(fieldnamesL{i});
                end
                % old bug; delete later
                try
                    Hcktobj.vectorAll = propArgs.vectorL;
                catch 
                    
                end
            end

            Hcktobj.Port2PortDist = containers.Map(0,0);
            Hcktobj.Port2PortDistForVectorise = containers.Map(0,0);
            Hcktobj.Port2VectorDist = containers.Map(0,zeros(1,Hcktobj.dim));
            Hcktobj.Port2PortDist_BinBin = Hcktobj.Port2PortDist;
        end
    end
    %% set hop
    methods
        function Hcktobj = set_home(Hcktobj,Subcktobj,PortInL,PortOutL,DescriptionL)
            arguments
                Hcktobj Hckt;
                Subcktobj Subckt;
                PortInL = 1;
                PortOutL = 2;
                DescriptionL = Subcktobj.Description;
            end
            Hcktobj = Hcktobj.set_hop(zeros(1,Hcktobj.dim),Subcktobj,PortInL,PortOutL,DescriptionL);
        end
        function Hcktobj = set_hop(Hcktobj,vector,Subcktobj,PortInL,PortOutL,DescriptionL,options)
            arguments
                Hcktobj Hckt;
                vector;
                Subcktobj Subckt;
                PortInL = 1;
                PortOutL = 2;
                DescriptionL = Subcktobj.Description;
                options.Port2PortMode = false;
            end
            [~,vectorlabel] = ismember(vector,Hcktobj.vectorAll,'rows');
            if vectorlabel == 0
                Hcktobj = add_empty_one(Hcktobj,vector);
                vectorlabel = size(Hcktobj.vectorAll,1);
            end
            %
            if vectorlabel == 1 % homecell
                seq =1;
            else
                seq = Hcktobj.NRPTS+1;
            end
            %disp(seq)
            Hcktobj.ScktL(seq,1) = Subcktobj;
            Hcktobj.vectorL(seq) = vectorlabel;
            Hcktobj.PortInCell{seq,1} = PortInL;
            Hcktobj.PortOutCell{seq,1}= PortOutL;
            Hcktobj.DescriptionL{seq,1} = char(DescriptionL);
            % temp list mode
            % temp one port 2 another
            % bugs here!
            if options.Port2PortMode
                checkvectorL = Hckt.dim2vectorL(Hcktobj.dim);
                if ~all(vector==0)
                    if ismember(vector,checkvectorL,'rows')
                        %
                        Hcktobj.ScktL(seq,1).ReverseConnection = false;
                        %
                        Hcktobj.Port2VectorDist(PortInL) = vector;
                        Hcktobj.Port2PortDist(PortInL) = PortOutL;
                        Hcktobj.Port2PortDistForVectorise(PortInL) = PortOutL;
                        % BinBin convention
                        Hcktobj.Port2PortDist_BinBin(PortInL) = PortInL;
                    else
                        %
                        Hcktobj.ScktL(seq,1).ReverseConnection = true;
                        %
                        Hcktobj.Port2VectorDist(PortInL) =  zeros(1,Hcktobj.dim);
                        Hcktobj.Port2PortDist(PortInL) = PortInL;
                        Hcktobj.Port2PortDistForVectorise(PortInL) = PortOutL;
                        % BinBin convention
                        Hcktobj.Port2PortDist_BinBin(PortInL) = PortOutL;
                    end
                end
            end
        end
        function Hcktobj = add_empty_one(Hcktobj,vector)
            Hcktobj.vectorAll = [Hcktobj.vectorAll ;vector];
        end
    end
    %% get
    methods
        function OPT = get.OPT(HcktObj)
            switch  HcktObj.magnitude
                case 'p'
                    TRAN =".tran 1ns 10us";
                    DEFAULT_PARM = "%.PARAM InitV =0V CA = 100p C_hopping = 100p C_v = 100p";
                case 'n'
                    TRAN =".tran 50ns 500us";
                    DEFAULT_PARM = "%.PARAM InitV =0V CA = 100n C_hopping = 100n C_v = 100n";
                case 'u'
                    TRAN =".tran 1us 10ms";
                    DEFAULT_PARM = "%.PARAM InitV =0V CA = 100u C_hopping = 100u C_v = 100u";
            end
            OPT = [TRAN;HcktObj.Options;DEFAULT_PARM];
        end
        function HomeCell = get.HomeCell(HcktObj)
            HomeCell = HcktObj.ScktL(1);
        end
        function Components = get.Components(HcktObj)
            Components = HcktObj.NRPTS;
        end
        function NRPTS = get.NRPTS(HcktObj)
            NRPTS = length(HcktObj.ScktL); % Todo: use size
        end
        function dim = get.dim(HcktObj)
            dim = size(HcktObj.vectorAll,2);
        end
        function innerNetlist = get.innerNetlist(HcktObj)

        end
        function NetlistStr = get.NetlistStr(HcktObj)
            if options.pinned
                NetlistStr = HcktObj.pinnedNetlistStr;
            else
                NetlistStr = mat2str(HcktObj.Netlist);
            end
        end
        function pinnedNetlist = get.pinnedNetlist(HcktObj)
            pinnedNetlist = park.char2list(HcktObj.pinnedNetlistStr);
        end
        function Netlist = get.Netlist(HcktObj)
            if HcktObj.pinned
                Netlist = HcktObj.pinnedNetlist;
            else
                
                Netlist = HcktObj.innerNetlist;
            end
        end
    end
    %% I/O
    methods(Static)
        function [simulation_result]=read_hspice_ac(filename,options)
            % run 'hspiceAC' first on linux
            % need case.ac0.info & case.ac0
            arguments
                filename = "Hckt.ac0";
                options.fast = true;
                options.filenameInformation = '';
            end
            optionscell = namedargs2cell(options);
            [simulation_result]=Hckt.read_hspice_tr(filename,optionscell{:});
            
        end
        function [simulation_result]=read_hspice_tr(filename,options)
            % run 'hspiceTR' first on linux
            % need case.tr0.info & case.tr0
            arguments
                filename = "Hckt.tr0";
                options.fast = true;
                options.filenameInformation = '';
            end
            if options.fast
                [~,~,ext] = fileparts(filename);%clear filepath name
                file_extension=ext(1:3); %clear ext
                if strcmp(options.filenameInformation,'')
                    filenameInformation = char(filename);
                    filenameInformation = [filenameInformation,'.info'];
                else
                    filenameInformation = options.filenameInformation;
                end
                % info
                fileID = fopen(filenameInformation);
                content_file = textscan(fileID,'%s','HeaderLines',3);
                fclose(fileID);
                % find number of the variables
                content_str=string(content_file{1,1});
                data_start_ind=find(content_str==('$&%#'))+1;
                num_var=find(isnan(str2double(content_str(1:data_start_ind-2))),1)-1;
                simulation_result = Hckt.save_signal_names(num_var,data_start_ind,content_str,file_extension);
                % data read
                fileID = fopen(filename);
                formatSpec = '%13f%13f%13f%13f%13f%13f%[^\n\r]';
                %formatSpec = '%13c%13c%13c%13c%13c%13c%[^\n\r]';
                dataArray = textscan(fileID, formatSpec);
                % check dataarray
                checkEndRow = 6;
                for i = 1:6
                    if isnan(dataArray{i}(end))
                        checkEndRow = i-1;
                        break;
                    end
                end
                iTOT = length(dataArray{1});
                TOTnum = iTOT*6;
                %DATA = char(TOTnum,13);
                DATA = zeros(iTOT,6);
                count = 1;
                for i = 1:6
                    if i < checkEndRow
                        DATA(1:iTOT,i) = (dataArray{i}(1:end));
                        count = count+iTOT;
                        %disp('full');
                        %disp(count);
                    else
                        DATA(1:(iTOT-1),i) = (dataArray{i}(1:end-1));
                        count = count+iTOT-1;
                        %disp(count);
                    end
                end
                %count = count;
                DATA = DATA.';
                DATA = DATA(:);
                DATA(count:end)= [];
                %DATA = double(DATA(data,[data_frmt 'E[+-]..'],'match'));
                %NDATA = numel(DATA);
                NStruct= length(simulation_result);
                DATA = reshape(DATA,NStruct,numel(DATA)/NStruct).';
                for i = 1:length(simulation_result)
                    simulation_result(i).val = DATA(:,i);
                end
                fclose(fileID);
            else
                [simulation_result]= Hckt.read_hspice_tr_sw_ac(filename);

            end
        end
        function [simulation_result]=read_hspice_tr_sw_ac(filename)
            % Author
            %Mohammad Abu Raihan Miah
            %University of California, San Diego
            %ver 4.0.0, 11/23/21
            %This function reads a ASCII formatted HSPICE output file (.option post=2)
            %'filename.tr#' or 'filename.sw#' or 'filename.ac#' and saves all the
            %signals in the variable 'simulation_result' as a structure.
            %
            %simulation_result=a structure (2 fields) that contains the contents from
            %                  the HSPICE file.
            %simulation_result(#).var_name=name of the signals present in the file.
            %simulation_result(#).val=a vector contains the values of the signal with
            %                         the name simulation_result(#).var_name
            % determine the file extension
            [filepath,filetitle,ext] = fileparts(filename);%clear filepath name
            file_extension=ext(1:3); %clear ext
            % read data from the file
            file_handle = fopen(filename);


            content_file = textscan(file_handle,'%s','HeaderLines',3);
            fclose(file_handle);

            % find number of the variables and start point of the data
            content_str=string(content_file{1,1});clear content_file
            data_start_ind=find(content_str==('$&%#'))+1;
            num_var=find(isnan(str2double(content_str(1:data_start_ind-2))),1)-1;

            % save name of the signals in the .var_name field
            simulation_result = Hckt.save_signal_names(num_var,data_start_ind,content_str,file_extension);

            % save signal data in the .val field
            data=replace(strjoin(content_str(data_start_ind:end)),' ','');
            data_frmt=char(ones(1,find(char(content_str(data_start_ind))=='E',1)-1)*'.');
            data_separate=double(regexp(data,[data_frmt 'E[+-]..'],'match'));clear data
            for ii=1:length(simulation_result)
                simulation_result(ii).val=data_separate(ii:length(simulation_result):length(data_separate)-1)';
            end
        end
    end
    methods(Static)
        function simulation_result = save_signal_names(num_var,data_start_ind,content_str,file_extension)
            var_name_raw_ind=num_var+2:data_start_ind-2;
            temp=find(max(max((char(content_str(var_name_raw_ind))=='(')*2,1),[],2)>1);
            var_name_1st_part=content_str(num_var+1+temp);
            var_name_special=setdiff(var_name_raw_ind,var_name_raw_ind(temp));
            var_name_last_part(1:length(temp),1)="";
            for iii=1:length(var_name_special)
                tt=find(~(var_name_special(iii)>var_name_raw_ind(temp)),1)-1;
                var_name_last_part(tt)=content_str(var_name_raw_ind(tt+1));
            end
            var_name_raw=[content_str(num_var+1) var_name_1st_part'+var_name_last_part'+")"]';

            % add _real and _imag with the relevant variable names from .ac#
            % files
            if file_extension == '.tr' | file_extension == '.sw'
                simulation_result=cell2struct(cellstr(var_name_raw'),'var_name',1);
            elseif file_extension == '.ac'
                var_real_imag=find(rem(str2double(content_str(1:num_var)),7)==1);
                if isempty(var_real_imag)
                    var_name=var_name_raw;
                else
                    count=1;
                    for iii=1:length(var_name_raw)
                        if any(iii==var_real_imag)
                            var_name(count)=var_name_raw(iii)+"_real";count=count+1;
                            var_name(count)=var_name_raw(iii)+"_imag";
                        else
                            var_name(count)=var_name_raw(iii);
                        end
                        count=count+1;
                    end
                end
                simulation_result=cell2struct(cellstr(var_name'),'var_name',1);
            else
                error(['The toolbox cannot handle' file_extension '# files']);
            end
        end
    end
    %% DataClean
    methods(Static)
        function [VectorList,ObservationsMat,Vstruct,SpectrumX,TimeL] = extractObservations(simulation_result,Observations,options)
            arguments
                simulation_result struct;
                Observations = 'v';
                options.analysis {mustBeMember(options.analysis,{'tran','ac'})} = 'tran';
                options.mode = 'simple'; %Todo  temp
                options.x = true;
                options.fft = true;
                options.patten = "_"|".";
                options.Vstrcut = false;
                options.SelectNode = [];
            end
            chooseL = ones(1,length(simulation_result));
            for i = 1:numel(simulation_result)
                if ~strcmp(simulation_result(i).var_name(1),Observations)
                    chooseL(i) = false;
                elseif strcmp(simulation_result(i).var_name ,[Observations,'(0)'])
                    chooseL(i) = false;
                elseif ~options.x
                    if contains(simulation_result(i).var_name ,'x')
                        chooseL(i) = false;
                    end
                else

                end
            end
            switch  options.analysis
                case 'tran'
                    TimeL =  simulation_result(1).val;
                    DurationTime = max(TimeL);
                    TimeListN = round(TimeL*(1/DurationTime)*length(TimeL));
                    SpectrumX = (1:length(TimeL))/(DurationTime);
                    simulation_result_Observations = simulation_result(logical(chooseL));
                    VectorList = [];
                    for i = 1:size(simulation_result_Observations,1)
                        tmpChar = simulation_result_Observations(i).var_name;
                        tmpChar(1:2) = [];
                        tmpChar(end) = [];
                        tmlStrL = split(tmpChar,options.patten).';
                        tmlStrL(1) = [];
                        try
                            VectorList = [VectorList;double([string(tmlStrL)])];
                        catch
                            VectorList = [VectorList;zeros(size(VectorList(1,:)))];
                        end
                    end
                    % select node
                    if isempty(options.SelectNode)
                    else
                        SelectL = VectorList(:,end) == options.SelectNode;
                        VectorList = VectorList(SelectL,:);
                        simulation_result_Observations = simulation_result_Observations(SelectL);
                    end
                    %
                    if options.Vstrcut
                        ObservationsMat = [];
                        for i = 1:size(VectorList,1)
                            Vstruct(i).R = VectorList(i,1:end-1);
                            Vstruct(i).port = VectorList(i,end);
                            Vstruct(i).val = simulation_result_Observations(i).val;
                            Vstruct(i).fftval = (nufft(simulation_result_Observations(i).val.',TimeListN));% abs ?
                            if options.fft
                                ObservationsMat = [ObservationsMat;Vstruct(i).fftval];
                            else
                                ObservationsMat = [ObservationsMat;Vstruct(i).val.' ];
                            end
                        end
                    else
                        Vstruct = [];
                        ObservationsMat = zeros(size(VectorList,1),length(TimeL));
                        if options.fft
                            for i = 1:size(VectorList,1)
                                ObservationsMat(i,:) = (nufft(simulation_result_Observations(i).val.',TimeListN));
                            end
                        else
                            for i = 1:size(VectorList,1)
                                ObservationsMat(i,:) = simulation_result_Observations(i).val.' ;
                            end
                        end
                    end
                case 'ac'
                    SpectrumX =  simulation_result(1).val;
                    TimeL = [];
                    simulation_result_Observations = simulation_result(logical(chooseL));
                    VectorList = [];
                    NPROBE = size(simulation_result_Observations,1);
                    IsRealL = ones(NPROBE,1);
                    for i = 1:NPROBE
                        tmpChar = simulation_result_Observations(i).var_name;
                        if tmpChar(2) == 'i'
                            IsRealL(i) = 1i;
                        end
                        tmpChar([1:3]) = [];
                        tmpChar(end) = [];
                        tmlStrL = split(tmpChar,options.patten).';
                        tmlStrL(1) = [];
                        try
                            VectorList = [VectorList;double([string(tmlStrL)])];
                        catch
                            VectorList = [VectorList;zeros(size(VectorList(1,:)))];
                        end
                    end
                    % select node
                    if isempty(options.SelectNode)
                    else
                        SelectL = VectorList(:,end) == options.SelectNode;
                        VectorList = VectorList(SelectL,:);
                        IsRealL = IsRealL(SelectL);
                        simulation_result_Observations = simulation_result_Observations(SelectL);
                    end
                    %
                    Vstruct = [];
                    ObservationsMat = zeros(size(VectorList,1),length(SpectrumX));
                    for i = 1:size(VectorList,1)
                        ObservationsMat(i,:) = IsRealL(i) * simulation_result_Observations(i).val.' ;
                    end
                    [VectorList,ObservationsMat] = HollowKnight.generalcontractrow2(VectorList,ObservationsMat);
                otherwise

            end
        end
        function EIGENCAR = ProjectDOSCAR(DOSCAR,options)
            arguments
                DOSCAR;
                options.POSCAR = 'POSCAR';
                options.KPOINTS = 'KPOINTS';
                options.dir_seq = [1,2,3];
            end
            %
            if exist(options.KPOINTS,'file') && exist(options.POSCAR,'file')
                Rm=POSCAR_readin(options.POSCAR);
                [~,klist_l,klist_s,~,~]=kpathgen3D(Rm,options.KPOINTS);
                %[~,klist_l,~,kpoints_l,kpoints_name]=kpathgen3D(Rm,options.KPOINTS);
            else
                error('No KPOINTS POSCAR');
            end
            sizemesh = size(DOSCAR);
            Dim = length(sizemesh)-1;
            klist_s = mod(klist_s,1);
            if Dim == 1
                EIGENCAR = zeros(sizemesh(1),length(klist_l));
                klist1_s = linspace(0,1,sizemesh(2));
                %[K1grid,K2grid] = meshgrid(klist1_s,klist2_s);
                for i = 1:sizemesh(1)
                    EIGENCAR(i,:) = interp1(klist1_s,DOSCAR(i,:),klist_s(:,options.dir_seq(1)));
                end
            elseif Dim == 2
                EIGENCAR = zeros(sizemesh(3),length(klist_l));
                klist1_s = linspace(0,1,sizemesh(1));
                klist2_s = linspace(0,1,sizemesh(2));
                [K1grid,K2grid] = meshgrid(klist1_s,klist2_s);
                for i = 1:sizemesh(3)
                    EIGENCAR(i,:) = interp2(K1grid,K2grid,DOSCAR(:,:,i),klist_s(:,options.dir_seq(1)),klist_s(:,options.dir_seq(2)));
                end
            elseif Dim == 3
                EIGENCAR = zeros(sizemesh(4),length(klist_l));
                klist1_s = linspace(0,1,sizemesh(1));
                klist2_s = linspace(0,1,sizemesh(2));
                klist3_s = linspace(0,1,sizemesh(3));
                [K1grid,K2grid,K3grid] = meshgrid(klist1_s,klist2_s,klist3_s);
                %P = [2 1 3];
                %K1grid = permute(K1grid, P);
                %K2grid = permute(K2grid, P);
                %K3grid = permute(K3grid, P);
                for i = 1:sizemesh(4)
                    %V = permute(DOSCAR(:,:,:,i), P);
                    %V = permute(DOSCAR(:,:,:,i),[1,2,3]);
                    V = DOSCAR(:,:,:,i);
                    EIGENCAR(i,:) = interp3(klist1_s,klist2_s,klist3_s,V,klist_s(:,options.dir_seq(1)),klist_s(:,options.dir_seq(2)),klist_s(:,options.dir_seq(3))).';
                    %EIGENCAR(i,:) = interpn(K1grid,K2grid,K3grid,V,klist_s(:,1),klist_s(:,2),klist_s(:,3));
                end
            else

            end
        end
        function [EIGENCAR] = CollectVstruct1D(ObservationsMat)
            % need fix
            %[VectorList,ObservationsMat]=HollowKnight.generalcontractrow2(VectorList(:,1:end-1),ObservationsMat);
            %DurationTime = max(TimeL);
            %redundancy = 1;
            %TimeListN = round(TimeL*(1/DurationTime)*length(TimeL)*redundancy);
            %SpectrumX = (1:length(length(TimeL)))/(DurationTime);
            %Xlist = 1:length(SelectL);
            %TkCar = fft(ObservationsMat(SelectL,:).',[],2);
            OmegaCar = ObservationsMat.';
            %EIGNECAR = abs(nufftn(ObservationsMat(SelectL,:).',{TimeListN,Xlist}));
            %EIGNECAR = reshape(EIGNECAR,[],length(SelectL));
            EIGENCAR = abs(fft((OmegaCar),[],2));
            %EIGNECAR = abs(fft(real(OmegaCar),[],2) +ifft(imag(OmegaCar),[],2));
            %EIGNECAR = abs(OmegaCar);
            %EIGNECAR = abs(fft2(ObservationsMat(SelectL,:))).';
            EIGENCAR = EIGENCAR(:,[1:end,1]);
        end
        function [DOSCAR_3D,klist1,klist2,OmgL] = CollectVstruct2D(VectorList,ObservationsMat,OmegaCut,SpectrumX)
            SelectOmegaL = find(SpectrumX >= OmegaCut(1) & SpectrumX <= OmegaCut(2));
            % mesh X & Y
            %NSizeMesh = size(VectorList,1);
            Rvector = VectorList(:,1:end-1);
            NX = max(Rvector(:,1));
            NY = max(Rvector(:,2));
            sizemesh = [NX,NY];
            % Todo: compatible with Cell
            klist1 = 1:NX+1;
            klist2 = 1:NY+1;
            OmgL = SpectrumX(SelectOmegaL);
            ObservationsMat = ObservationsMat(:,SelectOmegaL);
            %Xmesh = reshape(Rvector(:,1),sizemesh(1),sizemesh(2));
            %Ymesh = reshape(Rvector(:,2),sizemesh(1),sizemesh(2));
            DOSCAR_3D = zeros(sizemesh(1),sizemesh(2),length(OmgL));
            for i = 1:numel(SelectOmegaL)
                DOSCAR_3D(:,:,i) = abs(fft2(reshape(ObservationsMat(:,i),NX,NY)));
            end
            DOSCAR_3D(sizemesh(1)+1,:,:) = DOSCAR_3D(1,:,:);
            DOSCAR_3D(:,sizemesh(2)+1,:) = DOSCAR_3D(:,1,:);
            %
            % R = [Vstruct
        end
        function [DOSCAR_4D,klist1,klist2,klist3,OmgL] = CollectVstruct3D(VectorList,ObservationsMat,OmegaCut,SpectrumX,options)
            arguments
                VectorList
                ObservationsMat
                OmegaCut
                SpectrumX
                options.fin_dir = [0,0,0];
            end
            SelectOmegaL = find(SpectrumX >= OmegaCut(1) & SpectrumX <= OmegaCut(2));
            % mesh X & Y
            %NSizeMesh = size(VectorList,1);
            Rvector = VectorList(:,1:end-1);
            if options.fin_dir(1)
                NX = 1;
                klist1 = 1;
            else
                NX = max(Rvector(:,1));
                klist1 = 1:NX+1;
            end
            if options.fin_dir(2)
                NY = 1;
                klist2 = 1
            else
                NY = max(Rvector(:,2));
                klist2 = 1:NY+1;
            end
            if options.fin_dir(2)
                NZ = 1;
                klist3 = 1;
            else
                NZ = max(Rvector(:,3));
                klist3 = 1:NZ+1;
            end
            sizemesh = [NX,NY,NZ];
            % Todo: compatible with Cell
            OmgL = SpectrumX(SelectOmegaL);
            ObservationsMat = ObservationsMat(:,SelectOmegaL);
            %Xmesh = reshape(Rvector(:,1),sizemesh(1),sizemesh(2));
            %Ymesh = reshape(Rvector(:,2),sizemesh(1),sizemesh(2));
            DOSCAR_4D = zeros(sizemesh(1),sizemesh(2),sizemesh(3),length(OmgL));
            for i = 1:numel(SelectOmegaL)
                tmpDATA = abs(fftn(reshape(ObservationsMat(:,i),NZ,NY,NX)));
                DOSCAR_4D(:,:,:,i) = permute(tmpDATA,[3,2,1]);
            end
            if ~options.fin_dir(1)
                DOSCAR_4D(sizemesh(1)+1,:,:,:) = DOSCAR_4D(1,:,:,:);
            end
            if ~options.fin_dir(2)
                DOSCAR_4D(:,sizemesh(2)+1,:,:) = DOSCAR_4D(:,1,:,:);
            end
            if ~options.fin_dir(3)
                DOSCAR_4D(:,:,sizemesh(3)+1,:) = DOSCAR_4D(:,:,1,:);
            end
            %
            % R = [Vstruct
        end
    end
    %% plot
    methods(Static)
        function [fig,ax] = Frequencyplot(ObservationsMat,SpectrumL,OmegaCut,options)
            arguments
                ObservationsMat;
                SpectrumL;
                OmegaCut = [0,2e4];
                options.fig =  handle([]);
                options.ax =  handle([]);
                options.FontName = 'Helvetica';
                options.FontSize = 24;
                options.Units = 'pixels';
                options.Position = [];
                options.Color =  [rand rand rand];
                options.title = '';
                options.xlabel='Frequency (Hz)';
                options.ylabel='';
            end
            if isempty(options.fig) && isempty(options.ax)
                [fig,ax] = vasplib_tool.creat_figure('Units',options.Units,'Position',options.Position);
            else
                fig = options.fig;
                ax = options.ax;
            end
            for i = 1:size(ObservationsMat,1)
                plot(SpectrumL,ObservationsMat(i,:),'DisplayName',num2str(i),'Color',options.Color);
                hold on
            end
            %--------  title  -------
            title(ax,options.title);
            %
            xlabel(ax,options.xlabel);
            ylabel(ax,options.ylabel);
            xlim(OmegaCut);
        end
        function [ax,WaveFunc] = waveplot(ObservationsMat,VectorList,SpectrumL,orbL,options_select,options)
            arguments
                ObservationsMat = [];
                VectorList = [];
                SpectrumL = [];
                orbL = [];
                options_select.Frequency = -1;
                options_select.Width = 0;
                options_select.scale = 1;
                options.ax =  handle([]);
                options.Rm = [];
                options.POSCAR = 'POSCAR';
                options.WaveMin = 1e-3;
                options.WaveColor = 'r';
                options.WaveSize = 1;
                options.OrbColor = 'k';
                options.OrbSize = 1;
            end
            %
            optionscell = namedargs2cell(options);
            %
            if options_select.Frequency == -1
                options_select.Frequency = SpectrumL(end/2);
            end
            if options_select.Width == 0
                Scale = round(log(options_select.Frequency)/log(10));
                options_select.Width = 10^(Scale-3);
            end
            %
            OmegaCut = [options_select.Frequency - options_select.Width,options_select.Frequency + options_select.Width];
            ChooseL = SpectrumL >= OmegaCut(1) & SpectrumL <= OmegaCut(2);
            % For ABS
            [VectorList,ObservationsMat] = HollowKnight.generalcontractrow2(VectorList,ObservationsMat);
            ObservationsMat = abs(ObservationsMat);
            % Rvector enforce three
            NVectorList = size(VectorList,1);
            RvectorL = VectorList(:,1:end-1);
            NodeL = VectorList(:,end);
            Dim = size(RvectorL,2);
            if Dim <3
                RvectorL = [RvectorL,zeros(NVectorList,3-Dim)];
            end
            % True orbL
            littleorbL = orbL(NodeL,:);
            ORBL = RvectorL+littleorbL;
            % WaveFunc
            WaveFunc = normalize(sum(ObservationsMat(:,ChooseL),2),'range',[0,1])*options_select.scale;
            % 
            ax = vasplib_plot.waveplot(ORBL,WaveFunc,optionscell{:});
        end
        function [ax] = bandplot(EIGNECAR,OmegaCut,SpectrumL,options)
            arguments
                EIGNECAR = [];
                OmegaCut = [0 20000];
                SpectrumL = [];
                options.ax =  handle([]);
                options.FontName = 'Helvetica';
                options.FontSize = 24;
                options.Units = 'pixels';
                options.Position = [];
                options.Color =  parula;
                options.title = '';
                options.xlabel='';
                options.ylabel='Frequency (Hz)';
                options.alpha = 1;
                options.KPOINTS = 'KPOINTS';
                options.POSCAR = 'POSCAR';
                options.shift = false;
            end
            if isempty(options.ax)
                Fig = vasplib_plot.create_figure('Position',options.Position);
                ax = Fig.axes;
            else
                ax = options.ax;
            end
            %
            if exist(options.KPOINTS,'file') && exist(options.POSCAR,'file')
                Rm=POSCAR_readin(options.POSCAR);
                [~,klist_l,~,kpoints_l,kpoints_name]=kpathgen3D(Rm,options.KPOINTS);
                UseKpath = true;
                if options.shift
                    X = linspace(min(klist_l),max(klist_l),size(EIGNECAR,2)+1);
                else
                    X = klist_l;
                end
            else
                UseKpath = false;
                if options.shift
                    X = [1:size(EIGNECAR,2),size(EIGNECAR,2)+1];
                else
                    X = [1:size(EIGNECAR,2),size(EIGNECAR,2)];
                end
            end
            % 
            ChooseL = SpectrumL>OmegaCut(1) & SpectrumL<OmegaCut(2);
            Y = SpectrumL(ChooseL);
            if options.shift
                Z = EIGNECAR(ChooseL,[1:size(EIGNECAR,2),1]);
            else
                Z = EIGNECAR(ChooseL,:);
            end
            [Xgrid,Ygrid] = meshgrid(X,Y);
            surf(ax,Xgrid,Ygrid,Z,'EdgeColor','none','FaceAlpha',options.alpha);
            axis(ax,'tight');
            view(ax,2);
            colormap(options.Color);
            %--------  title  -------
            title(ax,options.title);
            %
            if UseKpath
                Xcut = [min(X),max(X)];
                Ycut = [min(Y),max(Y)];
                %--------reference -------
                ax= vasplib_plot.set_reference(kpoints_l,kpoints_name,Xcut,Ycut,...
                    'ax',ax,...
                    'xlabel',options.xlabel,...
                    'ylabel',options.ylabel ...
                    );
            else
                xlabel(ax,options.xlabel);
                ylabel(ax,options.ylabel);
            end
        end
    end
    %% Gen
    methods
        function HcktObj = vectorize(HcktObj,fin_dir)
            arguments
                HcktObj
                fin_dir = zeros(HcktObj.dim,1);
            end
            basePriVname = 'PriV';
            for i = 1:HcktObj.dim
                basePriVname = [basePriVname,'_1'];
            end
            count = 0;
            checkvectorlist = HcktObj.dim2vectorL(HcktObj.dim);
            % here we dont search zero
            checkvectorlist(1,:) = [];
            % homecell
            count = count +1;
            netlistTmp{count} = strcat(HcktObj.ScktL(1).type,num2str(count));
            for i = 1:length(HcktObj.PortInCell{1})
                netlistTmp{count} = [...
                    netlistTmp{count},...
                    string(HcktObj.PortInCell{1}(i)),...
                    ];
            end
            for i = 1:length(HcktObj.PortInCell{1})
                netlistTmp{count} = [...
                    netlistTmp{count},...
                    string((HcktObj.PortOutCell{1}(i))),...
                    ];
            end
            if HcktObj.ScktL(1).ToGND
                netlistTmp{count} = [netlistTmp{count},'TOGND'];
            end
            netlistTmp{count} = [netlistTmp{count},string(HcktObj.ScktL(1).name)];
            % other
            for i = 2:size(HcktObj.ScktL,1)
                % check the port of vector
                tmpVector = HcktObj.Port2VectorDist(HcktObj.PortInCell{i});
                [~,seq] = ismember(tmpVector,checkvectorlist,'rows');
                if seq ~= 0
                    count = count +1;
                    netlistTmp{count} = strcat(HcktObj.ScktL(i).type,num2str(count));
                    for j = 1:length(HcktObj.PortInCell{i})
                        netlistTmp{count} = [...
                            netlistTmp{count},...
                            "n"+string(HcktObj.Port2PortDistForVectorise(HcktObj.PortOutCell{i}(j))),...
                            string(HcktObj.PortInCell{i}(j))];
                    end
                    if HcktObj.ScktL(i).ToGND
                        netlistTmp{count} = [netlistTmp{count},'TOGND'];
                    end
                    netlistTmp{count} = [netlistTmp{count},string(HcktObj.ScktL(i).name)];
                end
            end
            % incell and outcell should be the same
            HcktObj.VectorizeScktL = Subckt('X',basePriVname,...
                ["n"+([string(HcktObj.PortInCell{1})]),([string(HcktObj.PortOutCell{1})]),...
                'TOGND'],...
                HcktObj.HomeCell.Description,netlistTmp);
            % <----- For open boundary -----> Giveup Directly check general
            % method
            if sum(fin_dir>0)
                AddingGNDScktL =sum(fin_dir>0);
            end
        end
        function [HcktObj,Basename] = hspice_gen(HcktObj,filename,ops,options,options2,options3,options4,options5)
            arguments
                HcktObj Hckt;
                filename = "";
                ops.checking = true;
                options.mesh = repmat(10,[HcktObj.dim,1]);
                options.fin_dir = zeros(HcktObj.dim,1);
                options2.probenode {mustBeMember(options2.probenode,{'Allnode','Selectnode'})}= 'Allnode';
                options2.prefix = "+ ";
                options2.father = "";
                options2.nodelist = [];
                options3.fft = false;
                options3.libmode = false;
                options4.ExciteNodes = 1;
                options4.mode {mustBeMember(options4.mode,{'general','vectorized','BinBin'})}= 'general';% finished -> general
                options5.analysis {mustBeMember(options5.analysis,{'tran','ac'})} = 'tran';
            end
            optionscell = namedargs2cell(options);
            options2cell = namedargs2cell(options2);
            options3cell = namedargs2cell(options3);
            options5cell = namedargs2cell(options5);
            % special node
            NodeStr  = 'n';
            Basename = HcktObj.title;
            for  i = 1:length(options.mesh)
                NodeStr = [NodeStr,'_',num2str(options.mesh(i))];
                Basename = [Basename,'_',num2str(options.mesh(i))];
            end
            if size(options4.ExciteNodes,2) == 1
                for i = 1:size(options4.ExciteNodes,1)
                    NodeStrList{i} = [NodeStr,'_',num2str(options4.ExciteNodes(i))];
                end
            elseif size(options4.ExciteNodes,2) >1
                for i = 1:size(options4.ExciteNodes,1)
                    NodeStr  = 'n';
                    for j = 1: size(options4.ExciteNodes,2)-1
                        NodeStr = [NodeStr,'_',num2str(options4.ExciteNodes(i,j))];
                    end
                    NodeStr = [NodeStr,'_',num2str(options4.ExciteNodes(i,j+1))];
                    NodeStrList{i} = NodeStr;
                end
            end
            % IC
            ICSTRING = ['.ic v(',NodeStrList{1},') = 1'];
            % Ipulse
            switch HcktObj.magnitude
                case 'p'
                    IpulseSTRING = ['Ipulse ',NodeStrList{1},' GND PU 0 1 5n 5n 50u'];
                case 'u'
                    IpulseSTRING = ['Ipulse ',NodeStrList{1},' GND PU 0 1 5n 5n 50m'];
                case 'm'
                    IpulseSTRING = ['Ipulse ',NodeStrList{1},' GND PU 0 1 5n 5n 500m'];
            end
            % Vac
            VacSTRING = "Vac source GND AC 1 0";
            VacSTRING = [VacSTRING;string(['R_for_ac ',NodeStrList{1},' source 100'])];
            for i = 2:length(NodeStrList)
                % IC
                ICSTRING = [ICSTRING;string(['.ic v(',NodeStrList{i},') = 1'])];
                % Ipulse
                switch HcktObj.magnitude
                    case 'p'
                        IpulseSTRING = [IpulseSTRING;string(['Ipulse ',NodeStrList{i},' GND PU 0 1 5n 5n 50u'])];
                    case 'u'
                        IpulseSTRING = [IpulseSTRING;string(['Ipulse ',NodeStrList{i},' GND PU 0 1 5n 5n 50m'])];
                    case 'm'
                        IpulseSTRING = [IpulseSTRING;string(['Ipulse ',NodeStrList{i},' GND PU 0 1 5n 5n 500m'])];
                end
                % Vac
                VacSTRING = "Vac source GND AC 1 0";
                VacSTRING = [VacSTRING;string(['R_for_ac ',NodeStrList{i},' source 100'])];
            end
            % 
            if strcmp(HcktObj.Vac,"")
                HcktObj.Vac = VacSTRING;
            end
            if strcmp(HcktObj.IC,"")
                HcktObj.IC = string(ICSTRING);
            end
            if strcmp(HcktObj.Ipulse,"")
                HcktObj.Ipulse = string(IpulseSTRING);
            end
            % filename
            if strcmp(filename,"")
                switch options5.analysis
                    case 'ac'
                        filename = strcat(Basename,'_AC.sp');
                    otherwise
                        filename = strcat(Basename,'.sp');
                end
            end
            %
            switch options4.mode 
                case 'vectorized'
                    HcktObj = hspice_gen_vectorized(HcktObj,filename,optionscell{:},options2cell{:},options3cell{:});
                case 'BinBin'
                    HcktObj = hspice_gen_vectorized(HcktObj,filename,optionscell{:},options2cell{:},options3cell{:},'BinBin',true);
                case 'general'
                    if ops.checking
                        HcktObj = HcktObj.autohermi();
                    end
                    HcktObj = HcktObj.half();
                    HcktObj = hspice_gen_general(HcktObj,filename,optionscell{:},options2cell{:},options3cell{:},options5cell{:});
            end
        end
        function optionslib_write(HcktObj,fid,options3,options5)
            arguments
                HcktObj
                fid
                options3.libmode = false;
                options5.analysis {mustBeMember(options5.analysis,{'tran','ac'})} = 'tran';
            end
            % options
            fprintf(fid,"* -------- Initial --------\n");
            for i = 1:size(HcktObj.IC,1)
                fprintf(fid,'* ');
                fprintf(fid,(HcktObj.IC(i,:)));fprintf(fid,"\n");
            end
            fprintf(fid,"* -------- OPTIONS --------\n");
            switch  HcktObj.magnitude
                case 'p'
                    TRAN =".tran 1ns 10us";
                    AC = ".AC LIN 1000 0.5e04 2e07";
                    DEFAULT_PARM = "%.PARAM InitV =0V CA = 100p C_hopping = 100p C_v = 100p";
                case 'n'
                    TRAN =".tran 50ns 500us";
                    AC = ".AC LIN 1000 0.5e04 1e06";
                    DEFAULT_PARM = "%.PARAM InitV =0V CA = 100n C_hopping = 100n C_v = 100n";
                case 'u'
                    TRAN =".tran 1us 10ms";
                    AC = ".AC LIN 1000 0.5e04 2e04";
                    DEFAULT_PARM = "%.PARAM InitV =0V CA = 100u C_hopping = 100u C_v = 100u";
            end
            switch options5.analysis
                case 'ac'
                    OPTcase = [AC;HcktObj.Options;DEFAULT_PARM];
                case 'tran'
                    OPTcase = [TRAN;HcktObj.Options;DEFAULT_PARM];
                otherwise
                    OPTcase = [TRAN;HcktObj.Options;DEFAULT_PARM];
            end
            for i = 1:size(OPTcase,1)
                fprintf(fid,(OPTcase(i,:)));fprintf(fid,"\n");
            end
            fprintf(fid,"* -------- -------- --------\n");
            % lib
            fprintf(fid,"* -------- lib --------\n");
            for i = 1:size(HcktObj.Lib,1)
                fprintf(fid,(HcktObj.Lib(i,:)));fprintf(fid,"\n");
            end
            fprintf(fid,"\n");
            fprintf(fid,"* -------- -------- --------\n");
            % Subckt
            fprintf(fid,"* -------- Subckt for Homecell --------\n");
            fprintf(fid,"*\n");
            fprintf(fid,HcktObj.HomeCell.OutputStr);
            fprintf(fid,"\n");
            fprintf(fid,"*\n");
            fprintf(fid,"* -------- Subckt for Hopping --------\n");
            fprintf(fid,"*\n");
            if options3.libmode 
                for i = 2: HcktObj.NRPTS
                    if HcktObj.ScktL(i).type == 'X'
                        fprintf(fid,"* %s\n",mat2str(HcktObj.SktL(i).name));
                        fprintf(fid,HcktObj.ScktL(i).OutputStr);
                        fprintf(fid,"\n");
                        fprintf(fid,"* --------");
                    end
                end
            end
            fprintf(fid,"*\n");
            % 
            % Ipulse
            switch options5.analysis
                case 'tran'
                    fprintf(fid,"* -------- Ipulse --------\n");
                    fprintf(fid,"%s\n",HcktObj.Ipulse); 
                case 'ac'
                    fprintf(fid,"* -------- V ac --------\n");
                    fprintf(fid,"%s\n",HcktObj.Vac);
                otherwise

            end
        end
        function HcktObj = hspice_gen_vectorized(HcktObj,filename,options,options2,options3)
            arguments
                HcktObj Hckt;
                filename = strcat("Hckt_test",".sp");
                options.mesh = repmat(10,[HcktObj.dim,1]);
                options.fin_dir = zeros(HcktObj.dim,1);
                options.BinBin = false;
                options2.probenode {mustBeMember(options2.probenode,{'Allnode','Selectnode'})}= 'Allnode';
                options2.prefix = "+ ";
                options2.father = "";
                options2.nodelist = [];
                options3.fft = false;
                options3.libmode = false;
            end
            %
            optionscell = namedargs2cell(options);
            %
            HcktObj = HcktObj.vectorize(options.fin_dir);
            fid = fopen(filename,'w');
            % title
            fprintf(fid,HcktObj.title);
            fprintf(fid,"\n");
            % options & lib
            HcktObj.optionslib_write(fid,'libmode',options3.libmode);
            % Vectorize,
            fprintf(fid,"* -------- Subckt for Vectorize --------\n");
            fprintf(fid,"*\n");
            fprintf(fid,"*%s\n",HcktObj.VectorizeScktL(1).name);
            fprintf(fid,HcktObj.VectorizeScktL(1).OutputStr);
            fprintf(fid,"\n");
            fprintf(fid,"*\n");
            fprintf(fid,"* -------- Subckt for Vectorize (OpenBoundary) --------\n");
            fprintf(fid,"*\n");
            for i = 2: length(HcktObj.VectorizeScktL)
                fprintf(fid,"* %s\n",HcktObj.VectorizeScktL(i).name);
                fprintf(fid,HcktObj.VectorizeScktL(i).OutputStr);
                fprintf(fid,"\n");
                fprintf(fid,"* --------");
            end
            fprintf(fid,"* -------- -------- --------\n");
            % Netlist
            fprintf(fid,"* -------- Netlist --------\n");
            fprintf(fid,"*\n");
            [HcktObj,ScktnameL,AllPortStrL] = HcktObj.netlist_gen(fid,optionscell{:},'convention','II');
            fprintf(fid,"*\n");
            fprintf(fid,"* -------- Netlist end --------\n");
            fprintf(fid,"*\n");
            % Probe List
            fprintf(fid,"* -------- Probe List --------\n");
            fprintf(fid,"*\n");
            fprintf(fid,".PROBE TRAN\n");
            optionsProbe = options2;
            optionsProbe.comment = false;
            optionsProbeCell = namedargs2cell(optionsProbe);
            switch options2.probenode
                case 'Allnode'
                    Hckt.Portlist_gen(fid,AllPortStrL,optionsProbeCell{:});
                case 'Selectnode'
                    Hckt.Portlist_gen(fid,ScktnameL,optionsProbeCell{:});
            end
            fprintf(fid,"* -------- Probe List end --------\n");
            fprintf(fid,"*\n");
            % End
            % fft List
            fprintf(fid,"* -------- FFT List --------\n");
            fprintf(fid,"*\n");
            optionsFFT = options2;
            optionsFFT.comment = ~options3.fft;
            optionsFFT.prefix = ".FFT ";
            optionsFFTcell = namedargs2cell(optionsFFT);
            switch options2.probenode
                case 'Allnode'
                    Hckt.Portlist_gen(fid,AllPortStrL,optionsFFTcell{:});
                case 'Selectnode'
                    Hckt.Portlist_gen(fid,ScktnameL,optionsFFTcell{:});
            end
            fprintf(fid,"* -------- FFT List end --------\n");
            fprintf(fid,"*\n");
            % End
            fprintf(fid,".END\n");
            fclose(fid);
        end
        function HcktObj = hspice_gen_general(HcktObj,filename,options,options2,options3,options5)
            arguments
                HcktObj Hckt;
                filename = strcat("Hckt_test",".sp");
                options.mesh = repmat(10,[HcktObj.dim,1]);
                options.fin_dir = zeros(HcktObj.dim,1);
                options2.probenode {mustBeMember(options2.probenode,{'Allnode','Selectnode'})}= 'Allnode';
                options2.prefix = "+ ";
                options2.father = "";
                options2.nodelist = [];
                options3.fft = false;
                options3.libmode = false;
                options5.analysis {mustBeMember(options5.analysis,{'tran','ac'})} = 'tran';
            end
            %
            optionscell = namedargs2cell(options);
            %
            fid = fopen(filename,'w');
            % title
            fprintf(fid,HcktObj.title);
            fprintf(fid,"\n");
            % options & lib
            HcktObj.optionslib_write(fid,"analysis",options5.analysis );
            % Netlist Pri
            fprintf(fid,"* -------- Netlist Homecell --------\n");
            fprintf(fid,"*\n");
            [HcktObj,ScktnameL,AllPortStrL,ijkL] = netlist_gen(HcktObj,fid,optionscell{:});
            fprintf(fid,"*\n");
            fprintf(fid,"* -------- Netlist Homecell end --------\n");
            fprintf(fid,"*\n");
            % Netlist Hopping
            fprintf(fid,"* -------- Netlist Hopping --------\n");
            fprintf(fid,"*\n");
            for i = 2:HcktObj.Components
                fprintf(fid,"* \\ %s|%s <-- %s \\ %s\n",...
                    mat2str(HcktObj.vectorAll(HcktObj.vectorL(i),:)),...
                    mat2str(HcktObj.PortOutCell{i,:}),...
                    mat2str(HcktObj.PortInCell{i,:}),...
                    HcktObj.ScktL(i).name(1) ...
                    );
                hoppinglist_gen(HcktObj,fid,i,ijkL,ScktnameL,optionscell{:});
            end
            fprintf(fid,"*\n");
            fprintf(fid,"* -------- Netlist Hopping end --------\n");
            fprintf(fid,"*\n");
            % Probe List
            switch options5.analysis
                case 'tran'
                    fprintf(fid,"* -------- Probe List --------\n");
                    fprintf(fid,"*\n");
                    fprintf(fid,".PROBE TRAN\n");
                    optionsProbe = options2;
                    optionsProbe.comment = false;
                    optionsProbeCell = namedargs2cell(optionsProbe);
                    switch options2.probenode
                        case 'Allnode'
                            Hckt.Portlist_gen(fid,AllPortStrL,optionsProbeCell{:});
                        case 'Selectnode'
                            Hckt.Portlist_gen(fid,ScktnameL,optionsProbeCell{:});
                    end
                    fprintf(fid,"* -------- Probe List end --------\n");
                    fprintf(fid,"*\n");
                    % End
                    % fft List
                    fprintf(fid,"* -------- FFT List --------\n");
                    fprintf(fid,"*\n");
                    optionsFFT = options2;
                    optionsFFT.comment = ~options3.fft;
                    optionsFFT.prefix = ".FFT ";
                    optionsFFTcell = namedargs2cell(optionsFFT);
                    switch options2.probenode
                        case 'Allnode'
                            Hckt.Portlist_gen(fid,AllPortStrL,optionsFFTcell{:});
                        case 'Selectnode'
                            Hckt.Portlist_gen(fid,ScktnameL,optionsFFTcell{:});
                    end
                    fprintf(fid,"* -------- FFT List end --------\n");
                    fprintf(fid,"*\n");
                case 'ac'
                    fprintf(fid,"* -------- Probe List --------\n");
                    fprintf(fid,"*\n");
                    fprintf(fid,".PROBE AC\n");
                    optionsProbe = options2;
                    optionsProbe.comment = false;
                    optionsProbe.type = 'VR';
                    optionsProbeCell = namedargs2cell(optionsProbe);
                    switch options2.probenode
                        case 'Allnode'
                            Hckt.Portlist_gen(fid,AllPortStrL,optionsProbeCell{:});
                        case 'Selectnode'
                            Hckt.Portlist_gen(fid,ScktnameL,optionsProbeCell{:});
                    end
                    optionsProbe.type = 'VI';
                    optionsProbeCell = namedargs2cell(optionsProbe);
                    switch options2.probenode
                        case 'Allnode'
                            Hckt.Portlist_gen(fid,AllPortStrL,optionsProbeCell{:});
                        case 'Selectnode'
                            Hckt.Portlist_gen(fid,ScktnameL,optionsProbeCell{:});
                    end
                    fprintf(fid,"* -------- Probe List end --------\n");
                    fprintf(fid,"*\n");
                    % End
            end
            % End
            fprintf(fid,".END\n");
            fclose(fid);
        end
        function [ScktnameL] = hoppinglist_gen(HcktObj,fid,seq,ijkL,ScktnameL,options)
            arguments
                HcktObj Hckt;
                fid = 1;
                seq = 2;
                ijkL = [];
                ScktnameL = [];
                options.BinBin = false;
                options.mesh = repmat(10,[HcktObj.dim,1]);
                options.fin_dir = zeros(HcktObj.dim,1);
                options.convention = 'I';
            end
            switch class(fid)
                case {'char','string'}
                    fid = fopen(fid,'w');
                case 'double'
            end
            if isempty(HcktObj.mesh)
                HcktObj.mesh = options.mesh ;
            end
            if length(HcktObj.mesh) ~= HcktObj.dim
                error('wrong input: mesh! The dim of HcktObj is %d',HcktObj.dim);
            end
            if isempty(HcktObj.fin_dir)
                HcktObj.fin_dir = options.fin_dir ;
            end
            if sum(HcktObj.fin_dir)
                fin_dir_mode=true;
            else
                fin_dir_mode=false;
            end
            %
            % first Periodic Boundary Conditions
            UsefulMesh = HcktObj.mesh(HcktObj.mesh>1);
            Ntotmesh = prod(UsefulMesh);
            % SnameL encoding
            ScktnameL = strrep(ScktnameL,'X',...
            [HcktObj.ScktL(seq).type,num2str(seq)]);
            underscore_L=repmat('_',[Ntotmesh,1]);
            nodepreL = repmat('n',[Ntotmesh,1]);
            NPortInCell = length(HcktObj.PortInCell{seq});
            NPortOutCell = length(HcktObj.PortOutCell{seq});
            NPort = NPortInCell+NPortOutCell;
            portL{NPort} = ijkL;
            GeneralRvector =  HcktObj.vectorAll(HcktObj.vectorL(seq),:);
            %AllPortStrL = [];
            if strcmp(options.convention,'I')
                % firstly InCell Consider InCell belong them thelf!
                for i = 1:NPortInCell
                    ConnectionPort{i} = HcktObj.PortInCell{seq}(i);% conect to myself
                    % Be alerted for node1 conect to R[,,,,]node2
                    portL{i} = ijkL;
                    % nevercheck exceed l
                    %portL{i}(portL{i}>HcktObj.mesh) = 1;
                    % nevercheck exceed r
                    %for j = 1:HcktObj.dim
                    %    portL{i}((portL{i}(:,j) ==0),j)= HcktObj.mesh(j);
                    %end
                    %
                    %AllPortStrL = [AllPortStrL,string(tmpPortStrL)];
                end
                if fin_dir_mode
                    % prepare adding node
                    % Duplicate portL
                    DupliteportL = portL;
                    DupliteportLabel = false(length(ScktnameL),1);
                    exceed_logical{NPort} = false(length(ScktnameL),1);
                end
                % Secondly InCell Consider InCell belong other selves!
                for i = (NPortInCell+1):NPort
                    ConnectionPort{i} = HcktObj.PortOutCell{seq}(i-NPortInCell);% conect to R[???] cell
                    portL{i} = ijkL + GeneralRvector;
                    if fin_dir_mode
                        DupliteportL{i} = portL{i};
                    end
                    exceed_logical{i} = false(length(ScktnameL),1);
                    % exceed r
                    for j = 1:HcktObj.dim
                        exceed_r_logical = portL{i}(:,j)>HcktObj.mesh(j);
                        exceed_r_label = find(exceed_r_logical);
                        portL{i}(exceed_r_label,j) = mod(portL{i}(exceed_r_label,j),HcktObj.mesh(j));
                        if fin_dir_mode
                            if HcktObj.fin_dir(j)
                                % add another
                                exceed_logical{i} = exceed_logical{i} + exceed_r_logical;
                                DupliteportL{i}(exceed_r_label,j) = mod(portL{i}(exceed_r_label,j),HcktObj.mesh(j));
                                DupliteportL{i-NPortInCell}(exceed_r_label,j) = 0;
                            else
                                DupliteportL{i}(exceed_r_label,j) = portL{i}(exceed_r_label,j);
                            end
                        end
                    end
                    % exceed l
                    for j = 1:HcktObj.dim
                        exceed_l_logical = (portL{i}(:,j)-1) < 0 ;
                        exceed_l_label = find(exceed_l_logical);
                        portL{i}(exceed_l_label,j)= mod(portL{i}(exceed_l_label,j)-1,HcktObj.mesh(j));
                        if fin_dir_mode
                            if HcktObj.fin_dir(j)
                                % add another
                                exceed_logical{i} = exceed_logical{i} + exceed_l_logical;
                                DupliteportL{i}(exceed_l_label,j) = mod(portL{i}(exceed_l_label,j)-1,HcktObj.mesh(j));
                                DupliteportL{i-NPortInCell}(exceed_l_label,j) = 0;
                            else
                                DupliteportL{i}(exceed_l_label,j) = portL{i}(exceed_l_label,j);
                            end
                        end
                    end
                    if fin_dir_mode
                        % if exceed
                        portL{i}(logical(exceed_logical{i}) ,j) = 0 ;
                        DupliteportLabel = DupliteportLabel + exceed_logical{i};
                    end
                    %AllPortStrL = [AllPortStrL,string(tmpPortStrL)];
                end
                for i = 1:NPort                 %
                    tmpPortStrL= nodepreL;
                    for j = 1:HcktObj.dim
                        tmpcharL = park.num2strwithzero(portL{i}(:,j));
                        tmpPortStrL =[tmpPortStrL,underscore_L,tmpcharL];
                    end
                    tmpPortStrL = [tmpPortStrL,underscore_L,park.num2strwithzero(ConnectionPort{i}*ones(Ntotmesh,1))];
                    PortStrL{i}  = string(tmpPortStrL);
                    if i > NPortInCell
                        PortStrL{i}(logical(exceed_logical{i})) = "GND";
                    else

                    end
                end
                % addtional
                if fin_dir_mode
                    % prepare adding components
                    DupliteportLabel = logical(DupliteportLabel);
                    AddtionalScktnameL = ScktnameL(DupliteportLabel,:);
                    Naddtional = size(AddtionalScktnameL,1);
                    AddtionalScktnameL = strcat(AddtionalScktnameL,repmat('_extra',[Naddtional,1]));
                    % prepare adding node
                    AddtionalnodepreL = repmat('n',[Naddtional,1]);
                    Addtionalunderscore_L =repmat('_',[Naddtional,1]);
                    for i = 1:NPort
                        tmpPortStrL= AddtionalnodepreL;
                        for j = 1:HcktObj.dim
                            tmpcharL = park.num2strwithzero([DupliteportL{i}(DupliteportLabel,j);portL{i}(:,j)]);
                            tmpcharL = tmpcharL(1:Naddtional,:);
                            tmpPortStrL =[tmpPortStrL,Addtionalunderscore_L,tmpcharL];
                        end
                        ConnectionPortStrL = park.num2strwithzero(ConnectionPort{i}*ones(Naddtional,1));
                        tmpPortStrL = [tmpPortStrL,Addtionalunderscore_L,ConnectionPortStrL];
                        AddtionalPortStrL{i}  = string(tmpPortStrL);
                        if i <= NPortInCell
                            AddtionalPortStrL{i}(1:Naddtional) = "GND";
                        else

                        end
                    end
                end
            else
                % check reverseconnection!!!!
            end
            %
            %AllPortStrL = unique(AllPortStrL);
            % Suppose every ScktL has GND
            %GND_nameL = repmat('GND',[Ntotmesh,1]);
            % NameL
            if strcmp(HcktObj.ScktL(seq).type,'X')
                basePriVname = strcat("GND ",HcktObj.ScktL(seq).name);
            else
                basePriVname = strcat(HcktObj.ScktL(seq).name,HcktObj.ScktL(seq).description);
            end
            % basePriVname = ['PriV_',num2str(ones(1,HcktObj.dim))];
            PriVnameL = repmat(basePriVname,[Ntotmesh,1]);% To be extended
            if fin_dir_mode
                AddtionalPriVnameL = repmat(basePriVname,[Naddtional,1]);% To be extended
            end
            % if need, save this variable
            % Gen NetlistMesh
            for i = 1:Ntotmesh
                % Scktcode
                fprintf(fid,"%s ",ScktnameL(i,:));
                % node
                for n = 1:NPort
                    fprintf(fid,"%s ",PortStrL{n}(i,:));
                end
                fprintf(fid,"%s %s\n",PriVnameL(i,:),HcktObj.DescriptionL{seq});
            end
            if sum(HcktObj.fin_dir)
                % Gen Addtional NetlistMesh
                for i = 1:Naddtional
                    % Scktcode
                    fprintf(fid,"%s ",AddtionalScktnameL(i,:));
                    % node
                    for n = 1:NPort
                        fprintf(fid,"%s ",AddtionalPortStrL{n}(i,:));
                    end
                    fprintf(fid,"%s %s\n",AddtionalPriVnameL(i,:),HcktObj.DescriptionL{seq});
                end
            end
        end
        function [HcktObj,ScktnameL,AllPortStrL,ijkL] = netlist_gen(HcktObj,fid,options)
            arguments
                HcktObj Hckt;
                fid = 1;
                options.BinBin = false;
                options.mesh = repmat(10,[HcktObj.dim,1]);
                options.fin_dir = zeros(HcktObj.dim,1);
                options.convention = 'I';
            end
            switch class(fid)
                case {'char','string'}
                    fid = fopen(fid,'w');
                case 'double'
            end
            if isempty(HcktObj.mesh)
                HcktObj.mesh = options.mesh ;
            end
            if length(HcktObj.mesh) ~= HcktObj.dim
                error('wrong input: mesh! The dim of HcktObj is %d',HcktObj.dim);
            end
            if isempty(HcktObj.fin_dir)
                HcktObj.fin_dir = options.fin_dir ;
            end
            %
            for t = 1:HcktObj.dim
                dimMesh{t}  = (1:HcktObj.mesh(t)).';
            end
            if options.BinBin 
                PPdist = HcktObj.Port2PortDist_BinBin;
            else
                PPdist = HcktObj.Port2PortDist;
            end
            % refer
            % https://ww2.mathworks.cn/help/matlab/ref/bsxfun.html
            %
            % first Periodic Boundary Conditions
            UsefulMesh = HcktObj.mesh(HcktObj.mesh>1);
            Ntotmesh = prod(UsefulMesh);
            ijkL = zeros(Ntotmesh,HcktObj.dim);
            count = 0;
            for t = find(HcktObj.mesh>1)
                count = count +1;
                if count >1
                    oldMesh = prod(UsefulMesh(1:count-1));
                else
                    oldMesh = 1;
                end
                ijkL(:,t) = repmat(kron(dimMesh{t},ones(oldMesh,1)),...
                    [Ntotmesh/(HcktObj.mesh(t)*oldMesh),1]);
            end
            % work!
            % SnameL encoding
            ScktnameL = repmat('X',[Ntotmesh,1]);
            underscore_L=repmat('_',[Ntotmesh,1]);
            nodepreL = repmat('n',[Ntotmesh,1]);
            for i = 1:size(ijkL,2)
                tmpcharL = park.num2strwithzero(ijkL(:,i));
                ScktnameL = [ScktnameL,underscore_L,tmpcharL];
            end
            ScktnameL = string(ScktnameL);
            % when print use 80 rules!;
            portL{HcktObj.Nports} = ijkL;
            AllPortStrL = [];
            if strcmp(options.convention,'II')
                for i = 1:HcktObj.Nports
                    ConnectionPort = PPdist(i);% conect 2 myself
                    % Be alerted for node1 conect to R[,,,,]node2

                    portL{i} = portL{HcktObj.Nports} + HcktObj.Port2VectorDist(i);
                    % exceed l
                    portL{i}(portL{i}>HcktObj.mesh) = 1;
                    % exceed r
                    for j = 1:HcktObj.dim
                        portL{i}((portL{i}(:,j) ==0),j)= HcktObj.mesh(j);
                    end
                    %
                    tmpPortStrL= nodepreL;
                    for j = 1:HcktObj.dim
                        tmpcharL = park.num2strwithzero(portL{i}(:,j));
                        tmpPortStrL =[tmpPortStrL,underscore_L,tmpcharL];
                    end
                    tmpPortStrL = [tmpPortStrL,underscore_L,park.num2strwithzero(ConnectionPort*ones(Ntotmesh,1))];
                    PortStrL{i}  = tmpPortStrL;
                    AllPortStrL = [AllPortStrL,string(tmpPortStrL)];
                end
            else
                for i = 1:HcktObj.Nports
                    ConnectionPort = i;% conect to myself
                    % Be alerted for node1 conect to R[,,,,]node2
                    portL{i} = ijkL;
                    % nevercheck exceed l
                    %portL{i}(portL{i}>HcktObj.mesh) = 1;
                    % nevercheck exceed r
                    %for j = 1:HcktObj.dim
                    %    portL{i}((portL{i}(:,j) ==0),j)= HcktObj.mesh(j);
                    %end
                    %
                    tmpPortStrL= nodepreL;
                    for j = 1:HcktObj.dim
                        tmpcharL = park.num2strwithzero(portL{i}(:,j));
                        tmpPortStrL =[tmpPortStrL,underscore_L,tmpcharL];
                    end
                    tmpPortStrL = [tmpPortStrL,underscore_L,park.num2strwithzero(ConnectionPort*ones(Ntotmesh,1))];
                    PortStrL{i}  = tmpPortStrL;
                    AllPortStrL = [AllPortStrL,string(tmpPortStrL)];
                end
            end
            %
            AllPortStrL = unique(AllPortStrL);
            % Suppose every Homecell has GND
            %GND_nameL = repmat('GND',[Ntotmesh,1]);
            % NameL
            if strcmp(options.convention,'II')
                basePriVname = 'PriV';
                for i = 1:HcktObj.dim
                    basePriVname = [basePriVname,'_1'];
                end
            else
                basePriVname = 'Pri';
            end
            % basePriVname = ['PriV_',num2str(ones(1,HcktObj.dim))];
            PriVnameL = repmat(basePriVname,[Ntotmesh,1]);
            % if need, save this variable
            % Gen NetlistMesh
            for i = 1:Ntotmesh
                % Scktcode
                fprintf(fid,"%s ",ScktnameL(i,:));
                % node
                for n = 1:HcktObj.Nports
                    fprintf(fid,"%s ",PortStrL{n}(i,:));
                end
                fprintf(fid,"%s %s %s\n",'GND',PriVnameL(i,:),HcktObj.AddtionInformation);
            end
        end
    end
    methods(Static)
        function WriteVendorComponentLibrary

        end
        function WritePlugIns(PlugIns,fid,magnitude)
            arguments
                PlugIns char;
                fid;
                magnitude = 'p';
            end
            switch magnitude
                case 'p'
                    Cmagnitude = 'p';
                    Lmagnitude = 'u';
                    Rmagnitude = '';
                case 'u'
                    Cmagnitude = 'p';
                    Lmagnitude = 'u';
                    Rmagnitude = '';
                case 'm'
                    Cmagnitude = 'p';
                    Lmagnitude = 'u';
                    Rmagnitude = '';
            end
            fprintf(fid,"* ---------------\n");
            switch PlugIns
                case 'Na-(Mb+Lc)' % checked
                    fprintf(fid,"* Na-(Mb+Lc) \n");
                    fprintf(fid,"*\n");
                    fprintf(fid,"*Rb = Rf/M; Rc = Rf/L Ra = (2+M+L-N)Rf RNf = N Rf \n");
                    fprintf(fid,".SubCkt AdderSubtractor va vb vc vo TOGND " + ...
                        "VarRf=100%s VarRa=100%s VarRb=100%s VarRc=100%s VarRNf=300%s\n", ...
                        Rmagnitude,Rmagnitude,Rmagnitude,Rmagnitude,Rmagnitude);
                    fprintf(fid,"Ra va v_plus VarRa \n");
                    fprintf(fid,"Rb vb v_minus VarRb \n");
                    fprintf(fid,"Rc vc v_minus VarRc \n");
                    fprintf(fid,"Rf1 TOGND v_minus VarRf \n");
                    fprintf(fid,"Rf2 vo v_minus VarRf \n");
                    fprintf(fid,"RNf TOGND v_plus VarRNf \n");
                    fprintf(fid,"E_opamp vo TOGND v_plus v_minus  level=1\n");
                    fprintf(fid,".ends AdderSubtractor\n");
                case 'VoltageFollower'
                    fprintf(fid,"* VoltageFollower \n");
                    fprintf(fid,"*\n");
                    fprintf(fid,".SubCkt VoltageFollower ui uo TOGND\n");
                    fprintf(fid,"E_opamp uo TOGND ui uo  level=1\n");
                    fprintf(fid,".ends VoltageFollower\n");
            end
        end
        function WriteModules(Modules,fid,magnitude)
            arguments
                Modules char;
                fid;
                magnitude = 'p';
            end
            switch magnitude
                case 'p'
                    Cmagnitude = 'p';
                    Lmagnitude = 'u';
                    Rmagnitude = 'k';
                case 'u'
                    Cmagnitude = 'p';
                    Lmagnitude = 'u';
                    Rmagnitude = 'k';
                case 'm'
                    Cmagnitude = 'p';
                    Lmagnitude = 'u';
                    Rmagnitude = 'k';
            end
            fprintf(fid,"* ---------------\n");
            switch Modules
                case 'Basis' % checked
                    fprintf(fid,"* BasisC3_origin \n");
                    fprintf(fid,"*\n");
                    fprintf(fid,".SubCkt BasisC3_origin n1 n2 n3 TOGND VarL0=1%s InitV=0V R_L=1u \n",Lmagnitude);
                    fprintf(fid,"L1 n1 n2 VarL0 R=R_L \n");
                    fprintf(fid,"L2 n2 n3 VarL0 R=R_L \n");
                    fprintf(fid,"L3 n3 n2 VarL0 R=R_L \n");
                    fprintf(fid,".ends BasisC3_origin\n");
                case 'Basis_SOI' % checked
                    fprintf(fid,"* BasisC3_origin \n");
                    fprintf(fid,"*\n");
                    fprintf(fid,"* n1 n2 n3 n1_prime n2_prime n3_prime nR1_prime nR2_prime nR3_prime  \n");
                    fprintf(fid,".SubCkt BasisC3_origin n1 n2 n3 TOGND " + ...
                        "VarL0=1%s VarCg=100%s VarRh=0.1%s R_L=1u InitV=0V \n",Lmagnitude,Cmagnitude,Rmagnitude);
                    fprintf(fid,"L1 n1 n2 VarL0 R=R_L \n");
                    fprintf(fid,"L2 n2 n3 VarL0 R=R_L \n");
                    fprintf(fid,"L3 n3 n2 VarL0 R=R_L \n");
                    fprintf(fid,"X1_1prime n1 n1_prime TOGND VoltageFollower \n");
                    fprintf(fid,"X2_2prime n2 n2_prime TOGND VoltageFollower \n");
                    fprintf(fid,"X3_3prime n3 n3_prime TOGND VoltageFollower \n");
                    fprintf(fid,"X1prime n2_prime n1_prime n3_prime NR1_prime TOGND AdderSubtractor \n");
                    fprintf(fid,"X2prime n3_prime n2_prime n1_prime NR2_prime TOGND AdderSubtractor \n");
                    fprintf(fid,"X3prime n1_prime n3_prime n2_prime NR3_prime TOGND AdderSubtractor \n");
                    fprintf(fid,"R1prime n2 NR1_prime VarRh \n");
                    fprintf(fid,"R2prime n3 NR2_prime VarRh \n");
                    fprintf(fid,"R3prime n1 NR3_prime VarRh \n");
                    fprintf(fid,"C1 n1 TOGND  VarCg \n");
                    fprintf(fid,"C2 n2 TOGND  VarCg \n");
                    fprintf(fid,"C3 n3 TOGND  VarCg \n");
                    fprintf(fid,".ends BasisC3_origin\n");
                case '+sigma_0' % checked
                    fprintf(fid,"* +sigma_0 \n");
                    fprintf(fid,"*\n");
                    fprintf(fid,".SubCkt PlusSigma0 L_n1 L_n2 L_n3 R_n1 R_n2 R_n3 TOGND InitV=0V " + ...
                        "VarC0=100%s Var2C0=200%s VarR0=100%s VarR0_2=50%s\n",Cmagnitude,Cmagnitude,Rmagnitude,Rmagnitude);
                    fprintf(fid,"C1 L_n1 R_n1  VarC0 IC=InitV\n");
                    fprintf(fid,"C2 L_n2 R_n2 VarC0 IC=InitV\n");
                    fprintf(fid,"C3 L_n3 R_n3 VarC0 IC=InitV\n");
                    fprintf(fid,".ends PlusSigma0\n");
                case '-sigma_0' % checked
                    fprintf(fid,"* -sigma_0 \n");
                    fprintf(fid,"*\n");
                    fprintf(fid,".SubCkt MinusSigma0 L_n1 L_n2 L_n3 R_n1 R_n2 R_n3 TOGND InitV=0V " + ...
                        "VarC0=100%s Var2C0=200%s VarR0=100%s VarR0_2=50%s\n",Cmagnitude,Cmagnitude,Rmagnitude,Rmagnitude);
                    fprintf(fid,"C1 L_n1  R_n2 VarC0 IC=InitV\n");
                    fprintf(fid,"C2 L_n2 R_n1  VarC0 IC=InitV\n");
                    fprintf(fid,"C3 L_n3 R_n1  VarC0 IC=InitV\n");
                    fprintf(fid,"C4 L_n1  R_n3 VarC0 IC=InitV\n");
                    fprintf(fid,"C5 L_n2 R_n3 VarC0 IC=InitV\n");
                    fprintf(fid,"C6 L_n3 R_n2 VarC0 IC=InitV\n");
                    fprintf(fid,".ends MinusSigma0\n");
                case '+isigma_0' % checked
                    fprintf(fid,"* +isigma_0 \n");
                    fprintf(fid,"*\n");
                    fprintf(fid,".SubCkt PlusiSigma0 L_n1 L_n2 L_n3 R_n1 R_n2 R_n3 TOGND InitV=0V " + ...
                        "VarC0=100%s Var2C0=200%s VarR0=100%s VarR0_2=50%s\n",Cmagnitude,Cmagnitude,Rmagnitude,Rmagnitude);
                    fprintf(fid,"R1 L_n1  R_n2 VarR0\n");
                    fprintf(fid,"R2 L_n2 R_n1  VarR0\n");
                    fprintf(fid,"R3 L_n3 R_n1  VarR0\n");
                    fprintf(fid,"R4 L_n1  R_n3 VarR0\n");
                    fprintf(fid,"R5 L_n2 R_n3 VarR0\n");
                    fprintf(fid,"R6 L_n3 R_n2 VarR0\n");
                    fprintf(fid,".ends PlusiSigma0\n");
                case '-isigma_0' % checked
                    fprintf(fid,"* -isigma_0 \n");
                    fprintf(fid,"*\n");
                    fprintf(fid,".SubCkt MinusiSigma0 L_n1 L_n2 L_n3 R_n1 R_n2 R_n3 TOGND InitV=0V " + ...
                        "VarC0=100%s Var2C0=200%s VarR0=100%s VarR0_2=50%s\n",Cmagnitude,Cmagnitude,Rmagnitude,Rmagnitude);
                    fprintf(fid,"R1 L_n1  R_n1  VarR0\n");
                    fprintf(fid,"R2 L_n2 R_n2 VarR0\n");
                    fprintf(fid,"R3 L_n3 R_n3 VarR0\n");
                    fprintf(fid,".ends MinusiSigma0\n");
                case '+sigma_1' % checked
                    fprintf(fid,"* +sigma_1 \n");
                    fprintf(fid,"*\n");
                    fprintf(fid,".SubCkt PlusSigma1 L_n1 L_n2 L_n3 R_n1 R_n2 R_n3 TOGND InitV=0V " + ...
                        "VarC0=100%s Var2C0=200%s VarR0=100%s VarR0_2=50%s\n",Cmagnitude,Cmagnitude,Rmagnitude,Rmagnitude);
                    fprintf(fid,"C1 L_n1  R_n2 VarC0 IC=InitV\n");
                    fprintf(fid,"C2 L_n2 R_n1  VarC0 IC=InitV\n");
                    fprintf(fid,"C3 L_n3 R_n3 VarC0 IC=InitV\n");
                    fprintf(fid,".ends PlusSigma1\n");
                case '-sigma_1' % checked
                    fprintf(fid,"* -sigma_1 \n");
                    fprintf(fid,"*\n");
                    fprintf(fid,".SubCkt MinusSigma1 L_n1 L_n2 L_n3 R_n1 R_n2 R_n3 TOGND InitV=0V " + ...
                        "VarC0=100%s Var2C0=200%s VarR0=100%s VarR0_2=50%s\n",Cmagnitude,Cmagnitude,Rmagnitude,Rmagnitude);
                    fprintf(fid,"C1 L_n1  R_n1  VarC0 IC=InitV\n");
                    fprintf(fid,"C2 L_n2 R_n2 VarC0 IC=InitV\n");
                    fprintf(fid,"C3 L_n3 R_n1  VarC0 IC=InitV\n");
                    fprintf(fid,"C4 L_n1  R_n3 VarC0 IC=InitV\n");
                    fprintf(fid,"C5 L_n2 R_n3 VarC0 IC=InitV\n");
                    fprintf(fid,"C6 L_n3 R_n2  VarC0 IC=InitV\n");
                    fprintf(fid,".ends MinusSigam1\n");
                case '+isigma_1' % checked
                    fprintf(fid,"* +isigma_1 \n");
                    fprintf(fid,"*\n");
                    fprintf(fid,".SubCkt PlusiSigma1 L_n1 L_n2 L_n3 R_n1 R_n2 R_n3 TOGND InitV=0V " + ...
                        "VarC0=100%s Var2C0=200%s VarR0=100%s VarR0_2=50%s\n",Cmagnitude,Cmagnitude,Rmagnitude,Rmagnitude);
                    fprintf(fid,"R1 L_n1  R_n1  VarR0\n");
                    fprintf(fid,"R2 L_n2 R_n2 VarR0\n");
                    fprintf(fid,"R3 L_n3 R_n1  VarR0\n");
                    fprintf(fid,"R4 L_n1  R_n3 VarR0\n");
                    fprintf(fid,"R5 L_n2 R_n3 VarR0\n");
                    fprintf(fid,"R6 L_n3 R_n2 VarR0\n");
                    fprintf(fid,".ends PlusiSigma1\n");
                case '-isigma_1' % checked
                    fprintf(fid,"* -isigma_1 \n");
                    fprintf(fid,"*\n");
                    fprintf(fid,"SubCkt MinusiSigma1 L_n1 L_n2 L_n3 R_n1 R_n2 R_n3 TOGND InitV=0V " + ...
                        "VarC0=100%s Var2C0=200%s VarR0=100%s VarR0_2=50%s\n",Cmagnitude,Cmagnitude,Rmagnitude,Rmagnitude);
                    fprintf(fid,"R1 L_n1  R_n2 VarR0\n");
                    fprintf(fid,"R2 L_n2 R_n1  VarR0\n");
                    fprintf(fid,"R3 L_n3 R_n3 VarR0\n");
                    fprintf(fid,".ends MinusiSigma1\n");
                case '+gen3sigma_2' % checked
                    fprintf(fid,"* +gen3sigma_2 \n");
                    fprintf(fid,"*\n");
                    fprintf(fid,".SubCkt PlusGen3Sigma2 L_n1 L_n2 L_n3 R_n1 R_n2 R_n3 TOGND InitV=0V " + ...
                        "VarC0=100%s Var2C0=200%s VarR0=100%s VarR0_2=50%s\n",Cmagnitude,Cmagnitude,Rmagnitude,Rmagnitude);
                    fprintf(fid,"C1 L_n1  R_n2 VarC0  IC=InitV\n");
                    fprintf(fid,"C2 L_n2 R_n1  VarC0  IC=InitV\n");
                    fprintf(fid,"C3 L_n3 R_n3 VarC0  IC=InitV\n");
                    fprintf(fid,"C4 L_n1  R_n3 Var2C0 IC=InitV\n");
                    fprintf(fid,"C5 L_n2 R_n2 Var2C0 IC=InitV\n");
                    fprintf(fid,"C6 L_n3 R_n1  Var2C0 IC=InitV\n");
                    fprintf(fid,".ends PlusGen3Sigma2\n");
                case '-gen3sigma_2' % checked
                    fprintf(fid,"* -gen3sigma_2 \n");
                    fprintf(fid,"*\n");
                    fprintf(fid,".SubCkt MinusGen3Sigma2 L_n1 L_n2 L_n3 R_n1 R_n2 R_n3 TOGND InitV=0V " + ...
                        "VarC0=100%s Var2C0=200%s VarR0=100%s VarR0_2=50%s\n",Cmagnitude,Cmagnitude,Rmagnitude,Rmagnitude);
                    fprintf(fid,"C1 L_n1  R_n2 VarC0  IC=InitV\n");
                    fprintf(fid,"C2 L_n2 R_n1  VarC0  IC=InitV\n");
                    fprintf(fid,"C3 L_n3 R_n3 VarC0  IC=InitV\n");
                    fprintf(fid,"C4 L_n1  R_n1  Var2C0 IC=InitV\n");
                    fprintf(fid,"C5 L_n2 R_n3 Var2C0 IC=InitV\n");
                    fprintf(fid,"C6 L_n3 R_n2 Var2C0 IC=InitV\n");
                    fprintf(fid,".ends MinusGen3Sigma2\n");
                case '+igen3sigma_2' % checked
                    fprintf(fid,"* +igen3sigma_2 \n");
                    fprintf(fid,"*\n");
                    fprintf(fid,".SubCkt PlusiGen3Sigma2 L_n1 L_n2 L_n3 R_n1 R_n2 R_n3 TOGND InitV=0V " + ...
                        "VarC0=100%s Var2C0=200%s VarR0=100%s VarR0_2=50%s\n",Cmagnitude,Cmagnitude,Rmagnitude,Rmagnitude);
                    fprintf(fid,"R1 L_n1 R_n2  VarR0  \n");
                    fprintf(fid,"R2 L_n2 R_n1  VarR0  \n");
                    fprintf(fid,"R3 L_n3 R_n3 VarR0  \n");
                    fprintf(fid,"R4 L_n1 R_n1   VarR0_2\n");
                    fprintf(fid,"R5 L_n2 R_n3 VarR0_2\n");
                    fprintf(fid,"R6 L_n3 R_n2 VarR0_2\n");
                    fprintf(fid,".ends PlusiGen3Sigma2\n");
                case '-igen3sigma_2'
                    fprintf(fid,"* -igen3sigma_2 \n");
                    fprintf(fid,"*\n");
                    fprintf(fid,".SubCkt MinusiGen3Sigma2 L_n1 L_n2 L_n3 R_n1 R_n2 R_n3 TOGND InitV=0V " + ...
                        "VarC0=100%s Var2C0=100%s   VarR0=100%s VarR0_2=100%s\n",Cmagnitude,Cmagnitude,Rmagnitude,Rmagnitude);
                    fprintf(fid,"R1 L_n1 R_n2  VarR0  \n");
                    fprintf(fid,"R2 L_n2 R_n1  VarR0  \n");
                    fprintf(fid,"R3 L_n3 R_n3 VarR0  \n");
                    fprintf(fid,"R4 L_n1 R_n3  VarR0_2\n");
                    fprintf(fid,"R5 L_n2 R_n2 VarR0_2\n");
                    fprintf(fid,"R6 L_n3 R_n1  VarR0_2\n");
                    fprintf(fid,".ends MinusiGen3Sigma2\n");
                case '+gen3sigma_3'
                    fprintf(fid,"* +gen3sigma_3 \n");
                    fprintf(fid,"*\n");
                    fprintf(fid,".SubCkt PlusGen3Sigma3 L_n1 L_n2 L_n3 R_n1 R_n2 R_n3 TOGND InitV=0V " + ...
                        "VarC0=100%s Var2C0=200%s VarR0=100%s VarR0_2=50%s\n",Cmagnitude,Cmagnitude,Rmagnitude,Rmagnitude);
                    fprintf(fid,"R1 L_n1  R_n1  VarR0  \n");
                    fprintf(fid,"R2 L_n2 R_n2 VarR0  \n");
                    fprintf(fid,"R3 L_n3 R_n3 VarR0  \n");
                    fprintf(fid,"R4 L_n1  R_n2 VarR0_2\n");
                    fprintf(fid,"R5 L_n2 R_n3 VarR0_2\n");
                    fprintf(fid,"R6 L_n3 R_n1  VarR0_2\n");
                    fprintf(fid,".ends PlusGen3Sigma3\n");
                case '-gen3sigma_3'
                    fprintf(fid,"* -gen3sigma_3 \n");
                    fprintf(fid,"*\n");
                    fprintf(fid,".SubCkt MinusGen3Sigma3 L_n1 L_n2 L_n3 R_n1 R_n2 R_n3 TOGND InitV=0V " + ...
                        "VarC0=100%s Var2C0=200%s VarR0=100%s VarR0_2=50%s\n",Cmagnitude,Cmagnitude,Rmagnitude,Rmagnitude);
                    fprintf(fid,"R1 L_n1 R_n1   VarR0  \n");
                    fprintf(fid,"R2 L_n2 R_n2 VarR0  \n");
                    fprintf(fid,"R3 L_n3 R_n3 VarR0  \n");
                    fprintf(fid,"R4 L_n1 R_n3  VarR0_2\n");
                    fprintf(fid,"R5 L_n2 R_n1  VarR0_2\n");
                    fprintf(fid,"R6 L_n3 R_n2 VarR0_2\n");
                    fprintf(fid,".ends MinusGen3Sigma3\n");
                case '+igen3sigma_3'
                    fprintf(fid,"* +igen3sigma_3 \n");
                    fprintf(fid,"*\n");
                    fprintf(fid,".SubCkt PlusiGen3Sigma3 L_n1 L_n2 L_n3 R_n1 R_n2 R_n3 TOGND InitV=0V " + ...
                        "VarC0=100%s Var2C0=200%s VarR0=100%s VarR0_2=50%s\n",Cmagnitude,Cmagnitude,Rmagnitude,Rmagnitude);
                    fprintf(fid,"C1 L_n1 R_n1 VarC0 \n");
                    fprintf(fid,"C2 L_n2 R_n2 VarC0 \n");
                    fprintf(fid,"C3 L_n3 R_n3 VarC0 \n");
                    fprintf(fid,"C4 L_n1 R_n2 Var2C0\n");
                    fprintf(fid,"C5 L_n2 R_n3 Var2C0\n");
                    fprintf(fid,"C6 L_n3 R_n1 Var2C0\n");
                    fprintf(fid,".ends PlusiGen3Sigma3\n");
                case '-igen3sigma_3'
                    fprintf(fid,"* -igen3sigma_3 \n");
                    fprintf(fid,"*\n");
                    fprintf(fid,".SubCkt MinusiGen3Sigma3 L_n1 L_n2 L_n3 R_n1 R_n2 R_n3 TOGND InitV=0V " + ...
                        "VarC0=100%s Var2C0=200%s VarR0=100%s VarR0_2=50%s\n",Cmagnitude,Cmagnitude,Rmagnitude,Rmagnitude);
                    fprintf(fid,"C1 L_n1 R_n1 VarC0 \n");
                    fprintf(fid,"C2 L_n2 R_n2 VarC0 \n");
                    fprintf(fid,"C3 L_n3 R_n3 VarC0 \n");
                    fprintf(fid,"C4 L_n1 R_n3 Var2C0\n");
                    fprintf(fid,"C5 L_n2 R_n1 Var2C0\n");
                    fprintf(fid,"C6 L_n3 R_n2 Var2C0\n");
                    fprintf(fid,".ends MinusiGen3Sigma3\n");
                case '+isigma_1_SOI' % checked
                    fprintf(fid,"* +isigma_1_SOI \n");
                    fprintf(fid,"*\n");
                    fprintf(fid,".SubCkt PlusiSigma1 L_n1_prime L_n2_prime L_n3_prime R_n1 R_n2 R_n3 TOGND InitV=0V " + ...
                        "VarC0=100%s Var2C0=200%s VarR0=100%s VarR0_2=50%s\n",Cmagnitude,Cmagnitude,Rmagnitude,Rmagnitude);
                    fprintf(fid,"R1 L_n1 R_n1  VarR0\n");
                    fprintf(fid,"R2 L_n2 R_n2 VarR0\n");
                    fprintf(fid,"R3 L_n3 R_n1  VarR0\n");
                    fprintf(fid,"R4 L_n1 R_n3 VarR0\n");
                    fprintf(fid,"R5 L_n2 R_n3 VarR0\n");
                    fprintf(fid,"R6 L_n3 R_n2 VarR0\n");
                    fprintf(fid,"X1_1prime L_n1_prime L_n1 TOGND VoltageFollower \n");
                    fprintf(fid,"X2_2prime L_n2_prime L_n2 TOGND VoltageFollower \n");
                    fprintf(fid,"X3_3prime L_n3_prime L_n3 TOGND VoltageFollower \n");
                    fprintf(fid,".ends PlusiSigma1\n");
                case '-isigma_1_SOI' % checked
                    fprintf(fid,"* -isigma_1_SOI \n");
                    fprintf(fid,"*\n");
                    fprintf(fid,".SubCkt MinusiSigma1 L_n1 L_n2 L_n3 R_n1_prime R_n2_prime R_n3_prime TOGND InitV=0V " + ...
                        "VarC0=100%s Var2C0=200%s VarR0=100%s VarR0_2=50%s\n",Cmagnitude,Cmagnitude,Rmagnitude,Rmagnitude);
                    fprintf(fid,"R1 L_n1 R_n2 VarR0\n");
                    fprintf(fid,"R2 L_n2 R_n1 VarR0\n");
                    fprintf(fid,"R3 L_n3 R_n3 VarR0\n");
                    fprintf(fid,"X1_1prime R_n1_prime R_n1 TOGND VoltageFollower \n");
                    fprintf(fid,"X2_2prime R_n2_prime R_n2 TOGND VoltageFollower \n");
                    fprintf(fid,"X3_3prime R_n3_prime R_n3 TOGND VoltageFollower \n");
                    fprintf(fid,".ends MinusiSigma1\n");
                case '+igen3sigma_2_SOI' % checked
                    fprintf(fid,"* +igen3sigma_2_SOI \n");
                    fprintf(fid,"*\n");
                    fprintf(fid,".SubCkt PlusiGen3Sigma2 L_n1_prime L_n2_prime L_n3_prime R_n1 R_n2 R_n3 TOGND InitV=0V " + ...
                        "VarC0=100%s Var2C0=200%s VarR0=100%s VarR0_2=50%s\n",Cmagnitude,Cmagnitude,Rmagnitude,Rmagnitude);
                    fprintf(fid,"R1 L_n1 R_n2 VarR0  \n");
                    fprintf(fid,"R2 L_n2 R_n1 VarR0  \n");
                    fprintf(fid,"R3 L_n3 R_n3 VarR0  \n");
                    fprintf(fid,"R4 L_n1 R_n1 VarR0_2\n");
                    fprintf(fid,"R5 L_n2 R_n3 VarR0_2\n");
                    fprintf(fid,"R6 L_n3 R_n2 VarR0_2\n");
                    fprintf(fid,"X1_1prime L_n1_prime L_n1 TOGND VoltageFollower \n");
                    fprintf(fid,"X2_2prime L_n2_prime L_n2 TOGND VoltageFollower \n");
                    fprintf(fid,"X3_3prime L_n3_prime L_n3 TOGND VoltageFollower \n");
                    fprintf(fid,".ends PlusiGen3Sigma2\n");
                case '-igen3sigma_2_SOI'
                    fprintf(fid,"* -igen3sigma_2_SOI \n");
                    fprintf(fid,"*\n");
                    fprintf(fid,".SubCkt MinusiGen3Sigma2 L_n1 L_n2 L_n3 R_n1_prime R_n2_prime R_n3_prime TOGND InitV=0V " + ...
                        "VarC0=100%s Var2C0=100%s   VarR0=100%s VarR0_2=100%s\n",Cmagnitude,Cmagnitude,Rmagnitude,Rmagnitude);
                    fprintf(fid,"R1 L_n1 R_n2  VarR0  \n");
                    fprintf(fid,"R2 L_n2 R_n1  VarR0  \n");
                    fprintf(fid,"R3 L_n3 R_n3 VarR0  \n");
                    fprintf(fid,"R4 L_n1 R_n3  VarR0_2\n");
                    fprintf(fid,"R5 L_n2 R_n2 VarR0_2\n");
                    fprintf(fid,"R6 L_n3 R_n1  VarR0_2\n");
                    fprintf(fid,"X1_1prime R_n1_prime R_n1 TOGND VoltageFollower \n");
                    fprintf(fid,"X2_2prime R_n2_prime R_n2 TOGND VoltageFollower \n");
                    fprintf(fid,"X3_3prime R_n3_prime R_n3 TOGND VoltageFollower \n");
                    fprintf(fid,".ends MinusiGen3Sigma2\n");
                case '+igen3sigma_3_SOI'
                    fprintf(fid,"* +igen3sigma_3_SOI \n");
                    fprintf(fid,"*\n");
                    fprintf(fid,".SubCkt PlusiGen3Sigma3 L_n1_prime L_n2_prime L_n3_prime R_n1 R_n2 R_n3 TOGND InitV=0V " + ...
                        "VarC0=100%s Var2C0=200%s VarR0=100%s VarR0_2=50%s\n",Cmagnitude,Cmagnitude,Rmagnitude,Rmagnitude);
                    fprintf(fid,"C1 L_n1 R_n1 VarC0 \n");
                    fprintf(fid,"C2 L_n2 R_n2 VarC0 \n");
                    fprintf(fid,"C3 L_n3 R_n3 VarC0 \n");
                    fprintf(fid,"C4 L_n1 R_n2 Var2C0\n");
                    fprintf(fid,"C5 L_n2 R_n3 Var2C0\n");
                    fprintf(fid,"C6 L_n3 R_n1 Var2C0\n");
                    fprintf(fid,"X1_1prime L_n1_prime L_n1 TOGND VoltageFollower \n");
                    fprintf(fid,"X2_2prime L_n2_prime L_n2 TOGND VoltageFollower \n");
                    fprintf(fid,"X3_3prime L_n3_prime L_n3 TOGND VoltageFollower \n");
                    fprintf(fid,".ends PlusiGen3Sigma3\n");
                case '-igen3sigma_3_SOI'
                    fprintf(fid,"* -igen3sigma_3_SOI \n");
                    fprintf(fid,"*\n");
                    fprintf(fid,".SubCkt MinusiGen3Sigma3 L_n1 L_n2 L_n3 R_n1_prime R_n2_prime R_n3_prime TOGND InitV=0V " + ...
                        "VarC0=100%s Var2C0=200%s VarR0=100%s VarR0_2=50%s\n",Cmagnitude,Cmagnitude,Rmagnitude,Rmagnitude);
                    fprintf(fid,"C1 L_n1 R_n1 VarC0 \n");
                    fprintf(fid,"C2 L_n2 R_n2 VarC0 \n");
                    fprintf(fid,"C3 L_n3 R_n3 VarC0 \n");
                    fprintf(fid,"C4 L_n1 R_n3 Var2C0\n");
                    fprintf(fid,"C5 L_n2 R_n1 Var2C0\n");
                    fprintf(fid,"C6 L_n3 R_n2 Var2C0\n");
                    fprintf(fid,"X1_1prime R_n1_prime R_n1 TOGND VoltageFollower \n");
                    fprintf(fid,"X2_2prime R_n2_prime R_n2 TOGND VoltageFollower \n");
                    fprintf(fid,"X3_3prime R_n3_prime R_n3 TOGND VoltageFollower \n");
                    fprintf(fid,".ends MinusiGen3Sigma3\n");
            end
            fprintf(fid,"* ---------------\n");
            fprintf(fid,"\n");
        end
        function WriteComponent(Component,fid,magnitude)
            arguments
                Component char;
                fid
                magnitude;
            end
            fprintf(fid,"* ---------------\n");
            switch Component
                case 'E_A_LC'
                    fprintf(fid,"* onsite term\n");
                    fprintf(fid,"*\n");
                    fprintf(fid,".SubCkt E_A TONET TOGND  VarC0=100%s VarL0=1u InitV=0V R_L=1u\n",magnitude);
                    fprintf(fid,"C0 TONET TOGND VarC0 IC=InitV\n");
                    fprintf(fid,"L0 TONET TOGND VarL0 R=R_L\n");
                    fprintf(fid,".ends E_A\n");
                case 'hopping_C'
                    fprintf(fid,"* C \n");
                    fprintf(fid,"*\n");
                    fprintf(fid,".SubCkt C n+ n- TOGND  C_hopping=100%s \n",magnitude);
                    fprintf(fid,"Cp n+ n- C_hopping \n");
                    fprintf(fid,".ends C\n");
                case 'hopping_minusC'
                    fprintf(fid,"* minusC \n");
                    fprintf(fid,"*\n");
                    fprintf(fid,".SubCKt minusC n+ n- TOGND C_hopping = 100%s R_0 = 1 R_L=1u\n",magnitude);
                    fprintf(fid,"E_opamp1 Op1out TOGND   n+ Op1-  1E6 max=+100 min=-100\n");
                    fprintf(fid,"E_opamp2 Op2-   TOGND   n+ Op2-  1E6 max=+100 min=-100\n");
                    fprintf(fid,"E_opamp3 Op3out TOGND   n- Op3-  1E6 max=+100 min=-100\n");
                    fprintf(fid,"E_opamp4 Op4-   TOGND   n- Op4-  1E6 max=+100 min=-100\n");
                    fprintf(fid,"C_Op1 n+ Op1out C_hopping\n");
                    fprintf(fid,"R_Op1 Op1out Op1- R_0\n");
                    fprintf(fid,"C_Op3 Op3out n- C_hopping\n");
                    fprintf(fid,"R_Op3 Op3- Op3out R_0\n");
                    fprintf(fid,"R_Op4 Op1- Op4- R_0\n");
                    fprintf(fid,"R_Op2 Op2- Op3- R_0\n");
                    fprintf(fid,".ends minusC\n");
                case 'hopping_minusC_DM'
                    fprintf(fid,"* minusC 2 \n");
                    fprintf(fid,"*\n");
                    fprintf(fid,".SubCKt minusC_DM n+ n- TOGND C_hopping = 100%s C_h = 'SQRT(C_hopping)' L_DM = 1u C_DM = 1\n",magnitude);
                    fprintf(fid,"Xonsite n1 TOGND E_A VarC0 =  C_DM  VarL0 = L_DM\n");
                    fprintf(fid,"C1 n+ n1 C_h\n");
                    fprintf(fid,"C2 n1 n- C_h\n");
                    fprintf(fid,".ends minusC_DM\n");
                case 'E_A_RC'
                    fprintf(fid,"* onsite term\n");
                    fprintf(fid,"*\n");
                    fprintf(fid,".SubCkt E_A TONET TOGND  VarR0=1 VarC0=100%s InitV=0V \n",magnitude);
                    fprintf(fid,"R0 TONET TOGND VarR0\n");
                    fprintf(fid,"C0 TONET TOGND VarC0 IC=InitV\n");
                    fprintf(fid,".ends E_A\n");
                case 'hopping_L'
                    fprintf(fid,"* L \n");
                    fprintf(fid,"*\n");
                    fprintf(fid,".SubCkt L n+ n- TOGND  L_hopping=1%s R_L=1u\n",magnitude);
                    fprintf(fid,"Lp n+ n- L_hopping R=R_L\n");
                    fprintf(fid,".ends L\n");
                case 'hopping_AL'
                    fprintf(fid,"* hopping_AL \n");
                    fprintf(fid,"*\n");
                    fprintf(fid,".SubCKt AL n+ n- TOGND L_hopping2 = 2%s R_L=1u\n",magnitude);
                    fprintf(fid,"E_opamp1 Op1-   TOGND   n+ Op1-  level=1\n");
                    fprintf(fid,"Lp2 Op1- n- L_hopping2 R=R_L\n");
                    fprintf(fid,".ends AL\n");
                case 'hopping_LA'
                    fprintf(fid,"* hopping_LA \n");
                    fprintf(fid,"*\n");
                    fprintf(fid,".SubCKt LA n+ n- TOGND L_hopping2 = 2%s R_L=1u\n",magnitude);
                    fprintf(fid,"E_opamp1 Op1-   TOGND   n- Op1-  level=1 \n");
                    fprintf(fid,"Lp2 Op1- n+ L_hopping2 R=R_L\n");
                    fprintf(fid,".ends LA\n");
                case 'Port_L3'
                    fprintf(fid,"* Lport3 \n");
                    fprintf(fid,"*\n");
                    fprintf(fid,".SubCkt Lport3 TONET TOGND  L_hopping=1%s L_hopping2=2%s L_hopping3=3%s R_L=1u \n",magnitude,magnitude,magnitude);
                    fprintf(fid,"Lp1 TONET TOGND L_hopping R=R_L\n");
                    fprintf(fid,"Lp2 TONET TOGND L_hopping2 R=R_L\n");
                    fprintf(fid,"Lp3 TONET TOGND L_hopping3 R=R_L\n");
                    fprintf(fid,".ends Lport3\n");
                case 'Port_L1'
                    fprintf(fid,"* Lport1 \n");
                    fprintf(fid,"*\n");
                    fprintf(fid,".SubCkt Lport1 TONET TOGND  L_hopping=1%s R_L=1u\n",magnitude);
                    fprintf(fid,"Lp TONET TOGND L_hopping R=R_L \n");
                    fprintf(fid,".ends Lport1\n");
            end
            fprintf(fid,"* ---------------\n");
            fprintf(fid,"\n");
        end
        function Genlib(filename,options)
            arguments
                filename = 'Hckt.lib';
                options.magnitude = 'p';
                options.ComponentLib = ["E_A_LC","hopping_C","hopping_minusC","hopping_minusC_DM"];
                options.ModulesLib = [];
                options.PlugIns = [];
                options.WorkingArea {mustBeMember(options.WorkingArea,{'realTB','nonH','Sigma_SOI'})}= 'realTB';
            end
            magnitude = options.magnitude;
            switch  options.WorkingArea
                case 'realTB'
                    ComponentLib  = ["E_A_LC","hopping_C","hopping_minusC","hopping_minusC_DM"];
                    ModulesLib = options.ModulesLib ;
                    PlugIns = options.PlugIns ;
                case 'nonH'
                    ComponentLib  = ["E_A_RC","hopping_L","hopping_AL","hopping_LA","Port_L3","Port_L1"];
                    ModulesLib = options.ModulesLib ;
                    PlugIns = options.PlugIns ;
                case 'Sigma_SOI'
                    ComponentLib  = [];
                    PlugIns = ["Na-(Mb+Lc)","VoltageFollower"];
                    ModulesLib = [
                        "Basis_SOI",...
                        "+sigma_0" ,...
                        "-sigma_0" ,...
                        "+isigma_0",...
                        "-isigma_0",...
                        "+sigma_1" ,...
                        "-sigma_1" ,...
                        "+gen3sigma_2" ,...
                        "-gen3sigma_2" ,...
                        "+gen3sigma_3",...
                        "-gen3sigma_3",...
                        "+isigma_1_SOI" ,...
                        "-isigma_1_SOI" ,...
                        "+igen3sigma_2_SOI" ,...
                        "-igen3sigma_2_SOI",...
                        "+igen3sigma_3_SOI",...
                        "-igen3sigma_3_SOI" ...
                        ];
                otherwise
                    ComponentLib = options.ComponentLib;
                    ModulesLib = options.ModulesLib ;
                    PlugIns = options.PlugIns ;
            end
            fid = fopen(filename,'w');
            for PlugIn = PlugIns
                Hckt.WritePlugIns(char(PlugIn),fid,magnitude);
            end
            for Modules = ModulesLib
                Hckt.WriteModules(char(Modules),fid,magnitude);
            end
            for Component = ComponentLib
                Hckt.WriteComponent(char(Component),fid,magnitude);
            end
        end
        function Portlist_gen(fid,AllPortStrL,options,options2)
            arguments
                fid = 1;
                AllPortStrL = [""];
                options.type = 'v';
                options2.comment = false;
                options2.probenode = 'Allnode';
                options2.prefix = "+ ";
                options2.father = "";
                options2.nodelist = [];
            end
            if options2.comment
                switch options2.probenode
                    case 'Allnode'
                        for i = 1:numel(AllPortStrL)
                            fprintf(fid,"*"+options2.prefix+options.type+"(%s)\n",AllPortStrL(i));
                        end
                    case 'Selectnode'
                        for i = 1:numel(AllPortStrL)
                            for j = 1:numel(options2.nodelist)
                                fprintf(fid,"*"+options2.prefix+options.type+"(%s"+options2.father+"%s)\n",AllPortStrL(i),options2.nodelist(j));
                            end
                        end
                end
            else
                switch options2.probenode
                    case 'Allnode'
                        for i = 1:numel(AllPortStrL)
                            fprintf(fid,options2.prefix+options.type+"(%s)\n",AllPortStrL(i));
                        end
                    case 'Selectnode'
                        for i = 1:numel(AllPortStrL)
                            for j = 1:numel(options2.nodelist)
                                fprintf(fid,options2.prefix+options.type+"(%s"+options2.father+"%s)\n",AllPortStrL(i),options2.nodelist(j));
                            end
                        end
                end
            end
        end
    end
    %% Modify
    methods
        function HcktObj = autohermi(HcktObj)
            arguments
                HcktObj Hckt;
            end
            HcktObj_tmp = HcktObj;
            for i = 2:HcktObj.NRPTS
                % Duality_vector_dist();
                vector_tmp = HcktObj.vectorAll(HcktObj.vectorL(i),:);
                portin_tmp = HcktObj.PortInCell{i,1};
                portout_tmp = HcktObj.PortOutCell{i,1};
                vector_tmp_oppo = -vector_tmp;
                portin_tmp_oppo = portout_tmp;
                portout_tmp_oppo = portin_tmp;
                % even vector 
                [~,vector_tmp_oppo_label]= ismember(vector_tmp_oppo,HcktObj.vectorAll,'rows');
                if vector_tmp_oppo_label == 0
                    fprintf('The opposite vector hamilton does not exist, build it : %s!\n', ...
                       string(HcktObj.ScktL(i).name));
                    HcktObj_tmp = HcktObj_tmp.set_hop(vector_tmp_oppo,HcktObj.ScktL(i),portin_tmp_oppo,portout_tmp_oppo);
                    continue;
                end
                % check vector
                vector1 = vector_tmp_oppo_label == HcktObj.vectorL;
                i1 = park.checkcell(portin_tmp_oppo,HcktObj.PortInCell);
                j1 = park.checkcell(portout_tmp_oppo,HcktObj.PortInCell);
                j = find(i1&j1&vector1);
                % homecell
                if i == j

                end
                % 
                if j == 0
                    fprintf('The opposite vector hamilton does not exist, build it : %s!\n', ...
                        string(HcktObj.ScktL(i).name));
                    HcktObj_tmp = HcktObj_tmp.set_hop(vector_tmp_oppo,HcktObj.HcktObj.ScktL(i),portin_tmp_oppo,portout_tmp_oppo);
                    continue;
                end
                %
                if ~(HcktObj.ScktL(i) == HcktObj.ScktL(j))
                    disp([vector_tmp vector_tmp_oppo]);
                   fprintf('The opposite vector Subckt does not exist, build it with : %s vs %s!\n',HcktObj.ScktL(i).name,HcktObj.ScktL(j).name);
                end
            end
            HcktObj = HcktObj_tmp;
        end
        function HcktObj = half(HcktObj)
            vectorLcheck = HcktObj.dim2vectorL(HcktObj.dim);
            [whichexist,~] = ismember(HcktObj.vectorAll ,vectorLcheck,'rows') ;
            ToKeepRef = find(whichexist);
            ToKeepL = ismember(HcktObj.vectorL,ToKeepRef);
            HcktObj = HcktObj.reseq(ToKeepL);
        end
        function HcktObj = reseq(HcktObj,seqL)
            HcktObj.ScktL = HcktObj.ScktL(seqL,:);
            HcktObj.vectorL = HcktObj.vectorL(seqL);
            HcktObj.PortInCell = HcktObj.PortInCell(seqL,:);
            HcktObj.PortOutCell = HcktObj.PortOutCell(seqL,:);
            HcktObj.DescriptionL = HcktObj.DescriptionL(seqL,:);
        end
        function HcktObj = clean(HcktObj)
            % hoppinglist
            hollowlabel = find(HcktObj.ScktL == 0);
            HcktObj = HcktObj.reseq(hollowlabel);
        end
    end
    %%
    methods(Static)
        function vectorL = dim2vectorL(dim)
            switch dim
                case 1
                    vectorL= [0;1];
                case 2
                    vectorL= [0,0;1,0;0,1;1,1;1,-1];
                case 3
                    vectorL= [...
                        0,0,0;...
                        1,0,0;...
                        0,1,0;...
                        0,0,1;...
                        1,1,0;...
                        0,1,1;...
                        1,0,1;...
                        -1,1,0;...
                        0,-1,1;...
                        1,0,-1;...
                        1,1,1;...
                        ];
                case 4
                    % Todo sb. fix it
                    vectorL= [...
                        0,0,0;...
                        1,0,0;...
                        0,1,0;...
                        0,0,1;...
                        1,1,0;...
                        0,1,1;...
                        1,0,1;...
                        -1,1,0;...
                        0,-1,1;...
                        1,0,-1;...
                        1,1,1;...
                        ];
                case 5
            end
        end
    end
end