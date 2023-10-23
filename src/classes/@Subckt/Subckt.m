classdef Subckt < matlab.mixin.CustomDisplay
    %SUBCKT 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
        portin = [];
        portout = [];
        portname = [""];
        name = '';
        ReverseConnection = false; 
        netlist cell  = {["X1","+","-","Sckt"]};
    end
    properties
        pinned = false;
        hollow = false;
        commute = true;
    end
    properties
        deviceL {mustBeMember(deviceL,{'C','L','R','X','K','D','Q','M','V','I','G',''})} = '';% 
        %         R
        %         C
        %         L
        %         X
        %         K
        %         D
        %         Q
        %         M
        %         V
        %         I
        %         G
        %         E
    end
    properties % for set
        OnsiteLSubscktL;
        HoppingSubscktL;
    end
    properties
        devicenameL ;
        type ;
        varL  ;
        paraL ;
    end
    properties
        description;
    end
    properties
        Options = '';
        Lib = '';
    end
    properties(Dependent)
        innerOutput;
        OutputStr;
        pinnedOutput;
        Output;
        Description;
        ToGND;
    end
    properties
        toGND = true;
        pinnedOutputStr;
    end
    %% Define which properties show
    methods (Access = protected)
        function propgrp = getPropertyGroups(~)
            proplist = {'name','pinned','portname','deviceL','Output'};
            propgrp = matlab.mixin.util.PropertyGroup(proplist);
        end
    end
    %% Construction
    methods
        function ScktObj = Subckt(device,name,portname,description,netlist,options)
            arguments
                device 
                name = 'Sckt';
                portname string= ["+","-"];
                description string ="";
                netlist cell = {["X1","+","-","Sckt"]};
                options.pinned = true;
                options.enforce_pinned = false;
                options.magicnumber = 4;
            end
            if isempty(device)
                ScktObj.hollow= true;
                return;
            end
            if options.enforce_pinned
                ScktObj.pinned = options.pinned;
                ScktObj.pinnedOutputStr = device;
                return;
            elseif  options.pinned && nargin == 1
                ScktObj.pinned = options.pinned;
                ScktObj.pinnedOutputStr = device;
            end
            switch class(device)
                case 'char'
                    inputchar = strtrim(device);
                    if length(inputchar) == 1
                        
                    else
                        StrMtmp = park.char2mat(inputchar);
                        tmpchar = char(StrMtmp(1));
                        if size(StrMtmp,1) > 1
                            ScktObj.type = 'X';
                            if tmpchar(1) == '.'
                                % remove title
                                Information = StrMtmp(1,:);
                                StrMtmp(1,:) = [];
                                % remove ends
                                tmpchar2 = char(StrMtmp(end,1));
                                if tmpchar2(1) == '.'
                                    StrMtmp(end,:) = [];
                                end
                            else
                                Information = '';
                            end
                            % node
                            nodeList= '';
                            % name
                            if ~strcmp(Information,'')
                                InformationL = split(Information,' ');
                                [description,InformationL1] = park.ExtractContainPat(InformationL);
                                name = InformationL1(2);
                                portname = InformationL1(3:end);
                            else
                                %name = 'Subckt0';
                                %portname =unique(nodeList);
                                %description = [''];
                            end
                            % device
                            [device,ScktObj.devicenameL] = park.ExtractCharFromStr(StrMtmp(:,1));
                        else
                            device = tmpchar(1);
                            ScktObj.type = device;
                            magicnumber =options.magicnumber;
                            name = StrMtmp(magicnumber:end);
                            portname = StrMtmp(2:magicnumber-1);
                            portin = 1:floor(length(portname)/2);
                            portout = floor(length(portname)/2)+1:length(portname);
                            ScktObj.portin = portin;
                            ScktObj.portout = portout;
                        end
                        if size(StrMtmp,2) > 1
                            for i = 1:size(StrMtmp,1)
                                netlist{i} = StrMtmp(i,:);
                            end
                        end
                    end
                case 'string'
                    inputchar = char(device);
                    if length(inputchar) == 1

                    else
                        StrMtmp = park.char2mat(inputchar);
                        tmpchar = char(StrMtmp(1));
                        if size(StrMtmp,1) > 1
                            ScktObj.type = 'X';
                            if tmpchar(1) == '.'
                                % remove title
                                Information = StrMtmp(1,:);
                                StrMtmp(1,:) = [];
                                % remove ends
                                tmpchar2 = char(StrMtmp(end,1));
                                if tmpchar2(1) == '.'
                                    StrMtmp(end,:) = [];
                                end
                            else
                                Information = '';
                            end
                            % node
                            nodeList= '';
                            % name
                            if ~strcmp(Information,'')
                                InformationL = split(Information,' ');
                                [description,InformationL1] = park.ExtractContainPat(InformationL);
                                name = InformationL1(2);
                                portname = InformationL1(3:end);
                            else
                                %name = 'Subckt0';
                                %portname =unique(nodeList);
                                %description = [''];
                            end
                            % device
                            [device,ScktObj.devicenameL] = park.ExtractCharFromStr(StrMtmp(:,1));
                        else
                            device = tmpchar(1);
                            ScktObj.type = device;
                            magicnumber =options.magicnumber;
                            name = StrMtmp(magicnumber:end);
                        end
                        if size(StrMtmp,2) > 1
                            for i = 1:size(StrMtmp,1)
                                netlist{i} = StrMtmp(i,:);
                            end
                        end
                    end
                otherwise 
            end
            ScktObj.deviceL = device;
            ScktObj.name = name;
            ScktObj.portname = portname;
            ScktObj.description = description;
            ScktObj.netlist = netlist;
        end
        function outputArg = method1(obj,inputArg)
            %METHOD1 此处显示有关此方法的摘要
            %   此处显示详细说明
            outputArg = obj.port + inputArg;
        end
    end
    methods(Static)
        function ScktObj = FromHomecellList(HnumL,vectorL,WAN_NUM,options)
            arguments
                HnumL
                vectorL = [];
                WAN_NUM = 4;
                options.magnitude = 'p';
                options.parameters = [];
                options.name = 'Pri';
                options.Dim = 3;
                options.prefix = 100;
            end
            if isempty(options.parameters)
                options.parameters = ['VarL0 = 1u C0 = 100',options.magnitude,' InitV = 0V'];
            end
            optionsCell = namedargs2cell(options);
            switch class(HnumL)
                case 'double'
                    sym_mode = false;
                case 'sym'
                    sym_mode = true;
                    DefaultC = '100p';
                    error('Not support yet');
                case 'HR'
                    H_hr = HnumL;
                    options.Dim = H_hr.Dim;
                    optionsCell = namedargs2cell(options);
                    ScktObj = Subckt.FromHomecellList(H_hr.HnumL,H_hr.vectorL,H_hr.WAN_NUM,optionsCell{:});
                    sym_mode = false;
                    return;
                otherwise
                    sym_mode = false;
            end
            % 
            ScktObjDevice = 'X';
            ScktObjName = options.name;
            ScktObjNode = [string(1:WAN_NUM),'TOGND'];
            ScktObjNetlist = Subckt.GenNetlistFromHList(HnumL,vectorL,'sym_mode',sym_mode,'magnitude',options.magnitude,'Dim',options.Dim,'prefix',options.prefix);
            ScktObjDescription = options.parameters;
            ScktObj = Subckt(ScktObjDevice,ScktObjName,ScktObjNode,ScktObjDescription,ScktObjNetlist);% 
        end
        function Netlist = GenNetlistFromHList(HnumL,vectorL,options)
            arguments
                HnumL
                vectorL
                options.magnitude = 'p';
                options.sym_mode = false;
                options.Dim = 3;
                options.vectorL = [0,0,0];
                options.prefix = 100;
            end
            DIM = options.Dim;
            if size(options.vectorL ,2) ~= DIM
                options.vectorL = zeros(1,DIM);
            end
            % Select vector default is homecell
            [SelectL] = ismember(double(vectorL(:,1:DIM)),options.vectorL,'rows');
            vectorL = vectorL(SelectL,:);
            HnumL = HnumL(SelectL);
            nComponents = size(vectorL,1);
            Netlist{nComponents} = '';
            selectL = logical(1:nComponents);
            for n = 1:nComponents
                i = vectorL(n,DIM+1);
                j = vectorL(n,DIM+2);
                Name = string(['X',num2str(i),num2str(j)]);
                if i == j
                    PortL = string(i);
                    TOGND = "TOGND";
                    Ref = "E_A";
                    Description = "VarC0 = "+ num2str(HnumL(n)*options.prefix) + options.magnitude;
                elseif i<j
                    PortL = string([i,j]);
                    TOGND = "TOGND";
                    if HnumL(n) < 0 
                        Ref = "C";
                    else
                        Ref = "MinusC";
                    end
                    Description = "C_hopping = "+ num2str(abs(HnumL(n))*options.prefix) + options.magnitude;
                else
                    %Hermition!
                    selectL(n) = false;
                    continue;
                end
                Netlist{n} = [Name,PortL,TOGND,Ref,Description];
            end
            Netlist = Netlist(selectL);
            %                 'X1 n+_1 TOGND E_A VarC0=Cc\n',...
            %                 'X12 n+_1 n+_2 TOGND C\n',...
            %                 'X2 n+_2 TOGND E_A VarC0=C0\n' ...
            %                 'X23 n+_2 n+_3 TOGND C\n',...
            %                 'X3 n+_3 TOGND E_A VarC0=C0\n',...
            %                 'X34 n+_3 n+_4 TOGND C\n',...
            %                 'X4 n+_4 TOGND E_A VarC0=Cc\n' ...
            %                 'X41 n+_4 n+_1 TOGND minusC\n',...
        end
    end

    % get
    methods 
        function ToGND = get.ToGND(ScktObj)
            if ScktObj.type == 'X'
                ToGND = ScktObj.toGND;
            else
                ToGND = false;
            end

        end
        function Description = get.Description(ScktObj)
            Description = ScktObj.description;
        end
        function innerOutput = get.innerOutput(ScktObj)
            innerOutput = park.char2list(ScktObj.OutputStr);
        end
        function OutputStr = get.OutputStr(ScktObj)
            if ScktObj.pinned
                OutputStr = ScktObj.pinnedOutputStr;
            else
                PortName = "";
                for i = 1:length(ScktObj.portname)
                    PortName = strcat(PortName,string(ScktObj.portname(i))," ");
                end
                OutputStr = strcat('.SubCkt'," ",ScktObj.name," ",...
                    PortName,ScktObj.description,'\n');
                endline = strcat('.ends'," ",ScktObj.name,"\n");
                OutputStr = strcat(OutputStr,park.cellcatwithspice(ScktObj.netlist));
                OutputStr =  strcat(OutputStr,endline);
            end
        end
        function pinnedOutput = get.pinnedOutput(ScktObj)
            pinnedOutput = park.char2list(ScktObj.pinnedOutputStr);
        end
        function Output = get.Output(ScktObj)
            if ScktObj.pinned
                Output = ScktObj.pinnedOutput;
            else
                Output = ScktObj.innerOutput;
            end
        end
    end
    %% 重载
    methods
        function C = eq(A,B)
            if length(A) ==1 && length(B) >1
                C = logical(1:length(B));
                for i = 1:numel(B)
                    C(i) = eq(B(i),A);
                end
                return;
            elseif length(A) >1 && length(B) ==1
                C = logical(1:length(A));
                for i = 1:numel(A)
                    C(i) = eq(A(i),B);
                end
                return;
            elseif length(A) >1 && length(B) >1
                error('Todo')
            else
                if xor(isa(A,'Subckt'),isa(B,'Subckt'))
                    if (A.hollow && B ==0)||(B.hollow && A ==0)
                        C = true;
                    end
                    return;
                elseif isa(A,'Subckt') && isa(B,'Subckt')
                    C = true;
                    if ~strcmp(A.name,B.name)
                        C = false;
                        return;
                    end
                    if ~strcmp(A.type,B.type)
                        C = false;
                        return;
                    end
                else
                    C = false;
                end
            end
        end
    end
    %% tools
    methods
        function ScktObjStr = fprintf(ScktObj,fileID)
            arguments
                ScktObj Subckt;
                fileID  double=1; % 从 fopen 获取的文件标识符。1 表示标准输出（屏幕）。2 表示标准错误
            end
            ScktObjStr = '';
            for i = 1:numel(ScktObj.Output)
                ScktObjStrTmp = ScktObj.Output{i};
                fprintf(fileID,'%s\n',ScktObjStrTmp);
                ScktObjStr = [ScktObjStr,ScktObjStrTmp,'\n'];
            end
        end
    end

end

