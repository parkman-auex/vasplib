classdef HR <vasplib & matlab.mixin.CustomDisplay
    %% HR - a powerful TB tools
    %
    % * Label: class TB model
    %
    %% Description of the dlass:
    % Implements the Model class, which describes a tight-binding model.
    %
    %% Usage:
    %
    % *
    %
    %% example:
    %
    %% Note:
    % For the number part, it is a direct copy of tbmodels
    % <https://github.com/Z2PackDev/TBmodels>
    % * (c) 2015-2018, ETH Zurich, Institut fuer Theoretische Physik
    % * Author: Dominik Gresch <greschd@gmx.ch>
    % For the symbolic part, add the feature consistent with vasplib
    %
    %% Change log
    %
    % * Document Date: 2021/02/23
    % * Creation Date: 2021/01/04
    % * Last updated : 2021/04/07
    %
    %% Copyright
    %
    % * parkman
    % * <parkman@buaa.edu.cn>
    %
    
    %% public properties
    properties
        vectorL  ; % int8 -128 ~ 128 % the Rvector list ;int64 list; int32 will be used more widely
        HnumL   ; % Hnum_list
        HcoeL   ; % Hcoe_list
    end
    %% private properties
    properties %(GetAccess = protected,Hidden = true)
        Duality_vector_dist; % the opposite vector sequnce dictionary.
    end
    properties %(GetAccess = protected)
        Type;            % the type of this TB : mat(default)  sparse(for large system) list(for symmetry)
        overlap logical= false; % overlap(for Atom orbital)
        num logical= false;        % num(numeric).
        coe logical= true;        % coe(symbolic).
        soc logical= false;     % spinful model?
        AvectorL;
        BvectorL;
        CvectorL;  %test
        vectorhopping = false;  % use for symmetrize for speed up symbolic;
    end
    %% temp properties
    properties (Transient,Hidden = true)
        R_vector_dist;  % the R vector sequnce dictionary.
    end
    %% dynamic properties dependent
    properties(Dependent = true)
        % Basis_num;
        NRPTS    ; % the total number of H(Rn)
        WAN_NUM  ; % the number of wannier like orbs
        Line_000 ; % the sequence of homecell
        homecell ; 
    end
    properties(Dependent = true,Hidden = true)
        
    end
    %% hidden properties
    properties (Hidden = true)
        % for spacegroup
        nn_store_smart   ; % nn_store
        nn_sparse_n      ;
        Atom_store_smart ;
        Rnn_map          ;
        ScoeL            ;
        SnumL            ;
        vectorL_overlap  ;
    end
    %% Define which properties show
    methods (Access = protected)
        function propgrp = getPropertyGroups(~)
            proplist = {'WAN_NUM','NRPTS','Type','HcoeL','HnumL','vectorL'};
            propgrp = matlab.mixin.util.PropertyGroup(proplist);
        end
    end
    %     events
    %       InsufficientFunds;
    %     end
    %% construction method
    % --------------------- construction method -----------------------
    methods
        function H_hr = HR(WAN_NUM,vectorL,options)
            % HR H_hr = HR(WAN_NUM,vectorL) construct a empty TB obj
            %   H_hr = HR() construct a 4-orbs empty TB obj
            %   H_hr = HR(WAN_NUM) construct a WAN_NUM-orbs empty TB obj
            %   H_hr = HR(WAN_NUM,vectorL) construct a WAN_NUM-orbs with vectorL H_R empty TB obj
            %   H_hr = HR(WAN_NUM,vectorL,HnumL,HcoeL) construct a TB obj
            %   with full information
            %   H_hr = HR(WAN_NUM,vectorL,HnumL,HcoeL,Type) indicate the
            %   Type of this TB obj.
            %   See also FROM_HSTRUCT, FROM_HSPARSE, FROM_HDF5, FROM_WANNIER90.
            
            % ----------- nargin ----------
            arguments
                WAN_NUM double{mustBeInteger} =4;
                vectorL = int32([0 ,0 ,0]);
                options.HnumL double=[]   ;
                options.HcoeL sym=sym([]) ;
                options.Type char = 'mat' ;
                options.overlap logical = false;
                options.SnumL double=[]   ;
                options.ScoeL sym=sym([]) ;
                options.vectorL_overlap = int32([0 ,0 ,0]);
                options.sym = true;
            end
            H_hr = H_hr@vasplib();
            H_hr.Basis_num = WAN_NUM;
            % -------------check---------------
            Type = options.Type;
            % vectorL
            if strcmp(Type,'list')
                if nargin < 2
                    vectorL = int32([0 ,0 ,0,WAN_NUM,WAN_NUM]);
                end
            else
                if size(vectorL,2) == 1
                    switch vectorL
                        case 1
                            tmp_vectorL = [0,0,0];
                        case 3
                            tmp_vectorL = [-1,0,0;0 0 0;1 0 0 ];
                        case 5
                            tmp_vectorL = [1,0,0;0 1 0;0 0 0;0 -1 0;-1 0 0];
                        case 7
                            tmp_vectorL = [0 0 1;1,0,0;0 1 0;0 0 0;0 -1 0;-1 0 0;0 0 -1];
                        case 9
                            tmp_vectorL = [1,0,0;0 1 0;0 0 0;0 -1 0;-1 0 0;...
                                1 1 0;1 -1 0;-1 1 0;-1 -1 0;];
                        case 11
                            tmp_vectorL = [0,0,0];
                        case 27
                            tmp_vectorL = [0,0,0];
                    end
                    vectorL = int32(tmp_vectorL);
                end
            end
            % HnumL
            if isempty(options.HnumL)
                [NRPTS,~]  = size(vectorL);
                if strcmp(Type,'sparse')
                    HnumL{NRPTS} = sparse(WAN_NUM,WAN_NUM);
                    for i =1:NRPTS-1
                        HnumL{i} = sparse(WAN_NUM,WAN_NUM);
                    end
                elseif strcmp(Type,'list')
                    HnumL = zeros(1,1);
                else
                    HnumL = zeros(WAN_NUM,WAN_NUM,NRPTS);
                end
            else
                HnumL = options.HnumL;
            end
            % HcoeL
            if options.sym
                if isempty(options.HcoeL)
                    if strcmp(Type,'sparse')
                        HcoeL = sym([]);
                    else
                        HcoeL = sym(HnumL);
                    end
                else
                    HcoeL = options.HcoeL;
                end
            else
                HcoeL = options.HcoeL;
            end
            if options.overlap
                if isempty(options.SnumL )
                    if strcmp(Type,'sparse')
                        SnumL{NRPTS} = sparse(WAN_NUM,WAN_NUM);
                        for i =1:NRPTS-1
                            SnumL{i} = sparse(WAN_NUM,WAN_NUM);
                        end
                    else
                        SnumL = zeros(WAN_NUM,WAN_NUM,NRPTS);
                    end
                    H_hr.SnumL = SnumL;
                else
                    H_hr.SnumL = options.SnumL  ; % Snum_list
                end
                if isempty(options.ScoeL ) && isempty(options.SnumL )
                    if strcmp(Type,'sparse')
                        ScoeL = [];
                    else
                        ScoeL = sym(SnumL);
                    end
                    H_hr.ScoeL = ScoeL;
                else
                    H_hr.ScoeL = options.ScoeL  ; % Scoe_list
                end
                H_hr.vectorL_overlap =  int32(options.vectorL_overlap);
            end
            %
            % H_hr.NRPTS   =  NRPTS; % the total number of H(Rn)
            % H_hr.WAN_NUM =  WAN_NUM ; % the num of wannier like orbs
            H_hr.vectorL = vectorL; % the Rvector list
            H_hr.HnumL = HnumL  ; % Hnum_list
            H_hr.HcoeL = HcoeL  ; % Hcoe_list
            H_hr.Type  = Type   ;
            %             if strcmp(Type,'sparse')
            %                 H_hr.orbL  = sparse(WAN_NUM,3);
            %             else
            %                 H_hr.orbL  = zeros(WAN_NUM,3);
            %             end
            H_hr.orbL  = zeros(WAN_NUM,3);
            H_hr.overlap = options.overlap;
            %H_hr.Duality_vector_dist = containers.Map('KeyType','double','ValueType','double');
        end
    end
    % --------------------- external construction method -----------------------
    methods (Static)
        % --------------------- construction method -----------------------
        function H_hr = from_Hstruct(Hstruct)
            [NRPTS,~]=size(Hstruct);
            try
                WAN_NUM = length(Hstruct(1).Hnum);
            catch
                WAN_NUM = length(Hstruct(1).Hcoe);
            end
            V = [Hstruct.vector];
            vectorL =int32(reshape(V,3,length(V)/3)');
            HnumL = zeros(WAN_NUM,WAN_NUM,NRPTS);
            HcoeL = sym(zeros(WAN_NUM,WAN_NUM,NRPTS));
            for i = 1:NRPTS
                try
                    HnumL(:,:,i) = Hstruct(i).Hnum/Hstruct(i).Degen;
                catch
                    HnumL(:,:,i) = zeros(WAN_NUM)  ;
                end
                try
                    HcoeL(:,:,i) = Hstruct(i).Hcoe;
                catch
                    
                end
            end
            H_hr = HR(WAN_NUM,vectorL,'HnumL',HnumL,'HcoeL',HcoeL);
        end
        function H_hr = from_Hsparse(Hsparse)
            vectorL = Hsparse.vectorL;
            [NRPTS,~] = size(vectorL);
            WAN_NUM = length(Hsparse.HnumL{1});
            HcoeL = [];
            HnumL = zeros(WAN_NUM,WAN_NUM,NRPTS);
            for i = 1:NRPTS
                HnumL(:,:,i) = full(Hsparse.HnumL{i});
            end
            H_hr = HR(WAN_NUM,vectorL,'HnumL',HnumL,'HcoeL',HcoeL);
        end
        function H_hr = from_POSCAR_SE(POSCAR_file,options)
            %POSCAR_file,r_max_search,level_cut,onsite_mode,Accuracy,search_range,mode)
            % ---------------------   nargin   ------------------------
            arguments
                POSCAR_file char= 'POSCAR';
                options.WAN_NUM double{mustBeInteger}= -1;
                options.r_max = 3.3;
                options.level_cut = -1;
                options.per_dir double = [1,1,1];
                options.chiral logical = false;
                options.spin logical = false;
                options.deltarule  double{mustBeMember(options.deltarule,[0,1,2])}= 0;
                options.alpharule  double{mustBeMember(options.alpharule,[0,1,2])}= 0;
                options.onsite logical = true;
                options.Accuracy = 1e-3;
                options.search_range = [2 2 2];
                options.Type char{mustBeMember(options.Type,{'mat','list','sparse'})}= 'mat';
                options.overlap = false;
                options.symbolic = false;
                options.simplify = false;
                options.Rd  = -1;
                options.para struct=struct();
                options.silence = true;
                options.E0 = 0;
            end
            r_max_search =options.r_max ;
            level_cut  = options.level_cut;
            per_dir = options.per_dir;
            onsite = options.onsite;
            Accuracy = options.Accuracy ;
            search_range = options.search_range ;
            Type = options.Type;
            %
            [~,sites,~,~,~]=vasplib.POSCAR_readin(POSCAR_file,'tbsk');
            if options.WAN_NUM == -1
                H_hr = HR(length(sites),'overlap',options.overlap,'Type',Type);
            else
                H_hr = HR(options.WAN_NUM,'overlap',options.overlap,'Type',Type);
            end
            H_hr.Basis_num = H_hr.WAN_NUM;
            
            if options.symbolic
                H_hr =  POSCAR_file > H_hr;
            else
                H_hr = H_hr < POSCAR_file;
            end
            %H_hr.overlap = options.overlap;
            switch Type
                case 'mat'
                    H_hr.Type = 'mat';
                    H_hr = H_hr.nn(search_range,Accuracy,r_max_search);
                    H_hr = H_hr.H_TBSK_gen(...
                        'chiral', options.chiral,...
                        'spin',options.spin ,...
                        'onsite',logical(onsite),...
                        'per_dir',per_dir ,...
                        'level_cut',level_cut ...
                        );
                    if isequal(sym(options.E0),(0))
                    else
                        H_hr = H_hr.set_hop_mat(eye(H_hr.WAN_NUM)*options.E0,[0,0,0],'symadd');
                    end
                case 'list'
                    H_hr.Type = 'list';
                    %H_hr = H_hr.rewrite();
                    H_hr = H_hr.nn(search_range,Accuracy,r_max_search);
                    % H_hr.Type = 'list';
                    H_hr = H_hr.H_TBSK_gen(...
                        'chiral', options.chiral,...
                        'spin',options.spin ,...
                        'onsite',logical(onsite),...
                        'per_dir',per_dir ,...
                        'level_cut',level_cut ...
                        );
                    if isequal(sym(options.E0),(0))
                    else
                        H_hr = H_hr.set_hop_mat(eye(H_hr.WAN_NUM)*options.E0,[0,0,0],'symadd');
                    end
                case 'sparse'
                    H_hr = H_hr.nn(search_range,Accuracy,r_max_search);
                    H_hr = H_hr.H_TBSK_gen_sparse(...
                        'chiral', options.chiral,...
                        'spin',options.spin ,...
                        'onsite',logical(onsite),...
                        'per_dir',per_dir ,...
                        'level_cut',level_cut,...
                        'Rd',options.Rd,...
                        'para',options.para,...
                        'deltarule',options.deltarule,...
                        'alpharule',options.alpharule ...
                        );
                    if isequal(sym(options.E0),(0))
                    else
                        H_hr = H_hr.set_hop_mat(eye(H_hr.WAN_NUM)*options.E0,[0,0,0],'add');
                    end
            end
            if ~strcmp(Type,'sparse')
                if options.deltarule
                    fprintf('applying delta rule ...\n');
                    H_hr = H_hr.deltarule(level_cut,options.deltarule,'Rd',options.Rd);
                end
                if options.alpharule
                    fprintf('applying alpha rule ...\n');
                    H_hr = H_hr.alpharule(level_cut,options.alpharule,'Rd',options.Rd,'silence',options.silence);
                end
                if options.simplify
                    fprintf('simplify the Hamitoian...\n');
                    H_hr =H_hr.simplify();
                end
                disp(H_hr.symvar_list);
                if H_hr.overlap
                    disp(symvar(H_hr.ScoeL));
                end
            end
        end
        function H_hr = from_hdf5(filename)
        end
        function H_hr = from_wannier90(filename,Type,options)
            arguments
                filename  = 'wannier90_hr.dat';
                Type  char = 'mat';
                options.Accuracy = 1e-6;
                options.overlap = false
            end
            %% read hopping terms
            if options.overlap
                [dataArray,NRPT_list,NRPTS,NUM_WAN]=HR.hrdat_read(filename{1});
                [dataArray2,NRPT_list_S,~,~]=HR.hrdat_read(filename{2});
                if strcmp(Type ,'mat')
                    %%
                    Vec_Fir = dataArray{:, 1};      % a1 direction
                    Vec_Sec = dataArray{:, 2};      % a2 direction
                    Vec_Thi = dataArray{:, 3};      % a3 direction
                    %Orb_fir = dataArray{:, 4};
                    %Orb_sec = dataArray{:, 5};
                    h_real = dataArray{:, 6};
                    h_imag = dataArray{:, 7};
                    %% Clear temporary variables
                    % read
                    V_f = reshape(Vec_Fir, NUM_WAN*NUM_WAN, NRPTS);
                    V_s = reshape(Vec_Sec, NUM_WAN*NUM_WAN, NRPTS);
                    V_t = reshape(Vec_Thi, NUM_WAN*NUM_WAN, NRPTS);
                    vectorL=[V_f(1,:)',V_s(1,:)',V_t(1,:)'];
                    HnumL_real = reshape(h_real, NUM_WAN,NUM_WAN, NRPTS);
                    HnumL_imag = reshape(h_imag, NUM_WAN,NUM_WAN, NRPTS);
                    HnumL = HnumL_real+1i*HnumL_imag;
                    if strcmp(Type,'n')
                        for i = 1 : NRPTS
                            HnumL(:,:,i) = HnumL(:,:,i)/NRPT_list(i);
                            %From here , get hints
                        end
                    else
                    end
                    %% S
                    Vec_Fir = dataArray2{:, 1};      % a1 direction
                    Vec_Sec = dataArray2{:, 2};      % a2 direction
                    Vec_Thi = dataArray2{:, 3};      % a3 direction
                    %Orb_fir = dataArray{:, 4};
                    %Orb_sec = dataArray{:, 5};
                    h_real = dataArray2{:, 6};
                    h_imag = dataArray2{:, 7};
                    %% Clear temporary variables
                    % read
                    V_f = reshape(Vec_Fir, NUM_WAN*NUM_WAN, NRPTS);
                    V_s = reshape(Vec_Sec, NUM_WAN*NUM_WAN, NRPTS);
                    V_t = reshape(Vec_Thi, NUM_WAN*NUM_WAN, NRPTS);
                    vectorL_overlap=[V_f(1,:)',V_s(1,:)',V_t(1,:)'];
                    SnumL_real = reshape(h_real, NUM_WAN,NUM_WAN, NRPTS);
                    SnumL_imag = reshape(h_imag, NUM_WAN,NUM_WAN, NRPTS);
                    SnumL = SnumL_real+1i*SnumL_imag;
                    if strcmp(Type,'n')
                        for i = 1 : NRPTS
                            SnumL(:,:,i) = SnumL(:,:,i)/NRPT_list(i);
                            %From here , get hints
                        end
                    else
                    end
                    HcoeL = sym([]);
                    ScoeL = sym([]);
                    %HcoeL = sym(zeros(size(HnumL)));
                    %ScoeL = sym(zeros(size(SnumL)));
                elseif strcmp(Type ,'list')
                    DATAARRAY = cell2mat(dataArray(1:7));
                    HOPARRAY = DATAARRAY(:,6)+1i*DATAARRAY(:,7);
                    % NRPT_lists = kron(NRPT_list,zeros(NUM_WAN*NUM_WAN,1));
                    HnumL_select = abs(HOPARRAY)>options.Accuracy;
                    vectorL = DATAARRAY(HnumL_select,1:5);
                    [~,~,original_label] = unique(vectorL(:,1:3),'rows');
                    HnumL = HOPARRAY(HnumL_select)./NRPT_list(original_label);
                    HcoeL = sym([]);
                    DATAARRAY = cell2mat(dataArray2(1:7));
                    OVERLAPARRAY = DATAARRAY(:,6)+1i*DATAARRAY(:,7);
                    % NRPT_lists = kron(NRPT_list,zeros(NUM_WAN*NUM_WAN,1));
                    SnumL_select = abs(OVERLAPARRAY)>options.Accuracy;
                    vectorL_overlap = DATAARRAY(SnumL_select,1:5);
                    [~,~,original_label] = unique(vectorL_overlap(:,1:3),'rows');
                    SnumL = OVERLAPARRAY(SnumL_select)./NRPT_list_S(original_label);
                    ScoeL = sym([]);
                end
                H_hr = HR(NUM_WAN,vectorL,'HnumL',HnumL,'HcoeL',HcoeL,'SnumL',SnumL,'ScoeL',ScoeL,'Type',Type,'overlap',true,'sym','false');
                H_hr.num = true;
                H_hr.coe = false;
                H_hr.vectorL_overlap = vectorL_overlap;
                H_hr.overlap = true;
            else
                [dataArray,NRPT_list,NRPTS,NUM_WAN]=HR.hrdat_read(filename);
                if strcmp(Type ,'mat')
                    %%
                    Vec_Fir = dataArray{:, 1};      % a1 direction
                    Vec_Sec = dataArray{:, 2};      % a2 direction
                    Vec_Thi = dataArray{:, 3};      % a3 direction
                    %Orb_fir = dataArray{:, 4};
                    %Orb_sec = dataArray{:, 5};
                    h_real = dataArray{:, 6};
                    h_imag = dataArray{:, 7};
                    %% Clear temporary variables
                    % read
                    V_f = reshape(Vec_Fir, NUM_WAN*NUM_WAN, NRPTS);
                    V_s = reshape(Vec_Sec, NUM_WAN*NUM_WAN, NRPTS);
                    V_t = reshape(Vec_Thi, NUM_WAN*NUM_WAN, NRPTS);
                    vectorL=[V_f(1,:)',V_s(1,:)',V_t(1,:)'];
                    HnumL_real = reshape(h_real, NUM_WAN,NUM_WAN, NRPTS);
                    HnumL_imag = reshape(h_imag, NUM_WAN,NUM_WAN, NRPTS);
                    HnumL = HnumL_real+1i*HnumL_imag;
                    if strcmp(Type,'n')
                        for i = 1 : NRPTS
                            HnumL(:,:,i) = HnumL(:,:,i)/NRPT_list(i);
                            %From here , get hints
                        end
                    else
                    end
                    HcoeL = sym([]);
                elseif strcmp(Type ,'list')
                    DATAARRAY = cell2mat(dataArray(1:7));
                    HOPARRAY = DATAARRAY(:,6)+1i*DATAARRAY(:,7);
                    % NRPT_lists = kron(NRPT_list,zeros(NUM_WAN*NUM_WAN,1));
                    HnumL_select = abs(HOPARRAY)>options.Accuracy;
                    vectorL = DATAARRAY(HnumL_select,1:5);
                    [~,~,original_label] = unique(vectorL(:,1:3),'rows');
                    HnumL = HOPARRAY(HnumL_select)./NRPT_list(original_label);
                    HcoeL = sym([]);
                end
                H_hr = HR(NUM_WAN,vectorL,'HnumL',HnumL,'HcoeL',HcoeL,'Type',Type,'sym',false);
                H_hr.num = true;H_hr.coe = false;
            end
            H_hr.Basis_num = H_hr.WAN_NUM;
        end
    end
    methods(Static,Hidden,Access= protected)
        function [dataArray,NRPT_list,NRPTS,NUM_WAN]=hrdat_read(filename)
            if nargin < 1
                filename = 'wannier90_hr.dat';
            end
            %% read Nbands and NRPTS
            delimiter = ' ';
            startRow = 2;
            endRow = 3;
            formatSpec = '%d%*s%[^\n\r]';
            fileID = fopen(filename,'r');
            dataArray = textscan(fileID, formatSpec, endRow-startRow+1, 'Delimiter', delimiter,...
                'MultipleDelimsAsOne', true, 'TextType', 'string',...
                'HeaderLines', startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
            NBR = dataArray{:, 1};
            NUM_WAN = NBR(1);
            NRPTS = NBR(2);
            NRPTS_num1=fix(double(NRPTS)/15);
            NRPTS_num2=mod(NRPTS,15);
            fclose(fileID);
            %% read NRTT_list
            fileID = fopen(filename,'r');
            startRow = 4;
            endRow = startRow+NRPTS_num1;% calculate here
            formatSpec = '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';
            dataArray = textscan(fileID, formatSpec,NRPTS_num1+1 , 'Delimiter', delimiter,...
                'MultipleDelimsAsOne', true, 'TextType', 'string', 'EmptyValue', 0,...
                'HeaderLines', startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
            NRPT_list = [dataArray{1:end-1}];
            NRPT_list = reshape(NRPT_list',NRPTS_num1*15+15,1);
            fclose(fileID);
            %% read hopping terms
            fileID = fopen(filename,'r');
            if NRPTS_num2==0
                startRow = endRow;
            else
                startRow = endRow+1;
            end
            formatSpec = '%f%f%f%f%f%f%f%*s%*s%*s%*s%*s%*s%*s%*s%[^\n\r]';
            %%
            dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
            fclose(fileID);
            %disp(h_imag);
            %clearvars filename startRow formatSpec fileID dataArray ans;
        end
    end
    % ---------------  indirect set property method -------------------
    methods
        % add a empty NRPT
        function H_hr = add_empty_one(H_hr,vector)
            if H_hr.vectorhopping
                H_hr.vectorL = [H_hr.vectorL;vector];
                nvector  = size(vector,1);
                H_hr.AvectorL = blkdiag(H_hr.AvectorL ,eye(nvector));
                H_hr.BvectorL = blkdiag(H_hr.BvectorL ,eye(nvector));
                H_hr.CvectorL = [blkdiag(H_hr.CvectorL(1:end/2,:),1*eye(nvector));...
                    blkdiag(H_hr.CvectorL(end/2+1:end,:),1*eye(nvector))];
                return;
            end
            for i = 1:size(vector,1)
                vector_single = int32(vector(i,:));
                try
                    if (ismember(vector_single,H_hr.vectorL,'rows') && ~H_hr.overlap) || ...
                            (ismember(vector_single,H_hr.vectorL,'rows') ...
                            && ismember(vector_single,H_hr.vectorL_overlap,'rows') && H_hr.overlap)
                        continue;
                    end
                catch
                    % ugly
                end
                NRPTS_new = H_hr.NRPTS +1;
                %             H_hr.NRPTS = NRPTS_new;
                H_hr.vectorL(NRPTS_new,:) = (vector_single);
                if strcmp(H_hr.Type ,'mat')
                    if H_hr.coe
                        H_hr.HcoeL(:,:,NRPTS_new) = sym(zeros(H_hr.WAN_NUM));
                    end
                    if H_hr.num
                        H_hr.HnumL(:,:,NRPTS_new) = zeros(H_hr.WAN_NUM);
                    end
                    if H_hr.overlap
                        H_hr.vectorL_overlap(NRPTS_new,:) = vector_single;
                        if H_hr.coe
                            H_hr.ScoeL(:,:,NRPTS_new) = sym(zeros(H_hr.WAN_NUM));
                        end
                        if H_hr.num
                            H_hr.SnumL(:,:,NRPTS_new) = zeros(H_hr.WAN_NUM);
                        end
                    end
                elseif strcmp(H_hr.Type ,'sparse')
                    if H_hr.coe
                        H_hr.HcoeL{NRPTS_new} = sym(sparse(H_hr.WAN_NUM));
                    end
                    if H_hr.num
                        H_hr.HnumL{NRPTS_new} = sparse(H_hr.WAN_NUM);
                    end
                    if H_hr.overlap
                        H_hr.vectorL_overlap(NRPTS_new,:) = vector_single;
                        if H_hr.coe
                            H_hr.ScoeL{NRPTS_new} = sym(sparse(H_hr.WAN_NUM));
                        end
                        if H_hr.num
                            H_hr.SnumL{NRPTS_new} = sparse(H_hr.WAN_NUM);
                        end
                    end
                elseif strcmp(H_hr.Type ,'list')
                    H_hr.HcoeL(NRPTS_new,1) = sym(0);
                    H_hr.HnumL(NRPTS_new,1) = 0;
                    if H_hr.overlap
                        H_hr.ScoeL(NRPTS_new,1) = sym(0);
                        H_hr.SnumL(NRPTS_new,1) = 0;
                    end
                end
            end
        end
        function H_hr = expand_empty_one(H_hr,orbOne,QuantumOne,elementOne)
            if nargin < 2
                orbOne = [0,0,0];
            end
            if nargin < 3
                QuantumOne = [1,0,0,1];
            end
            if nargin < 4
                elementOne = 1;
            end
            H_hr.orbL = [H_hr.orbL;orbOne];
            H_hr.quantumL = [H_hr.quantumL;QuantumOne];
            H_hr.elementL = [H_hr.elementL;elementOne];
            if strcmp(H_hr.Type ,'list')
                NRPTS_new = H_hr.NRPTS +size(orbOne,1);
                H_hr.HcoeL(NRPTS_new,1) = sym(0);
                H_hr.HnumL(NRPTS_new,1) = 0;
            else
                WANNUM = H_hr.WAN_NUM+size(orbOne,1);
                H_hr.HcoeL(WANNUM,WANNUM,:) = sym(0);
                H_hr.HnumL(WANNUM,WANNUM,:) = 0;
            end
        end
        % set onsite needed
        % set hop
        function H_hr = set_hop(H_hr,amp,hi,hj,vector_list,mode)
            % -------- nargin -------------
            if nargin <6
                mode = 'set';
            end
            % -------- init -------------
            WAN_NUM_half = H_hr.WAN_NUM/2;
            [n_vector,~] = size(vector_list);
            %V = H_hr.vectorL;
            if strcmp(mode,'sym')||strcmp(mode,'symadd')
                H_hr.num = false;
                H_hr.coe = true;
            end
            for i = 1:n_vector
                vector = vector_list(i,:);
                %[~,seq]=ismember(vector,V,'rows');
                %% test speed up
                if length(hi) == length(hj) && length(hi)>1 && strcmp(H_hr.Type,'mat')
                    WANNUM = H_hr.WAN_NUM;
                    switch mode
                        case {'set','add'}
                            amp_mat =zeros(WANNUM);
                        case {'sym','symadd'}
                            amp_mat =sym(zeros(WANNUM));
                    end
                    amp_mat(sub2ind([WANNUM,WANNUM],hi,hj)) = amp;
                    H_hr = H_hr.set_hop_mat(amp_mat,vector,mode);
                elseif  length(hi) == length(hj) && strcmp(H_hr.Type,'sparse') % && length(hi)>1
                    WANNUM = H_hr.WAN_NUM;
                    amp_mat = sparse(hi,hj,amp,WANNUM,WANNUM);
                    H_hr = H_hr.set_hop_mat(amp_mat,vector,mode);
                elseif length(hi) == length(hj) && length(hi)>1 && strcmp(H_hr.Type,'list')
                    for j = 1:length(hi)
                        H_hr = H_hr.set_hop_single(amp(j),hi(j),hj(j),vector,mode);
                    end
                else
                    length_amp = length(amp);
                    switch length_amp
                        case 1
                            %                     % single mode
                            H_hr = H_hr.set_hop_single(amp,hi,hj,vector,mode);
                        case 2
                            %                     % spin mode
                            H_hr = H_hr.set_hop_single(amp(1,1),hi,hj,vector,mode);
                            H_hr = H_hr.set_hop_single(amp(1,2),hi,hj+WAN_NUM_half,vector,mode);
                            H_hr = H_hr.set_hop_single(amp(2,1),hi+WAN_NUM_half,hj,vector,mode);
                            H_hr = H_hr.set_hop_single(amp(2,2),hi+WAN_NUM_half,hj+WAN_NUM_half,vector,mode);
                        case 4
                            % four-band mode
                            disp('if needed');
                        case 8
                            % eight-band mode
                            disp('if needed');
                        otherwise
                            H_hr = H_hr.set_hop_mat(amp,vector,mode);
                    end
                end
            end
        end
        % set hop mat
        function H_hr = set_hop_mat(H_hr,amp,vector,mode)
            V = H_hr.vectorL;
            switch H_hr.Type
                case 'list'
                    sizeamp = size(amp);
                    for n = 1:numel(amp)
                        [hi,hj] = ind2sub(sizeamp,n);
                        if ~isequal(sym(amp(n)),sym(0))
                            H_hr = set_hop_single(H_hr,amp(n),hi,hj,vector,mode);
                        end
                    end
                case 'sparse'
                    [~,seq]=ismember(int32(vector),V,'rows');
                    % for new block
                    if seq == 0
                        seq = H_hr.NRPTS +1;
                        H_hr = H_hr.add_empty_one(vector);
                    end
                    switch mode
                        case 'set'
                            %zeros_num_mat = sparse(H_hr.WAN_NUM);
                            %if H_hr.HnumL{seq} ~= zeros_num_mat
                            %    warning('May be you should use add mode on this NRPT');
                            %    fprintf('%d %d %d %d %d \n',vector(1),vector(2),vector(3));
                            %end
                            H_hr.HnumL{seq} = amp ;
                        case 'add'
                            H_hr.HnumL{seq} = H_hr.HnumL{seq} + amp  ;
                        case 'sym'
                            error('not be implemented yet');
                            zeros_coe_mat = sym(zeros(H_hr.WAN_NUM));
                            if H_hr.HcoeL(:,:,seq) ~= zeros_coe_mat
                                warning('May be you should use symadd mode on this NRPT ');
                                fprintf('%d %d %d %d %d \n',vector(1),vector(2),vector(3));
                            end
                            H_hr.HcoeL(:,:,seq) =  amp  ;
                        case 'symadd'
                            error('not be implemented yet');
                            H_hr.HcoeL(:,:,seq) = H_hr.HcoeL(:,:,seq) + amp  ;
                    end
                otherwise
                    [~,seq]=ismember(int32(vector),V,'rows');
                    % for new block
                    if seq == 0
                        seq = H_hr.NRPTS +1;
                        H_hr = H_hr.add_empty_one(vector);
                    end
                    switch mode
                        case 'set'
                            zeros_num_mat = zeros(H_hr.WAN_NUM);
                            if H_hr.HnumL(:,:,seq) ~= zeros_num_mat
                                warning('May be you should use add mode on this NRPT');
                                fprintf('%d %d %d %d %d \n',vector(1),vector(2),vector(3));
                            end
                            H_hr.HnumL(:,:,seq) = amp ;
                        case 'add'
                            H_hr.HnumL(:,:,seq) = H_hr.HnumL(:,:,seq) + amp  ;
                        case 'sym'
                            zeros_coe_mat = sym(zeros(H_hr.WAN_NUM));
                            if H_hr.HcoeL(:,:,seq) ~= zeros_coe_mat
                                warning('May be you should use symadd mode on this NRPT ');
                                fprintf('%d %d %d %d %d \n',vector(1),vector(2),vector(3));
                            end
                            H_hr.HcoeL(:,:,seq) =  amp  ;
                        case 'symadd'
                            H_hr.HcoeL(:,:,seq) = H_hr.HcoeL(:,:,seq) + amp  ;
                    end
            end
        end
        % set hop single
        function H_hr = set_hop_single(H_hr,amp,hi,hj,vector,mode)
            V = H_hr.vectorL;
            if strcmp(H_hr.Type,'list')
                vector = [vector,hi,hj];
            else
            end
            if isempty(V)
                seq = H_hr.NRPTS +1;
                H_hr = H_hr.add_empty_one(vector);
            else
                [~,seq]=ismember(int32(vector),V,'rows');
                % for new block
                if seq == 0
                    seq = H_hr.NRPTS +1;
                    H_hr = H_hr.add_empty_one(vector);
                end
            end
            switch H_hr.Type
                case 'list'
                    switch mode
                        case 'set'
                            if H_hr.HnumL(seq) ~= 0
                                warning('May be you should use add mode on this NRPT and hi hj');
                                fprintf('%d %d %d %d %d \n',vector(1),vector(2),vector(3),vector(4),vector(5));
                            end
                            H_hr.HnumL(seq) = amp ;
                        case 'add'
                            H_hr.HnumL(seq) = H_hr.HnumL(seq) + amp  ;
                        case 'sym'
                            if H_hr.HcoeL(seq) ~= sym(0)
                                warning('May be you should use symadd mode on this NRPT and hi hj');
                                fprintf('%d %d %d %d %d \n',vector(1),vector(2),vector(3),vector(4),vector(5));
                            end
                            H_hr.HcoeL(seq) =  amp  ;
                        case 'symadd'
                            H_hr.HcoeL(seq) = H_hr.HcoeL(seq) + amp  ;
                    end
                case 'sparse'
                    switch mode
                        case 'set'
                            if H_hr.HnumL{seq}(hi,hj) ~= 0
                                warning('May be you should use add mode on this NRPT and hi hj');
                                fprintf('%d %d %d %d %d \n',vector(1),vector(2),vector(3),hi,hj);
                            end
                            H_hr.HnumL{seq}(hi,hj) = amp ;
                        case 'add'
                            H_hr.HnumL{seq}(hi,hj) = H_hr.HnumL(hi,hj,seq) + amp  ;
                        case 'sym'
                            if H_hr.HcoeL{seq}(hi,hj) ~= sym(0)
                                warning('May be you should use symadd mode on this NRPT and hi hj');
                                fprintf('%d %d %d %d %d \n',vector(1),vector(2),vector(3),hi,hj);
                            end
                            H_hr.HcoeL{seq}(hi,hj) =  amp  ;
                        case 'symadd'
                            H_hr.HcoeL{seq}(hi,hj) = H_hr.HcoeL{seq}(hi,hj) + amp  ;
                    end
                otherwise
                    switch mode
                        case 'set'
                            if H_hr.HnumL(hi,hj,seq) ~= 0
                                warning('May be you should use add mode on this NRPT and hi hj');
                                fprintf('%d %d %d %d %d \n',vector(1),vector(2),vector(3),hi,hj);
                            end
                            H_hr.HnumL(hi,hj,seq) = amp ;
                        case 'add'
                            H_hr.HnumL(hi,hj,seq) = H_hr.HnumL(hi,hj,seq) + amp  ;
                        case 'sym'
                            if H_hr.HcoeL(hi,hj,seq) ~= sym(0)
                                warning('May be you should use symadd mode on this NRPT and hi hj');
                                fprintf('%d %d %d %d %d \n',vector(1),vector(2),vector(3),hi,hj);
                            end
                            H_hr.HcoeL(hi,hj,seq) =  amp  ;
                        case 'symadd'
                            H_hr.HcoeL(hi,hj,seq) = H_hr.HcoeL(hi,hj,seq) + amp  ;
                    end
            end
        end
        % set overlap
        function H_hr = set_overlap(H_hr,amp,si,sj,vector_list,mode)
            % -------- nargin -------------
            if nargin <6
                mode = 'set';
            end
            % -------- init -------------
            WAN_NUM_half = H_hr.WAN_NUM/2;
            [n_vector,~] = size(vector_list);
            %V = H_hr.vectorL;
            for i = 1:n_vector
                vector = vector_list(i,:);
                %[~,seq]=ismember(vector,V,'rows');
                %% test amp
                if length(si) == length(sj) && length(si)>1 && strcmp(H_hr.Type,'mat')
                    WANNUM = H_hr.WAN_NUM;
                    switch mode
                        case {'set','add'}
                            amp_mat =zeros(WANNUM);
                        case {'sym','symadd'}
                            amp_mat =sym(zeros(WANNUM));
                    end
                    amp_mat(sub2ind([WANNUM,WANNUM],si,sj)) = amp;
                    H_hr = H_hr.set_overlap_mat(amp_mat,vector,mode);
                elseif length(si) == length(sj) && length(si)>1 && strcmp(H_hr.Type,'list')
                    for j = 1:length(si)
                        H_hr = H_hr.set_overlap_single(amp(j),si(j),sj(j),vector,mode);
                    end
                else
                    length_amp = length(amp);
                    switch length_amp
                        case 1
                            %                     % single mode
                            H_hr = H_hr.set_overlap_single(amp,si,sj,vector,mode);
                        case 2
                            %                     % spin mode
                            H_hr = H_hr.set_overlap_single(amp(1,1),si,sj,vector,mode);
                            H_hr = H_hr.set_overlap_single(amp(1,2),si,sj+WAN_NUM_half,vector,mode);
                            H_hr = H_hr.set_overlap_single(amp(2,1),si+WAN_NUM_half,sj,vector,mode);
                            H_hr = H_hr.set_overlap_single(amp(2,2),si+WAN_NUM_half,sj+WAN_NUM_half,vector,mode);
                        case 4
                            % four-band mode
                            disp('if needed');
                        case 8
                            % eight-band mode
                            disp('if needed');
                        otherwise
                            H_hr = H_hr.set_overlap_mat(amp,vector,mode);
                    end
                end
            end
        end
        % set overlap mat
        function H_hr = set_overlap_mat(H_hr,amp,vector,mode)
            V = H_hr.vectorL_overlap;
            switch H_hr.Type
                case 'list'
                    sizeamp = size(amp);
                    for n = 1:numel(amp)
                        [si,sj] = ind2sub(sizeamp,n);
                        if ~isequal(sym(amp(n)),sym(0))
                            H_hr = set_overlap_single(H_hr,amp(n),si,sj,vector,mode);
                        end
                    end
                otherwise
                    [~,seq]=ismember(int32(vector),V,'rows');
                    % for new block
                    if seq == 0
                        seq = H_hr.NRPTS +1;
                        H_hr = H_hr.add_empty_one(vector);
                    end
                    switch mode
                        case 'set'
                            zeros_num_mat = zeros(H_hr.WAN_NUM);
                            if H_hr.HnumL(:,:,seq) ~= zeros_num_mat
                                warning('May be you should use add mode on this NRPT');
                                fprintf('%d %d %d %d %d \n',vector(1),vector(2),vector(3));
                            end
                            H_hr.SnumL(:,:,seq) = amp ;
                        case 'add'
                            H_hr.SnumL(:,:,seq) = H_hr.SnumL(:,:,seq) + amp  ;
                        case 'sym'
                            zeros_coe_mat = sym(zeros(H_hr.WAN_NUM));
                            if H_hr.ScoeL(:,:,seq) ~= zeros_coe_mat
                                warning('May be you should use symadd mode on this NRPT ');
                                fprintf('%d %d %d %d %d \n',vector(1),vector(2),vector(3));
                            end
                            H_hr.ScoeL(:,:,seq) =  amp  ;
                        case 'symadd'
                            H_hr.ScoeL(:,:,seq) = H_hr.ScoeL(:,:,seq) + amp  ;
                    end
            end
        end
        % set overlap single
        function H_hr = set_overlap_single(H_hr,amp,si,sj,vector,mode)
            V = H_hr.vectorL_overlap;
            if strcmp(H_hr.Type,'list')
                vector = [vector,si,sj];
            else
                
            end
            if isempty(V)
                seq = H_hr.NRPTS +1;
                H_hr = H_hr.add_empty_one(vector);
            else
                [~,seq]=ismember(int32(vector),V,'rows');
                % for new block
                if seq == 0
                    seq = H_hr.NRPTS +1;
                    H_hr = H_hr.add_empty_one(vector);
                end
            end
            switch H_hr.Type
                case 'list'
                    switch mode
                        case 'set'
                            if H_hr.SnumL(seq) ~= 0
                                warning('May be you should use add mode on this NRPT and hi hj');
                                fprintf('%d %d %d %d %d \n',vector(1),vector(2),vector(3),vector(4),vector(5));
                            end
                            H_hr.SnumL(seq) = amp ;
                        case 'add'
                            H_hr.SnumL(seq) = H_hr.SnumL(seq) + amp  ;
                        case 'sym'
                            if H_hr.ScoeL(seq) ~= sym(0)
                                warning('May be you should use symadd mode on this NRPT and hi hj');
                                fprintf('%d %d %d %d %d \n',vector(1),vector(2),vector(3),vector(4),vector(5));
                            end
                            H_hr.ScoeL(seq) =  amp  ;
                        case 'symadd'
                            H_hr.ScoeL(seq) = H_hr.ScoeL(seq) + amp  ;
                    end
                otherwise
                    switch mode
                        case 'set'
                            if H_hr.SnumL(si,sj,seq) ~= 0
                                warning('May be you should use add mode on this NRPT and si sj');
                                fprintf('%d %d %d %d %d \n',vector(1),vector(2),vector(3),si,sj);
                            end
                            H_hr.SnumL(si,sj,seq) = amp ;
                        case 'add'
                            H_hr.SnumL(si,sj,seq) = H_hr.SnumL(si,sj,seq) + amp  ;
                        case 'sym'
                            if H_hr.ScoeL(si,sj,seq) ~= sym(0)
                                warning('May be you should use symadd mode on this NRPT and si sj');
                                fprintf('%d %d %d %d %d \n',vector(1),vector(2),vector(3),si,sj);
                            end
                            H_hr.ScoeL(si,sj,seq) = amp  ;
                        case 'symadd'
                            H_hr.ScoeL(si,sj,seq) = H_hr.ScoeL(si,sj,seq) + amp  ;
                    end
            end
        end
    end
    %% set
    % ----------------  set property method --------------------
    methods
        function H_hr = set.NRPTS( H_hr, ~ )
            fprintf ("%s%d\n", "NRPTS is ", H_hr.NRPTS);
            error ( "You dont need set NRPTS explicitly.");
        end
        function H_hr = set.WAN_NUM(H_hr,~)
            fprintf ("%s%d\n", "WAN_NUM is ", H_hr.WAN_NUM);
            error ( "You dont need set WAN_NUM explicitly.");
        end
    end
    %% get
    % ----------------  get property method --------------------
    methods
        function NRPTS = get.NRPTS(H_hr)
            NRPTS = size(H_hr.vectorL,1);
        end
        function WAN_NUM = get.WAN_NUM(H_hr)
            if strcmp(H_hr.Type,'sparse')
                WAN_NUM = size(H_hr.HnumL{1},1);
            elseif strcmp(H_hr.Type,'list')
                try
                    [Sparse_vector,~,~] = unique(H_hr.vectorL(:,4:5),'rows');
                    WAN_NUM = double(max(Sparse_vector(:)));
                catch
                    WAN_NUM = 0;
                end
            else
                %WAN_NUM = max(size(H_hr.HnumL,1),size(H_hr.HcoeL,1));
                if H_hr.num && ~H_hr.coe
                    WAN_NUM = size(H_hr.HnumL,2);
                elseif ~H_hr.num && H_hr.coe
                    WAN_NUM = size(H_hr.HcoeL,2);
                elseif H_hr.num && H_hr.coe
                    WAN_NUM = max(size(H_hr.HnumL,1),size(H_hr.HcoeL,1));
                else
                    WAN_NUM = H_hr.Basis_num;
                end
            end
        end
        %         function Basis_num = get.Basis_num(H_hr)
        %             Basis_num = H_hr.WAN_NUM;
        %         end
        function Line_000 = get.Line_000(H_hr)
            vector = [0,0,0];
            [~,Line_000]=ismember(int32(vector),H_hr.vectorL(:,1:3),'rows');
        end
        function homecell = get.homecell(H_hr)
            if strcmp(H_hr.Type,'list')
                homecell = H_hr.list([0,0,0]);
            else
                if H_hr.coe
                    homecell = H_hr.HcoeL(:,:,H_hr.Line_000);
                else
                    homecell = H_hr.HnumL(:,:,H_hr.Line_000);
                end
            end
        end
        function Type = type(H_hr)
            Type = H_hr.Type;
        end
        %         function AvectorL = get.AvectorL(H_hr)
        %             AvectorL = H_hr.HnumL;
        %         end
        %         function BvectorL = get.BvectorL(H_hr)
        %             BvectorL = H_hr.HcoeL;
        %         end
    end
    %% unknown(Not Classify yet)
    methods
        %
        function vectorSeq = Getvector(H_hr,vector)
            [~,vectorSeq]=ismember(int32(vector),H_hr.vectorL,'rows');
        end
        % autohermi
        function H_hr = autohermi(H_hr,mode,options)
            arguments
                H_hr HR;
                mode =  'sym';
                options.enforce_list = true;
            end
            if nargin<2
                mode = 'sym';
            end
            if options.enforce_list
                if ~strcmp(H_hr.type, 'list')
                    H_hr = H_hr.rewrite();
                    giveback  = true;
                else
                    giveback  = false;
                end
            else
                giveback  = true;
            end
            H_hr_tmp = H_hr;
            switch mode
                case 'sym'
                    fprintf('For sym Hr, The requirement is strick!');
                    fprintf('The sym vasiable is real or ');
                    fprintf('in the original Hr, use conj()\n');
                    switch H_hr.type
                        case 'mat'
                            for i = 1:H_hr.NRPTS
                                %Duality_vector_dist();
                                vector_tmp = H_hr.vectorL(i,:);
                                vector_tmp_oppo = -vector_tmp;
                                [~,j]=ismember(vector_tmp_oppo,H_hr_tmp.vectorL,'rows');
                                fprintf('check %3d th : NRPT,\n the vector is %3d %3d %3d,\n the opposite vector is %3d,%3d,%3d the %3d th NRPT\n',i,...
                                    vector_tmp(1),vector_tmp(2),vector_tmp(3),...
                                    vector_tmp_oppo(1),vector_tmp_oppo(2),vector_tmp_oppo(3),j);
                                if i == j
                                    if ~isequal(H_hr.HcoeL(:,:,i) ,H_hr_tmp.HcoeL(:,:,j)')
                                        fprintf('The homecell hamilton is not hermi, never mind, we will hermi it enforcely!\n');
                                        disp(H_hr.HcoeL(:,:,i));
                                        H_hr_tmp = H_hr_tmp.set_hop_mat(H_hr.HcoeL(:,:,i)'/2+H_hr_tmp.HcoeL(:,:,j)/2,[0,0,0],'sym');
                                        fprintf('change into\n');
                                        disp(H_hr_tmp.HcoeL(:,:,i));
                                    end
                                    continue;
                                end
                                if j == 0
                                    fprintf('The opposite vector hamilton does not exist, build it!\n');
                                    H_hr_tmp = H_hr_tmp.set_hop_mat(H_hr.HcoeL(:,:,i)',-vector_tmp,'sym');
                                    disp(H_hr_tmp.HcoeL(:,:,i)');
                                    continue;
                                elseif ~isequal(H_hr.HcoeL(:,:,i) ,H_hr_tmp.HcoeL(:,:,j)')
                                    fprintf('The opposite vector hamilton is not hermi, replace it by strong mat! \n');
                                    N1 = nnz(H_hr.HcoeL(:,:,i));
                                    N2 = nnz(H_hr_tmp.HcoeL(:,:,j)');
                                    %disp([N1,N2]);
                                    if N1 >= N2
                                        fprintf('The %3d th NRPT is stonger!\n',i);
                                        disp(H_hr.HcoeL(:,:,i));
                                        H_hr_tmp = H_hr_tmp.set_hop_mat(H_hr.HcoeL(:,:,i)',-vector_tmp,'sym');
                                    elseif N1 < N2
                                        fprintf('The %3d th NRPT is stonger!\n',j);
                                        disp(H_hr_tmp.HcoeL(:,:,j)');
                                        H_hr_tmp = H_hr_tmp.set_hop_mat(H_hr_tmp.HcoeL(:,:,j)',vector_tmp,'sym');
                                    else
                                    end
                                elseif isequal(H_hr.HcoeL(:,:,i),H_hr_tmp.HcoeL(:,:,j)')
                                    disp('hermi test pasts');
                                else
                                    warning('!!!!!');
                                end
                            end
                        case 'list'
                            for i = 1:H_hr.NRPTS
                                % Duality_vector_dist();
                                vector_tmp = H_hr.vectorL(i,:);
                                vector_tmp_oppo(1:3) = -vector_tmp(1:3);
                                vector_tmp_oppo(4) = vector_tmp(5);
                                vector_tmp_oppo(5) = vector_tmp(4);
                                [~,j]=ismember(vector_tmp_oppo,H_hr_tmp.vectorL,'rows');
                                % homecell
                                if i == j
                                    if ~isequal(H_hr.HcoeL(i) ,H_hr_tmp.HcoeL(j)')
                                        if H_hr.HcoeL(i)== sym(0)
                                            fprintf('The testing homecell hamilton %3d,%3d,%3d[i:%3d j:%3d] does not exist, build it with : %s!\n', ...
                                                vector_tmp(1),vector_tmp(2),vector_tmp(3),vector_tmp(4),vector_tmp(5), ...
                                                string(H_hr.HcoeL(j)'));
                                            H_hr_tmp = H_hr_tmp.set_hop(H_hr.HcoeL(j)',vector_tmp(4),vector_tmp(5),vector_tmp(1:3),'sym');
                                        elseif H_hr.HcoeL(j)== sym(0)
                                            fprintf('The opposite homecell hamilton %3d,%3d,%3d[i:%3d j:%3d] does not exist, build it : %s!\n', ...
                                                vector_tmp_oppo(1),vector_tmp_oppo(2),vector_tmp_oppo(3), vector_tmp_oppo(4),vector_tmp_oppo(5), ...
                                                string(H_hr.HcoeL(i)'));
                                            H_hr_tmp = H_hr_tmp.set_hop(H_hr.HcoeL(i)',vector_tmp_oppo(4),vector_tmp_oppo(5),vector_tmp_oppo(1:3),'sym');
                                        else
                                            fprintf('The homecell hamilton is not hermi, never mind, we will not hermi it enforcely!\n');
                                        end
                                    end
                                    continue;
                                end
                                %
                                if j == 0
                                    fprintf('The opposite vector hamilton %3d,%3d,%3d[i:%3d j:%3d] does not exist, build it : %s!\n', ...
                                        vector_tmp_oppo(1),vector_tmp_oppo(2),vector_tmp_oppo(3), vector_tmp_oppo(4),vector_tmp_oppo(5), ...
                                        string(H_hr.HcoeL(i)'));
                                    H_hr_tmp = H_hr_tmp.set_hop(H_hr.HcoeL(i)',vector_tmp_oppo(4),vector_tmp_oppo(5),vector_tmp_oppo(1:3),'sym');
                                    continue;
                                end
                                %
                                if ~isequal(H_hr.HcoeL(i) ,H_hr_tmp.HcoeL(j)')
                                    if H_hr.HcoeL(i)== sym(0)
                                        fprintf('The testing vector hamilton %3d,%3d,%3d[i:%3d j:%3d] does not exist, build it with : %s!\n', ...
                                            vector_tmp(1),vector_tmp(2),vector_tmp(3),vector_tmp(4),vector_tmp(5), ...
                                            string(H_hr.HcoeL(j)'));
                                        H_hr_tmp = H_hr_tmp.set_hop(H_hr.HcoeL(j)',vector_tmp(4),vector_tmp(5),vector_tmp(1:3),'sym');
                                    elseif H_hr.HcoeL(j)== sym(0)
                                        fprintf('The opposite vector hamilton %3d,%3d,%3d[i:%3d j:%3d] does not exist, build it : %s!\n', ...
                                            vector_tmp_oppo(1),vector_tmp_oppo(2),vector_tmp_oppo(3), vector_tmp_oppo(4),vector_tmp_oppo(5), ...
                                            string(H_hr.HcoeL(i)'));
                                        H_hr_tmp = H_hr_tmp.set_hop(H_hr.HcoeL(i)',vector_tmp_oppo(4),vector_tmp_oppo(5),vector_tmp_oppo(1:3),'sym');
                                    else
                                        fprintf(['check on the %3d th NRPT,\n' ...
                                            '---------the vector is %3d %3d %3d[i:%d j:%d],\n' ...
                                            'the opposite vector is %3d,%3d,%3d[i:%d j:%d],\n' ...
                                            'find in the %3d th NRPT is not hermi, average them\n'], ...
                                            i,...
                                            vector_tmp(1),vector_tmp(2),vector_tmp(3),vector_tmp(4),vector_tmp(5),...
                                            vector_tmp_oppo(1),vector_tmp_oppo(2),vector_tmp_oppo(3),vector_tmp_oppo(4),vector_tmp_oppo(5), ...
                                            j);
                                        tmpsym = (H_hr.HcoeL(i)+H_hr.HcoeL(j)')/2;
                                        H_hr_tmp = H_hr_tmp.set_hop(tmpsym,vector_tmp(4),vector_tmp(5),vector_tmp(1:3),'sym');
                                        H_hr_tmp = H_hr_tmp.set_hop(tmpsym',vector_tmp_oppo(4),vector_tmp_oppo(5),vector_tmp_oppo(1:3)','sym');
                                    end
                                end
                            end
                        otherwise
                            
                    end
                case 'num'
                    fprintf('For num Hr, The requirement is less strick!');
                    for i = 1:H_hr.NRPTS
                        %Duality_vector_dist();
                        vector_tmp = H_hr.vectorL(i,:);
                        vector_tmp_oppo = -vector_tmp;
                        [~,j]=ismember(vector_tmp_oppo,H_hr_tmp.vectorL,'rows');
                        fprintf('check %d th : NRPT,\n the vector is %d %d %d,\n the opposite vector is %d,%d,%d the %d th NRPT\n',i,...
                            vector_tmp(1),vector_tmp(2),vector_tmp(3),...
                            vector_tmp_oppo(1),vector_tmp_oppo(2),vector_tmp_oppo(3),j);
                        if i == j
                            if ~isequal(H_hr.HnumL(:,:,i) ,H_hr_tmp.HnumL(:,:,j)')
                                fprintf('The homecell hamilton is not hermi, never mind, we will hermi it enforcely!\n');
                                disp(H_hr.HnumL(:,:,i));
                                H_hr_tmp = H_hr_tmp.set_hop_mat(H_hr.HnumL(:,:,i)'/2+H_hr_tmp.HnumL(:,:,j)/2,[0,0,0],'set');
                                fprintf('change into\n');
                                disp(H_hr_tmp.HnumL(:,:,i));
                            end
                            continue;
                        end
                        if j == 0
                            fprintf('The opposite vector hamilton does not exist, build it!\n');
                            H_hr_tmp = H_hr_tmp.set_hop_mat(H_hr.HnumL(:,:,i)',-vector_tmp,'set');
                            disp(H_hr_tmp.HnumL(:,:,i)');
                            continue;
                        elseif ~isequal(H_hr.HnumL(:,:,i) ,H_hr_tmp.HnumL(:,:,j)')
                            fprintf('The opposite vector hamilton is not hermi, replace it by strong mat! \n');
                            N1 = nnz(H_hr.HnumL(:,:,i));
                            N2 = nnz(H_hr_tmp.HnumL(:,:,j)');
                            %disp([N1,N2]);
                            if N1 >= N2
                                fprintf('The %d th NRPT is stonger!\n',i);
                                disp(H_hr.HnumL(:,:,i));
                                H_hr_tmp = H_hr_tmp.set_hop_mat(H_hr.HnumL(:,:,i)',-vector_tmp,'set');
                            elseif N1 < N2
                                fprintf('The %d th NRPT is stonger!\n',j);
                                disp(H_hr_tmp.HnumL(:,:,j)');
                                H_hr_tmp = H_hr_tmp.set_hop_mat(H_hr_tmp.HnumL(:,:,j)',vector_tmp,'set');
                            else
                            end
                        elseif isequal(H_hr.HnumL(:,:,i),H_hr_tmp.HnumL(:,:,j)')
                            disp('hermi test pasts');
                        else
                            warning('!!!!!');
                        end
                    end
            end
            if giveback
                H_hr = H_hr_tmp.rewind();
            else
                H_hr = H_hr_tmp;
            end
        end
    end
    %% overload
    methods
        % overload +
        function H_hr = plus(A,B)
            if isa(A,'HR') && isa(B,'HR')
                % --------vectorhopping--------
                H_hr1 = A;
                H_hr2 = B;
                % ---------check----------
                if H_hr1.WAN_NUM ~= H_hr2.WAN_NUM
                    error('WAN_NUM different');
                end
                
                if ~strcmp(A.Type,'list') && ~strcmp(B.Type,'list')
                    H_hr = H_hr1;
                    for i = 1:H_hr2.NRPTS
                        vector = H_hr2.vectorL(i,:);
                        if H_hr.coe
                            amp_sym = H_hr2.HcoeL(:,:,i);
                            H_hr = H_hr.set_hop_mat(amp_sym,vector,'symadd');
                        end
                        if H_hr.num
                            amp = H_hr2.HnumL(:,:,i);
                            H_hr = H_hr.set_hop_mat(amp,vector,'add');
                        end
                    end
                elseif strcmp(H_hr1.Type,'list') && strcmp(H_hr2.Type,'list')
                    [vectorList,~,~] = unique([H_hr1.vectorL;H_hr2.vectorL],'rows');
                    C = setdiff(vectorList,H_hr1.vectorL,'rows');
                    D = setdiff(vectorList,H_hr2.vectorL,'rows');
                    H_hr1 = H_hr1.add_empty_one(C);
                    H_hr2 = H_hr2.add_empty_one(D);
                    [~,seqA] = ismember(vectorList,H_hr1.vectorL,'rows');
                    [~,seqB] = ismember(vectorList,H_hr2.vectorL,'rows');
                    H_hr = H_hr1.reseq(':',seqA);
                    if H_hr.vectorhopping && B.vectorhopping
                        H_hr.AvectorL = H_hr.AvectorL + H_hr2.AvectorL(seqB,:);
                        H_hr.BvectorL = H_hr.BvectorL + H_hr2.BvectorL(seqB,:);
                        CL1 = H_hr2.CvectorL(1:end/2,:);
                        CL2 = H_hr2.CvectorL(end/2+1:end,:);
                        H_hr.CvectorL = H_hr.CvectorL + [CL1(seqB,:);CL2(seqB,:)];
                        return;
                    end
                    if H_hr.coe
                        H_hr.HcoeL = H_hr.HcoeL + H_hr2.HcoeL(seqB) ;
                    end
                    if H_hr.num
                        H_hr.HnumL = H_hr.HnumL + H_hr2.HnumL(seqB) ;
                    end
                elseif ~strcmp(H_hr1.Type,H_hr2.Type)
                    H_hr1 = H_hr1.rewrite();
                    H_hr2 = H_hr2.rewrite();
                    H_hr = H_hr1+H_hr2;
                    return;
                end
            elseif isa(A,'HR') && ~isa(B,'HR')
                H_hr1 = A;
                H_hr = H_hr1;
                %H_hr = H_hr.line_000_gen();
                i = H_hr.Line_000;
                %                 H_hr2 = B;
                % ---------check----------
                %                 if H_hr1.WAN_NUM ~= H_hr2.WAN_NUM
                %                     error('WAN_NUM different');
                %                 end
                if isa(B,'sym')
                    H_hr.HcoeL(:,:,i) = H_hr.HcoeL(:,:,i) + B;
                elseif isa(B,'Term')
                    H_hk_tmp = HK(H_hr.WAN_NUM,2);
                    H_hk_tmp.Rm = H_hr.Rm;
                    H_hk_tmp = H_hk_tmp + B;
                    H_hr = H_hr.plus(H_hk_tmp.kp2TB());
                elseif isa(B,'numeric')
                    H_hr.HnumL(:,:,i) = H_hr.HnumL(:,:,i) + B;
                else
                    
                end
            elseif ~isa(A,'HR') && isa(B,'HR')
                H_hr1 = B;
                H_hr = H_hr1;
                %                 H_hr2 = B;
                % ---------check----------
                %                 if H_hr1.WAN_NUM ~= H_hr2.WAN_NUM
                %                     error('WAN_NUM different');
                %                 end
                if isa(A,'sym')
                    for i = 1:H_hr.NRPTS
                        H_hr.HcoeL(:,:,i) = H_hr.HcoeL(:,:,i) + A;
                    end
                elseif isa(A,'numeric')
                    for i = 1:H_hr.NRPTS
                        H_hr.HnumL(:,:,i) = H_hr.HnumL(:,:,i) + A;
                    end
                else
                    
                end
            end
        end
        % overload -
        function H_hr = minus(A,B)
            if isa(A,'HR') && isa(B,'HR')
                H_hr1 = A;
                H_hr2 = B;
                % ---------check----------
                if H_hr1.WAN_NUM ~= H_hr2.WAN_NUM
                    error('WAN_NUM different');
                end
                H_hr = H_hr1;
                for i = 1:H_hr2.NRPTS
                    vector = H_hr2.vectorL(i,:);
                    amp = -H_hr2.HnumL(:,:,i);
                    amp_sym = -H_hr2.HcoeL(:,:,i);
                    H_hr = H_hr.set_hop_mat(amp,vector,'add');
                    H_hr = H_hr.set_hop_mat(amp_sym,vector,'symadd');
                end
            elseif isa(A,'HR') && ~isa(B,'HR')
                H_hr1 = A;
                H_hr = H_hr1;
                i = H_hr.Line_000;
                %                 H_hr2 = B;
                % ---------check----------
                %                 if H_hr1.WAN_NUM ~= H_hr2.WAN_NUM
                %                     error('WAN_NUM different');
                %                 end
                if isa(B,'sym')
                    H_hr.HcoeL(:,:,i) = H_hr.HcoeL(:,:,i) - B;
                elseif isa(B,'numeric')
                    H_hr.HnumL(:,:,i) = H_hr.HnumL(:,:,i) - B;
                else
                    
                end
            elseif ~isa(A,'HR') && isa(B,'HR')
                H_hr1 = B;
                H_hr = H_hr1;
                %                 H_hr2 = B;
                % ---------check----------
                %                 if H_hr1.WAN_NUM ~= H_hr2.WAN_NUM
                %                     error('WAN_NUM different');
                %                 end
                if isa(A,'sym')
                    for i = 1:H_hr.NRPTS
                        H_hr.HcoeL(:,:,i) = A + -H_hr.HcoeL(:,:,i) ;
                    end
                elseif isa(A,'numeric')
                    for i = 1:H_hr.NRPTS
                        H_hr.HnumL(:,:,i) = A + -H_hr.HnumL(:,:,i) ;
                    end
                else
                end
                
            end
        end
        % overload -self
        function H_hr = uminus(H_hr)
            for i = 1:H_hr.NRPTS
                H_hr.HcoeL(:,:,i) = -H_hr.HcoeL(:,:,i);
                H_hr.HnumL(:,:,i) = -H_hr.HnumL(:,:,i);
            end
        end
        % overload .*
        function H_hr = times(A,B)
            if isa(A,'HR') && isa(B,'HR')
                H_hr1 = A;
                H_hr2 = B;
                % ---------check----------
                if H_hr1.WAN_NUM ~= H_hr2.WAN_NUM
                    error('WAN_NUM different');
                end
                % ---------init-----------
                H_hr =  HR(H_hr2.WAN_NUM,...
                    unique([H_hr1.vectorL;H_hr2.vectorL],'rows'));
                %             zeros_num_mat = zeros(H_hr.WAN_NUM);
                %             zeros_coe_mat = sym(zeros(H_hr.WAN_NUM));
                for i = 1:H_hr.NRPTS
                    vector = H_hr.vectorL(i,:);
                    [~,seq1]=ismember(vector,H_hr1.vectorL,'rows');
                    [~,seq2]=ismember(vector,H_hr2.vectorL,'rows');
                    if seq1 ~= 0 && seq2 ~=0
                        amp = H_hr1.HnumL(:,:,seq1) .* H_hr2.HnumL(:,:,seq2);
                        amp_sym = H_hr1.HcoeL(:,:,seq1) .* H_hr2.HcoeL(:,:,seq2);
                        H_hr = H_hr.set_hop_mat(amp,vector,'set');
                        H_hr = H_hr.set_hop_mat(amp_sym,vector,'sym');
                    end
                end
                
            elseif isa(A,'HR') && ~isa(B,'HR')
                H_hr1 = A;
                H_hr = H_hr1;
                %                 H_hr2 = B;
                % ---------check----------
                if H_hr.WAN_NUM ~= length(B)
                    error('WAN_NUM different');
                end
                if isa(B,'sym')
                    for i = 1:H_hr.NRPTS
                        H_hr.HcoeL(:,:,i) = H_hr.HcoeL(:,:,i) .* B;
                    end
                elseif isa(B,'numeric')
                    for i = 1:H_hr.NRPTS
                        H_hr.HnumL(:,:,i) = H_hr.HnumL(:,:,i) .* B;
                    end
                else
                    
                end
                
            elseif ~isa(A,'HR') && isa(B,'HR')
                H_hr1 = B;
                H_hr = H_hr1;
                H_hr = H_hr.line_000_gen();
                i = H_hr.Line_000;
                % ---------check----------
                if H_hr.WAN_NUM ~= length(B)
                    error('WAN_NUM different');
                end
                if isa(A,'sym')
                    H_hr.HcoeL(:,:,i) = A .* H_hr.HcoeL(:,:,i) ;
                elseif isa(A,'numeric')
                    H_hr.HnumL(:,:,i) = A .* H_hr.HnumL(:,:,i) ;
                else
                end
            end
        end
        % overload /
        function H_hr = mrdivide(A,B)
            if isa(A,'HR') && isa(B,'HR')
                H_hr1 = A;
                H_hr2 = B;
                % ---------check----------
                if H_hr1.WAN_NUM ~= H_hr2.WAN_NUM
                    error('WAN_NUM different');
                end
                % ---------init-----------
                H_hr =  HR(H_hr2.WAN_NUM,...
                    unique([H_hr1.vectorL;H_hr2.vectorL],'rows'));
                %             zeros_num_mat = zeros(H_hr.WAN_NUM);
                %             zeros_coe_mat = sym(zeros(H_hr.WAN_NUM));
                for i = 1:H_hr.NRPTS
                    vector = H_hr.vectorL(i,:);
                    [~,seq1]=ismember(vector,H_hr1.vectorL,'rows');
                    [~,seq2]=ismember(vector,H_hr2.vectorL,'rows');
                    if seq1 ~= 0 && seq2 ~=0
                        amp = H_hr1.HnumL(:,:,seq1) / H_hr2.HnumL(:,:,seq2);
                        amp_sym = H_hr1.HcoeL(:,:,seq1) / H_hr2.HcoeL(:,:,seq2);
                        H_hr = H_hr.set_hop_mat(amp,vector,'set');
                        H_hr = H_hr.set_hop_mat(amp_sym,vector,'sym');
                    end
                end
                
            elseif isa(A,'HR') && ~isa(B,'HR')
                H_hr = A;
                
                %                 H_hr2 = B;
                % ---------check----------
                %                 if H_hr.WAN_NUM ~= length(B)
                %                     error('WAN_NUM different');
                %                 end
                if A.Type == 'list'
                    H_hr.HcoeL = H_hr.HcoeL/B;
                    H_hr.HnumL = H_hr.HnumL/B;
                else
                    if isa(B,'sym')
                        for i = 1:H_hr.NRPTS
                            H_hr.HcoeL(:,:,i) = H_hr.HcoeL(:,:,i) / B;
                        end
                    elseif isa(B,'numeric')
                        for i = 1:H_hr.NRPTS
                            H_hr.HnumL(:,:,i) = H_hr.HnumL(:,:,i) / B;
                        end
                    else
                        
                    end
                end
                
            elseif ~isa(A,'HR') && isa(B,'HR')
                H_hr1 = B;
                H_hr = H_hr1;
                H_hr = H_hr.line_000_gen();
                i = H_hr.Line_000;
                % ---------check----------
                if H_hr.WAN_NUM ~= length(B)
                    error('WAN_NUM different');
                end
                if isa(A,'sym')
                    H_hr.HcoeL(:,:,i) = A / H_hr.HcoeL(:,:,i) ;
                elseif isa(A,'numeric')
                    H_hr.HnumL(:,:,i) = A / H_hr.HnumL(:,:,i) ;
                else
                end
            end
        end
        % overload * mtimes(A,H_hr)
        function H_hr = premtimes(A,B)
            if isa(A,'HR') && ~isa(B,'HR')
                H_hr1 = A;
                H_hr = H_hr1;
                %                 H_hr2 = B;
                % ---------check----------
                if isa(B,'sym')
                    for i = 1:H_hr.NRPTS
                        H_hr.HcoeL(:,:,i) = B*H_hr.HcoeL(:,:,i);
                    end
                elseif isa(B,'numeric')
                    for i = 1:H_hr.NRPTS
                        H_hr.HnumL(:,:,i) = B*H_hr.HnumL(:,:,i) ;
                    end
                else
                    
                end
            end
        end
        function H_hr = mtimes(A,B)
            if isa(A,'HR') && isa(B,'HR')
                H_hr1 = A;
                H_hr2 = B;
                % ---------check----------
                if H_hr1.WAN_NUM ~= H_hr2.WAN_NUM
                    error('WAN_NUM different');
                end
                % ---------init-----------
                H_hr =  HR(H_hr2.WAN_NUM,...
                    unique([H_hr1.vectorL;H_hr2.vectorL],'rows'));
                %             zeros_num_mat = zeros(H_hr.WAN_NUM);
                %             zeros_coe_mat = sym(zeros(H_hr.WAN_NUM));
                for i = 1:H_hr.NRPTS
                    vector = H_hr.vectorL(i,:);
                    [~,seq1]=ismember(vector,H_hr1.vectorL,'rows');
                    [~,seq2]=ismember(vector,H_hr2.vectorL,'rows');
                    if seq1 ~= 0 && seq2 ~=0
                        amp = H_hr1.HnumL(:,:,seq1) * H_hr2.HnumL(:,:,seq2);
                        amp_sym = H_hr1.HcoeL(:,:,seq1) * H_hr2.HcoeL(:,:,seq2);
                        H_hr = H_hr.set_hop_mat(amp,vector,'set');
                        H_hr = H_hr.set_hop_mat(amp_sym,vector,'sym');
                    end
                end
            elseif isa(A,'HR') && ~isa(B,'HR')
                H_hr1 = A;
                H_hr = H_hr1;
                %                 H_hr2 = B;
                % ---------check----------
                if isa(B,'sym')
                    for i = 1:H_hr.NRPTS
                        H_hr.HcoeL(:,:,i) = H_hr.HcoeL(:,:,i) * B;
                    end
                elseif isa(B,'numeric')
                    for i = 1:H_hr.NRPTS
                        H_hr.HnumL(:,:,i) = H_hr.HnumL(:,:,i) * B;
                    end
                else
                    
                end
            elseif ~isa(A,'HR') && isa(B,'HR')
                H_hr1 = B;
                H_hr = H_hr1;
                %H_hr = H_hr.line_000_gen();
                %i = H_hr.Line_000;
                % ---------check----------
                if isa(A,'sym')
                    for i = 1:H_hr.NRPTS
                        H_hr.HcoeL(:,:,i) = A * H_hr.HcoeL(:,:,i);
                    end
                elseif isa(A,'numeric')
                    for i = 1:H_hr.NRPTS
                        H_hr.HnumL(:,:,i) = A * H_hr.HnumL(:,:,i);
                    end
                else
                    
                end
            end
        end
        % overload .^ power(H_hr,b)
        function H_hr = power(H_hr,b)
            for i = 1:H_hr.NRPTS
                H_hr.HcoeL(:,:,i) = H_hr.HcoeL(:,:,i).^b;
                H_hr.HnumL(:,:,i) = H_hr.HnumL(:,:,i).^b;
            end
        end
        % overload ^ power(H_hr,H_hr)
        function H_hr = mpower(A,B)
            if isa(A,'HR') && isa(B,'numeric')
                if length(B) >1
                    error('only support a number');
                end
                H_hr = A;
                c = B; % power num
                for i = 1:H_hr.NRPTS
                    H_hr.HcoeL(:,:,i) = H_hr.HcoeL(:,:,i)^c;
                    H_hr.HnumL(:,:,i) = H_hr.HnumL(:,:,i)^c;
                end
            elseif isa(B,'HR') && isa(A,'numeric')
                if length(A) >1
                    error('only support a number');
                end
                H_hr = B;
                c = A; % power num
                for i = 1:H_hr.NRPTS
                    H_hr.HcoeL(:,:,i) = c^H_hr.HcoeL(:,:,i);
                    H_hr.HnumL(:,:,i) = c^H_hr.HnumL(:,:,i);
                end
            else
                error('wrong input');
            end
            
        end
        % overload a == b eq(a,b)
        function varargout = eq(H_hr1,H_hr2)
            logical_num = true;
            if H_hr1.NRPTS ~= H_hr2.NRPTS
                logical_num = false;
                varargout{1} = logical_num;
                return;
            end
            if H_hr1.WAN_NUM ~= H_hr2.WAN_NUM
                logical_num = false;
                varargout{1} = logical_num;
                return;
            end
            if H_hr1.NRPTS ~= H_hr2.NRPTS
                logical_num = false;
                varargout{1} = logical_num;
                return;
            end
            [issame_list,reshape_list] = ismember(H_hr1.vectorL,H_hr2.vectorL,'rows');
            if sum(issame_list) ~= H_hr1.NRPTS
                logical_num = false;
                varargout{1} = logical_num;
                return;
            else
                reshape_list = reshape_list';
                H_hr2 = H_hr2.reseq(':',reshape_list);
                varargout{2} = H_hr2 ;
            end
            
            if H_hr1.HcoeL ~= H_hr2.HcoeL
                logical_num = false;
                varargout{1} = logical_num;
                varargout{2} = (H_hr1.HcoeL == H_hr2.HcoeL);
                return;
            end
            if ~isequal(H_hr1.HnumL ,H_hr2.HnumL)
                logical_num = false;
                varargout{1} = logical_num;
                return;
            end
            varargout{1} = logical_num;
            %             if H_hr1.el ~= H_hr2.NRPTS
            %                 logicalal_num = logical(0);
            %                 return;
            %             end
            %             if H_hr1.NRPTS ~= H_hr2.NRPTS
            %                 logicalal_num = logical(0);
            %                 return;
            %             end
        end
        % overload a ~= b ne(a,b)
        function logicalal_num = ne(H_hr1,H_hr2)
            logicalal_num = ~(H_hr1 == H_hr2);
        end
        % overload a' ctranspose(a)
        function H_hr = ctranspose(H_hr)
            H_hr = H_hr.dualize();
            if H_hr.vectorhopping
                H_hr.AvectorL = H_hr.AvectorL(H_hr.Duality_vector_dist,:);
                H_hr.BvectorL = -H_hr.BvectorL(H_hr.Duality_vector_dist,:);
                CLr = H_hr.CvectorL(1:end/2,:);
                CLi = H_hr.CvectorL(end/2+1:end,:);
                H_hr.CvectorL = [CLr(H_hr.Duality_vector_dist,:);-CLi(H_hr.Duality_vector_dist,:)];
            else
                
                if strcmp(H_hr.Type,'list')
                    if H_hr.num
                        H_hr.HnumL = conj(H_hr.HnumL(H_hr.Duality_vector_dist));
                    end
                    if H_hr.coe
                        H_hr.HcoeL = conj(H_hr.HcoeL(H_hr.Duality_vector_dist));
                    end
                else
                    if H_hr.coe
                        H_hr.HcoeL = pagectranspose(H_hr.HcoeL(:,:,H_hr.Duality_vector_dist));
                    end
                    if H_hr.num
                        H_hr.HnumL = pagectranspose(H_hr.HnumL(:,:,H_hr.Duality_vector_dist));
                    end
                end
            end
            %             for i = 1:H_hr.NRPTS
            %                 % one to one
            %                 H_hr.HcoeL(:,:,i) = H_hr.HcoeL(:,:,H_hr.Duality_vector_dist(i))';
            %                 H_hr.HnumL(:,:,i) = H_hr.HnumL(:,:,H_hr.Duality_vector_dist(i))';
            %             end
        end
        % overload a.' ctranspose(a)
        function H_hr = transpose(H_hr)
            H_hr.HcoeL = pagetranspose(H_hr.HcoeL);
            H_hr.HnumL = pagetranspose(H_hr.HnumL);
            %             for i = 1:H_hr.NRPTS
            %                 H_hr.HcoeL(:,:,i) = H_hr.HcoeL(:,:,i).';
            %                 H_hr.HnumL(:,:,i) = H_hr.HnumL(:,:,i).';
            %             end
        end
        % overload 	horzcat(a,b,...)
        function H_hr = horzcat(A,B)
            if isa(A,'HR') && isa(B,'HR')
                H_hr1 = A;
                H_hr2 = B;
                % ---------init-----------
                %                 H_hr =  HR(H_hr1.WAN_NUM+H_hr2.WAN_NUM,...
                %                     unique([H_hr1.vectorL;H_hr2.vectorL],'rows'));
                H_hr = A;
                H_hr.vectorL = unique([H_hr1.vectorL;H_hr2.vectorL],'rows');
                H_hr.orbL = [H_hr1.orbL;H_hr2.orbL];
                H_hr.quantumL = [H_hr1.quantumL;H_hr2.quantumL];
                H_hr.elementL = [H_hr1.elementL;H_hr2.elementL];
                H_hr.orb_symL = [H_hr1.orb_symL;H_hr2.orb_symL];
                H_hr.HnumL = zeros(H_hr1.WAN_NUM+H_hr2.WAN_NUM,H_hr1.WAN_NUM+H_hr2.WAN_NUM,H_hr.NRPTS);
                H_hr.HcoeL = sym(H_hr.HnumL);
                
                %             zeros_num_mat = zeros(H_hr.WAN_NUM);
                %             zeros_coe_mat = sym(zeros(H_hr.WAN_NUM));
                for i = 1:H_hr.NRPTS
                    vector = H_hr.vectorL(i,:);
                    [~,seq1]=ismember(vector,H_hr1.vectorL,'rows');
                    [~,seq2]=ismember(vector,H_hr2.vectorL,'rows');
                    if seq1 ~= 0 && seq2 ~=0
                        if H_hr1.num
                            amp = blkdiag(H_hr1.HnumL(:,:,seq1) ,...
                                H_hr2.HnumL(:,:,seq2));
                            H_hr = H_hr.set_hop_mat(amp,vector,'set');
                        end
                        if H_hr1.coe
                            amp_sym = blkdiag(H_hr1.HcoeL(:,:,seq1) ,...
                                H_hr2.HcoeL(:,:,seq2));
                            H_hr = H_hr.set_hop_mat(amp_sym,vector,'sym');
                        end
                    end
                end
            elseif isa(A,'HR') && ~isa(B,'HR')
                H_hr1 = A;
                H_hr = H_hr1;
                %                 H_hr2 = B;
                % ---------check----------
                if isa(B,'sym')
                    for i = 1:H_hr.NRPTS
                        H_hr.HcoeL(:,:,i) = blkdiag(H_hr.HcoeL(:,:,i),B);
                    end
                elseif isa(B,'numeric')
                    for i = 1:H_hr.NRPTS
                        H_hr.HnumL(:,:,i) = blkdiag(H_hr.HnumL(:,:,i),B);
                    end
                else
                    
                end
            elseif ~isa(A,'HR') && isa(B,'HR')
                H_hr1 = B;
                H_hr = H_hr1;
                % ---------check----------
                if isa(A,'sym')
                    for i = 1:H_hr.NRPTS
                        H_hr.HcoeL(:,:,i) = blkdiag( A , H_hr.HcoeL(:,:,i));
                    end
                elseif isa(A,'numeric')
                    for i = 1:H_hr.NRPTS
                        H_hr.HnumL(:,:,i) = blkdiag( A , H_hr.HnumL(:,:,i));
                    end
                else
                    
                end
            end
        end
        % overload 	vertcat(a,b,...)
        function H_hr = vertcat(A,B)
            H_hr = horzcat(A,B);
        end
        %         function subsref(a,s)
        % overload 	subsref(a,S)
        %%         function  varargout  = subsref(H_hr,S)
        %             if length(S) == 1
        %                 switch S(1).type
        %                     %                 case '.'
        %                     %                     builtin('subsref',S) % call builtin
        %                     case '{}'
        %                         switch length(S.subs)
        %                             case 1
        %                                 slist = S.subs{1};
        %                                 %                                 if ~issame(slist,':')
        %                                 %                                     H_hr = H_hr.set_NRPTS(length(slist));
        %                                 %                                 end
        %                                 H_hr.vectorL  =H_hr.vectorL(slist,:);
        %                                 H_hr.HnumL = H_hr.HnumL(:,:,slist);
        %                                 H_hr.HcoeL = H_hr.HcoeL(:,:,slist);
        %                             case 2 % extract mat
        %                                 slist = S.subs{1};
        %                                 numorsym_label =  S.subs{2};
        %                                 if numorsym_label == 0 % coe
        %                                     H_hr = H_hr.HcoeL(:,:,slist);
        %                                 elseif numorsym_label == 1 % num
        %                                     H_hr = H_hr.HnumL(:,:,slist);
        %                                 elseif isequal(numorsym_label,':')
        %                                     H_hr_tmp.HcoeL = H_hr.HcoeL(:,:,slist);
        %                                     H_hr_tmp.HnumL = H_hr.HnumL(:,:,slist);
        %                                     disp(H_hr.vectorL(slist,:));
        %                                     disp(H_hr_tmp.HcoeL );
        %                                     disp(H_hr_tmp.HnumL);
        %                                     H_hr = H_hr_tmp;
        %                                 end
        %                         end
        %                         varargout{1} =  H_hr;
        %                     case '()'
        %                         switch length(S.subs)
        %                             case 1
        %                                 slist = S.subs{1};
        %                                 if ~isequal(slist,':')
        %                                     %                                     H_hr = H_hr.set_WAN_NUM(length(slist));
        %                                     H_hr.HnumL=H_hr.HnumL(slist,slist,:);
        %                                     H_hr.HcoeL=H_hr.HcoeL(slist,slist,:);
        %                                 else
        %
        %                                 end
        %                             case 2 % reshape data
        %                                 wan_list = S.subs{1};
        %                                 nrpt_list =  S.subs{2};
        %                                 H_hr = H_hr.reseq(wan_list,nrpt_list);
        %                             case 3 % H_hr
        %                                 V = H_hr.vectorL;
        %                                 if isequal(S.subs{1},':')
        %                                     xlist =min(V(:,1)):max(V(:,1));
        %                                 else
        %                                     xlist = S.subs{1};
        %                                 end
        %                                 if isequal(S.subs{2},':')
        %                                     ylist =min(V(:,2)):max(V(:,2));
        %                                 else
        %                                     ylist = S.subs{2};
        %                                 end
        %                                 if isequal(S.subs{3},':')
        %                                     zlist =min(V(:,3)):max(V(:,3));
        %                                 else
        %                                     zlist = S.subs{3};
        %                                 end
        %                                 count = 0;
        %                                 for k = 1: length(zlist)
        %                                     for j = 1: length(ylist)
        %                                         for i = 1: length(xlist)
        %                                             count = count +1;
        %                                             vector_list(count,:) = [xlist(i),...
        %                                                 ylist(j),...
        %                                                 zlist(k)];
        %                                         end
        %                                     end
        %                                 end
        %                                 [~,label_list] = ismember(vector_list,V,'row');
        %                                 label_list_true = label_list(label_list>0);
        %                                 H_hr = H_hr.reseq(':',label_list_true);
        %                             case 4 % extract mat
        %                                 S_temp = S;
        %                                 S_temp.subs(4) = [];
        %                                 H_hr_temp = H_hr.subsref(S_temp);
        %                                 numorsym_label =  S.subs{4};
        %                                 if numorsym_label == 0 % coe
        %                                     H_hr = H_hr_temp.HcoeL;
        %                                 elseif numorsym_label == 1 % num
        %                                     H_hr = H_hr_temp.HnumL;
        %                                 elseif isequal(numorsym_label,':')
        %                                     disp(H_hr_temp.vectorL);
        %                                     %disp(H_hr.HcoeL );
        %                                     %disp(H_hr.HnumL);
        %                                     H_hr_tmp.HcoeL = H_hr_temp.HcoeL;
        %                                     H_hr_tmp.HnumL = H_hr_temp.HnumL;
        %                                     H_hr = H_hr_tmp;
        %                                 end
        %                         end
        %                         varargout{1} = H_hr;
        %                     otherwise
        %                         varargout = builtin('subsref',H_hr,S);% return inner func
        %                         nOutputs = nargout;
        %                         if nOutputs > 1
        %                             disp('diappointed');
        %                         elseif isa(varargout,'cell')
        %                             disp('diappointed');
        %                         else
        %
        %                             temppout{1} = varargout;
        %                             varargout = temppout;
        %
        %                         end
        %                 end
        %                 %             elseif length(S) == 2
        %                 %                 disp('?')
        %             else
        %
        %                 nOutputs = nargout;
        %                 if nOutputs > 1
        %                     %disp('diappointed');
        %                     switch nOutputs
        %                         case 2
        %                             [varargout{1},varargout{2}]= builtin('subsref',H_hr,S);% return inner func
        %                         case 3
        %                             [varargout{1},varargout{2},varargout{3}]= builtin('subsref',H_hr,S);% return inner func
        %                         case 4
        %                             [varargout{1},varargout{2},varargout{3},varargout{4}]= builtin('subsref',H_hr,S);% return inner func
        %                         case 5
        %                             [varargout{1},varargout{2},varargout{3},varargout{4},varargout{5}]= builtin('subsref',H_hr,S);% return inner func
        %                         case 6
        %                             [varargout{1},varargout{2},varargout{3},varargout{4},varargout{5},varargout{6}]= builtin('subsref',H_hr,S);% return inner func
        %                         case 7
        %                             [varargout{1},varargout{2},varargout{3},varargout{4},varargout{5},varargout{6},varargout{7}]= builtin('subsref',H_hr,S);%
        %                         otherwise
        %
        %                     end
        %                 else
        %
        %                     varargout = builtin('subsref',H_hr,S);% return inner func
        %                     temppout{1} = varargout;
        %                     varargout = temppout;
        %
        %                 end
        %             end
        %
        %         end
        %%
        % overload 	gt(A,B)
        function C = gt(B,A)
            if isa(A,'HR') && isa(B,'HR')
                H_hr1 = A;
                H_hr2 = B;
                % ---------check----------
                if H_hr1.WAN_NUM ~= H_hr2.WAN_NUM
                    error('WAN_NUM different');
                end
                error('not support at present.')
            elseif isa(B,'HR') && ~isa(A,'HR')
                switch class(A)
                    case 'char'
                        switch A(1)
                            case {'P','p'}
                                
                            case {'w','W'}
                                C = B.Gen_hr(A,'hr_dat');
                            otherwise
                        end
                    otherwise
                end
            elseif isa(A,'HR') && ~isa(B,'HR')
                switch class(B)
                    case 'char'
                        switch B(1)
                            case {'P','p'}
                                C = A.input_orb_struct(B,'tbsk','symbolic',true);
                            case {'w','W'}
                                
                            case {'k','K'}
                                C = A.kpathgen3D(B);
                            otherwise
                        end
                    case 'double'
                        switch size(B,1)
                            case A.WAN_NUN
                                C = A.input_orb_init(B);
                            otherwise
                        end
                    otherwise
                end
            end
            
        end
        % overload 	lt(A,B)
        function C = lt(A,B)
            if isa(A,'HR') && isa(B,'HR')
                H_hr1 = A;
                H_hr2 = B;
                % ---------check----------
                if H_hr1.WAN_NUM ~= H_hr2.WAN_NUM
                    error('WAN_NUM different');
                end
                error('not support at present.')
            elseif isa(A,'HR') && ~isa(B,'HR')
                switch class(B)
                    case 'char'
                        switch B(1)
                            case {'P','p'}
                                C = A.input_orb_struct(B,'tbsk');
                            case {'w','W'}
                                
                            case {'k','K'}
                                C = A.kpathgen3D(B);
                            otherwise
                        end
                    case 'double'
                        switch size(B,1)
                            case A.WAN_NUN
                                C = A.input_orb_init(B);
                            otherwise
                        end
                    otherwise
                end
            elseif ~isa(A,'HR') && isa(B,'HR')
                error('not support at present.');
            end
            
        end
        % overload le
        function C = le(A,B)
            if isa(A,'HR') && isa(B,'HR')
                H_hr1 = A;
                H_hr2 = B;
                % ---------check----------
                if H_hr1.WAN_NUM ~= H_hr2.WAN_NUM
                    error('WAN_NUM different');
                end
                error('not support at present.')
            elseif isa(A,'HR') && ~isa(B,'HR')
                switch class(B)
                    case 'char'
                        switch B(1)
                            case {'P','p'}
                                C = A.input_orb_struct(B,'vasp');
                            case {'w','W'}
                                
                            case {'k','K'}
                                C = A.kpathgen3D(B);
                            otherwise
                        end
                    case 'double'
                        switch size(B,1)
                            case A.WAN_NUN
                                C = A.input_orb_init(B);
                            otherwise
                        end
                    otherwise
                end
            elseif ~isa(A,'HR') && isa(B,'HR')
                error('not support at present.');
            end
            
        end
        % ---------------------  overload function   ----------------------
        function H_hr = kron(A,B)
            if isa(A,'HR') && isa(B,'HR')
                H_hr1 = A;
                H_hr2 = B;
                % ---------init-----------
                H_hr =  HR(H_hr1.WAN_NUM * H_hr2.WAN_NUM,...
                    unique([H_hr1.vectorL;H_hr2.vectorL],'rows'));
                for i = 1:H_hr.NRPTS
                    vector = H_hr.vectorL(i,:);
                    [~,seq1]=ismember(vector,H_hr1.vectorL,'rows');
                    [~,seq2]=ismember(vector,H_hr2.vectorL,'rows');
                    if seq1 ~= 0 && seq2 ~=0
                        amp = kron(H_hr1.HnumL(:,:,seq1) ,...
                            H_hr2.HnumL(:,:,seq2));
                        amp_sym = kron(H_hr1.HcoeL(:,:,seq1) ,...
                            H_hr2.HcoeL(:,:,seq2));
                        H_hr = H_hr.set_hop_mat(amp,vector,'set');
                        H_hr = H_hr.set_hop_mat(amp_sym,vector,'sym');
                    end
                end
            elseif isa(A,'HR') && ~isa(B,'HR')
                H_hr1 = A;
                H_hr = H_hr1;
                %                 H_hr = H_hr.set_WAN_NUM(H_hr.WAN_NUM*length(B));
                %                 H_hr2 = B;
                % ---------check----------
                if isa(B,'sym')
                    for i = 1:H_hr.NRPTS
                        H_hr.HcoeL(:,:,i) = kron(H_hr.HcoeL(:,:,i),B);
                    end
                elseif isa(B,'numeric')
                    for i = 1:H_hr.NRPTS
                        H_hr.HnumL(:,:,i) = kron(H_hr.HnumL(:,:,i),B);
                    end
                else
                    
                end
            elseif ~isa(A,'HR') && isa(B,'HR')
                H_hr1 = B;
                H_hr = H_hr1;
                %                 H_hr = H_hr.set_WAN_NUM(H_hr.WAN_NUM*length(A));
                % ---------check----------
                if isa(A,'sym')
                    for i = 1:H_hr.NRPTS
                        H_hr.HcoeL(:,:,i) = kron( A , H_hr.HcoeL(:,:,i));
                    end
                elseif isa(A,'numeric')
                    for i = 1:H_hr.NRPTS
                        H_hr.HnumL(:,:,i) = kron( A , H_hr.HnumL(:,:,i));
                    end
                else
                    
                end
            end
        end
        function H_hr = conj(H_hr)
            H_hr = H_hr.dualize();
            H_hr.HcoeL = conj(H_hr.HcoeL(:,:,H_hr.Duality_vector_dist));
            H_hr.HnumL = conj(H_hr.HnumL(:,:,H_hr.Duality_vector_dist));
            %             for i = 1:H_hr.NRPTS
            %                 H_hr.HcoeL(:,:,i) = conj(H_hr.HcoeL(:,:,H_hr.Duality_vector_dist(i)));
            %                 H_hr.HnumL(:,:,i) = conj(H_hr.HnumL(:,:,H_hr.Duality_vector_dist(i)));
            %             end
        end
        function H_hr = full(H_hr)
            if ~strcmp(H_hr.Type,'sparse')
                return;  %  error('The type is not the sparse ');
            end
            HnumL_temp = zeros(H_hr.WAN_NUM,H_hr.WAN_NUM,H_hr.NRPTS);
            for i = 1:H_hr.NRPTS
                HnumL_temp(:,:,i) = full(H_hr.HnumL{i});
            end
            H_hr.HnumL = HnumL_temp;
            H_hr.Type = 'mat';
        end
        function H_hr = sparse(H_hr)
            if strcmp(H_hr.Type,'sparse')
                return;
            end
            HnumL_temp{H_hr.NRPTS} = zeros(H_hr.WAN_NUM,H_hr.WAN_NUM);
            for i = 1:H_hr.NRPTS
                HnumL_temp{i}= sparse(H_hr.HnumL(:,:,i));
            end
            H_hr.HnumL = HnumL_temp;
            H_hr.HcoeL = [];
            H_hr.Type = 'sparse';
        end
        function [H_hr,EQL] = subs(H_hr,varargin)
            if strcmp(varargin{end},'all')
                all_mode = true;
                nargin_check = nargin-1;
            else
                all_mode = false;
                nargin_check = nargin;
            end
            switch nargin_check
                case 1
                    H_hr.HcoeL = expand(simplify(subs(H_hr.HcoeL)));
                case 2
                    SymVarL = H_hr.symvar_list;
                    if all_mode
                        H_hr.HcoeL = subs(H_hr.HcoeL,SymVarL,varargin{1});
                    else
                        if strcmp(H_hr.Type , 'list')
                            H_hr.HcoeL = expand(simplify(subs(H_hr.HcoeL,SymVarL,varargin{1})));
                        else
                            H_hr.HcoeL = (simplify(subs(H_hr.HcoeL,SymVarL,varargin{1})));
                        end
                    end
                    EQL =(SymVarL==vpa(varargin{1}));
                case 3
                    if all_mode
                        H_hr.HcoeL = subs(H_hr.HcoeL,varargin{1},varargin{2});
                    else
                        if strcmp(H_hr.Type , 'list')
                            H_hr.HcoeL = expand(simplify(subs(H_hr.HcoeL,varargin{1},varargin{2})));
                        else
                            H_hr.HcoeL = simplify(subs(H_hr.HcoeL,varargin{1},varargin{2}));
                        end
                    end
                    EQL =(varargin{1}==vpa(varargin{2}));
                case 5
                    if strcmp(H_hr.Type , 'list')
                        H_hr.HcoeL = expand(simplify(subs(H_hr.HcoeL,varargin{1},varargin{2})));
                        H_hr.ScoeL = expand(simplify(subs(H_hr.ScoeL,varargin{3},varargin{4})));
                    else
                        H_hr.HcoeL = simplify(subs(H_hr.HcoeL,varargin{1},varargin{2}));
                        H_hr.ScoeL = simplify(subs(H_hr.ScoeL,varargin{3},varargin{4}));
                    end
                    EQL{1} =(varargin{1}==vpa(varargin{2}));
                    EQL{2} =(varargin{3}==vpa(varargin{3}));
            end
            if  isempty(H_hr.symvar_list) && ~H_hr.overlap
                H_hr = H_hr.Subsall();
            elseif isempty(H_hr.symvar_list) &&  H_hr.overlap && isempty(symvar(H_hr.ScoeL))
                H_hr = H_hr.Subsall();
            else
                
            end
            % waiting to add ...
        end
        function H_hr = simplify(H_hr,Accuracy)
            if nargin < 2
                Accuracy = 1e-6;
            end
            if H_hr.vectorhopping % becareful! numerical error!
                %                 AL = orth(H_hr.AvectorL);
                %                 BL = orth(H_hr.BvectorL);
                AL = rref(H_hr.AvectorL.').';
                BL = rref(H_hr.BvectorL.').';
                CL = rref(H_hr.CvectorL.').';
                AL = AL(:,1:rank(AL));
                BL = BL(:,1:rank(BL));
                CL = CL(:,1:rank(CL));
                rAL = real(AL);rAL(abs(rAL)<Accuracy) = 0;
                %iAL = imag(AL);iAL(abs(iAL)<Accuracy) = 0;
                rBL = real(BL);rBL(abs(rBL)<Accuracy) = 0;
                %iBL = imag(BL);iBL(abs(iBL)<Accuracy) = 0;
                rCL = real(CL);rCL(abs(rCL)<Accuracy) = 0;
                %iCL = imag(CL);iCL(abs(iCL)<Accuracy) = 0;
                H_hr.AvectorL = rAL;%+1i*iAL;
                H_hr.BvectorL = rBL;%+1i*iBL;
                % test
                H_hr.CvectorL = rCL;%+1i*iCL;
                %NRPTS_list = false(H_hr.NRPTS,1);
                %                 NRPTS_list = sum(abs(H_hr.AvectorL)>Accuracy^2,2)+...
                %                     sum(abs(H_hr.BvectorL)<Accuracy^2,2);
                %                 H_hr = H_hr.reseq(':',logical(NRPTS_list));
                return;
            end
            if H_hr.coe
                H_coeL_tmp = simplify(H_hr.HcoeL);
                H_hr.HcoeL = H_coeL_tmp;
                if strcmp(H_hr.Type,'list')
                    NRPTS_list = find(H_coeL_tmp ~=sym(0));
                    H_hr = H_hr.reseq(':',NRPTS_list);
                    %[unique_HcoeL,ia,ic] = H_hr.HcoeL;
                elseif strcmp(H_hr.Type,'mat')
                    
                end
            end
            if H_hr.num
                H_numL_tmp = H_hr.HnumL;
                %H_hr.HnumL = H_numL_tmp;
                if strcmp(H_hr.Type,'list')
                    NRPTS_list = find(abs(H_numL_tmp) > Accuracy);
                    H_hr = H_hr.reseq(':',NRPTS_list);
                    %[unique_HcoeL,ia,ic] = H_hr.HcoeL;
                elseif strcmp(H_hr.Type,'mat')
                    zerosMat = ones(size(H_numL_tmp(:,:,1)))*Accuracy;
                    NRPTS_list = true(H_hr.NRPTS,1);
                    for i = 1:H_hr.NRPTS
                        if sum(abs(H_numL_tmp(:,:,i)) > zerosMat,'all')
                            %NRPTS_list(i) = true;
                        else
                            NRPTS_list(i) = false;
                        end
                    end
                    H_hr = H_hr.reseq(':',NRPTS_list);
                end
            end
        end
        function H_hr = filter(H_hr,Accuracy)
            if nargin < 2
                Accuracy = 1e-6;
            end
            switch H_hr.Type
                case {'list','mat'}
                    HnumList = H_hr.HnumL;
                    HnumList_real = real(HnumList);
                    HnumList_imag = imag(HnumList);
                    HnumList_real(abs(HnumList_real) < Accuracy) = 0;
                    HnumList_imag(abs(HnumList_imag) < Accuracy) = 0;
                    H_hr.HnumL = HnumList_real + 1i*HnumList_imag;
                otherwise
            end
            
        end
        function Hsym = sym(H_hr,options)
            arguments
                H_hr HR;
                options.cartesian = true;
                options.simple = false;
            end
            if strcmp(H_hr.Type,'list')
                H_hr = H_hr.rewind();
            elseif strcmp(H_hr.Type,'sparse')
                H_hr = H_hr.full();
            end
            Hsym = zeros(H_hr.WAN_NUM,'sym');
            H_hr = H_hr.tjmti_gen('sym');
            vectorList = double(H_hr.vectorL);
            %tij_mat_r = H_hr.tjmti{1};
            [H_hr.num,~] = H_hr.NumOrCoe;
            if H_hr.num
                H_hr.HcoeL = sym(H_hr.HnumL);
            end
            if options.cartesian
                tij_mat_k = H_hr.tjmti{3};
                syms k_x k_y k_z real;
                exp_pre = exp(1i... % 1i possible
                    *[k_x k_y k_z]*(...
                    vectorList*H_hr.Rm)');
                for i = 1:H_hr.NRPTS
                    Hsym = Hsym+H_hr.HcoeL(:,:,i)*exp_pre(i);
                end
                if options.simple
                else
                    Hsym = Hsym.*tij_mat_k;
                end
            else
                tij_mat_s = H_hr.tjmti{4};
                syms k_1 k_2 k_3 real;
                exp_pre = exp(1i*2*pi... % 1i possible
                    *[k_1 k_2 k_3]*(...
                    vectorList)');
                for i = 1:H_hr.NRPTS
                    Hsym = Hsym+H_hr.HcoeL(:,:,i)*exp_pre(i);
                end
                if options.simple
                else
                    Hsym = Hsym.*tij_mat_s;
                end
            end
        end
        function H_hr = sum(H_hr_list)
            H_hr = H_hr_list(1);
            for i = 2:length(H_hr_list)
                H_hr = H_hr + H_hr_list(i);
            end
        end
        function [H_hr,Sublist,Unique_term] = unique(H_hr,seed,checklist,Accuracy)
            arguments
                H_hr HR;
                seed char = 'gamma';
                checklist =sym([sqrt(3)]);
                Accuracy = (1e-6);
            end
            if strcmp(H_hr.Type,'list')
                HcoeL_tmp = [real(H_hr.HcoeL);imag(H_hr.HcoeL)];
                % Magnification divide
                [CoeffsList,factor_list_2] = vasplib.factorAll(HcoeL_tmp);
                % unique
                [Unique_term,~,ic] = unique(factor_list_2);
                %
                SymVar = sym(seed,[length(Unique_term),1],'real');
                % zero test
                SymVar(find(Unique_term==sym(0))) = sym(0);
                %SymVar(find(abs(CoeffsList)<sym(Accuracy))) = sym(0); % it
                % is a track
                %
                %
                % if we use factor all it seems like we need not apply the
                % check routine!
                % minus test sqrt3 test 2 test
                for coeff_for_check = checklist
                    [Lia,Locb] = ismember(expand(Unique_term),expand(coeff_for_check*Unique_term));
                    Equationlist = HR.isolateAll((SymVar(Lia) == coeff_for_check*SymVar(Locb(Locb>0)) ));
                    SymVar = subs(SymVar,lhs(Equationlist),rhs(Equationlist));
                end
                % waiting
                Sublist = SymVar == Unique_term;
                %Sublist_i = SymVar_i == Unique_term_i;
                H_hr.HcoeL = SymVar(ic(1:end/2)) .* CoeffsList(1:end/2) ...
                    +1i*SymVar(ic(end/2+1:end)) .* CoeffsList(end/2+1:end) ;
            end
        end
    end
    %% modify function
    % ----------------  change --------------------
    methods
        function H_hr = reseq(H_hr,wan_list,nrpt_list,nrpt_list_S)
            if nargin < 3
                nrpt_list = ':';
            end
            if nargin < 4
                nrpt_list_S = ':';
            end
            % wan first
            if ~isequal(wan_list,':')
                %                 H_hr = H_hr.set_WAN_NUM(length(wan_list));
                if strcmp(H_hr.Type,'sparse')
                    for i = 1:H_hr.NRPTS
                        H_hr.HnumL{i}=H_hr.HnumL{i}(wan_list,wan_list);
                    end
                elseif strcmp(H_hr.Type,'list')
                    wan_list =  all(ismember(H_hr.vectorL(:,4:5),int32(wan_list)),2);
                    H_hr.HnumL=H_hr.HnumL(wan_list,:);
                    H_hr.HcoeL=H_hr.HcoeL(wan_list,:);
                    H_hr.vectorL=H_hr.vectorL(wan_list,:);
                    if H_hr.overlap
                        H_hr.SnumL=H_hr.SnumL(wan_list,:);
                        H_hr.ScoeL=H_hr.ScoeL(wan_list,:);
                        H_hr.vectorL_overlap=H_hr.vectorL_overlap(wan_list,:);
                    end
                elseif strcmp(H_hr.Type,'mat')
                    if H_hr.num
                        H_hr.HnumL=H_hr.HnumL(wan_list,wan_list,:);
                        if H_hr.overlap
                            H_hr.SnumL=H_hr.SnumL(wan_list,wan_list,:);
                        end
                    else
                        H_hr.HnumL = [];
                        if H_hr.overlap
                            H_hr.SnumL=[];
                        end
                    end
                    if H_hr.coe
                        H_hr.HcoeL=H_hr.HcoeL(wan_list,wan_list,:);
                        if H_hr.overlap
                            H_hr.ScoeL=H_hr.ScoeL(wan_list,wan_list,:);
                        end
                    else
                        H_hr.HcoeL = sym([]);
                        if H_hr.overlap
                            H_hr.ScoeL=sym([]);
                        end
                    end
                end
                if ~isempty(H_hr.sites)
                    try
                        H_hr.sites = H_hr.sites(wan_list);
                    catch
                        % bug
                    end
                end
                if ~isempty( H_hr.orbL )
                    H_hr.orbL = H_hr.orbL(wan_list,:);
                end
                if ~isempty( H_hr.elementL )
                    H_hr.elementL = H_hr.elementL(wan_list,:);
                end
                if ~isempty( H_hr.quantumL)
                    H_hr.quantumL = H_hr.quantumL(wan_list,:);
                end
                % bug here sym_orbL waiting
            end
            % nrpt first
            if ~isequal(nrpt_list,':')
                %                 H_hr = H_hr.set_NRPTS(length(nrpt_list));
                H_hr.vectorL=H_hr.vectorL(nrpt_list,:);
                if H_hr.overlap
                    H_hr.vectorL_overlap=H_hr.vectorL_overlap(nrpt_list_S,:);
                end
                if strcmp(H_hr.Type,'sparse')
                    H_hr.HnumL=H_hr.HnumL(nrpt_list);
                elseif strcmp(H_hr.Type,'list')
                    if H_hr.vectorhopping
                        H_hr.AvectorL=H_hr.AvectorL(nrpt_list,:);
                        H_hr.BvectorL=H_hr.BvectorL(nrpt_list,:);
                        CL1 = H_hr.CvectorL(1:end/2,:);
                        CL2 = H_hr.CvectorL(end/2+1:end,:);
                        H_hr.CvectorL = [CL1(nrpt_list,:);CL2(nrpt_list,:)];
                        return;
                    end
                    if H_hr.num
                        H_hr.HnumL=H_hr.HnumL(nrpt_list);
                    end
                    if H_hr.coe
                        H_hr.HcoeL=H_hr.HcoeL(nrpt_list);
                    end
                    if H_hr.overlap
                        if H_hr.num
                            H_hr.SnumL=H_hr.SnumL(nrpt_list_S,:);
                        end
                        if H_hr.coe
                            H_hr.ScoeL=H_hr.ScoeL(nrpt_list_S,:);
                        end
                    end
                elseif strcmp(H_hr.Type,'mat')
                    if H_hr.num
                        H_hr.HnumL=H_hr.HnumL(:,:,nrpt_list);
                        if H_hr.overlap
                            H_hr.SnumL=H_hr.SnumL(:,:,nrpt_list);
                        end
                    end
                    if H_hr.coe
                        H_hr.HcoeL=H_hr.HcoeL(:,:,nrpt_list);
                        if H_hr.overlap
                            H_hr.ScoeL=H_hr.ScoeL(:,:,nrpt_list);
                        end
                    end
                end
            end
        end
        function H_hr = cut_orb(H_hr,rm_list,options)
            arguments
                H_hr HR;
                rm_list = [];
                options.rmfunc  function_handle=@()(1);
            end
            orb_tmp = H_hr.orbL;
            if isempty(rm_list) && ~strcmp(functions(options.rmfunc).function , '@()(1)')
                rm_list = options.rmfunc(orb_tmp(:,1),orb_tmp(:,2),orb_tmp(:,3));
            elseif isempty(rm_list)
                rm_list = zeros(size(H_hr.orbL,1),1);
            else
                rmlist = zeros(size(H_hr.orbL,1),1);
                rmlist(rm_list) = 1;
                rm_list = rmlist;
            end
            wan_list = ~rm_list;
            H_hr = H_hr.reseq(wan_list);
        end
        function H_hr = clean(H_hr,WANNUM)
            if nargin < 2
                WANNUM = H_hr.WAN_NUM;
            end
            if strcmp(H_hr.Type,'Sparse')
            elseif strcmp(H_hr.Type,'mat')
                H_hr.HnumL = zeros(WANNUM,WANNUM,H_hr.NRPTS);
                if H_hr.coe
                    H_hr.HcoeL = sym(H_hr.HnumL);
                else
                    H_hr.HcoeL =sym([]);
                end
                if H_hr.overlap
                    H_hr.SnumL = zeros(WANNUM,WANNUM,H_hr.NRPTS);
                    if H_hr.coe
                        H_hr.ScoeL = sym(H_hr.SnumL);
                    else
                        H_hr.ScoeL =sym([]);
                    end
                    
                end
            end
        end
        function H_hr = project(H_hr,BASIS_MAT)
            new_WAN_NUM = size(BASIS_MAT,1);
            BASIS_MAT_prime = BASIS_MAT';
            switch H_hr.Type
                case 'sparse'
                    disp('waiting');
                case 'list'
                    disp('waiting');
                case 'mat'
                    %disp('testing')
                    if H_hr.num
                        H_hr.HcoeL = zeros(new_WAN_NUM,new_WAN_NUM,H_hr.NRPTS);
                        tempHnumL = zeros(new_WAN_NUM,new_WAN_NUM,H_hr.NRPTS);
                        for i = 1:H_hr.NRPTS
                            tempHnumL(:,:,i) = BASIS_MAT*H_hr.HnumL(:,:,i)*BASIS_MAT_prime ;
                        end
                        H_hr.HnumL = tempHnumL;
                    end
                    if H_hr.coe
                        H_hr.HnumL = zeros(new_WAN_NUM,new_WAN_NUM,H_hr.NRPTS);
                        tempHcoeL = sym(zeros(new_WAN_NUM,new_WAN_NUM,H_hr.NRPTS));
                        for i = 1:H_hr.NRPTS
                            tempHcoeL(:,:,i) = BASIS_MAT*H_hr.HcoeL(:,:,i)*BASIS_MAT_prime ;
                        end
                        H_hr.HcoeL = tempHcoeL;
                    end
            end
        end
        function H_hr = charalize(H_hr)
            for i = 1:H_hr.WAN_NUM
                for j = 1:H_hr.WAN_NUM
                    if H_hr.elementL(i) == H_hr.elementL(j) && i~=j
                        H_hr.HnumL(i,j,:) = 0;
                    end
                end
            end
        end
        function H_hr = ForceToMat(H_hr)
            switch H_hr.type
                case 'sparse'
                    H_hr =H_hr.full();
                case 'mat'

                case 'list'
                    H_hr = H_hr.rewind;
                otherwise
                    error('Not support yet.');
            end
        end
        function H_hr = ForceTosparse(H_hr)
            switch H_hr.type
                case 'sparse'
                case 'mat'
                    H_hr = H_hr.sparse;
                case 'list'
                    H_hr = H_hr.rewind;
                    H_hr = H_hr.rewind;
                otherwise
                    error('Not support yet.');
            end
        end
        function H_hr = ForceTolist(H_hr)
            switch H_hr.type
                case 'sparse'
                    H_hr =H_hr.full();
                    H_hr = H_hr.rewrite;
                case 'mat'
                    H_hr = H_hr.rewrite;
                case 'list'

                otherwise
                    error('Not support yet.');
            end
        end
    end
    % ----------------  add something --------------------
    methods
        function H_hr = enlarge(H_hr,dir,amp)
            if isa(dir,'char')
                H_hr = H_hr.rewrite;
                vectorList = H_hr.vectorL;
                orbList= H_hr.orbL;
                orbListdiff = orbList(vectorList(:,4),:) - orbList(vectorList(:,5),:);
                orbListdiff_r = orbListdiff*H_hr.Rm;
                switch dir
                    case 'x'
                        selectL = logical(orbListdiff_r(:,1));
                    case 'y'
                        selectL = logical(orbListdiff_r(:,2));
                    case 'z'
                        selectL = logical(orbListdiff_r(:,3));
                end
                H_hr.HnumL(selectL)=H_hr.HnumL(selectL)*amp;
            else
                
                
                switch size(dir,2)
                    case 1
                        selectL = H_hr.vectorL(:,dir) ;
                    case 2
                    case 3
                    case 5
                end
            end
        end
        function H_hr = add_soc(H_hr)
            if ~H_hr.soc
                quantumList_up = H_hr.quantumL;
                quantumList_up(:,4) = 1;
                quantumList_dn = quantumList_up;
                quantumList_dn(:,4) = -1;
                % H_hr = [H_hr,conj(H_hr)]; % check !
                H_hr = [H_hr,(H_hr)];
            end
            H_hr.quantumL = [quantumList_up;quantumList_dn];
            H_hr.HcoeL = sym(H_hr.HnumL);
            H_hr.coe = true;
            H_hr.num = false;% check!
            H_hr = H_hr + H_hr.H_atom_soc();
        end
        % ? why duplicate funciton here?
        function H_hr = addorb(H_hr,orblist,options)
            arguments
                H_hr HR;
                orblist;
                options.inner  = true;
            end
            WANNUM = H_hr.WAN_NUM;
            norb = length(orblist);
            HnumLtmp = H_hr.HnumL;
            HnumLtmp2 = zeros(size(H_hr.HnumL)+[norb,norb,0]);
            HnumLtmp2(1:WANNUM,1:WANNUM,:) = HnumLtmp;
            HnumLtmp2(WANNUM+1:WANNUM+norb,1:WANNUM,:) = HnumLtmp(orblist,:,:);
            HnumLtmp2(1:WANNUM,WANNUM+1:WANNUM+norb,:) = HnumLtmp(:,orblist,:);
            H_hr.HnumL = HnumLtmp2;
        end
        function H_hr = add_orb(H_hr,hop_struct,orbOne,QuantumOne,elementOne)
            nstruct = length( hop_struct);
            if nargin < 3
                orbOne = repmat([0,0,0],[nstruct,1]);
            end
            if nargin < 4
                QuantumOne = repmat([1,0,0,1],[nstruct,1]);
            end
            if nargin < 5
                elementOne = ones(nstruct,1);
            end
            WANNUM = H_hr.WAN_NUM;
            
            H_hr = H_hr.expand_empty_one(orbOne,QuantumOne,elementOne);
            for j = 1:nstruct
                
                for i = 1:length(hop_struct(j).hop)
                    H_hr = H_hr.set_hop(hop_struct(j).hop(i),hop_struct(j).hi(i),WANNUM+j,hop_struct(j).vector,'set');
                    H_hr = H_hr.set_hop(conj(hop_struct(j).hop(i)),WANNUM+j,hop_struct(j).hi(i),-hop_struct(j).vector,'set');
                end
            end
        end
        function H_hr = addsoc(H_hr,quantumL)
            if nargin >1
                % set
                H_hr.quantumL = quantumL;
            end
            H_hr = H_hr + H_hr.H_atom_soc;
        end
        function H_hr = deltarule(H_hr,level_cut,mode,options)
            arguments
                H_hr HR;
                level_cut double{mustBeInteger} = 1;
                mode double = 0;
                options.Rd = -1;
            end
            import park.*;
            if mode == 0
                return;
            end
            if level_cut == -1
                level_cut = length(H_hr.Rnn);
            end
            % init
            Rnn = H_hr.nn_information();
            base_string = ["VssS","VspS","VsdS","VppS","VpdS","VppP","VpdP","VddS","VddP","VddD"];
            strvar_list = string(H_hr.symvar_list);
            base_symvar_min = sym([]);
            base_num_min = ones(1,length(base_string));
            count = 0;
            for ibase_string = base_string
                for i = 1:level_cut
                    if contains(ibase_string+"_"+string(i),strvar_list)
                        count = count +1;
                        base_symvar_min(count) = sym(ibase_string+"_"+string(i),'real');
                        base_num_min(count) = i;
                        base_string(count) = ibase_string;
                        break;
                    end
                end
            end
            base_string(count+1:end) = [];
            base_num_min(count+1:end) = [];
            if H_hr.overlap
                base_string_S = ["SssS","SspS","SsdS","SppS","SpdS","SppP","SpdP","SddS","SddP","SddD"];
                symvar_list_S = symvar(vasplibobj.ScoeL);
                strvar_list_S  = [string(symvar_list_S)];
                base_symvar_min_S = sym([]);
                base_num_min_S = ones(1,length(base_string));
                count_S = 0;
                for ibase_string_S = base_string_S
                    for i = 1:level_cut
                        if contains(ibase_string_S+"_"+string(i),strvar_list_S)
                            count_S = count_S +1;
                            base_symvar_min_S(count_S) = sym(ibase_string_S+"_"+string(i),'real');
                            base_num_min_S(count_S) = i;
                            base_string_S(count_S) = ibase_string_S;
                            break;
                        end
                    end
                end
                base_string_S(count+1:end) = [];
                base_num_min_S(count+1:end) = [];
            end
            % Rd
            if options.Rd == -1
                RdL = Rnn;
            else
                disp(base_string);
                if length(options.Rd) > 1
                    RdL = repmat(options.Rd,[1 count]);
                else
                    RdL = options.Rd;
                end
                base_num_min = 1:count;
                if H_hr.overlap
                    base_num_min_S = 1:count_S;
                end
            end
            %
            switch mode
                case 0
                    return;
                case 1
                    for j = 2:level_cut
                        for i= 1:length(base_string)
                            ibase_string = base_string(i);
                            V_n = ibase_string+"_"+string(j);
                            if strcontain(V_n ,strvar_list)
                                delta = sym('delta','real');
                                Coeffs_tmp =  (Rnn(j) - RdL(base_num_min(i)));
                                V_subs = base_symvar_min(i)*exp(-Coeffs_tmp/delta);
                                H_hr.HcoeL = subs(H_hr.HcoeL,sym(V_n),V_subs);
                                if H_hr.overlap
                                    H_hr.ScoeL = subs(H_hr.ScoeL,sym(V_n),V_subs);
                                end
                            end
                        end
                    end
                    if H_hr.overlap
                        for j = 2:level_cut
                            for i= 1:length(base_string_S)
                                ibase_string_S = base_string_S(i);
                                S_n = ibase_string_S+"_"+string(j);
                                if strcontain(S_n ,strvar_list_S)
                                    delta = sym('delta__2','real');
                                    Coeffs_tmp =  (Rnn(j) - RdL(base_num_min_S(i)));
                                    S_subs = base_num_min_S(i)*exp(-Coeffs_tmp/delta);
                                    H_hr.ScoeL = subs(H_hr.ScoeL,sym(S_n),S_subs);
                                end
                            end
                        end
                    end
                case 2
                    for j = 2:level_cut
                        for i= 1:length(base_string)
                            ibase_string = base_string(i);
                            V_n = ibase_string+"_"+string(j);
                            if strcontain(V_n ,strvar_list)
                                delta = sym(['delta_',num2str(i)],'real');
                                Coeffs_tmp =  (Rnn(j) - RdL(base_num_min(i)));
                                V_subs = base_symvar_min(i)*exp(-Coeffs_tmp/delta);
                                H_hr.HcoeL = subs(H_hr.HcoeL,sym(V_n),V_subs);
                            end
                        end
                    end
                    if H_hr.overlap
                        for j = 2:level_cut
                            for i= 1:length(base_string_S)
                                ibase_string_S = base_string_S(i);
                                S_n = ibase_string_S+"_"+string(j);
                                if strcontain(S_n ,strvar_list_S)
                                    delta = sym(['delta__2_',num2str(i)],'real');
                                    Coeffs_tmp =  (Rnn(j) - RdL(base_num_min_S(i)));
                                    S_subs = base_num_min_S(i)*exp(-Coeffs_tmp/delta);
                                    H_hr.ScoeL = subs(H_hr.ScoeL,sym(S_n),S_subs);
                                end
                            end
                        end
                    end
            end
        end
        function H_hr = alpharule(H_hr,level_cut,mode,options)
            arguments
                H_hr HR;
                level_cut double{mustBeInteger} = -1;
                mode double = 0;
                options.Rd double = -1;
                options.silence = true;
            end
            if options.Rd == -1
                Auto_Rd = true;
            else
                Auto_Rd = false;
            end
            if mode == 0
                return;
            end
            if level_cut == -1
                level_cut = length(H_hr.Rnn);
            end
            % init
            Rnn = H_hr.nn_information(options.silence);
            base_string = ["VssS","VspS","VsdS","VppS","VpdS","VppP","VpdP","VddS","VddP","VddD"];
            strvar_list = string(H_hr.symvar_list);
            base_symvar_min = sym([]);
            base_num_min = ones(1,length(base_string));
            count = 0;
            for ibase_string = base_string
                for i = 1:level_cut
                    if contains(ibase_string+"_"+string(i),strvar_list)
                        count = count +1;
                        base_symvar_min(count) = sym(ibase_string+"_"+string(i),'real');
                        base_num_min(count) = i;
                        base_string(count) = ibase_string;
                        break;
                    end
                end
            end
            base_string(count+1:end) = [];
            base_num_min(count+1:end) = [];
            if H_hr.overlap
                base_string_S = ["SssS","SspS","SsdS","SppS","SpdS","SppP","SpdP","SddS","SddP","SddD"];
                symvar_list_S = symvar(H_hr.ScoeL);
                strvar_list_S  = [string(symvar_list_S)];
                base_symvar_min_S = sym([]);
                base_num_min_S = ones(1,length(base_string));
                count_S = 0;
                for ibase_string_S = base_string_S
                    for i = 1:level_cut
                        if contains(ibase_string_S+"_"+string(i),strvar_list_S)
                            count_S = count_S +1;
                            base_symvar_min_S(count_S) = sym(ibase_string_S+"_"+string(i),'real');
                            base_num_min_S(count_S) = i;
                            base_string_S(count_S) = ibase_string_S;
                            break;
                        end
                    end
                end
                base_string_S(count_S+1:end) = [];
                base_num_min_S(count_S+1:end) = [];
            end
            %             debug =true;
            %             if debug
            %             disp(base_num_min);
            %             disp(base_num_min_S);
            %             end
            N_base_string = length(base_string);
            V_n_str_L = repmat(base_string.',[level_cut-1,1]);
            level_cut_str_L = string(kron((2:level_cut).',ones(N_base_string,1)));
            V_n_list = sym(strcat(V_n_str_L,repmat('_',[(level_cut-1)*N_base_string,1]),level_cut_str_L),'real');
            V_subs_list = sym(zeros((level_cut-1)*N_base_string,1));
            if H_hr.overlap
                N_base_string_S = length(base_string_S);
                S_n_str_L = repmat(base_string_S.',[level_cut-1,1]);
                S_n_list = sym(strcat(S_n_str_L,repmat('_',[(level_cut-1)*N_base_string_S,1]),level_cut_str_L),'real');
                S_subs_list = sym(zeros((level_cut-1)*N_base_string_S,1));
            end
            switch mode
                case 0
                    return;
                case 1
                    alpha = sym('alpha');
                    if Auto_Rd
                        for n = 1:(level_cut-1)*N_base_string
                            [j,i] = ind2sub([level_cut-1,N_base_string],n);
                            j = j+1;
                            Coeffs_tmp =  (Rnn(j) - Rnn(base_num_min(i)))/Rnn(base_num_min(i));
                            V_subs_list(n) = base_symvar_min(i)*exp(-Coeffs_tmp*alpha);
                        end
                        if H_hr.overlap
                            alpha = sym('alpha__2');
                            for n = 1:(level_cut-1)*N_base_string
                                [j,i] = ind2sub([level_cut-1,N_base_string],n);
                                j = j+1;
                                Coeffs_tmp =  (Rnn(j) - Rnn(base_num_min_S(i)))/Rnn(base_num_min_S(i));
                                S_subs_list(n) = base_symvar_min_S(i)*exp(-Coeffs_tmp*alpha);
                            end
                        end
                    else
                        for n = 1:(level_cut-1)*N_base_string
                            [j,i] = ind2sub([level_cut-1,N_base_string],n);
                            j = j+1;
                            Coeffs_tmp =  (Rnn(j) - options.Rd)/options.Rd;
                            V_subs_list(n) = base_symvar_min(i)*exp(-Coeffs_tmp*alpha);
                        end
                        if H_hr.overlap
                            alpha = sym('alpha__2');
                            for n = 1:(level_cut-1)*N_base_string
                                [j,i] = ind2sub([level_cut-1,N_base_string],n);
                                j = j+1;
                                Coeffs_tmp =  (Rnn(j) - options.Rd)/options.Rd;
                                S_subs_list(n) = base_symvar_min_S(i)*exp(-Coeffs_tmp*alpha);
                            end
                        end
                    end
                case 2
                    if Auto_Rd
                        Rd = Rnn(base_num_min);
                        for k = 1:length(base_symvar_min)
                            if isa(Rd(k),'sym')
                                fprintf('%s : %s \\AA\n',base_symvar_min(k),string(Rd(k)));
                            else
                                fprintf('%s : %5.3f \\AA\n',base_symvar_min(k),Rd(k));
                            end
                        end
                        for n = 1:(level_cut-1)*N_base_string
                            [j,i] = ind2sub([level_cut-1,N_base_string],n);
                            j = j+1;
                            alpha = sym(['alpha_',num2str(i),'_', char(base_symvar_min(i))]);
                            Coeffs_tmp =  (Rnn(j) - Rnn(base_num_min(i)))/Rnn(base_num_min(i));
                            %Coeffs_tmp =  (Rnn(j)/Rnn(base_num_min(i))-1);
                            %disp(Coeffs_tmp)
                            V_subs_list(n) = base_symvar_min(i)*exp(-Coeffs_tmp*alpha);
                        end
                        if H_hr.overlap
                            Rd = Rnn(base_num_min_S);
                            for k = 1:length(base_symvar_min_S)
                                if isa(Rd(k),'sym')
                                    fprintf('%s : %s \\AA\n',base_symvar_min_S(k),string(Rd(k)));
                                else
                                    fprintf('%s : %5.3f \\AA\n',base_symvar_min_S(k),Rd(k));
                                end
                            end
                            for n = 1:(level_cut-1)*N_base_string
                                [j,i] = ind2sub([level_cut-1,N_base_string],n);
                                j = j+1;
                                alpha = sym(['alpha__2_',num2str(i),'_', char(base_symvar_min_S(i))]);
                                Coeffs_tmp =  (Rnn(j) - Rnn(base_num_min_S(i)))/Rnn(base_num_min_S(i));
                                S_subs_list(n) = base_symvar_min_S(i)*exp(-Coeffs_tmp*alpha);
                            end
                        end
                    else
                        Rd = options.Rd;
                        if length(Rd) == 1
                            Rd = repmat(Rd,[N_base_string 1]);
                        end
                        for k = 1:length(base_symvar_min)
                            if isa(Rd(k),'sym')
                                fprintf('%s : %s \\AA\n',base_symvar_min(k),string(Rd(k)));
                            else
                                fprintf('%s : %5.3f \\AA\n',base_symvar_min(k),Rd(k));
                            end
                        end
                        for n = 1:(level_cut-1)*N_base_string
                            [j,i] = ind2sub([level_cut-1,N_base_string],n);
                            j = j+1;
                            alpha = sym(['alpha_',num2str(i),'_', char(base_symvar_min(i))]);
                            Coeffs_tmp =  (Rnn(j) - Rd(i))/Rd(i);
                            V_subs_list(n) = base_symvar_min(i)*exp(-Coeffs_tmp*alpha);
                        end
                        if H_hr.overlap
                            for k = 1:length(base_symvar_min_S)
                                if isa(Rd(k),'sym')
                                    fprintf('%s : %s \\AA\n',base_symvar_min_S(k),string(Rd(k)));
                                else
                                    fprintf('%s : %5.3f \\AA\n',base_symvar_min_S(k),Rd(k));
                                end
                            end
                            for n = 1:(level_cut-1)*N_base_string
                                [j,i] = ind2sub([level_cut-1,N_base_string],n);
                                j = j+1;
                                alpha = sym(['alpha__2_',num2str(i),'_', char(base_symvar_min_S(i))]);
                                Coeffs_tmp =  (Rnn(j) -Rd(i))/Rd(i);
                                S_subs_list(n) = base_symvar_min_S(i)*exp(-Coeffs_tmp*alpha);
                            end
                        end
                    end
            end
            H_hr.HcoeL = subs(H_hr.HcoeL,V_n_list,V_subs_list);
            if H_hr.overlap
                H_hr.ScoeL = subs(H_hr.ScoeL,S_n_list,S_subs_list);
            end
        end
    end
    % ----------------  expand --------------------
    methods
        function H_hr = cut_piece(H_hr,repeatnum,fin_dir,glue_edges,vacuum_mode)
            %         Constructs a (d-1)-dimensional tight-binding model out of a d-dimensional one by repeating the unit cell a given number of times along one of the periodic lattice vectors.
            %          The real-space lattice vectors of the returned model are the same as those of
            %          the original model; only the dimensionality of reciprocal space is reduced.
            %
            %         :param repeatnum: How many times to repeat the unit cell.
            %
            %         :param fin_dir: Index of the real space lattice vector along
            %           which you no longer wish to maintain periodicity.
            
            %         :param orbital_init
            
            %         :param glue_edges: Optional boolean parameter specifying whether to
            %           allow hoppings from one edge to the other of a cut model.
            %
            %         :returns:
            %           * **fin_model** -- Object of type
            %             Out_H_xyz representing a cutout
            %             tight-binding model.
            %             Orbitals in *fin_model* are
            %             numbered so that the i-th orbital of the n-th unit
            %             cell has index i+norb*n (here norb is the number of
            %             orbitals in the original model).
            %
            %         Example usage::
            %           Out_H_xyz = cut_piece(H_xyz,repeatnum,fin_dir,glue_edgs,orbital_init)
            
            %           A = tb_model(3, 3, ...)
            %           % Construct two-dimensional model B out of three-dimensional
            %           % model A by repeating model along second lattice vector ten times
            %           B = A.cut_piece(10, 1)
            %           % Further cut two-dimensional model B into one-dimensional model
            %           % A by repeating unit cell twenty times along third lattice
            %           % vector and allow hoppings from one edge to the other
            %           C = B.cut_piece(20, 2, glue_edgs=True)
            
            %  ---- nargin ----
            arguments
                H_hr HR;
                repeatnum double{mustBeInteger} =  10;
                fin_dir double{mustBeMember(fin_dir,[1,2,3])} = 3;
                glue_edges logical = false;
                vacuum_mode logical = false;
            end
            Ns = [1 0 0;0 1 0;0 0 1];
            Ns(fin_dir,:) = Ns(fin_dir,:) * repeatnum;
            fin_dir_list = [0 0 0];
            %[sc_orb,~] = H_hr.supercell_orb(Ns);
            [sc_orb,~,sc_elementL,sc_quantumL] = H_hr.supercell_orb(Ns);
            sc_orb = double(sc_orb);
            % -- check POSCAR --
            if vacuum_mode
                fin_dir_list(fin_dir) = 1;
            end
            if isempty(H_hr.sites)
                if exist('POSCAR','file')
                    [Rm_tmp,sites_tmp,Atom_name_tmp,Atom_num_tmp,~]=HR.POSCAR_readin('POSCAR','vasp');
                    H_hr = H_hr.supercell(Ns,'POSCAR_super_fin',Rm_tmp,sites_tmp,Atom_name_tmp,Atom_num_tmp,fin_dir_list);
                else
                    % interaction
                end
            else
                try
                    H_hr = H_hr.supercell(Ns,'POSCAR_super_fin',H_hr.Rm,H_hr.sites,H_hr.Atom_name,H_hr.Atom_num,fin_dir_list);
                catch
                end
            end
            % check value of orbital_init
            norb =size(H_hr.orbL,1);
            if norb ~= H_hr.WAN_NUM
                error("\n\nOribital_init is wrong,please give a right orbital init or just forget this parm!");
            end
            % check value of num
            if repeatnum<1
                error("\n\nArgument num must be positive!");
            end
            if repeatnum == 1 && glue_edges
                error("\n\nCan't have num==1 and glueing of the edges!");
            end
            if vacuum_mode
                % rebuild fin_orb
                Rmlength1 = norm (H_hr.Rm(1,:));
                Rmlength2 = norm (H_hr.Rm(2,:));
                Rmlength3 = norm (H_hr.Rm(3,:));
                Rm_s_fin_add = [10*H_hr.Rm(1,:)*fin_dir_list(1)/Rmlength1;...
                    10*H_hr.Rm(2,:)*fin_dir_list(2)/Rmlength2;...
                    10*H_hr.Rm(3,:)*fin_dir_list(3)/Rmlength3];
                Rm_s_fin = H_hr.Rm + Rm_s_fin_add ;
                Rc_s_fin_add = [1/2, 1/2 ,1/2] ;
                Rr_s_fin_add = Rc_s_fin_add * Rm_s_fin_add;
                [nfinorb,~ ]= size(sc_orb);
                for  i = 1:nfinorb
                    Rr_orb = sc_orb(i,:)*H_hr.Rm;
                    Rr_s_fin = Rr_orb + Rr_s_fin_add;
                    Rc_s_fin = Rr_s_fin /Rm_s_fin;
                    sc_orb(i,:) = Rc_s_fin ;
                end
            end
            % generate periodic directions of a finite model
            OUT_WAN_NUM = H_hr.WAN_NUM*repeatnum ;
            %OUT_WAN_NUM_2 = OUT_WAN_NUM^2 ;
            % init
            vectorList = double(H_hr.vectorL);
            NRPT_seq = false(H_hr.NRPTS,1);
            if H_hr.overlap
                vectorList_overlap = double(H_hr.vectorL_overlap);
                NRPTS_S = size(vectorList_overlap,1);
                NRPT_seq_S = false(NRPTS_S,1);
            end
            WANNUM = H_hr.WAN_NUM;
            if strcmp(H_hr.Type,'sparse')
                OUT_HnumL{H_hr.NRPTS} = sparse(OUT_WAN_NUM,OUT_WAN_NUM);
                for i =1: H_hr.NRPTS-1
                    OUT_HnumL{i} = sparse(OUT_WAN_NUM,OUT_WAN_NUM);
                end
                pb = vasplib_tool_outer.CmdLineProgressBar('NRPT : ');
                for iN = 1:H_hr.NRPTS
                    % lattice vector of the hopping
                    ind_R = vectorList(iN,:);
                    jump_fin=ind_R(fin_dir);
                    % store by how many cells is the hopping in finite direction
                    %temp_Hnum =sparse(OUT_WAN_NUM,OUT_WAN_NUM);
                    pb.print(iN,H_hr.NRPTS);
                    % speed up more and more !!
                    % go over all NRPTS
                    [ilist,jlist,amplist] = find(H_hr.HnumL{iN});
                    %
                    nhopping = length(amplist);
                    for ih = 1:nhopping
                        i = ilist(ih);
                        j = jlist(ih);
                        amp = amplist(ih);
                        for icur_sc_vec = 1:repeatnum % go over all cells in finite direction mini
                            hi= i + (icur_sc_vec-1)*WANNUM ;
                            %disp(hi);
                            hj= j + (icur_sc_vec+jump_fin-1)*WANNUM ;
                            %disp(hj);
                            % decide whether this hopping should be added or not
                            to_add=1;
                            %disp(hj);
                            % if edges are not glued then neglect all jumps that spill out
                            if ~glue_edges
                                if hj <= 0 || hj > OUT_WAN_NUM
                                    to_add=0;
                                end
                                % if edges are glued then do mod division to wrap up the hopping
                            else
                                hj= mod(hj,OUT_WAN_NUM);
                                if hj ==0
                                    hj = OUT_WAN_NUM;
                                end
                            end
                            if to_add == 1
                                OUT_HnumL{iN}(hi,hj) = amp;
                                %OUT_Hnum_list(hi,hj,ih) = amp;
                                %OUT_Hnum_list(hi,hj,i) = amp;
                            end
                        end
                    end
                    %OUT_Hnum_list(:,ih) = temp_Hnum;
                end
                pb.delete();
                %OUT_Hnum_list(:,:,ih) = temp_Hnum;
                %
                H_hr.HnumL = OUT_HnumL;
                H_hr.orbL = sc_orb    ;
                H_hr.Type = 'sparse'     ;
            elseif strcmp(H_hr.Type,'mat')
                OUT_HnumL = zeros(OUT_WAN_NUM,OUT_WAN_NUM,H_hr.NRPTS);
                NRPTS_record = H_hr.NRPTS;
                if H_hr.overlap
                    OUT_SnumL = zeros(OUT_WAN_NUM,OUT_WAN_NUM,H_hr.NRPTS);
                    % wrong
                end
                pb = vasplib_tool_outer.CmdLineProgressBar('NRPT : ');
                for ih = 1:H_hr.NRPTS % go over all NRPTS
                    % lattice vector of the hopping
                    ind_R = vectorList(ih,:);
                    jump_fin=ind_R(fin_dir);
                    tmpHnum = H_hr.HnumL(:,:,ih);
                    Nonzero_list = find(tmpHnum ~= 0);
                    [hi,hj] = ind2sub([WANNUM,WANNUM],Nonzero_list);
                    pb.print(ih,H_hr.NRPTS);
                    HI_list = zeros(1,OUT_WAN_NUM*OUT_WAN_NUM);
                    HJ_list = zeros(1,OUT_WAN_NUM*OUT_WAN_NUM);
                    CONTAIN_list = zeros(1,OUT_WAN_NUM*OUT_WAN_NUM);
                    rm_count  = 1;
                    for icur_sc_vec = 1:repeatnum %
                        ind_R_tmp = ind_R;
                        hi_list = hi + (icur_sc_vec-1)*WANNUM ;
                        hj_list = hj + (icur_sc_vec+jump_fin-1)*WANNUM ;
                        if ~glue_edges
                            Contain_list = find(hj_list > 0 & hj_list <= OUT_WAN_NUM);
                            hj_list = hj_list(Contain_list);
                            hi_list = hi_list(Contain_list);
                            ind_R_tmp(fin_dir) = 0;
                        else
                            % if edges are glued then do mod division to wrap up the hopping
                            hj_list= mod(hj_list,OUT_WAN_NUM);
                            ind_R_tmp(fin_dir) = floor(hj./OUT_WAN_NUM);
                            hj_list(hj_list ==0) = OUT_WAN_NUM;
                            error('There is a bug, contact me: parkman@buaa.edu.cn');
                        end
                        [~,IH] = ismember(ind_R_tmp,vectorList,'rows');
                        % store by how many cells is the hopping in finite direction
                        %temp_Hnum =sparse(OUT_WAN_NUM,OUT_WAN_NUM);
                        nContain_list = length(Contain_list) ;
                        HI_list(rm_count:rm_count +nContain_list-1) = hi_list ;
                        HJ_list(rm_count:rm_count +nContain_list-1 ) = hj_list;
                        CONTAIN_list(rm_count:rm_count +nContain_list-1 ) = Contain_list;
                        rm_count = rm_count + nContain_list;
                    end
                    
                    HI_list(rm_count:end) = [];
                    HJ_list(rm_count:end) = [];
                    CONTAIN_list(rm_count:end) = [];
                    if IH == 0
                        IH = NRPTS_record+1;
                        NRPTS_record = NRPTS_record+1;
                        vectorList(IH,:) = ind_R_tmp;
                        OUT_HnumL(:,:,IH) = full(sparse(HI_list,HJ_list,tmpHnum(Nonzero_list(CONTAIN_list)),OUT_WAN_NUM,OUT_WAN_NUM));
                    else
                        OUT_HnumL(:,:,IH) = OUT_HnumL(:,:,IH) + full(sparse(HI_list,HJ_list,tmpHnum(Nonzero_list(CONTAIN_list)),OUT_WAN_NUM,OUT_WAN_NUM));
                    end
                    
                    NRPT_seq(IH)  = true;
                    %OUT_Hnum_list(:,ih) = temp_Hnum;
                end
                
                pb.delete();
                
                if H_hr.overlap
                    pb = vasplib_tool_outer.CmdLineProgressBar('NRPT_S : ');
                    for is = 1:NRPTS_S % go over all NRPT_S
                        % lattice vector of the hopping
                        ind_R_S = vectorList_overlap(is,:);
                        jump_fin=ind_R_S(fin_dir);
                        tmpSnum = H_hr.SnumL(:,:,is);
                        Nonzero_list = find(tmpSnum ~= 0);
                        [si,sj] = ind2sub([WANNUM,WANNUM],Nonzero_list);
                        % store by how many cells is the hopping in finite direction
                        %temp_Hnum =sparse(OUT_WAN_NUM,OUT_WAN_NUM);
                        pb.print(is,NRPTS_S);
                        SI_list = zeros(1,OUT_WAN_NUM*OUT_WAN_NUM);
                        SJ_list = zeros(1,OUT_WAN_NUM*OUT_WAN_NUM);
                        CONTAIN_list = zeros(1,OUT_WAN_NUM*OUT_WAN_NUM);
                        rm_count  = 1;
                        % speed up more and more !!
                        for icur_sc_vec = 1:repeatnum % go over all cells in finite direction mini
                            ind_R_tmp  = ind_R_S;
                            si_list = si + (icur_sc_vec-1)*WANNUM ;
                            sj_list = sj + (icur_sc_vec+jump_fin-1)*WANNUM ;
                            % decide whether this hopping should be added or not
                            if ~glue_edges
                                Contain_list = find(sj_list > 0 & sj_list <= OUT_WAN_NUM);
                                sj_list = sj_list(Contain_list) ;
                                si_list = si_list(Contain_list) ;
                                ind_R_tmp(fin_dir) = 0;
                                % if edges are glued then do mod division to wrap up the hopping
                            else
                                sj_list= mod(sj_list,OUT_WAN_NUM);
                                ind_R_tmp(fin_dir) = floor(hj./OUT_WAN_NUM);
                                sj_list(sj_list ==0) = OUT_WAN_NUM;
                                error('There is a bug, contact me: parkman@buaa.edu.cn');
                            end
                            [~,IS] = ismember(ind_R_tmp,vectorList_overlap,'rows');
                            NRPT_seq_S(IS) = true;
                            nContain_list = length(Contain_list) ;
                            SI_list(rm_count:rm_count +nContain_list-1) = si_list ;
                            SJ_list(rm_count:rm_count +nContain_list-1 ) = sj_list;
                            CONTAIN_list(rm_count:rm_count +nContain_list-1 ) = Contain_list;
                            rm_count = rm_count + nContain_list;
                        end
                        SI_list(rm_count:end) = [];
                        SJ_list(rm_count:end) = [];
                        CONTAIN_list(rm_count:end) = [];
                        OUT_SnumL(:,:,IS)= OUT_SnumL(:,:,IS)+full(sparse(SI_list,SJ_list,tmpSnum(Nonzero_list(CONTAIN_list)),OUT_WAN_NUM,OUT_WAN_NUM));
                    end
                end
                
                %OUT_Hnum_list(:,:,ih) = temp_Hnum;
                %
                H_hr.vectorL = int8(vectorList);
                H_hr.HnumL = OUT_HnumL;
                if H_hr.overlap
                    H_hr.SnumL = OUT_SnumL;
                end
                if H_hr.overlap
                    H_hr = H_hr.reseq(':',NRPT_seq,NRPT_seq_S);
                else
                    H_hr = H_hr.reseq(':',NRPT_seq);
                end
            elseif strcmp(H_hr.Type,'list')
                OUT_HnumL = zeros(H_hr.NRPTS*repeatnum,1);
                vector_list = H_hr.vectorL;
                OUT_vectorList = zeros(H_hr.NRPTS*repeatnum,5);
                pb = vasplib_tool_outer.CmdLineProgressBar('H:NRPT : ');
                count = 0;
                for ih = 1:H_hr.NRPTS % go over all NRPTS
                    % lattice vector of the hopping
                    ind_R = vectorList(ih,1:3);
                    jump_fin=ind_R(fin_dir);
                    % store by how many cells is the hopping in finite direction
                    %temp_Hnum =sparse(OUT_WAN_NUM,OUT_WAN_NUM);
                    pb.print(ih,H_hr.NRPTS);
                    % speed up more and more !!
                    i = vector_list(ih,4);
                    j = vector_list(ih,5);
                    % amplitude of the hop is the same
                    amp = H_hr.HnumL(ih);
                    if norm(amp) > 0
                        for icur_sc_vec = 1:repeatnum % go over all cells in finite direction mini
                            hi= i + (icur_sc_vec-1)*WANNUM ;
                            %disp(hi);
                            hj= j + (icur_sc_vec+jump_fin-1)*WANNUM ;
                            %disp(hj);
                            % decide whether this hopping should be added or not
                            to_add=1;
                            %disp(hj);
                            % if edges are not glued then neglect all jumps that spill out
                            if ~glue_edges
                                if hj <= 0 || hj > OUT_WAN_NUM
                                    to_add=0;
                                end
                                % if edges are glued then do mod division to wrap up the hopping
                            else
                                hj= mod(hj,OUT_WAN_NUM);
                                if hj ==0
                                    hj = OUT_WAN_NUM;
                                end
                            end
                            if to_add == 1
                                count = count+1;
                                OUT_HnumL(count) = amp;
                                OUT_vectorList(count,:) = [ind_R,hi,hj];
                                %OUT_Hnum_list(hi,hj,ih) = amp;
                                %                        OUT_Hnum_list(hi,hj,i) = amp;
                            end
                        end
                    end
                end
                pb.delete();
                fprintf('remove %d empty orbs',H_hr.NRPTS*repeatnum- count);
                OUT_HnumL(count+1:end,:) = [];
                OUT_vectorList(count+1:end,:) = [];
                %OUT_Hnum_list(:,:,ih) = temp_Hnum;
                %
                H_hr.HnumL = OUT_HnumL;
                OUT_vectorList(:,fin_dir) = 0;
                H_hr.vectorL = OUT_vectorList;
                if H_hr.overlap
                    NRPTS_S = size(H_hr.vectorL_overlap,2);
                    vectorList_overlap = H_hr.vectorL_overlap;
                    OUT_vectorList_overlap = zeros(NRPTS_S*repeatnum,5);
                    OUT_SnumL = zeros(NRPTS_S*repeatnum,1);
                    pb = vasplib_tool_outer.CmdLineProgressBar('S:NRPT : ');
                    count = 0;
                    for ih = 1:NRPTS_S % go over all NRPTS
                        % lattice vector of the hopping
                        ind_R = vectorList_overlap(ih,1:3);
                        jump_fin=ind_R(fin_dir);
                        % store by how many cells is the hopping in finite direction
                        %temp_Hnum =sparse(OUT_WAN_NUM,OUT_WAN_NUM);
                        pb.print(ih,H_hr.NRPTS);
                        % speed up more and more !!
                        i = vectorList_overlap(ih,4);
                        j = vectorList_overlap(ih,5);
                        % amplitude of the hop is the same
                        amp = H_hr.SnumL(ih);
                        if norm(amp) > 0
                            for icur_sc_vec = 1:repeatnum % go over all cells in finite direction mini
                                si= i + (icur_sc_vec-1)*WANNUM ;
                                %disp(hi);
                                sj= j + (icur_sc_vec+jump_fin-1)*WANNUM ;
                                %disp(hj);
                                % decide whether this hopping should be added or not
                                to_add=1;
                                %disp(hj);
                                % if edges are not glued then neglect all jumps that spill out
                                if ~glue_edges
                                    if sj <= 0 || sj > OUT_WAN_NUM
                                        to_add=0;
                                    end
                                    % if edges are glued then do mod division to wrap up the hopping
                                else
                                    sj= mod(sj,OUT_WAN_NUM);
                                    if sj ==0
                                        sj = OUT_WAN_NUM;
                                    end
                                end
                                if to_add == 1
                                    count = count+1;
                                    OUT_SnumL(count) = amp;
                                    OUT_vectorList_overlap(count,:) = [ind_R,si,sj];
                                    %OUT_Hnum_list(hi,hj,ih) = amp;
                                    %                        OUT_Hnum_list(hi,hj,i) = amp;
                                end
                            end
                        end
                    end
                    pb.delete();
                    H_hr.SnumL = OUT_SnumL;
                    OUT_vectorList_overlap(:,fin_dir) = 0;
                    H_hr.vectorL_overlap = OUT_vectorList_overlap;
                end
            end
            H_hr.num = true    ;
            H_hr.coe = false    ;
            H_hr.orbL = sc_orb  ;
            H_hr.quantumL = sc_quantumL; % update orbL
            H_hr.elementL = sc_elementL; % update orbL
        end
        function H_hr = supercell_hr(H_hr,Ns,options)
            % Returns tight-binding model representing a super-cell of a current object.
            % This function can be used together with *cut_piece* in order to create slabs with arbitrary surfaces.
            %
            % * Label: play_hr
            %
            %--------  init  --------
            %import vasplib_tool.*
            arguments
                H_hr HR;
                Ns double = eye(3);
                options.Accuracy double = 1e-6;
                options.force_list = false;
            end
            %--------  init  --------
            V= abs(round(det(Ns)));% intger forcely
            OUT_WAN_NUM = H_hr.WAN_NUM*V ;
            WANNUM= H_hr.WAN_NUM;
            Accuracy = options.Accuracy;
            pc_orb = H_hr.orbL;
            %--------  check  --------
            % checks on super-lattice Ns
            [sc_orb,sc_vec,sc_elementL,sc_quantumL] = H_hr.supercell_orb(Ns,Accuracy);
            % create super-cell tb_model object to be returned
            %OUT_H_hr_ = struct('seq',[],'vector',[],'Degen',[],'key',[],'nokey',[],'Hstr',[],'Hsym',[],'Hcoe',[],'Hnum',[]);
            % OUT_H_xyz = repmat(H_xyz_ ,[NRPTS,1]);    % Hamiltonian of every cell;
            OUT_H_hr = H_hr;
            OUT_H_hr = OUT_H_hr.clean(OUT_WAN_NUM);
            OUT_H_hr.orbL = sc_orb; % update orbL
            OUT_H_hr.quantumL = sc_quantumL; % update orbL
            OUT_H_hr.elementL = sc_elementL; % update orbL
            OUT_H_hr.Rm = Ns * OUT_H_hr.Rm;
            NRPTS_  = H_hr.NRPTS;
            if H_hr.overlap
                NRPTS_S = size(H_hr.vectorL_overlap,1);
            end
            fprintf('Search done; begin set hoppings\n');
            fprintf('We can improve the perfomance later\n');
            %             t1=clock;
            num_sc = size(sc_vec,1);
            %set hopping terms
            Accuracy_roundn = round(log(Accuracy)/log(10));
            % define orb_sc_L
            pc_orb_super = repmat(pc_orb,[num_sc,1]);
            orb_sc_vevL = sc_orb*Ns - pc_orb_super;
            Leak_L = ~ismember(orb_sc_vevL,sc_vec,'rows');
            orb_sc_labelL = find(Leak_L,1);
            %Leak_orb_sc_vecL = orb_sc_vevL(Leak_L,:);
            if isempty(orb_sc_labelL) && strcmp(H_hr.Type,'mat') ||options.force_list
                fprintf("supercell has leak sites, enforce list mode!");
                H_hr = H_hr.rewrite();
                OUT_H_hr = OUT_H_hr.rewrite();
                %options.force_list = true;
            end
            switch H_hr.Type
                case 'mat'
                    % bug fix
                    % hj need double check
                    % bug still exist, there is one situation not
                    % considered!
                    for icur_sc_vec = 1:num_sc % go over all super-cell vectors
                        cur_sc_vec = double(sc_vec(icur_sc_vec,:));
                        pb = vasplib_tool_outer.CmdLineProgressBar(...
                            ['Generate process: SUPERCELL(',...
                            num2str(icur_sc_vec),',',num2str(num_sc),') NRPT:']);
                        for ih = 1:NRPTS_ % go over all hopping terms of the original model
                            % lattice vector of the hopping
                            ind_R = double(H_hr.vectorL(ih,:));
                            % super-cell component of hopping lattice vector
                            % shift also by current super cell vector
                            indR_in_supercell=double(ind_R+cur_sc_vec)/Ns;
                            % fix num bug
                            indR_in_supercell = roundn(indR_in_supercell,Accuracy_roundn);
                            sc_part=floor(indR_in_supercell); % round down!
                            %disp(sc_part);
                            %find remaining vector in the original reduced coordinates
                            orig_part=ind_R+cur_sc_vec-double(sc_part*Ns);
                            %disp(orig_part);
                            % remaining vector must equal one of the super-cell vectors
                            pair_ind = 9999;
                            for jcur_sc_vec = 1:num_sc
                                pair_sc_vec = sc_vec(jcur_sc_vec,:);
                                if pair_sc_vec==orig_part
                                    if pair_ind ~= 9999
                                        error("\n\nFound duplicate super cell vector!");
                                    end
                                    pair_ind=jcur_sc_vec;
                                end
                            end
                            if pair_ind==9999
                                disp(orig_part);
                                disp('Cant find sc in ');
                                disp(sc_vec);
                                disp('orig_part=ind_R+cur_sc_vec-sc_part*Ns;');
                                disp(ind_R);
                                disp(cur_sc_vec);
                                %continue;
                                % The bug from numerical issue! we will use sym instead
                                error("\n\nDid not find super cell vector!");
                            end
                            %disp(pair_ind);
                            % index of "from" and "to" hopping indices
                            if H_hr.num
                                tmpHnum = H_hr.HnumL(:,:,ih);
                                Nonzero_list = find(tmpHnum ~= 0);
                                %N_Nonzero_list = numel(Nonzero_list);
                                [hi,hj] = ind2sub([WANNUM,WANNUM],Nonzero_list);
                                hi= hi + (icur_sc_vec-1)*H_hr.WAN_NUM ;
                                hj= hj + (pair_ind-1)*H_hr.WAN_NUM ;
                                tmp_mat = full(sparse(hi,hj,tmpHnum(Nonzero_list),OUT_WAN_NUM,OUT_WAN_NUM));
                                % deal with leak orbitals
                                %                                 for i = 1:length(orb_sc_labelL)
                                %                                     tmp_mat_leak = zeros(size(tmp_mat));
                                %                                     tmp_mat_leak(orb_sc_labelL(i),:) =  tmp_mat_leak(orb_sc_labelL(i),:);
                                %                                     tmp_mat_leak(:,orb_sc_labelL(i)) =  tmp_mat_leak(:,orb_sc_labelL(i));
                                %                                     indR_in_supercell_leak=double(ind_R+Leak_orb_sc_vecL(i,:))/Ns;
                                %                                     % fix num bug
                                %                                     indR_in_supercell_leak = roundn(indR_in_supercell_leak,Accuracy_roundn);
                                %                                     sc_part_leak=floor(indR_in_supercell_leak); % round down!
                                %                                     OUT_H_hr = OUT_H_hr.set_hop_mat(tmp_mat_leak,sc_part_leak,'add');
                                %                                 end
                                %                                 tmp_mat(orb_sc_labelL,:) = 0;
                                %
                                OUT_H_hr = OUT_H_hr.set_hop_mat(tmp_mat,sc_part,'add');
                            end
                            if H_hr.coe
                                tmpHcoe = H_hr.HcoeL(:,:,ih);
                                Nonzero_list = find(tmpHcoe ~= sym(0));
                                [hi,hj] = ind2sub([WANNUM,WANNUM],Nonzero_list);
                                hi= hi + (icur_sc_vec-1)*H_hr.WAN_NUM ;
                                hj= hj + (pair_ind-1)*H_hr.WAN_NUM ;
                                tmp_mat = sym(zeros(OUT_WAN_NUM));
                                tmp_mat(sub2ind([OUT_WAN_NUM,OUT_WAN_NUM],hi,hj)) = tmpHcoe(Nonzero_list);
                                % deal with leak orbitals
                                %                                 for i = 1:length(orb_sc_labelL)
                                %                                     tmp_mat_coe = tmp_mat;
                                %                                     checkHcoeL = zeros(size(tmp_mat_coe),class(tmp_mat_coe));
                                %                                     pre_tmp_mat_coe = tmp_mat_coe(:,orb_sc_labelL(i));
                                %                                     tmp_mat_coe(orb_sc_labelL(i),orb_sc_labelL(i)) = sym(0);
                                %                                     checkHcoeL(orb_sc_labelL(i),:) = tmp_mat_coe(orb_sc_labelL(i),:);
                                %                                     checkHcoeL(:,orb_sc_labelL(i)) = pre_tmp_mat_coe;
                                %                                     Nonzero_list_leak = find(checkHcoeL ~= sym(0));
                                %                                     if ~isempty(Nonzero_list_leak)
                                %                                         tmpHcoeLeak = checkHcoeL(Nonzero_list_leak);
                                %                                         [hi_leakL,hj_leakL] = ind2sub([OUT_WAN_NUM,OUT_WAN_NUM],Nonzero_list);
                                %                                         % ----- find Rvector --------
                                %                                         hi_originL = mod((hi_leakL-1),WANNUM)+1;
                                %                                         hj_originL = mod((hj_leakL-1),WANNUM)+1;
                                %                                         indRtiL = pc_orb(hi_originL,:);
                                %                                         indRti_in_supercellL = sc_orb(hi_leakL,:);
                                %                                         real_sc_vecL = indRti_in_supercellL*Ns - indRtiL;
                                %                                         indRtjL = real_sc_vecL+ind_R+pc_orb(hj_originL,:);
                                %                                         sc_part_leakL = (double(indRtjL)/Ns);
                                %                                         sc_part_leakL  = floor(sc_part_leakL);
                                %                                         for ileak = 1:length(Nonzero_list_leak)
                                %                                             OUT_H_hr = OUT_H_hr.set_hop(tmpHcoeLeak(ileak),...
                                %                                                 hi_leakL(ileak),hj_leakL(ileak),sc_part_leakL(ileak,:),'symadd');
                                %                                         end
                                %                                     end
                                %                                 end
                                %tmp_mat(orb_sc_labelL,:) = sym(0);
                                %tmp_mat(:,orb_sc_labelL) = sym(0);
                                %
                                OUT_H_hr = OUT_H_hr.set_hop_mat(tmp_mat,sc_part,'symadd');
                            end
                            pb.print(ih,NRPTS_,' Hopping ...');
                            %                     fprintf("Generate process: SUPERCELL(%d,%d) NRPT(%d,%d) RUNINGTIME: %f s.\n",...
                            %                         icur_sc_vec,num_sc,ih,H_hr.NRPTS,etime(clock,t1));
                        end
                        pb.delete();
                        % overlap
                        if H_hr.overlap
                            pb = vasplib_tool_outer.CmdLineProgressBar(...
                                ['Generate process: SUPERCELL(S mat)(',...
                                num2str(icur_sc_vec),',',num2str(num_sc),') NRPT:']);
                            for is = 1:NRPTS_S
                                % lattice vector of the hopping
                                ind_R_S = double(H_hr.vectorL_overlap(is,:));
                                % super-cell component of hopping lattice vector
                                % shift also by current super cell vector
                                indR_in_supercell=double(ind_R_S+cur_sc_vec)/Ns;
                                % fix num bug
                                indR_in_supercell = roundn(indR_in_supercell,Accuracy_roundn);
                                sc_part=floor(indR_in_supercell); % round down!
                                %disp(sc_part);
                                %find remaining vector in the original reduced coordinates
                                orig_part=ind_R_S+cur_sc_vec-double(sc_part*Ns);
                                %disp(orig_part);
                                % remaining vector must equal one of the super-cell vectors
                                pair_ind = 9999;
                                for jcur_sc_vec = 1:num_sc
                                    pair_sc_vec = sc_vec(jcur_sc_vec,:);
                                    if pair_sc_vec==orig_part
                                        if pair_ind ~= 9999
                                            error("\n\nFound duplicate super cell vector!");
                                        end
                                        pair_ind=jcur_sc_vec;
                                    end
                                end
                                if pair_ind==9999
                                    disp(orig_part);
                                    disp('Cant find sc in ');
                                    disp(sc_vec);
                                    disp('orig_part=ind_R+cur_sc_vec-sc_part*Ns;');
                                    disp(ind_R);
                                    disp(cur_sc_vec);
                                    %continue;
                                    % The bug from numerical issue! we will use sym instead
                                    error("\n\nDid not find super cell vector!");
                                end
                                %disp(pair_ind);
                                % index of "from" and "to" hopping indices
                                if H_hr.num
                                    tmpSnum = H_hr.SnumL(:,:,is);
                                    Nonzero_list = find(tmpSnum ~= 0);
                                    %N_Nonzero_list = numel(Nonzero_list);
                                    [si,sj] = ind2sub([WANNUM,WANNUM],Nonzero_list);
                                    si= si + (icur_sc_vec-1)*H_hr.WAN_NUM ;
                                    sj= sj + (pair_ind-1)*H_hr.WAN_NUM ;
                                    tmp_mat = full(sparse(si,sj,tmpSnum(Nonzero_list),OUT_WAN_NUM,OUT_WAN_NUM));
                                    OUT_H_hr = OUT_H_hr.set_overlap_mat(tmp_mat,sc_part,'add');
                                end
                                if H_hr.coe
                                    tmpScoe = H_hr.ScoeL(:,:,is);
                                    Nonzero_list = find(tmpScoe ~= sym(0));
                                    [si,sj] = ind2sub([WANNUM,WANNUM],Nonzero_list);
                                    si= si + (icur_sc_vec-1)*H_hr.WAN_NUM ;
                                    sj= sj + (pair_ind-1)*H_hr.WAN_NUM ;
                                    tmp_mat = sym(zeros(OUT_WAN_NUM));
                                    tmp_mat(sub2ind([OUT_WAN_NUM,OUT_WAN_NUM],si,sj)) = tmpScoe(Nonzero_list);
                                    OUT_H_hr = OUT_H_hr.set_overlap_mat(tmp_mat,sc_part,'symadd');
                                end
                                pb.print(is,NRPTS_S,' Overlap ...');
                            end
                            pb.delete();
                            %                     fprintf("Generate process: SUPERCELL(%d,%d) NRPT(%d,%d) RUNINGTIME: %f s.\n",...
                            %                         icur_sc_vec,num_sc,ih,H_hr.NRPTS,etime(clock,t1));
                        end
                    end
                case 'list'
                    if H_hr.num
                        HnumList = repmat(H_hr.HnumL,[num_sc,1]);
                        nHopping = length(H_hr.HnumL);
                    end
                    if H_hr.coe
                        HcoeList = repmat(H_hr.HcoeL,[num_sc,1]);
                        nHopping = length(H_hr.HcoeL);
                    end
                    VectorList = double(H_hr.vectorL);
                    OutVectorList = repmat(VectorList,[num_sc,1]);
                    pb = vasplib_tool_outer.CmdLineProgressBar(...
                        'Generate process: SUPERCELL:');
                    for icur_sc_vec = 1:num_sc % go over all super-cell vectors
                        cur_sc_vec = double(sc_vec(icur_sc_vec,:));
                        hiL = VectorList(:,4);
                        hjL = VectorList(:,5);
                        ind_RL = double(VectorList(:,1:3));
                        % ----- find hj in orgin or sc_vec_list --------
                        ind_R_in_supercellL = double(ind_RL+cur_sc_vec)/Ns;
                        ind_R_in_supercellL = roundn(ind_R_in_supercellL,Accuracy_roundn);
                        sc_partL=floor(ind_R_in_supercellL); % round down!
                        orig_partL=ind_RL+cur_sc_vec-double(sc_partL*Ns);
                        [~,pair_indL] = ismember(orig_partL,sc_vec,'rows');
                        % Below line relies on the right sc_orbL with right
                        % sequence.We are lucky, this procedure goes well
                        sc_hjL = hjL+(pair_indL-1)*WANNUM;
                        % ----- ******* --------
                        % ----- find Rvector --------
                        sc_hiL = hiL + (icur_sc_vec-1)*WANNUM;
                        indRtiL = pc_orb(hiL,:);
                        indRti_in_supercellL = sc_orb(sc_hiL,:);
                        real_sc_vecL = indRti_in_supercellL*Ns - indRtiL;
                        real_sc_vecL = round(real_sc_vecL);% ? this numerical error may cause bug!
                        indRtjL = real_sc_vecL+ind_RL+pc_orb(hjL,:);
                        indRtj_in_supercellL = roundn(double(indRtjL)/Ns,Accuracy_roundn);% ? this numerical error may cause bug!
                        indR_in_supercellL  = floor(indRtj_in_supercellL);
                        OutVectorList((icur_sc_vec-1)*nHopping+1:(icur_sc_vec)*nHopping,:) = ...
                            [indR_in_supercellL,sc_hiL,sc_hjL];
                        pb.print(icur_sc_vec,num_sc,' ...');
                    end
                    pb.delete();
                    if H_hr.num
                        OUT_H_hr.HnumL = HnumList;
                    end
                    if H_hr.coe
                        OUT_H_hr.HcoeL = HcoeList;
                    end
                    OUT_H_hr.vectorL = OutVectorList;
                    
            end
            H_hr = OUT_H_hr;
            H_hr.Basis_num = OUT_WAN_NUM;
            if options.force_list
                if strcmp(H_hr.Type,'mat')
                    H_hr = H_hr.rewrite();
                end
            else
                if ~strcmp(H_hr.Type,'mat')
                    H_hr = H_hr.rewind();
                end
                
            end
        end
        function H_hr = unfold_hr(H_hr,Ns,options)
            % Returns tight-binding model representing a primitive cell of a current object.
            %
            % * Label: play_hr
            %
            %--------  init  --------
            %import vasplib_tool.*
            arguments
                H_hr HR;
                Ns double = eye(3);
                options.Accuracy double = 1e-6;
                options.force_list = false;
                options.orb_idL = [];
            end
            %--------  init  --------
            V= 1/abs(round(det(Ns)));% intger forcely
            OUT_WAN_NUM = H_hr.WAN_NUM*V;
            if  rem(OUT_WAN_NUM,1)
                error('The supercell matrix for primitive cell is not right.');
            else
                
            end
            WANNUM= H_hr.WAN_NUM;
            Accuracy = options.Accuracy;
            sc_orbL = H_hr.orbL;
            %--------  check  --------
            % checks on super-lattice Ns
            [pc_orb,pc_orbL_full,pc_elementL,pc_quantumL,sc_orb_idL,~,pc_orb_selectL] = H_hr.unfold_orb(Ns,Accuracy,options.orb_idL);
            %H_hr2 = H_hr.reseq(find(pc_orb_selectL));
            % create super-cell tb_model object to be returned
            %OUT_H_hr_ = struct('seq',[],'vector',[],'Degen',[],'key',[],'nokey',[],'Hstr',[],'Hsym',[],'Hcoe',[],'Hnum',[]);
            % OUT_H_xyz = repmat(H_xyz_ ,[NRPTS,1]);    % Hamiltonian of every cell;
            OUT_H_hr = H_hr;
            OUT_H_hr = OUT_H_hr.clean(WANNUM);
            OUT_H_hr.orbL = pc_orb; % update orbL
            OUT_H_hr.quantumL = pc_quantumL; % update orbL
            OUT_H_hr.elementL = pc_elementL; % update orbL
            OUT_H_hr.Rm = Ns\OUT_H_hr.Rm;
            NRPTS_  = H_hr.NRPTS;
            if strcmp(H_hr.Type,'sparse')
                H_hr = H_hr.full();
            end
            fprintf('Search done; begin set hoppings\n');
            fprintf('We can improve the perfomance later\n');
            %             t1=clock;
            %num_pc = size(pc_vec,1);
            %num_sc = size(sc_vec,1);
            %set hopping terms
            Accuracy_roundn = round(log(Accuracy)/log(10));
            if (strcmp(H_hr.Type,'mat')) ||options.force_list
                fprintf("Attention enforce list mode in the folding process !");
                H_hr = H_hr.rewrite();
                OUT_H_hr = OUT_H_hr.rewrite();
                %options.force_list = true;
            end
            switch H_hr.Type
                case 'list'
                    VectorList = double(H_hr.vectorL);
                    %OutVectorList = VectorList;
                    pb = vasplib_tool_outer.CmdLineProgressBar(...
                        'Generate process: UNFOLDING:');
                    sc_hiL = VectorList(:,4);
                    sc_hjL = VectorList(:,5);
                    hjL_orb_id_in_primitiveL = sc_orb_idL(sc_hjL);
                    hiL_orb_id_in_primitiveL = sc_orb_idL(sc_hiL);
                    Npc_orb_selectL = find(pc_orb_selectL);
                    % pc_orb_selectL(Npc_orb_selectL,:);
                    SelectedL = ismember(sc_hiL,Npc_orb_selectL);
                    if H_hr.num
                        HnumList = H_hr.HnumL(SelectedL,:);
                        %nHopping = length(H_hr.HnumL);
                    end
                    if H_hr.coe
                        HcoeList = H_hr.HcoeL(SelectedL,:);
                        %nHopping = length(H_hr.HcoeL);
                    end
                    hjL_orb_id_in_primitiveL = hjL_orb_id_in_primitiveL(SelectedL);
                    hiL_orb_id_in_primitiveL = hiL_orb_id_in_primitiveL(SelectedL);
                    Selected_vectorL = VectorList(SelectedL,:);
                    ind_R_in_supercellL = double(Selected_vectorL(:,1:3));
                    Selected_sc_hiL = Selected_vectorL(:,4);
                    Selected_sc_hjL = Selected_vectorL(:,5);
                    TijL_in_supercellL = ind_R_in_supercellL + ...
                        sc_orbL(Selected_sc_hjL,:) - sc_orbL(Selected_sc_hiL,:);
                    TijL_in_primitiveL = TijL_in_supercellL*Ns;
                    hiL_orbL_in_primitiveL = pc_orbL_full(Selected_sc_hiL,:);
                    hjL_plusR_orbL_in_primitiveL = hiL_orbL_in_primitiveL + TijL_in_primitiveL;
                    RvectorL_in_primitiveL = floor(hjL_plusR_orbL_in_primitiveL);
                    %hjL_orbL_in_primitiveL = hjL_plusR_orbL_in_primitiveL - RvectorL_in_primitiveL;
                    OutVectorList = [RvectorL_in_primitiveL,hiL_orb_id_in_primitiveL.',hjL_orb_id_in_primitiveL.'];
                    pb.delete();
                    if H_hr.num
                        OUT_H_hr.HnumL = HnumList;
                    end
                    if H_hr.coe
                        OUT_H_hr.HcoeL = HcoeList;
                    end
                    OUT_H_hr.vectorL = OutVectorList;
                otherwise
            end
            H_hr = OUT_H_hr;
            H_hr.Basis_num = OUT_WAN_NUM;
            if options.force_list
                if strcmp(H_hr.Type,'mat')
                    H_hr = H_hr.rewrite();
                end
            else
                if ~strcmp(H_hr.Type,'mat')
                    H_hr = H_hr.rewind();
                end
                
            end
        end
        function H_hr = translation(H_hr,translation_vector,options)
            arguments
                H_hr HR;
                translation_vector double = [0,0,0];
                options.Accuracy double = 1e-6;
                options.force_list = false;
            end
            %Accuracy_roundn = log10(options.Accuracy );
            %WANNUM = H_hr.WAN_NUM;
            OUT_H_hr = H_hr;
            pc_orb = H_hr.orbL;
            OUT_H_hr.orbL = mod(H_hr.orbL + translation_vector,1);
            sc_orb = OUT_H_hr.orbL;
            if strcmp(H_hr.Type,'sqarse')
                H_hr = H_hr.full();
            end
            if strcmp(H_hr.Type,'mat')
                fprintf("Translation must use list mode!");
                H_hr = H_hr.rewrite();
                OUT_H_hr = OUT_H_hr.rewrite();
                GiveBackMat =true;
                %options.force_list = true;
            else
                GiveBackMat =false;
            end
            if H_hr.num
                HnumList = H_hr.HnumL;
                %nHopping = length(H_hr.HnumL);
            end
            if H_hr.coe
                HcoeList = H_hr.HcoeL;
                %nHopping = length(H_hr.HcoeL);
            end
            VectorList = double(H_hr.vectorL);
            %OutVectorList = VectorList;
            %             pb = vasplib_tool_outer.CmdLineProgressBar(...
            %                 'Generate process: TRANSLATION:');
            %;
            hiL = VectorList(:,4);
            hjL = VectorList(:,5);
            ind_RL = double(VectorList(:,1:3));
            % ----- find hj --------
            %ind_R_in_supercellL = double(ind_RL+cur_sc_vec)/Ns;
            %ind_R_in_supercellL = roundn(ind_R_in_supercellL,Accuracy_roundn);
            %sc_partL=floor(ind_R_in_supercellL); % round down!
            %orig_partL=ind_RL+cur_sc_vec-double(sc_partL);
            %[~,pair_indL] = ismember(orig_partL,sc_vec,'rows');
            % ----- ******* --------
            % ----- find Rvector --------
            % cur_sc_vec = [0,0,0]
            sc_hjL = hjL;
            sc_hiL = hiL;
            indRtiL = pc_orb(hiL,:);
            indRti_in_supercellL = sc_orb(sc_hiL,:);
            real_sc_vecL = indRti_in_supercellL - indRtiL;
            indRtjL = real_sc_vecL+ind_RL+pc_orb(hjL,:);
            indRtj_in_supercellL = (double(indRtjL));
            indR_in_supercellL  = floor(indRtj_in_supercellL);
            OutVectorList= ...
                [indR_in_supercellL,sc_hiL,sc_hjL];
            %             pb.print(icur_sc_vec,num_sc,' ...');
            %             pb.delete();
            if H_hr.num
                OUT_H_hr.HnumL = HnumList;
            end
            if H_hr.coe
                OUT_H_hr.HcoeL = HcoeList;
            end
            OUT_H_hr.vectorL = OutVectorList;
            H_hr = OUT_H_hr;
            H_hr.Basis_num = OUT_H_hr.WAN_NUM;
            if options.force_list
                if strcmp(H_hr.Type,'mat')
                    H_hr = H_hr.rewrite();
                end
            elseif GiveBackMat
                if ~strcmp(H_hr.Type,'mat')
                    H_hr = H_hr.rewind();
                end
            else
            end
        end
        function [sc_orb,sc_vec,sc_elementL,sc_quantumL] = supercell_orb(H_hr,Ns,Accuracy)
            %--------  init  --------
            %import vasplib_tool.*
            arguments
                H_hr HR;
                Ns double = eye(3);
                Accuracy double = 1e-6;
            end
            %--------  init  --------
            %             use_sc_red_lat=Ns;
            %--------  check  --------
            % checks on super-lattice Ns
            if ~(Ns == round(Ns))
                error("sc_red_lat array elements must be integers");
            end
            if det(Ns) < Accuracy
                error("Super-cell lattice vectors length/area/volume too close to zero, or zero.");
            end
            if det(Ns)<0.0
                error("Super-cell lattice vectors need to form right handed system.");
            end
            nAccuracy = round(log10(Accuracy));
            %--------  init  --------
            % Ns = round(Ns);
            orb_init = H_hr.orbL;
            %OUT_WAN_NUM = H_hr.WAN_NUM*V ;
            % conservative estimate on range of search for super-cell vectors
            max_R=max(abs(Ns))*3;
            % candidates for super-cell vectors
            % this is hard-coded and can be improved!
            sc_cands = zeros( (2*max_R(1)+1) * (2*max_R(2)+1) * (2*max_R(3)+1) ,3 );
            count_tmp = 1;
            for i = -max_R(1):max_R(1)
                for j = -max_R(2):max_R(2)
                    for k = -max_R(3):max_R(3)
                        sc_cands(count_tmp,:)=[i,j,k];
                        count_tmp = count_tmp +1;
                    end
                end
            end
            % find all vectors inside super-cell
            % store them here
            sc_vec=([]);
            eps_shift=sqrt(2.0)*1.0E-8;% shift of the grid, so to avoid double counting
            %
            for ivec = 1:length(sc_cands)
                % compute reduced coordinates of this candidate vector in the super-cell frame
                tmp_red=HR.to_red_sc(sc_cands(ivec,:),Ns);
                % check if in the interior
                inside = 1;
                for it = 1:length(tmp_red)
                    t = tmp_red(it);
                    if t <= -1.0*eps_shift || t>1.0-eps_shift
                        inside=0;
                    end
                end
                if inside == 1
                    %                 disp(tmp_red);
                    %         disp(sc_cands(ivec,:));
                    sc_vec=[sc_vec;sc_cands(ivec,:)];
                end
            end
            % number of times unit cell is repeated in the super-cell
            [num_sc,~] = size(sc_vec);
            %disp(num_sc);
            % check that found enough super-cell vectors
            if round(round(abs(det(Ns)))) ~= num_sc
                error("\n\nSuper-cell generation failed! Wrong number of super-cell vectors found.");
            end
            % -------- bug-fix ----------
            % sym instead
            % Ns = sym(Ns);
            % sc_vec = sym(sc_vec);
            WANNUM = H_hr.WAN_NUM;
            % cartesian vectors of the super lattice
            % sc_cart_lat=Ns*Rm; %must Ns Rm
            % orbitals of the super-cell tight-binding model
            sc_orb=zeros(num_sc*WANNUM,3);
            sc_elementL = repmat(H_hr.elementL,[num_sc,1]);
            sc_quantumL = repmat(H_hr.quantumL,[num_sc,1]);
            count_tmp = 1;
            for icur_sc_vec = 1:num_sc % go over all super-cell vectors
                cur_sc_vec = sc_vec(icur_sc_vec,:);
                for iorb = 1:WANNUM % go over all orbitals
                    %                     orb = orbital_init(iorb,:);
                    % shift orbital and compute coordinates in
                    % reduced coordinates of super-cell
                    sc_orb(count_tmp,:) = roundn((orb_init(iorb,:)+cur_sc_vec)/Ns,nAccuracy);
                    if sc_orb(count_tmp,1) ==1
                        sc_orb(count_tmp,1)=0;
                    end
                    if sc_orb(count_tmp,2) ==1
                        sc_orb(count_tmp,2)=0;
                    end
                    if sc_orb(count_tmp,3) ==1
                        sc_orb(count_tmp,3)=0;
                    end
                    count_tmp =  count_tmp +1;
                end
            end
            sc_vec = int32(sc_vec);
            sc_orb = mod(sc_orb,1);% attention this may wrongly set!
        end
        function [pc_orb,pc_orbL_full,pc_elementL,pc_quantumL,orb_id_L,pc_orb_id_L,pc_orb_selectL] = unfold_orb(H_hr,Ns,Accuracy,orb_id_L)
            %--------  init  --------
            %import vasplib_tool.*
            arguments
                H_hr HR;
                Ns double = eye(3);
                Accuracy double = 1e-6;
                orb_id_L = [];
            end
            %--------  init  --------
            %             use_sc_red_lat=Ns;
            %--------  check  --------
            % checks on super-lattice Ns
            if ~(Ns == round(Ns))
                error("sc_red_lat array elements must be integers");
            end
            if det(Ns) < Accuracy
                error("Super-cell lattice vectors length/area/volume too close to zero, or zero.");
            end
            if det(Ns)<0.0
                error("Super-cell lattice vectors need to form right handed system.");
            end
            nAccuracy = round(log10(Accuracy));
            %--------  init  --------
            % Ns = round(Ns);
            orb_init = H_hr.orbL;
            %OUT_WAN_NUM = H_hr.WAN_NUM*V ;
            % -------- bug-fix ----------
            % sym instead
            % Ns = sym(Ns);
            % sc_vec = sym(sc_vec);
            WANNUM = H_hr.WAN_NUM;
            % orbitals of the super-cell tight-binding model
            pc_orb=zeros(WANNUM,3);
            translation_vector=zeros(WANNUM,3);
            % Suppose we know nothing about the supercell sequence

            count_tmp = 1;
            cur_sc_vec = [0 0 0];
            for iorb = 1:WANNUM % go over all orbitals
                % shift orbital and compute coordinates
                % reduced coordinates of primitive-cell
                pc_orb_tmp = roundn((orb_init(iorb,:)*Ns-cur_sc_vec),nAccuracy);
                [pc_orb_incell,translation_vector(count_tmp,:)] = vasplib.translation_orb(pc_orb_tmp);
                pc_orb(count_tmp,:) = pc_orb_incell;
                count_tmp =  count_tmp +1;
            end
            pc_orbL_full = pc_orb;
            %try
            %    PTL = [pc_orb,H_hr.elementL,H_hr.quantumL];
            %catch
            PTL = [pc_orb,translation_vector];
            %end
            [uniquePTL,uniqueAll,TrackingLater] = unique(PTL,'rows','stable');
            %
            translation_vector_mini = translation_vector(uniqueAll,:);
            pc_orb_mini = pc_orb(uniqueAll,:);
            % [~,uniqueTV] = unique(translation_vector,'rows','stable');
            [~,uniqueOrb] = uniquetol(pc_orb_mini,10^(nAccuracy+1),'ByRows',true);
            % Unique_translation_vector = translation_vector_mini(uniqueOrb,:);
            [pc_orb_selectL,~] = ismember(PTL,uniquePTL(uniqueOrb,:),'rows');
            pc_orb = pc_orbL_full(pc_orb_selectL,:);
            % pc_vec = int32(translation_vector(pc_orb_selectL,:));
            % pc_orb = mod(pc_orb,1);% attention this may wrongly set!
            pc_elementL = H_hr.elementL(pc_orb_selectL,:);
            pc_quantumL = H_hr.quantumL(pc_orb_selectL,:);
            pc_orb_id_L = 1:size(pc_orb,1);
            if isempty(orb_id_L)
                pc_peqL = [pc_orb,pc_elementL,pc_quantumL];
                sc_peqL = [pc_orbL_full,H_hr.elementL,H_hr.quantumL];
            else
                pc_peqL = [pc_orb,orb_id_L(pc_orb_selectL).'];
                sc_peqL = [pc_orbL_full,orb_id_L.'];
            end
            % change nan
            pc_peqL(isnan(pc_peqL)) = 0;
            sc_peqL(isnan(sc_peqL)) = 0;
            %
            if size(unique(pc_peqL,"rows"),1) == size(pc_peqL,1)
            else
                warning('check duplicate orbital in primitive cell, the unfolding process is unreliable')
            end
            [~,sc_orb_selectL] = ismembertol(sc_peqL,pc_peqL,10^(nAccuracy+1),'ByRows',true);
            orb_id_L = pc_orb_id_L(sc_orb_selectL);
        end
        function H_hr = Hnanowire_gen(H_hr,Nslab,np,vacuum_mode,options)
            %---------------------------
            arguments
                H_hr HR;
                Nslab double= [1 1 1];
                np double{mustBeInteger} =1
                vacuum_mode double = 0;
                options.symbolic logical = false;
                options.fast logical = true;
            end
            if options.symbolic
                % not recommend
                H_hr.coe = true;
            else
                H_hr.coe = false;
                H_hr.num = true;
            end
            %--------  reshape  --------
            sort_dir = [3,2,1];
            vectorList = double(H_hr.vectorL);
            [vector_list_init,sort_label] = sortrows(vectorList,sort_dir) ;% sort fin_dir
            H_hr = H_hr.reseq(':',sort_label);
            WANNUM =  H_hr.WAN_NUM;
            % --------  init  --------
            WAN_NUM_x = WANNUM*Nslab(1);
            WAN_NUM_y = WAN_NUM_x*Nslab(2);
            % --------  3D begining  --------
            [unique_z,unique_label_z]= unique(vector_list_init(:,3),'rows');
            cutlist_z = HR.unique_label2cutlist(unique_label_z,H_hr.NRPTS);
            NRPTS_z = length(unique_z);
            vector_list_wire = int32(zeros(NRPTS_z,3));
            vector_list_wire(:,3) = int32(unique_z);
            % init
            vertor_list_xy{NRPTS_z} = vector_list_init(cutlist_z(NRPTS_z,1):cutlist_z(NRPTS_z,2),:);
            Hnum_list_wire{NRPTS_z} = sparse(WAN_NUM_y,WAN_NUM_y);
            if strcmp(H_hr.Type,'sparse')
                Hnum_list_xy{NRPTS_z} = H_hr.HnumL(cutlist_z(NRPTS_z,1):cutlist_z(NRPTS_z,2));
            else
                Hnum_list_xy{NRPTS_z} = H_hr.HnumL(:,:,cutlist_z(NRPTS_z,1):cutlist_z(NRPTS_z,2));
            end
            % --------  vector  --------
            if  strcmp(H_hr.Type,'sparse')
                for iz = 1:NRPTS_z-1
                    %                     vector_list_wire(iz,3) = unique_z(iz);
                    vertor_list_xy{iz} = vector_list_init(cutlist_z(iz,1):cutlist_z(iz,2),:);
                    Hnum_list_xy{iz} = H_hr.HnumL(cutlist_z(iz,1):cutlist_z(iz,2));
                    Hnum_list_wire{iz} = sparse(WAN_NUM_y,WAN_NUM_y);
                end
            else
                for iz = 1:NRPTS_z-1
                    %                     vector_list_wire(iz,3) = unique_z(iz);
                    vertor_list_xy{iz} = vector_list_init(cutlist_z(iz,1):cutlist_z(iz,2),:);
                    Hnum_list_xy{iz}   = H_hr.HnumL(:,:,cutlist_z(iz,1):cutlist_z(iz,2));
                    Hnum_list_wire{iz} = sparse(WAN_NUM_y,WAN_NUM_y);
                end
            end
            
            TYPE = H_hr.Type;
            if np >1
                disp('parallel mode, we will use local settings, please set before.');
                np_handle = parpool('local',np);
                parfor iz = 1:NRPTS_z
                    fprintf('Gen (%d/%d) NRPT z \n',iz,NRPTS_z);
                    Hnum_list_xy_iz = Hnum_list_xy{iz};
                    vertor_list_xy_iz = vertor_list_xy{iz};
                    Hnum_list_wire{iz} = HR.Hnum_list_wire_iz_gen(Hnum_list_xy_iz,vertor_list_xy_iz,iz,WANNUM,WAN_NUM_x,WAN_NUM_y,Nslab,TYPE);
                end
                delete(np_handle);
            else
                for iz = 1:NRPTS_z
                    fprintf('Gen (%d/%d) NRPT z \n',iz,NRPTS_z);
                    Hnum_list_xy_iz = Hnum_list_xy{iz};
                    vertor_list_xy_iz = vertor_list_xy{iz};
                    Hnum_list_wire{iz} = HR.Hnum_list_wire_iz_gen(Hnum_list_xy_iz,vertor_list_xy_iz,iz,WANNUM,WAN_NUM_x,WAN_NUM_y,Nslab, TYPE);
                end
            end
            if strcmp(H_hr.Type,'sparse')
                H_hr.HnumL = Hnum_list_wire;
            else
                HnumL_temp = zeros(WAN_NUM_y,WAN_NUM_y,NRPTS_z);
                for iz = 1:NRPTS_z
                    HnumL_temp(:,:,iz) = full(Hnum_list_wire{iz});
                end
                H_hr.HnumL = HnumL_temp ;
            end
            H_hr.vectorL = int32(vector_list_wire);
            H_hr.orbL = H_hr.nanowire_orb(Nslab,vacuum_mode,'fast',options.fast);
            for i = 1:3
                if  Nslab(i)<1
                    Nslab(i) = 1;
                end
            end
            num_sc = Nslab(1)* Nslab(2)* Nslab(3);
            H_hr.elementL = repmat(H_hr.elementL,[num_sc,1]);
            H_hr.quantumL = repmat(H_hr.quantumL,[num_sc,1]);
        end
        function H_hr = descritize(H_hr,Nslab,options)
            arguments
                H_hr HR;
                Nslab double = [0,10,0];
                options.rmfunc function_handle=@()(1);
                options.Rotation  = sym(eye(3));
                options.np = 1;
                options.vacuum_mode = 0;
            end
            if strcmp(functions(options.rmfunc).function , '@()(1)')
                rm_mode =false;
            else
                rm_mode =true;
            end
            if isequal(options.Rotation ,sym(eye(3)))
                rotate_mode = false;
            else
                rotate_mode = true;
            end
        end
        function orbital_out = nanowire_orb(H_hr,fin_dir,vacuum_mode,options)
            arguments
                H_hr HR;
                fin_dir = [1 1 0];
                vacuum_mode = false;
                options.fast = true;
            end
            if isempty(H_hr.orbL)
                orbital_init = zeros(H_hr.WAN_NUM,3);
            else
                orbital_init = H_hr.orbL;
            end
            
            switch vacuum_mode
                case 0
                    orbital_out = orbital_init;
                    % orb
                    for i = 1:3
                        Nslab = fin_dir(i);
                        if  Nslab == 0
                            Nslab = 1;
                        end
                        count = 0;
                        WAN_NUM_tmp = size(orbital_out,1);
                        fin_orb = zeros(WAN_NUM_tmp*Nslab,3);
                        for inum = 1:Nslab     % go over all cells in finite direction
                            for j = 1:WAN_NUM_tmp  % go over all orbitals in one cell
                                count =count +1;
                                % make a copy of j-th orbital
                                orb_tmp=orbital_out(j,:);
                                % change coordinate along finite direction ; fractional
                                orb_tmp(i)= (orb_tmp(i) + double(inum-1)) / Nslab;
                                % add to the list
                                fin_orb(count,:)=orb_tmp;
                                % do the onsite energies at the same time
                            end
                        end
                        orbital_out = fin_orb;
                    end
                    % POSCAR
                    Ns = [1 0 0;0 1 0;0 0 1];
                    Ns = Ns.*fin_dir;
                    fin_dir_list = [0 0 0];
                    %disp(fin_dir_list);
                    if ~options.fast
                        [Rm_tmp,sites_tmp,Atom_name_tmp,Atom_num_tmp]=HR.POSCAR_readin();
                        H_hr.supercell(Ns,'POSCAR_super_fin',Rm_tmp,sites_tmp,Atom_name_tmp,Atom_num_tmp,fin_dir_list);
                    end
                case 1
                    orbital_out = orbital_init;
                    % orb
                    for i = 1:3
                        Nslab = fin_dir(i);
                        if  Nslab == 0
                            Nslab = 1;
                        end
                        count = 0;
                        WAN_NUM_tmp = size(orbital_out,1);
                        fin_orb = zeros(WAN_NUM_tmp*Nslab,3);
                        for inum = 1:Nslab % go over all cells in finite direction
                            for j = 1:WAN_NUM_tmp  % go over all orbitals in one cell
                                count =count +1;
                                % make a copy of j-th orbital
                                orb_tmp=orbital_out(j,:);
                                % change coordinate along finite direction ; fractional
                                orb_tmp(i)= (orb_tmp(i) + double(inum-1))/Nslab;
                                % add to the list
                                fin_orb(count,:)=orb_tmp;
                                % do the onsite energies at the same time
                            end
                        end
                        orbital_out = fin_orb;
                    end
                    
                    Ns = [1 0 0;0 1 0;0 0 1];
                    Ns = Ns.*fin_dir;
                    fin_dir_list = double(fin_dir>1);
                    if ~options.fast
                        %disp(fin_dir_list);
                        [Rm_tmp,sites_tmp,Atom_name_tmp,Atom_num_tmp]=HR.POSCAR_readin();
                        % gen POSCAR
                        H_hr.supercell(Ns,'POSCAR_super_fin',Rm_tmp,sites_tmp,Atom_name_tmp,Atom_num_tmp,fin_dir_list);
                    else
                        [Rm_tmp,sites_tmp,Atom_name_tmp,Atom_num_tmp]=HR.POSCAR_readin();
                    end
                    % rebuild fin_orb
                    Rm_tmp = Ns*Rm_tmp;
                    Rmlength1 = norm (Rm_tmp(1,:));
                    Rmlength2 = norm (Rm_tmp(2,:));
                    Rmlength3 = norm (Rm_tmp(3,:));
                    Rm_s_fin_add = [10*Rm_tmp(1,:)*fin_dir_list(1)/Rmlength1;...
                        10*Rm_tmp(2,:)*fin_dir_list(2)/Rmlength2;...
                        10*Rm_tmp(3,:)*fin_dir_list(3)/Rmlength3];
                    Rm_s_fin = Rm_tmp + Rm_s_fin_add ;
                    Rc_s_fin_add = [1/2, 1/2 ,1/2] ;
                    Rr_s_fin_add = Rc_s_fin_add * Rm_s_fin_add;
                    [nfinorb,~ ]= size(fin_orb);
                    for  i = 1:nfinorb
                        Rr_orb = fin_orb(i,:)*Rm_tmp;
                        Rr_s_fin = Rr_orb + Rr_s_fin_add;
                        Rc_s_fin = Rr_s_fin / Rm_s_fin;
                        fin_orb(i,:) = Rc_s_fin ;
                    end
                    orbital_out = fin_orb;
            end
            
        end
        function H_hr = supercell(H_hr,Ns,filename,Rm,sites,Atom_name,Atom_num,findir)
            % usage: [sites]=supercell(Ns)
            % Ns =
            %
            % [ N_11, N_12, N_13]
            % [ N_21, N_22, N_23]
            % [ N_31, N_32, N_33]
            % sites
            % note : site=struct('seq',[],'inseq',[],'rc1',[],'rc2',[],'rc3',[],'name',[],'nameseq',[]);
            %   Rm
            %  a1 a1x a1y a1z
            %  a2 a2x a2y a2z
            %  a3 a3x a3y a3z
            %   unit A?e-6
            Accuracy=8;
            
            % nargin
            if nargin < 4
                Rm = H_hr.Rm;
            end
            if nargin < 5
                sites = H_hr.sites;
            end
            if nargin < 6
                Atom_name = H_hr.Atom_name;
            end
            if nargin < 7
                Atom_num = H_hr.Atom_num;
            end
            if nargin < 8
                findir = [0,0,0];
            end
            if nargin == 4
                Origin_POSCAR = Rm;
                [Rm,sites,Atom_name,Atom_num,~]=HR.POSCAR_readin(Origin_POSCAR,'vasp');
                findir = [0,0,0];
            end
            if nargin < 3
                filename  = 'POSCAR_super';
            end
            if nargin < 2
                Ns = eye(3);
            end
            
            format long;
            %Ns_=inv(Ns);
            V_Ns=dot(Ns(1,:),cross(Ns(2,:),Ns(3,:)));
            if V_Ns < 0
                Ns(3,:) = -Ns(3,:);
            end
            V_Ns = abs(V_Ns);
            if V_Ns < 1
                error('check Ns');
            end
            
            % refresh data first
            [~,nsites] = size(sites);
            for i =1:nsites
                sites(i).rc1 = HR.plusrc(sites(i).rc1);
                if sites(i).rc1  < 0
                    disp('bugtest')
                end
                sites(i).rc2 = HR.plusrc(sites(i).rc2);
                if sites(i).rc2  < 0
                    disp('bugtest')
                end
                sites(i).rc3 = HR.plusrc(sites(i).rc3);
                if sites(i).rc3  < 0
                    disp('bugtest')
                end
            end
            
            if V_Ns >= 1
                % conservative estimate on range of search for super-cell vectors
                max_R=max(abs(Ns))*3;
                % candidates for super-cell vectors
                % this is hard-coded and can be improved!
                sc_cands = [];
                for i = -max_R(1):max_R(1)
                    for j = -max_R(2):max_R(2)
                        for k = -max_R(3):max_R(3)
                            %for k = 1:1
                            sc_cands=[sc_cands;[i,j,k]];
                        end
                    end
                end
                
                % find all vectors inside super-cell
                % store them here
                sc_vec=[];
                eps_shift=sqrt(2.0)*1.0E-8;% shift of the grid, so to avoid double counting
                
                %
                for ivec = 1:length(sc_cands)
                    % compute reduced coordinates of this candidate vector in the super-cell frame
                    tmp_red=HR.to_red_sc(sc_cands(ivec,:),Ns);
                    % check if in the interior
                    inside = 1;
                    for it = 1:length(tmp_red)
                        t = tmp_red(it);
                        if t <= -1.0*eps_shift || t>1.0-eps_shift
                            inside=0;
                        end
                    end
                    if inside == 1
                        %disp(tmp_red);
                        %disp(sc_cands(ivec,:));
                        sc_vec=[sc_vec;sc_cands(ivec,:)];
                        % number of times unit cell is repeated in the super-cell
                    end
                end
                
                [num_sc,~] = size(sc_vec);
                %disp(num_sc);
                % check that found enough super-cell vectors
                if round(round(abs(det(Ns)))) ~= num_sc
                    error("\n\nSuper-cell generation failed! Wrong number of super-cell vectors found.");
                end
                
                % cartesian vectors of the super lattice
                % sc_cart_lat=Ns*Rm; %must Ns Rm
                countnum=0;
                countnum2=0;
                for elementseq=1:length(Atom_num)
                    % for each Atom pattern
                    for n=1:Atom_num(elementseq)
                        % for each atom in Atom pattern
                        countnum2=countnum2+1;
                        Rpc1=sites(countnum2).rc1;
                        Rpc2=sites(countnum2).rc2;
                        Rpc3=sites(countnum2).rc3;
                        Rc=[Rpc1 Rpc2 Rpc3];
                        
                        for icur_sc_vec = 1:num_sc % go over all super-cell vectors
                            cur_sc_vec = sc_vec(icur_sc_vec,:);
                            Rsc = HR.to_red_sc(Rc+cur_sc_vec,Ns);
                            Rsc = round(Rsc.*10^Accuracy)./10^Accuracy;
                            rsc1=Rsc(1);
                            rsc2=Rsc(2);
                            rsc3=Rsc(3);
                            if rsc1 ==1
                                rsc1 = 0;
                            end
                            if rsc2 ==1
                                rsc2 = 0;
                            end
                            if rsc3 ==1
                                rsc3 = 0;
                            end
                            countnum=countnum+1;
                            sites_s(countnum).rc1=rsc1;%along a1
                            sites_s(countnum).rc2=rsc2;%along a2
                            sites_s(countnum).rc3=rsc3;%along a3
                            sites_s(countnum).name=Atom_name(elementseq);
                        end
                        
                    end
                end
                sites=sites_s;
            end
            % Atom_name Keep
            % Atom_number [Rm,sites]=supercell(Ns,Rm,sites,Atom_name,Atom_num,fin_dir,filename);
            Atom_num_s=round(Atom_num*V_Ns);
            Rm_s=Ns*Rm;
            [~,nsites] = size(sites);
            
            if findir(1) ~=0 || findir(2) ~=0 || findir(3) ~=0
                Rm = Ns*Rm;
                disp('findir_mode');
                %init
                
                Rmlength1 = abs(norm (Rm(1,:)));
                Rmlength2 = abs(norm (Rm(2,:)));
                Rmlength3 = abs(norm (Rm(3,:)));
                Rm_s_fin_add = [10*Rm(1,:)*findir(1)/Rmlength1;...
                    10*Rm(2,:)*findir(2)/Rmlength2;...
                    10*Rm(3,:)*findir(3)/Rmlength3];
                Rm_s_fin = Rm + Rm_s_fin_add ;
                Rc_s_fin_add = [1/2, 1/2 ,1/2] ;
                Rr_s_fin_add = Rc_s_fin_add * Rm_s_fin_add;
                for  i = 1:nsites
                    Rr = [sites(i).rc1,sites(i).rc2,sites(i).rc3]*Rm;
                    Rr_s_fin = Rr + Rr_s_fin_add;
                    Rc_s_fin = Rr_s_fin /Rm_s_fin;
                    sites_s(i).rc1 = Rc_s_fin(1);
                    sites_s(i).rc2 = Rc_s_fin(2);
                    sites_s(i).rc3 = Rc_s_fin(3);
                end
                Rm_s=Rm_s_fin;
            end
            
            % ------------------------------------------------------------------------------
            H_hr.Atom_num = Atom_num_s;
            H_hr.Atom_name = Atom_name;
            H_hr.sites=sites_s;
            H_hr.Rm=Rm_s;
            H_hr.POSCAR_gen(filename);
        end
        function [H_hr_out,H_hr_pi_plus,H_hr_pi_minus] = realmap(H_hr)
            switch H_hr.Type
                case 'sparse'
                    H_hr = H_hr.full;
                case 'list'
                    H_hr = H_hr.rewind;
                case 'mat'
            end
            H_hr_out = H_hr;
            H_hr_out.quantumL = [H_hr_out.quantumL;H_hr_out.quantumL];
            H_hr_out.elementL = [H_hr_out.elementL;H_hr_out.elementL];
            H_hr_out.orbL = [H_hr_out.orbL;H_hr_out.orbL];
            new_WAN_NUM = 2*H_hr_out.WAN_NUM;
            WANNUM = H_hr.WAN_NUM;
            NRPTS_ = H_hr.NRPTS;
            PI_plus = 1/2 * [eye(WANNUM), -1i*eye(WANNUM);1i*eye(WANNUM),eye(WANNUM)];
            PI_minus = 1/2 * [eye(WANNUM), 1i*eye(WANNUM);-1i*eye(WANNUM),eye(WANNUM)]; 
            if H_hr.coe
                HcoeList = H_hr.HcoeL;
                realHcoeList = real(HcoeList);
                imagHcoeList = imag(HcoeList);
                tempHcoeL = sym(zeros(new_WAN_NUM,new_WAN_NUM,H_hr.NRPTS));
                tempHcoeL_pi_plus = tempHcoeL;
                tempHcoeL_pi_minus = tempHcoeL;
                for i = 1:NRPTS_
                    tempHcoeL(:,:,i) = [realHcoeList(:,:,i),imagHcoeList(:,:,i);-imagHcoeList(:,:,i),realHcoeList(:,:,i)] ;
                    tempHcoeL_pi_plus(:,:,i)  = PI_plus*tempHcoeL(:,:,i)*PI_plus;
                    tempHcoeL_pi_minus(:,:,i)  = PI_minus*tempHcoeL(:,:,i)*PI_minus;
                end
                H_hr_out.HcoeL = tempHcoeL;
                H_hr_pi_plus = H_hr_out;
                H_hr_pi_plus.HcoeL = tempHcoeL_pi_plus;
                H_hr_pi_minus = H_hr_out;
                H_hr_pi_minus.HcoeL = tempHcoeL_pi_minus;
            elseif H_hr.num
                HnumList = H_hr.HnumL;
                realHnumList = real(HnumList);
                imagHnumList = imag(HnumList);
                tempHnumL = sym(zeros(new_WAN_NUM,new_WAN_NUM,H_hr.NRPTS));
                PI_plus_page = repmat(PI_plus,[1 1 NRPTS_]);
                PI_minus_page = repmat(PI_minus,[1 1 NRPTS_]);
                for i = 1:HNRPTS_
                    tempHnumL(:,:,i) = [realHnumList(:,:,i),imagHnumList(:,:,i);-imagHnumList(:,:,i),realHnumList(:,:,i)] ;
                end
                H_hr_out.HnumL = tempHnumL;
                H_hr_pi_plus = H_hr_out;
                H_hr_pi_plus.HnumL = pagemtimes(pagemtimes(PI_plus_page,tempHnumL),PI_plus_page);
                H_hr_pi_minus = H_hr_out;
                H_hr_pi_minus.HnumL = pagemtimes(pagemtimes(PI_minus_page,tempHnumL),PI_minus_page);
            end

        end
    end
    % ----------------  modify tool --------------------
    methods (Static,Hidden,Access= protected)
        function Hnum_list_wire_iz =  Hnum_list_wire_iz_gen(Hnum_list_xy_iz,vertor_list_xy_iz,iz,WAN_NUM,WAN_NUM_x,WAN_NUM_y,Nslab,Type)
            Hnum_list_wire_iz = zeros(WAN_NUM_y,WAN_NUM_y);
            [NRPTS_xy,~] = size(vertor_list_xy_iz);
            % --------  2D begining  --------
            [unique_y,unique_label_y]= unique(vertor_list_xy_iz(:,2),'rows');
            cutlist_y = HR.unique_label2cutlist(unique_label_y,NRPTS_xy);
            NRPTS_y = length(unique_y);
            % vertor_list_y = repmat([0 0 unique_z(iz)],3,1);
            %Hnum_list_iy_iz{NRPTS_y} = sparses(WAN_NUM_x,WAN_NUM_x);
            for iy = 1:NRPTS_y
                fprintf('%d th NRPT z ---- Gen (%d/%d) NRPT y \n',iz,iy,NRPTS_y);
                %         vertor_list_y(iy,2) = unique_y(iy);
                vertor_list_x_iy_iz = vertor_list_xy_iz(cutlist_y(iy,1):cutlist_y(iy,2),:);
                if strcmp(Type,'sparse')
                    Hnum_list_x_iy_iz   = Hnum_list_xy_iz(cutlist_y(iy,1):cutlist_y(iy,2));
                else
                    Hnum_list_x_iy_iz   = Hnum_list_xy_iz(:,:,cutlist_y(iy,1):cutlist_y(iy,2));
                end
                [NRPTS_x,~] = size(vertor_list_x_iy_iz);
                % --------  1D begining  --------
                Hnum_list_iy_iz =zeros(WAN_NUM_x,WAN_NUM_x);
                if strcmp(Type,'sparse')
                    for ix = 1:NRPTS_x % go over all NRPTS
                        [ilist,jlist,amplist] = find(Hnum_list_x_iy_iz{ix});
                        %
                        nhopping = length(amplist);
                        
                        fprintf('%d th NRPT z ---- %d th NRPT y ---- Gen (%d/%d) NRPT x \n',iz,iy,ix,NRPTS_x);
                        % lattice vector of the hopping
                        ind_R_x = vertor_list_x_iy_iz(ix,:);
                        jump_fin_x=ind_R_x(1);
                        % store by how many cells is the hopping in finite direction
                        %temp_Hnum =sparse(OUT_WAN_NUM,OUT_WAN_NUM);
                        % speed up more and more !!
                        for ih = 1:nhopping
                            i = ilist(ih);
                            j = jlist(ih);
                            amp = amplist(ih);
                            for icur_sc_vec = 1:Nslab(1) % go over all cells in finite direction mini
                                hi= i + (icur_sc_vec-1)*WAN_NUM ;
                                %disp(hi);
                                hj= j + (icur_sc_vec+jump_fin_x-1)*WAN_NUM ;
                                %disp(hj);
                                % decide whether this hopping should be added or not
                                to_add=1;
                                %disp(hj);
                                % if edges are not glued then neglect all jumps that spill out
                                if hj <= 0 || hj > WAN_NUM_x
                                    to_add=0;
                                end
                                if to_add == 1
                                    %Hnum_list_iy_iz(hi,hj) = Hnum_list_iy_iz(hi,hj)+ amp;
                                    Hnum_list_iy_iz(hi,hj) = amp;
                                end
                            end
                            
                        end
                        %OUT_Hnum_list(:,ih) = temp_Hnum;
                    end
                    %Hnum_list_iy_iz =sparse(WAN_NUM_x,WAN_NUM_x);
                else
                    pb = vasplib_tool_outer.CmdLineProgressBar('Gen NRPT x: ');
                    for ix = 1:NRPTS_x % go over all NRPTS
                        pb.print(ix,NRPTS_x);
                        % lattice vector of the hopping
                        ind_R_x = vertor_list_x_iy_iz(ix,:);
                        jump_fin_x=ind_R_x(1);
                        % store by how many cells is the hopping in finite direction
                        %temp_Hnum =sparse(OUT_WAN_NUM,OUT_WAN_NUM);
                        % speed up more and more !!
                        for  i = 1:WAN_NUM
                            for j = 1:WAN_NUM
                                % amplitude of the hop is the same
                                amp = Hnum_list_x_iy_iz(i,j,ix);
                                if norm(amp) > 0
                                    for icur_sc_vec = 1:Nslab(1) % go over all cells in finite direction mini
                                        hi= i + (icur_sc_vec-1)*WAN_NUM ;
                                        %disp(hi);
                                        hj= j + (icur_sc_vec+jump_fin_x-1)*WAN_NUM ;
                                        %disp(hj);
                                        % decide whether this hopping should be added or not
                                        to_add=1;
                                        %disp(hj);
                                        % if edges are not glued then neglect all jumps that spill out
                                        if hj <= 0 || hj > WAN_NUM_x
                                            to_add=0;
                                        end
                                        if to_add == 1
                                            %Hnum_list_iy_iz(hi,hj) = Hnum_list_iy_iz(hi,hj)+ amp;
                                            Hnum_list_iy_iz(hi,hj) = amp;
                                        end
                                    end
                                end
                            end
                        end
                        %OUT_Hnum_list(:,ih) = temp_Hnum;
                    end
                    pb.delete
                end
                jump_fin_y =  unique_y(iy);
                [ilist,jlist,amplist] = find(Hnum_list_iy_iz);
                nhopping = length(amplist);
                for ih = 1:nhopping
                    i = ilist(ih);
                    j = jlist(ih);
                    amp = amplist(ih);
                    for icur_sc_vec = 1:Nslab(2) % go over all cells in finite direction mini
                        hi= i + (icur_sc_vec-1)*WAN_NUM_x ;
                        %disp(hi);
                        hj= j + (icur_sc_vec+jump_fin_y-1)*WAN_NUM_x ;
                        %disp(hj);
                        % decide whether this hopping should be added or not
                        to_add=1;
                        %disp(hj);
                        % if edges are not glued then neglect all jumps that spill out
                        if hj <= 0 || hj > WAN_NUM_y
                            to_add=0;
                        end
                        if to_add == 1
                            %Hnum_list_wire_iz(hi,hj) = Hnum_list_wire_iz(hi,hj)+ amp;
                            Hnum_list_wire_iz(hi,hj) =  amp;
                        end
                    end
                end
            end
            Hnum_list_wire_iz = sparse(Hnum_list_wire_iz);
        end
        function Poly_priciplayer_mat = Poly_priciplayer_mat_gen(principle_layer)
            Poly_priciplayer_mat(:,:,principle_layer) = eye(principle_layer);
            for i = 1:principle_layer-1
                base_mat = zeros(principle_layer);
                base_mat(1:principle_layer-i,i+1:principle_layer) = eye(principle_layer-i);
                Poly_priciplayer_mat(:,:,principle_layer+i) = base_mat;
                Poly_priciplayer_mat(:,:,principle_layer-i) = base_mat';
            end
        end
        function cutlist= unique_label2cutlist(unique_label,NRPTS)
            cutlist(:,1)= unique_label;
            cutlist(1:end-1,2)= unique_label(2:end)-1;
            cutlist(end,2) = NRPTS;
        end
    end
    %% Green_prepare
    methods
        function H_hr= hrz_gen(H_hr,kpoints_f,fin_dir,mode)
            % change a HR_class into H_hrz by a certain Kpoints and direction
            % useful for cut_piece and GreenFunction
            %
            % * Label: HR class
            %
            %% Description of the Function:
            %%
            %% Usage:
            %
            % * H_hr_new = H_hr_old.hrz_gen(kpoints_f,fin_dir)
            % * H_hr_new = H_hr_old.hrz_gen(kpoints_f,fin_dir,mode)
            %
            %% Note:
            %
            %  Take advantage of the scope of application of the function.
            %
            %% Change log
            %
            % * Document Date: 2020/12/16
            % * Creation Date: 2020/12/16
            % * Last updated : 2021/04/16
            %
            %--------  narg  --------
            if nargin < 4
                mode = '1D';
            end
            if nargin <3
                fin_dir = 3;
            end
            if nargin < 2
                kpoints_f = [0 ,0 ,0];
            end
            %--------  init  --------
            import vasplib_tool.*
            %nkpoints_f = size(kpoints_f,1);
            factor_list = exp(1i*2*pi*(double(H_hr.vectorL)*kpoints_f.'));
            [Num,Coe] = H_hr.NumOrCoe;
            switch H_hr.Type
                case {'sparse'}
                    Hnum_kf_list = reshape(full(cell2mat(H_hr.HnumL)),H_hr.WAN_NUM,H_hr.WAN_NUM,H_hr.NRPTS);
                case 'mat'
                    if Num
                        Hnum_kf_list = H_hr.HnumL;
                    end
                    if Coe
                        Hcoe_kf_list = H_hr.HcoeL;
                    end
                case 'list'
                otherwise
                    
            end
            switch H_hr.Type
                case {'sparse'}
                    for i = 1 : H_hr.NRPTS
                        Hnum_kf_list(:,:,i) = Hnum_kf_list(:,:,i)*factor_list(i);
                    end
                case 'mat'
                    if Num
                        for i = 1 : H_hr.NRPTS
                            Hnum_kf_list(:,:,i) = Hnum_kf_list(:,:,i)*factor_list(i);
                        end
                    end
                    if Coe
                        for i = 1 : H_hr.NRPTS
                            Hcoe_kf_list(:,:,i) = Hcoe_kf_list(:,:,i)*factor_list(i);
                        end
                    end
            end
            %--------  main  --------
            switch mode
                case '0D'
                    switch H_hr.Type
                        case {'sparse'}
                            H_hr.HnumL = sum(Hnum_kf_list,3);
                        case {'mat'}
                            if Num
                                H_hr.HnumL = sum(Hnum_kf_list,3);
                            end
                            if Coe
                                H_hr.HcoeL = sum(Hcoe_kf_list,3);
                            end
                    end
                case '1D'
                    [vector_list_new,sort_label] = sortrows(H_hr.vectorL,fin_dir) ;% sort fin_dir
                    switch H_hr.Type
                        case {'sparse'}
                            Hnum_kf_list_new = Hnum_kf_list(:,:,sort_label); % sort hr
                        case {'mat'}
                            if Num
                                Hnum_kf_list_new = Hnum_kf_list(:,:,sort_label); % sort hr
                            end
                            if Coe
                                Hcoe_kf_list_new = Hcoe_kf_list(:,:,sort_label); % sort hr;
                            end
                    end
                    [unique_dir,unique_label]= unique(vector_list_new(:,fin_dir));
                    %cutlist
                    cutlist(:,1)= unique_label;
                    cutlist(1:end-1,2)= unique_label(2:end)-1;
                    cutlist(end,2) = H_hr.NRPTS;
                    %--------  H_hrz  --------
                    vector_list_hrz = zeros(length(unique_dir),3);
                    NRPTS_hrz = length(unique_label);
                    switch H_hr.Type
                        case {'sparse'}
                            Hnum_list_hrz{NRPTS_hrz} = sparse(H_hr.WAN_NUM,H_hr.WAN_NUM);
                            for i = 1:NRPTS_hrz
                                vector_list_hrz(i,fin_dir) = unique_dir(i,:);
                                Hnum_list_hrz{i} = sum(Hnum_kf_list_new(:,:,cutlist(i,1):cutlist(i,2)),3);
                            end
                        case {'mat'}
                            if Num
                                Hnum_list_hrz = zeros(H_hr.WAN_NUM,H_hr.WAN_NUM,NRPTS_hrz);
                                for i = 1:NRPTS_hrz
                                    vector_list_hrz(i,fin_dir) = unique_dir(i,:);
                                    Hnum_list_hrz(:,:,i) = sum(Hnum_kf_list_new(:,:,cutlist(i,1):cutlist(i,2)),3);
                                    Hcoe_list_hrz = sym(zeros(H_hr.WAN_NUM,H_hr.WAN_NUM,NRPTS_hrz));
                                end
                            end
                            if Coe
                                Hcoe_list_hrz = sym(zeros(H_hr.WAN_NUM,H_hr.WAN_NUM,NRPTS_hrz));
                                for i = 1:NRPTS_hrz
                                    vector_list_hrz(i,fin_dir) = unique_dir(i,:);
                                    Hcoe_list_hrz(:,:,i) = sum(Hcoe_kf_list_new(:,:,cutlist(i,1):cutlist(i,2)),3);
                                end
                            end
                    end
                    switch H_hr.Type
                        case {'sparse'}
                            H_hr.HnumL = Hnum_list_hrz;
                        case {'mat'}
                            if Num
                                H_hr.HnumL = Hnum_list_hrz;
                            end
                            if Coe
                                H_hr.HcoeL = Hcoe_list_hrz;
                            end
                    end
                    H_hr.vectorL = int32(vector_list_hrz);
                case '2D' % need check!!!
                    [vector_list_new,sort_label] = sortrows(H_hr.vectorL,fin_dir) ;% sort fin_dir
                    [unique_dir,unique_label]= unique(vector_list_new(:,fin_dir),'rows');
                    switch H_hr.Type
                        case {'sparse'}
                            Hnum_kf_list_new = Hnum_kf_list(:,:,sort_label); % sort hr
                        case {'mat'}
                            if Num
                                Hnum_kf_list_new = Hnum_kf_list(:,:,sort_label); % sort hr
                            end
                            if Coe
                                Hcoe_kf_list_new = Hcoe_kf_list(:,:,sort_label); % sort hr;
                            end
                    end
                    %cutlist
                    cutlist(:,1)= unique_label;
                    cutlist(1:end-1,2)= unique_label(2:end)-1;
                    cutlist(end,2) = H_hr.NRPTS;
                    %--------  H_hrz  --------
                    vector_list_hrz = zeros(length(unique_dir),3);
                    NRPTS_hrz = length(unique_label);
                    switch H_hr.Type
                        case {'sparse','num'}
                            Hnum_list_hrz = zeros(H_hr.WAN_NUM,H_hr.WAN_NUM,NRPTS_hrz);
                            for i = 1:NRPTS_hrz
                                vector_list_hrz(i,fin_dir) = unique_dir(i,:);
                                Hnum_list_hrz(:,:,i) = sum(Hnum_kf_list_new(:,:,cutlist(i,1):cutlist(i,2)),3);
                            end
                        case {'mat'}
                            Hnum_list_hrz = zeros(H_hr.WAN_NUM,H_hr.WAN_NUM,NRPTS_hrz);
                            for i = 1:NRPTS_hrz
                                vector_list_hrz(i,fin_dir) = unique_dir(i,:);
                                if Num
                                    Hnum_list_hrz(:,:,i) = sum(Hnum_kf_list_new(:,:,cutlist(i,1):cutlist(i,2)),3);
                                end
                                if Coe
                                    Hcoe_list_hrz = sym(zeros(H_hr.WAN_NUM,H_hr.WAN_NUM,NRPTS_hrz));
                                end
                            end
                    end
                    switch H_hr.Type
                        case {'sparse'}
                            H_hr.HnumL = Hnum_list_hrz;
                        case {'mat'}
                            if Num
                                H_hr.HnumL = Hnum_list_hrz;
                            end
                            if Coe
                                H_hr.HcoeL = Hcoe_list_hrz;
                            end
                    end
                    H_hr.vectorL = int32(vector_list_hrz);
            end
        end
        function [H00_H11_cell_list_1,H00_H11_cell_list_2] = H00_H11_cell_list_gen(H_hr,fin_dir,principle_layer)
            %             disp('hr_dat Hxyz wannierTOOls(wt) ');
            if nargin <3
                principle_layer = 3;
            end
            if nargin <2
                fin_dir = 2;
            end
            k_n=length(H_hr.klist_s(:,1));
            H00_H11_cell_list_1 = zeros(H_hr.WAN_NUM,H_hr.WAN_NUM,k_n); % H00
            H00_H11_cell_list_2 = zeros(H_hr.WAN_NUM,H_hr.WAN_NUM,k_n); % H01
            %hr_sparse = hr2hr_sparse(H_hr);
            pb = vasplib_tool_outer.CmdLineProgressBar('Hamiltionian Generating ');
            for j = 1:k_n
                pb.print(j,k_n);
                kpoints_f = H_hr.klist_s(j,:);
                H_hrz = H_hr.hrz_gen(kpoints_f,fin_dir);
                [H00,H01,~] = H_hrz.Green_prepare(principle_layer,fin_dir);
                H00_H11_cell_list_1(:,:,j) = H00;
                H00_H11_cell_list_2(:,:,j) = H01;
            end
            pb.delete();
        end
        function varargout = Green_prepare(H_hrz,principle_layer,fin_dir)
            if nargin <3
                fin_dir = 2;
            end
            if nargin <2
                principle_layer = 3;
            end
            % init
            WANNUM = H_hrz.WAN_NUM;
            vector_list = H_hrz.vectorL ;
            label_list = vector_list(:,fin_dir);
            %lin_000 = find(vector_list(:,fin_dir) == 0);
            Poly_priciplayer_mat = HR.Poly_priciplayer_mat_gen(principle_layer);
            %main
            switch H_hrz.Type
                case 'mat'
                    H00 = zeros(WANNUM*principle_layer);
                    for i = -(principle_layer-1):(principle_layer-1)
                        z = i;
                        label = find(label_list == z);
                        if label >0
                            H00=H00+kron(Poly_priciplayer_mat(:,:,i+principle_layer),H_hrz.HnumL(:,:,label));
                        end
                    end
                    H10 = zeros(WANNUM*principle_layer);
                    for i = 1:principle_layer
                        z = i;
                        label = find(label_list == z);
                        if label >0
                            H10=H10+kron(Poly_priciplayer_mat(:,:,i),H_hrz.HnumL(:,:,label));
                        end
                    end
                case 'sparse'
                    H00 = sparse(WANNUM*principle_layer);
                    for i = -(principle_layer-1):(principle_layer-1)
                        z = i;
                        label = find(label_list == z);
                        if label >0
                            H00=H00+kron(Poly_priciplayer_mat(:,:,i+principle_layer),H_hrz.HnumL{label});
                        end
                    end
                    H10 = sparse(WANNUM*principle_layer);
                    for i = 1:principle_layer
                        z = i;
                        label = find(label_list == z);
                        if label >0
                            H10=H10+kron(Poly_priciplayer_mat(:,:,i),H_hrz.HnumL{label});
                        end
                    end
            end
            H01 = H10' ;
            if nargout == 2
                varargout{1} = H00;
                varargout{2} = H01;
                return;
            else
                varargout{1} = H00;
                varargout{2} = H01;
            end
            switch H_hrz.Type
                case 'mat'
                    H_hrz_green.HnumL(:,:,1) = H10;H_hrz_green.vectorL(1,:) = int32([0 ,0, -1]);
                    H_hrz_green.HnumL(:,:,2) = H00;H_hrz_green.vectorL(2,:) = int32([0 ,0, 0]);
                    H_hrz_green.HnumL(:,:,3) = H01;H_hrz_green.vectorL(3,:) = int32([0 ,0, 1]);
                case 'sparse'
                    H_hrz_green.HnumL{1} = H10;H_hrz_green.vectorL(1,:) = int32([0 ,0, -1]);
                    H_hrz_green.HnumL{2}= H00;H_hrz_green.vectorL(2,:) = int32([0 ,0, 0]);
                    H_hrz_green.HnumL{3} = H01;H_hrz_green.vectorL(3,:) = int32([0 ,0, 1]);
            end
            
            varargout{3} = H_hrz_green;
        end
        % ----------------------     data update    ------------------------
        function [klist_car,klist_f,gap_list,fig] = findnodes2(H_hr,kzf_list,Noccupy,tolerance)
            if nargin < 4
                tolerance = 0.01;
            end
            if nargin <3
                Noccupy = H_hr.WAN_NUM/2;
            end
            if nargin <2
                kzf_list = 0:0.01:1;
            end
            import vasplib_tool.*;
            [fig,ax] = creat_figure();
            xlabel(ax,'kz_f');
            ylabel(ax,'min gap in k_z plane');
            title(ax,'findnode k_z plane mode');
            Gk_ = (eye(3)*2*pi)/H_hr.Rm;
            klist_car = [];
            klist_f = [];
            gap_list =[];
            for i = 1:length(kzf_list)
                H_hr.klist_s(:,3) = kzf_list(i);
                EIGENCAR = H_hr.EIGENCAR_gen();
                [bandgap,label] = min(EIGENCAR(Noccupy+1,:)- EIGENCAR(Noccupy,:));
                scatter(ax,kzf_list(i),EIGENCAR(Noccupy,label));
                scatter(ax,kzf_list(i),EIGENCAR(Noccupy+1,label));
                drawnow limitrate;
                if bandgap < tolerance
                    fprintf('find it: %7.4f eV (%6.3f, %6.3f,%6.3f)\n',bandgap,H_hr.klist_s(label,1),H_hr.klist_s(label,2),kzf_list(i));
                    klist_car = [ klist_car ;H_hr.klist_r(label,1),H_hr.klist_r(label,2),kzf_list(i)*(Gk_(3,3))];
                    klist_f = [klist_f ;H_hr.klist_s(label,1),H_hr.klist_s(label,2),kzf_list(i)];
                    gap_list =[gap_list;bandgap];
                else
                    fprintf('min gap: %7.4f eV (%6.3f, %6.3f,%6.3f)\n',bandgap,H_hr.klist_s(label,1),H_hr.klist_s(label,2),kzf_list(i));
                end
                
            end
        end
        function H_hr = Subsall(H_hr,mode)
            if nargin <2
                mode = 'num';
                %                 if strcmp(H_hr.Trpe,'sym')
                %
                %                 end
            end
            % ------ init ------- load para form inner or outer
            if exist('para.mat','file') && strcmp(mode,'file')
                %disp('The para in para.mat is first to be considered?');
                load('para.mat');
            else
                
            end
            switch mode
                case {'num','file'}
                    HcoeL_temp = subs(H_hr.HcoeL);
                    symname = symvar(HcoeL_temp);
                    if ~isempty(symname)
                        for i = 1:length(symname)
                            warning('this varible: %s is not defind',string(symname(i)));
                        end
                        disp(symname);
                        error('please introduce the varible value above or consider use sym mode');
                    else
                        if isempty(HcoeL_temp)
                            
                        else
                            H_hr.HnumL = double(HcoeL_temp);
                            H_hr.HcoeL = sym([]); % test
                        end
                        H_hr.num = true;
                        H_hr.coe = false;  % test
                    end
                    if H_hr.overlap
                        ScoeL_temp = subs(H_hr.ScoeL);
                        symname = symvar(ScoeL_temp);
                        if ~isempty(symname)
                            for i = 1:length(symname)
                                warning('this varible: %s is not defind',string(symname(i)));
                            end
                            disp(symname);
                            error('please introduce the varible value above or consider use sym mode');
                        else
                            H_hr.SnumL = double(ScoeL_temp);
                            if ~strcmp(H_hr.Type,'list')
                                H_hr.num = true;
                                H_hr.coe = false;  % test
                            end
                        end
                    end
                case 'sym'
                    H_hr.HcoeL = subs(H_hr.HcoeL);
                otherwise
                    
            end
        end
        function H_hr = input_orb_struct(H_hr,filename,mode,options)
            arguments
                H_hr HR;
                filename string ='POSCAR';
                mode char {mustBeMember(mode,{'vasp','tbsk','sym'})} = 'vasp';
                options.symbolic logical = false;
                options.Operation logical = false;
                options.warning logical = true
                options.spin  char {mustBeMember(options.spin,{'spinless','wannier','block'})}= 'spinless';
            end
            optionsCell = namedargs2cell(options);
            H_hr.Basis_num = H_hr.WAN_NUM;
            H_hr = input_orb_struct@vasplib(H_hr,filename,mode,...
                optionsCell{:});
        end
    end
    %% Gen
    methods
        function H_hr = GenfromOrth(H_hr,seed_r,seed_i,Accuracy,options)
            arguments
                H_hr HR;
                seed_r = 'gamma__r_';
                seed_i = 'gamma__i_';
                Accuracy = 1e-6;
                options.fromCvectorL = false;
            end
            nAccuracy = floor(-log(Accuracy)/log(10));
            %H_hr = H_hr.simplify(Accuracy);
            switch H_hr.Type
                case 'list'
                    if H_hr.vectorhopping
                        if  options.fromCvectorL % test
                            CL= vpa(H_hr.CvectorL(:,1:rank(H_hr.CvectorL)),nAccuracy);
                            SymVar_r = sym(seed_r,[size(CL,2),1],'real');
                            H_hr.HcoeL = CL(1:end/2,:)*SymVar_r + 1i*CL(end/2+1:end,:)*SymVar_r;
                            H_hr.vectorhopping = false;
                            return;
                        end
                        AL = vpa(H_hr.AvectorL,nAccuracy);
                        % ugly
                        %AL = vpa(real(rref(AL.').'),nAccuracy);
                        %AL(abs(AL)>Accuracy^(-1)) = 0;
                        %[U,S,D ] = svd(AL);
                        %                         nA = size(AL,2);
                        %                         [uAL,~,ic] = unique(AL,'rows');
                        %                         line0 = ismember(zeros(1,nA,class(AL)),uAL,'rows');
                        %                         LineNoneZero = ~line0;
                        %                         try
                        %                             uAL(LineNoneZero,:) = eye(nA);
                        %                         catch
                        %
                        %                         end
                        %                         AL = uAL(ic,:);
                        BL = vpa(H_hr.BvectorL,nAccuracy);
                        %BL = vpa(real(rref(BL.').'),nAccuracy);
                        %                         nB = size(BL,2);
                        %                         [uBL,~,ic] = unique(BL,'rows');
                        %                         line0 = ismember(zeros(1,nB,class(BL)),uBL,'rows');
                        %                         LineNoneZero = ~line0;
                        %                         try
                        %                             uBL(LineNoneZero,:) = eye(nB);
                        %                         catch
                        %
                        %                         end
                        %                         BL = uBL(ic,:);
                        %                        AL = AL.*conj(AL);
                        %                        BL = BL.*conj(BL);
                        SymVar_r = sym(seed_r,[size(AL,2),1],'real');
                        SymVar_i = sym(seed_i,[size(BL,2),1],'real');
                        H_hr.HcoeL = AL*SymVar_r + 1i*BL*SymVar_i;
                        H_hr.vectorhopping = false;
                        %H_hr.HcoeL = H_hr.HcoeL
                    end
            end
        end
        function hrdat = Gen_hr(H_hr,filename,mode)
            if nargin < 3
                mode = 'hr_dat';
            end
            if nargin <2
                filename = 'wannier90_hr.dat';
            end
            if strcmp(H_hr.Type,'sparse')
                H_hr  = H_hr.full();
            end
            if strcmp(H_hr.Type,'list')
                H_hr  = H_hr.rewind();
            end
            hrdat=fopen(filename ,'w');
            pb = vasplib_tool_outer.CmdLineProgressBar('Writing hr.dat - NRPT:');
            switch mode
                case 'hr_dat'
                    % write title Nbands and NRPTS
                    date_=date;
                    fprintf(hrdat,"%s\n",date_);
                    %1
                    %2
                    fprintf(hrdat,"         %d\n",H_hr.WAN_NUM);
                    %2
                    %3
                    fprintf(hrdat,"         %d\n",H_hr.NRPTS);
                    %3
                    %4
                    % write NRPT_list
                    % NRPTS_num1=fix(NRPTS/15);
                    % NRPTS_num2=rem(NRPTS,15);
                    COUNT_FLAG=0;
                    for i =1:H_hr.NRPTS
                        fprintf(hrdat,"    %d",1);
                        COUNT_FLAG=COUNT_FLAG+1;
                        if COUNT_FLAG == 15
                            COUNT_FLAG = 0;
                            fprintf(hrdat,"\n");
                        end
                    end
                    %another \n
                    fprintf(hrdat,"\n");
                    % write hopping terms
                    
                    for i=1:H_hr.NRPTS
                        pb.print(i,H_hr.NRPTS,' ...');
                        %fprintf('Wrinting (%d/%d)  NRPT\n',i,H_hr.NRPTS);
                        for k=1:H_hr.WAN_NUM
                            %fprintf('Wrinting %d th WAN_ORB \n',k);
                            for j=1:H_hr.WAN_NUM
                                fprintf(hrdat,"    %d    %d    %d",H_hr.vectorL(i,1),H_hr.vectorL(i,2),H_hr.vectorL(i,3));
                                fprintf(hrdat,"    %d    %d",j,k);
                                real_part=real(H_hr.HnumL(j,k,i));
                                imag_part=imag(H_hr.HnumL(j,k,i));
                                fprintf(hrdat,"    %.7f    %.7f\n",real_part,imag_part);
                            end
                        end
                    end
                    fclose(hrdat);
                otherwise
            end
        end
        function [H_hr_bk] = POSCAR_gen(H_hr,filename,options)
            arguments
                H_hr HR;
                filename = 'POSCAR_gen.vasp';
                options.title = "PARCHG Gen by vasplib. " + string(date);
                options.a_crystal_constance = 1;
                options.vacuum = false;
                options.fin_dir_list = [1 1 0];
                options.vacuum_length = 10;
            end
            if options.vacuum
                [H_hr.orbL,H_hr.Rm] = vasplib.AddVacuumLayer(H_hr.orbL,H_hr.Rm,options.fin_dir_list,'vacuum_length',options.vacuum_length);
            end
            H_hr_bk = H_hr;
            quantumList = H_hr.quantumL;
            H_hr = H_hr.reseq(quantumList(:,end)>=0 |isnan(quantumList(:,end)));
            elementList = H_hr.elementL;
            [elementList,sortL] = sort(elementList);
            H_hr = H_hr.reseq(sortL);
            orbList = H_hr.orbL;
            %--------  init  --------
            try
                elements = readtable('elements.txt');
            catch
                elements = [];
            end
            if ~isempty(elements)
                %elements.Properties.RowNames = elements.atom_number;
            end
            unique_elementL = unique(elementList);
            AtomTypes = length(unique_elementL);
            Atom_num = zeros(1,AtomTypes);
            Atom_name = "";
            Atom_name = repmat(Atom_name,size(Atom_num));
            Rm = H_hr.Rm;
            for i = 1:AtomTypes
                atom_number = unique_elementL(i);
                Atom_num(i) = sum(elementList == atom_number);
                rows = elements.atom_number ==atom_number;
                Atom_name(i) = string((elements.atom_symbol(rows)));
            end
            %  write POSCAR
            % Initialize variables
            fileID = fopen(filename,'w');
            %
            fprintf(fileID,"%s\n",options.title);
            fprintf(fileID,"%d\n",options.a_crystal_constance);
            %fprintf(fileID,"  ",Rm(i,j));
            for i=1:3
                for j=1:3
                    fprintf(fileID,"  %f",Rm(i,j));
                end
                fprintf(fileID,"\n");
            end
            for i=1:length(Atom_name)
                fprintf(fileID,"%s ",Atom_name(i));
            end
            fprintf(fileID,"\n  ");
            for i=1:length(Atom_num)
                fprintf(fileID,"%d ",Atom_num(i));
            end
            fprintf(fileID,"\n");
            fprintf(fileID,"Direct\n  ");
            % sites
            for i=1:size(orbList,1)
                fprintf(fileID,"%f  ",mod(orbList(i,1),1));
                fprintf(fileID,"%f  ",mod(orbList(i,2),1));
                fprintf(fileID,"%f  ",mod(orbList(i,3),1));
                %         if ~strcmp(string(sites(i).name),"")
                %             %fprintf(fileID,"%s\n  ",sites(i).name);
                %             fprintf(fileID,"\n  ");
                %         else
                %             fprintf(fileID,"\n  ");
                %         end
                fprintf(fileID,"\n  ");
            end
            fclose(fileID);
            [H_hr_bk.Rm,H_hr_bk.sites,H_hr_bk.Atom_name,H_hr_bk.Atom_num,~]=POSCAR_readin(filename);
        end
        function H_hr = tbbox_in_gen(H_hr,options)
            arguments
                H_hr HR;
                options.Accuracy = 1e-6;
                options.OperL = [];
                options.hr = true;
                options.wan = false;
                options.lda = true;
            end
            import spglib_matlab.*;
            Accuracy = -log10(options.Accuracy);
            %fprintf('Reseq the orbital by your self');
            [unique_orbL,unique_seqL,unique_seqL_inv] = unique(H_hr.orbL,'stable','rows');
            unique_orbL_counts = accumarray(unique_seqL_inv,1);
            lattice = H_hr.Rm;
            position = unique_orbL.';
            types = H_hr.quantumL(unique_seqL,1);
            if isempty(options.OperL )
                try
                    SpglibDataset  = spglib_matlab.spg_get_dataset_from_sites(lattice,position,types);
                    rotations = double(SpglibDataset.rotations);
                    translations = double(SpglibDataset.translations);
                    n_operations = SpglibDataset.n_operations;
                catch
                    rotations = eye(3);
                    translations = zeros(1,3);
                    n_operations = 1;
                end
            else
                n_operations = length(options.OperL);
                rotations = zeros(3,3,n_operations);
                translations  = zeros(1,3,n_operations);
                for i = 1:n_operations
                    rotations(:,:,i) = options.OperL(i).R;
                    translations(:,:,i) = options.OperL(i).t;
                end
            end
            sym_oper_list = zeros(n_operations,9);
            for i = 1:n_operations
                sym_oper_list(i,1) = i;
                sym_oper_list(i,2) = det(rotations(:,:,i));
                [n,theta]= spglib_matlab.Rotation2nTheta(rotations(:,:,i),H_hr.Rm);
                sym_oper_list(i,3) = roundn(theta,-Accuracy);
                sym_oper_list(i,4) = n(1);
                sym_oper_list(i,5) = n(2);
                sym_oper_list(i,6) = n(3);
                sym_oper_list(i,7) =roundn(translations(1,1,i),-Accuracy);
                sym_oper_list(i,8) =roundn(translations(1,2,i),-Accuracy);
                sym_oper_list(i,9) =roundn(translations(1,3,i),-Accuracy);
            end
            %
            if sum(H_hr.quantumL(:,4)) > 0 || options.lda
                casename  = 'lda';
                SOCcounts = 1;
            else
                casename  = 'soc';
                SOCcounts = 2;
            end
            filename_hr = casename+"_hr.dat";
            ntau = length(unique_seqL);
            %% begin write
            fileID = fopen('tbbox.in','w');
            %
            fprintf(fileID,' case = %s !lda or soc ! spinless of spinful ！dont use soc or consider\n\n',casename);
            %
            fprintf(fileID,' proj:\n');
            %fprintf(fileID,'orbt = 1  !Orbital convertion 1: s, px, py, pz, xy, yz, zx, x2-y2, 3z2-r2 !!!\n');
            fprintf(fileID,' orbt = 1\n');
            %fprintf(fileID,'ntau = %d !Orbital convertion 2: s, pz, px, py, 3z2-r2, xz, yz, x2-y2, xy (Wannier90) !!!\n',ntau);
            fprintf(fileID,' ntau = %d\n',ntau);
            for i  = 1:ntau
                fprintf(fileID,' %11.8f %11.8f %11.8f %2d %2d',...
                    unique_orbL(i,1),unique_orbL(i,2),unique_orbL(i,3),...
                    H_hr.quantumL(unique_seqL(i),1),unique_orbL_counts(i)/SOCcounts);
                %                 switch i
                switch i
                    case 1
                        fprintf(fileID,' !! x1, x2, x3, itau, iorbit \n');
                    case 2
                        fprintf(fileID,' !! ! itau 指标第几个原子,iorbit 指标什么类型的轨道 \n');
                    case 3
                        fprintf(fileID,' ! 1 s (ir2tb_pz 1 pz)  \n');
                    case 4
                        fprintf(fileID,' ! 2 (ir2tb_pz 2 s pz)  \n');
                    case 5
                        fprintf(fileID,' ! 3 px py pz  \n');
                    case 6
                        fprintf(fileID,' ! 4 s px py pz  \n');
                    case 7
                        fprintf(fileID,' ! 5 d  \n');
                    case 8
                        fprintf(fileID,' ! 6 s d \n');
                    case 9
                        fprintf(fileID,' ! 7 f  \n');
                    case 10
                        fprintf(fileID,' ! 8 p d  \n');
                    case 11
                        fprintf(fileID,' ! 9 s p d  \n');
                    otherwise
                        fprintf(fileID,'\n');
                end
            end
            fprintf(fileID,' end projections\n\n');
            %%
            fprintf(fileID,' kpoint:\n');
            
            fprintf(fileID,' kmesh = 1\n');
            fprintf(fileID,' Nk = 8  !8 lines\n');
            for i =1:8
                switch i
                    case 1
                        b1 =0;b2 =0;b3 =0;
                    case 2
                        b1 =0.5;b2 =0;b3 =0;
                    case 3
                        b1 =0;b2 =0.5;b3 =0;
                    case 4
                        b1 =0.5;b2 =0.5;b3 =0;
                    case 5
                        b1 =0;b2 =0;b3 =0.5;
                    case 6
                        b1 =0.5;b2 =0;b3 =0.5;
                    case 7
                        b1 =0;b2 =0.5;b3 =0.5;
                    case 8
                        b1 =0.5;b2 =0.5;b3 =0.5;
                end
                fprintf(fileID,' %11.8f %11.8f %11.8f !k%d :b1 b2 b3 \n',b1,b2,b3,i);
            end
            fprintf(fileID,' end kpoint_path\n\n');
            %% unitcell
            Gk = eye(3)/H_hr.Rm;
            fprintf(fileID,' unit_cell:\n');
            for i =1:3
                fprintf(fileID,'    %12.9f %12.9f %12.9f    %12.9f %12.9f %12.9f \n',...
                    H_hr.Rm(i,1),H_hr.Rm(i,2),H_hr.Rm(i,3),Gk(1,i),Gk(2,i),Gk(3,i));
            end
            %% spacegroup operators
            for i =1 :size(sym_oper_list,1)
                fprintf(fileID,' %4d ',sym_oper_list(i,1));
                %fprintf('%2d ',sym_oper_list(i,1));
                %11.7f %11.7f %11.7f  %11.7f %11.7f %11.7 %11.7f %11.7f\n',...
                for j =2:9
                    fprintf(fileID,' %11.8f ',sym_oper_list(i,j));
                    %fprintf('%11.7f ',sym_oper_list(i,j));
                end
                fprintf(fileID,'\n');
            end
            fprintf(fileID,' end unit_cell_cart\n');
            if options.hr
                if options.wan
                    Norb = H_hr.WAN_NUM;
                    H_hr = H_hr.reseq([(1:Norb/2)*2-1,(1:Norb/2)*2]);
                end
                H_hr.Gen_hr(filename_hr);
            end
            fclose(fileID);
        end
        function H_hr = pythtb_gen(H_hr,filename,options)
            arguments
                H_hr HR ;
                filename char='HRTB.py'
                options.sym = true;
                options.band logical= true;
                options.python_env = '/usr/bin/python';
                options.KPOINTS = 'KPOINTS';
                options.WilsonLoop = false;
            end
            % check
            H_hr = H_hr.rewrite();
            % write
            fileID = fopen(filename,'w');
            %
            fprintf(fileID,['#!',options.python_env ,'\n']);
            fprintf(fileID,'from pythtb import *\n');
            % lattice information
            fprintf(fileID,'# lattice vectors and orbital positions\n');
            % Rm
            latstring = "lat=[";
            for i = 1:3
                latstring = latstring+vasplib.mat2str_python(H_hr.Rm(i,:))+", ";
            end
            latstring = latstring + "]";
            fprintf(fileID,latstring+'\n');
            % orb
            orbstring = "orb=[\n";
            for i =1:H_hr.WAN_NUM
                orbstring = orbstring + vasplib.mat2str_python(H_hr.orbL(i,:))+", \n";
            end
            orbstring = orbstring + "]";
            fprintf(fileID,orbstring+'\n');
            fprintf(fileID,'\n');
            %
            fprintf(fileID,'# three-dimensional tight-binding model from HR\n');
            fprintf(fileID,['HRTB=tb_model(',num2str(3),',',num2str(3),', lat, orb)\n']);
            fprintf(fileID,'\n');
            fprintf(fileID,'# define hopping between orbitals\n');
            % hopping homecell
            SelectLabel_home = all([H_hr.vectorL(:,1:3) == 0,H_hr.vectorL(:,5)>=H_hr.vectorL(:,4)],2);
            SelectLabel_other = all(H_hr.vectorL(:,1:3) >=0,2) & ~all(H_hr.vectorL(:,1:3) == 0,2);
            SelectLabel = SelectLabel_home|SelectLabel_other;
            SelectVector=H_hr.vectorL(SelectLabel,:);
            if options.sym
                SelectHop=H_hr.HcoeL(SelectLabel);
                for i = 1:length(H_hr.symvar_list)
                    fprintf(fileID,string(SelectHop(i))+"= 1 \n");
                end
            else
                SelectHop=H_hr.HnumL(SelectLabel);
            end
            for i = 1:numel(SelectHop)
                fprintf(fileID,'HRTB.set_hop('+string(SelectHop(i))+","+...
                    num2str(SelectVector(i,4)-1)+","+num2str(SelectVector(i,5)-1)+","+...
                    vasplib.mat2str_python(SelectVector(i,1:3))+')\n');
            end
            fprintf(fileID,'\n');
            % band
            if options.band
                [kpoints,nodes,kpoints_name_tmp] = vasplib.KPOINTS_read(options.KPOINTS);
                nkline = length(kpoints_name_tmp)-1;
                kpoints_f = kpoints(1+(0:nkline-1)*2,:);
                kpoints_f = [kpoints_f;kpoints(end,:)];
                fprintf(fileID,'# solve model on a path in k-space\n');
                fprintf(fileID,'k=[');
                for i = 1:size(kpoints_f,1)
                    fprintf(fileID,[vasplib.mat2str_python(kpoints_f(i,:)),',']);
                end
                fprintf(fileID,']\n');
                fprintf(fileID,['(k_vec,k_dist,k_node)=HRTB.k_path(k,',num2str(nodes),')\n']);
                fprintf(fileID,'evals=HRTB.solve_all(k_vec)\n');
                
                fprintf(fileID,'# plot bandstructure\n');
                fprintf(fileID,'import matplotlib.pyplot as plt\n');
                fprintf(fileID,'fig, ax = plt.subplots()\n');
                fprintf(fileID,['for i in range(1,',num2str(H_hr.WAN_NUM+1),'):\n']);
                fprintf(fileID,'    ax.plot(k_dist,evals[i-1,:])\n');
                fprintf(fileID,'ax.set_xticks(k_node)\n');
                fprintf(fileID,['ax.set_xticklabels(',vasplib.mat2str_python(kpoints_name_tmp),')\n']);
                fprintf(fileID,'ax.set_xlim(k_node[0],k_node[-1])\n');
                fprintf(fileID,'fig.savefig("band.png")\n');
            end
            % WilsonLoop
            if options.WilsonLoop
                fprintf(fileID,'#calculate my-array\n');
                fprintf(fileID,'my_array=wf_array(HRTB,[41,41,41])\n');
                fprintf(fileID,'# solve model on a regular grid, and put origin of\n');
                fprintf(fileID,'# Brillouin zone at [0,0,0]  point\n');
                fprintf(fileID,'my_array.solve_on_grid([0,0,0])\n');
                fprintf(fileID,'\n');
                fprintf(fileID,'# calculate Berry phases around the BZ in the k_x direction\n');
                fprintf(fileID,'# (which can be interpreted as the 1D hybrid Wannier centers\n');
                fprintf(fileID,'# in the x direction) and plot results as a function of k_y\n');
                fprintf(fileID,'#\n');
                fprintf(fileID,'# Following the ideas in\n');
                fprintf(fileID,'#   A.A. Soluyanov and D. Vanderbilt, PRB 83, 235401 (2011)\n');
                fprintf(fileID,'#   R. Yu, X.L. Qi, A. Bernevig, Z. Fang and X. Dai, PRB 84, 075119 (2011)\n');
                fprintf(fileID,'# the connectivity of these curves determines the Z2 index\n');
                fprintf(fileID,'#\n');
                fprintf(fileID,['wan_cent = my_array.berry_phase(',...
                    vasplib.mat2str_python(1:H_hr.WAN_NUM/2),...
                    ',dir=1,contin=False,berry_evals=True)\n']);
                fprintf(fileID,'wan_cent/=(2.0*np.pi)\n');
                fprintf(fileID,'\n');
                fprintf(fileID,'nky=wan_cent.shape[0]\n');
                fprintf(fileID,'ky=np.linspace(0.,1.,nky)\n');
                fprintf(fileID,'# draw Wannier center positions\n');
                fprintf(fileID,'fig, ax2 = plt.subplots()\n');
                fprintf(fileID,['for iband in range(0,',num2str(H_hr.WAN_NUM/2),'):\n']);
                fprintf(fileID,'    ax2.plot(ky,wan_cent[iband,:,0],"k.")\n');
                fprintf(fileID,'ax2.set_ylim(-1.0,1.0)\n');
                fprintf(fileID,'ax2.set_ylabel(''Wannier center along x'')\n');
                fprintf(fileID,'ax2.set_xlabel(''k_y'')\n');
                fprintf(fileID,'ax2.set_xticks([0.0,0.5,1.0])\n');
                fprintf(fileID,'ax2.set_xlim(0.0,1.0)\n');
                fprintf(fileID,'ax2.set_xticklabels([r"$0$",r"$\\pi$", r"$2\\pi$"])\n');
                fprintf(fileID,'ax2.axvline(x=.5,linewidth=0.5, color=''k'')\n');
                fprintf(fileID,'ax2.set_title("1D Wannier centers phase")\n');
                fprintf(fileID,'\n');
                fprintf(fileID,'fig.tight_layout()\n');
                fprintf(fileID,'fig.savefig("wcc.pdf")\n');
            end
        end
    end
    %% ChangeClass
    methods
        function HcktObj = HR2Hckt(H_hr,options,options_homecell)
            arguments
                H_hr
                options.title = 'HcktFromHR';
                options.dim = 2;
                options.autominus = false;
                options.magnitude = 'p';
                options.mode {mustBeMember(options.mode,{'real','sigma'})}= 'real';
                options_homecell.homecell = "normal";
            end
            %
            VarC0=1;
            Var2C0=2;
            VarR0=1;
            VarR0_2=0.5;
            switch options.magnitude
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
            % Declare the primitive of Octagraphene
            TITLE = options.title;
            Dim = options.dim;
            % 
            if options.autominus
                warning off;
                for i = 1:H_hr.WAN_NUM
                    BaseOnsite = sum(H_hr.HnumL(H_hr.vectorL(:,4)==i));%H_hr.HnumL(ismember(H_hr.vectorL,[0,0,0,i,i],'rows')) +
                    H_hr = H_hr.set_hop(BaseOnsite,...
                        i,i,[0,0,0],'set');
                end
            end
            %  declare a Hckt
            HomeVector = zeros(1,Dim);
            NBAND= double(H_hr.WAN_NUM);
            switch options.mode
                case 'real'
                    H_hr = H_hr.ForceTolist();
                    HcktObj = Hckt('title',TITLE,'Nports',NBAND,'vectorL',HomeVector,'magnitude',options.magnitude);
                    fprintf('Limitations: Only Numerical HR, real hopping support.\n');
                    %% set_homecell
                    HomeCell = Subckt.FromHomecellList(H_hr,'magnitude',options.magnitude);
                    HcktObj = HcktObj.set_home(HomeCell,1:round(NBAND/2),round(NBAND/2)+1:H_hr.WAN_NUM);
                    %% set_hop
                    Cplus  = Subckt('Xplus 1 2  C') ;
                    Cminus  = Subckt('Xminus 1 2  minusC') ;
                    vectorList = double(H_hr.vectorL);
                    for n = 1:H_hr.NRPTS
                        Rvector = vectorList(n,1:Dim);
                        if isequal(Rvector,HomeVector)
                            continue;
                        end
                        i = vectorList(n,4);
                        j = vectorList(n,5);
                        % vector,Subcktobj,PortInL,PortOutL,DescriptionL
                        if H_hr.HnumL(n) > 0
                            HoppingSckt = Cminus;
                        elseif H_hr.HnumL(n) < 0
                            HoppingSckt = Cplus;
                        end
                        HcktObj = HcktObj.set_hop(Rvector,HoppingSckt,i,j,['C_hopping = ',num2str(abs(H_hr.HnumL(n)*100)),options.magnitude]);
                    end
                case 'sigma'
                    H_hr = H_hr.ForceToMat();
                    % check
                    VarC0 = 1;Var
                    BasisC3_origin   =Subckt('XBasisC3_origin   l1 l2 l3 r1 r2 r3 TOGND BasisC3_origin  ','magicnumber',8);
                    HoppingDist{1}   =Subckt('XPlusSigma0       l1 l2 l3 r1 r2 r3 TOGND PlusSigma0      ','magicnumber',8);
                    HoppingDist{2}   =Subckt('XMinusSigma0      l1 l2 l3 r1 r2 r3 TOGND MinusSigma0     ','magicnumber',8);
                    HoppingDist{3}   =Subckt('XPlusiSigma0      l1 l2 l3 r1 r2 r3 TOGND PlusiSigma0     ','magicnumber',8);
                    HoppingDist{4}   =Subckt('XMinusiSigma0     l1 l2 l3 r1 r2 r3 TOGND MinusiSigma0    ','magicnumber',8);
                    HoppingDist{5}   =Subckt('XPlusSigma1       l1 l2 l3 r1 r2 r3 TOGND PlusSigma1      ','magicnumber',8);
                    HoppingDist{6}   =Subckt('XMinusSigam1      l1 l2 l3 r1 r2 r3 TOGND MinusSigam1     ','magicnumber',8);
                    HoppingDist{7}   =Subckt('XPlusiSigma1      l1 l2 l3 r1 r2 r3 TOGND PlusiSigma1     ','magicnumber',8);
                    HoppingDist{8}   =Subckt('XMinusiSigma1     l1 l2 l3 r1 r2 r3 TOGND MinusiSigma1    ','magicnumber',8);
                    HoppingDist{9}   =Subckt('XPlusGen3Sigma2   l1 l2 l3 r1 r2 r3 TOGND PlusGen3Sigma2  ','magicnumber',8);
                    HoppingDist{10}  =Subckt('XMinusGen3Sigma2  l1 l2 l3 r1 r2 r3 TOGND MinusGen3Sigma2 ','magicnumber',8);
                    HoppingDist{11}  =Subckt('XPlusiGen3Sigma2  l1 l2 l3 r1 r2 r3 TOGND PlusiGen3Sigma2 ','magicnumber',8);
                    HoppingDist{12}  =Subckt('XMinusiGen3Sigma2 l1 l2 l3 r1 r2 r3 TOGND MinusiGen3Sigma2','magicnumber',8);
                    HoppingDist{13}  =Subckt('XPlusGen3Sigma3   l1 l2 l3 r1 r2 r3 TOGND PlusGen3Sigma3  ','magicnumber',8);
                    HoppingDist{14}  =Subckt('XMinusGen3Sigma3  l1 l2 l3 r1 r2 r3 TOGND MinusGen3Sigma3 ','magicnumber',8);
                    HoppingDist{15}  =Subckt('XPlusiGen3Sigma3  l1 l2 l3 r1 r2 r3 TOGND PlusiGen3Sigma3 ','magicnumber',8);
                    HoppingDist{16}  =Subckt('XMinusiGen3Sigma3 l1 l2 l3 r1 r2 r3 TOGND MinusiGen3Sigma3','magicnumber',8);
                    if NBAND ~= 2
                        fprintf('Choosing wrong mode!\n');
                        error('Nband ~= 2!!!!!');
                    end
                    HcktObj = Hckt('title',TITLE,'Nports',NBAND+1,'vectorL',HomeVector,'magnitude',options.magnitude);
                    fprintf('Limitations: Only Numerical HR, Two band model support.\n');
                    %% set_homecell
                    switch options_homecell.homecell
                        case 'normal'
                            HomeCell = BasisC3_origin;
                        otherwise
                            HomeCell = BasisC3_origin;
                    end
                    HcktObj = HcktObj.set_home(HomeCell,1:3,4:6);
                    %% set_hop
                    vectorList = double(H_hr.vectorL);
                    HnumList = double(H_hr.HnumL);
                    for n = 1:H_hr.NRPTS
                        Rvector = vectorList(n,1:Dim);
                        if isequal(Rvector,HomeVector)
                            continue;
                        end
                        [CoeForPauli] = pauliDecompositionNumerial(HnumList(:,:,n));
                        reaLCoeForPauli = real(CoeForPauli);
                        imagCoeForPauli = imag(CoeForPauli);
                        Coe16 = zeros(ones(1,16));
                        PlusReal  = [1,5,9,13];
                        MinusReal = [2,6,10,14];
                        PlusImag  = [3,7,11,15];
                        MinusImag = [4,8,12,16];
                        Coe16( PlusReal(reaLCoeForPauli>0)) =    reaLCoeForPauli(reaLCoeForPauli>0);
                        Coe16(MinusReal(reaLCoeForPauli<0)) =   -reaLCoeForPauli(reaLCoeForPauli<0);
                        Coe16( PlusImag(imagCoeForPauli>0)) =    imagCoeForPauli(imagCoeForPauli>0);
                        Coe16(MinusImag(imagCoeForPauli<0)) =   -imagCoeForPauli(imagCoeForPauli<0);
                        % vector,Subcktobj,PortInL,PortOutL,DescriptionL
                        for i = 1:16
                            if Coe16(i) == 0
                            else
                                HcktObj = HcktObj.set_hop(Rvector,HoppingDist{i},[1 2 3],[4 5 6],...
                                    [' VarC0 = ',num2str(abs(Coe16(i)*VarC0)),Cmagnitude,...
                                     ' Var2C0 = ',num2str(abs(Coe16(i)*Var2C0)),Cmagnitude,...
                                     ' VarR0 = ',num2str(abs(Coe16(i)*VarR0)),Rmagnitude,...
                                     ' VarR0_2 = ',num2str(abs(Coe16(i)*VarR0_2)),Rmagnitude ...
                                        ]);
                            end
                            
                        end
                    end
            end

        end
        function H_hr_forHckt = HRforHckt(H_hr,options)
            arguments
                H_hr HR;
                options.fast = false
                options.direct = false;
                options.Accuracy = 1e-6;
                options.C_0 = 1;
            end
            if strcmp(H_hr.Type,'mat')
                H_hr = H_hr.rewrite();
            end
            H_hr_forHckt = H_hr;
            if H_hr.coe
                H_hr_forHckt.HcoeL = -H_hr_forHckt.HcoeL;
                % make symbolic term > 0
                %SYMVAR = H_hr.symvar_list;
                %for i = 1:length(SYMVAR)
                %    assume(SYMVAR(i)>0);
                %end
                C_0 = sym('C_0','real');
                assume(C_0>0);
                BaseOnsiteL = repmat(C_0,[H_hr_forHckt.WAN_NUM,1]);
                maxOnsite = C_0;
                % set onsite
                for i = 1:H_hr_forHckt.WAN_NUM
                    maxOnsite = max(BaseOnsiteL(i)-sum(H_hr_forHckt.HcoeL(H_hr_forHckt.vectorL(:,4)==i)),maxOnsite);
                end
                for i = 1:H_hr_forHckt.WAN_NUM
                    BaseOnsiteL(i) = maxOnsite + sum(H_hr_forHckt.HcoeL(H_hr_forHckt.vectorL(:,4)==i));
                    H_hr_forHckt = H_hr_forHckt.set_hop(maxOnsite,...
                        i,i,[0,0,0],'symadd');
                end
            else
                H_hr_forHckt.HnumL = -H_hr_forHckt.HnumL;
                % make symbolic term > 0
                C_0 = options.C_0 ;
                BaseOnsiteL = repmat(C_0,[H_hr_forHckt.WAN_NUM,1]);
                maxOnsite = C_0;
                % set onsite
                for i = 1:H_hr_forHckt.WAN_NUM
                    maxOnsite = max(BaseOnsiteL(i)-sum(H_hr_forHckt.HnumL(H_hr_forHckt.vectorL(:,4)==i)),maxOnsite);
                end
                for i = 1:H_hr_forHckt.WAN_NUM
                    BaseOnsiteL(i) = maxOnsite + sum(H_hr_forHckt.HnumL(H_hr_forHckt.vectorL(:,4)==i));
                    H_hr_forHckt = H_hr_forHckt.set_hop(maxOnsite,...
                        i,i,[0,0,0],'add');
                end
            end
        end
        function H_htrig = HR2Htrig(H_hr,options)
            arguments
                H_hr HR;
                options.sym =true;
                options.num = false;
                options.fast = false
                options.direct = false;
                options.Accuracy = 1e-6;
                options.Type {mustBeMember(options.Type,{'exp','sincos','mat','list'})}='sincos' ;
            end
            if H_hr.num
                options.num = true;
                options.sym = false;
            end
            if  options.sym && ~options.fast
                H_hr.Rm = sym(H_hr.Rm);
                H_hr.orbL = sym(H_hr.orbL);
                fprintf('please check whether the Rm and orbL are properly set.\n');
                disp(H_hr.Rm);
                disp(H_hr.orbL);
            end
            if options.direct
                Hsym = sym(H_hr);
                H_htrig = Htrig(Hsym);
                H_htrig = H_htrig.vasplibCopy(H_hr);
            elseif options.fast
                if ~strcmp(H_hr.Type,'list')
                    H_hr = H_hr.rewrite('Accuracy',options.Accuracy);
                    %H_hr = H_hr.simplify();
                else
                    H_hr = H_hr.simplify(options.Accuracy);
                end
                H_htrig = Htrig(H_hr.WAN_NUM,'Type','list');
                H_htrig = H_htrig.vasplibCopy(H_hr);
                %H_htrig = H_htrig.tjmti_gen();
                tiL = H_htrig.orbL(H_hr.vectorL(:,4),:);
                tjL = H_htrig.orbL(H_hr.vectorL(:,5),:);
                RtjmtiL = (double(H_hr.vectorL(:,1:3))+ tjL-tiL); %
                % H_htrig.vectorL(:,1:3)+ tjL-tiL use this convention
                ijL = double(H_hr.vectorL(:,[4,5]));
                if options.num
                    H_htrig.HnumL = H_hr.HnumL;
                    H_htrig.HsymL_numL = [RtjmtiL*double(H_htrig.Rm),ijL] ;
                    H_htrig.num = true;
                    H_htrig.coe = false;
                end
                if options.sym
                    H_htrig.HcoeL = H_hr.HcoeL;
                    H_htrig.HsymL_coeL = [RtjmtiL*H_htrig.Rm,ijL] ;
                    H_htrig.coe = true;
                    H_htrig.num = false;
                end
                if ~strcmp(options.Type,'list')
                    H_htrig = H_htrig.rewrite();
                end
            else
                H_htrig = Htrig(H_hr.WAN_NUM,'Type',options.Type );
                H_htrig = H_htrig.vasplibCopy(H_hr);
                if options.num
                    H_hr.HcoeL = sym(H_hr.HnumL);
                end
                H_htrig = H_htrig.tjmti_gen('sym');
                syms k_x k_y k_z real;
                vectorList = double(H_hr.vectorL);
                pb = vasplib_tool_outer.CmdLineProgressBar('Transforming ');
                if strcmp(H_hr.Type,'list')
                    tji_mat_cart = H_htrig.tjmti{1};
                    NRPTS_ = H_hr.NRPTS;
                    for k = 1:NRPTS_
                        pb.print(k,NRPTS_,' th hopping into Htrig obj ...');
                        SymHopping = H_hr.HcoeL(k);
                        if SymHopping ~= sym(0)
                            i = vectorList(k,4);
                            j = vectorList(k,5);
                            % $H_{i j}^{\mathbf{k}}=\sum_{\mathbf{R}} e^{i \mathbf{k} \cdot\left(\mathbf{R}+\mathbf{t}_{j}-\mathbf{t}_{i}\right)} H_{i j}(\mathbf{R})$
                            SymVar = exp(...
                                1i...
                                *[k_x k_y k_z]*((...
                                vectorList(k,1:3)*H_hr.Rm)'+reshape(tji_mat_cart(i,j,:),[3,1])...
                                ));
                            %SymVar = simplify(rewrite(SymVar,'sincos'));
                            H_htrig = H_htrig.set_hop(SymHopping,SymVar,[i,j]);
                            %zeros_mat = zeros(H_htrig.Basis_num);
                            %zeros_mat(j,i) = 1; % an enforced transpose % bug!!!!
                            %H_htrig = H_htrig.setup_rough(expand(SymHopping*SymVar),zeros_mat,true);
                            %H_htrig = H_htrig.set_hop([i,j]);
                        end
                    end
                else
                    tji_mat_kcart = H_htrig.tjmti{3};
                    NRPTS_ = H_hr.NRPTS;
                    exp_preL = exp(1i...
                        *[k_x k_y k_z]*((...
                        vectorList*H_hr.Rm)'));
                    WAN_WAN = numel(tji_mat_kcart);
                    sizeWAN = size(tji_mat_kcart);
                    sym0 = sym(0);
                    [iL,jL] = ind2sub(sizeWAN,1:WAN_WAN);
                    for i = 1:NRPTS_
                        pb.print(i,NRPTS_,' th H(NRPT) into Htrig obj ...');
                        SymHoppingMat = H_hr.HcoeL(:,:,i);
                        exp_preL_i  = exp_preL(i);
                        for j = 1:WAN_WAN
                            SymHopping = SymHoppingMat(j);
                            if SymHopping ~= sym0
                                exp_pre = exp_preL_i*tji_mat_kcart(j);
                                SymVar = simplify(rewrite(exp_pre,'sincos'));
                                H_htrig = H_htrig.set_hop(SymHopping,SymVar,[iL(j),jL(j)]);
                            end
                        end
                    end
                end
                pb.delete();
            end
        end
        function H_hk = HR2HK(H_hr,kpoints_frac,options)
            arguments
                H_hr HR;
                kpoints_frac = [0,0,0];
                options.sym =false;
                options.Order = 2;
                options.convention = 'I';
            end
            if strcmp(H_hr.Type,'list')
                % waiting
                H_hr = H_hr.rewrite('rewind',true);
                H_hk = HR2HK(H_hr,kpoints_frac,...
                    'sym',options.sym,...
                    'Order',options.Order,...
                    'convention',options.convention);
            else
                syms k_x k_y k_z real;
                H_hk = HK(H_hr.WAN_NUM,options.Order);
                H_hk = H_hk.vasplibCopy(H_hr);
                %WANNUM = H_hr.WAN_NUM;
                if  options.sym
                    H_hr.Rm = sym(H_hr.Rm);
                    H_hr.orbL = sym(H_hr.orbL);
                    fprintf('please check whether the Rm and orbL are properly set.');
                    disp(H_hr.Rm);
                    disp(H_hr.orbL);
                end
                kpoints_r = kpoints_frac *H_hr.Gk;
                vectorList = double(H_hr.vectorL);
                if strcmp(options.convention,'II')
                    for n = 1:H_hr.NRPTS
                        kpoints_phase = exp(1i*2*pi*(vectorList(n,:))*kpoints_frac.');
                        tmp_HsymL_trig = exp(1i*(vectorList(n,:))*H_hk.Rm*[k_x k_y k_z].');
                        symbolic_polynomial = taylor(tmp_HsymL_trig,[k_x k_y k_z],'Order',options.Order);
                        % debug
                        H_hk = H_hk.setup_rough(symbolic_polynomial,H_hr.HcoeL(:,:,n)*kpoints_phase);
                    end
                elseif strcmp(options.convention,'I')
                    H_hr = H_hr.tjmti_gen('sym');
                    if options.sym
                        tji_mat =  H_hr.tjmti{1};
                    else
                        tji_mat = double(H_hr.tjmti{1});
                    end
                    for n = 1:H_hr.NRPTS
                        for i = 1:H_hr.WAN_NUM
                            for j = 1:H_hr.WAN_NUM
                                tjitmp = reshape(tji_mat(i,j,:),[1 3]);
                                R_add_t_vector = (vectorList(n,1:3)*H_hr.Rm+tjitmp).';
                                kpoints_phase = exp(1i*kpoints_r*(R_add_t_vector));
                                tmp_HsymL_trig = exp(1i*[k_x k_y k_z]*(R_add_t_vector));
                                symbolic_polynomial = taylor(tmp_HsymL_trig,[k_x k_y k_z],'Order',options.Order);
                                % debug
                                H_hk = H_hk.setup_single(simplify(symbolic_polynomial*H_hr.HcoeL(i,j,n)*kpoints_phase),i,j);
                            end
                        end
                    end
                else
                    % waiting
                    H_htrig = HR2Htrig(H_hr,'sym',options.sym);
                    H_hk = H_htrig.Htrig2HK(kpoints_frac,'sym',options.sym,'Order',options.Order);
                end
            end
            H_hk.HcoeL = simplify(H_hk.HcoeL);
        end
    end
    %% Eigen
    methods
        varargout = EIGENCAR_gen(H_hr,options)
        EIGENCARout = EIGENCAR_gen_sparse(H_hr,fermi,norb_enforce,klist_s_tmp)
    end
    %% symmetry
    methods
        function H_hr = rewrite(H_hr,options)
            arguments
                H_hr HR;
                options.rewind = false;
                options.Accuracy = 1e-6;
                options.type = '';
            end
            WANNUM = H_hr.WAN_NUM;
            if ~strcmp(H_hr.Type ,'list')
                
                if  H_hr.num
                    if isvector(H_hr.HnumL)
                        warning('May not need to rewrite. Do nothing');
                        return;
                    end
                    NRPTS_ = numel(H_hr.HnumL);
                    %HnumLtmp = zeros(1,NRPTS_);
                    sizeHcoeL = size(H_hr.HnumL);
                    HnumLtmp = reshape(H_hr.HnumL,[NRPTS_,1]);
                    % vector program
                    [iL,jL,kL]= ind2sub(sizeHcoeL,1:NRPTS_);
                    vectorList = [H_hr.vectorL(kL,:),iL.',jL.'];
                    H_hr.HnumL = HnumLtmp;
                    H_hr.vectorL = vectorList;
                    if H_hr.overlap
                        NRPTS_S = numel(H_hr.SnumL);
                        sizeScoeL = size(H_hr.SnumL);
                        SnumLtmp = reshape(H_hr.SnumL,[NRPTS_S,1]);
                        % vector program
                        [iL,jL,kL]= ind2sub(sizeScoeL,1:NRPTS_S);
                        vectorList_overlap = [H_hr.vectorL_overlap(kL,:),iL.',jL.'];
                        H_hr.SnumL = SnumLtmp;
                        H_hr.vectorL_overlap = vectorList_overlap;
                    end
                elseif H_hr.coe && ~H_hr.num
                    if isvector(H_hr.HcoeL)
                        warning('May not need to rewrite. Do nothing');
                        return;
                    end
                    NRPTS_ = numel(H_hr.HcoeL);
                    %HnumLtmp = zeros(1,NRPTS_);
                    sizeHcoeL = size(H_hr.HcoeL);
                    HcoeLtmp = reshape(H_hr.HcoeL,[NRPTS_,1]);
                    % vector program
                    [iL,jL,kL]= ind2sub(sizeHcoeL,1:NRPTS_);
                    vectorList = [H_hr.vectorL(kL,:),iL.',jL.'];
                    H_hr.HcoeL = HcoeLtmp;
                    %H_hr.HnumL = zeros(size(HcoeLtmp));
                    H_hr.vectorL = vectorList;
                    if H_hr.overlap
                        NRPTS_S = numel(H_hr.ScoeL);
                        sizeScoeL = size(H_hr.ScoeL);
                        ScoeLtmp = reshape(H_hr.ScoeL,[NRPTS_S,1]);
                        % vector program
                        [iL,jL,kL]= ind2sub(sizeScoeL,1:NRPTS_S);
                        vectorList_overlap = [H_hr.vectorL_overlap(kL,:),iL.',jL.'];
                        H_hr.ScoeL = ScoeLtmp;
                        H_hr.vectorL_overlap = vectorList_overlap;
                    end
                else
                    H_hr.HnumL = [];
                    H_hr.HcoeL = sym([]);
                    H_hr.vectorL = [];
                end
                H_hr.Type = 'list';
                H_hr = H_hr.simplify(options.Accuracy);
            elseif  strcmp(H_hr.Type ,'list') && options.rewind
                %vectorList = int32([0,0,0]);
                if H_hr.overlap
                    [vectorList_overlap,~,ic_S] = unique(H_hr.vectorL_overlap(:,1:3),'rows');
                    NRPTS_S= size(vectorList_overlap,1);
                    SnumLtmp = zeros(WANNUM,WANNUM,NRPTS_S);
                    ScoeLtmp = sym(zeros(WANNUM,WANNUM,NRPTS_S));
                end
                if H_hr.num
                    [vectorList,~,icL] = unique(H_hr.vectorL(:,1:3),'rows');
                    NRPTS_= size(vectorList,1);
                    HnumLtmp = zeros(WANNUM,WANNUM,NRPTS_);
                    sizemesh = [WANNUM,WANNUM,NRPTS_];
                    if H_hr.overlap
                        for n = 1:size(H_hr.vectorL_overlap,1)
                            SnumLtmp(H_hr.vectorL_overlap(n,4),H_hr.vectorL_overlap(n,5),ic_S(n)) = H_hr.SnumL(n);
                        end
                    end
                else
                    [vectorList,~,icL] = unique(H_hr.vectorL(:,1:3),'rows');
                    NRPTS_= size(vectorList,1);
                    HcoeLtmp = sym(zeros(WANNUM,WANNUM,NRPTS_));
                    sizemesh = [WANNUM,WANNUM,NRPTS_];
                    if H_hr.overlap
                        for n = 1:NRPTS_S
                            ScoeLtmp(H_hr.vectorL_overlap(n,4),H_hr.vectorL_overlap(n,5),ic_S(n)) = H_hr.ScoeL(n);
                        end
                    end
                end
                if H_hr.num
                    iL = double(H_hr.vectorL(:,4));
                    jL = double(H_hr.vectorL(:,5));
                    indL = sub2ind(sizemesh,iL,jL,icL);
                    HnumLtmp(indL) = H_hr.HnumL;
                    H_hr.HnumL = HnumLtmp;
                    H_hr.vectorL = vectorList;
                    if H_hr.overlap
                        H_hr.SnumL = SnumLtmp;
                        H_hr.vectorL_overlap = int32(vectorList_overlap);
                    end
                end
                if H_hr.coe
                    iL = double(H_hr.vectorL(:,4));
                    jL = double(H_hr.vectorL(:,5));
                    indL = sub2ind(sizemesh,iL,jL,icL);
                    HcoeLtmp(indL) = H_hr.HcoeL;
                    H_hr.HcoeL = HcoeLtmp;
                    H_hr.vectorL = vectorList;
                    if H_hr.overlap
                        H_hr.ScoeL = ScoeLtmp;
                        H_hr.vectorL_overlap = int32(vectorList_overlap);
                    end
                end
                H_hr.Type = 'mat';
            else
                
            end
            
        end
        function H_hr = rewind(H_hr)
            H_hr = H_hr.rewrite('rewind',true);
        end
        function H_hr = init(H_hr,options)
            arguments
                H_hr HR;
                options.level_cut {mustBeInteger} = 1;
                options.onsite logical = false;
                options.per_dir double = [1,1,1];
                options.chiral logical = false;
                options.spin logical = false;
                options.method  = 'nn_sparse';
                options.rough logical= false;
                options.hermitize = true;
                options.vectorL = [];
                options.fast logical= false;
            end
            % rough init also
            if ~isempty(options.vectorL)
                switch size(options.vectorL,2)
                    case 3
                        
                    case 5
                        
                end
                if options.hermitize
                    H_hr = H_hr.hermitize();
                end
                return;
            end
            % check nn_store
            if isempty(H_hr.nn_store)
                if isempty(H_hr.nn_store_smart)
                    error('You have not run the function: nn/nn_smart');
                else
                    if strcmp(options.method,'nn_sparse')
                        error('You should use the key value: ''method'',''nn_smart'' as input.');
                    end
                end
            end
            %
            %             if ~options.rough
            %                 if isempty(H_hr.elementL)
            %                     detailed(1) = false;
            %                 else
            %                     detailed(1) = true;
            %                 end
            %                 if isempty(H_hr.quantumL)
            %                     detailed(2) = false;
            %                 else
            %                     detailed(2) = true;
            %                 end
            %             else
            %                 detailed(1) = false;
            %                 detailed(2) = false;
            %             end
            % init
            N_orbit = H_hr.WAN_NUM;
            N_tot = H_hr.WAN_NUM;
            level_cut = options.level_cut;
            if ~strcmp(H_hr.Type,'list')
                H_hr = H_hr.rewrite();
            end
            % attention
            H_hr.coe = true;
            H_hr.num = false;
            if options.fast
                H_hr.vectorhopping = true;
                select_nn_store = H_hr.nn_store(H_hr.nn_store(:,10)<=options.level_cut,:);
                % try to speed up
                if options.chiral
                    element1L =  H_hr.elementL(select_nn_store(:,1));
                    element2L =  H_hr.elementL(select_nn_store(:,2));
                    select_nn_store(element1L == element2L,:) = [];
                end
                if options.spin
                    spin1L = H_hr.quantumL(select_nn_store(:,1),4);
                    spin2L = H_hr.quantumL(select_nn_store(:,2),4);
                    select_nn_store(spin1L == spin2L,:) = [];
                end
                if options.onsite
                    onsite_nnL = [(1:N_orbit)',(1:N_orbit)',zeros(N_orbit,8)];
                    select_nn_store = [select_nn_store;onsite_nnL];
                end
                nselect_nn_store = size(select_nn_store,1);
                H_hr.vectorL = select_nn_store(:,[6,7,8,1,2]);
                H_hr.AvectorL = eye(nselect_nn_store);
                H_hr.BvectorL = eye(nselect_nn_store);
                H_hr.CvectorL = eye(nselect_nn_store*2);
            else
                % list mode enforced
                if options.onsite
                    if strcmp(options.method,'nn_smart')
                        for i=1:N_orbit
                            A = sym(['A','_',num2str(i),'_',num2str(i),'_0_ubar'],'real');
                            B = sym(['B','_',num2str(i),'_',num2str(i),'_0_ubar'],'real');
                            H_hr = H_hr.set_hop_single( A+1i*B,i,i,[0 0 0],'sym');
                        end
                    else
                        onsite_nnL = [(1:N_orbit)',(1:N_orbit)',zeros(N_orbit,8)];
                        H_hr.nn_store = [H_hr.nn_store;onsite_nnL];
                    end
                end
                %   on-site
                %   hop
                if strcmp(options.method,'nn_smart')
                    pb = vasplib_tool_outer.CmdLineProgressBar('Setting ');
                    for i=1:N_orbit
                        %fprintf("setting (%4d/%4d) th orbital ... \n",i,N_orbit);
                        for j=1:N_tot
                            set_or_not = true;
                            if options.chiral
                                element1 =  H_hr.elementL(i);
                                element2 =  H_hr.elementL(j);
                                if element1 == element2
                                    set_or_not = false;
                                end
                            end
                            if options.spin
                                spin1 = H_hr.quantumL(i,4);
                                spin2 = H_hr.quantumL(j,4);
                                if spin1 == spin2
                                    set_or_not = false;
                                end
                            end
                            if set_or_not
                                nn = H_hr.nn_store_smart(i,j).nn;
                                for k = 1:size(nn,1)
                                    if nn(k).nn_level <= level_cut
                                        vector = nn(k).R_vector.*per_dir;
                                        A = HR.SymbolicVarible("A",vector,[i,j],nn(k).nn_level);
                                        B = HR.SymbolicVarible("B",vector,[i,j],nn(k).nn_level);
                                        SymHopping = A+1i*B;
                                        H_hr = H_hr.set_hop(SymHopping,i,j,vector,'sym');
                                    end
                                end
                            end
                            pb.print([i,j],[N_orbit,N_tot],' th orbital ...');
                        end
                    end
                    pb.delete();
                else
                    select_nn_store = H_hr.nn_store(H_hr.nn_store(:,10)<=options.level_cut,:);
                    pb = vasplib_tool_outer.CmdLineProgressBar('Setting ');
                    % try to speed up
                    if options.chiral
                        element1L =  H_hr.elementL(select_nn_store(:,1));
                        element2L =  H_hr.elementL(select_nn_store(:,2));
                        select_nn_store(element1L == element2L,:) = [];
                    end
                    if options.spin
                        spin1L = H_hr.quantumL(select_nn_store(:,1),4);
                        spin2L = H_hr.quantumL(select_nn_store(:,1),4);
                        select_nn_store(spin1L == spin2L,:) = [];
                    end
                    %
                    nselect_nn_store = size(select_nn_store,1);
                    iL = select_nn_store(:,1);STRiL = string(iL);
                    jL = select_nn_store(:,2);STRjL = string(jL);
                    R1L = select_nn_store(:,6);minusR1L = R1L <0;STRR1L = string(abs(R1L));STRR1L(minusR1L) =STRR1L(minusR1L)+"_bar";
                    R2L = select_nn_store(:,7);minusR2L = R2L <0;STRR2L = string(abs(R2L));STRR2L(minusR2L) =STRR2L(minusR2L)+"_bar";
                    R3L = select_nn_store(:,8);minusR3L = R3L <0;STRR3L = string(abs(R3L));STRR3L(minusR3L) =STRR3L(minusR3L)+"_bar";
                    nn_levelL = select_nn_store(:,10);STRnnL = string(nn_levelL)+"_ubar";
                    SuperscriptL = repmat("__",[nselect_nn_store 1]);
                    SubscriptL = repmat("_",[nselect_nn_store 1]);
                    affix_L = SuperscriptL+STRR1L+SuperscriptL+STRR2L+SuperscriptL+STRR3L+...
                        SubscriptL+STRiL+SubscriptL+STRjL+SubscriptL+STRnnL;
                    AL = "A" + affix_L;
                    BL = "B" + affix_L;
                    AsymL = str2sym(AL);
                    BsymL = str2sym(BL);
                    assume(AsymL,'real');
                    assume(BsymL,'real');
                    H_hr.vectorL = select_nn_store(:,[6,7,8,1,2]);
                    H_hr.HcoeL = AsymL+1i*BsymL;
                    %                     for in = 1:nselect_nn_store
                    %                         % there is a notation bug in List mode i(the first
                    %                         % one) means in home cell, j(the second one) is the
                    %                         % other. however the nn_sparse seems give a
                    %                         % different result. cause nn_sparse is also an
                    %                         % important function, we change here temperally
                    %                         % we fix the bug
                    %                         i = select_nn_store(in,1);
                    %                         j = select_nn_store(in,2);
                    %                         set_or_not = true;
                    %                         %                         if options.chiral
                    %                         %                             element1 =  H_hr.elementL(i);
                    %                         %                             element2 =  H_hr.elementL(j);
                    %                         %                             if element1 == element2
                    %                         %                                 set_or_not = false;
                    %                         %                             end
                    %                         %                         end
                    %                         %                         if options.spin
                    %                         %                             spin1 = H_hr.quantumL(i,4);
                    %                         %                             spin2 = H_hr.quantumL(j,4);
                    %                         %                             if spin1 == spin2
                    %                         %                                 set_or_not = false;
                    %                         %                             end
                    %                         %                         end
                    %                         if set_or_not
                    %                             nn_level = select_nn_store(in,10);
                    %                             vector = select_nn_store(in,[6,7,8]);
                    %                             A = HR.SymbolicVarible("A",vector,[i,j],nn_level);
                    %                             B = HR.SymbolicVarible("B",vector,[i,j],nn_level);
                    %                             SymHopping = A+1i*B;
                    %                             H_hr = H_hr.set_hop(SymHopping,i,j,vector,'sym');
                    %                         end
                    %                         pb.print(in,nselect_nn_store,' th orbital ...');
                    %                     end
                    pb.delete();
                end
            end
            if options.hermitize
                H_hr = H_hr.hermitize();
            end
        end
        function H_hr = applyOper(H_hr,SymOper,options)
            arguments
                H_hr HR;
                SymOper Oper = Oper();
                options.generator = false;
                options.HcoeL = false;
                options.Accuracy = 1e-6;
                options.fast = false;
                options.center = [0,0,0];
                options.Ugen = false;
            end
            options2 =options;
            nSymOper = length(SymOper);
            if options.Ugen
                try
                    BasisFunction = BasisFunc(H_hr);
                    % refresh SymOper
                    SymOper = SymOper.Ugen(BasisFunction,'Rm',H_hr.Rm,'center',options.center);
                catch
                    
                end
            end
            % use orth vector space
            if options.fast
                if isempty(H_hr.AvectorL) && isempty(H_hr.BvectorL)
                    H_hr = H_hr.init('fast',true);
                    H_hr = H_hr.hermitize();
                end
                if options.generator
                    options2.generator = false;
                    optionsCell = namedargs2cell(options2);
                    for i = 1:nSymOper
                        fprintf('******** apply (%d/%d)symmetry ********\n',i,nSymOper);
                        disp(SymOper(i));
                        SymOper_tmp = SymOper(i).generate_group();
                        H_hr =  applyOper(H_hr,SymOper_tmp,optionsCell{:});
                        fprintf('----------   SymVarNum: %d   ----------\n',...
                            rank(H_hr.CvectorL));
                    end
                    H_hr = H_hr.hermitize;
                    H_hr = H_hr.simplify(options.Accuracy);
                else
                    H_hr_R = H_hr;
                    nSymOper = length(SymOper);
                    pb = vasplib_tool_outer.CmdLineProgressBar('Applying Symmetry ...');
                    for j = 1:nSymOper
                        [H_hr_R(j),H_hr] = applyRU(H_hr,SymOper(j));
                        pb.print(j,nSymOper);
                    end
                    pb.delete();
                    H_hr = sum(H_hr_R);
                    H_hr = H_hr.hermitize;
                    H_hr = H_hr.simplify(options.Accuracy);
                end
                return;
            end
            % use numerical of symbolic
            
            if ~H_hr.coe && ~H_hr.num
                H_hr = H_hr.init();
                H_hr = H_hr.hermitize();
            end
            if ~strcmp(H_hr.Type , 'list')
                H_hr = H_hr.rewrite();
            end
            if length(SymOper) == 1
                if ~SymOper.conjugate && ~SymOper.antisymmetry && isequal(SymOper.R,eye(3))
                    return;
                end
                nSymOper = length(SymOper);
                fprintf('******** apply (%d/%d)symmetry ********\n',1,nSymOper);
                disp(SymOper);
                if options.generator
                    SymOper_tmp = SymOper.generate_group();
                    nSymOper_tmp = length(SymOper_tmp);
                    pb = vasplib_tool_outer.CmdLineProgressBar('Applying Symmetry ...');
                    %   for k = 1 : 10
                    %       pb.print(k,10)
                    %       % do stuff
                    %   end
                    H_hr_R = H_hr;
                    %H_hr_R = repmat(H_hr_R,[]);
                    for j = 1:nSymOper_tmp
                        pb.print(j,nSymOper_tmp);
                        [H_hr_R(j),H_hr] = applyRU(H_hr,SymOper_tmp(j));
                    end
                    pb.delete();
                    H_hr = sum(H_hr_R)/nSymOper_tmp;
                    H_hr = H_hr.simplify();
                elseif options.HcoeL
                    % when apply U it will reseq    the HcoeL
                    [H_hr_R,H_hr] = applyRU(H_hr,SymOper);
                    if isequal(H_hr_R.HcoeL,H_hr.HcoeL)
                    else
                        Equationlist_r = (real(H_hr.HcoeL - H_hr_R.HcoeL) == 0);
                        Equationlist_i = (imag(H_hr.HcoeL - H_hr_R.HcoeL) == 0);
                        %Equationlist_r = Htrig.isolateAll(Equationlist_r,real(H_hr.HcoeL));
                        %Equationlist_i = Htrig.isolateAll(Equationlist_i,imag(H_hr.HcoeL));
                        Equationlist_r = HR.isolateAll(Equationlist_r);
                        Equationlist_i = HR.isolateAll(Equationlist_i);
                        HcoeLtmp = H_hr.HcoeL ;
                        HcoeLtmp_r = subs(real(HcoeLtmp),lhs(Equationlist_r),rhs(Equationlist_r));
                        HcoeLtmp_i = subs(imag(HcoeLtmp),lhs(Equationlist_i),rhs(Equationlist_i));
                        H_hr.HcoeL = HcoeLtmp_r + 1i*HcoeLtmp_i;
                    end
                    H_hr = H_hr.simplify();
                end
            else
                nSymOper = length(SymOper);
                for i = 1:nSymOper
                    fprintf('******** apply (%d/%d)symmetry ********\n',i,nSymOper);
                    disp(SymOper(i));
                    if options.generator
                        SymOper_tmp = SymOper(i).generate_group();
                        nSymOper_tmp = length(SymOper_tmp);
                        pb = vasplib_tool_outer.CmdLineProgressBar('Applying Symmetry ...');
                        %   for k = 1 : 10
                        %       pb.print(k,10)
                        %       % do stuff
                        %   end
                        H_hr_R = H_hr;
                        %H_hr_R = repmat(H_hr_R,[]);
                        for j = 1:nSymOper_tmp
                            pb.print(j,nSymOper_tmp);
                            [H_hr_R(j),H_hr] = applyRU(H_hr,SymOper_tmp(j));
                        end
                        pb.delete();
                        H_hr = sum(H_hr_R)/nSymOper_tmp;
                        H_hr = H_hr.simplify(options.Accuracy);
                    else
                        %fprintf('    ');
                        H_hr = H_hr.applyOper(SymOper(i),'generator','false');
                    end
                    fprintf('----------   SymVarNum: %d   ----------\n',length(H_hr.symvar_list));
                end
            end
        end
        function H_hr = symmetrize(H_hr,SymOper,options)
            arguments
                H_hr HR;
                SymOper Oper = Oper();
                options.generator = false;
            end
            
        end
        function H_hr = dualize(H_hr)
            NRPTS_ = H_hr.NRPTS;
            vectorList = H_hr.vectorL;
            vectorList_oppo(:,1:3) = -vectorList(:,1:3);
            if size(vectorList,2) == 5
                vectorList_oppo(:,4) = vectorList(:,5);
                vectorList_oppo(:,5) = vectorList(:,4);
            end
            for i = 1:NRPTS_
                vector_tmp_oppo = vectorList_oppo(i,:);
                [~,j]=ismember(vector_tmp_oppo,H_hr.vectorL,'rows');
                if j == 0
                    %H_hr.vectorL =  [H_hr.vectorL;vector_tmp_oppo];
                    H_hr = H_hr.add_empty_one(vector_tmp_oppo);
                    j = H_hr.NRPTS;
                    H_hr.Duality_vector_dist(j) = i ;
                end
                H_hr.Duality_vector_dist(i) = j  ;
            end
        end
        function [H_hr,R_vector_dist_] = dualizeR(H_hr,Rf)
            NRPTS_ = H_hr.NRPTS;
            ONESLIST = ones(NRPTS_,1);
            H_hr = H_hr.timtj_gen();
            timtj_mat = H_hr.timtj{2};
            Size_timtj_mat = size(timtj_mat);
            vectorList = double(H_hr.vectorL(:,1:3));
            IndList1 = sub2ind(Size_timtj_mat,H_hr.vectorL(:,4),H_hr.vectorL(:,5),ONESLIST);
            IndList2 = sub2ind(Size_timtj_mat,H_hr.vectorL(:,4),H_hr.vectorL(:,5),ONESLIST*2);
            IndList3 = sub2ind(Size_timtj_mat,H_hr.vectorL(:,4),H_hr.vectorL(:,5),ONESLIST*3);
            vectorL_addtional = [timtj_mat(IndList1),timtj_mat(IndList2),timtj_mat(IndList3)];
            vectorL_  = floor(H_hr.orbL(H_hr.vectorL(:,4),:)+(vectorList-vectorL_addtional)*double(Rf));
            for i = 1:NRPTS_
                vector_tmp_oppo = vectorL_(i,:);
                [~,j]=ismember(vector_tmp_oppo,vectorList,'rows');
                if j == 0
                    %H_hr.vectorL =  [H_hr.vectorL;vector_tmp_oppo];
                    H_hr = H_hr.add_empty_one([vector_tmp_oppo,H_hr.vectorL(i,14:5)]);
                    j = H_hr.NRPTS;
                    H_hr.R_vector_dist(j) = i ;
                end
                H_hr.R_vector_dist(i) = j  ;
                R_vector_dist_ = H_hr.R_vector_dist;
            end
        end
        function [ml_cell,ij_list] = Tij2lm(H_hr,Rf)
            %WAN_NUM = H_hr.WAN_NUM;
            NRPTS_ = H_hr.NRPTS;
            orblist = H_hr.orbL;
            Norb = size(orblist,1);
            ij_list= unique(double(H_hr.vectorL(:,4:5)),'rows');
            addtional_ij_list = [];
            count  = 0;
            ml_cell{size(ij_list,1)} = [];
            for n = 1:size(ij_list,1)
                t_i = orblist(ij_list(n,1),:);
                t_j = orblist(ij_list(n,2),:);
                % in PRM, they use Tml =Sg(tm −tl)−tj −ti
                % t_i_p_t_j = t_j + t_i;
                % maybe a typo
                t_j_m_t_i = t_j - t_i;
                ml_cell{n} = [];
                tmp_count = 0;
                % go over lm
                for m = 1:Norb
                    t_m = orblist(m,:);
                    for l = 1:Norb
                        t_l = orblist(l,:);
                        % test if a lattice vector
                        T_ij__ml = (t_m - t_l)*Rf - t_j_m_t_i;
                        if vasplib.LatticeVectorTest(T_ij__ml)
                            % addtional i j ?
                            [~,j]=ismember([l,m],ij_list,'rows');
                            if j == 0
                                count = count+1;
                                addtional_ij_list(count,:) = [l,m];
                            end
                            tmp_count = tmp_count +1;
                            ml_cell{n}(tmp_count,:) = [m,l,T_ij__ml];
                        end
                    end
                end
            end
            %
            while ~isempty(addtional_ij_list)
                ij_list = [ij_list;addtional_ij_list];
                addtional_ij_list = [];
                count  = 0;
                for n = n+1:size(ij_list,1)
                    t_i = orblist(ij_list(n,1),:);
                    t_j = orblist(ij_list(n,2),:);
                    t_j_m_t_i = t_j + t_i;
                    ml_cell{n} = [];
                    tmp_count = 0;
                    % go over lm
                    for l = 1:Norb
                        t_l = orblist(l,:);
                        for m = 1:Norb
                            t_m = orblist(m,:);
                            % test if a lattice vector
                            T_ij__ml =  (t_m - t_l)*Rf - t_j_m_t_i;
                            if vasplib.LatticeVectorTest(T_ij__ml)
                                % addtional i j ?
                                [~,j]=ismember([l,m],ij_list,'rows');
                                if j == 0
                                    count = count+1;
                                    addtional_ij_list(count,:) = [l,m];
                                end
                                tmp_count = tmp_count +1;
                                ml_cell{n}(tmp_count,:) = [l,m];
                            end
                        end
                    end
                end
            end
        end
        function [H_hr,VectorDistMat] = dualizeOper(H_hr,SymOper)
            % init
            Accuracy = 1e-6;
            NRPTS_ = H_hr.NRPTS;
            Rf = double(Oper.Rc2Rf(inv(SymOper.R),H_hr.Rm));
            U = SymOper.U;
            %ONESLIST = ones(NRPTS_,1);
            invRf = inv(double(Rf));
            invU = inv(U);
            [ml_cell,ij_list] = Tij2lm(H_hr,Rf);
            vectorList = double(H_hr.vectorL(:,1:3));
            ij_list_full = double(H_hr.vectorL(:,4:5));
            %
            VectorDistMat = zeros(NRPTS_,NRPTS_);
            % Rotated vector
            % https://doi.org/10.1103/PhysRevMaterials.2.103805
            % 𝐻̃^𝑖𝑗[𝐑']=1/(|𝐺|) ∑_(((𝑔∈𝐺)¦(𝑙,𝑚)))  𝐷_𝑖𝑙 (𝑔)𝐻^𝑙𝑚 [𝑆_𝑔^(−1) (𝐑'−𝐓_𝑖𝑗^𝑚𝑙 )] 𝐷_𝑚𝑗 (𝑔^(−1) )
            for in = 1:NRPTS_
                ij = ij_list_full(in,:);
                Rvector = vectorList(in,:);
                [~,k] = ismember(ij,ij_list,'rows');
                ml_list = ml_cell{k};
                i = ij(1);
                j = ij(2);
                for kn  = 1 : size(ml_list,1)
                    T_ij__ml = ml_list(kn,3:5);
                    m =  ml_list(kn,1);
                    l =  ml_list(kn,2);
                    vector_tmp_oppo = (Rvector - T_ij__ml)*invRf;
                    tmp_vector = int32([vector_tmp_oppo,[l,m]]);
                    [~,jn] = ismember(tmp_vector,H_hr.vectorL,'rows');
                    %[~,jn]=ismember(vector_tmp_oppo,vectorList,'rows');
                    if jn == 0
                        %H_hr.vectorL =  [H_hr.vectorL;vector_tmp_oppo];
                        H_hr = H_hr.add_empty_one(tmp_vector);
                        jn = H_hr.NRPTS;
                        U_tmp = U(j,m)*invU(l,i);
                        U_tmp_r = real(U_tmp);
                        U_tmp_i = imag(U_tmp);
                        if abs(U_tmp_r) < Accuracy
                            U_tmp_r = 0 ;
                        end
                        if abs(U_tmp_i) < Accuracy
                            U_tmp_i = 0 ;
                        end
                        VectorDistMat(jn,in) = U_tmp_r + 1i * U_tmp_i ;
                    end
                    U_tmp = U(i,l)*invU(m,j);
                    U_tmp_r = real(U_tmp);
                    U_tmp_i = imag(U_tmp);
                    if abs(double(U_tmp_r)) < Accuracy
                        U_tmp_r = 0 ;
                    end
                    if abs(double(U_tmp_i)) < Accuracy
                        U_tmp_i = 0 ;
                    end
                    VectorDistMat(in,jn) = U_tmp_r + 1i * U_tmp_i ;
                end
            end
            
        end
        function H_hr = hermitize(H_hr)
            H_hr_bk = H_hr';
            if H_hr.vectorhopping
                H_hr = H_hr+H_hr_bk;
                H_hr = H_hr.simplify;
            else
                
                if H_hr.coe
                    Equationlist_r = real(H_hr.HcoeL - H_hr_bk.HcoeL) == 0;
                    Equationlist_i = imag(H_hr.HcoeL - H_hr_bk.HcoeL) == 0;
                    %Equationlist_r = Htrig.isolateAll(Equationlist_r,real(H_hr.HcoeL));
                    %Equationlist_i = Htrig.isolateAll(Equationlist_i,imag(H_hr.HcoeL));
                    Equationlist_r = HR.isolateAll(Equationlist_r);
                    Equationlist_i = HR.isolateAll(Equationlist_i);
                    HcoeLtmp = subs(H_hr.HcoeL,lhs(Equationlist_r),rhs(Equationlist_r));
                    HcoeLtmp = subs(HcoeLtmp,lhs(Equationlist_i),rhs(Equationlist_i));
                    H_hr.HcoeL = HcoeLtmp;
                end
                if H_hr.num
                    H_hr.HnumL = (H_hr_bk.HnumL + H_hr.HnumL )/2;
                end
            end
        end
        function H_hr_bk = subsOper(H_hr,SymOper)
            arguments
                H_hr Htrig;
                SymOper Oper = Oper();
            end
            %             if isequal(zeros(size(H_hr.HnumL)),H_hr.HnumL)
            %                 H_hr.num = false;
            %             else
            %                 H_hr.num = true;
            %             end
            if isequal(sym(zeros(size(H_hr.HcoeL))),H_hr.HcoeL)
                H_hr.coe = false;
            else
                H_hr.coe = true;
            end
            
            if ~H_hr.coe
                H_hr = H_hr.init();
                H_hr = H_hr.hermitize();
            end
            
            if length(SymOper) == 1
                if ~SymOper.conjugate && ~SymOper.antisymmetry && isequal(SymOper.R,eye(3))
                    return;
                end
                [H_hr_bk,H_hr]  = H_hr.applyR(inv(SymOper.R));
                if isnan(SymOper.U)
                    %build U
                end
                H_hr_bk  = H_hr_bk.applyU(SymOper.U,SymOper.conjugate,SymOper.antisymmetry);
            end
        end
        function [H_hr_R,H_hr] = applyRU(H_hr,SymOper )
            arguments
                H_hr HR;
                SymOper     ;
            end
            [H_hr,VectorDistMat] = dualizeOper(H_hr,SymOper);
            % square Mat test
            if size(VectorDistMat,1)~=size(VectorDistMat,2)
                error('check why?');
            end
            %
            H_hr_R = H_hr;
            if H_hr.vectorhopping
                AvectorLtmp = H_hr_R.AvectorL;
                BvectorLtmp = H_hr_R.BvectorL;
                CvectorLtmp = H_hr_R.CvectorL;
                if SymOper.conjugate
                    BvectorLtmp = -(BvectorLtmp);
                    %CvectorLtmp = conj(CvectorLtmp);
                    CvectorLtmp(end/2+1:end,:) = -CvectorLtmp(end/2+1:end,:);
                end
                if SymOper.antisymmetry
                    AvectorLtmp = -AvectorLtmp;
                    BvectorLtmp = -BvectorLtmp;
                    CvectorLtmp = -CvectorLtmp;
                end
                H_hr_R.AvectorL = VectorDistMat*AvectorLtmp;
                H_hr_R.BvectorL = VectorDistMat*BvectorLtmp;
                %VectorDistMat*CvectorLtmp;
                CL1 = VectorDistMat*CvectorLtmp(1:end/2,:);
                CL2 = VectorDistMat*CvectorLtmp(end/2+1:end,:);
                H_hr_R.CvectorL = [real(CL1)-imag(CL2);imag(CL1)+real(CL2)];
                return;
            end
            
            if H_hr.coe
                HcoeLtmp = H_hr_R.HcoeL;
                if SymOper.conjugate
                    %H_hr = H_hr.applyR(diag([-1,-1,-1]));
                    HcoeLtmp = conj(HcoeLtmp);
                    %HcoeLtmp = Htrig.matrixtimespage(H_hr.factorlist_parity(),HcoeLtmp);
                end
                if SymOper.antisymmetry
                    HcoeLtmp = -HcoeLtmp;
                end
                H_hr_R.HcoeL = VectorDistMat*HcoeLtmp;
            end
            if H_hr.num == true
                HnumLtmp = H_hr_R.HnumL;
                if SymOper.conjugate
                    HnumLtmp = conj(HnumLtmp);
                end
                if SymOper.antisymmetry
                    HnumLtmp = -HnumLtmp;
                end
                H_hr_R.HnumL = VectorDistMat*HnumLtmp;
            end
        end
        function [H_hr_R,H_hr] = applyR(H_hr,R)
            arguments
                H_hr HR;
                R ;
            end
            
            Rf = Oper.Rc2Rf(inv(R),H_hr.Rm);
            % dualize R
            H_hr = H_hr.rewrite();
            [H_hr,R_vector_dist_] = dualizeR(H_hr,Rf);
            H_hr_R = H_hr;
            %apply for vectorL
            %[H_hr,Smat] = H_hr.Smatgen(R);
            
            if H_hr.num
                H_hr_R.HnumL = H_hr_R.HnumL(R_vector_dist_);
            end
            if H_hr.coe
                H_hr_R.HcoeL = H_hr_R.HcoeL(R_vector_dist_);
            end
            %apply R
            
        end
        function [H_hr,Smat] = Smatgen(H_hr,R,Accuracy)
            arguments
                H_hr Htrig;
                R  ;
                Accuracy double = 6;
            end
            if isa(R,'sym')
                H_hr.coe = true;
            else
                H_hr.coe = false;
            end
            BASIS_NUM = H_hr.Basis_num;
            %HsymC_bk = H_hr.HsymL_trig;
            HsymC = H_hr.HsymL_trig;
            syms k_x k_y k_z real;
            % add in property
            %             varlist = [sin(k_x),sin(k_y),sin(k_z),...
            %                 cos(k_x),cos(k_y),cos(k_z)...
            %                 ];
            %varlist = H_hr.HsymL_trig_bk;
            
            if length(R) == 3
                k = R*[k_x; k_y; k_z];
                HsymC = subs(HsymC,[k_x k_y k_z],k.');
            elseif length(R) == 2
                k = R*[k_x; k_y];
                HsymC = subs(HsymC,[k_x k_y],k.');
            elseif length(R) == 1
                k = R*[k_x];
                HsymC = subs(HsymC,[k_x],k.');
            end
            HsymC = expand(simplify(HsymC));
            % find new db
            for k = 1:length(HsymC)
                [coeff_trig,~] = coeffs(HsymC(k),H_hr.HsymL_trig_bk);
                for i = 1:numel(coeff_trig)
                    tmp_label = contains(string(coeff_trig),H_hr.seeds);
                    if sum(tmp_label)
                        [~,coeff_trig_list] = coeffs(coeff_trig(i));
                        for j = 1:numel(coeff_trig_list)
                            tmp_label2 = contains(string(coeff_trig_list(j)),H_hr.seeds);
                            if sum(tmp_label2)
                                H_hr = H_hr.find_HsymL_trig_bk(coeff_trig_list(j));
                            end
                        end
                    end
                end
            end
            % find new basis & update
            for i = 1:numel(HsymC)
                [~,B] = coeffs(HsymC(i),H_hr.HsymL_trig_bk);
                for k = 1:length(B)
                    Kind = H_hr.k_symbol2Kind(B(k));
                    if isempty(Kind)
                        Kind = H_hr.Kinds+1;
                        H_hr.HsymL_trig(Kind) = B(k);
                        H_hr.HcoeL(:,:,Kind) = sym(zeros(BASIS_NUM,BASIS_NUM,1));
                        H_hr.HnumL(:,:,Kind)  = (zeros(BASIS_NUM,BASIS_NUM,1));
                    end
                end
            end
            % redo
            HsymC_bk = H_hr.HsymL_trig;
            %             HsymC = HsymC_bk;
            %             if length(R) == 3
            %                 k = R*[k_x; k_y; k_z];
            %                 HsymC = subs(HsymC,[k_x k_y k_z],k.');
            %             elseif length(R) == 2
            %                 k = R*[k_x; k_y];
            %                 HsymC = subs(HsymC,[k_x k_y],k.');
            %             elseif length(R) == 1
            %                 k = R*[k_x];
            %                 HsymC = subs(HsymC,[k_x],k.');
            %             end
            %             HsymC = expand(HsymC);
            Smat =sym(zeros(numel(HsymC_bk)));
            for i = 1:numel(HsymC)
                [A,B] = coeffs(HsymC(i),H_hr.HsymL_trig_bk);
                for k = 1:length(B)
                    tempSym = B(k);
                    for l  = 1:numel(HsymC_bk)
                        if isequal(tempSym,HsymC_bk(l))
                            Smat(l,i)=A(k);
                            break;
                        end
                    end
                end
            end
            if ~H_hr.coe
                Smat = roundn(double(Smat),-Accuracy);
            end
            %                 if i == 3
            %                     error('ss');
            %                 end
        end
        function Factorlist_parity = factorlist_parity(H_hr)
            syms k_x k_y k_z real;
            HsymC = H_hr.HsymL_trig;
            HsymC = subs(HsymC,[k_x k_y k_z],-[k_x k_y k_z]);
            Factorlist_parity = simplify(HsymC./H_hr.HsymL_trig);
        end
        function H_hr = nn(H_hr,search_range,Accuracy,Rlength_cut,options)
            arguments
                H_hr HR;
                search_range double = [0 0 0];
                Accuracy double = 1e-4;
                Rlength_cut double = 5;
                options.MAX_NN_LENGTH =  10000000;
                options.onsite = false;
            end
            optionsCell = namedargs2cell(options);
            H_hr = nn@vasplib(H_hr,search_range,Accuracy,Rlength_cut,optionsCell{:});
            %H_hr = H_hr.rewrite();
        end
    end
    %% script
    methods
        function [EIGENCAR_3D,klist1,klist2,WEIGHTCAR_3D,WAVECAR_3D] = EIGENCAR_gen_3D(H_hr,kmesh,k3d,options)
            arguments
                H_hr HR
                kmesh = [100,100];
                k3d = [-0.5, -0.5, 0;... %begin kpoints
                    1.0 , 0,   0;... %  first vector
                    0   , 1.0,   0;... %  second vector
                    ];
                options.output ='raw_data';
                options.LWAVE = false;
                options.cartisian = true;
                options.fin_dir = 3;
                options.ProjectionMethod = 'slab';
                options.ProjectionStruct = struct('discrimination',0.1,'center',[0.5,0.5,0.5],'orientation',2,'sign',false);
                options.WEIGHTCAR = false;
                options.norb = -1;
                options.fermi = 0;
            end
            options2 = rmfield(options,{'output','cartisian','fin_dir'});
            optionsCell = namedargs2cell(options2);
            [H_hr.klist_s,klist1,klist2] = H_hr.kmesh3D(kmesh,k3d,'fermiarc');
            if options.WEIGHTCAR
                [EIGENCAR_3D,WAVECAR_3D,WEIGHTCAR_3D] = H_hr.EIGENCAR_gen(...
                    optionsCell{:});
                %                     'LWAVE',options.LWAVE,'WEIGHTCAR',options.WEIGHTCAR,...
                %                     'ProjectionMethod',options.ProjectionMethod,'ProjectionStruct',options.ProjectionStruct);
            else
                [EIGENCAR_3D,WAVECAR_3D] = H_hr.EIGENCAR_gen('LWAVE',options.LWAVE);
                WEIGHTCAR_3D = [];
            end
            if options.LWAVE
                
            else
                WAVECAR_3D = [];
            end
            if options.cartisian
                options.output = 'refined';
                klist = H_hr.klist_s*H_hr.Gk;
            end
            if options.LWAVE
                WEIGHTCAR_3D = WAVECAR_3D;
            end
            if strcmp(options.output ,'refined')
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
                EIGENCAR_3D = reshape(EIGENCAR_3D.',kmesh(1),kmesh(2),[]);
                if options.WEIGHTCAR
                    WEIGHTCAR_3D = reshape(WEIGHTCAR_3D.',kmesh(1),kmesh(2),[]);
                end
            end
        end
        function [EIGENCAR_slab,klist_l,kpoints_l,kpoints_name] = slab(H_hr,repeatnum,fin_dir,KPOINTS_slab,norb_enforce,fermi)
            % -------------- nargin ------------------
            if nargin < 6
                fermi = 0;
            end
            if nargin < 5
                norb_enforce  = -1;
            end
            if nargin <4
                KPOINTS_slab = 'KPOINTS_slab';
            end
            if nargin < 3
                fin_dir     =  2;
            end
            if nargin < 2
                repeatnum   = 10;
            end
            glue_edges  = false;
            vacuum_mode = 1;
            % Gen Slab
            H_hr_slab = H_hr.cut_piece(repeatnum,fin_dir,glue_edges,vacuum_mode);
            % load KPOINTS
            H_hr_slab = H_hr_slab < KPOINTS_slab;
            % slab band
            EIGENCAR_slab = H_hr_slab.EIGENCAR_gen('fermi',fermi,'norb',norb_enforce);
            [klist_l,kpoints_l,kpoints_name] = H_hr_slab.kpath_information();
        end
        function [EIGENCAR_slab,WEIGHTCAR_slab,klist1,klist2] = slab_fermiarc(H_hr,repeatnum,fin_dir,KPOINTS_slab,norb_enforce,fermi)
            
        end
        function [EIGENCAR,orb_list,WEIGHTCAR,klist_l,kpoints_l,kpoints_name] = EIGENCAR_gen_wire(H_hr,Nslab,fermi,norb_enforce,KPOINTS_wire,vacuum_mode,np)
            %--------  init  --------
            import vasplib_tool.*
            % -------------- nargin ------------------
            if nargin < 2
                Nslab = [0 0 0];
            end
            if nargin < 3
                fermi = 0;
            end
            if nargin < 4
                norb_enforce  = -1;
            end
            if nargin < 5
                KPOINTS_wire = 'KPOINTS_wire';
            end
            if nargin < 6
                vacuum_mode = 1;
            end
            if nargin < 7
                np = 0;
            end
            
            if isequal(Nslab, [0,0,0])
                H_hr_wire = H_hr;
            else
                H_hr_wire = H_hr.Hnanowire_gen(Nslab,np,vacuum_mode);
            end
            H_hr_wire = H_hr_wire < KPOINTS_wire;
            H_hr_wire = H_hr_wire.sparse();
            %     [NRPTS,~]=size(vectorlist);
            %--------  swich  --------
            if np >1
                disp('parallel mode, we will use local settings, please set before.');
                np_handle = parpool('local',np);
            else
                %disp('without parallel mode.');
            end
            %--------  transform  --------
            % %             disp('The Nslab is ');
            % %             disp('we suppose the nanowire along a3 direction, please rotate it before by supercell_hr.');
            Hnum_list_wire = H_hr_wire.HnumL;
            vector_list_wire = double(H_hr_wire.vectorL);
            [nz,~] = size(vector_list_wire);
            factor_list_wire = exp(1i*2*pi*H_hr_wire.klist_s*vector_list_wire');
            %--------  HSVCAR  --------
            % generate orbitals of a finite model
            % -------------- nargin ------------------
            orb_list  = H_hr_wire.orbL;
            %HSVCAR_hinge = vasplib.HSVCAR_gen(orb_list,'hinge',0.05,[0.5,0.5,0.5],-3);
            HSVCAR_hinge = vasplib.HSVCAR_gen(orb_list,'hinge');
            NWAVE = H_hr_wire.WAN_NUM;
            %--------  init  --------
            if norb_enforce <0
                NBANDS=NWAVE;
            elseif norb_enforce >0
                NBANDS=norb_enforce;
            else
            end
            kn = size(H_hr_wire.klist_s,1);
            EIGENCAR = zeros(NBANDS,kn);
            WEIGHTCAR = zeros(NBANDS,kn);
            %--------  calculate  --------
            if np >1 && kn >1
                fprintf('begining parfor loop, we have %d kpoints\n',kn);
                parfor ki =1:kn
                    A = zeros(NWAVE,NWAVE);
                    U = zeros(NWAVE,NWAVE);
                    factor_list = factor_list_wire(ki,:);
                    Hout = sparse(NWAVE,NWAVE);
                    for iz = 1:nz
                        Hout = Hout+Hnum_list_wire{iz}*factor_list(iz);
                    end
                    Hout = (Hout+Hout')/2;
                    % U = zeros(NBANDS);
                    % A = zeros(NWAVE,NBANDS);
                    if norb_enforce <0
                        [A, U]=eig(full(Hout));
                    elseif norb_enforce >0
                        [A, U]=eigs(Hout,NBANDS,fermi);
                        [A, U]= park.sorteig(U,A);
                    else
                    end
                    EIGENCAR(:,ki) = diag(U);
                    [~,WEIGHTCAR(:,ki)] = COLORCAR_gen(A,HSVCAR_hinge);
                    fprintf('%d th kpoints has been calculated in %d kpoints total\n',ki,kn);
                end
                %delete(np_handle);
            elseif  kn <2
                WEIGHTCAR = zeros(NWAVE,NBANDS);
                for ki =1:kn
                    factor_list = factor_list_wire(ki,:);
                    Hout = sparse(NWAVE,NWAVE);
                    for iz = 1:nz
                        Hout = Hout+Hnum_list_wire{iz}*factor_list(iz);
                    end
                    Hout = (Hout+Hout')/2;
                    if norb_enforce <0
                        [A, U]=eig(full(Hout));
                    elseif norb_enforce >0
                        [A, U]=eigs(Hout,NBANDS,fermi);
                        [A, U]= park.sorteig(U,A);
                    else
                    end
                    EIGENCAR(:,ki) = diag(U);
                    if kn >1
                        [~,WEIGHTCAR(:,ki)] = COLORCAR_gen(A,HSVCAR_hinge);
                    else
                        WEIGHTCAR= A;
                    end
                    fprintf('%d th kpoint(%7.4f, %7.4f, %7.4f) has been calculated in %d kpoints total\n',...
                        ki,H_hr_wire.klist_s(ki,1),H_hr_wire.klist_s(ki,2),H_hr_wire.klist_s(ki,3),kn);
                end
            else
                for ki =1:kn
                    factor_list = factor_list_wire(ki,:);
                    Hout = sparse(NWAVE,NWAVE);
                    for iz = 1:nz
                        Hout = Hout+Hnum_list_wire{iz}*factor_list(iz);
                    end
                    Hout = (Hout+Hout')/2;
                    if norb_enforce <0
                        [A, U]=eig(full(Hout));
                    elseif norb_enforce >0
                        [A, U]=eigs(Hout,NBANDS,fermi);
                        [A, U]= park.sorteig(U,A);
                    else
                    end
                    EIGENCAR(:,ki) = diag(U);
                    [~,WEIGHTCAR(:,ki)] = COLORCAR_gen(A,HSVCAR_hinge);
                    fprintf('%d th kpoint(%7.4f, %7.4f, %7.4f) has been calculated in %d kpoints total\n',...
                        ki,H_hr_wire.klist_s(ki,1),H_hr_wire.klist_s(ki,2),H_hr_wire.klist_s(ki,3),kn);
                end
            end
            
            if np >1
                delete(np_handle);
            end
            [klist_l,kpoints_l,kpoints_name] = H_hr_wire.kpath_information();
        end
        function [EIGENCAR,orb_list,WAVECAR] = EIGENCAR_gen_disk(H_hr,Nslab,fermi,norb_enforce,kpoints,vacuum_mode,np)
            %--------  init  --------
            import vasplib_tool.*
            % -------------- nargin ------------------
            if nargin < 2
                Nslab = [10 10 0];
            end
            if nargin < 3
                fermi = 0;
            end
            if nargin < 4
                norb_enforce  = -1;
            end
            if nargin < 5
                kpoints = [0 0 0];
            end
            if nargin < 6
                vacuum_mode = 1;
            end
            if nargin < 7
                np = 0;
            end
            if isequal(Nslab, [0,0,0])
                H_hr_wire = H_hr;
            else
                H_hr_wire = H_hr.Hnanowire_gen(Nslab,np,vacuum_mode);
            end
            % H_hr_wire = H_hr.Hnanowire_gen(Nslab,np,vacuum_mode);
            % H_hr_wire = H_hr_wire < KPOINTS_wire;
            H_hr_wire = H_hr_wire.sparse();
            %     [NRPTS,~]=size(vectorlist);
            %--------  transform  --------
            % %             disp('The Nslab is ');
            % %             disp('we suppose the nanowire along a3 direction, please rotate it before by supercell_hr.');
            Hnum_list_wire = H_hr_wire.HnumL;
            vector_list_wire = double(H_hr_wire.vectorL);
            [nz,~] = size(vector_list_wire);
            factor_list_wire = exp(1i*2*pi*kpoints*vector_list_wire');
            %--------  HSVCAR  --------
            % generate orbitals of a finite model
            % -------------- nargin ------------------
            orb_list  = H_hr_wire.orbL;
            NWAVE = H_hr_wire.WAN_NUM;
            %--------  init  --------
            if norb_enforce <0
                NBANDS=NWAVE;
            elseif norb_enforce >0
                NBANDS=norb_enforce;
            else
            end
            kn = size(kpoints,1);
            EIGENCAR = zeros(NBANDS,kn);
            %--------  calculate  --------
            
            WAVECAR = zeros(NWAVE,NBANDS);
            for ki =1:kn
                factor_list = factor_list_wire(ki,:);
                Hout = sparse(NWAVE,NWAVE);
                for iz = 1:nz
                    Hout = Hout+Hnum_list_wire{iz}*factor_list(iz);
                end
                Hout = (Hout+Hout')/2;
                if norb_enforce <0
                    [A, U]=eig(full(Hout));
                elseif norb_enforce >0
                    [A, U]=eigs(Hout,NBANDS,fermi);
                    [A, U]= HR.sorteig(U,A);
                else
                end
                EIGENCAR(:,ki) = diag(U);
                WAVECAR= A;
                fprintf('%d th kpoint(%7.4f, %7.4f, %7.4f) has been calculated in %d kpoints total\n',...
                    ki,kpoints(1),kpoints(2),kpoints(3),kn);
            end
            
            if np >1
                delete(np_handle);
            end
            
            
        end
        function [DOSCAR_l,DOSCAR_b,DOSCAR_r,w_list,klist_l,kpoints_l,kpoints_name] = surf(H_hr,w_range,fin_dir,KPOINTS_surf,principle_layer,eta,fermi,mode)
            %nargin
            if nargin < 8
                mode = 'Green_iter';
            end
            if nargin < 7
                fermi = 0;
            end
            if nargin < 6
                eta = 0.01;
            end
            if nargin < 5
                principle_layer = 2;
            end
            if nargin < 4
                KPOINTS_surf = 'KPOINTS_surf';
            end
            if nargin < 3
                fin_dir = 2 ;
            end
            if nargin < 2
                w_range = [-1,1,100];
            end
            w_list = linspace(w_range(1),w_range(2),w_range(3))+fermi;
            H_hr_surf = H_hr < KPOINTS_surf;
            [H00_H01_cell_list_1,H00_H01_cell_list_2] = H_hr_surf.H00_H11_cell_list_gen(fin_dir,principle_layer);
            GREENCAR = HR.GREENCAR_gen(w_list,eta,H00_H01_cell_list_1,H00_H01_cell_list_2,mode);
            DOSCAR_b = HR.DOSCAR_gen(GREENCAR.bulk,'green');
            DOSCAR_l = HR.DOSCAR_gen(GREENCAR.surf_l,'green');
            DOSCAR_r = HR.DOSCAR_gen(GREENCAR.surf_r,'green');
            [klist_l,kpoints_l,kpoints_name] = H_hr_surf.kpath_information();
            w_list = w_list - fermi;
        end
        function [DOSCAR_l,DOSCAR_b,DOSCAR_r,klist1,klist2] = fermiarc(H_hr,w_arc,fin_dir,kmesh,kfermiarc,principle_layer,eta,fermi,mode)
            %nargin
            if nargin < 9
                mode = 'Green_iter';
            end
            if nargin < 8
                fermi = 0;
            end
            if nargin < 7
                eta = 0.01;
            end
            if nargin < 6
                principle_layer = 2;
            end
            if nargin < 5
                kfermiarc = [-0.5, -0.5, 0;... %begin kpoints
                    1.0 , 0,   0;... %  first vector
                    0   , 1.0,   0;... %  second vector
                    ];
            end
            if nargin < 4
                kmesh = [100 100];
            end
            if nargin < 3
                fin_dir = 3 ;
            end
            if nargin < 2
                w_arc = 0;
            end
            w_arc = w_arc + fermi;
            H_hr_arc = H_hr;
            [H_hr_arc.klist_s,klist1,klist2] = H_hr_arc.kmesh3D(kmesh,kfermiarc,'fermiarc');
            switch fin_dir
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
            [H00_H01_cell_list_1,H00_H01_cell_list_2] = H_hr_arc.H00_H11_cell_list_gen(fin_dir,principle_layer);
            GREENCAR = HR.GREENCAR_gen(w_arc,eta,H00_H01_cell_list_1,H00_H01_cell_list_2,mode);
            DOSCAR_b = HR.DOSCAR_gen(GREENCAR.bulk,'green');
            DOSCAR_l = HR.DOSCAR_gen(GREENCAR.surf_l,'green');
            DOSCAR_r = HR.DOSCAR_gen(GREENCAR.surf_r,'green');
            % reshape it
            DOSCAR_b = reshape(DOSCAR_b,kmesh(2),kmesh(1));
            DOSCAR_l = reshape(DOSCAR_l,kmesh(2),kmesh(1));
            DOSCAR_r = reshape(DOSCAR_r,kmesh(2),kmesh(1));
        end
        function [DOSCAR_l,DOSCAR_b,DOSCAR_r,klist1,klist2,E_list] = fermiarc3D(H_hr,w_range,fin_dir,kmesh,kfermiarc,options)
            arguments
                H_hr HR;
                w_range double;
                fin_dir{mustBeMember(fin_dir,[1,2,3])} =3;
                kmesh double = [100 100];
                kfermiarc = [-0.5, -0.5, 0;... %begin kpoints
                    1.0 , 0,   0;... %  first vector
                    0   , 1.0,   0;... %  second vector
                    ];
                options.principle_layer = 1;
                options.eta = 0.01;
                options.fermi = 0;
                options.mode = 'Green_iter';
            end
            
            %nargin
            if length(w_range) == 3
                nw = w_range(3);
                w_range = linspace(w_range(1),w_range(2),w_range(3))+options.fermi;
            else
                w_range = w_range + +options.fermi;
                nw = length(w_range);
            end
            H_hr_arc = H_hr;
            [H_hr_arc.klist_s,klist1,klist2] = H_hr_arc.kmesh3D(kmesh,kfermiarc,'fermiarc');
            switch fin_dir
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
            [H00_H01_cell_list_1,H00_H01_cell_list_2] = H_hr_arc.H00_H11_cell_list_gen(fin_dir,options.principle_layer);
            GREENCAR = HR.GREENCAR_gen(w_range,options.eta,H00_H01_cell_list_1,H00_H01_cell_list_2,options.mode);
            DOSCAR_b = HR.DOSCAR_gen(GREENCAR.bulk,'green');
            DOSCAR_l = HR.DOSCAR_gen(GREENCAR.surf_l,'green');
            DOSCAR_r = HR.DOSCAR_gen(GREENCAR.surf_r,'green');
            % reshape it
            DOSCAR_b = reshape(DOSCAR_b.',kmesh(2),kmesh(1),nw);
            DOSCAR_l = reshape(DOSCAR_l.',kmesh(2),kmesh(1),nw);
            DOSCAR_r = reshape(DOSCAR_r.',kmesh(2),kmesh(1),nw);
            E_list = w_range;
        end
    end
    methods(Static)
        function H_hr = from_POSCAR_Script(varargin)
            H_hr = HR.from_POSCAR_SE(varargin{:});
        end
    end
    %% information
    methods
        function A = list(H_hr,vectorL,options)
            arguments
                H_hr;
                vectorL = [0,0,0];
                options.vpa =true;
                options.digits = 5;
                options.numeric = false;
                options.disp = true;
            end
            if H_hr.vectorhopping
                H_hr = H_hr.GenfromOrth();
            end
            
            H_hr = H_hr.rewrite();
            fprintf('# R1 R2 R3 i j real imag\n');
            if nargin <2
                if H_hr.coe && ~options.numeric
                    if options.vpa
                        A = [sym(H_hr.vectorL),...
                            vpa(real(H_hr.HcoeL), options.digits) ,...
                            vpa(imag(H_hr.HcoeL),options.digits)];
                    else
                        A = [sym(H_hr.vectorL),...
                            (real(H_hr.HcoeL)) ,...
                            imag(H_hr.HcoeL)];
                    end
                    disp(A);
                elseif H_hr.num
                    disp([double(H_hr.vectorL),real(H_hr.HnumL),imag(H_hr.HnumL)]);
                end
            else
                vectorList = double(H_hr.vectorL);
                % if more?
                if isa(vectorL,'double')
                    switch size(vectorL,2)
                        case 2
                            [seq] = find(all(vectorL == vectorList(:,4:5),2));
                        case 3
                            [seq] = find(all(vectorL == vectorList(:,1:3),2));
                        case 5
                            [seq] = find(all(vectorL == vectorList,2));
                    end
                elseif isa(vectorL,'sym')
                    [seq] = find(park.strcontain(string(H_hr.HcoeL),string(vectorL)));
                end
                if H_hr.coe  && ~options.numeric
                    if options.vpa
                        A = [sym(H_hr.vectorL(seq,:)),...
                            vpa(real(H_hr.HcoeL(seq,:)), options.digits) ,...
                            vpa(imag(H_hr.HcoeL(seq,:)),options.digits)];
                    else
                        A = [sym(H_hr.vectorL(seq,:)),real(H_hr.HcoeL(seq,:)) ,imag(H_hr.HcoeL(seq,:))];
                    end
                    if options.disp
                        disp(vpa(A));
                    end
                elseif H_hr.num
                    disp([double(H_hr.vectorL(seq,:)),real(H_hr.HnumL(seq,:)),imag(H_hr.HnumL(seq,:))]);
                end
            end
        end
        function [fig,ax] = show(H_hr,mode,options)
            arguments
                H_hr HR;
                mode char {mustBeMember(mode,{'POSCAR','NRPTS','HOPPING'})} = 'HOPPING';
                options.fig = figure('WindowState','maximized');
                options.scale = 1;
                options.atomscale = 0.3;
                options.TwoD = false;
                options.vectorList = [];
                options.fast = true;
                options.ax = [];
                options.Select = [];
            end
            import vasplib_plot.*;
            Rm_ = H_hr.Rm*options.scale;
            switch mode
                case 'POSCAR'
                    fig  = POSCAR_plot(Rm_,H_hr.sites,H_hr.Atom_name,H_hr.Atom_num);
                case 'NRPTS'
                    if isempty(options.vectorList)
                        vectorList = double(unique(H_hr.vectorL(:,1:3),'rows'));
                    else
                        vectorList = options.vectorList ;
                    end
                    NRPTS_ = size(vectorList,1);
                    positions = H_hr.orbL;
                    sites_ = [];
                    % elementL2 name and num
                    if isempty(H_hr.elementL)
                        elementList = 6*ones(size(positions,1),1);
                    else
                        elementList = H_hr.elementL;
                    end
                    Atom_name_ = repmat(elementList,[NRPTS_,1]);
                    %                     Atom_num_ = repmat(H_hr.elementL,[NRPTS_,1])
                    for i = 1:NRPTS_
                        sites_ = [sites_;positions+vectorList(i,:)] ;
                    end
                    [fig,ax]  = POSCAR_plot(Rm_,sites_,Atom_name_,...
                        'vectorL',vectorList,...
                        'fig',options.fig,...
                        'ax',options.ax,...
                        'scale',options.scale,...
                        'atomscale',options.atomscale,...
                        'TwoD',options.TwoD ,...
                        'fast',options.fast);
                case 'HOPPING'
                    H_hr = H_hr.rewrite();
                    [fig,ax] = H_hr.show('NRPTS',...
                        'scale',options.scale,...
                        'ax',options.ax,...
                        'atomscale',options.atomscale,...
                        'fig',options.fig,...
                        'TwoD',options.TwoD );
                    Rm_ = Rm_ * options.scale;
                    if H_hr.coe
                        if isempty(options.Select)
                            Plot_Hopping(H_hr.vectorL,H_hr.HcoeL,Rm_,H_hr.orbL,'ax',ax,'TwoD',options.TwoD );
                        else
                            SelectL = H_hr.HcoeL == options.Select;
                            Plot_Hopping(H_hr.vectorL(SelectL,:),H_hr.HcoeL(SelectL),Rm_,H_hr.orbL,'ax',ax,'TwoD',options.TwoD );
                        end
                    else
                        if isempty(options.Select)
                            Plot_Hopping(H_hr.vectorL,H_hr.HnumL,Rm_,H_hr.orbL,'ax',ax,'TwoD',options.TwoD );
                        else
                            SelectL = H_hr.HnumL == options.Select;
                            Plot_Hopping(H_hr.vectorL(SelectL,:),H_hr.HnumL(SelectL),Rm_,H_hr.orbL,'ax',ax,'TwoD',options.TwoD );
                        end
                    end
                    %                     if options.TwoD
                    %                         view(ax,0,90);
                    %                         %% lighting
                    %                         l = light(ax);
                    %                         l.Position = [1,1,3];
                    %                     end
                    if options.TwoD
                        view(ax,0,90);
                        %% lighting
                        l = light(ax);
                        l.Position = [0,0,1];
                    end
                otherwise
            end
        end
        function Hout = printout(H_hr,print_list,mode)
            if nargin < 2
                print_list = 1:H_hr.NRPTS;
            end
            if nargin < 3
                mode = 'HcoeL';
            end
            switch mode
                case 'HnumL-sum'
                    Hout = sum(H_hr.HnumL,3);
                case 'HnumL'
                    %
                    for i = 1:length(print_list)
                        disp([H_hr.vectorL(print_list(i),:)]);
                        disp(H_hr.HnumL(:,:,print_list(i)));
                    end
                    Hout = H_hr.HnumL(:,:,print_list);
                case 'HcoeL'
                    %
                    for i = 1:length(print_list)
                        disp([H_hr.vectorL(print_list(i),:)]);
                        disp(H_hr.HcoeL(:,:,print_list(i)));
                    end
                    Hout = H_hr.HcoeL(:,:,print_list);
                case 'HcoeL-sum'
                    Hout =  sum(H_hr.HcoeL,3);
                otherwise
                    fprintf("usage: H_hr.printout([1,3,4],'HcoeL')\n");
                    fprintf("usage: H_hr.printout([1,3,4],'HnumL')\n");
                    fprintf("Optional: HcoeL,HnumL,HnumL-sum,HcoeL-sum,\n");
            end
        end
        function [VarInit,EQL2] =GetInit(H_hr,H_hr2,vectorL)
            
            H_hr = H_hr.rewrite();
            H_hr = H_hr.simplify();
            H_hr2 =H_hr2.rewrite();
            %Choosing Vector
            if nargin < 3
                [ia,ic] = ismember(H_hr.vectorL,H_hr2.vectorL,'row');
                HcoeLtmp = H_hr.HcoeL(ia);
                HnumLtmp = H_hr2.HnumL(ic);
            else
                if size(vectorL,1) == 1
                    ic = all(int32(vectorL) == H_hr.vectorL(:,1:3),2);
                    HcoeLtmp = H_hr.HcoeL(ic);
                    vector1 = H_hr.vectorL(ic,:);
                    ic = all(int32(vectorL) == H_hr2.vectorL(:,1:3),2);
                    HnumLtmp = H_hr2.HnumL(ic);
                    vector2 = H_hr2.vectorL(ic,:);
                    [ia,ic] = ismember(vector1,vector2,'row');
                    HcoeLtmp = HcoeLtmp(ia);
                    HnumLtmp = HnumLtmp(ic);
                end
            end
            
            [Unique_term,ia,ic] = unique(HcoeLtmp);
            for i  = 1:length(ia)
                EQL(i,:) = Unique_term(i) == max(HnumLtmp(ic == i));
            end
            [Unique_term2,ia,ic] = unique(lhs(EQL));
            HcoeLtmp2 = rhs(EQL);
            for i  = 1:length(ia)
                EQL2(i,:) = Unique_term2(i) == max(HcoeLtmp2(ic == i));
            end
            %   try
            VarInit = solve(EQL2,H_hr.symvar_list);
            EQL2 = vpa(EQL2);
        end
        function H_atom_soc = H_atom_soc(H_hr)
            % check
            if isempty(H_hr.quantumL)
                error('you should provide quantum number list');
            end
            if size(H_hr.quantumL,1) ~= H_hr.WAN_NUM
                error('size quantum number list qrong');
            end
            if mod(size(H_hr.quantumL,1),2) ~= 0
                error('size quantum number list odd, cant be spinful');
            end
            if sum(H_hr.quantumL(:,4)) ~= 0
                error('spin up ~= spin dn');
            end
            H_atom_soc = sym(zeros(H_hr.WAN_NUM));
            for i = 1:H_hr.WAN_NUM
                l1 = H_hr.quantumL(i,2);
                m1 = H_hr.quantumL(i,3);
                s1 = H_hr.quantumL(i,4);
                element1 = H_hr.elementL(i);
                orb1 = H_hr.orbL(i,:);
                for j =  1:H_hr.WAN_NUM
                    l2 = H_hr.quantumL(j,2);
                    m2 = H_hr.quantumL(j,3);
                    s2 = H_hr.quantumL(j,4);
                    element2 = H_hr.elementL(j);
                    orb2 = H_hr.orbL(j,:);
                    if isequal(orb1,orb2) && element1 == element2 && l1 == l2
                        H_atom_soc(i,j) = soc_term_gen(l1,l2,m1,m2,s1,s2,element1);
                    end
                end
            end
        end
        function kloop = kloop_gen(H_hr,input_mat,mode)
            import spglib_matlab.*;
            if nargin < 3
                mode = 'kline';
            end
            if nargin < 2
                input_mat = ...
                    [ 0, 0,  0 ;... %kpoints_f
                    0, 0,  1 ;... %fin_dir
                    0, 1,101];%start end nodes
            end
            %%
            kpoints_f = input_mat(1,:);
            if strcmp(mode,'kline')
                fin_dir_list = input_mat(2,:);
                fin_dir = fin_dir_list ==1;
                nodes = input_mat(3,3);
                if nodes <2
                    error('nodes >= 2!!!!');
                end
                k_start = input_mat(3,1);k_end = input_mat(3,2);
                dk = (k_end-k_start)/(nodes-1);
                klists =repmat(kpoints_f,[nodes,1]);
                klists(:,fin_dir) = (k_start:dk:k_end).';
                kloop = klists;
            elseif  strcmp(mode,'kplane')
                kpoints_r = kpoints_f *H_hr.Gk;
                n_vector = input_mat(2,:);
                nx =n_vector(1);ny =n_vector(1);nz = n_vector(3);
                nodes = input_mat(3,3);
                dk = input_mat(3,1);
                dtheta = 360/(nodes-1);
                theta_init_d = input_mat(3,2);
                if norm(abs(n_vector)-[1,0,0])<1e-8
                    disp('rotation along x axis.');
                    dkx = 0;
                    dky = 0;
                    dkz = dk;
                    dk_default_init = [dkx,dky,dkz];
                elseif ny~=0 && nz~=0
                    dkx = 0;
                    dky = sqrt(nz^2/(ny^2+nz^2)) * dk;
                    dkz = nz*dky/ny;
                    dk_default_init = [dkx,dky,dkz];
                elseif ny==0 && nz~=0
                    dkx = 0;
                    dky = dk;
                    dkz = 0;
                    dk_default_init = [dkx,dky,dkz];
                elseif ny~=0 && nz==0
                    dkx = 0;
                    dky = 0;
                    dkz = dk;
                    dk_default_init = [dkx,dky,dkz];
                else
                    dk_default_init = [dk,0,0];
                end
                rotation_mat_init = spglib_matlab.nTheta2RotationMat(n_vector,theta_init_d);
                dk_init = ((rotation_mat_init*dk_default_init.').'*H_hr.Gk).';
                kloop = zeros(nodes,3);
                kloop(1,:) = dk_init.';
                for i =2:nodes
                    rotation_mat = spglib_matlab.nTheta2RotationMat(n_vector,theta_init_d+(i-1)*dtheta);
                    kloop(i,:) =(rotation_mat*dk_init).';
                end
                %kloopf = kloop/H_hr.Gk+kpoints_f;
                kloop =  kloop + kpoints_r;
            else
                
            end
        end
        function stucture
        end
        function [Rnn,nn_store_smart,Atom_store_smart,Rnn_map] = nn_information(H_hr,silence)
            if nargin < 2
                silence = true;
            end
            Rnn = H_hr.Rnn;
            if ~isempty(H_hr.nn_store_smart)
                nn_store_smart = H_hr.nn_store_smart;
            else
                nn_store_smart = H_hr.nn_store;
            end
            try
                Atom_store_smart = H_hr.Atom_store_smart;
            catch
                Atom_store_smart = [];
            end
            Rnn_map = H_hr.Rnn_map ;
            if ~silence
                for i = 1:length(Rnn)
                    fprintf('the %3d th Rnn vector: %7.4f Angstrom\n',i,Rnn(i));
                end
            end
        end
    end
    %% SKTB
    methods
        function H_hr = H_TBSK_gen(H_hr,options)
            arguments
                H_hr HR;
                options.level_cut {mustBeInteger} = -1;
                options.onsite logical = false;
                options.per_dir double = [1,1,1];
                options.chiral logical = false;
                options.spin logical = false;
                options.method  = 'nn_sparse';
                options.rough = false;
                options.vectorL = [];
                %options.overlap = false;
            end
            % test_mode = true;
            H_hr.num = false;
            H_hr.coe = true;
            % check nn_store
            %
            if isempty(H_hr.nn_store)
                if isempty(H_hr.nn_store_smart)
                    error('You have not run the function: nn/nn_smart');
                else
                    if strcmp(options.method,'nn_sparse')
                        error('You should use the key value: ''method'',''nn_smart'' as input.');
                    end
                end
            end
            % init
            N_orbit = H_hr.WAN_NUM;
            %H_hr = H_hr.rewrite();% list mode enforced
            if options.level_cut >0
                select_nn_store = H_hr.nn_store(H_hr.nn_store(:,10)<=options.level_cut,:);
            else
                select_nn_store = H_hr.nn_store;
            end
            % try to speed up
            if options.chiral
                i_L = select_nn_store(:,1);
                j_L = select_nn_store(:,2);
                element1_L = H_hr.elementL(i_L);
                element2_L =  H_hr.elementL(j_L);
                select_nn_store =select_nn_store(element1_L ~= element2_L,:);
            end
            if options.spin
                i_L = select_nn_store(:,1);
                j_L = select_nn_store(:,2);
                spin1_L = H_hr.quantumL(i_L,4);
                spin2_L = H_hr.quantumL(j_L,4);
                select_nn_store =select_nn_store(spin1_L == spin2_L,:);
            end
            %
            % nselect_nn_store = size(select_nn_store,1);
            % sort first
            Rvector_L = select_nn_store(:,6:8);
            [Rvector_L_unique,sorted_label,cut_slice] = vasplib.cut_tools(Rvector_L);
            select_nn_store = select_nn_store(sorted_label,:);
            %
            i_L = select_nn_store(:,1);
            j_L = select_nn_store(:,2);
            Rlength_L = select_nn_store(:,9);
            l_L = select_nn_store(:,3)./Rlength_L;
            m_L = select_nn_store(:,4)./Rlength_L;
            n_L = select_nn_store(:,5)./Rlength_L;
            nn_level_L = select_nn_store(:,10);
            L_1_L = H_hr.quantumL(i_L,2);
            L_2_L = H_hr.quantumL(j_L,2);
            m_1_L = H_hr.quantumL(i_L,3);
            m_2_L = H_hr.quantumL(j_L,3);
            NRPTS_ = size(cut_slice,1);
            fprintf('Generating Symbolic Hopping term ...')
            if ~H_hr.overlap
                SymHopping_L = HR.TBSK_Var_gen(L_1_L,L_2_L,m_1_L,m_2_L,nn_level_L,l_L,m_L,n_L);
            else
                [SymHopping_L,SymOverlap_L] = HR.TBSK_Var_gen(L_1_L,L_2_L,m_1_L,m_2_L,nn_level_L,l_L,m_L,n_L,'overlap',true);
            end
            fprintf('< ok >\n')
            %%
            if strcmp(H_hr.Type,'mat')
                H_hr = H_hr.add_empty_one(Rvector_L);
            end
            %%
            pb = vasplib_tool_outer.CmdLineProgressBar('Setting NRPT ');
            for i = 1:NRPTS_
                if ~H_hr.overlap
                    H_hr = H_hr.set_hop(...
                        SymHopping_L(cut_slice(i,1):cut_slice(i,2)),...
                        i_L(cut_slice(i,1):cut_slice(i,2)),...
                        j_L(cut_slice(i,1):cut_slice(i,2)),...
                        Rvector_L_unique(i,:),'sym');
                else
                    H_hr = H_hr.set_hop(...
                        SymHopping_L(cut_slice(i,1):cut_slice(i,2)),...
                        i_L(cut_slice(i,1):cut_slice(i,2)),...
                        j_L(cut_slice(i,1):cut_slice(i,2)),...
                        Rvector_L_unique(i,:),'sym');
                    H_hr = H_hr.set_overlap(...
                        SymOverlap_L(cut_slice(i,1):cut_slice(i,2)),...
                        i_L(cut_slice(i,1):cut_slice(i,2)),...
                        j_L(cut_slice(i,1):cut_slice(i,2)),...
                        Rvector_L_unique(i,:),'sym');
                end
                pb.print(i,NRPTS_,' ...');
            end
            if options.onsite
                for i=1:N_orbit
                    onsite_sym_name = "E__"+string(H_hr.elementL(i,1))+"_"...
                        +string(H_hr.quantumL(i,2));...
                        H_hr = H_hr.set_hop_single(sym(onsite_sym_name,'real'),i,i,[0 0 0],'symadd');
                end
            end
            if H_hr.overlap
                for i=1:N_orbit
                    H_hr = H_hr.set_overlap_single(1,i,i,[0 0 0],'symadd');
                end
            end
            pb.delete();
        end
        function H_hr = H_TBSK_gen_sparse(H_hr,options)
            arguments
                H_hr HR;
                options.level_cut {mustBeInteger} = -1;
                options.onsite logical = false;
                options.per_dir double = [1,1,1];
                options.chiral logical = false;
                options.spin logical = false;
                options.method  = 'nn_sparse';
                options.rough = false;
                options.vectorL = [];
                options.deltarule  double{mustBeMember(options.deltarule,[0,1,2])}= 0;
                options.alpharule  double{mustBeMember(options.alpharule,[0,1,2])}= 0;
                options.para = struct();
                options.Rd double = -1;
                %options.overlap = false;
            end
            % At present, we only consider numerically gen sparse model
            H_hr.num = true;
            H_hr.coe = false;
            % check nn_store
            if isempty(H_hr.nn_store)
                if isempty(H_hr.nn_store_smart)
                    error('You have not run the function: nn/nn_smart');
                else
                    if strcmp(options.method,'nn_sparse')
                        error('You should use the key value: ''method'',''nn_smart'' as input.');
                    end
                end
            end
            % init
            N_orbit = H_hr.WAN_NUM;
            %H_hr = H_hr.rewrite();% list mode enforced
            if options.level_cut >0
                select_nn_store = H_hr.nn_store(H_hr.nn_store(:,10)<=options.level_cut,:);
            else
                select_nn_store = H_hr.nn_store;
            end
            % try to speed up
            if options.chiral
                i_L = select_nn_store(:,1);
                j_L = select_nn_store(:,2);
                element1_L = H_hr.elementL(i_L);
                element2_L =  H_hr.elementL(j_L);
                select_nn_store =select_nn_store(element1_L ~= element2_L,:);
            end
            if options.spin
                i_L = select_nn_store(:,1);
                j_L = select_nn_store(:,2);
                spin1_L = H_hr.quantumL(i_L,4);
                spin2_L = H_hr.quantumL(j_L,4);
                select_nn_store =select_nn_store(spin1_L ~= spin2_L,:);
            end
            %
            % nselect_nn_store = size(select_nn_store,1);
            % sort first
            Rvector_L = select_nn_store(:,6:8);
            [Rvector_L_unique,sorted_label,cut_slice] = vasplib.cut_tools(Rvector_L);
            select_nn_store = select_nn_store(sorted_label,:);
            %
            i_L = select_nn_store(:,1);
            j_L = select_nn_store(:,2);
            Rlength_L = select_nn_store(:,9);
            l_L = select_nn_store(:,3)./Rlength_L;
            m_L = select_nn_store(:,4)./Rlength_L;
            n_L = select_nn_store(:,5)./Rlength_L;
            nn_level_L = select_nn_store(:,10);
            L_1_L = H_hr.quantumL(i_L,2);
            L_2_L = H_hr.quantumL(j_L,2);
            m_1_L = H_hr.quantumL(i_L,3);
            m_2_L = H_hr.quantumL(j_L,3);
            NRPTS_ = size(cut_slice,1);
            fprintf('Generating Numeric Hopping term ...')
            if ~H_hr.overlap
                NumHopping_L = HR.TBSK_Var_gen_sparse(L_1_L,L_2_L,m_1_L,m_2_L,Rlength_L,l_L,m_L,n_L,...
                    'para',options.para,'Rd',options.Rd,...
                    'deltarule',options.deltarule ,...
                    'alpharule',options.alpharule ...
                    );
            else
                [NumHopping_L,NumOverlap_L] = HR.TBSK_Var_gen(L_1_L,L_2_L,m_1_L,m_2_L,nn_level_L,l_L,m_L,n_L,'overlap',true);
                % error('Not support yet');
            end
            fprintf('< ok >\n')
            %%
            if strcmp(H_hr.Type,'mat')
                H_hr = H_hr.add_empty_one(Rvector_L);
            end
            %%
            pb = vasplib_tool_outer.CmdLineProgressBar('Setting NRPT ');
            for i = 1:NRPTS_
                if ~H_hr.overlap
                    H_hr = H_hr.set_hop(...
                        NumHopping_L(cut_slice(i,1):cut_slice(i,2)),...
                        i_L(cut_slice(i,1):cut_slice(i,2)),...
                        j_L(cut_slice(i,1):cut_slice(i,2)),...
                        Rvector_L_unique(i,:),'set');
                else
                    H_hr = H_hr.set_hop(...
                        NumHopping_L(cut_slice(i,1):cut_slice(i,2)),...
                        i_L(cut_slice(i,1):cut_slice(i,2)),...
                        j_L(cut_slice(i,1):cut_slice(i,2)),...
                        Rvector_L_unique(i,:),'set');
                    H_hr = H_hr.set_overlap(...
                        NumOverlap_L(cut_slice(i,1):cut_slice(i,2)),...
                        i_L(cut_slice(i,1):cut_slice(i,2)),...
                        j_L(cut_slice(i,1):cut_slice(i,2)),...
                        Rvector_L_unique(i,:),'set');
                end
                pb.print(i,NRPTS_,' ...');
            end
            if options.onsite
                for i=1:N_orbit
                    onsite_sym_name = "E__"+string(H_hr.elementL(i,1))+"_"...
                        +string(H_hr.quantumL(i,2));...
                        % find the value in para.symvar
                    % modify it later
                    H_hr = H_hr.set_hop_single(sym(onsite_sym_name,'real'),i,i,[0 0 0],'add');
                end
            end
            if H_hr.overlap
                for i=1:N_orbit
                    H_hr = H_hr.set_overlap_single(1,i,i,[0 0 0],'add');
                end
            end
            pb.delete();
        end
    end
    %% sewing HR
    methods
        function H_hr_n = Connect2D(H_hr_n,H_hr_unitcell,opt)
            %--------  init  --------
            arguments
                H_hr_n HR;
                H_hr_unitcell HR;
                opt.Sequ double = [1:20];
                opt.findir double = 2;
                opt.Subsequ {mustBeMember(opt.Subsequ,{'normal','reverse','combined'})} = 'normal';
                opt.reverse_sequ double = [H_hr_unitcell.WAN_NUM:-1:1];
                opt.combined_norm double = [2:19];
            end
            
            % prepare
            WANNUM = H_hr_n.WAN_NUM;
            WANNUM_unit = H_hr_unitcell.WAN_NUM;
            
            if opt.findir==2
                Nrep_y = length(opt.Sequ);
                Nrep_x = WANNUM/WANNUM_unit/Nrep_y;
                Nrep_tot = Nrep_x*Nrep_y;
                for k = -1:1:1 % to consider 3 NRTPs with -1 0 1
                    slct = ismember(H_hr_unitcell.vectorL(:,1:3),[1 k 0],'rows');
                    if all(slct==0)
                        continue;
                    end
                    if strcmp(opt.Subsequ,'normal')
                        [row,col,num] = find(H_hr_unitcell.HnumL(:,:,slct) ); %HR list gen
                    elseif strcmp(opt.Subsequ,'reverse')
                        %                         col = flip(col); % low
                        %                         disp('This reverse mode only work for reflection symmetric inter-cell hoppings with nocross hopping in square lattice');
                        HnumL_tmp = H_hr_unitcell.HnumL(:,:,slct) ; %flip(H_hr_unitcell.HnumL(:,:,slct) ); % flip to reverse orbitals by their spacial sequence
                        %                         [row_norm,col_norm,num_norm] = find(HnumL_tmp );
                        HnumL_reverse = HnumL_tmp(:,opt.reverse_sequ);
                        [row,col,num] = find(HnumL_reverse ); %HR list gen
                        disp('This reverse mode only work for reflection symmetric inter-cell hoppings, and may have bugs for double layer materials');
                        if num~=flip(num)
                            warning('This unitcell HR may not be suitable for the reverse mode. Check the relection symmetry for Hr_unitcell!');
                        end
                    elseif strcmp(opt.Subsequ,'combined')
                        HnumL_tmp = H_hr_unitcell.HnumL(:,:,slct) ; %flip(H_hr_unitcell.HnumL(:,:,slct) ); % flip to reverse orbitals by their spacial sequence
                        [row_norm,col_norm,num_norm] = find(HnumL_tmp );
                        HnumL_reverse = HnumL_tmp(:,opt.reverse_sequ);
                        [row_rev,col_rev,num_rev] = find(HnumL_reverse ); %HR list gen
                        disp('This combined mode only work for reflection symmetric inter-cell hoppings, and may have bugs for double layer materials');
                        if num_rev~=flip(num_rev)
                            warning('This unitcell HR may not be suitable for the combined mode. Check the relection symmetry for Hr_unitcell!');
                        end
                        if max(opt.combined_norm)>WANNUM || min(opt.combined_norm)<1
                            warning('Check the combined_norm range!');
                        end
                    end
                    % set hop
                    if strcmp(opt.Subsequ,'combined')
                        count = 1;norm_length = length(opt.combined_norm);
                        for i = 1:length(opt.Sequ)
                            labR = Nrep_y*(Nrep_x-1)+i;
                            labL =  opt.Sequ(i) +k;
                            if labL>0&&labL<=Nrep_y % discard those out of the range
                                if count<=norm_length && i == opt.combined_norm(count)
                                    for j =1: length(row_norm)
                                        H_hr_n = H_hr_n.set_hop(...
                                            num_norm(j),WANNUM_unit*(labR-1)+row_norm(j),WANNUM_unit*(labL-1)+col_norm(j),[1 0 0],'add');
                                        H_hr_n = H_hr_n.set_hop(...
                                            conj(num_norm(j)),WANNUM_unit*(labL-1)+col_norm(j),WANNUM_unit*(labR-1)+row_norm(j),[-1 0 0],'add');
                                    end
                                    count = count +1;
                                else
                                    for j =1: length(row_rev)
                                        H_hr_n = H_hr_n.set_hop(...
                                            num_rev(j),WANNUM_unit*(labR-1)+row_rev(j),WANNUM_unit*(labL-1)+col_rev(j),[1 0 0],'add');
                                        H_hr_n = H_hr_n.set_hop(...
                                            conj(num_rev(j)),WANNUM_unit*(labL-1)+col_rev(j),WANNUM_unit*(labR-1)+row_rev(j),[-1 0 0],'add');
                                    end
                                end
                            end
                        end
                    else
                        for i = 1:length(opt.Sequ)
                            % set hop
                            labR = Nrep_y*(Nrep_x-1)+i;
                            %                 labL = Nrep_x*mod( (opt.Sequ(i) -1+k),Nrep_y); % assume
                            %                 PBC in the y dir
                            labL =  opt.Sequ(i) +k;
                            if labL>0&&labL<=Nrep_y % discard those out of the range
                                for j =1: length(row)
                                    H_hr_n = H_hr_n.set_hop(...
                                        num(j),WANNUM_unit*(labR-1)+row(j),WANNUM_unit*(labL-1)+col(j),[1 0 0],'add');
                                    H_hr_n = H_hr_n.set_hop(...
                                        conj(num(j)),WANNUM_unit*(labL-1)+col(j),WANNUM_unit*(labR-1)+row(j),[-1 0 0],'add');
                                end
                            end
                        end
                    end
                end
            elseif opt.findir==1
                % prepare
                Nrep_x = length(opt.Sequ);
                Nrep_y = WANNUM/WANNUM_unit/Nrep_x;
                for k = -1:1:1 % to consider 3 NRTPs with -1 0 1
                    slct = ismember(H_hr_unitcell.vectorL(:,1:3),[k 1 0],'rows');
                    if all(slct==0)
                        continue;
                    end                          
                    if strcmp(opt.Subsequ,'normal')
                        [row,col,num] = find(H_hr_unitcell.HnumL(:,:,slct) ); %HR list gen
                    elseif strcmp(opt.Subsequ,'reverse')
                        %                         col = flip(col); % low
                        HnumL_tmp = H_hr_unitcell.HnumL(:,:,slct) ; %flip(H_hr_unitcell.HnumL(:,:,slct) ); % flip to reverse orbitals by their spacial sequence
                        HnumL_reverse = HnumL_tmp(:,opt.reverse_sequ);
                        [row,col,num] = find(HnumL_reverse ); %HR list gen        
                        disp('This reverse mode only work for reflection symmetric inter-cell hoppings, and may have bugs for double layer materials');
                        if num~=flip(num)
                            warning('This unitcell HR may not be suitable for the reverse mode. Check the relection symmetry for Hr_unitcell!');
                        end
                    end
                    for i = 1:length(opt.Sequ)
                        % set hop
                        %                 labR = Nrep_x*(Nrep_y-1) + i;
                        %                 labL = opt.Sequ(i);
                        labR = i*(Nrep_y);
                        tmp = (opt.Sequ(i)-1+k);
                        labL = Nrep_y*tmp + 1;
                        if tmp>=0&&tmp<Nrep_x
                            for j =1: length(row)
                                H_hr_n = H_hr_n.set_hop(...
                                    num(j),WANNUM_unit*(labR-1)+row(j),WANNUM_unit*(labL-1)+col(j),[0 1 0],'add');
                                H_hr_n = H_hr_n.set_hop(...
                                    conj(num(j)),WANNUM_unit*(labL-1)+col(j),WANNUM_unit*(labR-1)+row(j),[0 -1 0],'add');
                            end
                        end
                    end
                end
            end
            %         H_hr_n.Rm = H_hr_n.Rm*diag([Nrep_x Nrep_y 1]);
        end
    end
    % -------------------------------
    methods(Static,Hidden,Access= protected)
        function varargout = TBSK_Var_gen(L_1_L,L_2_L,m_1_L,m_2_L,nn_level_L,l_L,m_L,n_L,options)
            arguments
                L_1_L double{mustBeInteger};
                L_2_L double{mustBeInteger};
                m_1_L double{mustBeInteger};
                m_2_L double{mustBeInteger};
                nn_level_L double{mustBeInteger} = -1;
                l_L  = 0;
                m_L  = 0;
                n_L  = 0;
                options.overlap logical = false;
                %                 l double{mustBeInRange(l,0,1)} = 0;
                %                 m double{mustBeInRange(m,0,1)} = 0;
                %                 n double{mustBeInRange(n,0,1)} = 0;
            end
            nhopping = length(L_1_L);
            if nhopping == 1
                switch nargout
                    case 3
                        [varargout{1},varargout{2},varargout{3}] = ...
                            HR.TBSK_Var_gen_single(L_1_L,L_2_L,m_1_L,m_2_L,nn_level_L,l_L,m_L,n_L,'overlap', options.overlap );
                    case 2
                        [varargout{1},varargout{2}] = ...
                            HR.TBSK_Var_gen_single(L_1_L,L_2_L,m_1_L,m_2_L,nn_level_L,l_L,m_L,n_L,'overlap', options.overlap );
                    case 1
                        [varargout{1}] = ...
                            HR.TBSK_Var_gen_single(L_1_L,L_2_L,m_1_L,m_2_L,nn_level_L,l_L,m_L,n_L,'overlap', options.overlap );
                end
                return;
            end
            % <https://en.wikipedia.org/wiki/Tight_binding>
            % L_1,L_2 : orbital angular momentum: s 0 p 1 d 2 f 3
            % m_1,m_2 : magnetic quantum number:
            %           - s 0
            %           - p y:-1 z:0 x:1 % because py = i/sart(2)*(Y-1
            %           +Y1), we label this i as a minus number
            %           - d xy:-2 yz:-1 z2:0 xz:1 x2-y2:2
            %           - f
            %           y(3x2-y2):-3 xyz:-2  yz2:-1
            %           z3:0
            %           x(x2-3y2):3  z(x2-y2):2  xz2:1
            
            % Coeff = zeros(nhopping,3);
            % when generating a batch of Hoppings, we must code vectorized
            exchange_list = L_1_L > L_2_L;
            L_3_L = L_2_L(exchange_list);
            L_2_L(exchange_list) = L_1_L(exchange_list);
            L_1_L(exchange_list) = L_3_L;
            l_L(exchange_list) = - l_L(exchange_list);
            m_L(exchange_list) = - m_L(exchange_list) ;
            n_L(exchange_list) = - n_L(exchange_list);
            % Coeff first
            Coeff = zeros(nhopping,3);
            for i = 1:nhopping
                Coeff(i,:) = HR.TBSK_Coeff_gen(L_1_L(i),L_2_L(i),m_1_L(i),m_2_L(i),l_L(i),m_L(i),n_L(i));
            end
            nn_level_str_L= strcat(repmat('_',[nhopping,1]),string(nn_level_L));
            if nargout == 1 && options.overlap
                varargout{1} = Coeff;
                return;
            end
            % make a dist
            %             E_hop_sym_out = repmat(sym(0),[nhopping,1]);
            %             E_integral_sym_out = repmat(sym(0),[nhopping,1]);
            l2str_map=['s';'p';'d';'f'];
            str_base_L = [l2str_map(L_1_L+1),l2str_map(L_2_L+1)];
            V_str_base_L = [repmat('V',[nhopping,1]),str_base_L];
            V_str_1_L = strcat(V_str_base_L,repmat('S',[nhopping,1]),nn_level_str_L);
            V_str_2_L = strcat(V_str_base_L,repmat('P',[nhopping,1]),nn_level_str_L);
            V_str_3_L = strcat(V_str_base_L,repmat('D',[nhopping,1]),nn_level_str_L);
            if options.overlap
                S_str_base_L = [repmat('S',[nhopping,1]),str_base_L];
                S_str_1_L = strcat(S_str_base_L,repmat('S',[nhopping,1]),nn_level_str_L);
                S_str_2_L = strcat(S_str_base_L,repmat('P',[nhopping,1]),nn_level_str_L);
                S_str_3_L = strcat(S_str_base_L,repmat('D',[nhopping,1]),nn_level_str_L);
            end
            if nargout == 3
                varargout{1} = Coeff;
                varargout{2} = [V_str_1_L,V_str_2_L,V_str_3_L];
                varargout{3} = [S_str_1_L,S_str_2_L,S_str_3_L];
                return;
            end
            E_hop_sym_out = sym(V_str_1_L,'real').*Coeff(:,1) +...
                sym(V_str_2_L,'real').*Coeff(:,2) +...
                sym(V_str_3_L,'real').*Coeff(:,3) ;
            if options.overlap
                E_integral_sym_out = sym(S_str_1_L,'real').*Coeff(:,1) +...
                    sym(S_str_2_L,'real').*Coeff(:,2) +...
                    sym(S_str_3_L,'real').*Coeff(:,3) ;
            end
            switch nargout
                case 1
                    varargout{1} = E_hop_sym_out;
                case 2
                    varargout{1} = E_hop_sym_out;
                    varargout{2} = E_integral_sym_out;
                otherwise
                    error('!!');
            end
        end
        function varargout = TBSK_Var_gen_sparse(L_1_L,L_2_L,m_1_L,m_2_L,Rlength_L,l_L,m_L,n_L,options)
            arguments
                L_1_L double{mustBeInteger};
                L_2_L double{mustBeInteger};
                m_1_L double{mustBeInteger};
                m_2_L double{mustBeInteger};
                Rlength_L double{mustBeNonnegative};
                l_L  = 0;
                m_L  = 0;
                n_L  = 0;
                options.overlap logical = false;
                options.deltarule  double{mustBeMember(options.deltarule,[0,1,2])}= 0;
                options.alpharule  double{mustBeMember(options.alpharule,[0,1,2])}= 0;
                options.para = struct();
                options.Rd double = -1;
            end
            nhopping = length(L_1_L);
            % when generating a batch of Hoppings, we must code vectorized
            exchange_list = L_1_L > L_2_L;
            L_3_L = L_2_L(exchange_list);
            L_2_L(exchange_list) = L_1_L(exchange_list);
            L_1_L(exchange_list) = L_3_L;
            l_L(exchange_list) = - l_L(exchange_list);
            m_L(exchange_list) = - m_L(exchange_list) ;
            n_L(exchange_list) = - n_L(exchange_list);
            % Coeff first
            Coeff = zeros(nhopping,3);
            for i = 1:nhopping
                Coeff(i,:) = HR.TBSK_Coeff_gen(L_1_L(i),L_2_L(i),m_1_L(i),m_2_L(i),l_L(i),m_L(i),n_L(i));
            end
            % Base Second
            l2str_map=['s';'p';'d';'f'];
            str_base_L = [l2str_map(L_1_L+1),l2str_map(L_2_L+1)];
            V_str_base_L = [repmat('V',[nhopping,1]),str_base_L];
            V_str_1_L = string(strcat(V_str_base_L,repmat('S',[nhopping,1])));
            V_str_2_L = string(strcat(V_str_base_L,repmat('P',[nhopping,1])));
            V_str_3_L = string(strcat(V_str_base_L,repmat('D',[nhopping,1])));
            strvar = options.para.strvar;
            numvar = options.para.numvar;
            [V_seq_1_L_contain,V_seq_1_L] = ismember(V_str_1_L,strvar);
            [V_seq_2_L_contain,V_seq_2_L] = ismember(V_str_2_L,strvar);
            [V_seq_3_L_contain,V_seq_3_L] = ismember(V_str_3_L,strvar);
            V_seq_1_L = V_seq_1_L(V_seq_1_L_contain);
            V_seq_2_L = V_seq_2_L(V_seq_2_L_contain);
            V_seq_3_L = V_seq_3_L(V_seq_3_L_contain);
            %
            V_num_1_L = zeros(nhopping,1);V_num_2_L = V_num_1_L;V_num_3_L = V_num_1_L;
            %
            V_num_1_L(V_seq_1_L_contain) = numvar(V_seq_1_L);
            V_num_2_L(V_seq_2_L_contain) = numvar(V_seq_2_L);
            V_num_3_L(V_seq_3_L_contain) = numvar(V_seq_3_L);
            % Scale third
            V_scale1_L = ones(nhopping,1);V_scale2_L = V_scale1_L;V_scale3_L = V_scale1_L;
            if options.alpharule > 0
                if options.alpharule == 1 && length(options.para.delta) == 1
                    V_scale1_L(V_seq_1_L_contain) = exp(-(Rlength_L(V_seq_1_L_contain) - options.Rd(V_seq_1_L).')./options.Rd(V_seq_1_L)*options.para.delta);
                    V_scale2_L(V_seq_2_L_contain) = exp(-(Rlength_L(V_seq_2_L_contain) - options.Rd(V_seq_2_L).')./options.Rd(V_seq_2_L)*options.para.delta);
                    V_scale3_L(V_seq_3_L_contain) = exp(-(Rlength_L(V_seq_3_L_contain) - options.Rd(V_seq_3_L).')./options.Rd(V_seq_3_L)*options.para.delta);
                elseif options.alpharule == 2 && length(options.para.delta) > 1
                    V_scale1_L(V_seq_1_L_contain) = exp(-(Rlength_L(V_seq_1_L_contain) - options.Rd(V_seq_1_L).')./options.Rd(V_seq_1_L).*options.para.delta(V_seq_1_L).');
                    V_scale2_L(V_seq_2_L_contain) = exp(-(Rlength_L(V_seq_2_L_contain) - options.Rd(V_seq_2_L).')./options.Rd(V_seq_2_L).*options.para.delta(V_seq_2_L).');
                    V_scale3_L(V_seq_3_L_contain) = exp(-(Rlength_L(V_seq_3_L_contain) - options.Rd(V_seq_3_L).')./options.Rd(V_seq_3_L).*options.para.delta(V_seq_3_L).');
                end
            end
            if options.deltarule > 0
                if options.deltarule == 1 && length(options.para.delta) == 1
                    V_scale1_L(V_seq_1_L_contain) = exp(-(Rlength_L(V_seq_1_L_contain) - options.Rd(V_seq_1_L).')/options.para.delta);
                    V_scale2_L(V_seq_2_L_contain) = exp(-(Rlength_L(V_seq_2_L_contain) - options.Rd(V_seq_2_L).')/options.para.delta);
                    V_scale3_L(V_seq_3_L_contain) = exp(-(Rlength_L(V_seq_3_L_contain) - options.Rd(V_seq_3_L).')/options.para.delta);
                elseif options.deltarule == 2 && length(options.para.delta) > 1
                    V_scale1_L(V_seq_1_L_contain) = exp(-(Rlength_L(V_seq_1_L_contain) - options.Rd(V_seq_1_L).')./options.para.delta(V_seq_1_L).');
                    V_scale2_L(V_seq_2_L_contain) = exp(-(Rlength_L(V_seq_2_L_contain) - options.Rd(V_seq_2_L).')./options.para.delta(V_seq_2_L).');
                    V_scale3_L(V_seq_3_L_contain) = exp(-(Rlength_L(V_seq_3_L_contain) - options.Rd(V_seq_3_L).')./options.para.delta(V_seq_3_L).');
                end
            end
            HoppingL = Coeff(:,1).*V_num_1_L.*V_scale1_L + ...
                Coeff(:,2).*V_num_2_L.*V_scale2_L + ...
                Coeff(:,3).*V_num_3_L.*V_scale3_L;
            if nargout ==1
                varargout{1} = HoppingL;
            end
        end
        function varargout = TBSK_Var_gen_single(L_1,L_2,m_1,m_2,nn_level,l,m,n,options)
            arguments
                L_1 double{mustBeInteger};
                L_2 double{mustBeInteger};
                m_1 double{mustBeInteger};
                m_2 double{mustBeInteger};
                nn_level double{mustBeInteger} = -1;
                l  = 0;
                m  = 0;
                n  = 0;
                options.overlap logical = false;
                %                 l double{mustBeInRange(l,0,1)} = 0;
                %                 m double{mustBeInRange(m,0,1)} = 0;
                %                 n double{mustBeInRange(n,0,1)} = 0;
            end
            % <https://en.wikipedia.org/wiki/Tight_binding>
            % L_1,L_2 : orbital angular momentum: s 0 p 1 d 2 f 3
            % m_1,m_2 : magnetic quantum number:
            %           - s 0
            %           - p y:-1 z:0 x:1 % because py = i/sart(2)*(Y-1
            %           +Y1), we label this i as a minus number
            %           - d xy:-2 yz:-1 z2:0 xz:1 x2-y2:2
            %           - f
            %           y(3x2-y2):-3 xyz:-2  yz2:-1
            %           z3:0
            %           x(x2-3y2):3  z(x2-y2):2  xz2:1
            if nargin < 6 || (isa(l,'sym') || isa(m,'sym') ||isa(n,'sym')  )
                sym_mode = true;
                if ~isa(l,'sym')
                    syms l real;
                end
                if ~isa(m,'sym')
                    syms m real;
                end
                if ~isa(n,'sym')
                    syms n real;
                end
                
                %syms l m n real;% debug
            else
                sym_mode = false;
            end
            if L_1 > L_2
                if nargout == 1
                    Coeff = HR.TBSK_Var_gen_single(L_2,L_1,m_2,m_1,nn_level,-l,-m,-n,'overlap',options.overlap);
                    varargout{1} = Coeff;
                    return;
                elseif nargout == 2
                    [E_hop,E_int] = ...
                        HR.TBSK_Var_gen_single(L_2,L_1,m_2,m_1,nn_level,-l,-m,-n,'overlap',options.overlap);
                    varargout{1} = E_hop;
                    varargout{2} = E_int;
                    return;
                elseif nargout == 3
                    [Coeff,varargout{2},varargout{3}] = ...
                        HR.TBSK_Var_gen_single(L_2,L_1,m_2,m_1,nn_level,-l,-m,-n,'overlap',options.overlap);
                    varargout{1} = Coeff;
                    return;
                end
                
            end
            if sym_mode
                Coeff = HR.TBSK_Coeff_gen(L_1,L_2,m_1,m_2,l,m,n,'sym_mode',true);
            else
                Coeff = HR.TBSK_Coeff_gen(L_1,L_2,m_1,m_2,l,m,n);
            end
            %
            if nargout == 1 &&  options.overlap
                varargout{1} = Coeff;
                return
            else
                switch L_1
                    case 0
                        Char1 = 's';
                    case 1
                        Char1 = 'p';
                    case 2
                        Char1 = 'd';
                    case 3
                        Char1 = 'f';
                end
                switch L_2
                    case 0
                        Char2 = 's';
                    case 1
                        Char2 = 'p';
                    case 2
                        Char2 = 'd';
                    case 3
                        Char2 = 'f';
                end
                E_hop{1} = ['V',Char1,Char2,'S'] ;
                E_hop{2} = ['V',Char1,Char2,'P'] ;
                E_hop{3} = ['V',Char1,Char2,'D'] ;
                if options.overlap || nargout ~= 1
                    E_integral{1} = ['S',Char1,Char2,'S'] ;
                    E_integral{2} = ['S',Char1,Char2,'P'] ;
                    E_integral{3} = ['S',Char1,Char2,'D'] ;
                else
                    E_integral = [];
                end
                if nn_level > 0
                    E_hop{1} = [E_hop{1},'_',num2str(nn_level)];
                    E_hop{2} = [E_hop{2},'_',num2str(nn_level)];
                    E_hop{3} = [E_hop{3},'_',num2str(nn_level)];
                    if options.overlap || nargout ~= 1
                        E_integral{1} = [E_integral{1},'_',num2str(nn_level)];
                        E_integral{2} = [E_integral{2},'_',num2str(nn_level)];
                        E_integral{3} = [E_integral{3},'_',num2str(nn_level)];
                    end
                end
            end
            if nargout == 3
                varargout{1} = Coeff;
                varargout{2} = E_hop;
                varargout{3} = E_integral;
                return
            elseif nargout == 2
                E_hop_sym = [sym(E_hop{1},'real'),...
                    sym(E_hop{2},'real'),...
                    sym(E_hop{3},'real') ]*Coeff.';
                E_integral_sym = [sym(E_integral{1},'real'),...
                    sym(E_integral{2},'real'),...
                    sym(E_integral{3},'real') ]*Coeff.';
                varargout{1} = E_hop_sym;
                varargout{2} = E_integral_sym;
            elseif nargout == 1 && ~options.overlap
                E_hop_sym = [sym(E_hop{1},'real'),...
                    sym(E_hop{2},'real'),...
                    sym(E_hop{3},'real') ]*Coeff.';
                varargout{1} = E_hop_sym;
            end
            
        end
        function Coeff = TBSK_Coeff_gen(L_1,L_2,m_1,m_2,l,m,n,options)
            arguments
                L_1 double{mustBeInteger};
                L_2 double{mustBeInteger};
                m_1 double{mustBeInteger};
                m_2 double{mustBeInteger};
                l  = 0;
                m  = 0;
                n  = 0;
                options.sym_mode = false;
            end
            if options.sym_mode
                Coeff = sym(zeros(1,3));
            else
                Coeff = zeros(1,3);
            end
            switch L_1
                case 0
                    switch L_2
                        case 0
                            Coeff(1) = 1;
                        case 1 % s - p
                            switch m_2
                                case -1
                                    Coeff(1) = m;
                                case 0
                                    Coeff(1) = n;
                                case 1
                                    Coeff(1) = l;
                            end
                        case 2 % s - d
                            SQRT3 = sqrt(3) ;
                            switch m_2
                                case -2
                                    Coeff(1) = SQRT3*l*m;
                                case -1
                                    Coeff(1) = SQRT3*m*n;
                                case 0
                                    Coeff(1) = (n^2 - (l^2 + m^2)/2);
                                case 1
                                    Coeff(1) = SQRT3*n*l;
                                case 2
                                    Coeff(1) = SQRT3/2 * (l^2- m^2);
                            end
                        case 3
                            %
                    end
                case 1
                    switch L_2
                        case 1 % p - p
                            if m_1 == m_2
                                switch m_1
                                    case -1
                                        Coeff(1) = m^2 ; Coeff(2) = 1 - Coeff(1);
                                    case 0
                                        Coeff(1) = n^2 ; Coeff(2) = 1 - Coeff(1);
                                    case 1
                                        Coeff(1) = l^2 ; Coeff(2) = 1 - Coeff(1);
                                end
                            else
                                switch m_1
                                    case -1
                                        switch m_2
                                            case 0
                                                Coeff(1) = m*n ; Coeff(2) = - Coeff(1) ;
                                            case 1
                                                Coeff(1) = m*l ; Coeff(2) = - Coeff(1);
                                        end
                                    case 0
                                        switch m_2
                                            case 1
                                                Coeff(1) = n*l ; Coeff(2) = - Coeff(1) ;
                                            case -1
                                                Coeff(1) = n*m ; Coeff(2) = - Coeff(1) ;
                                        end
                                    case 1
                                        switch m_2
                                            case -1
                                                Coeff(1) = l*m ; Coeff(2) = - Coeff(1) ;
                                            case 0
                                                Coeff(1) = l*n ; Coeff(2) = - Coeff(1) ;
                                        end
                                end
                            end
                        case 2 % p - d            - d xy:-2 yz:-1 z2:0 xz:1 x2-y2:2
                            SQRT3 = sqrt(3) ;
                            switch m_1
                                case -1 % py
                                    switch m_2
                                        case -2
                                            lM2 = m^2*l;
                                            Coeff(1) = SQRT3*lM2 ; Coeff(2) = l-2*lM2;
                                        case -1
                                            nM2 = n*m^2 ;
                                            Coeff(1) = SQRT3*nM2 ; Coeff(2) = n-2*nM2;
                                        case 1
                                            LMN = l*m*n ;
                                            Coeff(1) = SQRT3*LMN ; Coeff(2) = -2*LMN;
                                        case 0
                                            mN2 = m*n^2 ;
                                            Coeff(1) = mN2-m^3/2-m*l^2/2 ; Coeff(2) = -SQRT3*mN2;
                                        case 2
                                            M3 = m^3;mL2 = l*m^2;
                                            Coeff(1) = SQRT3/2 * (mL2 - M3);Coeff(2) = -m + M3+mL2;
                                    end
                                case 0  % pz
                                    switch m_2
                                        case -2
                                            LMN = l*m*n ;
                                            Coeff(1) = SQRT3*LMN ; Coeff(2) = -2*LMN;
                                        case -1
                                            mN2 = m*n^2 ;
                                            Coeff(1) = SQRT3*mN2 ; Coeff(2) = m-mN2;
                                        case 1
                                            lN2 = l*n^2 ;
                                            Coeff(1) = SQRT3*lN2 ; Coeff(2) = l-2*lN2;
                                        case 0
                                            nL2 = n*l^2;nM2 = n*m^2;
                                            Coeff(1) = n^3 - nL2/2 - nM2/2 ; Coeff(2) = SQRT3*(nM2+nL2);
                                        case 2
                                            nL2 = n*l^2;nM2 = n*m^2;
                                            Coeff(1) = SQRT3/2 * (nL2 - nM2);Coeff(2) = nM2 - nL2;
                                    end
                                case 1  % px
                                    switch m_2
                                        case -2
                                            mL2 = l^2*m;
                                            Coeff(1) = SQRT3*mL2 ; Coeff(2) = m-2*mL2;
                                        case -1
                                            LMN = l*m*n ;
                                            Coeff(1) = SQRT3*LMN ; Coeff(2) = -2*LMN;
                                        case 1
                                            nL2 = n*l^2 ;
                                            Coeff(1) = SQRT3*nL2 ; Coeff(2) = n-2*nL2;
                                        case 0
                                            lN2 = l*n^2 ;
                                            Coeff(1) = lN2-l^3/2-l*m^2/2 ; Coeff(2) = -SQRT3*lN2;
                                        case 2
                                            L3 = l^3;lM2 = l*m^2;
                                            Coeff(1) = SQRT3/2 * (L3 - lM2);Coeff(2) = l-L3+lM2;
                                    end
                            end
                        case 3
                            %
                    end
                case 2
                    switch L_2
                        case 2
                            % -2: l*m -1:m*n 1:n*l
                            switch m_1
                                case {-2,-1,1}
                                    switch m_2
                                        case {-2,-1,1}
                                            L2 = l^2;M2 = m^2;N2 = n^2;
                                            if m_1 == -2 && m_2 == -2
                                                L2M2 = L2*M2;
                                                Coeff(1) = 3*L2M2 ;
                                                Coeff(2) = L2 + M2 - 4 * L2M2;
                                                Coeff(3) = N2 + L2M2 ;
                                            elseif m_1 == -1 && m_2 == -1
                                                M2N2 = M2*N2 ;
                                                Coeff(1) = 3*M2N2 ;
                                                Coeff(2) = M2 + N2 - 4 * M2N2;
                                                Coeff(3) = L2 + M2N2;
                                            elseif m_1 == 1 && m_2 == 1
                                                N2L2 = N2*L2;
                                                Coeff(1) = 3*N2L2 ;
                                                Coeff(2) = N2 + L2 - 4 * N2L2;
                                                Coeff(3) = M2 + N2L2;
                                            elseif m_1 == -2 && m_2 == -1
                                                LN = l*n;LNM2 = LN*M2;
                                                Coeff(1) = 3*LNM2 ;
                                                Coeff(2) = LN*(1-4*M2);
                                                Coeff(3) = LNM2 - LN;
                                            elseif m_1 == -2 && m_2 == 1
                                                MN = m*n;MNL2 = MN*L2;
                                                Coeff(1) = 3*MNL2 ;
                                                Coeff(2) = MN*(1-4*L2);
                                                Coeff(3) = MNL2 - MN;
                                            elseif m_1 == -1 && m_2 == 1
                                                LM = l*m;LMN2 = LM*N2;
                                                Coeff(1) = 3*LMN2 ;
                                                Coeff(2) = LM*(1-4*N2);
                                                Coeff(3) = LMN2 - LM;
                                            elseif m_1 == -1 && m_2 == -2
                                                LN = l*n;LNM2 = LN*M2; % same with  m_1 == -2 && m_2 == -1
                                                Coeff(1) = 3*LNM2 ;
                                                Coeff(2) = LN*(1-4*M2);
                                                Coeff(3) = LNM2 - LN;
                                            elseif m_1 == 1 && m_2 == -2
                                                MN = m*n;MNL2 = MN*L2; % same with  m_1 == -2 && m_2 == 1
                                                Coeff(1) = 3*MNL2 ;
                                                Coeff(2) = MN*(1-4*L2);
                                                Coeff(3) = MNL2 - MN;
                                            elseif m_1 == 1 && m_2 == -1
                                                LM = l*m;LMN2 = LM*N2; % m_1 == -1 && m_2 == 1
                                                Coeff(1) = 3*LMN2 ;
                                                Coeff(2) = LM*(1-4*N2);
                                                Coeff(3) = LMN2 - LM;
                                            else
                                                error('unexpected error!');
                                            end
                                        case {0,2}
                                            if m_1 == -2 && m_2 ==0
                                                SQRT3 =sqrt(3); LM =l*m;L2 = l^2;M2 = m^2;N2 = n^2;
                                                L2pM2 = L2+M2;
                                                Coeff(1) = SQRT3/2 *LM*(N2-L2pM2);
                                                Coeff(2) = -2*LM*N2 ;
                                                Coeff(3) = LM*(1+N2)/2;
                                            elseif m_1 == -1 && m_2 ==0
                                                SQRT3 =sqrt(3); MN =m*n;L2 = l^2;M2 = m^2;N2 = n^2;
                                                L2pM2 = L2+M2;
                                                Coeff(1) = SQRT3/2 *MN*(N2-L2pM2);
                                                Coeff(2) = MN*(L2pM2-N2) ;
                                                Coeff(3) = -MN*(1+L2pM2)/2;
                                            elseif m_1 == 1 && m_2 ==0
                                                SQRT3 =sqrt(3); NL =n*l;L2 = l^2;M2 = m^2;N2 = n^2;
                                                L2pM2 = L2+M2;
                                                Coeff(1) = SQRT3/2 *NL*(N2-L2pM2);
                                                Coeff(2) = NL*(L2pM2-N2) ;
                                                Coeff(3) = NL*(1+L2pM2)/2;
                                            elseif m_1 == -2 && m_2 ==2
                                                LM =l*m;L2 = l^2;M2 = m^2;N2 = n^2;
                                                L2mM2 = L2-M2;
                                                Coeff(1) = 1.5*LM*L2mM2;
                                                Coeff(2) = -2*LM*L2mM2 ;
                                                Coeff(3) = LM*L2mM2/2;
                                            elseif m_1 == -1 && m_2 ==2
                                                MN =m*n;L2 = l^2;M2 = m^2;N2 = n^2;
                                                L2mM2 = L2-M2;
                                                Coeff(1) = 1.5*MN*L2mM2;
                                                Coeff(2) = -2*MN*L2mM2 - MN ;
                                                Coeff(3) = MN*L2mM2/2 + MN;
                                            elseif m_1 == 1 && m_2 ==2
                                                NL =l*m;L2 = l^2;M2 = m^2;N2 = n^2;
                                                L2mM2 = L2-M2;
                                                Coeff(1) = 1.5*NL*L2mM2;
                                                Coeff(2) = -2*NL*L2mM2 + NL ;
                                                Coeff(3) = NL*L2mM2/2 - NL;
                                            end
                                    end
                                case {0,2}
                                    switch m_2
                                        case {0,2}
                                            L2 = l^2;M2 = m^2;N2 = n^2;
                                            if m_1 == 2 && m_2 == 2
                                                L2mM2_2 = (L2-M2)^2;
                                                Coeff(1) = 0.75*L2mM2_2;
                                                Coeff(2) = L2 + M2 - L2mM2_2;
                                                Coeff(3) = N2 + L2mM2_2/4;
                                            elseif m_1 == 0 && m_2 == 0
                                                L2pM2 = (L2+M2);
                                                Coeff(1) = (N2 - L2pM2/2)^2;
                                                Coeff(2) = 3*N2*L2pM2;
                                                Coeff(3) = 0.75*L2pM2^2;
                                            else
                                                SQRT3 = sqrt(3); L2mM2 = (L2-M2); L2pM2 = (L2+M2);
                                                Coeff(1) = SQRT3/2 * L2mM2*(N2-L2pM2/2 );
                                                Coeff(2) = -N2*L2mM2;
                                                Coeff(3) = (1+N2)*L2mM2/4;
                                            end
                                        case {-2,-1,1}
                                            if m_1 == 0 && m_2 ==-2
                                                SQRT3 =sqrt(3); LM =l*m;L2 = l^2;M2 = m^2;N2 = n^2;
                                                L2pM2 = L2+M2;
                                                Coeff(1) = SQRT3/2 *LM*(N2-L2pM2);
                                                Coeff(2) = -2*LM*N2 ;
                                                Coeff(3) = LM*(1+N2)/2;
                                            elseif m_1 == 0 && m_2 ==-1
                                                SQRT3 =sqrt(3); MN =m*n;L2 = l^2;M2 = m^2;N2 = n^2;
                                                L2pM2 = L2+M2;
                                                Coeff(1) = SQRT3/2 *MN*(N2-L2pM2);
                                                Coeff(2) = MN*(L2pM2-N2) ;
                                                Coeff(3) = -MN*(1+L2pM2)/2;
                                            elseif m_1 == 0 && m_2 ==1
                                                SQRT3 =sqrt(3); NL =n*l;L2 = l^2;M2 = m^2;N2 = n^2;
                                                L2pM2 = L2+M2;
                                                Coeff(1) = SQRT3/2 *NL*(N2-L2pM2);
                                                Coeff(2) = NL*(L2pM2-N2) ;
                                                Coeff(3) = NL*(1+L2pM2)/2;
                                            elseif m_1 == 2 && m_2 ==-2
                                                LM =l*m;L2 = l^2;M2 = m^2;N2 = n^2;
                                                L2mM2 = L2-M2;
                                                Coeff(1) = 1.5*LM*L2mM2;
                                                Coeff(2) = -2*LM*L2mM2 ;
                                                Coeff(3) = LM*L2mM2/2;
                                            elseif m_1 == 2 && m_2 ==-1
                                                MN =m*n;L2 = l^2;M2 = m^2;N2 = n^2;
                                                L2mM2 = L2-M2;
                                                Coeff(1) = 1.5*MN*L2mM2;
                                                Coeff(2) = -2*MN*L2mM2 - MN ;
                                                Coeff(3) = MN*L2mM2/2 + MN;
                                            elseif m_1 == 2 && m_2 ==1
                                                NL =l*m;L2 = l^2;M2 = m^2;N2 = n^2;
                                                L2mM2 = L2-M2;
                                                Coeff(1) = 1.5*NL*L2mM2;
                                                Coeff(2) = -2*NL*L2mM2 + NL ;
                                                Coeff(3) = NL*L2mM2/2 - NL;
                                            end
                                    end
                            end
                        case 3
                            %
                    end
                case 3
                    switch L_2
                        case 3
                    end
                    % not support
            end
        end
    end
    %% Debug method
    methods
        function A = GetProperty(H_hr,name)
            A = H_hr.(name);
        end
    end
    methods (Static)
        EQ_list = Test_TBSK_Var_gen(testmode)
    end
    %% protected method
    % ----------------  old function remove later --------------------
    methods(Hidden,Access= protected)
        function H_hr = nn_sk_smart(H_hr,search_range,Accuracy,Rlength_cut)
            %%
            % Caculate the nn_SK for a primitive cell
            % input : Rm sites
            % output : Atom_store nn_store,Rnn
            % usage : [Atom_store,nn_store,Rnn]=nn(Rm,sites)
            %         [Atom_store,nn_store,Rnn]=nn(Rm,sites,search_range)
            % note :
            % Rm
            %a1 a1x a1y a1z
            %a2 a2x a2y a2z
            %a3 a3x a3y a3z
            %     a1=Rm(1,:);
            %     a2=Rm(2,:);
            % note : site=struct('seq',[],'inseq',[],'rc1',[],'rc2',[],'rc3',[],'name',[],'nameseq',[]);
            %       Atom=struct('totseq',[],'elementseq',[],'Re',[],'seq_in',[],'Rc',[],'name',[]);
            %       nn_t=struct('totseq',[],'orbit_on',[],'orbit_in',[],'Rlength',[],'Rc',[],'name',[],'nn_level',[]);
            %       Rnn sort the Ri -> ti in a primitivecell
            %       the Accuracy takes matter
            
            %% -------- nargin --------
            if nargin <4
                Rlength_cut = 15;
            end
            if nargin <3
                Accuracy = 4;
            end
            if nargin <2
                search_range = [0 0 0];
            end
            %% -------- init --------
            sites_num = size(H_hr.sites,2);
            search_rangex=search_range(1);
            search_rangey=search_range(2);
            search_rangez=search_range(3);
            Atom_smart_t = struct('R_fractional_from',[],...
                'R_fractional_to',[],'R_fractional_diff',[],...
                'seq_from',[],'seq_to',[],...
                'orb_sym_from',[],'orb_sym_to',[],...
                'l_name_from',[],'l_name_to',[],...
                'handyname',[]);
            nn_smart_t = struct('seq_from',[],'seq_to',[],'nn',[]);
            H_hr.Atom_store_smart = repmat(Atom_smart_t,[sites_num sites_num]);
            H_hr.nn_store_smart = repmat(nn_smart_t,[sites_num sites_num]);
            Rnn_list = [];
            %% -------- save ------------
            pb = vasplib_tool_outer.CmdLineProgressBar('Seaching ');
            for j  = 1:sites_num
                site2 = H_hr.sites(j); % homecell
                for i = 1:sites_num
                    site1 = H_hr.sites(i);
                    H_hr.Atom_store_smart(i,j) = HR.Atom_smart_t_gen(site1,site2);
                    [Rnn_list_temp,H_hr.nn_store_smart(i,j)] = ...
                        HR.nn_smart_t_gen(H_hr.Atom_store_smart(i,j),H_hr.Rm,search_rangex,search_rangey,search_rangez,Accuracy,Rlength_cut);
                    Rnn_list = [Rnn_list;Rnn_list_temp];
                end
                pb.print(j,sites_num,' th orb nn information ...');
            end
            pb.delete();
            %% -------- calculate ------------
            H_hr.Rnn = sort(unique(Rnn_list,'row'));
            %% -------- make map ------------
            H_hr.Rnn_map = containers.Map('KeyType','double','ValueType','double');
            for i = 1:length(H_hr.Rnn)
                H_hr.Rnn_map(H_hr.Rnn(i)) = i;
            end
            %%  give level
            pb = vasplib_tool_outer.CmdLineProgressBar('Giving ');
            for j  = 1:sites_num
                pb.print(j,sites_num,' th orb nn_level ...');
                for i = 1:sites_num
                    for k = 1:length(H_hr.nn_store_smart(i,j).nn)
                        nn_level = H_hr.Rnn_map(H_hr.nn_store_smart(i,j).nn(k).Rlength);
                        H_hr.nn_store_smart(i,j).nn(k).nn_level = nn_level ;
                        H_hr.nn_store_smart(i,j).nn(k).hop = HR.hop_gen(H_hr.nn_store_smart(i,j).nn(k).hop_pre,nn_level);
                        if H_hr.overlap
                            H_hr.nn_store_smart(i,j).nn(k).overlap = sym(strrep(string(H_hr.nn_store_smart(i,j).nn(k).hop),'V','S'));
                        end
                    end
                end
                
            end
            pb.delete();
            
        end
        function H_hr = nn_sk_sparse(H_hr,search_range,Accuracy,Rlength_cut)
            MAX_NN_LENGTH = 5000*5000*9;
            %% -------- nargin --------
            if nargin <4
                Rlength_cut = 5;
            end
            if nargin <3
                Accuracy = 4;
            end
            if nargin <2
                search_range = [0 0 0];
            end
            %% -------- init --------
            sites_num = size(H_hr.sites,2);
            search_rangex=search_range(1);
            search_rangey=search_range(2);
            search_rangez=search_range(3);
            % for sparse mode, space is the only requirement
            max_nn_length = min(H_hr.WAN_NUM*H_hr.WAN_NUM*H_hr.NRPTS,MAX_NN_LENGTH);
            nn_sparse = zeros(max_nn_length,10);
            % i j r1 r2 r3 Rlength nn_level coff1(sigma) coff2(pi) coff3(delta)
            Rnn_list = [];
            count  = 1;
            %% -------- save ------------
            for j  = 1:sites_num
                site2 = H_hr.sites(j); % homecell
                for i = 1:sites_num
                    site1 = H_hr.sites(i);
                    [nn_sparse_temp,Rnn_list_temp] = ...
                        HR.nn_sparse_t_gen(site1,site2,H_hr.Rm,search_rangex,search_rangey,search_rangez,Accuracy,Rlength_cut);
                    countadd = size(nn_sparse_temp,1);
                    nn_sparse_temp(:,1) = i;
                    nn_sparse_temp(:,2) = j;
                    nn_sparse(count:count+countadd-1,:) = nn_sparse_temp;
                    count = count+countadd;
                    Rnn_list = [Rnn_list;Rnn_list_temp];
                end
                fprintf('The (%4d/%4d) th orb nn information has been generated.\n',j,sites_num);
            end
            %% -------- clean --------------
            nn_sparse(count:max_nn_length,:) = [];
            %% -------- calculate ------------
            H_hr.Rnn = sort(unique(Rnn_list,'row'));
            %% -------- make map ------------
            H_hr.Rnn_map = containers.Map('KeyType','double','ValueType','double');
            for i = 1:length(H_hr.Rnn)
                H_hr.Rnn_map(H_hr.Rnn(i)) = i;
            end
            %%  give level
            for i  = 1:count-1
                fprintf('Giving (%4d/%4d) th hopping  ... \n',i,count);
                nn_sparse(i,7) = H_hr.Rnn_map(nn_sparse(i,6));
            end
            H_hr.nn_store_smart = nn_sparse;
        end
        function H_hr = H_TB_gen_SK(H_hr,options)
            arguments
                H_hr HR;
                options.level_cut {mustBeInteger} = 1;
                options.onsite logical = false;
                options.per_dir double = [1,1,1];
                options.chiral logical = false;
                options.spin logical = false;
                options.method  = 'nn_sk_smart';
                options.rough = false;
                options.vectorL = [];
            end
            if isempty(H_hr.nn_store_smart)
                if isempty(H_hr.nn_store_smart)
                    error('You have not run the function: nn_sk_smart');
                else
                    if strcmp(options.method,'nn_sparse')
                        error('You should use the key value: ''method'',''nn_smart'' as input.');
                    end
                end
            end
            % init
            N_orbit = H_hr.WAN_NUM;
            N_tot = H_hr.WAN_NUM;
            %orb_list(N_orbit,:) = [0  0  0];
            H_hr.HcoeL(:,:,1)  = sym(zeros(N_orbit));
            H_hr.HnumL(:,:,1)  = zeros(N_orbit);
            level_cut = options.level_cut;
            per_dir = options.per_dir;
            %   on-site
            if options.onsite
                for i=1:N_orbit
                    onsite_sym_name = "E__"+string(H_hr.elementL(i,1))+"_"...
                        +string(H_hr.quantumL(i,2));...
                        H_hr = H_hr.set_hop_single(sym(onsite_sym_name,'real'),i,i,[0 0 0],'sym');
                end
            end
            %             if onsite_mode >0
            %                 for j=1:N_orbit
            %                     onsite_sym_name = "E_onsite_"+string(H_hr.quantumL(j,1))+"_"...
            %                         +string(H_hr.quantumL(j,2));...
            %                         %                         +"_"+string(H_hr.quantumL(j,3))...
            %
            %                     EVALstring="syms "+ onsite_sym_name ...
            %                         +" real;";
            %                     eval(EVALstring);
            %                     % enforce onsite = 0
            %                     %H_xyz = set_hop(-sym(H(j,j)),j,j,[0 0 0],H_xyz,'sym');
            %                     H_hr = H_hr.set_hop(str2sym(onsite_sym_name),j,j,[0 0 0],'sym');
            %                 end
            %             end
            %   hop
            pb = vasplib_tool_outer.CmdLineProgressBar('Setting ');
            for j=1:N_orbit
                pb.print(j,N_orbit,' th orbital ...');
                spin1 = H_hr.quantumL(j,4);
                element1 =  H_hr.elementL(j);
                for i=1:N_tot
                    spin2 = H_hr.quantumL(i,4);
                    element2 =  H_hr.elementL(i);
                    set_or_not = true;
                    if options.chiral
                        if element1 == element2
                            set_or_not = false;
                        end
                    end
                    if options.spin
                        if spin1 == spin2
                            set_or_not = false;
                        end
                    end
                    if set_or_not
                        nn = H_hr.nn_store_smart(i,j).nn;
                        for k = 1:size(nn,1)
                            if nn(k).nn_level <= level_cut
                                vector_init = nn(k).R_vector.*per_dir;
                                H_hr = H_hr.set_hop(nn(k).hop,i,j,vector_init,'sym');
                                if H_hr.overlap
                                    H_hr = H_hr.set_hop(nn(k).overlap,i,j,vector_init,'sym');
                                end
                            end
                        end
                    end
                end
            end
            pb.delete();
            H_hr.coe = true;
            H_hr.num = false;
        end
        function H_hr = H_TB_gen_SK_sparse(H_hr,level_cut,para_filename,onsite_mode)
            % input: level-cut  (nn hopping max cut)
            % warning : The POSCAR must be properly set, no extra atoms!!!
            if nargin <4
                onsite_mode = 0;
            end
            if nargin <3
                para_filename = 'para.mat';
            end
            if nargin <2
                level_cut = 1;
            end
            % load first
            if exist(para_filename,'file')
                load(para_filename);
            else
                error('in sparse mode, we give value for the parameter first!');
            end
            %% init
            N_orbit = H_hr.WAN_NUM;
            %N_tot = H_hr.WAN_NUM;
            % vectorL
            %orb_list(N_orbit,:) = [0  0  0];
            % uodata nn_sparse
            nn_sparse_num = zeros(size(H_hr.nn_store_smart,1),5);
            H_hr.vectorL = unique(H_hr.nn_store_smart(:,3:5),'rows');
            Nhopping = size(H_hr.nn_store_smart,1);
            for i = 1: Nhopping
                nn_sparse_num(i,1) = H_hr.nn_store_smart(i,1);
                nn_sparse_num(i,2) = H_hr.nn_store_smart(i,2);
                [~,nn_sparse_num(i,3)] = ismember(...
                    [H_hr.nn_store_smart(i,3),...
                    H_hr.nn_store_smart(i,4),...
                    H_hr.nn_store_smart(i,5)],...
                    H_hr.vectorL,'rows');
                % hopping
                orb1  = H_hr.sites(nn_sparse_num(i,1)).orb;
                orb2  = H_hr.sites(nn_sparse_num(i,2)).orb;
                if strcmp(orb1,'p')
                    if strcmp(orb2,'s')
                        orb1 = 's';
                        orb2 = 'p';
                    end
                end
                nn_sparse_num(i,4) = H_hr.nn_store_smart(i,7);
                if nn_sparse_num(i,4)  <= level_cut
                    nn_sparse_num(i,5)=nn_sparse_num(i,5)+...
                        H_hr.nn_store_smart(i,8)*subs(str2sym("V"+orb1+orb2+'S_'+string(H_hr.nn_store_smart(i,7))))+...
                        H_hr.nn_store_smart(i,9)*subs(str2sym("V"+orb1+orb2+'P_'+string(H_hr.nn_store_smart(i,7))));...
                        %                     H_hr.nn_store_smart(i,10)*subs(str2sym("V"+orb1+orb2+'D_'+string(H_hr.nn_store_smart(i,7))))+...
                else
                    nn_sparse_num(i,5)=0;
                end
                
            end
            H_hr.HnumL = [];
            H_hr.HnumL{H_hr.NRPTS}  = sparse(N_orbit,N_orbit);
            for i =1:H_hr.NRPTS -1
                H_hr.HnumL{i}  = sparse(N_orbit,N_orbit);
            end
            %   on-site
            if onsite_mode == 1
                for j=1:N_orbit
                    onsite_sym_name = "E_onsite_"+string(H_hr.quantumL(j,1))+"_"...
                        +string(H_hr.quantumL(j,2));...
                        %                         +"_"+string(H_hr.quantumL(j,3))...
                    %                     EVALstring="syms "+ onsite_sym_name ...
                    %                         +" real;";
                    %                     eval(EVALstring);
                    % enforce onsite = 0
                    %H_xyz = set_hop(-sym(H(j,j)),j,j,[0 0 0],H_xyz,'sym');
                    H_hr.HnumL{i}(j,j) = subs(str2sym(onsite_sym_name));
                end
            end
            %   hop
            for Nh=1: Nhopping
                fprintf("setting (%4d/%4d) th hopping ... \n",Nh, Nhopping);
                if nn_sparse_num(Nh,4)  <= level_cut
                    H_hr.HnumL{nn_sparse_num(Nh,3)}(nn_sparse_num(Nh,1),nn_sparse_num(Nh,2)) = ...
                        nn_sparse_num(Nh,5);
                end
            end
            H_hr.nn_sparse_n = nn_sparse_num;
        end
        % ------ old function rm later
        function Atom_smart_t = Atom_smart_t_gen(site1,site2) % othercell -> homecell
            Rc1 = [site1.rc1,site1.rc2,site1.rc3];
            Rc2 = [site2.rc1,site2.rc2,site2.rc3];
            Atom_smart_t.R_fractional_from = Rc1 ;
            Atom_smart_t.R_fractional_to   = Rc2; %fix
            Atom_smart_t.R_fractional_diff = -(Rc1 - Rc2);
            Atom_smart_t.seq_from = site1.seq;
            Atom_smart_t.seq_to = site2.seq;
            Atom_smart_t.l_name_from = site1.orb  ;
            Atom_smart_t.l_name_to   = site2.orb  ;
            Atom_smart_t.orb_sym_from = site1.orb_sym;
            Atom_smart_t.orb_sym_to = site2.orb_sym;
            Atom_smart_t.handyname = strcat(site1.name,' -> ',site2.name);
        end
        function [nn_sparse_temp,Rnn_list] = nn_sparse_t_gen(site1,site2,Rm,search_rangex,search_rangey,search_rangez,Accuracy,Rlength_cut)
            Rc1 = [site1.rc1,site1.rc2,site1.rc3];
            Rc2 = [site2.rc1,site2.rc2,site2.rc3];
            R_fractional_diff = -(Rc1 - Rc2);
            %     nn_smart_t.R_cartesian_to = Atom_smart_t.R_fractional_to*Rm;
            %     nn_smart_t.R_cartesian_from = Atom_smart_t.R_fractional_from*Rm;
            count = 1;
            reducible_num=(2*search_rangex+1)*(2*search_rangey+1)*(2*search_rangez+1);
            Rnn_list = zeros(reducible_num,1);
            nn_sparse_temp= zeros(reducible_num,10);
            % nn_t = struct('R_vector',[],'R_fractional_diff',[],'Rlength',[],'nn_level',[],'hop_pre',[],'hop',[]);
            % nn = repmat(nn_t,[reducible_num 1]);
            for Rf_a1=-search_rangex:search_rangex
                for Rf_a2=-search_rangey:search_rangey
                    for Rf_a3=-search_rangez:search_rangez
                        R_vector = [Rf_a1 Rf_a2 Rf_a3];
                        Rij_cart = (R_vector + R_fractional_diff)*Rm ;
                        Rlength = norm(Rij_cart);
                        Rlength = roundn(Rlength,-Accuracy);
                        Rlmn = Rij_cart/Rlength;
                        if  0 < Rlength && Rlength < Rlength_cut
                            orb1 = site1.orb   ;
                            orb2 = site2.orb   ;
                            orb_sym1 =  site1.orb_sym;
                            orb_sym2 =  site2.orb_sym;
                            orbsym1_n =  HR.subs_xyz(orb_sym1 ,Rlmn);
                            orbsym2_n =  HR.subs_xyz(orb_sym2 ,Rlmn);
                            nn_sparse_temp(count,3) = Rf_a1;
                            nn_sparse_temp(count,4) = Rf_a2;
                            nn_sparse_temp(count,5) = Rf_a3;
                            nn_sparse_temp(count,6) =  Rlength;
                            [~,Coff] = HR.TBSK_hop_gen(orb1,orb2,orbsym1_n,orbsym2_n,orb_sym1,orb_sym2);
                            nn_sparse_temp(count,8) =   Coff(1);
                            nn_sparse_temp(count,9) =   Coff(2);
                            nn_sparse_temp(count,10) =  Coff(3);
                            
                            Rnn_list(count,:) = Rlength;
                            count = count +1;
                        end
                    end
                end
            end
            if count <= reducible_num
                Rnn_list(count:reducible_num,:) = [];
                nn_sparse_temp(count:reducible_num,:) = [];
            end
        end
        function [Rnn_list,nn_smart_t] = nn_smart_t_gen(Atom_smart_t,Rm,search_rangex,search_rangey,search_rangez,Accuracy,Rlength_cut)
            nn_smart_t.seq_from = Atom_smart_t.seq_from;
            nn_smart_t.seq_to = Atom_smart_t.seq_to;
            %     nn_smart_t.R_cartesian_to = Atom_smart_t.R_fractional_to*Rm;
            %     nn_smart_t.R_cartesian_from = Atom_smart_t.R_fractional_from*Rm;
            count = 1;
            reducible_num=(2*search_rangex+1)*(2*search_rangey+1)*(2*search_rangez+1);
            Rnn_list = zeros(reducible_num,1);
            nn_t = struct('R_vector',[],'R_fractional_diff',[],'Rlength',[],'nn_level',[],'hop_pre',[],'hop',[]);
            nn = repmat(nn_t,[reducible_num 1]);
            for Rf_a1=-search_rangex:search_rangex
                for Rf_a2=-search_rangey:search_rangey
                    for Rf_a3=-search_rangez:search_rangez
                        R_vector = [Rf_a1 Rf_a2 Rf_a3];
                        Rij_cart = (R_vector + Atom_smart_t.R_fractional_diff)*Rm ;
                        Rlength = norm(Rij_cart);
                        Rlength = roundn(Rlength,-Accuracy);
                        Rlmn = Rij_cart/Rlength;
                        if  0 < Rlength && Rlength < Rlength_cut
                            orb1 = Atom_smart_t.l_name_from ;
                            orb2 = Atom_smart_t.l_name_to   ;
                            orb_sym1 =  Atom_smart_t.orb_sym_from;
                            orb_sym2 =  Atom_smart_t.orb_sym_to;
                            orbsym1_n =  HR.subs_xyz(orb_sym1 ,Rlmn);
                            orbsym2_n =  HR.subs_xyz(orb_sym2 ,Rlmn);
                            TBSK_hop = HR.TBSK_hop_gen(orb1,orb2,orbsym1_n,orbsym2_n,orb_sym1,orb_sym2);
                            nn(count,1).hop_pre = TBSK_hop;
                            nn(count,1).R_vector = R_vector;
                            nn(count,1).R_fractional_diff = Atom_smart_t.R_fractional_diff;
                            nn(count,1).Rlength = Rlength;
                            Rnn_list(count,:) = Rlength;
                            count = count +1;
                        end
                    end
                end
            end
            if count <= reducible_num
                Rnn_list(count:reducible_num,:) = [];
                nn(count:reducible_num,:) = [];
            end
            nn_smart_t.nn = nn;
        end
        function hop = hop_gen(hop_pre,nn_level)
            tempstr = string(simplify(hop_pre));
            tempstr = strcat(tempstr,"_",string(nn_level));
            hop = sym(tempstr,'real');
        end
        function [TBSK_hop,Coff] = TBSK_hop_gen(orb1,orb2,orbsym1_n,orbsym2_n,orbsym1,orbsym2)
            %disp([orbsym1_n,orbsym2_n]);
            
            if strcmp(orb1,'p')
                if strcmp(orb2,'s')
                    orb1 = 's';
                    orb2 = 'p';
                    orbsym1_n = -orbsym1_n;
                end
            end
            Coff(1) = orbsym1_n*orbsym2_n;
            Coff(2) = HR.delta_orb(orb1,orb2)*(HR.delta_orb_sym(orbsym1,orbsym2)-orbsym1_n*orbsym2_n);
            Coff(3) = 0;
            TBSK_hop = Coff(1) *str2sym("V"+orb1+orb2+'S')+...
                Coff(2)*str2sym("V"+orb1+orb2+'P');
            
        end
        function [TBSK_hop,Coff] =  TBSK_hop_gen_sparse(site1,site2,Rij_cart,Rlength,nn_level)
            Rlmn = Rij_cart/Rlength;
            orb1 = site1.orb   ;
            orb2 = site2.orb   ;
            orb_sym1 =  site1.orb_sym;
            orb_sym2 =  site2.orb_sym;
            orbsym1_n =  HR.subs_xyz(orb_sym1 ,Rlmn);
            orbsym2_n =  HR.subs_xyz(orb_sym2 ,Rlmn);
            if strcmp(orb1,'p')
                if strcmp(orb2,'s')
                    orb1 = 's';
                    orb2 = 'p';
                    orbsym1_n = -orbsym1_n;
                end
            end
            Coff(1) = orbsym1_n*orbsym2_n;
            Coff(2) = HR.delta_orb(orb1,orb2)*(HR.delta_orb_sym(orb_sym1,orb_sym2)-orbsym1_n*orbsym2_n);
            if abs(Coff(1)) < 1e-10
                Coff(1) = 0;
            end
            if abs(Coff(2)) < 1e-10
                Coff(2) = 0;
            end
            Coff(3) = 0;
            VP = "V"+orb1+orb2+'S_'+num2str(nn_level);
            VS = "V"+orb1+orb2+'P_'+num2str(nn_level);
            VD = "";
            TBSK_hop =  Coff(1) * sym(VP,'real')+...
                Coff(2) * sym(VS,'real');
        end
    end
end



