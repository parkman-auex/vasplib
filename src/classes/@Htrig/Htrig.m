classdef Htrig < vasplib & matlab.mixin.CustomDisplay
    properties
        HcoeL;
        HnumL;
        HsymL_trig = sym(1);
    end
    properties % (Hidden = true)
        HsymL_coeL ;
        HsymL_numL ;
        Htrig_num   ;
        Nslab = [0,0,0];
        Htrig_latex ;
        HsymL_trig_bk ;
        seeds = [];
        seedsvar = sym([]);
        % for spacegroup
        Sigmas        = []     ;
        rm_list                ;
        Trig_list =    Trig()  ;
        Type {mustBeMember(Type,{'mat','list','sincos','exp','slab','sparse'})}= 'sincos';
        Hmat_pre;
        symvarL;
    end
    properties(Dependent = true)
        Htrig_sym   ;
        HsymL;
        Kinds;
    end
    properties %(GetAccess = protected,Hidden = true)
        num = false;        % num(numeric).
        coe = true;        % coe(symbolic).
        Duality_vector_dist; % the opposite vector sequnce dictionary.
    end
    %% Define which properties show
    methods (Access = protected)
        function propgrp = getPropertyGroups(~)
            proplist = {'Basis_num','Kinds','Type','HsymL','symvar_list','Dim'};
            propgrp = matlab.mixin.util.PropertyGroup(proplist);
        end
    end
    %% constuction
    methods
        function H_htrig = Htrig(BASIS_NUM,Trig_list,options,propArgs)
            arguments
                BASIS_NUM = 4;
                Trig_list =[];
                options.Type = 'sincos';
                propArgs.?vasplib;
            end
            %
            propArgsCell = namedargs2cell(propArgs);                                              
            H_htrig = H_htrig@vasplib(propArgsCell{:});
            VarUsing = H_htrig.VarsSeqLcart(1:H_htrig.Dim);
            H_htrig.seedsvar = VarUsing;
            H_htrig.seeds = string(H_htrig.seedsvar);
            %syms k_x k_y k_z real;
            switch nargin
                case 1
                    if isnumeric(BASIS_NUM)
                        H_htrig.Basis_num = BASIS_NUM;
                        H_htrig.HcoeL = sym(zeros(BASIS_NUM,BASIS_NUM,1));
                        H_htrig.HnumL = (zeros(BASIS_NUM,BASIS_NUM,1));
                        H_htrig.Trig_list = Trig()  ;
                        switch options.Type
                            case 'sincos'
                                H_htrig.HsymL_trig_bk = [sin(VarUsing) cos(VarUsing) ];
                            case 'exp'
                                H_htrig.HsymL_trig_bk = [exp(1i*VarUsing) exp(-1i*VarUsing)];
                            case 'mat'
                                H_htrig.HsymL_coeL = sym(zeros(1,H_htrig.Dim)); % k_x k_y k_z ...
                                H_htrig.HsymL_numL = zeros(1,H_htrig.Dim);
                            case 'list'
                                H_htrig.HsymL_coeL = [zeros(1,H_htrig.Dim,'sym'),ones(1,2,'sym')];
                                H_htrig.HcoeL = sym(0);
                                H_htrig.HsymL_numL = [zeros(1,H_htrig.Dim,'double'),ones(1,2,'double')];
                                H_htrig.HnumL = 0;
                        end
                        H_htrig.Type = options.Type;
                    else
                        if isa(BASIS_NUM,'sym')
                            Hsym = BASIS_NUM;
                            if size(Hsym,1) ~= size(Hsym,2)
                                error('square sym mat required!')
                            end
                            BASIS_NUM = size(Hsym,1);
                            H_htrig.Basis_num = BASIS_NUM;
                            H_htrig.HcoeL = sym(zeros(BASIS_NUM,BASIS_NUM,1));
                            H_htrig.HnumL = (zeros(BASIS_NUM,BASIS_NUM,1));
                            H_htrig.Trig_list = Trig()  ;
                            H_htrig = H_htrig + Hsym;
                        end
                    end
                case 2
                    opitonsCell = namedargs2cell(options);
                    H_htrig = Htrig(BASIS_NUM,opitonsCell{:},propArgsCell{:});
                    H_htrig = H_htrig + Trig_list;
            end
        end
    end
    %% set construction
    methods
        function H_htrig = setup_rough(H_htrig,symbolic_polynomial,pauli_mat,silence)
            if nargin < 4
                silence = true;
            end
            [coeff_trig,symvar_list_trig,H_htrig] = H_htrig.split_sym_eq((simplify(symbolic_polynomial,'IgnoreAnalyticConstraints',true)));
            nc = length(coeff_trig);
            if nc >0
                switch H_htrig.Type
                    case 'mat'
                        for i =1:nc
                            H_htrig = H_htrig.set_hop(coeff_trig(i)*sum(double(pauli_mat),3),symvar_list_trig(i,:));
                        end
                    case 'list'
                        Mtmp  = sum(double(pauli_mat),3);
                        LabelnonZero = find(Mtmp~=0);
                        nLabelnonZero = length(LabelnonZero);
                        [iL,jL] = ind2sub(size(Mtmp),LabelnonZero);
                        Mtmp = Mtmp(LabelnonZero);
                        Mlist = Mtmp(:);
                        Mtmp_list = kron(Mlist,ones(nc,1));
                        iL = kron(iL,ones(nc,1,'sym'));
                        jL = kron(jL,ones(nc,1,'sym'));
                        symvar_list_trig = repmat(symvar_list_trig,[nLabelnonZero,1]);
                        coeff_trig = repmat(coeff_trig(:),[nLabelnonZero,1]);
                        H_htrig = H_htrig.set_hop(coeff_trig.*Mtmp_list,[symvar_list_trig,iL,jL]);
                    otherwise
                        for i =1:nc
                            k_cell{i} = symvar_list_trig(i);
                            mat_cell{i} = sum(double(pauli_mat),3);
                            Var_cell{i} = coeff_trig(i);
                        end
                        H_htrig = setup(H_htrig,Var_cell,k_cell,mat_cell,silence);
                end
            end
        end
        function H_htrig = setup(H_htrig,Var_cell,k_cell,mat_cell,silence)
            if nargin < 5
                silence = false;
            end
            BASIS_NUM = H_htrig.Basis_num;
            if length(Var_cell)~=length(k_cell) && length(k_cell)~=length(mat_cell)
                error('error!');
            end
            if ~silence
                pb = vasplib_tool_outer.CmdLineProgressBar('Setting ');
            end
            nVar_cell = length(Var_cell);
            for i =1:nVar_cell
                Var = Var_cell{i};
                k_symbol = (k_cell{i});
                %disp(k_symbol);
                matcell = mat_cell{i};
                Kind = H_htrig.k_symbol2Kind(k_symbol);
                if isempty(Kind)
                    Kind = H_htrig.Kinds+1;
                    H_htrig.HsymL_trig(Kind) = k_symbol;
                    H_htrig.HcoeL(:,:,Kind) = sym(zeros(BASIS_NUM,BASIS_NUM,1));
                    H_htrig.HnumL(:,:,Kind)  = (zeros(BASIS_NUM,BASIS_NUM,1));
                end
                if ~silence
                    pb.print(i,nVar_cell,'Htrig ...');
                end
                switch class(Var)
                    case 'sym'
                        H_htrig.HcoeL(:,:,Kind) = H_htrig.HcoeL(:,:,Kind)+ matcell*Var;
                    case 'string'
                        
                    case 'double'
                        H_htrig.HnumL(:,:,Kind) = H_htrig.HnumL(:,:,Kind)+ matcell*Var;
                end
                
            end
            if ~silence
                pb.delete();
            end
        end
        function H_htrig = set_hop(H_htrig,SymHopping,SymVar,Indij)
            if isempty(H_htrig.num) || isempty(H_htrig.coe)
                [~,~,H_htrig] = H_htrig.NumOrCoe();
            end
            if nargin < 4
                %[coeff_trig,symvar_list_trig,H_htrig] = split_sym_eq(H_htrig,simplify(SymVar));
                nSymVar = size(SymVar,1);
                if size(SymVar,1) >1 && size(SymVar,2) ==1
                    SymVar = SymVar.';
                end
                if isempty(H_htrig.num) || isempty(H_htrig.coe)
                    [~,~,H_htrig] = H_htrig.NumOrCoe;
                end
                if H_htrig.num
                    SymHopping = double(SymHopping);
                    %SymVar = double(SymHopping)
                else
                    SymHopping = sym(SymHopping);
                    %SymVar = sym(SymHopping);
                end
                if isvector(SymHopping)
                    matmode = false;%listmode = true;
                elseif size(SymHopping,3) == 1 && ismatrix(SymHopping)
                    SymHopping = repmat(SymHopping,[1 1 nSymVar]);
                    matmode = true;%listmode = false;
                else
                    matmode = true;%listmode = false;
                end
                if matmode
                    for i =1:nSymVar
                        Kind = H_htrig.k_symbol2Kind(SymVar(i,:));
                        if isempty(Kind)
                            H_htrig = add_empty_one(H_htrig,SymVar(i));
                            Kind = H_htrig.Kinds;
                        else
                            try
                                if Kind == 0
                                    H_htrig = add_empty_one(H_htrig,SymVar(i));
                                    Kind = H_htrig.Kinds;
                                end
                            catch

                            end
                        end
                        if H_htrig.coe
                            H_htrig.HcoeL(:,:,Kind) = H_htrig.HcoeL(:,:,Kind) + SymHopping(:,:,i);
                        else
                            H_htrig.HnumL(:,:,Kind) = H_htrig.HnumL(:,:,Kind) + SymHopping(:,:,i);
                        end
                    end
                else
                    if H_htrig.coe
                        HsymL_coeLtmp= [H_htrig.HsymL_coeL;SymVar];
                        Hcoetmp = [H_htrig.HcoeL;SymHopping];
                        [H_htrig.HsymL_coeL,H_htrig.HcoeL] = ...
                            HollowKnight.generalcontractrow(HsymL_coeLtmp,Hcoetmp);
                    else
                        HsymL_numLtmp= [H_htrig.HsymL_numL;SymVar];
                        Hnumtmp = [H_htrig.HnumL;SymHopping];
                        [H_htrig.HsymL_numL,H_htrig.HnumL] = ...
                            HollowKnight.generalcontractrow(HsymL_numLtmp,Hnumtmp);
                    end
                end
            else
                % only support mat model now
                [coeff_trig,symvar_list_trig,H_htrig] = split_sym_eq(H_htrig,simplify(SymVar));
                nc = length(coeff_trig);
                for i =1:nc
                    SymVar = symvar_list_trig(i);
                    Kind = H_htrig.k_symbol2Kind(SymVar);
                    if isempty(Kind)
                        H_htrig = add_empty_one(H_htrig,SymVar);
                        Kind = H_htrig.Kinds;
                    end
                    H_htrig.HcoeL(Indij(1),Indij(2),Kind) = H_htrig.HcoeL(Indij(1),Indij(2),Kind) +  SymHopping*coeff_trig(i);
                end              
            end
        end
        function [coeff_trig,symvar_list_trig,H_htrig] = split_sym_eq(H_htrig,symbolic_polynomial)
            % key function!
            switch H_htrig.Type
                case {'exp'}
                    k_symbol = combine(rewrite(symbolic_polynomial,'exp'));
                case {'mat','list'}% for test
                    k_symbol = expand(simplify(rewrite(symbolic_polynomial,'exp'),'IgnoreAnalyticConstraints',true));
                case 'sincos'
                    k_symbol = expand(simplify(rewrite(symbolic_polynomial,'sincos'),'IgnoreAnalyticConstraints',true));
                otherwise
                    k_symbol = expand(rewrite(simplify(symbolic_polynomial),'sincos'));
            end
            k_symbol_str = string(k_symbol);
            symvar_list = symvar(k_symbol);
            % Standerd input
            for i = 1:length(symvar_list)
                str_tmp = string(symvar_list(i));
                switch str_tmp
                    case {'k_x','k_X','K_X','K_x','kx','kX','KX','Kx'}
                        k_symbol_str = strrep(k_symbol_str,str_tmp,'k_x');
                    case {'k_y','k_Y','K_y','K_Y','ky','kY','Ky','KY'}
                        k_symbol_str = strrep(k_symbol_str,str_tmp,'k_y');
                    case {'k_z','k_Z','K_z','K_Z','kz','kZ','Kz','KZ'}
                        k_symbol_str = strrep(k_symbol_str,str_tmp,'k_z');
                    case {'k_w','k_W','K_w','K_W','kw','kW','Kw','KW'}
                        k_symbol_str = strrep(k_symbol_str,str_tmp,'k_w');
                end
            end
            % select
            switch H_htrig.Type
                case {'exp','mat','list'}
                    k_symbol = combine(str2sym(k_symbol_str),'exp');
                otherwise
                    k_symbol = expand(combine(str2sym(k_symbol_str),'sincos'));
            end
            % use children
            % an unexpected bug; when only one term;
            % check if only oneterm
            k_symbol_children1 = children(k_symbol);
            if isequal(simplify(fold(@mtimes,[k_symbol_children1{:}])),k_symbol)
                % only one term
                k_symbol_children{1} = k_symbol;
            elseif isequal(simplify(fold(@plus,[k_symbol_children1{:}])),k_symbol)
                 k_symbol_children = k_symbol_children1;
            else
                k_symbol_children{1} = k_symbol;
            end
            %
            switch H_htrig.Type
                case {'mat','list'}
                    % check son first
                    %
                    nSon = length(k_symbol_children);
                    symvar_list_trig = zeros(nSon,3,'sym');
                    coeff_trig = zeros(nSon,1,'sym');
                    % fix bug
                    if k_symbol == sym(0)

                    else
                        for k = 1:nSon
                            [TmpCoe,TmpVar] = coeffs(k_symbol_children{k});
                            ActualVar = simplify(log(TmpVar),'IgnoreAnalyticConstraints',true)/1i;
                            [ActualCoe,ktype] = coeffs(ActualVar,H_htrig.seedsvar);
                            if ismember(sym(1),ktype)
                                n1 = find(sym(1)==ktype);
                                TmpCoe = simplify(exp(ActualCoe(n1)*(1i)),'IgnoreAnalyticConstraints',true)*TmpCoe;
                                ActualCoe(n1) = [];
                                ktype(n1) = [];
                            end
                            %[ActualCoe,ktype] = coeffs(ActualVar,H_htrig.seedsvar);
                            % temp support k_x k_y k_z
                            HsymL_coeLtmp = zeros(1,3,'sym');
                            for i = 1:3
                                ik = find(H_htrig.seedsvar(i) == ktype);
                                if ~isempty(ik)
                                    HsymL_coeLtmp(i) = ActualCoe(ik);
                                end
                            end
                            symvar_list_trig(k,:) = HsymL_coeLtmp;
                            coeff_trig(k) = TmpCoe;
                        end
                    end
                    if strcmp(H_htrig.Type,'mat')
                        % list do not need empty one
                        H_htrig = H_htrig.add_empty_one(symvar_list_trig);
                    end
                    return;
                otherwise
                    for k = 1:length(k_symbol_children)
                        try % ugly fix bug
                            [coeff_trig,~] = coeffs(k_symbol_children{k},H_htrig.HsymL_trig_bk);
                        catch
                            coeff_trig = k_symbol_children{k};
                        end
                        for i = 1:numel(coeff_trig)
                            tmp_label = contains(string(coeff_trig),H_htrig.seeds);
                            if sum(tmp_label)
                                [~,coeff_trig_list] = coeffs(coeff_trig(i));
                                for j = 1:numel(coeff_trig_list)
                                    tmp_label2 = contains(string(coeff_trig_list(j)),H_htrig.seeds);
                                    if sum(tmp_label2)
                                        H_htrig = H_htrig.find_HsymL_trig_bk(coeff_trig_list(j));
                                    end
                                end
                            end
                        end
                    end
            end
            coeff_trig = sym([]);
            symvar_list_trig= sym([]);
            for k = 1:length(k_symbol_children)
                [coeff_trig_tmp,symvar_list_trig_tmp] = coeffs(k_symbol_children{k},H_htrig.HsymL_trig_bk);
                coeff_trig = [coeff_trig,coeff_trig_tmp];
                symvar_list_trig = [symvar_list_trig,symvar_list_trig_tmp];
            end
        end
        function H_htrig = find_HsymL_trig_bk(H_htrig,coeff_trig)
            VarUsing = H_htrig.VarsSeqLcart(1:H_htrig.Dim);
            for j =1:length(coeff_trig)
                coeff_trig_str = string(sym('A')*coeff_trig(j));
                if isequal(H_htrig.seeds,string(VarUsing)) && strcmp(H_htrig.Type,'sincos')
                    %pat2 = "*"+("sin"|"cos");
                    pat1 = "*"+("sin"|"cos");
                    % "^"+wildcardPattern+
                    pat2 = "^"+digitsPattern(1) + "*"+("sin"|"cos");
                    % we cant solve the sin/cos(k_x)^n problem
                    %coeff_trig_str_list = split(coeff_trig_str,["*sin","*cos"]);
                    Type = 'sincos';
                elseif sum(contains(H_htrig.seeds,"Sigma_x"|"Sigma_y"|"Sigma_z"|"Sigma_w")) && strcmp(H_htrig.Type,'slab')
                    pat1 = "*"+"Sigma_";
                    
                    Type = 'slab';
                elseif isequal(H_htrig.seeds,string(VarUsing)) && strcmp(H_htrig.Type,'exp')
                    pat1 = "*"+"exp";
                    Type = 'exp';
                else
                    
                end
                if strcmp(Type,'sincos')
                    coeff_trig_str_list = split(coeff_trig_str,[pat1,pat2]);
                    %coeff_trig_str_list = split(coeff_trig_str,["*sin","*cos"]);
                    for i = 2:numel(coeff_trig_str_list)
                        if coeff_trig_str_list(i) ~= ""
                            sym1 = str2sym(strcat('cos',coeff_trig_str_list(i)));
                            if ~ismember(sym1,H_htrig.HsymL_trig_bk)
                                sym2 = str2sym(strcat('sin',coeff_trig_str_list(i)));
                                H_htrig.HsymL_trig_bk = [H_htrig.HsymL_trig_bk,sym1,sym2];
                            end
                        end
                    end
                elseif strcmp(Type,'slab')
                    coeff_trig_str_list = split(coeff_trig_str,[pat1]);
                    %coeff_trig_str_list = split(coeff_trig_str,["*Sigma","*Sigma"]);
                    for i = 2:numel(coeff_trig_str_list)
                        if coeff_trig_str_list(i) ~= ""
                            if ~contains(string(H_htrig.HsymL_trig_bk),coeff_trig_str_list(i))
                                sym1 = str2sym(strcat('Sigma_',coeff_trig_str_list(i)));
                                %                           sym2 = str2sym(strcat('sin',coeff_trig_str_list(i)));
                                H_htrig.HsymL_trig_bk = [H_htrig.HsymL_trig_bk,sym1];
                            end
                        end
                    end
                elseif strcmp(Type,'exp')
                    coeff_trig_str_list = split(coeff_trig_str,[pat1]);
                    %coeff_trig_str_list = split(coeff_trig_str,["*Sigma","*Sigma"]);
                    for i = 2:numel(coeff_trig_str_list)
                        if coeff_trig_str_list(i) ~= ""
                            if ~contains(string(H_htrig.HsymL_trig_bk),coeff_trig_str_list(i))
                                sym1 = str2sym(strcat('exp',coeff_trig_str_list(i)));
                                sym2 = str2sym(strcat('exp(-',coeff_trig_str_list(i),")"));
                                H_htrig.HsymL_trig_bk = [H_htrig.HsymL_trig_bk,sym1,sym2];
                            end
                        end
                    end
                else
                end
            end
            %disp('debug');
        end
        function Kind = k_symbol2Kind(H_htrig,k_symbol)
            switch H_htrig.Type
                case {'exp','sincos'}
                    % or use sym to find
                    str_2_compare = string(k_symbol);
                    Kind = find(string(H_htrig.HsymL_trig) == str_2_compare);
                case {'mat','list'}
                    if isa(k_symbol,'sym')
                        [~,Kind]=ismember(k_symbol,H_htrig.HsymL_coeL,'rows');
                    else
                        [~,Kind]=ismember(k_symbol,H_htrig.HsymL_numL,'rows');
                    end
                otherwise
                    % or use sym to find
                    str_2_compare = string(k_symbol);
                    Kind = find(string(H_htrig.HsymL_trig) == str_2_compare);
            end
        end
        function H_htrig = reseq(H_htrig,wan_list,kinds_list)
            if nargin < 3
                kinds_list = ':';
            end
            % wan first
            if ~isequal(wan_list,':')
                %                 H_htrig = H_htrig.set_WAN_NUM(length(wan_list));
                switch H_htrig.Type
                    case 'sparse'
                        for i = 1:H_htrig.Kinds
                            H_htrig.HnumL{i}=H_htrig.HnumL{i}(wan_list,wan_list);
                        end
                    case 'list'
                        wan_list =  all(ismember(H_htrig.vectorL(:,H_htrig.Dim+1:H_htrig.Dim+2),int32(wan_list)),2);
                        if H_htrig.num
                            H_htrig.HnumL=H_htrig.HnumL(wan_list,:);
                            H_htrig.HsymL_numL=H_htrig.H_htrig.HsymL_numL(wan_list,:);
                            H_htrig.HsymL_coeL = [];
                            H_htrig.HcoeL= [];
                        end
                        if H_htrig.coe
                            H_htrig.HcoeL=H_htrig.HcoeL(wan_list,:);
                            H_htrig.HsymL_coeL=H_htrig.H_htrig.HsymL_coeL(wan_list,:);
                            H_htrig.HsymL_numL = [];
                            H_htrig.HnumL= [];
                        end
                    case 'mat'
                        if H_htrig.num
                            H_htrig.HnumL=H_htrig.HnumL(wan_list,wan_list,:);
                        else
                            H_htrig.HnumL = [];
                            H_htrig.HsymL_numL = [];
                        end
                        if H_htrig.coe
                            H_htrig.HcoeL=H_htrig.HcoeL(wan_list,wan_list,:);
                        else
                            H_htrig.HcoeL= [];
                            H_htrig.HcoeL = sym([]);
                        end
                    otherwise
                       if H_htrig.num
                            H_htrig.HnumL=H_htrig.HnumL(wan_list,wan_list,:);
                        else
                            H_htrig.HnumL = [];
                        end
                        if H_htrig.coe
                            H_htrig.HcoeL=H_htrig.HcoeL(wan_list,wan_list,:);
                        else
                            H_htrig.HcoeL = sym([]);
                        end
                end
                if ~isempty(H_htrig.sites)
                    try
                        H_htrig.sites = H_htrig.sites(wan_list);
                    catch
                        % bug
                    end
                end
                if ~isempty( H_htrig.orbL )
                    H_htrig.orbL = H_htrig.orbL(wan_list,:);
                end
                if ~isempty( H_htrig.elementL )
                    H_htrig.elementL = H_htrig.elementL(wan_list,:);
                end
                if ~isempty( H_htrig.quantumL)
                    H_htrig.quantumL = H_htrig.quantumL(wan_list,:);
                end
                % bug here sym_orbL waiting
            end
            % kinds_list
            if ~isequal(kinds_list,':')
                %                 H_htrig = H_htrig.set_NRPTS(length(nrpt_list));
                switch H_htrig.Type
                    case 'sparse'
                        H_htrig.HnumL=H_htrig.HnumL(kinds_list);
                    case 'list'
                    if H_htrig.num
                        H_htrig.HsymL_numL=H_htrig.HsymL_numL(kinds_list,:);
                        H_htrig.HnumL=H_htrig.HnumL(kinds_list);
                        if ~isempty(H_htrig.HcoeL)
                            try
                                H_htrig.HcoeL=H_htrig.HcoeL(:,:,kinds_list);
                            catch
                            end
                        end
                        if ~isempty(H_htrig.HsymL_coeL)
                            try
                                H_htrig.HsymL_coeL=H_htrig.HsymL_coeL(kinds_list,:);
                            catch
                            end
                        end
                    end
                    if H_htrig.coe
                        H_htrig.HsymL_coeL=H_htrig.HsymL_coeL(kinds_list,:);
                        H_htrig.HcoeL=H_htrig.HcoeL(kinds_list);
                        if ~isempty(H_htrig.HnumL)
                            try
                                H_htrig.HnumL=H_htrig.HnumL(:,:,kinds_list);
                            catch
                            end
                        end
                        if ~isempty(H_htrig.HsymL_numL)
                            try
                                H_htrig.HsymL_numL=H_htrig.HsymL_numL(kinds_list,:);
                            catch
                            end
                        end
                    end
                    case  'mat'
                    if H_htrig.num
                        H_htrig.HsymL_numL=H_htrig.HsymL_numL(kinds_list,:);
                        H_htrig.HnumL=H_htrig.HnumL(:,:,kinds_list);
                        if ~isempty(H_htrig.HcoeL)
                            try
                                H_htrig.HcoeL=H_htrig.HcoeL(:,:,kinds_list);
                            catch
                            end
                        end
                        if ~isempty(H_htrig.HsymL_coeL)
                            try
                                H_htrig.HsymL_coeL=H_htrig.HsymL_coeL(kinds_list,:);
                            catch
                            end
                        end
                    end
                    if H_htrig.coe
                        H_htrig.HsymL_coeL=H_htrig.HsymL_coeL(kinds_list,:);
                        H_htrig.HcoeL=H_htrig.HcoeL(:,:,kinds_list);
                        if ~isempty(H_htrig.HnumL)
                            try
                                H_htrig.HnumL=H_htrig.HnumL(:,:,kinds_list);
                            catch
                            end
                        end
                        if ~isempty(H_htrig.HsymL_numL)
                            try
                                H_htrig.HsymL_numL=H_htrig.HsymL_numL(kinds_list,:);
                            catch
                            end
                        end
                    end
                    otherwise
                        %[num_label,coe_label] = H_htrig.NumOrCoe();
                        H_htrig.HsymL_trig = H_htrig.HsymL_trig(kinds_list);
                        if H_htrig.coe
                            H_htrig.HcoeL = H_htrig.HcoeL(:,:,kinds_list);
                        end
                        if H_htrig.num
                            H_htrig.HnumL = H_htrig.HnumL(:,:,kinds_list);
                        end
                end
            end
        end
        % add a empty one
        function H_htrig = add_empty_one(H_htrig,vector)
            if isempty(H_htrig.coe)||isempty(H_htrig.num)
                [~,~,H_htrig] = H_htrig.NumOrCoe();
            end
            %vectorClass = class(vector);
            switch H_htrig.Type
                case {'mat','list'}
                    if H_htrig.num
                        vector = double(vector);
                        for i = 1:size(vector,1)
                            vector_single = (vector(i,:));
                            try
                                if (ismember(vector_single,H_htrig.HsymL_numL,'rows'))
                                    continue; % 进入下一个
                                end
                            catch
                                % ugly
                            end
                            Kind = H_htrig.Kinds +1;
                            H_htrig.HsymL_numL(Kind,:) = (vector_single);
                            H_htrig.HnumL(:,:,Kind) = zeros(H_htrig.Basis_num);
                            %H_htrig.Kinds = Kind;
                        end
                    end
                    if H_htrig.coe
                        vector = sym(vector);
                        for i = 1:size(vector,1)
                            vector_single = (vector(i,:));
                            try
                                if (ismember(vector_single,H_htrig.HsymL_coeL,'rows'))
                                    continue; % 进入下一个
                                end
                            catch
                                % ugly
                            end
                            Kind = H_htrig.Kinds +1;
                            H_htrig.HsymL_coeL(Kind,:) = (vector_single);
                            H_htrig.HcoeL(:,:,Kind) = zeros(H_htrig.Basis_num,'sym');
                            %H_htrig.Kinds = Kind;
                        end
                    end
                otherwise
                    Kind = H_htrig.Kinds+1;
                    BASIS_NUM = H_htrig.Basis_num;
                    H_htrig.HsymL_trig(Kind) = vector;
                    H_htrig.HcoeL(:,:,Kind) = zeros(BASIS_NUM,BASIS_NUM,'sym');
                    H_htrig.HnumL(:,:,Kind)  = (zeros(BASIS_NUM,BASIS_NUM));
            end

   
        end
    end
    %% modify
    methods
        function H_htrig2 = rewrite(H_htrig,mode)
            arguments
                H_htrig Htrig;
                mode {mustBeMember(mode,{'sincos','exp','mat','list',''})}= '';
            end
            if nargin < 2 || strcmp(mode , '')
                if strcmp(H_htrig.Type,'sincos')
                    mode = 'exp';
                elseif strcmp(H_htrig.Type,'exp')
                    mode = 'sincos';
                elseif strcmp(H_htrig.Type,'mat')
                    mode = 'list';
                elseif strcmp(H_htrig.Type,'list')
                    mode = 'mat';
                end
            end
            H_htrig2 = H_htrig;
            VarUsing = H_htrig.VarsSeqLcart(1:H_htrig.Dim);
            switch mode
                case 'exp'
                    HsymL_trig_tmp  = combine(rewrite(H_htrig.HsymL_trig,'exp'));
                    H_htrig2.HsymL_trig = sym([]);
                    H_htrig2.HsymL_trig_bk = [exp(1i*VarUsing) exp(-1i*VarUsing)];
                    H_htrig2.HcoeL = sym([]);
                    H_htrig2.Type = 'exp';
                    count = 0;
                    for i = 1:H_htrig.Kinds
                        [coeff_trig,symvar_list_trig,H_htrig2] = split_sym_eq(H_htrig2,HsymL_trig_tmp(i));
                        for j =1:numel(coeff_trig)
                            count = count+1;
                            k_cell{count} = symvar_list_trig(j);
                            mat_cell{count} = H_htrig.HcoeL(:,:,i);%
                            Var_cell{count} = coeff_trig(j);
                        end
                    end
                    H_htrig2 = H_htrig2.setup(Var_cell,k_cell,mat_cell);
                case 'sincos'
                case 'mat'
                    [~,~,H_htrig2] = H_htrig2.NumOrCoe();
                    WANNUM = H_htrig2.Basis_num;
                    %vectorList = int32([0,0,0]);
                    if H_htrig2.num
                        [vectorList,~,icL] = unique(H_htrig2.HsymL_numL(:,1:H_htrig.Dim),'rows');
                        KINDS= size(vectorList,1);
                        HnumLtmp = zeros(WANNUM,WANNUM,KINDS);
                        sizemesh = [WANNUM,WANNUM,KINDS];
                    else
                        [vectorList,~,icL] = unique(H_htrig2.HsymL_coeL(:,1:H_htrig.Dim),'rows');
                        KINDS= size(vectorList,1);
                        HcoeLtmp = sym(zeros(WANNUM,WANNUM,KINDS));
                        sizemesh = [WANNUM,WANNUM,KINDS];
                    end
                    if H_htrig2.num
                        iL = double(H_htrig2.HsymL_numL(:,H_htrig.Dim+1));
                        jL = double(H_htrig2.HsymL_numL(:,H_htrig.Dim+2));
                        indL = sub2ind(sizemesh,iL,jL,icL);
                        HnumLtmp(indL) = H_htrig2.HnumL;
                        H_htrig2.HnumL = HnumLtmp;
                        H_htrig2.HsymL_numL = vectorList;
                    else
                        iL = double(H_htrig2.HsymL_coeL(:,H_htrig.Dim+1));
                        jL = double(H_htrig2.HsymL_coeL(:,H_htrig.Dim+2));
                        indL = sub2ind(sizemesh,iL,jL,icL);
                        HcoeLtmp(indL) = H_htrig2.HcoeL;
                        H_htrig2.HcoeL = HcoeLtmp;
                        H_htrig2.HsymL_coeL = vectorList;
                    end
                    H_htrig2.Type = 'mat';
                case 'list'
                    if strcmp(H_htrig2.Type,'sincos')|strcmp(H_htrig2.Type,'exp')
                        [num_label,coe_label,H_htrig2] = H_htrig2.NumOrCoe();
                        if H_htrig2.num
                            
                        end
                    else
                        % WANNUM = H_htrig2.Basis_num;
                        [num_label,coe_label,H_htrig2] = H_htrig2.NumOrCoe();
                        if H_htrig2.num
                            if isvector(H_htrig2.HnumL)
                                warning('May not need to rewrite. Do nothing');
                                return;
                            end
                            NRPTS_ = numel(H_htrig2.HnumL);
                            %HnumLtmp = zeros(1,NRPTS_);
                            sizeHcoeL = size(H_htrig2.HnumL);
                            HnumLtmp = reshape(H_htrig2.HnumL,[NRPTS_,1]);
                            % vector program
                            [iL,jL,kL]= ind2sub(sizeHcoeL,1:NRPTS_);
                            vectorList = [H_htrig2.HsymL_numL(kL,:),iL.',jL.'];
                            H_htrig2.HnumL = HnumLtmp;
                            H_htrig2.HsymL_numL = vectorList;
                        elseif coe_label && ~num_label
                            if isvector(H_htrig2.HcoeL)
                                warning('May not need to rewrite. Do nothing');
                                return;
                            end
                            NRPTS_ = numel(H_htrig2.HcoeL);
                            %HnumLtmp = zeros(1,NRPTS_);
                            sizeHcoeL = size(H_htrig2.HcoeL);
                            HcoeLtmp = reshape(H_htrig2.HcoeL,[NRPTS_,1]);
                            % vector program
                            [iL,jL,kL]= ind2sub(sizeHcoeL,1:NRPTS_);
                            vectorList = [H_htrig2.HsymL_coeL(kL,:),iL.',jL.'];
                            H_htrig2.HcoeL = HcoeLtmp;
                            H_htrig2.HnumL = zeros(size(HcoeLtmp));
                            H_htrig2.HsymL_coeL = vectorList;
                        end

                    end
                     H_htrig2.Type = 'list';
                otherwise
            end

        end
        function H_htrig = diff(H_htrig,dir,options)
            arguments
                H_htrig Htrig;
                dir = 1;
                options.Accuracy = 1e-6;
            end
            switch H_htrig.Type
                case 'sincos'
                    VarUsing = H_htrig.VarsSeqLcart(1:H_htrig.Dim);
                    K = VarUsing;
                    % sin cos diff make sin 2 cos cos 2 sin
                    % the HsymL_trig_bk has all sin cos pair！
                    % actually we do not need clasify if HsymL has another
                    % symbolic term like: a, for HsymL has a,b,c, diff
                    % works failed
                    HsymL_trig_tmp = diff(H_htrig.HsymL,K(dir));
                    AdditionCoe = zeros(H_htrig.Kinds,1,'sym');
                    for i  = 1:H_htrig.Kinds
                        [AdditionCoeTmp,HsymL_trigTmp] = coeffs(HsymL_trig_tmp(i));
                        if isempty(AdditionCoeTmp) && isempty(HsymL_trigTmp)
                            AdditionCoe(i) = 0;
                            HsymL_trig_tmp(i) = 1;
                        else
                            AdditionCoe(i) = fold(@times,AdditionCoeTmp);
                            HsymL_trig_tmp(i) =  fold(@times,HsymL_trigTmp);
                        end
                    end
                    % mat mode
                    if H_htrig.num
                        H_htrig.HsymL_trig = HsymL_trig_tmp;
                        H_htrig.HnumL = vasplib.matrixtimespage(double(AdditionCoe),H_htrig.HnumL);
                    else
                        H_htrig.HsymL_trig = HsymL_trig_tmp;
                        H_htrig.HcoeL = vasplib.matrixtimespage(AdditionCoe,H_htrig.HcoeL);
                    end
                case 'exp'
                    syms k_x k_y k_z real;
                    K = [k_x;k_y;k_z];
                    % exp diff dont change anything
                    % mat mode
                    AdditionCoe = diff(H_htrig.HsymL,K(dir))./H_htrig.HsymL;
                    if H_htrig.num
                        H_htrig.HnumL = vasplib.matrixtimespage(double(AdditionCoe),H_htrig.HnumL);
                    else
                        H_htrig.HcoeL = vasplib.matrixtimespage(AdditionCoe,H_htrig.HcoeL);
                    end
                case 'mat'
                    % mat mode
                    if H_htrig.num
                        H_htrig.HnumL = 1i*vasplib.matrixtimespage(H_htrig.HsymL_numL(:,dir),H_htrig.HnumL);
                    else
                        H_htrig.HcoeL = 1i*vasplib.matrixtimespage(H_htrig.HsymL_coeL(:,dir),H_htrig.HcoeL);
                    end
                case 'list'
                    % constrain list mode
                    if H_htrig.num
                        H_htrig.HnumL = 1i*H_htrig.HnumL.*H_htrig.HsymL_numL(:,dir);
                    else
                        H_htrig.HcoeL = 1i*H_htrig.HcoeL.*H_htrig.HsymL_coeL(:,dir);
                    end
            end
            H_htrig =  simplify(H_htrig,options.Accuracy);
        end
        function H_htrig = simplify(H_htrig,Accuracy,options)
            arguments
                H_htrig Htrig;
                Accuracy = 1e-6;
                options.reduce = false;
            end
            if options.reduce
                for i = 1:numel(H_htrig.HcoeL)
                    if H_htrig.HcoeL(i)~=sym(0)
                        H_htrig.HcoeL(i) = vasplib.cleanVar(H_htrig.HcoeL,log10(Accuracy));
                    end
                end
            end
            [~,~,H_htrig] = NumOrCoe(H_htrig);
            if H_htrig.coe
                H_coeL_tmp = simplify(H_htrig.HcoeL);
                H_htrig.HcoeL = H_coeL_tmp;
                switch H_htrig.Type
                    case 'list'
                        Kinds_list = find(H_coeL_tmp ~=sym(0));
                        H_htrig = H_htrig.reseq(':',Kinds_list);
                        %[unique_HcoeL,ia,ic] = H_htrig.HcoeL;
                    case {'mat','exp','sincos'}
                        zerosMat = zeros(size(H_coeL_tmp(:,:,1)),'sym');
                        Kinds_list = true(H_htrig.Kinds,1);
                        for i = 1:H_htrig.Kinds
                            if isequal(zerosMat,H_coeL_tmp(:,:,i))
                                Kinds_list(i) = false;
                            end
                        end
                        H_htrig = H_htrig.reseq(':',Kinds_list);
                    otherwise
                end
            end
            if H_htrig.num
                H_numL_tmp = H_htrig.HnumL;
                %H_htrig.HnumL = H_numL_tmp;
                switch H_htrig.Type
                    case 'list'
                        Kinds_list = find(abs(H_numL_tmp) > Accuracy);
                        H_htrig = H_htrig.reseq(':',Kinds_list);
                        %[unique_HcoeL,ia,ic] = H_htrig.HcoeL;
                    case {'mat','exp','sincos'}
                        zerosMat = ones(H_htrig.Basis_num)*Accuracy;
                        Kinds_list = true(H_htrig.Kinds,1);
                        for i = 1:H_htrig.Kinds
                            if sum(sum(abs(H_numL_tmp(:,:,i)) > zerosMat))
                                %NRPTS_list(i) = true;
                            else
                                Kinds_list(i) = false;
                            end
                        end
                        H_htrig = H_htrig.reseq(':',Kinds_list);
                    otherwise
                end
            end
        end
    end
    methods
        function [HoutL] = HCAR_gen(H_htrig,klist,options)
            arguments
                H_htrig;
                klist;
                options.Hermi = true;
            end
            kn = size(klist,1);
            switch H_htrig.Type
                case 'list'
                    if isempty(H_htrig.Sparse_vector) && isempty(H_htrig.CutList)
                        H_htrig = SliceGen(H_htrig);
                    end
                    HoutL = zeros(H_htrig.Basis_num,H_htrig.Basis_num,kn);
                    Hnum_list_k = H_htrig.HnumL.*exp(1i*H_htrig.HsymL_numL(:,1:3)*klist.');
                    for ki = 1:kn
                        Hout = zeros(H_htrig.Basis_num,H_htrig.Basis_num);
                        Hnum_list_ktmp = Hnum_list_k(:,ki);
                        for i=1:H_htrig.N_Sparse_vector
                            Hout(H_htrig.Sparse_vector(i,1),H_htrig.Sparse_vector(i,2)) = sum(Hnum_list_ktmp(H_htrig.CutList(i,1):H_htrig.CutList(i,2)));
                        end
                        %Factorlist = H_htrig.HnumL.*exp(1i*H_htrig.HsymL_numL(:,1:3)*klist_cart_tmp(ki,:).');
                        %[ij_unique,sumFactorlist] = HollowKnight.generalcontractrow(H_htrig.HsymL_numL(:,4:5),Factorlist);
                        %indL = sub2ind(sizemesh,ij_unique(:,1),ij_unique(:,2));
                        %Hout(indL) = sumFactorlist;
                        if options.Hermi
                            HoutL(:,:,ki) = (Hout+Hout')/2;
                        else
                            HoutL(:,:,ki) = HoutL;
                        end
                    end
                case 'sparse'
                    % this type has not be coded
                    for ki =1:kn
                        Htemp=sparse(H_htrig.Basis_num ,H_htrig.Basis_num);
                        for ki=1:H_htrig.Kinds
                            Htemp = Htemp +Hnum_list{ki}*double(H_htrig.HsymL_trig(ki));
                        end
                        HoutL{ki} = Htemp;
                    end
                case 'mat'
                    HoutL = zeros(H_htrig.Basis_num,H_htrig.Basis_num,kn);
                    for ki =1:kn
                        Factorlist = exp(1i*H_htrig.HsymL_numL*klist(ki,:).');
                        Hout = sum(vasplib.matrixtimespage(Factorlist,H_htrig.HnumL),3);
                        if options.Hermi
                            HoutL(:,:,ki) = (Hout+Hout')/2;
                        else
                            HoutL(:,:,ki) = HoutL;
                        end
                    end
                otherwise
                    H_htrig.Htrig_num = subs(H_htrig.Htrig_sym);
                    H_htrig.Hfun = matlabFunction(H_htrig.Htrig_num,'Vars',[sym('k_x'),sym('k_y'),sym('k_z')]);
                    HoutL = zeros(H_htrig.Basis_num,H_htrig.Basis_num,kn);
                    for ki =1:kn
                        Hout = H_htrig.Hfun(klist(ki,1),klist(ki,2),klist(ki,3));
                        if options.Hermi
                            HoutL(:,:,ki) = (Hout+Hout')/2;
                        else
                            HoutL(:,:,ki) = HoutL;
                        end
                    end

            end
        end
        function varargout = diff_klist(H_htrig,dir,klist,options)
            arguments
                H_htrig
                dir = 1;
                klist = [];
                options.Accuracy = 1e-6;
            end
            if ~strcmp(H_htrig.Type,'list')
                %error('list enforced!!');
            end
            optionscell = namedargs2cell(options);
            H_htrig_tmp = H_htrig;
            varargout{nargout} = H_htrig_tmp.diff(dir(numel(dir)),optionscell{:}).HCAR_gen(klist);
            % constrain list mode
            for i = 1:numel(dir)-1
                %H_htrig2 = H_htrig.diff(dir(i),optionscell{:});
                H_htrig_tmp = H_htrig;
                varargout{i} = H_htrig_tmp.diff(dir(i),optionscell{:}).HCAR_gen(klist);
            end
        end
    end
    %% symmetry
    methods
%         function H_htrig = seeds(H_htrig,varlist)
%             if nargin < 2
%                 syms k_x k_y k_z real;
%                 varlist = [k_x k_y k_z];
%             end
%             HsymL_trig_in = sym([]);
%             for i = 1:length(varlist)
%                 HsymL_trig_in = [HsymL_trig_in,sin(varlist(i)),cos(varlist(i))];
%             end
%             H_htrig.HsymL_trig = HsymL_trig_in;
%             [sizeX,sizeY] = size(H_htrig.HcoeL);
%             sizeH = [sizeX,sizeY,H_htrig.Kinds];
%             HcoeL_tmp = sym('A',sizeH,'real')+1i*sym('B',sizeH,'real');
%             for i =1:H_htrig.Kinds
%                 HcoeL_tmp(:,:,i) = triu(HcoeL_tmp(:,:,i))+triu(HcoeL_tmp(:,:,i),1)';
%             end
%             H_htrig.HcoeL = HcoeL_tmp;
%             H_htrig = H_htrig.hermitize();
%         end
        function H_htrig = init(H_htrig,HsymL_trig_in,options)
            arguments
                H_htrig Htrig;
                HsymL_trig_in sym=sym([]);
                options.level_cut double = 1;
                options.per_dir double = [1,1,1];
                options.onsite_mode double = 0;
            end
            syms k_x k_y k_z real;
            k = [k_x;k_y;k_z];
            if nargin < 2
                %syms k_x k_y k_z real;
                if strcmp(H_htrig.Type,'sincos')
                    HsymL_trig_in = [sym(1),sin(k_x),cos(k_x),sin(k_y),cos(k_y),sin(k_z),cos(k_z)];
                elseif strcmp(H_htrig.Type,'exp')
                    
                end
            end
            if strcmp(H_htrig.Type,'sincos')
                H_htrig.HsymL_trig = HsymL_trig_in;
                [sizeX,sizeY] = size(H_htrig.HcoeL);
                sizeH = [sizeX,sizeY,H_htrig.Kinds];
                HcoeL_tmp = sym('A',sizeH,'real')+1i*sym('B',sizeH,'real');
                for i =1:H_htrig.Kinds
                    HcoeL_tmp(:,:,i) = triu(HcoeL_tmp(:,:,i))+triu(HcoeL_tmp(:,:,i),1)';
                end
                H_htrig.HcoeL = HcoeL_tmp;
            elseif strcmp(H_htrig.Type,'exp')
                select_nn_store = H_htrig.nn_store(H_htrig.nn_store(:,10)<=options.level_cut,:);
                nselect_nn_store = size(select_nn_store,1);
                if options.onsite_mode == 1
                    for i = 1:H_htrig.Basis_num
                        A = sym(['A','_',num2str(i),'_',num2str(i),'_0_ubar'],'real');
                        B = sym(['B','_',num2str(i),'_',num2str(i),'_0_ubar'],'real');
                        SymHopping = A+1i*B;
                        SymVar = sym(1);
                        H_htrig = H_htrig.set_hop(SymHopping,SymVar,[i,i]);
                    end
                end
                for i = 1:nselect_nn_store
                    A = Htrig.SymbolicVarible("A",select_nn_store(i,[6,7,8]),select_nn_store(i,1:2),select_nn_store(i,10));
                    B = Htrig.SymbolicVarible("B",select_nn_store(i,[6,7,8]),select_nn_store(i,1:2),select_nn_store(i,10));
                    SymHopping = A+1i*B;
                    SymVar = expand(exp(1i*select_nn_store(i,[3,4,5])*k));
                    H_htrig = H_htrig.set_hop(SymHopping,SymVar,select_nn_store(i,1:2));
                end
                H_htrig =H_htrig.hermitize();
                
            else
                
            end

        end
        function H_htrig = applyOper(H_htrig,SymOper,options)
            arguments
                H_htrig Htrig;
                SymOper Oper = Oper();
                options.fast = true; % will not check basis completity 
            end
            % [num_label,coe_label] = H_htrig.NumOrCoe();
            %                 num_label = false;
            %             else
            %                 num_label = true;
            %             end
            [~,coe_label] = H_htrig.NumOrCoe();
            
            if ~coe_label
                H_htrig = H_htrig.init();
                H_htrig = H_htrig.hermitize();
            end
            
            if length(SymOper) == 1
                if ~SymOper.conjugate && ~SymOper.antisymmetry && isequal(SymOper.R,eye(3))
                    return;
                end
                [H_htrig_bk,H_htrig]  = H_htrig.applyR(inv(SymOper.R));
                if isnan(SymOper.U)
                    %build U
                end
                % when apply U it will reseq    the HcoeL
                H_htrig_bk  = H_htrig_bk.applyU(SymOper.U,SymOper.conjugate,SymOper.antisymmetry);
                Equationlist_r = (real(H_htrig.HcoeL - H_htrig_bk.HcoeL) == 0);
                Equationlist_i = (imag(H_htrig.HcoeL - H_htrig_bk.HcoeL) == 0);
                %Equationlist_r = Htrig.isolateAll(Equationlist_r,real(H_htrig.HcoeL));
                %Equationlist_i = Htrig.isolateAll(Equationlist_i,imag(H_htrig.HcoeL));
                Equationlist_r = Htrig.isolateAll(Equationlist_r);
                Equationlist_i = Htrig.isolateAll(Equationlist_i);
                HcoeLtmp = H_htrig.HcoeL ;
                HcoeLtmp_r = subs(real(HcoeLtmp),lhs(Equationlist_r),rhs(Equationlist_r));
                HcoeLtmp_i = subs(imag(HcoeLtmp),lhs(Equationlist_i),rhs(Equationlist_i));
                H_htrig.HcoeL = HcoeLtmp_r + 1i*HcoeLtmp_i;
            else
                for i = 1:length(SymOper)
                    H_htrig = H_htrig.applyOper(SymOper(i));
                end
            end
        end
        function H_htrig = dualize(H_htrig)
            if strcmp(H_htrig.Type,'exp')
                BASIS_NUM = H_htrig.Basis_num;
                for i=1:length(H_htrig.HsymL_trig)
                    k_symbol = conj(H_htrig.HsymL_trig(i));
                    Kind = H_htrig.k_symbol2Kind(k_symbol);
                    if isempty(Kind)
                        Kind = H_htrig.Kinds+1;
                        H_htrig.HsymL_trig(Kind) = k_symbol;
                        H_htrig.HcoeL(:,:,Kind) = sym(zeros(BASIS_NUM,BASIS_NUM,1));
                        H_htrig.HnumL(:,:,Kind) = (zeros(BASIS_NUM,BASIS_NUM,1));
                    end
                end
                [~,H_htrig.Duality_vector_dist] = ismember(conj(H_htrig.HsymL_trig),H_htrig.HsymL_trig);
            end
        end
        function H_htrig = hermitize(H_htrig)
            [num_label,coe_label] = H_htrig.NumOrCoe();
            
            H_htrig_bk = H_htrig';
            if coe_label 
                Equationlist_r = real(H_htrig.HcoeL - H_htrig_bk.HcoeL) == 0;
                Equationlist_i = imag(H_htrig.HcoeL - H_htrig_bk.HcoeL) == 0;
                %Equationlist_r = Htrig.isolateAll(Equationlist_r,real(H_htrig.HcoeL));
                %Equationlist_i = Htrig.isolateAll(Equationlist_i,imag(H_htrig.HcoeL));
                Equationlist_r = Htrig.isolateAll(Equationlist_r);
                Equationlist_i = Htrig.isolateAll(Equationlist_i);
                HcoeLtmp = subs(H_htrig.HcoeL,lhs(Equationlist_r),rhs(Equationlist_r));
                HcoeLtmp = subs(HcoeLtmp,lhs(Equationlist_i),rhs(Equationlist_i));
                H_htrig.HcoeL = HcoeLtmp;
            end
            if num_label
                H_htrig.HnumL = (H_htrig_bk.HnumL + H_htrig.HnumL )/2;
            end
        end
        function H_htrig_bk = subsOper(H_htrig,SymOper)
            arguments
                H_htrig Htrig;
                SymOper Oper = Oper();
            end
            % [num_label,coe_label] = H_htrig.NumOrCoe();
            %                 num_label = false;
            %             else
            %                 num_label = true;
            %             end
            if isequal(sym(zeros(size(H_htrig.HcoeL))),H_htrig.HcoeL)
                coe_label = false;
            else
                coe_label = true;
            end
            
            if ~coe_label
                H_htrig = H_htrig.init();
                H_htrig = H_htrig.hermitize();
            end
            
            if length(SymOper) == 1
                if ~SymOper.conjugate && ~SymOper.antisymmetry && isequal(SymOper.R,eye(3))
                    return;
                end
                [H_htrig_bk,H_htrig]  = H_htrig.applyR(inv(SymOper.R));
                if isnan(SymOper.U)
                    %build U
                end
                H_htrig_bk  = H_htrig_bk.applyU(SymOper.U,SymOper.conjugate,SymOper.antisymmetry);
            end
        end
        function H_htrig = applyU(H_htrig,U,conjugate ,antisymmetry )
            arguments
                H_htrig Htrig;
                U =nan;
                conjugate logical =false;
                antisymmetry logical = false;
            end
            if isa(U,'Oper')
                conjugate = U.conjugate;
                antisymmetry = U.antisymmetry;
                U = U.U;
            end
            [num_label,coe_label] = H_htrig.NumOrCoe();
            
            U_inv = inv(U);
            if coe_label == true
                HcoeLtmp = H_htrig.HcoeL;
                if conjugate
                    if strcmp(H_htrig.Type,'sincos')
                        H_htrig = H_htrig.applyR(diag([-1,-1,-1]));
                        HcoeLtmp = H_htrig.HcoeL;
                    end
                    %H_htrig = H_htrig.applyR(diag([-1,-1,-1]));
                    HcoeLtmp = conj(HcoeLtmp);
                    %HcoeLtmp = Htrig.matrixtimespage(H_htrig.factorlist_parity(),HcoeLtmp);
                end
                if antisymmetry
                    HcoeLtmp = -HcoeLtmp;
                end
                HcoeLtmp = Htrig.page_mtimes_matrix(Htrig.matrix_mtimes_page(U,HcoeLtmp),U_inv);
                H_htrig.HcoeL = HcoeLtmp;
            end
            if num_label == true
                U_page = repmat(U,[1 1 size(H_htrig.HnumL,3)]);
                U_inv_page = repmat(U_inv,[1 1 size(H_htrig.HnumL,3)]);
                HnumLtmp = H_htrig.HnumL;
                if conjugate
                    HnumLtmp = conj(HnumLtmp);
                    HnumLtmp = Htrig.matrixtimespage(H_htrig.factorlist_parity(),HnumLtmp);
                end
                if antisymmetry
                    HnumLtmp = -HnumLtmp;
                end
                HnumLtmp = pagemtimes(pagemtimes(U_page,HnumLtmp),U_inv_page);
                H_htrig.HnumL = HnumLtmp;
            end
        end
        function [H_htrig,H_htrig2] = applyR(H_htrig,R)
            arguments
                H_htrig Htrig;
                R ;
            end
            [num_label,coe_label] = H_htrig.NumOrCoe();
            

            %applyk
            [H_htrig,Smat] = H_htrig.Smatgen(R);
            H_htrig2 = H_htrig;
            if num_label
                H_htrig.HnumL = Htrig.matrixtimespage(Smat,H_htrig.HnumL);
            end
            if coe_label
                H_htrig.HcoeL = Htrig.matrixtimespage(Smat,H_htrig.HcoeL);
            end
            %apply R
            
        end
        function [H_htrig,H_htrig2] = applyRt(H_htrig,Rt)
        end
        function [H_htrig,Smat] = Smatgen(H_htrig,R,Accuracy)
            arguments
                H_htrig Htrig;
                R  ;
                Accuracy double = 6;
            end
            if isa(R,'sym')
                coe_label = true;
            else
                coe_label = false;
            end
            BASIS_NUM = H_htrig.Basis_num;
            %HsymC_bk = H_htrig.HsymL_trig;
            HsymC = H_htrig.HsymL_trig;
            syms k_x k_y k_z real;
            % add in property
            %             varlist = [sin(k_x),sin(k_y),sin(k_z),...
            %                 cos(k_x),cos(k_y),cos(k_z)...
            %                 ];
            %varlist = H_htrig.HsymL_trig_bk;
            
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
                [coeff_trig,~] = coeffs(HsymC(k),H_htrig.HsymL_trig_bk);
                for i = 1:numel(coeff_trig)
                    tmp_label = contains(string(coeff_trig),H_htrig.seeds);
                    if sum(tmp_label)
                        [~,coeff_trig_list] = coeffs(coeff_trig(i));
                        for j = 1:numel(coeff_trig_list)
                            tmp_label2 = contains(string(coeff_trig_list(j)),H_htrig.seeds);
                            if sum(tmp_label2)
                                H_htrig = H_htrig.find_HsymL_trig_bk(coeff_trig_list(j));
                            end
                        end
                    end
                end
            end
            % find new basis & update
            for i = 1:numel(HsymC)
                [~,B] = coeffs(HsymC(i),H_htrig.HsymL_trig_bk);
                for k = 1:length(B)
                   Kind = H_htrig.k_symbol2Kind(B(k));
                   if isempty(Kind)
                       Kind = H_htrig.Kinds+1;
                       H_htrig.HsymL_trig(Kind) = B(k);
                       H_htrig.HcoeL(:,:,Kind) = sym(zeros(BASIS_NUM,BASIS_NUM,1));
                       H_htrig.HnumL(:,:,Kind)  = (zeros(BASIS_NUM,BASIS_NUM,1));
                   end
                end
            end
            % redo
            HsymC_bk = H_htrig.HsymL_trig;
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
                [A,B] = coeffs(HsymC(i),H_htrig.HsymL_trig_bk);
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
            if ~coe_label
                Smat = roundn(double(Smat),-Accuracy);
            end
            %                 if i == 3
            %                     error('ss');
            %                 end
        end
        function Factorlist_parity = factorlist_parity(H_htrig)
            syms k_x k_y k_z real;
            HsymC = H_htrig.HsymL_trig;
            HsymC = subs(HsymC,[k_x k_y k_z],-[k_x k_y k_z]);
            Factorlist_parity = simplify(HsymC./H_htrig.HsymL_trig);
        end
        function H_htrig = nn(H_htrig,search_range,Accuracy,Rlength_cut)
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
            H_htrig = nn@vasplib(H_htrig,search_range,Accuracy,Rlength_cut);
            H_htrig.Type = 'exp';
        end
    end
    %% get methods
    methods
        function symvarL = get.symvarL(H_htrig)
            symvarL = [H_htrig.symvar_list,symvar(H_htrig.HsymL)];
        end
        function Htrig_sym = get.Htrig_sym(H_htrig)
            Htrig_sym = sym(zeros(H_htrig.Basis_num,H_htrig.Basis_num));
            if H_htrig.num
                H_htrig.HcoeL = sym(H_htrig.HnumL);
                H_htrig.HsymL_coeL = sym(H_htrig.HsymL_numL);
                %H_htrig.HsymL_trig = sym(H_htrig.HsymL_trig);
            end
            switch H_htrig.Type
                case {'exp','sincos'}
                    try
                        for i =1:H_htrig.Kinds
                            Htrig_sym = Htrig_sym + H_htrig.HcoeL(:,:,i)*H_htrig.HsymL_trig(i);
                        end
                    catch
                        Htrig_sym =sym([]);
                    end
                case {'slab'}
                    try
                        for i =1:H_htrig.Kinds
                            Htrig_sym = Htrig_sym + H_htrig.HcoeL(:,:,i)*H_htrig.HsymL_trig(i);
                        end
                    catch
                        Htrig_sym =sym([]);
                    end
                    if strcmp(H_htrig.Type , 'slab' )
                        tmpHsymL_trig = H_htrig.HsymL_trig;
                        for i = 1:length(tmpHsymL_trig)
                            A = char(tmpHsymL_trig(i));
                            A(end)= A(end-3);
                            A(end-3)='N';
                            tmpHsymL_trig(i) = str2sym(A);
                        end
                        Htrig_sym = tril(Htrig_sym)+subs(triu(Htrig_sym,1),H_htrig.HsymL_trig,tmpHsymL_trig);
                    end
                case 'mat'
                    syms k_x k_y k_z;
                    HvarL = exp(1i*H_htrig.HsymL_coeL*[k_x;k_y;k_z]);
                    for i =1:H_htrig.Kinds
                        Htrig_sym = Htrig_sym + H_htrig.HcoeL(:,:,i)*HvarL(i);
                    end
                case 'list'
                    syms k_x k_y k_z;
                    Factorlist = H_htrig.HcoeL.*exp(1i*H_htrig.HsymL_coeL(:,1:3)*[k_x;k_y;k_z]);
                    [ij_unique,sumFactorlist] = HollowKnight.generalcontractrow(double(H_htrig.HsymL_coeL(:,4:5)),Factorlist);
                    Htrig_sym = zeros(H_htrig.Basis_num,'sym');
                    indL = sub2ind([H_htrig.Basis_num,H_htrig.Basis_num],ij_unique(:,1),ij_unique(:,2));
                    Htrig_sym(indL) = sumFactorlist;
            end

            
        end
        function Htrig_latex = get.Htrig_latex(H_htrig)
            Htrig_latex = latex(H_htrig.Htrig_sym);
        end
        function HsymL = get.HsymL(H_htrig)
            switch H_htrig.Type
                case {'exp','sincos'}
                    HsymL = H_htrig.HsymL_trig;
                case {'mat','list'}
                    if H_htrig.coe
                        HsymL = H_htrig.HsymL_coeL;
                    else
                        HsymL = H_htrig.HsymL_numL;
                    end
                otherwise
                    HsymL = H_htrig.HsymL_trig;
            end
        end
        function Kinds = get.Kinds(H_htrig)
            switch H_htrig.Type
                case {'exp','sincos'}
                    Kinds = length(H_htrig.HsymL_trig);
                case {'mat','list'}
                    Kinds = max(size(H_htrig.HsymL_numL,1),size(H_htrig.HsymL_coeL,1));
                otherwise
                    Kinds = length(H_htrig.HsymL_trig);
            end
        end
    end
    %% methods reload
    methods
        function C = plus(A,B)
            if isa(A,'Htrig') && isa(B,'Htrig')
                C = A;
                mat_cell = mat2cell(reshape(B.HcoeL,B.Basis_num,B.Basis_num*B.Kinds).',repmat(B.Basis_num,[1 B.Kinds]));
                Var_cell = mat2cell(sym(ones(1,B.Kinds)),1);
                k_cell = mat2cell(B.HsymL_trig,1);
                C = C.setup(Var_cell,k_cell,mat_cell);
            elseif isa(A,'Htrig') && ~isa(B,'Htrig')
                if isa(B,'Term')
                    % disp('gg');
                    C = A;
                    % term2trig
                    %                     for i = 1:length(B)
                    %                         C = C.setup_rough(B(i).symbolic_polynomial,B(i).pauli_mat);
                    %                     end
                elseif isa(B,'Trig')
                    C = A;
                    C.Trig_list = C.Trig_list+B;
                    for i =1:length(B)
                        C = C.setup_rough(B(i).symbolic_polynomial,B(i).pauli_mat);
                    end
                elseif isa(B,'sym')
                    basis_num = A.Basis_num;
                    if basis_num~= length(B) || size(B,2)  ~=  size(B,2)
                        error('Basis wrong!')
                    end
                    for i = 1:basis_num
                        for j = 1: basis_num
                            if B(i,j)~=sym(0)
                                tempmat = zeros( basis_num);
                                tempmat(i,j) =1 ;
                                %                                 switch A.Type
                                %                                     case 'exp'
                                %                                         SymPoly = (rewrite(B(i,j),'exp'));
                                %                                     case 'sincos'
                                %                                         SymPoly = (rewrite(B(i,j),'sincos'));
                                %                                     otherwise
                                %                                         SymPoly = B(i,j);
                                %                                 end
                                SymPoly = B(i,j);
                                A = A.setup_rough(SymPoly,tempmat);
                            end
                        end
                    end
                    C = A;
                else
                    C = A;
                end
            elseif ~isa(A,'Htrig') && isa(B,'Htrig')
                if isa(A,'Term')
                    C = B;
                    %                     for i = 1:length(A)
                    %                         C = C.setup_rough(A(i).symbolic_polynomial,A(i).pauli_mat);
                    %                     end
                elseif isa(A,'Trig')
                    C = B;
                    C.Trig_list = C.Trig_list+A;
                    for i = 1:length(A)
                        C = C.setup_rough(A(i).symbolic_polynomial,A(i).pauli_mat);
                    end
                else
                    C = B;
                end
            else
                C = 0;
            end
            %disp(C.Term_to_save);
        end
        function H_htrig = uminus(H_htrig)
            H_htrig.HnumL = -H_htrig.HnumL ;
            H_htrig.HcoeL = -H_htrig.HcoeL ;
            %
            error('not be implmented');
        end
        function C = minus(A,B)
            C = plus(A,-B);
        end
        function C = mtimes(A,B)
            if isa(A,'Htrig') && isa(B,'Htrig')
                error('not be implmented');
            elseif isa(A,'Htrig') && ~isa(B,'Htrig')
                C =A;
                if isa(B,'double')
                    % disp('gg');
                    for i = 1:C.Kinds
                        C.HcoeL(:,:,i) = C.HcoeL(:,:,i)*B;
                        C.HnumL(:,:,i) = C.HnumL(:,:,i)*B;
                    end
                elseif isa(B,'sym')
                    for i = 1:C.Kinds
                        C.HcoeL(:,:,i) = C.HcoeL(:,:,i)*B;
                        %                         A.HnumL(:,:,i) = A.HcoeL(:,:,i)*B;
                    end
                else
                    
                end
            elseif ~isa(A,'Htrig') && isa(B,'Htrig')
                C = B;
                if isa(A,'double')
                    % disp('gg');
                    for i = 1:C.Kinds
                        C.HcoeL(:,:,i) = A*C.HcoeL(:,:,i);
                        C.HnumL(:,:,i) = A*C.HnumL(:,:,i);
                    end
                elseif isa(A,'sym')
                    for i = 1:C.Kinds
                        C.HcoeL(:,:,i) = A*C.HcoeL(:,:,i);
                        %                         A.HnumL(:,:,i) = A.HcoeL(:,:,i)*B;
                    end
                else
                    
                end
            else
                C = 0;
            end
        end
        % ?
        function C = mtimes_inv(B,A)
            if ~isa(A,'Htrig') && isa(B,'Htrig')
                C = B;
                if isa(A,'double')
                    % disp('gg');
                    for i = 1:C.Kinds
                        C.HcoeL(:,:,i) = A*C.HcoeL(:,:,i);
                        C.HnumL(:,:,i) = A*C.HnumL(:,:,i);
                    end
                elseif isa(A,'sym')
                    for i = 1:C.Kinds
                        C.HcoeL(:,:,i) = A*C.HcoeL(:,:,i);
                        %                         A.HnumL(:,:,i) = A.HcoeL(:,:,i)*B;
                    end
                else
                    
                end
            end
        end
        % overload 	lt(A,B)
        function C = gt(B,A)
            if isa(A,'Htrig') && isa(B,'Htrig')
                H_htrig1 = A;
                H_htrig2 = B;
                % ---------check----------
                if H_htrig1.Bassis_num ~= H_htrig2.Bassis_num
                    error('Bassis_num different');
                end
                error('not support at present.')
            elseif isa(A,'Htrig') && ~isa(B,'Htrig')
                switch class(B)
                    case 'char'
                        switch B(1)
                            case {'P','p'}
                                C = A.input_orb_struct(B,'sym');
                                C.Rm = sym(C.Rm );
                                C.orbL = sym(C.orbL );
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
            elseif ~isa(A,'Htrig') && ~isa(B,'Htrig')
                error('not support at present.');
            end
        end
        function C = lt(A,B)
            if isa(A,'Htrig') && isa(B,'Htrig')
                H_htrig1 = A;
                H_htrig2 = B;
                % ---------check----------
                if H_htrig1.Bassis_num ~= H_htrig2.Bassis_num
                    error('Bassis_num different');
                end
                error('not support at present.')
            elseif isa(A,'Htrig') && ~isa(B,'Htrig')
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
            elseif ~isa(A,'Htrig') && ~isa(B,'Htrig')
                error('not support at present.');
            end
            
        end
        function C = le(A,B)
            if isa(A,'Htrig') && isa(B,'Htrig')
                H_htrig1 = A;
                H_htrig2 = B;
                % ---------check----------
                if H_htrig1.Bassis_num ~= H_htrig2.Bassis_num
                    error('Bassis_num different');
                end
                error('not support at present.')
            elseif isa(A,'Htrig') && ~isa(B,'Htrig')
                switch class(B)
                    case 'char'
                        switch B(1)
                            case {'P','p'}
                                C = A.input_orb_struct(B,'sym');
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
            elseif ~isa(A,'Htrig') && ~isa(B,'Htrig')
                error('not support at present.');
            end
        end
        function C = horzcat(A,B)
            if isa(A,'Htrig') && isa(B,'Htrig')
                H_htrig1 = A;
                H_htrig2 = B;
                
                % ---------init-----------
                %                 H_htrig =  HR(H_htrig1.WAN_NUM+H_htrig2.WAN_NUM,...
                %                     unique([H_htrig1.vectorL;H_htrig2.vectorL],'rows'));
                H_htrig = A;
                H_htrig.Basis = [H_htrig1.Basis;H_htrig2.Basis];
                H_htrig.Basis_num = H_htrig1.Basis_num+H_htrig2.Basis_num;
                H_htrig.HnumL = zeros(H_htrig.Basis_num,H_htrig.Basis_num,H_htrig.Kinds);
                H_htrig.HcoeL = sym(H_htrig.HnumL);
                %             zeros_num_mat = zeros(H_htrig.WAN_NUM);
                %             zeros_coe_mat = sym(zeros(H_htrig.WAN_NUM));
                % need select here
                for i = 1:H_htrig.Kinds
                    
                    H_htrig.HnumL(:,:,i) = blkdiag(H_htrig1.HnumL(:,:,i) ,...
                        H_htrig2.HnumL(:,:,i));
                    H_htrig.HcoeL(:,:,i) = blkdiag(H_htrig1.HcoeL(:,:,i) ,...
                        H_htrig2.HcoeL(:,:,i));
                    
                end
                H_htrig.Trig_to_save =sym(zeros(H_htrig.Basis_num,H_htrig.Basis_num));
                %H_htrig.Term_to_save =blkdiag(H_htrig1.Term_to_save,H_htrig2.Term_to_save);
            else
            end
            C = H_htrig;
        end
        function H_htrig = ctranspose(H_htrig)
            [num_label,coe_label] = H_htrig.NumOrCoe();
            
            for i =1:H_htrig.Kinds
                if num_label
                    H_htrig.HnumL(:,:,i) = H_htrig.HnumL(:,:,i)';
                end
                if coe_label
                    H_htrig.HcoeL(:,:,i) = H_htrig.HcoeL(:,:,i)';
                end
            end
            if strcmp(H_htrig.Type,'exp')
                H_htrig =H_htrig.dualize();
                if num_label
                    H_htrig.HnumL(:,:,:) = H_htrig.HnumL(:,:,H_htrig.Duality_vector_dist);
                end
                if coe_label
                    H_htrig.HcoeL(:,:,:) = H_htrig.HcoeL(:,:,H_htrig.Duality_vector_dist);
                end
            end
        end
        function H_htrig = conj(H_htrig)
            [num_label,coe_label] = H_htrig.NumOrCoe();
            for i =1:H_htrig.Kinds
                if num_label
                    H_htrig.HnumL(:,:,i) = conj(H_htrig.HnumL(:,:,i));
                end
                if coe_label
                    H_htrig.HcoeL(:,:,i) = conj(H_htrig.HcoeL(:,:,i)');
                end
            end
            H_htrig.HsymL_trig = conj(H_htrig.HsymL_trig);
        end
        % disp
        %         function Htrig_sym= disp(H_htrig)
        %             % Htrig_sym = ;
        %             %
        %             % disp(H_htrig.Htrig_latex);
        %             builtin('disp',H_htrig);
        %             %disp(H_htrig.Htrig_latex);
        %             Htrig_sym = H_htrig.Htrig_sym;
        %         end
        function Htrig_sym = sym(H_htrig,options)
            arguments
                H_htrig Htrig;
                options.simple = false;
                options.Type {mustBeMember(options.Type,{'exp','sincos','sin','cos','cot','tan','sinh','cosh','cosh','tanh',''})}= '';
            end
            Htrig_sym = H_htrig.Htrig_sym;
            if strcmp(options.Type,'')
                switch H_htrig.Type
                    case {'sincos'}
                        type = 'sincos';
                    case {'exp','list','mat'}
                        type = 'exp';
                    otherwise
                        type = 'sincos';
                end
            else
                type = options.Type;
            end
            if options.simple
                H_htrig = H_htrig.timtj_gen('sym');
                Htrig_sym = simplify(Htrig_sym./H_htrig.timtj{3}.');
            end
            Htrig_sym = simplify(rewrite(Htrig_sym,type),'IgnoreAnalyticConstraints',true);
        end
        function Htrig_latex = latex(H_htrig)
            Htrig_latex = H_htrig.Htrig_latex;
        end
    end
    %% manipulate
    methods
        function H_htrig = translate(H_htrig,U)
            U_inv = inv(U);
            if isa(U,'double')
                % disp('gg');
                for i = 1:C.Kinds
                    H_htrig.HcoeL(:,:,i) = U_inv*H_htrig.HcoeL(:,:,i)*U;
                    H_htrig.HnumL(:,:,i) = U_inv*H_htrig.HnumL(:,:,i)*U;
                end
            elseif isa(U,'sym')
                for i = 1:H_htrig.Kinds
                    H_htrig.HcoeL(:,:,i) = U_inv*H_htrig.HcoeL(:,:,i)*U;
                    %                         A.HnumL(:,:,i) = A.HcoeL(:,:,i)*B;
                end
            else
                
            end
        end
        function [H_htrig,EQL] = subsVar(H_htrig,varargin)
            if strcmp(varargin{end},'all')
                all_mode = true;
                nargin_check = nargin-1;
            else
                all_mode = false;
                nargin_check = nargin;
            end
            switch nargin_check
                case 1
                    H_htrig.HcoeL = expand(simplify(subs(H_htrig.HcoeL)));
                case 2
                    SymVarL = H_htrig.symvar_list;
                    if all_mode
                        H_htrig.HcoeL = subs(H_htrig.HcoeL,SymVarL,varargin{1});
                    else
                        H_htrig.HcoeL = (simplify(subs(H_htrig.HcoeL,SymVarL,varargin{1})));
                    end
                    EQL =(SymVarL==vpa(varargin{1}));
                case 3
                    if all_mode
                        H_htrig.HcoeL = subs(H_htrig.HcoeL,varargin{1},varargin{2});
                    else
                        H_htrig.HcoeL = simplify(subs(H_htrig.HcoeL,varargin{1},varargin{2}));
                    end
                    EQL =(varargin{1}==vpa(varargin{2}));
                case 5
                    H_htrig.HcoeL = simplify(subs(H_htrig.HcoeL,varargin{1},varargin{2}));
                    H_htrig.ScoeL = simplify(subs(H_htrig.ScoeL,varargin{3},varargin{4}));
                    EQL{1} =(varargin{1}==vpa(varargin{2}));
                    EQL{2} =(varargin{3}==vpa(varargin{3}));
            end
            if  isempty(H_htrig.symvar_list) 
                H_htrig = H_htrig.Subsall();
            else

            end
            % waiting to add ...
        end
        function H_htrig2 = subs(H_htrig,varargin)
            HsymL_trig_tmp = H_htrig.HsymL_trig;
            switch length(varargin)
                case 0
                    HsymL_trig_tmp = expand(simplify(subs(HsymL_trig_tmp)));
                case 1
                    HsymL_trig_tmp = expand(simplify(subs(HsymL_trig_tmp,varargin{1})));
                case 2
                    HsymL_trig_tmp = expand(simplify(subs(HsymL_trig_tmp,varargin{1},varargin{2})));
            end
            % reset HcoeL
            H_htrig2 = H_htrig;
            H_htrig2.HcoeL = sym([]);
            H_htrig2.HnumL = [];
            H_htrig2.HsymL_trig = sym([]);
            count = 0;
            for i = 1:H_htrig.Kinds
                [coeff_trig,symvar_list_trig,H_htrig2] = split_sym_eq(H_htrig2,HsymL_trig_tmp(i));
                for j =1:numel(coeff_trig)
                    count = count+1;
                    k_cell{count} = symvar_list_trig(j);
                    mat_cell{count} = H_htrig.HcoeL(:,:,i);%
                    Var_cell{count} = coeff_trig(j);
                end
            end
            H_htrig2 = H_htrig2.setup(Var_cell,k_cell,mat_cell);
        end
        function H_htrig = rotation(H_htrig,Rotation)
            if nargin <2
                Rotation = inv(sym(H_htrig.Rm));
            end
            if ischar(Rotation)
                if strcmp(Rotation,'auto')
                    Rotation = inv(sym(H_htrig.Rm));
                end
            end
            % syms k_x k_y k_z real;
            VarUsing = H_htrig.VarsSeqLcart(1:H_htrig.Dim);
            k = Rotation * VarUsing.';
            for i = 1:H_htrig.Dim
                VarUsingCell{i} = VarUsing(i);
                kCell{i} = k(i);
            end
            H_htrig = H_htrig.subs(VarUsingCell,kCell);
        end
        function H_htrig2 = discretize(H_htrig,Nslab,options)
            arguments
                H_htrig Htrig;
                Nslab double = [0,10,0];
                options.rmfunc function_handle=@()(1);
                options.Rotation  = sym(eye(3));
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
            % rotate the mode first
            if rotate_mode
                H_htrig = H_htrig.rotation(options.Rotation);
            end
            % replace seeds as Sigma_x Sigma_y Sigma_z
            H_htrig2 = Htrig(H_htrig.Basis_num,'Dim',H_htrig.Dim);
            H_htrig2.seeds = ["Sigma_x","Sigma_y","Sigma_z","Sigma_w"];
            H_htrig2.Type = 'slab';
            H_htrig2.Nslab = Nslab;
            NSLAB = (H_htrig2.Nslab ==0) + H_htrig2.Nslab;
            NS = fold(@times,NSLAB);
            % subs in any(1-4) direction
            case_d = Nslab>1;
            % expand seeds for more dimension
            seeds = ["x","y","z","w"];
            StrSinUsing = "sin(k_" + seeds + ")";
            StrCosUsing = "cos(k_" + seeds + ")";
            StrSigma_1__N  = "Sigma_" + seeds + "_1__N";
            StrSigma_N__N  = "Sigma_" + seeds + "_N__N";
            pat_d_pre   = "Sigma_" + seeds + "_";

            SymSinUsing = str2sym(StrSinUsing);
            SymCosUsing = str2sym(StrCosUsing);
            SymSigma_1__N = str2sym(StrSigma_1__N);
            SymSigma_N__N = str2sym(StrSigma_N__N);
            %syms Sigma_x_1__N Sigma_x_N__N real;
            %syms Sigma_y_1__N Sigma_y_N__N real;
            %syms Sigma_z_1__N Sigma_z_N__N real;
            %syms Sigma_w_1__N Sigma_w_N__N real;
            for d = 1:numel(case_d)
                if case_d(d)
                    H_htrig.HsymL_trig = subs(H_htrig.HsymL_trig,SymSinUsing(d),1i/2 * SymSigma_1__N(d));
                    H_htrig.HsymL_trig = subs(H_htrig.HsymL_trig,SymCosUsing(d),1 /2 * SymSigma_1__N(d));
                    pat_d = pat_d_pre(d)+digitsPattern(1)+"__N";
                    for i =1:length(H_htrig.HsymL_trig)
                        if ~contains(string(H_htrig.HsymL_trig(i)),pat_d)
                            H_htrig.HsymL_trig(i) = H_htrig.HsymL_trig(i)*SymSigma_N__N(d);
                        end
                    end
                    H_htrig.Sigmas =[H_htrig.Sigmas;[SymSigma_N__N(d),SymSigma_1__N(d)]];
                end
            end
            %             if case_x
            %                 H_htrig.HsymL_trig = subs(H_htrig.HsymL_trig,str2sym('sin(k_x)'),str2sym('1i*Sigma_x_1__N/2'));
            %                 H_htrig.HsymL_trig = subs(H_htrig.HsymL_trig,str2sym('cos(k_x)'),str2sym('Sigma_x_1__N/2'));
            %                 pat_x = "Sigma_x_"+digitsPattern(1)+"__N";
            %                 for i =1:length(H_htrig.HsymL_trig)
            %                     if ~contains(string(H_htrig.HsymL_trig(i)),pat_x)
            %                         H_htrig.HsymL_trig(i) = H_htrig.HsymL_trig(i)*str2sym('Sigma_x_N__N');
            %                     end
            %                 end
            %                 H_htrig.Sigmas =[H_htrig.Sigmas;[str2sym('Sigma_x_N__N'),str2sym('Sigma_x_1__N')]];
            %             end
            %             if case_y
            %                 H_htrig.HsymL_trig = subs(H_htrig.HsymL_trig,str2sym('sin(k_y)'),str2sym('1i*Sigma_y_1__N/2'));
            %                 H_htrig.HsymL_trig = subs(H_htrig.HsymL_trig,str2sym('cos(k_y)'),str2sym('Sigma_y_1__N/2'));
            %                 pat_y = "Sigma_y_"+digitsPattern(1)+"__N";
            %                 for i =1:length(H_htrig.HsymL_trig)
            %                     if ~contains(string(H_htrig.HsymL_trig(i)),pat_y)
            %                         H_htrig.HsymL_trig(i) = H_htrig.HsymL_trig(i)*str2sym('Sigma_y_N__N');
            %                     end
            %                 end
            %                 H_htrig.Sigmas =[H_htrig.Sigmas;[str2sym('Sigma_y_N__N'),str2sym('Sigma_y_1__N')]];
            %             end
            %             if case_z
            %                 H_htrig.HsymL_trig = subs(H_htrig.HsymL_trig,str2sym('sin(k_z)'),str2sym('1i*Sigma_z_1__N/2'));
            %                 H_htrig.HsymL_trig = subs(H_htrig.HsymL_trig,str2sym('cos(k_z)'),str2sym('Sigma_z_1__N/2'));
            %                 pat_z = "Sigma_z_"+digitsPattern(1)+"__N";
            %                 for i =1:length(H_htrig.HsymL_trig)
            %                     if ~contains(string(H_htrig.HsymL_trig(i)),pat_z)
            %                         H_htrig.HsymL_trig(i) = H_htrig.HsymL_trig(i)*str2sym('Sigma_z_N__N');
            %                     end
            %                 end
            %                 H_htrig.Sigmas =[H_htrig.Sigmas;[str2sym('Sigma_z_N__N'),str2sym('Sigma_z_1__N')]];
            %             end
            % reset HcoeL
            H_htrig2.HsymL_trig_bk = [SymSigma_1__N SymSigma_N__N];
            H_htrig2.HsymL_trig = sym([]);
            count = 0;
            for i = 1:H_htrig.Kinds
                [coeff_trig,symvar_list_trig,H_htrig2] = split_sym_eq(H_htrig2,H_htrig.HsymL_trig(i));
                for j =1:numel(coeff_trig)
                    count = count+1;
                    k_cell{count} = symvar_list_trig(j);
                    mat_cell{count} = H_htrig.HcoeL(:,:,i);%
                    Var_cell{count} = coeff_trig(j);
                end
            end
            H_htrig2 = H_htrig2.setup(Var_cell,k_cell,mat_cell);
            % supercell orb
            orb_tmp = zeros(NS*H_htrig.Basis_num,H_htrig.Dim);
            NWAVE  = NS*H_htrig.Basis_num;
            if isempty(H_htrig.orbL)
                H_htrig.orbL = zeros(H_htrig.Basis_num,H_htrig.Dim);
            end
            % General form need!!!
            switch H_htrig.Dim
                case 1
                    [i1L            ] = ind2sub(NSLAB,1:NS);
                    OrbAddL = i1L.' -1;
                case 2
                    [i1L,i2L        ] = ind2sub(NSLAB,1:NS);
                    OrbAddL = [i1L.' i2L.'] -1;
                case 3
                    [i1L,i2L,i3L    ] = ind2sub(NSLAB,1:NS);
                    OrbAddL = [i1L.' i2L.' i3L.'] -1;
                case 4
                    [i1L,i2L,i3L,i4L] = ind2sub(NSLAB,1:NS);
                    OrbAddL = [i1L.' i2L.' i3L.' i4L.'] -1;
            end
            for i=1:NS
                orb_tmp((i-1)*H_htrig.Basis_num+1:i*H_htrig.Basis_num,:) = H_htrig.orbL+OrbAddL(i,:);
            end
            orb_tmp = orb_tmp./NSLAB;
            SizeOrb_tmp = size(orb_tmp);
            % if cut the final mode
            if rm_mode
                try
                    for d = 1:SizeOrb_tmp(2)
                        Input{d} = orb_tmp(:,d);
                    end
                    H_htrig2.rm_list = options.rmfunc(Input{:});
                catch
                    H_htrig2.rm_list = false(1,NWAVE);
                    for i =1:NWAVE
                        for d = 1:SizeOrb_tmp(2)
                            Input{d} = orb_tmp(i,d);
                        end
                        H_htrig2.rm_list(i) = options.rmfunc(Input{:});
                    end
                end
            else
               
            end
            orb_tmp( H_htrig2.rm_list,:) = [];
            H_htrig2.orbL = orb_tmp;
            % Hmatpre
            H_htrig2.Hmat_pre{numel(H_htrig2.HsymL_trig)} = sparse(NS,NS);
            for i = 1:numel(H_htrig2.HsymL_trig)
                H_htrig2.Hmat_pre{i} = H_htrig2.HsymL_trig2mat(H_htrig2.HsymL_trig(i));
            end
        end
        
        function H_hr = Htrig2HR(H_htrig,options)
            arguments
                H_htrig Htrig;
                options.POSCAR = 'POSCAR';
            end
            switch H_htrig.Type
                case 'sincos'
                    Htrig_exp = H_htrig.rewrite();
                case 'exp'
                    Htrig_exp = H_htrig;
                otherwise

            end
            %             if Htrig_sincos.Type ~= "sincos"
            %                 warning('wrong type of trig class')
            %             end
            %%
            Hexp = Htrig_exp.HcoeL;
            WAN_NUM = Htrig_exp.Basis_num;
            %% empty HR class
            H_hr = HR(WAN_NUM,'Dim',Htrig_exp.Dim);
            H_hr = H_hr.vasplibCopy(Htrig_exp);            
            if isempty(H_hr.orbL)
                warning(['use ',options.POSCAR,'enforcely, please check it']);
                H_hr = H_hr < options.POSCAR;
            end
            H_hr = H_hr.tjmti_gen();
            tji_mat_r = H_hr.tjmti{1};
            %%
            VarsUsing = H_htrig.VarsSeqLcart(1:H_htrig.Dim);
            %syms k_x k_y k_z real
            hsym = Htrig_exp.HsymL;
            RM = H_hr.Rm;
            DIM = H_hr.Dim;
            for n = 1:length(hsym)
                if isequal(hsym(n), sym(1))
                    kd_num = zeros(1,DIM);
                else
                    %ikd = children(hsym(n));
                    [ChirdrenCell,Type] = vasplib.fixedchildren(hsym(n),'exp_inner');
                    if strcmp(Type,'sum') || strcmp(Type,'prod') % combine
                        error('? something wrong? You should Debug here!');
                    elseif strcmp(Type,'inner')
                        kd = ChirdrenCell{1}/1i;
                    end
                    for d = 1:DIM
                        cDim{d} = subs(kd,VarsUsing(d),1) - subs(kd,VarsUsing(d),0);
                    end
                    kd_num = double(fold(@horzcat,cDim));
                end
                % Factorlist_R $H_{i j}^{\mathbf{k}}=\left\langle\chi_{i}^{\mathbf{k}}|H| \chi_{j}^{\mathbf{k}}\right\rangle=\sum_{\mathbf{R}} e^{i \mathbf{k} \cdot\left(\mathbf{R}+\mathbf{t}_{j}-\mathbf{t}_{i}\right)} H_{i j}(\mathbf{R})$
                % exp(i(R+tj-ti))
                for i = 1:WAN_NUM
                    for j = 1:WAN_NUM
                        if isequal(Hexp(i,j,n), sym(0))
                            continue               
                        end
                        kji_num = reshape(tji_mat_r(i,j,:),[1,DIM]);
                        kr =   kd_num-kji_num;
                        vector = round(kr/RM);
                        H_hr = H_hr.set_hop((Hexp(i,j,n)),i,j,vector,'symadd');
                    end
                end
            end
        end
        
        function H_hk = Htrig2HK(H_htrig,kpoints_f,options)
            arguments
                H_htrig Htrig;
                kpoints_f = [0,0,0];
                options.sym =false;
                options.Order = 1;
            end
            syms k_x k_y k_z real;
            H_hk = HK(H_htrig.Basis_num,options.Order);
            H_hk = H_hk.vasplibCopy(H_htrig);
            if  options.sym
                H_hk.Rm = sym(H_hk.Rm);
                kpoints_r = kpoints_f * H_hk.Gk;
                fprintf('please check whether the kpoint(cartesian) is properly set.');
                disp(kpoints_r);
            else
                kpoints_r = kpoints_f * H_hk.Gk;
            end
            pb = vasplib_tool_outer.CmdLineProgressBar('Transforming ');
            for i = 1:H_htrig.Kinds
                pb.print(i,H_htrig.Kinds,' th HsymL_trig into HK obj ...');
                tmp_HsymL_trig = subs(H_htrig.HsymL_trig(i),...
                    [k_x k_y k_z],...
                    [k_x-kpoints_r(1),k_y-kpoints_r(2),k_z-kpoints_r(3)]);
                symbolic_polynomial = taylor(tmp_HsymL_trig,[k_x k_y k_z],'Order',options.Order+1);
                H_hk = H_hk.setup_rough(symbolic_polynomial,H_htrig.HcoeL(:,:,i),true);
            end
            pb.delete();
        end
    end
    
    %% input and output like function   
    methods
        function varargout = EIGENCAR_gen(H_htrig,options)
            arguments
                H_htrig Htrig;
                options.fermi double = 0;
                options.norb double = -1;
                options.klist double = H_htrig.klist_cart;
                options.para  = [];
                options.paraname ;
                options.show = false;
                options.ax = handle([]);
                options.printmode = true;
            end
            %
            switch H_htrig.Type
                case 'slab'
                    fprintf('use slab eigencar gen or ?\n');
                    return;
                otherwise
            end
            % -------------- nargin ------------------
            fermi = options.fermi;
            norb_enforce  = options.norb;
            if isempty(options.klist)
                H_htrig = H_htrig.kpathgen3D('KPOINTS');
                klist_cart_tmp = H_htrig.klist_cart;
            else
                klist_cart_tmp = options.klist;
            end
            if options.show
                if isempty(options.ax)
                    %[fig,ax] = create_figure();
                    ax = vasplib.BZplot(H_htrig.Rm,'color','r');
                else
                    ax = options.ax;
                end
            end
            % -------------- nargin ------------------
            %disp("EIGENCAR gen for H_xyz(wt TB) Type: HR class ");
            Hnum_list = H_htrig.HnumL ;
            if isempty(options.para)
                %             Bassis_mat = H_htrig.Bassis_mat ;
                if H_htrig.Basis_num > 500
                    print_mode = 1;
                else
                    print_mode = 0;
                end
                [kn,~] = size(klist_cart_tmp);
                %--------  check  --------
                if norb_enforce <0
                    NBANDS=H_htrig.Basis_num;
                elseif norb_enforce >0
                    NBANDS=norb_enforce;
                else
                    
                end
                if strcmp(H_htrig.Type,'list')
                    H_htrig=H_htrig.SliceGen();
                end
                WAVECAR  = zeros(H_htrig.Basis_num,NBANDS,kn);
                EIGENCAR = zeros(NBANDS,kn);
                sizemesh = [H_htrig.Basis_num,H_htrig.Basis_num];
                if options.printmode
                    pb = vasplib_tool_outer.CmdLineProgressBar('BAND calculating ');
                end
                try
                    %syms k_x k_y k_z real;
                    HsymL_fun = (matlabFunction( H_htrig.HsymL,'Vars',H_htrig.VarsSeqLcart(1:H_htrig.Dim)));
                catch 
                    
                end
                for ki =1:kn
                    %k_x=klist_cart_tmp(ki,1);
                    %k_y=klist_cart_tmp(ki,2);
                    %k_z=klist_cart_tmp(ki,3);
                    switch H_htrig.Type
                        case 'sparse'
                            Htemp=sparse(H_htrig.Basis_num ,H_htrig.Basis_num);
                            for i=1:H_htrig.Kinds
                                Htemp = Htemp +Hnum_list{i}*double(H_htrig.HsymL_trig(i));
                            end
                            Hout = Htemp;
                        case {'mat'}
                            Factorlist = exp(1i*H_htrig.HsymL_numL*klist_cart_tmp(ki,:).');
                            Hout = sum(vasplib.matrixtimespage(Factorlist,H_htrig.HnumL),3);
                            Hout = (Hout+Hout')/2;
                        case 'list'
                            Hout = zeros(sizemesh);
                            Hnum_list_k = H_htrig.HnumL.*exp(1i*H_htrig.HsymL_numL(:,1:H_htrig.Dim)*klist_cart_tmp(ki,:).');
                            for i=1:H_htrig.N_Sparse_vector
                                Hout(H_htrig.Sparse_vector(i,1),H_htrig.Sparse_vector(i,2)) = sum(Hnum_list_k(H_htrig.CutList(i,1):H_htrig.CutList(i,2)));
                            end
                            %Factorlist = H_htrig.HnumL.*exp(1i*H_htrig.HsymL_numL(:,1:3)*klist_cart_tmp(ki,:).');
                            %[ij_unique,sumFactorlist] = HollowKnight.generalcontractrow(H_htrig.HsymL_numL(:,4:5),Factorlist);
                            %indL = sub2ind(sizemesh,ij_unique(:,1),ij_unique(:,2));
                            %Hout(indL) = sumFactorlist;
                            Hout = (Hout+Hout')/2;
                        case 'sincos'
                            Input = num2cell(klist_cart_tmp(ki,:)); 
                            kL = HsymL_fun(Input{:});
                            Hout = H_htrig.HnumL;
                            for i =1:H_htrig.Kinds
                                Hout(:,:,i) = Hout(:,:,i).*kL(:,i);
                            end
                            Hout = sum(Hout,3);
                            %                 end
                            Hout = (Hout+Hout')/2;
                        otherwise
                            Input = num3cell(klist_cart_tmp(ki,:)); 
                            Hout = H_htrig.Hfun(Input{:});
                            Hout = (Hout+Hout')/2;
                    end
                    if norb_enforce <0
                        try
                            [A, U]=eig(full(Hout));
                        catch
                            disp([ki,klist_cart_tmp(ki,:)] );
                            disp(Hout);
                            disp(H_htrig.Hfun);
                            error('check this k point');
                        end
                    elseif norb_enforce >0
                        [A, U]=eigs(Hout,NBANDS,fermi);
                        [A, U]=park.sorteig(U,A);
                    else
                    end
                    EIGENCAR(:,ki) = diag(U);
                    WAVECAR(:,:,ki) = A;
                    if options.printmode
                        pb.print(ki,kn,' ...');
                    end
                end
                if options.printmode
                    pb.delete;
                end
            else
                Npara = size(options.para ,1);
                paraN = size(options.para ,2);
                %             Bassis_mat = H_htrig.Bassis_mat ;
                if H_htrig.Basis_num > 500
                    print_mode = 1;
                else
                    print_mode = 0;
                end
                [kn,~] = size(klist_cart_tmp);
                %--------  check  --------
                if norb_enforce <0
                    NBANDS=H_htrig.Basis_num;
                elseif norb_enforce >0
                    NBANDS=norb_enforce;
                else
                    
                end
                WAVECAR  = [];
                EIGENCAR{Npara} = zeros(NBANDS,kn);
                Htrig_sym = H_htrig.Hsym;
                for j = 1:Npara
                    fprintf('**************************************************************************************\n');
                    for i = 1:paraN
                        fprintf('%s :',mat2str(string(sym(options.paraname(i)))));
                        fprintf('%f\n',options.para(j,i));
                    end
                   % fprintf('**************************************************************************************\n');
                    EIGENCAR_tmp = zeros(NBANDS,kn);
                    %                     for n = 1:paraN
                    %
                    %                     end
                    if strcmp(H_htrig.Type,'sparse')
                        for i = 1:H_htrig.Kinds
                            Hnum_list{i} = subs(H_htrig.HnumL{i},sym(options.paraname),options.para(j,:));
                        end
                    else
                        H_fun_tmp = matlabFunction(subs(Htrig_sym,sym(options.paraname),options.para(j,:)),'Vars',sym(["k_x","k_y","k_z"]));
                    end
                    
                    for ki =1:kn
                        k_x=klist_cart_tmp(ki,1);
                        k_y=klist_cart_tmp(ki,2);
                        k_z=klist_cart_tmp(ki,3);
                        if strcmp(H_htrig.Type,'sparse')
                            Htemp=sparse(H_htrig.Basis_num ,H_htrig.Basis_num);
                            for i=1:H_htrig.Kinds
                                Htemp = Htemp +Hnum_list{i}*double(H_htrig.HsymL_trig(i));
                            end
                            Hout = Htemp;
                        else
                            Hout = H_fun_tmp(k_x,k_y,k_z);
                            Hout = (Hout+Hout')/2;
                        end
                        if norb_enforce <0
                            try
                                [A, U]=eig(full(Hout));
                            catch
                                disp([ki,k_x,k_y,k_z] );
                                disp(Hout);
                                disp(H_htrig.Hfun);
                                error('check this k point');
                            end
                        elseif norb_enforce >0
                            [A, U]=eigs(Hout,NBANDS,fermi);
                            [A, U]=park.sorteig(U,A);
                        else
                        end
                        EIGENCAR_tmp(:,ki) = diag(U);
                        if print_mode ==1
                            fprintf('%d th kpoints has been caculated in %d kpoints total\n',ki,kn);
                        end
                    end
                    EIGENCAR{j} = EIGENCAR_tmp;
                end
            end
            %             if kn >1
            %                 EIGENCARout = EIGENCAR;
            %                 %EIGENCARout{2} = WAVECAR;
            %             else
            %                 EIGENCARout.EIGENCAR= EIGENCAR;
            %                 EIGENCARout.WAVECAR = WAVECAR;
            %             end
            varargout{1} = EIGENCAR;
            varargout{2} = WAVECAR;
            if options.show
                [varargout{3}] = vasplib.klist_show(...
                    'klist',klist_cart_tmp,...
                    'ax',ax);
            end
        end
        function [EIGENCAR,WAVECAR,WEIGHTCAR] = EIGENCAR_gen_slab(H_htrig,options)
            arguments
                H_htrig Htrig;
                options.fermi double = 0;
                options.norb double = -1;
                options.klist  = H_htrig.klist_cart;
                options.para  = [];
                options.paraname ;
            end
            % 
            switch H_htrig.Type
                case 'sincos'
                    fprintf('use normal eigencar gen or ?\n');
                    return;
                otherwise
            end
            % -------------- nargin ------------------
            fermi = options.fermi;
            norb_enforce  = options.norb;
            if isempty(H_htrig.klist_cart)
                if isstr(options.klist)
                    H_htrig = H_htrig.kpathgen3D(options.klist);
                else
                    H_htrig = H_htrig.kpathgen3D('KPOINTS_slab');
                end 
                options.klist = H_htrig.klist_cart;
            end
            klist_cart_tmp = options.klist;
            % -------------- nargin ------------------
            %disp("EIGENCAR gen for H_xyz(wt TB) Type: HR class ");
            HcoeList = H_htrig.HcoeL ;
            VarUsing = H_htrig.VarsSeqLcart(1:H_htrig.Dim);
            NSLAB = (H_htrig.Nslab ==0) + H_htrig.Nslab;
            NS = NSLAB(1)* NSLAB(2)* NSLAB(3);
            NWAVE = H_htrig.Basis_num* NS;
            orb_list = H_htrig.orbL;
            
            if H_htrig.Nslab(1) > 1
                signlist = sign((orb_list(:,1)-0.5));
                dir = 1;
            elseif H_htrig.Nslab(2) > 1
                dir = 2;
                signlist = sign((orb_list(:,2)-0.5));
            elseif H_htrig.Nslab(3) > 1
                dir = 3;
                signlist = sign((orb_list(:,3)-0.5));
            elseif H_htrig.Nslab(4) > 1
                dir = 4;
                signlist = sign((orb_list(:,4)-0.5));
            end
            HSVCAR_slab = vasplib.HSVCAR_gen(orb_list,'slab',0.05,[0.5,0.5,0.5],dir);
            signlist(HSVCAR_slab(:,1) == 0) = 0;
            if size(H_htrig.orbL,1) > 1000
                print_mode = true;
            else
                print_mode = false;
            end
            if isempty(options.para)
                %             Bassis_mat = H_htrig.Bassis_mat ;

                [kn,~] = size(klist_cart_tmp);
                %--------  check  --------
                if norb_enforce <0
                    NBANDS=NWAVE;
                elseif norb_enforce >0
                    NBANDS=norb_enforce;
                else
                    
                end
                
                WAVECAR  = zeros(NWAVE,NBANDS,kn);
                EIGENCAR = zeros(NBANDS,kn);
                WEIGHTCAR = zeros(NBANDS,kn);
                for ki =1:kn
                    k_d=klist_cart_tmp(ki,:);
                    Input = num2cell(k_d);
                    Hmat = zeros(NWAVE,NWAVE,numel(H_htrig.HsymL_trig));
                    for i = 1:numel(H_htrig.HsymL_trig)
                        try
                            Hfuntemp = matlabFunction(HcoeList(:,:,i),'Vars',VarUsing); %?
                        catch
                            error('You have not subs the function');
                        end
                        %disp(Hfuntemp)
                        %disp(Hfuntemp(k_x,k_y,k_z))
                        %disp(Hmat_pre(:,:,i));
                        Hmat(:,:,i) = kron(H_htrig.Hmat_pre{i},Hfuntemp(Input{:}));
                    end
                    Hout = sum(Hmat,3);
                    %disp(tril(Hout,-1));
                    Hout = tril(Hout,-1)+diag(real(diag((Hout))))+tril(Hout,-1)';

                    if norb_enforce <0
                        try
                            [A, U]=eig(full(Hout));
                        catch
                            disp([ki,Input] );
                            disp(Hout);
                            disp(H_htrig.Hfun);
                            error('check this k point');
                        end
                    elseif norb_enforce >0
                        [A, U]=eigs(Hout,NBANDS,fermi);
                        [A, U]=park.sorteig(U,A);
                    else
                    end
                    EIGENCAR(:,ki) = diag(U);
                    WAVECAR(:,:,ki) = A;
                    [~,WEIGHTCAR(:,ki)] = Htrig.COLORCAR_gen(A,HSVCAR_slab,signlist);
                    if print_mode ==1
                        fprintf('%d th kpoints has been caculated in %d kpoints total\n',ki,kn);
                    end
                end
                
            else
                Npara = size(options.para ,1);
                paraN = size(options.para ,2);
                %             Bassis_mat = H_htrig.Bassis_mat ;
                [kn,~] = size(klist_cart_tmp);
                %--------  check  --------
                if norb_enforce <0
                    NBANDS=NWAVE;
                elseif norb_enforce >0
                    NBANDS=norb_enforce;
                else
                    
                end
                WEIGHTCAR = [];
%                 Hmat_pre = zeros(NS,NS,numel(H_htrig.HsymL_trig));
%                 for i = 1:numel(H_htrig.HsymL_trig)
%                     Hmat_pre(:,:,i) = H_htrig.HsymL_trig2mat(H_htrig.HsymL_trig(i));
%                 end
                WAVECAR  = [];
                EIGENCAR{Npara} = zeros(NBANDS,kn);
                %Htrig_num_tmp = H_htrig.Htrig_num;
                for j = 1:Npara
                    fprintf('**************************************************************************************\n');
                    for i = 1:paraN
                        fprintf('%s :',mat2str(string(sym(options.paraname(i)))));
                        fprintf('%f\n',options.para(j,i));
                    end
                   % fprintf('**************************************************************************************\n');
                    EIGENCAR_tmp = zeros(NBANDS,kn);
                    %                     for n = 1:paraN
                    %
                    %                     end
                    if strcmp(H_htrig.Type,'sparse')
                        for i = 1:H_htrig.Kinds
                            HcoeList{i} = subs(H_htrig.HnumL{i},sym(options.paraname),options.para(j,:));
                        end
                    else
                        
                        for i =1:numel( H_htrig.HsymL_trig)
                            
                            % sorry i dont know how to improve
                              %H_fun_t{i} =@(k_x,k_y,k_z) H_htrig.Hfun{i}(k_x,k_y,k_z,([options.para(j,:)]));
                             temp_str = ["H_fun_t{i} =@(k_x,k_y,k_z) H_htrig.Hfun{i}(k_x,k_y,k_z",string(options.para(j,:))];
                             temp_str = strjoin(temp_str,',');
                             temp_str = temp_str+");";
                             eval(temp_str);
                        end
%                         H_fun_tmp =  H_htrig.Hfun;
%                         H_fun_tmp = matlabFunction(subs(Htrig_num_tmp,sym(options.paraname),options.para(j,:)),'Vars',sym(["k_x","k_y","k_z"]));
                    end
                    
                    for ki =1:kn
                        k_x=klist_cart_tmp(ki,1);
                        k_y=klist_cart_tmp(ki,2);
                        k_z=klist_cart_tmp(ki,3);
                        
                        if strcmp(H_htrig.Type,'sparse')
                            Htemp=sparse(H_htrig.Basis_num ,H_htrig.Basis_num);
                            for i=1:H_htrig.Kinds
                                Htemp = Htemp +HcoeList{i}*double(H_htrig.HsymL_trig(i));
                            end
                            Hout = Htemp;
                        else
                            Hmat = zeros(NWAVE,NWAVE,numel(H_htrig.HsymL_trig));
                            for i = 1:numel(H_htrig.HsymL_trig)
                                try
                                    Hfuntemp = H_fun_t{i};
                                catch
                                    error('You have not subs the function');
                                end
                                %disp(Hfuntemp)
                                %disp(Hfuntemp(k_x,k_y,k_z))
                                %disp(Hmat_pre(:,:,i));
                                Hmat(:,:,i) = kron(H_htrig.Hmat_pre{i},Hfuntemp(k_x,k_y,k_z));
                            end
                            Hout = sum(Hmat,3);
                            %disp(tril(Hout,-1));
                            Hout = tril(Hout,-1)+diag(real(diag((Hout))))+tril(Hout,-1)';
                        end
                        if norb_enforce <0
                            try
                                [A, U]=eig(full(Hout));
                            catch
                                disp([ki,k_x,k_y,k_z] );
                                disp(Hout);
                                disp(H_htrig.Hfun);
                                error('check this k point');
                            end
                        elseif norb_enforce >0
                            try
                                [A, U]=eigs(Hout,NBANDS,fermi);
                                [A, U]=park.sorteig(U,A);
                            catch
                                [A, U]=eig(Hout,'vector');
                                U = diag(U(NWAVE/2-NBANDS/2+1:NWAVE/2+NBANDS/2));
                            end
                        else
                        end
                        EIGENCAR_tmp(:,ki) = diag(U);
                        %WAVECAR(:,:,ki) = A;
                        if print_mode ==1
                            fprintf('%d th kpoints has been caculated in %d kpoints total\n',ki,kn);
                        end
                    end
                    EIGENCAR{j} = EIGENCAR_tmp;
                    
                end
            end
            %             if kn >1
            %                 EIGENCARout = EIGENCAR;
            %                 %EIGENCARout{2} = WAVECAR;
            %             else
            %                 EIGENCARout.EIGENCAR= EIGENCAR;
            %                 EIGENCARout.WAVECAR = WAVECAR;
            %             end
        end
        function [EIGENCAR,WAVECAR,WEIGHTCAR] = EIGENCAR_gen_wire(H_htrig,options)
            arguments
                H_htrig Htrig;
                options.fermi double = 0;
                options.norb double = -1;
                options.klist  = H_htrig.klist_cart;
                options.para  = [];
                options.paraname ;
                options.issym  logical =false;
                options.spdB double = 0;
            end
            % 
            switch H_htrig.Type
                case 'sincos'
                    fprintf('use normal eigencar gen or ?\n');
                    return;
                otherwise
            end
            % -------------- nargin ------------------
            fermi = options.fermi;
            norb_enforce  = options.norb;
            if isempty(H_htrig.klist_cart) && ~isnumeric(options.klist)
                if isstr(options.klist)
                    H_htrig = H_htrig.kpathgen3D(options.klist);
                else
                    try
                        H_htrig = H_htrig.kpathgen3D('KPOINTS_wire');
                    catch
                        
                    end
                end 
                options.klist = H_htrig.klist_cart;
            end
            klist_cart_tmp = options.klist;
            % -------------- nargin ------------------
            %disp("EIGENCAR gen for H_xyz(wt TB) Type: HR class ");
            Hnum_list = H_htrig.HnumL ;
            NSLAB = (H_htrig.Nslab ==0) + H_htrig.Nslab;
            NS = prod(NSLAB);
            NWAVE_origin =  H_htrig.Basis_num* NS;
            NWAVE = H_htrig.Basis_num* NS-sum(H_htrig.rm_list);
            orb_list  = H_htrig.orbL;
            HcoeList = H_htrig.HcoeL ;
            VarUsing = H_htrig.VarsSeqLcart(1:H_htrig.Dim);
            if isempty(H_htrig.rm_list)
                HSVCAR_hinge = vasplib.HSVCAR_gen(orb_list,'hinge',0.05,[0.5,0.5,0.5],3);
            else
                HSVCAR_hinge = vasplib.HSVCAR_gen(orb_list,'surf',0.05,[0.5,0.5,0.5],3);
            end
            signlist = sign((orb_list(:,1)-0.5).*(orb_list(:,2)-0.5));
            if isempty(options.para)
                %             Bassis_mat = H_htrig.Bassis_mat ;
                if H_htrig.Basis_num > 1000
                    print_mode = 1;
                else
                    print_mode = 0;
                end
                [kn,~] = size(klist_cart_tmp);
                %--------  check  --------
                if norb_enforce <0
                    NBANDS=NWAVE;
                elseif norb_enforce >0
                    NBANDS=norb_enforce;
                else
                    
                end
                
                WAVECAR  = zeros(NWAVE,NBANDS,kn);
                EIGENCAR = zeros(NBANDS,kn);
                WEIGHTCAR = zeros(NBANDS,kn);
                for ki =1:kn
                    k_d=klist_cart_tmp(ki,:);
                    Input = num2cell(k_d);
                    Hmat{numel(H_htrig.HsymL_trig)} = sparse(NWAVE_origin,NWAVE_origin);
                    for i = 1:numel(H_htrig.HsymL_trig)
                        try
                            Hfuntemp = matlabFunction(HcoeList(:,:,i),'Vars',VarUsing); %?
                        catch
                            error('You have not subs the function');
                        end
                        Hmat{i} = kron(H_htrig.Hmat_pre{i},Hfuntemp(Input{:}));
                    end
                    Hout = fold(@plus,Hmat);
                    %disp(tril(Hout,-1));
                    Hout = tril(Hout,-1)+diag(real(diag((Hout))))+tril(Hout,-1)';

                    if ~isempty(H_htrig.rm_list)
                        Hout(H_htrig.rm_list,:) = [];
                        Hout(:,H_htrig.rm_list) = [];
                    end
                    if norb_enforce <0
                        try
                            [A, U]=eig(full(Hout));
                        catch
                            disp([ki,Input] );
                            disp(Hout);
                            disp(H_htrig.Hfun);
                            error('check this k point');
                        end
                    elseif norb_enforce >0
                        if options.issym
                            Hout = Hout+eye()*options.spdB;
                        end
                        [A, U]=eigs(Hout,NBANDS,options.spdB+fermi+1e-6,'IsSymmetricDefinite',options.issym);
                        [A, U]=park.sorteig(U,A);
                    else
                    end
                    EIGENCAR(:,ki) = diag(U);
                    WAVECAR(:,:,ki) = A;
                    [~,WEIGHTCAR(:,ki)] = Htrig.COLORCAR_gen(A,HSVCAR_hinge,signlist);
                    if print_mode ==1
                        fprintf('%d th kpoints has been caculated in %d kpoints total\n',ki,kn);
                    end
                end
                %WEIGHTCAR = normalize(WEIGHTCAR,2,'range');
                
            else
                Npara = size(options.para ,1);
                paraN = size(options.para ,2);
                %             Bassis_mat = H_htrig.Bassis_mat ;
                if H_htrig.Basis_num > 500
                    print_mode = 1;
                else
                    print_mode = 0;
                end
                [kn,~] = size(klist_cart_tmp);
                if kn ==1
                    tmpmode = 'para';
                else
                    tmpmode = 'wire';
                end
                %--------  check  --------
                if norb_enforce <0
                    NBANDS=NWAVE;
                elseif norb_enforce >0
                    NBANDS=norb_enforce;
                else
                    
                end
                

                
                %Htrig_num_tmp = H_htrig.Htrig_num;
                if strcmp(tmpmode,'wire')
                    % abandon WAVECAR
                    WAVECAR  = [];
                    % 
                    EIGENCAR{Npara} = zeros(NBANDS,kn);
                    WEIGHTCAR{Npara} = zeros(NBANDS,kn);
                    for j = 1:Npara
                        fprintf('**************************************************************************************\n');
                        for i = 1:paraN
                            fprintf('%s :',mat2str(string(sym(options.paraname(i)))));
                            fprintf('%f\n',options.para(j,i));
                        end
                        % fprintf('**************************************************************************************\n');
                        EIGENCAR_tmp = zeros(NBANDS,kn);
                        WEIGHTCAR_tmp = zeros(NBANDS,kn);
                        %                     for n = 1:paraN
                        %
                        %                     end
                        if strcmp(H_htrig.Type,'sparse')
                            for i = 1:H_htrig.Kinds
                                Hnum_list{i} = subs(H_htrig.HnumL{i},sym(options.paraname),options.para(j,:));
                            end
                        else
                            for i =1:numel( H_htrig.HsymL_trig)
                                % sorry i dont know how to improve
                                %H_fun_t{i} =@(k_x,k_y,k_z) H_htrig.Hfun{i}(k_x,k_y,k_z,([options.para(j,:)]));
                                temp_str = ["H_fun_t{i} =@(k_x,k_y,k_z) H_htrig.HcoeL(:,:,i)(k_x,k_y,k_z",string(options.para(j,:))];
                                temp_str = strjoin(temp_str,',');
                                temp_str = temp_str+");";
                                eval(temp_str);
                            end
                            %                         H_fun_tmp =  H_htrig.Hfun;
                            %                         H_fun_tmp = matlabFunction(subs(Htrig_num_tmp,sym(options.paraname),options.para(j,:)),'Vars',sym(["k_x","k_y","k_z"]));
                        end
                        
                        for ki =1:kn
                            k_x=klist_cart_tmp(ki,1);
                            k_y=klist_cart_tmp(ki,2);
                            k_z=klist_cart_tmp(ki,3);
                            
                            if strcmp(H_htrig.Type,'sparse')
                                Htemp=sparse(H_htrig.Basis_num ,H_htrig.Basis_num);
                                for i=1:H_htrig.Kinds
                                    Htemp = Htemp +Hnum_list{i}*double(H_htrig.HsymL_trig(i));
                                end
                                Hout = Htemp;
                            else
                                Hmat{numel(H_htrig.HsymL_trig)} = sparse(NWAVE_origin,NWAVE_origin);
                                for i = 1:numel(H_htrig.HsymL_trig)
                                    try
                                        Hfuntemp = H_fun_t{i};
                                    catch
                                        error('You have not subs the function');
                                    end
                                    %disp(Hfuntemp)
                                    %disp(Hfuntemp(k_x,k_y,k_z))
                                    %disp(Hmat_pre(:,:,i));
                                    Hmat{i} = kron(H_htrig.Hmat_pre{i},Hfuntemp(k_x,k_y,k_z));
                                end
                                Hout = Hmat{1};
                                for i = 2:numel(H_htrig.HsymL_trig)
                                    Hout =Hout +Hmat{i};
                                end
                                %disp(tril(Hout,-1));
                                Hout = tril(Hout,-1)+diag(real(diag((Hout))))+tril(Hout,-1)';
                            end
                            if ~isempty(H_htrig.rm_list)
                                Hout(H_htrig.rm_list,:) = [];
                                Hout(:,H_htrig.rm_list) = [];
                            end
                            if norb_enforce <0
                                try
                                    [A, U]=eig(full(Hout));
                                catch
                                    disp([ki,k_x,k_y,k_z] );
                                    disp(Hout);
                                    disp(H_htrig.Hfun);
                                    error('check this k point');
                                end
                            elseif norb_enforce >0
                                try
                                    if options.issym
                                        Hout = Hout+eye()*options.spdB;
                                    end
                                    [A, U]=eigs(Hout,NBANDS,options.spdB+fermi+1e-6,'IsSymmetricDefinite',options.issym);
                                    [A, U]=park.sorteig(U,A);
                                catch
                                    [A, U]=eig(Hout,'vector');
                                    U = diag(U(NWAVE/2-NBANDS/2+1:NWAVE/2+NBANDS/2));
                                end
                            else
                            end
                            EIGENCAR_tmp(:,ki) = diag(U);
                            [~,WEIGHTCAR_tmp(:,ki)] = Htrig.COLORCAR_gen(A,HSVCAR_hinge,signlist);
                            %WAVECAR(:,:,ki) = A;
                            if print_mode ==1
                                fprintf('%d th kpoints has been caculated in %d kpoints total\n',ki,kn);
                            end
                        end
                        EIGENCAR{j} = EIGENCAR_tmp;
                        %WEIGHTCAR_tmp = normalize(WEIGHTCAR_tmp,2,'range');
                        WEIGHTCAR{para} = WEIGHTCAR_tmp;
                    end
                else
                    WAVECAR  = zeros(NWAVE,NBANDS,Npara);
                    EIGENCAR = zeros(NBANDS,Npara);
                    WEIGHTCAR = zeros(NBANDS,Npara);
                    k_x=klist_cart_tmp(1);
                    k_y=klist_cart_tmp(2);
                    k_z=klist_cart_tmp(3);
                    for j = 1:Npara
                        for i = 1:paraN
                            fprintf('%s :',mat2str(string(sym(options.paraname(i)))));
                            fprintf('%f\n',options.para(j,i));
                        end
                        TheHcoeL = H_htrig.HcoeL;
                        TheHcoeL =  subs(TheHcoeL,sym(options.paraname),options.para(j,:));
                        for i =1:numel( H_htrig.HsymL_trig)
                            H_fun_t{i} = matlabFunction(TheHcoeL(:,:,i),'Vars',H_htrig.VarsSeqLcart(1:H_htrig.Dim));
                        end
                        if strcmp(H_htrig.Type,'sparse')
                            Htemp=sparse(H_htrig.Basis_num ,H_htrig.Basis_num);
                            for i=1:H_htrig.Kinds
                                Htemp = Htemp +Hnum_list{i}*double(H_htrig.HsymL_trig(i));
                            end
                            Hout = Htemp;
                        else
                            Hmat{numel(H_htrig.HsymL_trig)} = sparse(NWAVE_origin,NWAVE_origin);
                            for i = 1:numel(H_htrig.HsymL_trig)
                                try
                                    Hfuntemp = H_fun_t{i};
                                catch
                                    error('You have not subs the function');
                                end
                                %disp(Hfuntemp)
                                %disp(Hfuntemp(k_x,k_y,k_z))
                                %disp(Hmat_pre(:,:,i));
                                Hmat{i} = kron(H_htrig.Hmat_pre{i},Hfuntemp(k_x,k_y,k_z));
                            end
                            
                            Hout = Hmat{1};
                            for i = 2:numel(H_htrig.HsymL_trig)
                                Hout =Hout +Hmat{i};
                            end
                            %disp(tril(Hout,-1));
                            Hout = tril(Hout,-1)+diag(real(diag((Hout))))+tril(Hout,-1)';
                        end
                        if ~isempty(H_htrig.rm_list)
                            Hout(H_htrig.rm_list,:) = [];
                            Hout(:,H_htrig.rm_list) = [];
                        end
                        if norb_enforce <0
                            try
                                [A, U]=eig(full(Hout));
                            catch
                                disp([k_x,k_y,k_z] );
                                disp(Hout);
                                disp(H_htrig.Hfun);
                                error('check this k point');
                            end
                        elseif norb_enforce >0
                            if options.issym
                                Hout = Hout+eye()*options.spdB;
                            end
                            [A, U]=eigs(Hout,NBANDS,options.spdB+fermi+1e-6,'IsSymmetricDefinite',options.issym);
                            [A, U]= park.sorteig(U,A);
                        else
                        end
                        EIGENCAR(:,j) = diag(U);
                        [~,WEIGHTCAR(:,j)] = Htrig.COLORCAR_gen(A,HSVCAR_hinge,signlist);
                        WAVECAR(:,:,j) = A;
                        if print_mode ==1
                            fprintf('%d th para has been caculated in %d Npara total\n',j,Npara);
                        end
                    end
                    %WEIGHTCAR = normalize(WEIGHTCAR,2,'range');
                    %WEIGHTCAR = WEIGHTCAR;
                end
            end
            %             if kn >1
            %                 EIGENCARout = EIGENCAR;
            %                 %EIGENCARout{2} = WAVECAR;
            %             else
            %                 EIGENCARout.EIGENCAR= EIGENCAR;
            %                 EIGENCARout.WAVECAR = WAVECAR;
            %             end
        end
        function H_htrig = Subsall(H_htrig,mode,AddtionVar)
            arguments
                H_htrig
                mode {mustBeMember(mode,{'num','file','gen','para','sym','slab','disk'})}= 'num';
                AddtionVar = sym([]);
            end
            if nargin <2
                mode = 'num';
            end
            if nargin < 3
                AddtionVar = sym([]);
            end
            % ------ init ------- load para form inner or outer
            if exist('para.mat','file') && strcmp(mode,'file')
                %disp('The para in para.mat is first to be considered?');
                load('para.mat');
            else

            end

            if ~isempty(AddtionVar)
                varlist = [ H_htrig.VarsSeqLcart(1:H_htrig.Dim),AddtionVar];%sym(AddtionVar);%
            else
                varlist = H_htrig.VarsSeqLcart(1:H_htrig.Dim);%sym([]);
            end
            
            H_htrig.HcoeL = subs( H_htrig.HcoeL);
            H_htrig.HsymL_trig = subs(H_htrig.HsymL_trig);
            H_htrig.HsymL_coeL  = subs(H_htrig.HsymL_coeL);
            H_htrig.HsymL_trig_bk = subs(H_htrig.HsymL_trig_bk);
            H_htrig.Htrig_num  = subs(H_htrig.Htrig_sym);
            H_htrig.Hsym = H_htrig.Htrig_num;
            SymvarListInHcoeL = symvar(H_htrig.HcoeL);
            switch mode
                case {'num','file','gen','para'}
                    HcoeL_temp = H_htrig.HcoeL;
                    if ~isempty(SymvarListInHcoeL)
                        for i = 1:length(SymvarListInHcoeL)
                            fprintf('this varible: %s is the extra parameter\n',string(SymvarListInHcoeL(i)));
                            if isempty(find(varlist == SymvarListInHcoeL(i), 1))
                                try
                                    HcoeL_temp = subs(HcoeL_temp,SymvarListInHcoeL(i),evalin('base',string(SymvarListInHcoeL(i))));
                                catch
                                    warning('please subs this varible, Subsall must return numerical obj!\n');
                                end
                            end
                        end
                    else
                        H_htrig.HnumL = double(HcoeL_temp);
                    end 
                    if strcmp(mode,'para')
                        H_htrig.HcoeL = subs(H_htrig.HcoeL);
                        return;
                    else
                        H_htrig.HnumL = double(HcoeL_temp);
                        H_htrig.HcoeL = sym([]);
                    end
                    H_htrig.num = true;
                    H_htrig.coe = false;
                    if strcmp(H_htrig.Type,'mat') ||strcmp(H_htrig.Type,'list')
                        H_htrig.HsymL_numL  = double(H_htrig.HsymL_coeL);
                        H_htrig.HsymL_coeL = sym([]);
                        return;
                    end
                    %                     H_htrig.Htrig_sym = sym(zeros(H_htrig.Basis_num,H_htrig.Basis_num));
                    %                     for i =1:H_htrig.Kinds
                    %                         H_htrig.Htrig_sym = H_htrig.Htrig_sym + H_htrig.HcoeL(:,:,i)*H_htrig.HsymL_k(i);
                    %                     end
                    if strcmp(mode,'gen')
                        % H_htrig.Htrig_num  = subs(H_htrig.Htrig_sym);
                        % H_htrig.Hfun = matlabFunction(H_htrig.Htrig_num,'Vars',varlist,'File','H_htrig.m');
                    elseif ~strcmp(H_htrig.Type,'slab')
                        H_htrig.Htrig_num  = subs(H_htrig.Htrig_sym);
                        if H_htrig.Basis_num >500
                            H_htrig = sparse(H_htrig);
                        else
                            % H_htrig.Hfun = matlabFunction(H_htrig.Htrig_num,'Vars',varlist);
                        end
                        H_htrig.HcoeL =[];
                    else
                        if H_htrig.Basis_num >500
                            %H_htrig = sparse(H_htrig);
                        else

                        end
                    end
                case {'sym','slab','disk'}
                    H_htrig.HcoeL = subs(H_htrig.HcoeL);
            end
        end
        function mat = HsymL_trig2mat(H_htrig,HsymL_trig)
            %? When shall we use this function? I forgot!
            %Reply : For slab  
            NSLAB = (H_htrig.Nslab ==0)+H_htrig.Nslab;
            NS = prod(NSLAB);
            %NWAVE = H_htrig.Basis_num* NS;
            % We need to expand this function!
            Delta_ijk = Htrig.Delta_Oper(Htrig.coeff_extract(HsymL_trig));
            tmpmat = meshgrid((1:NS).',(1:NS).');
            NS1 = reshape(tmpmat,NS^2,1);
            NS2 = reshape(tmpmat.',NS^2,1);
            [i1,i2,i3] = ind2sub(NSLAB,NS1);
            [j1,j2,j3] = ind2sub(NSLAB,NS2);
            mat = double(reshape(Delta_ijk(i1,i2,i3,j1,j2,j3),NS,NS)).';
            %mat = double(reshape(Delta_ijk(i1,i2,i3,j1,j2,j3).',NS,NS)).';
            %             mat = zeros(NS);
            %             for i = 1:NS
            %                 [i1,i2,i3] = ind2sub(NSLAB,i);
            %                 for j = 1:NS
            %                     [j1,j2,j3] = ind2sub(NSLAB,j);
            %                     if Delta_ijk(i1,i2,i3,j1,j2,j3)
            %                         mat(i,j) = 1;
            %                     end
            %                 end
            %             end
        end
        function [num_label,coe_label,H_htrig] = NumOrCoe(H_htrig)
            if isempty(H_htrig.num) ||  isempty(H_htrig.coe)
                [num_label,coe_label] = NumOrCoe@vasplib(H_htrig);
                H_htrig.num = num_label;
                H_htrig.coe = coe_label;
            else
                num_label = H_htrig.num;
                coe_label = H_htrig.coe;
            end
        end
    end
    methods (Static)
        function H_htrig = HR2Htrig(H_hr)
            H_htrig =H_hr.HR2Htrig();
        end
        
        function Delta_Oper = Delta_Oper(Oper_list)
            tmpeq_all = "1";
            syms i1 i2 i3 j1 j2 j3 integer;
            tmpeqstr_list  =["(i1 == j1)","(i2 == j2)","(i3 == j3)"];
            for i = 1:length(Oper_list)
                tmpeq ="";
                strtmp = strsplit(string(Oper_list(i)),{'_'},'CollapseDelimiters',true);
                switch strtmp{2}
                    case 'x'
                        tmpeq = tmpeq+'i1 == (j1';
                    case 'y'
                        tmpeq = tmpeq+'i2 == (j2';
                    case 'z'
                        tmpeq = tmpeq+'i3 == (j3';
                    otherwise
                end
                if strcmp(strtmp{3},'N')
                    tmpeq = tmpeq + ")";
                else
                    tmpeq = tmpeq + " + "+strtmp{3}+")";
                end
                tmpeq = "(" +tmpeq+")";
                switch strtmp{2}
                    case 'x'
                        tmpeqstr_list(1) = tmpeq;
                    case 'y'
                        tmpeqstr_list(2) = tmpeq;
                    case 'z'
                        tmpeqstr_list(3) = tmpeq;
                    otherwise
                end
            end
            tmpeq_all = tmpeq_all + "&" +...
                tmpeqstr_list(1)+'&'+...
                tmpeqstr_list(2)+'&'+...
                tmpeqstr_list(3);
            Delta_Oper = matlabFunction(str2sym(tmpeq_all),'Vars',[i1 i2 i3 j1 j2 j3]);
        end
        % modify if 
        function sym_coe_list = coeff_extract(sym_term)
            str_list_tmp = strsplit(string(sym_term),'*');
            sym_coe = sym(1);
            for i = 1:length(str_list_tmp)
                if Htrig.strcontain(str_list_tmp{i},['x','y','z'])
                    sym_coe_list(i) = sym_coe * str2sym(str_list_tmp{i});
                end
            end
        end
    end

%% 
end