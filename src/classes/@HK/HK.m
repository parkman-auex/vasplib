%HK a powerful kp tools
classdef HK < vasplib & matlab.mixin.CustomDisplay
    properties
        Degree = 1;
        Kinds;
        HcoeL;
        HnumL;
        HstrL;
        HsymL = sym([]);
        Hk_num   ;
    end
    properties %(Hidden = true)
        % for spacegroup
        HsymL_xyz = sym([]);
        nn_store_smart   ; % nn_store
        nn_sparse_n      ;
        Atom_store_smart ;
        Rnn_map          ;
        Term_to_save     ;
        Trig_to_save     ;
        num = [];        % num(numeric).
        coe = [];        % coe(symbolic).
        Type = 'empty';
        %PlaneWave = false;
        %pqoL_G =[];
        %PlaneWaveN = 0;
        %PlaneWaveExpandDirection = [0 0 0];
    end
    properties (Dependent=true,Transient = true)
        HsymL_k;
        Hk_sym ;
        Hk_latex ;
    end
    %% Define which properties show
    methods (Access = protected)
        function propgrp = getPropertyGroups(~)
            proplist = {'Degree','Kinds','HsymL','HcoeL','HnumL'};
            propgrp = matlab.mixin.util.PropertyGroup(proplist);
        end
    end
    %% constuction
    methods
        function H_hk = HK(BASIS_NUM,Degree,Term_list,propArgs)
            arguments
                BASIS_NUM = 4;
                Degree = 2;
                Term_list = [];
                propArgs.?vasplib;
            end
            %
            propArgsCell = namedargs2cell(propArgs);      
            H_hk = H_hk@vasplib(propArgsCell{:});
            switch nargin
                case 1
                    if isnumeric(BASIS_NUM)
                        H_hk.Degree = 1;
                        H_hk = H_hk.Degree2Kinds();
                        H_hk.Basis_num = BASIS_NUM;
                        H_hk.HcoeL = sym(zeros(BASIS_NUM,BASIS_NUM,H_hk.Kinds));
                        H_hk.HnumL = (zeros(BASIS_NUM,BASIS_NUM,H_hk.Kinds));
                        %H_hk.Hk_sym = sym(zeros(BASIS_NUM,BASIS_NUM));
                        H_hk.Term_to_save =sym(zeros(BASIS_NUM,BASIS_NUM));
                        H_hk.Trig_to_save =sym(zeros(BASIS_NUM,BASIS_NUM));
                    else
                        if isa(BASIS_NUM,'sym')
                            Hsym = BASIS_NUM;
                            if size(Hsym,1) ~= size(Hsym,2)
                                error('square sym mat required!')
                            end
                            try
                                H_hk.Degree = HK.checkDegree(Hsym);
                            catch
                                %H_hk.Degree = 1;
                            end
                            H_hk =H_hk.Degree2Kinds();
                            BASIS_NUM = length(Hsym);
                            H_hk.Basis_num = BASIS_NUM;
                            H_hk.HcoeL = sym(zeros(BASIS_NUM,BASIS_NUM,H_hk.Kinds));
                            H_hk.HnumL = (zeros(BASIS_NUM,BASIS_NUM,H_hk.Kinds));
                            %H_hk.Hk_sym = sym(zeros(BASIS_NUM,BASIS_NUM));
                            H_hk.Term_to_save =sym(zeros(BASIS_NUM,BASIS_NUM));
                            H_hk.Trig_to_save =sym(zeros(BASIS_NUM,BASIS_NUM));
                            H_hk = H_hk + Hsym;
                        end
                    end
                case 2
                    if isnumeric(BASIS_NUM)
                        H_hk.Degree = Degree;
                        H_hk =H_hk.Degree2Kinds();
                        H_hk.Basis_num = BASIS_NUM;
                        H_hk.HcoeL = sym(zeros(BASIS_NUM,BASIS_NUM,H_hk.Kinds));
                        H_hk.HnumL = (zeros(BASIS_NUM,BASIS_NUM,H_hk.Kinds));
                        %H_hk.Hk_sym = sym(zeros(BASIS_NUM,BASIS_NUM));
                        H_hk.Term_to_save = sym(zeros(BASIS_NUM,BASIS_NUM));
                        H_hk.Trig_to_save = sym(zeros(BASIS_NUM,BASIS_NUM));
                    else

                    end
                case 3
                    if isnumeric(BASIS_NUM)
                        H_hk = HK(BASIS_NUM,Degree);
                        H_hk = H_hk + Term_list;
                    else

                    end

            end
            %   此处显示详细说明

        end
        function H_hk = setup_single(H_hk,symbolic_polynomial,i,j,silence)
            if nargin <5
                silence = false;
            end
            tempmat = zeros( H_hk.Basis_num);
            tempmat(i,j) =1 ;
            H_hk = H_hk.setup_rough(symbolic_polynomial,tempmat,silence);

        end
        function H_hk = setup_rough(H_hk,symbolic_polynomial,pauli_mat,silence)
            if nargin < 4
                silence = false;
            end
            [coeffs_list,symvar_monomial_list] = coeffs(expand(symbolic_polynomial));
            nc = length(coeffs_list);
            for i =1:nc
                [Var_cell{i},k_cell{i},~] = HK.coeff_extract(HK.standardize_sym(symvar_monomial_list(i)));
                mat_cell{i} = pauli_mat;
                Var_cell{i} = Var_cell{i}*coeffs_list(i);
            end
            if nc >0
                H_hk = setup(H_hk,Var_cell,k_cell,mat_cell,silence);
            end
        end
        function H_hk = setup(H_hk,Var_cell,k_cell,mat_cell,silence)
            if nargin < 5
                silence = false;
            end
            if length(Var_cell)~=length(k_cell) && length(k_cell)~=length(mat_cell)
                error('error!');
            end
            if ~silence
                pb = vasplib_tool_outer.CmdLineProgressBar('Setting ');
            end
            nVar_cell = length(Var_cell);
            for i =1:nVar_cell
                Var = Var_cell{i};
                k_symbol = k_cell{i};
                %disp(k_symbol);
                if isa(mat_cell{i},'sym')
                    matcell = sum(sym(mat_cell{i}),3);
                else
                    matcell = sum(double(mat_cell{i}),3);
                end
                Kind = H_hk.k_symbol2Kind(k_symbol);
                if ~silence
                    pb.print(i,nVar_cell,' term into HK ...');
                    %                      fprintf('setting %dth (%s) Hk\n',Kind,string(H_hk.HsymL(Kind)));
                end
                switch class(Var)
                    case 'sym'
                        H_hk.HcoeL(:,:,Kind) = H_hk.HcoeL(:,:,Kind)+ matcell*Var;
                    case 'string'

                    case 'double'
                        H_hk.HnumL(:,:,Kind) = H_hk.HcoeL(:,:,Kind)+ matcell*Var;
                end

            end
            if ~silence
                pb.delete();
            end
        end
    end
    %% de constuction
    methods
        function [H_sym_Gamma,H_latex_Gamma] = GammaDecomposition(H_hk)
            H_sym = H_hk.Hk_sym;
            %
            [H_sym_Gamma,~,H_latex_Gamma] = vasplib.GammaDecomposition(H_sym);
        end
        function [H_sym_pauli,H_latex_pauli] = pauliDecomposition(H_hk)
            H_sym = H_hk.Hk_sym;
            %
            [H_sym_pauli,~,H_latex_pauli]= vasplib.pauliDecomposition(H_sym);
        end
    end
    %% get
    methods
        function HsymL_k = get.HsymL_k(H_hk)
            HsymL_k = H_hk.HsymL;
        end
        function Hk_sym = get.Hk_sym(H_hk)
            Hk_sym = sym(zeros(H_hk.Basis_num,H_hk.Basis_num));
            for i =1:H_hk.Kinds
                Hk_sym = Hk_sym + H_hk.HcoeL(:,:,i)*H_hk.HsymL_k(i);
            end
            try
                Hk_sym = Hk_sym + H_hk.Trig_to_save; %temp
            catch
            end
        end
        function Hk_latex = get.Hk_latex(H_hk)
            Hk_latex = latex(H_hk.Hk_sym);
        end
        function Type = get.Type(H_hk)
            if isequal(H_hk.Trig_to_save,sym(zeros(size(H_hk.Trig_to_save))))
                TB = 0;
            else
                TB = 1;
            end
            if strcmp(class(H_hk.Term_to_save), class(sym(zeros(H_hk.Basis_num))))
                KP = 0 ;
            else
                KP =1;
            end
            if ~TB && ~KP
                Type = 'empty';
            elseif ~TB && KP
                Type = 'kp';
            elseif TB && ~KP
                Type = 'tb';
            elseif TB && KP
                Type = 'kp&tb';
            end
        end
        function [num_label,coe_label,H_hk] = NumOrCoe(H_hk)
            if isempty(H_hk.num) ||  isempty(H_hk.coe)
                [num_label,coe_label] = NumOrCoe@vasplib(H_hk);
                H_hk.num = num_label;
                H_hk.coe = coe_label;
            else
                num_label = H_hk.num;
                coe_label = H_hk.coe;
            end
        end
    end
    %% modify
    methods
        function H_hk_out = PlaneWaveExpand(H_hk,N,ExpandDirection)
            % The Hamiltonian is written in componenets of <mnl|H|pqo>,
            % where |mnl> corresponds to state with wave vector
            % k + m*G_1 + n*G_2 +l*G_3 (G^M_i is Moire reciprocal vectors)
            % m,n,l, p,q,o changes in the range [-N,N]
            arguments
                H_hk HK;
                N double{mustBeInteger} = 5; %integernum for plane wave expanding
                ExpandDirection = [1 1 0];
            end
            % check
            if H_hk.num
                H_hk.HcoeL = sym(H_hk.HnumL);
                H_hk.num = false; H_hk.coe = true;
            else

            end
            %if strcmp(H_hk.Type,'kp')
            %else
            %    error('!!!');
            %end
            %
            ExpandNum =(2*N+1)^sum(ExpandDirection) ;
            BasisNum = ExpandNum * H_hk.Basis_num;
            NL = (-N:N).';
            if ExpandDirection(1)
                mL = NL;
            else
                mL = 0;
            end
            if ExpandDirection(2)
                nL = NL;
            else
                nL = 0;
            end
            if ExpandDirection(3)
                lL = NL;
            else
                lL = 0;
            end
            mnlL = [...
                kron(mL,ones(length(nL)*length(lL),1)),...
                kron(ones(length(mL),1),kron(nL,ones(length(lL),1))),...
                kron(ones(length(mL)*length(nL),1),lL)];
            pqoL = mnlL;
            syms G_0_0_0 k_x k_y k_z real;
            H_hk.HcoeL = H_hk.HcoeL;
            HcoeListpre = H_hk.HcoeL;
            pqoL_G = kron(ones(H_hk.Basis_num ,1 ),pqoL * H_hk.Gk);
            Symvarlist = symvar(HcoeListpre);
            SymvarlistStrL = string(Symvarlist);
            G_Symvarlist = Symvarlist(contains(SymvarlistStrL,"G_"));
            nG_Symvarlist = length(G_Symvarlist);
            G_SymvarCell{nG_Symvarlist} = sym('0','real');
            G_SymvarCellSubs{nG_Symvarlist} = zeros(ExpandNum);
            for i = 1:nG_Symvarlist
                G_SymvarCell{i} = G_Symvarlist(i);
                VectorTmp = park.Variable2Vector(G_Symvarlist(i));
                G_SymvarCellSubs{i} = park.DeltaList(mnlL,VectorTmp+pqoL);
            end
            G_SymvarCell{nG_Symvarlist+1} = G_0_0_0;
            G_SymvarCellSubs{nG_Symvarlist+1} = eye(ExpandNum);
            %
            HcoeListpre = expand(simplify(HcoeListpre*G_0_0_0));
            % 
            HcoeListpreStr = string(HcoeListpre);
            HcoeListpreStr = strrep(HcoeListpreStr,"G_0_0_0*G","G");
            HcoeListpreModify = str2sym(HcoeListpreStr);
            %for i = 1:numel(HcoeListpreStr)
            %    HcoeListpreModify(i) = sym(0);
            %    iHcoeListStr = strrep(HcoeListpreStr(i),'-',' + -');
            %    iHcoeListStrSplit = split(iHcoeListStr,{' + '});
            %    for j = 1:length(iHcoeListStrSplit)
            %        if contains(iHcoeListStrSplit{j},"G_")
            %            iHcoeListStrSplit{j} = strcat("G_0_0_0 * ",iHcoeListStrSplit{j});
            %        end
            %        HcoeListpreModify(i) = HcoeListpreModify(i)+str2sym(iHcoeListStrSplit{j});
            %    end
            %end
            %
            HcoeList = subs(HcoeListpreModify,G_SymvarCell,G_SymvarCellSubs);
            %Hksym = sym(zeros(BasisNum,BasisNum));
            %H_hk_HsymL_k = H_hk.HsymL_k;
            %for i =1:H_hk.Kinds
            %    Hksym = Hksym + HcoeList(:,:,i) * expand(subs(H_hk_HsymL_k(i),...
            %        {k_x,k_y,k_z},{k_x + diag(pqoL_G(:,1)),k_y + + diag(pqoL_G(:,2)),k_z + diag(pqoL_G(:,3)) }));
            %end
            H_hk_out = H_hk;
            H_hk_out.HcoeL = HcoeList;
            H_hk_out.Basis_num = BasisNum;
            %H_hk_out.PlaneWave = true;
            H_hk_symL = H_hk.HsymL;
            H_hk_symL_2 = sym([ones(BasisNum,H_hk.Kinds)]);
            %syms Kx Ky Kz real;
            for i = 2:H_hk.Kinds
                for j = 1:BasisNum
                    H_hk_symL_2(j,i) = subs(H_hk_symL(i),[k_x,k_y,k_z],[...
                        k_x+pqoL_G(j,1),...
                        k_y+pqoL_G(j,2),...
                        k_z+pqoL_G(j,3)...
                        ]);
                end
            end
            H_hk_out.HsymL = H_hk_symL_2;
            %H_hk_out.PlaneWaveN = N;
            %H_hk_out.PlaneWaveExpandDirection = ExpandDirection;
            %H_hk_out.pqoL_G = H_hk.pqoL_G;
        end
        function H_hk = reseq(H_hk,basis_list)
            if isempty(H_hk.coe)||isempty(H_hk.num)
                [~,~,H_hk] = H_hk.NumOrCoe();
            end
            % wan first
            if ~isequal(basis_list,':')
                if H_hk.num
                    H_hk.HnumL=H_hk.HnumL(basis_list,basis_list,:);
                else
                    H_hk.HnumL = [];
                end
                if H_hk.coe
                    H_hk.HcoeL=H_hk.HcoeL(basis_list,basis_list,:);
                else
                    H_hk.HcoeL = sym([]);
                end
            end
            % fix bug
            H_hk.Basis_num = length(basis_list);
            H_hk.Trig_to_save = sym(zeros(H_hk.Basis_num));
            if ~isempty(H_hk.sites)
                try
                    H_hk.sites = H_hk.sites(basis_list);
                catch
                    % bug
                end
            end
            if ~isempty( H_hk.orbL )
                H_hk.orbL = H_hk.orbL(basis_list,:);
            end
            if ~isempty( H_hk.elementL )
                H_hk.elementL = H_hk.elementL(basis_list,:);
            end
            if ~isempty( H_hk.quantumL)
                H_hk.quantumL = H_hk.quantumL(basis_list,:);
            end
            % bug here sym_orbL waiting
        end
    end
    %% methods reload
    methods
        function C = plus(A,B)
            if isa(A,'HK') && isa(B,'HK')
                H_hk1 = A;
                H_hk2 = B;
                % ---------check----------
                if H_hk1.Basis_num ~= H_hk2.Basis_num
                    error('basis num differ');
                end
                if H_hk1.Degree ~= H_hk2.Degree
                    error('Degree differ');
                end
                C = H_hk1;
                C.HcoeL = C.HcoeL + H_hk2.HcoeL;
            elseif isa(A,'HK') && ~isa(B,'HK')
                if isa(B,'Term')
                    % disp('gg');
                    C = A;
                    if contains(C.Type,'kp')
                        C.Term_to_save = C.Term_to_save + B;
                    else
                        C.Term_to_save = B;
                    end
                    for i = 1:length(B)
                        C = C.setup_rough(B(i).symbolic_polynomial,B(i).pauli_mat);
                    end
                elseif isa(B,'Trig')
                    C = A;
                    C.Trig_to_save = C.Trig_to_save+B.symbolic_polynomial*double(B.pauli_mat);
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
                                A = A.setup_rough(B(i,j),tempmat);
                            end
                        end
                    end
                    C =A;
                else
                    C = A;
                end
            elseif ~isa(A,'HK') && isa(B,'HK')
                if isa(A,'Term')
                    C = B;
                    if contains(C.Type,'kp')
                        C.Term_to_save = C.Term_to_save + A;
                    else
                        C.Term_to_save = A;
                    end
                    for i = 1:length(A)
                        C = C.setup_rough(A(i).symbolic_polynomial,A(i).pauli_mat);
                    end
                elseif isa(A,'Trig')
                    C = B;
                     
                    C.Trig_to_save = C.Trig_to_save + A.symbolic_polynomial*double(A.pauli_mat);
                else
                    C = B;
                end
            else
                C = 0;
            end
            %disp(C.Term_to_save);
        end
        function C = minus(A,B)
            if isa(A,'HK') && isa(B,'HK')
                C = A;
            elseif isa(A,'HK') && ~isa(B,'HK')
                if isa(B,'Term')
                    % disp('gg');
                    for i = 1:length(B)
                        C = A.setup_rough(B(i).symbolic_polynomial,B(i).pauli_mat);
                    end
                else
                    C = A;
                end
            elseif ~isa(A,'HK') && isa(B,'HK')
                if isa(A,'Term')
                    if length(B) == 2
                        for i = 1:length(A)
                            C = B.setup_rough(A(i).symbolic_polynomial,A(i).pauli_mat);
                        end
                    else
                        C = B;
                    end
                else
                    C = B;
                end
            else
                C = 0;
            end
        end
        function C = mrdivide(A,B)
            if isa(A,'HK') && isa(B,'HK')
            elseif isa(A,'HK') && ~isa(B,'HK')
                C = A;
                
                %                 H_hr2 = B;
                % ---------check----------
                %                 if H_hr.WAN_NUM ~= length(B)
                %                     error('WAN_NUM different');
                %                 end
                C.HcoeL =  C.HcoeL/B;
            else
            end
        end
        function C = mtimes(A,B)
            if isa(A,'HK') && isa(B,'HK')
                C = A;
            elseif isa(A,'HK') && ~isa(B,'HK')
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
            elseif ~isa(A,'HK') && isa(B,'HK')
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
        function C = mtimes_inv(B,A)
            if ~isa(A,'HK') && isa(B,'HK')
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
        function C = lt(A,B)% overload 	lt(A,B)
            if isa(A,'HK') && isa(B,'HK')
                H_hk1 = A;
                H_hk2 = B;
                % ---------check----------
                if H_hk1.Bassis_num ~= H_hk2.Bassis_num
                    error('Bassis_num different');
                end
                error('not support at present.')
            elseif isa(A,'HK') && ~isa(B,'HK')
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
            elseif ~isa(A,'HK') && ~isa(B,'HK')
                error('not support at present.');
            end
            
        end
        function C = le(A,B)
            if isa(A,'HK') && isa(B,'HK')
                H_hk1 = A;
                H_hk2 = B;
                % ---------check----------
                if H_hk1.Bassis_num ~= H_hk2.Bassis_num
                    error('Bassis_num different');
                end
                error('not support at present.')
            elseif isa(A,'HK') && ~isa(B,'HK')
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
            elseif ~isa(A,'HK') && ~isa(B,'HK')
                error('not support at present.');
            end
        end
        function C = horzcat(A,B)
            if isa(A,'HK') && isa(B,'HK')
                H_hk1 = A;
                H_hk2 = B;
                % ---------init-----------
                %                 H_hk =  HR(H_hk1.WAN_NUM+H_hk2.WAN_NUM,...
                %                     unique([H_hk1.vectorL;H_hk2.vectorL],'rows'));
                H_hk = A;
                H_hk.Basis = [H_hk1.Basis;H_hk2.Basis];
                H_hk.Basis_num = H_hk1.Basis_num+H_hk2.Basis_num;
                H_hk.HnumL = zeros(H_hk.Basis_num,H_hk.Basis_num,H_hk.Kinds);
                H_hk.HcoeL = sym(H_hk.HnumL);
                %             zeros_num_mat = zeros(H_hk.WAN_NUM);
                %             zeros_coe_mat = sym(zeros(H_hk.WAN_NUM));
                for i = 1:H_hk.Kinds
                    
                    H_hk.HnumL(:,:,i) = blkdiag(H_hk1.HnumL(:,:,i) ,...
                        H_hk2.HnumL(:,:,i));
                    H_hk.HcoeL(:,:,i) = blkdiag(H_hk1.HcoeL(:,:,i) ,...
                        H_hk2.HcoeL(:,:,i));
                    
                end
                H_hk.Trig_to_save =sym(zeros(H_hk.Basis_num,H_hk.Basis_num));
                %H_hk.Term_to_save =blkdiag(H_hk1.Term_to_save,H_hk2.Term_to_save);
            else
            end
            C = H_hk;
        end
        function H_hk = conj(H_hk)
            for i =1:H_hk.Kinds
                H_hk.HnumL(:,:,i) = conj(H_hk.HnumL(:,:,i));
                H_hk.HcoeL(:,:,i) = conj(H_hk.HcoeL(:,:,i));
            end
        end
        function H_hk = ctranspose(H_hk)
            for i =1:H_hk.Kinds
                H_hk.HnumL(:,:,i) = H_hk.HnumL(:,:,i)';
                H_hk.HcoeL(:,:,i) = H_hk.HcoeL(:,:,i)';
            end
        end
        function H_hk_sym =sym(H_hk)
            H_hk_sym = H_hk.Hk_sym;
        end
        %         function H_hk = disp(H_hk)% disp
        %             % Hk_sym = ;
        %
        %             %
        %             % disp(H_hk.Hk_latex);
        %             builtin('disp',H_hk);
        %             disp(H_hk.Hk_latex);
        %         end
        function [H_hk,Sublist,Unique_term] =  unique(H_hk,seed,checklist,options)
            arguments
                H_hk HK;
                seed char = 'eta';
                checklist =sym([-1,sqrt(3),-sqrt(3),2,-2,3,-3,6,-6]);
                options.Accuracy = 1e-6;
                options.simplify logical=true;
            end
            HcoeL_tmp_1 = H_hk.HcoeL(:);
            HcoeL_tmp = [real(HcoeL_tmp_1);imag(HcoeL_tmp_1)];
            % Magnification divide give up
            if  options.simplify
                HcoeL_tmp = vasplib.cleanVar(HcoeL_tmp,options.Accuracy);
            end
            % HcoeL_tmp(find(abs(CoeffsList)<sym(Accuracy))) = sym(0);
            % unique
            [Unique_term,~,ic] = unique(HcoeL_tmp);
            %[CoeffsList,~] = HK.factorAll(Unique_term);
            %
            SymVar = sym(seed,[length(Unique_term),1],'real');
            % zero test
            SymVar(find(Unique_term==sym(0))) = sym(0);
            % minus test sqrt3 test 2 test
            for coeff_for_check = checklist
                [Lia,Locb] = ismember(expand(Unique_term),expand(coeff_for_check*Unique_term));
                Equationlist = HR.isolateAll((SymVar(Lia) == coeff_for_check*SymVar(Locb(Locb>0)) ));
                SymVar = subs(SymVar,lhs(Equationlist),rhs(Equationlist));
            end
            % 
            % waiting
            Sublist = SymVar == Unique_term;
            %Sublist_i = SymVar_i == Unique_term_i;
            HcoeL_tmp_1 = SymVar(ic(1:end/2)) ...
                +1i*SymVar(ic(end/2+1:end));
            %             HcoeL_tmp_1 = SymVar(ic(1:end/2)) .* CoeffsList(1:end/2) ...
            %                 +1i*SymVar(ic(end/2+1:end)) .* CoeffsList(end/2+1:end) ;
            H_hk.HcoeL = reshape(HcoeL_tmp_1,H_hk.Basis_num,H_hk.Basis_num,H_hk.Kinds);
        end
        function H_hk = sum(H_hk_list)
            H_hk = H_hk_list(1);
            for i = 2:length(H_hk_list)
                H_hk = H_hk + H_hk_list(i);
            end
        end
        function H_hk = simplify(H_hk,Accuracy)
            H_hk.HcoeL = simplify(H_hk.HcoeL);
        end
    end
    %% conventional
    methods       
%         function H_hk = setup_rough_trig(H_hk,symbolic_polynomial_trig,pauli_mat)
%             [coeffs_list,symvar_monomial_list] = coeffs(expand(symbolic_polynomial));
%             for i =1:length(coeffs_list)
%                 [Var_cell{i},k_cell{i},~] = HK.coeff_extract(HK.standardize_sym(symvar_monomial_list(i)));
%                 mat_cell{i} = pauli_mat;
%                 Var_cell{i} = Var_cell{i}*coeffs_list(i);
%             end
%         end 
        function H_hk = Degree2Kinds(H_hk)
            %syms k_x k_y k_z real;
            VarUsing = H_hk.VarsSeqLcart(1:H_hk.Dim);
            str_tmp = string(expand((1+fold(@plus,VarUsing))^H_hk.Degree));      
            str_cell = strsplit(str_tmp,'+');
            H_hk.Kinds = length(str_cell);
            coeff_list = zeros(H_hk.Kinds,1);
            Degree_list = zeros(H_hk.Kinds,1);
            for i = 1:H_hk.Kinds
                tmpsym = str2sym(str_cell{i});
                Degree_list(i) = polynomialDegree(tmpsym);
                coeff_list(i) = coeffs(tmpsym);
            end
            [~,sort_label] = sort(Degree_list);
            coeff_list = coeff_list(sort_label);
            %[coeff_list,sort_label] = sort(sort_init_list,'descend');
            %sort_label = [sort_label(H_hk.Kinds);sort_label(1:H_hk.Kinds-1)];
            %coeff_list = [coeff_list(H_hk.Kinds);coeff_list(1:H_hk.Kinds-1)];
            str_cell = str_cell(sort_label);
            for i = 1:H_hk.Kinds
                H_hk.HsymL(i) = str2sym(str_cell{i})/coeff_list(i);
            end
            H_hk.HstrL = string(H_hk.HsymL);
            temp_strL = H_hk.HstrL;
             temp_strL = strrep(temp_strL,'k_x','x');
             temp_strL = strrep(temp_strL,'k_y','y');
             temp_strL = strrep(temp_strL,'k_z','z');
            for i = 1:H_hk.Kinds
                H_hk.HstrL(i) = HK.quadk2K(temp_strL(i));
            end
            H_hk.HsymL_xyz = str2sym(temp_strL);
        end
    end
    %% symmetry
    methods
        function H_hk  = init(H_hk)
            sizeH = size(H_hk.HcoeL);
            %For compatible with Yang Result
            %HcoeL_tmp = sym('A',sizeH,'real') - 1i*sym('B',sizeH,'real');
            %
            HcoeL_tmp = sym('A',sizeH,'real')+1i*sym('B',sizeH,'real');
            % HcoeL_tmp = sym('C',sizeH,'complex');
            for i =1:H_hk.Kinds
                HcoeL_tmp(:,:,i) = triu(HcoeL_tmp(:,:,i))+triu(HcoeL_tmp(:,:,i),1)';
            end
            H_hk.HcoeL = HcoeL_tmp;
            %     HcoeL_tmp
            %     H_hk.
        end
        function H_hk = applyU(H_hk,U,conjugate ,antisymmetry )
            arguments
                H_hk HK;
                U = nan;
                conjugate logical =false;
                antisymmetry logical = false;
            end
            if isa(U,'Oper')
                conjugate = U.conjugate;
                antisymmetry = U.antisymmetry;
                U = U.U;
            end
            if isequal(zeros(size(H_hk.HnumL)),H_hk.HnumL)
                num_label = false;
            else
                num_label = true;
            end
            if isequal(sym(zeros(size(H_hk.HcoeL))),H_hk.HcoeL)
                coe_label = false;
            else
                coe_label = true;
            end
            
            U_inv = inv(U);
            if coe_label == true
                HcoeLtmp = H_hk.HcoeL;
                if conjugate
                    HcoeLtmp = conj(HcoeLtmp);
                    HcoeLtmp = HK.matrixtimespage(HK.factorlist_parity(H_hk.Degree),HcoeLtmp);
                end
                if antisymmetry
                    HcoeLtmp = -HcoeLtmp;
                end
                HcoeLtmp = HK.page_mtimes_matrix(HK.matrix_mtimes_page(U,HcoeLtmp),U_inv);
                H_hk.HcoeL = HcoeLtmp;
            end
            if num_label == true
                U_page = repmat(U,[1 1 size(H_hk.HnumL,3)]);
                U_inv_page = repmat(U_inv,[1 1 size(H_hk.HnumL,3)]);
                HnumLtmp = H_hk.HnumL;
                if conjugate
                    HnumLtmp = conj(HnumLtmp);
                    HnumLtmp = HK.matrixtimepage(HK.factorlist_parity(H_hk.Degree),HnumLtmp);
                end
                if antisymmetry
                    HnumLtmp = -HnumLtmp;
                end
                HnumLtmp = pagemtimes(pagemtimes(U_page,HnumLtmp),U_inv_page);
                H_hk.HnumL = HnumLtmp;
            end
        end
        % ---------------------   read function   ------------------------
        function matcell = matgen(H_hk,R,Accuracy)
            arguments
                H_hk HK;
                R  ;
                Accuracy double = 6;
            end
            if isa(R,'sym')
                coe_label = true;
            else
                coe_label = false;
            end
            VarUsing = H_hk.VarsSeqLcart(1:H_hk.Dim);
            %
            for i = 1:H_hk.Degree+1
                Orderlist = HK.orderlist(i-1);
                nOrderlist = length(Orderlist);
                if coe_label
                    matcell{i} = sym(zeros(nOrderlist));
                else
                    matcell{i} = zeros(nOrderlist);
                end
                HsymC{i}(1:nOrderlist) = H_hk.HsymL_k(Orderlist);
            end  
            HsymC_bk = HsymC;
            k_orgin = VarUsing.';
            k_R = R*k_orgin;
            for i = 1:H_hk.Degree+1
                HsymC{i} = subs(HsymC{i},k_orgin,k_R);
            end
            for i = 1:H_hk.Degree+1
                H_symL_i = simplify(HsymC{i}) ;
                nOrderlist = length(H_symL_i);
                for j = 1:nOrderlist
                    [A,B] = coeffs(H_symL_i(j),VarUsing);
                    for k = 1:length(B)
                        tempSym = B(k);
                        for l  = 1:nOrderlist
                            if isequal(tempSym,HsymC_bk{i}(l))
                                %matcell{i}(l,j)=A(k);% for the need of symmetry constrain
                                matcell{i}(j,l)=A(k);
                                break;
                            end
                        end
                    end
                end
                if ~coe_label
                    matcell{i} = roundn(matcell{i},-Accuracy);
                end
%                 if i == 3
%                     error('ss');
%                 end
            end
            
        end
        function H_hk = applyR(H_hk,R)
            arguments
                H_hk HK;
                R ;
            end
            if isequal(zeros(size(H_hk.HnumL)),H_hk.HnumL)
                num_label = false;
            else
                num_label = true;
            end
            if isequal(sym(zeros(size(H_hk.HcoeL))),H_hk.HcoeL)
                coe_label = false;
            else
                coe_label = true;
            end
            matcell = H_hk.matgen(R);
            if num_label
                for i = 0:H_hk.Degree
                    Orderlist = HK.orderlist(i);
                    H_hk.HnumL(:,:,Orderlist) = HK.matrixtimespage(matcell{i+1},H_hk.HnumL(:,:,Orderlist));
                end
            end
            if coe_label
                for i = 0:H_hk.Degree
                    Orderlist = HK.orderlist(i);
                    H_hk.HcoeL(:,:,Orderlist) = HK.matrixtimespage(matcell{i+1},H_hk.HcoeL(:,:,Orderlist));
                end
            end
        end
        function [H_hk_R,H_hk] = applyRU(H_hk,SymOper)
            if ~SymOper.conjugate && ~SymOper.antisymmetry && isequal(SymOper.R,eye(3))
                H_hk_R = H_hk;
                return;
            end
            if isequal(zeros(size(H_hk.HnumL)),H_hk.HnumL)
                num_label = false;
            else
                num_label = true;
            end
            if isequal(sym(zeros(size(H_hk.HcoeL))),H_hk.HcoeL)
                coe_label = false;
            else
                coe_label = true;
            end
            % left ? right
            if SymOper.conjugate
                R_inv = -inv(SymOper.R);
            else
                R_inv = inv(SymOper.R);
            end
            %
%             if SymOper.conjugate
%                 R_inv = -(SymOper.R);
%             else
%                 R_inv = (SymOper.R);
%             end
            H_hk_R = H_hk;
            matcell = H_hk.matgen(R_inv);
            if num_label
                for i = 0:H_hk.Degree
                    Orderlist = HK.orderlist(i);
                    % !
                    H_hk_R.HnumL(:,:,Orderlist) = HK.matrixtimespage(matcell{i+1}.',H_hk_R.HnumL(:,:,Orderlist));
                end
            end
            if coe_label
                for i = 0:H_hk.Degree
                    Orderlist = HK.orderlist(i);
                    H_hk_R.HcoeL(:,:,Orderlist) = HK.matrixtimespage(matcell{i+1}.',H_hk_R.HcoeL(:,:,Orderlist));
                end
            end
            
            if isnan(SymOper.U)
                %build U
            else
                U = SymOper.U;
            end
            U_inv = inv(U);
            if coe_label == true
                HcoeLtmp = H_hk_R.HcoeL;
                if SymOper.conjugate
                    HcoeLtmp = conj(HcoeLtmp);
                    %HcoeLtmp = HK.matrixtimespage(HK.factorlist_parity(H_hk_R.Degree),HcoeLtmp);
                end
                if SymOper.antisymmetry
                    HcoeLtmp = -HcoeLtmp;
                end
                HcoeLtmp = HK.page_mtimes_matrix(HK.matrix_mtimes_page(U,HcoeLtmp),U_inv);
                H_hk_R.HcoeL = (HcoeLtmp);
            end
            if num_label == true
                U_page = repmat(U,[1 1 size(H_hk_R.HnumL,3)]);
                U_inv_page = repmat(U_inv,[1 1 size(H_hk_R.HnumL,3)]);
                HnumLtmp = H_hk_R.HnumL;
                if conjugate
                    HnumLtmp = conj(HnumLtmp);
                    HnumLtmp = HK.matrixtimepage(HK.factorlist_parity(H_hk_R.Degree),HnumLtmp);
                end
                if antisymmetry
                    HnumLtmp = -HnumLtmp;
                end
                HnumLtmp = pagemtimes(pagemtimes(U_page,HnumLtmp),U_inv_page);
                H_hk_R.HnumL = HnumLtmp;
            end
            
        end
        function [H_hk] = applyOper(H_hk,SymOper,options)
            arguments
                H_hk HK;
                SymOper Oper = Oper();
                options.generator  logical = true;
                options.Accuracy = 1e-3;
                options.oneshot = false;
                options.center = [0,0,0];
            end
            %             if isequal(zeros(size(H_hk.HnumL)),H_hk.HnumL)
            %                 num_label = false;
            %             else
            %                 num_label = true;
            %             end
            if isequal(sym(zeros(size(H_hk.HcoeL))),H_hk.HcoeL)
                coe_label = false;
            else
                coe_label = true;
            end
            
            if ~coe_label
                H_hk = H_hk.init();
                H_hk = H_hk.hermitize();
            end
            % refresh SymOper
            try
            BasisFunction = BasisFunc(H_hk);
            SymOper = SymOper.Ugen(BasisFunction,'Rm',H_hk.Rm,'center',options.center);
            catch
                error('Oper U is nan');
            end
            if length(SymOper) == 1
                nSymOper = length(SymOper);
                i = 1;
                fprintf('******** apply (%d/%d)symmetry ********\n',i,nSymOper);
                disp(SymOper(i));
                if options.generator
                    SymOper_tmp = SymOper(i).generate_group();
                    nSymOper_tmp = length(SymOper_tmp);
                    pb = vasplib_tool_outer.CmdLineProgressBar('Applying Symmetry ...');
                    H_hk_R = H_hk;
                    %H_hr_R = repmat(H_hr_R,[]);
                    for j = 1:nSymOper_tmp
                        pb.print(j,nSymOper_tmp);
                        [H_hk_R(j),H_hk] = applyRU(H_hk,SymOper_tmp(j));
                    end
                    pb.delete();
                    H_hk = sum(H_hk_R)/nSymOper_tmp;
                    H_hk = H_hk.simplify(options.Accuracy);
                else
                    [H_hk,H_hk_bk] = applyRU(H_hk,SymOper);
                    Equationlist_r = real(H_hk.HcoeL - H_hk_bk.HcoeL) == 0;
                    Equationlist_i = imag(H_hk.HcoeL - H_hk_bk.HcoeL) == 0;
                    %Equationlist_r = HK.isolateAll(Equationlist_r,real(H_hk.HcoeL));
                    %Equationlist_i = HK.isolateAll(Equationlist_i,imag(H_hk.HcoeL));
                    Equationlist_r = HK.isolateAll(Equationlist_r);
                    Equationlist_i = HK.isolateAll(Equationlist_i);
                    HcoeLtmp = H_hk.HcoeL ;
                    HcoeLtmp_r = subs(real(HcoeLtmp),lhs(Equationlist_r),rhs(Equationlist_r));
                    HcoeLtmp_i = subs(imag(HcoeLtmp),lhs(Equationlist_i),rhs(Equationlist_i));
                    H_hk.HcoeL = HcoeLtmp_r + 1i*HcoeLtmp_i;
                end
                %H_hk = H_hk.applyOper(SymOper(i),'generator','false');
                fprintf('----------   SymVarNum: %d   ----------\n',length(symvar(H_hk.HcoeL)));
            else
                if options.oneshot
                    SymOper_tmp = SymOper.generate_group();
                    nSymOper_tmp = length(SymOper_tmp);
                    pb = vasplib_tool_outer.CmdLineProgressBar('Applying Symmetry ...');
                    H_hk_R = H_hk;
                    %H_hr_R = repmat(H_hr_R,[]);
                    for j = 1:nSymOper_tmp
                        pb.print(j,nSymOper_tmp);
                        [H_hk_R(j),H_hk] = applyRU(H_hk,SymOper_tmp(j));
                    end
                    pb.delete();
                    H_hk = sum(H_hk_R)/nSymOper_tmp;
                    H_hk = H_hk.simplify(options.Accuracy);
                    return;
                end
                nSymOper = length(SymOper);
                for i = 1:nSymOper
                    fprintf('******** apply (%d/%d)symmetry ********\n',i,nSymOper);
                    disp(SymOper(i));
                    if options.generator
                        SymOper_tmp = SymOper(i).generate_group();
                        nSymOper_tmp = length(SymOper_tmp);
                        pb = vasplib_tool_outer.CmdLineProgressBar('Applying Symmetry ...');
                        H_hk_R = H_hk;
                        %H_hr_R = repmat(H_hr_R,[]);
                        for j = 1:nSymOper_tmp
                            pb.print(j,nSymOper_tmp);
                            [H_hk_R(j),H_hk] = applyRU(H_hk,SymOper_tmp(j));
                        end
                        pb.delete();
                        H_hk = sum(H_hk_R)/nSymOper_tmp;
                        H_hk = H_hk.simplify(options.Accuracy);
                    else
                        %fprintf('    ');
                        H_hk = H_hk.applyOper(SymOper(i),'generator','false');
                    end
                    fprintf('----------   SymVarNum: %d   ----------\n',length(symvar(H_hk.HcoeL)));
                end
            end
        end
        function H_hk = hermitize(H_hk)
            if isequal(zeros(size(H_hk.HnumL)),H_hk.HnumL)
                num_label = false;
            else
                num_label = true;
            end
            if isequal(sym(zeros(size(H_hk.HcoeL))),H_hk.HcoeL)
                coe_label = false;
            else
                coe_label = true;
            end
            H_hk_bk = H_hk';
            if coe_label
                H_hk = (H_hk + H_hk_bk)/2;
                %                 Equationlist_r = real(H_hk.HcoeL - H_hk_bk.HcoeL) == 0;
                %                 Equationlist_i = imag(H_hk.HcoeL - H_hk_bk.HcoeL) == 0;
                %                 %Equationlist_r = HK.isolateAll(Equationlist_r,real(H_hk.HcoeL));
                %                 %Equationlist_i = HK.isolateAll(Equationlist_i,imag(H_hk.HcoeL));
                %                 Equationlist_r = HK.isolateAll(Equationlist_r);
                %                 Equationlist_i = HK.isolateAll(Equationlist_i);
                %                 HcoeLtmp = subs(H_hk.HcoeL,lhs(Equationlist_r),rhs(Equationlist_r));
                %                 HcoeLtmp = subs(HcoeLtmp,lhs(Equationlist_i),rhs(Equationlist_i));

            end
            if num_label
                H_hk.HnumL = (H_hk_bk.HnumL + H_hk.HnumL )/2;
            end
        end
        function H_hk_bk = subsOper(H_hk,SymOper)
            [~,H_hk_bk] = H_hk.applyOper(SymOper);
%             Equationlist_r = real(H_hk.HcoeL - H_hk_bk.HcoeL) == 0;
%             Equationlist_i = imag(H_hk.HcoeL - H_hk_bk.HcoeL) == 0;
%             Equationlist_r = HK.isolateAll(Equationlist_r,real(H_hk.HcoeL));
%             Equationlist_i = HK.isolateAll(Equationlist_i,imag(H_hk.HcoeL));
%             HcoeLtmp = HK.subs(H_hk.HcoeL,lhs(Equationlist_r),rhs(Equationlist_r));
%             HcoeLtmp = HK.subs(HcoeLtmp,lhs(Equationlist_i),rhs(Equationlist_i));
%             H_hk.HcoeL = HcoeLtmp;
        end
        function J_mat = Jmat_gen(H_hk,QuantumL,options)
            arguments
                H_hk  HK;
                QuantumL double = H_hk.quantumL;
                options.include_0 logical = false;
                options.strict logical = false;
            end
            J_mat = Oper.spin_matrices_from_orb(QuantumL,options.include_0,options.strict);
        end
        %%
        function H_hk = Subsall(H_hk,mode)
            if nargin <2
                mode = 'num';
            end
            % ------ init ------- load para form inner or outer
            if exist('para.mat','file') && strcmp(mode,'file')
                %disp('The para in para.mat is first to be considered?');
                load('para.mat');
            else
                
            end
            switch mode
                case {'num','file','gen'}
                    HcoeL_temp = subs(H_hk.HcoeL);
                    symname = symvar(HcoeL_temp);
                    if ~isempty(symname)
                        for i = 1:length(symname)
                            warning('this varible: %s is not defind',string(symname(i)));
                        end
                        disp(symname);
                        error('please introduce the varible value above');
                    else
                        H_hk.HnumL = double(HcoeL_temp);
                    end
%                     H_hk.Hk_sym = sym(zeros(H_hk.Basis_num,H_hk.Basis_num));
%                     for i =1:H_hk.Kinds
%                         H_hk.Hk_sym = H_hk.Hk_sym + H_hk.HcoeL(:,:,i)*H_hk.HsymL_k(i);
%                     end
                    H_hk.Hk_num  = subs(H_hk.Hk_sym);
                    H_hk.Hsym  = H_hk.Hk_num ;
                case 'sym'
                    H_hk.HcoeL = subs(H_hk.HcoeL);
                    H_hk.Trig_to_save = subs(H_hk.Trig_to_save);
            end
        end
        function H_hk  = subs(H_hk,varargin)
            switch length(varargin)
                case 1
                    H_hk.HcoeL = subs(H_hk.HcoeL,varargin{1});
                    H_hk.Trig_to_save = subs(H_hk.Trig_to_save,varargin{1});
                case 2
                    H_hk.HcoeL = subs(H_hk.HcoeL,varargin{1},varargin{2});
                    H_hk.Trig_to_save = subs(H_hk.Trig_to_save,varargin{1},varargin{2});
                case 3
                    H_hk.HcoeL = subs(H_hk.HcoeL,varargin{1},varargin{2},varargin{3});
                    H_hk.Trig_to_save = subs(H_hk.Trig_to_save,varargin{1},varargin{2},varargin{3});
            end
        end
       
    end
    %% tools function   
    methods
        function varargout= kp2TB(H_hk,kpoints_f,groups,options)
            %
            arguments
                H_hk HK;
                kpoints_f double= [0,0,0];
                groups = [];
                options.Accuracy = 1e-4;
                options.search_range = [1 1 1];
                options.R_max = -1;
                options.level_cut = 1;
                options.mini = false;
            end
            if isempty(groups)
                if isequal(kpoints_f,[0,0,0])
                    mode = 'Gamma';
                else
                    mode = 'kpoints';
                end
            else
                mode = 'TBsubs';
            end
            %
            R_struct = H_hk.Rm2abc(H_hk.Rm);
            if ~strcmp(mode,'TBsubs')
                if abs(abs(cosd(R_struct.gamma))) < 1e-4
                    lattice_mode = 'T';
                    %disp('wow')
                elseif abs((cosd(R_struct.gamma)) - -1/2) < 1e-4
                    lattice_mode = 'H';
                elseif abs((cosd(R_struct.gamma)) - 1/2) < 1e-4
                    lattice_mode = 'H2';
                else
                    error('unsupport yet! only tetragonal of hexagonal\n');
                end
                if strcmp(mode , 'Gamma')
                    
                else
                    warning('This function is not well-written, only support Term mode now');
                    kpoints = kpoints_f * H_hk.Gk;
                    H_hk = HK(H_hk.Basis_num,H_hk.Degree,H_hk.Term_to_save.subsk(kpoints));
                end
            end
            %%%%%%%%%%%%%%%%%%%%
             H_hr = HR(H_hk.Basis_num,'sym',true);
             H_hr = H_hr.vasplibCopy(H_hk);
             H_hr.coe = true;
             H_hr.num = false;
%             H_hr.Rm = H_hk.Rm;
%             H_hr.orbL = H_hk.orbL; % the orbital list fractional
%             H_hr.elementL = H_hk.elementL; % the element of each orbital
%             H_hr.quantumL = H_hk.quantumL; % [n,l,m,s] the three quantum number combined with spin
%             H_hr.orb_symL = H_hk.orb_symL; % orbital function, makes your own
%             % for spacegroup
%             H_hr.sgn       = H_hk.sgn    ;
%             H_hr.Atom_name = H_hk.Atom_name   ;
%             H_hr.Atom_num  = H_hk.Atom_num   ;
%             H_hr.sites = H_hk.sites           ;
%             H_hr.symmetry_operation = H_hk.symmetry_operation;
%             H_hr.klist_cart = H_hk.klist_cart          ;
%             H_hr.klist_l =  H_hk.klist_l       ;
%             H_hr.klist_frac = H_hk.klist_frac          ;
%             H_hr.kpoints_l = H_hk.kpoints_l         ;
%             H_hr.kpoints_name = H_hk.kpoints_name     ;
%             H_hr.Rnn = H_hk.Rnn                  ;
%             H_hr.nn_store = H_hk.nn_store        ;
%             H_hr = H_hr.timtj_gen('sym')      ;
            % bug fix
            
            if strcmp(mode,'TBsubs')
                if options.R_max == -1
                    R_a = norm(H_hr.Rm(1,:));
                    R_b = norm(H_hr.Rm(2,:));
                    R_c = norm(H_hr.Rm(3,:));
                    R_max = max(options.search_range)*max([R_a,R_b,R_c])+options.Accuracy; % need fix
                else
                    R_max = options.R_max;
                end
                Accuracy = options.Accuracy;
                NAccuracy = round(-log(Accuracy)/log(10));
                H_hr= H_hr.nn(options.search_range,Accuracy,R_max);
                %
                H_hr = H_hr.init('level_cut',options.level_cut,"onsite",1); %最近邻
                %
                %H_hk_bk = H_hk;
                basis_num = H_hk.Basis_num;
                
                if options.mini
                    nonzeromat = zeros(H_hk.Basis_num);
                    for i = 1:length(groups)
                        H_hk_bk = H_hk.applyU(groups(i)) ;
                        nonzeromat =nonzeromat+sum(logical((H_hk_bk.HcoeL~=sym(0))),3);
                    end
                    if strcmp(H_hr.type,'list')
                        H_hr = H_hr.rewind();
                    end
                    zeromat = (nonzeromat == 0);
                    for i = 1:basis_num
                        for j = 1:basis_num
                            if zeromat(i,j)
                                H_hr.HcoeL(i,j,:) = sym(0);
                            end
                        end
                    end
                    H_hr = H_hr.rewrite();
                end

                %
                H_hr = H_hr.applyOper(groups,'generator',true);
                [H_hr,Sublist,Unique_term] = H_hr.unique();% defeat with numerical error
                fprintf('******** Simplify Hamiontonian ********\n');
                fprintf('----------   SymVarNum: %d   ----------\n',length(H_hr.symvar_list));
                %
                H_htrig = H_hr.HR2Htrig();
                H_hr_hk = H_htrig.Htrig2HK(kpoints_f,'Order',H_hk.Degree);
                % find non zero element in H_hk
                label_list = find(H_hr_hk.HcoeL~=sym(0));
                HcoeL_tmp_1 = (H_hr_hk.HcoeL(label_list));
                HcoeL_tmp_2 = (H_hk.HcoeL(label_list));
                Symvar_for_sub = symvar(HcoeL_tmp_1);
                Equationlist_r = real(HcoeL_tmp_1)==real(HcoeL_tmp_2);
                Equationlist_i = imag(HcoeL_tmp_1)==imag(HcoeL_tmp_2);
                Equationlist = vpa([Equationlist_r;Equationlist_i],NAccuracy);
                Equationlist = HK.isolateAll(Equationlist,Symvar_for_sub);
                %
                Equationlist = unique(simplify(Equationlist));
                %Equationlist_force = Equationlist(find(rhs(Equationlist)~=sym(0)));
                %Equationlist_i = unique(Equationlist_i);
                % solve liner Equation
                Solved = solve(Equationlist,Symvar_for_sub);
                Solved_value = sym(ones(1,length(Symvar_for_sub)));
                count = 0;
                for symvartmp = Symvar_for_sub
                    count = count+1;
                    if isempty(Solved.(char(symvartmp)))
                        Solved_value(count) = sym(0);
                    else
                        Solved_value(count) = Solved.(char(symvartmp));
                    end
                end 
                HcoeLtmp = subs(H_hr.HcoeL,Symvar_for_sub,Solved_value);
                Allzero =length(symvar(HcoeLtmp));
                if Allzero == 0
                    warning('No subs result, check by yourselt')
                    Problem_hr.HR = H_hr;
                    Problem_hr.Equationlist = Equationlist;
                    Problem_hr.Sublist = Sublist;
                    Problem_hr.Unique_term = Unique_term;
                    %H_hr.HcoeL = subs(H_hr.HcoeL,lhs(Equationlist_force),rhs(Equationlist_force));
                    %Problem_hr.HR_force = H_hr;
                    varargout{1} = Problem_hr;
                    if nargout == 2
                        varargout{2} = H_htrig;
                    end
                    return;
                else
                    H_hr.HcoeL = HcoeLtmp;
                end
                H_hr = H_hr.rewrite();
                H_hr.HcoeL = vasplib.cleanVar(H_hr.HcoeL,Accuracy);
                varargout{1} = H_hr;
                if nargout == 2
                    H_htrig.HcoeL  = subs(H_htrig.HcoeL,Symvar_for_sub,Solved_value);
                    varargout{2} = H_htrig;
                end
            else
                %%%%%%%%%%%%%%%%%%%%
                for i = 1:H_hk.Kinds
                    %                         if sum(sum(H_hk.HcoeL(:,:,i))) ~= sym(0)
                    if  sum(all(H_hk.HcoeL(:,:,i) == sym(0) ))
                        [vector_list,Coeffs_list] = HK.HstrL_classify(H_hk.HstrL(i),R_struct,lattice_mode);
                        %disp(vector_list);
                        for j =1:length(Coeffs_list)
                            H_hr = H_hr.set_hop_mat(H_hk.HcoeL(:,:,i)*Coeffs_list(j),vector_list(j,:),'symadd');
                        end
                    end
                end
                varargout{1} = H_hr;
            end
        end
        
        function varargout = EIGENCAR_gen(H_hk,options)
            arguments
                H_hk HK;
                options.fermi double = 0;
                options.norb double = -1;
                options.klist double = H_hk.klist_cart;
                options.para  = [];
                options.paraname ;
                options.show = false;
                options.OriginGk = [0 0 0];
                options.ax = handle([]);
                options.printmode = false;
            end
            % -------------- nargin ------------------
            fermi = options.fermi;
            norb_enforce  = options.norb;
            if isempty(options.klist)
                H_hk = H_hk.kpathgen3D('KPOINTS');
                klist_cart_tmp = H_hk.klist_cart;
            else
                klist_cart_tmp = options.klist;
            end
            klist_cart_tmp = klist_cart_tmp + options.OriginGk;
            if options.show
                if iisempty(options.ax)
                   ax = vasplib.BZplot(H_hk.Rm,'color','r');
                else
                    ax = options.ax;
                end
            end
            % -------------- nargin ------------------
            %disp("EIGENCAR gen for H_xyz(wt TB) type: HR class ");
            % Hnum_list = H_hk.HnumL ;
            %             Bassis_mat = H_hk.Bassis_mat ;
            
            if H_hk.Basis_num > 500
                print_mode = 1;
            else
                print_mode = 0;
            end
            [kn,~] = size(klist_cart_tmp);
            %--------  check  --------
            if norb_enforce <0
                NBANDS=H_hk.Basis_num;
            elseif norb_enforce >0
                NBANDS=norb_enforce;
            else
                
            end
            WAVECAR  = zeros(H_hk.Basis_num,NBANDS,kn);
            EIGENCAR = zeros(NBANDS,kn);
            %syms k_x k_y k_z real;
            HsymL_fun = matlabFunction( H_hk.HsymL,'Vars',H_hk.VarsSeqLcart(1:H_hk.Dim));
            %numL = reshape(H_hk.HnumL, NBANDS^2, []);
            for ki =1:kn
                Input = num2cell(klist_cart_tmp(ki,:));
                kL = HsymL_fun(Input{:});
                Hout = H_hk.HnumL;
                %Hout = reshape(numL*kL, NBANDS, NBANDS);
                for i =1:H_hk.Kinds
                    Hout(:,:,i) = Hout(:,:,i).*kL(:,i);
                end
                Hout = sum(Hout,3);
                %                 end
                Hout = (Hout+Hout')/2;
                if norb_enforce <0
                    try
                        [A, U]=eig(full(Hout));
                    catch
                        disp([ki,x,y,z] );
                        disp(Hout);
                        disp(H_hk.Hfun);
                        error('check this k point');
                    end
                elseif norb_enforce >0
                    [A, U]=eigs(Hout,NBANDS,fermi);
                    [A, U]= park.sorteig(U,A);
                else
                end
                EIGENCAR(:,ki) = diag(U);
                WAVECAR(:,:,ki) = A;
                if print_mode ==1
                    fprintf('%d th kpoints has been caculated in %d kpoints total\n',ki,kn);
                end
            end
            varargout{1} = EIGENCAR;
            varargout{2} = WAVECAR;
            if options.show
                [varargout{3}] = vasplib.klist_show(...
                    'klist',klist_cart_tmp,...
                    'ax',ax);
            end
            %             if kn >1
%                 EIGENCARout = EIGENCAR;
%                 %EIGENCARout{2} = WAVECAR;
%             else
%                 EIGENCARout.EIGENCAR= EIGENCAR;
%                 EIGENCARout.WAVECAR = WAVECAR;
%             end
        end
        function Kind = k_symbol2Kind(H_hk,k_symbol)
            k_symbol = string(k_symbol);
            symvar_list = symvar(str2sym(k_symbol));
            %             syminit = sym(1);
            for i = 1:length(symvar_list)
                str_tmp = string(symvar_list(i));
                switch str_tmp
                    case {'k_x','k_X','K_X','K_x','kx','kX','KX','Kx','x','X'}
                        k_symbol = strrep(k_symbol,str_tmp,'k_x');
                    case {'k_y','k_Y','K_y','K_Y','ky','kY','Ky','KY','y','Y'}
                        k_symbol = strrep(k_symbol,str_tmp,'k_y');
                    case {'k_z','k_Z','K_z','K_Z','kz','kZ','Kz','KZ','z','Z'}
                        k_symbol = strrep(k_symbol,str_tmp,'k_z');
                    case {'k_w','k_W','K_w','K_W','kw','kW','Kw','KW'}
                        k_symbol_str = strrep(k_symbol_str,str_tmp,'k_w');
                end
            end
            coeff_tmp = coeffs(str2sym(k_symbol));
            % or use sym to find
            str_2_compare = string((str2sym(k_symbol)/coeff_tmp));
            Kind = find(string(H_hk.HsymL) == str_2_compare);
        end
    end
    methods(Static)
        function Degree = checkDegree(Hsym,Dim)
            arguments
                Hsym
                Dim = 3;
            end
            VarsSeqLcart = [sym('k_x'),sym('k_y'),sym('k_z'),sym('k_w')];
            Degree_mat = polynomialDegree(Hsym,VarsSeqLcart(1:Dim));
            Degree = max(Degree_mat,[],'all');
        end
        function Orderlist  = orderlist(order)
            if order == 0
                Orderlist = 1;
            else
                Orderlist = nchoosek(4+order-2,order-1)+1:nchoosek(4+order-1,order);
            end
            
        end
        function Factorlist_parity = factorlist_parity(Degree)
            Nterm = nchoosek(4+Degree-1,Degree);
            Factorlist_parity = zeros(1,Nterm);
            for i =0:Degree
                if mod(i, 2) == 0
                    % number is even
                    Factorlist_parity(HK.orderlist(i)) = 1;
                else
                    % number is odd
                    Factorlist_parity(HK.orderlist(i)) = -1;
                end
                
            end
        end
        function sym_term = standardize_sym(sym_term)
            str_tmp = string(sym_term);
            str_tmp = replace(str_tmp,{'k_x','k_X','K_X','K_x','kx','kX','KX','Kx'},'k_x');
            str_tmp = replace(str_tmp,{'k_y','k_Y','K_y','K_Y','ky','kY','Ky','KY'},'k_y');
            str_tmp = replace(str_tmp,{'k_z','k_Z','K_z','K_Z','kz','kZ','Kz','KZ'},'k_z');
            sym_term = str2sym(str_tmp);
        end
        function [sym_coe,sym_single,str_single] = coeff_extract(sym_term,Dim)
%             str_list_tmp = strsplit(string(sym_term),'*');
%             sym_coe = sym(1);
%             for i = 1:length(str_list_tmp)
%                 if ~HK.strcontain(str_list_tmp{i},["k_x","k_y","k_z"])
%                     sym_coe = sym_coe * str2sym(str_list_tmp{i});
%                 end
%             end
%             sym_single = sym_term/sym_coe;
%             str_single = string(sym_single);
            % test use coeffs if bugs exist, try to umcomment above and comment below!
            arguments
                sym_term
                Dim = 3;
            end
            VarsSeqLcart = [sym('k_x'),sym('k_y'),sym('k_z'),sym('k_w')];
            [sym_coe,sym_single] = coeffs(sym_term,VarsSeqLcart(1:Dim));
            str_single = string(sym_single);
        end
    end
    %% protected method
    methods (Static,Hidden,Access= protected)
        function str_out = quadk2K(str_in)
            str_out = strrep(str_in,'x^2',"X");
            str_out = strrep(str_out,'y^2',"Y");
            str_out = strrep(str_out,'z^2',"Z");
            str_out = strrep(str_out,'x^3',"X*x");
            str_out = strrep(str_out,'y^3',"Y*y");
            str_out = strrep(str_out,'z^3',"Z*z");
            str_out = strrep(str_out,'x^4',"X^2");
            str_out = strrep(str_out,'y^4',"Y^2");
            str_out = strrep(str_out,'z^4',"Z^2");
            str_out = strrep(str_out,'x^5',"X^2*x");
            str_out = strrep(str_out,'y^5',"Y^2*y");
            str_out = strrep(str_out,'z^5',"Z^2*z");
            % can expand to n
        end
        function [vector_list,Coeffs_list] = HstrL_classify(H_strL_tmp,R_struct,mode)
            if nargin <3
                mode = 'T';
            end
            if strcmp(mode,'T')
                symvar_list = symvar(str2sym(H_strL_tmp));
                if isempty(symvar_list)
                    vector_list = [0,0,0];
                    Coeffs_list = 1;
                elseif length(symvar_list) == 1
                    [vector_list,Coeffs_list] = HK.Hstr_mapping(H_strL_tmp,R_struct,mode);
                else
                    [vector_list,Coeffs_list] = HK.Hstr_mapping(string(symvar_list(1)),R_struct,mode);
                    for i = 2:length(symvar_list)
                        [VL2,CL2] = HK.Hstr_mapping(string(symvar_list(i)),R_struct,mode);
                        [vector_list,Coeffs_list] = HK.VLCL_ltimes(vector_list,Coeffs_list,VL2,CL2);
                    end
                end
            elseif strcmp(mode,'H')
                [vector_list,Coeffs_list] = HK.Hstr_mapping(H_strL_tmp,R_struct,mode);
            elseif strcmp(mode,'H2')
                [vector_list,Coeffs_list] = HK.Hstr_mapping(H_strL_tmp,R_struct,mode);
                
            end
        end
        function [vector_list,Coeffs_list] = Hstr_mapping(H_str,R_struct,mode)
            if nargin <3
                mode = 'T';
            end
            if strcmp(mode,'T')
                switch H_str
                    case '1'
                        vector_list =  [0 0 0];
                        Coeffs_list = 1;
                    case 'x'
                        vector_list =  [1 0 0;-1 0 0];
                        Coeffs_list = (1/(2*R_struct.a*1i))*[1;-1];                        
                    case 'y'
                        vector_list =  [0 1 0;0 -1 0];
                        Coeffs_list = (1/(2*R_struct.b*1i))*[1;-1];
                    case 'z'
                        vector_list =  [0 0 1;0 0 -1];
                        Coeffs_list = (1/(2*R_struct.c*1i))*[1;-1];  
                    case 'X'
                        vector_list =  [0 0 0;1 0 0;-1 0 0];
                        Coeffs_list = (1/R_struct.a^2)*[2;-1;-1];                        
                    case 'Y'
                        vector_list =  [0 0 0;0 1 0;0 -1 0];
                        Coeffs_list = (1/R_struct.b^2)*[2;-1;-1];      
                    case 'Z'
                        vector_list =  [0 0 0;0 0 1;0 0 -1];
                        Coeffs_list = (1/R_struct.c^2)*[2;-1;-1];   
                    otherwise 
                        disp(H_str);
                        warning('the above Term cant be kp2tb, the result is wrong, please check!');
                end
            elseif strcmp(mode,'H')
                % A bug fixing
                %
                switch H_str
                    case '1'
                        vector_list =  [0 0 0];
                        Coeffs_list = 1;
                    case 'x'
                        vector_list =  [...
                            1 0 0 ;-1  0 0;...
                            0 1 0 ;0 -1 0;...
                            -1 -1 0;1  1 0;...
                            ];
                        Coeffs_list = (1/(6*R_struct.a*1i))*[...
                            2;-2;
                            (exp(-pi*2i/3)+exp(pi*2i/3));-(exp(-pi*2i/3)+exp(pi*2i/3));
                            (exp(-pi*4i/3)+exp(pi*4i/3));-(exp(-pi*4i/3)+exp(pi*4i/3));
                            ];
                    case 'y'
                        vector_list =  [...
                            0 1 0 ;0 -1 0;...
                            -1 -1 0;1  1 0;...
                            ];
                        Coeffs_list = (1/(6*R_struct.b*1i))*[...
                             -(exp(pi*1i/6)+exp(-pi*1i/6));(exp(pi*1i/6)+exp(-pi*1i/6));
                             -(exp(pi*5i/6)+exp(-pi*5i/6));(exp(pi*5i/6)+exp(-pi*5i/6));
                            ];
                    case 'z'
                        vector_list =  [0 0 1;0 0 -1];
                        Coeffs_list = (1/2*R_struct.c*1i)*[1;-1];
                    case 'X'
                        vector_list =  [...
                            0 0 0;1 0 0 ;-1  0 0;...
                            0 0 0;1  1 0;-1 -1 0;...
                            0 0 0;0 1 0 ;0 -1 0;...
                            ];
                        Coeffs_list = (1/(6*R_struct.a^2))*[...
                            12;-6;-6;...
                            0;0;0;...
                            0;0;0;...
                            ];
                    case 'Y'
                        vector_list =  [...
                            0 0 0;1 0 0 ;-1  0 0;...
                            0 0 0;1  1 0;-1 -1 0;...
                            0 0 0;0 1 0 ;0 -1 0;...
                            ];
                        Coeffs_list = (1/(6*R_struct.b^2))*[...
                            -4;2;2;...
                            8;-4;-4;...
                            8;-4;-4;...
                            ];
                    case 'x*y'
                        vector_list =  [...
                            1  1 0;-1 -1 0;...
                            0 1 0 ;0 -1 0;...
                            ];
                        Coeffs_list = (1/(3*R_struct.a*R_struct.a))*[
                            -sqrt(3);-sqrt(3);...
                            sqrt(3);sqrt(3);...
%                             -1/2-1i*sqrt(3)/2;-1/2-1i*sqrt(3)/2 ;...
%                             -1/2+1i*sqrt(3)/2;-1/2+1i*sqrt(3)/2;...
                            ];
%                             -(exp(pi*1i/6)) *[2 ;-1;-1];...
%                             -(exp(pi*5i/6)) *[2 ;-1;-1]];
                    case 'Z'
                        vector_list =  [0 0 0;0 0 1;0 0 -1];
                        Coeffs_list = (1/R_struct.c^2)*[2;-1;-1];
                    otherwise  
                        if HK.strcontain(H_str,"x*y") && ~strcmp(H_str,'x*y')
                            H_str_tmp = strrep(H_str,'x*y','1');
                            [vector_list,Coeffs_list] = HK.Hstr_mapping("x*y",R_struct,mode);
                            [VL2,CL2] = HK.Hstr_mapping(H_str_tmp,R_struct,mode);
                            [vector_list,Coeffs_list] =  HK.VLCL_ltimes(vector_list,Coeffs_list,VL2,CL2);
                            return;
                            %disp('i')
                        else
                            H_str_tmp = string(simplify(str2sym(H_str)));
                            symvar_list = symvar(str2sym(H_str));
                            if isempty(symvar_list)
                                vector_list = [0,0,0];
                                Coeffs_list = 1;
                            elseif length(symvar_list) == 1
                                [vector_list,Coeffs_list] = HK.Hstr_mapping(H_str_tmp,R_struct,mode);
                            else
                                [vector_list,Coeffs_list] = HK.Hstr_mapping(string(symvar_list(1)),R_struct,mode);
                                for i = 2:length(symvar_list)
                                    [VL2,CL2] = HK.Hstr_mapping(string(symvar_list(i)),R_struct,mode);
                                    [vector_list,Coeffs_list] = HK.VLCL_ltimes(vector_list,Coeffs_list,VL2,CL2);
                                end
                                
                            end
                        end
                end
            else
                disp('figure out the arbitrary kp2tb!!!');
            end
        end
        function [vector_list,Coeffs_list] = VLCL_ltimes(VL1,CL1,VL2,CL2)
            vector_list = zeros(size(CL1,1)*size(CL2,1),3);
            Coeffs_list = zeros(size(CL1,1)*size(CL2,1),1);
            count = 0;
            for i = 1:size(VL1,1)
                for j = 1:size(VL2,1)
                    count = count +1;
                    vector_list(count,:) = VL1(i,:)+VL2(j,:);
                    Coeffs_list(count) = CL1(i)*CL2(j);
                end
            end
            
        end
    end
end


