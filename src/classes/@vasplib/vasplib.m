classdef vasplib < matlab.mixin.CustomDisplay
    % Tools function Storage in vasplib
    %   此处显示详细说明

    %% public properties
    properties
        Dim = 3;
        Basis_num ;
        Rm = [];
        Hermitian= true;
    end
    %% Hsym and Hfun
    properties
        VarsSeqLcart = [];
        VarsSeqLfrac = [];
        Hsym;
    end
    %% private properties
    properties (GetAccess = protected,Hidden = true)
        
    end
    properties (GetAccess = protected)
        
    end
    %% temp properties
    properties (Transient,Hidden = true)


    end
    %% dynamic properties dependent
    properties(Dependent = true)
        Gk;
        symvar_list ; % the symvar in HcoeL
        Nbands ;
        Hfun;
    end
    properties(Dependent = true,Hidden = true)
         
    end
    %% hidden properties
    properties 
        Basis     ;
        % for orb
        orbL     = []; % the orbital list fractional
        elementL = []; % the element of each orbital
        quantumL = []; % [n,l,m,s] the three quantum number combined with spin
        orb_symL = sym([]); % orbital function, makes your own
        % for spacegroup
        sgn       = 1    ;
        Atom_name = []   ;
        Atom_num  = []   ;
        sites            ;
        symmetry_operation = [];
        klist_cart          ;
        klist_l          ;
        klist_frac          ;
        kpoints_frac        ;
        kpoints_l        ;
        kpoints_name     ;
        Rnn              ;
        nn_store         ;
        timtj            ;
        tjmti            ;
        Sparse_vector = [];
        N_Sparse_vector =[];
        CutList =[]      ;
    end
    %% Define which properties show
    methods (Access = protected)
        function propgrp = getPropertyGroups(~)
            proplist = {'Basis_num','Rm','Gk','Dim'};
            propgrp = matlab.mixin.util.PropertyGroup(proplist);
        end
    end
    %% method Constuction
    methods
        function vasplibobj =vasplib(propArgs)
            arguments
                propArgs.?vasplib;
            end
            Fieldnames = fieldnames(propArgs);
            if isempty(Fieldnames)
            else
                for iFieldname = 1:numel(Fieldnames)
                    vasplibobj.(Fieldnames{iFieldname}) = propArgs.(Fieldnames{iFieldname});
                end
            end
            vasplibobj = vasplibobj.vasplib_init();
        end
        function vasplibobj = vasplib_init(vasplibobj)
            if isempty(vasplibobj.Rm)
                vasplibobj.Rm = eye(vasplibobj.Dim);
            end
            if isempty(vasplibobj.VarsSeqLcart)
                syms k_x k_y k_z k_w real;
                vasplibobj.VarsSeqLcart = [k_x k_y k_z k_w];
            end
            if isempty(vasplibobj.VarsSeqLfrac)
                syms k_1 k_2 k_3 k_4 real;
                vasplibobj.VarsSeqLfrac = [k_1 k_2 k_3 k_4];
            end
        end
        function vasplibobj = vasplibCopy(vasplibobj,vasplibobj_in)
            CopyItem = [...
                "Rm","orbL","Dim","elementL","quantumL","orb_symL","sgn","Atom_name","Atom_num",...
                "sites","symmetry_operation","klist_cart","klist_l","klist_frac","kpoints_l","kpoints_name",...
                "Rnn","nn_store"...
                ];
            for i = 1:numel(CopyItem)
                vasplibobj.(CopyItem(i)) = vasplibobj_in.(CopyItem(i));
            end
        end
    end
    %% get
    methods
        function Hfun = get.Hfun(vasplibobj)
            Hfun = matlabFunction(vasplibobj.Hsym,'Vars',vasplibobj.VarsSeqLcart(1:vasplibobj.Dim));
        end
        function Gk = get.Gk(vasplibobj)
            Gk = (eye(length(vasplibobj.Rm))*2*pi/(vasplibobj.Rm)).';
        end
        function symvar_list = get.symvar_list(vasplibobj)
            symvar_list = symvar(vasplibobj.HcoeL);
        end
        function Nbands = get.Nbands(vasplibobj)
            if isfield(vasplibobj,"WAN_NUM")
                Nbands = vasplibobj.WAN_NUM;
            else
                Nbands = vasplibobj.Basis_num;
            end
        end
    end
    %% (ti tj)
    methods
        function vasplibobj = timtj_gen(vasplibobj,mode)
            if nargin < 2
                mode = 'num';
            end
            % {1} cartesian
            % {2} fractional
            % {3} symbolic term cartesian
            % {4} symbolic term fractional
            % Prepare ti - tj
            vasplibobj.timtj{1} = zeros([size(vasplibobj.orbL,1),size(vasplibobj.orbL,1),vasplibobj.Dim]);
            if strcmp(mode,'sym')
                vasplibobj.timtj{1} = sym( vasplibobj.timtj{1});
                vasplibobj.timtj{3} = sym(zeros(size(vasplibobj.orbL,1)));
                %syms k_x k_y k_z real;
                %syms k_1 k_2 k_3 real;
            end
            vasplibobj.timtj{2} = vasplibobj.timtj{1};
            for i = 1:size(vasplibobj.orbL,1)
                ti= vasplibobj.orbL(i,:);
                for j = 1:size(vasplibobj.orbL,1)
                    tj= vasplibobj.orbL(j,:);
                    tij = ti-tj;
                    tij_cart = tij*vasplibobj.Rm;
                    for d = 1:vasplibobj.Dim
                        vasplibobj.timtj{1}(i,j,d) = tij_cart(d);
                        vasplibobj.timtj{2}(i,j,d) = tij(d);
                    end
                end
            end
            if strcmp(mode,'sym')
                ExpInnerTerm = vasplib.matrixtimespage(vasplibobj.VarsSeqLcart(1:vasplibobj.Dim),vasplibobj.timtj{1});
                vasplibobj.timtj{3} = exp(1i*(sum(ExpInnerTerm,3)));
                ExpInnerTermFrac = vasplib.matrixtimespage(vasplibobj.VarsSeqLcart(1:vasplibobj.Dim),vasplibobj.timtj{2});
                vasplibobj.timtj{4} =  exp(1i*(sum(ExpInnerTermFrac,3)));
            else
                vasplibobj.timtj{3} = [];
                vasplibobj.timtj{4} = [];
            end
        end
        function vasplibobj = tjmti_gen(vasplibobj,mode)
            if nargin < 2
                mode = 'num';
            end
            vasplibobj = vasplibobj.timtj_gen(mode);
            % Prepare tj - ti
            for i = 1:2
                vasplibobj.tjmti{i} = - vasplibobj.timtj{i};
            end
            if strcmp(mode,'sym')
                ExpInnerTerm = vasplib.matrixtimespage(vasplibobj.VarsSeqLcart(1:vasplibobj.Dim),vasplibobj.tjmti{1});
                vasplibobj.tjmti{3} = exp(1i*(sum(ExpInnerTerm,3)));
                ExpInnerTermFrac = vasplib.matrixtimespage(vasplibobj.VarsSeqLcart(1:vasplibobj.Dim),vasplibobj.tjmti{2});
                vasplibobj.tjmti{4} =  exp(1i*(sum(ExpInnerTermFrac,3)));
            else
                vasplibobj.tjmti{3} = [];
                vasplibobj.tjmti{4} = [];
            end
        end
        function vasplibobj = SliceGen(vasplibobj)
            % $H_{i j}^{\mathbf{k}}=\left\langle\chi_{i}^{\mathbf{k}}|H| \chi_{j}^{\mathbf{k}}\right\rangle=\sum_{\mathbf{R}} e^{i \mathbf{k} \cdot\left(\mathbf{R}+\mathbf{t}_{j}-\mathbf{t}_{i}\right)} H_{i j}(\mathbf{R})$
            DIM = vasplibobj.Dim;
            if strcmp(vasplibobj.Type,'list')
                switch class(vasplibobj)
                    case 'Htrig'
                        Kind = size(vasplibobj.HsymL_numL,1);
                        [ij_list,index_row] = sortrows(vasplibobj.HsymL_numL(:,DIM+1:DIM+2));
                    case 'HR'
                        Kind = vasplibobj.NRPTS;
                        [ij_list,index_row] = sortrows(vasplibobj.vectorL(:,DIM+1:DIM+2));
                end
                vasplibobj = vasplibobj.reseq(':',index_row);% ?
                [vasplibobj.Sparse_vector,SliceList] = unique(ij_list,'rows');
                vasplibobj.N_Sparse_vector = size(vasplibobj.Sparse_vector,1);
                vasplibobj.CutList = [SliceList,[SliceList(2:end)-1;Kind]];
            end
        end
    end
    %% information
    methods
        function vasplibobj = input_Rm(vasplibobj,Rm)
            if nargin <2 && exist('POSCAR','file')
                [Rm,~,~,~,~] = vasplib.POSCAR_read('POSCAR','vasp');
            elseif nargin <2 && ~exist('POSCAR','file')
                Rm = eye(vasplibobj.Dim);
                %warning('POSCAR or Rm needed');
            else
                if isa(Rm,'string') || isa(Rm,'char')
                    [Rm,~,~,~,~] = vasplibobj.POSCAR_read(Rm,'vasp');
                end
            end
            Rm = (Rm);
            vasplibobj.Rm = Rm;
        end
        function vasplibobj = input_orb_struct(vasplibobj,filename,mode,options)
            arguments
                vasplibobj vasplib;
                filename string ='POSCAR';
                mode char {mustBeMember(mode,{'vasp','tbsk','sym'})} = 'vasp';
                options.symbolic logical = false;
                options.Operation logical = false;
                options.warning logical = true;
                options.spin char {mustBeMember(options.spin,{'spinless','wannier','block'})}= 'spinless';% spinless;block;wannier;
            end
            if options.Operation
                import spglib_matlab.*;
            end
            [vasplibobj.Rm,tmpsites,vasplibobj.Atom_name,vasplibobj.Atom_num,elements]=vasplibobj.POSCAR_read(filename,mode);

            % error check
            if options.warning
                if sum(vasplibobj.Atom_num) ~= vasplibobj.Basis_num
                    if isa(vasplibobj,'HR')
                        warning('The sum of POSCAR 7th line is different with WAN_NUM you defined.\n We use WAN_NUM enforcedly');
                    else
                        warning('The sum of POSCAR 7th line is different with Basis_num you defined.\n We use Basis_num enforcedly');
                    end
                end
                if length(tmpsites) ~= vasplibobj.Basis_num
                    if isa(vasplibobj,'HR')
                        warning('The sum of POSCAR 9~end lines is different with WAN_NUM you defined.\n We use WAN_NUM enforcedly');
                    else
                        warning('The sum of POSCAR 9~end lines is different with Basis_num you defined.\n We use Basis_num enforcedly');
                    end
                end
            end
            % reset tmpsites
            switch mode
                case 'vasp'

                case {'tbsk','sym'}
                    %if length(tmpsites) < vasplibobj.Basis_num
                    tmpsites = vasplib.MakeUpSites(tmpsites,options.spin);
                    %end
            end
            vasplibobj.orbL = [[tmpsites.rc1].',[tmpsites.rc2].',[tmpsites.rc3].'];
            if options.Operation
                % position = [[vasplibobj.sites.rc1].',[vasplibobj.sites.rc2].',[vasplibobj.sites.rc3].'].';
                types = [];
                for i = 1:length(vasplibobj.Atom_num)
                    types = [types,ones(1,vasplibobj.Atom_num(i)*i)];
                end
                spglib_database =  spglib_matlab.spg_get_dataset_from_sites(vasplibobj.Rm,vasplibobj.orbL.',types);
                vasplibobj.sgn = spglib_database.spacegroup_number;
                vasplibobj.SymmetryOperation = [spglib_database.rotations;spglib_database.translations];
            end
            %              = tmp_orbL(1:vasplibobj.Bassis_num,:);
            if ~isempty(elements)
                elements.Properties.RowNames = elements.atom_symbol;
            end
            switch mode
                case 'vasp'
                    vasplibobj.sites = tmpsites;
                case 'tbsk'
                    tmp_orb_symL = [tmpsites.orb_sym].';
                    if isempty(vasplibobj.Basis_num)
                        vasplibobj.Basis_num = size(vasplibobj.orbL,1);
                    end
                    for i = 1:vasplibobj.Basis_num
                        tmpdata = table2array(elements(char(tmpsites(i).element),{'atom_number','n'}));
                        element_seq = tmpdata(1);
                        orb_symbk = tmp_orb_symL(i,:);
                        n = tmpdata(2);
                        l = vasplib.orb2l(tmpsites(i).orb);
                        m =  vasplib.orb_sym2m(orb_symbk);
                        vasplibobj.elementL(i,:) = element_seq;
                        vasplibobj.quantumL(i,1) = n; % element
                        vasplibobj.quantumL(i,2) = l; % L
                        vasplibobj.quantumL(i,3) = m; % lz
                        vasplibobj.quantumL(i,4) = tmpsites(i).spin;
                        try
                            vasplibobj.orb_symL(i,:) = vasplib.Ymlsym(l,m,orb_symbk);
                        catch
                        end
                    end
                    vasplibobj.sites = tmpsites;
                case 'sym'
                    tmp_orb_symL = [tmpsites.orb_sym].';
                    for i = 1:vasplibobj.Basis_num
                        tmpdata = table2array(elements(char(tmpsites(i).element),{'atom_number','n'}));
                        element_seq = tmpdata(1);
                        orb_symbk = tmp_orb_symL(i,:);
                        n = tmpdata(2);
                        l = vasplib.orb2l(tmpsites(i).orb);
                        m =  vasplib.orb_sym2m(orb_symbk);
                        vasplibobj.elementL(i,:) = element_seq;
                        vasplibobj.quantumL(i,1) = n; % element
                        vasplibobj.quantumL(i,2) = l; % L
                        vasplibobj.quantumL(i,3) = m; % lz
                        vasplibobj.quantumL(i,4) = tmpsites(i).spin; % sz/Jz
                        vasplibobj.quantumL(i,5) = tmpsites(i).J;   % J
                        vasplibobj.orb_symL(i,:) = vasplib.Ymlsym(l,m,orb_symbk);
                    end

            end
            if options.symbolic
                recommend_accuracy = false;
                Rm_tmp = sym(vasplibobj.Rm);
                for i = 1:numel(Rm_tmp)
                    tmpstr =char(string(Rm_tmp(i)));
                    if length(tmpstr)>15
                        [k,j] = ind2sub([3,3],i);
                        fprintf('Rm(%d,%d) %f seems not reach the accuracy :%s\n',k,j,vasplibobj.Rm(k,j),tmpstr);
                        recommend_accuracy = true;
                    end
                end
                orb_tmp = sym(vasplibobj.orbL);
                sizeorbL = size(orb_tmp);
                for i = 1:numel(orb_tmp)
                    tmpstr =char(string(orb_tmp(i)));
                    if length(tmpstr)>15
                        [k,j] = ind2sub(sizeorbL,i);
                        fprintf('orbL(%d,%d) %f seems not reach the accuracy :%s\n',k,j,orb_tmp(k,j),tmpstr);
                        recommend_accuracy = true;
                    end
                end
                if recommend_accuracy
                    fprintf('************* warning *****************\n');
                    fprintf('the accuracy of the digits should  reach to 1e-13. eg.\n');
                    fprintf('0.333333333333333 -> %s\n',string(sym(0.333333333333333)));
                    fprintf('0.866025403784439 -> %s\n',string(sym(0.866025403784439)));
                    fprintf('You can use ''format long'', test all the digits in input struct.\n');
                end
                vasplibobj.Rm = Rm_tmp;
                vasplibobj.orbL = orb_tmp;
            end
        end
        function vasplibobj = nn(vasplibobj,search_range,Accuracy,Rlength_cut,options)
            %%
            % Caculate the nn for a primitive cell
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
            %       the Accuracy takes weight
            % Debug: Unsupport For Dim
            arguments
                vasplibobj ;
                search_range = zeros(1,vasplibobj.Dim);
                Accuracy = 1e-4;
                Rlength_cut = 5;
                options.MAX_NN_LENGTH =  10000000;
                options.onsite = false;
            end
            MAX_NN_LENGTH = options.MAX_NN_LENGTH;
            %% -------- init --------
            if isempty(vasplibobj.orbL)
                error('no orbL, please load orbL.');
            end
            if Accuracy > 1
                warning('Attention，your ACCURACY may set wrongly!')
            end
            sites_num = size(vasplibobj.orbL,1);
            search_rangex=search_range(1);
            search_rangey=search_range(2);
            search_rangez=search_range(3);
            % for sparse mode, space is the only requirement
            MAX_NN_LENGTH2 = sites_num^2*(search_rangex+1)^2*(search_rangey+1)^2*(search_rangez+1)^2;
            max_nn_length = min(MAX_NN_LENGTH2,MAX_NN_LENGTH);
            nn_sparse = zeros(max_nn_length,10);
            % i j rcart1 rcart2 rcart3 R1 R2 R3 Rlength nn_level
            Rnn_list = [];
            count  = 1;
            if isa(vasplibobj.Rm,'sym')
                sym_mode = true;
            else
                sym_mode = false;
            end
            if sym_mode
                Accuracy = -1;
                Rm_ = (vasplibobj.Rm);
            else
                Rm_ = double(vasplibobj.Rm);
            end

            % -------- save ------------
            pb = vasplib_tool_outer.CmdLineProgressBar('Generating ');
            for j  = 1:sites_num
                orb2 = double(vasplibobj.orbL(j,:));
                for i = 1:sites_num
                    orb1 = double(vasplibobj.orbL(i,:));
                    [nn_sparse_temp,Rnn_list_temp] = ...
                        vasplibobj.nn_sparse_gen(orb1,orb2,Rm_,search_rangex,search_rangey,search_rangez,Accuracy,Rlength_cut,options.onsite);
                    countadd = size(nn_sparse_temp,1);
                    nn_sparse_temp(:,1) = i;
                    nn_sparse_temp(:,2) = j;
                    nn_sparse(count:count+countadd-1,:) = nn_sparse_temp;
                    count = count+countadd;
                    Rnn_list = [Rnn_list;Rnn_list_temp];
                end
                pb.print(j,sites_num,' th orb nn information ...');
            end
            pb.delete();
            %% -------- clean --------------
            nn_sparse(count:max_nn_length,:) = [];
            %% -------- caculate ------------
            vasplibobj.Rnn = sort(unique(Rnn_list,'row'));
            %% -------- make map ------------
            %             H_hr.Rnn_map = containers.Map('KeyType','double','ValueType','double');
            %             for i = 1:length(vasplibobj.Rnn)
            %                 H_hr.Rnn_map(vasplibobj.Rnn(i)) = i;
            %             end
            %%  give level
            if options.onsite
                level_bias = 1;
            else
                level_bias = 0;
            end
            pb = vasplib_tool_outer.CmdLineProgressBar('Giving nn_level for ');
            for i  = 1:count-1
                pb.print(i,count-1,' th hopping  ...');
                nn_sparse(i,10) = find(double(nn_sparse(i,9)) == vasplibobj.Rnn);
            end
            pb.delete();
            nn_sparse(:,10) = nn_sparse(:,10)-level_bias;
            vasplibobj.nn_store = nn_sparse;
            %             if sym_mode
            %                 vasplibobj.nn_store = sym(nn_sparse);
            %             else
            %                 vasplibobj.nn_store = nn_sparse;
            %             end
        end
        function [klist_l,kpoints_l,kpoints_name] = kpath_information(vasplibobj)
            klist_l = vasplibobj.klist_l;
            kpoints_l = vasplibobj.kpoints_l;
            kpoints_name = vasplibobj.kpoints_name;
        end
        function vasplibobj = kpathgen3D(vasplibobj,KPOINTS_name,nodes,kpoints_name_tmp)
            %% To generate klist for caculation and plot
            %
            % creat kpath from KPOINTS or otherfile
            %
            %% Description of the Function:
            %%
            %% Usage:
            %
            % * [klist_cart,klist_l,klist_frac,kpoints_l,kpoints_name]=kpathgen3D(Rm,kpoints,nodes,kpoints_name)
            % * [klist_cart,klist_l,klist_frac,kpoints_l,kpoints_name]=kpathgen3D(Rm,kpoints,nodes)
            % * [klist_cart,klist_l,klist_frac,kpoints_l,kpoints_name]=kpathgen3D(Rm,kpoints)
            % * [klist_cart,klist_l,klist_frac,kpoints_l,kpoints_name]=kpathgen3D(Rm)
            %
            %% Input:
            %
            % # input1:
            % # input2:
            % # kpoints:
            %kpoints=[ k1x k1y k1z;...
            %          k2x k2y k2z;...
            %          k2x k2y k2z;...
            %          k3x k3y k3z;...
            %                .
            %                .
            %                .
            %          kn-1x kn-1y kn-1z;]
            %          kn-1x kn-1y kn-1z;]
            %          knx kny knz;]
            % # nodes: (nodes along each k-path)
            % # kpoints_name: [K1;K2;K3]
            %
            %% Output:
            %
            % # klist_l      (for plot : the x-> will be klist_l for each band ,one-dimension)
            % # kpoints_l    (for plot : To gen high symmetry points along k-path,one-dimension)
            % # kpoints_name (for plot : To gen high symmetry points name along k-path,one-dimension)
            % # klist_cart   (for caculation: cati real)
            % # klist_frac   (for other function you need )
            %
            %% example:
            %   commmad
            %
            %
            %   result
            %
            %% Note:
            % if [klist_cart,klist_l,klist_frac,kpoints_l,kpoints_name]=kpathgen3D()
            %          the kpoints , nodes ,kpoints_name will be given by KPOINTS POSCAR file
            %
            % if [klist_cart,klist_l,klist_frac,kpoints_l,kpoints_name]=kpathgen3D(Rm)
            %          the kpoints , nodes ,kpoints_name will be given by KPOINTS file
            % note : ['Line-Mode' ],['Reciprocal' ]}
            %
            %% Change log
            %
            % * Document Date: 2020/12/03
            % * Creation Date: 2020/12/03
            % * Last updated : 2022/11/06
            %
            %% Copyright
            %
            % * parkman
            % * <parkman@buaa.edu.cn>
            %
            arguments
                vasplibobj
                KPOINTS_name = 'KPOINTS';
                nodes = [];
                kpoints_name_tmp = [];
            end
            if isempty(vasplibobj.Rm)
                vasplibobj = vasplibobj.input_Rm();
            end
            if isempty(nodes)
                [kpoints,nodes,kpoints_name_tmp] = KPOINTS_read(KPOINTS_name);
            else
                kpoints = KPOINTS_name;
            end
            [vasplibobj.klist_cart,vasplibobj.klist_frac,vasplibobj.klist_l,vasplibobj.kpoints_l,vasplibobj.kpoints_frac] =...
                vasplib.kpathgen(kpoints,nodes,vasplibobj.Gk);
            %high symmetry points store
            vasplibobj.kpoints_name = kpoints_name_tmp;
        end
        function [Rm,sites,Atom_name,Atom_num] = POSCAR_gen(vasplibobj,filename,Rm,sites,Atom_name,Atom_num)
            %% note: POSCAR is Direct mode
            % input : Stucture information
            %                             a_crystal_constance
            %                             a1,a2,a3,Rm=[a1;a2;a3]
            %                             atom_name;atom_num
            %                             sites
            %          element_information:
            %                             sites
            % Output : POSCAR:
            %                             a_crystal_constance
            %                             a1,a2,a3,Rm=[a1;a2;a3]
            %                             atom_name;atom_num
            %                             coordinates pattern

            % note : site=struct('seq',[],'inseq',[],'rc1',[],'rc2',[],'rc3',[],'name',[],'nameseq',[]);
            % Rm
            %a1 a1x a1y a1z
            %a2 a2x a2y a2z
            %a3 a3x a3y a3z
            % unit A?e-6
            % usage : POSCAR_gen(Rm,sites,Atom_name,Atom_num,filename);
            arguments
                vasplibobj
                filename = 'POSCAR';
                Rm = vasplibobj.Rm;
                sites = vasplibobj.sites;
                Atom_name = vasplibobj.Atom_name;
                Atom_num = vasplibobj.Atom_num;
            end
            title = "POSCAR Gen by vasplib";
            %% write POSCAR
            % Initialize variables
            %filename = "POSCAR_"+num2str(term2)+".vasp";
            fileID = fopen(filename,'w');
            %% Rm
            a_crystal_constance=1;
            %%
            fprintf(fileID,"%s\n",title);
            fprintf(fileID,"%d\n",a_crystal_constance);
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
            [~,sites_num]=size(sites);
            for i=1:sites_num
                fprintf(fileID,"%f  ",mod(sites(i).rc1,1));
                fprintf(fileID,"%f  ",mod(sites(i).rc2,1));
                fprintf(fileID,"%f  ",mod(sites(i).rc3,1));
                fprintf(fileID,"\n  ");
            end
            fclose(fileID);
        end
        function [klist_dos,klist1,klist2]=kmesh3D(vasplibobj,mesh,kz,mode,options)
            %% To generate k-mesh for DOSplot or 2D-Bandplot
            % usage [klist_dos]=kmesh3D(Rm,mesh)
            % input  : mesh [meshx meshy meshz]
            %          Rm
            %a1 a1x a1y a1z
            %a2 a2x a2y a2z
            %a3 a3x a3y a3z
            % output : mesh [meshx meshy meshz]
            arguments
                vasplibobj ;
                mesh = [21,21];
                kz = [];
                mode = 'bulk-DOS';
                options.KPOINTS = '';
            end
            %--------  nargin  --------
            if isempty(kz) && strcmp(mode,'bulk-DOS')
            elseif strcmp(mode,'bulk-DOS')
                mode = '3Dbandplot';
            end
            %--------  init  --------
            %             vasplibobj.Gk  = (eye(3)*2*pi/vasplibobj.Rm)' ;% Reciprocal vector
            %--------  distri  --------
            switch mode
                case 'bulk-DOS'
                    b_x=sum(vasplibobj.Gk(:,1));
                    b_y=sum(vasplibobj.Gk(:,2));
                    b_z=sum(vasplibobj.Gk(:,3));
                    xnodes=mesh(1);
                    ynodes=mesh(2);
                    znodes=mesh(3);
                    xlist=linspace(-b_x,b_x,xnodes);
                    ylist=linspace(-b_y,b_y,ynodes);
                    zlist=linspace(-b_z,b_z,znodes);
                    klist_dos=[];
                    for i=1:xnodes
                        for j=1:ynodes
                            for k=1:znodes
                                klist_temp=[xlist(i) ylist(j) zlist(k)];
                                klist_dos=[klist_dos;klist_temp];
                            end
                        end
                    end
                    klist1 = xlist;
                    klist2 = ylist;
                case '3Dbandplot'
                    disp('Primitive cell')
                    disp(vasplibobj.Rm);
                    disp('Reciprocal Lattice')
                    disp(vasplibobj.Gk);
                    disp("Rm*Gk' ");
                    disp(sym(vasplibobj.Rm*vasplibobj.Gk'));
                    xnodes=mesh(1);
                    ynodes=mesh(2);
                    klist1=linspace(-0.5,0.5,xnodes);
                    klist2=linspace(-0.5,0.5,ynodes);
                    %     kz_r = kz*Gk(3,:);
                    klist_dos=[];
                    for i=1:xnodes
                        for j=1:ynodes
                            klist_temp=[klist1(i) klist2(j) kz]*vasplibobj.Gk;
                            klist_dos=[klist_dos;klist_temp];
                        end
                    end
                case 'fermiarc'
                    %disp('fermi-arc');
                    xnodes = mesh(1);
                    ynodes = mesh(2);
                    kfermiarc = kz;
                    klist1 = [linspace(0,kfermiarc(2,1),xnodes)',...
                        linspace(0,kfermiarc(2,2),xnodes)',...
                        linspace(0,kfermiarc(2,3),xnodes)',...
                        ];
                    klist2 = [linspace(0,kfermiarc(3,1),ynodes)',...
                        linspace(0,kfermiarc(3,2),ynodes)',...
                        linspace(0,kfermiarc(3,3),ynodes)',...
                        ];
                    klist_1 = klist1 ;
                    klist_2 = klist2 ;
                    %         disp( klist_1 );
                    %         disp( klist_2 );
                    klist_dos = zeros(xnodes*ynodes,3);
                    count =0;
                    for i=1:xnodes
                        for j=1:ynodes

                            count = count+1;
                            klist_dos(count,:) = klist_1(i,:)+klist_2(j,:)+ kfermiarc(1,:);
                        end
                    end
                    klist1  = klist_1 + kfermiarc(1,:);
                    klist2  = klist_2 + kfermiarc(1,:);
                    %         klist1 = linspace(kfermiarc(1,1),max(kfermiarc(2,:)),xnodes)';
                    %         klist2 = linspace(kfermiarc(1,2),max(kfermiarc(3,:)),xnodes)';
            end
            if strcmp(options.KPOINTS,'')
            else
                fileID = fopen(options.KPOINTS,'w');
                for i = 1:size(klist1,1)
                    fprintf(fileID,'%7.5f %7.5f %7.5f\n',klist2(1,1)+klist1(i,1),klist2(1,2)+klist1(i,2),klist2(1,3)+klist1(i,3));
                    fprintf(fileID,'%7.5f %7.5f %7.5f\n',klist2(end,1)+klist1(i,1),klist2(end,2)+klist1(i,2),klist2(end,3)+klist1(i,3));
                    fprintf(fileID,'\n');
                end
                fclose(fileID);
            end
        end
        function [num_label,coe_label] = NumOrCoe(vasplibobj)
            if isequal(zeros(size(vasplibobj.HnumL)),vasplibobj.HnumL)
                num_label = false;
            else
                num_label = true;
            end
            if isequal(sym(zeros(size(vasplibobj.HcoeL))),vasplibobj.HcoeL)
                coe_label = false;
            else
                coe_label = true;
            end
        end
        function varargout = klist_show(vasplibobj,options)
            arguments
                vasplibobj vasplib;
                options.klist double = vasplibobj.klist_cart;
                options.ax = handle([]);
                options.show = true;
                options.color = [rand rand rand];
                options.LineWidth = 2;
            end
            if options.show
                if isempty(options.ax)
                    [~,ax] = vasplib.BZplot(vasplibobj.Rm);
                else
                    ax = options.ax;
                end
            end
            if options.show
                plot3(options.klist(:,1),options.klist(:,2),options.klist(:,3),...
                    'Color',options.color,...
                    'DisplayName','k-path',...
                    'LineWidth',options.LineWidth);
                try % test
                    for i = 1:length(vasplibobj.kpoints_name)-1
                        kpoints_r = vasplibobj.kpoints_frac(2*i-1,:)*vasplibobj.Gk;
                        text(ax,kpoints_r(1),kpoints_r(2),kpoints_r(3),vasplibobj.kpoints_name(i),'FontSize',24);
                    end
                    kpoints_r = vasplibobj.kpoints_frac(end,:)*vasplibobj.Gk;
                    text(ax,kpoints_r(1),kpoints_r(2),kpoints_r(3),vasplibobj.kpoints_name(end),'FontSize',24);
                catch

                end
                varargout{1} = ax;
            end
        end
        function varargout = PARCHG_gen(varargin)
            varargout{:} = PARCHG_gen(varargin{:});
        end
    end
    methods(Static)
        function [kloop1_frac,kloop2_frac,kloop1_cart,kloop2_cart,klist_l,kstart_frac,kstart_cart] = kloop2D(Rm,options)
            arguments
                Rm
                options.knum_evol   = 51;
                options.knum_int   = 51;
                options.kstart  = [-0.5,-0.5,0];
                options.kevolution   = [1,0,0];
                options.kintegral   = [0,1,0];
                options.cartesian = false;
                options.dir_seq = [1,2,3];
                options.dir_start = 'kcar';
                options.shift = true;
            end
            if ~isa(Rm,'double')
                Rm = Rm.Rm;
            end
            Gk = (2*pi*eye(size(Rm,1))/Rm).';
            if size(Gk,1) < 3 || size(Gk,2) <3
                Gk(3,3) = 1;
            end
            if options.cartesian
                Gk_ = vasplib.CartisianMat(Gk,options.dir_seq,options.dir_start);
                kstart_frac  = options.kstart * Gk_ /Gk;
            else
                Gk_ = Gk;
                kstart_frac = options.kstart;
            end
            %
            knum1 = options.knum_evol;
            knum2 = options.knum_int;
            %
            [kloop1_cart,kloop1_frac,~,~] =...
                vasplib.kpathgen([[0,0,0];options.kevolution],knum1,Gk_,Gk);
            [kloop2_cart,kloop2_frac,~,~] =...
                vasplib.kpathgen([[0,0,0];options.kintegral],knum2,Gk_,Gk);
            %
            klist_l = zeros(knum1,1);
            normklist_l = norm(options.kevolution)/norm(kloop1_frac(end,:));
            for i = 1:knum1
                klist_l(i) = norm(kloop1_frac(i,:)*(eye(3)*2*pi))*normklist_l;
            end
            klist_l = klist_l + sum(sign(kstart_frac))* norm(kstart_frac*(eye(3)*2*pi))*normklist_l;
            kstart_cart = options.kstart*Gk_;
        end
        function [klist_cart,klist_frac,klist_cart_plot,sizemesh,Gk_,Grid] = kmesh2D(Rm,options,optionsEdge)
            arguments
                Rm =POSCAR_read;
                options.knum1   = 51;
                options.knum2   = 51;
                options.kstart  = [-0.5,-0.5,0];
                options.kdir1   = [1,0,0];
                options.kdir2   = [0,1,0];
                options.cartesian = false;
                options.dir_seq = [1,2,3];
                options.dir_start = 'kcar';
                options.shift = true;
                optionsEdge.full = false;
            end
            if ~isa(Rm,'double') && ~isa(Rm,'sym')
                Rm = Rm.Rm;
            end
            Gk = (2*pi*eye(size(Rm,1))/Rm).';
            if size(Gk,1) < 3 || size(Gk,2) <3
                Gk(3,3) = 1;
            end
            %
            if options.cartesian
                Gk_ = vasplib.CartisianMat(Gk,options.dir_seq,options.dir_start);
                kstart_s  = options.kstart * Gk_ /Gk;
            else
                Gk_ = Gk;
                kstart_s = options.kstart;
            end
            %
            knum1 = options.knum1+1;
            knum2 = options.knum2+1;
            [~,klist_s_1,~,~] =...
                vasplib.kpathgen([[0,0,0];options.kdir1],knum1,Gk_,Gk);
            %             klist_l = zeros(size(klist_s_1,1),1);
            %             klist_l(1) = sum(sign(klist_s_1(1,:)))*norm(klist_s_1(1,:)*(eye(3)*2*pi));
            %klist_s_1_ = klist_s_1 ;
            %normklist_l = norm(options.kdir1)/norm(klist_s_1(end,:));
            %             for i = 1:size(klist_s_1_,1)
            %                 klist_l(i) = norm(klist_s_1_(i,:)*(eye(3)*2*pi))*normklist_l;
            %             end
            %             klist_l = klist_l + sum(sign(kstart_s))* norm(kstart_s*(eye(3)*2*pi))*normklist_l;
            [~,klist_s_2,~,~] =...
                vasplib.kpathgen([[0,0,0];options.kdir2],knum2,Gk_,Gk);
            % half-closed half-open
            dk_s_1 = options.kdir1/(options.knum1-1);
            dk_s_2 = options.kdir2/(options.knum2-1);
            dk_s_start = dk_s_1/2 + dk_s_2/2;
            if optionsEdge.full
                klist_s_1_ = klist_s_1(1:knum1-1,:);
                klist_s_2_ = klist_s_2(1:knum2-1,:);
                knum1_ = knum1-1;
                knum2_ = knum2-1;
                sizemesh = [knum1_,knum2_];
            else
                klist_s_1 = klist_s_1(1:knum1-1,:);
                klist_s_2 = klist_s_2(1:knum2-1,:);
                knum1 = knum1-1;
                knum2 = knum2-1;
                sizemesh = [knum1,knum2];
            end
            %BC_WAVECAR = zeros(vasplibobj.Basis_num,options.knum1*options.knum2);
            %BCCAR = zeros(options.knum1,options.knum_evol);
            %BF_WAVECAR = zeros(length(BAND_index),length(BAND_index),options.knum_evol);
            %kstart_r = options.kstart*Gk_;
            % make klist
            %
            klist_frac = repmat(klist_s_1,[ knum2 1 ])+kron(klist_s_2,ones( knum1, 1 )) + kstart_s;
            klist_cart = klist_frac*Gk_;
            if options.shift
                shift = dk_s_start*Gk_;
            else
                shift = [0,0,0];
            end
            if optionsEdge.full
                klist_frac_ = repmat(klist_s_1_,[ knum2_ 1 ])+kron(klist_s_2_,ones( knum1_, 1 )) + kstart_s;
                klist_cart_ = klist_frac_*Gk_;
                klist_cart_plot = klist_cart_ + shift;
            else
                klist_cart_plot = klist_cart + shift;
            end
            Grid(:,:,1)  = reshape(klist_cart_plot(:,1),sizemesh);
            Grid(:,:,2)  = reshape(klist_cart_plot(:,2),sizemesh);
            Grid(:,:,3)  = reshape(klist_cart_plot(:,3),sizemesh);
            
        end
        function [klist_cart,klist_frac] = kloop1D(kpoint_frac,Orientation,radius,opt)
            arguments
                kpoint_frac;
                Orientation = [0,0,1];
                radius = 0.01;
                opt.nk = 100;
                opt.Gk = [];
                opt.inputCar logical= false
                opt.enforceCar logical= false
                opt.theta_0 = 0;
                opt.anticlockwise = true;
            end
            if isempty(opt.Gk)
                Rm = vasplib.POSCAR_read;
                opt.Gk = (2*pi*eye(3)/Rm).';
            end
            if opt.inputCar
                kpoint_cart = kpoint_frac;
                Orientation_cart = Orientation;
                radius_cart = radius;
            else
                kpoint_cart = kpoint_frac*opt.Gk;
                Orientation_cart = Orientation*opt.Gk;
                radius_cart = norm(opt.Gk)*radius;
            end
            if opt.enforceCar
                %Orientation_r = Orientation;
            end
            %
            % centered on (x0, y0, z0), with a radius of r, and with
            % normal vector n
            n=Orientation_cart; %法向量n
            r=radius_cart; %圆的半径为1
            c=kpoint_cart; %圆心的坐标
            if opt.anticlockwise
                theta=linspace(opt.theta_0,opt.theta_0+2*pi,opt.nk)'; %theta角从0到2*pi
            else
                theta=linspace(opt.theta_0,opt.theta_0-2*pi,opt.nk)'; %theta角从0到2*pi
            end
            a=cross([0 1 0],n); %n与j叉乘，求取a向量
            if ~any(a) %如果a为零向量，将n与j叉乘
                a=cross(n,[1 0 0]);
            end
            b=cross(n,a); %求取b向量
            a=a/norm(a); %单位化a向量
            b=b/norm(b); %单位化b向量
            c1=c(1)*ones(size(theta,1),1);
            c2=c(2)*ones(size(theta,1),1);
            c3=c(3)*ones(size(theta,1),1);
            k_xL=c1+r*a(1)*cos(theta)+r*b(1)*sin(theta);%圆上各点的x坐标
            k_yL=c2+r*a(2)*cos(theta)+r*b(2)*sin(theta);%圆上各点的y坐标
            k_zL=c3+r*a(3)*cos(theta)+r*b(3)*sin(theta);%圆上各点的z坐标
            %h = plotg(k_xL,k_yL,k_zL);
            %set(h,'LineWidth',3);
            klist_cart = [k_xL,k_yL,k_zL];
            klist_frac = klist_cart/opt.Gk;
        end
    end
    %% wt
    methods
        function fid = kpath_card_gen(vasplibobj,mode,filename)
            arguments
                vasplibobj ;
                mode = 'wt';
                filename = 'wt.in';
            end
            % function string
            if strcmp(mode,'wt')
                fid = fopen(filename,'a');
                head_string="KPATH_BULK            ! k point path ";
                fprintf(fid,"%s\n",head_string);
                fprintf(fid,"%d\n",abs(length(vasplibobj.kpoints_name-1)));
                for i=1:abs(length(vasplibobj.kpoints_name-1))
                    fprintf(fid,"%1s ",strrep(vasplibobj.kpoints_name(i),'GAMMA','G'));
                    fprintf(fid,"%9f ",num2str(vasplibobj.kpoints_frac(2*i-1,1)));
                    fprintf(fid,"%9f ",num2str(vasplibobj.kpoints_frac(2*i-1,2)));
                    fprintf(fid,"%9f ",num2str(vasplibobj.kpoints_frac(2*i-1,3)));
                    fprintf(fid,"%1s ",strrep(vasplibobj.kpoints_name(i+1),'GAMMA','G'));
                    fprintf(fid,"%9s ",num2str(vasplibobj.kpoints_frac(2*i,1)));
                    fprintf(fid,"%9s ",num2str(vasplibobj.kpoints_frac(2*i,2)));
                    fprintf(fid,"%9s \n",num2str(vasplibobj.kpoints_frac(2*i,3)));
                end
                fclose(fid);
            end
        end
        function fid = write_pj(vasplibobj,mode,filename)
            arguments
                vasplibobj ;
                mode = 'wt';
                filename = 'wt.in';
            end
            % function string
            if strcmp(mode,'wt')
                fid = fopen(filename,'a');
                %% wt.in projectors card gen
                maprule3= containers.Map({-1,0,1,2,4,5,6,8},{"","s", "pz px py","dz2 dxz dyz dx2-y2 dxy","s pz px py","s dz2 dxz dyz dx2-y2 dxy","pz px py dz2 dxz dyz dx2-y2 dxy","s pz px py dz2 dxz dyz dx2-y2 dxy"});
                maprule4= containers.Map({-1,0,1,2,3,4,5,6,7,8,9},{0,1,3,5,7,4,6,8,8,9,16});
                % function string
                head_string="PROJECTORS ";
                fprintf(fid,"%s\n",head_string);
                num_wan=0;
                for i=1:vasplibobj.basis_num
                    tempnum=maprule4(Projector_list(i));
                    num_wan = num_wan+tempnum;
                    fprintf(fid,"%d ",tempnum);
                end
                fprintf(fid,"\n");
                for i=1:vasplibobj.basis_num
                    if Projector_list(i) ~= -1
                        fprintf(fid,"%s ",'X');
                        % fprintf(fid,"%s ",vasplibobj.elementL(i)); % here we can
                        fprintf(fid,"%s ",maprule3(Projector_list(i)));
                        fprintf(fid,"\n");
                    end
                end
                fclose(fid);
            end
        end
    end
    %% Abstract
    %% tools
    % ---------------------   VASP I/O   ------------------------
    methods(Static)
        function [orbital_out,Rm_s_fin] = AddVacuumLayer(orbital_init,POSCAR_file,fin_dir_list,options)
            arguments
                orbital_init
                POSCAR_file = 'POSCAR';
                fin_dir_list = [1 1 0];
                options.fast = true;
                options.vacuum_length = 10;
            end
            orbital_out = orbital_init;
            fin_orb = orbital_out;
            Ns = [1 0 0;0 1 0;0 0 1];
            if ~options.fast
                %disp(fin_dir_list);
                [Rm_tmp,sites_tmp,Atom_name_tmp,Atom_num_tmp]=vasplib.POSCAR_read(POSCAR_file);
                % gen POSCAR
                H_hr.supercell(Ns,'POSCAR_super_fin',Rm_tmp,sites_tmp,Atom_name_tmp,Atom_num_tmp,fin_dir_list);
            else
                switch class(POSCAR_file)
                    case {'string','char'}
                        [Rm_tmp,~,~,~]=vasplib.POSCAR_read(POSCAR_file);
                    case 'double'
                        Rm_tmp = POSCAR_file;
                end
            end
            % rebuild fin_orb
            Rm_tmp = Ns*Rm_tmp;
            Rmlength1 = norm (Rm_tmp(1,:));
            Rmlength2 = norm (Rm_tmp(2,:));
            Rmlength3 = norm (Rm_tmp(3,:));
            Rm_s_fin_add = [options.vacuum_length*Rm_tmp(1,:)*fin_dir_list(1)/Rmlength1;...
                options.vacuum_length*Rm_tmp(2,:)*fin_dir_list(2)/Rmlength2;...
                options.vacuum_length*Rm_tmp(3,:)*fin_dir_list(3)/Rmlength3];
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
            orbital_out = fin_orb;% change coordinate along finite direction ; fractional
        end
        function POSCAR = POSCAR_cell_read(filename,formatSpec)
            startRow = 2;
            delimiter = {'\t',' '};
            fileID = fopen(filename,'r');
            dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'TextType', 'string', 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
            fclose(fileID);
            raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
            for col=1:length(dataArray)-1
                raw(1:length(dataArray{col}),col) = mat2cell(dataArray{col}, ones(length(dataArray{col}), 1));
            end
            rawStringColumns = string(raw(:,:));
            POSCAR=rawStringColumns;
            clearvars filename delimiter startRow formatSpec fileID dataArray ans rawStringColumns raw;
        end
        function [Rm,sites,Atom_name,Atom_num,elements]=POSCAR_read(filename,mode,options)
            %
            % get poscar information from vasp and others
            % * Label:
            %
            % Description of the Function:
            %
            arguments
                filename string ='POSCAR';
                mode char {mustBeMember(mode,{'symble','vasp','tbsk','tbsym','sym','list'})} = 'vasp';
                options.digits = 15;
            end
            %--------  init  --------
            try
                elements = readtable('elements.txt');
            catch
                elements = [];
            end
            %--------  chek  --------
            if ~exist(filename,'file')
                fprintf('No such file: %s\n',filename);
                error('file not exist!');
            end
            %--------  juge  --------
            switch mode
                case {'vasp','symble'}
                    site=struct('seq',[],'inseq',[],'rc1',[],'rc2',[],'rc3',[],'name',[],'nameseq',[]);
                    formatSpec = '%s%s%s%s%s%[^\n\r]';
                case {'tbsk','tbsym'}
                    site=struct('seq',[],'inseq',[],'rc1',[],'rc2',[],'rc3',[],'name',[],'nameseq',[],'element',[],'orb',[],'orb_sym',[],'spin',[]);
                    formatSpec = '%s%s%s%s%s%s%s%[^\n\r]';
                case 'sym'
                    site=struct('seq',[],'inseq',[],'rc1',[],'rc2',[],'rc3',[],'name',[],'nameseq',[],'element',[],'orb',[],'orb_sym',[],'spin',[],'Jz',[]);
                    formatSpec = '%s%s%s%s%s%s%s%s%[^\n\r]';
                case 'list'
                    site=struct('seq',[],'rc1',[],'rc2',[],'rc3',[],'Hue',[],'surf_level',[],'hing_level',[]);

            end
            %--------  init  --------
            digits(options.digits); %take care
            POSCAR = vasplib.POSCAR_cell_read(filename,formatSpec);
            % a_crystal_constance=str2double(char(POSCAR(1,1)));
            % Rm
            a1=[str2double(char(POSCAR(2,1))) str2double(char(POSCAR(2,2))) str2double(char(POSCAR(2,3)))];
            a2=[str2double(char(POSCAR(3,1))) str2double(char(POSCAR(3,2))) str2double(char(POSCAR(3,3)))];
            a3=[str2double(char(POSCAR(4,1))) str2double(char(POSCAR(4,2))) str2double(char(POSCAR(4,3)))];
            Rm=[a1;a2;a3];
            %
            %Rm = (Rm)/det(Rm);
            % read atom
            n_atom=length(POSCAR(5,:));
            for i=1:n_atom
                if POSCAR(5,i) ~= ""
                    Atom_name(i) = POSCAR(5,i);
                end
            end
            % need
            for i=1:length(Atom_name)
                Atom_num(i)=str2double(char(POSCAR(6,i)));
            end
            %site_num
            sites_num=sum(Atom_num);
            sites=repmat(site,[1 sites_num]);
            % coordinate pattern
            % Coordinates_pattern=POSCAR(7,1);
            % temp : frac support
            % struct
            sequence=1;
            n_atom=length(Atom_name);
            labelcut_list = vasplib.labelcut_list_gen(Atom_num);

            %first
            if n_atom >= 1
                for i=1:n_atom
                    inseq =1;
                    for j=labelcut_list(i,1):labelcut_list(i,2)
                        if j > size(POSCAR,1)
                            error('Your POSCAR 7th line mismatches the fractional coordinates of each orbital  ')
                        end
                        %id
                        sites(sequence).seq=sequence;
                        %inid
                        sites(sequence).inseq=inseq ;
                        %name
                        sites(sequence).nameseq=i;
                        sites(sequence).name=Atom_name(i)+num2str(sites(sequence).inseq);
                        sites(sequence).ion_num=Atom_num(i);
                        %sites(sequence).ion_type_num=i;
                        %incoordinate
                        sites(sequence).rc1=str2double(char(POSCAR(j,1)));
                        sites(sequence).rc2=str2double(char(POSCAR(j,2)));
                        sites(sequence).rc3=str2double(char(POSCAR(j,3)));
                        if strcmp(mode,'tbsk')
                            sites(sequence).element = string(POSCAR(j,4));
                            sites(sequence).orb = string(POSCAR(j,5));
                            sites(sequence).orb_sym = string(POSCAR(j,6));
                            if length(POSCAR(j,:)) > 6
                                sites(sequence).spin = str2double(POSCAR(j,7));
                            else
                                sites(sequence).spin = 1;
                            end
                        elseif strcmp(mode,'tbsym')
                            sites(sequence).element = string(POSCAR(j,4));
                            sites(sequence).orb = string(POSCAR(j,5));
                            sites(sequence).orb_sym = string(POSCAR(j,6));
                            sites(sequence).spin = str2double(POSCAR(j,7));
                        elseif strcmp(mode,'sym')
                            sites(sequence).element = string(POSCAR(j,4));
                            sites(sequence).orb = string(POSCAR(j,5));
                            sites(sequence).orb_sym = string(POSCAR(j,6));
                            sites(sequence).spin = str2double(POSCAR(j,7));
                            sites(sequence).J = str2double(POSCAR(j,8));
                        end
                        %
                        sequence=sequence+1;
                        inseq =inseq + 1;

                    end
                end
            end

            if length(POSCAR)>j && strcmp(mode,'mag')
                %disp('mag_mode');
                sequence=1;
                beginline=j+1;
                %first
                for j=beginline:beginline+Atom_num(1)-1
                    %id
                    sites(sequence).mag=double(POSCAR(j,1:3));
                    sequence=sequence+1;
                end
                %other
                if n_atom >= 2
                    for i=2:n_atom
                        beginline=beginline+Atom_num(i-1);
                        for j=beginline:beginline+Atom_num(i)-1
                            %id
                            sites(sequence).mag=double(POSCAR(j,1:3));
                            sequence=sequence+1;
                        end
                    end
                end
            end
        end
        function tmpsites = MakeUpSites(sites,spintype)
            tmpsites= sites(1);
            tmpsites(1) = [];
            % review
            target_sites = sites(1);
            for i = 2:length(sites)
                if strcmp(sites(i).element,"") && strcmp(sites(i).orb,"")
                    sites(i).element = target_sites.element;
                    sites(i).orb = target_sites.orb;
                    sites(i).orb_sym = target_sites.orb_sym;
                    % sites(i).spin = target_sites.spin;
                elseif ~strcmp(sites(i).element,"") && ~strcmp(sites(i).orb,"")
                    target_sites = sites(i);
                else
                    error('Your POSCAR file format is not right!');
                end
            end
            % expand
            for i = 1:length(sites)
                if strcmp(sites(i).orb_sym,"")
                    switch char(sites(i).orb)
                        case 's'
                            tmpsites_slib = sites(i);
                            tmpsites_slib(1).orb_sym = "s";tmpsites_slib(1).orb = "s";
                        case 'sp'
                            tmpsites_slib = repmat(sites(i),[1,4]);
                            tmpsites_slib(1).orb_sym = "s";tmpsites_slib(1).orb = "s";
                            tmpsites_slib(2).orb_sym = "y";tmpsites_slib(2).orb = "p";
                            tmpsites_slib(3).orb_sym = "z";tmpsites_slib(3).orb = "p";
                            tmpsites_slib(4).orb_sym = "x";tmpsites_slib(4).orb = "p";
                        case 'spz'
                            tmpsites_slib = repmat(sites(i),[1,2]);
                            tmpsites_slib(1).orb_sym = "s";tmpsites_slib(1).orb = "s";
                            tmpsites_slib(2).orb_sym = "z";tmpsites_slib(2).orb = "p";
                        case 'pxpy'
                            tmpsites_slib = repmat(sites(i),[1,2]);
                            tmpsites_slib(1).orb_sym = "x";tmpsites_slib(1).orb = "p";
                            tmpsites_slib(2).orb_sym = "y";tmpsites_slib(2).orb = "p";
                        case 'p'
                            tmpsites_slib = repmat(sites(i),[1,3]);
                            tmpsites_slib(1).orb_sym = "z";tmpsites_slib(1).orb = "p";
                            tmpsites_slib(2).orb_sym = "x";tmpsites_slib(2).orb = "p";
                            tmpsites_slib(3).orb_sym = "y";tmpsites_slib(3).orb = "p";
                        case 'sd'
                            tmpsites_slib = repmat(sites(i),[1,6]);
                            tmpsites_slib(1).orb_sym = "s";tmpsites_slib(1).orb = "s";
                            tmpsites_slib(2).orb_sym = "z^2";tmpsites_slib(2).orb = "d";
                            tmpsites_slib(3).orb_sym = "xz";tmpsites_slib(3).orb = "d";
                            tmpsites_slib(4).orb_sym = "yz";tmpsites_slib(4).orb = "d";
                            tmpsites_slib(5).orb_sym = "x^2-y^2";tmpsites_slib(5).orb = "d";
                            tmpsites_slib(6).orb_sym = "xy";tmpsites_slib(6).orb = "d";
                        case 'spd'
                            tmpsites_slib = repmat(sites(i),[1,9]);
                            tmpsites_slib(1).orb_sym = "s";tmpsites_slib(1).orb = "s";
                            tmpsites_slib(2).orb_sym = "z";tmpsites_slib(2).orb = "p";
                            tmpsites_slib(3).orb_sym = "x";tmpsites_slib(3).orb = "p";
                            tmpsites_slib(4).orb_sym = "y";tmpsites_slib(4).orb = "p";
                            tmpsites_slib(5).orb_sym = "z^2";tmpsites_slib(5).orb = "d";
                            tmpsites_slib(6).orb_sym = "xz";tmpsites_slib(6).orb = "d";
                            tmpsites_slib(7).orb_sym = "yz";tmpsites_slib(7).orb = "d";
                            tmpsites_slib(8).orb_sym = "x^2-y^2";tmpsites_slib(8).orb = "d";
                            tmpsites_slib(9).orb_sym = "xy";tmpsites_slib(9).orb = "d";
                        case 'pd'
                            tmpsites_slib = repmat(sites(i),[1,8]);
                            tmpsites_slib(1).orb_sym = "z";tmpsites_slib(1).orb = "p";
                            tmpsites_slib(2).orb_sym = "x";tmpsites_slib(2).orb = "p";
                            tmpsites_slib(3).orb_sym = "y";tmpsites_slib(3).orb = "p";
                            tmpsites_slib(4).orb_sym = "z^2";tmpsites_slib(4).orb = "d";
                            tmpsites_slib(5).orb_sym = "xz";tmpsites_slib(5).orb = "d";
                            tmpsites_slib(6).orb_sym = "yz";tmpsites_slib(6).orb = "d";
                            tmpsites_slib(7).orb_sym = "x^2-y^2";tmpsites_slib(7).orb = "d";
                            tmpsites_slib(8).orb_sym = "xy";tmpsites_slib(8).orb = "d";
                        case 'd'
                            tmpsites_slib = repmat(sites(i),[1,5]);
                            tmpsites_slib(1).orb_sym = "z^2";tmpsites_slib(1).orb = "d";
                            tmpsites_slib(2).orb_sym = "xz";tmpsites_slib(2).orb = "d";
                            tmpsites_slib(3).orb_sym = "yz";tmpsites_slib(3).orb = "d";
                            tmpsites_slib(4).orb_sym = "x^2-y^2";tmpsites_slib(4).orb = "d";
                            tmpsites_slib(5).orb_sym = "xy";tmpsites_slib(5).orb = "d";
                        case 'f'
                            %
                    end
                else
                    tmpsites_slib = sites(i);
                end
                tmpsites = [tmpsites,tmpsites_slib];
            end
            switch spintype
                case 'spinless'
                    return
                case 'wannier'
                    ntmpsites = length(tmpsites);
                    tmpsites1 = tmpsites;
                    if isnan(tmpsites1(1).spin)
                        [tmpsites1(1:ntmpsites).spin] = deal(0.5);
                    end
                    tmpsites2 = tmpsites1;
                    % [A(1:100).age]=deal(0) use deal !
                    [tmpsites2(1:ntmpsites).spin] = deal(-tmpsites2(1).spin);
                    tmpsites = repmat(tmpsites,[1 2]);
                    tmpsites(1:2:ntmpsites*2-1) = tmpsites1;
                    tmpsites(2:2:ntmpsites*2) = tmpsites2;
                case 'block'
                    ntmpsites = length(tmpsites);
                    tmpsites1 = tmpsites;
                    if isnan(tmpsites1(1).spin)
                        [tmpsites1(1:ntmpsites).spin] = deal(0.5);
                    end
                    tmpsites2 = tmpsites1;
                    [tmpsites2(1:ntmpsites).spin] = deal(-tmpsites2(1).spin);
                    tmpsites = [tmpsites1,tmpsites2];
            end
        end
        
        function [klist_cart,klist_frac,klist_l,kpoints_l,kpoints_frac] = kpathgen(kpoints,nodes,Gk,Gk_,options)
            arguments
                kpoints;
                nodes = 60;
                Gk = [];
                Gk_ = [];
                options.Dim = 3;
            end
            %-------- nargin --------
            if nargin > 3 && ~isequal(Gk,Gk_)
                mode = 'return mode';
            else
                mode = 'norm';
            end
            Dim = options.Dim;
            %--------  init  --------
            kn = size(kpoints,1)/2;
            if length(nodes) == 1
                klist_frac=zeros(kn*nodes,Dim);
                nodes = ones(kn,1)*nodes;
            else
                klist_frac=zeros(sum(nodes),Dim);
            end
            for i = 1:kn
                for d = 1:Dim
                    klisttemp{d} = linspace(kpoints(2*i-1,d),kpoints(2*i,d),nodes(i)).';
                end
                start_seq = sum(nodes(1:i-1))+1;
                end_seq = sum(nodes(1:i));
                klist_frac(start_seq :end_seq ,:)   = fold(@horzcat,klisttemp);
            end
            klist_cart = klist_frac*Gk;
            %klist_liner & kpoints_liner
            kpoints_l = zeros(kn+1,1);
            kpoints_l(1) = 0;
            for i = 1:kn
                lenghth = norm((kpoints(2*i,:)-kpoints(2*i-1,:))*Gk);
                temp = kpoints_l(i)+lenghth;
                kpoints_l(i+1) = temp;
            end
            klist_l = zeros(1,sum(nodes));
            for i = 1:kn
                klisttemp = linspace(kpoints_l(i),kpoints_l(i+1),nodes(i));
                start_seq = sum(nodes(1:i-1))+1;
                end_seq = sum(nodes(1:i));
                klist_l(1,start_seq :end_seq) = klisttemp;
            end
            if strcmp(mode,'return mode')
                klist_frac = klist_cart/Gk_;
            end
            kpoints_frac = kpoints;
        end
    end
    % ---------------------   semi infinite Green function   ------------------------
    methods(Static)
        function DOSCAR = DOSCAR_gen(GREENCAR,mode)
            % nargin
            if nargin <2
                mode = 'green';
            end
            if strcmp(mode,'green')
                if iscell(GREENCAR)
                    [Nsize1,Nsize2] = size(GREENCAR);
                    DOSCAR =zeros(Nsize1,Nsize2);
                    for i = 1:Nsize1
                        for j = 1:Nsize2
                            DOSCAR(i,j) = -trace(imag(GREENCAR{i,j}));
                        end
                    end
                else
                    disp('temp support cell format for GREENCAR');
                end
            end
        end
        function GREENCAR = GREENCAR_gen(w_list,eta,H00_H01_cell_list_1,H00_H01_cell_list_2,mode)
            if nargin < 4
                mode = 'Green_iter';
            end
            switch mode
                case 'Green_iter'
                    %H00_H01_cell_list = H00_H01_cell_list;
                    WAN_NUM = size(H00_H01_cell_list_1,1);
                    kn =size( H00_H01_cell_list_2,3);
                    Nsize1 = length(w_list);
                    if Nsize1 == 1
                        % arc
                    end
                    Nsize2 = kn;
                    GREENCAR1{Nsize1,Nsize2} = zeros(WAN_NUM) ;
                    GREENCAR2{Nsize1,Nsize2} = zeros(WAN_NUM) ;
                    GREENCAR3{Nsize1,Nsize2} = zeros(WAN_NUM) ;
                    count = 0;
                    pb = vasplib_tool_outer.CmdLineProgressBar('Green Solving ');
                    Nsize = Nsize1*Nsize2;
                    for i = 1:Nsize1
                        for j= 1:Nsize2
                            count = count +1;
                            pb.print(count,Nsize);
                            w =w_list(i);
                            H00 = H00_H01_cell_list_1(:,:,j);
                            H01 = H00_H01_cell_list_2(:,:,j);
                            [GREENCAR3{i,j},GREENCAR1{i,j},GREENCAR2{i,j}] = HR.GW_iter(H00,H01,w,eta);
                        end
                    end
                    pb.delete()
                    GREENCAR.bulk = GREENCAR1;
                    GREENCAR.surf_l = GREENCAR2;
                    GREENCAR.surf_r = GREENCAR3;
                case 'Tmatrix_iter'
                    %H00_H01_cell_list = H00_H01_cell_list;
                    WAN_NUM = size(H00_H01_cell_list_1,1);
                    kn =size( H00_H01_cell_list_2,3);
                    Nsize1 = length(w_list);
                    Nsize2 = kn;
                    GREENCAR1{Nsize1,Nsize2} = zeros(WAN_NUM) ;
                    GREENCAR2{Nsize1,Nsize2} = zeros(WAN_NUM) ;
                    GREENCAR3{Nsize1,Nsize2} = zeros(WAN_NUM) ;
                    for i = 1:Nsize1
                        for j= 1:Nsize2
                            w =w_list(i);
                            H00 = H00_H01_cell_list_1(:,:,j);
                            H01 = H00_H01_cell_list_2(:,:,j);
                            [GREENCAR3{i,j},GREENCAR1{i,j},GREENCAR2{i,j}] = HR.Tmatrix_iter(H00,H01,w,eta);
                            % [GREENCAR1{i,j}, GREENCAR2{i,j},GREENCAR3{i,j}]=   Tmatrix2Green00(Tmatrix,H00,H01,w,eta);
                        end
                    end
                    GREENCAR.bulk = GREENCAR1;
                    GREENCAR.surf_l = GREENCAR2;
                    GREENCAR.surf_r = GREENCAR3;
                    %                 case  't'
                    %                     if iscell(H00_H01_cell_list)
                    %                         disp('do nothing');
                    %                     elseif isstruct(H00_H01_cell_list)
                    %                         vasplibobj = H00_H01_cell_list;
                    %                         H00_H01_cell_list =H_cell_list_gen(vasplibobj);
                    %                     end
                    %                     Nsize1 = length(w_list);
                    %                     Nsize2 = length(H00_H01_cell_list);
                    %                     %w_list= linspace(w_range(1),w_range(2),w_number);
                    %                     WAN_NUM = length(H00_H01_cell_list{1});
                    %                     GREENCAR{Nsize1,Nsize2} = zeros(WAN_NUM) ;
                    %                     for i = 1:Nsize1
                    %                         for j= 1:Nsize2
                    %                             GREENCAR{i,j} =  GREENCAR_one_gen(H00_H01_cell_list{j},w_list(i),eta,'TB');
                    %                         end
                    %                     end
                    %                 case 'Tmatrix_eig'
                    %                     %H00_H01_cell_list = H00_H01_cell_list;
                    %                     WAN_NUM = length(H00_H01_cell_list{1,1});
                    %                     [~,kn] =size( H00_H01_cell_list);
                    %                     Nsize1 = length(w_list);
                    %                     Nsize2 = kn;
                    %                     GREENCAR{Nsize1,Nsize2} = zeros(WAN_NUM) ;
                    %                     for i = 1:Nsize1
                    %                         for j= 1:Nsize2
                    %                             w =w_list(i);
                    %                             H00 = H00_H01_cell_list{1,j};
                    %                             H01 = H00_H01_cell_list{2,j};
                    %                             Tmatrix = HR.Tmatrix_gen(H00,H01,w,0);
                    %                             GREENCAR{i,j} = Tmatrix2pband(Tmatrix,H00,H01,w);
                    %                         end
                    %                     end
                case 'Tmatrix_inv'
                    %H00_H01_cell_list = H00_H01_cell_list;
                    WAN_NUM = size(H00_H01_cell_list_1,1);
                    kn =size( H00_H01_cell_list_2,3);
                    Nsize1 = length(w_list);
                    Nsize2 = kn;
                    GREENCAR1{Nsize1,Nsize2} = zeros(WAN_NUM) ;
                    GREENCAR2{Nsize1,Nsize2} = zeros(WAN_NUM) ;
                    GREENCAR3{Nsize1,Nsize2} = zeros(WAN_NUM) ;
                    for i = 1:Nsize1
                        for j= 1:Nsize2
                            w =w_list(i);
                            H00 = H00_H01_cell_list_1(:,:,j);
                            H01 = H00_H01_cell_list_2(:,:,j);
                            Tmatrix = HR.Tmatrix_gen(H00,H01,w,eta);
                            [GREENCAR1{i,j}, GREENCAR2{i,j},GREENCAR3{i,j}]=   HR.Tmatrix2Green00(Tmatrix,H00,H01,w,eta);
                        end
                    end
                    GREENCAR.bulk = GREENCAR1;
                    GREENCAR.surf_l = GREENCAR2;
                    GREENCAR.surf_r = GREENCAR3;


            end
        end
        function [Gwl,Gwb,Gwr] = Tmatrix_iter(H00,H01,w,eta,mu_max,infinity_small)
            if nargin <6
                infinity_small= 1e-10;
            end
            %%
            if nargin <5
                mu_max = 100;
            end
            oumega =  (w+1i*eta)*eye(length(H00));
            oumega_00 = oumega-H00;
            t0 =(oumega_00)\H01';
            t0_bar =(oumega_00)\H01;
            T = t0;
            T_bar = t0_bar;
            ti = t0;
            ti_bar = t0_bar;
            ti_times = t0;
            ti_bar_times = t0_bar;
            I = eye(length(H00));
            for i = 1:mu_max
                M1 = ti*ti_bar;
                M2 = ti_bar*ti;
                M3 = M1 + M2;
                ti = (I-M3)\ti^2;
                ti_bar = (I-M3)\ti_bar^2;
                T =T + ti_bar_times*ti;
                T_bar =T_bar + ti_times*ti_bar;
                %difference = ti +ti_bar;
                if norm(M3,'fro') < infinity_small
                    %disp('reach_accuracy')
                    break;
                end

                ti_times = ti_times*ti;
                ti_bar_times = ti_bar_times*ti_bar;

            end

            Gwl = inv(oumega_00-H01*T);
            Gwb = inv((oumega_00- H01*T-H01'*T_bar));
            Gwr = inv((oumega_00-H01'*T_bar));
        end
        function [Gwl,Gwb,Gwr] = GW_iter(H00,H01,w,eta,mu_max,infinity_small)
            %%
            if nargin <6
                infinity_small= 1e-10;
            end
            %%
            if nargin <5
                mu_max = 100;
            end
            epsilon0 = H00;
            alpha0 = H01;
            beta0 = H01';

            %oumega_bk = E_w;
            oumega =  (w+1i*eta)*eye(length(H00));
            %save('Green_prepare.mat');
            epsilon_s0 = epsilon0;
            epsilon_s_bar0 = epsilon0;
            %[epsilon,epsilon_s,epsilon_s_bar,~,~] = GW_iter_coff(mu,oumega,alpha0,beta0,epsilon0);

            for i = 1:mu_max
                Ml = (oumega - epsilon0)\alpha0;
                Mr = (oumega - epsilon0)\beta0;
                M1 = alpha0 * (Ml);
                M2 = beta0  * (Mr);
                M3 = alpha0 * (Mr);
                M4 = beta0  * (Ml);

                difference = M3+ M4;
                if norm(difference,'fro') < infinity_small
                    %disp('reach_accuracy')
                    epsilon0 = epsilon0 +M3+ M4;
                    epsilon_s0 = epsilon_s0 + M3;
                    epsilon_s_bar0 = epsilon_s_bar0 + M4;
                    break;
                end
                %     disp(max( difference(:)));
                alpha0 = M1;
                beta0 = M2;
                epsilon0 = epsilon0 +M3+ M4;
                epsilon_s0 = epsilon_s0 + M3;
                epsilon_s_bar0 = epsilon_s_bar0 + M4;

                %     [epsilon,epsilon_s,epsilon_s_bar,alpha,beta,difference] = GW_iter_coff_once(oumega,epsilon0,epsilon_s0,epsilon_s_bar0,alpha0,beta0);
                %     alpha0 = alpha;
                %     beta0 = beta;
                %     epsilon0 = epsilon;
                %     epsilon_s0 = epsilon_s;
                %     epsilon_s_bar0 = epsilon_s_bar;
                %disp(i);
                %disp(max(difference(:)));

            end
            % epsilon = epsilon0;
            % epsilon_s = epsilon_s0 ;
            % epsilon_s_bar = epsilon_s_bar0;
            Gwl = inv((oumega-epsilon_s0));
            Gwb = inv((oumega-epsilon0));
            Gwr = inv((oumega-epsilon_s_bar0));
        end
        function Tmatrix = Tmatrix_gen(H00,H01,w,eta)
            Dimi = length(H00);
            wc = w + 1i*eta;
            W = (wc*eye(Dimi)-H00);

            if rcond(H01) < 1e-10

                %temp_diag = diag(H01)==0;
                %disp('work');

                H01 = H01 +  eta*eye(Dimi);
                %isnan(H01)

                %H01 = (10^-nfin)*eye(Dimi);
                %disp('nor work');
                %Tmatrix = [pinv(H01)*W,-pinv(H01)*(H01');eye(Dimi) zeros(Dimi)];
                Tmatrix = [H01\W,-H01\(H01');eye(Dimi) zeros(Dimi)];
                if rcond(Tmatrix) < 1e-10
                    Tmatrix = [lsqminnorm(H01,W),-lsqminnorm(H01,(H01'));eye(Dimi) zeros(Dimi)];
                    return;
                end
            else
                Tmatrix = [H01\W,-H01\(H01');eye(Dimi) zeros(Dimi)];
            end


        end
        function [Green_00,Green_00s1,Green_00s2]= Tmatrix2Green00(Tmatrix,H00,H01,w,eta,n)
            if nargin < 6
                n=2;
            end
            Dimi = length(H00);
            wc = w + 1i*eta;
            %disp(Tmatrix);
            Tmatrix = Tmatrix^n ;
            eyemat = eye(length(Tmatrix));
            %             if abs(det(Tmatrix)) <1e-6
            %                 Tmatrix = Tmatrix + (1e-6)*eye(Dimi*2)*1i;
            %             end
            if rcond(Tmatrix) < 1e-10
                %[A,U] = eig(eyemat,Tmatrix);
                Tmatrix = Tmatrix + (1e-6)*eye(Dimi*2);
                [A,U] = eig(Tmatrix);
            else
                [A,U] = eig(Tmatrix);
            end
            U_abs =abs(U);
            [Asort,Usort] =park.sorteig(U_abs,A);
            Lambda= diag(Usort);
            S = Asort(:,Lambda<1);
            SS = Asort(:,Lambda>1);
            S2 = S(1:Dimi,:);
            S1 = S(Dimi+1:2*Dimi,:);
            S3 = SS(Dimi+1:2*Dimi,:);
            S4 = SS(1:Dimi,:);
            phi1 = H01*S2/S1;
            phi2 = H01'*S3/S4;
            Green_00 = inv(wc*eye(Dimi)-H00-phi1-phi2);
            Green_00s2 = inv(wc*eye(Dimi)-H00-phi1);
            Green_00s1 = inv(wc*eye(Dimi)-H00-phi2);
            %     Green00.Green_00 =Green_00;
            %     Green00.Green_00s1 =Green_00s1;
            %     Green00.Green_00s2 =Green_00s2;
        end
    end
    % ---------------------   site  Green function   ------------------------
    methods(Static)
        function [Hue,surf_level,hing_level] = orbone2hsv(orb_one,discrimination,center,orientation)
            %--------  nargin  --------
            if nargin < 2
                discrimination = 0.1;
            end
            if nargin <3
                center  = [0.5, 0.5,0.5];
            end
            if nargin <4
                orientation = 3;
            end

            orb_init = orb_one - center;
            %--------  init  --------
            switch orientation
                case {1,-1}
                    x = orb_init(2);
                    y = orb_init(3);
                case {2,-2}
                    x = orb_init(3);
                    y = orb_init(1);
                case {3,-3}
                    x = orb_init(1);
                    y = orb_init(2);
            end
            r = norm([x y]);
            z = x+1i*y;
            Hue = (angle(z)+pi)/(2*pi);
            G_hinge = ((abs(x)+1i*discrimination-0.5)*(abs(y)+1i*discrimination-0.5))^-1;
            %G_hinge = (r+1i*discrimination-0.5*1.414)^-1;
            G_surf = (r+1i*discrimination-0.5)^-1;
            if orientation > 0
                hing_level = -imag(G_hinge);
                surf_level = -imag(G_surf);
            else
                hing_level = -imag(G_hinge)*sign(x*y);
                surf_level = -imag(G_surf)*sign(x*y);
            end

        end
        function HSVCAR = HSVCAR_gen(orb_list,mode,discrimination,center,orientation)
            %--------  nargin  --------
            if nargin < 3
                discrimination = 0.1;
            end
            if nargin < 4
                center  = [0.5, 0.5,0.5];
            end

            if nargin < 5
                orientation = 3;
            end
            if nargin < 2
                mode = 'hinge';
            end
            [norb,~] = size(orb_list);
            HSVCAR  = zeros(norb,1);
            switch mode
                case 'hinge'
                    for i = 1:norb
                        orb_one = orb_list(i,:);
                        [~,~,HSVCAR(i)] = vasplib.orbone2hsv(orb_one,discrimination,center,orientation);
                    end
                case 'surf'
                    for i = 1:norb
                        orb_one = orb_list(i,:);
                        [~,HSVCAR(i),~] = vasplib.orbone2hsv(orb_one,discrimination,center,orientation);
                    end
                case 'select'
                    nCenter = size(center,1);
                    for i = 1:nCenter/2
                        %((abs(x)+1i*discrimination-0.5)*(abs(y)+1i*discrimination-0.5))^-1;
                        HSVCAR = HSVCAR + ...
                            imag(((abs(orb_list(:,1) -center(i,1))+1i*discrimination)...
                            .*(abs(orb_list(:,2) -center(i,2))+1i*discrimination)...
                            .*(abs(orb_list(:,3) -center(i,3))+1i*discrimination))...
                            .^-1);
                    end
                    for i = nCenter/2+1 :nCenter
                        %((abs(x)+1i*discrimination-0.5)*(abs(y)+1i*discrimination-0.5))^-1;
                        HSVCAR = HSVCAR - ...
                            imag(((abs(orb_list(:,1) -center(i,1))+1i*discrimination)...
                            .*(abs(orb_list(:,2) -center(i,2))+1i*discrimination)...
                            .*(abs(orb_list(:,3) -center(i,3))+1i*discrimination))...
                            .^-1);
                    end
                    HSVCAR = HSVCAR/nCenter;
                case 'select-points'
                    nCenter = size(center,1);
                    for i = 1:nCenter
                        %((abs(x)+1i*discrimination-0.5)*(abs(y)+1i*discrimination-0.5))^-1;
                        HSVCAR = HSVCAR + ...
                            imag(((abs(orb_list(:,1) -center(i,1))+1i*discrimination)...
                            .*(abs(orb_list(:,2) -center(i,2))+1i*discrimination)...
                            .*(abs(orb_list(:,3) -center(i,3))+1i*discrimination))...
                            .^-1);
                    end
                    HSVCAR = HSVCAR/nCenter;
                case 'orient'
                    for i = 1:norb
                        orb_one = orb_list(i,:);
                        [HSVCAR(i),~,~] = vasplib.orbone2hsv(orb_one,discrimination,center,orientation);
                    end
                case 'slab'
                    G_slab = abs((abs(orb_list(:,orientation)-center(orientation))-0.5)) < discrimination;
                    HSVCAR = ~(G_slab);
            end
            % 0-1
            WAN_NUM = norb;
            HSVCAR = (-normalize(HSVCAR,'range')+1)/2;
            HSVCAR(:,2:3) = ones(WAN_NUM,2);
        end
        function OberserveValue = Observecar_gen(WAVECAR,Oper)
            OberserveValue = zeros(size(WAVECAR,2),1);
            for i = 1:size(WAVECAR,2)
                OberserveValue(i) = WAVECAR(:,i)'*Oper*WAVECAR(:,i);
            end
        end
        function [COLORCAR,WEIGHTCAR] = COLORCAR_gen(WAVECAR,HSVCAR,signlist)
            if nargin < 3
                signlist = 1;
            end
            if all(signlist==1)
                SIGN_mode = false;
            else
                SIGN_mode = true;
            end
            [~,norb,kn] = size(WAVECAR);
            temp.rgb = [0 0 0];
            COLORCAR = repmat(temp,norb,kn);
            WEIGHTCAR = zeros(norb,kn);
            if SIGN_mode
                for ki = 1:kn
                    WAVECAR_one = WAVECAR(:,:,ki);
                    for orbi = 1:norb
                        WAVEFUNC = WAVECAR_one(:,orbi);
                        [WEIGHTCAR(orbi,ki),~]  = vasplib.COLOR_one_gen(WAVEFUNC,HSVCAR,signlist);
                    end
                end
            else
                for ki = 1:kn
                    WAVECAR_one = WAVECAR(:,:,ki);
                    for orbi = 1:norb
                        WAVEFUNC = WAVECAR_one(:,orbi);
                        [rgb,SIGN_one]  = vasplib.COLOR_one_gen(WAVEFUNC,HSVCAR);
                        COLORCAR(orbi,ki).rgb = rgb;
                        hsv_temp = rgb2hsv(rgb);
                        WEIGHTCAR(orbi,ki) = SIGN_one*(hsv_temp(1)*100+hsv_temp(3)*10+hsv_temp(2)*1);
                    end
                end
            end
        end
        function [COLOR_one,SIGN_one] = COLOR_one_gen(WF,HSVCAR,signlist)
            if nargin <3
                SIGN_one = 1;
                signlist = 1;
            else

            end

            n = abs(WF.*conj(WF));
            %[WAN_NUM,~] = size(HSVCAR);
            %HSVCAR(:,2:3) = ones(WAN_NUM,2);
            if nargin == 3
                COLOR_one = sum(signlist.*HSVCAR(:,1).*n,1);
                %COLOR_one = normalize(COLOR_one,'range',[-1,1]);
            else
                RGBCAR =  hsv2rgb(HSVCAR);
                COLOR_one = sum(RGBCAR.*n,1);
                COLOR_one = normalize(COLOR_one,'range');
            end
            if nargin == 3
                SIGN_one = sign(sum(signlist.*n,1));
            end

            if all(signlist==1)
                SIGN_one = 1;
                return;
            end

        end     
    end
    % ---------------------   handy   ------------------------
    methods(Static)
        function str = relaceSpaceByComma(str)
            str = strrep(str,' ',',');
        end
        function str = mat2str_python(mat)
            str = mat2str(mat);
            str = vasplib.relaceSpaceByComma(str);
        end
        function SIGN_one = SIGN_one_gen(WF,signlist)
            n = abs(WF.*conj(WF));
            %[WAN_NUM,~] = size(HSVCAR);
            %HSVCAR(:,2:3) = ones(WAN_NUM,2);

            SIGN_one = sign(sum(signlist.*n,1));
        end
        function Symble = SymbolicVarible(seeds,superscript,subscript,level)
            Superscript = "";
            Subscript = "";
            for i = 1:length(superscript)
                if superscript(i) < 0
                    Superscript = Superscript+"__"+num2str(round(abs(superscript(i))))+"_bar";
                else
                    Superscript = Superscript+"__"+num2str(round(abs(superscript(i))));
                end
            end
            for i = 1:length(subscript)
                if subscript(i) < 0
                    Subscript = Subscript+"_"+num2str(round(abs(subscript(i))))+"_bar";
                else
                    Subscript = Subscript+"_"+num2str(round(abs(subscript(i))));
                end
            end
            if nargin == 4
                Level = "_"+num2str(level)+"_ubar";
                Symble = sym(seeds+Subscript+Superscript+Level,'real');
            else
                Symble = sym(seeds+Subscript+Superscript,'real');
            end
        end
    end
    %% plot-tools
    methods
        function varargout = bandplot(vasplibobj,varargin)
            % -------------- plot ------------------
            if  vasplibobj.coe
                vasplibobj = vasplibobj.Subsall();
            end
            if  nargin-1 ==1 && isvector(varargin{1})
                EIGENCAR = vasplibobj.EIGENCAR_gen();
                Ecut = varargin{1};
                varargin = varargin(2:end);
                bandplot(EIGENCAR,Ecut,vasplibobj.klist_l,vasplibobj.kpoints_l,vasplibobj.kpoints_name,varargin{:});
            elseif nargin ==1
                EIGENCAR = vasplibobj.EIGENCAR_gen();
                Ecut=  [-3,3];
                bandplot(EIGENCAR,Ecut,vasplibobj.klist_l,vasplibobj.kpoints_l,vasplibobj.kpoints_name);
            elseif ismatrix(varargin{1}) && ~isnumeric(varargin{2}) && ~isvector(varargin{2})
                EIGENCAR = varargin{1};
                varargin = varargin(2:end);
                Ecut=  [-3,3];
                [varargout{:}]  = bandplot(EIGENCAR,Ecut,vasplibobj.klist_l,vasplibobj.kpoints_l,vasplibobj.kpoints_name,varargin{:});
            elseif  ismatrix(varargin{1}) && isnumeric(varargin{1}) && length(varargin{1}) ==2 && isvector(varargin{1})
                EIGENCAR = varargin{1};Ecut = varargin{2};
                varargin = varargin(3:end);
                varargout{:} = bandplot(EIGENCAR,Ecut,vasplibobj.klist_l,vasplibobj.kpoints_l,vasplibobj.kpoints_name,varargin{:});
            else
                [varargout{:}] = bandplot(varargin{:});
            end
        end
    end
    methods(Static)
        % https://ww2.mathworks.cn/help/matlab/ref/varargout.html
        function varargout = WilsonLoopPlot(varargin)
            [varargout{1:nargout}] = WilsonLoopPlot(varargin{:});
        end
        function varargout = BZplot(varargin)
            [varargout{1:nargout}] = BZplot(varargin{:});
        end
        function varargout= pbandplot(varargin)
            [varargout{1:nargout}] = pbandplot(varargin{:});
        end
        function varargout = bandcompare(varargin)
            [varargout{1:nargout}] = bandcompare(varargin{:});
        end
        function varargout = waveplot(varargin)
            [varargout{1:nargout}] = waveplot(varargin{:});
        end
        function varargout = bandplot_3d(varargin)
            [varargout{1:nargout}] = bandplot3d(varargin{:}); % Name Truncation
        end
        function varargout = heatplot(varargin)
            [varargout{1:nargout}] = heatplot(varargin{:});
        end
        function varargout = heatplot3D(varargin)
            [varargout{1:nargout}] = heatplot3D(varargin{:});
        end
        function varargout = BCplot2D(varargin)
            [varargout{1:nargout}] = BCplot2D(varargin{:});
        end
        function varargout = ShowSurfIntegral(varargin)
            [varargout{1:nargout}] = ShowSurfIntegral(varargin{:});
        end
    end
    %% Gemometric Phase
    methods
        function Chern_number = Chern(vasplibobj,options,options_Chern)
            arguments
                vasplibobj ;
                options.BAND_index = [];
                options.knum1   = 51;
                options.knum2   = 51;
                options.kstart  = [-0.5,-0.5,0];
                options.kdir1   = [1,0,0];
                options.kdir2   = [0,1,0];
                options.cartesian = false;
                options.dir_seq = [1,2,3];
                options.dir_start = 'k_z';
                options.fig = false;
                options.plot = false;
                options.Oper = [];
                options.subband = [];
                options.ProjectionMethod = 'sign';
                options.ProjectionStruct = struct('field','imag');
                options_Chern.Accuracy = 1e-6;
            end
            optionsCell = namedargs2cell(options);
            [BCCAR,~,~,~] = BC_2D(vasplibobj,optionsCell{:});
            Chern_number = roundn(sum(BCCAR,'all')/(2*pi),round(log10(options_Chern.Accuracy)));
        end
        function [BCCAR,Grid,BC_WAVECAR,klist_r_plot] = BC_2D(vasplibobj,optionsK,options,optionsPlot)
            arguments
                vasplibobj ;
                optionsK.knum1   = 51;
                optionsK.knum2   = 51;
                optionsK.kstart  = [-0.5,-0.5,0];
                optionsK.kdir1   = [1,0,0];
                optionsK.kdir2   = [0,1,0];
                optionsK.cartesian = false;
                optionsK.dir_seq = [1,2,3];
                optionsK.dir_start = 'kcar';
                options.BAND_index = [];
                options.fig = false;
                options.plot = false;
                options.Oper = [];
                options.subband = [];
                options.sum = true;
                options.ProjectionMethod {mustBeMember(options.ProjectionMethod,{'sign','dynamicsign'})}= 'sign';
                options.ProjectionStruct = struct('field','imag');
                optionsPlot.oneshot = true;
                optionsPlot.view = [0,90];
            end
            % --- nargin
            switch class(vasplibobj)
                case {'Htrig','HK'}
                    %
                case 'HR'
                    vasplibobj.Basis_num = vasplibobj.WAN_NUM;
            end
            optionsKcell = namedargs2cell(optionsK);
            [klist_cart,klist_frac,klist_r_plot,sizemesh,Gk_,Grid] = vasplib.kmesh2D(vasplibobj.Rm,optionsKcell{:},'full',true);
            % individial chern number 
            if isempty(options.Oper)
                project = false;
                set_divide = 2;
                AOperU = [];
            else
                project = true;
                set_divide = 4;
                OperObj = options.Oper; % only support one Oper now
                switch options.ProjectionMethod
                    case 'sign'
                        OperU = roundn(OperObj.U,-6);
                        [AOperU,UOperU] = eig(OperU);
                        if strcmp(options.ProjectionStruct.field,'imag')
                            [AOperU,Ucheck] = park.sorteig(imag(UOperU),AOperU);
                        elseif strcmp(options.ProjectionStruct.field,'real')
                            [AOperU,Ucheck] = park.sorteig(real(UOperU),AOperU);
                        end
                    case 'dynamicsign'
                        OperU = OperObj.U;
                        switch class(vasplibobj)
                            case {'Htrig','HK'}
                                klist = klist_cart;
                            case 'HR'
                                klist = klist_frac;
                        end
                        AOperU = zeros(vasplibobj.Basis_num,vasplibobj.Basis_num,size(klist,1));
                        for i = 1:size(klist,1)
                            [AOperUtmp,UOperU] = eig(OperU(klist(i,1),klist(i,2),klist(i,3)));
                            if strcmp(options.ProjectionStruct.field,'imag')
                                [AOperUtmp,Ucheck] = park.sorteig(imag(UOperU),AOperUtmp);
                            elseif strcmp(options.ProjectionStruct.field,'real')
                                [AOperUtmp,Ucheck] = park.sorteig(real(UOperU),AOperUtmp);
                            end
                            AOperU(:,:,i) = AOperUtmp;
                        end
                    otherwise

                end
            end
            if isempty(options.BAND_index)
                switch class(vasplibobj)
                    case {'Htrig','HK','HR'}
                        BAND_index = 1:(vasplibobj.Basis_num/set_divide);
                    case ''
                        BAND_index = 1:(vasplibobj.Basis_num/set_divide);
                end
            else
                BAND_index = options.BAND_index;
            end
            if project
                if isempty(options.subband)
                    subband = 1:vasplibobj.Basis_num/2;
                else
                    subband = options.subband;
                end
            else
                subband = options.subband;
            end
            %
            switch class(vasplibobj)
                case {'Htrig','HK'}
                    [~,BC_WAVECAR] = vasplibobj.EIGENCAR_gen(...
                        'klist',klist_cart);
                case 'HR'
                    [~,BC_WAVECAR] = vasplibobj.EIGENCAR_gen(...
                        'klist',klist_frac,...
                        'convention','I','printmode',false, ...
                        'Umat',AOperU,'subband',subband);
            end
            BC_WAVECAR = BC_WAVECAR(:,BAND_index,:);
            %reshape(BC_WAVECAR(:,BAND_index,:),[vasplibobj.Basis_num options.knum1*options.knum2]);
            if options.sum
                BCCAR = sum(vasplib.BerryCuvature_2D( BC_WAVECAR ,sizemesh),3);
            else
                BCCAR = vasplib.BerryCuvature_2D( BC_WAVECAR ,sizemesh,'sum',false);
            end
            % remove the edge data

            if options.plot
                knum1 = sizemesh(1);
                knum2 = sizemesh(2);
                dk_1 = (optionsK.kdir1)/knum1*Gk_;
                dk_2 = (optionsK.kdir2)/knum2*Gk_;
                optionsPlotcell = namedargs2cell(optionsPlot);
                if options.sum
                    dSumL = BCCAR(:)/(2*pi);
                else
                    dSumL = sum(BCCAR,3)/(2*pi);
                    dSumL = dSumL(:);
                end
                [~,~] = vasplib.ShowSurfIntegral(vasplibobj.Gk,klist_r_plot,dk_1,dk_2,dSumL,optionsPlotcell{:});
            end
        end
        function [BFCAR,BF_WAVECAR,klist_l,WAVELOOPCAR] = WilsonLoop(vasplibobj,optionsK,options)
            arguments
                vasplibobj ;
                optionsK.knum_int    = 31;
                optionsK.knum_evol   = 51;
                optionsK.kstart      = [0,0,0];
                optionsK.kintegral   = [0,1,0];
                optionsK.kevolution  = [1,0,0];
                optionsK.cartesian = false;
                optionsK.dir_seq = [1,2,3];
                optionsK.dir_start = 'kcar';
                options.ax = handle([]);
                options.plot = false;
                options.V = [];
                options.Accuracy = 1E-12;
                options.LWAVE = false;
                options.BAND_index = [];
                options.script {mustBeMember(options.script,...
                    {...
                    'nu_x(k_y)','nu_y(k_x)',....
                    'nu_y(k_z)','nu_z(k_y)',...
                    'nu_z(k_x)','nu_x(k_z)',...
                    'nu_1(k_2)','nu_2(k_1)',...
                    'nu_2(k_3)','nu_3(k_2)',...
                    'nu_3(k_1)','nu_1(k_3)',...
                    ''})}= '';
            end
            % --- nargin
            switch class(vasplibobj)
                case {'Htrig','HK'}
                    NBAND = vasplibobj.Basis_num;
                case 'HR'
                    NBAND = vasplibobj.WAN_NUM;
            end
            if isempty(options.BAND_index)
                BAND_index = 1:(NBAND/2);
            else
                BAND_index = options.BAND_index;
            end
            % use script
            switch options.script
                case {'nu_x(k_y)','nu_1(k_2)'}
                    optionsK.kintegral   = [1,0,0];
                    optionsK.kevolution  = [0,1,0];
                    optionsK.kstart      = [0,-0.5,0];
                case {'nu_y(k_z)','nu_2(k_3)'}
                    optionsK.kintegral   = [0,1,0];
                    optionsK.kevolution  = [0,0,1];
                    optionsK.kstart      = [0,0,-0.5];
                case {'nu_z(k_x)','nu_3(k_1)'}
                    optionsK.kintegral   = [0,0,1];
                    optionsK.kevolution  = [1,0,0];
                    optionsK.kstart      = [-0.5,0,0];
                case {'nu_y(k_x)','nu_2(k_1)'}
                    optionsK.kintegral   = [0,1,0];
                    optionsK.kevolution  = [1,0,0];
                    optionsK.kstart      = [-0.5,0,0];
                case {'nu_z(k_y)','nu_3(k_2)'}
                    optionsK.kintegral   = [0,0,1];
                    optionsK.kevolution  = [0,1,0];
                    optionsK.kstart      = [0,-0.5,0];
                case {'nu_x(k_z)','nu_1(k_3)'}
                    optionsK.kintegral   = [1,0,0];
                    optionsK.kevolution  = [0,0,1];
                    optionsK.kstart      = [0,0,-0.5];
                case {''}
            end
            %
            if contains(string(options.script),["x","y","z"])
                optionsK.cartesian = true;
                optionsK.dir_seq = [1,2,3];
                optionsK.dir_start = 'kcar';
            end
            optionsKcell = namedargs2cell(optionsK);
            % Gen kloop 2D
            [kloop1_frac,kloop2_frac,kloop1_cart,kloop2_cart,klist_l,kstart_frac,kstart_cart] = vasplib.kloop2D(vasplibobj.Rm,optionsKcell{:});
            % prepare
            BFCAR = zeros(length(BAND_index),optionsK.knum_evol);
            BF_WAVECAR = zeros(length(BAND_index),length(BAND_index),optionsK.knum_evol);
            if options.LWAVE
                if isempty(options.V )
                    V = diag((1:NBAND) *sqrt(options.Accuracy));
                else
                    V = options.V ;
                end
                %EIGENCAR = zeros(length(BAND_index),optionsK.knum_int,optionsK.knum_evol);
                WAVELOOPCAR = zeros(NBAND,length(BAND_index),optionsK.knum_int,optionsK.knum_evol);
            end
            switch class(vasplibobj)
                case {'Htrig','HK'}
                    if isempty(vasplibobj.orbL)
                        vasplibobj.orbL = zeros(vasplibobj.Basis_num,3);
                    end
                    for i = 1:optionsK.knum_evol
                        klist_tmp = kloop1_cart(i,:)+kloop2_cart+kstart_cart;
                        [E,WAVECAR_loop] = vasplibobj.EIGENCAR_gen(...
                            'klist',klist_tmp,'printmode',false);
                        WAVECAR_loop_tmp = WAVECAR_loop(:,BAND_index,:);
                        % The last bloch state is the same as the first up to a phase factor
                        WAVECAR_loop_tmp(:,:,end) = WAVECAR_loop_tmp(:,:,1).* exp(-1i*(vasplibobj.orbL*vasplibobj.Rm*(klist_tmp(end,:)-klist_tmp(1,:)).'));
                        if options.LWAVE
                            WAVECAR_loop_tmp = vasplib.cleanWAVECAR(WAVECAR_loop_tmp,E(BAND_index,:),V,options.Accuracy);
                            WAVELOOPCAR(:,:,:,i) = WAVECAR_loop_tmp;
                        end
                        [BFCAR(:,i),BF_WAVECAR(:,:,i)] = vasplib.wancenter_1D(WAVECAR_loop_tmp);
                    end
                case 'HR'
                    for i = 1:optionsK.knum_evol
                        klist_tmp = kloop1_frac(i,:)+kloop2_frac+kstart_frac;
                        [E,WAVECAR_loop] = vasplibobj.EIGENCAR_gen(...
                            'klist',klist_tmp,...
                            'convention','II','printmode',false);
                        % If we use convention II, each wavefactor should
                        % give back the factor
                        % C^{nk}_j = C^{nk}_j_tilde * e^{-ik·tj}.
                        WAVECAR_loop_tmp = WAVECAR_loop(:,BAND_index,:);
                        % normalize phases to get u instead of phi
                        for j =1:size(WAVECAR_loop_tmp,3)
                            WAVECAR_loop_tmp(:,:,j) = WAVECAR_loop_tmp(:,:,j).* exp(-2i*pi*(vasplibobj.orbL*klist_tmp(j,:).'));
                        end
                        % The last bloch state is the same as the first up to a phase factor
                        WAVECAR_loop_tmp(:,:,end) = WAVECAR_loop_tmp(:,:,1).* exp(-2i*pi*(vasplibobj.orbL*(klist_tmp(end,:)-klist_tmp(1,:)).'));
                        if options.LWAVE
                            WAVECAR_loop_tmp = vasplib.cleanWAVECAR(WAVECAR_loop_tmp,E(BAND_index,:),V,options.Accuracy);
                            WAVELOOPCAR(:,:,:,i) = WAVECAR_loop_tmp;
                        end
                        [BFCAR(:,i),BF_WAVECAR(:,:,i)] = vasplib.wancenter_1D(WAVECAR_loop_tmp);
                    end
            end
        end
        function [BFCAR,WEIGHTCAR,klist_l] = ProjectedWilsonLoop(vasplibobj,options)
            arguments
                vasplibobj ;
                options.BAND_index = [];
                options.knum_int    = 31;
                options.knum_evol   = 51;
                options.kstart      = [0,0,0];
                options.kintegral   = [0,1,0];
                options.kevolution  = [1,0,0];
                options.cartesian = false;
                options.dir_seq = [1,2,3];
                options.dir_start = 'k_z';
                options.fig = false;
                options.plot = false;
                options.Oper = [];
                options.ProjectionMethod = 'sign';
                options.ProjectionStruct = struct('field','imag');
            end
            % --- nargin
            if isempty(options.Oper)
                project = false;
                set_divide = 2;
            else
                project = true;
                set_divide = 4;
                OperObj = options.Oper; % only support one Oper now
                OperU = roundn(OperObj.U,-6);
                [AOperU,UOperU] = eig(OperU);
                switch options.ProjectionMethod
                    case 'sign'
                        if strcmp(options.ProjectionStruct.field,'imag')
                            [AOperU,Ucheck] = park.sorteig(imag(UOperU),AOperU);
                            %ProjectionFunction = @(WAVECAR) sign(imag(OperObj.EIGEN(WAVECAR)));
                        elseif strcmp(options.ProjectionStruct.field,'real')
                            [AOperU,Ucheck] = park.sorteig(real(UOperU),AOperU);
                            %ProjectionFunction = @(WAVECAR) sign(real(OperObj.EIGEN(WAVECAR)));
                        end
                    otherwise

                end
            end
            if isempty(options.BAND_index)
                switch class(vasplibobj)
                    case {'Htrig','HK'}
                        
                        BAND_index = 1:vasplibobj.Basis_num/set_divide;
                    case 'HR'
                        BAND_index = 1:vasplibobj.WAN_NUM/set_divide;
                end
            else
                BAND_index = options.BAND_index;
            end
            %
            if options.cartesian
                if options.plot
                    [fig,ax] = vasplib.BZplot(vasplibobj.Rm,'color','r','alpha',0.1);
                    view(ax,1.465409092282576e+02,16.15644804763928);
                end
                Gk_ = vasplib.CartisianMat(vasplibobj.Gk,options.dir_seq,options.dir_start);
                if options.plot
                    [fig,ax] = vasplib.BZplot(Gk_,'Gk',true,'color','y','alpha',0.3,'ax',ax,'fig',fig);
                    view(ax,1.465409092282576e+02,16.15644804763928);
                    disp(Gk_);
                end
                kstart_s  = options.kstart * Gk_ /vasplibobj.Gk;
            else
                Gk_ = vasplibobj.Gk;
                kstart_s = options.kstart;
            end
            %
            [klist_r_1,klist_s_1,~,~] =...
                vasplib.kpathgen([[0,0,0];options.kevolution],options.knum_evol,Gk_,vasplibobj.Gk);
            klist_l = zeros(size(klist_s_1,1),1);
            %             klist_l(1) = sum(sign(klist_s_1(1,:)))*norm(klist_s_1(1,:)*(eye(3)*2*pi));
            klist_s_1_ = klist_s_1 ;
            normklist_l = norm(options.kevolution)/norm(klist_s_1(end,:));
            for i = 1:size(klist_s_1_,1)
                klist_l(i) = norm(klist_s_1_(i,:)*(eye(3)*2*pi))*normklist_l;
            end
            klist_l = klist_l + sum(sign(kstart_s))* norm(kstart_s*(eye(3)*2*pi))*normklist_l;
            [klist_r_2,klist_s_2,~,~] =...
                vasplib.kpathgen([[0,0,0];options.kintegral],options.knum_int,Gk_,vasplibobj.Gk);
            BFCAR = zeros(length(BAND_index)*2,options.knum_evol);
            WEIGHTCAR = BFCAR;
            BF_WAVECAR = zeros(length(BAND_index),length(BAND_index),options.knum_evol);
            kstart_r = options.kstart*Gk_;
            switch class(vasplibobj)
                case {'Htrig','HK'}
                    for i = 1:options.knum_evol
                        klist_tmp = klist_r_1(i,:)+klist_r_2+kstart_r;
                        [~,WAVECAR_loop1] = vasplibobj.EIGENCAR_gen(...
                            'klist',klist_tmp,'show',options.plot);
                        WAVECAR_loop_tmp = WAVECAR_loop1(:,BAND_index,:);
                        if project
                            WEIGHTCAR_SYM = ProjectionFunction(WAVECAR_loop_tmp);
                            kn = size(WAVECAR_loop_tmp,3);
                            Norb = size(WAVECAR_loop_tmp,1);
                            switch options.ProjectionMethod
                                case 'sign'
                                    LABELCAR1 = logical(WEIGHTCAR_SYM > 0);
                                    LABELCAR2 = logocal(WEIGHTCAR_SYM < 0);
                                    PlusNum = sum(LABELCAR1(:,round(end/2)));
                                    MinusNum = sum(LABELCAR2(:,round(end/2)));
                                    WAVECAR_loop_tmp1 = zeros(Norb,PlusNum,kn);
                                    WAVECAR_loop_tmp2 = zeros(Norb,MinusNum,kn);
                                    for j = 1:kn
                                        WAVECAR_loop_tmp1 = WAVECAR_loop_tmp(:,LABELCAR1(:,j),j);
                                        WAVECAR_loop_tmp2 = WAVECAR_loop_tmp(:,LABELCAR2(:,j),j);
                                    end
                                    [BFCAR1,~] = vasplib.wancenter_1D(WAVECAR_loop_tmp1);% [BFCAR1,BFWAVECAR1]
                                    [BFCAR2,~] = vasplib.wancenter_1D(WAVECAR_loop_tmp2);% [BFCAR1,BFWAVECAR2]
                                    BFCAR(:,i)=[BFCAR1;BFCAR2];
                                    WEIGHTCAR(:,i) = [ones(PlusNum,1);-ones(MinusNum,1)];
                                otherwise
                            end
                        else
                            [BFCAR(:,i),BF_WAVECAR(:,:,i)] = vasplib.wancenter_1D(WAVECAR_loop_tmp);
                        end
                    end
                case 'HR'
                    LABELCAR1 = 1:vasplibobj.WAN_NUM/2;
                    LABELCAR2 = (vasplibobj.WAN_NUM/2+1):vasplibobj.WAN_NUM;
                    for i = 1:options.knum_evol
                        klist_tmp = klist_s_1(i,:)+klist_s_2+kstart_s;
                        if project
                            if options.plot
                                [~,WAVECAR_loop1,ax] = vasplibobj.EIGENCAR_gen(...
                                    'klist',klist_tmp,...
                                    'convention','I','show',options.plot,'ax',ax,'Umat',AOperU,'subband',LABELCAR1);
                                drawnow;
                            else
                                [~,WAVECAR_loop1] = vasplibobj.EIGENCAR_gen(...
                                    'klist',klist_tmp,...
                                    'convention','I','printmode',false,'Umat',AOperU,'subband',LABELCAR1);
                            end
                            %error('debug');
                            if options.plot
                                [~,WAVECAR_loop2,ax] = vasplibobj.EIGENCAR_gen(...
                                    'klist',klist_tmp,...
                                    'convention','I','show',options.plot,'ax',ax,'Umat',AOperU,'subband',LABELCAR2);
                                drawnow;
                            else
                                [~,WAVECAR_loop2] = vasplibobj.EIGENCAR_gen(...
                                    'klist',klist_tmp,...
                                    'convention','I','printmode',false,'Umat',AOperU,'subband',LABELCAR2);
                            end
                            % If we use convention II, each wavefactor should
                            % give back the factor
                            % C^{nk}_j = C^{nk}_j_tilde * e^{-ik·tj}.
                            % not selected band here?
                            WAVECAR_loop_tmp1 = WAVECAR_loop1(:,BAND_index,:);
                            WAVECAR_loop_tmp2 = WAVECAR_loop2(:,BAND_index,:);
                            % normalize phases to get u instead of phi
                            %for j =1:size(WAVECAR_loop_tmp,3)
                            %    WAVECAR_loop_tmp(:,:,j) = WAVECAR_loop_tmp(:,:,j).* exp(-2i*pi*(vasplibobj.orbL*klist_tmp(j,:).'));
                            %end
                            % The last bloch state is the same as the first up to a phase factor
                            WAVECAR_loop_tmp1(:,:,end) = WAVECAR_loop_tmp1(:,:,1).* exp(-2i*pi*(vasplibobj.orbL*(klist_tmp(end,:)-klist_tmp(1,:)).'));
                            WAVECAR_loop_tmp2(:,:,end) = WAVECAR_loop_tmp2(:,:,1).* exp(-2i*pi*(vasplibobj.orbL*(klist_tmp(end,:)-klist_tmp(1,:)).'));
                            %WAVECAR_loop_tmp(:,:,end) = WAVECAR_loop_tmp(:,:,1).* exp(-2i*pi*(vasplibobj.orbL*(klist_tmp(end,:)-klist_tmp(1,:)).'));

                            %LABELCAR1 = logical(WEIGHTCAR_SYM > 0);
                            %LABELCAR2 = logical(WEIGHTCAR_SYM < 0);
                            PlusNum = length(BAND_index);
                            MinusNum = length(BAND_index);
                            %PlusNum = sum(LABELCAR1(:,round(end/2)));
                            %MinusNum = sum(LABELCAR2(:,round(end/2)));
                            %for j = 1:kn
                            %    WAVECAR_loop_tmp1 = WAVECAR_loop_tmp(:,LABELCAR1(:,j),j);
                            %    WAVECAR_loop_tmp2 = WAVECAR_loop_tmp(:,LABELCAR2(:,j),j);
                            %end
                            WAVECAR_loop_tmp1 = WAVECAR_loop_tmp1(:,BAND_index,:);
                            WAVECAR_loop_tmp2 = WAVECAR_loop_tmp2(:,BAND_index,:);
                            [BFCAR1,~] = vasplib.wancenter_1D(WAVECAR_loop_tmp1);% [BFCAR1,BFWAVECAR1]
                            [BFCAR2,~] = vasplib.wancenter_1D(WAVECAR_loop_tmp2);% [BFCAR1,BFWAVECAR2]
                            BFCAR(:,i)=[BFCAR1;BFCAR2];
                            WEIGHTCAR(:,i) = [ones(PlusNum,1);-ones(MinusNum,1)];
                        else
                            klist_tmp = klist_s_1(i,:)+klist_s_2+kstart_s;
                            if options.plot
                                [~,WAVECAR_loop,ax] = vasplibobj.EIGENCAR_gen(...
                                    'klist',klist_tmp,...
                                    'convention','II','show',options.plot,'ax',ax);
                                drawnow;
                            else
                                [~,WAVECAR_loop] = vasplibobj.EIGENCAR_gen(...
                                    'klist',klist_tmp,...
                                    'convention','II','printmode',false);
                            end
                            % If we use convention II, each wavefactor should
                            % give back the factor
                            % C^{nk}_j = C^{nk}_j_tilde * e^{-ik·tj}.
                            WAVECAR_loop_tmp = WAVECAR_loop(:,BAND_index,:);
                            % normalize phases to get u instead of phi
                            for j =1:size(WAVECAR_loop_tmp,3)
                                WAVECAR_loop_tmp(:,:,j) = WAVECAR_loop_tmp(:,:,j).* exp(-2i*pi*(vasplibobj.orbL*klist_tmp(j,:).'));
                            end
                            % The last bloch state is the same as the first up to a phase factor
                            WAVECAR_loop_tmp(:,:,end) = WAVECAR_loop_tmp(:,:,1).* exp(-2i*pi*(vasplibobj.orbL*(klist_tmp(end,:)-klist_tmp(1,:)).'));
                            [BFCAR(:,i),BF_WAVECAR(:,:,i)] = vasplib.wancenter_1D(WAVECAR_loop_tmp);
                        end
                    end
            end
        end
        function [nested_BFCAR,nested_BF_ALL,klist_l] = nested_WilsonLoop(vasplibobj,optionsK,options,optionsNested)
            arguments
                vasplibobj ;
                optionsK.knum_int    = 31;
                optionsK.kintegral   = [0,1,0];
                optionsK.kevolution  = [1,0,0];
                optionsK.cartesian = false;
                optionsK.dir_seq = [1,2,3];
                optionsK.dir_start = 'kcar';
                options.BAND_index = [];
                options.ax = handle([]);
                options.plot = false;
                options.V = [];
                optionsNested.kstart      = [0,0,0];
                optionsNested.knum_evol   = 51;
                optionsNested.nested_kevolution  = [0,0,1];
                optionsNested.nested_BAND_index = [];
            end
            % --- nargin
            switch class(vasplibobj)
                case {'Htrig','HK'}
                    NBAND = vasplibobj.Basis_num;
                case 'HR'
                    NBAND = vasplibobj.WAN_NUM;
            end
            if isempty(options.BAND_index)
                BAND_index = 1:(NBAND/2);
            else
                BAND_index = options.BAND_index;
            end
            if isempty(optionsNested.nested_BAND_index)
                nested_BAND_index = 1:(NBAND/4);
            else
                nested_BAND_index = options.nested_BAND_index;
            end
            %
            %kstart_r = options.kstart*Gk_;
            %
            optionscell = namedargs2cell(options);
            optionsK.knum_evol = optionsK.knum_int;
            kstart_s = optionsNested.kstart;
            knum_evol = optionsNested.knum_evol;
            % nested-1 fix kz semi-fix kx integral ky
            % k_z
            [~,klist_s_3,~,~] =...
                vasplib.kpathgen([[0,0,0];optionsNested.nested_kevolution],knum_evol,vasplibobj.Gk,vasplibobj.Gk);
            klist_l = zeros(size(klist_s_3,1),1);
            normklist_l = norm(optionsNested.nested_kevolution)/norm(klist_s_3(end,:));
            for i = 1:size(klist_s_3,1)
                klist_l(i) = norm(klist_s_3(i,:)*(eye(3)*2*pi))*normklist_l;
            end
            klist_l = klist_l + sum(sign(kstart_s))* norm(kstart_s*(eye(3)*2*pi))*normklist_l;
            klist_s_3 = klist_s_3+kstart_s;
            %
            nested_BFCAR = zeros(length(nested_BAND_index),knum_evol);
            nested_BF_ALL = zeros(length(nested_BAND_index),optionsK.knum_int,knum_evol);
            %
            pb = vasplib_tool_outer.CmdLineProgressBar('Nested Berry Phase caculating ');
            for k = 1:knum_evol
                % first WAN(y)
                % 在固定kx的情况下，计算沿着ky方向的Wilson loop并得到对应的本征矢量
                % 将先沿着ky方向再沿着kx方向的所有哈密顿量本征波函数存储
                % 注意ky首尾相连 而kx没有
                optionsKbk = optionsK;
                optionsKbk.kstart = klist_s_3(k,:);
                optionsKbk.knum_evol = optionsK.knum_int;
                optionsKbkcell = namedargs2cell(optionsKbk);
                [~,BF_WAVECAR,~,WAVELOOPCAR] = vasplibobj.WilsonLoop(optionsKbkcell{:},optionscell{:},'V',options.V,'LWAVE',true,'LWAVE',true);
                %OccupyBand = size(BFCAR,1);
                NBAND = size(WAVELOOPCAR,1);
                knum_evol_nested = optionsK.knum_int;
                knum_int_nested  = optionsK.knum_int;
                pmulist = zeros(length(nested_BAND_index),knum_evol_nested);
                for ik_evolution = 1:knum_evol_nested %(k_y direction)
                    NuWAVELOOP = zeros(NBAND,length(nested_BAND_index),optionsK.knum_int);%(k_x direction)
                    for ik_int= 1:knum_int_nested %(k_x direction)
                        WAVELOOP_ike_ikint = WAVELOOPCAR(:,:,ik_evolution,ik_int);
                        WCCWAVECAR_ikint = BF_WAVECAR(:,nested_BAND_index,ik_int);
                        NuWAVELOOP(:,:,ik_int)=sum(WAVELOOP_ike_ikint * WCCWAVECAR_ikint,2);
                    end
                    % The last  is the same as the first 
                    NuWAVELOOP(:,:,end) = NuWAVELOOP(:,:,1);%.* exp(-1i*(vasplibobj.orbL*vasplibobj.Rm*(klist_tmp(end,:)-klist_tmp(1,:)).'));
                    pmulist(:,ik_evolution) = vasplib.wancenter_1D(NuWAVELOOP);
                end
                % second Construct WANhamiltionian
                pmulist(pmulist<0) = pmulist(pmulist<0)+2*pi;
                nested_BFCAR(:,k) = mean(pmulist,2);
                pb.print(k,knum_evol,' ...');
            end
            pb.delete();
        end
    end
    methods
        function [pmulist,klist_l] = nestedWilsonLoop(vasplibobj,optionsK,options,optionsNested)
            arguments
                vasplibobj ;
                optionsK.knum_int    = 31;
                optionsK.knum_evol   = 51;
                optionsK.kstart      = [0,0,0];
                optionsK.kintegral   = [0,1,0];
                optionsK.kevolution  = [1,0,0];
                optionsK.cartesian = false;
                optionsK.dir_seq = [1,2,3];
                optionsK.dir_start = 'kcar';
                options.BAND_index = [];
                options.V = [];
                optionsNested.nested_BAND_index = [];
                optionsNested.script {mustBeMember(optionsNested.script,...
                    {...
                    'p_x(nu_y)','p_y(nu_x)',....
                    'p_y(nu_z)','p_z(nu_y)',...
                    'p_z(nu_x)','p_x(nu_z)',...
                    'p_1(nu_2)','p_2(nu_1)',...
                    'p_2(nu_3)','p_3(nu_2)',...
                    'p_3(nu_1)','p_1(nu_3)',...
                    ''})}= '';
            end
            % --- nargin
            switch class(vasplibobj)
                case {'Htrig','HK'}
                    NBAND = vasplibobj.Basis_num;
                case 'HR'
                    NBAND = vasplibobj.WAN_NUM;
            end
            if isempty(options.BAND_index)
                BAND_index = 1:(NBAND/2);
            else
                BAND_index = options.BAND_index;
            end
            if isempty(optionsNested.nested_BAND_index)
                nested_BAND_index = 1:(NBAND/4);
            else
                nested_BAND_index = options.nested_BAND_index;
            end
            switch optionsNested.script
                case {'p_x(nu_y)','p_1(nu_2)'} % Calculate nu_y(k_x)
                    optionsK.kintegral   = [0,1,0];
                    optionsK.kevolution  = [1,0,0];
                    optionsK.kstart      = [-0.5,0,0]+optionsK.kstart;
                case {'p_y(nu_z)','p_2(nu_3)'}
                    optionsK.kintegral   = [0,0,1];
                    optionsK.kevolution  = [0,1,0];
                    optionsK.kstart      = [0,-0.5,0]+optionsK.kstart;
                case {'p_z(nu_x)','p_3(nu_1)'}
                    optionsK.kintegral   = [1,0,0];
                    optionsK.kevolution  = [0,0,1];
                    optionsK.kstart      = [0,0,-0.5]+optionsK.kstart;
                case {'p_y(nu_x)','p_2(nu_1)'}
                    optionsK.kintegral   = [1,0,0];
                    optionsK.kevolution  = [0,1,0];
                    optionsK.kstart      = [0,-0.5,0]+optionsK.kstart;
                case {'p_z(nu_y)','p_3(nu_2)'}
                    optionsK.kintegral   = [0,1,0];
                    optionsK.kevolution  = [0,0,1];
                    optionsK.kstart      = [0,0,-0.5]+optionsK.kstart;
                case {'p_x(nu_z)','p_1(nu_3)'}
                    optionsK.kintegral   = [0,0,1];
                    optionsK.kevolution  = [1,0,0];
                    optionsK.kstart      = [-0.5,0,0]+optionsK.kstart;
            end
            %
            if contains(string(optionsNested.script),["x","y","z"])
                optionsK.cartesian = true;
                optionsK.dir_seq = [1,2,3];
                optionsK.dir_start = 'kcar';
            end
            optionsKcell = namedargs2cell(optionsK);
            [~,WannierCenterWAVECAR,~,WAVELOOPCAR] = vasplibobj.WannierCenter(optionsKcell{:},'V',options.V,'LWAVE',true,'BAND_index',BAND_index);
            % WAVELOOPCAR format: NBAND Occupiband k_intgrel k_evolution
            % takecare WAVELOOPCAR(:,:,:,end) should be WAVELOOPCAR(:,:,:,1)
            % takes nu_y(k_x) as example, k_int = y dir k_evol = x dir ->
            % exchange two direction
            %[kloop1_s,kloop2_s,kloop1_r,kloop2_r,klist_l,kstart_s,kstart_r] = vasplib.kloop2D(vasplibobj.Rm,optionsKcell{:});
            % p_x(nu_y)(k_y)
            knum_evol_nested = optionsK.knum_int;
            knum_int_nested  = optionsK.knum_evol;
            %
            pmulist = zeros(length(nested_BAND_index),knum_evol_nested);
            for ik_evolution = 1:knum_evol_nested %(k_y direction)
                NuWAVELOOP = zeros(NBAND,length(nested_BAND_index),optionsK.knum_evol);%(k_x direction)
                for ik_int= 1:knum_int_nested %(k_x direction)
                    WAVELOOP_ike_ikint = WAVELOOPCAR(:,:,ik_evolution,ik_int);
                    WCCWAVECAR_ikint = WannierCenterWAVECAR(:,nested_BAND_index,ik_int);
                    NuWAVELOOP(:,:,ik_int)=sum(WAVELOOP_ike_ikint * WCCWAVECAR_ikint,2);
                end
                % The last  is the same as the first 
                NuWAVELOOP(:,:,end) = NuWAVELOOP(:,:,1);%.* exp(-1i*(vasplibobj.orbL*vasplibobj.Rm*(klist_tmp(end,:)-klist_tmp(1,:)).'));
                pmulist(:,ik_evolution) = vasplib.WannierCenter1D(NuWAVELOOP);
            end
            klist_l = linspace(-pi,pi,knum_evol_nested);
        end
        function [WannierCenterCAR,WannierCenterWAVECAR,klist_l,WAVELOOPCAR] = WannierCenter(vasplibobj,optionsK,options)
            arguments
                vasplibobj ;
                optionsK.knum_int    = 31;
                optionsK.knum_evol   = 51;
                optionsK.kstart      = [0,0,0];
                optionsK.kintegral   = [0,1,0];
                optionsK.kevolution  = [1,0,0];
                optionsK.cartesian = false;
                optionsK.dir_seq = [1,2,3];
                optionsK.dir_start = 'kcar';
                options.V = [];
                options.Accuracy = 1E-12;
                options.LWAVE = false;
                options.BAND_index = [];
                options.script {mustBeMember(options.script,...
                    {...
                    'nu_x(k_y)','nu_y(k_x)',....
                    'nu_y(k_z)','nu_z(k_y)',...
                    'nu_z(k_x)','nu_x(k_z)',...
                    'nu_1(k_2)','nu_2(k_1)',...
                    'nu_2(k_3)','nu_3(k_2)',...
                    'nu_3(k_1)','nu_1(k_3)',...
                    ''})}= '';
            end
            % --- nargin
            switch class(vasplibobj)
                case {'Htrig','HK'}
                    NBAND = vasplibobj.Basis_num;
                case 'HR'
                    NBAND = vasplibobj.WAN_NUM;
            end
            if isempty(options.BAND_index)
                BAND_index = 1:(NBAND/2);
            else
                BAND_index = options.BAND_index;
            end
            % use script
            switch options.script
                case {'nu_x(k_y)','nu_1(k_2)'}
                    optionsK.kintegral   = [1,0,0];
                    optionsK.kevolution  = [0,1,0];
                    optionsK.kstart      = [0,-0.5,0];
                case {'nu_y(k_z)','nu_2(k_3)'}
                    optionsK.kintegral   = [0,1,0];
                    optionsK.kevolution  = [0,0,1];
                    optionsK.kstart      = [0,0,-0.5];
                case {'nu_z(k_x)','nu_3(k_1)'}
                    optionsK.kintegral   = [0,0,1];
                    optionsK.kevolution  = [1,0,0];
                    optionsK.kstart      = [-0.5,0,0];
                case {'nu_y(k_x)','nu_2(k_1)'}
                    optionsK.kintegral   = [0,1,0];
                    optionsK.kevolution  = [1,0,0];
                    optionsK.kstart      = [-0.5,0,0];
                case {'nu_z(k_y)','nu_3(k_2)'}
                    optionsK.kintegral   = [0,0,1];
                    optionsK.kevolution  = [0,1,0];
                    optionsK.kstart      = [0,-0.5,0];
                case {'nu_x(k_z)','nu_1(k_3)'}
                    optionsK.kintegral   = [1,0,0];
                    optionsK.kevolution  = [0,0,1];
                    optionsK.kstart      = [0,0,-0.5];
                case {''}
            end
            %
            if contains(string(options.script),["x","y","z"])
                optionsK.cartesian = true;
                optionsK.dir_seq = [1,2,3];
                optionsK.dir_start = 'kcar';
            end
            optionsKcell = namedargs2cell(optionsK);
            % Gen kloop 2D
            [kloop1_frac,kloop2_frac,kloop1_cart,kloop2_cart,klist_l,kstart_frac,kstart_cart] = vasplib.kloop2D(vasplibobj.Rm,optionsKcell{:});
            % prepare
            WannierCenterCAR = zeros(length(BAND_index),optionsK.knum_evol);
            WannierCenterWAVECAR = zeros(length(BAND_index),length(BAND_index),optionsK.knum_evol);
            if options.LWAVE
                if isempty(options.V )
                    V = diag((1:NBAND) *sqrt(options.Accuracy));
                else
                    V = options.V ;
                end
                %EIGENCAR = zeros(length(BAND_index),optionsK.knum_int,optionsK.knum_evol);
                WAVELOOPCAR = zeros(NBAND,length(BAND_index),optionsK.knum_int,optionsK.knum_evol);
            end
            % EIGEN
            switch class(vasplibobj)
                case {'Htrig','HK'}
                    if isempty(vasplibobj.orbL)
                        vasplibobj.orbL = zeros(vasplibobj.Basis_num,3);
                    end
                    for i = 1:optionsK.knum_evol
                        klist_tmp = kloop1_cart(i,:)+kloop2_cart+kstart_cart;
                        [E,WAVECAR_loop] = vasplibobj.EIGENCAR_gen(...
                            'klist',klist_tmp,'printmode',false);
                        WAVECAR_loop_tmp = WAVECAR_loop(:,BAND_index,:);
                        % The last bloch state is the same as the first up to a phase factor
                        WAVECAR_loop_tmp(:,:,end) = WAVECAR_loop_tmp(:,:,1).* exp(-1i*(vasplibobj.orbL*vasplibobj.Rm*(klist_tmp(end,:)-klist_tmp(1,:)).'));
                        % WAVECAR_loop_tmp(:,:,end) = WAVECAR_loop_tmp(:,:,1);
                        if options.LWAVE
                            WAVECAR_loop_tmp = vasplib.cleanWAVECAR(WAVECAR_loop_tmp,E(BAND_index,:),V,options.Accuracy);
                            WAVELOOPCAR(:,:,:,i) = WAVECAR_loop_tmp;
                        end
                        [WannierCenterCAR(:,i),WannierCenterWAVECAR(:,:,i)] = vasplib.WannierCenter1D(WAVECAR_loop_tmp);
                    end
                case 'HR'
                    for i = 1:optionsK.knum_evol
                        klist_tmp = kloop1_frac(i,:)+kloop2_frac+kstart_frac;
                        [E,WAVECAR_loop] = vasplibobj.EIGENCAR_gen(...
                            'klist',klist_tmp,...
                            'convention','II','printmode',false);
                        % If we use convention II, each wavefactor should
                        % give back the factor
                        % C^{nk}_j = C^{nk}_j_tilde * e^{-ik·tj}.
                        WAVECAR_loop_tmp = WAVECAR_loop(:,BAND_index,:);
                        % normalize phases to get u instead of phi
                        for j =1:size(WAVECAR_loop_tmp,3)
                            WAVECAR_loop_tmp(:,:,j) = WAVECAR_loop_tmp(:,:,j).* exp(-2i*pi*(vasplibobj.orbL*klist_tmp(j,:).'));
                        end
                        % The last bloch state is the same as the first up to a phase factor
                        WAVECAR_loop_tmp(:,:,end) = WAVECAR_loop_tmp(:,:,1).* exp(-2i*pi*(vasplibobj.orbL*(klist_tmp(end,:)-klist_tmp(1,:)).'));
                        if options.LWAVE
                            WAVECAR_loop_tmp = vasplib.cleanWAVECAR(WAVECAR_loop_tmp,E(BAND_index,:),V,options.Accuracy);
                            WAVELOOPCAR(:,:,:,i) = WAVECAR_loop_tmp;
                        end
                        [WannierCenterCAR(:,i),WannierCenterWAVECAR(:,:,i)] = vasplib.WannierCenter1D(WAVECAR_loop_tmp);
                    end
            end
        end
    end
    methods
        function [WindingNumber,WL] = WindingNumber(Ham_obj,GammaOper,kloop,options)
            % \mathcal{N}=\frac{1}{4 \pi i} \oint_{C} \operatorname{Tr} \sigma_{z} \mathcal{H}_{\mathrm{eff}}^{-1}(\boldsymbol{q}) \nabla_{\boldsymbol{q}} \mathcal{H}_{\mathrm{eff}}(\boldsymbol{q}) \cdot d \boldsymbol{q}
            % https://journals.aps.org/prb/abstract/10.1103/PhysRevB.99.121106 (8)
            arguments
                Ham_obj;
                GammaOper =[]; % the chiral symmetry oper of Ham_obj
                kloop double = [];
                options.plot logical = false;
                options.dir = [1,2,3];
            end
            dir = options.dir ;
            % prepare dH & dH_dk
            switch class(Ham_obj)
                case "HR"
                    switch Ham_obj.Type
                        case {'sparse'}
                            HnumList = reshape(full(cell2mat(Ham_obj.HnumL)),Ham_obj.WAN_NUM,Ham_obj.WAN_NUM,Ham_obj.NRPTS);
                        case {'mat'}
                            HnumList = Ham_obj.HnumL;
                        case {'list'}
                            Ham_obj = Ham_obj.rewind();
                            HnumList = Ham_obj.HnumL;
                    end
                    NRPTS_ = Ham_obj.NRPTS;
                    vectorList = double(Ham_obj.vectorL);
                    vectorList_r = vectorList * Ham_obj.Rm;
                    % partial R
                    HnumLpx = 1i*pagemtimes(reshape(vectorList_r(:,dir(1)),[1 1 NRPTS_]),HnumList);
                    HnumLpy = 1i*pagemtimes(reshape(vectorList_r(:,dir(2)),[1 1 NRPTS_]),HnumList);
                    HnumLpz = 1i*pagemtimes(reshape(vectorList_r(:,dir(3)),[1 1 NRPTS_]),HnumList);
                    % partial titj  % we dont consider tij mat here check!
                    Ham_obj = Ham_obj.tjmti_gen();
                    tji_mat_frac = Ham_obj.tjmti{2};
                    HnumLpA_tji = tji_mat_frac(:,:,dir(1));
                    HnumLpB_tji = tji_mat_frac(:,:,dir(2));
                    HnumLpC_tji = tji_mat_frac(:,:,dir(3));
                    kloop_cart = kloop * Ham_obj.Gk;
                case {"Htrig"}
                    [dH_dkx_fun,dH_dky_fun,dH_dkz_fun] = Ham_diff(Ham_obj);
                    HfunTmp = Ham_obj.Hfun;
                    kloop_cart = kloop;
                case 'HK'
                    syms k_x  k_y k_z real;
                    TargetH =  Ham_obj.Hk_num ;
                    TargetH_Pk_x = diff(TargetH,k_x);
                    TargetH_Pk_y = diff(TargetH,k_y);
                    TargetH_Pk_z = diff(TargetH,k_z);
                    HfunTmp = matlabFunction(TargetH,'Vars',[k_x k_y k_z]);
                    dH_dkx_fun = matlabFunction(TargetH_Pk_x,'Vars',[k_x k_y k_z]);
                    dH_dky_fun = matlabFunction(TargetH_Pk_y,'Vars',[k_x k_y k_z]);
                    dH_dkz_fun = matlabFunction(TargetH_Pk_z,'Vars',[k_x k_y k_z]);
                    kloop_cart = kloop;
                otherwise
                    syms k_x k_y k_z real;
                    TargetH =  Ham_obj.Hsym ;
                    TargetH_Pk_x = diff(TargetH,k_x);
                    TargetH_Pk_y = diff(TargetH,k_y);
                    TargetH_Pk_z = diff(TargetH,k_z);
                    HfunTmp = matlabFunction(TargetH,'Vars',[k_x k_y k_z]);
                    dH_dkx_fun = matlabFunction(TargetH_Pk_x,'Vars',[k_x k_y k_z]);
                    dH_dky_fun = matlabFunction(TargetH_Pk_y,'Vars',[k_x k_y k_z]);
                    dH_dkz_fun = matlabFunction(TargetH_Pk_z,'Vars',[k_x k_y k_z]);
                    kloop_cart = kloop;
            end
            % diff k
            % check 1 == end
            if ~Oper.isclose(kloop_cart(1,:),kloop_cart(end,:))
                error('loop enforced!');
            end
            dkloop = (diff(kloop_cart,1,1) + diff(kloop_cart([end-1,1:end-1],:),1,1))/2;
            AL = zeros(size(dkloop));
            switch class(Ham_obj)
                case "HR"
                    for kn = 1:size(dkloop,1)
                        % efactor R
                        FactorListki = exp(1i*2*pi*vectorList*kloop(kn,:).');
                        % HRmat
                        HRmat = sum(pagemtimes(HnumList,reshape(FactorListki,[1 1 NRPTS_])),3);
                        % pHRmat
                        HRmatpA = sum(pagemtimes(HnumLpx,reshape(FactorListki,[1 1 NRPTS_])),3);
                        HRmatpB = sum(pagemtimes(HnumLpy,reshape(FactorListki,[1 1 NRPTS_])),3);
                        HRmatpC = sum(pagemtimes(HnumLpz,reshape(FactorListki,[1 1 NRPTS_])),3);
                        % efactor orb
                        kjiL_A =  tji_mat_frac(:,:,1).*kloop(kn,1);
                        kjiL_B =  tji_mat_frac(:,:,2).*kloop(kn,2);
                        kjiL_C =  tji_mat_frac(:,:,3).*kloop(kn,3);
                        Hmat_tji = exp(1i*2*pi*(kjiL_A+kjiL_B+kjiL_C));
                        %
                        Hmat_tjipA = Hmat_tji.* HnumLpA_tji;
                        Hmat_tjipB = Hmat_tji.* HnumLpB_tji;
                        Hmat_tjipC = Hmat_tji.* HnumLpC_tji;
                        %
                        vxk = HRmatpA.*Hmat_tji + HRmat.*Hmat_tjipA;% vx partial_A_tmp
                        vyk = HRmatpB.*Hmat_tji + HRmat.*Hmat_tjipB;% vy partial_B_tmp
                        vzk = HRmatpC.*Hmat_tji + HRmat.*Hmat_tjipC;% vx partial_A_tmp
                        H = HRmat.*Hmat_tji;
                        Ax = trace(GammaOper/H*vxk);
                        Ay = trace(GammaOper/H*vyk);
                        Az = trace(GammaOper/H*vzk);
                        AL(kn,:)= [Ax,Ay,Az];
                    end
                otherwise
                    for kn = 1:size(dkloop,1)
                        kx = kloop_cart(kn,1); ky = kloop_cart(kn,2); kz = kloop_cart(kn,3);
                        dH_dkx = dH_dkx_fun(kx,ky,kz);dH_dky = dH_dky_fun(kx,ky,kz);
                        dH_dkz = dH_dkz_fun(kx,ky,kz);H = HfunTmp(kx,ky,kz);
                        Ax = trace(GammaOper/H*dH_dkx);
                        Ay = trace(GammaOper/H*dH_dky);
                        Az = trace(GammaOper/H*dH_dkz);
                        AL(kn,:)= [Ax,Ay,Az];
                    end
            end
            WL = dot(AL,dkloop,2)/(4i*pi);
            WindingNumber = sum(WL);
        end
        function [HL,pHL] = Topo1DpreHL(vasplibobj,klist,options)
            
        end
        function WAVECAR_loop = Topo1DpreWAVECAR(vasplibobj,klist,options)
            arguments
                vasplibobj
                klist
                options.BAND_index = [];
                options.ax = handle([]);
                options.plot = false;
                options.LWAVE = false;
            end
            % --- nargin
            switch class(vasplibobj)
                case {'Htrig','HK'}
                    %
                case 'HR'
                    vasplibobj.Basis_num = vasplibobj.WAN_NUM;
            end
            if isempty(options.BAND_index)
                set_divide = 2;
                switch class(vasplibobj)
                    case {'Htrig','HK','HR'}
                        BAND_index = 1:(vasplibobj.Basis_num/set_divide);
                    case ''
                        BAND_index = 1:(vasplibobj.Basis_num/set_divide);
                    otherwise
                        BAND_index = 1:(vasplibobj.Basis_num/set_divide);
                end
            else
                BAND_index = options.BAND_index;
            end
            switch class(vasplibobj)
                case {'Htrig','HK'}
                    klist_tmp = klist;
                    [~,WAVECAR_loop] = vasplibobj.EIGENCAR_gen(...
                        'klist',klist_tmp,'printmode',false);
                    WAVECAR_loop = WAVECAR_loop(:,BAND_index,:);
                    % The last bloch state is the same as the first up to a phase factor
                    WAVECAR_loop(:,:,end) = WAVECAR_loop(:,:,1).* exp(-1i*(vasplibobj.orbL*vasplibobj.Rm*(klist_tmp(end,:)-klist_tmp(1,:)).'));
                case 'HR'
                    klist_tmp = klist;
                    [~,WAVECAR_loop] = vasplibobj.EIGENCAR_gen(...
                        'klist',klist_tmp,...
                        'convention','II','printmode',false);
                    % If we use convention II, each wavefactor should
                    % give back the factor
                    % C^{nk}_j = C^{nk}_j_tilde * e^{-ik·tj}.
                    WAVECAR_loop = WAVECAR_loop(:,BAND_index,:);
                    % normalize phases to get u instead of phi
                    for j =1:size(WAVECAR_loop,3)
                        WAVECAR_loop(:,:,j) = WAVECAR_loop(:,:,j).* exp(-2i*pi*(vasplibobj.orbL*klist_tmp(j,:).'));
                    end
                    % The last bloch state is the same as the first up to a phase factor
                    WAVECAR_loop(:,:,end) = WAVECAR_loop(:,:,1).* exp(-2i*pi*(vasplibobj.orbL*(klist_tmp(end,:)-klist_tmp(1,:)).'));
                otherwise
                    klist_tmp = klist;
                    [~,WAVECAR_loop] = vasplib.EIGENSOLVE(vasplibobj.Hfun,klist_tmp,vasplibobj.Basis_num);
                    WAVECAR_loop = WAVECAR_loop(:,BAND_index,:);
                    % The last bloch state is the same as the first up to a phase factor
                    WAVECAR_loop(:,:,end) = WAVECAR_loop(:,:,1).* exp(-1i*(vasplibobj.orbL*vasplibobj.Rm*(klist_tmp(end,:)-klist_tmp(1,:)).'));
            end
            
        end
        function [BP,WAVECAR_loop] = BP_1D(vasplibobj,klist,options,optionsWAVE)
            arguments
                vasplibobj
                klist
                options.BAND_index = [];
                options.ax = handle([]);
                options.plot = false;
                optionsWAVE.LWAVE = false;
            end
            optionsCell = namedargs2cell(options);
            WAVECAR_loop = Topo1DpreWAVECAR(vasplibobj,klist,optionsCell{:});
            if optionsWAVE.LWAVE

            end
            BP = sum(vasplib.wancenter_1D(WAVECAR_loop));
        end
    end
    methods(Static)
        function WAVECAR = cleanWAVECAR(WAVECAR,EIGENCAR,V,Accuracy)
            arguments 
                WAVECAR;
                EIGENCAR;
                V;
                Accuracy = 1e-12;
            end
            kn = size(EIGENCAR,2);
            %ki = 1;
            %EIGEN_ki = EIGENCAR(:,ki);
            %DegenPairBase = vasplib.checkDegen(EIGEN_ki,Accuracy);
            %for ki = 2:kn
            %    EIGEN_ki = EIGENCAR(:,ki);
            %    DegenPair = vasplib.checkDegen(EIGEN_ki,Accuracy);
            %    if ~isequal(DegenPairBase,DegenPair)
            %        error('!');
            %    end
            %end
            %WAVECAR = vasplib.smoothDegen2(WAVECAR);
            for ki = 1:kn
                EIGEN_ki = EIGENCAR(:,ki);
                DegenPair = vasplib.checkDegen(EIGEN_ki,Accuracy);
                WAVE_ki = WAVECAR(:,:,ki);
                for i =1:size(DegenPair,1)
                    % $\hat{H}^{(0)} \psi^{(1)}+\hat{H}^{\prime} \psi^{(0)}=E^{(0)} \psi^{(1)}+E^{(1)} \psi^{(0)}$
                    WAVE_ki(:,DegenPair(i,1):DegenPair(i,2),:) = vasplib.smoothDegen(WAVE_ki(:,DegenPair(i,1):DegenPair(i,2)),V); 
                end
                WAVECAR(:,:,ki) = WAVE_ki;
            end
        end
        function WAVEFuncL = smoothDegen2(WAVEFuncL)
            % Wave1 Todo
            WAVEFuncLabs =abs(WAVEFuncL(:,1,:));
            [~,maxWAVEFuncLabs] = max(WAVEFuncLabs,[],1);
            [~,SelectOrb] = max(sum(maxWAVEFuncLabs,2));
            for i = 1:size(WAVEFuncL,2)
                for k = 1:size(WAVEFuncL,3)
                    WAVEFuncLtmp = exp(-1i*angle(WAVEFuncL(SelectOrb,1,i)))*WAVEFuncL(:,:,i);
                    WAVEFuncL(:,:,i) = WAVEFuncLtmp;
                end
            end
            %WAVEFuncLimag =imag(WAVEFuncL);
        end
        function WAVEFuncL = smoothDegen(WAVEFuncL,V)
            % Degenerate perturbation theory
            [A,U] = eig(WAVEFuncL'*V*WAVEFuncL);
            %disp(WAVEFuncL)
            [A,~] = park.sorteig(U,A);
            WAVEFuncL = WAVEFuncL*A;
        end
        function DegenPair = checkDegen(EIGEN,Accuracy)
            arguments
                EIGEN {mustBeVector}
                Accuracy = 1e-8;
            end
            Nband= length(EIGEN);
            DEGEN_init = (EIGEN(2:Nband)- EIGEN(1:(Nband-1)))<Accuracy;
            tail = false;
            DegenPair= [];
            for i = 1:numel(DEGEN_init)
                if DEGEN_init(i)==1
                    if  ~tail
                        DEGEN_Pairlittle = i-1+1;
                        tail = true;
                    end
                elseif DEGEN_init(i)==0 
                    if tail
                        DEGEN_Pairlittle = [DEGEN_Pairlittle,i-1+1];
                        tail = false;
                        DegenPair= [DegenPair;DEGEN_Pairlittle];
                    end
                end
            end
            if tail
                DEGEN_Pairlittle = [DEGEN_Pairlittle,i+1];
                tail = false;
                DegenPair= [DegenPair;DEGEN_Pairlittle];
            end
            if ~tail
                return;
            else
                error('?');
            end
            %DEGEN_char = reshape(char(string(DEGEN_init)),1,[]);
        end
        function [WAVECAR_loop_] = modify_WAVECAR(WAVECAR_loop,BF_WAVECAR)
            WAVECAR_loop_ = zeros(size(WAVECAR_loop,1),size(BF_WAVECAR,2),size(BF_WAVECAR,3));
            for i = 1:size(WAVECAR_loop,3) % For k evolution 
                uk = WAVECAR_loop(:,:,i);
                vk = BF_WAVECAR(:,:,i);
                for j= 1:size(BF_WAVECAR,2) % for all occupied WAN band
                    WAVECAR_loop_(:,j,i) = sum(uk.*vk(:,j).',2);
                end
            end
        end
        function [WCC,WCCvec,HWan] = WannierCenter1D(WAVECAR_loop)
            %
            % W(C)=M^{k0,k1}·...·M^{kn−1,kn,}
            % M_{m, n}^{\mathbf{k}_{i}, \mathbf{k}_{j}}=\left\langle u_{m, \mathbf{k}_{i}} \mid u_{n, \mathbf{k}_{j}}\right\rangle
            % M_{m, n}^{ki,kj} = <u_{m,ki}|u_{n,kj}>
            % Take care!
            % Bloch states : Psi_{n,k} = e^{i*k*r}*u_{n,k}
            % So the u_{n,k} should be normalized phases instead of phi
            HWan = eye(size(WAVECAR_loop,2));
            Nloop = size(WAVECAR_loop,3);
            for kj = 1:Nloop-1
                HWan = HWan* vasplib.BerryConnection(WAVECAR_loop(:,:,kj),WAVECAR_loop(:,:,kj+1));
            end
            
            %HWan = eye(size(WAVECAR_loop,2));
            %Nloop = size(WAVECAR_loop,3);
            %WAVECAR_loop(:,:,Nloop+1) =WAVECAR_loop(:,:,1);
            %for kj = 1:Nloop-1
            %    if kj == Nloop
            %        F = vasplib.BerryConnection(WAVECAR_loop(:,:,kj),WAVECAR_loop(:,:,1));
            %    else
            %        F = vasplib.BerryConnection(WAVECAR_loop(:,:,kj),WAVECAR_loop(:,:,kj+1));
            %    end
            %    [U,~,V] = svd(F);
            %    HWan = HWan*(U*V');
            %    %Wan = Wan*F;
            %end
            %HWan = WAVECAR_loop(:,:,1)';
            %Nloop = size(WAVECAR_loop,3);
            %Wan = eye(size(WAVECAR_loop,1));
            %WAVECAR_loop(:,:,Nloop+1) =WAVECAR_loop(:,:,1);
            %for kj = 2:Nloop-1
            %    Wan = Wan*(WAVECAR_loop(:,:,kj)*WAVECAR_loop(:,:,kj)');
            %end
            %HWan = HWan*Wan*WAVECAR_loop(:,:,end);%HWan*Wan*WAVECAR_loop(:,:,1)
            %Wan
            [WCCvec,WCCU] = eig(HWan);
            %Ei = mod((angle(diag(WCCU))/(2*pi)),1);
            Ei = mod(real(log(diag(WCCU))/(2*pi*1i)),1);
            [WCCvec,WCC] = park.sorteig(Ei,WCCvec);
        end
        function [BF,BF_WAVE,Wan] = wancenter_1D(WAVECAR_loop,mode)
            %
            % W(C)=M^{k0,k1}·...·M^{kn−1,kn,}
            % M_{m, n}^{\mathbf{k}_{i}, \mathbf{k}_{j}}=\left\langle u_{m, \mathbf{k}_{i}} \mid u_{n, \mathbf{k}_{j}}\right\rangle
            % M_{m, n}^{ki,kj} = <u_{m,ki}|u_{n,kj}>
            % Take care!
            % Bloch states : Psi_{n,k} = e^{i*k*r}*u_{n,k}
            % So the u_{n,k} should be normalized phases instead of phi
            formula = 2;
            if nargin < 2
                switch formula
                    case 1
                        Nband = size(WAVECAR_loop,2);
                        Wan = eye(Nband);
                        Nloop = size(WAVECAR_loop,3);
                        %WAVECAR_loop(:,:,Nloop+1) =WAVECAR_loop(:,:,1);
                        for k =1:Nloop-1
                            WanTmp = zeros(Nband) ;
                            for i = 1:Nband
                                for j = 1:Nband
                                    WanTmp(i,j) = WAVECAR_loop(:,j,k+1).'*conj(WAVECAR_loop(:,i,k));
                                end
                            end
                            Wan = Wan*WanTmp;
                        end
                        %pause(1);
                    case 2
                        Wan = eye(size(WAVECAR_loop,2));
                        Nloop = size(WAVECAR_loop,3);
                        %WAVECAR_loop(:,:,Nloop+1) =WAVECAR_loop(:,:,1);
                        for kj = 1:Nloop-1
                            Wan = Wan* vasplib.BerryConnection(WAVECAR_loop(:,:,kj),WAVECAR_loop(:,:,kj+1));
                        end
                        %Wan = Wan;%* vasplib.BerryConnection(WAVECAR_loop(:,:,end),WAVECAR_loop(:,:,1));
                    case 3
                        Wan = eye(size(WAVECAR_loop,2));
                        Nloop = size(WAVECAR_loop,3);
                        %WAVECAR_loop(:,:,Nloop+1) =WAVECAR_loop(:,:,1);
                        for kj = 1:Nloop-1
                            if kj == Nloop
                                F = vasplib.BerryConnection(WAVECAR_loop(:,:,kj),WAVECAR_loop(:,:,1));
                            else
                                F = vasplib.BerryConnection(WAVECAR_loop(:,:,kj),WAVECAR_loop(:,:,kj+1));
                            end
                            [U,~,V] = svd(F);
                            Wan = Wan*(U*V');
                            %Wan = Wan*F;
                        end
                end
                switch formula
                    case 1
                        [BF_WAVE,Ei] = eig(Wan);
                        Ei = angle((Ei));
                        [BF_WAVE,Ei] = park.sorteig(Ei,BF_WAVE);
                        BF=diag(Ei);
                    case {2,3}
                        [BF_WAVE,Ei] = eig(Wan);
                        Ei = angle((Ei));
                        [BF_WAVE,Ei] = park.sorteig(Ei,BF_WAVE);
                        BF=diag(Ei);
                end

            elseif strcmp(mode,'nested')
                %                 if size(site_weight,2) ==size(WAVECAR_loop,3)
                %                     Wan = (WAVECAR_loop(:,:,1).*site_weight(:,1))';
                %                     for kj = 2:size(WAVECAR_loop,3)
                %                         Wan = Wan* (WAVECAR_loop(:,:,kj).*site_weight(:,kj))*(WAVECAR_loop(:,:,kj).*site_weight(:,kj))';
                %                     end
                %                     Wan = Wan* WAVECAR_loop(:,:,1).*site_weight(:,1);
                %                 else
                %                     Wan = (WAVECAR_loop(:,:,1).*site_weight(:,1))';
                %                     for kj = 2:size(WAVECAR_loop,3)
                %                         Wan = Wan* (WAVECAR_loop(:,:,kj).*site_weight(:,1))*(WAVECAR_loop(:,:,kj).*site_weight(:,1))';
                %                     end
                %                     Wan = Wan* (WAVECAR_loop(:,:,1).*site_weight(:,1));
                %                 end
                Wan = WAVECAR_loop(:,:,1)';
                for kj = 2:size(WAVECAR_loop,3)
                    Wan = Wan* WAVECAR_loop(:,:,kj)*WAVECAR_loop(:,:,kj)';
                end
                Wan = Wan* WAVECAR_loop(:,:,1);
                [BF_WAVE,Ei] = eig(Wan);
                BF=angle(sum(diag(Ei)));
            end
        end
        function F = BerryConnection(W1,W2)
            %             Norb1 = size(W1,2);
            %             Norb2 = size(W2,2);
            %             for n = 1:Norb1
            %                 for m = 1:Norb2
            %                     F(n,m) = W1(:,n)'* W2(:,m);
            %                 end
            %             end
            F = W1'*W2;
        end
        function F = Berryphase_2D()

        end
        % arrangement of klist

        % integral method
        function [BC_2D,BCL] = BerryCuvature_2D(WAVECAR,sizemesh,options)
            arguments
                WAVECAR
                sizemesh
                options.sum = true;
            end
            if options.sum 
                Nband = 1;
            else
                Nband = size(WAVECAR,2);
            end
            sizemesh_WAVECAR = sizemesh+1;% Contaion Edage % For intgrel method
            kn = size(WAVECAR,3);
            BC_2D = zeros([sizemesh,Nband]);
            [i_list,j_list] = ind2sub(sizemesh_WAVECAR,1:kn);
            % ind_VV_list = 1:kn;
            % Vk
            Vk_list = 1:kn;
            seqL = i_list <=sizemesh(1) & j_list <=sizemesh(2);
            Vk_list = Vk_list(seqL);
            i_list = i_list(seqL);
            j_list = j_list(seqL);
            % Vk1
            Vk1_i_list = i_list+1;%Vk1_i_list(Vk1_i_list>sizemesh(1)) = 1;
            ind_Vk1_list = sub2ind(sizemesh_WAVECAR,Vk1_i_list,j_list);
            % Vk2
            Vk1_j_list = j_list+1;%Vk1_j_list(Vk1_j_list>sizemesh(2)) = 1;
            ind_Vk2_list = sub2ind(sizemesh_WAVECAR,i_list,Vk1_j_list);
            % Vk1k2
            ind_Vk1k2_list = sub2ind(sizemesh_WAVECAR,Vk1_i_list,Vk1_j_list);
            %
            kn = length(Vk_list);
            [iL,jL]= ind2sub(sizemesh,1:kn);
            if options.sum
                for k = 1:kn
                    VV = WAVECAR(:,:,Vk_list(k));
                    Vk1 = WAVECAR(:,:,ind_Vk1_list(k));
                    Vk2 = WAVECAR(:,:,ind_Vk2_list(k));
                    Vk1k2 = WAVECAR(:,:,ind_Vk1k2_list(k));
                    BC_2D(iL(k),jL(k),:) = vasplib.BerryCuvature_Discrete_2D(VV,Vk1,Vk2,Vk1k2);
                end
            else
                for i =1:Nband
                    for k = 1:kn
                        VV = WAVECAR(:,i,Vk_list(k));
                        Vk1 = WAVECAR(:,i,ind_Vk1_list(k));
                        Vk2 = WAVECAR(:,i,ind_Vk2_list(k));
                        Vk1k2 = WAVECAR(:,i,ind_Vk1k2_list(k));
                        BC_2D(iL(k),jL(k),i) = vasplib.BerryCuvature_Discrete_2D(VV,Vk1,Vk2,Vk1k2);
                    end
                end
            end
            if nargout == 2
                BCL = [];
            end
        end
        function BC = BerryCuvature_Discrete_2D(VV,Vk1,Vk2,Vk1k2)
            % Consider Vector them for speed up
            % VV    W_d_0;
            % Vk1   W_d_1;  略偏离kx的波函数
            % Vk2   W_d_2;  略偏离ky的波函数
            % Vk1k2 W_d_12; 略偏离kx，ky的波函数
            % ----- approach 1 -----
            % U(1) like variable
            % 10.1143/JPSJ.74.1674
            % https://topocondmat.org/w4_haldane/ComputingChern.html
            Uk1 = det(VV'*Vk1);Uk1 = Uk1/norm(Uk1);
            Uk2 = det(VV'*Vk2);Uk2 = Uk2/norm(Uk2);
            Uk1_k2 = det(Vk2'*Vk1k2);Uk1_k2 = Uk1_k2/norm(Uk1_k2);
            Uk2_k1 = det(Vk1'*Vk1k2);Uk2_k1 = Uk2_k1/norm(Uk2_k1);
            % ----- approach 2 -----
            %             U = VV'*(Vk1*Vk1')*(Vk1k2*Vk1k2')*(Vk2*Vk2')*VV;
            %             BC = angle(eig(U));
            % berry curvature
            BC = log(Uk1*Uk2_k1/(Uk1_k2*Uk2));
            BC = imag(BC);
            %BC = sum(BC);
        end
        function BC = nBerryCuvature_Discrete_2D(VV,Vk1,Vk2,Vk1k2)
            % Consider Vector them for speed up
            % VV    W_d_0;
            % Vk1   W_d_1;  略偏离kx的波函数
            % Vk2   W_d_2;  略偏离ky的波函数
            % Vk1k2 W_d_12; 略偏离kx，ky的波函数
            % ----- approach 1 -----
            % ----- approach 2 -----
            U = VV'*(Vk1*Vk1')*(Vk1k2*Vk1k2')*(Vk2*Vk2')*VV;
            BC = angle(eig(U));
        end
        function BC = BerryCuvature_Discrete_3D(VV,Vk1,Vk2,Vk3,Vk1k2,Vk2k3,Vk3k1)
            BC(1) = vasplib.BerryCuvature_Discrete_2D(VV,Vk2,Vk3,Vk2k3);
            BC(2) = vasplib.BerryCuvature_Discrete_2D(VV,Vk3,Vk1,Vk3k1);
            BC(3) = vasplib.BerryCuvature_Discrete_2D(VV,Vk1,Vk2,Vk1k2);
        end
    end
    % analytical way, not for real calcultion
    methods(Static)
        % ------------ Berryphase Line  ------------
        % fun
        function F = BerryPhaseLine_fun(Fun,kloopr,options)
            arguments
                Fun function_handle;
                kloopr double;
                options.Rm = [];
                options.oneshot = true;
                options.plot = false;
                options.car = true;
                options.funtype = 'Wfun';
            end
            if strcmp(options.funtype,'Wfun')
                WfunType = true;
                Wfun = Fun;
                FfunType = false;
            elseif strcmp(options.funtype,'Ffun')
                FfunType = true;
                Ffun = Fun;
                WfunType=false;
            end
            if isempty(options.Rm)
                Rm = POSCAR_read;
            else
                Rm = options.Rm;
            end
            if options.car
            else
                Gk = 2*pi*eye(3)/Rm;
                kloopr = kloopr * Gk;
            end
            if size(kloopr,2)<3
                kloopr = [kloopr,zeros(size(kloopr,1),3-size(kloopr,2))];
            end
            kn = size(kloopr,1);
            if WfunType
                % $\phi \equiv-\operatorname{Im} \ln \left[\left\langle u_{0} \mid u_{1}\right\rangle\left\langle u_{1} \mid u_{2}\right\rangle \cdots\left\langle u_{N-1} \mid u_{0}\right\rangle\right]=-\sum_{j=0}^{N-1} \operatorname{Im} \ln \left\langle u_{j} \mid u_{j+1}\right\rangle$
                W1 = Wfun(kloopr(1,1),kloopr(1,2),kloopr(1,3));
                WL = zeros(size(W1,1),kn);
                for i = 1:kn
                    WL(:,i) =  Wfun(kloopr(i,1),kloopr(i,2),kloopr(i,3));
                end
                dF = -imag(log(vasplib.BerryConnection(WL(:,1),WL(:,2))));
                F = dF;
            elseif FfunType
                dkloopr = diff(kloopr);
                %dkloopr =[kloopr(1,:) - kloopr(end,:);dkloopr];
                dF = Ffun(kloopr(1,1),kloopr(1,2),kloopr(1,3),dkloopr(1,1),dkloopr(1,2),dkloopr(1,3));
                F = dF;
            end
            warning off;
            if options.plot
                dBF = sum(dF);
                fig = figure('PaperType','a4letter','PaperSize',[16 8],'Color','white','Units','normalized','Position',[0.1,0.1,0.8,0.6]);
                ax0 = subplot(1,3,1,'LineWidth',1.5,'FontSize',24,'FontName',"Helvetica",'Parent',fig);
                title(ax0,'integral path');
                box(ax0,'on');
                hold(ax0,'all');
                axis(ax0,'equal');
                [fig,ax0] = BZplot(Rm,'r',0.2,fig,ax0);
                ax1 = subplot(1,3,2,'LineWidth',1.5,'FontSize',24,'FontName',"Helvetica",'Parent',fig);
                hold(ax1,'all');
                xlabel(ax1,'kpath num');
                ylabel(ax1,'dBerryphase');
                ax2 = subplot(1,3,3,'LineWidth',1.5,'FontSize',24,'FontName',"Helvetica",'Parent',fig);
                hold(ax2,'all');
                xlabel(ax2,'kpath num');
                ylabel(ax2,'Berryphase');
                count = 0;
                dBf_old = dBF;
                K_last = kloopr(1,:);  
                axis(ax0,'equal');
            end
            for i = 2:kn-1
                k_x = kloopr(i,1);
                k_y = kloopr(i,2);
                k_z = kloopr(i,3);
                K = kloopr(i,:);
                dk = norm(K-K_last); % test
                if WfunType
                    dF = -imag(log(vasplib.BerryConnection(WL(:,i),WL(:,i+1))))*dk;
                elseif FfunType
                    dk_x = dkloopr(i,1);
                    dk_y = dkloopr(i,2);
                    dk_z = dkloopr(i,3);
                    dF = Ffun(k_x,k_y,k_z,dk_x,dk_y,dk_z)*dk;
                end
                F = F +dF;
                if options.plot
                    dBF = real(sum(dF));
                    count = count +1;
                    line(ax0,[K_last(1),k_x],[K_last(2),k_y],[K_last(3),k_z],'LineWidth',3);
                    %(ax0,[K_last(1),k_x],[K_last(2),k_y],[K_last(3),k_z],'LineWidth',3);
                    K_last = K;
                    xlabel(ax1,string(K(1))+","+string(K(2))+","+string(K(3)));
                    area(ax1,[count-1,count],[dBf_old,dBF]);
                    dBf_old= dBF;
                    stem(ax2,count,sum(F));
                    if options.oneshot
                    else
                        drawnow;
                    end
                end
            end
            if options.plot
                title(ax2,'ChernNumber = '+string((mod(real(sum(F)),(2*pi)))));
            end
            F = mod(real(F),2*pi);
        end
        function F = BerryPhaseLine_definition_sym_i(Hsym,options)
            arguments
                Hsym sym;
                options.BAND_index = [];
            end
            % temp 2d
            syms k_x k_y k_z delta_k_x delta_k_y delta_k_z real;
            [Norb,~] = size(Hsym);
            set_divide = 2;
            if isempty(options.BAND_index)
                BAND_index = 1:(Norb/set_divide);
            else
                BAND_index = options.BAND_index;
            end
            [W,~] = eig(Hsym);
            W = W(:,BAND_index);
            [~,Nbands] = size(W);
            F = zeros([Nbands,1],class(W));
            for i = 1:Nbands
                Eigenvector_sym = (vasplib.NomalizeEigenvector(W(:,i)));
                A_x = vasplib.BerryConnection_definition(Eigenvector_sym,k_x);
                A_y = vasplib.BerryConnection_definition(Eigenvector_sym,k_y);
                A_z = vasplib.BerryConnection_definition(Eigenvector_sym,k_z);
                F(i) = ([delta_k_x delta_k_y delta_k_z]*[A_x,A_y,A_x].');
                %WL(:,:,i) = HsymP1*u(:,i)*u(:,i)'*HsymP2;
            end
        end
        % ------------ ***************  ------------ 
        % ------------ wannier center  ------------       
        % fun WL direct
        function [BFCAR,BF_WAVECAR,klist_l] = WilsonLoop_fun(Hfun,optionsK,options)
            arguments
                Hfun ;
                optionsK.knum_int    = 31;
                optionsK.knum_evol   = 51;
                optionsK.kstart      = [0,0,0];
                optionsK.kintegral   = [0,1,0];
                optionsK.kevolution  = [1,0,0];
                optionsK.cartesian = false;
                optionsK.dir_seq = [1,2,3];
                optionsK.dir_start = 'kcar';
                options.WLfun = false;
                options.ax = handle([]);
                options.plot = false;
                options.V = [];
                options.Accuracy = 1E-12;
                options.LWAVE = false;
                options.BAND_index = [];
                options.script {mustBeMember(options.script,...
                    {...
                    'nu_x(k_y)','nu_y(k_x)',....
                    'nu_y(k_z)','nu_z(k_y)',...
                    'nu_z(k_x)','nu_x(k_z)',...
                    'nu_1(k_2)','nu_2(k_1)',...
                    'nu_2(k_3)','nu_3(k_2)',...
                    'nu_3(k_1)','nu_1(k_3)',...
                    ''})}= '';
            end
            % --- nargin
            switch class(Hfun)
                case {'sym'}
                    syms k_x k_y k_z real;
                    Hfun = matlabFunction(Hfun,'Vars',[k_x,k_y,k_z]);
                case 'handle'

                otherwise
            end
            if ~options.WLfun
                testH = Hfun(0,0,0);
                NBAND = size(testH,1);
                if isempty(options.BAND_index)
                    set_divide = 2;
                    BAND_index = 1:(Norb/set_divide);
                else
                    BAND_index = options.BAND_index;
                end
            else
            end
            % use script
            switch options.script
                case {'nu_x(k_y)','nu_1(k_2)'}
                    optionsK.kintegral   = [1,0,0];
                    optionsK.kevolution  = [0,1,0];
                    optionsK.kstart      = [0,-0.5,0];
                case {'nu_y(k_z)','nu_2(k_3)'}
                    optionsK.kintegral   = [0,1,0];
                    optionsK.kevolution  = [0,0,1];
                    optionsK.kstart      = [0,0,-0.5];
                case {'nu_z(k_x)','nu_3(k_1)'}
                    optionsK.kintegral   = [0,0,1];
                    optionsK.kevolution  = [1,0,0];
                    optionsK.kstart      = [-0.5,0,0];
                case {'nu_y(k_x)','nu_2(k_1)'}
                    optionsK.kintegral   = [0,1,0];
                    optionsK.kevolution  = [1,0,0];
                    optionsK.kstart      = [-0.5,0,0];
                case {'nu_z(k_y)','nu_3(k_2)'}
                    optionsK.kintegral   = [0,0,1];
                    optionsK.kevolution  = [0,1,0];
                    optionsK.kstart      = [0,-0.5,0];
                case {'nu_x(k_z)','nu_1(k_3)'}
                    optionsK.kintegral   = [1,0,0];
                    optionsK.kevolution  = [0,0,1];
                    optionsK.kstart      = [0,0,-0.5];
                case {''}
            end
            %
            if contains(string(options.script),["x","y","z"])
                optionsK.cartesian = true;
                optionsK.dir_seq = [1,2,3];
                optionsK.dir_start = 'kcar';
            end
            optionsKcell = namedargs2cell(optionsK);
            % Gen kloop 2D
            [~,~,kloop1_cart,kloop2_cart,klist_l,~,kstart_cart] = vasplib.kloop2D(vasplibobj.Rm,optionsKcell{:});
            % prepare
            BFCAR = zeros(length(BAND_index),optionsK.knum_evol);
            BF_WAVECAR = zeros(length(BAND_index),length(BAND_index),optionsK.knum_evol);
            if options.LWAVE
                if isempty(options.V )
                    V = diag((1:NBAND) *sqrt(options.Accuracy));
                else
                    V = options.V ;
                end
                %EIGENCAR = zeros(length(BAND_index),optionsK.knum_int,optionsK.knum_evol);
                WAVELOOPCAR = zeros(NBAND,length(BAND_index),optionsK.knum_int,optionsK.knum_evol);
            end
            for i = 1:optionsK.knum_evol
                klist_tmp = kloop1_cart(i,:)+kloop2_cart+kstart_cart;
                [E,WAVECAR_loop] = vasplib.EIGENSOLVE(Hfun, ...
                    klist_tmp,NBAND);
                WAVECAR_loop_tmp = WAVECAR_loop(:,BAND_index,:);
                % The last bloch state is the same as the first up to a phase factor
                WAVECAR_loop_tmp(:,:,end) = WAVECAR_loop_tmp(:,:,1);
                if options.LWAVE
                    WAVECAR_loop_tmp = vasplib.cleanWAVECAR(WAVECAR_loop_tmp,E(BAND_index,:),V,options.Accuracy);
                    WAVELOOPCAR(:,:,:,i) = WAVECAR_loop_tmp;
                end
                [BFCAR(:,i),BF_WAVECAR(:,:,i)] = vasplib.wancenter_1D(WAVECAR_loop_tmp);
            end

        end
        % ------------ Berryphase Curvature  ------------ 
        % fun BC direct
        function [BCCAR,Grid,klist_r_plot] = BerryCuvature_fun(Hfun,optionsK,options,optionsPlot) 
            arguments
                Hfun ;
                optionsK.knum1   = 51;
                optionsK.knum2   = 51;
                optionsK.kstart  = [-0.5,-0.5,0];
                optionsK.kdir1   = [1,0,0];
                optionsK.kdir2   = [0,1,0];
                optionsK.cartesian = false;
                optionsK.dir_seq = [1,2,3];
                optionsK.dir_start = 'kcar';
                options.knodes = 100;
                options.Rm = [];
                options.plot = false;
                options.car = true;
                options.BAND_index = [];
                options.BCfun = false;
                optionsPlot.oneshot = true;
                optionsPlot.view = [0,90];
            end
            if isempty(options.Rm)
                Rm = POSCAR_read;
            else
                Rm = options.Rm;
            end
            Gk = (2*pi*eye(3)/Rm).';
            %
            optionsKcell = namedargs2cell(optionsK);
            [klist_cart,~,klist_r_plot,sizemesh,Gk_,Grid] = vasplib.kmesh2D(Rm,optionsKcell{:},'full',true);
            if isa(Hfun,'sym')
                Hfun = matlabFunction(Hfun,'Vars',[sym('k_x'),sym('k_y'),sym('k_z')]);
            end
            if ~options.BCfun
                testH = Hfun(0,0,0);
                Norb = size(testH,1);
                if isempty(options.BAND_index)
                    set_divide = 2;
                    BAND_index = 1:(Norb/set_divide);
                else
                    BAND_index = options.BAND_index;
                end
                % integral method
                [~,BC_WAVECAR] = vasplib.EIGENSOLVE(Hfun,klist_cart,Norb);
                BC_WAVECAR = BC_WAVECAR(:,BAND_index,:);
                %reshape(BC_WAVECAR(:,BAND_index,:),[vasplibobj.Basis_num options.knum1*options.knum2]);
                BCCAR = vasplib.BerryCuvature_2D( BC_WAVECAR ,sizemesh);
            else
                BCfun = Hfun;
                BCCAR = zeros(sizemesh);
                for ki = 1:numel(BCCAR)
                    K = klist_r_plot(ki,:);
                    Bcki = BCfun(K(1),K(2),K(3));
                    BCCAR(ki) = sum(Bcki);
                end
            end
            if options.plot
                knum1 = sizemesh(1);
                knum2 = sizemesh(2);
                dk_1 = (optionsK.kdir1)/knum1*Gk_;
                dk_2 = (optionsK.kdir2)/knum2*Gk_;
                if ~options.BCfun
                    dS = 1;
                else
                    dS = cross(dk_1,dk_2);
                end
                optionsPlotcell = namedargs2cell(optionsPlot);
                [~,~] = vasplib.ShowSurfIntegral(Gk,klist_cart,dk_1,dk_2,BCCAR(:)*norm(dS)/(2*pi),optionsPlotcell{:});
            end
        end
        % origin definition for full symbolic obj
        function BC = BC_definition(Hsym,para1,para2,epsilon,options)
            arguments
                Hsym sym;
                para1 sym;
                para2 sym;
                epsilon = 1;
                options.BAND_index = [];
            end
            [Norb,~] = size(Hsym);
            set_divide = 2;
            if isempty(options.BAND_index)
                BAND_index = 1:(Norb/set_divide);
            else
                BAND_index = options.BAND_index;
            end
            %HsymP1 = diff(Hsym,para1);
            %HsymP2 = diff(Hsym,para2);
            [W,~] = eig(Hsym);
            W = W(:,BAND_index);
            [~,Nbands] = size(W);
            BC = zeros([Nbands,1],class(W));
            %u = zeros(size(W),class(W));
            for i = 1:Nbands
                Eigenvector_sym = (vasplib.NomalizeEigenvector(W(:,i)));
                A_1 = vasplib.BerryConnection_definition(Eigenvector_sym,para1);
                A_2 = vasplib.BerryConnection_definition(Eigenvector_sym,para2);
                BC(i) = vasplib.BerryCurvature_definition(A_1,A_2,para1,para2);
                %WL(:,:,i) = HsymP1*u(:,i)*u(:,i)'*HsymP2;
            end
            BC = BC*epsilon;
            %E = simplify(diag(U));
            % $\begin{aligned} \Omega_{\mu \nu}^{n}(\boldsymbol{R}) &=\frac{\partial}{\partial R_{\mu}} A_{\nu}^{n}(\boldsymbol{R})-\frac{\partial}{\partial R_{\nu}} A_{\mu}^{n}(\boldsymbol{R}) \\ &=i\left[\frac{\partial}{\partial R_{\mu}}\left\langle n(\boldsymbol{R})\left|\frac{\partial}{\partial R_{\nu}}\right| n(\boldsymbol{R})\right\rangle-(\nu \leftrightarrow \mu)\right] \end{aligned}$
            %Bc = simplify(diff(A_2,para1)-diff(A_1,para2));
        end
        % kubo formula for full symbolic obj
        function Bc = BC_kubo_sym(Hsym,para1,para2,epsilon,options)
            if nargin < 4
                epsilon = 1;
                options.BAND_index = [];
            end
            [Norb,~] = size(Hsym);
            set_divide = 2;
            if isempty(options.BAND_index)
                BAND_index = 1:(Norb/set_divide);
            else
                BAND_index = options.BAND_index;
            end
            HsymP1 = diff(Hsym,para1);
            HsymP2 = diff(Hsym,para2);
            [W,U] = eig(Hsym);
            E = simplify(diag(U));
            WL = zeros(Norb,Norb,Norb,class(W));
            %WL2 = WL;% full
            for i = 1:Norb
                W(:,i) = (vasplib.NomalizeEigenvector(W(:,i)));
                WL(:,:,i) = HsymP1*W(:,i)*W(:,i)'*HsymP2;
                %WL2(:,:,i) = HsymP2*W(:,i)*W(:,i)'*HsymP1;% full
            end
            u = W(:,BAND_index);
            [Norb,Nbands] = size(u);
            Bc = zeros([Nbands,1],class(W))+1i*zeros([Nbands,1],class(W));
            for i = 1:Nbands
                for j = 1:Norb
                    %$\Omega_{\mu \nu}^{n}(\boldsymbol{R})=i \sum_{n^{\prime} \neq n} \frac{\left\langle n\left|\frac{\partial H}{\partial R_{\mu}}\right| n^{\prime}\right\rangle\left\langle n^{\prime}\left|\frac{\partial H}{\partial R_{\nu}}\right| n\right\rangle-(\nu \leftrightarrow \mu)}{\left(E_{n}-E_{n^{\prime}}\right)^{2}}$
                    if j ~= i && simplify(E(i)-E(j)) ~= sym(0)
                        %                         BC(i) = BC(i) + (1i/(E(i)-E(j)))*((u(:,i)'*HsymP1*u(:,i))*(u(:,i)'*HsymP2*u(:,i))- ...
                        %                 (u(:,i)'*HsymP2*u(:,j))*(u(:,j)'*HsymP1*u(:,i)));
                        % if Hermite
                        Bc(i) = Bc(i) +(1/(E(i)-E(j))^2)*(u(:,i)'*WL(:,:,j)*u(:,i));
                        % full
                        %Bc(i) = Bc(i)+(1/(E(i)-E(j))^2)*((u(:,i)'*WL(:,:,j)*u(:,i)) -(u(:,i)'*WL2(:,:,j)*u(:,i)) ) ;
                    end
                end
            end
            Bc = (Bc*2i*epsilon);
            %Bc = (Bc*1i*epsilon);% full
        end
        % tool
        function A_1 = BerryConnection_definition(Eigenvector_sym,para1)
            A_1 = 1i*Eigenvector_sym'*diff(Eigenvector_sym,para1);
        end
        function Bc = BerryCurvature_definition(A_1,A_2,para1,para2)
            Bc = diff(A_2,para1)-diff(A_1,para2);
        end
        % symbolic diff 2
        function Bc = Berry_curvature_D2(Eigenvector_sym,para1,para2)
            A_1 = Eigenvector_sym'*1i*diff(Eigenvector_sym,para1);
            A_2= Eigenvector_sym'*1i*diff(Eigenvector_sym,para2);
            Bc = simplify(diff(A_2,para1)-diff(A_1,para2));
        end
    end
    %% eig
    methods(Static)
        function HoutL = HCAR_gen(Hfun,klist_cart,Norb)
            arguments
                Hfun function_handle;
                klist_cart double;
                Norb;
            end
            kn = size(klist_cart,1);
            HoutL = zeros(Norb,Norb,kn);
            for i = 1:kn
                Input = num2cell(klist_cart(i,:));
                HoutL(:,:,i) = Hfun(Input{:});
            end
        end
        function [EIGENCAR,WAVECAR,HoutL] = EIGENSOLVE(Hfun,klist_cart,Norb,opt)
            arguments
                Hfun
                klist_cart = [];
                Norb = -1;
                opt.Hermitian = true;
            end
            if isa(Hfun,'vasplib')
                if ~isempty(Hfun.klist_cart)
                    klist_cart = Hfun.klist_cart;
                end
                Norb = Hfun.Basis_num;
                Hfun = Hfun.Hfun;
            end
            if Norb < 0
                Norb = length(Hfun(0,0,0));
            end
            %HoutL = vasplib.HCAR_gen(Hfun,klist_cart,Norb);
            if isempty(klist_cart)
                warning('Empty klist!');
                EIGENCAR = [];
                WAVECAR = [];
                HoutL = [];
                return;
            end
            kn = size(klist_cart,1);
            EIGENCAR = zeros(Norb,kn);
            if nargout >1
                WAVECAR  = zeros(Norb,Norb,kn);
            end
            if nargout >2
                HoutL = zeros(Norb,Norb,kn);
            end
            for i = 1:kn
                Input = num2cell(klist_cart(i,:));
                Hout= Hfun(Input{:});
                if opt.Hermitian
                    Hout = (Hout+Hout')/2;
                end
                [A, U]=eig(Hout);
                [A, U]= park.sorteig(U,A);
                EIGENCAR(:,i) = diag(U);
                if nargout >1
                    WAVECAR(:,:,i)=A;
                end
                if nargout >2
                    HoutL(:,:,i) = Hout;
                end
            end
            
        end
        function EIGENCAR = arrangeEIGENCAR(EIGENCAR,REFCAR,method,opt)
            arguments
                EIGENCAR ;
                REFCAR = EIGENCAR; 
                method {mustBeMember(method,{'NonHermitian','PBAND_single'})} = 'NonHermitian';
                opt.SelectL = [];
                opt.Ecut = [-3,3];
                opt.Disp = false;
            end
            switch  method
                case 'NonHermitian'
                    %XList = real(REFCAR);
                    %YList = imag(REFCAR);
                    Nbands = size(EIGENCAR,1);
                    Nk = size(EIGENCAR,2);
                    EIGENCAR_OUT = EIGENCAR;
                    NormelSeq = (1:Nbands).';
                    % dsearchn
                    for i = 2:Nk
                        XListTemp = real(EIGENCAR_OUT(:,i-1:i));
                        YListTemp = imag(EIGENCAR_OUT(:,i-1:i));
                        
                        PQ = [XListTemp(:,1),YListTemp(:,1)];
                        P  = [XListTemp(:,2),YListTemp(:,2)];
                        [k,dist] = dsearchn(P,PQ);
                        
                        if ~isequal(k,NormelSeq)
                            if opt.Disp
                                disp(i);
                                disp(k);
                            end
                            for j = i:Nk
                                EIGENCAR_OUT(:,j) = EIGENCAR_OUT(k,j);
                            end
                        else

                        end
                    end
                    EIGENCAR = EIGENCAR_OUT;
                case 'PBAND_single'
                    %WEIGHTCAR  = REFCAR;
                    Nbands = size(EIGENCAR,1);
                    Nk = size(EIGENCAR,2);
                    EIGENCAR_OUT = EIGENCAR;
                    NormelSeq = (1:Nbands).';
                    [NBANDS,NK,NPROJECTION]= size(REFCAR);
                    if NBANDS ~= Nbands || Nk ~= NK
                        error('REFCAR and EIGENCAR mismatch!');
                    end
                    if isempty(opt.SelectL)
                        SelectL = 1:NPROJECTION;
                    else
                        SelectL = opt.SelectL;
                    end
                    Ecut =  opt.Ecut;
                    % dsearchn
                    for i = 2:Nk
                        %XListTemp = real(EIGENCAR_OUT(:,i-1:i));
                        %YListTemp = imag(EIGENCAR_OUT(:,i-1:i));
                        PQ = reshape(REFCAR(:,i-1,:),NBANDS,NPROJECTION);
                        P  = reshape(REFCAR(:,i,:),NBANDS,NPROJECTION);
                        Eselect1 = Ecut(1) < EIGENCAR_OUT(:,i-1) & EIGENCAR_OUT(:,i-1) < Ecut(2) ;
                        Eselect2 = Ecut(1) < EIGENCAR_OUT(:,i) & EIGENCAR_OUT(:,i) < Ecut(2) ;
                        Eselect = find(Eselect1 & Eselect2);
                        PQ = PQ(Eselect,SelectL);
                        P = P(Eselect,SelectL);
                        [k,~] = dsearchn(P,PQ);
                        realk = [(1:min(Eselect)-1),k'+min(Eselect)-1,(max(Eselect)+1):NBANDS];
                        if ~isequal(realk,NormelSeq)
                            %disp(i);
                            %disp(k);
                            %for j = i:Nk
                            %    EIGENCAR_OUT(:,j) = EIGENCAR_OUT(k,j);
                            %end
                             EIGENCAR_OUT(:,i:Nk) = EIGENCAR_OUT(realk,i:Nk);
                        else

                        end
                    end
                    EIGENCAR = EIGENCAR_OUT;
            end


        end
    end
    %% math tools
    methods(Static)
        function pagenew = page_mtimes_matrix(page,mat)
            if length(mat) ~= size(page,2)
                raise ValueError('the length of the mat must equals with the pages.')
            end
            pagenew = page;
            for i = 1:size(page,3)
                pagenew(:,:,i) = page(:,:,i)*mat;
            end
        end
        function pagenew = matrix_mtimes_page(mat,page)
            if length(mat) ~= size(page,2)
                raise ValueError('the length of the mat must equals with the pages.')
            end
            pagenew = page;
            for i = 1:size(page,3)
                pagenew(:,:,i) = mat*page(:,:,i);
            end
        end
        function pagenew = matrixtimespage(mat,page)
            if length(mat) ~= size(page,3)
                raise ValueError('the length of the mat must equals with the pages.')
            end
            pagenew = page;
            if isvector(mat)
                for i = 1:length(mat)
                    pagenew(:,:,i) = mat(i)*page(:,:,i);
                end
            else
                for i = 1:length(mat)
                    pagenew(:,:,i) = sum(vasplib.matrixtimespage(mat(i,:),page),3);
                end
            end
        end
        function [c,t] = coeffsAll(p,varargin)
            [Ni,Nj] = size(p);
            for i = 1:Ni
                for j = 1:Nj
                    switch length(varargin)
                        case 0
                            [c{i,j},t{i,j}] = coeffs(p(i,j));
                        otherwise
                            [c{i,j},t{i,j}] = coeffs(p(i,j),varargin);
                    end
                end
            end
        end
        function Equation_list = isolateAll(Equation_list,Symvar_list)
            if nargin < 2
                Symvar_list = Equation_list;
            end
            zeroeq = sym(0) ==sym(0);
            if isequal(size(Equation_list),size(Symvar_list))
                for i = 1:numel(Equation_list)
                    Equation_list_tmp  = simplify(Equation_list(i));
                    symvartmp = symvar(simplify(Symvar_list(i)));
                    if  ~isequal(zeroeq,Equation_list_tmp) && ~isempty(symvartmp)
                        try
                            Equation_list(i) = isolate(Equation_list_tmp,symvartmp(1));
                        catch
                            fprintf('cant find the solution of Equation_list(%d):\n',i);
                            disp(Equation_list(i));
                        end
                    end
                end
            elseif isequal(size(Equation_list(:,:,1)),size(Symvar_list))
                %Nsym = numel(Symvar_list);
                SizeE = Equation_list ;
                for i = 1:numel(Equation_list)
                    [row,col,~] = ind2sub(SizeE,i);
                    symvartmp = symvar(Symvar_list(row,col));
                    if Symvar_list(row,col)~=sym(0) && ~isequal(zeroeq,Equation_list(i))
                        Equation_list(i) = isolate(Equation_list(i),symvartmp(1));
                    end
                end
            elseif isvector(Symvar_list)
                % not be implemented

            else
                if Symvar_list~=sym(0)
                    Symvar_list = symvar(Symvar_list);
                    for i = 1:numel(Equation_list)
                        if ~isequal(zeroeq,Equation_list(i))
                            Equation_list(i) = isolate(Equation_list(i),Symvar_list);
                        end
                    end
                end

            end
        end
        function [factor_list_1,factor_list_2] = factorAll(SymListRow)
            factor_list_1 = sym(ones(size(SymListRow)));
            factor_list_2 = sym(zeros(size(SymListRow)));
            for i = 1:length(SymListRow)
                F = factor(SymListRow(i));
                if length(F)==2
                    factor_list_1(i) = F(1);
                    factor_list_2(i)  = F(2);
                else
                    factor_list_2(i)  = F(1);
                end
            end
        end
        function SymListRow = cleanVar(SymListRow,Accuracy)
            Accur  = round(-log(Accuracy)/log(10));
            nSymListRow = length(SymListRow);
            notChoose = true(ones(1,nSymListRow));
            for i = 1:nSymListRow
                [coeffs_list,symvar_L] = coeffs(SymListRow(i));
                jChoose = false;
                sumCoeffsTmp = 0;
                for j = 1:length(coeffs_list)
                    CoeffsTmp = coeffs(coeffs_list(j));
                    try
                        dCoeffsTmp = double(vpa(CoeffsTmp,Accur));
                        if dCoeffsTmp >1e30
                            warning('addtional numerical error!');
                        end
                    catch
                        jChoose = true;
                        break;
                    end
                    if ~isempty(dCoeffsTmp)
                        sumCoeffsTmp = sumCoeffsTmp +abs(dCoeffsTmp);
                    end
                end

                if ~jChoose
                    if sumCoeffsTmp>Accuracy
                        dcoeffs_list = roundn(double(coeffs_list),-Accur);
                        jChoose = true;
                        SymListRow(i) = sum(dcoeffs_list.*symvar_L);
                        %disp(sumCoeffsTmp);
                        %disp(vpa(SymListRow(i),6));
                    end
                end
                notChoose(i) = ~jChoose;
            end
            SymListRow(notChoose) = sym(0);
        end
        function Pmat = Pvector(A)
            if nargin < 2
                mode = 'left';
            end
            if strcmp(mode,'left')
                Pmat = A.'*A/(A*A.');
            else
                Pmat = A*A.'/(A.'*A);
            end
        end
        function Pmat = Pplane(A,mode)
            % project a 3D vector to a plane
            if nargin < 2
                mode = 'left';
            end
            if strcmp(mode,'left')
                Pmat = A.'*inv(A*A.')*A;
            else
                Pmat = A*inv(A.'*A)*A.';
            end
        end
        function Gknew = CartisianMat(Gk,dir_seq,kstart)
            if nargin < 2
                dir_seq = [1,2,3];
            end
            if nargin < 3
                kstart = 'k_x';
            end
            switch kstart
                case 'kcar'
                    Gknew = diag(diag(Gk));
                case 'k_x'
                    k_x = Gk(dir_seq(1),:);
                    Pk_x = vasplib.Pvector(k_x);
                    if  dir_seq(2) == 0
                        k_y = Gk(2,:);
                    else
                        k_y = Gk(dir_seq(2),:) - Gk(dir_seq(2),:)*Pk_x;
                    end
                    Pk_xy = vasplib.Pplane([k_x;k_y]);
                    k_z = Gk(dir_seq(3),:) - Gk(dir_seq(3),:)*Pk_xy;
                    Gknew =  [k_x;k_y;k_z];
                case 'k_y'
                case 'k_z'
                    k_z = Gk(dir_seq(3),:);
                    Pk_z = vasplib.Pvector(k_z);
                    k_x = Gk(dir_seq(1),:) - Gk(dir_seq(1),:)*Pk_z;
                    Pk_zx = vasplib.Pplane([k_z;k_x]);
                    k_y = Gk(dir_seq(2),:) - Gk(dir_seq(2),:)*Pk_zx;
                    Gknew =  [k_x;k_y;k_z];
                    Gknew = Gknew(dir_seq,:);
            end
            % debug
            %             V_Gk=dot(Gk(1,:),cross(Gk(2,:),Gk(3,:)));
            %             V_Gknew=dot(Gknew(1,:),cross(Gknew(2,:),Gknew(3,:)));
            %             disp([V_Gk,V_Gknew]);
        end
        function M  = P2M(P)
            P_index = abs(P);
            P_phase = P./P_index;%exp(1i*angle(P));
            M  = zeros(length(P));
            for i =1:numel(P_index)
                M(i,P_index(i)) = P_phase(i);
            end
        end
        function Eigenvetor=NomalizeEigenvector(Eigenvetor)
            Eigenvetor = ...
                Eigenvetor/((Eigenvetor'*Eigenvetor)^(1/length(Eigenvetor)));
        end
    end
    %% Decomposition
    methods(Static)
        function [H_sym_Gamma,H_sym_Gamma_L,H_latex_Gamma] = GammaDecomposition(H_sym)
            if ~isequal(   size(H_sym) , [4,4])
                error('Gamma Decomposition requires 4*4 Ham!');
            end
            Smat_inv = gamma_matric.S();
            Gamma_L = gamma_matric.L();
            H_sym_L = sym(zeros(1,16));
            H_sym_Gamma_L = sym(zeros(1,16));
            %
            tmp_mat_r = real(H_sym);
            tmp_mat_i = imag(H_sym);
            H_sym_L(1:4) = diag(tmp_mat_r);
            H_sym_L(5:7) = diag(tmp_mat_r,1);
            H_sym_L(8:9) = diag(tmp_mat_r,2);
            H_sym_L(10)  = diag(tmp_mat_r,3);
            H_sym_L(11:13) = diag(tmp_mat_i,1);
            H_sym_L(14:15) = diag(tmp_mat_i,2);
            H_sym_L(16)  = diag(tmp_mat_i,3);
            %
            for i = 1:16
                Label_tmp = find(Smat_inv(i,:));
                H_sym_Gamma_L(Label_tmp) = H_sym_Gamma_L(Label_tmp)+H_sym_L(i)*Smat_inv(i,Label_tmp);
            end
            H_sym_Gamma_L = simplify(H_sym_Gamma_L);
            H_sym_Gamma = sym(0);
            H_latex_Gamma = 'H = ';
            count = 0;
            for i = 1:16
                if H_sym_Gamma_L(i)~=sym(0)
                    count = count+1;
                    H_sym_Gamma = H_sym_Gamma + H_sym_Gamma_L(i)*Gamma_L(i);
                    str1 = latex(H_sym_Gamma_L(i));
                    str2 = latex(Gamma_L(i));
                    if count > 1
                        H_latex_Gamma = [H_latex_Gamma,'+','\left(',str1,'\right)',str2];
                    else
                        H_latex_Gamma = ['\left(',str1,'\right)',str2];
                    end
                end
            end
            %
        end
        function [CoeForPauli] = pauliDecompositionNumerial(H_double)
            if ~isequal(   size(H_double) , [2,2])
                error('Pauli Decomposition requires 2*2 Ham!');
            end
            % sigma_0+z sigma_0-z sigma_+ sigma_-
            % 1/2*(sigma_0+sigma_z)
            % 1/2*(sigma_0-sigma_z)
            % 1/2*(sigma_x+1i*sigma_y)
            % 1/2*(sigma_x-1i*sigma_y)
            tmp_mat_r = real(H_double);
            tmp_mat_i = imag(H_double);
            S = 1/2 * [[1;0;0;1],[1;0;0;-1],[0;1;1i;0],[0;1;-1i;0]];
            S = [S,S*1i];
            Phi_8 = (zeros(8,1));
            Phi_8(1) = tmp_mat_r(1,1);
            Phi_8(2) = tmp_mat_r(2,2);
            Phi_8(3) = tmp_mat_r(1,2);
            Phi_8(4) = tmp_mat_r(2,1);
            Phi_8(5) = tmp_mat_i(1,1);
            Phi_8(6) = tmp_mat_i(2,2);
            Phi_8(7) = tmp_mat_i(1,2);
            Phi_8(8) = tmp_mat_i(2,1);
            H_sym_Gamma_L = S * Phi_8;
            count = 0;
            if isa(H_double,'double')
                CoeForPauli = H_sym_Gamma_L;
            else
                for i = 1:4
                    if H_sym_Gamma_L(i)~=(0)
                        count = count+1;
                        H_sym_Gamma = H_sym_Gamma + H_sym_Gamma_L(i)*Gamma_L(i);
                        str1 = latex(H_sym_Gamma_L(i));
                        str2 = latex(Gamma_L(i));
                        if count > 1
                            H_latex_Gamma = [H_latex_Gamma,'+','\left(',str1,'\right)',str2];
                        else
                            H_latex_Gamma = ['\left(',str1,'\right)',str2];
                        end
                    end
                end
            end
        end
        function [H_sym_pauli,H_sym_pauli_L,H_latex_pauli] = pauliDecomposition(H_sym)
            if isa(H_sym,'vasplib') || isa(H_sym,'HR') || isa(H_sym,'Htrig') || isa(H_sym,'HK')
                H_sym = H_sym.sym();
            elseif isa(H_sym,'double')
                H_sym = sym(H_sym);
            elseif isa(H_sym,'sym')
            else
                try
                    H_sym = sym(H_sym);
                catch ME
                    error('cant convert into symbolic');
                end
            end
            if ~isequal(   size(H_sym) , [4,4]) && ~isequal(   size(H_sym) , [2,2])
                error('Pauli Decomposition requires 4*4/2*2 Ham!');
            end
            if isequal(   size(H_sym) , [4,4])
                [~,H_sym_pauli_L] = vasplib.GammaDecomposition(H_sym);
                Pauli_L = gamma_matric.pauli_L();
                
                H_sym_pauli = sym(0);
                H_latex_pauli = 'H = ';
                count = 0;
                ChooseL = logical(1:16);
                for i = 1:16
                    if H_sym_pauli_L(i)~=sym(0)
                        count = count+1;
                        H_sym_pauli = H_sym_pauli + H_sym_pauli_L(i)*Pauli_L(i);
                        str1 = latex(H_sym_pauli_L(i));
                        str2 = latex(Pauli_L(i));
                        if count >1
                            H_latex_pauli = [H_latex_pauli,'+','\left(',str1,'\right)',str2];
                        else
                            H_latex_pauli = ['\left(',str1,'\right)',str2];
                        end
                    else
                        ChooseL(i) = false;
                    end
                end
                H_sym_pauli_L = [H_sym_pauli_L;Pauli_L];
                H_sym_pauli_L = H_sym_pauli_L(:,ChooseL);
            end
            if isequal(   size(H_sym) , [2,2])
                tmp_mat_r = real(H_sym);
                tmp_mat_i = imag(H_sym);
                H_sym_L = sym(zeros(1,4));
                H_sym_pauli_L = sym(zeros(1,4));
%                 sigma_L = pauli_matric();
                %
                Smat = [
                    1 0 0 1;
                    1 0 0 -1;
                    0 1 1i 0;
                    0 1 -1i 0];
%                 Smat = sym([1/2,1/2,0,0; ...
%                                 0,0,1,0; ...
%                                 0,0,0,-1; ...
%                                 1/2,-1/2,0,0]);
                Smat_inv = inv(Smat);
                %
                syms sigma_0 sigma_x sigma_y sigma_z real;
                Pauli_L = [sigma_0 sigma_x sigma_y sigma_z];
                H_sym_L(1:2) = diag(H_sym);
                H_sym_L(3) = diag(H_sym,1);
                H_sym_L(4) = diag(H_sym,-1);
                for i = 1:numel(H_sym_L)
                    Label_tmp = find(Smat_inv(i,:));
                    H_sym_pauli_L(i) = H_sym_L(Label_tmp)*Smat_inv(i,Label_tmp)';
                end
%                 H_sym_pauli_L_temp = Smat*H_sym_L';
%                 H_sym_pauli_L = H_sym_pauli_L_temp';
                H_sym_pauli_L = simplify(H_sym_pauli_L);
                H_sym_pauli = sym(0);
                H_latex_pauli = 'H = ';
                count = 0;
                ChooseL = logical(1:4);
                %Pauli_L = sym(Pauli_L);
                for i = 1:4
                    if H_sym_pauli_L(i)~=sym(0)
                        count = count+1;
                        H_sym_pauli = H_sym_pauli + H_sym_pauli_L(i)*sym(Pauli_L(i));
                        str1 = latex(H_sym_pauli_L(i));
                        str2 = latex(sym(Pauli_L(i)));
                        if count >1
                            H_latex_pauli = [H_latex_pauli,'+','\left(',str1,'\right)',str2];
                        else
                            H_latex_pauli = ['\left(',str1,'\right)',str2];
                        end
                    else
                        ChooseL(i) = false;
                    end
                end
                H_sym_pauli_L = [H_sym_pauli_L;Pauli_L];
                H_sym_pauli_L = H_sym_pauli_L(:,ChooseL);
            end
        end
    end
    %% fit tools
    methods(Static)
        function Value = loss_func(parameters,extra_parm,options)
            arguments
                parameters;
                extra_parm double =0;
                options.mode = 'extra';
                options.algorithm  = 'pure_comparison';
                options.DFTBAND = 'EIGENCAR_DFT';
                options.FITobj = 'H_vasplib';
                options.Varlist = 'Varlist';
                options.Varlock = 'Varlock';
                options.namespace = 'base';
                options.extra = 'options_extra';
                options.show = false;
                options.fig = handle([]);
                options.ax = handle([]);
                options.KPOINTS = 'KPOINTS';
            end
            %%%%%%%%%%%%%%
            EIGENCAR_DFT = evalin(options.namespace,options.DFTBAND);
            FITobj = evalin(options.namespace,options.FITobj);
            try
                Varlist = evalin(options.namespace,options.Varlist);
            catch
                try
                    Varlist = FITobj.symvar_list;
                catch
                    Varlist = [];
                end
            end
            %%%%%%%%%%%%%%%
            if isa(parameters,'table')
                fitmethod = 'Bayes';
            elseif isa(parameters,'double')
                if length(parameters) > 1
                    fitmethod = 'NM';
                else
                    fitmethod = 'Single';
                end
            end
            %%%%%%%%%%%%%%%
            switch fitmethod
                case {'Single','NM'}
                    if isa(FITobj,'function_handle')
                        FITobj_n = FITobj(parameters);
                        FITobj_n = FITobj_n <options.KPOINTS;
                    else
                        try
                            FITobj = FITobj.subs(Varlist,parameters);
                        catch
                            Varlist = FITobj.symvar_list;
                            FITobj = FITobj.subs(Varlist,parameters);
                        end
                        FITobj_n = FITobj.Subsall();
                    end
                    EIGENCAR_vasplib = FITobj_n.EIGENCAR_gen('printmode',false);
                case 'Bayes'
                    %%%%%%%%%%%%%%%%%
                    try
                        Varlock = evalin(options.namespace,options.Varlock);
                    catch

                    end
                    %%%%%%%%%%%%%%%%%
                    for i  = 1:length(Varlist)
                        try  % Generate Field Names from Variables  dynamic fieldnames, or sometimes dynamic field names.
                            FITobj = FITobj.subs(Varlist(i),parameters.(string(Varlist(i))));
                        catch
                            FITobj = FITobj.subs(Varlist(i),Varlock(i));
                        end
                    end
                    FITobj_n = FITobj.Subsall();
                    EIGENCAR_vasplib = FITobj_n.EIGENCAR_gen('printmode',false);
                otherwise
                    error('not be implemented');
            end
            Value = vasplib.EIGENCAR_Value(EIGENCAR_DFT,EIGENCAR_vasplib,extra_parm,...
                'mode',options.mode,...
                'algorithm',options.algorithm,...
                'namespace',options.namespace,...
                'extra',options.extra,...
                'show',options.show,...
                'fig',options.fig,...
                'ax',options.ax ...
                );
        end
        function Varbayes = VarBayes(Varlist,VarGuess,VarWidth,VarFix)
            if nargin < 4
                VarFix = sym([]);
            end
            for i = 1:length(Varlist)
                if isnan(VarGuess(i))
                    VarGuess(i) = 0;
                end
                if ismember(Varlist(i),VarFix)
                    Varbayes(i) = optimizableVariable(string(Varlist(i)),...
                        [-VarWidth(i)+VarGuess(i),VarGuess(i)+VarWidth(i)],...
                        'Optimize',true);
                else
                    Varbayes(i) = optimizableVariable(string(Varlist(i)),...
                        [-VarWidth(i)+ VarGuess(i),VarGuess(i)+VarWidth(i)]);
                end
            end
        end
        function ValueTotal = EIGENCAR_Value(EIGENCAR_DFT,EIGENCAR,extra_parm,options)
            arguments
                EIGENCAR_DFT double;
                EIGENCAR double;
                extra_parm double = 0;
                options.mode = 'extra';
                options.namespace = 'base';
                options.extra = 'options_extra'
                options.algorithm  = 'pure_comparison';
                options.show = false;
                options.fig = handle([]);
                options.ax = handle([]);
            end
            %%%%%%%%%%%%
            if  ~strcmp(options.mode,'extra')
                DATA1 = EIGENCAR_DFT;
                DATA2 = EIGENCAR;
                weight_list = ones(2,1);
                Echoose = ':';
            elseif strcmp(options.mode,'extra')
                options_extra = evalin(options.namespace,options.extra);
                NBAND_range_DFT = options_extra.NBAND_range_DFT;
                NBAND_range = options_extra.NBAND_range;
                klist_range = options_extra.klist_range;
                try
                    weight_list = options_extra.weight_list;
                catch
                    weight_list = ones(2,1);
                end
                DATA1 = EIGENCAR_DFT(NBAND_range_DFT,klist_range);
                DATA2 = EIGENCAR(NBAND_range,klist_range);
                if isfield(options_extra, 'E_range')
                    if ~isempty(options_extra.E_range)
                        E_range =  options_extra.E_range;
                    else
                        E_range = [min(DATA1,[],'all'),max(DATA1,[],'all')];
                    end
                else
                    E_range = [min(DATA1,[],'all'),max(DATA1,[],'all')];
                end
                Echoose =DATA1 >= E_range(1)  & DATA1 <= E_range(2);
            end
            %%%%%%%%%%%%%%%%%%%%
            % Alldata
            count = 0;
            % Direct minus
            count = count +1 ;

            Value(count) = sqrt(mean(mean(abs(DATA1(Echoose)-DATA2(Echoose))).^2));
            % Direct diff minus
            count = count +1 ;
            Value(count) = sqrt(mean(mean(abs(diff(DATA1(Echoose).')-diff(DATA2(Echoose).')))));
            if strcmp(options.algorithm,'pure_comparison')
                ValueTotal = sum(Value.*weight_list);
                return;
            else
                Noccu = options_extra.Noccu;
                GapThreshold = options_extra.GapThreshold;
                E_range =  options_extra.E_range;
                highK = options_extra.highK;
                GapL = (abs(EIGENCAR(Noccu +1,:)-EIGENCAR(Noccu ,:)));...
                    minGap = abs((min(EIGENCAR_DFT(Noccu +1,:)-EIGENCAR_DFT(Noccu ,:)) - min(GapL)));
                Value_Gap_Label =min(GapL)> GapThreshold;
                % Node Structure
                Node_punishment = options_extra.Node_punishment;
                Node_insect_punishment = options_extra.Node_insect_punishment;
                [DEGENCAR1,NODEINSECT1] = vasplib.Degeneracy_EIGENCAR(EIGENCAR_DFT(1:Noccu+1,:),highK,GapThreshold/2);
                [DEGENCAR2,NODEINSECT2] = vasplib.Degeneracy_EIGENCAR(EIGENCAR(1:Noccu+1,:),highK,GapThreshold/2);
                Value_Node = Node_punishment*sum(sum(abs(DEGENCAR1-DEGENCAR2))) ...
                    + Node_insect_punishment*sum(sum(abs(NODEINSECT1-NODEINSECT2)));
                % Flash
                IM1 = vasplib.EIGENCAR2IMG(EIGENCAR_DFT,GapThreshold/2,E_range);
                IM2 = vasplib.EIGENCAR2IMG(EIGENCAR,GapThreshold/2,E_range);
                Value_ssimval = 1/ssim(IM1,IM2)-1;
                % High symmetry line energy
                E_highK_DFT = EIGENCAR_DFT(NBAND_range_DFT,highK);
                E_highK_TB = EIGENCAR(NBAND_range_DFT,highK);
                Value_E_HighK = mean(mean(abs(E_highK_DFT-E_highK_TB)));
            end
            if strcmp(options.algorithm,'dirac')
                Dirac_area = options_extra.Dirac_area;
                Dirac_punishment = options_extra.Dirac_punishment;
                Gapless_points = GapL<0.3;
                Value_Dirac = 0;
                for i = 1:length(Gapless_points )
                    if Gapless_points(i) == 1 && ~ismember(i,Dirac_area)
                        Value_Dirac = Value_Dirac+ Dirac_punishment;
                    end
                end
            end
            switch options.algorithm
                case 'dirac'
                    count = count +1 ;
                    Value(count) = Value_E_HighK;
                    ValueTotal = (1-extra_parm)*(Value_ssimval^2)*sum(Value.*weight_list(1:3));
                    %
                    ValueTotal = ValueTotal...
                        +extra_parm * Value_Gap_Label * weight_list(4)...
                        +extra_parm *(Value_ssimval^3)  *(Value_Node+1) * weight_list(5)...
                        +extra_parm * Value_Dirac * weight_list(6);
                case 'insulator'
                    Gapless_punishment = options_extra.Gapless_punishment;
                    count = count +1 ;
                    Value(count) = Value_E_HighK;
                    ValueTotal = (1-extra_parm) * (Value_ssimval^2) * sum(Value.*weight_list(1:3));
                    %
                    ValueTotal = ValueTotal...
                        +extra_parm * minGap * weight_list(4)...
                        +extra_parm * (Value_ssimval^3) * (Value_Node+1) * weight_list(5)...
                        +extra_parm * Gapless_punishment * ~Value_Gap_Label * weight_list(6);
                case 'metal'
                    Gap_punishment = options_extra.Gap_punishment;
                    count = count +1 ;
                    Value(count) = Value_E_HighK;
                    ValueTotal = (1-extra_parm)*(Value_ssimval^2)*sum(Value.*weight_list(1:3));
                    %
                    ValueTotal = ValueTotal...
                        +extra_parm * maxGap * weight_list(4)...
                        +extra_parm * (Value_ssimval^3) * (Value_Node+1) * weight_list(5)...
                        +extra_parm * Gap_punishment * Value_Gap_Label * weight_list(6);
                otherwise
            end

            %ValueTotal = ValueTotal^2;

        end
        function options_extra = FitOptionHelper(EIGENCAR_DFT,algorithm,options)
            arguments
                EIGENCAR_DFT;
                algorithm {mustBeMember(algorithm,['pure_comparison','dirac','insulator','metal'])} = pure_comparison ;
                options.show logical = false;
                options.WAN_NUM {mustBeInteger} = -1;
                options.Noccu {mustBeInteger} = 1;
                options.klist_range  = ':';
                options.NBAND_range_DFT = 1:size(EIGENCAR_DFT,1)/2;
                options.NBAND_range = [];
                options.weight_list = [];
                options.KPOINTS = 'KPOINTS';
                options.highK = [1,size(EIGENCAR_DFT,2)];
                options.E_range = [];
                options.GapThreshold = 0.1;
                options.Node_punishment = 3;
                options.Node_insect_punishment = 0.5;
                options.Dirac_area = [];
                options.Dirac_width = 3;
                options.Dirac_punishment = 5;
                options.Gap_punishment = 5;
                options.Gapless_punishment = 5;
            end
            if options.WAN_NUM  == -1
                options_extra.WAN_NUM = options.Noccu*2;
            else
                options_extra.WAN_NUM = options.WAN_NUM;
            end
            %
            if isempty(options.NBAND_range )
                options_extra.NBAND_range = 1:options.Noccu;
            else
                options_extra.NBAND_range = options.NBAND_range;
            end
            %
            if isempty(options.E_range)
                options_extra.E_range  = [min(min(EIGENCAR_DFT(options_extra.NBAND_range,:))),...
                    max(max(EIGENCAR_DFT(options_extra.NBAND_range,:)))];
            else
                options_extra.E_range  = options.E_range;
            end
            %
            %设定拟合范围：
            options_extra.NBAND_range_DFT = options.NBAND_range_DFT;
            options_extra.klist_range = options.klist_range;
            %
            options_extra.Noccu = options.Noccu;
            options_extra.GapThreshold = options.GapThreshold;
            options_extra.highK = options.highK;
            %
            options_extra.Node_punishment = options.Node_punishment;
            options_extra.Node_insect_punishment = options.Node_insect_punishment;
            options_extra.Dirac_punishment = options.Dirac_punishment;
            options_extra.Gap_punishment = options.Gap_punishment;
            options_extra.Gapless_punishment = options.Gapless_punishment;
            %
            if strcmp(algorithm,'dirac')
                if isempty(options.Dirac_area)
                    Dirac_area = find(abs(EIGENCAR_DFT(options_extra.Noccu+1,:) - EIGENCAR_DFT(options_extra.Noccu,:))...
                        <options_extra.GapThreshold);
                    DiracA = [];
                    for i = Dirac_area
                        DiracA = [DiracA,i-options.Dirac_width : i+options.Dirac_width];
                    end
                    options_extra.Dirac_area = unique(DiracA);
                else
                    options_extra.Dirac_area = options.Dirac_area;
                end
            end
            %默认大小和斜率等权
            if isempty(options.weight_list)
                switch algorithm
                    case 'pure_comparison'
                        options_extra.weight_list = [1,1];
                    case 'dirac'
                        options_extra.weight_list = [1,1,1,1,1,1];
                    case 'insulator'
                        options_extra.weight_list = [1,1,1,1,1,1];
                    case 'metal'
                        options_extra.weight_list = [1,1,1,1,1,1];
                    otherwise
                        options_extra.weight_list = [1,1,1,1,1,1];
                end
            else
                options_extra.weight_list = options.weight_list;
                % disp([';;']);
            end
        end
        function IMG = EIGENCAR2IMG(EIGENCAR,dE,Erange,kselect)

            % nargin
            if nargin < 2
                dE  = 0.1;
            end
            if nargin <3
                Erange = [min(min(EIGENCAR)),max(max(EIGENCAR))];
            end
            if nargin <4
                kselect = 1:size(EIGENCAR,2);
            end
            NBANDS = size(EIGENCAR,1);
            DATA = EIGENCAR(:,kselect);

            X_nodes = size(DATA,2);
            Y_nodes = round(abs(Erange(2)-Erange(1))/dE);
            Emin = min(Erange);
            IMG_sparse = sparse(Y_nodes,X_nodes);
            for i =1:X_nodes
                for j = 1:NBANDS
                    NE = ceil((DATA(j,i)-Emin)/dE);
                    if NE>0 && NE <= Y_nodes
                        IMG_sparse(Y_nodes-NE+1,i) = IMG_sparse(Y_nodes-NE+1,i) +0.1;
                    end
                end
            end
            %UNTITLED2 此处显示有关此函数的摘要
            %   此处显示详细说明
            IMG = full(IMG_sparse);
        end
        function [DEGENCAR, NODEINSECT]= Degeneracy_EIGENCAR(EIGENCAR,highK,dE)
            if nargin <3
                dE =0.1;
            end
            [NBANDS,~] = size(EIGENCAR);
            NhighK = length(highK);
            DEGENCAR = zeros(NBANDS,NhighK);
            NODEINSECT = zeros(NBANDS-1,NhighK-1);
            % NhighK
            for i = highK
                DEGENCAR(1,i) = 1;
                count = 1;
                for j = 2:NBANDS
                    if (EIGENCAR(j,i)-EIGENCAR(j-1,i)) < dE
                        DEGENCAR(count,i) = DEGENCAR(count,i)+1;
                    else
                        count = count+1;
                        DEGENCAR(count,i) = 1;
                    end
                end
            end
            %
            for i = 1:NBANDS-1
                for j = 1:NhighK-1
                    NODEINSECT(i,j) = sum(abs(EIGENCAR(i+1,highK(j)+1:highK(j+1)-1)-...
                        EIGENCAR(i,highK(j)+1:highK(j+1)-1)) < dE);
                end
            end
        end

    end
    %% tools 大杂烩
    methods(Static,Hidden)
        function [ChirdrenCell,Type] = fixedchildren(SymVar,mode)
            arguments
                SymVar sym;
                mode = 'exp_inner';
            end
            if strcmp(mode,'exp_inner')
                subexpr = children(combine(SymVar,'exp'));
                if isequal(simplify(fold(@plus,subexpr)-SymVar),sym(0))
                    ChirdrenCell = subexpr;
                    Type = 'sum';
                    return;
                end
                if isequal(simplify(fold(@times,subexpr)-SymVar),sym(0))
                    ChirdrenCell{1} = SymVar;
                    Type = 'prod';
                    return;
                end
                ChirdrenCell = subexpr;
                Type = 'inner';
            else
                subexpr = children(SymVar);
                if isequal(simplify(fold(@plus,subexpr)-SymVar),sym(0))
                    ChirdrenCell = subexpr;
                    Type = 'sum';
                    return;
                end
                if isequal(simplify(fold(@times,subexpr)-SymVar),sym(0))
                    ChirdrenCell{1} = SymVar;
                    Type = 'prod';
                    return;
                end
                ChirdrenCell{1} = SymVar;
                Type = 'inner';
            end
        end
        function [Asort,Usort] = sorteig(U,A)
            if nargin <2
                mode = 'eigenval';
                NUM_WAN = length(U);
            else
                mode = 'whole';
                NUM_WAN = length(A);
            end
            NBANDS = length(U);
            if ~isvector(U)
                SortTmp=diag(U);%抽取特征值
                vec = false;
            else
                SortTmp = U;
                vec = true;
            end
            if strcmp(mode,'whole')
                % 按从大到小的特征值顺序排序重新组合对应的特征向量
                Asort=zeros(NUM_WAN ,NBANDS );
                [Usort,IJ]=sort(SortTmp,1,'ComparisonMethod','real');
                for jj=1:NBANDS
                    Asort(:,jj)=A(:,IJ(jj));%取特征向量的列向量
                end
                if ~vec 
                    Usort = diag(Usort);
                end
            elseif strcmp(mode,'eigenval')
                % 按从大到小的特征值顺序排序重新组合对应的特征向量
                SortTmp=diag(U);%抽取特征值
                [Usort,~]=sort(SortTmp,1,'ComparisonMethod','real');
                if ~vec
                    Usort = diag(Usort);
                end
                Asort = [];
            end
        end
        function SymVar = shelling(SymVar)
            % cos sin exp
            StrVar = char(SymVar);
            StrVar(end) = [];
            StrVar(1:4) = [];
            SymVar = str2sym(StrVar);
        end
        function R = Sph2Cart(S)
            % R =S;
            r= S(:,1);
            theta = S(:,2);
            phi = S(:,3);
            z = r .*cos(theta);
            rcoselev = r .* sin(theta);
            x = rcoselev .* cos(phi);
            y = rcoselev .* sin(phi);
            R = [x,y,z];
        end
    end
    methods (Static,Hidden,Access= protected)
        function C = isclose(A,B)
            %             if (A-B).^2<1e-6
            %                 C = true;
            %             else
            %                 C = false;
            %             end
            C = (A-B).^2<1e-12 ;
        end
        function C = allclose(A,B)
            if all(all(Oper.isclose(A,B)))
                C = true;
            else
                C = false;
            end
        end
        function TrueOrFalse = LatticeVectorTest(Vector,Accuracy)
            if nargin < 2
                Accuracy = 1e-6;
            end
            TrueOrFalse = true;
            for i = 1:numel(Vector)
                if abs(rem(Vector(i),1)) > Accuracy
                    TrueOrFalse = false;
                    return;
                end
            end
        end
        function angle = name_angle(theta, Latex)
            arguments
                theta double ;
                Latex logical = false;
            end
            frac = rat(theta / pi, 1e-4);
            frac = simplify(str2sym(frac))*pi;
            if Latex
                angle = latex(frac);
            else
                angle = string(frac);
            end
        end
        function symPage = subspage(symPage,lhs,rhs)
            for i = 1:numel(lhs)
                symPage = subs(symPage,lhs(i),rhs(i));
            end
        end
        function [nn_sparse_temp,Rnn_list] = nn_sparse_gen(orb1,orb2,Rm,search_rangex,search_rangey,search_rangez,Accuracy,Rlength_cut,onsite)
            arguments
                orb1 ;
                orb2 ;
                Rm   ;
                search_rangex;
                search_rangey;
                search_rangez;
                Accuracy;
                Rlength_cut;
                onsite = false;
            end
            Rc1 = orb1;
            Rc2 = orb2;
            R_fractional_diff = -(Rc1 - Rc2);
            %     nn_smart_t.R_cartesian_to = Atom_smart_t.R_fractional_to*Rm;
            %     nn_smart_t.R_cartesian_from = Atom_smart_t.R_fractional_from*Rm;
            count = 1;
            reducible_num=(2*search_rangex+1)*(2*search_rangey+1)*(2*search_rangez+1);
            Rnn_list = zeros(reducible_num,1);
            if Accuracy == -1
                nn_sparse_temp= sym(zeros(reducible_num,10));
                symmode = true;
            else
                nn_sparse_temp= zeros(reducible_num,10);
                symmode = false;
            end
            % nn_t = struct('R_vector',[],'R_fractional_diff',[],'Rlength',[],'nn_level',[],'hop_pre',[],'hop',[]);
            % nn = repmat(nn_t,[reducible_num 1]);
            for Rf_a1=-search_rangex:search_rangex
                for Rf_a2=-search_rangey:search_rangey
                    for Rf_a3=-search_rangez:search_rangez
                        R_vector = [Rf_a1 Rf_a2 Rf_a3];
                        Rij_cart = (R_vector + R_fractional_diff)*Rm ;
                        Rlength = norm(Rij_cart);
                        if ~symmode
                            Rlength = roundn(Rlength,round(log(Accuracy)/log(10)));
                            %Rlength = vpa
                            if  ((0 < Rlength)||(0 == Rlength && onsite)) && Rlength < Rlength_cut
                                nn_sparse_temp(count,6) = Rf_a1;
                                nn_sparse_temp(count,7) = Rf_a2;
                                nn_sparse_temp(count,8) = Rf_a3;
                                nn_sparse_temp(count,3) = Rij_cart(1);
                                nn_sparse_temp(count,4) = Rij_cart(2);
                                nn_sparse_temp(count,5) = Rij_cart(3);
                                nn_sparse_temp(count,9) =  Rlength;
                                Rnn_list(count,:) = Rlength;
                                count = count +1;
                            end
                        else
                            if  (sym(0) < Rlength||(sym(0) == Rlength && onsite)) && Rlength < sym(Rlength_cut)
                                nn_sparse_temp(count,6) = Rf_a1;
                                nn_sparse_temp(count,7) = Rf_a2;
                                nn_sparse_temp(count,8) = Rf_a3;
                                nn_sparse_temp(count,3) = Rij_cart(1);
                                nn_sparse_temp(count,4) = Rij_cart(2);
                                nn_sparse_temp(count,5) = Rij_cart(3);
                                nn_sparse_temp(count,9) =  Rlength;
                                Rnn_list(count,:) = Rlength;
                                count = count +1;
                            end
                        end
                    end
                end
            end
            if count <= reducible_num
                Rnn_list(count:reducible_num,:) = [];
                nn_sparse_temp(count:reducible_num,:) = [];
            end
        end
        function labelcut_list = labelcut_list_gen(Atom_num)
            n = length(Atom_num);
            beginline =8 ;
            sum_n =beginline;
            for i =1:n
                sum_n2 = sum_n+Atom_num(i)-1;
                labelcut_list(i,:) = [sum_n sum_n2];
                sum_n = sum_n2+1;
            end

        end
        function l_num = orb2l(input)
            input = string(input);
            switch input
                case {'0','s','S'}
                    l_num = 0;
                case {'1','p','P','px','py','pz','Px','Py','Pz','p_x','p_y','p_z'}
                    l_num = 1;
                case {'2','d','D','dx2-y2','dz2','dxy','dyx','dxz','dyz','dzx','dzy','dx2y2'}
                    l_num = 2;
                case {'3','f'}
                    l_num = 3;
                case {'-1','sp'}
                    l_num = -1;
                case {'-2','sp2'}
                    l_num = -2;
                case {'-3','sp3'}
                    l_num = -3;
                case {'-4','sp3d'}
                    l_num = -4;
                case {'-5','sp3d2'}
                    l_num = -5;
                otherwise
                    l_num = 1i;
            end


        end
        function m_num = orb_sym2m(input)
            input = string(input);
            switch input
                case {'I','s','S','z','pz','Pz','p_z','z^2','dz^2','d_z^2','z^3','fz_3','0'}
                    m_num = 0;
                case {'-1','y','py','p_y','P_y','yz','zy','dyz','d_yz','dzy','d_zy',...
                        'yz^2','z^2y','f_yz^2','fyz^2','f_z^2y','fz^2y'}
                    m_num = -1;
                case {'1','x','px','p_x','P_x','xz','zx','dxz','d_xz','dzx','d_zx',...
                        'xz^2','z^2x','f_xz^2','fxz^2','f_z^2x','fz^2x'}
                    m_num = 1;
                case {'-2','xy','yx','dxy','d_xy','dyx','d_yx',...
                        'xyz','zxy','yzx','xzy','zyx','yxz',...
                        'fxyz','fzxy','fyzx','fxzy','fzyx','fyxz',...
                        'f_xyz','f_zxy','f_yzx','f_xzy','f_zyx','f_yxz'}
                    m_num = -2;
                case {'2','(x^2-y^2)','-(y^2-x^2)','d(x^2-y^2)','d_(x^2-y^2)','x^2-y^2',...
                        'dx2y2','d_x2y2','dx2-y2','d_x2-y2',...
                        '(x^2-y^2)z','z(x^2-y^2)','f_(x^2-y^2)z','f_z(x^2-y^2)',...
                        'f(x^2-y^2)z','fz(x^2-y^2)','fzx2y2','f_zx2y2','fx2y2z','f_x2y2z'}
                    m_num = 2;
                case {'3','y(3x^2-y2)','(3x^2-y^2)','y3x2-y2','3x2-y2y','y3x2y2',...
                        'fy(3x^2-y2)','f(3x^2-y^2)','fy3x2-y2','f3x2-y2y','fy3x2y2',...
                        'f_y(3x^2-y2)','f_(3x^2-y^2)','f_y3x2-y2','f_3x2-y2y','f_y3x2y2',...
                        }
                    m_num = 3;
                case {'-3','x(3y^2-x2)','(3y^2-x^2)','x3y2-x2','3y2-x2x','x3y2x2',...
                        'fx(3y^2-x2)','f(3y^2-x^2)','fx3y2-x2','f3y2-x2x','fx3y2x2',...
                        'f_x(3y^2-x2)','f_(3y^2-x^2)','f_x3y2-x2','f_3y2-x2x','f_x3y2x2',...
                        }
                    m_num = -3;
                case {'4'}
                    m_num = 4;
                case {'-4'}
                    m_num = -4;
                case {'5'}
                    m_num = 5;
                case {'-5'}
                    m_num = -5;
                otherwise
                    m_num = 1i;
            end


        end
        function orb_sym = Ymlsym(l,m,orb_symbk)
            if l == 1i && m == 1i
                orb_sym = orb_symbk;
                return;
            end
            switch l
                case 0
                    switch m
                        case 0
                            orb_sym = str2sym('1');
                        otherwise
                            warning('check your input POSCAR');
                    end
                case 1
                    switch m
                        case 0
                            orb_sym = str2sym('z');
                        case -1
                            orb_sym = str2sym('y');
                        case 1
                            orb_sym = str2sym('x');
                        otherwise
                            warning('check your input POSCAR');
                    end
                case 2
                    switch m
                        case 0
                            orb_sym = str2sym('z^2');
                        case -1
                            orb_sym = str2sym('y*z');
                        case 1
                            orb_sym = str2sym('x*z');
                        case -2
                            orb_sym = str2sym('x*y');
                        case 2
                            orb_sym = str2sym('x^2-y^2');
                        otherwise
                            warning('check your input POSCAR');
                    end
                case 3
                    switch m
                        case 0
                            orb_sym = str2sym('z^3');
                        case -1
                            orb_sym = str2sym('y*z^2');
                        case 1
                            orb_sym = str2sym('x*z^2');
                        case -2
                            orb_sym = str2sym('x*y*z');
                        case 2
                            orb_sym = str2sym('(x^2-y^2)*z');
                        case -3
                            orb_sym = str2sym('y*(3*x^2-y^2)');
                        case 3
                            orb_sym = str2sym('x*(x^2-3*y^2)');
                        otherwise
                            warning('check your input POSCAR');
                    end
                case -1
                    switch m
                        case 1
                            orb_sym = str2sym('1+x');
                        case 2
                            orb_sym = str2sym('1-x');
                        otherwise
                            warning('check your input POSCAR');
                    end
                case -2
                    switch m
                        case 1
                            orb_sym = str2sym('3^(-1/2)-6^(-1/2)*x+2^(-1/2)*y');
                        case 2
                            orb_sym = str2sym('3^(-1/2)-6^(-1/2)*x-2^(-1/2)*y');
                        case 3
                            orb_sym = str2sym('3^(-1/2)+6^(-1/2)*x+6^(-1/2)*x');
                        otherwise
                            warning('check your input POSCAR');
                    end
                case -3
                    switch m
                        case 1
                            orb_sym = str2sym('1+x+y+z');
                        case 2
                            orb_sym = str2sym('1+x-y-z');
                        case 3
                            orb_sym = str2sym('1-x+y-z');
                        case 4
                            orb_sym = str2sym('1-x-y+z');
                        otherwise
                            warning('check your input POSCAR');
                    end
                case -4 %
                    switch m
                        case 1
                            orb_sym = str2sym('3^(-1/2)-6^(-1/2)*x+2^(-1/2)*y');
                        case 2
                            orb_sym = str2sym('3^(-1/2)-6^(-1/2)*x-2^(-1/2)*y');
                        case 3
                            orb_sym = str2sym('3^(-1/2)+6^(-1/2)*x+6^(-1/2)*x');
                        case 4
                            orb_sym = str2sym('2^(-1/2)*z+2^(-1/2)*z^2');
                        case 5
                            orb_sym = str2sym('-2^(-1/2)*z+2^(-1/2)*z^2');
                        otherwise
                            warning('check your input POSCAR');
                    end
                case -5
                    switch m
                        case 1
                            orb_sym = str2sym('6^(-1/2)-2^(-1/2)*x-12^(-1/2)*z^2+2^(-1)*(x^2-y^2)');
                        case 2
                            orb_sym = str2sym('6^(-1/2)+2^(-1/2)*x-12^(-1/2)*z^2+2^(-1)*(x^2-y^2)');
                        case 3
                            orb_sym = str2sym('6^(-1/2)-2^(-1/2)*x-12^(-1/2)*z^2-2^(-1)*(x^2-y^2)');
                        case 4
                            orb_sym = str2sym('6^(-1/2)+2^(-1/2)*x-12^(-1/2)*z^2-2^(-1)*(x^2-y^2)');
                        case 5
                            orb_sym = str2sym('6^(-1/2)-2^(-1/2)*z+3^(-1)*(z^2)');
                        case 6
                            orb_sym = str2sym('6^(-1/2)+2^(-1/2)*z+3^(-1)*(z^2)');
                        otherwise
                            warning('check your input POSCAR');
                    end
                otherwise
                    warning('check your input POSCAR');
            end
        end
        function orbsym_n = subs_xyz(orbsym,Rlmn)
            % subs ok
            %disp(Rlmn);
            if orbsym == sym(1)
                orbsym_n = 1;
            elseif orbsym == sym('x')
                orbsym_n = Rlmn(1);

            elseif orbsym == sym('y')
                orbsym_n = Rlmn(2);

            elseif orbsym == sym('z')
                orbsym_n = Rlmn(3);

            else
                orbsym_n = 1;
            end
            % l m n fixed ?
            % orbsym_n =abs(orbsym_n);
            %disp(orbsym_n);
        end
        function out = delta_orb(orb1,orb2)
            if strcmp(orb1, orb2)
                out=1;
            else
                out =0;
            end
        end
        function out = delta_orb_sym(orb_sym1,orb_sym2)
            if orb_sym1 == orb_sym2
                out=1;
            else
                out =0;
            end
        end
        function R_struct = Rm2abc(Rm)
            R_struct.a = norm(Rm(1,:));
            R_struct.b = norm(Rm(2,:));
            R_struct.c = norm(Rm(3,:));
            R_struct.alpha = acos(dot(Rm(2,:),Rm(3,:))/(norm(R_struct.b)*norm(R_struct.c)))/pi*180;
            R_struct.beta  = acos(dot(Rm(3,:),Rm(1,:))/(norm(R_struct.c)*norm(R_struct.a)))/pi*180;
            R_struct.gamma = acos(dot(Rm(1,:),Rm(2,:))/(norm(R_struct.a)*norm(R_struct.b)))/pi*180;
        end
        function flag = strcontain(str,containlist)
            %% Retun 1 or 0 to show if there a string contains any string in string list
            N=length(containlist);
            flagtemp=0;
            for i =1:N
                temp_add =  contains(str,string(containlist(i)));
                if strcmp(containlist(i),"")
                    temp_add  = 1;
                end
                flagtemp = temp_add+flagtemp;
            end

            flag =flagtemp;
        end
        function rc_plus = plusrc(rc)
            if rc<0
                %disp('bugtest');
                rc_plus = rc +1;

            elseif rc>=1
                rc_plus =mod(rc,1);
            else
                rc_plus =rc;
            end
        end
        function To_red_sc=to_red_sc(red_vec_orig ,Ns)
            To_red_sc=red_vec_orig/Ns;
            %To_red_sc= To_red_sc';
        end
        function To_red_pc=to_red_pc(red_vec_sc ,Ns)
            To_red_pc = red_vec_sc*Ns;
        end
        function [orb_one_incell,translation_vector]=translation_orb(orb_one)
            orb_one_incell = mod(orb_one,1);
            if orb_one_incell(1) ==1
                orb_one_incell(1)=0;
            end
            if orb_one_incell(2) ==1
                orb_one_incell(2)=0;
            end
            if orb_one_incell(3) ==1
                orb_one_incell(3)=0;
            end
            translation_vector = round(orb_one_incell - orb_one);
        end
        function [list_obj_unique,sorted_label,cut_slice] = cut_tools(list_obj)
            [list_obj_sorted,sorted_label] = sortrows(list_obj);
            [list_obj_unique,unique_label] = unique(list_obj_sorted,'rows');
            cut_slice = [(unique_label),[(unique_label(2:end)-1);size(list_obj_sorted,1)]];
        end
    end
end

