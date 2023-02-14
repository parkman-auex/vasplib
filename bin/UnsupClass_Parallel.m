function [t_list, topo_of_Ham, topo] = UnsupClass_Parallel(Ham_obj,dim,t_range,t_num,opts,prec_opts)
% unsupervised classification of topological gapped system 
% 10.1103/PhysRevLett.130.036601
arguments
    Ham_obj
    dim int16;
    t_range (:,2) double = [0, 2; 0, 2]
    t_num (1,:) double = 20
    % opts.perturbation = []
    opts.ncore int16 = 2
    prec_opts.noccu int16 = 0
    prec_opts.metal_threshold double = 1e-4
    prec_opts.exist_metal_phase logical = false
end
%% check inputs
symvarL = Ham_obj.symvar_list;
t_type = length(symvarL);
if t_type ~= size(t_range,1)
    error("The inputs are not compatible! Please check the Ham_obj.symvar_list and the t_range")
end
switch length(t_num)
    case 1
        t_num = ones(1,t_type).*t_num;
    case t_type
        
    otherwise
        error("The inputs are not compatible! Please check the t_num and the t_range")
end
%% half occupation
nbands = Ham_obj.Nbands;
if prec_opts.noccu == 0
    prec_opts.noccu = nbands/2;
end
%% mesh of parameters
t_list = para_mesh_gen(t_range, t_num, 'shift', 'uni-rand');
nsamples = length(t_list);
nparts = nsamples/opts.ncore;
if round(nparts) ~= nparts
    error("nsamples="+nsamples+" can not be divided evenly by opts.ncore="+ opts.ncore)
end
%% allocate memory space
t_list_cell = cell(opts.ncore,1);
topo_of_Ham_cell = cell(opts.ncore,1);
topo_cell = cell(opts.ncore,1);
%% slice the t_list to avoid broadcast
for ci = 1:opts.ncore
    t_list_cell{ci} = t_list((ci-1)*nparts+1:ci*nparts,:);
end
%% start parallel pool
poolobj = parpool(opts.ncore);
tic
parfor (ci = 1:opts.ncore,opts.ncore)
    [topo_of_Ham_cell{ci}, topo_cell{ci}] = UnsupClass_inner(Ham_obj, dim, t_list_cell{ci},...
        prec_opts.noccu, prec_opts.metal_threshold, prec_opts.exist_metal_phase);
end
toc
%% end parallel pool
delete(poolobj)
%% clustering method
topo = topo_cell{1};
topo_of_Ham = topo_of_Ham_cell{1};

for ci = 2:opts.ncore    
    topoB = topo_cell{ci};   
    topo_of_HamB = topo_of_Ham_cell{ci};
    
    ai_skip = [];
    for bi = 1:length(topoB)        
        LA = length(topo);
        is_included = false;
        
        for ai = 1:LA
            if ismember(ai, ai_skip)
                continue
            end
            if isTopoEquv(topoB{bi}, topo{ai}, dim, prec_opts.noccu, prec_opts.metal_threshold)
                is_included = true;
                break
            end
        end

        if is_included
            topo_of_HamB = subs(topo_of_HamB,bi,-ai);
            ai_skip = [ai_skip, ai];
        else       
            topo{LA+1} = topoB{bi};
            topo_of_HamB = subs(topo_of_HamB,bi,-(LA+1));
            ai_skip = [ai_skip, LA+1];
        end
    end
    topo_of_Ham = [topo_of_Ham, -int16(topo_of_HamB)];
end
end

function [topo_of_Ham, topo] = UnsupClass_inner(Ham_obj, dim, t_list,...
    noccu, metal_threshold, exist_metal_phase)
% unsupervised classification of topological gapped system 
% 10.1103/PhysRevLett.130.036601
arguments
    Ham_obj
    dim int16;
    t_list double
    noccu int16
    metal_threshold double = 1e-4
    exist_metal_phase logical = false
    % perturbation = []
end
%% check inputs
symvarL = Ham_obj.symvar_list;
nsamples = length(t_list);
%% init
topo = {};
topo(1) = {Ham_obj.subs(symvarL, t_list(1,:))};
topo_of_Ham = zeros(1,nsamples);
topo_of_Ham(1) = 1;
%% loop
for sample_i = 2:nsamples
    H2 = Ham_obj.subs(symvarL, t_list(sample_i,:));
    
    %% exclude metal phase and set 'topo_of_Ham' to be zeros
    if exist_metal_phase
        if isMetal_inner(H2, dim, noccu, metal_threshold)
            % topo_of_Ham(sample_i) = 0;            
            break
        end               
    end
    
    %%
    nclass = length(topo);
    is_included = false;  
    for class_i = 1:nclass
        if isTopoEquv(topo{class_i}, H2, dim, noccu, metal_threshold)
            is_included = true;            
            break
        end
    end
    
%     %% add pertubations to overcome accidental crossing point
%     if ~is_included && ~isempty(perturbation)
%         H2p = H2 + perturbation;
%         for class_i = 1:nclass
%             H1p = topo{class_i} + perturbation;       
%             if isTopoEquv(H1p,H2p,dim,noccu)
%                 is_included = true;            
%                 break
%             end
%         end
%     end          
    
    if is_included
        topo_of_Ham(sample_i) = class_i;
    else       
        topo(nclass+1) = {H2};
        topo_of_Ham(sample_i) = nclass+1;
    end
end
end

function ismetal = isMetal_inner(Ham_obj, dim, noccu, metal_threshold)
% compare two topological GAPPED systems
arguments
    Ham_obj
    dim int16
    noccu int16
    metal_threshold double
end
%% precision control
%% very important !!! nk<=4 may miss the gapless points at rare cases
nk = 5;
options = optimset('TolFun',metal_threshold/10,'TolX',metal_threshold/10,'Display','off');
%% generate initial k-mesh
switch dim
    case 1
        [klist_s,klist_r] = kmesh_gen(Ham_obj,[],'nk',[nk, 1, 1]);        
    case 2
        [klist_s,klist_r] = kmesh_gen(Ham_obj,[],'nk',[nk, nk, 1]);
    case 3
        [klist_s,klist_r] = kmesh_gen(Ham_obj,[],'nk',[nk, nk, nk]);
end
nkpts = size(klist_s,1);
%% handle of band gap search
get_gap_kpoint = @(kpoint) get_gap_inner(Ham_obj, kpoint, noccu);
%%
ismetal = false;
switch class(Ham_obj)
    case "HR"
        for i = 1:nkpts
            [~,gapout]= fminsearch(get_gap_kpoint,klist_s(i,:),options);
            if gapout < metal_threshold
                ismetal = true;
                break
            end
        end
    case {"Htrig","HK"}
        for i = 1:nkpts
            [~,gapout]= fminsearch(get_gap_kpoint,klist_r(i,:),options);
            if gapout < metal_threshold
                ismetal = true;
                break
            end
        end
end
end

function gap = get_gap_inner(Ham_obj, kpoint, noccu)
EIGENCAR = Ham_obj.EIGENCAR_gen('klist',kpoint,'printmode',false);
gap = EIGENCAR(noccu+1,:) - EIGENCAR(noccu,:);
end