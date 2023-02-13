function [t_list, topo_of_Ham, topo] = UnsupClass_Parallel(Ham_obj,dim,t_range,t_num,opts)
% unsupervised classification of topological gapped system 
% 10.1103/PhysRevLett.130.036601
arguments
    Ham_obj
    dim int16;
    t_range (:,2) double = [0, 2; 0, 2]
    t_num (1,:) double = 20
    opts.perturbation = []
    opts.ncore = 10
end
%% check inputs
dH = opts.perturbation;
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
    [topo_of_Ham_cell{ci}, topo_cell{ci}] = UnsupClass_single(...
        Ham_obj, dim, t_list_cell{ci}, dH);
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
            if isTopoEquv(topoB{bi},topo{ai},dim)
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

function [topo_of_Ham, topo] = UnsupClass_single(Ham_obj,dim,t_list, perturbation)
% unsupervised classification of topological gapped system 
% 10.1103/PhysRevLett.130.036601
arguments
    Ham_obj
    dim int16;
    t_list double
    perturbation = []
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
    
    nclass = length(topo);
    is_included = false;  
    for class_i = 1:nclass
        if isTopoEquv(topo{class_i},H2,dim)
            is_included = true;            
            break
        end
    end
    
    %% add pertubations to overcome accidental crossing point
    if ~is_included && ~isempty(perturbation)
        H2p = H2 + perturbation;
        for class_i = 1:nclass
            H1p = topo{class_i} + perturbation;       
            if isTopoEquv(H1p,H2p,dim,noccu)
                is_included = true;            
                break
            end
        end
    end          
    
    if is_included
        topo_of_Ham(sample_i) = class_i;
    else       
        topo(nclass+1) = {H2};
        topo_of_Ham(sample_i) = nclass+1;
    end
end
end