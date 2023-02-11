function [t_list,topo_of_Ham, topo] = UnsupClass(Ham_obj,dim,t_range,t_num,opts)
% unsupervised classification of topological gapped system 
% 10.1103/PhysRevLett.130.036601
arguments
    Ham_obj
    dim int16;
    t_range (:,2) double = [0, 2; 0, 2]
    t_num (1,:) double = 20
    opts.perturbation = []
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
%% mesh of parameters
t_list = para_mesh_gen(t_range , t_num);
nsamples = length(t_list);
%% init
tic
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
    if ~is_included && ~isempty(opts.perturbation)
        H2p = H2 + opts.perturbation;
        for class_i = 1:nclass
            H1p = topo{class_i} + opts.perturbation;       
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
toc
end