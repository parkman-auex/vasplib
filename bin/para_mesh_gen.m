%para_mesh_gen generate a mesh in form of list, in arbitrary
% dimensional complex parameter space。
% It can be seen as an encapsulation of ndgrid 
%
% mesh_list = para_mesh_gen(AtoB, np)
% AtoB = [p1A, p1B; p2A, p2B; ... ]
%   the value of parameters ranges from A to B
% np = [np1, np2, ... ] or np
%   control the mesh density of each parameter by inputing a full list,
%   or use the same density by inputing a single number

function mesh_list = para_mesh_gen(AtoB, np, opts)
arguments
    AtoB (:,2) double = [0, 1; 0, 1; 0, 1]
    np double = 10
    opts.mode {mustBeMember(opts.mode,{'uniform','random'})} = 'uniform';
end
%% check
nparas = size(AtoB,1);
switch length(np)
    case 1
        np_list = ones(1,nparas).*np;
    case nparas
        np_list = np;
    otherwise
        error("The dimensions of inputs are not compatible!")
end
%% input vectors
vecs = cell(1,nparas);
switch opts.mode
    case 'uniform'        
        for i = 1:nparas
            vecs{i} = linspace(AtoB(i,1), AtoB(i,2), np_list(i));
        end
    case 'random'
        for i = 1:nparas
            vecsi = AtoB(i,1) + (AtoB(i,2)-AtoB(i,1)).*rand(1, np_list(i));
            vecs{i} = sort(vecsi);
        end
end
%% EVAL can handle arbitrary inputs and outputs, but it is not recommended by MATLAB
%% So we use SWITCH with fixed inputs and outputs instead.
switch nparas
    case 2
        [p1, p2] = ndgrid(vecs{1}, vecs{2});
        mesh_list = [p1(:), p2(:)];
    case 3
        [p1, p2, p3] = ndgrid(vecs{1}, vecs{2}, vecs{3});
        mesh_list = [p1(:), p2(:), p3(:)];
    case 4
        [p1, p2, p3, p4] = ndgrid(vecs{1}, vecs{2}, vecs{3}, vecs{4});
        mesh_list = [p1(:), p2(:), p3(:), p4(:)];
    otherwise
        error("dim="+nparas+" not support, sorry")
end
end