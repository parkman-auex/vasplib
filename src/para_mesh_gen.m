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
% opts.shift - 'none','uni-rand','tot-rand'
%   'none' - no shift
%   'uni-rand' - uniform" random shift
%   'tot-rand' - total random shift

function mesh_list = para_mesh_gen(AtoB, np_list, opts)
arguments
    AtoB (:,2) double = [0, 1; 0, 1; 0, 1]
    np_list double = [10, 10, 10]
    opts.shift {mustBeMember(opts.shift,{'none','uni-rand','tot-rand'})} = 'none';
end
%% input vectors
nparas = length(np_list);
vecs = cell(1,nparas);
for i = 1:nparas
    vecs{i} = linspace(AtoB(i,1), AtoB(i,2), np_list(i));
end
%% shift
step = (AtoB(:,2) - AtoB(:,1))./np_list';
switch opts.shift
    case 'none'        

    case 'uni-rand'
        for i = 1:nparas
            vecs{i} = vecs{i} + 1e-2*step(i)*(1+rand());
        end
    case 'tot-rand'
        for i = 1:nparas
            vecs{i} = vecs{i} + 1e-2*step(i)*(1+rand(1,np_list(i)));
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
        error("para_dim="+nparas+" not support, sorry")
end
end