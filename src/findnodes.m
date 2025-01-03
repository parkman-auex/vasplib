function [nodes_s, nodes_r] = findnodes(Ham_obj,opts,kopts)
arguments
    Ham_obj vasplib
    opts.Num_Occupied int8 = 0
    opts.Gap_Threshold double = 0.0001
    opts.debug_mode logical = false
    
    kopts.nk int8 = [10 10 10]; % [nk1 nk2 nk3]
    kopts.vk = [1 0 0; 0 1 0; 0 0 1]; % [vk1; vk2; vk3]
    kopts.original_point = [-0.5 -0.5 -0.5];
    kopts.mode {mustBeMember(kopts.mode,{'corner','center'})} = 'corner'
    kopts.edge {mustBeMember(kopts.edge,{'half','full'})} = "half";
end
koptscell = namedargs2cell(kopts);
if opts.Num_Occupied == 0
    nbands = Ham_obj.Nbands;
    occu = nbands/2;
else
    occu = opts.Num_Occupied;
end
%%
[klist_s, klist_r] = kmesh_gen(Ham_obj,[],koptscell{:});
switch class(Ham_obj)
    case 'HR'
        get_gap_kpoint = @(kpoint) get_gap(Ham_obj,kpoint,occu);
        klist = klist_s;
    case {'Htrig','HK'}
        Hfun = Ham_obj.Hfun;
        get_gap_kpoint = @(kpoint) get_gap_Hfun(Hfun,kpoint,occu);
        klist = klist_r;
end
nkpts = size(klist_s,1);
nodes_tmp = [];

if ~opts.debug_mode
    pb = vasplib_tool_outer.CmdLineProgressBar('FindNode: ');
    count = 0;
    msgtail = [' hit:',num2str(count),'/',num2str(nkpts)];
    for i = 1:nkpts
        [kout,gapout]= fminsearch(get_gap_kpoint,klist(i,:));
        if gapout < opts.Gap_Threshold
            nodes_tmp = [nodes_tmp;kout];
            count = count +1;
            msgtail = [' hit:',num2str(count),'/',num2str(nkpts)];
        end
        pb.print(i,nkpts,msgtail);
    end
    pb.delete();
else
    optm_opts = optimset('Display','off');
    for i = 1:nkpts
        [kout,gapout]= fminsearch(get_gap_kpoint,klist(i,:),optm_opts);
        if gapout < opts.Gap_Threshold
            nodes_tmp = kout;
            break
        end
    end
end

for i = 1:3
    if kopts.nk(i) == 1
        nodes_tmp(:,i) = 0;
    end
end
%% to fractional nodes
switch class(Ham_obj)
    case "HR"
        if isempty(nodes_tmp)
            nodes_s = [];
            nodes_r = [];
            return;
        else
            nodes_s = nodes_tmp;
        end
    case {"Htrig","HK"}
        if isempty(nodes_tmp)
            nodes_s = [];
            nodes_r = [];
            return;
        else
            nodes_s = nodes_tmp / Ham_obj.Gk;
        end
end
%% remove to a given user defined block
nodes_s = kshift(nodes_s,[kopts.original_point;kopts.vk]);
nodes_s = uniquetol(nodes_s,opts.Gap_Threshold*10,'ByRows',true);
nodes_r = nodes_s * Ham_obj.Gk;
end
function gap = get_gap_Hfun(Hfun,kpoint,occu)
H = Hfun(kpoint(1),kpoint(2),kpoint(3));
EIGENCAR = (eig((H + H')/2));
gap = EIGENCAR(occu+1,:) - EIGENCAR(occu,:);
end

function gap = get_gap(Ham_obj,kpoint,occu)
EIGENCAR = Ham_obj.EIGENCAR_gen('klist',kpoint,'printmode',false);
gap = EIGENCAR(occu+1,:) - EIGENCAR(occu,:);
end