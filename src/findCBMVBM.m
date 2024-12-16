function [vbm_E, cbm_E, vbm_position, cbm_position] = findCBMVBM(Ham_obj,opts,kopts)
arguments
    Ham_obj vasplib
    opts.Num_Occupied int8 = 0
    
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
        klist = klist_s;
    case {'Htrig','HK'}
        klist = klist_r;
end

get_cbm_kpoint = @(kpoint) get_cbm(Ham_obj,kpoint,occu);
get_vbm_kpoint = @(kpoint) get_vbm(Ham_obj,kpoint,occu);
cbm_E2 = 1000; % very large value for init
vbm_E2 = 1000;

nkpts = size(klist_s,1);

pb = vasplib_tool_outer.CmdLineProgressBar('Find_CBM: '); 
for i = 1:nkpts  
    [cbm_k, cbm_E1]= fminsearch(get_cbm_kpoint,klist(i,:));
    if cbm_E1 < cbm_E2
        cbm_position = cbm_k;
        cbm_E2 = cbm_E1;
    end
    msgtail = [];%[num2str(i),'/',num2str(nkpts)];
    pb.print(i, nkpts, msgtail);
end
pb.delete();

pb = vasplib_tool_outer.CmdLineProgressBar('Find_VBM: '); 
for i = 1:nkpts
    [vbm_kout, vbm_E1]= fminsearch(get_vbm_kpoint,klist(i,:));
    if vbm_E1 < vbm_E2
        vbm_position = vbm_kout;
        vbm_E2 = vbm_E1;
    end
    msgtail = [];%[num2str(i),'/',num2str(nkpts)];
    pb.print(i, nkpts, msgtail); 
end
pb.delete();

%%
vbm_E = - vbm_E2;
cbm_E = cbm_E2;
end


function cbm = get_cbm(Ham_obj, kpoint, occu)
EIGENCAR = Ham_obj.EIGENCAR_gen('klist',kpoint,'printmode',false);
cbm = EIGENCAR(occu+1);
end

function vbm = get_vbm(Ham_obj, kpoint, occu) % 加一个负号转化为最小值问题
EIGENCAR = Ham_obj.EIGENCAR_gen('klist',kpoint,'printmode',false);
vbm = -EIGENCAR(occu);
end