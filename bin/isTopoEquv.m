function logi = isTopoEquv(H1, H2, dim, noccu, metal_threshold)
% compare two topological GAPPED systems
arguments
    H1
    H2
    dim int16
    noccu int16
    metal_threshold double = 1e-4
end
%% precision control
%% very important !!! nk<=4 may miss the gapless points at rare cases
nbands = H1.Nbands;
nk = 5;
options = optimset('TolFun',metal_threshold/10,'TolX',metal_threshold/10,'Display','off');
%% generate initial k-mesh
switch dim
    case 1
        [klist_s,klist_r] = kmesh_gen(H1,[],'nk',[nk, 1, 1]);        
    case 2
        [klist_s,klist_r] = kmesh_gen(H1,[],'nk',[nk, nk, 1]);
    case 3
        [klist_s,klist_r] = kmesh_gen(H1,[],'nk',[nk, nk, nk]);
end
nkpts = size(klist_s,1);
%% handle of band gap search
get_gap_kpoint = @(kpoint) get_gap_flatten(H1,H2, kpoint, nbands, noccu);
%%
logi = true;
switch class(H1)
    case "HR"
        for i = 1:nkpts
            [~,gapout]= fminsearch(get_gap_kpoint,klist_s(i,:),options);
            if gapout < metal_threshold
                logi = false;
                break
            end
        end
    case {"Htrig","HK"}
        for i = 1:nkpts
            [~,gapout]= fminsearch(get_gap_kpoint,klist_r(i,:),options);
            if gapout < metal_threshold
                logi = false;
                break
            end
        end
end
end

function gap = get_gap_flatten(H1,H2, kp, nbands,noccu)
[~, WAV1] = H1.EIGENCAR_gen('klist',kp,'printmode',false);
[~, WAV2] = H2.EIGENCAR_gen('klist',kp,'printmode',false);
psi1 = WAV1(:,1:noccu);
Q1 = eye(nbands) - 2*(psi1*psi1');
psi2 = WAV2(:,1:noccu);
Q2 = eye(nbands) - 2*(psi2*psi2');
gap = min(abs(eig(Q1+Q2)));
end