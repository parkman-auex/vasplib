function logi = isTopoEquv_M(H1, H2, ori_dim, noccu)
% short name for 'is Topological Equivalent Metal'
% compare two topological GAPLESS systems by
% adding an extra dimension to a N-dimentional gapless model,
% to construct a (N+1)-dimensional gapped model
arguments
    H1
    H2
    ori_dim int16 % original dimension
    noccu int16    
end
%% precision control
metal_threshold = 1e-4; % after flattening, so the energy level of H1 does not matter
nk = 5; % nk<=4 may miss the gapless points at rare cases
options = optimset('TolFun',metal_threshold/10,'TolX',metal_threshold/10,'Display','off');
nbands = H1.Nbands;
%% generate initial k-mesh
switch ori_dim + 1
    case 1
        [klist_s,klist_r] = kmesh_gen(H1,[],'nk',[nk, 1, 1]);        
    case 2
        [klist_s,klist_r] = kmesh_gen(H1,[],'nk',[nk, nk, 1]);
    case 3
        [klist_s,klist_r] = kmesh_gen(H1,[],'nk',[nk, nk, nk]);
end
nkpts = size(klist_s,1);
%% handle of band gap search
get_gap_kpoint = @(kpoint) get_gap_flatten_M(H1,H2, kpoint, nbands, noccu);
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

function gap = get_gap_flatten_M(H1,H2, kp, nbands,noccu)
[~, WAV1] = H1.EIGENCAR_gen('klist',kp,'printmode',false);
[~, WAV2] = H2.EIGENCAR_gen('klist',kp,'printmode',false);
psi1 = WAV1(:,1:noccu);
Q1 = eye(nbands) - 2*(psi1*psi1');
psi2 = WAV2(:,1:noccu);
Q2 = eye(nbands) - 2*(psi2*psi2');
gap = min(abs(eig(Q1+Q2)));
end