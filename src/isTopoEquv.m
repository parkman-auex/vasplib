function logi = isTopoEquv(H1, H2, klist, noccu)
% compare two topological GAPPED systems
arguments
    H1
    H2
    klist (:,3) double
    noccu int16    
end
%% precision control
metal_threshold = 1e-6; % after flattening, so the energy level of H1 does not matter
options = optimset('TolFun',metal_threshold/10,'TolX',metal_threshold/10,'Display','off');

nbands = H1.Nbands;
nkpts = size(klist,1);
%% handle of band gap search
get_gap_kpoint = @(kpoint) get_gap_flatten(H1,H2, kpoint, nbands, noccu);
%%
logi = true;
for i = 1:nkpts
    [~,gapout]= fminsearch(get_gap_kpoint,klist(i,:),options);
    if gapout < metal_threshold
        logi = false;
        break
    end            
end
end

function gap = get_gap_flatten(H1, H2, kp, nbands, noccu)
[~, WAV1] = H1.EIGENCAR_gen('klist',kp,'printmode',false);
[~, WAV2] = H2.EIGENCAR_gen('klist',kp,'printmode',false);
psi1 = WAV1(:,1:noccu);
Q1 = eye(nbands) - 2*(psi1*psi1');
psi2 = WAV2(:,1:noccu);
Q2 = eye(nbands) - 2*(psi2*psi2');

Q_alpha = Q1+Q2;
gap = min(abs(eig(Q_alpha)));
end