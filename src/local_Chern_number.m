function [C_local, C_orb_spinful, C_orb_spinless] = local_Chern_number(Ham_obj,opts)
% for 2D slab model
arguments
    Ham_obj;
    opts.Num_Occupied = 0;
    opts.Nz double = 1; % the number of layers in z-direction
    opts.nk double = 100;
    opts.plot logical = true
end
%% basic info
[klist_s, klist_r] = kmesh_gen(Ham_obj,"nk",[opts.nk opts.nk 1]);
if sum(abs(Ham_obj.vectorL)) == 0
    klist_s = [0 0 0];
    klist_r = [0 0 0];
    opts.plot = false;
end
nkpts = size(klist_s,1);
Area = cross(Ham_obj.Rm(:,1),Ham_obj.Rm(:,2));
Area = Area(3);

nbands = Ham_obj.Nbands;
if opts.Num_Occupied == 0
    Nc = nbands/2;
    disp('Num_Occupied not give, suppose half occupied')
else
    Nc = opts.Num_Occupied;
end
%% prepare dH_dk
switch class(Ham_obj)
    case "HR"
        HnumList = Ham_obj.HnumL;
        NRPTS_ = Ham_obj.NRPTS;
        vectorList = double(Ham_obj.vectorL);
        Ham_obj = Ham_obj.tjmti_gen();

        tji_mat_cart = Ham_obj.tjmti{1};
        tji_mat_frac = Ham_obj.tjmti{2};

        % vectorList_cart = vectorList*Ham_obj.Rm;
        vectorList_r = vectorList * Ham_obj.Rm;
        % partial R
        HnumLpx = 1i*pagemtimes(reshape(vectorList_r(:,1),[1 1 NRPTS_]),HnumList);
        HnumLpy = 1i*pagemtimes(reshape(vectorList_r(:,2),[1 1 NRPTS_]),HnumList);
        % partial titj
        HnumLpx_tji = 1i*tji_mat_cart(:,:,1);
        HnumLpy_tji = 1i*tji_mat_cart(:,:,2);
    case {"Htrig","HK"}
        [dH_dkx_fun,dH_dky_fun,~] = Ham_diff(Ham_obj);
        [EIGENCAR, WAVECAR] = Ham_obj.EIGENCAR_gen('klist',klist_r);
end
%% real space kubo loop
tic
C_orb = zeros(nbands,1);
for ki = 1:nkpts
    switch class(Ham_obj)
        case "HR"
            FactorListki = exp(1i*2*pi*vectorList*klist_s(ki,:).'); 
            % HRmat here
            HRmat = sum(pagemtimes(HnumList,reshape(FactorListki,[1 1 NRPTS_])),3); % sum
            kjiL_x =  tji_mat_frac(:,:,1).*klist_s(ki,1);
            kjiL_y =  tji_mat_frac(:,:,2).*klist_s(ki,2);           
            Hmat_tji = exp(1i*2*pi*(kjiL_x+kjiL_y));

            HRmat = HRmat .* Hmat_tji;
            HRmat = (HRmat+HRmat')/2;

            HRmatpA = sum(pagemtimes(HnumLpx,reshape(FactorListki,[1 1 NRPTS_])),3);
            HRmatpB = sum(pagemtimes(HnumLpy,reshape(FactorListki,[1 1 NRPTS_])),3);

            Hmat_tijpA = Hmat_tji.* HnumLpx_tji;
            Hmat_tijpB = Hmat_tji.* HnumLpy_tji;
            % (AB)' = A'B + AB'
            dH_dkx = HRmatpA.*Hmat_tji + HRmat.*Hmat_tijpA;% vx partial_A_tmp
            dH_dky = HRmatpB.*Hmat_tji + HRmat.*Hmat_tijpB;% vy partial_B_tmp

            % vx_oper = HRmatpA;
            % vy_oper = HRmatpB;
            
            [WAV,EIG] = eig(HRmat);
            [WAV,EIG] = park.sorteig(EIG,WAV);
            EIG = diag(EIG);
            WAV_v = WAV(:, 1:Nc);
            WAV_c = WAV(:, Nc+1:nbands);
        case {"HK","Htrig"}
            kx = klist_r(ki,1); ky = klist_r(ki,2); kz = klist_r(ki,3);   
            dH_dkx = dH_dkx_fun(kx,ky,kz);
            dH_dky = dH_dky_fun(kx,ky,kz);
            
            EIG = EIGENCAR(:,ki);
            WAV_v = WAVECAR(:, 1:Nc,       ki);
            WAV_c = WAVECAR(:, Nc+1:nbands,ki);
    end
     
    Xvc = (WAV_v' * 1i*dH_dkx * WAV_c);
    Yvc = (WAV_v' * 1i*dH_dky * WAV_c);
    
    for vi = 1:Nc
        for ci = 1:(nbands-Nc)
            dEcv = EIG(ci+Nc) - EIG(vi);
            if abs(dEcv) < 1e-6
                Xvc(vi,ci) = 0;
                Yvc(vi,ci) = 0;
            else
                Xvc(vi,ci) = Xvc(vi,ci) / dEcv;
                Yvc(vi,ci) = Yvc(vi,ci) / dEcv;
            end
        end
    end      
    
    % Xvc1 = (WAV_v' * vx1 * WAV_c) ./ dEcv;
    % Yvc1 = (WAV_v' * vy1 * WAV_c) ./ dEcv;
        
    C_orb = C_orb + imag(diag(WAV_v * Xvc * Yvc' * WAV_v'));  
end
toc
C_orb = -4*pi/nkpts * C_orb;
%% reduce spin freedom
C_orb_spinful = C_orb;
C_orb_spinless = zeros(nbands/2,1);
% specfic for four-band model in 1st BZ
for i = 1:nbands/4
    C_orb_spinless(2*i-1) = C_orb(4*i-3) + C_orb(4*i-2);
    C_orb_spinless(2*i  ) = C_orb(4*i-1) + C_orb(4*i  );
end
%% orb degeneracy
% xyz = Ham_obj.orbL;
% for bi = 1:nbands
%     orb_Nd = sum(ismember(xyz,xyz(bi,1:3),'rows'));
%     C_orb(bi) = C_orb(bi)*orb_Nd;
% end
%%
% nsites = nbands/opts.Nz;
C_local = zeros(opts.Nz, 1);
for iz = 1:opts.Nz
    C_local_index = Ham_obj.orbL(:,3) >= (iz - 1)/opts.Nz &...
        Ham_obj.orbL(:,3) < iz/opts.Nz;
    C_local(iz) = sum(C_orb(C_local_index))/Area;
end
C_local = round(C_local,5);
%% plot C_local
if opts.plot == true
    fig = Figs(1,1);
    C_sum = zeros(opts.Nz,1);
    for i = 1:opts.Nz
        C_sum(i) = sum(C_local(1:i));
    end

    C_sum = [0;C_sum];

    hold on
    area(fig.axes(1,1),0:opts.Nz, C_sum,...
        'LineWidth',2,...
        'EdgeColor','#E6AE41',...
        'FaceColor','#FAEED1');
    plot(fig.axes(1,1),1:opts.Nz, C_local,...
        '-o','LineWidth',3,'MarkerSize',10,'Color','#0072BD');
    hold off

    yticks([-1 -0.5 0 0.5 1]);
    xlabel('Layer Index');
    ylabel('Layer Chern Marker');
end
end