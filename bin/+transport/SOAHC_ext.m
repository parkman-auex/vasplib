function D_bd_fermi = SOAHC_ext(Ham_obj,kstruct,sub_index,opts)
% extrinsic 2nd order Anomalous Hall effect, or the so-called nonliner Hall effect(NHE).
% ref: 10.1103/PhysRevLett.115.216806 
% this version only calculates the BCD !!!
arguments
    Ham_obj;
    kstruct struct;
    sub_index string = "xz";
    opts.T double = 20; % temperature (Kelvin)
    opts.Ef double = 0;
    opts.Ef_num double = 100;
    opts.Ef_range double = [-1,1];
    opts.plot logical = true;
    opts.plotBZ logical = true;
    opts.plotBZ_bands double = 0; 
end
%% kmesh info
klist_s = kstruct.klist_s;
klist_r = kstruct.klist_r;
dk_s = kstruct.dk_s;
dk_r = kstruct.dk_r;
nkpts = size(klist_s,1);
%% index convert
sub_index = char(sub_index);
b_index = string(sub_index(1));
d_index = string(sub_index(2));
%% Fermi info
if length(opts.Ef_range) == 1
    Ef_list = opts.Ef_range;
    opts.plot = false;
else
    Ef_list = linspace(opts.Ef_range(1),opts.Ef_range(2),opts.Ef_num);
end
Ef_list_real = Ef_list + opts.Ef;
nef = length(Ef_list_real);
E_min = min(Ef_list_real) - 0.05;
E_max = max(Ef_list_real) + 0.05;

nbands = Ham_obj.Nbands;
%% prepare dH_dk
switch class(Ham_obj)
    case "HR"
        HnumList = Ham_obj.HnumL;
        NRPTS_ = Ham_obj.NRPTS;
        vectorList = double(Ham_obj.vectorL);
        vectorList_r = vectorList * Ham_obj.Rm;

        % partial R
        HnumLpx = 1i*pagemtimes(reshape(vectorList_r(:,1),[1 1 NRPTS_]),HnumList);
        HnumLpy = 1i*pagemtimes(reshape(vectorList_r(:,2),[1 1 NRPTS_]),HnumList);
        HnumLpz = 1i*pagemtimes(reshape(vectorList_r(:,3),[1 1 NRPTS_]),HnumList);
        [EIGENCAR, WAVECAR] = Ham_obj.EIGENCAR_gen('klist',klist_s,'printmode',false);
    case {"Htrig","HK"}
        [dH_dkx_fun,dH_dky_fun,dH_dkz_fun] = Ham_diff(Ham_obj);
        [EIGENCAR, WAVECAR] = Ham_obj.EIGENCAR_gen('klist',klist_r,'printmode',false);
end
%% kubo loop
D_bd = zeros(nkpts,nbands);
D_bd_fermi = zeros(nef,1);
hclass = class(Ham_obj);
count = 0;

tic
for kn = 1:nkpts
    WAV_kn = WAVECAR(:,:,kn);
    EIG_kn = EIGENCAR(:,kn);
    
    % band index in Fermi energy neighborhood
    nbands_in = find(EIG_kn<E_max & EIG_kn>E_min).'; 
    if sum(nbands_in) == 0       
        continue
    end
    count = count + 1;
    
    switch hclass        
        case "HR"
            % efactor R
            FactorListki = exp(1i*2*pi*vectorList*klist_s(kn,:).');
            dH_dkx = sum(pagemtimes(HnumLpx,reshape(FactorListki,[1 1 NRPTS_])),3);
            dH_dky = sum(pagemtimes(HnumLpy,reshape(FactorListki,[1 1 NRPTS_])),3);
            dH_dkz = sum(pagemtimes(HnumLpz,reshape(FactorListki,[1 1 NRPTS_])),3);
        case {"HK","Htrig"}
            kx = klist_r(kn,1); ky = klist_r(kn,2); kz = klist_r(kn,3);   
            dH_dkx = dH_dkx_fun(kx,ky,kz);
            dH_dky = dH_dky_fun(kx,ky,kz);
            dH_dkz = dH_dkz_fun(kx,ky,kz);
    end
    
    vx = WAV_kn' * dH_dkx * WAV_kn;
    vy = WAV_kn' * dH_dky * WAV_kn;
    vz = WAV_kn' * dH_dkz * WAV_kn;
    switch b_index
        case "x"
            vb = vx;
        case "y"
            vb = vy;
    end
    switch d_index
        case "x"
            v1 = vy;
            v2 = vz;
        case "y"
            v1 = vz;
            v2 = vx;
        case "z"
            v1 = vx;
            v2 = vy;
    end
    
    dEdb = real(diag(vb));    
    omega_d_tmp = zeros(nbands,1);

    for n = nbands_in
        for p = 1:nbands
            dEnp = EIG_kn(n) - EIG_kn(p);
            if abs(dEnp) < 1e-6
                continue
            end                   
            omega_d_tmp(n) = omega_d_tmp(n) + imag(v1(n,p)*v2(p,n))/dEnp^2;
        end
    end
    D_bd_tmp = omega_d_tmp.*dEdb;
    D_bd(kn,1:nbands) = D_bd_tmp;
    
    for fi = 1:nef
        dfdE   = transport.Fermi_1(EIG_kn,Ef_list_real(fi),opts.T);
        D_bd_fermi(fi) = D_bd_fermi(fi) + sum(D_bd_tmp.*dfdE);
    end
end     
toc
%%
D_bd_fermi = dk_r*(2*pi)^-2 * (-1)*(-2) .* D_bd_fermi;
%%
disp("integration domain: "+count*dk_s*100+" % of 1st BZ")
disp("please check the convergence of integration domain yourself")
%% plot
f12 = Figs(1,2);
if opts.plot  
    hold on
    plot(f12.axes(1,1),Ef_list,D_bd_fermi,'linewidth',3);
    plot(f12.axes(1,1),Ef_list,zeros(nef,1),'--','linewidth',1)
    hold off
    ylabel(f12.axes(1,1),"D_{"+b_index+d_index+"}")
    xlabel(f12.axes(1,1),"Ef(eV)")
    title("2nd AHC")
end
if opts.plotBZ
    if opts.plotBZ_bands == 0
        ib = 1:nbands/2;
    else
        ib = opts.plotBZ_bands;
    end  
    D_bd_reshape = reshape(sum(D_bd(:,ib),2),kstruct.nk);  
    contourf(f12.axes(1,2),D_bd_reshape);
    title(f12.axes(1,2),"D_{"+b_index+d_index+"}")   
end
end