function [Chi_abc_fermi, G_ac_reshape, G_bc_reshape, Chi_abc_reshape]= SOAHC_int(Ham_obj,klist,sub_index,opts)
% intrinsic 2nd order Anomalous Hall effect
% ref: 10.1103/PhysRevLett.127.277202
% in SI unit, Ampere*Volt^-2*meter for 2D case
% klist == klist_s for HR, klist_r for HK and Htrig
arguments
    Ham_obj {mustBeA(Ham_obj,{'HR','HK','Htrig'})};
    klist double; % klist_s for HR, klist_r for HK and Htrig
    sub_index string = "xyy";
    opts.T double = 20; % temperature (Kelvin)
    opts.Ef double = 0;
    opts.Ef_num double = 100;
    opts.Ef_range double = [-1,1];
    opts.plot logical = true;
    opts.plotBZ logical = false;
    opts.plotBZ_bands double = 0;
    opts.nk double = [];
end
%% kmesh info
nkpts = size(klist,1);
%% index convert
sub_index = char(sub_index);
a_index = string(sub_index(1));
b_index = string(sub_index(2));
c_index = string(sub_index(3));
%% constant
charge = 1.6e-19; % C
hbar = 6.5821e-16; % eV.s
%% Fermi info
if length(opts.Ef_range) == 1
    Ef_list = opts.Ef_range;
    opts.plot = false;
else
    Ef_list = linspace(opts.Ef_range(1),opts.Ef_range(2),opts.Ef_num);
end
Ef_list_real = Ef_list + opts.Ef;
nef = length(Ef_list_real);
E_min = min(Ef_list_real) - 0.1;
E_max = max(Ef_list_real) + 0.1;

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
        [EIGENCAR, WAVECAR] = Ham_obj.EIGENCAR_gen('klist',klist,'printmode',false);
    case {"Htrig","HK"}
        [dH_dkx_fun,dH_dky_fun,~] = Ham_diff(Ham_obj);
        [EIGENCAR, WAVECAR] = Ham_obj.EIGENCAR_gen('klist',klist,'printmode',false);
end
%% kubo loop
G_ac = zeros(nkpts,nbands);
G_bc = zeros(nkpts,nbands);
Chi_abc = zeros(nkpts,nbands);
Chi_abc_fermi = zeros(nef,1);
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
            FactorListki = exp(1i*2*pi*vectorList*klist(kn,:).');
%             % HRmat
%             HRmat = sum(pagemtimes(HnumList,reshape(FactorListki,[1 1 NRPTS_])),3);
%             HRmat = (HRmat+HRmat')/2;
%             [WAV_kn,EIG_kn] = eig(HRmat);
%             [WAV_kn,EIG_kn] = vasplib.sorteig(EIG_kn,WAV_kn);
%             EIG_kn = diag(EIG_kn);

            dH_dkx = sum(pagemtimes(HnumLpx,reshape(FactorListki,[1 1 NRPTS_])),3);
            dH_dky = sum(pagemtimes(HnumLpy,reshape(FactorListki,[1 1 NRPTS_])),3);   
        case {"HK","Htrig"}
            kx = klist(kn,1); ky = klist(kn,2); kz = klist(kn,3);   
            dH_dkx = dH_dkx_fun(kx,ky,kz);
            dH_dky = dH_dky_fun(kx,ky,kz);
    end
    
    vx = WAV_kn' * dH_dkx * WAV_kn;
    vy = WAV_kn' * dH_dky * WAV_kn;
    switch a_index
        case "x"
            va = vx;
        case "y"
            va = vy;
    end
    switch b_index
        case "x"
            vb = vx;
        case "y"
            vb = vy;
    end
    switch c_index
        case "x"
            vc = vx;
        case "y"
            vc = vy;
    end

    dEda = real(diag(va));
    dEdb = real(diag(vb));

    G_ac_tmp = zeros(nbands,1);
    G_bc_tmp = zeros(nbands,1);

    for n = nbands_in
        for p = 1:nbands
            dEnp = EIG_kn(n) - EIG_kn(p);
            if abs(dEnp) < 1e-10
                continue
            end
            G_ac_tmp(n) = G_ac_tmp(n) + real(va(n,p)*vc(p,n))/dEnp^3;
            G_bc_tmp(n) = G_bc_tmp(n) + real(vb(n,p)*vc(p,n))/dEnp^3;
        end
    end
    Chi_abc_tmp = G_bc_tmp.*dEda-G_ac_tmp.*dEdb;
    
    dfdE_tmp = transport.Fermi_1(EIG_kn,opts.Ef,opts.T);
    
    G_ac(kn,1:nbands) = G_ac_tmp .*dfdE_tmp;
    G_bc(kn,1:nbands) = G_bc_tmp .*dfdE_tmp;
    Chi_abc(kn,1:nbands) = Chi_abc_tmp .*dfdE_tmp;
    
    for fi = 1:nef
        dfdE   = transport.Fermi_1(EIG_kn,Ef_list_real(fi),opts.T);
        Chi_abc_fermi(fi) = Chi_abc_fermi(fi) + sum(Chi_abc_tmp.*dfdE);            
    end
end
toc
%%
Area = cross(Ham_obj.Rm(1,:),Ham_obj.Rm(2,:));
Area = Area(3);

Chi_abc_fermi = -2 * charge/hbar * (1e-10 /Area /nkpts) * Chi_abc_fermi;
%%
% disp("integration domain: "+count*dk_s*100+" % of 1st BZ")
% disp("please check the convergence of integration domain yourself")
%% plot
if opts.plot
    f11 = Figs(1,1);
    hold on
    plot(f11.axes(1,1),Chi_abc_fermi,Ef_list,'linewidth',3);
    plot(f11.axes(1,1),zeros(nef,1),Ef_list,'--','linewidth',1)
    hold off
    xlabel("\sigma_{"+a_index+b_index+c_index+"}")
    ylabel("E(eV)")
    title("2nd AHC")
end


if opts.plotBZ
    if isempty(opts.nk)
        nk_guess = sqrt(nkpts);
        if nk_guess == round(nk_guess)
            nk = [nk_guess nk_guess 1];
        else
            error('opts.nk is needed to plot distributions in BZ');
        end
    else
        nk = opts.nk;
    end
    
    if opts.plotBZ_bands == 0
        ib = 1:nbands/2;
    else
        ib = opts.plotBZ_bands;
    end
    
    G_ac_reshape = reshape(sum(G_ac(:,ib),2),nk);
    G_bc_reshape = reshape(sum(G_bc(:,ib),2),nk);    
    Chi_abc_reshape = reshape(sum(Chi_abc(:,ib),2),nk);
    
    G_ac_reshape_2 = [G_ac_reshape, -G_ac_reshape];
    G_bc_reshape_2 = [G_bc_reshape, -G_bc_reshape];
    Chi_abc_reshape_2 = [Chi_abc_reshape, -Chi_abc_reshape];
    
    f13 = Figs(1,3);
    colormap(ColorMap.redblue(101))
    
    surf(f13.axes(1,1),G_ac_reshape_2);
    surf(f13.axes(1,2),G_bc_reshape_2);
    surf(f13.axes(1,3),Chi_abc_reshape_2);
    title(f13.axes(1,1),"G_{"+a_index+c_index+"}")
    title(f13.axes(1,2),"G_{"+b_index+c_index+"}")
    title(f13.axes(1,3),"\Lambda_{"+a_index+b_index+c_index+"}")
    
    xlabel(f13.axes(1,1),"k_y")
    xlabel(f13.axes(1,2),"k_y")
    xlabel(f13.axes(1,3),"k_y")
    
    ylabel(f13.axes(1,1),"k_x")
    ylabel(f13.axes(1,2),"k_x")
    ylabel(f13.axes(1,3),"k_x")
    
    for i = 1:3
        xlim(f13.axes(1,i),[1 nk(2)])
        ylim(f13.axes(1,i),[1 nk(1)])
        shading(f13.axes(1,i),'interp');
        xticks(f13.axes(1,i),[]);
        yticks(f13.axes(1,i),[]);
        % xlabel(f13.axes(1,i),"k_y");
        % ylabel(f13.axes(1,i),"k_x");
        % xticks(f13.axes(1,i),[0 nk(1)/2 nk(1)])
        % yticks(f13.axes(1,i),[0 nk(2)/2 nk(2)])
        % xticklabels(f13.axes(1,i),["" "k_y" ""])
        % yticklabels(f13.axes(1,i),["" "k_x" ""])
    end
end
end