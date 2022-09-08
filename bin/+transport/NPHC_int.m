function Chi_ijlq = NPHC_int(Ham_obj,kstruct,sub_index,options)
arguments
    Ham_obj;
    kstruct struct;
    sub_index string;
    options.T double = 20; % temperature Kelvin
    options.Ef double = 0;
    options.Ef_range double = [-1,1];
    options.plot logical = true;
    options.spin_mode = "new";
end
%% kmesh info
klist_s = kstruct.klist_s;
klist_r = kstruct.klist_r;
dk_s = kstruct.dk_s;
dk_r = kstruct.dk_r;
nkpts = size(klist_s,1);
%% index convert
sub_index = char(sub_index);
i_index = string(sub_index(1)); % three electric field indexes
j_index = string(sub_index(2));
l_index = string(sub_index(3));
q_index = string(sub_index(4)); % magnetic field index
%% constant
hbar = 6.5821e-16; % eV.s
mu_B = 5.7884e-5; % eV/Tesla
Lande_g = 2;
%% basic info
if length(options.Ef_range) == 1
    Ef_list = options.Ef_range;
    options.plot = false;
else
    Ef_list = linspace(options.Ef_range(1),options.Ef_range(2),100);
end
Ef_list_real = Ef_list + options.Ef;
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
    case {"Htrig","HK"}
        [dH_dkx_fun,dH_dky_fun,~] = Ham_diff(Ham_obj);
        [EIGENCAR, WAVECAR] = Ham_obj.EIGENCAR_gen('klist',klist_r,'printmode',false);
end
%% pauli matrix
switch q_index
    case "x"
        B = [1 0 0];
    case "y"
        B = [0 1 0];
end
s_dot_B = add_zeeman_field(Ham_obj,B,'spin_mode',options.spin_mode,'natural_unit',true);
%% kubo loop
Chi_ijlq = zeros(nef,1);
count = 0;

tic
for kn = 1:nkpts
    switch class(Ham_obj)
        case "HR"
            % efactor R
            FactorListki = exp(1i*2*pi*vectorList*klist_s(kn,:).');
            % HRmat
            HRmat = sum(pagemtimes(HnumList,reshape(FactorListki,[1 1 NRPTS_])),3);
            HRmat = (HRmat+HRmat')/2;
            [WAV_kn,EIG_kn] = eig(HRmat);
            [WAV_kn,EIG_kn] = vasplib.sorteig(EIG_kn,WAV_kn);
            EIG_kn = diag(EIG_kn);

            % band index in Fermi energy neighborhood
            nbands_in = find(EIG_kn<E_max & EIG_kn>E_min).'; 
            if sum(nbands_in) == 0       
                continue
            end
            count = count + 1;

            HRmatpx = sum(pagemtimes(HnumLpx,reshape(FactorListki,[1 1 NRPTS_])),3);
            HRmatpy = sum(pagemtimes(HnumLpy,reshape(FactorListki,[1 1 NRPTS_])),3);   
            vx = WAV_kn' * HRmatpx * WAV_kn;
            vy = WAV_kn' * HRmatpy * WAV_kn;
        case {"HK","Htrig"}
            WAV_kn = WAVECAR(:,:,kn);
            EIG_kn = EIGENCAR(:,kn);

            % band index in Fermi energy neighborhood
            nbands_in = find(EIG_kn<E_max & EIG_kn>E_min).'; 
            if sum(nbands_in) == 0       
                continue
            end
            count = count + 1;

            kx = klist_r(kn,1); ky = klist_r(kn,2); kz = klist_r(kn,3);   
            vx = WAV_kn' * dH_dkx_fun(kx,ky,kz) * WAV_kn;
            vy = WAV_kn' * dH_dky_fun(kx,ky,kz) * WAV_kn;
    end

    switch i_index
        case "x"
            vi = vx;
        case "y"
            vi = vy;
    end
    switch j_index
        case "x"
            vj = vx;
        case "y"
            vj = vy;
    end
    switch l_index
        case "x"
            vl = vx;
        case "y"
            vl = vy;
    end
    sq = WAV_kn' * s_dot_B * WAV_kn;

    dvi = real(diag(vi));
    dvj = real(diag(vj));
    dsq = real(diag(sq));

    G_il = zeros(nbands,1);
    G_jl = zeros(nbands,1);
    Lambda_ilq = zeros(nbands,1);
    Lambda_jlq = zeros(nbands,1);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              

    for n = nbands_in
        for p = 1:nbands
            dEnp = EIG_kn(n) - EIG_kn(p);
            if abs(dEnp) < 1e-10
                continue
            end

            G_il(n) = G_il(n) + real(vi(n,p)*vl(p,n))/dEnp^3;
            G_jl(n) = G_jl(n) + real(vj(n,p)*vl(p,n))/dEnp^3;
            Lambda_ilq(n) = Lambda_ilq(n) + real(-3*(sq(n,n)-sq(p,p))*vi(n,p)*vl(p,n))/dEnp^4;
            Lambda_jlq(n) = Lambda_jlq(n) + real(-3*(sq(n,n)-sq(p,p))*vj(n,p)*vl(p,n))/dEnp^4;

            for q = 1:nbands
                dEnq = EIG_kn(n) - EIG_kn(q);
                dEpq = EIG_kn(p) - EIG_kn(q);

                if abs(dEnq) > 1e-10
                    Lambda_ilq(n) = Lambda_ilq(n) + ...
                        real(sq(n,q)*(vi(q,p)*vl(p,n)+vi(p,n)*vl(q,p)))/dEnp^3*dEnq;
                    Lambda_jlq(n) = Lambda_jlq(n) + ...
                        real(sq(n,q)*(vj(q,p)*vl(p,n)+vj(p,n)*vl(q,p)))/dEnp^3*dEnq;
                end  


                if abs(dEpq) > 1e-10
                    Lambda_ilq(n) = Lambda_ilq(n) + ...
                        real(sq(p,q)*(vi(q,n)*vl(n,p)+vl(q,n)*vi(n,p)))/dEnp^3*dEpq;
                    Lambda_jlq(n) = Lambda_jlq(n) + ...
                        real(sq(p,q)*(vj(q,n)*vl(n,p)+vl(q,n)*vj(n,p)))/dEnp^3*dEpq;
                end
            end

        end
    end

    for fi = 1:nef
        dfdE_num   = transport.Fermi_1( EIG_kn,Ef_list_real(fi),options.T);
        df2dE2_num = transport.Fermi_2(EIG_kn,Ef_list_real(fi),options.T);
        Chi_ijlq(fi) = Chi_ijlq(fi)...
            + sum((G_jl.*dvi      -G_il.*dvj      ).*df2dE2_num.*dsq)...
            + sum((Lambda_jlq.*dvi-Lambda_ilq.*dvj).*dfdE_num       );
    end
end
toc
%%
factor = Lande_g/2 * mu_B / hbar;
Chi_ijlq = -2* factor * real(Chi_ijlq);
Chi_ijlq = (1.6e-19 * 1e-10 * dk_r * (2*pi)^-2).*Chi_ijlq;
%%
disp("integration domain: "+count*dk_s*100+" % of 1st BZ")
disp("please check the convergence of integration domain yourself")
%% plot
if options.plot
    figure()
    hold on
    plot(Chi_ijlq,Ef_list,'linewidth',3);
    plot(zeros(nef,1),Ef_list,'--','linewidth',1)
    hold off
    xlabel("\chi_{"+i_index+j_index+l_index+q_index+"}")
    ylabel("Ef(eV)")
    title("Nonlinear Planar Hall Conductivity")
end
end