function [Rm,sites]=supercell(Ns,Rm,sites,Atom_name,Atom_num,findir,filename)
%% usage: [sites]=supercell(Ns)
% Ns =
%
% [ N_11, N_12, N_13]
% [ N_21, N_22, N_23]
% [ N_31, N_32, N_33]
% sites
% note : site=struct('seq',[],'inseq',[],'rc1',[],'rc2',[],'rc3',[],'name',[],'nameseq',[]);
%   Rm
%  a1 a1x a1y a1z
%  a2 a2x a2y a2z
%  a3 a3x a3y a3z
%   unit A?e-6
  Accuracy=8;
%% nargin
if nargin < 2
    [Rm,sites,Atom_name,Atom_num,elements,a_crystal_constance]=POSCAR_read;
    findir = [0,0,0];
    filename  = 'POSCAR_super';
end
if nargin < 6
    findir = [0,0,0];
    filename  = 'POSCAR_super';
end

if nargin < 7
    filename  = 'POSCAR_super';
end

    format long;
    %Ns_=inv(Ns);
    V_Ns=dot(Ns(1,:),cross(Ns(2,:),Ns(3,:)));
    if V_Ns < 0
        Ns(3,:) = -Ns(3,:);
    end
    V_Ns = abs(V_Ns);
    if V_Ns < 1
        error('check Ns');
    end

    % refresh data first
    [~,nsites] = size(sites);
    for i =1:nsites
        sites(i).rc1 = plusrc(sites(i).rc1);
        if sites(i).rc1  < 0 
            disp('bugtest')
        end
        sites(i).rc2 = plusrc(sites(i).rc2);
        if sites(i).rc2  < 0 
            disp('bugtest')
        end
        sites(i).rc3 = plusrc(sites(i).rc3);
        if sites(i).rc3  < 0 
            disp('bugtest')
        end
    end


    % -------------------------------------------------------------------------
    % old algorihm is not so good, we just copy from pythonTB
    % if V_Ns > 1
    %     Range1=[floor(min(Ns(:,1)))-1,ceil(max(Ns(:,1)))+1];
    %     Range2=[floor(min(Ns(:,2)))-1,ceil(max(Ns(:,2)))+1];
    %     Range3=[floor(min(Ns(:,3)))-1,ceil(max(Ns(:,3)))+1];
    %
    %     countnum=0;
    %     countnum2=0;
    %     for elementseq=1:length(Atom_num)
    %         % for each Atom pattern
    %             for n=1:Atom_num(elementseq)
    %              % for each atom in Atom pattern
    %                 countnum2=countnum2+1;
    %                 for i=Range1(1):Range1(2)
    %                     for j=Range2(1):Range2(2)
    %                         for k=Range3(1):Range3(2)
    %                             %countnum
    %
    %     %                         Tranlation_vector=[i j k];
    %     %                         Tx=Tranlation_vector(1);
    %     %                         Ty=Tranlation_vector(2);
    %     %                         Tz=Tranlation_vector(3);
    %                             Rpc1=sites(countnum2).rc1+i;
    %                             Rpc2=sites(countnum2).rc2+j;
    %                             Rpc3=sites(countnum2).rc3+k;
    %                             Rsc=[Rpc1 Rpc2 Rpc3]*inv(Ns);% use A\b for inv(A)*b or ???
    %                             rsc1=Rsc(1);
    %                             rsc2=Rsc(2);
    %                             rsc3=Rsc(3);
    %                             if rangein(Rsc)
    %                                 countnum=countnum+1;
    %                                 sites_s(countnum).rc1=rsc1;%along a1
    %                                 sites_s(countnum).rc2=rsc2;%along a2
    %                                 sites_s(countnum).rc3=rsc3;%along a3
    %                                 sites_s(countnum).name=Atom_name(elementseq);
    %                             end
    %                         end
    %                     end
    %                 end
    %
    %
    %             end
    %     end
    %
    % sites=sites_s;
    % end
    % -------------------------------------------------------------------------
    if V_Ns >= 1
        % conservative estimate on range of search for super-cell vectors
        max_R=max(abs(Ns))*3;
        % candidates for super-cell vectors
        % this is hard-coded and can be improved!
        sc_cands = [];
        for i = -max_R(1):max_R(1)
            for j = -max_R(2):max_R(2)
                for k = -max_R(3):max_R(3)
                    sc_cands=[sc_cands;[i,j,k]];
                end
            end
        end

        % find all vectors inside super-cell
        % store them here
        sc_vec=[];
        eps_shift=sqrt(2.0)*1.0E-8;% shift of the grid, so to avoid double counting

        %
        for ivec = 1:length(sc_cands)
            % compute reduced coordinates of this candidate vector in the super-cell frame
            tmp_red = park.to_red_sc(sc_cands(ivec,:),Ns);
            % check if in the interior
            inside = 1;
            for it = 1:length(tmp_red)
                t = tmp_red(it);
                if t <= -1.0*eps_shift | t>1.0-eps_shift
                    inside=0;
                end
            end
            if inside == 1
                %disp(tmp_red);
                %disp(sc_cands(ivec,:));
                sc_vec=[sc_vec;sc_cands(ivec,:)];
                % number of times unit cell is repeated in the super-cell
            end
        end

        [num_sc,~] = size(sc_vec);
        %disp(num_sc);
        % check that found enough super-cell vectors
        if round(round(abs(det(Ns)))) ~= num_sc
            error("\n\nSuper-cell generation failed! Wrong number of super-cell vectors found.");
        end

        % cartesian vectors of the super lattice
        % sc_cart_lat=Ns*Rm; %must Ns Rm
        countnum=0;
        countnum2=0;
        for elementseq=1:length(Atom_num)
            % for each Atom pattern
            for n=1:Atom_num(elementseq)
                % for each atom in Atom pattern
                countnum2=countnum2+1;
                Rpc1=sites(countnum2).rc1;
                Rpc2=sites(countnum2).rc2;
                Rpc3=sites(countnum2).rc3;
                Rc=[Rpc1 Rpc2 Rpc3];

                for icur_sc_vec = 1:num_sc % go over all super-cell vectors
                    cur_sc_vec = sc_vec(icur_sc_vec,:);
                    Rsc = park.to_red_sc(Rc+cur_sc_vec,Ns);
                    Rsc = round(Rsc.*10^Accuracy)./10^Accuracy;
                    rsc1=Rsc(1);
                    rsc2=Rsc(2);
                    rsc3=Rsc(3);
                    if rsc1 ==1
                        rsc1 = 0; 
                    end
                    if rsc2 ==1
                        rsc2 = 0;
                    end
                    if rsc3 ==1
                        rsc3 = 0;
                    end
                    countnum=countnum+1;
                    sites_s(countnum).rc1=rsc1;%along a1
                    sites_s(countnum).rc2=rsc2;%along a2
                    sites_s(countnum).rc3=rsc3;%along a3
                    sites_s(countnum).name=Atom_name(elementseq);
                end

            end
        end
        sites=sites_s;
    end
    % Atom_name Keep
    % Atom_number [Rm,sites]=supercell(Ns,Rm,sites,Atom_name,Atom_num,fin_dir,filename);
    Atom_num_s=round(Atom_num*V_Ns);
    Rm_s=Ns*Rm;
    Rm=Rm_s;
    [~,nsites] = size(sites);

    if findir(1) ~=0 | findir(2) ~=0 | findir(3) ~=0
        disp('findir_mode');
        %init

        Rmlength1 = abs(norm (Rm(1,:)));
        Rmlength2 = abs(norm (Rm(2,:)));
        Rmlength3 = abs(norm (Rm(3,:)));
        Rm_s_fin_add = [10*Rm(1,:)*findir(1)/Rmlength1;...
            10*Rm(2,:)*findir(2)/Rmlength2;...
            10*Rm(3,:)*findir(3)/Rmlength3];
        Rm_s_fin = Rm + Rm_s_fin_add ;
        Rc_s_fin_add = [1/2, 1/2 ,1/2] ;
        Rr_s_fin_add = Rc_s_fin_add * Rm_s_fin_add;
        for  i = 1:nsites
            Rr = [sites(i).rc1,sites(i).rc2,sites(i).rc3]*Rm;
            Rr_s_fin = Rr + Rr_s_fin_add;
            Rc_s_fin = Rr_s_fin /Rm_s_fin;
            sites_s(i).rc1 = Rc_s_fin(1);
            sites_s(i).rc2 = Rc_s_fin(2);
            sites_s(i).rc3 = Rc_s_fin(3);
        end
        Rm_s=Rm_s_fin;
        sites=sites_s;
        Rm=Rm_s;
    end

    % ------------------------------------------------------------------------------
    POSCAR_gen(Rm_s,sites_s,Atom_name,Atom_num_s,filename);


end



function rc_plus = plusrc(rc)
if rc<0
    %disp('bugtest');
    rc_plus = rc +1;

elseif rc>=1
     rc_plus =mod(rc,1);
else
    rc_plus =rc;
end


end

