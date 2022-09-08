%% supercell_hr
% Returns tight-binding model representing a super-cell of a current object.
% This function can be used together with *cut_piece* in order to create slabs with arbitrary surfaces.
%
% * Label: play_hr
%
%% Description of the Function:
%%
%% Usage:
%
% * [hrTB_super,Rm_super]=supercell_hr(Ns,Rm,H_xyz,orbital_init)
% * [hrTB_super,Rm_super]=supercell_hr(Ns,Rm,H_xyz)
%% Input:
%
% # input1:
% # input2:
% # input3:
%
%% Output:
%
% # fig:
% # ax:
%
%% example:
%
%% Note:
%
%  Take advantage of the scope of application of the function.
%
%% Change log
%
% * Document Date: 2020/12/08
% * Creation Date: 2020/12/08
% * Last updated : 2020/12/08
%
%% Copyright
%
% * parkman
% * <parkman@buaa.edu.cn>
%
%% Source Code :
%
function [hrTB_super,sc_orb]=supercell_hr(Ns,orbital_init,Rm,H_hr)
    %--------  init  --------
    import vasplib_tool.*
    %--------  narg  --------
    if nargin < 4
        H_hr  = Readhr('wannier90_hr.dat',Rm);
    end
    if nargin  < 3
        POSCAR_read;
    end
    if nargin < 2
        % if exist wannier90_xyz.dat
        if ~exist(filename,'wannier90_xyz.dat')
            fprint('Wannier center from : %s\n',filename);
        else
            inputchar = input('Copy POSCAR site(y) or Set to Zero point(n)');
            if strcmp('y',inputchar)

            else
                WAN_NUM=length(H_hr(1).Hnum);
                orbital_init =zeros(WAN_NUM,3);
            end

        end
    end
    if nargin < 1
        Ns = eye(3);
    end
    %--------  init  --------
    use_sc_red_lat=Ns;
    %--------  check  --------
    % checks on super-lattice Ns
    if ~(Ns == round(Ns))
    	error("\n\nsc_red_lat array elements must be integers");
    end
    if det(Ns)<1.0E-6
        error("\n\nSuper-cell lattice vectors length/area/volume too close to zero, or zero.");
    end
    if det(Ns)<0.0
        error("\n\nSuper-cell lattice vectors need to form right handed system.");
    end
    %--------  init  --------
    [NRPTS,~] =  size(H_hr);
    Ns = round(Ns);
    V= abs(det(Ns));% intger forcely
    WAN_NUM=length(H_hr(1).Hnum);
    OUT_WAN_NUM = WAN_NUM*V ;
    
    % conservative estimate on range of search for super-cell vectors
    max_R=max(abs(Ns))*3;
    % candidates for super-cell vectors
    % this is hard-coded and can be improved!
    sc_cands = [];
    for i = -max_R(1):max_R(1)+1
        for j = -max_R(2):max_R(2)+1
            for k = -max_R(3):max_R(3)+1
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
        tmp_red=to_red_sc(sc_cands(ivec,:),Ns);
        % check if in the interior
        inside = 1;
        for it = 1:length(tmp_red)
            t = tmp_red(it);
            if t <= -1.0*eps_shift | t>1.0-eps_shift
                inside=0;
            end
        end
        if inside == 1
            %                 disp(tmp_red);
            %         disp(sc_cands(ivec,:));
            sc_vec=[sc_vec;sc_cands(ivec,:)];
            
        end
    end
    % number of times unit cell is repeated in the super-cell
    [num_sc,~] = size(sc_vec);
    %disp(num_sc);
    % check that found enough super-cell vectors
    if round(round(abs(det(Ns)))) ~= num_sc
        error("\n\nSuper-cell generation failed! Wrong number of super-cell vectors found.");
    end
    % -------- bug-fix ----------
    % sym instead
    Ns = sym(Ns);
    sc_vec = sym(sc_vec);
    % cartesian vectors of the super lattice
    % sc_cart_lat=Ns*Rm; %must Ns Rm
    % orbitals of the super-cell tight-binding model
    sc_orb=[];
    for icur_sc_vec = 1:num_sc % go over all super-cell vectors
        cur_sc_vec = sc_vec(icur_sc_vec,:);
        for iorb = 1:length(orbital_init) % go over all orbitals
            orb = orbital_init(iorb,:);            
            % shift orbital and compute coordinates in
            % reduced coordinates of super-cell
            sc_orb=[sc_orb;to_red_sc(orb+cur_sc_vec,Ns)];
        end
    end
    % create super-cell tb_model object to be returned
    %OUT_H_hr_ = struct('seq',[],'vector',[],'Degen',[],'key',[],'nokey',[],'Hstr',[],'Hsym',[],'Hcoe',[],'Hnum',[]);
    % OUT_H_xyz = repmat(H_xyz_ ,[NRPTS,1]);    % Hamiltonian of every cell;
    OUT_H_hr = zero_H_xyz(H_hr,OUT_WAN_NUM);
    fprintf('Search done; begin set hoppings\n');
    fprintf('We can improve the perfomance later\n');
    t1=clock;
    %set hopping terms
    for icur_sc_vec = 1:num_sc % go over all super-cell vectors
        cur_sc_vec = sc_vec(icur_sc_vec,:);
        for ih = 1:NRPTS % go over all hopping terms of the original model
            temp_i_H_xyz = H_hr(ih);
            % lattice vector of the hopping
            ind_R = temp_i_H_xyz.vector;            
            % super-cell component of hopping lattice vector
            % shift also by current super cell vector
            sc_part=floor(to_red_sc(ind_R+cur_sc_vec,Ns)); % round down!
            %disp(sc_part);
            %find remaining vector in the original reduced coordinates
            orig_part=ind_R+cur_sc_vec-sc_part*use_sc_red_lat;
            %disp(orig_part);
            
            % remaining vector must equal one of the super-cell vectors
            pair_ind = 9999;
            for jcur_sc_vec = 1:num_sc
                pair_sc_vec = sc_vec(jcur_sc_vec,:);
                if pair_sc_vec==orig_part
                    if pair_ind ~= 9999
                        error("\n\nFound duplicate super cell vector!");
                    end
                    pair_ind=jcur_sc_vec;
                end
            end
            if pair_ind==9999
                disp(ind_R);
                disp(orig_part);
                disp('Cant find sc');
                %continue;
                % The bug from numerical issue! we will use sym instead
                error("\n\nDid not find super cell vector!");
            end
            %disp(pair_ind);
            % index of "from" and "to" hopping indices
            for  i = 1:WAN_NUM
                for j = 1:WAN_NUM
                    amp = H_hr(ih).Hnum(i,j);
                    hi= i + (icur_sc_vec-1)*WAN_NUM ;
                    %disp(hi);
                    hj= j + (pair_ind-1)*WAN_NUM ;
                    %disp(hj);
                    OUT_H_hr = set_hop(amp,hi,hj,double(sc_part),OUT_H_hr,'add');

                end
            end
            fprintf("Generate process: SUPERCELL(%d,%d) NRPT(%d,%d) RUNINGTIME: %f s.\n",...
                icur_sc_vec,num_sc,ih,NRPTS,etime(clock,t1));
        end
    end
    hrTB_super = OUT_H_hr;
end


function To_red_sc=to_red_sc(red_vec_orig ,Ns)
    To_red_sc=red_vec_orig*inv(Ns);
    %To_red_sc= To_red_sc';
end