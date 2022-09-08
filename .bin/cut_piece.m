function [OUT_H_xyz,fin_orb] = cut_piece(H_xyz,repeatnum,fin_dir,glue_edges,orbital_init,mode)
%%         Constructs a (d-1)-dimensional tight-binding model out of a d-dimensional one by repeating the unit cell a given number of times along one of the periodic lattice vectors.
%          The real-space lattice vectors of the returned model are the same as those of
%          the original model; only the dimensionality of reciprocal space is reduced.
%
%         :param repeatnum: How many times to repeat the unit cell.
%
%         :param fin_dir: Index of the real space lattice vector along
%           which you no longer wish to maintain periodicity.

%         :param orbital_init

%         :param glue_edges: Optional boolean parameter specifying whether to
%           allow hoppings from one edge to the other of a cut model.
%
%         :returns:
%           * **fin_model** -- Object of type
%             Out_H_xyz representing a cutout
%             tight-binding model.
%             Orbitals in *fin_model* are
%             numbered so that the i-th orbital of the n-th unit
%             cell has index i+norb*n (here norb is the number of
%             orbitals in the original model).
%
%         Example usage::
%           Out_H_xyz = cut_piece(H_xyz,repeatnum,fin_dir,glue_edgs,orbital_init)

%           A = tb_model(3, 3, ...)
%           % Construct two-dimensional model B out of three-dimensional
%           % model A by repeating model along second lattice vector ten times
%           B = A.cut_piece(10, 1)
%           % Further cut two-dimensional model B into one-dimensional model
%           % A by repeating unit cell twenty times along third lattice
%           % vector and allow hoppings from one edge to the other
%           C = B.cut_piece(20, 2, glue_edgs=True)
%
if nargin < 6
    mode = 'norm-mat';
    [NRPTS,~]=size(H_xyz);
    WAN_NUM=length(H_xyz(1).Hnum);
end
mode_mode = mode(1:4);
mode_cal = mode(6:8);
%%
if exist('POSCAR','file')
    Ns = [1 0 0;0 1 0;0 0 1];
    Ns(fin_dir,:) = Ns(fin_dir,:) * repeatnum;
    POSCAR_read;
    fin_dir_list = [0 0 0];
    if strcmp(mode_mode,'disk')
        fin_dir_list(fin_dir) = 1;
    end
    [Rm,~] = supercell(Ns,Rm,sites,Atom_name,Atom_num,fin_dir_list,'POSCAR_super_fin');
end

%% init
if strcmp(mode_mode,'disk')
    [NRPTS,~]=size(H_xyz);
    WAN_NUM=length(H_xyz(1).Hnum);
elseif strcmp(mode_mode,'hrln')
    Hnum_list = H_xyz.HnumL ;
    vector_list = H_xyz.vectorL ;
    [NRPTS,~]=size(vector_list);
    %    WAN_NUM=length(Hnum_list(:,:,1));
    WAN_NUM_2 = length(Hnum_list(:,:,1));
    WAN_NUM=sqrt(WAN_NUM_2);
elseif strcmp(mode_mode,'hrsn')
    Hnum_list = H_xyz.HnumL ;
    vector_list = H_xyz.vectorL ;
    [NRPTS,~]=size(vector_list);
    WAN_NUM=length(Hnum_list{1});
end
%% nargin
if nargin < 5
    orbital_init =zeros(WAN_NUM,3);
end

if nargin < 4
    glue_edges = 'true';
end


% check value of orbital_init
[norb ,~] =size(orbital_init);
if norb ~= WAN_NUM
    error("\n\nOribital_init is wrong,please give a right orbital init or just forget this parm!");
end
% check value of num
if repeatnum<1
    error("\n\nArgument num must be positive!");
end
if repeatnum == 1 & strcmp(glue_edges,'true')
    error("\n\nCan't have num==1 and glueing of the edges!");
end
num  = round(repeatnum);

% generate orbitals of a finite model
fin_orb = [];
%onsite = [] % store also onsite energies
for inum = 1:num % go over all cells in finite direction
    for j = 1:WAN_NUM  % go over all orbitals in one cell
        % make a copy of j-th orbital
        orb_tmp=orbital_init(j,:);
        % change coordinate along finite direction ; fractional
        orb_tmp(fin_dir)= (orb_tmp(fin_dir) + double(inum-1))/num;
        % add to the list
        fin_orb=[fin_orb;orb_tmp];
        % do the onsite energies at the same time
    end
end

% generate periodic directions of a finite model
OUT_WAN_NUM = WAN_NUM*num ;
%OUT_WAN_NUM_2 = OUT_WAN_NUM^2 ;
% init
if strcmp(mode_mode,'disk')
    % generate object of tb_model type that will correspond to a cutout
    %OUT_H_xyz_ = struct('seq',[],'vector',[],'Degen',[],'key',[],'nokey',[],'Hstr',[],'Hsym',[],'Hcoe',[],'Hnum',[]);
    OUT_H_xyz = zero_H_xyz(H_xyz,OUT_WAN_NUM);
elseif strcmp(mode_mode,'hrsn')
    for ih =1:NRPTS
        OUT_Hnum_list{ih} =sparse(OUT_WAN_NUM,OUT_WAN_NUM);
    end
    %OUT_Hnum_list = zeros(OUT_WAN_NUM,OUT_WAN_NUM,NRPTS);
    OUT_vector_list = vector_list ;
    OUT_H_xyz.vectorL  = OUT_vector_list ;
end

%% three type of hr
    if strcmp(mode_mode,'hrsn')
        disp('sparse mode');
        for ih = 1:NRPTS % go over all NRPTS
            % lattice vector of the hopping
            ind_R = OUT_vector_list(ih,:);
            jump_fin=ind_R(fin_dir);
            % store by how many cells is the hopping in finite direction
            %temp_Hnum =sparse(OUT_WAN_NUM,OUT_WAN_NUM);
            fprintf('NPRT: %d s \n',ih);
            % speed up more and more !!
            for  i = 1:WAN_NUM
                for j = 1:WAN_NUM
                    % amplitude of the hop is the same
                    amp = Hnum_list{ih}(i,j);
                    if norm(amp) > 0
                        for icur_sc_vec = 1:num % go over all cells in finite direction mini
                            hi= i + (icur_sc_vec-1)*WAN_NUM ;
                            %disp(hi);
                            hj= j + (icur_sc_vec+jump_fin-1)*WAN_NUM ;
                            %disp(hj);
                            % decide whether this hopping should be added or not
                            to_add=1;
                            %disp(hj);
                            % if edges are not glued then neglect all jumps that spill out
                            if strcmp(glue_edges,'false')
                                if hj <= 0 | hj > OUT_WAN_NUM
                                    to_add=0;
                                end
                                % if edges are glued then do mod division to wrap up the hopping
                            else
                                hj= mod(hj,OUT_WAN_NUM);
                                if hj ==0
                                    hj = OUT_WAN_NUM;
                                end
                            end
                            if to_add == 1
                                OUT_Hnum_list{ih}(hi,hj) = amp;
                                %OUT_Hnum_list(hi,hj,ih) = amp;
                                %                        OUT_Hnum_list(hi,hj,i) = amp;
                            end
                        end
                    end
                end
            end
            %OUT_Hnum_list(:,ih) = temp_Hnum;
        end
        %OUT_Hnum_list(:,:,ih) = temp_Hnum;
        %
        OUT_H_xyz.HnumL = OUT_Hnum_list;
    else
    end




if exist('POSCAR','file')
    % rebuild fin_orb
    Rmlength1 = norm (Rm(1,:));
    Rmlength2 = norm (Rm(2,:));
    Rmlength3 = norm (Rm(3,:));
    Rm_s_fin_add = [10*Rm(1,:)*fin_dir_list(1)/Rmlength1;...
        10*Rm(2,:)*fin_dir_list(2)/Rmlength2;...
        10*Rm(3,:)*fin_dir_list(3)/Rmlength3];
    Rm_s_fin = Rm + Rm_s_fin_add ;
    Rc_s_fin_add = [1/2, 1/2 ,1/2] ;
    Rr_s_fin_add = Rc_s_fin_add * Rm_s_fin_add;
    [nfinorb,~ ]= size(fin_orb);
    for  i = 1:nfinorb
        Rr_orb = fin_orb(i,:)*Rm;
        Rr_s_fin = Rr_orb + Rr_s_fin_add;
        Rc_s_fin = Rr_s_fin * inv(Rm_s_fin);
        fin_orb(i,:) = Rc_s_fin ;
    end
end

end

% function Mat_for =test(icur_sc_vec,jump_fin,WAN_NUM,OUT_WAN_NUM,glue_edges)
%     v  =1 ;
%     end