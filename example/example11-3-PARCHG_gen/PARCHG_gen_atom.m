function [ PARCHGCAR ,WF_COLORCAR ]= PARCHG_gen_atom(orb_list,WaveFunc,PARCHG_file,HIST,Rm,element_table)
% PARCHG gen
% usage: PARCHG_gen(orb_list,WaveFunc,PARCHG_file,curve_R,Rm,sites,Atom_name,Atom_num )
%        PARCHG_gen(orb_list,WaveFunc,PARCHG_file,curve_R)
%    a    PARCHG_gen(orb_list,WaveFunc,PARCHG_file)
%        PARCHG_gen(orb_list,WaveFunc)
% PARCHG contains the partial charge densities.
% This file contains
%   the lattice vectors,
%	atomic coordinates,
%	the total charge density multiplied by the volume { \rho (r)*V_{{{\rm {cell}}}}}\rho (r)*V_{{{\rm {cell}}}} on the fine FFT-grid (NG(X,Y,Z)F),
% and the PAW one-center occupancies.
% WRITE(IU,FORM) (((C(NX,NY,NZ),NX=1,NGXF),NY=1,NGYF),NZ=1,NGZF)
% Note that the real-space mesh (NX,NY,NZ) is uniform and is spanned by the primitive lattice vectors
% (a,b,c)defined in the POSCAR and can read explicitly
% \left(N_{x}, N_{y}, N_{z}\right) \triangleq \frac{N_{x}-1}{N_{G X P}} \vec{a}+\frac{N_{y}-1}{N_{G Y F}} \vec{b}+\frac{N_{z}-1}{N_{G Z F}} \vec{c}
%% naigin
if nargin < 3
    disp('glance mode')
    %PARCHG_file = "PARCHG";
    mode = 'glance';
else
    mode = 'PARCHG.cif';
end

if nargin <4 
    HIST = -1;
end

if nargin < 5
    if exist('POSCAR','file')
        [Rm,~,~,~,element_table]=POSCAR_readin('POSCAR','vasp');
    else
        error('POSCAR need！');
    end
end

[Nwave,Nlist] = size(WaveFunc);
Norb = size(orb_list,1);
if  Norb ~= Nwave
    error('Orbital list length is not equal to WaveFunc');
end
%% First map atom color to hue
rgb_ = [element_table.r,element_table.g,element_table.b];
Hsv_ = rgb2hsv(rgb_);
Hue_ = Hsv_(:,1);
%% Second reshape WaveFunc to Hue
WF_square_list = zeros(Norb,1);
for i = 1:Norb
    temp_C = 0;
    for j =1:Nlist
        temp_amp = WaveFunc(i,j);
        if temp_amp*temp_amp' < 1e-3
            temp_amp =0;
        end
        temp_C=temp_C+temp_amp*temp_amp';
    end
    WF_square_list(i) = temp_C;
end

WF_COLORCAR = normalize(WF_square_list,'range');
% hist ?
if HIST > 0
    HIST_frac = 1 /  HIST;
    WF_COLORCAR  = ceil(WF_COLORCAR/HIST_frac);
    WF_COLORCAR  = normalize(WF_COLORCAR ,'range')/2;
end

%% Third Gen WF plot information

namelist{Norb} = "";
% bug fix 
[~,orb_list(1,4)] = min(abs(Hue_-WF_COLORCAR(1)));
Recommand_list = [            orb_list(1,4),...
    WF_COLORCAR(1)];
%
for i = 1:Norb
    [~,orb_list(i,4)] = min(abs(Hue_-WF_COLORCAR(i)));
    if ~ismember(orb_list(i,4),Recommand_list(:,1))
        Recommand_list = [Recommand_list;...
            orb_list(i,4),...
            WF_COLORCAR(i)] ;
    end
    namelist{i} = element_table.atom_symbol(orb_list(i,4));
end
PARCHGCAR = orb_list;

%% Gen
if strcmp(mode,'glance')
    disp(Norb);
    %[Nwave,Nlist] = size(WaveFunc);
    if  Norb ~= Nwave
        error('Orbital list length is not equal to WaveFunc');
    end
    figure();
    for i = 1 : Norb
        %plot(PS(j,1),PS(j,2),'ro','MarkerSize',(Z1(j)'*Z1(j) + Z2(j)'*Z2(j) + Z3(j)'*Z3(j) + Z4(j)'*Z4(j) + Z5(j)'*Z5(j) + Z6(j)'*Z6(j) +0.001)*200,'MarkerFaceColor','r');
        if WF_COLORCAR(i)>0
            plot3(orb_list(i,1),orb_list(i,2),orb_list(i,3),'ro','MarkerSize',WF_COLORCAR(i)*100,'MarkerFaceColor',rgb_(orb_list(i,4),:));
        else
            plot3(orb_list(i,1),orb_list(i,2),orb_list(i,3),'ko','MarkerSize',1,'MarkerFaceColor','k');
        end
        hold on
    end
    

    %%
elseif strcmp(mode,'PARCHG.cif')
    %% check
    disp(PARCHG_file);
    %% init
    title = "# CIF Gen by MATLAB";
    filename = PARCHG_file;
    fileID = fopen(filename,'w');
    fprintf(fileID,"%s\n",title);
    fprintf(fileID,'#email: parkman@buaa.edu.cn\n');
    fprintf(fileID,'## Date: %s\n\n',date);
    %
    fprintf(fileID,'# %s\n\n',filename);
    
    % CRYSTAL INFOR
    R_struct = Rm2abc(Rm);
    fprintf(fileID,'%-20s %13.6f\n',"_cell_length_a",    R_struct.a);
    fprintf(fileID,'%-20s %13.6f\n',"_cell_length_b",   R_struct.b);
    fprintf(fileID,'%-20s %13.6f\n',"_cell_length_c", R_struct.c);
    fprintf(fileID,'%-20s %13.6f\n',"_cell_angle_alpha", R_struct.alpha);
    fprintf(fileID,'%-20s %13.6f\n',"_cell_angle_beta",  R_struct.beta);
    fprintf(fileID,'%-20s %13.6f\n\n',"_cell_angle_gamma ",R_struct.gamma);
    
    fprintf(fileID,'loop_\n');
    fprintf(fileID,'_space_group_symop.magn_id\n');
    fprintf(fileID,'_space_group_symop.magn_operation_xyz\n');
    fprintf(fileID,'1   x,y,z,+1\n\n');
    % SITES INFORMATION
    fprintf(fileID,'loop_\n');
    fprintf(fileID,'_atom_site_label\n');
    fprintf(fileID,'_atom_site_type_symbol\n');
    fprintf(fileID,'_atom_site_fract_x\n');
    fprintf(fileID,'_atom_site_fract_y\n');
    fprintf(fileID,'_atom_site_fract_z\n');
    fprintf(fileID,'_atom_site_occupancy\n');
    
    for i =1:Norb
        fprintf(fileID,'%4s %3s %7.5f %7.5f %7.5f %6.4f\n',...
            namelist{i}{1},namelist{i}{1},...
            orb_list(i,1),orb_list(i,2),orb_list(i,3),...
            1.0000);
    end
    fprintf(fileID,'\n');
    fclose(fileID);
    %% color radius needed 
    Recommand_list = sortrows(Recommand_list ,2,'descend');
    for i =  1: size(Recommand_list,1)
        tempcell = element_table.atom_symbol(Recommand_list(i,1));
        fprintf("The %d st componet: %s Radius_suggest: %5.3f \n",...
            i,tempcell{1},Recommand_list(i,2)*8+1);
    end
end
end

function R_struct = Rm2abc(Rm)
    R_struct.a = norm(Rm(1,:));
    R_struct.b = norm(Rm(2,:));
    R_struct.c = norm(Rm(3,:));
    R_struct.alpha = acos(dot(Rm(2,:),Rm(3,:))/(norm(R_struct.b)*norm(R_struct.c)))/pi*180;
    R_struct.beta  = acos(dot(Rm(3,:),Rm(1,:))/(norm(R_struct.c)*norm(R_struct.a)))/pi*180;
    R_struct.gamma = acos(dot(Rm(1,:),Rm(2,:))/(norm(R_struct.a)*norm(R_struct.b)))/pi*180;
end 


