%% Mapping Method : J1 VASP输入文件自动生成
%% 定义环境变量

clear;
% 外部程序地址

%% System environment
vaspkit_run = '/Users/parkman/Documents/soft/vaspkit.1.2.1/bin/vaspkit';
phonopy_run = '/Users/parkman/Documents/TOOLs/miniconda2/envs/my_pymatgen/bin/phonopy';
% 自定义变量初始化

Accuracy = 3 ;
level_cut  = 1;  % only J1

U_list = [1,3];
General_MagAmp  = 7;
mag_element_seq_list =  56:72; % rare_earth add more!!
search_range = [level_cut+2, level_cut+2, level_cut+2];

neckname = 'TMGO';
POSCAR_name = 'POSCAR_origin';

WORKSPACE=pwd;
% define U_struct
%U_list = [0,1,3,5,7,9]; % which U list
%% POSCAR处理
% 磁性ion的处理
% 分离 mag_ion 的POSCAR

[Rm,sites,Atom_name,Atom_num,elements]=POSCAR_readin(POSCAR_name,'vasp');
[mag_ion_seq,sites,Atom_name,Atom_num] = Magion_detect(mag_element_seq_list,Rm,sites,Atom_name,Atom_num,elements);
% 读取 mag_ion 的POSCAR

[Rm,sites_mag,Atom_name_mag,Mag_ion_num,~]=POSCAR_readin('POSCAR_mag');

[~,nn_store_smart,~,~]=nn_smart(Rm,sites_mag,search_range,Accuracy,10);
% 给出合适地supercell 矩阵
%% 
% * 选定命运之子

ion_Rf_list = [[sites_mag.rc1]'.^2+[sites_mag.rc2]'.^2+[sites_mag.rc3]'.^2];
[~,The_God_ion ]= min(ion_Rf_list);
%% 
% * 查看其 mag nn

The_God_ion_nn = nn_store_smart(The_God_ion,:).nn;
[~,sort_seq] = sort([The_God_ion_nn.nn_level]);
The_God_ion_nn = The_God_ion_nn(sort_seq); 
%% 
% * 只需要在 level_cut 的这一层 统计 R_vector 即可确定 最小supercell值

select_nn = [The_God_ion_nn.nn_level] == level_cut  ; 
R_vector_list = [The_God_ion_nn(select_nn).R_vector];
R_vector_list = abs(reshape(R_vector_list,3,length(R_vector_list)/3).');
[~,R_vector_label] = min(sum(R_vector_list,2));
Nslab = R_vector_list(R_vector_label,:)+1;
Ns = diag(Nslab);
V= abs(det(Ns));
%% 
% * 生成POSCAR_super

supercell(Ns,Rm,sites,Atom_name,Atom_num);
supercell(Ns,Rm,sites_mag,Atom_name_mag,Mag_ion_num,[0,0,0],'POSCAR_super_mag');
%% FM - AFM 法两个对位原子的选择
% 读取POSCAR_super 和 POSCAR_super_mag

[Rm,sites,Atom_name,Atom_num,~]=POSCAR_readin('POSCAR_super','vasp');
[~,sites_mag,~,~,~]=POSCAR_readin('POSCAR_super_mag','vasp');
[~,nn_store_smart,~,~]=nn_smart(Rm,sites_mag,[1,1,1],Accuracy,10);
% [mag_ion_seq] = Magion_detect(mag_element_seq_list,Rm,sites,Atom_name,Atom_num,elements);
% 生成 MAGMOM_list
%% 
% * 选定命运之子

ion_Rf_list = [[sites_mag.rc1]'.^2+[sites_mag.rc2]'.^2+[sites_mag.rc3]'.^2];
[~,The_God_ion ]= min(ion_Rf_list);
The_God_ion_nn = nn_store_smart(The_God_ion,:);
%% 
% * 构造 for 循环

for i = 1:level_cut
    
%% 
% * 选定撒旦之子

    R_vector_list = [];
    for j = 1:length(The_God_ion_nn)
        select_nn = [The_God_ion_nn(j).nn.nn_level] == i;
        temp_R_list = [The_God_ion_nn(j).nn(select_nn).R_vector];
        R_vector_list = [R_vector_list;...
            reshape( temp_R_list ,3,...
            length( temp_R_list)/3).';];
    end
    %R_vector_list = reshape(R_vector_list,3,length(R_vector_list)/3).';
    R_vector_norm_list =R_vector_list (:,1).^2 +R_vector_list (:,2).^2 +R_vector_list (:,3).^2 ;
    [~,The_Satan_ion ]= max(R_vector_norm_list );
%% 
% # 生成

    MAGMOM_list{i} = ones(4,length(sites_mag));
    % state 1
    MAGMOM_list{i}(1,The_God_ion) =  1; MAGMOM_list{i}(1,The_Satan_ion ) =   1;
    % state 2
    MAGMOM_list{i}(2,The_God_ion) =  1; MAGMOM_list{i}(2,The_Satan_ion ) =  -1;
end
%% 生成 INCAR
% 定义INCAR 变量

NBANDS = -1;
ion_kinds = length(Atom_name);
ions_num = sum(Atom_num) ;
U_struct = zeros(ion_kinds ,3);
for i = 1:ion_kinds 
    U_struct(i,1) = -1;
    U_struct(i,2) =  0;
    U_struct(i,3) =  0;
end
U_struct(mag_ion_seq,1) = 3;%f
% 生成 INCAR
% 创建INCAR_collection

%% Gen INCAR
mkdir('INCAR_collection');
cd('INCAR_collection');
% 遍历U

for U = U_list
    U_struct(mag_ion_seq,2) = U;
    U_dir = "U_EQ_"+string(U);
    mkdir(U_dir);
    cd(U_dir );
    
% 遍历Jn

    for Jn =  1:level_cut
        mkdir("J"+string(Jn));
        cd("J"+string(Jn));
        
% 遍历 w / wo SOC nCOL/COL

        for BIG_dir_CELL  = ["nosoc_col","nosoc_ncol","soc_ncol"]
            BIG_dir = BIG_dir_CELL{1};
            switch BIG_dir
                case 'nosoc_col'
                    LNOCOL = '.FALSE.';
                    LSOC = '.FALSE.'  ;
                    ICHARG = -1    ;
                    ISTART = 1     ;
                    ISPIN  = 2     ;
                case 'nosoc_ncol'
                    LNOCOL = '.TRUE.';
                    LSOC = '.FALSE.'  ;
                    ICHARG = -1    ;
                    ISTART = 1     ;
                    ISPIN  = 2     ;
                case 'soc_ncol'
                    LNOCOL = '.TRUE.';
                    LSOC = '.TRUE.'  ;
                    ICHARG = 11    ;
                    ISTART = 1     ;
                    ISPIN  = 1     ;
            end
            mkdir(BIG_dir);
            cd(BIG_dir);
            switch BIG_dir
                case 'nosoc_col'
                    CS_COUNT = 0;
                    %% State 1
                    CS_COUNT = CS_COUNT +1;
                    COL_SYSTEM(CS_COUNT).systemname = 'DirectMap-State1';
                    mag_list = MAGMOM_list{Jn}(1,:);
                    MagAmp   = General_MagAmp;
                    magtype  = 'AFM';
                    COL_SYSTEM(CS_COUNT).MAGMOM = MAGMOM_gen(mag_list,MagAmp,magtype,ions_num);
                    %% State 2
                    CS_COUNT = CS_COUNT +1;
                    COL_SYSTEM(CS_COUNT).systemname = 'DirectMap-State2';
                    mag_list = MAGMOM_list{Jn}(2,:);
                    MagAmp   = General_MagAmp;
                    magtype  = 'AFM';
                    COL_SYSTEM(CS_COUNT).MAGMOM = MAGMOM_gen(mag_list,MagAmp,magtype,ions_num);
                    %% GEN
                    for i = 1:length(COL_SYSTEM)
                        INCARsystem = neckname+"-"+BIG_dir+"-"+U_dir+"-"+COL_SYSTEM(i).systemname;
                        INCAR_gen(COL_SYSTEM(i).systemname,INCARsystem,COL_SYSTEM(i).MAGMOM,...
                            U_struct,ISPIN,ISTART,ICHARG,LSOC,LNOCOL,NBANDS);
                    end
                    clear('COL_SYSTEM');
                case {'nosoc_ncol','soc_ncol'}
                    NCS_COUNT = 0;
                    % FMz
                    NCS_COUNT = NCS_COUNT +1;
                    NCOL_SYSTEM(NCS_COUNT).systemname = 'DirectMap-MAE-z';
                    mag_list = ones(ions_num/V,1);
                    MagAmp   = General_MagAmp;
                    magtype  = 'FMz';
                    NCOL_SYSTEM(NCS_COUNT).MAGMOM = MAGMOM_gen(mag_list,MagAmp,magtype,ions_num*3/V);
                    %% FMx
                    NCS_COUNT = NCS_COUNT +1;
                    NCOL_SYSTEM(NCS_COUNT).systemname = 'DirectMap-MAE-x';
                    mag_list = ones(ions_num/V,1);
                    MagAmp   = General_MagAmp;
                    magtype  = 'FMx';
                    NCOL_SYSTEM(NCS_COUNT).MAGMOM = MAGMOM_gen(mag_list,MagAmp,magtype,ions_num*3/V);
                    %% z
                    % State 1
                    NCS_COUNT = NCS_COUNT +1;
                    NCOL_SYSTEM(NCS_COUNT).systemname = 'DirectMap-State1-z';
                    mag_list = MAGMOM_list{Jn}(1,:);
                    MagAmp   = General_MagAmp;
                    magtype  = 'AFMz';
                    NCOL_SYSTEM(NCS_COUNT).MAGMOM = MAGMOM_gen(mag_list,MagAmp,magtype,ions_num*3);
                    % State 2
                    NCS_COUNT = NCS_COUNT +1;
                    NCOL_SYSTEM(NCS_COUNT).systemname = 'DirectMap-State2-z';
                    mag_list = MAGMOM_list{Jn}(2,:);
                    MagAmp   = General_MagAmp;
                    magtype  = 'AFMz';
                    NCOL_SYSTEM(NCS_COUNT).MAGMOM = MAGMOM_gen(mag_list,MagAmp,magtype,ions_num*3);
                    %% x
                    % State 1
                    NCS_COUNT = NCS_COUNT +1;
                    NCOL_SYSTEM(NCS_COUNT).systemname = 'DirectMap-State1-x';
                    mag_list = MAGMOM_list{Jn}(1,:);
                    MagAmp   = General_MagAmp;
                    magtype  = 'AFMx';
                    NCOL_SYSTEM(NCS_COUNT).MAGMOM = MAGMOM_gen(mag_list,MagAmp,magtype,ions_num*3);
                    % State 2
                    NCS_COUNT = NCS_COUNT +1;
                    NCOL_SYSTEM(NCS_COUNT).systemname = 'DirectMap-State2-x';
                    mag_list = MAGMOM_list{Jn}(2,:);
                    MagAmp   = General_MagAmp;
                    magtype  = 'AFMx';
                    NCOL_SYSTEM(NCS_COUNT).MAGMOM = MAGMOM_gen(mag_list,MagAmp,magtype,ions_num*3);
                    %% GEN
                    for i = 1:length(NCOL_SYSTEM)
                        INCARsystem = neckname+"-"+BIG_dir+"-"+U_dir+"-"+NCOL_SYSTEM(i).systemname;
                        INCAR_gen(NCOL_SYSTEM(i).systemname,INCARsystem,NCOL_SYSTEM(i).MAGMOM,...
                            U_struct,ISPIN,ISTART,ICHARG,LSOC,LNOCOL,NBANDS);
                    end
            end
            cd('..');
        end
        cd('..');
    end
    cd(WORKSPACE+"/"+'INCAR_collection');
end
cd(WORKSPACE);
%% 准备POTCAR 和 POSCAR
% 生成POTCAR

%% POTCAR_gen
copyfile('POSCAR_mag','POSCAR');
if exist('POTCAR','file')
    delete('POTCAR');
end
system(strcat("(echo 104;echo ",Atom_name_mag,')|',vaspkit_run,'>>vaspkit.log'));
movefile('POTCAR',Atom_name_mag);

copyfile('POSCAR_others','POSCAR');
% delete('POTCAR');
system(strcat("(echo 103)|",vaspkit_run,'>>vaspkit.log'));
movefile('POTCAR','POT_others');
system(strcat("cat ",Atom_name_mag,' POT_others > POTCAR'));
[~,tmp_str] = system('grep VRHFIN POTCAR');
fprintf('%s\n',tmp_str);
% 复制 POSCAR_super -> POSCAR

copyfile('POSCAR_super','POSCAR');
% POSCAR_name = 'POSCAR_origin';
%% 生成WORKDIR

%% Gen workdir
mkdir('WORKDIR');
system('cp -r INCAR_collection/* WORKDIR/');
cd('WORKDIR');
% 遍历U

for U = U_list
    U_dir = "U_EQ_"+string(U);
    cd(U_dir);
    %copyfile(WORKSPACE+"/"+'RUN_one_scf.m','RUN_one_scf.m');
% 遍历Jn

    for Jn =  1:level_cut
        cd("J"+string(Jn));        
% 遍历 w / wo SOC nCOL/COL

        for BIG_dir_CELL  = ["nosoc_col","nosoc_ncol","soc_ncol"]
            BIG_dir = BIG_dir_CELL{1};
            cd(BIG_dir);
            %
            INCAR_system = dir();
            for i = 1:length(INCAR_system)
                system_neckname = INCAR_system(i).name;
                if isfile(system_neckname)
                    system_neckname_dir = system_neckname + "_scf";
                    if strcontain(system_neckname,"DirectMap")
                        fprintf('The system is : %s J%d %s %s\n',U_dir,Jn,BIG_dir,system_neckname);
                        mkdir(system_neckname_dir);
                        cd(system_neckname_dir);
                        copyfile(WORKSPACE+"/POTCAR",'POTCAR');
                        copyfile("../"+system_neckname,'INCAR');
                        if strcontain(system_neckname,"State")
                            copyfile(WORKSPACE+"/POSCAR_super",'POSCAR');
                        elseif  strcontain(system_neckname,"MAE")
                            copyfile(WORKSPACE+"/POSCAR_origin",'POSCAR');
                        end
                        %                 ls;
%                         fprintf('\n');
                        cd('..')
                    end
                end
            end
            %
            fprintf('_______________________________________________________________\n');
            fprintf('***************************************************************\n');
            cd('..');
        end
        fprintf('***************************************************************\n');
        cd('..');
    end
    cd(WORKSPACE+"/"+'WORKDIR');
end
cd(WORKSPACE);
fprintf('All done!\n');