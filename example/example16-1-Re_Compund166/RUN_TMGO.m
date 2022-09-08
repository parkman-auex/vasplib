clear;
%% RUN TMGO
WORKSPACE=pwd;
%% System environment
vaspkit_run = '/home/soft/vaspkit.1.2.0/bin/vaspkit';
phonopy_run = '/home/soft/anaconda2/bin/phonopy';

%% init
POSCAR_init = 'POSCAR_init';
neckname = "TMGO";
General_MagAmp  = 7;
%% Reshape POSCAR
eval("!"+phonopy_run+' --tolerance 0.01 --symmetry -c '+POSCAR_init+...
    '>'+'SymInfor_'+POSCAR_init);
copyfile('PPOSCAR','POSCAR');
% read POSCAR
[Rm,sites,Atom_name,Atom_num,elements]=POSCAR_readin('POSCAR','vasp');
% make peridoic table
elements.Properties.RowNames = elements.atom_symbol;
%
fprintf('We have %d ions:\n',sum(Atom_num));
rare_earth_ion_seq  = -1;
for i = 1:length(Atom_name)
    tmpdata = table2array(elements(Atom_name(i),{'atom_number','n'}));
    element_seq = tmpdata(1);
    fprintf('    %d %s(%d) ions \n',Atom_num(i),Atom_name(i),element_seq);
    if element_seq > 56 && element_seq < 72
        if rare_earth_ion_seq >0
            error('Dont support multiple rare-earth element at present');
        end
        if Atom_num(i) > 1
            error('Dont support more rare-earth ions in unit cell  at present');
        end
        rare_earth_ion_seq = i;
        rare_earth_ion_site_seq = sum(Atom_num(1:rare_earth_ion_seq));
    end
end
if rare_earth_ion_seq == -1
    fprintf('do not detect re-earth element');
    error('!!');
else
    fprintf('Rare earth element: %s\n',Atom_name(rare_earth_ion_seq));
end
% distri
Atom_name_rare_earth = Atom_name(rare_earth_ion_seq);
Atom_num_rare_earth  = 1;
sites_rare_earth     = sites(rare_earth_ion_site_seq);
POSCAR_gen(Rm,sites_rare_earth,Atom_name_rare_earth,Atom_num_rare_earth,'POSCAR_rarearth');
Atom_name(rare_earth_ion_seq) = [];
Atom_num(rare_earth_ion_seq) = [];
sites(rare_earth_ion_site_seq) = [];
POSCAR_gen(Rm,sites,Atom_name,Atom_num,'POSCAR_others');
%POSCAR_gen(Rm,sites,Atom_name,Atom_num,filename);
%% POTCAR_gen
copyfile('POSCAR_rarearth','POSCAR');
if exist('POTCAR','file')
    delete('POTCAR');
end
eval("!"+'(echo 104;echo '+Atom_name_rare_earth+')|'+vaspkit_run+'>>vaspkit.log');
movefile('POTCAR',Atom_name_rare_earth);

copyfile('POSCAR_others','POSCAR');
% delete('POTCAR');
eval("!"+'(echo 103)|'+vaspkit_run+'>>vaspkit.log');
movefile('POTCAR','POT_others');
eval("!"+'cat '+Atom_name_rare_earth+' POT_others > POTCAR');
eval("!"+'grep VRHFIN POTCAR');
%% POSCAR final gen
Atom_name = [Atom_name_rare_earth,Atom_name];
Atom_num = [1,Atom_num];
sites = [sites_rare_earth,sites];
POSCAR_gen(Rm,sites,Atom_name,Atom_num,'POSCAR');
eval("!"+phonopy_run+' --tolerance 0.01 --symmetry -c '+'POSCAR'+...
    '>'+'SymInfor_'+'POSCAR');
copyfile('PPOSCAR','POSCAR');
copyfile('PPOSCAR','POSCAR_origin');
%% flat POSCAR
findir = [0,0,0];
[Rm,sites,Atom_name,Atom_num,~]=POSCAR_readin('POSCAR','vasp');
Ns_flat = [1 -1 0;0 1 -1;1 0 0];
supercell(Ns_flat,Rm,sites,Atom_name,Atom_num,findir,'POSCAR');
[Rm,sites,Atom_name,Atom_num,elements]=POSCAR_readin('POSCAR','vasp');
%% supercell_POSCAR
Ns_221 = [2 0 0;0 2 0;0 0 1];
Ns_112 = [1 0 0;0 1 0;0 0 2];
Ns_gen3gen31 = [2 1 0;1 2 0;0 0 1];
supercell(Ns_221,Rm,sites,Atom_name,Atom_num,findir,'POSCAR_221');
supercell(Ns_112,Rm,sites,Atom_name,Atom_num,findir,'POSCAR_112');
supercell(Ns_gen3gen31,Rm,sites,Atom_name,Atom_num,findir,'POSCAR_gen3gen31');
%% Gen INCAR
U_struct = zeros(length(Atom_name),3);
for i = 1:length(Atom_name)
    U_struct(i,1) = -1;
    U_struct(i,2) = 0;
    U_struct(i,3) = 0;
end
U_struct(1,1) = 3;
NBANDS = -1;

mkdir('INCAR_collection');
cd('INCAR_collection');
for U = [0,1,3,5,7,9]
    U_struct(1,2) = U;
    U_dir = "U="+string(U);
    mkdir(U_dir);
    cd(U_dir );
    for BIG_dir_CELL  = {"nosoc_col","nosoc_ncol","soc_ncol"}
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
                %% PM
                systemname = 'PM-origin';
                magtype = 'PM';
                ions_num = sum(Atom_num);
                mag_list = zeros(ions_num,1);
                MagAmp = 0;
                INCARsystem = neckname+"-"+BIG_dir+"-"+U_dir+"-"+systemname;
                MAGMOM = MAGMOM_gen(mag_list,MagAmp,magtype,ions_num);
                INCAR_gen(systemname,INCARsystem,MAGMOM,U_struct,ISPIN,ISTART,ICHARG,LSOC,LNOCOL,NBANDS);
                %% FM origin 112 221 gen3gen31
                % origin
                systemname = 'FM-origin';
                magtype = 'FM';
                ions_num = sum(Atom_num);
                mag_list = 1;
                MagAmp = General_MagAmp;
                INCARsystem = neckname+"-"+BIG_dir+"-"+U_dir+"-"+systemname;
                MAGMOM = MAGMOM_gen(mag_list,MagAmp,magtype,ions_num);
                INCAR_gen(systemname,INCARsystem,MAGMOM,U_struct,ISPIN,ISTART,ICHARG,LSOC,LNOCOL,NBANDS);
                % 112
                systemname = 'FM-112';
                magtype = 'FM';
                ions_num = sum(Atom_num)*2;
                mag_list = [1 1];
                MagAmp = General_MagAmp;
                INCARsystem = neckname+"-"+BIG_dir+"-"+U_dir+"-"+systemname;
                MAGMOM = MAGMOM_gen(mag_list,MagAmp,magtype,ions_num);
                INCAR_gen(systemname,INCARsystem,MAGMOM,U_struct,ISPIN,ISTART,ICHARG,LSOC,LNOCOL,NBANDS);
                % gen3gen31
                systemname = 'FM-gen3gen31';
                magtype = 'FM';
                ions_num = sum(Atom_num)*3;
                mag_list = [1 1 1];
                MagAmp = General_MagAmp;
                INCARsystem = neckname+"-"+BIG_dir+"-"+U_dir+"-"+systemname;
                MAGMOM = MAGMOM_gen(mag_list,MagAmp,magtype,ions_num);
                INCAR_gen(systemname,INCARsystem,MAGMOM,U_struct,ISPIN,ISTART,ICHARG,LSOC,LNOCOL,NBANDS);
                % 221
                systemname = 'FM-221';
                magtype = 'FM';
                ions_num = sum(Atom_num)*4;
                mag_list = [1 1 1 1];
                MagAmp = General_MagAmp;
                INCARsystem = neckname+"-"+BIG_dir+"-"+U_dir+"-"+systemname;
                MAGMOM = MAGMOM_gen(mag_list,MagAmp,magtype,ions_num);
                INCAR_gen(systemname,INCARsystem,MAGMOM,U_struct,ISPIN,ISTART,ICHARG,LSOC,LNOCOL,NBANDS);
                %% AFM 112 221 gen3gen31
                % 112
                systemname = 'AFM-112';
                magtype = 'AFM';
                ions_num = sum(Atom_num)*2;
                mag_list = [1 -1];
                MagAmp = General_MagAmp;
                INCARsystem = neckname+"-"+BIG_dir+"-"+U_dir+"-"+systemname;
                MAGMOM = MAGMOM_gen(mag_list,MagAmp,magtype,ions_num);
                INCAR_gen(systemname,INCARsystem,MAGMOM,U_struct,ISPIN,ISTART,ICHARG,LSOC,LNOCOL,NBANDS);
                % gen3gen31
                systemname = 'AFM-gen3gen31-upupdown';
                magtype = 'AFM';
                ions_num = sum(Atom_num)*3;
                mag_list = [1 1 -1];
                MagAmp = General_MagAmp;
                INCARsystem = neckname+"-"+BIG_dir+"-"+U_dir+"-"+systemname;
                MAGMOM = MAGMOM_gen(mag_list,MagAmp,magtype,ions_num);
                INCAR_gen(systemname,INCARsystem,MAGMOM,U_struct,ISPIN,ISTART,ICHARG,LSOC,LNOCOL,NBANDS);
                % 221-stripe1
                systemname = 'AFM-221-stripe1';
                magtype = 'AFM';
                ions_num = sum(Atom_num)*4;
                mag_list = [1 -1 -1 1];
                MagAmp = General_MagAmp;
                INCARsystem = neckname+"-"+BIG_dir+"-"+U_dir+"-"+systemname;
                MAGMOM = MAGMOM_gen(mag_list,MagAmp,magtype,ions_num);
                INCAR_gen(systemname,INCARsystem,MAGMOM,U_struct,ISPIN,ISTART,ICHARG,LSOC,LNOCOL,NBANDS);
                % 221-stripe2
                systemname = 'AFM-221-stripe2';
                magtype = 'AFM';
                ions_num = sum(Atom_num)*4;
                mag_list = [1 1 -1 -1];
                MagAmp = General_MagAmp;
                INCARsystem = neckname+"-"+BIG_dir+"-"+U_dir+"-"+systemname;
                MAGMOM = MAGMOM_gen(mag_list,MagAmp,magtype,ions_num);
                INCAR_gen(systemname,INCARsystem,MAGMOM,U_struct,ISPIN,ISTART,ICHARG,LSOC,LNOCOL,NBANDS);
            case {'nosoc_ncol','soc_ncol'}
                %% PM
                systemname = 'PM-origin';
                magtype = 'PM';
                ions_num = sum(Atom_num)*3;
                mag_list = zeros(ions_num*3,1);
                MagAmp = 0;
                INCARsystem = neckname+"-"+BIG_dir+"-"+U_dir+"-"+systemname;
                MAGMOM = MAGMOM_gen(mag_list,MagAmp,magtype,ions_num);
                INCAR_gen(systemname,INCARsystem,MAGMOM,U_struct,ISPIN,ISTART,ICHARG,LSOC,LNOCOL,NBANDS);
                %% FM origin(z,x) 112(z,x) 221(z,x) gen3gen31(z,x)
                % origin z
                systemname = 'FM-origin-z';
                magtype = 'FMz';
                ions_num = sum(Atom_num);
                mag_list = 1;
                MagAmp = General_MagAmp;
                INCARsystem = neckname+"-"+BIG_dir+"-"+U_dir+"-"+systemname;
                MAGMOM = MAGMOM_gen(mag_list,MagAmp,magtype,ions_num);
                INCAR_gen(systemname,INCARsystem,MAGMOM,U_struct,ISPIN,ISTART,ICHARG,LSOC,LNOCOL,NBANDS);
                % origin x
                systemname = 'FM-origin-x';
                magtype = 'FMx';
                ions_num = sum(Atom_num);
                mag_list = 1;
                MagAmp = General_MagAmp;
                INCARsystem = neckname+"-"+BIG_dir+"-"+U_dir+"-"+systemname;
                MAGMOM = MAGMOM_gen(mag_list,MagAmp,magtype,ions_num);
                INCAR_gen(systemname,INCARsystem,MAGMOM,U_struct,ISPIN,ISTART,ICHARG,LSOC,LNOCOL,NBANDS);
                % 112 z
                systemname = 'FM-112-z';
                magtype = 'FMz';
                ions_num = sum(Atom_num)*2;
                mag_list = [1 1];
                MagAmp = General_MagAmp;
                INCARsystem = neckname+"-"+BIG_dir+"-"+U_dir+"-"+systemname;
                MAGMOM = MAGMOM_gen(mag_list,MagAmp,magtype,ions_num);
                INCAR_gen(systemname,INCARsystem,MAGMOM,U_struct,ISPIN,ISTART,ICHARG,LSOC,LNOCOL,NBANDS);
                % 112 x
                systemname = 'FM-112-x';
                magtype = 'FMx';
                ions_num = sum(Atom_num)*2;
                mag_list = [1 1];
                MagAmp = General_MagAmp;
                INCARsystem = neckname+"-"+BIG_dir+"-"+U_dir+"-"+systemname;
                MAGMOM = MAGMOM_gen(mag_list,MagAmp,magtype,ions_num);
                INCAR_gen(systemname,INCARsystem,MAGMOM,U_struct,ISPIN,ISTART,ICHARG,LSOC,LNOCOL,NBANDS);
                % gen3gen31 z
                systemname = 'FM-gen3gen31-z';
                magtype = 'FMz';
                ions_num = sum(Atom_num)*3;
                mag_list = [1 1 1];
                MagAmp = General_MagAmp;
                INCARsystem = neckname+"-"+BIG_dir+"-"+U_dir+"-"+systemname;
                MAGMOM = MAGMOM_gen(mag_list,MagAmp,magtype,ions_num);
                INCAR_gen(systemname,INCARsystem,MAGMOM,U_struct,ISPIN,ISTART,ICHARG,LSOC,LNOCOL,NBANDS);
                % gen3gen31 x
                systemname = 'FM-gen3gen31-x';
                magtype = 'FMx';
                ions_num = sum(Atom_num)*3;
                mag_list = [1 1 1];
                MagAmp = General_MagAmp;
                INCARsystem = neckname+"-"+BIG_dir+"-"+U_dir+"-"+systemname;
                MAGMOM = MAGMOM_gen(mag_list,MagAmp,magtype,ions_num);
                INCAR_gen(systemname,INCARsystem,MAGMOM,U_struct,ISPIN,ISTART,ICHARG,LSOC,LNOCOL,NBANDS);
                % 221 z
                systemname = 'FM-221-z';
                magtype = 'FMz';
                ions_num = sum(Atom_num)*4;
                mag_list = [1 1 1 1];
                MagAmp = General_MagAmp;
                INCARsystem = neckname+"-"+BIG_dir+"-"+U_dir+"-"+systemname;
                MAGMOM = MAGMOM_gen(mag_list,MagAmp,magtype,ions_num);
                INCAR_gen(systemname,INCARsystem,MAGMOM,U_struct,ISPIN,ISTART,ICHARG,LSOC,LNOCOL,NBANDS);
                % 221 x
                systemname = 'FM-221-x';
                magtype = 'FMx';
                ions_num = sum(Atom_num)*4;
                mag_list = [1 1 1 1];
                MagAmp = General_MagAmp;
                INCARsystem = neckname+"-"+BIG_dir+"-"+U_dir+"-"+systemname;
                MAGMOM = MAGMOM_gen(mag_list,MagAmp,magtype,ions_num);
                INCAR_gen(systemname,INCARsystem,MAGMOM,U_struct,ISPIN,ISTART,ICHARG,LSOC,LNOCOL,NBANDS);
                %% AFM 112 221 gen3gen31
                % 112 z
                systemname = 'AFM-112-z';
                magtype = 'AFMz';
                ions_num = sum(Atom_num)*2;
                mag_list = [1 -1];
                MagAmp = General_MagAmp;
                INCARsystem = neckname+"-"+BIG_dir+"-"+U_dir+"-"+systemname;
                MAGMOM = MAGMOM_gen(mag_list,MagAmp,magtype,ions_num);
                INCAR_gen(systemname,INCARsystem,MAGMOM,U_struct,ISPIN,ISTART,ICHARG,LSOC,LNOCOL,NBANDS);
                % 112 x
                systemname = 'AFM-112-x';
                magtype = 'AFMx';
                ions_num = sum(Atom_num)*2;
                mag_list = [1 -1];
                MagAmp = General_MagAmp;
                INCARsystem = neckname+"-"+BIG_dir+"-"+U_dir+"-"+systemname;
                MAGMOM = MAGMOM_gen(mag_list,MagAmp,magtype,ions_num);
                INCAR_gen(systemname,INCARsystem,MAGMOM,U_struct,ISPIN,ISTART,ICHARG,LSOC,LNOCOL,NBANDS);
                % gen3gen31 upupdown z
                systemname = 'AFM-gen3gen31-upupdown_z';
                magtype = 'AFMz';
                ions_num = sum(Atom_num)*3;
                mag_list = [1 1 -1];
                MagAmp = General_MagAmp;
                INCARsystem = neckname+"-"+BIG_dir+"-"+U_dir+"-"+systemname;
                MAGMOM = MAGMOM_gen(mag_list,MagAmp,magtype,ions_num);
                INCAR_gen(systemname,INCARsystem,MAGMOM,U_struct,ISPIN,ISTART,ICHARG,LSOC,LNOCOL,NBANDS);
                % gen3gen31 120
                systemname = 'AFM-gen3gen31-120';
                magtype = '120';
                ions_num = sum(Atom_num)*3;
                mag_list = [1 2 3];
                MagAmp = General_MagAmp;
                INCARsystem = neckname+"-"+BIG_dir+"-"+U_dir+"-"+systemname;
                MAGMOM = MAGMOM_gen(mag_list,MagAmp,magtype,ions_num);
                INCAR_gen(systemname,INCARsystem,MAGMOM,U_struct,ISPIN,ISTART,ICHARG,LSOC,LNOCOL,NBANDS);
                % 221-stripe1 z
                systemname = 'AFM-221-stripe1-z';
                magtype = 'AFMz';
                ions_num = sum(Atom_num)*4;
                mag_list = [1 -1 -1 1];
                MagAmp = General_MagAmp;
                INCARsystem = neckname+"-"+BIG_dir+"-"+U_dir+"-"+systemname;
                MAGMOM = MAGMOM_gen(mag_list,MagAmp,magtype,ions_num);
                INCAR_gen(systemname,INCARsystem,MAGMOM,U_struct,ISPIN,ISTART,ICHARG,LSOC,LNOCOL,NBANDS);
                % 221-stripe1 x
                systemname = 'AFM-221-stripe1-x';
                magtype = 'AFMx';
                ions_num = sum(Atom_num)*4;
                mag_list = [1 -1 -1 1];
                MagAmp = General_MagAmp;
                INCARsystem = neckname+"-"+BIG_dir+"-"+U_dir+"-"+systemname;
                MAGMOM = MAGMOM_gen(mag_list,MagAmp,magtype,ions_num);
                INCAR_gen(systemname,INCARsystem,MAGMOM,U_struct,ISPIN,ISTART,ICHARG,LSOC,LNOCOL,NBANDS);
                % 221-stripe2 z
                systemname = 'AFM-221-stripe2-z';
                magtype = 'AFMz';
                ions_num = sum(Atom_num)*4;
                mag_list = [1 1 -1 -1];
                MagAmp = General_MagAmp;
                INCARsystem = neckname+"-"+BIG_dir+"-"+U_dir+"-"+systemname;
                MAGMOM = MAGMOM_gen(mag_list,MagAmp,magtype,ions_num);
                INCAR_gen(systemname,INCARsystem,MAGMOM,U_struct,ISPIN,ISTART,ICHARG,LSOC,LNOCOL,NBANDS);
                % 221-stripe2 x
                systemname = 'AFM-221-stripe2-x';
                magtype = 'AFMx';
                ions_num = sum(Atom_num)*4;
                mag_list = [1 1 -1 -1];
                MagAmp = General_MagAmp;
                INCARsystem = neckname+"-"+BIG_dir+"-"+U_dir+"-"+systemname;
                MAGMOM = MAGMOM_gen(mag_list,MagAmp,magtype,ions_num);
                INCAR_gen(systemname,INCARsystem,MAGMOM,U_struct,ISPIN,ISTART,ICHARG,LSOC,LNOCOL,NBANDS);
        end
        cd('..');
    end
    cd(WORKSPACE+"/"+'INCAR_collection');
end
cd(WORKSPACE);
%% Gen workdir
!cp -r INCAR_collection WORKDIR
cd('WORKDIR');
for U = [0,1,3,5,7,9]
    U_dir = "U="+string(U);
    cd(U_dir);
    copyfile(WORKSPACE+"/"+'RUN_one_scf.m','RUN_one_scf.m');
    for BIG_dir_CELL  = {"nosoc_col","nosoc_ncol","soc_ncol"}
        BIG_dir = BIG_dir_CELL{1};
        cd(BIG_dir);
        %
        INCAR_system = dir();
        for i = 1:length(INCAR_system)
            system_neckname = INCAR_system(i).name;
		if isfile(system_neckname)
            system_neckname_dir = system_neckname + "_scf";
            if strcontain(system_neckname,"M")
                fprintf('The system is : \n%s\n%s\n%s\n',BIG_dir,U_dir,system_neckname);
                mkdir(system_neckname_dir);
                cd(system_neckname_dir);
                copyfile(WORKSPACE+"/POTCAR",'POTCAR');
                copyfile("../"+system_neckname,'INCAR');
                if strcontain(system_neckname,"origin")
                    copyfile(WORKSPACE+"/POSCAR_origin",'POSCAR');
                elseif strcontain(system_neckname,"112")
                    copyfile(WORKSPACE+"/POSCAR_112",'POSCAR');
                elseif strcontain(system_neckname,"221")
                    copyfile(WORKSPACE+"/POSCAR_221",'POSCAR');
                elseif strcontain(system_neckname,"gen3gen31")
                    copyfile(WORKSPACE+"/POSCAR_gen3gen31",'POSCAR');
	          else
                end
%                 ls;
                fprintf('Deploy done\n');
                cd('..')
            end
            end
        end
        %
        cd('..');
    end
    cd(WORKSPACE+"/"+'WORKDIR');
end
cd(WORKSPACE);
fprintf('All done!\n');

% eval("!"+' --tolerance 0.01 --symmetry -c '+POSCAR_init+...
%     '>'+'SymInfor_'+POSCAR_init);
