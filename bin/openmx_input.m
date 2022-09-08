function openmx_input(soc_flag)
% Input: 0 = non_soc; 1 = soc (default);
disp('此函数基于POSCAR、KPOINTS生成openmx3.9输入文件，与3.8不兼容');
disp('提醒：对于二维材料，默认的原子截断半径一般都不够大');
if ~exist('POSCAR','file')||~exist('KPOINTS','file')
    error('Missing files, POSCAR and KPOINTS needed');
end

if nargin == 0
    soc_flag = 1;
end
%%
load 'PAO_table_2019.mat' PAO_table;
[Rm,sites,Atom_name,Atom_num] = POSCAR_readin('POSCAR');
[kpoints,~,kpoints_name] = KPOINTS_read('KPOINTS');

Element = length(Atom_name);
Atom_index = zeros(Element,83,'logical');
for i = 1:Element
    softhard = ["Fe","Co","Ni","Cu","Zn"];
    if find(Atom_name(i) == softhard) ~= 0
        Atom_index(i,:) = (PAO_table.VPS==(Atom_name(i)+"_PBE19H"));      
    else
        Atom_index(i,:) = (PAO_table.VPS==(Atom_name(i)+"_PBE19"));
    end
end

%% 任务名，赝势库
fid = fopen('openmx.dat','w');
fprintf(fid, ...
    "System.Name               "+"openmx"+"\n"+...
    "DATA.PATH                 /home/soft/openmx3.9/DFT_DATA19\n\n");

%% 自洽参数，晶格弛豫，塞曼场
if soc_flag == 1
    fprintf(fid, ...
    "scf.SpinPolarization       nc         # On|Off|NC\n"+...
    "scf.SpinOrbit.Coupling     on         # On|Off, default=off\n");
else
    fprintf(fid, ...
    "scf.SpinPolarization       on         # On|Off|NC\n"+...
    "scf.SpinOrbit.Coupling     off        # On|Off, default=off\n");
end

fprintf(fid, ...
    "scf.restart                off         # on|off,default=off\n"+...
    "scf.XcType                 GGA-PBE    # LDA|LSDA-CA|LSDA-PW|GGA-PBE\n"+...
    "scf.maxIter                100        # default=40\n"+...
    "scf.EigenvalueSolver       band       # DC|GDC|Cluster|Band\n"+...
    "scf.Kgrid                  12 12 1    # means n1 x n2 x n3\n"+...
    "scf.criterion              1.0e-7     # default=1.0e-6 (Hartree)\n"+... 
    "scf.Mixing.Type            rmm-diisk  # Simple|Rmm-Diis|Gr-Pulay|Kerker|Rmm-Diisk\n"+...
    "scf.Hubbard.U              off\n"+...  
    "MD.Type                    nomd       # Nomd|Constant_Energy_MD|Opt\n"+...
    "HS.fileout                 ON         # on|off, default=Off\n\n");
fprintf(fid, ...    
    "scf.NC.Zeeman.Spin          off       # on|off, default=off\n"+...
    "scf.NC.Mag.Field.Spin       10.0      # default=0.0(Tesla)\n"+...
    "scf.NC.Zeeman.Orbital       off       # on|off, default=off\n"+...
    "scf.NC.Mag.Field.Orbital    10.0      # default=0.0(Tesla)\n\n");

%% 能带路径
lines = length(kpoints_name)-1;
fprintf(fid, ...
"Band.dispersion              on\n"+...
"  Band.Nkpath                "+lines+"\n"+...
"  <Band.kpath\n");
for i = 1:lines
    fprintf(fid,"  30 ");
    fprintf(fid,'%9.5f',kpoints(2*i-1,1),kpoints(2*i-1,2),kpoints(2*i-1,3));
    fprintf(fid,'%9.5f',kpoints(2*i,1),kpoints(2*i,2),kpoints(2*i,3));
    fprintf(fid,"  "+kpoints_name(i)+" "+kpoints_name(i+1)+"\n");
end
fprintf(fid,"  Band.kpath>\n\n");

%% 晶格常数
fprintf(fid,"Atoms.UnitVectors.Unit             Ang # Ang|AU\n"+...
            "<Atoms.UnitVectors\n"); 
for i = 1:3
    for j = 1:3
        fprintf(fid,'%14.10f',Rm(i,j));
    end
    fprintf(fid,"\n");
end
fprintf(fid,"  Atoms.UnitVectors>\n\n"); 

%% 原子基函数
fprintf(fid, ...
"Species.Number   "+Element+"\n"+...
"<Definition.of.Atomic.Species\n"); 
for i = 1:Element
    Atom_data = PAO_table(Atom_index(i,:),:);
    fprintf(fid,'%-4.2s',Atom_name(i));
    fprintf(fid,Atom_data.Standard+" "+Atom_data.VPS+"\n");
end
fprintf(fid,"Definition.of.Atomic.Species>\n\n"); 

%% 原子坐标
poscar = readcell('POSCAR');
format = char(poscar{8,1});
if format(1) == 'C'
    fprintf(fid,"Atoms.SpeciesAndCoordinates.Unit   Ang # Ang|AU\n");
elseif format(1) == 'D'
    fprintf(fid,"Atoms.SpeciesAndCoordinates.Unit   FRAC # Ang|AU\n");
end
fprintf(fid,...
"Atoms.Number   "+sum(Atom_num)+"\n"+...
"<Atoms.SpeciesAndCoordinates\n");
order = 0;
for i = 1:Element
    Atom_data = PAO_table(Atom_index(i,:),:);
    Atom_elec = Atom_data.Valenceelectrons/2;
    for j = 1:Atom_num(i)
        order = order+1;
        fprintf(fid,'%-5.f',order);
        fprintf(fid,'%-3.2s',Atom_name(i));
        fprintf(fid,'%9.5f',sites(order).rc1,sites(order).rc2,sites(order).rc3);
        
        if soc_flag == 1
            fprintf(fid," "+Atom_elec+" "+Atom_elec+" 0.0 0.0 0.0 0.0 1 off\n");
        else
            fprintf(fid," "+Atom_elec+" "+Atom_elec+" 1 off\n");
        end
    end
end
fprintf(fid,"Atoms.SpeciesAndCoordinates>\n\n");

%% DFT+U
fprintf(fid,"<Hubbard.U.values\n");
for i = 1:Element
    Atom_data = PAO_table(Atom_index(i,:),:);
    fprintf(fid,'%-4.2s',Atom_name(i));
    tot = strsplit(Atom_data.Standard,'-');
    tot = char(tot(2));
    len = length(tot)/2;
    for m = 1:len
        orbit = tot(2*m-1);
        orbit_num =str2double(tot(2*m));
        for L = 1:orbit_num
            fprintf(fid," "+num2str(L)+orbit+" 0.0");
        end
    end
    fprintf(fid,"\n");
end
fprintf(fid,"Hubbard.U.values>\n\n");

%%
fclose(fid);
