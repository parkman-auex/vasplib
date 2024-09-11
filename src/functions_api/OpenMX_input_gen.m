function OpenMX_input_gen(options)
arguments
    options.soc logical = true
    options.bandplot logical = true
    options.wannier logical = false
end
disp('此函数基于POSCAR、KPOINTS生成openmx3.9输入文件，与3.8不兼容');
disp('提醒：对于二维材料，默认的原子截断半径一般都不够大');
if ~exist('POSCAR','file')
    error('Missing POSCAR !');
end
%%
load 'PAO_table_2019.mat' PAO_table;
[Rm,sites,Atom_name,Atom_num] = POSCAR_read('POSCAR');

nElement = length(Atom_name);
Atom_index = zeros(nElement,83,'logical');
for i = 1:nElement
    softhard = ["Fe","Co","Ni","Cu","Zn"];
    if find(Atom_name(i) == softhard) ~= 0
        Atom_index(i,:) = (PAO_table.VPS==(Atom_name(i)+"_PBE19H"));      
    else
        Atom_index(i,:) = (PAO_table.VPS==(Atom_name(i)+"_PBE19"));
    end
end
%% path to the PAO lib
fid = fopen('openmx.dat','w');
fprintf(fid, ...
    "System.Name               "+"openmx"+"\n"+...
    "DATA.PATH                 /home/soft/openmx3.9/DFT_DATA19\n\n");

%% soc
if options.soc
    fprintf(fid, ...
    "scf.SpinPolarization       nc         # On|Off|NC\n"+...
    "scf.SpinOrbit.Coupling     on         # On|Off, default=off\n");
else
    fprintf(fid, ...
    "scf.SpinPolarization       off        # On|Off|NC\n"+...
    "scf.SpinOrbit.Coupling     off        # On|Off, default=off\n");
end
%% scf control parameters
fprintf(fid, ...
    "scf.EigenvalueSolver       band       # DC-LNO | Band\n"+...
    "scf.restart                off        # on|off,default=off\n"+...
    "scf.XcType                 GGA-PBE    # LDA|LSDA-CA|LSDA-PW|GGA-PBE\n"+...
    "scf.maxIter                100        # default=40\n"+...
    "scf.energycutoff           150.0      # default=150(Ry)\n"+...  
    "scf.criterion              1.0e-7     # default=1.0e-6 (Hartree)\n"+... 
    "scf.Mixing.Type            rmm-diisk  # Simple|Rmm-Diis|Gr-Pulay|Kerker|Rmm-Diisk\n"+...
    "scf.Mixing.History         30\n"+...
    "scf.Mixing.EveryPulay      1\n"+...
    "scf.Hubbard.U              off\n"+...  
    "MD.Type                    nomd       # Nomd|Constant_Energy_MD|Opt\n"+...
    "HS.fileout                 ON         # on|off, default=Off\n\n");
%% zeeman field
fprintf(fid, ...    
    "scf.NC.Zeeman.Spin         off        # on|off, default=off\n"+...
    "scf.NC.Mag.Field.Spin      10.0       # default=0.0(Tesla)\n"+...
    "scf.NC.Zeeman.Orbital      off        # on|off, default=off\n"+...
    "scf.NC.Mag.Field.Orbital   10.0       # default=0.0(Tesla)\n\n");
%% lattice parameters
fprintf(fid,"Atoms.UnitVectors.Unit             Ang # Ang|AU\n"+...
            "<Atoms.UnitVectors\n"); 
for i = 1:3
    for j = 1:3
        fprintf(fid,'%14.10f',Rm(i,j));
    end
    fprintf(fid,"\n");
end
fprintf(fid,"Atoms.UnitVectors>\n\n"); 
%% scf kpoints
a = norm(Rm(1,:));
b = norm(Rm(2,:));
c = norm(Rm(3,:));
fprintf(fid, ...
    "scf.Kgrid                  "+round(30/a)+" "+round(30/b)+" "+round(30/c)+"      # means n1 x n2 x n3\n\n");
%% PAO basis
if options.wannier
    fprintf(fid, ...
    "Species.Number   "+nElement*2+"\n"+...
    "<Definition.of.Atomic.Species\n");
else
    fprintf(fid, ...
    "Species.Number   "+nElement+"\n"+...
    "<Definition.of.Atomic.Species\n");
end
for i = 1:nElement
    Atom_data = PAO_table(Atom_index(i,:),:);
    fprintf(fid,'%-4.2s',Atom_name(i));
    fprintf(fid,Atom_data.Standard+" "+Atom_data.VPS+"\n");
    if options.wannier
        fprintf(fid,'%-5.3s',"p"+Atom_name(i));
        tmp = Atom_data.Standard.split('-');
        fprintf(fid,tmp(1)+"-s1p1d1f1 "+Atom_data.VPS+"\n");
    end       
end
fprintf(fid,"Definition.of.Atomic.Species>\n\n"); 

%% atomic sites
fprintf(fid,"Atoms.SpeciesAndCoordinates.Unit   Ang # Ang|AU\n");
fprintf(fid,...
"Atoms.Number   "+sum(Atom_num)+"\n"+...
"<Atoms.SpeciesAndCoordinates\n");
order = 0;
for i = 1:nElement
    Atom_data = PAO_table(Atom_index(i,:),:);
    Atom_elec = Atom_data.Valenceelectrons/2;
    for j = 1:Atom_num(i)
        order = order+1;
        fprintf(fid,'%-5.f',order);
        fprintf(fid,'%-3.2s',Atom_name(i));
        fprintf(fid,'%9.5f',sites(order).rc1,sites(order).rc2,sites(order).rc3);
        
        if options.soc
            fprintf(fid," "+Atom_elec+" "+Atom_elec+" 0.0 0.0 0.0 0.0 1 off\n");
        else
            fprintf(fid," "+Atom_elec+" "+Atom_elec+" 1 off\n");
        end
    end
end
fprintf(fid,"Atoms.SpeciesAndCoordinates>\n\n");

%% DFT+U
fprintf(fid,"<Hubbard.U.values\n");
for i = 1:nElement
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

%% band high-symmetry kpath
if options.bandplot
    if ~exist('KPOINTS','file')
        error('Missing KPOINTS !');
    end
    [High_Symm_Kpoints,nk_perpath,kpoints_name] = KPOINTS_read('KPOINTS');

    lines = length(kpoints_name)-1;
    fprintf(fid, ...
    "Band.dispersion             on\n"+...
    "Band.Nkpath                 "+lines+"\n"+...
    "<Band.kpath\n");
    for i = 1:lines
        fprintf(fid,"  %d ", nk_perpath);
        fprintf(fid,'%9.5f',High_Symm_Kpoints(2*i-1,1),High_Symm_Kpoints(2*i-1,2),High_Symm_Kpoints(2*i-1,3));
        fprintf(fid,'%9.5f',High_Symm_Kpoints(2*i,1),High_Symm_Kpoints(2*i,2),High_Symm_Kpoints(2*i,3));
        fprintf(fid,"  "+kpoints_name(i)+" "+kpoints_name(i+1)+"\n");
    end
    fprintf(fid,"Band.kpath>\n\n");
end

%% wannier 
if options.wannier
    fprintf(fid, ...
    "Wannier.Func.Calc        on           \n"+...
    "Wannier90.fileout        on           \n"+...
    "Wannier.MaxShells        12           # default value is 12.\n"+... 
    "Wannier.Kgrid           "+round(30/a)+" "+round(30/b)+" "+round(30/c)+"\n\n");

    fprintf(fid, ...
    "Wannier.Initial.Projectors.Unit     FRAC\n"+...
    "Wannier.Func.Num              0     #no default\n");

    fprintf(fid,...
    "<Wannier.Initial.Projectors\n");
    order = 0;
    for i = 1:nElement
        for j = 1:Atom_num(i)
            order = order + 1;
            fprintf(fid, "p"+Atom_name(i)+"-proj"+i);
            fprintf(fid,'%9.5f',sites(order).rc1,sites(order).rc2,sites(order).rc3);
            
            fprintf(fid,"  0.0 0.0 1.0    1.0  0.0  0.0\n");
        end
    end
    fprintf(fid,"Wannier.Initial.Projectors>\n\n");

    fprintf(fid, ...
    "Wannier.Outer.Window.Bottom   -10.0   #eV, relative to the Fermi level\n"+...
    "Wannier.Outer.Window.Top       10.0   #eV, relative to the Fermi level\n"+...
    "Wannier.Inner.Window.Bottom    -1.0   #eV, relative to the Fermi level\n"+...
    "Wannier.Inner.Window.Top        1.0   #eV, relative to the Fermi level\n");
end
%%
fclose(fid);
