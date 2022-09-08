% matlab bash for run vasp
warning off;
tic;
% 
Func = importdata('Func');Func = Func{1};
NP = importdata('NP');
VaspRunner = importdata('vasprun');VaspRunner = VaspRunner{1};
vasprun = ['!mpirun -np ',num2str(NP),' ',VaspRunner,'>log 2>err'];
format long;
[Rm,sites,Atom_name,Atom_num,~,factor]=POSCAR_readin('POSCAR');
natom = length(sites);
V = abs(dot(Rm(:,3),cross(Rm(:,1),Rm(:,2))))*factor^3;
fid = fopen('matlabrun.log','w');
fprintf(fid,'Volume is :%6.5f \\AA ^3 \n',V);
copyfile('POSCAR','POSCAR.bk0');
copyfile('POSCAR','POSCAR.bk');
copyfile('POSCAR','POSCAR.bk_converge');
% INCAR.RELAX INCAR.static

% mkdir
mkdir('RELAX');
mkdir('STATIC');
WORKDIR = pwd;

DecimalDepth = 3;
Vnum = 11; 
BeginingCrange = [-0.1,0.1];
Width = (BeginingCrange(2)-BeginingCrange(1))/2; 
DepthBase = (Vnum-1)/2;
lastVolume = 1;
lastVrange = BeginingCrange+lastVolume;
TOTVlist = zeros(DecimalDepth*Vnum,1);
TOTElist = zeros(DecimalDepth*Vnum,1);

fid2 =fopen(['EvsV_',Func,'.dat'],'w');

for j = 1:DecimalDepth 
    Vrangej = lastVrange;
    TOTElistj = zeros(Vnum,1);
    TOTclistj = zeros(Vnum,1);
    count = 0;
    Vrange = linspace(Vrangej(1),Vrangej(2),Vnum);
    copyfile('POSCAR.bk_converge','POSCAR.bk');
    Emin = [];
    for v_i = Vrange
        count = count + 1;
        copyfile('POSCAR.bk','POSCAR');
        evalstring = ['!sed -i "2c ',num2str(v_i),'" POSCAR'];
        [Rm,~,~,~,~,~] = POSCAR_readin('POSCAR','vasp','digits',10);
        Rstruct = park.Rm2abc(Rm*v_i);
        eval(evalstring);
        % cd RELAX
        cd('RELAX');
        !cp ../POSCAR ../POTCAR ../KPOINTS .
        !cp ../INCAR.RELAX ./INCAR
        eval(vasprun);
        !grep 'free  energy' OUTCAR |tail -n 1|awk '{print $5}' >tmpE.dat
        E_RELAX = importdata('tmpE.dat')/natom;
        [Rm_optimize,~,~,~,~,~]= POSCAR_readin('CONTCAR','vasp','digits',10);
        Rstruct_optimize = park.Rm2abc(Rm_optimize*v_i);
        
        % cd back
        cd(WORKDIR);
        % cd STATIC
        cd('STATIC');
        !cp ../RELAX/CONTCAR POSCAR 
        !cp ../POTCAR ../KPOINTS .
        !cp ../INCAR.STATIC ./INCAR
        eval(vasprun);
        !grep 'free  energy' OUTCAR |tail -n 1|awk '{print $5}' >tmpE.dat
        E_STATIC = importdata('tmpE.dat')/natom;
        if isempty(Emin)
            Emin = E_STATIC;
            !cp POSCAR ../POSCAR.bk_converge
        else
            if E_STATIC < Emin
                Emin = E_STATIC;
                !cp POSCAR ../POSCAR.bk_converge
            end
        end
        TOTElistj(count) = E_STATIC;

        % cd back
        cd(WORKDIR);
        fprintf(fid,'--- %s ---%2d-%2d [V:%7.5f E: %13.9f]\n',datetime,j,count,v_i^3,E_STATIC);
        fprintf(fid,'LatticeC:   a       b       c      alpha     beta     gamma \n');
        fprintf(fid,'%s %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f\n','    init',Rstruct.a,Rstruct.b,Rstruct.c,Rstruct.alpha,Rstruct.alpha,Rstruct.gamma);
        fprintf(fid,'%s %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f\n','converge',Rstruct_optimize.a,Rstruct_optimize.b,Rstruct_optimize.c,Rstruct_optimize.alpha,Rstruct_optimize.alpha,Rstruct_optimize.gamma);
        % write data
        fprintf(fid2,'%7.5f %11.7f %13.9f ',v_i,v_i^3*V,E_STATIC);
        fprintf(fid2,'%7.5f %7.5f %7.5f %7.5f %7.5f %7.5f\n',Rstruct_optimize.a,Rstruct_optimize.b,Rstruct_optimize.c,Rstruct_optimize.alpha,Rstruct_optimize.alpha,Rstruct_optimize.gamma);
    end
    [Eminj,Eminseq] = min(TOTElistj);
    Vmin = Vrange(Eminseq);
    lastVolume = Vmin;
    VWidth = Width/DepthBase^j;
    lastVrange = [Vmin-VWidth,Vmin+VWidth]; 
    TOTVlist(((j-1)*Vnum+1):j*Vnum) = Vrange;
    TOTElist(((j-1)*Vnum+1):j*Vnum) = TOTElistj;
end
fclose(fid);
fclose(fid2);
 TOTVlist_r = (TOTVlist).^3*V/natom;
save(['EvsV_',Func,'.mat']);
toc;
