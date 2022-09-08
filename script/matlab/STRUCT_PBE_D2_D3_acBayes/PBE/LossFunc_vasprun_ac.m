function E_STATIC = LossFunc_vasprun_ac(para,lattice,vasprun,WORKDIR,fid,fid2)
    a = para.a;
    c = para.c;
    Rm = ac2Rm(a,c,lattice);
    copyfile('POSCAR.bk','POSCAR');
    evalstring = ['!sed -i "3c ',num2str(Rm(1,1)),' ',num2str(Rm(1,2)),' ',num2str(Rm(1,3)),' ','" POSCAR'];
    eval(evalstring);
    evalstring = ['!sed -i "4c ',num2str(Rm(2,1)),' ',num2str(Rm(2,2)),' ',num2str(Rm(2,3)),' ','" POSCAR'];
    eval(evalstring);
    evalstring = ['!sed -i "5c ',num2str(Rm(3,1)),' ',num2str(Rm(3,2)),' ',num2str(Rm(3,3)),' ','" POSCAR'];
    eval(evalstring);
    [Rm,sites,~,~,~,v_i] = POSCAR_readin('POSCAR','vasp','digits',10);
    %     Rstruct = park.Rm2abc(Rm*v_i);
    natom = length(sites);
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
    %     if isempty(Emin)
    %         Emin = E_STATIC;
    %         !cp POSCAR ../POSCAR.bk_converge
    %     else
    %         if E_STATIC < Emin
    %             Emin = E_STATIC;
    %             !cp POSCAR ../POSCAR.bk_converge
    %         end
    %     end

    % cd back
    cd(WORKDIR);
    fprintf(fid,'--- %s --- [A:%7.5f C:%7.5f]E_RELAX: %13.9f E_STATIC: %13.9f\n',datetime,a,c,E_RELAX,E_STATIC);
    fprintf(fid,'LatticeC:   a       b       c      alpha     beta     gamma \n');
    fprintf(fid,'%s %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f\n','converge',Rstruct_optimize.a,Rstruct_optimize.b,Rstruct_optimize.c,Rstruct_optimize.alpha,Rstruct_optimize.alpha,Rstruct_optimize.gamma);
    % write data
    fprintf(fid2,'%7.5f %7.5f %13.9f ',a,c,E_STATIC);
    fprintf(fid2,'%7.5f %7.5f %7.5f %7.5f %7.5f %7.5f\n',Rstruct_optimize.a,Rstruct_optimize.b,Rstruct_optimize.c,Rstruct_optimize.alpha,Rstruct_optimize.alpha,Rstruct_optimize.gamma);
end