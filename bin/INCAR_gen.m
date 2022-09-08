function INCAR_file = INCAR_gen(filename,systemname,MAGMOM,U_struct,ISPIN,ISTART,ICHARG,LSOC,LNOCOL,NBANDS)
    INCAR_file = fopen(filename,'w');
    % init parameter
    fprintf(INCAR_file,'SYSTEM = %s\n',systemname);
    fprintf(INCAR_file,'ISTART = %d\n',ISTART);
    if ICHARG >= 0
        fprintf(INCAR_file,'ICHARG = %d\n',ICHARG);
    else
        fprintf(INCAR_file,'#ICHARG = %d\n',ICHARG);
    end
    fprintf(INCAR_file,'ISPIN  = %d\n', ISPIN);
    fprintf(INCAR_file,'LSORBIT  = %s\n', LSOC);
    fprintf(INCAR_file,'LNONCOLLINEAR  = %s\n\n', LNOCOL);
    
    % Static parmeter
    fprintf(INCAR_file,'LREAL  = .FALSE.          (Projection operators: automatic)\n');
    fprintf(INCAR_file,'ENCUT  =  500          (Cut-off energy for plane wave basis set, in eV,set it manually)\n');
    fprintf(INCAR_file,'PREC   =  Normal       (Precision level,PREC = Low | Normal(ENOUGH) | Single | Accurate )\n');
    fprintf(INCAR_file,'LWAVE  = .TRUE.        (Write WAVECAR or not)\n');
    fprintf(INCAR_file,'LCHARG = .TRUE.        (Write CHGCAR or not)\n');
    fprintf(INCAR_file,'ADDGRID= .TRUE.        (Increase grid; helps GGA convergence)\n');
    fprintf(INCAR_file,'KSPACING = 0.5\n');    
    fprintf(INCAR_file,'NCORE  = 4\n\n');
    
    fprintf(INCAR_file,'ISMEAR =  0            (gaussian smearing method)\n');
    fprintf(INCAR_file,'SIGMA  =  0.05         (please check the width of the smearing)\n');
    fprintf(INCAR_file,'LORBIT =  11           (PAW radii for projected DOS)\n');
    fprintf(INCAR_file,'NELM   =  500           (Max electronic SCF steps) \n');  
    fprintf(INCAR_file,'EDIFF  =  1E-06        (SCF energy convergence; in eV)\n');
    fprintf(INCAR_file,'GGA_COMPAT = .FALSE.   (Apply spherical cutoff on gradient field)\n');
    fprintf(INCAR_file,'VOSKOWN    =  1        (Enhances the magnetic moments and the magnetic energies)\n\n');
    
    % DFT + U parameter
    fprintf(INCAR_file,'LDAU   = .TRUE.        (Activate DFT+U)\n');
    fprintf(INCAR_file,'LDATYPE=  2            (Dudarev; only U-J matters)\n');
    fprintf(INCAR_file,'LDAUL  = ');
    for i = 1:size(U_struct,1)
    fprintf(INCAR_file,' %d',U_struct(i,1));         
    end
    fprintf(INCAR_file,'            (Orbitals for each species)\n');
    fprintf(INCAR_file,'LDAUU  = ');
    for i = 1:size(U_struct,1)
    fprintf(INCAR_file,' %d',U_struct(i,2));         
    end
    fprintf(INCAR_file,'            (U for each species)\n');
    fprintf(INCAR_file,'LDAUJ  = ');
    for i = 1:size(U_struct,1)
    fprintf(INCAR_file,' %d',U_struct(i,3));         
    end
    fprintf(INCAR_file,'            (J for each species)\n');
    fprintf(INCAR_file,'LMAXMIX=  6            (Mixing cut-off; 4-d, 6-f)\n\n');
    % MAGMOM
    fprintf(INCAR_file,'MAGMOM = %s\n',MAGMOM);
    
    % NBANDS 
    if NBANDS > 0 
        fprintf(INCAR_file,'NBANDS = %d\n',NBANDS);
    else
        fprintf(INCAR_file,'#NBANDS = \n');
    end
    
end