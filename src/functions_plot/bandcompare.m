function ax = bandcompare(EIGENCAR,EIGENCAR_DFT,Ecut,klist_l,kpoints_l,kpoints_name,options)
arguments
    EIGENCAR ;
    EIGENCAR_DFT = EIGENVAL_read();
    Ecut double= [-3,3];
    klist_l double= [];
    kpoints_l double= [];
    kpoints_name string= [];
    options.ax =  handle([]);
    options.fontname = 'Helvetica';
    options.KPOINTS = 'KPOINTS';
    options.POSCAR = 'POSCAR';
    options.Color =  @jet;
    options.title = '';
    options.klist_l = [];
    options.kpoints_l = [];
    options.kpoints_name = [];
    options.xlabel='';
    options.ylabel='E(eV)';
    options.LineSpec = '-';
    options.LineWidth = 1;
    options.MarkerSize = 3;
    options.MarkerEdgeColor = [];
    options.MarkerFaceColor = [];
end
options.LineSpec = '-';
options.Color = 'b';
propertyCell1 = namedargs2cell(options);
options.ax = bandplot(EIGENCAR_DFT,Ecut,klist_l,kpoints_l,kpoints_name,propertyCell1{:});
options.LineSpec = '.';
options.Color = 'r';
propertyCell2 = namedargs2cell(options);
ax = bandplot(EIGENCAR,Ecut,klist_l,kpoints_l,kpoints_name,propertyCell2{:});
end
