%% PBAND_DAT_gen 
% usage =PBAND_DAT_gen(outputfile,EIGENCAR,WEIGHTCAR,klist_l)
function test=PBAND_DAT_gen(outputfile,EIGENCAR,WEIGHTCAR,klist_l)

if nargin < 2
    %band_mode!
    disp('Atom_mode');
    disp("need more inputfiles : POSCAR EIGENVAL KPOINTS(vaspkit-type) DOSCAR(for fermi) ")
end

if nargin < 3
    %band_mode!
    disp('band_mode');
else
    disp('pband_mode');
end

if nargin < 4
    %band_mode!
    disp('KPOINTS_read');
    POSCAR_read;
    [~,klist_l,~,~,~]=kpathgen3D(Rm);
end

PBANDDATAFILE=fopen(outputfile,'w');

[bandsnum,kpointsnum]=size(EIGENCAR);

for i=1:bandsnum
    %
    fprintf(PBANDDATAFILE,"#%9s %10s %10s   ------BAND: %d\n","KPOINTS","EIGEN","WEIGHT",i);
    for j=1:kpointsnum
        fprintf(PBANDDATAFILE,"%10.6f %10.6f %10.6f\n",klist_l(j),EIGENCAR(i,j),WEIGHTCAR(i,j));
    end
    %
    %fprintf(PBANDDATAFILE,"\n");
end
test=0;



end
