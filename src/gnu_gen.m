function gnu_gen(title,E_cut,POSCAR_file,KPOINTS_file,fermi,mode)
% nargin
if nargin < 1
    title = 'BAND';
end
if nargin < 2
    E_cut = [-3,3];
end
if nargin < 3
    POSCAR_file = 'POSCAR';
end
if nargin < 4
    KPOINTS_file = 'KPOINTS';
end
if nargin < 5
    fermi = "DOSCAR";
end
if nargin < 6
    mode = 'vasp'; 
end
%
if strcmp(class(fermi),'double')
    
elseif isa(fermi,'string') ||  ischar(fermi)
    fermi = GetFermi(mode,fermi);
else
    fermi = 0;
end
[Rm,~,~,~]= POSCAR_readin(POSCAR_file,mode);
[kpoints,nodes,kpoints_name] = KPOINTS_read(KPOINTS_file,mode);
[~,~,~,kpoints_l,kpoints_name]=kpathgen3D(Rm,kpoints,nodes,kpoints_name);

fid = fopen(title+".gnu",'w');
fprintf(fid,'%s\n','set terminal  postscript enhanced color font ",20"');
fprintf(fid,'%s\n','set palette defined ( -1  "blue", 0 "grey", 1 "red" )');
fprintf(fid,'%s\n','set colorbox vertical  user origin .1,.35 size .02,.4');
fprintf(fid,'%s\n','#set term eps size 1920,1080');
fprintf(fid,'set output "%s.eps"\n',title);
fprintf(fid,'%s\n','set style data linespoints');
fprintf(fid,'%s\n','unset ztics');
fprintf(fid,'%s\n','unset key');
fprintf(fid,'%s\n','set pointsize 0.8');
fprintf(fid,'%s\n','set view 0,0');

fprintf(fid,'%s\n','set xtics font ",24"');
fprintf(fid,'%s\n','set ytics font ",24"');
fprintf(fid,'%s\n','set ylabel font ",24"');
fprintf(fid,'%s\n','set ylabel "Energy (eV)"');
fprintf(fid,'%s\n','set ylabel offset 1.5,0');
fprintf(fid,'%s%f%s%f%s\n','set xrange [',kpoints_l(1),':   ', kpoints_l(end),']');
fprintf(fid,'%s%f%s%f%s\n','set yrange [  ',E_cut(1),':   ',E_cut(2),']');

fprintf(fid,'%s','set xtics (');
for i = 1:length(kpoints_l)-1
    fprintf(fid,'"%s" %f,',kpoints_name(i),kpoints_l(i));
end

fprintf(fid,'"%s" %f',kpoints_name(end),kpoints_l(end));
fprintf(fid,'%s\n',')');

for i = 2:length(kpoints_l)-1
    fprintf(fid,'set arrow from    %f, %f to     %f, %f nohead\n',kpoints_l(i),E_cut(1),kpoints_l(i),E_cut(2));
end

fprintf(fid,'%s\n','#plot ''BAND.dat'' u 1:2  w lp lw 2 pt 7  ps 0.2');
fprintf(fid,'splot "%s.dat" u 1:2:19 w lp lw 2 pt 7 ps 0.2 palette\n',title);
fprintf(fid,'%s\n','set cbrange [-1:1]');
fprintf(fid,'%s\n','unset colorbox ');
fclose(fid);
end