function print_PDOS_information_in_dir(Pdos_namelist,width)
fprintf('We have these pDOS dats,the map is :\n');
for i = 1:length(Pdos_namelist)
    fprintf(" %d  :  %s\n",i,Pdos_namelist(i,:));
end
fprintf("And the format is :");
if width >10
    fprintf('s py pz px dxy dyz dz2 dxz dx2-y2 fy3x2  fxyz  fyz2 fz3 fxz2  fzx2  fx3 tot');
    fprintf('1  2  3  4   5   6   7   8      9    10    11    12  13   14  	 15   16  17');
else
    fprintf('s py pz px dxy dyz dz2 dxz dx2-y2 tot');
    fprintf('1  2  3  4   5   6   7   8      9  10');
end
end