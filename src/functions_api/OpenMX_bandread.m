function EIGENVAL = OpenMX_bandread(filename)
arguments
    filename string = 'openmx.Band'
end
%% get nbands and nkpts
fid = fopen(filename,'r');
%%
line1 = str2num(fgetl(fid));
nbands = line1(1);
Ef = line1(3);
line2 = fgetl(fid);
%%
nkpath = str2double(fgetl(fid));
nkpoints = 0;
for i = 1:nkpath
    linek = fgetl(fid);
    linek_list = split(linek);
    nkpoints = nkpoints + str2double(linek_list(1));
end
%%
EIGENVAL = zeros(nkpoints,nbands);
for i = 1:nkpoints
    kpoint = fgetl(fid);
    eigval_k = str2num(fgetl(fid));
    EIGENVAL(i,:) = eigval_k;
end
fclose(fid);

EIGENVAL = (EIGENVAL-Ef) * 27.21138;
EIGENVAL = EIGENVAL';
end