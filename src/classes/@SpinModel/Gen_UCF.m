function Gen_UCF(H_spin, filename, opts)
arguments
    H_spin SpinModel
    filename = 'vampire.UCF';
    opts.num_materials int8 = 2;
end
%% check Rm
H_spin.check_exchange_value();
if ~isdiag(H_spin.Rm)
    error('Vampire 6.0 only support orthogonal lattice vectors! Please check your Rm');
end
orbL = H_spin.orbL;
%% change int8 to double, very important for using mod and floor!
mag_atom_num = double(H_spin.mag_atom_num);
num_materials = double(opts.num_materials);
%% check the number of materials
if mod(mag_atom_num, num_materials) ~= 0
    error("The input number of 'opts.num_materials' can not divide the atoms evenly."+...
        "If it is what you want, please set opts.num_materials = 1 and modify the UCF file manually.")
end
nmBlock = mag_atom_num/num_materials;
%% head
fileID = fopen(filename ,'w');
%
fprintf(fileID,"# Unit cell size (Angstrom):\n");
for i=1:3
        fprintf(fileID,"  %f",H_spin.Rm(i,i));
end
%
fprintf(fileID,"\n");
fprintf(fileID,"# Unit cell vectors:\n");
fprintf(fileID,"  1.0 0.0 0.0\n");
fprintf(fileID,"  0.0 1.0 0.0\n");
fprintf(fileID,"  0.0 0.0 1.0\n");
%%
fprintf(fileID,"# Atoms num_atoms num_materials; id cx cy cz (mat cat hcat)\n");
fprintf(fileID,"%d  %d\n", mag_atom_num, num_materials);

% orbLs
for i = 1:size(orbL,1)
    fprintf(fileID,"  %d", i - 1);
    fprintf(fileID,"  %f", mod(orbL(i,1),1));
    fprintf(fileID,"  %f", mod(orbL(i,2),1));
    fprintf(fileID,"  %f", mod(orbL(i,3),1));
    fprintf(fileID,"  %d", floor((i-1)/nmBlock));
    fprintf(fileID,"\n");
end
%% body
pb = vasplib_tool_outer.CmdLineProgressBar('Writing UCF - NRPT:');
hr = H_spin.HR_obj;

J2meV = 6.2415e18 * 1e3;
hr = hr.subs( H_spin.exchange_label, H_spin.exchange_value./J2meV );

fprintf(fileID,"# Interactions; id i j dx dy dz Jij\n");
fprintf(fileID,"%d  ", nnz(hr.HnumL));
fprintf(fileID, H_spin.exchange_type + "\n");

iid = 0;
for i = 1:hr.NRPTS
    pb.print(i, hr.NRPTS, ' ...');
    for k = 1:hr.WAN_NUM
        for j = 1:hr.WAN_NUM
            if hr.HnumL(j,k,i) ~= 0            
                fprintf(fileID,"  %d", iid);
                iid = iid + 1;
                fprintf(fileID,"    %d    %d",j-1,k-1);
                fprintf(fileID,"    %d    %d    %d",...
                    hr.vectorL(i,1),hr.vectorL(i,2),hr.vectorL(i,3));         
                fprintf(fileID,"    %.3e\n", hr.HnumL(j,k,i));
            end
        end
    end
end
fclose(fileID);
end      
