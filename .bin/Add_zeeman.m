function H_xyz_new = Add_zeeman(H_xyz,H_zeeman_n,raw)

[NRPTS,~]=size(H_xyz);
%WAN_NUM=length(H_xyz(1).Hnum);
    H_xyz_new_ =struct('seq',[],'vector',[],'Degen',[],'key',[],'nokey',[],'Hstr',[],'Hsym',[],'Hcoe',[],'Hnum',[]);
    H_xyz_new = repmat(H_xyz_new_ ,[NRPTS,1]);    % Hamiltonian of every cell;
for n=1:NRPTS
       H_xyz_new(n).seq   = H_xyz(n).seq;
       H_xyz_new(n).Degen = H_xyz(n).Degen;
       H_xyz_new(n).Hnum  = H_xyz(n).Hnum;
       H_xyz_new(n).Hcoe  = sym(H_xyz(n).Hnum);
       H_xyz_new(n).vector= H_xyz(n).vector;
end
   
 H_xyz_new(raw).Hnum= H_xyz_new(raw).Hnum + H_zeeman_n;
 H_xyz_new(raw).Hcoe = H_xyz_new(raw).Hnum;
end