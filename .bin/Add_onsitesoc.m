function H_xyz_new = Add_onsitesoc(H_xyz,H_so_n,raw,mode)
if nargin <4
    mode = 'num';
end
if strcmp(mode,'num')
    [NRPTS,~]=size(H_xyz);
    %WAN_NUM=length(H_xyz(1).Hnum);
    H_xyz_new_ =struct('seq',[],'vector',[],'Degen',[],'key',[],'nokey',[],'Hstr',[],'Hsym',[],'Hcoe',[],'Hnum',[]);
    H_xyz_new = repmat(H_xyz_new_ ,[NRPTS,1]);    % Hamiltonian of every cell;
    for n=1:NRPTS
        H_xyz_new(n).seq=H_xyz(n).seq;
        H_xyz_new(n).Degen=H_xyz(n).Degen;
        H_xyz_new(n).Hnum = kron(eye(2),H_xyz(n).Hnum);
        H_xyz_new(n).Hcoe = kron(eye(2),H_xyz(n).Hcoe);
        H_xyz_new(n).vector =  H_xyz(n).vector;
    end
    
    H_xyz_new(raw).Hnum= H_xyz_new(raw).Hnum + H_so_n;
    H_xyz_new(raw).Hcoe = H_xyz_new(raw).Hcoe;
elseif strcmp(mode,'sym')
    [NRPTS,~]=size(H_xyz);
    %WAN_NUM=length(H_xyz(1).Hnum);
    H_xyz_new_ =struct('seq',[],'vector',[],'Degen',[],'key',[],'nokey',[],'Hstr',[],'Hsym',[],'Hcoe',[],'Hnum',[]);
    H_xyz_new = repmat(H_xyz_new_ ,[NRPTS,1]);    % Hamiltonian of every cell;
    for n=1:NRPTS
        H_xyz_new(n).seq=H_xyz(n).seq;
        H_xyz_new(n).Degen=H_xyz(n).Degen;
        H_xyz_new(n).Hnum = kron(eye(2),H_xyz(n).Hnum);
        H_xyz_new(n).Hcoe = kron(eye(2),H_xyz(n).Hcoe);
        H_xyz_new(n).vector =  H_xyz(n).vector;
    end
    
    H_xyz_new(raw).Hnum= H_xyz_new(raw).Hnum ;
    H_xyz_new(raw).Hcoe = H_xyz_new(raw).Hcoe + H_so_n;
end


end