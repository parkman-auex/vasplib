%% To generate klist for caculation and plot
%
% creat kpath from KPOINTS or otherfile 
%
%% Description of the Function:
%%
%% Usage: 
%
% * [klist_r,klist_l,klist_s,kpoints_l,kpoints_name]=kpathgen3D(Rm,kpoints,nodes,kpoints_name)
% * [klist_r,klist_l,klist_s,kpoints_l,kpoints_name]=kpathgen3D(Rm,kpoints,nodes)
% * [klist_r,klist_l,klist_s,kpoints_l,kpoints_name]=kpathgen3D(Rm,kpoints)
% * [klist_r,klist_l,klist_s,kpoints_l,kpoints_name]=kpathgen3D(Rm)
%
%% Input:
%  
% # input1: 
% # input2:
% # kpoints:
              %kpoints=[ k1x k1y k1z;...
              %          k2x k2y k2z;...
              %          k2x k2y k2z;...
              %          k3x k3y k3z;...
              %                .
              %                .
              %                .
              %          kn-1x kn-1y kn-1z;]
              %          kn-1x kn-1y kn-1z;]              
              %          knx kny knz;]
% # nodes: (nodes along each k-path)
% # kpoints_name: [K1;K2;K3]
%
%% Output:
%
% # klist_l      (for plot : the x-> will be klist_l for each band ,one-dimension)
% # kpoints_l    (for plot : To gen high symmetry points along k-path,one-dimension)
% # kpoints_name (for plot : To gen high symmetry points name along k-path,one-dimension)
% # klist_r   (for caculation: cati real)
% # klist_s   (for other function you need )       
%          
%% example:
%   commmad
% 
% 
%   result
%   
%% Note: 
% if [klist_r,klist_l,klist_s,kpoints_l,kpoints_name]=kpathgen3D() 
%          the kpoints , nodes ,kpoints_name will be given by KPOINTS POSCAR file
%
% if [klist_r,klist_l,klist_s,kpoints_l,kpoints_name]=kpathgen3D(Rm) 
%          the kpoints , nodes ,kpoints_name will be given by KPOINTS file
% note : ['Line-Mode' ],['Reciprocal' ]}
%
%% Change log
%
% * Document Date: 2020/12/03
% * Creation Date: 2020/12/03
% * Last updated : 2020/12/03
%
%% Copyright
%
% * parkman
% * <parkman@buaa.edu.cn>
%
%% Source code : 
%

function [klist_r,klist_l,klist_s,kpoints_l,kpoints_name]=kpathgen3D(Rm,kpoints,nodes,kpoints_name)   
%--------  narg  --------
    if nargin < 1
        POSCAR_read;
    end
    if nargin < 2
        [kpoints,nodes,kpoints_name] = KPOINTS_read();
    elseif nargin < 3
        [kpoints,nodes,kpoints_name] = KPOINTS_read(kpoints);
    end
    %--------  init  --------
    Gk  = (eye(length(Rm))*2*pi/Rm)' ;
    n = size(kpoints,1)/2;
    %--------  main  --------
    if length(nodes) == 1
        nodes = ones(n,1)*nodes;
    end
    %klist_symbolic
    if length(Rm) == 3
        klist_s=[];
        for i = 1:n
            klisttempX = linspace(kpoints(2*i-1,1),kpoints(2*i,1),nodes(i));
            klisttempY = linspace(kpoints(2*i-1,2),kpoints(2*i,2),nodes(i));
            klisttempZ = linspace(kpoints(2*i-1,3),kpoints(2*i,3),nodes(i));
            klist_s   = [klist_s ;klisttempX' klisttempY' klisttempZ'];
        end
        %klist_liner & kpoints_liner
        kpoints_l=[0];
        for i = 1:n
            lenghth = norm((kpoints(2*i,:)-kpoints(2*i-1,:))*Gk);
            temp = kpoints_l(i)+lenghth;
            kpoints_l = [kpoints_l;temp];
        end
        klist_l = [];
        for i = 1:n
            klisttemp = linspace(kpoints_l(i),kpoints_l(i+1),nodes(i));
            klist_l = [klist_l klisttemp];
        end
        %klist_real for caculation
        klist_r=klist_s*Gk;
        %high symmetry points store
    elseif length(Rm) == 2
        klist_s=[];
        for i = 1:n
            klisttempX = linspace(kpoints(2*i-1,1),kpoints(2*i,1),nodes(i));
            klisttempY = linspace(kpoints(2*i-1,2),kpoints(2*i,2),nodes(i));
%             klisttempZ = linspace(kpoints(2*i-1,3),kpoints(2*i,3),nodes);
            klist_s   = [klist_s ;klisttempX' klisttempY' ];
        end
        %klist_liner & kpoints_liner
        kpoints_l=[0];
        for i = 1:n
            lenghth = norm((kpoints(2*i,1:2)-kpoints(2*i-1,1:2))*Gk);
            temp = kpoints_l(i)+lenghth;
            kpoints_l = [kpoints_l;temp];
        end
        klist_l = [];
        for i = 1:n
            klisttemp = linspace(kpoints_l(i),kpoints_l(i+1),nodes(i));
            klist_l = [klist_l klisttemp];
        end
        %klist_real for caculation
        klist_r=klist_s*Gk;
        %high symmetry points store
    end
end




