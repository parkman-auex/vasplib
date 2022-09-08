%% 
% usage : [Energy_list,node_index_list,klist_s_list]=find_node(occupi_band_index,EIGENCAR,node_cut,klist_s)
%           [Energy_list,node_index_list,klist_s_list]=find_node(occupi_band_index,EIGENCAR,node_cut)
%               [Energy_list,node_index_list,klist_s_list]=find_node(occupi_band_index,EIGENCAR)
%               	[Energy_list,node_index_list,klist_s_list]=find_node(occupi_band_index)
function [Energy_list,node_index_list,klist_s_list]=find_node(occupi_band_index,EIGENCAR,node_cut,klist_s)
[norb,~] = size(EIGENCAR);
if nargin < 4
    POSCAR_read;
     [~,~,klist_s,~,~]=kpathgen3D(Rm);
end
if nargin <3
    node_cut  =0.0001;
end
if nargin <2
    EIGENCAR = EIGENVAl_read();
end
if nargin <1
    occupi_band_index = norb/2;
end

    
E_dif_list = EIGENCAR(occupi_band_index+1,:)-EIGENCAR(occupi_band_index,:);
node_index_list = [];
Energy_list = [];
klist_s_list = [];
for i = 1:length(E_dif_list)
    if E_dif_list(i) < node_cut
        node_index_list = [node_index_list,i];
        Energy_list = [Energy_list,(EIGENCAR(occupi_band_index+1,i)+EIGENCAR(occupi_band_index,i))/2 ]; 
        klist_s_list = [klist_s_list; klist_s(i,:)];
    end    
end
end