%% Green_prepare
function    [H00,H01,H_hrz_green] = Green_prepare(H_hrz,principle_layer,fin_dir)

if nargin <3
    fin_dir = 3;
end
if nargin <2
    principle_layer = 3;
end
% init
vector_list = H_hrz.vectorL ;
HnumL = H_hrz.HnumL ;
WAN_NUM = length(HnumL(:,:,1));
label_list = vector_list(:,fin_dir);
%lin_000 = find(vector_list(:,fin_dir) == 0);
Poly_priciplayer_mat = Poly_priciplayer_mat_gen(principle_layer);
%main
H00 = zeros(WAN_NUM*principle_layer);
for i = -(principle_layer-1):(principle_layer-1)
    z = i;
    label = find(label_list == z); 
    if label >0
        H00=H00+kron(Poly_priciplayer_mat(:,:,i+principle_layer),HnumL(:,:,label));
    end
end

H10 = zeros(WAN_NUM*principle_layer);
for i = 1:principle_layer
    z = i;
    label = find(label_list == z); 
    if label >0
        H10=H10+kron(Poly_priciplayer_mat(:,:,i),HnumL(:,:,label));
    end
end
H01 = H10' ;
H_hrz_green.HnumL(:,:,1) = H10;H_hrz_green.vectorL(1,:) = [0 ,0, -1];
H_hrz_green.HnumL(:,:,2) = H00;H_hrz_green.vectorL(2,:) = [0 ,0, 0];
H_hrz_green.HnumL(:,:,3) = H01;H_hrz_green.vectorL(3,:) = [0 ,0, 1];

end


function Poly_priciplayer_mat = Poly_priciplayer_mat_gen(principle_layer)
    Poly_priciplayer_mat(:,:,principle_layer) = eye(principle_layer);
    
    for i = 1:principle_layer-1
        base_mat = zeros(principle_layer);
        base_mat(1:principle_layer-i,i+1:principle_layer) = eye(principle_layer-i);
        Poly_priciplayer_mat(:,:,principle_layer+i) = base_mat;
        Poly_priciplayer_mat(:,:,principle_layer-i) = base_mat';
    end
end