function  fermi= fermi_get(EIGENCAR)
[num_wan,~] = size(EIGENCAR);
[~,label]=sort(abs(EIGENCAR(num_wan/2+1,:) - EIGENCAR(num_wan/2,:)));
label_list=sort(label(1:4));
fermi=mean(mean(EIGENCAR(num_wan/2:num_wan/2+1,label_list(3:4))));
end