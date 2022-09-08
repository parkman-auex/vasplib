function H_hr = clean_hr(H_hr,mode,inifinity_small)
if nargin < 2
    mode ='H_xyz';
end
if nargin <3
    inifinity_small = 1e-6;
end

Precision = 10/inifinity_small;
%disp(Precision);
switch mode
    case 'H_xyz'
    case 'hr_concise'
    case 'hr_sparse'
    case 'H_hrz'
        Hnum_list = H_hr.HnumL ;
        vector_list = H_hr.vectorL ;
        [NRPTS,~]=size(vector_list);
        %WAN_NUM = length(Hnum_list(:,:,1));
        for i = 1:NRPTS
            %label_mat = Hnum_list(:,:,i)<inifinity_small;
            %Hnum_list(Hnum_list(:,:,i)<inifinity_small,i)=0;
            Hnum_list(:,:,i) = round(Hnum_list(:,:,i)*Precision)/Precision;
        end
        H_hr.HnumL = Hnum_list;
        H_hr.vectorL  = vector_list ;
end
end