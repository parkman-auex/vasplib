%% half_hr
% we always take the upper half (R >= 0) of H_hr

function H_hr = half_hr(H_hr,mode,half_dir,direction)
    if nargin <2
        mode = 'H_xyz';
    end
    if nargin <3
        half_dir = 3;
    end
    if nargin <4
        direction = 1;
    end
    %%%%%
switch mode
    case 'H_xyz'
        H_hr = simple_hr(H_hr);
    case 'hr_concise'
    case 'hr_sparse'
    case 'H_hrz'
        H_hrz = clean_hr(H_hr,'H_hrz');
        Hnum_list = H_hrz.HnumL ;
        vector_list = H_hrz.vectorL ;
        label_list = find(vector_list(:,half_dir) >= 0);
        vector_list_new = vector_list(label_list ,:);
        Hnum_list_new = Hnum_list(:,:,label_list );
        H_hr.HnumL = Hnum_list_new;
        H_hr.vectorL  = vector_list_new ;
end
end