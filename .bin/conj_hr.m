%% 
function H_hr = conj_hr(H_hr)
    [NRPTS,~]=size(H_hr);
    for i = 1:NRPTS
        try
            H_hr(i).Hcoe = conj(H_hr(i).Hcoe);
        catch
            
        end
        try
            H_hr(i).Hnum = conj(H_xyz(i).Hnum);
        catch
            
        end
    end
end