%% for H_zeeman gen
% usage: H_zeeman = H_zeeman_gen(NUM_WANN,phase_list)

function H_zeeman = H_zeeman_gen(NUM_WANN,phase_list)
NUM_ORBI = NUM_WANN/2;
if NUM_ORBI == length(phase_list)
    syms m_z m_xy phi real

    Hm1=[m_z];
    Hm4=[-m_z];
    Hm2=[m_xy * exp(-1i*phi)];
    Hm3=[m_xy * exp(1i*phi)];



    Hm{1}=repmat(Hm1,1,NUM_ORBI);
    Hm{2}=repmat(Hm2,1,NUM_ORBI);
    Hm{3}=repmat(Hm3,1,NUM_ORBI);
    Hm{4}=repmat(Hm4,1,NUM_ORBI);

    %% For AFM
    for i =1:4
        for j = 1:NUM_ORBI
            Hm{i}(j) = Hm{i}(j)*phase_list(j);
        end
    end


    %% make list
    for i=1:4
        Hmag{i}=diag(Hm{i}*1);
    end


    H_zeeman=[Hmag{1},Hmag{2};Hmag{3},Hmag{4}];

    H_zeeman_latex=latex(H_zeeman);

    %diary latex_of_MODEL_BaMnPb_oneshot_AFM
    %disp(H_zeeman_latex)
else
    error("the phase list is wrong !")
end

end