function Value = loss_func_HK(parmeter,extra_parm,options)
arguments
    parmeter;
    extra_parm double =0;
    options.mode = 'default';
    options.algorithm  = 'pure_comparision';
end

%%%%%%%%%%%%%%
EIGENCAR_DFT = evalin('base','EIGENCAR_DFT');
H_hk_n = evalin('base','H_hk_n');
%%%%%%%%%%%%%%%5


if isa(parmeter,'struct')
    mode = 'Bayes';
elseif isa(parmeter,'double')
    if length(parmeter) > 1
        mode = 'NM';
    else
        mode = 'Single';
    end
end
switch mode
    case 'Single'

    case 'NM'
        %%%%%%%%%%%%%%%%%
        Varlist = evalin('base','Varlist');
        %%%%%%%%%%%%%%%%%
        H_hk_n = H_hk_n.subs(Varlist,parmeter);
        H_hk_n = H_hk_n.Subsall();
        EIGENCAR_HK = H_hk_n.EIGENCAR_gen();
        Value = EIGENCAR_Value(EIGENCAR_DFT,EIGENCAR_HK,extra_parm,'mode',options.mode,'algorithm',options.algorithm);
    case 'Bayes'
        %%%%%%%%%%%%%%%%%
        Varlist = evalin('base','Varlist');
        %%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%
        Varlock = evalin('base','Varlock');
        %%%%%%%%%%%%%%%%%
        for i  = 1:length(Varlist)
            try  % Generate Field Names from Variables  dynamic fieldnames, or sometimes dynamic field names.
                H_hk_n = H_hk_n.subs(Varlist(i),parmeter.(string(Varlist(i))));
            catch
                H_hk_n = H_hk_n.subs(Varlist(i),Varlock(i));
            end
        end
        H_hk_n = H_hk_n.Subsall();
        EIGENCAR_HK = H_hk_n.EIGENCAR_gen();
        Value = EIGENCAR_Value(EIGENCAR_DFT,EIGENCAR_HK,extra_parm,'mode',options.mode,'algorithm',options.algorithm);
end

end