function ValueTotal = EIGENCAR_Value(EIGENCAR_DFT,EIGENCAR,extra_parm,options)
    arguments
        EIGENCAR_DFT double;
        EIGENCAR double;
        extra_parm double = [1,1];
        options.mode = 'default';
        options.algorithm  = 'pure_comparision';
    end
    if strcmp(options.mode,'extra')
        options_extra = evalin('base','options_extra');
    end
    if strcmp(options.algorithm,'pure_comparision') && ~strcmp(options.mode,'extra')
        DATA1 = EIGENCAR_DFT;
        DATA2 = EIGENCAR;
    elseif strcmp(options.algorithm,'pure_comparision') && strcmp(options.mode,'extra')
        NBAND_range = options_extra.NBAND_range;
        klist_range = options_extra.klist_range;
        DATA1 = EIGENCAR_DFT(NBAND_range,klist_range);
        DATA2 = EIGENCAR(NBAND_range,klist_range);
    end
    count = 0;
    %% Alldata
    % Direct minus
    count = count +1 ;
    Value(count) = extra_parm(count)*sqrt(mean(mean(abs(DATA1-DATA2)).^2));
    % Direct diff minus
    count = count +1 ;
    Value(count) = extra_parm(count)*sqrt(mean(mean(abs(diff(DATA1.')-diff(DATA2.')))));
    %% ...
    %%
    ValueTotal = sum(Value);
end
