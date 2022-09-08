function VsmMagDataStruct = VsmMagDataStruct_gen(Block_struct,VsmMagData,InputName,PpmsVsmData)
% nargin
    if nargin <3
        InputName = input('Please give your simple a name\n','s');
    end
    if nargin <4
        PpmsVsmData = [];
        mode  = 'focus';
    else
        mode  = 'full';
    end
    VsmMagDataStruct = Block_struct;
    if isa(VsmMagData,'table')
        fprintf('test past');
    end
    [nLabel_block_struct,~] = size(Block_struct);
    for i  = 1:nLabel_block_struct
        ablock = [VsmMagDataStruct(i).block(1),VsmMagDataStruct(i).block(2)];
        VsmMagDataStruct(i).datatable = ...
            VsmMagData(ablock(1):ablock(2),:);
        title = string(InputName)+'-'+num2type(VsmMagDataStruct(i).type)+'-'+VsmMagDataStruct(i).name;
        rstr = '[\/\\\:\*\?\"\<\>\|]';
        title = regexprep(title,rstr,'_');
        disp(title);
        savedata_experiment('table',VsmMagDataStruct(i,:).datatable,title);
        if strcmp(mode,'full')
            savedata_experiment('table',PpmsVsmData(ablock(1):ablock(2),:),title+'_full');
        end
    end
end 
function savedata_experiment(mode,datasource,title)
switch mode
    case 'table'
        title = title + '.dat';
        writetable(datasource,title);
end
end
function type = num2type(num)
switch num
    case 1
        type = 'ZFC';
    case 2
        type = 'FC';
    case 3
        type = 'MH';
end   
end