function sym_mat = sym_mat_gen(filename,startRow, endRow)


delimiter = {',','[',']','#'};
if nargin <=2
    startRow = 6;
    endRow = inf;
end



formatSpec = '%s%s%s%s%[^\n\r]';

fileID = fopen(filename,'r');


dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines', startRow(1)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines', startRow(block)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end


fclose(fileID);

raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = mat2cell(dataArray{col}, ones(length(dataArray{col}), 1));
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));

for col=[2,3,4]

    rawData = dataArray{col};
    for row=1:size(rawData, 1)

        regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
        try
            result = regexp(rawData(row), regexstr, 'names');
            numbers = result.numbers;
            

            invalidThousandsSeparator = false;
            if numbers.contains(',')
                thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
                if isempty(regexp(numbers, thousandsRegExp, 'once'))
                    numbers = NaN;
                    invalidThousandsSeparator = true;
                end
            end

            if ~invalidThousandsSeparator
                numbers = textscan(char(strrep(numbers, ',', '')), '%f');
                numericData(row, col) = numbers{1};
                raw{row, col} = numbers{1};
            end
        catch
            raw{row, col} = rawData{row};
        end
    end
end


rawNumericColumns = raw(:, [2,3,4]);
rawStringColumns = string(raw(:, 1));



R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),rawNumericColumns);
rawNumericColumns(R) = {NaN}; 


sym_operation = table;
sym_operation.Label_name = rawStringColumns(:, 1);
sym_operation.rc1 = cell2mat(rawNumericColumns(:, 1));
sym_operation.rc2 = cell2mat(rawNumericColumns(:, 2));
sym_operation.rc3 = cell2mat(rawNumericColumns(:, 3));

% sym_operation= readtable(filename, opts);






% 
cut_label = find([sym_operation.Label_name]'=='atom_mapping:');
%
rc1 = [sym_operation.rc1];
rc2 = [sym_operation.rc2];
rc3 = [sym_operation.rc3];

sym_operation_list = [rc1,rc2,rc3];
sym_operation_list(cut_label:end,:) = [];


operation_num = size(sym_operation_list,1)/5;

sym_mat = zeros(4,3,operation_num);

for i = 1:operation_num
    cut1 = (i-1)*5 +2;
    cut2 = (i-1)*5 +5;
    sym_mat(:,:,i) = sym_operation_list(cut1:cut2,:);
end


end