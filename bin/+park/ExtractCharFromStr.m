function [CharL,StrL] = ExtractCharFromStr(StrL,number)
arguments
    StrL string{mustBeVector};
    number = 1;
end
    CharCell = convertStringsToChars(StrL);
    if ischar(CharCell)
        CharL= CharCell(1:number);
        CharCell(1:number) = [];
        StrL = convertCharsToStrings(CharCell);
        return;
    end
    CharL = repmat("",[numel(CharCell), number]);
    for i = 1:numel(CharCell)
        CharL(i) = CharCell{i}(1:number);
        CharCell{i}(1:number) = [];
    end
    CharL = char(CharL);
    StrL = convertCharsToStrings(CharCell);
end