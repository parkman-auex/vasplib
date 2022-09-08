clear;
json_file = 'db.json'; % line split !!
fid=fopen(json_file ,'r');
count = 1;
Database = jsondecode(fgetl(fid));
if  Database.bandgap > 0
    countadd = 1;
else
    countadd = 0;
end
%Database.literature_doi = 'none';
Fnames = fieldnames(Database );
NFnames = length(Fnames);

while ~feof(fid)
    J=fgetl(fid);
    STRUCT = jsondecode(J); % change json 2 struct
    count = count+1;
    if STRUCT.bandgap > 0
        for i =1:NFnames
            if ~isfield(STRUCT,Fnames{i})
                %STRUCT.Fnames{i} = nan;
                STRUCT = setfield(STRUCT,Fnames{i},nan);
            end
        end
        if isfield(STRUCT,'literature_doi')
            STRUCT=rmfield(STRUCT,'literature_doi');
        end
        countadd = countadd + 1;
        fprintf('%d - %7.4f - in %d th DATABASE \n',count,STRUCT.bandgap,countadd);
        Database(countadd,1) = STRUCT;
    end
end
fclose(fid);
Database_2Dmatpedia = Database;
save('Database_2Dmatpedia','Database_2Dmatpedia');