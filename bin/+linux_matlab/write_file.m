%% write_file
function write_file(string_list,filename)
[nline,~]=size(string_list);
fid = fopen(filename, 'w');
for i = 1:nline
    fprintf(fid,'%s\n',string_list(i,:));
end
fclose(fid);
end