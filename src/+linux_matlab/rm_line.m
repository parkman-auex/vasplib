% function rm_line
function [outstring] = rm_line(filename,rm_number_list,outfilename)
import linux_matlab.*
if nargin <3
    mode = 'silence';
    %outfilename = 'none';
else
    mode = 'aloud';
end
     outstring = read_file(filename);
     % rm
     outstring(rm_number_list,:) = [];
     
    if ~strcmp(mode,'silence')
        write_file(outstring,outfilename);
    end
end