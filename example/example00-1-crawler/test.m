prefix = 'https://www.sacada.info';
for i = 1:523
    thename = strcat(String_url_list{i}{2},'.tar.gz');
    if ~exist(thename,'file')
        the_download_url = strcat(prefix,strrep(String_url_list{i}{1},'"',''));
        try
            websave(thename,the_download_url);
            fprintf('%s has been download\n',thename);
        catch
            fprintf('%s cant be arrived\n',the_download_url);
        end
        
        %!tar -xzvf test.tar.gz
        pause(rand*10);
    else
        fprintf('%s has already in\n',thename);
    end
end