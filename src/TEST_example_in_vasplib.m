function TEST_example_in_vasplib(dir,runfile,basedir)
    workshop=pwd;
    mkdir('test');
    cd test;
    if nargin < 3
        basedir = '../example/';
    end
    copyfile([basedir,dir]);
    if nargin <2
        % check runfile
        % begin with test
    end
    try
    run(runfile);
    catch
        fprintf('contact me to fix this issue: parkman@buaa.edu.cn');
    end
    cd(workshop);
    rmdir test/ s

end