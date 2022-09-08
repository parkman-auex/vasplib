%% RUN_NaTmO2_scf
% workdir = '/Users/parkman/Documents/MATLAB/vasplib/example/example16-1-Re_Compund166/WORKDIR/U=0';
%% init -------------------------------------------------------------------
vasprun = '/home/soft/vasp.5.4.4/bin/vasp_std';
vasprun2 = '/home/soft/vasp.5.4.4/bin/vasp_ncl';
mpirun = 'mpirun';
workdir = pwd;
target_keywords = "";
time = 0;
np = 4;
ls;
%% init -------------------------------------------------------------------
date;
import linux_matlab.*;
for BIG_dir_CELL  = {'nosoc_col','nosoc_ncol','soc_ncol'}
    BIG_dir = BIG_dir_CELL{1};
    cd(BIG_dir);
    pwd;
    ls;
    cd('..')
end
cd(workdir);
fprintf('---------------------------------------------------------\n');


%% nosoc_col first
fprintf('******************************* nosoc_col ************************************\n');
cd('nosoc_col');
INCAR_system = dir();
job_total = 0;
fprintf('.............................................................\n');
for i = 1:length(INCAR_system)
    system_neckname = INCAR_system(i).name;
    if strcontain(system_neckname,"scf")
        fprintf('system:%s run\n',system_neckname);
        cd(system_neckname);
        if strcontain(system_neckname,target_keywords)
            eval("!"+mpirun+' -np '+string(np)+' '+vasprun+'>log 2>err &');
        end
        job_total = job_total +1;
        cd('..');
	  pwd;
    end
end
fprintf('.............................................................\n');
cd(workdir);
while 1
    %pause(30*60);
    fprintf('______________________________________________________________\n');
    pause(checktime*60);%test
    time = time +checktime/60;
    fprintf("Total time cost : %f h\n",time);
    !grep '1 F=' nosoc_col/*/OSZICAR |cat -n |tee RUNING_STATE 
    [~, numstr] = system(['wc -l ', 'RUNING_STATE']);
    numstr = strsplit(numstr);
    finished_job_num=str2double(numstr{1});
    fprintf('(%d/%d) jobs has runed over ...\n',finished_job_num,job_total);
    fprintf('check no scf in 500 steps.\n')
    !grep 'DAV: 500' nosoc_col/*/OSZICAR
    fprintf('Last log for each job.\n')
    !tail -n 1 nosoc_col/*/log
    if finished_job_num == job_total
        fprintf('All done for nosoc_col! next\n');
        break;
    end
    fprintf('______________________________________________________________\n');
end
cd(workdir);
date;
%% nosoc_ncol second
fprintf('****************************** nosoc_ncol ********************************\n');
cd('nosoc_ncol');
INCAR_system = dir();
job_total = 0;
fprintf('.............................................................\n');
for i = 1:length(INCAR_system)
    system_neckname = INCAR_system(i).name;
    if strcontain(system_neckname,"scf")
        cd(system_neckname);
        fprintf('system:%s run\n',system_neckname);
        if exist('CHGCAR','file')
            movefile('CHGCAR','CHGCAR.bk');
        end
        if strcontain(system_neckname,target_keywords)
            eval("!"+vasprun2+'>log 2>err &');
        end
        job_total = job_total +1;
        cd('..');
    end
end
fprintf('.............................................................\n');
cd(workdir);
%one check point
while 1
    %pause(30*60);
    fprintf('______________________________________________________________\n');
    pause(checktime*60);%test
    time = time +checktime/60;
    fprintf("Total time cost : %f h\n",time);
    !grep '1 F=' nosoc_ncol/*/OSZICAR |cat -n |tee RUNING_STATE 
    [~, numstr] = system(['wc -l ', 'RUNING_STATE']);
    numstr = strsplit(numstr);
    finished_job_num=str2double(numstr{1});
    fprintf('(%d/%d) jobs has runed over ...\n',finished_job_num,job_total);
    fprintf('check no scf in 500 steps.\n')
    !grep 'DAV: 500' nosoc_ncol/*/OSZICAR
    fprintf('Last log for each job.\n')
    !tail -n 1 nosoc_ncol/*/log
    if finished_job_num == job_total
        fprintf('All done for nosoc_ncol! next\n');
        break;
    end
    fprintf('______________________________________________________________\n');
end
cd(workdir);
date;
fprintf('******************************** soc_ncol **********************************\n');
%% soc_ncol last
cd('soc_ncol');
INCAR_system = dir();
job_total = 0;
fprintf('.............................................................\n');
for i = 1:length(INCAR_system)
    system_neckname = INCAR_system(i).name;
    if strcontain(system_neckname,"scf")
        cd(system_neckname);
        fprintf('system:%s run\n',system_neckname);
        system_nosoc_name = system_neckname;
        % CHGCAR and WAVECAR
        if strcontain(system_neckname,target_keywords)
            copyfile(workdir+"/nosoc_ncol/"+system_nosoc_name+"/CHGCAR",'CHGCAR');
            copyfile(workdir+"/nosoc_ncol/"+system_nosoc_name+"/WAVECAR",'WAVECAR');
            eval("!"+vasprun2+'>log 2>err &');
        end
        job_total = job_total +1;
        cd('..');
    end
end
fprintf('.............................................................\n');
cd(workdir);
%one check point
while 1
    fprintf('______________________________________________________________\n');
    pause(checktime*60);
    time = time +checktime/60;
    fprintf("Total time cost : %f h\n",time);
    !grep '1 F=' soc_ncol/*/OSZICAR |cat -n |tee RUNING_STATE
    if (isunix) 
        [~, numstr] = system( ['wc -l ', 'RUNING_STATE'] );
        numstr = strsplit(numstr);
        finished_job_num=str2double(numstr{1});
    end
    fprintf('(%d/%d) jobs has runed over ...\n',finished_job_num,job_total);
    fprintf('check no scf in 500 steps.\n')
    !grep 'DAV: 500' soc_ncol/*/OSZICAR 
    fprintf('Last log for each job.\n')
    !tail -n 1 soc_ncol/*/log
    if finished_job_num == job_total
        fprintf('All done for soc_ncol! next\n');
        break;
    end
    fprintf('______________________________________________________________\n');
end
cd(workdir);
% !rm -f */*/WAVECAR
% !rm -f */*/vasprun.xml
fprintf('All done');
fprintf("Total time cost : %f h\n",time);
