%% rm_list_gen
% for this funciton work, you should write a judge function return 0 or 1 in the current dir
% the default judge function name is rm_or_not
% the function's input parm is :rc_onepoint output parm is 1 or 0
% usage: [orb_rm_index_list,orb_rm_list] = Gen_rmlist(orb_list,rm_or_not_functionname)
function [orb_rm_index_list,orb_rm_list] = Gen_rmlist(orb_list,rm_or_not_functionname)
%% nargin 
if nargin < 2
    rm_or_not_functionname = 'rm_or_not';
    % we can add a draft stencil here for this function, waiting ...
end

%% test if exist 
rm_or_not_functionname_m = 'rm_or_not'+".m";
if ~exist(rm_or_not_functionname_m,'file')
    error('file is not exist! please give the right rm_or_not_functionname.');
end
%% test if work
disp('test the function work or not? default input rc = [0 0 0]');
rc_test = [0,0,0];
test =  rm_or_not_eval(rm_or_not_functionname,rc_test);
if test == 0 | test == 1
    disp('test past!');
else
    disp('test error');
    error('The rm_or_not_funciton disfunction! try another one');
end

%%  Gen orb_rm_index_list,orb_rm_list
orb_rm_index_list =[];
orb_rm_list = [];
[norb,~] = size(orb_list);
for i = 1:norb
    temp_orb  = orb_list(i,:);
    logic_result = rm_or_not_eval(rm_or_not_functionname,temp_orb);
%     disp(temp_orb);
%     disp(logic_result);
    if logic_result == 1
        orb_rm_index_list =[orb_rm_index_list;i];
        orb_rm_list = [orb_rm_list ;temp_orb ];
    end
end
%% 
return

end

function logic_result = rm_or_not_eval(rm_or_not_functionname,rc)
    eval_string  = rm_or_not_functionname+"(["+num2str(rc(1))+","+...
        num2str(rc(2))+","+...
        num2str(rc(3))+"]);";
    logic_result = eval(eval_string);
    
end