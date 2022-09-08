%% Read wannier90_hr.dat
% Be careful of the struct
% usage [H_xyz,row] =Readhr(hrdat)
% usage [H_xyz,row] =Readhr(hrdat,mode) mode = 'n'
% H_xyz_ =struct('seq',[],'vector',[],'Degen',[],'key',[],'nokey',[],'Hstr',[],'Hsym',[],'Hcoe',[],'Hnum',[]);

function [H_xyz,row] = Readhr(hrdat,mode)
    [H_xyz,row] = Read_hr(hrdat,mode);
end









