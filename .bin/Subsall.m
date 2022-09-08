%% subsall
% !! the add Degen is temperally set here
function out_H_xyz=Subsall(H_xyz,mode)
if nargin<2
    mode = 'struct';
end
if strcmp(mode,'cell')
    [~,N_H_xyz]=size(H_xyz);
    out_H_xyz=H_xyz;
    for i=1:N_H_xyz
%         out_H_xyz{i}.Hsym=subs(H_xyz{i}.Hsym);
        out_H_xyz{i}.Hcoe=subs(H_xyz{i}.Hcoe);
    end
elseif strcmp(mode,'struct')
    [N_H_xyz,~]=size(H_xyz);
    out_H_xyz=H_xyz;
    for i=1:N_H_xyz
%         out_H_xyz(i).Hsym=subs(H_xyz(i).Hsym);
%  disp(i);
        out_H_xyz(i).Hcoe=subs(H_xyz(i).Hcoe);
        
        out_H_xyz(i).Hnum=double(out_H_xyz(i).Hcoe);
        %
        out_H_xyz(i).Degen=1;
    end
elseif strcmp(mode,'file')
    [N_H_xyz,~]=size(H_xyz);
    out_H_xyz=H_xyz;
    load('parm.mat');
    for i=1:N_H_xyz
        %         out_H_xyz(i).Hsym=subs(H_xyz(i).Hsym);
        %  disp(i);
        out_H_xyz(i).Hcoe=subs(H_xyz(i).Hcoe);
        
        out_H_xyz(i).Hnum=double(out_H_xyz(i).Hcoe);
        %
        out_H_xyz(i).Degen=1;
    end
elseif strcmp(mode,'sym')
    [N_H_xyz,~]=size(H_xyz);
    out_H_xyz=H_xyz;
    for i=1:N_H_xyz
        %         out_H_xyz(i).Hsym=subs(H_xyz(i).Hsym);
        %  disp(i);
        out_H_xyz(i).Hcoe=subs(H_xyz(i).Hcoe);
        
        %out_H_xyz(i).Hnum=double(out_H_xyz(i).Hcoe);
        %
        out_H_xyz(i).Degen=1;
    end
else
    %% need more
    [N_H_xyz,~]=size(H_xyz);
    out_H_xyz=H_xyz;
    syms k_x k_y k_z real;
    for i=1:N_H_xyz
%         out_H_xyz(i).Hsym=subs(H_xyz(i).Hsym);
        out_H_xyz(i).Hcoe=subs(H_xyz(i).Hcoe);
        out_H_xyz(i).Hfun=matlabFunction(out_H_xyz(i).Hsym,'Vars',{k_x,k_y,k_z});
        out_H_xyz(i).Hnum=double(out_H_xyz(i).Hcoe);
        %
        out_H_xyz(i).Degen=1;
    end
end

end
    