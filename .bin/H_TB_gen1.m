%% To gen H_TB.m(use syms) with Rnn (same atom same parm)
% input: level-cut  (nn hopping max cut)
%        nn_store
%        sites
% usage: [H,H_xyz,Rm_string]=H_TB_gen2(level_cut,nn_store,Rm)
% output: 
%        H (H_up_triangle String Matrix (n,n), here n is the orbitals*atoms in one site)
%        H_TB.m ( DIY it for need)
% note : site=struct('seq',[],'inseq',[],'rc1',[],'rc2',[],'rc3',[],'name',[],'nameseq',[]);  
%        nn_t=struct('totseq',[],'orbit_on',[],'orbit_in',[],'Rlength',[],'Rc',[],'name',[],'nn_level',[]);
% warning : The POSCAR must be properly set, no extra atoms!!!  
function [H,H_xyz,orb_list]=H_TB_gen1(level_cut,nn_store,Rm,per_dir)
%% file=
disp("H_TB")
%%

if nargin <3
    syms Rm a11 a12 a13 a21 a22 a23 a31 a32 a33 real;
    Rm = [a11 a12 a13;a21 a22 a23;a31 a32 a33];
  
    disp("Model Mode");
    % mode = 'M';
end    
if nargin <4
    per_dir = [1 1 1];
end
  Rm_string=sym2str(Rm);
  disp(Rm_string);
%% term_storage

    syms t_a t_b t_c t_d t_e t_f t_g t_i t_j t_k t_l t_m t_n t_o t_p t_q t_r t_s real;
    
    syms k k_x k_y k_z real
    syms E_onsite real
    t_store=["t_a","t_b","t_c","t_d","t_e","t_f","t_g","t_f","t_i","t_j","t_k","t_l","t_m","t_n","t_o","t_p","t_q","t_r","t_s"];
    
%% init
    [N_tot,N_orbit]=size(nn_store);
    orb_list(N_orbit,:) = [0 0  0];
% H_string
    zerostring="";
    H=repmat(zerostring,[N_orbit N_orbit]);
    %H=zeros(N_orbit);% H_kj
    %H_xyz = struct('seq',[],'vector',[],'Degen',[],'key',[],'nokey',[],'Hstr',[],'Hsym',[],'Hcoe',[],'Hnum',[]);
    H_xyz(1).seq = 1;
    H_xyz(1).vector = [0 0 0];
    H_xyz(1).Degen  = 1;
    H_xyz(1).Hcoe  = sym(zeros(N_orbit));
    H_xyz(1).Hnum  = zeros(N_orbit);
%   on-site
    for j=1:N_orbit
         H(j,j)="E_onsite_"+string(nn_store(j,j).nameseq);
         EVALstring="syms E_onsite_"+string(nn_store(j,j).nameseq)+" real;";
         eval(EVALstring);
         H_xyz = set_hop(-sym(H(j,j)),j,j,[0 0 0],H_xyz,'sym');
         orb_list(j,:) = nn_store(j,j).Rc;
    end
    
    

    
%   upper triangle (only nn -term)
    for j=1:N_orbit
        for i=1:N_tot
            k=nn_store(i,j).orbit_on;
            level=nn_store(i,j).nn_level;
            if level <= level_cut
                % note sign problem
                %ion_type_num
                %disp(i);
                %disp(j);
                % count = level_cut*(nn_store(i,j).ion_type-1)+level ;
                count = level ;
                Stringrc=string(nn_store(i,j).Rr);
                %testRr = nn_store(i,j).Rr*Rm;
                %fprintf('level:%d Rr1:%s Rr2:%s Rr3:%s\n',level,string(testRr(1)),string(testRr(2)),string(testRr(3)));
                %                 disp("level=")
                %                 disp(level);
                %                 disp("k,j=")
                %                 disp(k);disp(j);
                %                 disp("term")
                %                 disp(count)
                %                 disp(Stringrc);
                % disp(Stringrc(1));
                % disp(Stringrc(2));
                % disp(Stringrc(3));
                vector_init = nn_store(i,j).Re;
                vector_init(1) = vector_init(1)*per_dir(1);
                vector_init(2) = vector_init(2)*per_dir(2);
                vector_init(3) = vector_init(3)*per_dir(3);
                vector = floor(vector_init);
                %fix bug here
                H(j,k)=H(j,k)+"-"+t_store(count)+"*exp(1i*dot([k_x,k_y,k_z],["+Stringrc(1)+" "+Stringrc(2)+" "+Stringrc(3)+"]*"+Rm_string+"))";
                H_xyz = set_hop(-sym(t_store(count)),j,k,vector,H_xyz,'sym');
            end
        end
    end

%% split
%H=str2sym(H);
%% gen
%H_file=matlabFunction(H,'File','H_TB');
end

function [out] = sym2str(sy)
%updated:  02-03-2009
%author: Marty Lawson
%
%converts symbolic variables to a matlab equation string insuring that
%only array opps are used.  
%Symbolic arrays are converted to linear Cell arrays of strings. 
%This function is especially usefull when combined with the eval() command.  
%Also, converts maple atan function to matlab atan2 and converts 
%maple "array([[a,b],[c,d]])" notation to matlab "[a,b;c,d]" notation.  
%
%Note: eval() of a matrix of functions only works if all the input 
%variables have single values.  i.e. vectors and arrays won't work.
%
%Note2: eval() does not work on Cell arrays directly.  Use "Cell_array{index}"
%inside of the eval() to keep eval() happy
%
%	EXAMPLE:
%
% X_t = dsolve('5*D2x+6*Dx+3*x = 10*sin(10*t)','x(0)=0','Dx(0)=3'); %solution is a symbolic function, X(t)
% X_t_str = sym2str(X_t);							                %convert from symbolic type to char type using array operations
% t = 0:.01:20;									                    %make "t" an array containing the time range of interest for X(t)
% X_t_vec = eval(X_t_str);							                %see "help eval" for details.
% plot(t,X_t_vec)									                %plot the results
% grid on										                    %make the plot look nice
% xlabel('time [radians]')
% ylabel('amplitude [-]')
% 
    sy = sym(sy); %insure input is symbolic
    stingin =  string(char(sy));   %find the number of elements in "sy"
 
        stingin = strrep(stingin,'^','.^');%insure that all martix opps are array opps
        stingin = strrep(stingin,'*','.*');
        stingin = strrep(stingin,'/','./');
        stingin = strrep(stingin,'atan','atan2'); %fix the atan function
        stingin = strrep(stingin,'matrix([[','['); %fix the matrix label function
        stingin = strrep(stingin,'array([[','['); %clean up any maple array notation
        stingin = strrep(stingin,'], [',';');
        stingin = strrep(stingin,']])',']');


    out = string(stingin);
end