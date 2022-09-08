%% What you have is a fully syms Hamiton
%  Change it into a TB
%  is cubic important ? here try
%  Be careful the sruct H_xyz
%% H_xyz_ =struct('seq',[],'vector',[],'key',[],'nokey',[],'Hstr',[],'Hsym',[],'Hcoe',[],'Hnum',[]);

%% example
% model ="real_para";
% syms_para_load(model);
% KaneFourband;
% discretize_kp_model;
% Then run this code
% H=subs(H);
% disp(H);
% expand(H);
%_______________________________
%% expand for 
H_check =  expand(H);
disp(H_check);
%char(H_check)
H_string=string(H_check);
[Ni,Nj]=size(H_string);
H_string=strrep(H_string,' - ',' +- ');% for division

%% Note: for bug fix ;It is not an elegant way indeed

      %------------------fix bug--------------------
NOKEY=["k_1","k_2","k_3","k_12"];
      %------------------fix bug--------------------

%division
for i =1:Ni
    for j=1:Nj
        %------------------fix bug--------------------
        tempflag=strcontain(H_string(i,j),NOKEY);
        if  tempflag== 0
                    %disp(strobj);
                    H_string(i,j)=H_string(i,j)+" ";
                    %disp(tempstring);
        end
        %------------------fix bug--------------------
        H_str{i,j}=strsplit(H_string(i,j),'+');
        
    end
end


% typical show the struct
H_xyz_ =struct('seq',[],'vector',[],'Hstr',[],'Hcoe',[]);
NRPTS = 27;
% you can innounce the H_xyz{} here use repmap [I am not sure the name of that]
H_xyz = repmat(H_xyz_,[NRPTS 1]);





%% You should mix here into a function 

%% zero
% The first is Constant 
%  H_xyz(1).seq=1  ; H_xyz(1).vector=[ 0  0  0] ; H_xyz(1).key=" ";              H_xyz(1).nokey=["k_x","k_y","k_z"];
% remind key is a string ;while nokey is a string list 
% The others is flux 
%For more combination use vector and basic string we can gen the key and nokey automatically
%% single
%  H_xyz(2).seq= 2 ; H_xyz(2).vector=[ 1  0  0] ; H_xyz(2).key=["k_x"] ;H_xyz(2).nokey=["k_y","k_z"];
%  H_xyz(3).seq= 3 ; H_xyz(3).vector=[ 0  1  0] ; H_xyz(3).key="exp(b*k_y*1i)" ;H_xyz(3).nokey=["k_z","k_x"];
%  H_xyz(4).seq= 4 ; H_xyz(4).vector=[ 0  0  1] ; H_xyz(4).key="exp(c*k_z*1i)" ;H_xyz(4).nokey=["k_x","k_y"];
% 
%  H_xyz(5).seq= 5 ; H_xyz(5).vector=[-1  0  0] ; H_xyz(5).key="exp(-a*k_x*1i)";H_xyz(5).nokey=["k_y","k_z"];
%  H_xyz(6).seq= 6 ; H_xyz(6).vector=[ 0 -1  0] ; H_xyz(6).key="exp(-b*k_y*1i)";H_xyz(6).nokey=["k_z","k_x"];
%  H_xyz(7).seq= 7 ; H_xyz(7).vector=[ 0  0 -1] ; H_xyz(7).key="exp(-c*k_z*1i)";H_xyz(7).nokey=["k_x","k_y"];
% %% Qual
%  H_xyz(8).seq= 8 ; H_xyz(8).vector=[ 1  1  0] ; H_xyz(8).key="exp(a*k_x*1i)*exp(b*k_y*1i)"  ;H_xyz(8).nokey=["k_z"];
%  H_xyz(9).seq= 9 ; H_xyz(9).vector=[ 1 -1  0] ; H_xyz(9).key="exp(a*k_x*1i)*exp(-b*k_y*1i)" ;H_xyz(9).nokey=["k_z"];
% H_xyz(10).seq=10 ;H_xyz(10).vector=[-1  1  0] ;H_xyz(10).key="exp(-a*k_x*1i)*exp(b*k_y*1i)" ;H_xyz(10).nokey=["k_z"];
% H_xyz(11).seq=11 ;H_xyz(11).vector=[-1 -1  0] ;H_xyz(11).key="exp(-a*k_x*1i)*exp(-b*k_y*1i)";H_xyz(11).nokey=["k_z"];
% 
% H_xyz(12).seq=12 ;H_xyz(12).vector=[ 0  1  1] ;H_xyz(12).key="exp(b*k_y*1i)*exp(c*k_z*1i)"  ;H_xyz(12).nokey=["k_x"];
% H_xyz(13).seq=13 ;H_xyz(13).vector=[ 0  1 -1] ;H_xyz(13).key="exp(b*k_y*1i)*exp(-c*k_z*1i)" ;H_xyz(13).nokey=["k_x"];
% H_xyz(14).seq=14 ;H_xyz(14).vector=[ 0 -1  1] ;H_xyz(14).key="exp(-b*k_y*1i)*exp(c*k_z*1i)" ;H_xyz(14).nokey=["k_x"];
% H_xyz(15).seq=15 ;H_xyz(15).vector=[ 0 -1 -1] ;H_xyz(15).key="exp(-b*k_y*1i)*exp(-c*k_z*1i)";H_xyz(15).nokey=["k_x"];
% 
% H_xyz(16).seq=16 ;H_xyz(16).vector=[ 1  0  1] ;H_xyz(16).key="exp(a*k_x*1i)*exp(c*k_z*1i)"  ;H_xyz(16).nokey=["k_y"];
% H_xyz(17).seq=17 ;H_xyz(17).vector=[-1  0  1] ;H_xyz(17).key="exp(-a*k_x*1i)*exp(c*k_z*1i)" ;H_xyz(17).nokey=["k_y"];
% H_xyz(18).seq=18 ;H_xyz(18).vector=[ 1  0 -1] ;H_xyz(18).key="exp(a*k_x*1i)*exp(-c*k_z*1i)" ;H_xyz(18).nokey=["k_y"];
% H_xyz(19).seq=19 ;H_xyz(19).vector=[-1  0 -1] ;H_xyz(19).key="exp(-a*k_x*1i)*exp(-c*k_z*1i)";H_xyz(19).nokey=["k_y"];
% %% Tri
%  H_xyz(20).seq=20 ;H_xyz(20).vector=[ 1  1  1] ;H_xyz(1).key="k";
%  H_xyz(21).seq=21 ;H_xyz(21).vector=[ 1  1 -1] ;H_xyz(1).key="k";
%  H_xyz(22).seq=22 ;H_xyz(22).vector=[ 1 -1  1] ;H_xyz(1).key="k";
%  H_xyz(23).seq=23 ;H_xyz(23).vector=[ 1 -1 -1] ;H_xyz(1).key="k";
% 
%  H_xyz(24).seq=24 ;H_xyz(24).vector=[ 1  1  1] ;H_xyz(1).key="k";
%  H_xyz(25).seq=25 ;H_xyz(25).vector=[ 1 -1  1] ;H_xyz(1).key="k";
%  H_xyz(26).seq=26 ;H_xyz(26).vector=[-1  1  1] ;H_xyz(1).key="k";
%  H_xyz(27).seq=27 ;H_xyz(27).vector=[-1 -1  1] ;H_xyz(1).key="k";

%% Other



[N_H_xyz,~]=size(H_xyz);
% line_000=14;
%% learn REG
for n=1:N_H_xyz
    % init
    H_temp = repmat("",[Ni Nj]);
    H_coe =  repmat(sym(0),[Ni Nj]);
    H_xyz(n).Hcoe=H_coe;   
    H_xyz(n).Hstr=H_temp;
end

for n1=-1:1:1
    for n2 =-1:1:1
        for n3 =-1:1:1
            vector=[n1,n2,n3];
            v1=vector(1)+2;v2=vector(2)+2;v3=vector(3)+2;
            seq=(v1-1)*9+(v2-1)*3+v3;
            H_xyz(seq).seq=seq;
            H_xyz(seq).vector=vector;
        end
    end
end

for i=1:Ni       
        for j=1:Nj
            for k=1:length(H_str{i,j})
                % you need a tool function here
                strobj=H_str{i,j}(k);
                [seq,vector,Hcoe_single,strout]=split_string2tb(strobj);
                H_xyz(seq).Hcoe(i,j)=H_xyz(seq).Hcoe(i,j)+Hcoe_single;
                H_xyz(seq).Hstr(i,j)=H_xyz(seq).Hstr(i,j)+strout;
                %H_xyz(seq).vector=vector;
            end
            
        end
        
end


%% disp whole
%Disp_H_xyz(H_xyz,'whole');
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        