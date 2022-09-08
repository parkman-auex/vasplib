%% to generate hr.dat 

% parkman 2019.9.23 

%% properties

% input :
%       H{}      (Hr.dat data format R000 R001   )
%       A[]      (Rk sequence) 
%       wann_num (dimention of matrix)
%       mp_grid    (the num of cell mp-grids)
%       \\R_abc    (3 dimension lattice matrix)
% output :
%       hr.dat

%% format of hr.dat

%1 date and time 
%2 [tab]wann_num (same with .win file)
%3 [tab]nprt     ( in .win file : mp-grid : a b c ; mpgrid=(2a+1)*(2b+1)*(2c+1) )
%4 degeneracy of each cell;15 num per line!!!!! 
%(15mod$nprt+1+5) Rk(1 2 3) M_m M_n H_r(real_part) H_i(imag_part)

%% function may be used
% 
% each H attain M_H and Traverse every elements get its real and imaginary part
% formational output

%% main_function
function main=hrdat_gen(H,A,wann_num,nprt)
%%init
%1
    hrdat=fopen('hr.dat','w');
    date_=date;
    fprintf(hrdat,"%s\n",date_);
%1
%2
    fprintf(hrdat,"         %d\n",wann_num);
%2
%3
    fprintf(hrdat,"         %d\n",nprt);
%3
%4
    linenprt=fix(nprt/15);
    remnprt=rem(nprt,15)
    while linenprt>0
        fprintf(hrdat,"    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1\n",nprt);
        linenprt=linenprt-1;
    end
    if remnprt>0
        for i=1:remnprt
            fprintf(hrdat,"    1");
        end
        fprintf(hrdat,"\n");
    end
%initlinenprt+1+3
%initlinenprt+1+3
    cellnum=length(H);
    mpgridnum=length(A);
    if cellnum ~= mpgridnum
        error("the length of H set is different from A(Rk set) ");
    end
    if wann_num ~= length(H{1})
        error("the dimision of H0  is different from wann_num ");
    end
    for i=1:cellnum
        for k=1:wann_num
            for j=1:wann_num
                fprintf(hrdat,"    %d    %d    %d",A(i,1),A(i,2),A(i,3));
                fprintf(hrdat,"    %d    %d",j,k);
                real_part=real(H{i}(j,k));
                imag_part=imag(H{i}(j,k));
                fprintf(hrdat,"    %.7f    %.7f\n",real_part,imag_part);
            end
        end
    end
    main=0;
%end
end


%% sub_function1 R_abc->Gk_123
% useless in these program
% function Gk=inverse_lattce(Ra)
% %if possible check R is 3*3 if not error
%     a=Ra(1,:);
%     b=Ra(2,:);
%     c=Ra(3,:);
%     V=abs(cross(a,b)*c');
%     Gk1= 2*pi*(cross(b,c))/V;
%     Gk2= 2*pi*(cross(c,a))/V;
%     Gk3= 2*pi*(cross(a,b))/V;
%     Gk=[Gk1;Gk2;Gk3];    
% end
% 
% %% sub_function2 Ksets_gen
% %% sub_function1 R_abc->Gk_123
% % useless in these program
% function [mp_gridsets Kset]=Ksets_gen(mp_grid,Gk)
% 
% % check mp_grid is 3 dimision vertix
% %   init
%     mp_gridsets=[];
%     Kset=[];
%     X=mp_grid(1);
%     Y=mp_grid(2);
%     Z=mp_grid(3);
%     Curve=[1/(2*X) 1/(2*X) 1/(2*X)];
% % Gk/2 W-S cell
%     Gk_1=Gk(1,:)/2;
%     Gk_2=Gk(2,:)/2;
%     Gk_3=Gk(3,:)/2;
% % reshape x y z ??????????
%     Gk_=[Gk_1(1)+Gk_1(2)+Gk_1(3) Gk_2(1)+Gk_2(2)+Gk_2(3) Gk_3(1)+Gk_3(2)+Gk_3(3)];
%   
%     
% %    nprt=X*Y*Z;
%  % gen ksets and mp_gridsets
%    for x=X:-1:-X
%         for y=Y:-1:-Y
%             for z=Z:-1:-Z
%                 tempvar=[x y z];
%                 mp_gridsets=[mp_gridsets;tempvar];
%                 Kset=[Kset;tempvar.*Gk_.*Curve];
%             end
%         end
%    end
%     
% end