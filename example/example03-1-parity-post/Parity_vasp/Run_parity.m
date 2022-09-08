%% Run parity

opts = delimitedTextImportOptions("NumVariables", 1);
opts.DataLines = [1, Inf];
opts.Delimiter = " ";
opts.VariableNames = "VarName1";
opts.VariableTypes = "double";
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts.LeadingDelimitersRule = "ignore";
HOMO = readtable("HOMO", opts);
HOMO = table2array(HOMO);
clear opts

%% Parity
Parity = Parity_eight(HOMO);

%% TRIM.dat
f = fopen('TRIM.dat','w');


fprintf(f,'# TRIM TOT PLUS MINUS\n');

fprintf(f,'  G %4d %4d %4d\n',Parity(1).tot,Parity(1).plus,Parity(1).minus);
fprintf(f,'  S %4d %4d %4d\n',Parity(2).tot,Parity(2).plus,Parity(2).minus);
fprintf(f,'  X %4d %4d %4d\n',Parity(3).tot,Parity(3).plus,Parity(3).minus);
fprintf(f,'  Y %4d %4d %4d\n',Parity(4).tot,Parity(4).plus,Parity(4).minus);
fprintf(f,'  Z %4d %4d %4d\n',Parity(5).tot,Parity(5).plus,Parity(5).minus);
fprintf(f,'  R %4d %4d %4d\n',Parity(6).tot,Parity(6).plus,Parity(6).minus);
fprintf(f,'  U %4d %4d %4d\n',Parity(7).tot,Parity(7).plus,Parity(7).minus);
fprintf(f,'  T %4d %4d %4d\n',Parity(8).tot,Parity(8).plus,Parity(8).minus);

fclose(f);
%% Real chern
f2 = fopen('RealChern.dat','w');
fprintf(f2,'# TRIM RealChernforPlanes totminus \n');

IsRealChern = 0;
%% kz = 0
ParityFour = [Parity(1).minus,...
    Parity(2).minus,...
    Parity(3).minus,...
    Parity(4).minus,...
    ];
totminus_tmp = sum(ParityFour);
RealChern_tmp = mod(totminus_tmp/2,2);
fprintf(f2,'kz0 %4d %4d\n',RealChern_tmp,totminus_tmp);
if RealChern_tmp ==1
    IsRealChern = 1;
end
%% kz = pi
ParityFour = [Parity(5).minus,...
    Parity(6).minus,...
    Parity(7).minus,...
    Parity(8).minus,...
    ];
totminus_tmp = sum(ParityFour);
RealChern_tmp = mod(totminus_tmp/2,2);
fprintf(f2,'kzP %4d %4d\n',RealChern_tmp,totminus_tmp);
if RealChern_tmp ==1
    IsRealChern = 1;
end
%% kx = 0
ParityFour = [Parity(1).minus,...
    Parity(4).minus,...
    Parity(5).minus,...
    Parity(8).minus,...
    ];
totminus_tmp = sum(ParityFour);
RealChern_tmp = mod(totminus_tmp/2,2);
fprintf(f2,'kx0 %4d %4d\n',RealChern_tmp,totminus_tmp);
if RealChern_tmp ==1
    IsRealChern = 1;
end
%% kx = pi
ParityFour = [Parity(2).minus,...
    Parity(3).minus,...
    Parity(6).minus,...
    Parity(7).minus,...
    ];
totminus_tmp = sum(ParityFour);
RealChern_tmp = mod(totminus_tmp/2,2);
fprintf(f2,'kxP %4d %4d\n',RealChern_tmp,totminus_tmp);
if RealChern_tmp ==1
    IsRealChern = 1;
end
%% ky = 0
ParityFour = [Parity(1).minus,...
    Parity(3).minus,...
    Parity(5).minus,...
    Parity(7).minus,...
    ];
totminus_tmp = sum(ParityFour);
RealChern_tmp = mod(totminus_tmp/2,2);
fprintf(f2,'ky0 %4d %4d\n',RealChern_tmp,totminus_tmp);
if RealChern_tmp ==1
    IsRealChern = 1;
end
%% ky = pi
ParityFour = [Parity(2).minus,...
    Parity(4).minus,...
    Parity(6).minus,...
    Parity(8).minus,...
    ];
totminus_tmp = sum(ParityFour);
RealChern_tmp = mod(totminus_tmp/2,2);
fprintf(f2,'kyP %4d %4d\n',RealChern_tmp,totminus_tmp);
if RealChern_tmp ==1
    IsRealChern = 1;
end
fclose(f2);
%% 
f3 = fopen('IsRealChern','w');
fprintf(f3,'%d',IsRealChern);
fclose(f3);



