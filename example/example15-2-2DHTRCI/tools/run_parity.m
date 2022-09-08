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
Parity = Parity_eight(HOMO,[1,2,3,4]);

%% TRIM.dat
f = fopen('TRIM.dat','w');


fprintf(f,'# TRIM TOT PLUS MINUS\n');

fprintf(f,'  G %4d %4d %4d\n',Parity(1).tot,Parity(1).plus,Parity(1).minus);
fprintf(f,'  S %4d %4d %4d\n',Parity(2).tot,Parity(2).plus,Parity(2).minus);
fprintf(f,'  X %4d %4d %4d\n',Parity(3).tot,Parity(3).plus,Parity(3).minus);
fprintf(f,'  Y %4d %4d %4d\n',Parity(4).tot,Parity(4).plus,Parity(4).minus);

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
fclose(f2);
%% 
f3 = fopen('IsRealChern','w');
fprintf(f3,'%d',IsRealChern);
fclose(f3);



