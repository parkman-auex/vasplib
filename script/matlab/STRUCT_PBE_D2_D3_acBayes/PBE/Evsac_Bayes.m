lattice = importdata('lattice');
lattice = lattice{1};
% matlab bash for run vasp
% https://www.mathworks.com/help/stats/bayesian-optimization-workflow.html
warning off;
tic;
% 
Func = importdata('Func');Func = Func{1};
NP = importdata('NP');
VaspRunner = importdata('vasprun');VaspRunner = VaspRunner{1};
vasprun = ['!mpirun -np ',num2str(NP),' ',VaspRunner,'>log 2>err'];
format long;
[Rm,sites,Atom_name,Atom_num,~,factor]=POSCAR_readin('POSCAR');
natom = length(sites);
[abase,cbase] = Rm2ac(Rm,lattice);
V = abs(dot(Rm(:,3),cross(Rm(:,1),Rm(:,2))))*factor^3;
fid = fopen('matlabrun.log','w');
fprintf(fid,'Volume is :%6.5f \\AA ^3 \n',V);
copyfile('POSCAR','POSCAR.bk0');
copyfile('POSCAR','POSCAR.bk');
copyfile('POSCAR','POSCAR.bk_converge');
% INCAR.RELAX INCAR.static

% mkdir
mkdir('RELAX');
mkdir('STATIC');
WORKDIR = pwd;

DecimalDepth = 3;
Vnum = 11; 
BeginingArange = [-0.1,0.1];
BeginingCrange = [-0.5,0.5];
Width = (BeginingCrange(2)-BeginingCrange(1))/2; 
DepthBase = (Vnum-1)/2;
lastVolume = 1;
lastVrange = BeginingCrange+lastVolume;
TOTVlist = zeros(DecimalDepth*Vnum,1);
TOTElist = zeros(DecimalDepth*Vnum,1);
fid2 =fopen(['EvsV_',Func,'_bayes.dat'],'w');
LossFunc  = @(para) LossFunc_vasprun_ac(para,lattice,vasprun,WORKDIR,fid,fid2);
a = optimizableVariable('a',BeginingArange + abase);
c = optimizableVariable('c',BeginingCrange + cbase);
results = bayesopt(LossFunc,[a,c],'Verbose',1,...
    'AcquisitionFunctionName','expected-improvement-plus');

% for j = 1:DecimalDepth
%     Vrangej = lastVrange;
%     TOTElistj = zeros(Vnum,1);
%     TOTclistj = zeros(Vnum,1);
%     count = 0;
%     Vrange = linspace(Vrangej(1),Vrangej(2),Vnum);
%     copyfile('POSCAR.bk_converge','POSCAR.bk');
%     Emin = [];
%     LossFunc_vasprun_ac
% 
%     [Eminj,Eminseq] = min(TOTElistj);
%     Vmin = Vrange(Eminseq);
%     lastVolume = Vmin;
%     VWidth = Width/DepthBase^j;
%     lastVrange = [Vmin-VWidth,Vmin+VWidth];
%     TOTVlist(((j-1)*Vnum+1):j*Vnum) = Vrange;
%     TOTElist(((j-1)*Vnum+1):j*Vnum) = TOTElistj;
% end
fclose(fid);
fclose(fid2);
% TOTVlist_r = (TOTVlist).^3*V/natom;
save(['EvsV_',Func,'_bayes.mat']);
toc;
