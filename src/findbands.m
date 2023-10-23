


function [allk_min,anyk_max]=findbands(EIGENCAR,Erange)
[Nbands,Ktotals]=size(EIGENCAR);
fprintf("%d bands in %d Kpoints ; (Efermi set to zero)\n\n ",Nbands,Ktotals);

if nargin < 2
    Erange(1)=input('please input Emin(eV)');
    Erange(2)=input('please input Emax(eV)');
end

Emin=Erange(1);
Emax=Erange(2);

findbands=[];
for i=1:Ktotals
    Minflag=1;
    Maxflag=1;
    for j=1:Nbands
        if EIGENCAR(j,i)>Emin & Minflag==1
            Minflag = 0;
            findbands(i,1)=j;
        end
        if EIGENCAR(j,i)>Emax & Maxflag==1
            Maxflag = 0;
            findbands(i,2)=j-1;
            findbands(i,3)=j-findbands(i,1);
            continue;
        end
        
    end
    if Maxflag==1
       findbands(i,2)=j;
       findbands(i,3)=j-findbands(i,1);
    end
end
if nargin < 2
disp(findbands);
end
[allk_min.nbands,allk_min.label] = min(findbands(:,3));
allk_min.n1 = findbands(allk_min.label,1);allk_min.n2 = findbands(allk_min.label,2);
[anyk_max.nbands,anyk_max.label] = max(findbands(:,3));
anyk_max.n1 = findbands(anyk_max.label,1);anyk_max.n2 = findbands(anyk_max.label,2);
fprintf("for all k : %d %d %d \n",allk_min.n1,allk_min.n2,allk_min.nbands);
fprintf("for any k : %d %d %d \n",anyk_max.n1,anyk_max.n2,anyk_max.nbands);
end

