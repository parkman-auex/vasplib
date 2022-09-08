function chern = chern_number(Ham,BAND_index)
if nargin == 1
    Size = size(Ham.HnumL);
    BAND_index = 1:Size(1)/2;
end
knum = 40;
gk = Ham.Gk;
klinex = linspace(0,gk(1,1),knum);
kliney = linspace(0,gk(2,2),knum);

WAVE_BZ = cell(knum,knum);
for ki = 1:knum
    for kj = 1:knum
        if class(Ham) == "Htrig"
            Ham.klist_r = [klinex(ki),kliney(kj),0];
            [~,WAVECAR] = Ham.EIGENCAR_gen();
        elseif class(Ham) == "HR"
            Ham.klist_s = [klinex(ki),kliney(kj),0]/gk;
            [~,WAVECAR] = Ham.EIGENCAR_gen('convention','I','printmode',false);
        end
        WAVE_BZ{ki,kj} = WAVECAR(:,BAND_index);
    end
end

F = zeros(knum-1,knum-1);
for ki = 1:knum-1
    for kj = 1:knum-1
        W = WAVE_BZ{ki  ,kj  };
        W_dkx = WAVE_BZ{ki+1,kj  };
        W_dky = WAVE_BZ{ki  ,kj+1};
        W_dkxy = WAVE_BZ{ki+1,kj+1};
        
        Ux = det(W'*W_dkx);
        Ux = Ux/norm(Ux);        
        Uy = det(W'*W_dky);
        Uy = Uy/norm(Uy);        
        Uxy = det(W_dky'*W_dkxy);
        Uxy = Uxy/norm(Uxy);        
        Uyx = det(W_dkx'*W_dkxy);
        Uyx = Uyx/norm(Uyx);
        
        F(ki,kj) = log(Ux*Uyx*Uy^-1*Uxy^-1);
    end
end
chern = round(sum(F,'all')/(2i*pi),3);
%%
% figure();
% contourf(F,'linewidth',0.5);
% % colormap('jet');
% colorbar;
% xticks([]);
% yticks([]);
disp("Chern number = "+chern);


