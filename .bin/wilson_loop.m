function wilson_loop(Ham,direction)
% wilson_loop(Ham, direction), direction= "kx" or "ky"

Size = size(Ham.HnumL);
BAND_index = 1:Size(1)/2;
knum1 = 200;
knum2 = 20;
gk = Ham.Gk;
if direction == "kx"
    klinex = linspace(0,gk(1,1),knum1);
    kliney = linspace(0,gk(2,2),knum2);
elseif direction == "ky"
    klinex = linspace(0,gk(1,1),knum2);
    kliney = linspace(0,gk(2,2),knum1);
end

WAVE_BZ = cell(knum1,knum2);
for ki = 1:knum1
    for kj = 1:knum2
        if direction == "kx"
            ktmp = [klinex(ki),kliney(kj),0];
        elseif direction == "ky"
            ktmp = [klinex(kj),kliney(ki),0];
        end
        if class(Ham) == "HR"
            Ham.klist_s = ktmp/gk;
            [~,WAVECAR] = Ham.EIGENCAR_gen('printmode',false,'convention','I');
        elseif class(Ham) == "Htrig"
            Ham.klist_r = ktmp;
            [~,WAVECAR] = Ham.EIGENCAR_gen('printmode',false);
        else
            [~,WAVECAR] = Ham.EIGENCAR_gen('printmode',false);
        end
        
        WAVE_BZ{ki,kj} = WAVECAR(:,BAND_index);
    end
end

theta = zeros(knum1,length(BAND_index));
for ki = 1:knum1
    Wan = WAVE_BZ{ki,1}';
    for kj = 2:knum2-1
        Wan = Wan* WAVE_BZ{ki,kj}*WAVE_BZ{ki,kj}';
    end
    Wan = Wan* WAVE_BZ{ki,1};
    [~,Ei] = eig(Wan);
    theta(ki,:)=angle(diag(Ei));
end

for i = BAND_index
    scatter(1:knum1,theta(:,1),'r.');
    hold on
end
%%
title("wilson loop for valence bands")
xlim([1,knum1]);
xticks([knum1/2,knum1]);
xticklabels(["\pi","2\pi"]);
xlabel(direction);
ylim([-pi,pi])
yticks([-pi,0,pi])
yticklabels(["-\pi",0,"\pi"]);
ylabel("\psi");