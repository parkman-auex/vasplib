function DOSCAR = DOSCAR_gen(GREENCAR,mode)
% nargin 
if nargin <2
    mode = 'green';
end

if strcmp(mode,'green')
    if iscell(GREENCAR)
        [Nsize1,Nsize2] = size(GREENCAR);
        DOSCAR =zeros(Nsize1,Nsize2);
        for i = 1:Nsize1
            for j = 1:Nsize2
                DOSCAR(i,j) = -trace(imag(GREENCAR{i,j}));
            end
        end
    else
        disp('temp support cell format for GREENCAR');
    end
end
end

