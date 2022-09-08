function [RCN,Plus,Minus]=RCN(PARITYCAR,WAN_NUM,TRIM_NUM)

if nargin < 3
    TRIM_NUM = 4;
end

if TRIM_NUM==4
Plus=zeros(4,1);
Minus=zeros(4,1);
for i=1:WAN_NUM/2
    if PARITYCAR(i,1,1) >0

    Plus(1) = Plus(1) +PARITYCAR(i,1,1) ;
    
    elseif PARITYCAR(i,1,1) < 0
    Minus(1) = Minus(1) + PARITYCAR(i,1,1);
    end
end

for i=1:WAN_NUM/2
    if PARITYCAR(i,1,2) >0

    Plus(2) = Plus(2) +PARITYCAR(i,1,2) ;
    
    elseif PARITYCAR(i,1,2) < 0
    Minus(2) = Minus(2) + PARITYCAR(i,1,2);
    end
end

for i=1:WAN_NUM/2
    if PARITYCAR(i,1,3) >0

    Plus(3) = Plus(3) +PARITYCAR(i,1,3) ;
    
    elseif PARITYCAR(i,1,3) < 0
    Minus(3) = Minus(3) + PARITYCAR(i,1,3);
    end
end

for i=1:WAN_NUM/2
    if PARITYCAR(i,1,4) >0

    Plus(4) = Plus(4) +PARITYCAR(i,1,4) ;
    
    elseif PARITYCAR(i,1,4) < 0
    Minus(4) = Minus(4) + PARITYCAR(i,1,4);
    end
end

RCN = mod( (Minus(1)+Minus(2)+Minus(3)+Minus(4))/2,2);

disp('Calculation of RCN finished.');

else if TRIM_NUM==2
        Plus=zeros(4,1);
Minus=zeros(4,1);
for i=1:WAN_NUM/2
    if PARITYCAR(i,1,1) >0

    Plus(1) = Plus(1) +PARITYCAR(i,1,1) ;
    
    elseif PARITYCAR(i,1,1) < 0
    Minus(1) = Minus(1) + PARITYCAR(i,1,1);
    end
end

for i=1:WAN_NUM/2
    if PARITYCAR(i,1,2) >0

    Plus(2) = Plus(2) +PARITYCAR(i,1,2) ;
    
    elseif PARITYCAR(i,1,2) < 0
    Minus(2) = Minus(2) + PARITYCAR(i,1,2);
    end
end

RCN = mod( (Minus(1)+Minus(2)*3)/2,2);
disp('Calculation of RCN finished.');
    end
end

end


