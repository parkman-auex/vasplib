function  [n,theta]= Rotation2nTheta(rotation,Rm)
    rotation = rotation*det(rotation);
    rotation_cart_inv = Rm.' * rotation / Rm.';
    rotation_cart = inv(rotation_cart_inv );
    R =  rotation_cart;
    theta = round(real(acosd((trace(rotation_cart)-1)/2)));
    if abs(theta- 0) < 1e-8
        n = [1 0 0];
       
    elseif abs(theta - 180) < 1e-8
        n = [sqrt(rotation_cart(1,1)/2+0.5) ,...
            ((((R(2,1)) >= 0)-0.5)/0.5)*sqrt(rotation_cart(2,2)/2+0.5) ,...
            ((((R(3,1)) >= 0)-0.5)/0.5)*sqrt(rotation_cart(3,3)/2+0.5)];
    else
        n = [R(3,2)-R(2,3),R(1,3)-R(3,1),R(2,1)-R(1,2)]/(2*sind(theta));
    end
end