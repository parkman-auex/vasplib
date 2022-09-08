function R_struct = Rm2abc(Rm)
    R_struct.a = norm(Rm(1,:));
    R_struct.b = norm(Rm(2,:));
    R_struct.c = norm(Rm(3,:));
    R_struct.alpha = acos(dot(Rm(2,:),Rm(3,:))/(norm(R_struct.b)*norm(R_struct.c)))/pi*180;
    R_struct.beta  = acos(dot(Rm(3,:),Rm(1,:))/(norm(R_struct.c)*norm(R_struct.a)))/pi*180;
    R_struct.gamma = acos(dot(Rm(1,:),Rm(2,:))/(norm(R_struct.a)*norm(R_struct.b)))/pi*180;
end 