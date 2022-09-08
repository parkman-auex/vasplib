%%
Mx = Oper.mirror([1,0,0]);
My = Oper.mirror([0,1,0]);
Mz = Oper.mirror([0,0,1]);
C6z = Oper.rotation(1/6,[0,0,1]);
C4z = Oper.rotation(1/4,[0,0,1]);
C3z = Oper.rotation(1/3,[0,0,1]);
C2z = Oper.rotation(1/2,[0,0,1]);
C2x = Oper.rotation(1/2,[1,0,0]);
C2x.R = round(C2x.R);
C2z.R = round(C2z.R);
C4z.R = round(C4z.R);
%%
rank = 4;
oper_real_space = Mz.R;
before = sym('X_',ones(1,rank)*3,'real');
%% anti symmetry
before(1,1,:,:) = sym(0);
before(2,2,:,:) = sym(0);
before(3,3,:,:) = sym(0);
before(3,1,:,:) = -before(1,3,:,:);
before(3,2,:,:) = -before(2,3,:,:);

before(2,1,2,1) = -before(1,2,2,1);
before(2,1,2,2) = -before(1,2,2,2);
before(1,2,1,1) = -before(2,1,1,1);
before(1,2,1,2) = -before(2,1,1,2);
%%
after = tensor_symm(before,oper_real_space);
%%
b4 = [before(2,1,1,2),before(1,2,2,2),before(2,1,1,1),before(1,2,2,1)];
a4 = [after(2,1,1,2),after(1,2,2,2),after(2,1,1,1),after(1,2,2,1)];
%% for 1/2 1/4 rotation, Mx My etc. These can be compared easily
disp(["yxxy" "xyyy" "yxxx" "xyyx"])
disp(b4);
disp(a4);
%% for 1/3 and 1/6 rotation, subs a certain num to compare
syms X_2_1_1_2 X_1_2_2_2 X_2_1_1_1 X_1_2_2_1
disp(["yxxy" "xyyy" "yxxx" "xyyx"])
b4num = subs(b4,[X_2_1_1_2, X_1_2_2_2, X_2_1_1_1, X_1_2_2_1],[1 2 -2 1]);
disp(b4num);
a4num = subs(a4,[X_2_1_1_2, X_1_2_2_2, X_2_1_1_1, X_1_2_2_1],[1 2 -2 1]);
disp(simplify(a4num));