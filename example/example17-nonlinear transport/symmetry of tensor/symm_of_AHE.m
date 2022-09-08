%%
% reference: 10.1103/PhysRevLett.127.277202
%%
P = Oper.inversion();
Mx = Oper.mirror([1,0,0]);
My = Oper.mirror([0,1,0]);
Mz = Oper.mirror([0,0,1]);
C6x = Oper.rotation(1/6,[1,0,0]);
C6z = Oper.rotation(1/6,[0,0,1]);
C3x = Oper.rotation(1/3,[1,0,0]);
C3z = Oper.rotation(1/3,[0,0,1]);
S6x = Mx*C6x;
S6z = Mz*C6z;
%%
rank = 3;
oper_real_space = My.R;
before = sym('X_',ones(1,rank)*3,'real');
%%
before(1,1,:) = sym(0);
before(2,2,:) = sym(0);
before(3,3,:) = sym(0);
before(3,1,:) = -before(1,3,:);
before(3,2,:) = -before(2,3,:);

before(1,2,1) = -before(2,1,1); % yxx
before(2,1,2) = -before(1,2,2); % xyy
%%
after = tensor_symm(before,oper_real_space);
%%
b2 = [before(2,1,1),before(1,2,2)];
a2 = [after(2,1,1), after(1,2,2) ];
%% for 1/2 1/4 rotation, Mx My etc. These can be compared easily
disp(["yxx" "xyy"])
disp(b2);
disp(a2);
%% for 1/3 and 1/6 rotation, subs a certain num to compare
syms X_2_1_1 X_1_2_2
disp(["yxx" "xyy"])
b2num = subs(b2,[X_2_1_1 X_1_2_2],[1 2]);
disp(b2num);
a2num = subs(a2,[X_2_1_1 X_1_2_2],[1 2]);
disp(simplify(a2num));