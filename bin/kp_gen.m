%KP_GEN kp_init (delete soon)
% init kp model 
%
%
% usage: 
% input:  LO:LinearOperation_collection
%         ALO:AntilinearOperation_collection
%         k_order
% outplut: 
%         H_kp
%         A file with name MODEL_kp
%         A file with name PARA_kp
%         A file with name RUN_kp
% 

%% Main function
function H_kp  = kp_gen(LO_set,ALO_set,k_order,mat_stencil,mode)

%% nargin 
if nargin < 3
    k_order = 2;
end

if nargin < 4
    mat_stencil = eye(4);
end

if nargin <5
    disp('full kp Hamitonian ');
    mode = 'fullkp';
end


%% check error
for i = 1:1
    if k_order >3 | k_order<0
        error('we only support k_order equals 0, 1, 2, 3');
    end
    [wLO,nLO] =  size(LO_set);
    [wALO,nALO] =  size(ALO_set);
    if ~(wLO == 0 | wLO == 2 | wALO == 0 | wALO == 2 )
        error('LO_set,ALO_set error, check completence');
    end
    if nLO == 0 & nALO == 0 
        error('LO_set,ALO_set error, empty');
    end  
    
    if nLO>0
        D_kp = length(LO_set{1,1});
        d_k  = length(LO_set{2,1});
        % check each 
        for i = 1:nLO
            if D_kp ~= length(LO_set{1,i})
               fprintf('check %d LO_set' ,i);
               error('LO_set,ALO_set error, dimision not equal');
            end
            if d_k ~= length(LO_set{2,i})
               fprintf('check %d LO_set(k)' ,i);
               error('LO_set,ALO_set error, dimision not equal');
            end
        end
    end     
    
    if ~isempty(ALO_set)
        if ~isempty(LO_set)
            if length(ALO_set{1,1}) ~= length(LO_set{1,1})
                error('LO_set,ALO_set error, dimision not equal');
            end
            if length(ALO_set{2,1}) ~= length(LO_set{2,1})
                error('LO_set,ALO_set error, dimision not equal');
            end
        end
        
        D_kp = length(ALO_set{1,1});
        d_k  = length(ALO_set{2,1});
        % check each 
        for i = 1:nALO
            if D_kp ~= length(ALO_set{1,i})
               fprintf('check %d ALO_set' ,i);
               error('LO_set,ALO_set error, dimision not equal');
            end
            if d_k ~= length(ALO_set{2,i})
               fprintf('check %d ALO_set(k)' ,i);
               error('LO_set,ALO_set error, dimision not equal');
            end
        end
    end
end   
    
    
    %% first init 
    syms K_x K_y K_z real;
    syms K_a K_b K_c real;
    if strcmp(mode,'fullkp')
        % for real_opereation
        Basis_set_Rm= k_sym_list_gen(k_order,[K_x,K_y,K_z]);
        % for k opereation
        %     Basis_set_Gk= Basis_set_Rm2Gk(Basis_set_Rm,[K_x,K_y,K_z],[K_a,K_b,K_c]);
        
        Basis_set_Hkp_init_Rm = kp_init_gen(D_kp,Basis_set_Rm);
        %     Basis_set_Hkp_init_Gk = kp_init_gen(D_kp,Basis_set_Gk);
        %     H_kp_init_Rm = H_kp_temp_gen(Basis_set_Hkp_init_Rm,D_kp);
        %     H_kp_init_temp_Rm_conj = H_kp_init_Rm';
        %     H_kp_init_Gk = H_kp_temp_gen(Basis_set_Hkp_init_Gk,D_kp);
        %% zero H' = H
        
    elseif strcmp(mode,'halfkp')
        Basis_set_Rm= k_sym_list_gen(k_order,[K_x,K_y,K_z]);
        Basis_set_Hkp_init_Rm = kp_init_gen(D_kp,Basis_set_Rm,mat_stencil);
    elseif strcmp(mode,'?')
    else
    end
    
    
%%  main functional   
    Basis_set_Hkp_Rm_temp = Basis_set_Hkp_init_Rm;
    Basis_set_Hkp_Gk_temp = Basis_set_Hkp_Rm2Gk(Basis_set_Hkp_Rm_temp,[K_x,K_y,K_z],[K_a,K_b,K_c]);
%     for i=1:1
%         %% -----------*-----init-----*-----------------
%         H_kp_init_Rm = H_kp_temp_gen(Basis_set_Hkp_Rm_temp,D_kp);
%         H_kp_init_temp_Rm_hermi = H_kp_init_Rm'; % neverneed
%         H_kp_init_Gk = H_kp_temp_gen(Basis_set_Hkp_Gk_temp,D_kp);
%         %  -----------*-----init-----*-----------------
%         %% Operation on Realspace
%         H_kp_temp = H_kp_init_temp_Rm_hermi;
%         Basis_set_Hkp_temp =  H_kp_temp2Basis_set_Hkp(H_kp_temp,Basis_set_Rm);
%         % Operation on Realspace kspace 
%         [K_a,K_b,K_c]=K_xyzintoK_123(eye(3),[K_x,K_y,K_z]);
%         H_kp_temp_2 = subs(H_kp_init_Gk);
%          Basis_set_Hkp_temp2 =  H_kp_temp2Basis_set_Hkp(H_kp_temp_2,Basis_set_Rm);
%         % solve and subs
%         Basis_set_Hkp_out = Basis_set_Hkp_compare(Basis_set_Hkp_temp,Basis_set_Hkp_temp2);
%         %% -----------*-----back-----*----------------- 
%         Basis_set_Hkp_Rm_temp = Basis_set_Hkp_out;
%         Basis_set_Hkp_Gk_temp = Basis_set_Hkp_Rm2Gk(Basis_set_Hkp_Rm_temp,[K_x,K_y,K_z],[K_a,K_b,K_c]);
%         %  -----------*-----back-----*-----------------
%     end
    %% first check ALO
%     Basis_set_Hkp_Rm_temp = Basis_set_Hkp_init_Rm;
%     Basis_set_Hkp_Gk_temp = Basis_set_Hkp_init_Gk;
    for i=1:nALO
        fprintf('************************************* The %d th ALO operation *************************************\n',i);
        fprintf('Operation of H \n');
        disp(ALO_set{1,i});
        fprintf('Operation of k \n');
        disp(ALO_set{2,i});
        %% -----------*-----init-----*-----------------
        H_kp_init_Rm = H_kp_temp_gen(Basis_set_Hkp_Rm_temp,D_kp);
        H_kp_init_temp_Rm_conj = conj(H_kp_init_Rm);
        H_kp_init_Gk = H_kp_temp_gen(Basis_set_Hkp_Gk_temp,D_kp);
        %  -----------*-----init-----*-----------------
        %% Operation on Realspace
%         disp(ALO_set{1,i});
%         disp(H_kp_init_temp_Rm_conj);
%         disp(ALO_set{1,i}');
        H_kp_temp = ALO_set{1,i}*H_kp_init_temp_Rm_conj*ALO_set{1,i}';
        disp('distrbute H ********************************************************');
        disp(H_kp_temp);
        Basis_set_Hkp_temp =  H_kp_temp2Basis_set_Hkp(H_kp_temp,Basis_set_Rm);
        % Operation on Realspace kspace 
        [K_a,K_b,K_c]=K_xyzintoK_123(ALO_set{2,i},[K_x,K_y,K_z]);
%         disp([K_a,K_b,K_c]);
%         disp(H_kp_init_Gk)
        H_kp_temp_2 = subs(H_kp_init_Gk);
        disp('distrbute k ********************************************************');
        disp(H_kp_temp_2);
         
         Basis_set_Hkp_temp2 =  H_kp_temp2Basis_set_Hkp(H_kp_temp_2,Basis_set_Rm);
         disp('caculating >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>');
        % solve and subs
        Basis_set_Hkp_out = Basis_set_Hkp_compare(Basis_set_Hkp_temp,Basis_set_Hkp_temp2);
        %% -----------*-----back-----*----------------- 
        clear('K_x', 'K_y',' K_z', 'K_a','K_b', 'K_c');
        syms K_x K_y K_z real;
        syms K_a K_b K_c real;
        Basis_set_Hkp_Rm_temp = Basis_set_Hkp_out;
        Basis_set_Hkp_Gk_temp = Basis_set_Hkp_Rm2Gk(Basis_set_Hkp_Rm_temp,[K_x,K_y,K_z],[K_a,K_b,K_c]);
        
        %  -----------*-----back-----*-----------------
    end
    %% next check LO
%     Basis_set_Hkp_Rm_temp = Basis_set_Hkp_init_Rm;
%     Basis_set_Hkp_Gk_temp = Basis_set_Hkp_init_Gk;
    for i=1:nLO
        fprintf('************************************* The %d th LO operation *************************************\n',i);
        fprintf('Operation of H \n');
        disp(LO_set{1,i});
        fprintf('Operation of k \n');
        disp(LO_set{2,i});
        %% -----------*-----init-----*-----------------
        H_kp_init_Rm = H_kp_temp_gen(Basis_set_Hkp_Rm_temp,D_kp);
        % H_kp_init_temp_Rm_conj = H_kp_init_Rm'; % neverneed
        H_kp_init_Gk = H_kp_temp_gen(Basis_set_Hkp_Gk_temp,D_kp);
        %  -----------*-----init-----*-----------------
        %% Operation on Realspace
        H_kp_temp = LO_set{1,i}*H_kp_init_Rm*LO_set{1,i}';
        %disp(H_kp_temp);
        disp('distrbute H ********************************************************');
        Basis_set_Hkp_temp =  H_kp_temp2Basis_set_Hkp(H_kp_temp,Basis_set_Rm);
        %disp(Basis_set_Hkp_temp(4).coe_i(2,1));
        % Operation on kspace 
        [K_a,K_b,K_c]=K_xyzintoK_123(LO_set{2,i},[K_x,K_y,K_z]);
%         disp([K_a,K_b,K_c]);
%         disp(H_kp_init_Gk);
        H_kp_temp_2 = subs(H_kp_init_Gk);
%         disp(H_kp_temp_2);
        disp(H_kp_temp);
        %disp(H_kp_temp_2);
        disp('distrbute k ********************************************************');
        disp(H_kp_temp_2);
        Basis_set_Hkp_temp2 =  H_kp_temp2Basis_set_Hkp(H_kp_temp_2,Basis_set_Rm);
        %disp(Basis_set_Hkp_temp2(4).coe_i(2,1));
        % solve and subs
         disp('caculating >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>');
        Basis_set_Hkp_out = Basis_set_Hkp_compare(Basis_set_Hkp_temp,Basis_set_Hkp_temp2);
        %% -----------*-----back-----*----------------- 
        clear('K_x', 'K_y',' K_z', 'K_a','K_b', 'K_c');
        syms K_x K_y K_z real;
        syms K_a K_b K_c real;
       
        Basis_set_Hkp_Rm_temp = Basis_set_Hkp_out;
        Basis_set_Hkp_Gk_temp = Basis_set_Hkp_Rm2Gk(Basis_set_Hkp_Rm_temp,[K_x,K_y,K_z],[K_a,K_b,K_c]);
        %  -----------*-----back-----*-----------------
    end


    %% print and gen a bitch of files    
    disp('************************************* Cogratulations *************************************');
    H_kp = simplify(H_kp_temp_gen(Basis_set_Hkp_out,D_kp));
    disp(H_kp)
end
 

%% Output/print function

%% Translation function 
function Basis_set_Hkp_temp =  H_kp_temp2Basis_set_Hkp(H_kp_temp,Basis_set_Rm)
%% a empty Basis_set_Hkp_temp
    [D_kp,~] =  size(H_kp_temp);
    K_type_number = length(Basis_set_Rm);
    Basis_list = [Basis_set_Rm.exp];
%     disp(Basis_list);
    for i = 1:K_type_number
        Basis_set_Hkp_empty(i).exp = Basis_set_Rm(i).exp;
        Basis_set_Hkp_empty(i).str  = Basis_set_Rm(i).str;
        Basis_set_Hkp_empty(i).coe_r = sym( zeros(D_kp,D_kp));
        Basis_set_Hkp_empty(i).coe_i = sym(zeros(D_kp,D_kp));
        Basis_set_Hkp_empty(i).coe = Basis_set_Hkp_empty(i).coe_r + 1i * Basis_set_Hkp_empty(i).coe_i;
    end     
    Basis_set_Hkp_temp = Basis_set_Hkp_empty;
    
    
    H_kp_temp_real = real(expand(simplify(H_kp_temp)));
    H_kp_temp_imag = imag(expand(simplify(H_kp_temp)));
    %disp(H_kp_temp_real);
    %disp(H_kp_temp_imag );
    H_kp_temp_real_str = string(expand(H_kp_temp_real));H_kp_temp_real_str =strrep(H_kp_temp_real_str ,' - ',' + -');
    H_kp_temp_imag_str = string(expand(H_kp_temp_imag));H_kp_temp_imag_str =strrep(H_kp_temp_imag_str ,' - ',' + -');
    for i=1:D_kp
        for j =1:D_kp
            H_kp_temp_real_str_set{i,j} = strsplit(H_kp_temp_real_str{i,j},' + ');
            H_kp_temp_imag_str_set{i,j} = strsplit(H_kp_temp_imag_str{i,j},' + ');
            H_kp_temp_real_exp_set{i,j} = str2sym(H_kp_temp_real_str_set{i,j});
            H_kp_temp_imag_exp_set{i,j} = str2sym(H_kp_temp_imag_str_set{i,j});
            % need a function to check
            %% real part
            for k = 1:length(H_kp_temp_real_exp_set{i,j})
                %disp(H_kp_temp_real_exp_set{i,j}(k)) ;
                [seq,exp_out] = check_in_basis(H_kp_temp_real_exp_set{i,j}(k),Basis_list);
                %disp(seq);
%                 disp(exp_out);
                Basis_set_Hkp_temp(seq).coe_r(i,j) =Basis_set_Hkp_temp(seq).coe_r(i,j) + exp_out;
            end
            
            %% imag part
            for k = 1:length(H_kp_temp_imag_exp_set{i,j})
                %disp(H_kp_temp_imag_exp_set{i,j}(k)) ;
                %disp(H_kp_temp_imag_exp_set{i,j}(k)) ;
                [seq,exp_out] = check_in_basis(H_kp_temp_imag_exp_set{i,j}(k),Basis_list);
                %disp(seq);
                %disp(exp_out);
                Basis_set_Hkp_temp(seq).coe_i(i,j) = Basis_set_Hkp_temp(seq).coe_i(i,j) + exp_out;
            end
        end
    end
    
    for i = 1:K_type_number
        Basis_set_Hkp_temp(i).coe = Basis_set_Hkp_temp(i).coe_r + 1i * Basis_set_Hkp_temp(i).coe_i;
    end
    
    %Basis_set_Hkp_temp =[];
    
%     H_kp_temp = sym(zeros(D_kp));
%     for i = 1:K_type_number
%         H_kp_temp = H_kp_temp + Basis_set_Hkp(i).coe * Basis_set_Hkp(i).exp;
%     end
end
function Basis_set_Gk = Basis_set_Rm2Gk(Basis_set_Rm,K_sym_list_Rm,K_sym_list_Gk)
    K_x = string(K_sym_list_Rm(1));
    K_y = string(K_sym_list_Rm(2));
    K_z = string(K_sym_list_Rm(3)); 
    K_a = string(K_sym_list_Gk(1));
    K_b = string(K_sym_list_Gk(2));
    K_c = string(K_sym_list_Gk(3)); 
    
    Basis_set_Gk = Basis_set_Rm;
    K_type_number = length(Basis_set_Rm);
    for i = 1:K_type_number
        Basis_set_Gk(i).str = strrep(Basis_set_Gk(i).str,K_x,K_a);
        Basis_set_Gk(i).str = strrep(Basis_set_Gk(i).str,K_y,K_b);
        Basis_set_Gk(i).str = strrep(Basis_set_Gk(i).str,K_z,K_c);
        Basis_set_Gk(i).str = regexprep(Basis_set_Gk(i).str,'\d[*]','');
        Basis_set_Gk(i).exp = str2sym(Basis_set_Gk(i).str);      
    end
end
function Basis_set_Hkp_Gk = Basis_set_Hkp_Rm2Gk(Basis_set_Hkp_Rm,K_sym_list_Rm,K_sym_list_Gk)
    K_x = string(K_sym_list_Rm(1));
    K_y = string(K_sym_list_Rm(2));
    K_z = string(K_sym_list_Rm(3)); 
    K_a = string(K_sym_list_Gk(1));
    K_b = string(K_sym_list_Gk(2));
    K_c = string(K_sym_list_Gk(3)); 
    
    Basis_set_Hkp_Gk = Basis_set_Hkp_Rm;
    K_type_number = length(Basis_set_Hkp_Rm);
    for i = 1:K_type_number
        Basis_set_Hkp_Gk(i).str = strrep(Basis_set_Hkp_Gk(i).str,K_x,K_a);
        Basis_set_Hkp_Gk(i).str = strrep(Basis_set_Hkp_Gk(i).str,K_y,K_b);
        Basis_set_Hkp_Gk(i).str = strrep(Basis_set_Hkp_Gk(i).str,K_z,K_c); 
        Basis_set_Hkp_Gk(i).str = regexprep(Basis_set_Hkp_Gk(i).str,'\d[*]','');
        Basis_set_Hkp_Gk(i).exp = str2sym(Basis_set_Hkp_Gk(i).str);      
    end
end
%% Construction function
function Basis_set_Hkp_init = kp_init_gen(D_kp,Basis_set_Rm,mat_stencil)
if nargin<3
    mat_stencil = ones(D_kp);
    mode = 'full';
else
    mode = 'half';
end
    K_type_number = length(Basis_set_Rm);
    
    for i = 1:K_type_number
        Basis_set_Hkp(i).exp = Basis_set_Rm(i).exp;
        Basis_set_Hkp(i).str  = Basis_set_Rm(i).str;
%         disp(Basis_set(i).coestr);
        temp_str = char(strrep(Basis_set_Rm(i).coestr,'_',''));
%         disp(temp_str);
        isvarname(temp_str);
        if strcmp(mode,'full')
            for j = 1:D_kp
                for k = 1:j
                    if j == k
                        Basis_set_Hkp(i).coe_r(j,j) = sym( char(temp_str+"_real_"+string(j)+"_"+string(k)),'real');

                        Basis_set_Hkp(i).coe_i(j,j) = sym(0);
                        Basis_set_Hkp(i).coe(j,j) = Basis_set_Hkp(i).coe_r(j,j) + 1i * Basis_set_Hkp(i).coe_i(j,j);
                    else
                        Basis_set_Hkp(i).coe_r(j,k) = sym( char(temp_str+"_real_"+string(j)+"_"+string(k)),'real');

                        Basis_set_Hkp(i).coe_i(j,k) = sym(char( temp_str+"_imag_"+string(j)+"_"+string(k)),'real');

                        Basis_set_Hkp(i).coe(j,k) = Basis_set_Hkp(i).coe_r(j,k) + 1i * Basis_set_Hkp(i).coe_i(j,k);
                        Basis_set_Hkp(i).coe_r(k,j) = sym( char(temp_str+"_real_"+string(j)+"_"+string(k)),'real');
                        Basis_set_Hkp(i).coe_i(k,j) = -sym(char( temp_str+"_imag_"+string(j)+"_"+string(k)),'real');
                        Basis_set_Hkp(i).coe(k,j) = Basis_set_Hkp(i).coe_r(k,j) + 1i * Basis_set_Hkp(i).coe_i(k,j);
                    end
                end
            end
        elseif strcmp(mode,'half')
            for j = 1:D_kp
                for k = 1:j
                    if j == k
                        if mat_stencil(j,k) == 1
                            Basis_set_Hkp(i).coe_r(j,j) = sym( char(temp_str+"_real_"+string(j)+"_"+string(k)),'real');
                            Basis_set_Hkp(i).coe_i(j,j) = sym(0);                        
                        else
                            Basis_set_Hkp(i).coe_r(j,j) = sym(0);
                            Basis_set_Hkp(i).coe_i(j,j) = sym(0);
                        end
                        Basis_set_Hkp(i).coe(j,j) = Basis_set_Hkp(i).coe_r(j,j) + 1i * Basis_set_Hkp(i).coe_i(j,j);
                    else
                        if mat_stencil(j,k) == 1
                            Basis_set_Hkp(i).coe_r(j,k) = sym( char(temp_str+"_real_"+string(j)+"_"+string(k)),'real');
                            Basis_set_Hkp(i).coe_i(j,k) = sym(char( temp_str+"_imag_"+string(j)+"_"+string(k)),'real');                           
                            Basis_set_Hkp(i).coe_r(k,j) = sym( char(temp_str+"_real_"+string(j)+"_"+string(k)),'real');
                            Basis_set_Hkp(i).coe_i(k,j) = -sym(char( temp_str+"_imag_"+string(j)+"_"+string(k)),'real');
                        else
                            Basis_set_Hkp(i).coe_r(j,k) = sym(0);
                            Basis_set_Hkp(i).coe_i(j,k) = sym(0);                        
                            Basis_set_Hkp(i).coe_r(k,j) = sym(0);
                            Basis_set_Hkp(i).coe_i(k,j) = sym(0);
                        end
                        Basis_set_Hkp(i).coe(j,k) = Basis_set_Hkp(i).coe_r(j,k) + 1i * Basis_set_Hkp(i).coe_i(j,k);
                        Basis_set_Hkp(i).coe(k,j) = Basis_set_Hkp(i).coe_r(k,j) + 1i * Basis_set_Hkp(i).coe_i(k,j);
                    end
                end
            end
        end
            
    end     
    Basis_set_Hkp_init = Basis_set_Hkp;
end
function Basis_set= k_sym_list_gen(k_order,K_sym_list)
K_x = K_sym_list(1);
K_y = K_sym_list(2);
K_z = K_sym_list(3);
% syms Const real;
%% basic sym
Basis  = expand((1+K_x+K_y+K_z)^k_order);
Basis_str = string(Basis);
Basis_str =strrep(Basis_str ,' - ',' +- ');

Basis_str_set = strsplit(Basis_str,'+');
K_type_number = length(Basis_str_set);
for i = 1:K_type_number 
Basis_set(i).str = Basis_str_set(i);
newStr = regexprep(Basis_set(i).str,'\d[*]','');
Basis_set(i).exp = str2sym(newStr);
% if Basis_set(i).exp == 1
%     Basis_set(i).exp = Const;
% end
newStr = regexprep(Basis_set(i).str,'\d[*]','');
newStr = strrep(strrep( strrep( newStr,'1','Const'),'^','o'),'*','t');
%disp(newStr);
Basis_set(i).coestr = strrep( newStr,' ','');
end

end
function H_kp_temp =  H_kp_temp_gen(Basis_set_Hkp,D_kp)
    K_type_number = length(Basis_set_Hkp);
    H_kp_temp = sym(zeros(D_kp));
    for i = 1:K_type_number
        H_kp_temp = H_kp_temp + Basis_set_Hkp(i).coe * Basis_set_Hkp(i).exp;
    end
end
%% Tools function
function Basis_set_Hkp_out = Basis_set_Hkp_compare(Basis_set_Hkp_temp,Basis_set_Hkp_temp2)
    Basis_set_Hkp_out = Basis_set_Hkp_temp;
    K_type_number = length(Basis_set_Hkp_out);
    [D_kp,~] =  size(Basis_set_Hkp_out(1).coe);
    
    %eval_mat(D_kp,D_kp,K_type_number,2)  = ''; 
    effect_run = 0;
    for k = 1:K_type_number
        for i = 1:D_kp
            for j = 1:D_kp
                %% real 
                
                temp_exp = Basis_set_Hkp_temp(k).coe_r(i,j) - Basis_set_Hkp_temp2(k).coe_r(i,j);
%                 disp(Basis_set_Hkp_temp(k).coe_r(i,j));
%                 disp(Basis_set_Hkp_temp2(k).coe_r(i,j));
%                 disp(temp_exp );
                if temp_exp ~=0
                    [evalstr1,eval_mat(i,j,k,1) ] = solve_one(temp_exp);
                    eval(evalstr1);
                    effect_run = 1;
                    %eval(evalstr2);
                end
                %% imag
                temp_exp = Basis_set_Hkp_temp(k).coe_i(i,j) - Basis_set_Hkp_temp2(k).coe_i(i,j);
                %disp(temp_exp );
%                  disp([Basis_set_Hkp_temp(k).exp i j]);
%                  disp(Basis_set_Hkp_temp(k).coe_i(i,j));
%                  disp(Basis_set_Hkp_temp2(k).coe_i(i,j));
%                 disp(temp_exp );
                if temp_exp ~=0
                    [evalstr1,eval_mat(i,j,k,2) ] = solve_one(temp_exp);
                    eval(evalstr1);
                    effect_run = 1;
                    %eval(evalstr2);                   
                end
            end
        end
    end
    %% fixbug
    if effect_run == 1
        %disp(eval_mat);
        eval_list = rmmissing(reshape(eval_mat,[],1));
        %disp(eval_list);
        
        eval_mat = strsplit_list(eval_list,'=');
        clean_label = find(eval_mat(:,2) == "0");
        eval_mat(clean_label,:) = [];
        eval_list2= eval_list_after_gen(eval_mat );
        eval_list = [eval_list;eval_list2];
        eval_list = zero_first(eval_list );
        
        %% eval list
        for i =1 :length(eval_list)
            %disp(eval_list(i,:));
            eval(strcat(eval_list(i,:),";"));
        end
    else
        disp('The operation does nothing');
    end
    %% subs
    save('subs_temp.mat');
    Basis_set_Hkp_out = subsall_Basis_set_Hkp(Basis_set_Hkp_out);
end
function eval_list_out = zero_first(eval_list )
    %disp(eval_list);
    eval_list_out = [];
    for i = 1:length(eval_list )
         if contains(eval_list(i),"=")
            %disp(eval_list(i));
                eval_mat = strsplit(eval_list(i),'=');
                %disp(eval_mat(2))
            if strcmp(eval_mat(2),"0")
                eval_list_out = [eval_list(i);eval_list_out ];
            else
                eval_list_out = [eval_list_out;eval_list(i) ];
            end
         end
    end
end
function eval_list_after= eval_list_after_gen(str_mat )
    basis_str = unique(str_mat(:,1));
    eval_list_after = '';
    for i =1:length(basis_str)
        seqlist = find(str_mat(:,1)==basis_str(i));
        if length(seqlist) == 1
%             disp('do nothing');
            continue;
        elseif length(seqlist) == 2
%             disp(str_mat);
%             disp(basis_str(i));
%             if contains(str_mat(seqlist(1),2),'-')
%                 temp_str = strrep(str_mat(seqlist(1),2),'-','')+"= - "+str_mat(seqlist(2),2);
%             else
%                 
%             end
            if strcmp(str_mat(seqlist(1),2),str_mat(seqlist(2),2))
                continue;
            end
            temp_str = str_mat(seqlist(1),2)+"-"+str_mat(seqlist(2),2);
            [~,temp_str  ] = solve_one(str2sym(temp_str));
            eval_list_after = [eval_list_after;temp_str];
            continue;
        elseif length(seqlist) > 2
            for j = 2:length(seqlist)
                for k = 1:j-1
                    if strcmp(str_mat(seqlist(k),2),str_mat(seqlist(j),2))
                        continue;
                    end
                    temp_str = str_mat(seqlist(k),2)+"-"+str_mat(seqlist(j),2);
                    [~,temp_str  ] = solve_one(str2sym(temp_str));
                    eval_list_after = [eval_list_after;temp_str];
                end
            end
        else
           disp('sonething wrong'); 
        end
        
    end
end
function str_list = strsplit_list(str_list_origin,strlabel)
    [size1,size2] = size(str_list_origin);
    str_list = [];
    for i = 1:size1
        tempstrlist1 = []; 
        for j = 1:size2
            tempstrlist = strsplit(str_list_origin(i,j),strlabel);
            tempstrlist1 = [tempstrlist1,tempstrlist];
        end
        str_list = [str_list; tempstrlist1  ]; 
    end
    
end
function Basis_set_Hkp_out = subsall_Basis_set_Hkp(Basis_set_Hkp_out)
    load('subs_temp.mat');
    K_type_number = length(Basis_set_Hkp_out);
    for k = 1:K_type_number

                Basis_set_Hkp_out(k).coe_r = subs(Basis_set_Hkp_out(k).coe_r);
                Basis_set_Hkp_out(k).coe_i = subs(Basis_set_Hkp_out(k).coe_i);
                Basis_set_Hkp_out(k).coe = subs(Basis_set_Hkp_out(k).coe);

    end
end
function [evalstr1,evalstr2] = solve_one(temp_exp)
    solvevar = symvar(temp_exp);
    if length(temp_exp) >2
        disp(solvevar);
        warning('cant solve this, sorry');
    end
    solveout = solve(temp_exp);
    
%     disp([solvevar solveout]);
    if length (solvevar) == 2
        evalstr1 = "syms "+string(solvevar(1))+" "+string(solvevar(2))+" real";
    elseif length (solvevar) == 1
%         disp('test it');
%         disp(temp_exp);
        evalstr1 = "syms "+string(solvevar)+" real";
    else
        disp(check);
        evalstr1 = "syms "+ " ";
        for i =1:length(temp_exp)
             evalstr1 =  evalstr1 + string(solvevar(i))+" ";
        end
        evalstr1 =  evalstr1 + " real";
    end
    evalstr2 = string(solvevar(1))+'='+string(solveout);
    %disp(evalstr2);
end

function [seq,exp_out] = check_in_basis(temp_exp,Basis_list)
    seq= 0;
    Basis_list_str = string(Basis_list);
    Basis_list_str = strrep(Basis_list_str,'1','Constant');
    %disp(Basis_list_str );
    for i =1:length(Basis_list_str)
        temp_exp_after = temp_exp/Basis_list(i);
        temp_exp_after_str= string(simplify(temp_exp_after));
        %disp(temp_exp_after_str);
        
        temp_contain_label = strlist_find_in_str(Basis_list_str,temp_exp_after_str);
        if temp_contain_label == 1
            seq = i;
            exp_out = temp_exp_after;
%             break
        end
    end
    if seq == 0
        error('error!!!');
    end
end
function label = strlist_find_in_str(temp_str_list,temp_str)
    label = 1;
%     disp([temp_str_list temp_str]);

    for i = 1:length(temp_str_list)
%         disp('right?');
        if contains(temp_str,temp_str_list(i))
            label = 0;
%             warning('label = 0');
            break
        end
    end
end

function [K1,K2,K3]=K_xyzintoK_123(Operation,K_sym_list_Rm)
    K = Operation*K_sym_list_Rm';
    K1= K(1);
    K2= K(2);
    K3= K(3);
end
