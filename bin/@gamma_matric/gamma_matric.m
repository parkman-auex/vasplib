classdef gamma_matric < pauli_matric 
    %pauli matric 
    %   此处显示详细说明
    
    properties
       
    end
    
    methods
        function GM = gamma_matric( label1 ,label2,options )
            arguments
               label1 double =0;
               label2 double =0;
               options.rep = 'ZYX';
            end
            GM = GM@pauli_matric(0);
            paulis =  pauli_matric();
            if nargin  ==1
                switch options.rep
                    case 'ZYX'
                        switch label1
                            case 0
                                pauli = paulis(1)*paulis(1);
                            case 1
                                pauli = paulis(1)*paulis(2);
                            case 2
                                pauli = paulis(3)*paulis(3);
                            case 3
                                pauli = paulis(1)*paulis(4);
                            case 4
                                pauli = paulis(2)*paulis(3);
                            case 5
                                pauli = paulis(4)*paulis(3);
                        end
                    case 'BJY'
                        switch label1
                            case 0
                                pauli = paulis(1)*paulis(1);
                            case 1
                                pauli = paulis(1)*paulis(2);
                            case 2
                                pauli = paulis(3)*paulis(3);
                            case 3
                                pauli = paulis(1)*paulis(4);
                            case 4
                                pauli = paulis(2)*paulis(3);
                            case 5
                                pauli = paulis(4)*paulis(3);
                        end
                    case 'ZSC'
                        switch label1
                            case 0
                                pauli = paulis(1)*paulis(1);
                            case 1
                                pauli = paulis(2)*paulis(2);
                            case 2
                                pauli = paulis(3)*paulis(2);
                            case 3
                                pauli = paulis(4)*paulis(2);
                            case 4
                                pauli = paulis(1)*paulis(3);
                            case 5
                                pauli = paulis(1)*paulis(4);
                        end
                    case 'KM'
                        switch label1
                            case 0
                                pauli = paulis(1)*paulis(1);
                            case 1
                                pauli = paulis(1)*paulis(2);
                            case 2
                                pauli = paulis(1)*paulis(4);
                            case 3
                                pauli = paulis(2)*paulis(3);
                            case 4
                                pauli = paulis(3)*paulis(3);
                            case 5
                                pauli = paulis(4)*paulis(3);
                        end
                    case 'Dirac'
                        switch label1
                            case 0
                                pauli = paulis(4)*paulis(1);
                            case 1
                                pauli = 1i*paulis(3)*paulis(2);
                            case 2
                                pauli = 1i*paulis(3)*paulis(3);
                            case 3
                                pauli = 1i*paulis(3)*paulis(4);
                            case 4
                                pauli = paulis(1)*paulis(1);
                            case 5
                                pauli = paulis(2)*paulis(1);
                        end
                    case 'Weyl'
                        switch label1
                            case 0
                                pauli = paulis(2)*paulis(1);
                            case 1
                                pauli = 1i*paulis(3)*paulis(2);
                            case 2
                                pauli = 1i*paulis(3)*paulis(3);
                            case 3
                                pauli = 1i*paulis(3)*paulis(4);
                            case 4
                                pauli = paulis(1)*paulis(1);
                            case 5
                                pauli = -paulis(4)*paulis(1);
                        end
                    otherwise
                        switch label1
                            case 0
                                pauli = paulis(1)*paulis(1);
                            case 1
                                pauli = paulis(1)*paulis(2);
                            case 2
                                pauli = paulis(1)*paulis(4);
                            case 3
                                pauli = paulis(2)*paulis(3);
                            case 4
                                pauli = paulis(3)*paulis(3);
                            case 5
                                pauli = paulis(4)*paulis(3);
                        end
                end
                GM.mat = pauli.mat;
            elseif nargin == 0
                count = 0;
                for i = 0:5
                    count = count +1;
                    GM(count)=gamma_matric(i,'rep',options.rep);
                end
                switch options.rep
                    case {'Dirac','Weyl'}
                        for i = [0,1,2,3,5]
                            for j =[0,1,2,3,5]
                                if i < j
                                    count = count +1;
                                    GM(count)=gamma_matric(i,j,'rep',options.rep);
                                end
                            end
                        end
                        
                    otherwise
                        for i = 2:5
                            for j =1:i-1
                                count = count +1;
                                GM(count)=gamma_matric(i,j,'rep',options.rep);
                            end
                        end
                end
            else
                GM = gamma_matric(label1,'rep',options.rep);
                GM =  GM.commute(gamma_matric(label2,'rep',options.rep));
                GM.mat = GM.mat/(2*1i);
            end
        end
    end
    methods(Static)
        function Smat_inv = S()
            Smat = zeros(16);
            GM_L = sym(gamma_matric);
            for i = 1:16
                tmp_mat = GM_L(:,:,i);
                tmp_mat_r = real(tmp_mat);
                tmp_mat_i = imag(tmp_mat);
                Smat(i,1:4) = diag(tmp_mat_r);
                Smat(i,5:7) = diag(tmp_mat_r,1);
                Smat(i,8:9) = diag(tmp_mat_r,2);
                Smat(i,10)  = diag(tmp_mat_r,3);
                Smat(i,11:13) = diag(tmp_mat_i,1);
                Smat(i,14:15) = diag(tmp_mat_i,2);
                Smat(i,16)  = diag(tmp_mat_i,3);
            end
            Smat_inv = inv(Smat);
        end
        function Gamma_L = L()
            count = 0;
            Gamma_L = sym(zeros(1,16));
            for i = 0:5
                count = count +1;
                Gamma_L(count)=sym(['Gamma_',num2str(i)]);
            end
            
            for i = 2:5
                for j =1:i-1
                    count = count +1;
                    Gamma_L(count)=sym(['Gamma_',num2str(j),'_',num2str(i)]);
                end
            end
            
        end
        function Pauli_L = pauli_L()
            Pauli_L = sym(zeros(1,16));
            syms sigma_0 sigma_x sigma_z tau_0 tau_x tau_z real;
            syms sigma_y tau_y ;
            %
            Pauli_L(1) = sigma_0*tau_0;
            Pauli_L(2) = sigma_0*tau_x;
            Pauli_L(3) = sigma_y*tau_y;
            Pauli_L(4) = sigma_0*tau_z;
            Pauli_L(5) = sigma_x*tau_y;
            Pauli_L(6) = sigma_z*tau_y;
            %
            Pauli_L(7)  = -sigma_y*tau_z;
            Pauli_L(8)  =  sigma_0*tau_y;
            Pauli_L(9)  = -sigma_y*tau_x;
            Pauli_L(10) = -sigma_x*tau_z;
            %
            Pauli_L(11) = sigma_z*tau_0;
            Pauli_L(12) = sigma_x*tau_x;
            %
            Pauli_L(13) = -sigma_z*tau_z;
            %
            Pauli_L(14) = -sigma_x*tau_0;
            Pauli_L(15) =  sigma_z*tau_x;
            %
            Pauli_L(16) =  sigma_y*tau_0;
        end
        %
%         function GM = Gamma( varargin )
% 
%         end
    end
    methods
        %         function mat = disp(GM)
        %             %METHOD1 此处显示有关此方法的摘要
%             %   此处显示详细说明
%             for i=1:length(GM)
%                 mat(:,:,i) = GM(i).mat;
%                 disp(GM(i).mat);
%             end
%         end
%         function GM = anticommute(GM1,GM2)
%             GM = GM1;
%             GM.mat = (GM1.mat *  GM2.mat + ...
%                 GM2.mat *  GM1.mat);
%         end
%         function GM = commute(GM1,GM2)
%             GM = GM1;
%             GM.mat = (GM1.mat *  GM2.mat - ...
%                 GM2.mat *  GM1.mat);
%         end
%         function GM =  colon(GM1,GM2)
%             GM = anticommute(GM1,GM2);
%         end
%         function GM = horzcat(GM1,GM2)
%             GM = commute(GM1,GM2);
%         end
%         function GM = plus(GM1,GM2)
%             GM =GM1;
%              GM.mat = GM1.mat +  GM2.mat ;
%         end
%       function GM = vertcat(GM1,GM2)
%          GM = GM1;
%          zeromat =zeros(length(GM.mat));
%          GM.mat = [GM.mat,zeromat ;zeromat ,GM2.mat];
%       end
%       function GM = minus(GM1,GM2)
%          GM = GM1;
%          GM.mat = GM.mat - GM2.mat;
%       end
%       function GM = times(GM1,GM2)
%           GM = GM1;
%           GM.mat = GM.mat * GM2.mat;
%       end
%       function GM = mtimes(GM1,GM2)
%           GM = GM1;
%           GM.mat = kron(GM.mat , GM2.mat);
%       end
%       function mat = double(GM)
%           for i=1:length(GM)
%               mat(:,:,i) = GM(i).mat;
%           end
%       end
    end
end

