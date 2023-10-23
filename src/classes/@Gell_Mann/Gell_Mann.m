classdef Gell_Mann < pauli_matric 
    %pauli matric 
    %   此处显示详细说明
    
    properties
       
    end
    
    methods
        function GM = Gell_Mann( label1 ,label2,options )
            arguments
               label1 double =0;
               label2 double =0;
               options.rep = 'norm';
            end
            GM = GM@pauli_matric(0);
            paulis =  pauli_matric();
            pauli_mat = zeros(3);
            if nargin  ==1
                switch options.rep
                    case 'norm'
                        switch label1
                            case 0
                                pauli_mat = eye(3);
                            case 1
                                pauli_mat(1:2,1:2) = paulis(2).mat;
                            case 2
                                pauli_mat(1:2,1:2) = paulis(3).mat;
                            case 3
                                pauli_mat(1:2,1:2) = paulis(4).mat;
                            case 4
                                pauli_mat([1,3],[1,3]) = paulis(2).mat;
                            case 5
                                pauli_mat([1,3],[1,3]) = paulis(3).mat;
                            case 6
                                pauli_mat(2:3,2:3) = paulis(2).mat;
                            case 7
                                pauli_mat(2:3,2:3) = paulis(3).mat;
                            case 8
                                pauli_mat = diag([1 1 -2]./sqrt(3));
                        end
                    otherwise
                        switch label1
                            case 1
                                pauli_mat(1:2,1:2) = paulis(2).mat;
                            case 2
                                pauli_mat(1:2,1:2) = paulis(3).mat;
                            case 3
                                pauli_mat(1:2,1:2) = paulis(4).mat;
                            case 4
                                pauli_mat([1,3],[1,3]) = paulis(2).mat;
                            case 5
                                pauli_mat([1,3],[1,3]) = paulis(3).mat;
                            case 6
                                pauli_mat(2:3,2:3) = paulis(2).mat;
                            case 7
                                pauli_mat(2:3,2:3) = paulis(3).mat;
                            case 8
                                pauli_mat = diag([1 1 -2]./sqrt(3));
                        end
                end
                GM.mat = pauli_mat;
            elseif nargin == 0
                count = 0;
                for i = 1:8
                    count = count +1;
                    GM(count)=Gell_Mann(i,'rep',options.rep);
                end
            else
                %wait
                GM = gamma_matric(label1,'rep',options.rep);
                GM =  GM.commute(Gell_Mann(label2,'rep',options.rep));
                GM.mat = GM.mat/(2*1i);
            end
        end
    end
    methods(Static)
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

