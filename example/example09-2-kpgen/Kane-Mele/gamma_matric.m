classdef gamma_matric
    %pauli matric 
    %   此处显示详细说明
    
    properties
       mat;
    end
    
    methods
        function GM = gamma_matric(label1,label2)
            if nargin <1
                count = 0;
                for i = 0:5
                    count = count +1;
                    GM(count)=gamma_matric(i);
                end
                for i = 2:5
                    for j =1:i-1
                        count = count +1;
                        GM(i)=gamma_matric(i,j);
                    end
                end
            elseif nargin <2
                switch label1
                    case 0
                        GM.mat = eye(4);
                    case 1
                        GM.mat = kron(eye(2),[1,0;0,-1]);
                    case 2
                        GM.mat = kron([0,1;1,0],[0,1;1,0]);
                    case 3
                        GM.mat = kron(eye(2),[0,1;1,0]);
                    case 4
                        GM.mat = kron([0,1;1,0],[0,-1i;1i,0]);
                    case 5
                        GM.mat = kron([1,0;0,-1],[0,-1i;1i,0]);
                end
            else
                GM = gamma_matric(label1);
                GM =  GM.commute(gamma_matric(label2));
                GM.mat = GM.mat/(2*1i);
            end
        end
        
        function mat = disp(GM)
            %METHOD1 此处显示有关此方法的摘要
            %   此处显示详细说明
            for i=1:length(GM)
                mat(:,:,i) = GM(i).mat;
                disp(GM(i).mat);
            end
        end
        function GM = anticommute(GM1,GM2)
            GM = GM1;
            GM.mat = (GM1.mat *  GM2.mat + ...
                GM2.mat *  GM1.mat);
        end
        function GM = commute(GM1,GM2)
            GM = GM1;
            GM.mat = (GM1.mat *  GM2.mat - ...
                GM2.mat *  GM1.mat);
        end
      function GM = plus(GM1,GM2)
         GM = GM1;
         zeromat =zeros(length(GM.mat));
         GM.mat = [GM.mat,zeromat ;zeromat ,GM2.mat];
      end
      function GM = minus(GM1,GM2)
         GM = GM1;
         GM.mat = GM.mat - GM2.mat;
      end
      function GM = times(GM1,GM2)
          GM = GM1;
          GM.mat = GM.mat * GM2.mat;
      end
      function GM = mtimes(GM1,GM2)
          GM = GM1;
          GM.mat = kron(GM.mat , GM2.mat);
      end
      function mat = double(GM)
          for i=1:length(GM)
              mat(:,:,i) = GM(i).mat;
          end
      end
    end
end

