classdef pauli_matric
    %pauli matric 
    %   此处显示详细说明
    
    properties
       mat;
    end
    
    methods
        function PM = pauli_matric(label)
            if nargin <1
                for i = 1:4
                    PM(i)=pauli_matric(i-1);
                end
            else
                switch label
                    case {'0','I',0}
                        PM.mat = [1,0;0,1];
                    case {'x','1','X',1}
                        PM.mat = [0,1;1,0];
                    case {'y','2','Y',2}
                        PM.mat = [0,-1i;1i,0];
                    case {'z','3','Z',3}
                        PM.mat = [1,0;0,-1];
                end
            end
        end
    end
    %% overload
    methods
        function mat = disp(PM)
            %METHOD1 此处显示有关此方法的摘要
            %   此处显示详细说明
            for i=1:length(PM)
                mat(:,:,i) = PM(i).mat;
                disp(PM(i).mat);
            end
        end
      function PM = plus(PM1,PM2)
          if isa(PM1,'pauli_matric') && isa(PM2,'pauli_matric')
              PM = PM1;
              PM.mat = PM.mat+PM2.mat;
          elseif isa(PM1,'pauli_matric') && ~isa(PM2,'pauli_matric')
              PM = PM1;
              if length(PM1.mat) == length(PM2)
                PM.mat = PM.mat + PM2;
              elseif length(PM2) == 1 
                PM.mat = PM.mat + PM2*eye(length(PM.mat));
              else
                error('+ obj error');
              end
          elseif ~isa(PM1,'pauli_matric') && isa(PM2,'pauli_matric')
              PM = PM2;
              if length(PM2.mat) == length(PM1)
                  PM.mat = PM.mat + PM1;
              elseif length(PM1) == 1
                  PM.mat = PM.mat + PM1*eye(length(PM.mat));
              else
                  error('+ obj error');
              end
          else
              error('obj + obj error');
          end
      end
      function PM = minus(PM1,PM2)
          if isa(PM1,'pauli_matric') && isa(PM2,'pauli_matric')
              PM = PM1;
              PM.mat = PM.mat - PM2.mat;
          elseif isa(PM1,'pauli_matric') && ~isa(PM2,'pauli_matric')
              PM = PM1;
              if length(PM1.mat) == length(PM2)
                PM.mat = PM.mat - PM2;
              elseif length(PM2) == 1 
                PM.mat = PM.mat -  PM2*eye(length(PM.mat));
              else
                error('+ obj error');
              end
          elseif ~isa(PM1,'pauli_matric') && isa(PM2,'pauli_matric')
              PM = PM2;
              if length(PM2.mat) == length(PM1)
                  PM.mat =  PM1 - PM.mat ;
              elseif length(PM1) == 1
                  PM.mat =  PM1*eye(length(PM.mat))- PM.mat ;
              else
                  error('obj + error');
              end
          else
              error('obj + obj error');
          end
      end
      function PM = uminus(PM)
          PM.mat = -PM.mat;
      end
      function PM = times(PM1,PM2)
          if isa(PM1,'pauli_matric') && isa(PM2,'pauli_matric')
              PM = PM1;
              PM.mat = PM.mat * PM2.mat;
          elseif isa(PM1,'pauli_matric') && ~isa(PM2,'pauli_matric')
              PM = PM1;
              if length(PM2.mat) == length(PM1)
                  PM.mat = kron(PM.mat , PM2);
              elseif length(PM1) == 1
                  PM.mat =  kron(PM.mat , PM2*eye(length(PM.mat)));
              else
                  error('.* obj error');
              end
          elseif ~isa(PM1,'pauli_matric') && isa(PM2,'pauli_matric')
              PM = PM2;
              if length(PM2.mat) == length(PM1)
                  PM.mat = kron(P<1 ,PM.mat);
              elseif length(PM1) == 1
                  PM.mat =  kron( PM1*eye(length(PM.mat)),PM.mat );
              else
                  error('.* obj error');
              end
          end
          
      end
      function PM = mtimes(PM1,PM2)
          if isa(PM1,'pauli_matric') && isa(PM2,'pauli_matric')
              PM = PM1;
              PM.mat = kron(PM.mat , PM2.mat);
          elseif isa(PM1,'pauli_matric') && ~isa(PM2,'pauli_matric')
              PM = PM1;
              if length(PM2) == length(PM1.mat)
                  PM.mat =PM.mat * PM2;
              elseif length(PM2) == 1
                  PM.mat =  PM.mat * PM2;
              else
                  error('* obj error');
              end
          elseif ~isa(PM1,'pauli_matric') && isa(PM2,'pauli_matric')
              PM = PM2;
              if length(PM2.mat) == length(PM1)
                  PM.mat = PM1 *PM.mat;
              elseif length(PM1) == 1
                  PM.mat =  PM1*PM.mat ;
              else
                  error(' obj * error');
              end
          end
          
      end
      function PM = mldivide(PM1,PM2)
          if isa(PM1,'pauli_matric') && isa(PM2,'pauli_matric')
              PM = PM1;
              PM.mat = PM.mat\PM2.mat;
          elseif isa(PM1,'pauli_matric') && ~isa(PM2,'pauli_matric')
              PM = PM1;
              if length(PM2.mat) == length(PM1)
                  PM.mat =PM.mat \ PM2;
              elseif length(PM1) == 1
                  PM.mat =  PM.mat \ PM2*eye(length(PM.mat));
              else
                  error('\ obj error');
              end
          elseif ~isa(PM1,'pauli_matric') && isa(PM2,'pauli_matric')
              PM = PM2;
              if length(PM2.mat) == length(PM1)
                  PM.mat = PM1 \PM.mat;
              elseif length(PM1) == 1
                  PM.mat =  PM1\eye(length(PM.mat))*PM.mat ;
              else
                  error(' obj \ error');
              end
          end
          
      end
      function PM = mrdivide(PM1,PM2)
          if isa(PM1,'pauli_matric') && isa(PM2,'pauli_matric')
              PM = PM1;
              PM.mat = PM.mat/PM2.mat;
          elseif isa(PM1,'pauli_matric') && ~isa(PM2,'pauli_matric')
              PM = PM1;
              if length(PM2.mat) == length(PM1)
                  PM.mat =PM.mat / PM2;
              elseif length(PM1) == 1
                  PM.mat =  PM.mat* (1/ PM2);
              else
                  error('/ obj error');
              end
          elseif ~isa(PM1,'pauli_matric') && isa(PM2,'pauli_matric')
              PM = PM2;
              if length(PM2.mat) == length(PM1)
                  PM.mat = PM1 /PM.mat;
              elseif length(PM1) == 1
                  PM.mat =  PM1/eye(length(PM.mat))*PM.mat ;
              else
                  error(' obj / error');
              end
          end
          
      end
      function mat = double(PM)
          if length(PM) == 1
              mat = PM.mat;
          else  
              for i=1:length(PM)
                  mat(:,:,i) = PM(i).mat;
              end
          end
      end
      function symmat = sym(PM)
          if length(PM) == 1
              mat = PM.mat;
          else
              for i=1:length(PM)
                  mat(:,:,i) = PM(i).mat;
              end
          end
          symmat = sym(mat);
      end
      function PM = inv(PM)
          PM.mat = inv(PM.mat);
      end
      function PM = expm(PM)
          PM.mat = expm(PM.mat);
      end
      function PM = ctranspose(PM)
          PM.mat = PM.mat';
      end
      function PM = transpose(PM)
          PM.mat = PM.mat.';
      end
      function PM = anticommute(PM1,PM2)
          PM = PM1;
          PM.mat = (PM1.mat *  PM2.mat + ...
              PM2.mat *  PM1.mat);
      end
      function PM = commute(PM1,PM2)
          PM = PM1;
          PM.mat = (PM1.mat *  PM2.mat - ...
              PM2.mat *  PM1.mat);
      end
      function PM =  colon(PM1,PM2)
          PM = anticommute(PM1,PM2);
      end
      function PM = vertcat(PM1,PM2)
         PM = PM1;
         zeromat =zeros(length(PM.mat));
         PM.mat = [PM.mat,zeromat ;zeromat ,PM2.mat];
      end
      function PM = horzcat(PM1,PM2)
          PM = commute(PM1,PM2);
      end
      function result = iscommute(PM1,PM2)
          mattmp = PM2*PM1*PM2'/PM1;
          tf = isdiag(mattmp);
          switch tf
              case 1
                  eigentmp = diag(mattmp);
                  if all (~(diff(eigentmp))) %判断a中元素是否全部相等 若相等，则输出1，若不相等，则输出0
                      result = eigentmp(1);
                  else
                      result = eigentmp;
                      disp('error')
                  end
                  
              case 0
                  result = nan;
                  disp('error')
          end


      end
    end
    %%
    methods(Static)
        function Smat_inv = S()
            Smat = zeros(4);
            Pauli_L = sym(pauli_matric);
            for i = 1:4
                tmp_mat = Pauli_L(:,:,i);
                tmp_mat_r = real(tmp_mat);
                tmp_mat_i = imag(tmp_mat);
                Smat(i,1:2) = diag(tmp_mat_r);
                Smat(i,3) = diag(tmp_mat_r,1);
                Smat(i,4) = diag(tmp_mat_i,1);
            end
            Smat_inv = inv(Smat);
        end
        function Pauli_L = L()
            count = 0;
            Pauli_L = sym(zeros(1,4));
            for i = 0:3
                count = count +1;
                Pauli_L(count)=sym(['Pauli_',num2str(i)]);
            end

        end
    end
end

