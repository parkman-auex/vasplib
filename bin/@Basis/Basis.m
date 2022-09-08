classdef Basis < handle
    % the Basis to describe an orbital
    % Basis fucntion are linerly combined of these function:
    % Ylm Yl^m X LS j mz
    %   此处提供详细说明

    properties
        BasisL {mustBeA(BasisL,["BasisFunc","cell","double","string","char","sym","Spin","Y_lm","Y_l__m","jm","lm"])};
        BcoeL = sym([]);
        BnumL = [];
        orbL = [];
    end
    properties(Dependent)
        Basis_num;
    end
    properties(Hidden)
        Rm = [1 0 0;0 1 0;0 0 1];
        num = false;
        coe = true;
    end
    methods % construction
        function Basisobj = Basis(BasisL,BnumL,BcoeL)
            %BASIS Basis = Basis(BasisFunc(Spin(1/2)) 
            %   此处提供详细说明
            arguments
                BasisL {mustBeA(BasisL,["BasisFunc","cell","double","string","char","sym","Spin","Y_lm","Y_l__m","jm","lm"])};
                BnumL double = [];
                BcoeL sym = sym([]);
            end
            % 
            switch class(BasisL)
                case 'double'
                    if length(BasisL) == 1
                        Basisobj.BasisL = BasisFunc(BasisL);
                    else 
                        % quantumL
                    end
                case {'string','char'}
                    %innerbuilt!
                    Basisobj = Basis.BasisDatabase(BasisL);
                    % not now
                case 'sym'
                    Basisobj.BasisL = BasisFunc(BasisL); % temp
                case {'Spin','jm','lm'}
                    %Basisobj.BasisL = BasisFunc(BasisL);
                    Basisobj.BasisL = (BasisL);
                case {'Y_lm','Y_l__m'}
                    %Basisobj.BasisL = BasisFunc(BasisL);
                    Basisobj.BasisL = (BasisL);
            end
            Basis_num = size(Basisobj.BasisL,1);
            if isempty(BnumL)
                Basisobj.BnumL = ones([1,Basis_num]);
            else
                Basisobj.BnumL = BnumL;
            end
            if isempty(BcoeL)
                Basisobj.BcoeL = sym(ones([1,Basis_num]));
            else
                Basisobj.BcoeL = BcoeL;
            end  
        end
    end
    methods % modify
        
    end
    methods(Static)
        function Basisobj = BasisDatabase(BasisName)
            switch BasisName
                case 'BHZ'
                    Basisobj = Basis([Spin(1/2,1/2);Spin(3/2,3/2,'parity',-1);Spin(1/2,-1/2);Spin(3/2,-3/2,'parity',-1)]);
                case 'Graphene'
                case 'WSM'
            end
        end
    end
    methods % get
        function Basis_num = get.Basis_num(Basisobj)
            Basis_num = length(Basisobj.BasisL);
        end
    end
    methods % disp
        function disp(Basisobj)
            BasisFunction = Basisobj.BasisL;
            builtin('disp',Basisobj) % call builtin
            if Basisobj.coe
                for i =1:size(BasisFunction,1)
                    if  Basisobj.BcoeL(i) == (1)
                        CoeStr(i) = "";
                    elseif  Basisobj.BcoeL(i) == (-1)
                        CoeStr(i) = "-";
                    else
                        CoeStr(i) = string(Basisobj.BcoeL(i));
                    end
                end
            elseif Basisobj.num
                for i =1:size(BasisFunction,1)
                    if  Basisobj.BcoeL(i) == (1)
                        CoeStr(i) = "";
                    elseif  Basisobj.BcoeL(i) == (-1)
                        CoeStr(i) = "-";
                    else
                        CoeStr(i) = string(Basisobj.BcoeL(i));
                    end
                end
            end
            switch class(BasisFunction)
                case 'BasisFunc'
                    for i =1:size(BasisFunction,1)
                        %fprintf('============================================\n')
                        %fprintf('The %d / %d th BasisFunction:\n',i,numel(BasisFunction));
                        fprintf(CoeStr(i) +"{ "+string(BasisFunction(i).BFuncL{1}));
                        for j = 2:length(BasisFunction(i).BFuncL)
                            fprintf(' + ');
                            fprintf(string(BasisFunction(i).BFuncL{j}));
                        end
                        fprintf('}\n')
                    end
                case 'Spin'
                    for i =1:size(BasisFunction,1)
                        fprintf(CoeStr(i) +"{ "+string(BasisFunction(i,:)));
                        fprintf('}\n')
                    end
            end

        end
    end
    % Ugen
    methods
        function Umat = U(Basisobj,rotm,t,rightorleft)
            arguments
                Basisobj Basis;
                rotm ;
                t = [0,0,0];
                rightorleft = 'right';
            end
            if isequal(size(rotm),[3,3])
                [n,theta]= Oper.Rotation2nTheta(rotm,Basisobj.Rm);
                axang = [n theta];
            elseif isequal(size(rotm),[4,1])
                axang = rotm;
            else

            end
            if isequal(t,[0,0,0])
                % pure rotate
                Umat = rotation(Basisobj.BasisL,axang,rightorleft);
            else

            end
        end

    end
    methods(Static)
        function U = Ugen(OperObj)

        end
    end
end