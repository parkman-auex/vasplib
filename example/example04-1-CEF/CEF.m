classdef CEF < vasplib
    %For doing point-charge calculations
    %   LS basis of J basis
    % learn EVENT
    properties
        a
        b
        c
        alpha
        beta
        gamma
    end
    properties
        ion ;
        ionS;
        ionL;
        ionJ;
    end
    properties
        H ;
        O ;
        B ;
    end
    
    %% construction method
    methods
        % --------------------- construction method -----------------------
        function H_CEF = CEF(ion,QuantumNumber,options)
            % HR H_hr = HR(WAN_NUM,vectorL) construct a empty TB obj
            %   H_hr = HR() construct a 4-orbs empty TB obj
            %   H_hr = HR(WAN_NUM) construct a WAN_NUM-orbs empty TB obj
            %   H_hr = HR(WAN_NUM,vectorL) construct a WAN_NUM-orbs with vectorL H_R empty TB obj
            %   H_hr = HR(WAN_NUM,vectorL,HnumL,HcoeL) construct a TB obj
            %   with full information
            %   H_hr = HR(WAN_NUM,vectorL,HnumL,HcoeL,Type) indicate the
            %   Type of this TB obj.
            %   See also FROM_HSTRUCT, FROM_HSPARSE, FROM_HDF5, FROM_WANNIER90.
            
            % ----------- nargin ----------
            arguments
                ion double{mustBeInteger} = 4;
                QuantumNumber = int32([0 ,0 ,0]);
                options.HnumL double=[];
                options.HcoeL sym=sym([]);
                options.Type char = 'mat';
            end
            H_CEF = H_CEF@vasplib();
            % -------------check---------------
        end
    end
    %% Hamitonian method
    methods
        function H_CEF =PointChargeModel(H_CEF)
            
        end
        function H_CEF =StenvenOperater(H_CEF)
            
        end
    end
    %% reload
    methods
        function display()
            
        end
    end
    %% Value
    methods
        function [EIGENCAR,WAVECAR] = EIGENCAR_gen()
        end
    end
    %% DATABASE
    methods
        
    end
    %% Script
    methods
        function CEF_PCM()
            
        end
    end
    %%
    methods
        function obj = untitled4(inputArg1,inputArg2)
            %UNTITLED4 构造此类的实例
            %   此处显示详细说明
            obj.Property1 = inputArg1 + inputArg2;
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 此处显示有关此方法的摘要
            %   此处显示详细说明
            outputArg = obj.Property1 + inputArg;
        end
    end
end

