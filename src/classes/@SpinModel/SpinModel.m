classdef SpinModel
    %% SpinModel - a toolbox of construct atomistic spin model
    % 2022/06/
    % 
    % At present, it focus on interface with Vampire (https://vampire.york.ac.uk/)
    % 
    %% public properties
    properties
        Rm (3,3) double
        sites
        mag_atom_num int8
        orbL (:,3) double
        % vectorL
        
        exchange_value double % The values of each exchange terms
        exchange_type
    end
    
    properties (Dependent = true)
        exchange_label sym % The individual labels of each exchange terms
        exchange_num int8 % how many exchange terms are used
    end
    %% hidden properties
    properties (Hidden = true)
        HR_obj_base HR % a base HR with only site informations
        HR_obj HR % a complete HR with exchange terms
    end
    %%  
    %% construction method
    methods
        function H_spin = SpinModel(POSCAR_name, options)        
            arguments
                POSCAR_name string = 'POSCAR'               
                options.exchange_type ...
                    {mustBeMember(options.exchange_type,...
                    {'isotropic','vectorial','tensorial'})} = 'isotropic'
            end            
            %% pre check
            [Rm, sites, ~, Atom_num] = POSCAR_read(POSCAR_name); 
            if length(Atom_num) ~= 1
                error('At present, you should provide a POSCAR file with only magnetic sites');
            end   
            H_spin.mag_atom_num = Atom_num;
            %% save properties
            H_spin.Rm = Rm;
            H_spin.sites = sites;                 
            H_spin.orbL = [[sites.rc1].',[sites.rc2].',[sites.rc3].'];
            % H_spin.vectorL = [0 0 0];   % init
            
            H_spin.exchange_type = options.exchange_type;              
            %% take innar HR objs to hidden unused informations
            HR_obj_base = HR(Atom_num);
            % HR_obj_base = HR_obj_base.input_Rm();
            HR_obj_base = HR_obj_base.input_orb_struct();
            HR_obj_base.quantumL = zeros(Atom_num, 4);
            
            H_spin.HR_obj      = HR_obj_base; % init
            H_spin.HR_obj_base = HR_obj_base;
        end        
    end
    %% dynamic properties
    methods
        function value = get.exchange_num(H_spin)
            value = length(H_spin.HR_obj.symvar_list);
        end
        
        function labels = get.exchange_label(H_spin)
            labels = H_spin.HR_obj.symvar_list;
        end
    end
    %% add exchange terms
    methods
        H_spin = add_nn(H_spin, exchange_num, exchange_label)
            
        H_spin = add_interlayer(H_spin, exchange_num, exchange_label, directions)
        
    end
    %% check the exchange strength
    methods
        check_exchange_value(H_spin)               
    end
    %% interface to Vampire
    methods
        Gen_UCF(H_spin, filename, opts)
    end
end




