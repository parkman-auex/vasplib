classdef jm <Spin
    % We choose the component Jz and denote the common
    % eigenstate of the operators J^2
    % and Jz by |j, m>. 
    % We know
    % \begin{aligned}
    % \boldsymbol{J}^{2}|j, m\rangle &=j(j+1) \hbar^{2}|j, m\rangle, j=0, \frac{1}{2}, 1, \frac{3}{2}, \ldots \\
    % J_{z}|j, m\rangle &=m \hbar|j, m\rangle, m=-j,-j+1, \ldots, j-1, j
    % \end{aligned}
    properties
        j
        m
    end
    methods
        function jm_basis = jm(j,jz,coe,options)
            %JM jm_basis = jm(3/2,1/2)
            %   此处提供详细说明
            arguments
                j = 3/2;
                jz = [];
                coe = 1;
                options.orientation = [0 0 1];
            end
            optionsCell = namedargs2cell(options);
            if isa(j,'Spin')
                Spin_list= j;
                j = [Spin_list.J];jz = [Spin_list.Jz];coe = [Spin_list.coe];
            end
            if isempty(jz)
                jz = j:-1:-j;
            else
                %jz = jz;
            end
            jm_basis = jm_basis@Spin(j,jz,coe,optionsCell{:});
        end
    end
    methods %get
        function j = get.j(lm_basis)
            j = lm_basis.J;
        end
        function m = get.m(lm_basis)
            m = lm_basis.Jz;
        end
    end
end