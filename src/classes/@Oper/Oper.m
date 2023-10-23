classdef Oper < group
    %% An operation of a continus or discontinus group.
    % a subclass of group class
    % * Label: sym group operation
    %
    %% Description of the class:
    %     R : sym or double Rotation
    %         Real space rotation action of the operator. Square matrix with size
    %         of the number of spatial dimensions.
    %     conjugate : boolean (default False)
    %         Whether the operation includes conplex conjugation (antiunitary operator)
    %     antisymmetry : boolean (default False)
    %         Whether the operator flips the sign of the Hamiltonian (antisymmetry)
    %
    %     U : double (optional)
    %         The unitary action on the Hilbert space.
    %         May be None, to be able to treat symmetry candidates
    %     _strict_eq : boolean (default False) ?
    %         Whether to test the equality of the unitary parts when comparing with
    %         other PointGroupElements. By default the unitary parts are ignored.
    %         If True, PointGroupElements are considered equal, if the unitary parts
    %         are proportional, an overall phase difference is still allowed.
    %% Warning:
    %
    %
    %% example:
    %
    %% Note:
    % This class is partially inspired by the python code 'Qsymm'
    % <https://github.com/quantum-tinkerer/qsymm>
    % Dániel Varjas, Tómas Ö Rosdahl, and Anton R Akhmerov
    % Qsymm: algorithmic symmetry finding and symmetric Hamiltonian generation
    % New J. Phys. 20 093026 (2018)
    %
    %
    %     Notes
    %     -----
    %     As U is floating point and has a phase ambiguity at least, hence
    %     it is ignored when comparing objects by default.
    %
    %     R is the real space rotation acion. Do not include minus sign for
    %     the k-space action of antiunitary operators, such as time reversal.
    %     This minus sign will be included automatically if 'conjugate=True'.
    %
    %     For most uses R can be provided as a floating point mat. ?
    %     It is necessary to use exact sym matrix representation if the PGE has
    %     to act on Models with complicated momentum dependence (not polynomial),
    %     as the function parts of models are compared exactly. If the momentum
    %     dependence is periodic (sine, cosine and exponential), use BlochModel,
    %     this works with floating point rotations.
    %     """
    %% Change log
    %
    % * Document Date: 2021/06/12
    % * Creation Date: 2021/06/12
    % * Last updated : 2022/01/02
    %
    %% Copyright
    %
    % * parkman
    % * <parkman@buaa.edu.cn>
    %
    %% public properties    
    properties
        R =[]; % rotation matrix
        % U =[]; % realspace action
        t = [0 0 0]; % translation vevtor
        %
        conjugate = false;
        antisymmetry = false;
        %BasisHandle = Basis([]);
        %order ;
        Rf = []; % attach a Rm matrix;
        tf = [0 0 0]; % translation vevtor attach a Rm matrix;
    end
    %% private properties
    properties (GetAccess = protected,Hidden = true)
        continuous = false;
        strict_eq = false;
    end
    properties (GetAccess = protected)
        
    end
    %% construction method 
    methods
        function SymOper = Oper(R,U,t,options)
            arguments
                % (,)size class {Functions} = default
                R  =[1 0 0;0 1 0;0 0 1];
                U  = nan;
                t  {mustBeVector} = zeros(1,length(R));
                options.R = [];
                options.U = [];
                options.t = [];
                options.Rlocal = [];
                options.conjugate logical = false;
                options.antisymmetry logical = false;
                options.strict_eq logical = false;
                %options.BasisHandle Basis = Basis([]);
            end
            %
            if isempty(options.U)
                SymOper_U = U; % must be integer?
            else
                SymOper_U = options.U; % must be integer?
            end
            SymOper = SymOper@group(SymOper_U);
            if isempty(options.R)
                SymOper.R = R; % must be integer?
            else
                SymOper.R = options.R; % must be integer?
            end
            SymOper.R = integer(SymOper.R);
            if isempty(options.t)
                SymOper.t = t; % must be integer?
            else
                SymOper.t = options.t; % must be integer?
            end
            SymOper.t = integer(SymOper.t);
            
            SymOper.Rf = options.Rlocal;
            SymOper.conjugate = options.conjugate;
            SymOper.antisymmetry = options.antisymmetry;
            SymOper.strict_eq = options.strict_eq;
            %SymOper.BasisHandle = options.BasisHandle ;
        end
        %
    end
    %% construction method 2
    methods(Static)
        function SymOper = identity(dim, shape,propArgs)
            %     Return identity operator with appropriate shape.
            %
            %     Parameters
            %     ----------
            %     dim : int
            %         Dimension of real space.
            %     shape : int (optional)
            %         Size of the unitary part of the operator.
            %         If not provided, U is set to None.
            %
            %     Returns
            %     -------
            %     id :SymOper
            arguments 
                dim  double {mustBeInteger} =3;
                shape = nan;
                propArgs.?Oper; % use public 
            end
            propertyCell = namedargs2cell(propArgs);
            if ~isnan(shape) 
                U_ = eye(shape);
            else
                U_ = nan;
            end
            SymOper = Oper(eye(dim),U_,propertyCell{:});
            SymOper.e = true;
        end
        function SymOper = time_reversal(realspace_dim, U, spin,propArgs)
            %     Return a time-reversal symmetry operator
            %
            %     parameters
            %     ----------
            %     realspace_dim : int
            %         Realspace dimension
            %     U: mat (optional)
            %         The unitary action on the Hilbert space.
            %         May be None, to be able to treat symmetry candidates.
            %     spin : float or sequence of arrays (optional)
            %         Spin representation to use for the unitary action of the time reversal
            %         operator. If float is provided, it should be integer or half-integer
            %         specifying the spin representation in the standard basis, see `spin_matrices`.
            %         Otherwise a sequence of 3 arrays of identical square size must be provided
            %         representing 3 components of the angular momentum operator. The unitary action
            %         of time-reversal operator is `U = exp(-i π s_y)`. Only one of `U` and `spin`
            %         may be provided.
            %
            %     Returns
            %     -------
            %     T : PointGroupElement
            arguments
                realspace_dim  double {mustBeInteger} =3;
                U double =nan;
                spin double = nan;
                propArgs.?Oper; % use public 
                propArgs.conjugate = true;
            end
            propertyCell = namedargs2cell(propArgs);
            if ~all(all(isnan(U))) && ~all(all(all(isnan(spin))))
                raise ValueError('Only one of `U` and `spin` may be provided.');
            end
            if ~isnan(spin)
                U_tr = Oper.spin_rotation(pi*[0, 1, 0], spin);
            else
                U_tr = U;
            end
            R_tr = eye(realspace_dim);
            SymOper = Oper(R_tr,U_tr,propertyCell{:});
        end
        function SymOper = particle_hole(realspace_dim, U,propArgs)
            %     Return a particle-hole symmetry operator
            %
            %     parameters
            %     ----------
            %     realspace_dim : int
            %         Realspace dimension
            %     U: ndarray (optional)
            %         The unitary action on the Hilbert space.
            %         May be None, to be able to treat symmetry candidates
            %
            %     Returns
            %     -------
            %     P : PointGroupElement
            %    
            arguments
                realspace_dim  double {mustBeInteger} =3;
                U double =nan;
                propArgs.?Oper; % use public 
                propArgs.conjugate = true;
                propArgs.antisymmetry = true;
            end
            propertyCell = namedargs2cell(propArgs);
            U_ph = U;
            R_ph = eye(realspace_dim);
            SymOper = Oper(R_ph,U_ph,propertyCell{:});
        end
        function SymOper = chiral(realspace_dim, U,propArgs)
            %     Return a chiral symmetry operator
            %
            %     parameters
            %     ----------
            %     realspace_dim : int
            %         Realspace dimension
            %     U: ndarray (optional)
            %         The unitary action on the Hilbert space.
            %         May be None, to be able to treat symmetry candidates
            %
            %     Returns
            %     -------
            %     P : PointGroupElement
            %     
            arguments
                realspace_dim  double {mustBeInteger} =3;
                U double =nan;
                propArgs.?Oper; % use public
                propArgs.antisymmetry = true;
            end
            propertyCell = namedargs2cell(propArgs);
            U_c = U;
            R_c = eye(realspace_dim);
            SymOper = Oper(R_c,U_c,propertyCell{:});
        end
        function SymOper = inversion(realspace_dim, U,quantumL,propArgs)
            %     Return an inversion operator
            %
            %     parameters
            %     ----------
            %     realspace_dim : int
            %         Realspace dimension
            %     U: ndarray (optional)
            %         The unitary action on the Hilbert space.
            %         May be None, to be able to treat symmetry candidates
            %
            %     Returns
            %     -------
            %     P : PointGroupElement
            %     
            arguments
                realspace_dim  double {mustBeInteger} =3;
                U double =nan;
                quantumL = nan;
                propArgs.?Oper; % use public
            end
            propertyCell = namedargs2cell(propArgs);
            U_inv = U;
            if ~isnan(quantumL)
                %
            end
            R_inv = -eye(realspace_dim);
            SymOper = Oper(R_inv,U_inv,propertyCell{:});
        end
        function SymOper = rotation(angle, axis, inversion, U, spin,options,propArgs)
            %     Return a rotation operator
            %
            %     parameters
            %     ----------
            %     angle : float
            %         Rotation angle in units of 2 pi.
            %     axis : ndarray or None (default)
            %         Rotation axis, optional. If not provided, a 2D rotation is generated
            %         around the axis normal to the plane. If a 3D vector is provided,
            %         a 3D rotation is generated around this axis. Does not need to be
            %         normalized to 1.
            %     inversion : bool (default False)
            %         Whether to generate a rotoinversion. By default a proper rotation
            %         is returned. Only valid in 3D.
            %     U: ndarray (optional)
            %         The unitary action on the Hilbert space.
            %         May be None, to be able to treat symmetry candidates
            %     spin : float or sequence of arrays (optional)
            %         Spin representation to use for the unitary action of the
            %         operator. If float is provided, it should be integer or half-integer
            %         specifying the spin representation in the standard basis, see `spin_matrices`.
            %         Otherwise a sequence of 3 arrays of identical square size must be provided
            %         representing 3 components of the angular momentum operator. The unitary action
            %         of rotation operator is `U = exp(-i n⋅s)`. In 2D the z axis is assumed to be
            %         the rotation axis. Only one of `U` and `spin` may be provided.
            %
            %     Returns
            %     -------
            %     P : PointGroupElement
            % 
            arguments
                angle  ;
                axis  = nan;
                inversion logical = false;
                U  = nan;
                spin  = nan;
                options.sym = false;
                options.rightorleft {mustBeMember(options.rightorleft,{'left','right'})}= 'right';
                propArgs.?Oper; % use public
            end
            propertyCell = namedargs2cell(propArgs);
            if isa(U,'double')
            if ~all(all(isnan(U))) && ~all(all(all(isnan(spin))))
                raise ValueError('Only one of `U` and `spin` may be provided.');
            end
            elseif isa(U,'pauli_matric')
                U = double(U);
            end
            if options.sym
                if isa(angle,'sym')
                else
                    angle = sym(2 * pi * angle);
                end
                U_PG = sym(U);
            else
                if isa(angle,'sym')
                else
                    angle = 2 * pi * angle;
                end
                U_PG = U;
            end
            if strcmp(options.rightorleft,'left')
                rightorleft = -1;
            else
                rightorleft = 1;
            end
            if isnan(axis)
                %2D
                R_PG = ([[cos(angle), sin(angle)];
                    [-sin(angle), cos(angle)]]);
                if ~isnan(spin)
                    U_PG = Oper.spin_rotation(angle * ([0, 0, 1]), spin);
                end
            elseif length(axis) == 3
                % 3D
                n = angle * axis / norm(axis);
                R_PG = Oper.spin_rotation(n, rightorleft*Oper.L_matrices(3, 1));
                if inversion
                    R_PG = -R_PG;
                else

                end

                if ~isnan(spin)
                    if ~isvector(spin)
                        U_PG = Oper.spin_rotation(n, spin);
                    else

                    end
                end
            else
                raise ValueError('`axis` needs to be `None` or a 3D vector.')
            end
            SymOper = Oper(R_PG,U_PG,propertyCell{:});          
        end
        function SymOper = spaceRotation(angle, axis,t, inversion, U, spin,options,propArgs)
            %     Return a space rotation operator
            %
            %     parameters
            %     ----------
            %     angle : float
            %         Rotation angle in units of 2 pi.
            %     axis : ndarray or None (default)
            %         Rotation axis, optional. If not provided, a 2D rotation is generated
            %         around the axis normal to the plane. If a 3D vector is provided,
            %         a 3D rotation is generated around this axis. Does not need to be
            %         normalized to 1.
            %     t:
            %
            %     inversion : bool (default False)
            %         Whether to generate a rotoinversion. By default a proper rotation
            %         is returned. Only valid in 3D.
            %     U: ndarray (optional)
            %         The unitary action on the Hilbert space.
            %         May be None, to be able to treat symmetry candidates
            %     spin : float or sequence of arrays (optional)
            %         Spin representation to use for the unitary action of the
            %         operator. If float is provided, it should be integer or half-integer
            %         specifying the spin representation in the standard basis, see `spin_matrices`.
            %         Otherwise a sequence of 3 arrays of identical square size must be provided
            %         representing 3 components of the angular momentum operator. The unitary action
            %         of rotation operator is `U = exp(-i n⋅s)`. In 2D the z axis is assumed to be
            %         the rotation axis. Only one of `U` and `spin` may be provided.
            %
            %     Returns
            %     -------
            %     P : Space Oper
            % 
            arguments
                angle double ;
                axis double = nan;
                t double =[0,0,0];
                inversion logical = false;
                U double =nan;
                spin double = nan;
                options.sym = false;
                propArgs.?Oper; % use public
            end
            propertyCell = namedargs2cell(propArgs);
            U_SG = U;
            t_SG = t;
            if ~all(all(isnan(U))) && ~all(all(all(isnan(spin))))
                raise ValueError('Only one of `U` and `spin` may be provided.');
            end
            angle = 2 * pi * angle;
            if isnan(axis)
                %2D
                R_SG = ([[cos(angle), -sin(angle)];
                    [sin(angle), cos(angle)]]);
                if ~isnan(spin)
                    U_SG = Oper.spin_rotation(angle * ([0, 0, 1]), spin);
                end
            elseif len(axis) == 3
                % 3D
                n = angle * axis / norm(axis);
                R_SG = Oper.spin_rotation(n, L_matrices(d=3, l=1));
                if inversion
                    R_SG = -R_SG;
                else
                end
                
                if ~isnan(spin)
                    U_SG = spin_rotation(n, spin);
                end
            else
                raise ValueError('`axis` needs to be `None` or a 3D vector.')
            end 
            SymOper = Oper(R_SG,U_SG,t_SG,propertyCell{:});
        end
        function SymOper = C3z(realspace_dim, inversion, U, spin)
        end
        function SymOper = C4z(realspace_dim, inversion, U, spin)
        end
        function SymOper = C6z(realspace_dim, inversion, U, spin)
        end
        function SymOper = mirror(axis, U, spin,options,propArgs)
            %     Return a mirror operator
            %
            %     Parameters
            %     ----------
            %     axis : ndarray
            %         Normal of the mirror. The dimensionality of the operator is the same
            %         as the length of `axis`.
            %     U: ndarray (optional)
            %         The unitary action on the Hilbert space.
            %         May be None, to be able to treat symmetry candidates
            %     spin : float or sequence of arrays (optional)
            %         Spin representation to use for the unitary action of the
            %         operator. If float is provided, it should be integer or half-integer
            %         specifying the spin representation in the standard basis, see `spin_matrices`.
            %         Otherwise a sequence of 3 arrays of identical square size must be provided
            %         representing 3 components of the angular momentum operator. The unitary action
            %         of mirror operator is `U = exp(-i π n⋅s)` where n is normalized to 1. In 2D the
            %         axis is treated as x and y coordinates. Only one of `U` and `spin` may be provided.
            %
            %     Returns
            %     -------
            %     P : PointGroupElement
            %
            %     Notes:
            %     ------
            %         Warning: in 2D the real space action of a mirror and and a 2-fold rotation
            %         around an axis in the plane is identical, however the action on angular momentum
            %         is different. Here we consider the action of the mirror, which is the same as the
            %         action of a 2-fold rotation around the mirror axis.
            % 
            arguments
                axis  = nan;
                U  =nan;
                spin double = nan;
                options.sym = false;
                propArgs.?Oper; % use public
            end
            propertyCell = namedargs2cell(propArgs);
            U_mirror = U;
            if ~all(all(isnan(U))) && ~all(all(all(isnan(spin))))
                raise ValueError('Only one of `U` and `spin` may be provided.');
            end
            
            axis = axis/norm(axis);
            axis = reshape(axis,length(axis),1);
            R_mirror = eye(length(axis)) - 2 *(axis*axis.');
            if ~isnan(spin)
                if length(axis) == 2
                    axis(3) = 0;
                end
                U_mirror = Oper.spin_rotation(pi * axis, spin);
            end
            if isa(R_mirror,'sym')
                R_mirror = simplify(R_mirror);
            end
            if isa(U_mirror,'sym')
                U_mirror = simplify(U_mirror);
            end
            SymOper = Oper(R_mirror,U_mirror,propertyCell{:});
        end
        function Mx
        end
        function My
        end
        function Mz
        end
        %  Predefined point groups
        function group = square(tr, ph, generators, spin,options,propArgs)
            %     
            %     Generate square point group in standard basis.
            %
            %     Parameters
            %     ----------
            %     tr, ph : bool (default True)
            %         Whether to include time-reversal and particle-hole
            %         symmetry.
            %     generators : bool (default false)
            %         Only return the group generators if True.
            %     spin : float or sequence of arrays (optional)
            %         Spin representation to use for the unitary action of the
            %         operator. If not provided, the PointGroupElements have the unitary
            %         action set to None. If float is provided, it should be integer or half-integer
            %         specifying the spin representation in the standard basis, see `spin_matrices`.
            %         Otherwise a sequence of 3 arrays of identical square size must be provided
            %         representing 3 components of the angular momentum operator. The unitary action
            %         of rotation operator is `U = exp(-i n⋅s)`. In 2D the z axis is assumed to be
            %         the rotation axis. If `ph` is True, `spin` may not be provided, as it is not
            %         possible to deduce the unitary representation of particle-hole symmetry from
            %         spin alone. In this case construct the particle-hole operator manually.
            %
            %     Returns
            %     -------
            %     set of PointGroupElement objects with integer rotations
            %
            %     Notes:
            %     ------
            %         Warning: in 2D the real space action of a mirror and and a 2-fold rotation
            %         around an axis in the plane is identical, however the action on angular  momentum
            %         is different. Here we consider the action of the mirror, which is the same as the
            %         action of a 2-fold rotation around the mirror axis, assuming inversion acts trivially.
            %
            arguments
                tr logical =true;
                ph logical =true;
                generators logical = false;
                spin double = nan;
                options.sym = false;
                propArgs.?Oper; % use public
            end
            optionsCell = namedargs2cell(options);
            propertyCell = namedargs2cell(propArgs);
            if ~all(all(isnan(ph))) && ~all(all(all(isnan(spin))))
                raise ValueError('If `ph` is True, `spin` may not be provided, as it is not '...
                         'possible to deduce the unitary representation of particle-hole symmetry '...
                         'from spin alone. In this case construct the particle-hole operator manually.');
            end
            Mx = Oper.mirror([1, 0],nan,spin,optionsCell{:},propertyCell{:});
            C4 = Oper.rotation(1/4, nan,false,nan, spin,optionsCell{:},propertyCell{:});
            gens = [Mx, C4];
            if tr
                propArgs.antisymmetry = true;
                propertyCell = namedargs2cell(propArgs);
                TR = Oper.time_reversal(2, spin,propertyCell{:});
                gens = [gens,TR];
            end
            if ph
                propArgs.antisymmetry = true;
                propArgs.conjugate = true;
                PH = Oper.particle_hole(2,propertyCell{:});
                gens = [gens,PH];
            end
            if generators
                group =  gens;
            else
                %fprintf('group muplicity: %d\n',gens.order);
                group =  gens.generate_group();
            end
        end
        function group = cubic(tr, ph, generators, spin,options,propArgs)
            %     Generate cubic point group in standard basis.
            %
            %     Parameters
            %     ----------
            %     tr, ph : bool (default True)
            %         Whether to include time-reversal and particle-hole
            %         symmetry.
            %     generators : bool (default false)
            %         Only return the group generators if True.
            %     spin : float or sequence of arrays (optional)
            %         Spin representation to use for the unitary action of the
            %         operator. If not provided, the PointGroupElements have the unitary
            %         action set to None. If float is provided, it should be integer or half-integer
            %         specifying the spin representation in the standard basis, see `spin_matrices`.
            %         Otherwise a sequence of 3 arrays of identical square size must be provided
            %         representing 3 components of the angular momentum operator. The unitary action
            %         of rotation operator is `U = exp(-i n⋅s)`. If `ph` is True, `spin` may not be
            %         provided, as it is not  possible to deduce the unitary representation of
            %         particle-hole symmetry from spin alone. In this case construct the
            %         particle-hole operator manually.
            %
            %     Returns
            %     -------
            %     set of PointGroupElement objects with integer rotations
            %
            %     Notes:
            %     ------
            %         We assume inversion acts trivially in spin space.
            arguments
                tr logical =true;
                ph logical =true;
                generators logical = false;
                spin double = nan;
                options.sym = false;
                propArgs.?Oper; % use public
            end
            optionsCell = namedargs2cell(options);
            propertyCell = namedargs2cell(propArgs);
            if ~isnan(ph) && ~isnan(spin)
                raise ValueError('If `ph` is True, `spin` may not be provided, as it is not '...
                    'possible to deduce the unitary representation of particle-hole symmetry '...
                    'from spin alone. In this case construct the particle-hole operator manually.');
            end
            if isnan(spin)
                I = Oper.inversion(3,propertyCell{:});
            else
                U = Oper.spin_rotation(zeros(3), spin);
                I = Oper.inversion(3, U,propertyCell{:});
            end
            C4 = Oper.rotation(1/4, [1, 0, 0],false,nan,spin,optionsCell{:},propertyCell{:});
            C3 = Oper.rotation(1/3, [1, 1, 1],false,nan,spin,optionsCell{:},propertyCell{:});
            gens = [I, C4, C3];
            if tr
                propArgs.antisymmetry = true;
                propertyCell = namedargs2cell(propArgs);
                TR = Oper.time_reversal(3, spin,propertyCell{:});
                gens = [gens,TR];
            end
            if ph
                propArgs.antisymmetry = true;
                propArgs.conjugate = true;
                propertyCell = namedargs2cell(propArgs);
                PH = Oper.particle_hole(3,propertyCell{:});
                gens = [gens,PH];
            end
            if generators
                group =  gens;
            else
                %fprintf('group muplicity: %d\n',gens.order);
                group =  gens.generate_group();
            end
        end
        function group = hexagonal_2D(tr, ph, generators, spin,options, propArgs)
            %     Generate hexagonal point group in standard basis in 2 dimensions.
            %     Mirror symmetries with the main coordinate axes as normals are included.
            %
            %     Parameters
            %     ----------
            %     tr, ph : bool (default True)
            %         Whether to include time-reversal and particle-hole
            %         symmetry.
            %     generators : bool (default True)
            %         Only return the group generators if True.
            %     spin : float or sequence of arrays (optional)
            %         Spin representation to use for the unitary action of the
            %         operator. If not provided, the PointGroupElements have the unitary
            %         action set to None. If float is provided, it should be integer or half-integer
            %         specifying the spin representation in the standard basis, see `spin_matrices`.
            %         Otherwise a sequence of 3 arrays of identical square size must be provided
            %         representing 3 components of the angular momentum operator. The unitary action
            %         of rotation operator is `U = exp(-i n⋅s)`. In 2D the z axis is assumed to be
            %         the rotation axis. If `ph` is True, `spin` may not be provided, as it is not
            %         possible to deduce the unitary representation of particle-hole symmetry from
            %         spin alone. In this case construct the particle-hole operator manually.
            %
            %     Returns
            %     -------
            %     Oper List
            %
            %     Notes:
            %     ------
            %         Warning: in 2D the real space action of a mirror and and a 2-fold rotation
            %         around an axis in the plane is identical, however the action on angular momentum
            %         is different. Here we consider the action of the mirror, which is the same as the
            %         action of a 2-fold rotation around the mirror axis, assuming inversion acts trivially.
            arguments
                tr logical =true;
                ph logical =true;
                generators logical = false;
                spin double = nan;
                options.sym = false;
                propArgs.?Oper; % use public
            end
            optionsCell = namedargs2cell(options);
            propertyCell = namedargs2cell(propArgs);
            if ~isnan(ph) && ~isnan(spin)
                raise ValueError('If `ph` is True, `spin` may not be provided, as it is not '...
                    'possible to deduce the unitary representation of particle-hole symmetry '...
                    'from spin alone. In this case construct the particle-hole operator manually.');
            end
%             if ~isnan(spin) 
%                 U6 = spin_rotation(pi / 3 * ([0, 0, 1]), spin);
%             else
%                 U6 = nan;
%             end
            Mx = Oper.mirror([1, 0], spin,optionsCell{:},propertyCell{:});
            C6 = Oper.rotation(1/6, nan,false,nan,spin,optionsCell{:},propertyCell{:});
            gens = [Mx, C6];
            if tr
                propArgs.antisymmetry = true;
                propertyCell = namedargs2cell(propArgs);
                TR = Oper.time_reversal(2, spin,propertyCell{:});
                gens = [gens,TR];
            end
            if ph
                propArgs.antisymmetry = true;
                propArgs.conjugate = true;
                PH = Oper.particle_hole(2,propertyCell{:});
                gens = [gens,PH];
            end
            if generators
                group =  gens;
            else
                %fprintf('group muplicity: %d\n',gens.order);
                group =  gens.generate_group();
            end
        end
        function group = hexagonal(tr, ph, generators, spin,options,propArgs)
            %     Generate hexagonal point group in standard basis in 3 dimensions.
            %     Mirror symmetries with the main coordinate axes as normals are included.
            %
            %     Parameters
            %     ----------
            %     tr, ph : bool (default True)
            %         Whether to include time-reversal and particle-hole
            %         symmetry.
            %     generators : bool (default True)
            %         Only return the group generators if True.
            %     spin : float or sequence of arrays (optional)
            %         Spin representation to use for the unitary action of the
            %         operator. If not provided, the PointGroupElements have the unitary
            %         action set to None. If float is provided, it should be integer or half-integer
            %         specifying the spin representation in the standard basis, see `spin_matrices`.
            %         Otherwise a sequence of 3 arrays of identical square size must be provided
            %         representing 3 components of the angular momentum operator. The unitary action
            %         of rotation operator is `U = exp(-i n⋅s)`. In 2D the z axis is assumed to be
            %         the rotation axis. If `ph` is True, `spin` may not be provided, as it is not
            %         possible to deduce the unitary representation of particle-hole symmetry from
            %         spin alone. In this case construct the particle-hole operator manually.
            %
            %     Returns
            %     -------
            %     Oper List
            %
            arguments
                tr logical =true;
                ph logical =true;
                generators logical = false;
                spin double = nan;
                options.sym = false;
                propArgs.?Oper; % use public
            end
            optionsCell = namedargs2cell(options);
            propertyCell = namedargs2cell(propArgs);
            if ~isnan(ph) && ~isnan(spin)
                raise ValueError('If `ph` is True, `spin` may not be provided, as it is not '...
                    'possible to deduce the unitary representation of particle-hole symmetry '...
                    'from spin alone. In this case construct the particle-hole operator manually.');
            end
            if isnan(spin)
                I = Oper.inversion(3,propertyCell{:});
            else
                U = Oper.spin_rotation(zeros(3), spin);
                I = Oper.inversion(3, U,propertyCell{:});
            end
            C2x = Oper.rotation(1/2, [1, 0, 0],false,nan,spin,optionsCell{:},propertyCell{:});
            C6 = Oper.rotation(1/6, [0, 0, 1],false,nan,spin,optionsCell{:},propertyCell{:});
            gens = [I ,C2x,C6];
            if tr
                propArgs.antisymmetry = true;
                propertyCell = namedargs2cell(propArgs);
                TR = Oper.time_reversal(3, spin,propertyCell{:});
                gens = [gens,TR];
            end
            if ph
                propArgs.antisymmetry = true;
                propArgs.conjugate = true;
                propertyCell = namedargs2cell(propArgs);
                PH = Oper.particle_hole(3,propertyCell{:});
                gens = [gens,PH];
            end
            if generators
                group =  gens;
            else
                %fprintf('group muplicity: %d\n',gens.order);
                group =  gens.generate_group();
            end
        end
    end
    %% Math operation
    methods 
        function [TrueOrNot,result] = commute(SymOper1,SymOper2)
            % commute is not right
            U1 = SymOper1.U;
            U2 = SymOper2.U;
            result = U1*U2-U2*U1;
            if isa(U1,'sym')||isa(U2,'sym')
                zeromat = sym(zeros(size(U1)));
            else
                zeromat = zeros(size(U1));
            end
            if isequal(result,zeromat)
                TrueOrNot = true;
            else
                TrueOrNot = false;
            end
        end
    end
    methods
        function SymOper = Ugen(SymOper,Basis,options)
            arguments
                SymOper Oper
                Basis
                options.Rm = POSCAR_read;
                options.sym = false;
                options.center = [0,0,0];
            end
            optionsCell = namedargs2cell(options);
            for i =1:numel(SymOper)
                if isempty(SymOper(i).Rf)
                    SymOper(i) = SymOper(i).attachRm(options.Rm);
                end
                if isnan(SymOper(i).U)
                    SymOper(i).U =Basis.rotation('Oper', SymOper(i),optionsCell{:});
                end
            end
        end
    end
    %% Character
    methods
        function [SYMCAR,OperObj] = Character_gen(OperObj,vasplibobj,klist,options)
            arguments
                OperObj Oper;
                vasplibobj vasplib;
                klist = [0,0,0;0.5,0,0;0,0.5,0;0.5,0.5,0;0,0,0.5;0.5,0,0.5;0,0.5,0.5;0.5,0.5,0.5;];
                options.BasisFunction = [];
                options.generate_group = false;
                options.classify = true;
                options.Accuracy = 1e-6;
                options.sum = false;
            end
            % alert: fractional klist!!
            % LG checking!
            Accuracy = round(log10(options.Accuracy));
            if isempty(options.BasisFunction)
                BasisFunction = BasisFunc(vasplibobj);
            else
                BasisFunction = options.BasisFunction;
            end
            % full group
            if options.generate_group
                OperObj = OperObj.generate_group();
            end
            % classify
            if options.classify
                % to be added
                % OperObj = OperObj.class();
            end
            nOper = length(OperObj);
            % refresh
            for i = 1:nOper
                if isnan(OperObj(i).U)
                    OperObj(i) = OperObj(i).attachRm(vasplibobj.Rm);
                    try
                        OperObj(i).U = BasisFunction.rotation('Oper',OperObj(i));
                    catch
                        error('fail to generate Oper matric!');
                    end
                end
            end
            kn = size(klist,1);
            switch class(vasplibobj)
                case 'HR'
                    Basis_num = vasplibobj.WAN_NUM;
                    [EIGENCAR,WAVECAR] = vasplibobj.EIGENCAR_gen('klist',klist,'convention','I','printmode',false);
                case {'Htrig','HK'}
                    Basis_num = vasplibobj.Basis_num;
                    [EIGENCAR,WAVECAR] = vasplibobj.EIGENCAR_gen('klist',klist*vasplibobj.Gk,'printmode',false);
            end
            % check convention here
            VCAR = zeros(kn,Basis_num,nOper);
            FactorL = zeros(kn,nOper);
            orbL = vasplibobj.orbL;
            for i = 1:nOper
                VCAR(:,:,i) = exp(1i*2*pi*(klist*OperObj(i).Rf-klist)*(orbL.'));
                FactorL(:,i) = exp(-1i*2*pi*(klist*OperObj(i).Rf)*OperObj(i).tf.');
            end
            SYMCAR = zeros(Basis_num,nOper,kn);
            for i =1:kn
                for j = 1:nOper
                    SYMCAR(:,j,i) = FactorL(i,j)*Oper.EIGEN_Kone(WAVECAR(:,:,i),OperObj(j).U,diag(VCAR(i,:,j)));
                end
            end
            SYMCAR = roundn(real(SYMCAR),Accuracy) + 1i*roundn(imag(SYMCAR),Accuracy);
            % sum 
            if options.sum
                % waiting
            end
        end
        function EIGENCAR_SYM = EIGEN(OperObj,WAVECAR,klist,orbL,options)
            arguments
                OperObj Oper;
                WAVECAR;
                klist;
                orbL;
                options.Accuracy = 1e-6;
                options.Obj = 'HR';
            end
            %
            switch options.Obj
                case 'HR'
                    % check convention here
                    VCAR = exp(1i*2*pi*(klist*OperObj.R-klist)*orbL.');
                    FactorL = exp(-1i*2*pi*klist*OperObj.R*OperObj.tf.');
                case {'Htrig','HK'}
            end
            %
            Accuracy = round(log10(options.Accuracy));
            %sizeWAVECAR = size(WAVECAR);
            kn = size(WAVECAR,3);
            NBANDS = size(WAVECAR,2);
            %Norb = size(WAVECAR,1);
            EIGENCAR_SYM = zeros(NBANDS,kn);
            for i =1:kn
                EIGENCAR_SYM(:,i) = FactorL(i)*Oper.EIGEN_Kone(WAVECAR(:,:,i),OperObj.U,diag(VCAR(i,:)));
            end
            EIGENCAR_SYM = roundn(real(EIGENCAR_SYM),Accuracy) + 1i*roundn(imag(EIGENCAR_SYM),Accuracy);
        end
    end
    methods(Static)
        function EIGENCAR_Kone = EIGEN_Kone(WAVECAR_one,D,V)
            %$\left\langle\psi_{n \mathbf{k}}\left|\mathcal{O}_{s}\right| \psi_{n \mathbf{k}}\right\rangle=e^{-i\left(R_{s} \mathbf{k} \cdot \mathbf{v}_{s}\right)}\left[\overline{C^{\dagger} V\left(R_{s} \mathbf{k}-\mathbf{k}\right) D C}\right]_{n n}$
            %with $\bar{V}(\mathbf{k})_{\mu^{\prime} \beta, \mu \alpha}=e^{i \mathbf{k} \cdot \tau_{\mu}} \delta_{\mu \mu^{\prime}} \delta_{\alpha \beta}, \bar{C}_{\mu \alpha, n}=C_{\mu \alpha}^{n}, \bar{D}_{\mu^{\prime} \beta, \mu \alpha}=\left\{\begin{array}{cl}D_{\beta \alpha}^{s, \mu} & \text { when } \mathbf{v}_{s}+R_{s} \tau_{\mu}=\mathbf{L}_{0}^{i}+\tau_{\mu^{\prime}} \\ 0 & \text { otherwise. }\end{array}\right.$
            % 10.1016/j.cpc.2020.107760
            %
            EIGENCAR_Kone = diag(WAVECAR_one'*V*D*WAVECAR_one);
            %
        end
        function EIGENCAR_Kone = EIGEN_Kone2(WAVECAR_one,U)
            NBANDS = size(WAVECAR_one,2);
            EIGENCAR_Kone = zeros(NBANDS,1);
            for i = 1:NBANDS
                EIGENCAR_Kone(i) = Oper.EIGEN_one(WAVECAR_one(:,i),U);
            end
        end
        function EIGEN = EIGEN_one(WAVEFUNC,U)
            EIGEN = (WAVEFUNC)'*(U*WAVEFUNC);
            %EIGEN = mean(U*WAVEFUNC ./WAVEFUNC);
        end
    end
    %% 
    methods(Static)
        function U = Umat(SymOper1,SymOper2,options)
            arguments
                SymOper1
                SymOper2
                options.disp = true;
            end
            if isa(SymOper1,'Oper')
                A =SymOper1.U;
            else
                A =SymOper1;
            end
            if isa(SymOper2,'Oper')
                B =SymOper2.U;
            else
                B =SymOper2;
            end
            % check eig
            [QA,DA]=eig(A);
            [QB,DB]=eig(B);
            k = DA/DB;
            U = sym((QB/QA)/sqrtm(k));
            if options.disp
                disp(('U * A * U^* = B'));
                sym(A)
                sym(B)
                sym(U)
            end
        end
        function P = similarity(SymOper1,SymOper2,options)
            arguments
                SymOper1
                SymOper2
                options.disp = true;
            end
            if isa(SymOper1,'Oper')
                A =SymOper1.U;
            else
                A =SymOper1;
            end
            if isa(SymOper2,'Oper')
                B =SymOper2.U;
            else
                B =SymOper2;
            end
            [V,~]=jordan(A);
            [U,~]=jordan(B);
            P = (V/(U));
            if options.disp
                disp(('P * A * P^-1 = B'));
                A
                B
                P
            end
        end
    end
    %% overload method alert SG oper still not be implemented waiting
    methods
        %disp is good
        function str = disp(SymOper,options)
            arguments
                SymOper Oper;
                options.full logical =false;
            end
            len = length(SymOper);
            if len == 1
                builtin('disp',SymOper);
                try
                    str = SymOper.pretty();
                    fprintf(str);
                catch
                end
                return;
            end
            if options.full
                for i =1:len
                    if len > 1
                        fprintf("=============== %dth Oper ===============\n",i);
                    end
                    builtin('disp',SymOper(i));
                    str = SymOper(i).pretty();
                    fprintf(str);
                end
            else
                %builtin('disp',SymOper);
                for i =1:len
                    str = SymOper(i).pretty('full',false);
                    fprintf("%3d : %s\n",i,str);
                end
            end
        end
        % sym 
        function [SymMat,SymR] = sym(SymOper)
            if numel(SymOper)>1
                for i = 1:numel(SymOper)
                    [SymOper(i).U,SymOper(i).R] = sym(SymOper(i));
                end
                SymMat = SymOper;SymR = [];
            else
                SymMat = sym(SymOper.U);
                SymR = sym(SymOper.R);
            end

        end
        % E
        function SymOper = E(SymOper)
            %    Return identity element with the same structure as self
            dim = length(SymOper.R);
            R_E = eye(dim) ;
            %             t = [0,0,0];
            if isnan(SymOper.U)
                U_E = SymOper.U;
            else
                U_E = eye(length(SymOper.U));
            end
            SymOper = Oper(R_E, U_E,SymOper.t,'conjugate',SymOper.conjugate,...
                'antisymmetry',SymOper.antisymmetry,...
                'strict_eq',SymOper.strict_eq...
                );
        end
    end
    % reload these functions carefully
    methods
        % ==
        function [basic_eq,U_eq]=eq(SymOper1,SymOper2)
            % need repair !!
            m = length(SymOper1);
            n = length(SymOper2);
            if m == 0 && n == 0
               basic_eq =true;
               U_eq = true;
            elseif m == 0 || n == 0
                basic_eq = false;
                U_eq = false;
                % a bug will occure when void == n?
            elseif m == 1 && n ==1
                R_eq = Oper.allclose(SymOper1.R, SymOper2.R);
                basic_eq = R_eq && SymOper1.conjugate ==SymOper2.conjugate && SymOper1.antisymmetry == SymOper2.antisymmetry && isequal(sym(SymOper1.t),  sym(SymOper2.t));
                if basic_eq && (SymOper1.strict_eq || SymOper2.strict_eq)
                    if isnan(SymOper1.U) && isnan(SymOper2.U)
                        U_eq = true;
                    elseif xor(isnan(SymOper1.U) , isnan(SymOper2.U))
                        U_eq = false;
                    else
                        SymOper3 = (SymOper1.inv()*SymOper2);
                        [prop, coeff] = Oper.prop_to_id(SymOper3.U);
                        U_eq = (prop && Oper.isclose(abs(coeff), 1));
                    end
                else
                    U_eq = true;
                end
            elseif m == 1
                basic_eq(n) = false;
                U_eq(n) = false;
                for i = 1:n
                    [basic_eq(i),U_eq(i)] = eq(SymOper1,SymOper2(i));
                end
            elseif n == 1
                basic_eq(m) = false;
                U_eq(m) = false;
                for i = 1:m
                    [basic_eq(i),U_eq(i)] = eq(SymOper2,SymOper1(i));
                end
            elseif m == n
                basic_eq = false(size(SymOper1));
                U_eq = false(size(SymOper1));
                SymOper1 = sort(SymOper1);
                SymOper2 = sort(SymOper2);
                for i =1:m
                    [basic_eq(i),U_eq(i)] = eq(SymOper2(i),SymOper1(i));
                end
%                 if isempty(~U_eq) 
%                     basic_eq = true;
%                     U_eq = true;
%                     return
%                 else
%                     SymOper3 = SymOper1(~U_eq);
%                     SymOper4 = SymOper2(~U_eq);
%                     for i = 1:length(SymOper3)
%                         if ~sum(SymOper3(i) == SymOper4)
%                             basic_eq = false;
%                             U_eq = false;
%                             return;
%                         end
%                         basic_eq = true;
%                         U_eq = true;
%                     end
%                 end
 
            else
                basic_eq = false;
                U_eq = false;
                
            end
        end
        % < 
        function result = lt(SymOper1,SymOper2)
            % Sort group elements:
            % First by conjugate and a, then R = identity, then the rest
            % lexicographically
            Rs = SymOper1.R;
            Ro = SymOper2.R;
            identity = eye(length(Rs));
            if SymOper1.conjugate ~= SymOper2.conjugate || ...
                    SymOper1.antisymmetry ~= SymOper2.antisymmetry
                
                if SymOper1.conjugate < SymOper2.conjugate
                    result =true;
                elseif SymOper1.conjugate == SymOper2.conjugate
                    if SymOper1.antisymmetry < SymOper2.antisymmetry
                        result =true;
                    else
                        result =false;
                    end
                else
                    result =false;
                end
            elseif xor(Oper.allclose(Rs,identity), Oper.allclose(Ro,identity))
                result = Oper.allclose(Rs,identity);
            else
                %                 tf = issymmetric(Ro-Rs);
                %                 if tf
                %                     result = all( eig(Ro-Rs) > 0 );
                %                 else
                %                     result = false;
                %                 end
                L = Rs(:) < Ro(:);
                B = ~Oper.isclose(Rs(:),Ro(:));
                for i =1:length(B)
                    if B(i)
                       result = L(i);
                       return;
                    end
                end
                result  = false;
            end
        end
        % .*
        function SymOper_out = times(SymOper1,SymOper2)
            m = length(SymOper1);
            n = length(SymOper2);
            if m == 1&& n ==1
                SymOper_out = SymOper1;
                if all(all((isnan(SymOper1.U)))) 
                    SymOper_out.U = nan;
                elseif SymOper1.conjugate
                    SymOper_out.U = SymOper1.U*conj(SymOper2.U);
                else
                    SymOper_out.U = SymOper1.U*(SymOper2.U);
                end
                SymOper_out.R =  SymOper1.R*(SymOper2.R);
                SymOper_out.t =  SymOper2.t*SymOper1.R.'+SymOper1.t;
                SymOper_out.conjugate = xor(SymOper1.conjugate,SymOper2.conjugate);
                SymOper_out.antisymmetry = xor(SymOper1.antisymmetry,SymOper2.antisymmetry);
                SymOper_out.strict_eq = SymOper1.strict_eq || SymOper2.strict_eq;
            else 
                error('.* is elemenary operator');
            end
        end
        % inv
        function SymOper = inv(SymOper)
            % Invert SymOper
            if isnan(SymOper.U) 
                Uinv = nan;
            elseif SymOper.conjugate
                Uinv = conj(inv(SymOper.U));
            else
                Uinv = inv(SymOper.U);
            end
            % Check if inverse is stored, if not, calculate it
            Rinv = inv(SymOper.R);
            tinv = SymOper.t*Rinv.';
            SymOper = Oper(Rinv, Uinv,tinv,'conjugate',SymOper.conjugate,...
                'antisymmetry',SymOper.antisymmetry,...
                'strict_eq',SymOper.strict_eq...
                );
        end
    end
    %% modify
    methods 
        function OperObj = attachRm(OperObj,Rm)
            for i =1:numel(OperObj)
                OperObj(i).Rf = double(Oper.Rc2Rf(OperObj(i).R,Rm));
                OperObj(i).tf = OperObj(i).t/Rm;
            end
        end
    end
    %% pretty print
    methods
        % string
        function SymOper_Str = string(SymOper)
            SymOper_Str = pretty(SymOper,'full',true);
        end
        % latex
        function SymOper_latex= latex(SymOper)
            SymOper_latex = pretty(SymOper,'full',true,'latex',true);
        end
        % repr pretty
        function SymOper_pretty= repr_pretty(SymOper,Cycle)
            SymOper_pretty = latex(SymOper,'full',false);
        end
        % repr latex
        function SymOper_latex= repr_latex(SymOper)
            SymOper_latex = pretty(SymOper,'full',false,'latex',true);
        end
        % pretty
        function name= pretty(SymOper,options)
            
            % Human readable group element names
            %
            %     Return a human readable string representation of PointGroupElement
            %
            %     Parameters
            %     ----------
            %
            %     SymOper : Oper class
            %         Symmetry Operation to be represented.
            %     full : bool (default False)
            %         Whether to return a full representation.
            %
            %         The default short representation only contains the real space action
            %         and the symbol of the Altland-Zirnbauer part of the symmetry (see below).
            %         The full representation presents the symmetry action on the Hamiltonian
            %         and the unitary Hilbert-space action if set.
            %     Latex : bool (default False)
            %         Whether to output LateX formatted string.
            %
            %     Returns
            %     -------
            %     name : string
            %         In the short representation it is a sting `rot_name + az_name`.
            %         In the long representation the first line is the action on the
            %         Hamiltonian, the second line is `rot_name` and the third line
            %         is the unitary action as a matrix, if set.
            %
            %         `rot_name` can be:
            %         - `1` for identity
            %         - `I` for inversion (in 1D mirror is the same as inversion)
            %         - `R(angle)` for 2D rotation
            %         - `R(angle, axis)` for 3D rotation (axis is not normalized)
            %         - `M(normal)` for mirror
            %         - `S(angle, axis)` for 3D rotoinversion (axis is not normalized)
            %
            %         `az_name` can be:
            %         - `T` for time-reversal (antiunitary symmetry)
            %         - `P` for particle-hole (antiunitary antisymmetry)
            %         - `C` for chiral (unitary antisymmetry)
            %         - missing if the symmetry is unitary
            %
            arguments
                SymOper Oper ;
                options.full = true;
                options.Latex = false;
            end
            diminsion = size(SymOper.R,1);
            switch diminsion
                case 1
                    if SymOper.R(1,1) == 1
                        rot_name = '1';
                    else
                        rot_name = 'I';
                    end
                case 2
                    if Oper.isclose(det(SymOper.R),1)
                        % pure rotation Proper Rotation
                        theta = atan2(SymOper.R(1,2),SymOper.R(1,1));
                        if Oper.isclose(theta,0)
                            rot_name = '1';
                        else
                            if options.Latex
                                rot_name = "R\left("+Oper.name_angle(theta,options.Latex)+"\right)";
                            else
                                rot_name = "R("+Oper.name_angle(theta)+")";
                            end
                        end
                    elseif Oper.isclose(det(SymOper.R),-1)
                        % Improper Rotation
                        [val, vec] = eig(SymOper.R);
                        [val,vec] =  Oper.sorteig(vec,val);
                        if SymOper.allclose(diag(vec).',[-1,1])
                            n = val(:,1).';
                            if options.Latex
                                rot_name = "M\left("+mat2str(Oper.round_axis(n))+"\right)";
                            else
                                rot_name = "M("+""+mat2str(Oper.round_axis(n))+")";
                            end
                        else
                            error('?');
                        end
                    else
                       error('det ~= \pm 1!');
                    end
         
                case 3
                    if Oper.isclose((det((SymOper.R))),1)
                        % pure rotation Proper Rotation
                        [n,theta] = Oper.Rotation2nTheta(SymOper.R);
                        if Oper.isclose(theta,0)
                            rot_name = '1';
                        else
                            if options.Latex
                                rot_name = "R\left("+Oper.name_angle(theta,options.Latex)+","+mat2str(Oper.round_axis(n))+"\right)";
                            else
                                rot_name = "R("+Oper.name_angle(theta)+","+mat2str(Oper.round_axis(n))+")";
                            end
                        end
                    elseif Oper.isclose((det((SymOper.R))),-1)
                        % rotoinversion
                        [n,theta] = Oper.Rotation2nTheta(-SymOper.R);
                        if Oper.isclose(theta,0)
                            % inversion
                            rot_name = 'I';
                        elseif Oper.isclose(theta,pi)
                            % Mirror
                            if options.Latex
                                rot_name = "M\left("+mat2str(Oper.round_axis(n))+"\right)";
                            else
                                rot_name = "M("+""+mat2str(Oper.round_axis(n))+")";
                            end
                        else
                            % Generic rotoinversion
                            if options.Latex
                                rot_name = "S\left("+Oper.name_angle(theta,latex)+","+mat2str(Oper.round_axis(n))+"\right)";
                            else
                                rot_name = "S("+Oper.name_angle(theta)+","+mat2str(Oper.round_axis(n))+")";
                            end
                        end
                    else
                        error('det ~= 1!');
                    end
                    
            end
            if options.full
                if options.Latex
                    name = "\begin{aligned} U H(\mathbf{{k}})";
                    if SymOper.conjugate
                        name = name + "^*";
                    end
                    name = name + " U^{{-1}} &= ";
                    if SymOper.antisymmetry
                        name = name + "-";
                    end  
                    name = name + "H(";
                    if SymOper.conjugate
                        name = name + "-";
                    end
                    name = name + "R\mathbf{{k}}) \\";
                else
                    name = "U·H(k)";
                    if SymOper.conjugate
                        name = name + "*";
                    end
                    name = name + "·U^-1 = ";
                    if SymOper.antisymmetry
                        name = name + "-";
                    end
                    name = name + "H(";
                    if SymOper.conjugate
                        name = name + "-";
                    end
                    name = name + "R·k)\n";
                end
                if options.Latex
                    name = name + 'R &= '+rot_name + '\\';
                else
                    name = name + 'R = '+rot_name + '\n';
                end
                if ~isnan(SymOper.U)
                    if isa(SymOper.U,'sym')
                        if options.Latex
                            Umat = latex(SymOper.U);
                            name = name + "U &= "+Umat+'\end{aligned}';
                        else
                            name = name + "U = "+mat2str(string(SymOper.U))+'\n';
                        end
                    else
                        if options.Latex
                            Umat = Oper.mat2latex(SymOper.U,6);
                            name = name + "U &= "+Umat+'\end{aligned}';
                        else
                            name = name + "U = "+mat2str(roundn(SymOper.U,-4),3)+'\n';
                        end
                    end
                else
                    if options.Latex
                        name = name + '\end{aligned}';
                    else
                        name = name + '\n';
                    end
                end
            else
                if SymOper.conjugate && ~SymOper.antisymmetry
                    if options.Latex
                        az_name = "\mathcal{T}";
                    else
                        az_name = "T";
                    end 
                elseif SymOper.conjugate && SymOper.antisymmetry
                    if options.Latex
                        az_name = "\mathcal{P}";
                    else
                        az_name = "P";
                    end
                elseif ~SymOper.conjugate && SymOper.antisymmetry
                    if options.Latex
                        az_name = "\mathcal{C}";
                    else
                        az_name = "C";
                    end
                else
                    az_name = "";
                end
                if strcmp(rot_name ,'1') && ~strcmp(az_name,"")
                    if ~strcmp(az_name,"")
                        name = az_name + " "+az_name;
                    else
                        name = az_name + ""+az_name;
                    end
                else
                    if ~strcmp(az_name,"")
                        name = rot_name + " "+az_name;
                    else
                        name = rot_name + ""+az_name;
                    end
                end
                
            end
%             if options.Latex
%                 name = "$"+name+"$";
%             else
%                 
%             end
        end
    end
    %% Group theory
    methods
        function symmetry_from_permutation()
            %     Construct symmetry operator for lattice systems with multiple sites.
            %
            %     Parameters
            %     ----------
            %     R : real space rotation
            %     perm : dict : {site: image_site}
            %         permutation of the sites under the symmetry action
            %     norbs : OrderedDict : {site : norbs_site} or tuple of tuples ((site, norbs_site), )
            %         sites are ordered in the order specified, with blocks of size norbs_site
            %         corresponding to each site.
            %     onsite_action : dict : {site: ndarray} or ndarray or None
            %         onsite symmetry action, such as spin rotation for each site. If only one
            %         array is specified, it is used for every site. If None (default), identity
            %         is used on every site. Size of the arrays must be consistent with norbs.
            %     antiunitary, antisymmetry : bool
            %
            %     Returns
            %     -------
            %     g : PointGroupElement
            %         PointGroupElement corresponding to the operation.
            %
            %     Notes:
            %     ------
            %     Sites can be indexed by any hashable identifiers, such as integers or stings.
            %     
            % if needed
        end
    end
    %% protected method
    methods (Access= protected)
    end
    %% Source method
    methods (Static)
        function S_mat = spin_matrices(s, include_0)
            %     Construct spin-s matrices for any half-integer spin.
            %
            %     Parameters
            %     ----------
            %
            %     s : float or int
            %         Spin representation to use, must be integer or half-integer.
            %     include_0 : bool (default False)
            %         If `include_0` is True, S[0] is the identity, indices 1, 2, 3
            %         correspond to x, y, z. Otherwise indices 0, 1, 2 are x, y, z.
            %
            %     Returns
            %     -------
            %
            %         Sequence of spin-s operators in the standard spin-z basis.
            %         mat of shape `( 2*s + 1, 2*s + 1,3)`, or if `include_0` is True
            %         `( 2*s + 1, 2*s + 1,4)`.
            arguments
                s double {Oper.mustBeHalfInteger(s)} =1/2;
                include_0 logical = false;
            end
            S = 2*s+1;
            Sz= 1/2 * diag(S-1:-2:-S);
            
            S_plus = zeros(S-1,1);
            % first diagonal for general s from en.wikipedia.org/wiki/Spin_(physics)
            for i  = 1:S-1
                S_plus(i) =  sqrt((s+1)*2*i-i*(i+1))/2 ;
            end
            Sx  = diag(S_plus,1) +diag(S_plus,-1);
%             disp(Sx);
            Sy  = -1i*diag(S_plus,1) +1i*diag(S_plus,-1); 
            if include_0
                S_mat(:,:,1) = eye(S);
                S_mat(:,:,2) = Sx;
                S_mat(:,:,3) = Sy;
                S_mat(:,:,4) = Sz;
            else
                S_mat(:,:,1) = Sx;
                S_mat(:,:,2) = Sy;
                S_mat(:,:,3) = Sz;
            end
        end
        function J_mat = spin_matrices_from_orb(quantumL,include_0,strict)
            arguments
                quantumL double =[1,0,0,0.5,0.5;1,0,0,0.5,-0.5];
                include_0 logical = false;
                strict logical = true;
            end
            Norb = size(quantumL,1);
%             S = 2*s+1;
%             Sz= 1/2 * diag(S-1:-2:-S);
            Sz =  diag(quantumL(:,5));
            S_plus = zeros(Norb);
            S_minus = zeros(Norb);
            % first diagonal for general s from en.wikipedia.org/wiki/Spin_(physics)
            for i  = 1:Norb
                for j = 1:Norb
                    J_l = quantumL(i,2)+quantumL(i,4);
                    jz_l = quantumL(i,5);
                    J_r = quantumL(j,2)+quantumL(j,4);
                    jz_r = quantumL(j,5);
                    if strict
                        S_plus(i,j) =  Oper.braket([J_l,jz_l ],'S+',[J_r,jz_r ]);
                        S_minus(i,j) =  Oper.braket([J_l,jz_l ],'S-',[J_r,jz_r ]);
                    else
                        S_plus(i,j) =  Oper.braket([J_l,jz_l ],'^S+',[J_r,jz_r ]);
                        S_minus(i,j) =  Oper.braket([J_l,jz_l ],'^S-',[J_r,jz_r ]);
                    end
                end
            end
            Sx  = 1/2 * (S_plus+S_minus);
            %             disp(Sx);
            Sy  = 1/(2*1i) * (S_plus- S_minus);
            if include_0
                J_mat(:,:,1) = eye(S);
                J_mat(:,:,2) = Sx;
                J_mat(:,:,3) = Sy;
                J_mat(:,:,4) = Sz;
            else
                J_mat(:,:,1) = Sx;
                J_mat(:,:,2) = Sy;
                J_mat(:,:,3) = Sz;
            end
        end
        function L_mat = L_matrices(d, l)
            %     Construct real space rotation generator matrices in d=2 or 3 dimensions.
            %     Can also be used to get angular momentum operators for real atomic orbitals
            %     in 3 dimensions, for p-orbitals use `l=1`, for d-orbitals `l=2`. The basis
            %     of p-orbitals is `p_x`, `p_y`, `p_z`, for d-orbitals `d_{x^2 - y^2}`,
            %     `d_{3 z^2 - r^2}`, `d_{xy}`, `d_{yz}`, `d_{zx}`. The matrices are all
            %     purely imaginary and antisymmetric.
            %     To generate finite rotations, use 'spin_rotation(n, L)'.

            arguments
                d double {mustBeMember(d,[2,3])} =3;
                l double {mustBeMember(l,[1,2])} =1;
            end
            
            if d == 2 && l==1
                L_mat= 1i * ([[0, -1];...
                    [1, 0]]);
            elseif d == 3 &&l==1
                L_mat(:,:,1)= 1i * ([[0, 0, 0];...
                    [0, 0, -1];...
                    [0, 1, 0]]);...
                    
                L_mat(:,:,2)= 1i * ([[0, 0, 1];...
                    [0, 0, 0];...
                    [-1, 0, 0]]);...
                    
                L_mat(:,:,3)= 1i * ([[0, -1, 0];...
                    [1, 0, 0];...
                    [0, 0, 0]]);
            elseif d == 3 && l==2
                s3 = sqrt(3);
                L_mat(:,:,1)= 1i * ([[0, 0, 0, -1, 0];...
                    [0, 0, 0, -s3, 0];...
                    [0, 0, 0, 0, 1];...
                    [1, s3, 0, 0, 0];...
                    [0, 0, -1, 0, 0]]);...
                    
                L_mat(:,:,2)= 1i * ([[0, 0, 0, 0, -1];...
                    [0, 0, 0, 0, s3];...
                    [0, 0, 0, -1, 0];...
                    [0, 0, 1, 0, 0];...
                    [1, -s3, 0, 0, 0]]);...
                    
                L_mat(:,:,3)= 1i * ([[0, 0, -2, 0, 0];...
                    [0, 0, 0, 0, 0];...
                    [2, 0, 0, 0, 0];...
                    [0, 0, 0, 0, 1];...
                    [0, 0, 0, -1, 0]]);
            else
                raise ValueError('Only 2 and 3 dimensions are supported.')
            end
        end
        function U = spin_rotation(n, s, roundint)
            %     Construct the unitary spin rotation matrix for rotation specified by the
            %     vector n (in radian units) with angular momentum `s`, given by
            %     `U = exp(-i n⋅s)`.
            %
            %     Parameters
            %     ----------
            %
            %     n : iterable
            %         Rotation vector. Its norm is the rotation angle in radians.
            %     s : float or sequence of arrays
            %         Spin representation to use for the unitary action of the
            %         operator. If float is provided, it should be integer or half-integer
            %         specifying the spin representation in the standard basis, see `spin_matrices`.
            %         Otherwise a sequence of 3 arrays of identical square size must be provided
            %         representing 3 components of the angular momentum operator.
            %     roundint : bool (default False)
            %         If roundint is True, result is converted to integer tinyarray if possible.
            %
            %     Returns
            %     -------
            %     U : mat
            %         Unitary spin rotation matrix of the same shape as the spin matrices or
            %         `(2*s + 1, 2*s + 1)`.
            arguments
                n ;
                s double  =3/2;
                roundint logical = false;
            end
            if length(s) == 1
                J_mat = Oper.spin_matrices(s);
            else
                J_mat = s;
            end
            % Make matrix exponential for rotation representation
            expmat = Oper.tensordot_naive(n,J_mat);
            U = expm(1i*expmat);
            if roundint
                Ur = round(real(U));
                if Oper.allclose(U,Ur)
                    U = round(Ur); 
                else
                    error('?');
                end   
            end
        end
        function phi = braket(phi1,Oper,phi2)
            switch Oper
                case 'S+'
                    S = phi2(1);
                    m = phi2(2);
                    phi = sqrt((S-m)*(S+m+1));
                    if m + 1 <= S 
                        phi2(2) = m+1;
                    else
                        phi2(2) = -S;
                    end
                    if isequal(phi1,phi2)
                        
                    else
                        phi= 0;
                    end
                case 'S-'
                    S = phi2(1);
                    m = phi2(2);
                    phi = sqrt((S+m)*(S-m+1));
                    if m - 1 >= -S
                        phi2(2) = m-1;
                    else
                        phi2(2) = S;
                    end
                    if isequal(phi1,phi2)
                        
                    else
                        phi= 0;
                    end
                case '^S+'
                    S = phi2(1);
                    m = phi2(2);
                    phi = sqrt((S-m)*(S+m+1));
                    if m + 1 <= S
                        phi2(2) = S;
                    else
                        phi2(2) = -S;
                    end
                    if isequal(phi1,phi2)
                        
                    else
                        phi= 0;
                    end
                case '^S-'
                    S = phi2(1);
                    m = phi2(2);
                    phi = sqrt((S+m)*(S-m+1));
                    if m - 1 >= -S
                        phi2(2) = -S;
                    else
                        phi2(2) = S;
                    end
                    if isequal(phi1,phi2)
                        
                    else
                        phi= 0;
                    end
            end

        end
    end
    %% Hidden method
    %     methods (Static,Hidden,Access= protected)
    methods(Static)
        % ---------------------   tools function   ------------------------
        function angle = name_angle(theta, Latex)
            arguments
                theta  ;
                Latex logical = false;
            end
            if isa(theta,'sym')
                frac = simplify(theta,'IgnoreAnalyticConstraints',true);
            else
                frac = rat(theta / pi, 1e-4);
                frac = simplify(str2sym(frac))*pi;
            end
            if Latex
                angle = latex(frac);
            else
                angle = string(frac);
            end
        end
        function str = mat2latex(mat,accuracy)
            mat = roundn(mat,-accuracy-1);
            Size = size(mat);
            str =  "\left(\begin{array}{";
            for j =1:Size(2)
                str = str+ 'c';
            end
            str = str+ '} ';
            for i = 1:Size(1)
                str = str + Oper.num2latex(mat(i,1),accuracy);
                for j = 2:Size(2)
                    str = str +'&'+Oper.num2latex(mat(i,j),accuracy);
                end
                str = str+ '\\';
            end
            str = str + '\end{array}\right)';
        end
        function str = num2latex(num,accuracy)
            arguments
                num double 
                accuracy double  =6;
            end
            digits(accuracy);
            theta = angle(num);
            r = abs(num);
            theta_name = Oper.name_angle(theta,1);
            theta_expr = "e^{i "+theta_name+"}";
            if Oper.isclose(imag(num),0)
                str = latex(sym(num));
            elseif Oper.isclose(real(num),0)
                str = latex(sym(num));
            else
                if Oper.isclose(r,1)
                    str = theta_expr;
                else
                    str = latex(sym(r))+theta_expr;
                end
            end
            
        end
        function C = isclose(A,B)
            if isa(A,'sym') ||  isa(B,'sym') 
               C = logical(simplify(sym(A))==simplify(sym(B)));
            else
               C = (A-B).^2<1e-12 ; 
            end
%             if (A-B).^2<1e-6 
%                 C = true;
%             else
%                 C = false;
%             end
            
        end
        function C = allclose(A,B)
            if all(all(Oper.isclose(A,B)))
                C = true;
            else
                C = false;
            end
        end
        function n = round_axis(n)
            % Try to find integer axis
            % if need
            if isa(n,'sym')
                n = simplify(n,'IgnoreAnalyticConstraints',true);
                n = string(n);
            else
                n = roundn(n, -2);
            end
        end
        function M = tensordot_naive(A,B,sizeLast)
            if nargin < 3
                sizeA = size(A);
                sizeB = size(B);
                sizeLast = sizeB(end);
            end
            if length(sizeA) ~= length(sizeB)
                A = reshape(A,sqrt(length(A(:))/sizeLast),sqrt(length(A(:))/sizeLast),sizeLast);
            end
            %B = reshape(B,sqrt(length(B(:))/sizeLast),[],sizeLast);
            M = A(:,:,1)* B(:,:,1);
            for i = 2:sizeLast
                M = M + A(:,:,i)* B(:,:,i);
            end
        end
        function RotationMat = nThetad2Rotation(thetad,n)
            if nargin < 2
                 n = [0,0,1];
            end
%             % Compute rotation matrices
%             cth = cos(theta);
%             sth = sin(theta);
%             vth = (1 - cth);
%             v = n;
%             % Preallocate input vectors
%             vx = zeros(1,1,numInputs,'like',axang);
%             vy = vx;
%             vz = vx;
%             
%             % Shape input vectors in depth dimension
%             vx(1,1,:) = v(:,1);
%             vy(1,1,:) = v(:,2);
%             vz(1,1,:) = v(:,3);
%             
%             % Explicitly specify concatenation dimension
%             tempR = cat(1, vx.*vx.*vth+cth,     vy.*vx.*vth-vz.*sth, vz.*vx.*vth+vy.*sth, ...
%                 vx.*vy.*vth+vz.*sth, vy.*vy.*vth+cth,     vz.*vy.*vth-vx.*sth, ...
%                 vx.*vz.*vth-vy.*sth, vy.*vz.*vth+vx.*sth, vz.*vz.*vth+cth);
%             
%             R = reshape(tempR, [3, 3, length(vx)]);
%             R = permute(R, [2 1 3]);
            n_x = n(1);
            n_y = n(2);
            n_z = n(3);
            %% left
            % RotationMat(1,1)=n_x^2*(1-cosd(theta))+cosd(theta);
            % RotationMat(1,2)=n_x*n_y*(1-cosd(theta))+n_z*sind(theta);
            % RotationMat(1,3)=n_x*n_z*(1-cosd(theta))-n_y*sind(theta);
            % RotationMat(2,1)=n_y*n_z*(1-cosd(theta))-n_z*sind(theta);
            % RotationMat(2,2)=n_y^2*(1-cosd(theta))+cosd(theta);
            % RotationMat(2,3)=n_y*n_z*(1-cosd(theta))+n_x*sind(theta);
            % RotationMat(3,1)=n_z*n_x*(1-cosd(theta))+n_y*sind(theta);
            % RotationMat(3,2)=n_z*n_y*(1-cosd(theta))-n_x*sind(theta);
            % RotationMat(3,3)=n_z^2*(1-cosd(theta))+cosd(theta);
            %% right
            RotationMat(1,1)=n_x^2*(1-cosd(thetad))+cosd(thetad);
            RotationMat(1,2)=n_x*n_y*(1-cosd(thetad))-n_z*sind(thetad);
            RotationMat(1,3)=n_x*n_z*(1-cosd(thetad))+n_y*sind(thetad);
            RotationMat(2,1)=n_y*n_x*(1-cosd(thetad))+n_z*sind(thetad);
            RotationMat(2,2)=n_y^2*(1-cosd(thetad))+cosd(thetad);
            RotationMat(2,3)=n_y*n_z*(1-cosd(thetad))-n_x*sind(thetad);
            RotationMat(3,1)=n_z*n_x*(1-cosd(thetad))-n_y*sind(thetad);
            RotationMat(3,2)=n_z*n_y*(1-cosd(thetad))+n_x*sind(thetad);
            RotationMat(3,3)=n_z^2*(1-cosd(thetad))+cosd(thetad);
        end
        function RotationMat = nTheta2Rotation(theta,n)
            % nTheta2Rotation Convert axis-angle rotation representation to rotation matrix
            %   RotationMat = nTheta2Rotation(theta,n) converts a 3D rotation given in axis-angle form,
            %   theta,n, to an orthonormal rotation matrix, R. theta is an N-by-4
            %   matrix of N axis-angle rotations. The first three elements of every
            %   row specify the rotation axis and the last element defines the rotation
            %   angle (in radians).
            %   The output, R, is an 3-by-3-by-N matrix containing N rotation matrices.
            %   Each rotation matrix has a size of 3-by-3 and is orthonormal.
            %
            %   Example:
            %      % Convert a rotation from axis-angle to rotation matrix
            %      axang = [0 1 0 pi/2];
            %      R = axang2rotm(axang)
            if nargin < 2
                n = [0,0,1];
            end
            n = normalize(n);
            thetad = rad2deg(theta);
            RotationMat = Oper.nThetad2Rotation(thetad,n);
        end
        function [n,theta]= Rotation2nTheta(R,Rm)
            % Convert 3D rotation matrix to axis and angle
            % by my own
%             rotation_cart =R;
%             theta = round(real(acosd((trace(rotation_cart)-1)/2)));
%             if abs(theta- 0) < 1e-8
%                 n = [1 0 0];
%                 
%             elseif abs(theta - 180) < 1e-8
%                 n = [sqrt(rotation_cart(1,1)/2+0.5) ,...
%                     ((((R(2,1)) >= 0)-0.5)/0.5)*sqrt(rotation_cart(2,2)/2+0.5) ,...
%                     ((((R(3,1)) >= 0)-0.5)/0.5)*sqrt(rotation_cart(3,3)/2+0.5)];
%             else
%                 n = [R(3,2)-R(2,3),R(1,3)-R(3,1),R(2,1)-R(1,2)]/(2*sind(theta));
%             end
            % from 
            arguments
                R  {Oper.mustBeSize(R,[3,3])} 
                Rm = [1 0 0;0 1 0;0 0 1];
            end
            if isequal(Rm,[1 0 0;0 1 0;0 0 1])
                %R = R.';
            else
                R = inv(Rm.' * R / Rm.');
                % need test
            end
            % 
            %rotation_cart = inv(rotation_cart_inv );
            if Oper.isclose(det(R),-1)
                R = -R;
                disp('det = -1');
            else
                R = R/det(R);
            end
            
            if Oper.allclose(R,eye(3))
                theta = 0;
                n = [1, 0 ,0];
                return
            end
            L_mat = Oper.L_matrices();
            % left ?
            % n = real(1i*[trace(L_mat(:,:,1)*R),trace(L_mat(:,:,2)*R),trace(L_mat(:,:,3)*R)]);
            % right
            n = real(-1i*[trace(L_mat(:,:,1)*R),trace(L_mat(:,:,2)*R),trace(L_mat(:,:,3)*R)]);
            
            % Choose direction to minimize number of minus signs
            absn = norm(n) * sign(sum(n));
            if Oper.isclose(absn,0)
                % n is zero for 2-fold rotation, find eigenvector with +1
                if isa(R,'sym')
                    [val, vec] = eig(simplify(R));
                else
                    [val, vec] = eig((R));
                end
                [val,vec] =  Oper.sorteig(vec,val);
                vec = diag(vec);
                if Oper.isclose(vec(end),1)
                    n = val(:, end).';
                    % Choose direction to minimize number of minus signs
                    n = n/ sign(sum(n));
                end
            else
                n = n/absn;
            end
            % adapt matlab robotic box ROTM2AXANG 
            try
                theta = real(acos(complex((1/2)*(R(1,1,:)+R(2,2,:)+R(3,3,:)-1))));
            catch
                theta = real(acos(((1/2)*(R(1,1,:)+R(2,2,:)+R(3,3,:)-1))));
            end
            %theta = atan2(absn,(real(trace(R)) - 1));
            if isa(n, 'sym')
                n = simplify(n,'IgnoreAnalyticConstraints',true);
                theta = simplify(theta,'IgnoreAnalyticConstraints',true);
            end
        end
        function [n,thetad]= Rotation2nThetad(R,Rm)
            arguments
                R double {mustBeReal,Oper.mustBeSize(R,[3,3])}
                Rm = [1 0 0;0 1 0;0 0 1];
            end
            [n,theta]= Oper.Rotation2nTheta(R,Rm);
            thetad = str2sym(Oper.name_angle(theta, false));
        end
        function Rf = Rc2Rf(Rc,Rm)
            Rf =int8( Rm*Rc/Rm);
        end
        function Rft = Rct2Rft(Rct,Rm)
            
        end
        function Rc = Rf2Rc(Rf,Rm)
            rotation_cart_inv = Rm.' * rotation / Rm.';
            Rc = inv(rotation_cart_inv );
        end
        function Rct = Rft2Rct(Rft,Rm)
        end
        %
        function [prop,Coeff] = prop_to_id(A)
            % Test if A is proportional to the identity matrix
            % and return the factor as well
            if Oper.isclose(A(1,1), 0)
                if matallclose(A, zeros(length(A)))
                    prop = true;
                    Coeff = 0;
                else
                    prop = false;
                    Coeff = 0;
                end
            else
                if Oper.allclose( A /A(1,1), eye(length(A)))
                    prop = true;
                    Coeff = A(1,1);
                else
                    prop = false;
                    Coeff = 0;
                end
            end
            
        end
        % Custom validator functions
        function mustBeOfClass(input,className)
            % Test for specific class name
            cname = class(input);
            if ~strcmp(cname,className)
                eid = 'Class:notCorrectClass';
                msg = ['Input must be of class ', cname];
                throwAsCaller(MException(eid,msg))
            end
        end
        function mustBeEqualSize(a,b)
            % Test for equal size
            if ~(isequal(size(a),size(b)))
                eid = 'Size:notEqual';
                msg = 'Inputs must have same size.';
                throwAsCaller(MException(eid,msg))
            end
        end
        function mustBeSize(a,b)
            % Test for size
            if sum(~(size(a)==b))
                eid = 'Size:notRequired';
                msg = 'Inputs must have size of '+mat2str(b);
                throwAsCaller(MException(eid,msg))
            end
        end
        function mustBeDims(input,numDims)
            % Test for number of dimensions
            if ~isequal(length(size(input)),numDims)
                eid = 'Size:wrongDimensions';
                msg = ['Input must have dimensions: ',num2str(numDims)];
                throwAsCaller(MException(eid,msg))
            end
        end
        function mustBeHalfInteger(A)
            if ~isnumeric(A) && ~islogical(A)
                throwAsCaller(createValidatorException('MATLAB:validators:mustBeNumericOrLogical'));
            end
            
            if ~isreal(A)
                throwAsCaller(createValidatorException('MATLAB:validators:mustBeReal'));
            end
            A = 2*A;
            if ~all(isfinite(A), 'all') || ~all(A == floor(A), 'all')
                throwAsCaller(createValidatorException('MATLAB:validators:mustBeHalfInteger'));
            end
        end
        % tool
        function [Asort,Usort] = sorteig(U,A)
            if nargin <2
                mode = 'eigenval';
                NUM_WAN = length(U);
            else
                mode = 'whole';
                NUM_WAN = length(A);
            end
            NBANDS = length(U);
            
            if strcmp(mode,'whole')
                % %--step9--:按从大到小的特征值顺序排序重新组合对应的特征向量
                
                if isa(U,'sym')
                    Asort=sym(zeros(NUM_WAN ,NBANDS ));
                    SortTmp=diag(simplify(U));%抽取特征值
                else
                    Asort=zeros(NUM_WAN ,NBANDS );
                    SortTmp=diag(U);%抽取特征值
                end
                [Usort,IJ]=sort(double(SortTmp),1,'ComparisonMethod','real');
                for jj=1:NBANDS
                    Asort(:,jj)=A(:,IJ(jj));%取特征向量的列向量
                end
                Usort = diag(Usort);
            elseif strcmp(mode,'eigenval')
                % %--step9--:按从大到小的特征值顺序排序重新组合对应的特征向量
                SortTmp=diag(U);%抽取特征值
                [Usort,~]=sort(SortTmp,1,'ComparisonMethod','real');
                Usort = diag(Usort);
                Asort = [];
            end
        end
    end
    methods(Static) % Copyby Robotics toolbox
        function eul = Rotation2eul( R, Rm )
            %calculateEulerAngles Calculate Euler angles from rotation matrix
            %   EUL = Rotation2eul(R, Rm) calculates the Euler angles, EUL,
            %   corresponding to the input rotation matrix, R. The Euler angles follow
            %   the axis order specified in SEQ.
            arguments
                R  {mustBeReal,Oper.mustBeSize(R,[3,3])}
                Rm = [1 0 0;0 1 0;0 0 1];
            end
            if isequal(Rm,[1 0 0;0 1 0;0 0 1])
                %R = R.';
            else
                R = inv(Rm.' * R / Rm.');
                % need test
            end
            if isa(R,'double')
                try
                    eul = rotm2eul(R,'ZYZ');
                catch
                    % Preallocate output
                    eul = zeros(1, 3, size(R,3), 'like', R);  %#ok<PREALL>
                    nextAxis = [2, 3, 1, 2];
                    % Pre-populate settings for different axis orderings
                    % Each setting has 4 values:
                    %   1. firstAxis : The right-most axis of the rotation order. Here, X=1,
                    %      Y=2, and Z=3.
                    %   2. repetition : If the first axis and the last axis are equal in
                    %      the sequence, then repetition = 1; otherwise repetition = 0.
                    %   3. parity : Parity is 0 if the right two axes in the sequence are
                    %      YX, ZY, or XZ. Otherwise, parity is 1.
                    %   4. movingFrame : movingFrame = 1 if the rotations are with
                    %      reference to a moving frame. Otherwise (in the case of a static
                    %      frame), movingFrame = 0.
                    seqSettings.ZYZ = [3, 1, 1, 1];
                    % Retrieve the settings for a particular axis sequence
                    setting = seqSettings.('ZYZ');
                    firstAxis = setting(1);
                    repetition = setting(2);
                    parity = setting(3);
                    movingFrame = setting(4);
                    % Calculate indices for accessing rotation matrix
                    i = firstAxis;
                    j = nextAxis(i+parity);
                    k = nextAxis(i-parity+1);
                    if repetition
                        % Find special cases of rotation matrix values that correspond to Euler
                        % angle singularities.
                        sy = sqrt(R(i,j,:).*R(i,j,:) + R(i,k,:).*R(i,k,:));
                        singular = sy < 10 * eps(class(R));
                        % Calculate Euler angles
                        eul = [atan2(R(i,j,:), R(i,k,:)), atan2(sy, R(i,i,:)), atan2(R(j,i,:), -R(k,i,:))];
                        % Singular matrices need special treatment
                        numSingular = sum(singular,3);
                        assert(numSingular <= length(singular));
                        if numSingular > 0
                            eul(:,:,singular) = [atan2(-R(j,k,singular), R(j,j,singular)), ...
                                atan2(sy(:,:,singular), R(i,i,singular)), zeros(1,1,numSingular,'like',R)];
                        end
                    else
                        % Find special cases of rotation matrix values that correspond to Euler
                        % angle singularities.
                        sy = sqrt(R(i,i,:).*R(i,i,:) + R(j,i,:).*R(j,i,:));
                        singular = sy < 10 * eps(class(R));
                        % Calculate Euler angles
                        eul = [atan2(R(k,j,:), R(k,k,:)), atan2(-R(k,i,:), sy), atan2(R(j,i,:), R(i,i,:))];
                        % Singular matrices need special treatment
                        numSingular = sum(singular,3);
                        assert(numSingular <= length(singular));
                        if numSingular > 0
                            eul(:,:,singular) = [atan2(-R(j,k,singular), R(j,j,singular)), ...
                                atan2(-R(k,i,singular), sy(:,:,singular)), zeros(1,1,numSingular,'like',R)];
                        end
                    end
                    if parity
                        % Invert the result
                        eul = -eul;
                    end
                    if movingFrame
                        % Swap the X and Z columns
                        eul(:,[1,3],:)=eul(:,[3,1],:);
                    end
                end
            elseif isa(R,'sym')
                firstAxis = 3;parity =1;movingFrame=1;
                %eul = sym(zeros(1, 3, size(R,3), 'like', R));
                nextAxis = [2, 3, 1, 2];
                % Calculate indices for accessing rotation matrix
                i = firstAxis;
                j = nextAxis(i+parity);
                k = nextAxis(i-parity+1);
                sy = simplify(sqrt(R(i,j,:).*R(i,j,:) + R(i,k,:).*R(i,k,:)));
                % Calculate Euler angles
                eul = [atan2(R(i,j,:), R(i,k,:)), atan2(sy, R(i,i,:)), atan2(R(j,i,:), -R(k,i,:))];
                singular = logical(sy ==sym(0));
                % Singular matrices need special treatment
                numSingular = sum(singular,3);
                assert(numSingular <= length(singular));
                if numSingular > 0
                    eul(:,:,singular) = [angle(-R(j,k,singular)*1i+R(j,j,singular)), ...
                        angle(1i*sy(:,:,singular)+R(i,i,singular)), zeros(1,1,numSingular,'like',R)];
                end
                if parity
                    % Invert the result
                    eul = -eul;
                end
                if movingFrame
                    % Swap the X and Z columns
                    eul(:,[1,3],:)=eul(:,[3,1],:);
                end
                eul = simplify(eul,'IgnoreAnalyticConstraints',true);
            end

        end
        function R = axang2rotm( axang )
            %AXANG2ROTM Convert axis-angle rotation representation to rotation matrix
            %   R = AXANG2ROTM(AXANG) converts a 3D rotation given in axis-angle form,
            %   AXANG, to an orthonormal rotation matrix, R. AXANG is an N-by-4
            %   matrix of N axis-angle rotations. The first three elements of every
            %   row specify the rotation axis and the last element defines the rotation
            %   angle (in radians).
            %   The output, R, is an 3-by-3-by-N matrix containing N rotation matrices.
            %   Each rotation matrix has a size of 3-by-3 and is orthonormal.
            %
            %   Example:
            %      % Convert a rotation from axis-angle to rotation matrix
            %      axang = [0 1 0 pi/2];
            %      R = axang2rotm(axang)
            %
            %   See also rotm2axang

            %   Copyright 2014-2019 The MathWorks, Inc.

            %#codegen

            %             robotics.internal.validation.validateNumericMatrix(axang, 'axang2rotm', 'axang', ...
            %                 'ncols', 4);

            % For a single axis-angle vector [ax ay az theta] the output rotation
            % matrix R can be computed as follows:
            % R =  [t*x*x + c	  t*x*y - z*s	   t*x*z + y*s
            %       t*x*y + z*s	  t*y*y + c	       t*y*z - x*s
            %       t*x*z - y*s	  t*y*z + x*s	   t*z*z + c]
            % where,
            % c = cos(theta)
            % s = sin(theta)
            % t = 1 - c
            % x = normalized axis ax coordinate
            % y = normalized axis ay coordinate
            % z = normalized axis az coordinate

            % Normalize the axis
            if isa(axang,'double')
                v = normalize(axang(:,1:3),'scale');
            else
                v =(axang(:,1:3));
            end

            % Extract the rotation angles and shape them in depth dimension
            numInputs = size(axang,1);
            theta = zeros(1,1,numInputs,class(axang));
            theta(1,1,:) = axang(:,4);

            % Compute rotation matrices
            cth = cos(theta);
            sth = sin(theta);
            vth = (1 - cth);

            % Preallocate input vectors
            vx = zeros(1,1,numInputs,'like',axang);
            vy = vx;
            vz = vx;

            % Shape input vectors in depth dimension
            vx(1,1,:) = v(:,1);
            vy(1,1,:) = v(:,2);
            vz(1,1,:) = v(:,3);

            % Explicitly specify concatenation dimension
            tempR = cat(1, vx.*vx.*vth+cth,     vy.*vx.*vth-vz.*sth, vz.*vx.*vth+vy.*sth, ...
                vx.*vy.*vth+vz.*sth, vy.*vy.*vth+cth,     vz.*vy.*vth-vx.*sth, ...
                vx.*vz.*vth-vy.*sth, vy.*vz.*vth+vx.*sth, vz.*vz.*vth+cth);

            R = reshape(tempR, [3, 3, length(vx)]);
            R = permute(R, [2 1 3]);
            if isa(R, 'sym')
                R = simplify(R);
            end
        end
        function eul = axang2eul(axang)
            R = Oper.axang2rotm( axang );
            eul = Oper.Rotation2eul(R);
        end
        function quat = rotm2quat( R )
            %ROTM2QUAT Convert rotation matrix to quaternion
            %   Q = ROTM2QUAT(R) converts a 3D rotation matrix, R, into the corresponding
            %   unit quaternion representation, Q. The input, R, is an 3-by-3-by-N matrix
            %   containing N orthonormal rotation matrices.
            %   The output, Q, is an N-by-4 matrix containing N quaternions. Each
            %   quaternion is of the form q = [w x y z], with w as the scalar number.
            %   Each element of Q must be a real number.
            %
            %   If the input matrices are not orthonormal, the function will
            %   return the quaternions that correspond to the orthonormal matrices
            %   closest to the imprecise matrix inputs.
            %
            %
            %   Example:a
            %      % Convert a rotation matrix to a quaternion
            %      R = [0 0 1; 0 1 0; -1 0 0];
            %      q = rotm2quat(R)
            %
            %   References:
            %   [1] I.Y. Bar-Itzhack, "New method for extracting the quaternion from a
            %       rotation matrix," Journal of Guidance, Control, and Dynamics,
            %       vol. 23, no. 6, pp. 1085-1087, 2000
            %
            %   See also quat2rotm

            %   Copyright 2014-2019 The MathWorks, Inc.

            %#codegen

            %robotics.internal.validation.validateRotationMatrix(R, 'rotm2quat', 'R');

            % Pre-allocate output
            if isa(R,'sym')
                quat = sym(zeros(size(R,3), 4, 'like', R));
            else
                quat = zeros(size(R,3), 4, 'like', R);
            end
            % Calculate all elements of symmetric K matrix
            K11 = R(1,1,:) - R(2,2,:) - R(3,3,:);
            K12 = R(1,2,:) + R(2,1,:);
            K13 = R(1,3,:) + R(3,1,:);
            K14 = R(3,2,:) - R(2,3,:);

            K22 = R(2,2,:) - R(1,1,:) - R(3,3,:);
            K23 = R(2,3,:) + R(3,2,:);
            K24 = R(1,3,:) - R(3,1,:);

            K33 = R(3,3,:) - R(1,1,:) - R(2,2,:);
            K34 = R(2,1,:) - R(1,2,:);

            K44 = R(1,1,:) + R(2,2,:) + R(3,3,:);

            % Construct K matrix according to paper
            K = [...
                K11,    K12,    K13,    K14;
                K12,    K22,    K23,    K24;
                K13,    K23,    K33,    K34;
                K14,    K24,    K34,    K44];

            K = K ./ 3;

            % For each input rotation matrix, calculate the corresponding eigenvalues
            % and eigenvectors. The eigenvector corresponding to the largest eigenvalue
            % is the unit quaternion representing the same rotation.
            for i = 1:size(R,3)
                [eigVec,eigVal] = eig(K(:,:,i));
                eigVal = diag(eigVal);
                [~,maxIdx] = max(real(eigVal));
                quat(i,:) = real([eigVec(4,maxIdx) eigVec(1,maxIdx) eigVec(2,maxIdx) eigVec(3,maxIdx)]);

                % By convention, always keep scalar quaternion element positive.
                % Note that this does not change the rotation that is represented
                % by the unit quaternion, since q and -q denote the same rotation.
                if quat(i,1) < 0
                    quat(i,:) = -quat(i,:);
                end
            end
        end
    end
end

function A = integer(A)
cycle = size(A);
for i = 1:cycle(1)
    for j = 1:cycle(2)
        delta = abs(A(i,j) - round(A(i,j)));
        if delta < 1e-6
            A(i,j) = round(A(i,j));
        end
    end
end
end