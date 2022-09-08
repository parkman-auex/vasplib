classdef wt
    % wt fileinput
    %   restore a wtfile

    %  TB_FILE namelist
    properties
        Particle = 'electron'
        Package = 'VASP'
        KPorTB = 'TB'
        Is_HrFile = .TRUE.
        Is_Sparse_Hr = .FALSE.
        Is_Sparse = .FALSE.
        Use_ELPA = .FALSE.
        vef = 0d0
        TB_FILE = 'wannier90_hr.dat'
    end

    %  CONTROL namelist
    properties
        BulkBand_calc = .FALSE.
        BulkBand_line_calc = .FALSE.
        BulkBand_unfold_line_calc = .FALSE.
        BulkBand_unfold_plane_calc = .FALSE.
        BulkFatBand_calc = .FALSE.
        BulkBand_plane_calc = .FALSE.
        BulkBand_cube_calc = .FALSE.
        BulkFS_calc = .FALSE.
        BulkFS_Plane_calc = .FALSE.
        BulkFS_plane_stack_calc = .FALSE.
        BulkGap_cube_calc = .FALSE.
        BulkGap_plane_calc = .FALSE.
        SlabBand_calc = .FALSE.
        SlabBandWaveFunc_calc = .FALSE.
        WireBand_calc = .FALSE.
        SlabSS_calc = .FALSE.
        SlabArc_calc = .FALSE.
        SlabQPI_calc = .FALSE.
        SlabQPI_kpath_calc = .FALSE.
        SlabQPI_kplane_calc = .FALSE.
        SlabSpintexture_calc = .FALSE.
        BulkSpintexture_calc = .FALSE.
        wanniercenter_calc = .FALSE.
        Z2_3D_calc = .FALSE.
        Chern_3D_calc = .FALSE.
        WeylChirality_calc = .FALSE.
        NLChirality_calc = .FALSE.
        BerryPhase_calc = .FALSE.
        BerryCurvature_calc = .FALSE.
        BerryCurvature_Cube_calc = .FALSE.
        BerryCurvature_slab_calc = .FALSE.
        MirrorChern_calc = .FALSE.
        Dos_calc = .FALSE.
        JDos_calc = .FALSE.
        EffectiveMass_calc = .FALSE.
        FindNodes_calc = .FALSE.
        LOTO_correction = .FALSE.
        Boltz_OHE_calc = .FALSE.
        AHC_calc = .FALSE.
        Hof_Butt_calc = .FALSE.
        LandauLevel_k_calc = .FALSE.
        LandauLevel_B_calc = .FALSE.
        LandauLevel_wavefunction_calc = .FALSE.
        OrbitalTexture_calc = .FALSE.
        OrbitalTexture_3D_calc = .FALSE.
        Fit_kp_calc = .FALSE.
        DMFT_MAG_calc = .FALSE.
        Symmetry_Import_calc = .FALSE.
        LanczosSeqDOS_calc = .FALSE.
        LandauLevel_kplane_calc = .FALSE.
        LandauLevel_k_dos_calc = .FALSE.
        LandauLevel_B_dos_calc = .FALSE.
        Translate_to_WS_calc = .FALSE.
        FermiLevel_calc = .FALSE.
    end

    % SYSTEM namelist
    properties
        % !> set system parameters by default
        Nslab = 10
        Nslab1 = 1
        Nslab2 = 1
        Numoccupied = 0
        Ntotch = 0
        SOC = 0
        SOC_in = 0
        E_FERMI = 0d0

        % !> By default magnetic field is zero
        Bx = 0d0
        By = 0d0
        Bz = 0d0

        Bmagnitude = 0d0
        Btheta = -99999d0
        Bphi = -99999d0
        surf_onsite = 0d0

        % !> By default, we don't add zeeman field
        Add_Zeeman_Field = .FALSE.

        % !> by default, g-factor is 2
        Effective_gfactor = 2d0
        Zeeman_energy_in_eV = 0d0

        % !> by default, Electric_field_in_eVpA=0
        Electric_field_in_eVpA = 0d0
        Symmetrical_Electric_field_in_eVpA = 0d0
        Inner_symmetrical_Electric_Field = .False.

        % !> by default, Vacuum_thickness_in_Angstrom= 12 Angstrom
        Vacuum_thickness_in_Angstrom = 12d0
    end

    % PARAMETERS namelist
    properties
        % !> set up parameters for calculation
        E_arc = 0.0d0
        Eta_Arc = 0.001d0
        OmegaNum = 100
        OmegaNum_unfold = 0
        OmegaMin = -1d0
        OmegaMax = 1d0
        Nk1 = 10
        Nk2 = 10
        Nk3 = 1
        NP = 2
        Gap_threshold = 0.01d0
        Tmin = 100. !in Kelvin
        Tmax = 100. !in Kelvin
        NumT = 1
        NBTau = 1
        BTauNum = 1
        Nslice_BTau_Max = 5000
        BTauMax = 0d0
        Rcut = 999999d0
        Magp = 1
        Magp_min = 0
        Magp_max = 0
        wcc_calc_tol = 0.08
        wcc_neighbour_tol = 0.3
        NumLCZVecs = 400
        NumRandomConfs = 1
        NumSelectedEigenVals = 0
        Beta = 100
        % !> by default, we only project on atoms for a given wave function
        projection_weight_mode = "NORMAL"
    end

    % LATTICE card
    properties
        WAN_NUM ;
        Rm =diag([1 1 1]);
        elementL;
        orbL = [];
        Gk
        quantumL
    end
    %  PROJECTORS card
    properties
        wanniercenter;
    end
    % MILLER_INDICES card
    properties
        MILLER_INDICES;
    end
    % SURFACE card
    properties
        Umatrix = eye(3); 
    end
    %  KPATH_BULK card
    properties
        %kpath information
    end
    %  KPATH_SLAB card
    properties
    end
    % KPLANE_BULK card
    properties
    end
    % KCUBE_BULK card
    properties
    end
    %  KPATH_BERRY card
    properties
    end
    % KPOINT_BULK card
    properties
    end  
    %  EFFECTIVE_MASS card
    properties
    end
    % KPOINTS_3D card
    properties
    end 
    % SURFACE_ATOMS card
    properties
    end 
    % SELECTED_ATOMS
    properties
    end
    
    methods

        function obj = untitled(inputArg1, inputArg2)
            %UNTITLED 构造此类的实例
            %   此处提供详细说明
            obj.Property1 = inputArg1 + inputArg2;
        end

        function outputArg = method1(obj, inputArg)
            %METHOD1 此处提供此方法的摘要
            %   此处提供详细说明
            outputArg = obj.Property1 + inputArg;
        end

    end

end
