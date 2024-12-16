%% Run this script to add all the classes and functions in vasplib to your matlab search path
%% classes
classes_path = "src/classes";
% HR
%   HR                    - - a powerful TB tools
%
% HK
%   HK                    - a powerful kp tools
%
% Htrig                  
%   Htrig
%
% SYMMETRY
%   Oper                  - An operation of a continus or discontinus group.
%
% OTHER CLASS
%   Trig                  - 
%   Term                  - 
%   pauli_matric          - pauli matric 
%   gamma_matric          - pauli matric 
%   MaterialAccom         - 
%% plot functions
plot_funcs_path = "src/functions_plot";
%   bandplot              - bandplot
%   bandcompare           - General
%   bandplot_3d           - An example for the visualization of 3D-bandstructure for your own by using MATLAB
%   pbandplot             - pbandplot
%   plot_3D               - 
%   Plotdata              - Plotdata(X1, Y1)
%   plotPolyhedron        -  
%   POSCAR_plot           - 
%   BZplot                - 
%   dosplot               - DOS_plot
%   heatplot              - 
%   POSCAR_play           - note: POSCAR is Direct mode
%% API functions with the files in VASP formats and etc
api_funcs_path = "src/functions_api";
%  READ
%   DOSCAR_read           - DOSCAR_read
%   KPOINTS_read          - KPOINTS_read
%   EIGENVAL_read         - To Get EIGENCAR from file EIGENVAL
%   POSCAR_read           - note: POSCAR is Direct mode
%   POSCAR_readin         - POSCAR_readin
%   PROCAR_read           - PROCAR_read
%   PROCAR_select         - PROCAR_select
%   MAGMOM_read           - 
%   OUTCAR_read           - 
%  OUT
%   INCAR_gen             - 
%   cif_gen               - 
%   COLOR_one_gen         - 
%   COLORCAR_gen          - 
%   PARCHG_gen            - PARCHG gen
%   EIGENCAR_gen          - EIGENCAR_gen
%   EIGENCAR_gen_wire     - 
%   kp_gen                - kp_init (delete soon)
%   kpath_card_gen        - KPOINT_path card
%   sym_mat_gen           - 
%   SymInfor_gen          - [POSCAR_file,POSCAR_syminfor_file,ID] = SymInfor_gen(POSCAR_name,phonopy_run)
%   WEIGHTCAR_gen         - WEIGHTCAR gen for VASPKIT -> dat
%   MAGMOM_gen            - for
%   gnu_gen               - nargin
%   PBAND_DAT_gen         - PBAND_DAT_gen 
%   POSCAR_gen            - note: POSCAR is Direct mode
%   wt_in_gen             - 
%   wannier90_win_gen     - To generate Different card for vasp2wannier
%   Parity                - headline
%   Parity_eight          - 
%  EXTRACT
%   fermi_get             - 
%   Degeneracy_EIGENCAR   - 
%   find_node             - usage : [Energy_list,node_index_list,klist_s_list]=find_node(occupi_band_index,EIGENCAR,node_cut,klist_s)
%   findbands             - 
%   GetFermi              - Get fermi
%   Magion_detect         - UNTITLED2 此处显示有关此函数的摘要
%   supercell             - usage: [sites]=supercell(Ns)
%  SCRIPT
%   vasp_matlab_prework   - make sure present work and save the var
%   openmx_input          - Input: 0 = non_soc; 1 = soc (default);
%   init_v2w              - make sure present work and save the var
%   AutoJ                 - 
%   nanodisk_vasp_prework - make sure present work and save the var
%   enforce_2D_POSCAR     - 
%

%
% TOOLS
%   kmesh3D               - To generate k-mesh for DOSplot or 2D-Bandplot
%   kpathgen3D            - To generate klist for caculation and plot
%   DOSCAR_gen            - nargin
%   EIGENCAR2IMG          - nargin
%   nn_smart              - Caculate the nn for a primitive cell
%
% OTHERS 
%   Direct_kp             - direct kp model
%   HSVCAR_gen            - --------  nargin  --------
%   Gen_rmlist            - rm_list_gen
%   orbone2hsv            - --------  nargin  --------
%   POSCAR_rmRc           - Remove a specfic rc or a rc list in POSCAR
%   PROCAR_data_location  - PROCAR_data_location
%   selectbands           - 
%   soc_term_gen          - 
%   sorteig               - 
%   str2double_mat        - 
%   strcontain            - Retun 1 or 0 to show if there a string contains any string in string list
%   supercell_orb         - 
%   ta_tb_effect          - main-fourband model
%   term_effect_plot2     - set default para
%   to_red_sc             - 
%   write_pj              - map_rule
%   Ymlsym                - 
%   Htrig                 - 
%% data files
datas_path = "src/datas";
%%
main_path = pwd;
addpath(main_path+"/src");
addpath(main_path+"/"+classes_path);
addpath(main_path+"/"+plot_funcs_path);
addpath(main_path+"/"+api_funcs_path);
addpath(main_path+"/"+datas_path);
addpath(main_path+"/lib/export_fig")
savepath;