%% example01 bandplot
% ****************************************
%             # vasp_nospin
% ****************************************
%
%% Label:
% * vasp
% * band
% * spinless
%
%% Description of the Script:
%
% plot a band stucture of standard vasp output
%
%% Usage: 
%
% * Modify or Understand this code
%
%% Input:
%  
% # Files: POSCAR KPOINTS 
%          EIGENVAL/BAND.dat
%          DOSCAR/Efermi
%
%% Output:
%
% # band-structure with good output
%
%% Note: 
%
%  Take advantage of the scope of application of the function.
%
%% Change log
%
% * Document Date: 2020/12/05
% * Creation Date: 2020/12/05
% * Last updated : 2020/12/05
%
%% Copyright
%
% * parkman
% * <parkman@buaa.edu.cn>
%
% Let's begin
%
%% make sure you have copy the POSCAR KPOINTS EIGENVAL DOSCAR in this directory
ls

%% for a default print,just type 
%   bandplot
% or 
%   bandplot()
bandplot;

%% for a user-defined bandplot
%
% we shoulf first load the EIGENCAR, 
%   EIGENCAR = EIGENVAL_read();
% which is a highly useful matric stucture
% with the format
%
% <html>
% <style type="text/css">
% .tg  {border-collapse:collapse;border-spacing:0;}
% .tg td{border-color:black;border-style:solid;border-width:1px;font-family:Arial, sans-serif;font-size:14px;
%   overflow:hidden;padding:10px 5px;word-break:normal;}
% .tg th{border-color:black;border-style:solid;border-width:1px;font-family:Arial, sans-serif;font-size:14px;
%   font-weight:normal;overflow:hidden;padding:10px 5px;word-break:normal;}
% .tg .tg-baqh{text-align:center;vertical-align:top}
% </style>
% <table class="tg" style="undefined;table-layout: fixed; width: 324px">
% <colgroup>
% <col style="width: 64px">
% <col style="width: 61px">
% <col style="width: 57px">
% <col style="width: 65px">
% <col style="width: 77px">
% </colgroup>
% <thead>
%   <tr>
%     <th class="tg-baqh"></th>
%     <th class="tg-baqh">kpoint 1</th>
%     <th class="tg-baqh">kpoint2</th>
%     <th class="tg-baqh">...</th>
%     <th class="tg-baqh">kpoint nk</th>
%   </tr>
% </thead>
% <tbody>
%   <tr>
%     <td class="tg-baqh">band 1</td>
%     <td class="tg-baqh">.</td>
%     <td class="tg-baqh">.</td>
%     <td class="tg-baqh"></td>
%     <td class="tg-baqh">.</td>
%   </tr>
%   <tr>
%     <td class="tg-baqh">band 2</td>
%     <td class="tg-baqh">.</td>
%     <td class="tg-baqh">.</td>
%     <td class="tg-baqh"></td>
%     <td class="tg-baqh">.</td>
%   </tr>
%   <tr>
%     <td class="tg-baqh">...</td>
%     <td class="tg-baqh"></td>
%     <td class="tg-baqh"></td>
%     <td class="tg-baqh"></td>
%     <td class="tg-baqh"></td>
%   </tr>
%   <tr>
%     <td class="tg-baqh">band nE</td>
%     <td class="tg-baqh">.</td>
%     <td class="tg-baqh">.</td>
%     <td class="tg-baqh"></td>
%     <td class="tg-baqh">.</td>
%   </tr>
% </tbody>
% </table>
% </html>
%
% <latex>
% \begin{table}[]
% \begin{tabular}{ccccc}
%         & kpoint 1 & kpoint2 & ... & kpoint nk \\
% band 1  & .        & .       &     & .         \\
% band 2  & .        & .       &     & .         \\
% ...     &          &         &     &           \\
% band nE & .        & .       &     & .        
% \end{tabular}
% \end{table}
% </latex>%
EIGENCAR = EIGENVAL_read();

%% bandplot parm
% The whole parm of bandplot is 
%   [fig,ax]=bandplot(EIGENCAR,Ecut,titlestring,color,klist_l,kpoints_l,...
%   kpoints_name,fontname,fig,ax)
% Here, we'd like check Ecut = [-6,6], with a title 'Band Stucture CuI w SOC'
% blue color, we just run and abtain a (title).eps file. You can open it with
% AI for publication
bandplot(EIGENCAR,[-6,6],'title','Band Stucture CuI w SOC','Color','b');

%% save plot data
% we can get handle and the plot data of a bandplot
[fig,ax,data_plot] = bandplot()