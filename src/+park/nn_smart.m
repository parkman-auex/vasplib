function [Atom_store_smart,nn_store_smart,Rnn,Rnn_map]=nn_smart(Rm,sites,search_range,Accuracy,Rlength_cut)
%% 
% Caculate the nn for a primitive cell
% input : Rm sites
% output : Atom_store nn_store,Rnn
% usage : [Atom_store,nn_store,Rnn]=nn(Rm,sites)
%         [Atom_store,nn_store,Rnn]=nn(Rm,sites)
% note : 
% Rm
    %a1 a1x a1y a1z
    %a2 a2x a2y a2z
    %a3 a3x a3y a3z
%     a1=Rm(1,:);
%     a2=Rm(2,:);
% note : site=struct('seq',[],'inseq',[],'rc1',[],'rc2',[],'rc3',[],'name',[],'nameseq',[]);  
%       Atom=struct('totseq',[],'elementseq',[],'Re',[],'seq_in',[],'Rc',[],'name',[]);
%       nn_t=struct('totseq',[],'orbit_on',[],'orbit_in',[],'Rlength',[],'Rc',[],'name',[],'nn_level',[]);
%       Rnn sort the Ri -> ti in a primitivecell 
%       the Accuracy should be paid attention


%% -------- nargin --------
if nargin <5
    Rlength_cut = 15;
end
if nargin <4
    Accuracy = 4;
end
if nargin <3
    search_range = [0 0 0];
end
%% -------- init --------
sites_num=size(sites,2);
search_rangex=search_range(1);
search_rangey=search_range(2);
search_rangez=search_range(3);   
Atom_smart_t = struct('R_fractional_from',[],'R_fractional_to',[],'R_fractional_diff',[],'seq_from',[],'seq_to',[],'handyname',[]);
nn_smart_t = struct('seq_from',[],'seq_to',[],'nn',[]);
Atom_store_smart = repmat(Atom_smart_t,[sites_num sites_num]);
nn_store_smart = repmat(nn_smart_t,[sites_num sites_num]);
Rnn_list = [];
%% -------- save ------------
for j  = 1:sites_num
    site2 = sites(j); % homecell
    for i = 1:sites_num
        site1 = sites(i);
        Atom_store_smart(i,j) = Atom_smart_t_gen(site1,site2);
        [Rnn_list_temp,nn_store_smart(i,j)] = ...
            nn_smart_t_gen(Atom_store_smart(i,j),Rm,search_rangex,search_rangey,search_rangez,Accuracy,Rlength_cut);
        Rnn_list = [Rnn_list;Rnn_list_temp];
    end
end
%% -------- caculate ------------
Rnn = sort(unique(Rnn_list,'row'));
%% -------- make map ------------
Rnn_map = containers.Map('KeyType','double','ValueType','double');
for i = 1:length(Rnn)
    Rnn_map(Rnn(i)) = i;
end
%%  give level
for j  = 1:sites_num
    for i = 1:sites_num
        for k = 1:length(nn_store_smart(i,j).nn)
            nn_store_smart(i,j).nn(k).nn_level = Rnn_map(nn_store_smart(i,j).nn(k).Rlength);
        end
    end
end
end

function Atom_smart_t = Atom_smart_t_gen(site1,site2) % othercell -> homecell
    Rc1 = [site1.rc1,site1.rc2,site1.rc3];
    Rc2 = [site2.rc1,site2.rc2,site2.rc3];
    Atom_smart_t.R_fractional_from = Rc1 ;
    Atom_smart_t.R_fractional_to   = Rc2; %fix
    Atom_smart_t.R_fractional_diff = -(Rc1 - Rc2);
    Atom_smart_t.seq_from = site1.seq;
    Atom_smart_t.seq_to = site2.seq;
    Atom_smart_t.handyname = strcat(site1.name,' -> ',site2.name);
end
function [Rnn_list,nn_smart_t] = nn_smart_t_gen(Atom_smart_t,Rm,search_rangex,search_rangey,search_rangez,Accuracy,Rlength_cut)
    nn_smart_t.seq_from = Atom_smart_t.seq_from;
    nn_smart_t.seq_to = Atom_smart_t.seq_to;
%     nn_smart_t.R_cartesian_to = Atom_smart_t.R_fractional_to*Rm;
%     nn_smart_t.R_cartesian_from = Atom_smart_t.R_fractional_from*Rm;
    count = 1;
    reducible_num=(2*search_rangex+1)*(2*search_rangey+1)*(2*search_rangez+1);
    Rnn_list = zeros(reducible_num,1);
    nn_t = struct('R_vector',[],'R_fractional_diff',[],'Rlength',[],'nn_level',[]);
    nn = repmat(nn_t,[reducible_num 1]);
    for Rf_a1=-search_rangex:search_rangex
        for Rf_a2=-search_rangey:search_rangey
            for Rf_a3=-search_rangez:search_rangez
                R_vector = [Rf_a1 Rf_a2 Rf_a3];
                Rlength = norm((R_vector + Atom_smart_t.R_fractional_diff)*Rm);
                Rlength = roundn(Rlength,-Accuracy);
                if  0 < Rlength && Rlength < Rlength_cut
                    nn(count,1).R_vector = R_vector;
                    nn(count,1).R_fractional_diff = Atom_smart_t.R_fractional_diff;
                    nn(count,1).Rlength = Rlength;
                    Rnn_list(count,:) = Rlength;
                    count = count +1;
                end
            end
        end
    end
    if count <= reducible_num
        Rnn_list(count:reducible_num,:) = [];
        nn(count:reducible_num,:) = [];
    end
    nn_smart_t.nn = nn;
end