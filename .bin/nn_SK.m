function [Atom_store,nn_store,Rnn]=nn_SK(Rm,sites,search_range, Accuracy)
%% 
% Caculate the nn_SK for a primitive cell
% input : Rm sites
% output : Atom_store nn_store,Rnn
% usage : [Atom_store,nn_store,Rnn]=nn(Rm,sites)
%         [Atom_store,nn_store,Rnn]=nn(Rm,sites,search_range)
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


%% init
if nargin <4
   Accuracy=4;
end
    sites_num=size(sites,2);

    Rlength_cut=10;
    if nargin <3
        search_range = [1 1 1];
        search_rangex=1;
        search_rangey=1;
        search_rangez=1;
    end
    search_rangex=search_range(1);
    search_rangey=search_range(2);
    search_rangez=search_range(3);

    reducible_num=(2*search_rangex+1)*(2*search_rangey+1)*(2*search_rangez+1);
    
    Atom_t=struct('totseq',[],'elementseq',[],'Re',[],'seq_in',[],'Rc',[],'name',[],'orb',[],'orb_sym',[]);
    nn_t=struct('totseq',[],'orbit_on',[],'orbit_in',[],'Rlength',[],'Rc',[],'Rr',[],'Re',[],'name',[],'nameseq',[],'nn_level',[],'ion_num',[]);
    Rlength_set=[];


%    Atom_in=struct('seq_in',[],'R_in',[],'name',[]);
    
%    element=repmat(Atom_in,[1 sites_num]);
     Atom_store=repmat(Atom_t,[1 sites_num*reducible_num]);
     nn_store=repmat(nn_t,[sites_num*reducible_num sites_num]);
   
% %% gen element---------------------------------------------------------------
%     for i=1:sites_num
%         element(i).seq_in=i;
%     end
%     clear i;
%     
%     element(1).name="C0";
%     element(2).name="C1";
%     
%     element(1).R_in=[0 0 0];
%     element(2).R_in=[0.25 0.25 0.25];
% 
% %------------------------------------------------------------------------------------
%% gen atom_store
    atomn=0;
    elementn=0;
    %
    for C_a1=-search_rangex:search_rangex
        for C_a2=-search_rangey:search_rangey
            for C_a3=-search_rangez:search_rangez
                 elementn=elementn+1;
                 Re=[C_a1 C_a2 C_a3];
                 for i=1:sites_num
                     %totseq
                      atomn=atomn+1;
                      Atom_store(atomn).totseq=atomn;
                     %siteseq
                      Atom_store(atomn).elementseq=elementn;
                     %Rc
                      Atom_store(atomn).Re=Re;
                     %seq_in
                      Atom_store(atomn).seq_in=i;
                     %Rr
                      R_in=[sites(i).rc1 sites(i).rc2 sites(i).rc3];
                      Atom_store(atomn).Rc=Re+R_in;
                     %name
                      Atom_store(atomn).name=num2str(Re)+sites(i).name;
                     % nn
                     
                     % 
                     Atom_store(atomn).orb = sites(i).orb;
                     Atom_store(atomn).orb_sym = sites(i).orb_sym;
                     
                     for j=1:sites_num
                         nn_store(atomn,j).totseq=atomn;
                         nn_store(atomn,j).orbit_in=j;
                         nn_store(atomn,j).orbit_on=Atom_store(atomn).seq_in;
                         R_r=(Re+R_in-[sites(j).rc1 sites(j).rc2 sites(j).rc3]);
                         R_c=[sites(j).rc1 sites(j).rc2 sites(j).rc3];
                         
                         Rij_cart =   R_r*Rm;
                         %[azimuth,elevation,Rlength] = cart2sph(Rij_cart(1),Rij_cart(2),Rij_cart(3));
                         Rlength=norm(Rij_cart);
                         Rlmn = Rij_cart/Rlength;
                         % 10^6 confine
                         Rlength=round(Rlength.*10^Accuracy)./10^Accuracy;
                         % Hopping_amplitude 
                         orb1 = sites(i).orb;
                         orb2 = sites(j).orb;
                         orb_sym1 = sites(i).orb_sym;
                         orb_sym2 = sites(j).orb_sym;                         
                         orbsym1_n =  subs_xyz(orb_sym1 ,Rlmn);
                         orbsym2_n =  subs_xyz(orb_sym2 ,Rlmn);
                         TBSK_hop =TBSK_hop_gen(orb1,orb2,orbsym1_n,orbsym2_n,orb_sym1,orb_sym2);
                        
                         
                         %
                         
                         nn_store(atomn,j).Rlength=Rlength;
                         %nn_store(atomn,j).Rlmn = Rlmn ;
                         nn_store(atomn,j).hop_pre =  TBSK_hop ;
                         nn_store(atomn,j).Rr=R_r;
                         nn_store(atomn,j).Re=Re;
                         nn_store(atomn,j).Rc=R_c;
                         nn_store(atomn,j).name=Atom_store(atomn).name+"<-"+sites(j).name;
                         nn_store(atomn,j).ion_num=sites(i).inseq;
                         nn_store(atomn,j).nameseq=sites(i).nameseq;
                         %nn_store(atomn,j).ion_type_num=e;
                         Rlength_set(atomn,j)=Rlength;
                     end

                 end
            end
        end
    end

%%  sort nn

    for i=1:sites_num
        %disp(sites(i).name);
        hi=tabulate(Rlength_set(:,i));
    end
    
    Rnnt=tabulate(Rlength_set(:));
    Rnn=Rnnt(2:end,1);
    
%%  give level
    for j=1:sites_num
        for i=1:atomn      
            signal=find(Rnn==nn_store(i,j).Rlength);
            nn_store(i,j).nn_level=signal;
            nn_store(i,j).hop = hop_gen(nn_store(i,j).hop_pre,signal);
        end
    end
    
   
end

function hop = hop_gen(hop_pre,nn_level)
    tempstr = string(simplify(hop_pre));
    tempstr = strrep(tempstr,'S',"S_"+string(nn_level));
    tempstr = strrep(tempstr,'P',"P_"+string(nn_level));
    hop = str2sym(tempstr);
end

function  TBSK_hop =TBSK_hop_gen(orb1,orb2,orbsym1_n,orbsym2_n,orbsym1,orbsym2)
    %disp([orbsym1_n,orbsym2_n]);
    TBSK_hop = orbsym1_n*orbsym2_n*sym("V"+orb1+orb2+'S')+...
                delta_orb(orb1,orb2)*(...
                delta_orb_sym(orbsym1,orbsym2)-orbsym1_n*orbsym2_n)*sym("V"+orb1+orb2+'P');    
end
 
function  orbsym_n = subs_xyz(orbsym,Rlmn)
    % subs ok   



    %disp(Rlmn);
    if orbsym == sym(1) 
        orbsym_n = 1;
    elseif orbsym == sym('x')
        orbsym_n = Rlmn(1); 
    elseif orbsym == sym('y')
        orbsym_n = Rlmn(2); 
    elseif orbsym == sym('z')
        orbsym_n = Rlmn(3);         
    end
    %disp(orbsym_n);
end
function out = delta_orb(orb1,orb2)
    if strcmp(orb1, orb2)
        out=1;
    else
        out =0;
    end
end

function out =delta_orb_sym(orb_sym1,orb_sym2)
    if orb_sym1 == orb_sym2
        out=1;
    else
        out =0;
    end
end

% bug on long range



