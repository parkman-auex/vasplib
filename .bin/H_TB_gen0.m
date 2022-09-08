%% To gen H_TB.m with Rnn 
% input: level-cut  (nn hopping max cut)
%        nn_store
%        sites
% usage: H=H_TB_gen(level_cut,nn_store,sites,pattern)
% output: 
%        H (H_up_triangle String Matrix (n,n), here n is the orbitals*atoms in one site)
%        H_TB.m ( DIY it for need)
% note : site=struct('seq',[],'inseq',[],'rc1',[],'rc2',[],'rc3',[],'name',[],'nameseq',[]);  
%        nn_t=struct('totseq',[],'orbit_on',[],'orbit_in',[],'Rlength',[],'Rc',[],'name',[],'nn_level',[]);

function H=H_TB_gen(level_cut,nn_store,sites,pattern)
    if nargin ==3
        pattern='bulk';
    end
%% file=
    fid = fopen('H_TB.m','w');
%% term_storage
    t_store=["ta","tb","tc","td","te","tf","tg","tf","ti","tj","tk","tl","tm"];
%% init
    [N_tot,N_orbit]=size(nn_store);
    
    
    if pattern == 'bulk'
%% bulk        
% function string
  head_string=" Hout=H_TB(k,Rm";
  for i=1:level_cut
      head_string=head_string+","+t_store(i);
  end
  head_string=head_string+")";
% H_string
    zerostring="";
    H=repmat(zerostring,[N_orbit N_orbit]);
    %H=zeros(N_orbit);% H_kj
%   on-site
    for j=1:N_orbit
        H(j,j)="E"+string(j);
    end
%   upper triangle (only nn -term)
    for j=1:N_orbit
        for i=1:N_tot
            k=nn_store(i,j).orbit_on;
            level=nn_store(i,j).nn_level;
            if level <= level_cut
                % note sign problem
                %
                Stringrc=string(nn_store(i,j).Rr);
                H(k,j)=H(k,j)+"-"+t_store(level)+"*exp(1i*dot(k,["+Stringrc(1)+" "+Stringrc(2)+" "+Stringrc(3)+"]*Rm))";
                
            end
        end
    end

%% file.out
   fprintf(fid,"%%%% This code is generated by H_TB_gen C. by parkman\n");  
   fprintf(fid,"%%  To describe a TB model with POSCAR \n");
   fprintf(fid,"%%  usage: %s\n",head_string);
   fprintf(fid,"%%  input: \n");
   fprintf(fid,"%%          k  (the k-vector in momont space)\n");
   fprintf(fid,"%%          Rm (the basis in real space)\n");
   fprintf(fid,"%%          nn_store (the k-vector in momont space)\n");
   for i=1:level_cut
    fprintf(fid,"%%          %s (the %d th nn hopping)\n",t_store(i),i);
   end
   fprintf(fid,"%% output: \n");
   fprintf(fid,"%%          H (H_nn ,where n is the orbitals*atoms in one site)\n");
   fprintf(fid,"%%            (It will be a Hermitian Matrix)\n");
   fprintf(fid,"%%  note: \n");
   %fprintf(fid,"%%       site=struct('seq',[],'inseq',[],'rc1',[],'rc2',[],'rc3',[],'name',[],'nameseq',[]);\n");
   %fprintf(fid,"%%       nn_t=struct('totseq',[],'orbit_on',[],'orbit_in',[],'Rlength',[],'Rc',[],'name',[],'nn_level',[]);\n");
   fprintf(fid,"%% \n");
   fprintf(fid,"%% \n");
   fprintf(fid,"%% \n");
   fprintf(fid,"function %s \n",head_string);
% end-head
% begin-init
   fprintf(fid,"\t%% init\n");
   fprintf(fid,"\n");
   fprintf(fid,"\t %%kx=k(1);\n");
   fprintf(fid,"\t %%ky=k(2);\n");
   fprintf(fid,"\t %%kz=k(3);\n");
   fprintf(fid,"\t H_up=zeros(%d);",N_orbit);
   fprintf(fid,"\t H_on=zeros(%d);",N_orbit);
   fprintf(fid,"\n");
   for j=1:N_orbit
        fprintf(fid,"\t E%d = 0;\t%% %s onsite\n",j,sites(j).name);     
   end
%% -TB-bulk
% H-on
%    fprintf(fid,"\t%% H-on\n");
% 
%    for k=1:N_orbit
%        fprintf(fid,"\t %% %s->%s\n",sites(k).name,sites(k).name); 
%        fprintf(fid,"\t H_on(%d,%d)=%s;\n",k,k,H(k,k));
%    end
% H-nn
   fprintf(fid,"\t%% H-nn\n");

   for k=1:N_orbit
       for j=1:N_orbit
           if H(k,j) ~= ""
                fprintf(fid,"\t %% %s->%s\n",sites(k).name,sites(j).name);  
                fprintf(fid,"\t H_up(%d,%d)=%s;\n",k,j,H(k,j));
           end
       end
   end
% Hout
   fprintf(fid,"\t%% Hout\n");
   fprintf(fid,"\t Hout=(H_up+H_up')/2;\n");
% end-init  
   fprintf(fid,"\nend  \n");
   
   
   
   
    elseif pattern =='nano'
        
        
        
%% nano        
% function string
  head_string=" Hout=H_TB(";
  for i=1:level_cut-1
      head_string=head_string+t_store(i)+",";
  end
  head_string=head_string+t_store(level_cut);
  head_string=head_string+")";
% H_string
    zerostring="";
    H=repmat(zerostring,[N_orbit N_orbit]);
    %H=zeros(N_orbit);% H_kj
%   upper triangle (only nn -term)
    for j=1:N_orbit
        for i=1:N_tot
            k=nn_store(i,j).orbit_on;
            level=nn_store(i,j).nn_level;
            if level <= level_cut
                % note sign problem
                %
                Stringrc=string(nn_store(i,j).Rr);
                H(k,j)=H(k,j)+"-"+t_store(level);
                
            end
        end
    end

%% file.out
   fprintf(fid,"%%%% This code is generated by H_TB_gen C. by parkman\n");  
   fprintf(fid,"%%  To describe a TB model with POSCAR \n");
   fprintf(fid,"%%  usage: %s\n",head_string);
   fprintf(fid,"%%  input: \n");
   %fprintf(fid,"%%          nn_store (the k-vector in momont space)\n");
   for i=1:level_cut
    fprintf(fid,"%%          %s (the %d th nn hopping)\n",t_store(i),i);
   end
   fprintf(fid,"%% output: \n");
   fprintf(fid,"%%          H (H_nn ,where n is the orbitals*atoms in one site)\n");
   fprintf(fid,"%%            (It will be a Hermitian Matrix)\n");
   fprintf(fid,"%%  note: \n");
   %fprintf(fid,"%%       site=struct('seq',[],'inseq',[],'rc1',[],'rc2',[],'rc3',[],'name',[],'nameseq',[]);\n");
   %fprintf(fid,"%%       nn_t=struct('totseq',[],'orbit_on',[],'orbit_in',[],'Rlength',[],'Rc',[],'name',[],'nn_level',[]);\n");
   fprintf(fid,"%% \n");
   fprintf(fid,"%% \n");
   fprintf(fid,"%% \n");
   fprintf(fid,"function %s \n",head_string);
% end-head
% begin-init
   fprintf(fid,"\t%% init\n");
   fprintf(fid,"\n");

   fprintf(fid,"\t H_up=zeros(%d);",N_orbit);

   fprintf(fid,"\n");

%% -TB-nano

   fprintf(fid,"\t%% H-nn\n");

   for k=1:N_orbit
       for j=1:N_orbit
           if H(k,j) ~= ""
                fprintf(fid,"\t %% %s->%s\n",sites(k).name,sites(j).name);  
                fprintf(fid,"\t H_up(%d,%d)=%s;\n",k,j,H(k,j));
           end
       end
   end
% Hout
   fprintf(fid,"\t%% Hout\n");
   fprintf(fid,"\t Hout=(H_up+H_up')/2;\n");
% end-init  
   fprintf(fid,"\nend  \n");
   
    end
end
