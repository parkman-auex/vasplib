%% map_rule
% non s p d f sp sd pd sf spd spdf 
%  -1 0 1 2 3  4  5  6  7   8    9
%% slecet Projection
% usage=num_wan=write_pj(Projector_list)

function num_wan=write_pj(Projector_list)

[~,sites,Atom_name,~]=POSCAR_readin('digits',16);

sites_num = length(sites);
maprule_message=["n s p d f sp sd pd sf spd spdf";"-1 0 1 2 3  4  5  6  7   8    9"];
maprule= containers.Map({'n','s','p','d','f','sp','sd','pd','sf','spd','spdf'},{-1,0,1,2,3,4,5,6,7,8,9});
%% wannier_type_input
if nargin<1
    sites_num = length(sites);
    for i=1:sites_num
        fprintf("%d ",i);
        fprintf("%4s ",sites(i).name);
        fprintf("%f ",sites(i).rc1);
        fprintf("%f ",sites(i).rc2);
        fprintf("%f ",sites(i).rc3);
        fprintf("\n")
    end
    Projector_list=linspace(-1,-1,sites_num);
    fprintf("%s\n%s\n",maprule_message(1),maprule_message(2));

    for i=1:sites_num
            fprintf("%d ",i);
            fprintf("%4s ",sites(i).name);
            fprintf("%f ",sites(i).rc1);
            fprintf("%f ",sites(i).rc2);
            fprintf("%f ",sites(i).rc3);
            fprintf("\n");
        tempchar=input('please input projector(n s p d f sp sd pd sf spd spdf) of this ion: ','s');
        Projector_list(i)=maprule(char(tempchar));
    end
end
%% case.win projectors card gen

maprule2= containers.Map({0,1,2,3,4,5,6,7,8,9},{"l=0","l=1","l=2","l=3","l=0;l=1","l=0;l=2","l=1;l=2","l=0;l=3","l=0;l=1;l=2","l=0;l=1;l=2;l=3"});

    fid = fopen('wannier90projector_card','w');
    % function string
    head_string="begin projections ";
    end_string ="end projections ";
    fprintf(fid,"%s\n",head_string);
    for i=1:sites_num
        if Projector_list(i) ~= -1
        fprintf(fid,"%s","f= ");
        fprintf(fid,"%f, ",sites(i).rc1);
        fprintf(fid,"%f, ",sites(i).rc2);
        fprintf(fid,"%f:",sites(i).rc3);
        fprintf(fid,"%s ",maprule2(Projector_list(i)));
        fprintf(fid,"\n");
        end
    end
    fprintf(fid,"%s\n",end_string);
    fclose(fid);
%% wt.in projectors card gen

maprule3= containers.Map({-1,0,1,2,4,5,6,8},["","s", "pz px py","dz2 dxz dyz dx2-y2 dxy","s pz px py","s dz2 dxz dyz dx2-y2 dxy","pz px py dz2 dxz dyz dx2-y2 dxy","s pz px py dz2 dxz dyz dx2-y2 dxy"]);
maprule4= containers.Map({-1,0,1,2,3,4,5,6,7,8,9},{0,1,3,5,7,4,6,8,8,9,16});

    fid = fopen('wt_in_projector_card','w');
    % function string
    head_string="PROJECTORS ";
    fprintf(fid,"%s\n",head_string);
    num_wan=0;
    for i=1:sites_num
        tempnum=maprule4(Projector_list(i));
        num_wan = num_wan+tempnum;
       fprintf(fid,"%d ",tempnum);     
    end
    fprintf(fid,"\n");
    for i=1:sites_num
        if Projector_list(i) ~= -1
        fprintf(fid,"%s ",Atom_name(sites(i).nameseq));

        fprintf(fid,"%s ",maprule3(Projector_list(i)));
        fprintf(fid,"\n");
        end
    end

    fclose(fid);
%% wannier90_orbital.dat gen
digits(16);
maprule5= containers.Map({-1,0,1,2,4,5,6,8},{[""],["s"],["p"],["d"],["s","p"],["s","d"],["p","d"],["s","p","d"]});
fid = fopen('wannier90_orbital.dat','w');
% function string
head_string="0 # spinless is 0; spinful is 1 ";
fprintf(fid,"%s\n",head_string);
fprintf(fid,"%d\n",num_wan);
for i=1:sites_num
    if Projector_list(i) ~= -1
        wan_stringL = maprule5(Projector_list(i));
        for j = wan_stringL
            fprintf(fid,"%s,",j);
            fprintf(fid,"%19.16f, ",sites(i).rc1);
            fprintf(fid,"%19.16f, ",sites(i).rc2);
            fprintf(fid,"%19.16f",sites(i).rc3);
            fprintf(fid,"\n");
        end
    end
end
fclose(fid);
end










