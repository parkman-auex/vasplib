function [freq, re_diag, im_diag] = OPTICS_read(opts)
%%
% Needed: INCAR OUTCAR, vasprun.xml
%
% ref: DOi: 10.1038/srep12285, Sec. Methods
%
% For metals, the frequency dependent dielectric function consists of 
% interband and Drude-like intraband contributions.
%
% The interband part is calculated with LOPTICS = T (IPA) or ALGO=Chi(RPA). 
% For IPA, the "optics.sh" script given by vaspwiki only collects the
% "density-density" response functions. Here we also collect the
% "current-current" response functions and take it in postprocessing.
% 
% The intraband part is calculated inside this matlab function,
% with the plasmon frequency tensors given by VASP in OUTCAR, 
% and the relaxation time set in INCAR using the tag RTIME (in unit of fs).
%
% WARNING: We are not sure how the tag RTIME takes effect in VASP
% calculations, because it can NOT be found in the category "INCAR tag" 
% in vaspwiki website. But a large RTIME does cause divergence near
% zero-point.
% And The default value of RTIME (-0.1 fs) is much smaller that common metals. 
% For robustness we recommend you to manually set RTIME in INCAR,
% and this function will read it from INCAR.
%%
arguments
    opts.dfhandle double = 0
    opts.mode {mustBeMember(opts.mode,{'ipa','rpa'})} = 'ipa'    
    opts.RTIME double = 0
    opts.count {mustBeMember(opts.count,{'all','inter','intra'})} = 'all'
    opts.ax handle = handle([])
    opts.title = ''
    opts.Xlim double = [0 6]
end
%%
%% draw all the blocks of frequency dependent data from vasprun.xml
% system("cp vasprun.xml vasprun.xml.bk");
[~, end_of_blocks] = system("grep '/imag' vasprun.xml");
endlines = splitlines(end_of_blocks);
nblock = length(endlines) - 1;
for ii = 1:nblock
    system("sed -i '0,/<imag>/{ s/<imag>/<I"+ii+">/;}' vasprun.xml");
    system("sed -i '0,/<real>/{ s/<real>/<R"+ii+">/;}' vasprun.xml");
    
    system("awk 'BEGIN{i=1} /I"+ii+"/,"+...
                    "/\/imag/ "+...
                    "{a[i]=$2 ; b[i]=$3 ; c[i]=$4; d[i]=$5 ; e[i]=$6 ; f[i]=$7; g[i]=$8; i=i+1} "+...
        "END{for (j=12;j<i-3;j++) print a[j],b[j],c[j],d[j],e[j],f[j],g[j]}' vasprun.xml > IMAG"+ii+".dat");

    system("awk 'BEGIN{i=1} /R"+ii+"/,"+...
                    "/\/real/ "+...
                    "{a[i]=$2 ; b[i]=$3 ; c[i]=$4; d[i]=$5 ; e[i]=$6 ; f[i]=$7; g[i]=$8; i=i+1} "+...
        "END{for (j=12;j<i-3;j++) print a[j],b[j],c[j],d[j],e[j],f[j],g[j]}' vasprun.xml > REAL"+ii+".dat");
end
%% interband
% switch opts.mode
%     case 'ipa'
%         disp("mode = IPA (indpendent particle approximation)")
%         
%         re_inter = readmatrix("REAL2.dat");
%         im_inter = readmatrix("IMAG2.dat");
%     case 'rpa'
%         disp("mode = RPA (random phase approximation)")
%         disp("Mention: LRPA=.FALSE. do not change the form of OUTCAR,");
%         disp("         so you can still use OPTICS_read('mode','rpa') to read it.");
%         
%         re_inter = readmatrix("REAL3.dat");
%         im_inter = readmatrix("IMAG3.dat");
% end
% % 3:end to avoid *** values
% freq=re_inter(3:end,1); % hbar * freq
% 
% re_inter_diag=re_inter(3:end,2:4);
% im_inter_diag=im_inter(3:end,2:4);
vaspxml = readstruct('vasprun.xml');
df = vaspxml.calculation.dielectricfunction;
if opts.dfhandle == 0
    switch length(df)
        case 2
            dfhandle = 1; % "density"
            % dfhandle = 2; % "current"
        case 4
            dfhandle = 1; % ipa
            % dfhandle = 3; % rpa
    end
else
    dfhandle = opts.dfhandle;
end
re_inter_tot = multi_split(df(dfhandle).real.array.set.r'); 
im_inter_tot = multi_split(df(dfhandle).imag.array.set.r'); 

freq          = re_inter_tot(3:end,1); % hbar * freq
re_inter_diag = re_inter_tot(3:end,2:4);
im_inter_diag = im_inter_tot(3:end,2:4);
%% intraband
if opts.RTIME == 0
    [~,tmp_char] = system("grep RTIME INCAR");
    if isempty(tmp_char)
        tau = -0.1; % fs, vasp default
    else
        tmp_str = split(string(tmp_char(1:end-1)));
        tmp_db  = split(tmp_str(end), "=");
        tau = str2double(tmp_db(end));
    end
    disp("The relaxation time is "+tau+" fs");
else
    tau = opts.RTIME;
    disp("The relaxation time is "+tau+" fs (manually set)");
end
tau = tau*1e-15;

h_eV_s = 4.1357e-15; % eV.s
gamma = h_eV_s * 1/tau;

plasma2 = zeros(1,3);
[~,tmp_char] = system("grep -A4 'intraband' OUTCAR");
if isempty(tmp_char)
    disp("Do not find 'plasma frequency squared' in OUTCAR! Auto changing opts.count");
    opts.count = 'inter';
else
    tmp_str = split(string(tmp_char));
    plasma2(1) = str2double(tmp_str(10));
    plasma2(2) = str2double(tmp_str(14));
    plasma2(3) = str2double(tmp_str(18));

    para_freq = (freq.^2 + gamma^2).^-1;
    re_intra_diag = 1 - para_freq * plasma2;
    im_intra_diag = (gamma./freq.*para_freq) * plasma2;
end
%% choose
switch opts.count
    case 'all'
        re_diag = re_inter_diag + re_intra_diag;
        im_diag = im_inter_diag + im_intra_diag;
    case 'inter'
        re_diag = re_inter_diag;
        im_diag = im_inter_diag;        
    case 'intra'
        re_diag = re_intra_diag;
        im_diag = im_intra_diag;        
end    
%% plot
if isempty(opts.ax)
    figure();
    ax = gca;
else
    if isvalid(opts.ax)
        ax = opts.ax;
    else
        figure();
        ax = gca;
    end
end

hold on
plot(ax, opts.Xlim, [0 0],':k','LineWidth',1.5,'HandleVisibility','off');

plot(ax, freq, re_diag(:,1),'-','LineWidth',1.5);
plot(ax, freq, re_diag(:,2),'-','LineWidth',1.5);
plot(ax, freq, re_diag(:,3),'-','LineWidth',1.5);

plot(ax, freq, im_diag(:,1),'--','LineWidth',1.5);
plot(ax, freq, im_diag(:,2),'--','LineWidth',1.5);
plot(ax, freq, im_diag(:,3),'--','LineWidth',1.5);
hold off

legend(ax, 'RE \epsilon_x','RE \epsilon_y','RE \epsilon_z',...
    'IM \epsilon_x','IM \epsilon_y','IM \epsilon_z');
xlabel(ax, 'Energy/eV');
ylabel(ax, 'Dielectric Constant');
title(ax, opts.title);

xlim(ax, opts.Xlim);
ylim(ax, [-20 20]);
set(ax, 'FontSize',21);
set(ax, 'LineWidth',1.5);

saveas(gcf,"optics_"+opts.mode+"_"+opts.count+".pdf");
end

function rdouble =  multi_split(str_array)
rstr = split(str_array," ");
rdouble = double(rstr);
end