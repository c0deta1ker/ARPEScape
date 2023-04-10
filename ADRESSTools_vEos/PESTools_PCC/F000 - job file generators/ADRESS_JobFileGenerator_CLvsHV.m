function ADRESS_JobFileGenerator_CLvsHV()
% ADRESS_JobFileGenerator_CLvsHV()
%   This is a function that generates a text-file that is compatible to use
%   in the 'Restore Tab' button on 'SmartGUI' at the ADRESS beamline. It
%   allows the user to create a list of scans to be taken as a function of
%   photon energy, where each scan can be referenced to a known Fermi-edge,
%   if desired. Very useful for XPS vs photon energy data sets.
%% (A) Defining the input parameters
FILENAME        = char("ADRESS_JobFile_CLvsHV_" + string(datestr(now,'yyyymmdd_HHMMSS')));
REPEATS   	    = 1;        % Total number of repeat acquisitions (or Sweeps) for high statistics (only applies to sample measurements, not Ef)
% -- Photon Energy cell array to file through
PHOTON_ENERGY{1} = [370,450,578,650,786,850,950,1050];  % Photon energy list #1
PHOTON_ENERGY{2} = [650,786,850,950,1050];              % Photon energy list #2
PHOTON_ENERGY{3} = [1050];                              % Photon energy list #3
%% (B) Set the SAMPLE measurement position
% - SAMPLE POSITION
x_samp      ='-3.60';
y_samp      ='-4.20';
z_samp      ='5.098';
tilt_samp   ='0.00'; 
theta_samp  ='-11.00';
phi_samp    ='-233.00';
ADef_samp   = '0.0';
% - SAMPLE ACQUISITION SETTINGS
slit_samp   = '1';
dt_samp     = '0.1';
pol_samp 	= 'LV';
mode_samp 	= 'MAM';
Ep_samp   	= '100'; 
%% (C) Set the Fermi-level REFERENCE position
use_ref     = 1;	% Set whether you want to use a reference; 1 = True, 0 = False
% - REFERENCE POSITION
x_ref       = x_samp;
y_ref       = y_samp;
z_ref       = z_samp;
tilt_ref    = tilt_samp;
theta_ref   = theta_samp;
phi_ref     = phi_samp;
ADef_ref    = ADef_samp;
% - REFERENCE ACQUISITION SETTINGS
name_ref    = 'EF';   % name_ref    = 'AuEF';   % name_ref    = 'Au4d';
Eb_ref      = '-2';     % Eb_ref      = '-88:0.02:-80';
slit_ref    = '20';     % slit_ref    = '1';
dt_ref      = '20';     % dt_ref      = '0.1';
pol_ref 	= pol_samp;
mode_ref	= mode_samp;
Ep_ref      = Ep_samp;
%% (D) Initialising the JobFile with the entries
POSITION = {};
%% - PHOTON ENERGIES #1
POSITION{1} = struct();
% -- In4d
ii=1;
POSITION{1}(ii).name='In4d'; POSITION{1}(ii).mode=mode_samp; POSITION{1}(ii).Ep=Ep_samp;
POSITION{1}(ii).pol=pol_samp; POSITION{1}(ii).slit=slit_samp; POSITION{1}(ii).dT=dt_samp;
POSITION{1}(ii).Eb='-33:0.02:-13';
POSITION{1}(ii).tilt=tilt_samp; POSITION{1}(ii).theta=theta_samp; POSITION{1}(ii).phi=phi_samp;
POSITION{1}(ii).x=x_samp; POSITION{1}(ii).y=y_samp; POSITION{1}(ii).z=z_samp;
% -- As3d
ii=ii+1;
POSITION{1}(ii).name='As3d'; POSITION{1}(ii).mode=mode_samp; POSITION{1}(ii).Ep=Ep_samp;
POSITION{1}(ii).pol=pol_samp; POSITION{1}(ii).slit=slit_samp; POSITION{1}(ii).dT=dt_samp;
POSITION{1}(ii).Eb='-51:0.02:-37';
POSITION{1}(ii).tilt=tilt_samp; POSITION{1}(ii).theta=theta_samp; POSITION{1}(ii).phi=phi_samp;
POSITION{1}(ii).x=x_samp; POSITION{1}(ii).y=y_samp; POSITION{1}(ii).z=z_samp;
% -- Al2p
ii=ii+1;
POSITION{1}(ii).name='Al2p'; POSITION{1}(ii).mode=mode_samp; POSITION{1}(ii).Ep=Ep_samp;
POSITION{1}(ii).pol=pol_samp; POSITION{1}(ii).slit=slit_samp; POSITION{1}(ii).dT=dt_samp;
POSITION{1}(ii).Eb='-82:0.02:-71';
POSITION{1}(ii).tilt=tilt_samp; POSITION{1}(ii).theta=theta_samp; POSITION{1}(ii).phi=phi_samp;
POSITION{1}(ii).x=x_samp; POSITION{1}(ii).y=y_samp; POSITION{1}(ii).z=z_samp;
% -- Pb4f
ii=ii+1;
POSITION{1}(ii).name='Pb4f'; POSITION{1}(ii).mode=mode_samp; POSITION{1}(ii).Ep=Ep_samp;
POSITION{1}(ii).pol=pol_samp; POSITION{1}(ii).slit=slit_samp; POSITION{1}(ii).dT=dt_samp;
POSITION{1}(ii).Eb='-147:0.02:-134';
POSITION{1}(ii).tilt=tilt_samp; POSITION{1}(ii).theta=theta_samp; POSITION{1}(ii).phi=phi_samp;
POSITION{1}(ii).x=x_samp; POSITION{1}(ii).y=y_samp; POSITION{1}(ii).z=z_samp;
% -- C1s
ii=ii+1;
POSITION{1}(ii).name='C1s'; POSITION{1}(ii).mode=mode_samp; POSITION{1}(ii).Ep=Ep_samp;
POSITION{1}(ii).pol=pol_samp; POSITION{1}(ii).slit=slit_samp; POSITION{1}(ii).dT=dt_samp;
POSITION{1}(ii).Eb='-290:0.02:-280';
POSITION{1}(ii).tilt=tilt_samp; POSITION{1}(ii).theta=theta_samp; POSITION{1}(ii).phi=phi_samp;
POSITION{1}(ii).x=x_samp; POSITION{1}(ii).y=y_samp; POSITION{1}(ii).z=z_samp;
%% - PHOTON ENERGIES #2
POSITION{2} = struct();
% -- In3d
kk=1;
POSITION{2}(kk).name='In3d';POSITION{2}(kk).mode=mode_samp;POSITION{2}(kk).Ep=Ep_samp;
POSITION{2}(kk).pol=pol_samp;POSITION{2}(kk).slit=slit_samp;POSITION{2}(kk).dT=dt_samp;
POSITION{2}(kk).Eb='-449:0.02:-440';
POSITION{2}(kk).tilt=tilt_samp;POSITION{2}(kk).theta=theta_samp;POSITION{2}(kk).phi=phi_samp;
POSITION{2}(kk).x=x_samp;POSITION{2}(kk).y=y_samp;POSITION{2}(kk).z=z_samp;
% -- O1s
kk=kk+1;
POSITION{2}(kk).name='O1s';POSITION{2}(kk).mode=mode_samp;POSITION{2}(kk).Ep=Ep_samp;
POSITION{2}(kk).pol=pol_samp;POSITION{2}(kk).slit=slit_samp;POSITION{2}(kk).dT=dt_samp;
POSITION{2}(kk).Eb='-535:0.02:-525';
POSITION{2}(kk).tilt=tilt_samp;POSITION{2}(kk).theta=theta_samp;POSITION{2}(kk).phi=phi_samp;
POSITION{2}(kk).x=x_samp;POSITION{2}(kk).y=y_samp;POSITION{2}(kk).z=z_samp;
%% - PHOTON ENERGIES #3
POSITION{3} = struct();
% -- Survey Scan
hh=1;
POSITION{3}(hh).name='Survey';POSITION{3}(hh).mode=mode_samp;POSITION{3}(hh).Ep=Ep_samp;
POSITION{3}(hh).pol=pol_samp;POSITION{3}(hh).slit=slit_samp;POSITION{3}(hh).dT='0.1';
POSITION{3}(hh).Eb='-555:0.2:5';
POSITION{3}(hh).tilt=tilt_samp;POSITION{3}(hh).theta=theta_samp;POSITION{3}(hh).phi=phi_samp;
POSITION{3}(hh).x=x_samp;POSITION{3}(hh).y=y_samp;POSITION{3}(hh).z=z_samp;

%% (E) Generate the JobFile Table
fid     = fopen([FILENAME '.txt'],'wt');
head    = '/*File|Comment|hv|Pol|Slit|Mode|Epass|Eb/Ek|Range|dt(s)|Sweeps|Theta|Tilt|Azimuth|X|Y|Z*/'; fprintf(fid,'%s\n',head);
index   = 1;
Energies = sort(unique(cell2mat(PHOTON_ENERGY)));
% - Filing through each one of the photon energies
for hv = Energies 
    % -- If necessary, acquire EF reference at the start
    if use_ref == 1
        line=strcat(num2str(index, '%02.f'),'_',name_ref,'_hv=',num2str(hv),'||',num2str(hv),'|',pol_ref,'|',slit_ref,'|',mode_ref,'|',Ep_ref,'|','Eb','|',Eb_ref,'|',ADef_ref,'|',dt_ref,'|1|',theta_ref,'|',tilt_ref,'|',phi_ref,'|',x_ref,'|',y_ref,'|',z_ref,'|');
        fprintf(fid,'%s\n',line);
        index=index+1;
    end
    % -- Filing through all the core-levels
    for n = 1:length(PHOTON_ENERGY)
        if ismember(hv,PHOTON_ENERGY{n})
            for kk=1:length(POSITION{n})
                line=strcat(num2str(index, '%02.f'),'_',POSITION{n}(kk).name,'_hv=',num2str(hv),'||',sprintf("ones(%i)*",REPEATS)+num2str(hv),'|',POSITION{n}(kk).pol,'|',POSITION{n}(kk).slit,'|',POSITION{n}(kk).mode,'|',POSITION{n}(kk).Ep,'|','Eb','|',POSITION{n}(kk).Eb,'|',ADef_samp,'|',POSITION{n}(kk).dT,'|1|',POSITION{n}(kk).theta,'|',POSITION{n}(kk).tilt,'|',POSITION{n}(kk).phi,'|',POSITION{n}(kk).x,'|',POSITION{n}(kk).y,'|',POSITION{n}(kk).z,'|');
                fprintf(fid,'%s\n',line);
                index=index+1;
            end
        end
    end
    % -- If necessary, acquire EF reference at the end
    if use_ref == 1
        line=strcat(num2str(index, '%02.f'),'_',name_ref,'_hv=',num2str(hv),'||',num2str(hv),'|',pol_ref,'|',slit_ref,'|',mode_ref,'|',Ep_ref,'|','Eb','|',Eb_ref,'|',ADef_ref,'|',dt_ref,'|1|',theta_ref,'|',tilt_ref,'|',phi_ref,'|',x_ref,'|',y_ref,'|',z_ref,'|');
        fprintf(fid,'%s\n',line);
        index=index+1;
    end
end
fclose(fid);
end