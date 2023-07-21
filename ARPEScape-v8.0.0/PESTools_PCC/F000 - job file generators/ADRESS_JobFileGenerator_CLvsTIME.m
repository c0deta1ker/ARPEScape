function ADRESS_JobFileGenerator_CLvsTIME()
% ADRESS_JobFileGenerator_CLvsTIME()
%   This is a function that generates a text-file that is compatible to use
%   in the 'Restore Tab' button on 'SmartGUI' at the ADRESS beamline. It
%   allows the user to create a list of scans to be taken as a function of
%   time, where each scan can be referenced to a known Fermi-edge, if 
%   desired. Very useful for ARPES / XPS vs time data sets, to build high
%   statistics for analysis.
%% (A) Defining the input parameters
FILENAME        = char("ADRESS_JobFile_CLvsTIME_" + string(datestr(now,'yyyymmdd_HHMMSS')));
CYCLES          = 12;       % Total number of periodic cycles to run
PHOTON_ENERGY   = 578;      % Photon energy of the XPS
REPEATS   	    = 1;        % Total number of repeat acquisitions (or Sweeps) for high statistics (only applies to sample measurements, not Ef)
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
%% (D) Initialising the JobFile with the first entries
POSITION = struct();
% - Regular Energy Jobs to be Cycled
% -- In4d
ii=1;
POSITION(ii).name='In4d'; POSITION(ii).mode=mode_samp; POSITION(ii).Ep=Ep_samp;
POSITION(ii).pol=pol_samp; POSITION(ii).slit=slit_samp; POSITION(ii).dT=dt_samp;
POSITION(ii).Eb='-26:0.02:-15';
POSITION(ii).tilt=tilt_samp; POSITION(ii).theta=theta_samp; POSITION(ii).phi=phi_samp;
POSITION(ii).x=x_samp; POSITION(ii).y=y_samp; POSITION(ii).z=z_samp;
% -- As3d
ii=ii+1;
POSITION(ii).name='As3d'; POSITION(ii).mode=mode_samp; POSITION(ii).Ep=Ep_samp;
POSITION(ii).pol=pol_samp; POSITION(ii).slit=slit_samp; POSITION(ii).dT=dt_samp;
POSITION(ii).Eb='-51:0.02:-39';
POSITION(ii).tilt=tilt_samp; POSITION(ii).theta=theta_samp; POSITION(ii).phi=phi_samp;
POSITION(ii).x=x_samp; POSITION(ii).y=y_samp; POSITION(ii).z=z_samp;
% -- In3d
ii=ii+1;
POSITION(ii).name='In3d';POSITION(ii).mode=mode_samp;POSITION(ii).Ep=Ep_samp;
POSITION(ii).pol=pol_samp;POSITION(ii).slit=slit_samp;POSITION(ii).dT=dt_samp;
POSITION(ii).Eb='-457:0.02:-450';
POSITION(ii).tilt=tilt_samp;POSITION(ii).theta=theta_samp;POSITION(ii).phi=phi_samp;
POSITION(ii).x=x_samp;POSITION(ii).y=y_samp;POSITION(ii).z=z_samp;
%% (E) Generate the JobFile Table
fid     = fopen([FILENAME '.txt'],'wt');
head    = '/*File|Comment|hv|Pol|Slit|Mode|Epass|Eb/Ek|Range|dt(s)|Sweeps|Theta|Tilt|Azimuth|X|Y|Z*/'; fprintf(fid,'%s\n',head);
index   = 1;
Energies   	= PHOTON_ENERGY .* ones(1,CYCLES);
% - Filing through each one of the photon energies
for hv = Energies 
    % -- If necessary, acquire EF reference first
    if use_ref == 1
        line=strcat(num2str(index, '%02.f'),'_',name_ref,'_hv=',num2str(hv),'||',num2str(hv),'|',pol_ref,'|',slit_ref,'|',mode_ref,'|',Ep_ref,'|','Eb','|',Eb_ref,'|',ADef_ref,'|',dt_ref,'|1|',theta_ref,'|',tilt_ref,'|',phi_ref,'|',x_ref,'|',y_ref,'|',z_ref,'|');
        fprintf(fid,'%s\n',line);
        index=index+1;
    end
    % -- Filing through core-levels
    for kk=1:length(POSITION)
        line=strcat(num2str(index, '%02.f'),'_',POSITION(kk).name,'_hv=',num2str(hv),'||',sprintf("ones(%i)*",REPEATS)+num2str(hv),'|',POSITION(kk).pol,'|',POSITION(kk).slit,'|',POSITION(kk).mode,'|',POSITION(kk).Ep,'|','Eb','|',POSITION(kk).Eb,'|',ADef_ref,'|',POSITION(kk).dT,'|1|',POSITION(kk).theta,'|',POSITION(kk).tilt,'|',POSITION(kk).phi,'|',POSITION(kk).x,'|',POSITION(kk).y,'|',POSITION(kk).z,'|');
        fprintf(fid,'%s\n',line);
        index=index+1;
    end
end
fclose(fid);
end