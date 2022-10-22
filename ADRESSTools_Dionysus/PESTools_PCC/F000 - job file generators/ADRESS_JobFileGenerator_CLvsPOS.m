function ADRESS_JobFileGenerator_CLvsPOS()
% ADRESS_JobFileGenerator_CLvsPOS()
%   This is a function that generates a text-file that is compatible to use
%   in the 'Restore Tab' button on 'SmartGUI' at the ADRESS beamline. It
%   allows the user to create a list of scans to be taken as a function of
%   position, where each scan can be referenced to a known Fermi-edge, if 
%   desired. Very useful for ARPES / XPS vs z-position data-sets, to
%   determine the homogeneity of a sample.
%% (A) Defining the input parameters
FILENAME        = char("ADRESS_JobFile_CLvsPOS_" + string(datestr(now,'yyyymmdd_HHMMSS')));
PHOTON_ENERGY   = 578;      % Photon energy of the XPS
REPEATS   	    = 1;        % Total number of repeat acquisitions (or Sweeps) for high statistics (only applies to sample measurements, not Ef)
%% (B) Set the SAMPLE measurement position
% - SAMPLE POSITION
x_samp      = -0.4:0.2:0.4; % can be a single value or vector
y_samp      ='-4.20';
z_samp      = -0.4:0.2:0.4; % can be a single value or vector
tilt_samp   ='0.00'; 
theta_samp  ='-11.00';
phi_samp    ='-233.00';
ADef_samp   = '0.0';
% - SAMPLE ACQUISITION SETTINGS
slit_samp   = '1';
dt_samp     = '0.1';
pol_samp 	= 'LV';
mode_samp 	= 'LAD';
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
name_ref    = 'AuEF';   % name_ref    = 'Au4d';
Eb_ref      = '-2';     % Eb_ref      = '-88:0.02:-80';
slit_ref    = '20';     % slit_ref    = '1';
dt_ref      = '30';     % dt_ref      = '0.1';
pol_ref 	= pol_samp;
mode_ref	= mode_samp;
Ep_ref      = Ep_samp;
%% (D) Initialising the JobFile with measurement information
POSITION = struct();
% -- aRPES / XPS info
ii=1;
POSITION(ii).name='As3d'; POSITION(ii).mode=mode_samp; POSITION(ii).Ep=Ep_samp;
POSITION(ii).pol=pol_samp; POSITION(ii).slit=slit_samp; POSITION(ii).dT=dt_samp;
POSITION(ii).Eb='-51:0.02:-39';
POSITION(ii).tilt=tilt_samp; POSITION(ii).theta=theta_samp; POSITION(ii).phi=phi_samp;
POSITION(ii).x=x_samp; POSITION(ii).y=y_samp; POSITION(ii).z=z_samp;
%% (E) Generate the JobFile Table
fid     = fopen([FILENAME '.txt'],'wt');
head    = '/*File|Comment|hv|Pol|Slit|Mode|Epass|Eb/Ek|Range|dt(s)|Sweeps|Theta|Tilt|Azimuth|X|Y|Z*/'; fprintf(fid,'%s\n',head);
index   = 1;
hv      = PHOTON_ENERGY;
% - Filing through each one of the photon energies
for ix = 1:length(x_samp)
    for iz = 1:length(z_samp)
        % -- Filing through core-levels
        for kk=1:length(POSITION)
            line=strcat(num2str(index, '%02.f'),'_',POSITION(kk).name,'_hv=',num2str(hv),'||',sprintf("ones(%i)*",REPEATS)+num2str(hv),'|',POSITION(kk).pol,'|',POSITION(kk).slit,'|',POSITION(kk).mode,'|',POSITION(kk).Ep,'|','Eb','|',POSITION(kk).Eb,'|',ADef_ref,'|',POSITION(kk).dT,'|1|',POSITION(kk).theta,'|',POSITION(kk).tilt,'|',POSITION(kk).phi,'|',num2str(POSITION(kk).x(ix)),'|',POSITION(kk).y,'|',num2str(POSITION(kk).z(iz)),'|');
            fprintf(fid,'%s\n',line);
            index=index+1;
        end
        % -- If necessary, acquire EF reference
        if use_ref == 1
            if length(string(x_ref)) == 1 && length(string(z_ref)) == 1
                line=strcat(num2str(index, '%02.f'),'_',name_ref,'_hv=',num2str(hv),'||',num2str(hv),'|',pol_ref,'|',slit_ref,'|',mode_ref,'|',Ep_ref,'|','Eb','|',Eb_ref,'|',ADef_ref,'|',dt_ref,'|1|',theta_ref,'|',tilt_ref,'|',phi_ref,'|',num2str(x_ref),'|',y_ref,'|',num2str(z_ref),'|');
            elseif length(string(x_ref)) > 1 && length(string(z_ref)) == 1
                line=strcat(num2str(index, '%02.f'),'_',name_ref,'_hv=',num2str(hv),'||',num2str(hv),'|',pol_ref,'|',slit_ref,'|',mode_ref,'|',Ep_ref,'|','Eb','|',Eb_ref,'|',ADef_ref,'|',dt_ref,'|1|',theta_ref,'|',tilt_ref,'|',phi_ref,'|',num2str(x_ref(ix)),'|',y_ref,'|',num2str(z_ref),'|');
            elseif length(string(x_ref)) == 1 && length(string(z_ref)) > 1
                line=strcat(num2str(index, '%02.f'),'_',name_ref,'_hv=',num2str(hv),'||',num2str(hv),'|',pol_ref,'|',slit_ref,'|',mode_ref,'|',Ep_ref,'|','Eb','|',Eb_ref,'|',ADef_ref,'|',dt_ref,'|1|',theta_ref,'|',tilt_ref,'|',phi_ref,'|',num2str(x_ref),'|',y_ref,'|',num2str(z_ref(iz)),'|');
            else
                line=strcat(num2str(index, '%02.f'),'_',name_ref,'_hv=',num2str(hv),'||',num2str(hv),'|',pol_ref,'|',slit_ref,'|',mode_ref,'|',Ep_ref,'|','Eb','|',Eb_ref,'|',ADef_ref,'|',dt_ref,'|1|',theta_ref,'|',tilt_ref,'|',phi_ref,'|',num2str(x_ref(ix)),'|',y_ref,'|',num2str(z_ref(iz)),'|');
            end
            fprintf(fid,'%s\n',line);
            index=index+1;
        end
    end
end
fclose(fid);
end