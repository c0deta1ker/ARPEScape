function JobFileGenerator_CLvsHV()
%% (A) Defining the consistent parameters
% (A1) Set the output file name to be saved as
sample_name = "Brugg";
filename    = char("JobFile_CLvsHV_" + sample_name + "_"+ string(datestr(now,'yyyymmdd_HHMMSS')));
% (A2) Set the list of photon energies to be used
Energies   	= [370,450,578,650,786,850,950,1045];
sEnergy1  	= Energies(3:end);      % For In3d, O1s,    > 550 eV energies only
sEnergy2 	= Energies(end);        % For XPS survey,   > 1000 eV only
% (A3) Define the total number of repeat acquisitions for high statistics (only applies to sample measurements, not Ef)
Repeats   	= 1;
sRepeats1  	= 1;
sRepeats2 	= 1;
% (A4) Define the list of relevant core-levels
cl_fm_In4d  = '-17.00';     cl_sm_In4d  = '-25:0.02:-10';
cl_fm_As3d  = '-41.00';     cl_sm_As3d  = '-51:0.02:-37';    % to include the oxide, use: '-51:0.02:-37';
cl_fm_Al2p  = '-74.00';     cl_sm_Al2p  = '-82:0.02:-71';
cl_fm_Au4d  = '-86.00';     cl_sm_Au4d  = '-91:0.02:-82';
cl_fm_C1s   = '-285.00';    cl_sm_C1s   = '-290:0.02:-280';
cl_fm_In3d  = '-444.00';    cl_sm_In3d  = '-449:0.02:-440';  % to include other SO component, use: '-456:0.05:-440';
cl_fm_O1s   = '-533.00';    cl_sm_O1s   = '-538:0.02:-529';
% (A5) Define the list of other useful binding energy ranges
vb_eb       = '-2.00';
cl_survey   = '-555:0.05:5';

%% (B) Set the SAMPLE measurement position
% - POSITION
x_samp      ='-1.704';
y_samp      ='-7.251';
z_samp      ='0.100';
tilt_samp   ='0.0'; 
theta_samp  ='0.0';
phi_samp    ='-57';
% - ACQUISITION SETTINGS
slitXPS_samp	= '1';      slitVB_samp		= '20'; 
dtXPS_samp    	= '0.2';    dtVB_samp    	= '60';
pol_samp 	= 'LV';
mode_samp 	= 'MAM';
Ep_samp   	= '100'; 

%% (C) Set the Fermi-level REFERENCE position
use_ref     = 1;	% Set whether you want to use a reference; 1 = True, 0 = False
% - POSITION
x_ref       = x_samp;
y_ref       = y_samp;
z_ref       = z_samp;
tilt_ref    = tilt_samp;
theta_ref   = theta_samp;
phi_ref     = phi_samp;
% - ACQUISITION SETTINGS
% slit_ref = '1'; dt_ref = '0.1'; Eb_ref = '-86';   % XPS, Au(4d) reference
slit_ref = '20'; dt_ref = '30'; Eb_ref = '-2';      % ARPES, Fermi-Edge reference for metallic samples
pol_ref 	= pol_samp;
mode_ref	= mode_samp;
Ep_ref      = Ep_samp;

%% (D) Initialising the JobFile with the first entries
% -- Regular Energy Jobs
% ii=1;
% Position(ii).name='VB'; Position(ii).mode=mode_samp; Position(ii).Ep=Ep_samp;
% Position(ii).pol=pol_samp; Position(ii).slit=slitVB_samp; Position(ii).dT=dtVB_samp;
% Position(ii).Eb='-2.00';
% Position(ii).tilt=tilt_samp; Position(ii).theta=theta_samp; Position(ii).phi=phi_samp;
% Position(ii).x=x_samp; Position(ii).y=y_samp; Position(ii).z=z_samp;

ii=1;
Position(ii).name='In4d'; Position(ii).mode=mode_samp; Position(ii).Ep=Ep_samp;
Position(ii).pol=pol_samp; Position(ii).slit=slitXPS_samp; Position(ii).dT=dtXPS_samp;
Position(ii).Eb='-25:0.02:-10';
Position(ii).tilt=tilt_samp; Position(ii).theta=theta_samp; Position(ii).phi=phi_samp;
Position(ii).x=x_samp; Position(ii).y=y_samp; Position(ii).z=z_samp;

ii=ii+1;
Position(ii).name='As3d'; Position(ii).mode=mode_samp; Position(ii).Ep=Ep_samp;
Position(ii).pol=pol_samp; Position(ii).slit=slitXPS_samp; Position(ii).dT=dtXPS_samp;
Position(ii).Eb='-51:0.02:-37';
Position(ii).tilt=tilt_samp; Position(ii).theta=theta_samp; Position(ii).phi=phi_samp;
Position(ii).x=x_samp; Position(ii).y=y_samp; Position(ii).z=z_samp;

ii=ii+1;
Position(ii).name='Al2p'; Position(ii).mode=mode_samp; Position(ii).Ep=Ep_samp;
Position(ii).pol=pol_samp; Position(ii).slit=slitXPS_samp; Position(ii).dT=dtXPS_samp;
Position(ii).Eb='-82:0.02:-71';
Position(ii).tilt=tilt_samp; Position(ii).theta=theta_samp; Position(ii).phi=phi_samp;
Position(ii).x=x_samp; Position(ii).y=y_samp; Position(ii).z=z_samp;

ii=ii+1;
Position(ii).name='C1s'; Position(ii).mode=mode_samp; Position(ii).Ep=Ep_samp;
Position(ii).pol=pol_samp; Position(ii).slit=slitXPS_samp; Position(ii).dT=dtXPS_samp;
Position(ii).Eb='-290:0.02:-280';
Position(ii).tilt=tilt_samp; Position(ii).theta=theta_samp; Position(ii).phi=phi_samp;
Position(ii).x=x_samp; Position(ii).y=y_samp; Position(ii).z=z_samp;

% -- Special Energy 1 (reserved for deep core-levels - In3d, O1s)
kk=1;
sPosition1(kk).name='In3d';sPosition1(kk).mode=mode_samp;sPosition1(kk).Ep=Ep_samp;
sPosition1(kk).pol=pol_samp;sPosition1(kk).slit=slitXPS_samp;sPosition1(kk).dT=dtXPS_samp;
sPosition1(kk).Eb='-449:0.02:-440';
sPosition1(kk).tilt=tilt_samp;sPosition1(kk).theta=theta_samp;sPosition1(kk).phi=phi_samp;
sPosition1(kk).x=x_samp;sPosition1(kk).y=y_samp;sPosition1(kk).z=z_samp;

kk=kk+1;
sPosition1(kk).name='O1s';sPosition1(kk).mode=mode_samp;sPosition1(kk).Ep=Ep_samp;
sPosition1(kk).pol=pol_samp;sPosition1(kk).slit=slitXPS_samp;sPosition1(kk).dT=dtXPS_samp;
sPosition1(kk).Eb='-538:0.02:-529';
sPosition1(kk).tilt=tilt_samp;sPosition1(kk).theta=theta_samp;sPosition1(kk).phi=phi_samp;
sPosition1(kk).x=x_samp;sPosition1(kk).y=y_samp;sPosition1(kk).z=z_samp;

% -- Special Energy 2 (reserved for Survey scan)
hh=1;
sPosition2(hh).name='Survey';sPosition2(hh).mode=mode_samp;sPosition2(hh).Ep=Ep_samp;
sPosition2(hh).pol=pol_samp;sPosition2(hh).slit=slitXPS_samp;sPosition2(hh).dT='0.1';
sPosition2(hh).Eb='-555:0.05:5';
sPosition2(hh).tilt=tilt_samp;sPosition2(hh).theta=theta_samp;sPosition2(hh).phi=phi_samp;
sPosition2(hh).x=x_samp;sPosition2(hh).y=y_samp;sPosition2(hh).z=z_samp;

%% (E) Generate the JobFile Table
fid     = fopen([filename '.txt'],'wt');
head    = '/*File|Comment|hv|Pol|Slit|Mode|Epass|Eb/Ek|Range|dt(s)|Theta|Tilt|Azimuth|X|Y|Z*/'; fprintf(fid,'%s\n',head);
index   = 1;
% - Filing through each one of the photon energies
for hv = Energies 
	% -- Filing through Position
    for kk=1:length(Position)
        if use_ref == 1
            line=strcat(num2str(index, '%02.f'),'_',Position(kk).name,'_hv=',num2str(hv),'_EF||',num2str(hv),'|',pol_ref,'|',slit_ref,'|',mode_ref,'|',Ep_ref,'|','Eb','|',Eb_ref,'|',dt_ref,'|1|',theta_ref,'|',tilt_ref,'|',phi_ref,'|',x_ref,'|',y_ref,'|',z_ref,'|');
            fprintf(fid,'%s\n',line);
            index=index+1;
        end
        line=strcat(num2str(index, '%02.f'),'_',Position(kk).name,'_hv=',num2str(hv),'||',sprintf("ones(%i)*",Repeats)+num2str(hv),'|',Position(kk).pol,'|',Position(kk).slit,'|',Position(kk).mode,'|',Position(kk).Ep,'|','Eb','|',Position(kk).Eb,'|',Position(kk).dT,'|1|',Position(kk).theta,'|',Position(kk).tilt,'|',Position(kk).phi,'|',Position(kk).x,'|',Position(kk).y,'|',Position(kk).z,'|');
        fprintf(fid,'%s\n',line);
        index=index+1;
    end
    % -- Filing through each Special Energy 1
    if ismember(hv,sEnergy1)
        for kk=1:length(sPosition1)    
            if use_ref == 1
                line=strcat(num2str(index, '%02.f'),'_',sPosition1(kk).name,'_hv=',num2str(hv),'_EF||',num2str(hv),'|',pol_ref,'|',slit_ref,'|',mode_ref,'|',Ep_ref,'|','Eb','|',Eb_ref,'|',dt_ref,'|1|',theta_ref,'|',tilt_ref,'|',phi_ref,'|',x_ref,'|',y_ref,'|',z_ref,'|');
                fprintf(fid,'%s\n',line);
                index=index+1;
            end
            line=strcat(num2str(index, '%02.f'),'_',sPosition1(kk).name,'_hv=',num2str(hv),'||',sprintf("ones(%i)*",sRepeats1)+num2str(hv),'|',sPosition1(kk).pol,'|',sPosition1(kk).slit,'|',sPosition1(kk).mode,'|',sPosition1(kk).Ep,'|','Eb','|',sPosition1(kk).Eb,'|',sPosition1(kk).dT,'|1|',sPosition1(kk).theta,'|',sPosition1(kk).tilt,'|',sPosition1(kk).phi,'|',sPosition1(kk).x,'|',sPosition1(kk).y,'|',sPosition1(kk).z,'|');
            fprintf(fid,'%s\n',line);
            index=index+1;
        end
    end
   % -- Filing through each Special Energy 2
    if ismember(hv,sEnergy2)
        for kk=1:length(sPosition2)    
            if use_ref == 1
                line=strcat(num2str(index, '%02.f'),'_',sPosition2(kk).name,'_hv=',num2str(hv),'_EF||',num2str(hv),'|',pol_ref,'|',slit_ref,'|',mode_ref,'|',Ep_ref,'|','Eb','|',Eb_ref,'|',dt_ref,'|1|',theta_ref,'|',tilt_ref,'|',phi_ref,'|',x_ref,'|',y_ref,'|',z_ref,'|');
                fprintf(fid,'%s\n',line);
                index=index+1;
            end
            line=strcat(num2str(index, '%02.f'),'_',sPosition2(kk).name,'_hv=',num2str(hv),'||',sprintf("ones(%i)*",sRepeats2)+num2str(hv),'|',sPosition2(kk).pol,'|',sPosition2(kk).slit,'|',sPosition2(kk).mode,'|',sPosition2(kk).Ep,'|','Eb','|',sPosition2(kk).Eb,'|',sPosition2(kk).dT,'|1|',sPosition2(kk).theta,'|',sPosition2(kk).tilt,'|',sPosition2(kk).phi,'|',sPosition2(kk).x,'|',sPosition2(kk).y,'|',sPosition2(kk).z,'|');
            fprintf(fid,'%s\n',line);
            index=index+1;
        end
    end
end

fclose(fid);
end