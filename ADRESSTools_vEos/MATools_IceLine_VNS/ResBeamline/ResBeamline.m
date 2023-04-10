function varargout = ResBeamline(varargin)
%RESBEAMLINE M-file for ResBeamline.fig
%      RESBEAMLINE, by itself, creates a new RESBEAMLINE or raises the existing
%      singleton*.
%      H = RESBEAMLINE returns the handle to a new RESBEAMLINE or the handle to
%      the existing singleton*.
%      RESBEAMLINE('Property','Value',...) creates a new RESBEAMLINE using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to ResBeamline_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%      RESBEAMLINE('CALLBACK') and RESBEAMLINE('CALLBACK',hObject,...) call the
%      local function named CALLBACK in RESBEAMLINE.M with the given input
%      arguments.
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ResBeamline
% Last Modified by GUIDE v2.5 14-Feb-2023 18:06:32
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ResBeamline_OpeningFcn, ...
                   'gui_OutputFcn',  @ResBeamline_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
   gui_State.gui_Callback = str2func(varargin{1});
end
if nargout 
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --- Executes just before ResBeamline is made visible.
function ResBeamline_OpeningFcn(hObject,~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% Choose default command line output for ResBeamline
handles.output = hObject;
% beamline parameters
handles.sigmaSizeV_ID=0.014; % rms vertical source size
handles.sigmaDivV_ID=0.014E-3; % rms vertical divergence of the coherent core
handles.fCM=17750;          % distance between the source and CM
handles.graLength=90;       % grating optical surface length
handles.fFM=9510;           % dispersion length
handles.Gratings=get(handles.popupmenuGrating,'String');    % available gratings
handles.Cff=[2.15 3 2.55];    % default Cff values
handles.Orders=get(handles.popupmenuOrder,'String');        % available orders
handles.sigmaSE_PM=1e-6*0.252; % rms slope error PM = 0.052 arcsec
handles.sigmaSE_Gr=1e-6*[0.07*4.848137 0.35 0.079]; % rms slope error of the grating
% handles.Cff=Res(:,1); handles.SigmaV=Res(:,2:end)/4;
% Update handles structure
guidata(hObject, handles);
% UIWAIT makes ResBeamline wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = ResBeamline_OutputFcn(~,~, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% Get default command line output from handles structure
varargout{1} = handles.output;

% % --- Executes on button press in pushbuttonUpdate.
% function pushbuttonUpdate_Callback(hObject, eventdata, handles)
% % hObject    handle to pushbuttonUpdate (see GCBO)
% set([handles.editSlitWidth handles.editResolution],'BackgroundColor','w')
% % retrieve the actual parameters
% [hEnergy,hGratingNo,hOrd,hCff,hSlitWidth]=mcaopen('X03MA-PGM:energy','X03MA-PGM-GRCH:GRATING',...
%                                              'X03MA-PGM:difforder0','X03MA-PGM:cff','X03MA-OP-SL:TRY_AP');
% [energy,gratingNo,ord,cff,slitWidth]=mcaget(hEnergy,hGratingNo,hOrd,hCff,hSlitWidth);
% mcaexit()
% % display the actual parameters
% set(handles.editEnergy,'String',num2str(energy));
% % gratingNo=2-gratingNo; % reverse grating order from EPICS
% set(handles.popupmenuGrating,'Value',gratingNo+1); 
% set(handles.popupmenuOrder,'Value',ord+1);
% set(handles.editCff,'String',num2str(cff));
% set(handles.editSlitWidth,'String',num2str(slitWidth)); slitWidth=slitWidth/1000;
% % try to calculate the resolution parameters
% [resLimit,resCoeff]=resParams(handles); if isempty(resLimit*resCoeff) return; end
% % slit limited resolution
% resSlit=slitWidth*energy^2/(resCoeff*1240*1e-6);
% % total resolution
% resTotal=sqrt(resSlit^2+resLimit^2);
% % display resolution
% set(handles.editResolution,'String',num2str(1000*resTotal))

% Cff solver
function alphabeta=cffSolver(cff,N,k,lambda)
% alphabeta=cffSolver(cff,N,k,lambda) finds the incident and exit angles 
% [alpha beta](rad) from the cff value, groove density N(l/mm), 
% diffraction order k and wavelength lambda(mm)
a=1-(1/cff)^2; b=-2*N*k*lambda; c=(1/cff)^2+(N*k*lambda)^2-1;
if cff==1
   disp('cffSolver => error: invalid cff=1'); beta=NaN;
else
   X=[-1*b-sqrt(b^2-4*a*c) -1*b+sqrt(b^2-4*a*c)]/(2*a);
   Beta=asin(X);
   [~,ind]=min( abs( imag(Beta) ) ); beta=Beta(ind);      % selecting the real value
end
% alpha
alpha=asin(N*k*lambda-sin(beta));
% output
alphabeta=[alpha beta];

function editEnergy_Callback(hObject,~, handles)
% hObject    handle to editEnergy (see GCBO)
% Hints: get(hObject,'String') returns contents of editEnergy as text
%        str2double(get(hObject,'String')) returns contents of editEnergy as a double
set([handles.editSlitWidth handles.editResolution],'String','')
% get energy value 
energy=str2num(get(hObject,'String'));
if ~isempty(energy) && (energy<160 || energy>1800); set(hObject,'String',''); return; end
% try to calculate the resolution parameters
resParams(handles);

% --- Executes on selection change in popupmenuGrating.
function popupmenuGrating_Callback(hObject,~, handles)
% hObject    handle to popupmenuGrating (see GCBO)
% Hints: contents = get(hObject,'String') returns popupmenuGrating contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenuGrating
set([handles.editSlitWidth handles.editResolution],'String','')
% default Cff
gratingNo=get(hObject,'Value'); cff=handles.Cff(gratingNo); set(handles.editCff,'String',num2str(cff))
% try to calculate the resolution parameters
resParams(handles);

% --- Executes on selection change in popupmenuOrder.
function popupmenuOrder_Callback(~,~, handles)
% hObject    handle to popupmenuOrder (see GCBO)
% Hints: contents = get(hObject,'String') returns popupmenuOrder contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenuOrder
set([handles.editSlitWidth handles.editResolution],'String','')
% try to calculate the resolution parameters
resParams(handles);

function editCff_Callback(hObject,~, handles)
% hObject    handle to editCff (see GCBO)
% Hints: get(hObject,'String') returns contents of editCff as text
%        str2double(get(hObject,'String')) returns contents of editCff as a double
set([handles.editSlitWidth handles.editResolution],'String','')
% get Cff value
cff=str2num(get(hObject,'String'));
if isempty(cff) || cff<1 || cff>10; set(hObject,'String',''); return; end
% try to calculate the resolution parameters
resParams(handles);

function editSlitWidth_Callback(hObject,~,handles)
% hObject    handle to editSlitWidth (see GCBO)
% Hints: get(hObject,'String') returns contents of editSlitWidth as text
%        str2double(get(hObject,'String')) returns contents of editSlitWidth as a double
set(handles.editResolution,'String','')
slitWidth=str2num(get(hObject,'String'))/1000; 
if isempty(slitWidth); set([hObject handles.editResolution],'String',''); return; end
energy=str2num(get(handles.editEnergy,'String'));
% try to calculate the resolution parameters
[resLimit,resCoeff]=resParams(handles); if isempty(resLimit*resCoeff) return; end
% slit limited resolution
resSlit=slitWidth*energy^2/(resCoeff*1240*1e-6);
% total resolution
resTotal=sqrt(resSlit^2+resLimit^2);
% display resolution
set(handles.editResolution,'String',num2str(round(1e6*resTotal)/1000))

function editResolution_Callback(hObject,~,handles)
% hObject    handle to editResolution (see GCBO)
% Hints: get(hObject,'String') returns contents of editResolution as text
%        str2double(get(hObject,'String')) returns contents of editResolution as a double
set(handles.editSlitWidth,'String','')
res=str2num(get(hObject,'String'))/1000; 
if isempty(res); set([hObject handles.editSlitWidth],'String',''); return; end
energy=str2num(get(handles.editEnergy,'String'));
% try to calculate the resolution parameters
[resLimit,resCoeff]=resParams(handles); if isempty(resLimit*resCoeff) return; end
% check the resolution limit
if res<resLimit res=resLimit; set(hObject,'String',num2str(res*1000)); end
% slit limited resolution
resSlit=sqrt(res^2-resLimit^2);
% slit width
slitWidth=resSlit*resCoeff*1240*1e-6/energy^2;
% display slit width
set(handles.editSlitWidth,'String',num2str(1000*slitWidth))

% --- Evaluates the resolution parameters from energy, groove density and cff
function [resLimit,resCoeff]=resParams(handles)
% input
% rms to fwhm
rms2fwhm=2*sqrt(2*log(2));
% - energy
energy=str2num(get(handles.editEnergy,'String'));
% - n
gratingNo=get(handles.popupmenuGrating,'Value');
n=str2num(handles.Gratings{gratingNo});
% - k
orderNo=get(handles.popupmenuOrder,'Value'); k=str2num(handles.Orders{orderNo});
% - Cff
cff=str2num(get(handles.editCff,'String')); if isempty(cff); resLimit=[]; return; end
% % sigmaV
% sigmaV=interp1(handles.Cff,handles.SigmaV(:,gratingNo),cff,'pchip');
% % display spot size
% set(handles.textSpotSize,'String',['Spot Size FWHM (um) = ' num2str(1000*2*sqrt(2*log(2))*sigmaV)])
% set(handles.textSpotSize,'String',['4-sigma Spot Size (um) = ' num2str(1000*4*sigmaV)])
% % if not zero order
% set([handles.editResolution handles.editSlitWidth],'Enable','On');
% wavelength
lambda=1e-6*1240/energy;
% alpha and beta
alphabeta=cffSolver(cff,n,k,lambda); alpha=alphabeta(1); beta=alphabeta(2);
% T-parameter
T=sqrt(lambda*k*n*(cff^2-1)/2);
% 4-sigma grating illumination
gratingIllumination=(4*handles.sigmaDivV_ID*handles.fCM/cos(alpha));
set(handles.textGratingIllumination,'String',['Grating Illumination (mm) = ' num2str(gratingIllumination)])
% resolving power
% - source contribution
handles.sigmaSizeV_ID*rms2fwhm;
R_source=handles.fCM*T/(handles.sigmaSizeV_ID*rms2fwhm);
% - PM-SE contribution
R_PM=T/(2*handles.sigmaSE_PM*rms2fwhm);
% - grating-SE contribution
R_Gr=T/(handles.sigmaSE_Gr(gratingNo)*rms2fwhm*(1+cff));
% - illuminated grooves
incAngleGra=pi/2-alpha; footGraV=handles.sigmaDivV_ID*rms2fwhm*handles.fCM/sin(incAngleGra);
R_N=min(handles.graLength,footGraV)*n;
% - total
R_tot=1/sqrt(R_source^-2+R_PM^-2+R_Gr^-2+R_N^-2);
resLimit=energy/R_tot;
% display resolution limit
set(handles.textResolutionLimit,'String',['Resolution Limit (meV) = ' num2str(1000*resLimit)])
% spot size at the slit
slitV=sqrt( (handles.sigmaSizeV_ID*rms2fwhm*handles.fFM/(cff*handles.fCM))^2 + ...
            (2*handles.sigmaSE_PM*rms2fwhm*handles.fFM/cff)^2 + ...
            (handles.sigmaSE_Gr(gratingNo)*rms2fwhm*(1+1/cff)*handles.fFM)^2 );
% display spot size at the slit
set(handles.textSpotSize,'String',['Focused spot size (um) = ' num2str(1000*slitV)])
% resolution coefficient
resCoeff=handles.fFM*n*k/cos(beta);
