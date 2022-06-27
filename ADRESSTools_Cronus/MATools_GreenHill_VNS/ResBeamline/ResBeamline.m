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
% Last Modified by GUIDE v2.5 29-Feb-2016 11:35:56
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
function ResBeamline_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% Choose default command line output for ResBeamline
handles.output = hObject;
% beamline parameters
handles.sigmaV_ID=0.014E-3; % vertical divergence of the coherent core
handles.fCM=17750;          % distance between the source and CM
handles.graLength=90;       % grating optical surface length
handles.fFM=9510;           % dispersion length
handles.Gratings=get(handles.popupmenuGrating,'String');    % available gratings
handles.Orders=get(handles.popupmenuOrder,'String');        % available orders
% results of ray tracing
Res=[1.00   0.050944241   0.050147058   0.048947233 ;...
     1.25   0.043866336   0.043192697   0.041464040 ;...
     1.50   0.039387180   0.038665573   0.036774471 ;...
     1.75   0.036796608   0.035641532   0.033910805 ;...
     2.00   0.034732560   0.033670998   0.031710985 ;...
     2.25   0.033312806   0.032375647   0.030001851 ;...
     2.50   0.032252078   0.031215571   0.028830829 ;...
     2.75   0.031627573   0.030258899   0.027933936 ;...
     3.00   0.030953200   0.029586359   0.027405348 ;...
     3.50   0.029808244   0.028834696   0.026129788 ;...
     4.00   0.029470368   0.028264421   0.025539756 ;...
     4.50   0.028950357   0.027646038   0.025166420 ;...
     5.00   0.028717878   0.027298163   0.024918597 ;...
     5.50   0.028479262   0.027152126   0.024656275 ;...
     6.00   0.028330398   0.027032641   0.024477859 ;...
     7.00   0.028072773   0.026834474   0.024244989 ;...
     8.00   0.028017667   0.026788041   0.024135201 ;...
     9.00   0.027826015   0.026671926   0.024099521 ;...
    10.00   0.027652850   0.026557109   0.024013327 ]; 
handles.Cff=Res(:,1); handles.SigmaV=Res(:,2:end)/4;
% Update handles structure
guidata(hObject, handles);
% UIWAIT makes ResBeamline wait for user response (see UIRESUME)
% uiwait(handles.figure1);

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
   [dum,ind]=min( abs( imag(Beta) ) ); beta=Beta(ind);      % selecting the real value
end
% alpha
alpha=asin(N*k*lambda-sin(beta));
% output
alphabeta=[alpha beta];

% --- Outputs from this function are returned to the command line.
function varargout = ResBeamline_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% Get default command line output from handles structure
varargout{1} = handles.output;

function editEnergy_Callback(hObject, eventdata, handles)
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
function popupmenuGrating_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuGrating (see GCBO)
% Hints: contents = get(hObject,'String') returns popupmenuGrating contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenuGrating
set([handles.editSlitWidth handles.editResolution],'String','')
% try to calculate the resolution parameters
resParams(handles);

% --- Executes on selection change in popupmenuOrder.
function popupmenuOrder_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuOrder (see GCBO)
% Hints: contents = get(hObject,'String') returns popupmenuOrder contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenuOrder
set([handles.editSlitWidth handles.editResolution],'String','')
% try to calculate the resolution parameters
resParams(handles);

function editCff_Callback(hObject, eventdata, handles)
% hObject    handle to editCff (see GCBO)
% Hints: get(hObject,'String') returns contents of editCff as text
%        str2double(get(hObject,'String')) returns contents of editCff as a double
set([handles.editSlitWidth handles.editResolution],'String','')
% get Cff value
cff=str2num(get(hObject,'String'));
if isempty(cff) || cff<1 || cff>10; set(hObject,'String',''); return; end
% try to calculate the resolution parameters
resParams(handles);

function editSlitWidth_Callback(hObject, eventdata, handles)
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
set(handles.editResolution,'String',num2str(1000*resTotal))

function editResolution_Callback(hObject, eventdata, handles)
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

% --- Executes on button press in pushbuttonUpdate.
function pushbuttonUpdate_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonUpdate (see GCBO)
set([handles.editSlitWidth handles.editResolution],'BackgroundColor','w')
% retrieve the actual parameters
[hEnergy,hGratingNo,hOrd,hCff,hSlitWidth]=mcaopen('X03MA-PGM:energy','X03MA-PGM-GRCH:GRATING',...
                                             'X03MA-PGM:difforder0','X03MA-PGM:cff','X03MA-OP-SL:TRY_AP');
[energy,gratingNo,ord,cff,slitWidth]=mcaget(hEnergy,hGratingNo,hOrd,hCff,hSlitWidth);
mcaexit()
% display the actual parameters
set(handles.editEnergy,'String',num2str(energy));
% gratingNo=2-gratingNo; % reverse grating order from EPICS
set(handles.popupmenuGrating,'Value',gratingNo+1); 
set(handles.popupmenuOrder,'Value',ord+1);
set(handles.editCff,'String',num2str(cff));
set(handles.editSlitWidth,'String',num2str(slitWidth)); slitWidth=slitWidth/1000;
% try to calculate the resolution parameters
[resLimit,resCoeff]=resParams(handles); if isempty(resLimit*resCoeff) return; end
% slit limited resolution
resSlit=slitWidth*energy^2/(resCoeff*1240*1e-6);
% total resolution
resTotal=sqrt(resSlit^2+resLimit^2);
% display resolution
set(handles.editResolution,'String',num2str(1000*resTotal))

% --- Evaluates the resolution parameters from energy, groove density and cff
function [resLimit,resCoeff]=resParams(handles)
% input
% - energy
energy=str2num(get(handles.editEnergy,'String'));
% - n
gratingNo=get(handles.popupmenuGrating,'Value');
n=str2num(handles.Gratings{gratingNo});
% - k
orderNo=get(handles.popupmenuOrder,'Value'); k=str2num(handles.Orders{orderNo});
% - Cff
cff=str2num(get(handles.editCff,'String')); if isempty(cff); resLimit=[]; return; end
% sigmaV
sigmaV=interp1(handles.Cff,handles.SigmaV(:,gratingNo),cff,'pchip');
% display spot size
set(handles.textSpotSize,'String',['Spot Size FWHM (um) = ' num2str(1000*2*sqrt(2*log(2))*sigmaV)])
% set(handles.textSpotSize,'String',['4-sigma Spot Size (um) = ' num2str(1000*4*sigmaV)])
% if zero order or no energy
if cff<1.01 || isempty(energy)
   set([handles.editResolution handles.editSlitWidth],'Enable','Off');
   set([handles.textGratingIllumination handles.textResolutionLimit],'String',''); 
   resLimit=[]; resCoeff=[]; return;
end
% if not zero order
set([handles.editResolution handles.editSlitWidth],'Enable','On');
% wavelength
lambda=1e-6*1240/energy;
% alpha and beta
alphabeta=cffSolver(cff,n,k,lambda); alpha=alphabeta(1); beta=alphabeta(2);
% 4-sigma grating illumination
gratingIllumination=(4*handles.sigmaV_ID*handles.fCM/cos(alpha))/handles.graLength;
set(handles.textGratingIllumination,'String',['Grating Illumination (%) = ' num2str(100*gratingIllumination)])
% diffraction limited resolution
illumLength=min(handles.graLength, 2*sqrt(2*log(2))*handles.sigmaV_ID*handles.fCM/cos(alpha));
resPower=illumLength*n*k; resDiffr=energy/resPower;
% spot size limited resolution
resCoeff=handles.fFM*n*k/cos(beta); resPower=resCoeff*lambda/(2*sqrt(2*log(2))*sigmaV);
resSpotSize=energy/resPower;
% total resolution limit
resLimit=sqrt(resDiffr^2+resSpotSize^2);
% display resolution limit
set(handles.textResolutionLimit,'String',['Resolution Limit (meV) = ' num2str(1000*resLimit)])
