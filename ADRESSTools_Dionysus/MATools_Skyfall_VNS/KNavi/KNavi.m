function varargout = KNavi(varargin)
% KNAVI MATLAB code for KNavi.fig
%      KNAVI, by itself, creates a new KNAVI or raises the existing
%      singleton*.
%      H = KNAVI returns the handle to a new KNAVI or the handle to
%      the existing singleton*.
%      KNAVI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in KNAVI.M with the given input arguments.
%      KNAVI('Property','Value',...) creates a new KNAVI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before KNavi_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to KNavi_OpeningFcn via varargin.
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
% See also: GUIDE, GUIDATA, GUIHANDLES
% Edit the above text to modify the response to help KNavi
% Last Modified by GUIDE v2.5 02-May-2020 12:45:34
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @KNavi_OpeningFcn, ...
                   'gui_OutputFcn',  @KNavi_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
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

% --- Executes just before KNavi is made visible.
function KNavi_OpeningFcn(hObject,~,handles,varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to KNavi (see VARARGIN)
% Choose default command line output for KNavi
handles.output = hObject;
% Update handles structure
guidata(hObject, handles);
% parameters
clear global
global ePhi alpha eB
ePhi=4.5; alpha=20; eB=0; % surfNx calculated for arbitrary eB, and all other vcalculations assume eB=0
% image the scheme
% - read in image
IS=imread('Scheme.png');
% - match the axis size to image
Origin=get(handles.axes1,'Position'); Origin=Origin(1:2);
WH=size(IS); WH=WH(1:2);
set(handles.axes1,'Position',[Origin WH]);
% - plot image
image(IS,'Parent',handles.axes1); axis image off

% UIWAIT makes KNavi wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = KNavi_OutputFcn(~,~,handles)
% Get default command line output from handles structure
varargout{1} = handles.output;

% ---------------- Origin section ----------------

% --- Executes on button press in pushbuttonOrgCalc.
function pushbuttonOrgCalc_Callback(~,~,handles)
global ePhi
% initial parameters
set(handles.textOrgSurfNormX,'String','')
% read in the parameters
hv=str2num(get(handles.editOrgHv,'String')); if hv<ePhi; hv=[]; end
eB=str2num(get(handles.editOrgEB,'String')); 
kx=str2num(get(handles.editOrgKx,'String')); 
thtM=str2num(get(handles.editOrgThtM,'String')); 
thtA=str2num(get(handles.editOrgThtA,'String'));
% surface normal calculations
surfNx = SurfNormX(hv,eB,kx,thtM,thtA);
% display results
set(handles.textOrgSurfNormX,'String',['surfNx = ' num2str(surfNx)])
set([handles.editKxSurfNx handles.editKySurfNx handles.editKzSurfNx], 'String',num2str(surfNx))

% ---------------- Kx section ----------------

% --- Executes on button press in pushbuttonKxCalc
function pushbuttonKxCalc_Callback(~,~,handles)
global ePhi alpha eB
% read in the parameters
surfNx=str2num(get(handles.editKxSurfNx,'String'));
hv=str2num(get(handles.editKxHv,'String')); if hv<ePhi; hv=[]; end
kx=str2num(get(handles.editKxKx,'String')); 
thtM=str2num(get(handles.editKxThtM,'String')); thtMPhys=thtM+surfNx;
thtA=str2num(get(handles.editKxThtA,'String'));
% calculation of kx/thtM/thtA
switch get(handles.popupmenuKx,'Value')
   case 1
      if isempty(surfNx*hv*thtM*thtA); msgbox('Some Inputs or Parameters Empty','modal'); return; end
      kx=0.5124*sqrt(hv-ePhi+eB)*sind(thtMPhys+thtA)-2*pi*hv*cosd(alpha+thtMPhys)/12400;
      set(handles.editKxKx,'String',num2str(kx))
   case 2
      if isempty(surfNx*hv*kx*thtA); msgbox('Some Inputs or Parameters Empty','modal'); return; end 
      [thtM,~,exitFlag]=fsolve(@(x) 0.5124*sqrt(hv-ePhi+eB)*sind(x+surfNx+thtA)-...
         2*pi*hv*cosd(alpha+x+surfNx)/12400-kx,0,optimset('Display','off','TolFun',1e-12));
      if exitFlag~=1; thtM=[]; msgbox('Error: Check Input Parameters','modal'); end
      set(handles.editKxThtM,'String',num2str(thtM))
   case 3
      if isempty(surfNx*hv*kx*thtM); msgbox('Some Inputs or Parameters Empty','modal'); return; end 
      [thtA,~,exitFlag]=fsolve(@(x) 0.5124*sqrt(hv-ePhi+eB)*sind(thtMPhys+x)-...
         2*pi*hv*cosd(alpha+thtMPhys)/12400-kx,0,optimset('Display','off','TolFun',1e-12));
      if exitFlag~=1; thtA=[]; msgbox('Error: Check Input Parameters','modal'); end
      set(handles.editKxThtA,'String',num2str(thtA))
end
% transfer hv value to the Ky and Kz sections
set([handles.editKyHv handles.editKzHv],'String',num2str(hv))
% transfer Kx to the Kz section
set(handles.editKzKx,'String',num2str(kx))
% transfer thtM value to the Ky and Kz sections
set([handles.editKyThtM handles.editKzThtM],'String',num2str(thtM))

% ---------------- Ky section ----------------

% --- Executes on button press in pushbuttonKyCalc
function pushbuttonKyCalc_Callback(~,~,handles)
global ePhi alpha eB
% read in the parameters
thtA=str2num(get(handles.editKxThtA,'String'));
surfNx=str2num(get(handles.editKySurfNx,'String'));
thtM=str2num(get(handles.editKyThtM,'String')); thtMPhys=thtM+surfNx;
surfNy=str2num(get(handles.editKySurfNy,'String'));
hv=str2num(get(handles.editKyHv,'String')); if hv<ePhi; hv=[]; end
ky=str2num(get(handles.editKyKy,'String'));
tlt=str2num(get(handles.editKyTlt,'String')); tlt=tlt-surfNy;
% calculation of ky/TLT
switch get(handles.popupmenuKy,'Value')
   case 1
      if isempty(thtM*surfNx*hv*tlt); msgbox('Some Inputs or Parameters Empty','modal'); return; end 
      ky=0.5124*sqrt(hv-ePhi+eB)*sind(tlt)*cosd(thtA+thtMPhys)+2*pi*hv*sind(alpha+thtMPhys)*sind(tlt)/12400;
      set(handles.editKyKy,'String',num2str(ky))
   case 2
      if isempty(thtM*surfNx*hv*ky); msgbox('Some Inputs or Parameters Empty','modal'); return; end  
      [tlt,~,exitFlag]=fsolve(@(x) 0.5124*sqrt(hv-ePhi+eB)*sind(x)*cosd(thtA+thtMPhys)+...
         2*pi*hv*sind(alpha+thtMPhys)*sind(x)/12400-ky,0,optimset('Display','off','TolFun',1e-12));
      if exitFlag~=1; tlt=[]; msgbox('Error: Check Input Parameters','modal'); end      
      set(handles.editKyTlt,'String',num2str(tlt+surfNy))
end
% transfer Ky to the Kz section
set(handles.editKzKy,'String',num2str(ky))

% ---------------- Kz section ----------------

% --- Executes on button press in pushbuttonKzCalc.
function pushbuttonKzCalc_Callback(~,~,handles)
global ePhi alpha eB
% inputs
surfNx=str2num(get(handles.editKzSurfNx,'String'));
thtM=str2num(get(handles.editKzThtM,'String')); thtMPhys=thtM+surfNx;
kx=str2num(get(handles.editKzKx,'String')); 
ky=str2num(get(handles.editKzKy,'String')); 
v000=str2num(get(handles.editKzV000,'String'));
hv=str2num(get(handles.editKzHv,'String')); if hv<ePhi; hv=[]; end
kz=str2num(get(handles.editKzKz,'String')); 
c=str2num(get(handles.editKzC,'String')); if ~isempty(c); tpc=2*pi/c; else tpc=[]; end
% calculation of kz/hv/V000
switch get(handles.popupmenuKz,'Value')
   case 1
      if isempty(thtM*surfNx*kx*ky*v000*hv); msgbox('Some Inputs or Parameters Empty','modal'); return; end  
      kz=0.5124*sqrt(hv+eB+v000-3.81*(kx^2+ky^2))+2*pi*hv*sind(alpha+thtMPhys)/12400;
      if get(handles.popupmenuKzUnits,'Value')==2; kz=kz/tpc; end   % kz units
      set(handles.editKzKz,'String',num2str(kz))
   case 2
      if isempty(thtM*surfNx*kx*ky*v000*kz); msgbox('Some Inputs or Parameters Empty','modal'); return; end 
      if get(handles.popupmenuKzUnits,'Value')==2; kz=kz*tpc; end   % kz units        
      [hv,~,exitFlag]=fsolve(@(x) 0.5124*sqrt(x+eB+v000-3.81*(kx^2+ky^2))+...
         2*pi*x*sind(alpha+thtMPhys)/12400-kz,0,optimset('Display','off','TolFun',1e-12));
      if exitFlag~=1; hv=[]; msgbox('Error: Check Input Parameters','modal'); end
      hv=real(hv);    % trim the imaginary residual error
      set(handles.editKzHv,'String',num2str(hv))
    case 3
      if isempty(thtM*surfNx*kx*ky*hv*kz); msgbox('Some Inputs or Parameters Empty','modal'); return; end 
      if get(handles.popupmenuKzUnits,'Value')==2; kz=kz*tpc; end   % kz units 
      v000=((kz-2*pi*hv*sind(alpha+thtMPhys)/12400)/0.5124)^2-hv-eB+3.81*(kx^2+ky^2);
      set(handles.editKzV000,'String',num2str(v000))   
end
% --- Executes on selection change in popupmenuKzUnits.
function popupmenuKzUnits_Callback(hObject,~,handles)
if get(hObject,'Value')==1
   set([handles.editKzC handles.textKzC],'Enable','off'); 
else
   set([handles.editKzC handles.textKzC],'Enable','on');
end
