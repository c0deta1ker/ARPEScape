function varargout = ViewSetup(varargin)
% VIEWSETUP MATLAB code for ViewSetup.fig
%      VIEWSETUP, by itself, creates a new VIEWSETUP or raises the existing
%      singleton*.
%      H = VIEWSETUP returns the handle to a new VIEWSETUP or the handle to
%      the existing singleton*.
%      VIEWSETUP('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VIEWSETUP.M with the given input arguments.
%      VIEWSETUP('Property','Value',...) creates a new VIEWSETUP or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ViewSetup_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ViewSetup_OpeningFcn via varargin.
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
% Last Modified by GUIDE v2.5 05-Nov-2020 12:36:57
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ViewSetup_OpeningFcn, ...
                   'gui_OutputFcn',  @ViewSetup_OutputFcn, ...
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

% --- Executes just before ViewSetup is made visible.
function ViewSetup_OpeningFcn(hObject,~,handles,varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% varargin   command line arguments to ViewSetup (see VARARGIN)
% Choose default command line output for ViewSetup
handles.output = hObject;
% Update handles structure
guidata(hObject, handles);
% figure properties
% set(gcf,'WindowStyle','modal')
% global parameters
% - initial values
global FracView FracViewIni XRange YRange viewMode gamma gammaIni render renderIni
if viewMode==2; set([handles.editFracView handles.editGamma handles.popupmenuRender],'Enable','Off');
else set([handles.editFracView handles.editGamma handles.popupmenuRender],'Enable','On'); end
if isempty(FracView); FracView=FracViewIni; end; set(handles.editFracView,'String',sprintf('%g:%g',FracView));
if isempty(XRange); set(handles.editXRange,'String','Full'); else; set(handles.editXRange,'String',sprintf('%g:%g',XRange)); end
if isempty(YRange); set(handles.editYRange,'String','Full'); else; set(handles.editYRange,'String',sprintf('%g:%g',YRange)); end
gammaIni=1; if isempty(gamma); gamma=gammaIni; end; set(handles.editGamma,'String',gamma);
renderIni='Flat'; if isempty(render); render=renderIni; end; if isequal(render,'Flat'); set(handles.popupmenuRender,'value',1); else set(handles.popupmenuRender,'value',2); end  
% % - figure position
% global xPos_View yPos_View
% Pos=get(gcf,'position'); offsetFigHeader=29;
% set(gcf,'position',[xPos_View yPos_View-Pos(4)-offsetFigHeader Pos(3) Pos(4)]);

% --- Outputs from this function are returned to the command line.
function varargout = ViewSetup_OutputFcn(~,~,handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in editFracView.
function editFracView_Callback(hObject,~,~)
global FracView FracViewIni
FracView=field2num(hObject); if length(FracView)~=2; FracView=FracViewIni; set(hObject,'String',sprintf('%g:%g',FracView)); end

% --- Executes on button press in editXRange.
function editXRange_Callback(hObject,~,~)
global XRange
XRange=field2num(hObject); if length(XRange)~=2; XRange=[]; set(hObject,'String','Full'); end

% --- Executes on button press in editYRange.
function editYRange_Callback(hObject,~,~)
global YRange
YRange=field2num(hObject); if length(YRange)~=2; YRange=[]; set(hObject,'String','Full'); end

function editGamma_Callback(hObject,~,~)
global gamma
gamma=field2num(hObject,1); set(hObject,'String',num2str(gamma))

% --- Executes on selection change in popupmenuRender.
function popupmenuRender_Callback(hObject,~,~)
global render
val=get(hObject,'Value'); if val==1; render='Flat'; else render='Interp'; end
