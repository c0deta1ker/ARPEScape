function varargout = LockSetup(varargin)
% LOCKSETUP MATLAB code for LockSetup.fig
%      LOCKSETUP, by itself, creates a new LOCKSETUP or raises the existing
%      singleton*.
%      H = LOCKSETUP returns the handle to a new LOCKSETUP or the handle to
%      the existing singleton*.
%      LOCKSETUP('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in LOCKSETUP.M with the given input arguments.
%      LOCKSETUP('Property','Value',...) creates a new LOCKSETUP or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before LockSetup_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to LockSetup_OpeningFcn via varargin.
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
% Last Modified by GUIDE v2.5 26-Apr-2020 14:27:52
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @LockSetup_OpeningFcn, ...
                   'gui_OutputFcn',  @LockSetup_OutputFcn, ...
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

% --- Executes just before LockSetup is made visible.
function LockSetup_OpeningFcn(hObject,~,handles,varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% varargin   command line arguments to LockSetup (see VARARGIN)
% Choose default command line output for LockSetup
handles.output = hObject;
% Update handles structure
guidata(hObject, handles);
% figure properties
% set(gcf,'WindowStyle','modal')
% global parameters
% - initial values
global feat_Lock setTo_Lock dESm_Lock dEWin_Lock 
if isequal(feat_Lock,'Peak'); set(handles.popupmenuFeat,'Value',2); else set(handles.popupmenuFeat,'Value',1); end 
set(handles.editSetTo,'String',num2str(setTo_Lock));
if isempty(dESm_Lock); set(handles.editSmooth,'String','Auto'); else; set(handles.editSmooth,'String',num2str(dESm_Lock)); end
if isempty(dEWin_Lock); set(handles.editSearch,'String','Auto'); else; set(handles.editSearch,'String',num2str(dEWin_Lock)); end
% % - figure position
% global xPos_Lock yPos_Lock
% Pos=get(gcf,'position'); offsetFigHeader=29;
% set(gcf,'position',[xPos_Lock yPos_Lock-Pos(4)-offsetFigHeader Pos(3) Pos(4)]);

% --- Outputs from this function are returned to the command line.
function varargout = LockSetup_OutputFcn(~,~,handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on selection change in popupmenuFeat.
function popupmenuFeat_Callback(hObject,~,~)
global feat_Lock
val=get(hObject,'Value'); if val==2; feat_Lock='Peak'; end

% --- Executes on button press in editSetTo.
function editSetTo_Callback(hObject,~,~)
global setTo_Lock
setTo_Lock=str2num(get(hObject,'String'));
if isempty(setTo_Lock); set(hObject,'String','0.0'); end

% --- Executes on button press in editSmooth.
function editSmooth_Callback(hObject,~,~)
global dESm_Lock
dESm_Lock=str2num(get(hObject,'String'));
if isempty(dESm_Lock); set(hObject,'String','Auto'); end

% --- Executes on button press in editSearch.
function editSearch_Callback(hObject,~,~)
global dEWin_Lock
dEWin_Lock=str2num(get(hObject,'String'));
if isempty(dEWin_Lock); set(hObject,'String','Auto'); end
