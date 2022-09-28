function varargout = NormSetup(varargin)
% NORMSETUP MATLAB code for NormSetup.fig
%      NORMSETUP, by itself, creates a new NORMSETUP or raises the existing
%      singleton*.
%      H = NORMSETUP returns the handle to a new NORMSETUP or the handle to
%      the existing singleton*.
%      NORMSETUP('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in NORMSETUP.M with the given input arguments.
%      NORMSETUP('Property','Value',...) creates a new NORMSETUP or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before NormSetup_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to NormSetup_OpeningFcn via varargin.
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
% Last Modified by GUIDE v2.5 04-Nov-2020 16:28:10
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @NormSetup_OpeningFcn, ...
                   'gui_OutputFcn',  @NormSetup_OutputFcn, ...
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

% --- Executes just before NormSetup is made visible.
function NormSetup_OpeningFcn(hObject,~,handles,varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% varargin   command line arguments to NormSetup (see VARARGIN)
% Choose default command line output for NormSetup
handles.output = hObject;
% Update handles structure
guidata(hObject, handles);
% figure properties
% set(gcf,'WindowStyle','modal')
% global parameters
% - initial values
global RangeNorm scCoeff scCoeffIni clipNormNeg
% if isempty(RangeNorm)||isempty(str2num(RangeNorm)); set(handles.editRange,'String','Full'); else set(handles.editRange,'String',sprintf('%g:%g',RangeNorm)); end
if isempty(RangeNorm); set(handles.editRange,'String','Full'); else set(handles.editRange,'String',sprintf('%g:%g',RangeNorm)); end
if isempty(scCoeff); scCoeff=scCoeffIni; end; set(handles.editSubtr,'String',num2str(scCoeff))
if isempty(clipNormNeg); clipNormNeg=1; end; set(handles.radiobuttonClipNormNeg,'Value',clipNormNeg)
% % - figure position
% global xPos_Norm yPos_Norm
% Pos=get(gcf,'position'); offsetFigHeader=29;
% set(gcf,'position',[xPos_Norm yPos_Norm-Pos(4)-offsetFigHeader Pos(3) Pos(4)]);

% --- Outputs from this function are returned to the command line.
function varargout = NormSetup_OutputFcn(~,~,handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in editRange.
function editRange_Callback(hObject,~,~)
% hObject    handle to editRange (see GCBO)
% clear normalization array
global RangeNorm IDiv ISubtr
IDiv=[]; ISubtr=[]; disp('- cleared IDiv ISubtr')
% energy range
% RangeStr=get(hObject,'String'); RangeStr=strrep(RangeStr,':',' '); Range=sort(str2num(RangeStr));
% if length(Range)==2; RangeNorm=Range; else RangeNorm=RangeIni; end
RangeNorm=field2num(hObject); if length(RangeNorm)~=2; RangeNorm=[]; set(hObject,'String','Full'); end
% set(hObject,'String',[num2str(RangeNorm(1)) ':' num2str(RangeNorm(2))]);

% --- Executes on button press in editSubtr.
function editSubtr_Callback(~,~,handles)
% hObject    handle to editSubtr (see GCBO)
global scCoeff scCoeffIni
% Scaling coefficient to subtract angle-integrated spectrum
scCoeff=field2num(handles.editSubtr,scCoeffIni);
% set(handles.editSubtr,'String',num2str(scCoeffIni));

% --- Executes on button press in radiobuttonClipNormNeg.
function radiobuttonClipNormNeg_Callback(hObject,~,~)
global clipNormNeg; clipNormNeg=get(hObject,'value');