function varargout = ImgProcSetup(varargin)
% IMGPROCSETUP MATLAB code for ImgProcSetup.fig
%      IMGPROCSETUP, by itself, creates a new IMGPROCSETUP or raises the existing
%      singleton*.
%      H = IMGPROCSETUP returns the handle to a new IMGPROCSETUP or the handle to
%      the existing singleton*.
%      IMGPROCSETUP('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in IMGPROCSETUP.M with the given input arguments.
%      IMGPROCSETUP('Property','Value',...) creates a new IMGPROCSETUP or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ImgProcSetup_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ImgProcSetup_OpeningFcn via varargin.
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
% Last Modified by GUIDE v2.5 04-Nov-2020 16:11:59
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ImgProcSetup_OpeningFcn, ...
                   'gui_OutputFcn',  @ImgProcSetup_OutputFcn, ...
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

% --- Executes just before ImgProcSetup is made visible.
function ImgProcSetup_OpeningFcn(hObject,~,handles,varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% varargin   command line arguments to ImgProcSetup (see VARARGIN)
% Choose default command line output for ImgProcSetup
handles.output = hObject;
% Update handles structure
guidata(hObject, handles);
% figure properties
% set(gcf,'WindowStyle','modal')
% global parameters and their initialization
global presmdX presmdY structMode sharp sharpIni postsmdX postsmdY clipImgProcNeg
if isempty(presmdX); presmdX=0; end; set(handles.editPresmdX,'String',sprintf('%g',presmdX)); 
if isempty(presmdY); presmdY=0; end; set(handles.editPresmdY,'String',sprintf('%g',presmdY));
if isempty(structMode); structMode=1; end; set(handles.menuStructMode,'Value',structMode)
if structMode==1; set([handles.textSharp handles.editSharp],'Visible','Off'); else set([handles.textSharp handles.editSharp],'Visible','On'); end
sharpIni=0; if isempty(sharp); sharp=sharpIni; end; set(handles.editSharp,'String',sprintf('%g',sharp));
if isempty(postsmdX); postsmdX=0; end; set(handles.editPostsmdX,'String',sprintf('%g',postsmdX)); 
if isempty(postsmdY); postsmdY=0; end; set(handles.editPostsmdY,'String',sprintf('%g',postsmdY));
if isempty(clipImgProcNeg); clipImgProcNeg=1; end; set(handles.radiobuttonClipImgProcNeg,'Value',clipImgProcNeg);
%if modeView==1; set(handles.menuStructMode,'String',{'None';'Curvature' ; '-d2I/dE2>0'}); else; set(handles.menuStructMode,'String',{'None';'Curvature'}); end
% % - figure position
% global xPos_View yPos_View
% Pos=get(gcf,'position'); offsetFigHeader=29;
% set(gcf,'position',[xPos_View yPos_View-Pos(4)-offsetFigHeader Pos(3) Pos(4)]);

% --- Outputs from this function are returned to the command line.
function varargout = ImgProcSetup_OutputFcn(~,~,handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in editPresmdX.
function editPresmdX_Callback(hObject,~,~)
global presmdX; presmdX=field2num(hObject,0);
% --- Executes on button press in editPresmdY.
function editPresmdY_Callback(hObject,~,~)
global presmdY; presmdY=field2num(hObject,0);

% --- Executes on selection change in menuStructMode.
function menuStructMode_Callback(hObject,~,handles)
global structMode; structMode=get(hObject,'Value');
if structMode==1; set([handles.textSharp handles.editSharp],'Visible','Off'); else set([handles.textSharp handles.editSharp],'Visible','On'); end

function editSharp_Callback(hObject, ~,~)
global sharp sharpIni; sharp=field2num(hObject,sharpIni);

% --- Executes on button press in editPostsmdX.
function editPostsmdX_Callback(hObject,~,~)
global postsmdX; postsmdX=field2num(hObject,0);
% --- Executes on button press in editPresmdY.
function editPostsmdY_Callback(hObject,~,~)
global postsmdY; postsmdY=field2num(hObject,0);

% --- Executes on button press in radiobuttonClipImgProcNeg.
function radiobuttonClipImgProcNeg_Callback(hObject,~,~)
global clipImgProcNeg; clipImgProcNeg=get(hObject,'value');

% Hints: get(hObject,'String') returns contents of editSharp as text
%        str2double(get(hObject,'String')) returns contents of editSharp as a double

% % --- Executes on button press in editXRange.
% function editXRange_Callback(hObject,~,~)
% global XRange
% XRange=field2num(hObject); if length(XRange)~=2; XRange=[]; set(hObject,'String','Full'); end
% 
% % --- Executes on button press in editYRange.
% function editYRange_Callback(hObject,~,~)
% global YRange
% YRange=field2num(hObject); if length(YRange)~=2; YRange=[]; set(hObject,'String','Full'); end

% Hint: get(hObject,'Value') returns toggle state of radiobuttonClipImgProcNeg
