function varargout = NaviSetup(varargin)
% Ver. 21Aug2022
% NAVISETUP MATLAB code for NaviSetup.fig
%      NAVISETUP, by itself, creates a new NAVISETUP or raises the existing
%      singleton*.
%      H = NAVISETUP returns the handle to a new NAVISETUP or the handle to
%      the existing singleton*.
%      NAVISETUP('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in NAVISETUP.M with the given input arguments.
%      NAVISETUP('Property','Value',...) creates a new NAVISETUP or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before NaviSetup_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to NaviSetup_OpeningFcn via varargin.
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
% Last Modified by GUIDE v2.5 21-May-2020 19:44:12
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @NaviSetup_OpeningFcn, ...
                   'gui_OutputFcn',  @NaviSetup_OutputFcn, ...
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

% --- Executes just before NaviSetup is made visible.
function NaviSetup_OpeningFcn(hObject,~,handles,varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% varargin   command line arguments to NaviSetup (see VARARGIN)
% Choose default command line output for NaviSetup
handles.output = hObject;
% Update handles structure
guidata(hObject, handles);
% figure properties
% set(gcf,'WindowStyle','modal')
% global parameters
global surfNx_Navi hv_N eB_N kx_N thtM thtM_N thtA_N alpha
set(handles.editHv,'String',sprintf('%0.3f',hv_N));
thtM_N=thtM; set(handles.editThtM,'String',sprintf('%0.3f',thtM_N));
set(handles.editEB,'String',sprintf('%0.3f',eB_N));
kx_N=0; % kx_N=[]; 
set(handles.editKx,'String',sprintf('%0.3f',kx_N));
set(handles.editThtA,'String',sprintf('%0.3f',thtA_N));
surfNx_Navi=SurfNormX(hv_N,eB_N,kx_N,thtM_N,thtA_N,alpha);
set(handles.textSurfNormX,'String',['surfNx = ' sprintf('%0.3f',surfNx_Navi)]);
% % - figure position
% global xPos_Navi yPos_Navi
% if ~isempty(xPos_Navi*yPos_Navi)
% Pos=get(gcf,'position'); offsetFigHeader=29; set(gcf,'position',[xPos_Navi yPos_Navi-Pos(4)-offsetFigHeader Pos(3) Pos(4)]);
% end

% --- Outputs from this function are returned to the command line.
function varargout = NaviSetup_OutputFcn(~,~,handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% Get default command line output from handles structure
varargout{1} = handles.output;

function editHv_Callback(hObject,~,handles)
global hv_N eB_N kx_N thtM_N thtA_N surfNx_Navi alpha
hv_N=field2num(hObject);
if length(hv_N)==1
   surfNx_Navi=SurfNormX(hv_N,eB_N,kx_N,thtM_N,thtA_N,alpha);
   set(handles.textSurfNormX,'String',['surfNx = ' sprintf('%0.3f',surfNx_Navi)]);
end
%set(hObject,'String',num2str(hv_N))

function editThtM_Callback(hObject,~,handles)
global hv_N eB_N kx_N thtM_N thtA_N surfNx_Navi alpha
thtM_N=field2num(hObject);
if length(thtM_N)==1
   surfNx_Navi=SurfNormX(hv_N,eB_N,kx_N,thtM_N,thtA_N,alpha);
   set(handles.textSurfNormX,'String',['surfNx = ' sprintf('%0.3f',surfNx_Navi)]);
end
%set(hObject,'String',num2str(thtM_N))

function editKx_Callback(hObject,~,handles)
global hv_N eB_N kx_N thtM_N thtA_N surfNx_Navi alpha
kx_N=field2num(hObject);
if length(kx_N)==1
   surfNx_Navi=SurfNormX(hv_N,eB_N,kx_N,thtM_N,thtA_N,alpha);
   set(handles.textSurfNormX,'String',['surfNx = ' sprintf('%0.3f',surfNx_Navi)]);
end
%set(hObject,'String',num2str(kx_N))

function editThtA_Callback(hObject,~,handles)
global hv_N eB_N kx_N thtM_N thtA_N surfNx_Navi alpha
thtA_N=field2num(hObject);
if length(thtA_N)==1 
   surfNx_Navi=SurfNormX(hv_N,eB_N,kx_N,thtM_N,thtA_N,alpha);
   set(handles.textSurfNormX,'String',['surfNx = ' sprintf('%0.3f',surfNx_Navi)]);
end
%set(hObject,'String',num2str(thtA_N))

function editEB_Callback(hObject,~,handles)
global hv_N eB_N kx_N thtM_N thtA_N surfNx_Navi alpha
eB_N=field2num(hObject);
if length(eB_N)==1
   surfNx_Navi=SurfNormX(hv_N,eB_N,kx_N,thtM_N,thtA_N,alpha);
   set(handles.textSurfNormX,'String',['surfNx= ' sprintf('%0.3f',surfNx_Navi)]);
end   
%set(hObject,'String',num2str(eB_N))
