function varargout = RawTableDisplay(varargin)
% Created by X. Wang, 02-Feb-2021

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @RawTableDisplay_OpeningFcn, ...
                   'gui_OutputFcn',  @RawTableDisplay_OutputFcn, ...
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

% --- Executes just before RawTableDisplay is made visible.
function RawTableDisplay_OpeningFcn(hObject,~, handles, varargin)
% This function has no output args, see OutputFcn.
% varargin   command line arguments to RawTableDisplay (see VARARGIN)
% Set table properties
emptyRow = TableColumns.emptyRow(); numberColumns = TableColumns.Activate - 1; set(handles.tableRaw, 'Data', emptyRow(1:numberColumns));
coledit  = true(1, numberColumns); set(handles.tableRaw, 'ColumnEditable', coledit);
colname = TableColumns.columnNames(); set(handles.tableRaw, 'ColumnName', colname(1:numberColumns));
colformat = TableColumns.columnFormats(); set(handles.tableRaw, 'ColumnFormat', colformat(1:numberColumns));
colwidth = TableColumns.columnWidth(); colwidth{TableColumns.File} = 0; colwidth{TableColumns.Comment} = 0; 
colwidth{TableColumns.Sweeps} = 50; colwidth{TableColumns.Y} = 'auto'; set(handles.tableRaw, 'ColumnWidth', colwidth);
% Find main GUI from input argument
rawData = varargin{1};
set(handles.tableRaw, 'Data', rawData);
handles.output = hObject;
% Update handles structure
guidata(hObject, handles); uiwait(hObject);

% --- Outputs from this function are returned to the command line.
function varargout = RawTableDisplay_OutputFcn(~,~,~) 
% varargout  cell array for returning output args (see VARARGOUT);
% Get default command line output from handles structure
%varargout{1} = handles.output;
