function [Angle,Energy,Scan,Data,ep,Note,RawTable]=ReadARPES(file)
% [Angle,Energy,Scan,Data,ep,Note]=ReadARPES(file) reads in HDF5 and SP2 files 
% to return Angle (row), Energy (column), Scan (empty for 2D arrays and row
% for 3D arrays, Data (2- or 3-dim array), ep (pass energy) and Note (info)
% Ver. 17.08.2022

disp(['- Reading in ' file])

RawTable = [];
if strcmp(file(end-1:end),'h5')
% read in HDF5
   [Data,Axes,Note]=ReaderHDF5(file);
   try; RawTable = load(file, 'SmartTable', '-mat').SmartTable.rawData; end
   Angle=Axes{1};
   Energy=(Axes{2})';
   if size(Axes,1)==3; Scan=Axes{3}; else Scan=[]; end
   if ndims(Data)==3; Data=double(permute(Data,[2 1 3])); else Data=double(Data'); end
%   pos=strfind(Note,'Epass'); ep=str2double(Note(pos+9:pos+12));
else
% read in SP2
   [data,NoteCell]=ReaderSP2(file);
   Data=data.image; Data=double(Data');
   Angle=linspace(data.arange(1),data.arange(2),size(Data,2));
   Energy=(linspace(data.erange(1),data.erange(2),size(Data,1)))';
   Scan=[];
%   ep=str2double(NoteCell{10}(22:27));
% convert Note from cell to string
   Note = ''; for i=1:size(NoteCell)
     Note = sprintf('%s\n%s', Note, NoteCell{i});
   end
end   
% ep evaluation
pos=strfind(Note,'Ep                 =');
% - sp2 or sp2/hdf5 converted file 
if ~isempty(pos)
   hash=strfind(Note(pos:end),'#');
   ep=str2double(Note(pos+20:pos+hash(1)-2));
else
% - hdf5 file   
   pos1=strfind(Note,'Epass   =');
   pos2=strfind(Note,'Ek/Eb   =');
   ep=str2double(Note(pos1+9:pos2-1));
end