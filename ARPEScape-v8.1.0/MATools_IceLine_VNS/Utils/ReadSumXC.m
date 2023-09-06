function [Angle,Energy,Scan,DataXC,ep,Note]=ReadSumXC(Files)
% ReadSumXC reads and adds separate sp2- or hdf5-files with x-correlation. 
% The files are assumed to have the same size and energy scale.
% The input parameter Files is a cell array like {'A','B','C'}.

% load files
DataXC=[];
for iFile=1:size(Files,2)
   try [Angle,Energy,Scan,Data,ep,Note]=ReadARPES(Files{iFile}); 
   catch msgbox(['Input Error in ' Files{iFile}],'modal'); return; end
   DataXC=cat(3,DataXC,Data);
end

% reduction
disp('- X-correlating data') 
maxLagE=0.25; [Energy,DataXC]=SumScanXC(Energy,DataXC,maxLagE);
