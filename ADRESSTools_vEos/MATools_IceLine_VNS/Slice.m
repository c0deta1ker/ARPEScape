function [DSlice,XSlice]=Slice(ACorr,ECorr,Data,xMode,Win)
% [DSlice,XSlice]=Slice(ACorr,ECorr,Data,xMode,Win) slices the 3D array Data defined as a 
% function of the [warping corrected] angles/parallel momenta ACorr and curvature corrected 
% [aligned] energies ECorr. The calculated 2D array DSlice is defined as a function of 2D 
% array XSlice(angles/momenta in Angle dimension, angles/momenta/hv in Scan dimension) for 
% xMode='IsoE' (case-insensitive) or XSlice(angles/momenta/hv in Scan dimension, energies 
% in Energy dimension) for xMode='isoK'. DSlice is integrated within the Win window 
% [win1 win2] of energies (xMode='isoE') or angles/momenta (xMode='isoK'). XSlice is 
% calculated in the middle of Win.
% Ver. 06.06.2020 (ver. 15.12.2017 + permanent ACorr/ECorr) 

disp('- Data slice formation')

% persistent DataSave ISum
persistent ACorrSave ECorrSave DataSave ISumE ISumK

% check supported modes
xMode=lower(xMode);
if ~isequal(xMode,'isoe') && ~isequal(xMode,'isok')
      XSlice=[]; DSlice=[]; disp('Error in Slice: Only ''IsoE'' and ''IsoK'' modes supported'); return;
end

% reset if new data (cumulative sum for another integration direction is persistent) 
if ~isequal(ACorr,ACorrSave) || ~isequal(ECorr,ECorrSave) || ~isequal(Data,DataSave)
   ISumE=[]; ISumK=[]; ACorrSave=ACorr; ECorrSave=ECorr; DataSave=Data;
end

% array permutation to reduce the 'isoK' mode to 'isoE'
if isequal(xMode,'isok')
   ACorr=permute(ACorr,[2 1 3]); ECorr=permute(ECorr,[2 1 3]); Data=permute(Data,[2 1 3]); 
   T=ACorr; ACorr=ECorr; ECorr=T; clear T;
end

% check integration window
Win=sort(Win);
% if Win(1)>max(max(max(ECorr))) || Win(2)<min(min(min(ECorr)))
if Win(1)<max(max(ECorr(1,:,:))) || Win(2)>min(min(ECorr(end,:,:)))
   XSlice=[]; DSlice=[]; disp('Error in Slice: Inconsistent integration window'); return;
end

% remove NaNs spoiling the cumulative sum
Data(isnan(Data))=0;

% expanding angle and energy arrays
% - 1D arrays
if size(ACorr,1)==1; ACorr=repmat(ACorr,size(Data,1),1); end
if size(ECorr,2)==1; ECorr=repmat(ECorr,1,size(Data,2)); end
% - 2D arrays
if size(ACorr,3)==1; ACorr=repmat(ACorr,[1 1 size(Data,3)]); end
if size(ECorr,3)==1; ECorr=repmat(ECorr,[1 1 size(Data,3)]); end

% shifted angle and energy arrays
ACorr1=ACorr(1:end-1,:,:); ECorr1=ECorr(1:end-1,:,:); ECorr2=ECorr(2:end,:,:);

% parameters
nA=size(Data,2); nS=size(Data,3); de=mean(mean(mean(diff(ECorr,1,1))));

% cumulative arrays
if isequal(xMode,'isoe')
   if isempty(ISumE); ISumE=cumsum(Data(1:end-1,:,:),1)-0.5*Data(1:end-1,:,:); end
   ISum=ISumE;
else
   if isempty(ISumK); ISumK=cumsum(Data(1:end-1,:,:),1)-0.5*Data(1:end-1,:,:); end
   ISum=ISumK;
end

% slice values
% point 1
Lin1=find(ECorr1<=Win(1)&ECorr2>Win(1));
% - subscript indices of the next point in energy
[I,J,K]=ind2sub(size(ECorr1),Lin1); Lin2=sub2ind(size(ECorr1),I+1,J,K);
% - interpolation of the linear matrices defined by linear indices
I1=ISum(Lin1)+(ISum(Lin2)-ISum(Lin1)).*(Win(1)-ECorr1(Lin1))/de; 
% point 2
Lin1=find(ECorr1<=Win(2)&ECorr2>Win(2));
[I,J,K]=ind2sub(size(ECorr1),Lin1); Lin2=sub2ind(size(ECorr1),I+1,J,K);
I2=ISum(Lin1)+(ISum(Lin2)-ISum(Lin1)).*(Win(2)-ECorr1(Lin1))/de;
% integral value
% size(ISum)
% size(I2)
% size(I1)
DSlice=I2-I1;
% normalization
DSlice=DSlice*de/diff(Win);
% reshape back to array
DSlice=(reshape(DSlice,nA,nS))';

% XSlice values interpolated in the middle of Win
% - this block runs ~25% of the integration block runtime
   midE=mean(Win);
   Lin1=find(ECorr1<=midE&ECorr2>midE);
   [I,J,K]=ind2sub(size(ECorr1),Lin1); Lin2=sub2ind(size(ECorr1),I+1,J,K);
   XSlice=ACorr1(Lin1)+(ACorr1(Lin2)-ACorr1(Lin1)).*(midE-ECorr1(Lin1))/de;
% - reshape back to array
   XSlice=(reshape(XSlice,nA,nS))';

% permutate back if the 'IsoK' mode
if isequal(xMode,'isok')
   XSlice=XSlice'; DSlice=DSlice';
end