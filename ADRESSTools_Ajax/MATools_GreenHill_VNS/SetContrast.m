function DataC=SetContrast(Data,minFrac,maxFrac,gamma)
% DataC=SetContrast(Data,minFrac,maxFrac [,gamma]) increases the image contrast 
% by saturating the lowest and highest data values within the minFrac and maxFrac
% fractions of all data points and (optionally) introducing the gamma-correction.
% The function ignores NaN elements. DataC is offset to zero in its minimum.
% Ver. 11.02.2018

% inputs
if nargin<4; gamma=1; end

% normalization to vary between 0 and 1 required for image processing functions
dMin=min(min(Data)); Data=Data-dMin; scale=max(max(Data)); Data=Data/scale;

% reference array without NaNs 
DataR=Data(:); DataR=DataR(~isnan(DataR));
if isempty(DataR); DataC=Data; return; end

% histogram enhancement; the output arrays varies between 0 and 1
Lims=stretchlim(DataR,[minFrac maxFrac]); 
DataC=imadjust(Data,Lims,[0 1],gamma);

% recovery of the NaNs removed by imadjust
DataC(isnan(Data))=NaN;

% scaling back
DataC=scale*DataC/diff(Lims)+dMin;