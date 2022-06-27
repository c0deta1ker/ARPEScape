function [XCut,DCut]=Cut(ACorr,ECorr,Data,xMode,Win)
% [XCut,DCut]=Cut(ACorr,ECorr,Data,xMode,Win) calculates 1D array DCut of intensity as
% a function of 1D array XCut(raw vector of angles/momenta) for xMode='mdc', or XCut
% (column vector of energies in Energy dimension) for xMode='edc'. The input 2D array 
% Data of intensity is defined on the 1D (row) or 2D array ACorr of [warping corrected] 
% angles/momenta and the 1D (column) or 2D array ECorr of curvature corrected [aligned] 
% energies. DCut is integrated within the Win window [win1 win2] in energies (xMode='mdc') 
% or angles/momenta (xMode='edc'). XCut is calculated in the middle of Win.
% Ver. 01.01.2021

% Revisions to Ver. 29.12.2017: Grossly reworked to correct the dimensions and linear/subscript indexing of the arrays meeting in the interpolation

% disp('- Data cut formation')

% check inputs
% - supported modes
xMode=lower(xMode);
if ~isequal(xMode,'mdc') && ~isequal(xMode,'edc')
      XCut=[]; DCut=[]; disp('Error: Only ''mdc'' and ''edc'' modes supported'); return;
end
% - integration window and its check
Win=sort(Win); 
Win(1)=Win(1)+1e-10*diff(Win); Win(2)=Win(2)-1e-10*diff(Win); % to avoid the edge effects in the code
if isequal(xMode,'mdc') Range=[max(ECorr(1,:)) min(ECorr(end,:))]; 
else Range=[max(ACorr(:,1)) min(ACorr(:,end))]; end
if Win(1)<Range(1)||Win(2)>Range(2) 
    XCut=[]; DCut=[]; disp('Error: Inconsistent integration window'); return; 
end

% shaping the arrays
% - remove NaNs spoiling the cumulative sum
Data(isnan(Data))=0;
% - expanding ACorr and ECorr if 1D arrays
if size(ACorr,1)==1; ACorr=repmat(ACorr,size(Data,1),1); end
if size(ECorr,2)==1; ECorr=repmat(ECorr,1,size(Data,2)); end
% - array permutation to reduce the 'edc' mode to 'mdc'
if isequal(xMode,'edc')
    ACorr=ACorr'; ECorr=ECorr'; Data=Data';
    T=ACorr; ACorr=ECorr; ECorr=T; clear T;
end
% - shifted energy arrays for locating the window
ECorr1=ECorr(1:end-1,:); ECorr2=ECorr(2:end,:);

% parameters
de=mean(mean(diff(ECorr,1)));

% array of intensity values DCut
% cumulative array
ISum=cumsum(Data(1:end,:),1)-0.5*Data(1:end,:)-0.5*Data(1,:);
% - point 1
% - - subscripts and linear indices of the grid point below Win(1)
[I,J]=find(ECorr1<=Win(1)&ECorr2>Win(1)); Lin1=sub2ind(size(ECorr),I,J); 
% - - linear indices of the grid point above Win(1)
Lin2=sub2ind(size(ECorr),I+1,J);
% - - interpolation of the linear matrices defined by linear indices
I1=ISum(Lin1)+(ISum(Lin2)-ISum(Lin1)).*(Win(1)-ECorr(Lin1))/de; 
% - point 2
[I,J]=find(ECorr1<=Win(2)&ECorr2>Win(2)); Lin1=sub2ind(size(ECorr),I,J);
Lin2=sub2ind(size(ECorr),I+1,J);
I2=ISum(Lin1)+(ISum(Lin2)-ISum(Lin1)).*(Win(2)-ECorr(Lin1))/de;
% integral value
DCut=I2-I1;
% normalization
DCut=DCut*de/(eps+diff(Win));

% array of coordinate values XCut interpolated in the middle of Win
% - this block runs ~25% of the integration block runtime
midE=mean(Win);
[I,J]=find(ECorr1<=midE&ECorr2>midE); 
Lin1=sub2ind(size(ECorr),I,J); Lin2=sub2ind(size(ECorr),I+1,J);
XCut=ACorr(Lin1)+(ACorr(Lin2)-ACorr(Lin1)).*(midE-ECorr(Lin1))/de;

% transpose if the MDC mode
if isequal(xMode,'mdc'); DCut=DCut'; XCut=XCut'; end