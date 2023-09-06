function AA=Gaco1(A,hw,hs)
% AA=Gaco1(A,hw [,hs]) convolutes vector or matrix A with Gaussian of 
% half-width hw and half-support hs(in inter-point distances). If skipped,
% the default hs=3*hw. The ends are padded symmetrically. If A is matrix, 
% the convolution proceeds along the columns. 

% transpose the row vector input
if size(A,1)==1; A=A'; tra=1; else tra=0; end

% default support
if nargin<3; hs=3*hw; end; hs=ceil(hs);

% convolution kernel
% - formation
K=(-1*hs:1:hs)'; K=exp(-1*(sqrt(log(2))*K/(hw+eps)).^2);
% - normalization
K=K/sum(K);

% padding A
% ALow=flipud(A(2:hs+1,:)); AHi=flipud(A(end-hs:end-1,:)); A=[ALow;A;AHi];
A=padarray(A,[hs 0],'replicate','both');

% convolution
AA=conv2(K,1,A,'valid');

% back transpose when the row vector input
if tra==1; AA=AA'; end