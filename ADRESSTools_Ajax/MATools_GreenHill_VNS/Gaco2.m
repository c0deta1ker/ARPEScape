function AA=Gaco2(A,hwX,hwY,hsX,hsY)
% AA=Gaco2(A,hwX,hwY [,hsX] [,hsY]) convolves vector or matrix A with Gaussian of half-width/  
% half-support along the rows hwX/hsX (in inter-point distances > 0) and along the columns hwY/hsY. 
% If empty or skipped, the default hsX=3*hwX and hsY=3*hwY. Setting hwX or hwY to zero disables 
% convolution along the corresponding direction. The input A is padded by the end values. 
% Gaco2 is equivalent to Gaco1 except that it convolves along both dimensions.
% Ver. 23.12.2020 

% Update to Ver. 17.06.2017: 2D convolution changed to two 1D convolutions to keep the dimensionality at ~50% increase of runtime 

% default supports
if nargin<4||isempty(hsX); hsX=3*hwX; end; hsX=round(hsX);
if nargin<5||isempty(hsY); hsY=3*hwY; end; hsY=round(hsY);

% convolution along rows
if hwX~=0
% - normalized kernel   
   KX=-1*hsX:hsX; KX=exp(-1*(sqrt(log(2))*KX/(hwX)).^2); KX=KX/sum(KX);
% - padding
   A=padarray(A,[0 hsX],'replicate','both');
% - convolution
   A=conv2(1,KX,A,'valid');
end

% convolution along columns
if hwY~=0
% - normalized kernel   
   KY=(-1*hsY:hsY)'; KY=exp(-1*(sqrt(log(2))*KY/(hwY)).^2); KY=KY/sum(KY);
% - padding
   A=padarray(A,[hsY 0],'replicate','both');
% - convolution   
   A=conv2(KY,1,A,'valid');
end

% output
AA=A;

% % Ver. 2017 (malfunctioning: changing dimensionality)
% % normalized convolution kernels
% % - along rows
% if hwX==0; hsX=0; KX=1; else
%    KX=(-1*hsX:hsX)'; KX=exp(-1*(sqrt(log(2))*KX/(hwX)).^2); 
%    KX=KX/sum(KX);
% end
% % - along columns
% if hwY==0; hsY=0; KY=1; else
%    KY=-1*hsY:hsY; KY=exp(-1*(sqrt(log(2))*KY/(hwY)).^2); KY=KY/sum(KY);
% end
% 
% % padding A
% A=padarray(A,[hsY hsX],'replicate','both');
% 
% % convolution
% AA=conv2(KY,KX,A,'valid');