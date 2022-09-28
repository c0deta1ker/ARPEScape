function AA=DiffFlt2(A,scheme,sharp,dX2dY)
% AA=DiffFlt2(A,scheme,sharp [,dX2dY]) applies a curvature filter to sharpen the structure in the 2D array A.
% Inputs: 
% scheme='1D' or '2D' selects 1D derivation along columns or 2D derivaton
% sharp from 0 to 1 tunes the sharpening essentially from smooth f'' to singular 1/f' 
% dX2dY (if scheme='2D') is the step-width ratio between rows and columns    
% Algorithm: Zhang et al, "A precise method for visualizing dispersive features in image plots," Revi. Sci. Instr. 82 (2011) 043712  

% check input
if ~isequal(scheme,'1D')&&~isequal(scheme,'2D'); disp('Input error in CurvFlt2: scheme must be ''1D'' or ''2D'); AA=[]; return; end
if sharp<0||sharp>1; disp('Input error in CurvFlt2: sharp must be between 0 and 1'); AA=[]; return; end

% sharpening constant
c=1/(sharp+eps)-1; 

% 1D derivation scheme
if isequal(scheme,'1D')
% - normalize the derivative   
   A1Sq=DifC(A,1).^2; A1Sq=A1Sq./mean(A1Sq(~isnan(A1Sq)));
% - calculate curvature      
   AA=DifC(A,2)./(c+A1Sq).^(3/2);   
% 2D derivation scheme (neglecting cross-derivatives in the exact curvature)
else
% - normalize the partial derivatives along X and Y
   A1XSq=DifC(A',1).^2; A1XSq=A1XSq./mean(A1XSq(~isnan(A1XSq)));
   A1YSq=DifC(A,1).^2; A1YSq=A1YSq./mean(A1YSq(~isnan(A1YSq)));
% - calculate curvature
   AA=(DifC(A',2)./(c+A1XSq).^(3/2))'+dX2dY^2*DifC(A,2)./(c+A1YSq).^(3/2);  
end
