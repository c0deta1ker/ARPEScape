% -------------------------------------------------------------------------
% FUNCTION CALL:    Img = CurvatureFlt2(Img [,order] [,CX] [,CY]) 
% -------------------------------------------------------------------------
% DESCRIPTION:
% This function computes the 1D or the 2D curvature as presented in [1] 
% of the image. 
% -------------------------------------------------------------------------
% INPUT:
% Img     inpute image
% order   dim parameter [String]: '1D', '2D' (1D or 2D curvature method)
% CX      weighting parameter x direction (in 2D: exact parameter defined in [1])
% CY      weighting parameter y direction (in 2D: exact parameter defined in [1])
% [1]     P. P. Zhang, "A precise method for visualizing dispersive features
%         in image plots," Review of Scientific Instruments, vol. 82, 
%         no. 4, 2011. 
% -------------------------------------------------------------------------
% OUTPUT:
% Img     filtered image (image is normalized to value range between (0,1))
% -------------------------------------------------------------------------
function [img] = CurvatureFlt2(img,varargin)
    if(isempty(varargin))
        img = curvature2D(img, 1,1);
    else
        switch varargin{1}
            case '1D'
                img = curvature1D(img,varargin{2},varargin{3});
            case '2D'
                img = curvature2D(img,varargin{2},varargin{3});     
        end
    end
end
function AA = curvature2D(A,CX,CY)
%compute 1st,2nd derivative along x,y

    AX  = DifC(A',1)';
    AXX = DifC(A',2)';
    AY  = (DifC(A,1));
    AYY = (DifC(A,2));
    AXY = (DifC(AX,1));

%compute curvature

    AA = subplus(-(...
        ((1+CX*AX.^2).*CY.*AYY)   -...
        (2*CX*CY*(AX.*AY.*AXY))    +...
        ((1+CY*AY.^2)*CX.*AXX))   ./...
        (...
        (1 + CX*AX.^2 +CY*AY.^2).^(1.5)));
    AA              = mat2gray(AA);
end

function AA = curvature1D(A,hw,vw)
    AX  = DifC(A',1)';
    AXX = DifC(A',2)';
    AY  = (DifC(A,1));
    AYY = (DifC(A,2));
    
    AA  = hw*AXX./(1+AX.^2).^(3/2) + vw*AYY./(1+AY.^2).^(3/2);
    AA  = subplus(-AA);
    
    AA  = mat2gray(AA);
end
