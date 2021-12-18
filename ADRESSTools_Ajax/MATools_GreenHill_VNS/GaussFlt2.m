% -------------------------------------------------------------------------
% FUNCTION CALL: GaussFlt2(Img,hwX,hwY,hsX,hsY)
% -------------------------------------------------------------------------
% DESCRIPTION:
% The image is convolved with a gaussian kernel. The gaussian kernel is a 
% kernel with [hsX,hsY] size and incorporates the horicontal, vertical
% half widths. The kernel is normalized.
% -------------------------------------------------------------------------
% HINT:
% For medianfilter use the MATLAB function medfilt2 (help medfilt2)
% -------------------------------------------------------------------------
% INPUT:
% Img    2D intensity image
% hwX    horizontal half width
% hwY    vertical half width
% hsX    horizontal kernel length
% hsY    vertical kernel length
% -------------------------------------------------------------------------
% OUTPUT:
% Img    image convolved with gaussian kernel
% -------------------------------------------------------------------------
function Img = GaussFlt2(Img,hwX,hwY,hsX,hsY)
    sigx    = hwX/(2*sqrt(2*log(2)));
    sigy    = hwY/(2*sqrt(2*log(2)));
    cov     = [sigx, 0; 0 sigy];
    kernel  = getGauss(hsX,hsY,cov);
    Img     = imfilter(Img,kernel);
end

function gauss = getGauss(hsX, hsY, cov)
        halfx   = (hsX-1)/2;
        halfy   = (hsY-1)/2;
        [x,y]   = meshgrid(-halfx:halfx,-halfy:halfy);
        A       = inv(cov);
        a1      = A(1,1)*x + A(1,2)*y;
        a2      = A(2,1)*x + A(2,2)*y;
        gauss   = exp(-0.5*(x.*a1 + y.*a2));
        gauss   = gauss./sum(gauss(:));
end

