% -------------------------------------------------------------------------
% FUNCTION CALL:      LaplaceFlt2(Img [,y2xRatio] [,order])
% -------------------------------------------------------------------------
% DESCRIPTION:
% For laplacian computation a 2nd and 4th order scheme is used, for details
% please see (https://en.wikipedia.org/wiki/Finite_difference_coefficient).
% The 2nd order scheme needs a x_(i-1) and x_(i+1) for computation of the
% 2nd derivative at x_i. While the 4th order scheme needs x_(i-2),...,
% x_(i+2) values for x_i computation. For boundary pixels padding is 
% therefore needed. Padding is done by mirroring boundary pixels.
% -------------------------------------------------------------------------
% EXAMPLE USE:
% img = imread('cameraman.tif');
% img = Laplacian2(img,1,'4th');
% imshow(img)
% -------------------------------------------------------------------------
% INPUT:
% Img       ---      image (NxM array)
% y2xRatio  ---      weigthing ratio between y and x laplacian
% order     ---      laplacian order [String]; '2nd','4th'
% -------------------------------------------------------------------------
% OUTPUT:
% Img       ---      output image (NxM array), image is normalized between
%                    (0,1)
% -------------------------------------------------------------------------

function [img] = LaplaceFlt2(img, varargin)
    if(isempty(varargin))
        img = laplacian2nd(img,1);       
    elseif(length(varargin)==2)
        y2x     = varargin{1};
        
        switch varargin{2}
            case '2nd'
                img = laplacian2nd(img, y2x);
                
            case '4th'
                img = laplacian4th(img, y2x);
        end 
    end

end

function [img] = laplacian2nd(img, y2x)
     a = 1;
     b = -2;
     
     %resize array
     
     [nsz, msz]              = size(img);
     AA                      = zeros(nsz+2,msz+2);
     AA(2:end-1,2:end-1)     = img;
     AA(1,2:end-1)           = (img(1,:));
     AA(end,2:end-1)         = (img(end,:));
     AA(:,1)                 = (AA(:,2));
     AA(:,end)               = (AA(:,end-1));
     
     img                       = AA;
     
     % horizontal
     
     Am = b*img(2:end-1,2:end-1);
     Al = a*img(2:end-1,1:end-2);
     Ar = a*img(2:end-1,3:end);
     
     Ah = subplus(-(Am + Al + Ar));
     % vertical
     
     Am = b*img(2:end-1,2:end-1);
     Al = a*img(1:end-2,2:end-1);
     Ar = a*img(3:end,2:end-1);
     
     Av = subplus(-(Am + Al + Ar));
     
%      % diag-1
%      
%      Am = b*img(2:end-1,2:end-1);
%      Al = a*img(1:end-2,1:end-2);
%      Ar = a*img(3:end,3:end);
%      
%      Ad = dw1*subplus(-(Am + Al + Ar)/2);
%      
%      
%      % diag-2
%      
%      Am = b*img(2:end-1,2:end-1);
%      Al = a*img(3:end,1:end-2);
%      Ar = a*img(1:end-2,3:end);
%      
%      Add = dw2*subplus(-(Am + Al + Ar)/2);
     
     Av  = y2x*mat2gray(Av);
     Ah  = mat2gray(Ah);
%      Add = mat2gray(Add);
%      Ad  = mat2gray(Ad);
%      img = max(max(Av,Ah),max(Add,Ad));
    img = max(Av,Ah);
end

function [img] = laplacian4th(A,y2x)
    %coeff 4th order

    a = -1/12;
    b = 4/3;
    c = -5/2;

    %enlarge array
    [nsz, msz]               = size(A);
    AA                      = zeros(nsz+4,msz+4);
    AA(3:end-2,3:end-2)     = A;
    AA(1:2,3:end-2)         = repmat(A(1,:),[2 1]);
    AA(end-1:end,3:end-2)   = repmat(A(end,:),[2 1]);
    AA(:,1:2)               = repmat(AA(:,3),[1 2]);
    AA(:,end-1:end)         = repmat(AA(:,end-2),[1 2]);

    A                       = AA;
    % horizontal
    Am      = c*A(3:end-2,3:end-2);
    All     = a*A(3:end-2,1:end-4);
    Al      = b*A(3:end-2,2:end-3);
    Arr     = a*A(3:end-2,5:end);
    Ar      = b*A(3:end-2,4:end-1);

    Axx     = subplus(-(All + Al + Am + Ar + Arr));

    % vertical
    Am      = c*A(3:end-2,3:end-2);
    All     = a*A(1:end-4,3:end-2);
    Al      = b*A(2:end-3,3:end-2);
    Arr     = a*A(5:end,3:end-2);
    Ar      = b*A(4:end-1,3:end-2);

    Ayy     = subplus(-(All + Al + Am + Ar + Arr));

%     % diag-1
%     Am      = c*A(3:end-2,3:end-2);
%     All     = a*A(1:end-4,1:end-4);
%     Al      = b*A(2:end-3,2:end-3);
%     Arr     = a*A(5:end,5:end);
%     Ar      = b*A(4:end-1,4:end-1);
% 
%     Aud     = dw1 * subplus(-(All + Al + Am + Ar + Arr)/2);
% 
%     % diag-2
%     Am      = c*A(3:end-2,3:end-2);
%     All     = a*A(5:end,1:end-4);
%     Al      = b*A(4:end-1,2:end-3);
%     Arr     = a*A(1:end-4,5:end);
%     Ar      = b*A(2:end-3,4:end-1);
% 
%     Abu     = dw2 * subplus(-(All + Al + Am + Ar + Arr)/2);

    Axx     = mat2gray(Axx);
    Ayy     = y2x*mat2gray(Ayy);
%     Aud     = mat2gray(Aud);
%     Abu     = mat2gray(Abu);
%     img     = max(max(Axx,Ayy),max(Aud,Abu));
    img     = max(Axx,Ayy);
end
