function zdat_fil = data_filter2D(zdat, filType, filXVal, filYVal)
% zdat_fil = data_filter2D(zdat, filType, filXVal, filYVal)
%   This function applies a filtering / smoothing operation on 2D data.
%
%   REQ. FUNCTIONS:
%   -   AA = Gaco2(A,hwX,hwY [,hsX] [,hsY])
%   -   AA = GaussFlt2(Img,hwX,hwY,hsX,hsY)
%   -   AA = CurvatureFlt2(Img [,order] [,CX] [,CY]) 
%   -   AA = LaplaceFlt2(Img [,y2xRatio] [,order])
%
%   IN:
%   -   zdat:           [N×M×O] array of the intensity data. Smoothes along the [N×M] axis.
%   -   filType:        string of the type of filtering to use. Default: "spline" ("none","Gaco2","GaussFlt2","LaplaceFlt2","CurvatureFlt2").
%   -   filVal:         scalar value of the x filtering / smoothing strength.
%   -   filVal:         scalar value of the y filtering / smoothing strength.
%
%   OUT:
%   -   zdat_fil:       [N×M×O] array of the filtered intensity data.

%% Default parameters
if nargin < 3; filType = "Gaco2"; end
if nargin < 4; filXVal = 1; end
if nargin < 5; filYVal = 1; end
if isempty(filType); filType = "Gaco2"; end
if isempty(filXVal); filXVal = 1; end
if isempty(filYVal); filYVal = 1; end
%% Validity checks on the input parameters
if filXVal < 0; filXVal = 0; end
if filYVal < 0; filYVal = 0; end

%% 1 - Performing the filtering operation over all scans
zdat_fil = [];
for i = 1:size(zdat, 3)
    if strcmpi(filType,"Gaco2")
        zdat_fil_01 = Gaco2(zdat(:,:,i), filXVal, filYVal); 
    elseif strcmpi(filType,"GaussFlt2")
        zdat_fil_01 = GaussFlt2(zdat(:,:,i), filXVal, filYVal, 40, 40); 
    elseif strcmpi(filType,"LaplaceFlt2")
        zdat_fil_01 = GaussFlt2(zdat(:,:,i), 5, 5, 40, 40);
        zdat_fil_01 = SetContrast(zdat_fil_01, 0.4, 0.999, 1.5);
        zdat_fil_01 = LaplaceFlt2(zdat_fil_01, (filYVal./filXVal));
    elseif strcmpi(filType,"CurvatureFlt2")
        zdat_fil_01 = GaussFlt2(zdat(:,:,i), 100, 100, 500, 500);
        zdat_fil_01 = SetContrast(zdat_fil_01, 0.4, 0.999, 1.5);
        zdat_fil_01 = CurvatureFlt2(zdat_fil_01, '2D', filXVal, filYVal);
        zdat_fil_01 = GaussFlt2(zdat_fil_01, 1, 1, 10, 10);
        zdat_fil_01 = SetContrast(zdat_fil_01, 0.2, 0.999);
    elseif strcmpi(filType,"none")
        zdat_fil_01 = zdat(:,:,i);
    end
    zdat_fil(:,:,i) = zdat_fil_01;
end
%% 2 - Setting NaN values to zero
zdat_fil(isnan(zdat_fil)) = 0;

end