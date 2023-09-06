function dataStr = data_filter3D(dataStr, filType, filArgs)
% dataStr = filter_data(dataStr, filter_args)
%   This function filters the 3D ARPES data given the filter arguments, which
%   includes the filter type ("Gaco2", "GaussFlt2", "LaplaceFlt2", "CurvatureFlt2") and filter parameters.
%
%   REQ. FUNCTIONS:
%   -   [xField, yField, zField, dField] = find_data_fields(dataStr);
%   -   AA = Gaco2(A,hwX,hwY [,hsX] [,hsY])
%   -   AA = GaussFlt2(Img,hwX,hwY,hsX,hsY)
%   -   AA = CurvatureFlt2(Img [,order] [,CX] [,CY]) 
%   -   AA = LaplaceFlt2(Img [,y2xRatio] [,order])
%
%   IN:
%   -   dataStr      	loaded MATLAB data structure.
%   -   filType:        string of the type of filtering to use. Default: "spline" ("none","Gaco2","GaussFlt2","LaplaceFlt2","CurvatureFlt2").
%   -   filArgs:       	cell-array of the filtering / smoothing arguments:
%                            -> gaco2:          2x1     {hwX, hwY}              : {1, 1}
%                            -> gaussflt2:      4x1     {hwX, hwY, hsX, hsY}    : {1, 1, 3, 3}
%                            -> laplaceflt2:    2x1     {y2xRatio, order}       : {1, '4th'}
%                            -> curvatureflt2:  1x1     {order, CX, CY}         : {'2D', 2, 2}
%
%   OUT:
%   -   dataStr         modified and filtered ARPES data structure.

%% Default parameters
if nargin < 2; filType = "gaco2"; end
if nargin < 3; filArgs = {1,1}; end
if isempty(filType); filType = "gaco2"; end
if isempty(filArgs);  filArgs = {1,1}; end
%% Validity checks on the input parameters
filType = string(filType);
% - Extracting the fields to be used with most recent processing
[~, ~, ~, dField] = find_data_fields(dataStr);
%% Appending to meta information
dataStr.meta.filType = filType;
dataStr.meta.filArgs = filArgs;

%% 1 - Performing the filtering operation over all scans
DATA = dataStr.(dField); DATA_fil = [];
for i = 1:size(DATA, 3)
    if strcmpi(filType,"none") || strcmpi(filType,"")
        filtered_data = DATA(:,:,i);
    elseif strcmpi(filType,"gaco2")
        filtered_data = Gaco2(DATA(:,:,i), filArgs{:}); 
    elseif strcmpi(filType,"gaussflt2") || strcmpi(filType,"gauss2")
        filtered_data = GaussFlt2(DATA(:,:,i), filArgs{:}); 
    elseif strcmpi(filType,"laplaceflt2") || strcmpi(filType,"laplace2") || strcmpi(filType,"l2")
        filtered_data = GaussFlt2(DATA(:,:,i), 5, 5, 20, 20);
        filtered_data = LaplaceFlt2(filtered_data, filArgs{:});
    elseif strcmpi(filType,"curvatureflt2") || strcmpi(filType,"curvature2") || strcmpi(filType,"c2")
        filtered_data = GaussFlt2(DATA(:,:,i), 5, 5, 20, 20);
        filtered_data = CurvatureFlt2(filtered_data, filArgs{:});
    end
    DATA_fil(:,:,i) = filtered_data;
end
%% 2 - Setting NaN values to zero
DATA_fil(isnan(DATA_fil)) = 0;
dataStr.(dField) = DATA_fil;

end
