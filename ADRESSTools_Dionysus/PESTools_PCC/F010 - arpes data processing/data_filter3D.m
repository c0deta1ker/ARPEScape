function dataStr = data_filter3D(dataStr, filter_args)
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
% 	-   filter_args: 	1x2 cell of {filter_type, [xFilterVal, yFilterVal]}.
%
%   OUT:
%   -   dataStr         modified and filtered ARPES data structure.

%% Default parameters
if nargin < 2; filter_args = {"Gaco2", [1, 1]}; end
if isempty(filter_args); filter_args = {"Gaco2", [1, 1]}; end
% disp('Data filtering...')
% wbar = waitbar(0, 'Executing data filtering...', 'Name', 'data_filter3D');

%% 1 - Initialising the filter parameters
dataStr.meta.filter_args = filter_args;
% - Extracting the fields to be used with most recent processing
[~, ~, ~, dField] = find_data_fields(dataStr);
% - Extracting filter parameters
filter_type = filter_args{1};
filter_val  =  filter_args{2};
if length(filter_val) == 1; filter_val = [1, 1] .*filter_val; end

%% 2 - Performing the filtering operation over all scans
for i = 1:size(dataStr.(dField), 3)
%     waitbar(i/size(dataStr.(dField), 3), wbar, 'Filtering ARPES data...', 'Name', 'data_filter3D');
    if filter_type == "Gaco2"
        filtered_data = Gaco2(dataStr.(dField)(:,:,i), filter_val(1), filter_val(2)); 
    elseif filter_type == "GaussFlt2"
        filtered_data = GaussFlt2(dataStr.(dField)(:,:,i), filter_val(1), filter_val(2), 40, 40); 
    elseif filter_type == "LaplaceFlt2"
        filtered_data = GaussFlt2(dataStr.(dField)(:,:,i), 5, 5, 40, 40);
        filtered_data = SetContrast(filtered_data, 0.4, 0.999, 1.5);
        filtered_data = LaplaceFlt2(filtered_data, filter_val(1));
    elseif filter_type == "CurvatureFlt2"
        filtered_data = GaussFlt2(dataStr.(dField)(:,:,i), 100, 100, 500, 500);
        filtered_data = SetContrast(filtered_data, 0.4, 0.999, 1.5);
        filtered_data = CurvatureFlt2(filtered_data, '2D', filter_val(1), filter_val(2));
        filtered_data = GaussFlt2(filtered_data, 1, 1, 10, 10);
        filtered_data = SetContrast(filtered_data, 0.2, 0.999);
    end
    dataStr.(dField)(:,:,i) = filtered_data;
end
%% 3 - Setting NaN values to zero
dataStr.(dField)(isnan(dataStr.(dField))) = 0;
%% Close wait-bar
% close(wbar);
end
