function dataStr = convert_to_kx_only(dataStr, kconv_args, plot_result)
% dataStr = convert_to_k(dataStr, kconv_args, plot_result)
%   - #3 - THIRD STEP IN PROCESSING ARPES DATA
%   This is a function that will convert the angles thtA into kx and the
%   scan parameters into either ky (for tltM) or kz (for hv).  Last updated in
%   October 2022, and includes input for new X-ray incidence angle.
%
%   REQ. FUNCTIONS:
%   -   [xField, yField, zField, dField] = find_data_fields(dataStr);
%   -   surfNormX = SurfNormX(hv, eB_ref, kx_ref, thtM, thtA_ref) calculates the surface normal angle in the MP.
%   -   Kxx = Kxx(HV, Eb, thtM, ThtA, surfNormX) calculates k// in the MP.
%   -   Kyy = Kyy(HV, Eb, thtM, TltM, surfNormX) calculates k// in the MP.
%   -   Kzz = Kzz(HV, Eb, thtM, ThtA, TltM, v000, surfNormX) calculates kz along the surface normal.
%
%   IN:
%   -   dataStr:        data structure of the ARPES data.
%   -   kconv_args:     {1×5} cell of {eB_ref, kx_ref, thtA_ref, v000, incAlpha}.
%                               eB_ref:    	single, constant value of the binding energy at the refernece point. 
%                               kx_ref:   	single, constant value of the wave-vector value at the reference point.
%                               thtA_ref: 	single, constant value of the tht angle value at the reference point.
%                               v000:       single, constant value of the inner potential. Tune to make sure spacing of Gamma points is correct. Default is 12.57.
%                               incAlpha:   nominal X-ray incidence angle, for the present geometry equal to 9deg. If omitted, old geometry with 20deg implied.
%   -   plot_result:        if 1, will plot figure of the alignment, otherwise it wont.
%   OUT:
%   dataStr - MATLAB data structure with new additional fields below;
%   -   .(surfNormX):   double or vector of surface normal vector.
%   -   .(kx):          2D or 3D array of kx from the Theta angle.
%	-   .(ky):          2D or 3D array of kx from the Tilt angle.
%	-   .(kz):          2D or 3D array of kx from the Photon Energy.

%% Default parameters
if nargin < 3; plot_result = 0; end
if nargin < 2; kconv_args = {0, 0, 0, 12.57, 20}; end
if isempty(kconv_args); kconv_args = {0, 0, 0, 12.57, 20}; end
if isempty(plot_result); plot_result = 0; end
%% Validity checks on the input parameters
if length(kconv_args) ~= 5
    error('error: kconv_args is not a {1x5} cell - make sure you defined the arguments correctly.');
else
    if isempty(kconv_args{1}) || ischar(kconv_args{2}) || isstring(kconv_args{3}); kconv_args{1} = 0; end
    if isempty(kconv_args{2}) || ischar(kconv_args{2}) || isstring(kconv_args{2}); kconv_args{2} = 0; end
    if isempty(kconv_args{3}) || ischar(kconv_args{3}) || isstring(kconv_args{3}); kconv_args{3} = 0; end
    if isempty(kconv_args{4}) || ischar(kconv_args{4}) || isstring(kconv_args{4}); kconv_args{4} = 12.57; end
    if isempty(kconv_args{5}) || ischar(kconv_args{5}) || isstring(kconv_args{5}); kconv_args{5} = 20; end
end

%% - 1 - Initialising the k conversion parameters
dataStr.meta.kconv_args = kconv_args;
% - Extracting conversion parameters
eB_ref      = kconv_args{1};
kx_ref      = kconv_args{2};
thtA_ref    = kconv_args{3};
v000        = kconv_args{4};
incAlpha    = kconv_args{5};
%% - 2 - Wave-vector conversions over all scans
if ~isfield(dataStr,'tltE'); dataStr.tltE = 0; end
if isnan(dataStr.tltE); dataStr.tltE = 0; end
% - 2.1 Wave-vector conversions for Eb(k)
if dataStr.Type == "Eb(k)" || dataStr.Type == "Eb(k,i)" 
    if length(thtA_ref) == 2; thtA_ref = thtA_ref(1); end
    if length(dataStr.index) == 1
        dataStr.surfNormX = SurfNormX(dataStr.hv, eB_ref, kx_ref, dataStr.thtM, thtA_ref, incAlpha);
        dataStr.kx = Kxx(dataStr.hv, dataStr.eb, dataStr.thtM, dataStr.tht, dataStr.surfNormX, incAlpha);
    else
        for i = 1:length(dataStr.index)
            dataStr.surfNormX = SurfNormX(dataStr.hv, eB_ref, kx_ref, dataStr.thtM, thtA_ref, incAlpha);
            dataStr.kx(:,:,i)  = Kxx(dataStr.hv, dataStr.eb(:,:,i), dataStr.thtM, dataStr.tht(:,:,i), dataStr.surfNormX, incAlpha);
        end
    end
% - 2.2 Wave-vector conversions for Eb(kx,ky)
elseif dataStr.Type == "Eb(kx,ky)"
    if length(thtA_ref) == 2; thtA_ref = linspace(thtA_ref(1), thtA_ref(2), size(dataStr.tltM, 2));
    else; thtA_ref = ones(size(dataStr.tltM, 2))*thtA_ref(1); 
    end
    % Mechanical Tilt Map
    if size(dataStr.tltM, 2) > 1
        for i = 1:size(dataStr.tltM, 2)
            dataStr.surfNormX(i) = SurfNormX(dataStr.hv, eB_ref, kx_ref, dataStr.thtM, thtA_ref(i), incAlpha);
            dataStr.kx(:,:,i) = Kxx(dataStr.hv, dataStr.eb(:,:,i), dataStr.thtM, dataStr.tht(:,:,i), dataStr.surfNormX(i), incAlpha);
        end
    % Electrostatic Tilt Map
    else
        for i = 1:size(dataStr.tltE, 2)
            dataStr.surfNormX(i) = SurfNormX(dataStr.hv, eB_ref, kx_ref, dataStr.thtM, thtA_ref(i), incAlpha);
            dataStr.kx(:,:,i) = Kxx(dataStr.hv, dataStr.eb(:,:,i), dataStr.thtM, dataStr.tht(:,:,i), dataStr.surfNormX(i), incAlpha);
        end
    end

% - 2.3 Wave-vector conversions for Eb(kx,kz)
elseif dataStr.Type == "Eb(kx,kz)"
    if length(thtA_ref) == 2; thtA_ref = linspace(thtA_ref(1), thtA_ref(2), size(dataStr.hv, 2)); 
    else; thtA_ref = ones(size(dataStr.hv, 2))*thtA_ref(1); 
    end
    for i = 1:size(dataStr.hv, 2)
        dataStr.surfNormX(i) = SurfNormX(dataStr.hv(i), eB_ref, kx_ref, dataStr.thtM, thtA_ref(i), incAlpha);
        dataStr.kx(:,:,i) = Kxx(dataStr.hv(i), dataStr.eb(:,:,i), dataStr.thtM, dataStr.tht(:,:,i), dataStr.surfNormX(i), incAlpha);
    end
end


[xField, yField, ~, dField]             = find_data_fields(dataStr);
HV = dataStr.hv;
% - 3.1 Remapping to 2D
% -- Remapping the THT / K domain to be 2D
if size(HV,1) == 1 || size(HV,2) == 1; HV=repmat(HV,[size(dataStr.(dField),1), 1]); end
% - 3.2 Remapping to 3D
% -- Remapping the THT / K domain to be 3D
if size(HV, 3) == 1; HV = repmat(HV,[1, 1, size(dataStr.(dField),2)]); end
HV = permute(HV, [1,3,2]);
dataStr.kz = HV;

%% -- For Debugging
if plot_result == 1
    fig = figure();
    fig.Position(3) = 2*450; 
    fig.Position(4) = 350;
    subplot(121); hold on;
    ImData(dataStr.tht(:,:,1), dataStr.eb(:,:,1), dataStr.data(:,:,1));
    crosshair(mean(thtA_ref(:)), 0, 'color', 'c', 'linewidth', 1.5, 'linestyle', '-');
    img_props(1, 'Color', [1 1 1]); colormap(hot); 
    axis([min(dataStr.tht(:)), max(dataStr.tht(:)), min(dataStr.eb(:)), max(dataStr.eb(:))]);
    colorbar; title('convert_to_k(): Before', 'Interpreter', 'none');
    subplot(122); hold on;
    ImData(dataStr.kx(:,:,1), dataStr.eb(:,:,1), dataStr.data(:,:,1));
    img_props(1, 'Color', [1 1 1]); colormap(hot); 
    axis([min(dataStr.kx(:)), max(dataStr.kx(:)), min(dataStr.eb(:)), max(dataStr.eb(:))]);
    colorbar; title('convert_to_k(): After', 'Interpreter', 'none');
end
end