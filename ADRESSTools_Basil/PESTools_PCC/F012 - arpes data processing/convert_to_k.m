function dataStr = convert_to_k(dataStr, kconv_args)
% dataStr = convert_to_k(dataStr, kconv_args)
%   This is a function that will convert the angles thtA into kx and the
%   scan parameters into either ky (for tltM) or kz (for hv). This performs
%   stage 3 processing - wavevector conversions
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
%   -   kconv_args:     1x4 cell of {eB_ref, kx_ref, thtA_ref, v000}.
%
%   OUT:
%   dataStr - MATLAB data structure with new additional fields below;
%   -   .(surfNormX):   double or vector of surface normal vector.
%   -   .(kx):          2D or 3D array of kx from the Theta angle.
%	-   .(ky):          2D or 3D array of kx from the Tilt angle.
%	-   .(kz):          2D or 3D array of kx from the Photon Energy.

%% Default parameters
if nargin < 2; kconv_args = {0, 0, 0, 12.57}; end
if isempty(kconv_args); kconv_args = {0, 0, 0, 12.57}; end

%% - 1 - Initialising the k conversion parameters
dataStr.meta.kconv_args = kconv_args;
% - Extracting conversion parameters
eB_ref      = kconv_args{1};
kx_ref      = kconv_args{2};
thtA_ref    = kconv_args{3};
v000        = kconv_args{4};
%% - 2 - Wave-vector conversions over all scans
% - 2.1 Wave-vector conversions for Eb(k)
if dataStr.Type == "Eb(k)" || dataStr.Type == "Eb(k,i)" 
    if length(thtA_ref) == 2; thtA_ref = thtA_ref(1); end
    dataStr.surfNormX = SurfNormX(dataStr.hv, eB_ref, kx_ref, dataStr.thtM, thtA_ref);
    dataStr.kx = Kxx(dataStr.hv, dataStr.eb, dataStr.thtM, dataStr.tht, dataStr.surfNormX);
    dataStr.ky = Kyy(dataStr.hv, dataStr.eb, dataStr.tht, dataStr.thtM, dataStr.tltM, dataStr.surfNormX);
    dataStr.kz = Kzz(dataStr.hv, dataStr.eb, dataStr.thtM, dataStr.tht, dataStr.tltM, v000, dataStr.surfNormX);
    % Finding the mean value of ky and kz
    dataStr.ky = mean(dataStr.ky(:));
    dataStr.kz = mean(dataStr.kz(:));
% - 2.2 Wave-vector conversions for Eb(kx,ky)
elseif dataStr.Type == "Eb(kx,ky)"
    if length(thtA_ref) == 2; thtA_ref = linspace(thtA_ref(1), thtA_ref(2), size(dataStr.tltM, 2));
    else; thtA_ref = ones(size(dataStr.tltM, 2))*thtA_ref(1); 
    end
    for i = 1:size(dataStr.tltM, 2)
        dataStr.surfNormX(i) = SurfNormX(dataStr.hv, eB_ref, kx_ref, dataStr.thtM, thtA_ref(i));
        dataStr.kx(:,:,i) = Kxx(dataStr.hv, dataStr.eb(:,:,i), dataStr.thtM, dataStr.tht(:,:,i), dataStr.surfNormX(i));
        dataStr.ky(:,:,i) = Kyy(dataStr.hv, dataStr.eb(:,:,i), dataStr.tht(:,:,i), dataStr.thtM, dataStr.tltM(i), dataStr.surfNormX(i));
        dataStr.kz(:,:,i) = Kzz(dataStr.hv, dataStr.eb(:,:,i), dataStr.thtM, dataStr.tht(:,:,i), dataStr.tltM(i), v000, dataStr.surfNormX(i));
    end
% - 2.3 Wave-vector conversions for Eb(kx,kz)
elseif dataStr.Type == "Eb(kx,kz)"
    if length(thtA_ref) == 2; thtA_ref = linspace(thtA_ref(1), thtA_ref(2), size(dataStr.hv, 2)); 
    else; thtA_ref = ones(size(dataStr.hv, 2))*thtA_ref(1); 
    end
    for i = 1:size(dataStr.hv, 2)
        dataStr.surfNormX(i) = SurfNormX(dataStr.hv(i), eB_ref, kx_ref, dataStr.thtM, thtA_ref(i));
        dataStr.kx(:,:,i) = Kxx(dataStr.hv(i), dataStr.eb(:,:,i), dataStr.thtM, dataStr.tht(:,:,i), dataStr.surfNormX(i));
        dataStr.ky(:,:,i) = Kyy(dataStr.hv(i), dataStr.eb(:,:,i), dataStr.tht(:,:,i), dataStr.thtM, dataStr.tltM, dataStr.surfNormX(i));
        dataStr.kz(:,:,i) = Kzz(dataStr.hv(i), dataStr.eb(:,:,i), dataStr.thtM, dataStr.tht(:,:,i), dataStr.tltM, v000, dataStr.surfNormX(i));
    end
end


end