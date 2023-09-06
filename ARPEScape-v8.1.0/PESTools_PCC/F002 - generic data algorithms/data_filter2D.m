function zdat_fil = data_filter2D(zdat, filType, filArgs, plot_result)
% zdat_fil = data_filter2D(zdat, filType, filArgs, plot_result)
%   This function applies a filtering / smoothing operation on 2D data.
%
%   REQ. FUNCTIONS:
%   -   AA = Gaco2(A,hwX,hwY [,hsX] [,hsY])
%   -   AA = GaussFlt2(Img,hwX,hwY,hsX,hsY)
%   -   AA = LaplaceFlt2(Img [,y2xRatio] [,order])
%   -   AA = CurvatureFlt2(Img [,order] [,CX] [,CY]) 
%
%   IN:
%   -   zdat:           [nX×nY] array of the intensity data. Smoothes along the [N×M] axis.
%   -   filType:        string of the type of filtering to use. Default: "spline" ("none","Gaco2","GaussFlt2","LaplaceFlt2","CurvatureFlt2").
%   -   filArgs:       	cell-array of the filtering / smoothing arguments:
%                            -> gaco2:          2x1     {hwX, hwY}              : {1, 1}
%                            -> gaussflt2:      4x1     {hwX, hwY, hsX, hsY}    : {1, 1, 3, 3}
%                            -> laplaceflt2:    2x1     {y2xRatio, order}       : {1, '4th'}
%                            -> curvatureflt2:  1x1     {order, CX, CY}         : {'2D', 2, 2}
%   -   plot_result:    if 1, will plot figure summary, otherwise it wont.
%
%   OUT:
%   -   zdat_fil:       [nX×nY] array of the filtered intensity data.

%% Default parameters
if nargin < 2; filType = "gaco2"; end
if nargin < 3; filArgs = {1,1}; end
if nargin < 4; plot_result=0; end
if isempty(filType); filType = "gaco2"; end
if isempty(filArgs);  filArgs = {1,1}; end
if isempty(plot_result); plot_result=0;  end
%% Validity checks on the input parameters
filType = string(filType);

%% 1 - Performing the filtering operation over all scans
zdat_fil = [];
for i = 1:size(zdat, 3)
    if strcmpi(filType,"none") || strcmpi(filType,"")
        zdat_fil_01 = zdat(:,:,i);
    elseif strcmpi(filType,"gaco2")
        zdat_fil_01 = Gaco2(zdat(:,:,i), filArgs{:}); 
    elseif strcmpi(filType,"gaussflt2") || strcmpi(filType,"gauss2")
        zdat_fil_01 = GaussFlt2(zdat(:,:,i), filArgs{:}); 
    elseif strcmpi(filType,"laplaceflt2") || strcmpi(filType,"laplace2") || strcmpi(filType,"l2")
        zdat_fil_01 = GaussFlt2(zdat(:,:,i), 5, 5, 20, 20);
        zdat_fil_01 = LaplaceFlt2(zdat_fil_01, filArgs{:});
    elseif strcmpi(filType,"curvatureflt2") || strcmpi(filType,"curvature2") || strcmpi(filType,"c2")
        zdat_fil_01 = GaussFlt2(zdat(:,:,i), 5, 5, 20, 20);
        zdat_fil_01 = CurvatureFlt2(zdat_fil_01, filArgs{:});
    end
    zdat_fil(:,:,i) = zdat_fil_01;
end
%% 2 - Setting NaN values to zero
zdat_fil(isnan(zdat_fil)) = 0;

%% -- For Debugging
if plot_result == 1
    [xdat,ydat] = meshgrid(1:size(zdat,2), 1:size(zdat,1));
    fig = figure(); 
    fig.Position(3) = 400*2; 
    fig.Position(4) = 400*0.9; 
    subplot(121); hold on;
    h1 = pcolor(xdat, ydat, zdat); set(h1,'EdgeColor','None','FaceColor','Flat');
    img_props();
    xlabel(' X ', 'fontweight', 'bold');
    ylabel(' Y ', 'fontweight', 'bold');
    axis([min(xdat(:)), max(xdat(:)), min(ydat(:)), max(ydat(:))]);
    subplot(122); hold on;
    h2 = pcolor(xdat, ydat, zdat_fil); set(h2,'EdgeColor','None','FaceColor','Flat');
    img_props(); cbar_props();
    xlabel(' X ', 'fontweight', 'bold');
    ylabel(' Y ', 'fontweight', 'bold');
    axis([min(xdat(:)), max(xdat(:)), min(ydat(:)), max(ydat(:))]);
end
end