function [roi_xdat, roi_ydat, roi_bgrnd] = PESBackground(xdat, ydat, bTYPE, LHS, RHS, BGR, ARGS, WIN, plot_result)
% [roi_xdat, roi_ydat, roi_bgrnd] = PESBackground(xdat, ydat, bTYPE, LHS, RHS, BGR, ARGS, WIN, plot_result)
%   Function that evaluates the background associated with 1D data or
%   photoelectron spectroscopy (PES) data acquired at the ADRESS beamline, 
%   SLS. The background type can be defined over a given region of interest
%   (ROI). 
%
%   REQ. FUNCTIONS: none
%
%   IN:
%   -   xdat:           N×1 column vector of the input domain (binding energy for PES).
%   -   ydat:           N×1 column vector of the intensity range (intensity for PES).
%   -   bTYPE:          string of the type of background to use for fitting. {Default: "Poly" ("none", "Poly", "Shir", "LinShir", "StepFDDGpL", "StepFDDGsL")}
%   -   LHS:            scalar of the LHS x-axis START position of the ROI. {Default: 25th percentile.}
%   -   RHS:            scalar of the RHS x-axis END position of the ROI. {Default: 75th percentile.}
%   -   BGR:            scalar for a constant y-axis background to be included. {Default: 0}
%   -   ARGS:           cell vector containing the arguments for the background type used;
%                           "none":         1×1 {}
%                           "Poly":         1×1 {polyOrder}
%                           "Shir":         1×1 {}
%                           "LinShir":      1×2 {lambda, delta}
%                           "StepFDDGpL":   1×3 {fdd_ef, fdd_T, fdd_fwhm}
%                           "StepFDDGsL":   1×3 {fdd_ef, fdd_T, fdd_fwhm}
%   -   WIN:            scalar of the window size around the LHS and RHS positions. {Default: 1% of ROI size.}
%   -   plot_result:    if 1, will plot figure summary, otherwise it wont.
%
%   OUT:
%   -   roi_xdat:       M×1 column vector of the ROI x-domain (binding energy for PES).
%   -   roi_ydat:       M×1 column vector of the ROI y-values (intensity for PES).
%   -   roi_bgrnd:      M×1 column vector of the best fit polynomial background to the ROI data.

%% Default parameters
% - Setting the default conditions
if nargin < 3;      bTYPE 	= "Poly"; end
if nargin < 4;    	LHS     = mean(xdat(:)) - abs(0.25*range(xdat(:))); end
if nargin < 5;   	RHS     = mean(xdat(:)) + abs(0.25*range(xdat(:))); end
if nargin < 6;   	BGR     = 0.0; end
if nargin < 7;   	ARGS    = {1}; end
if nargin < 8;    	WIN = abs(0.01*range(xdat(:))); end
if nargin < 9;    	plot_result = 0; end
if isempty(bTYPE);  bTYPE 	= "Poly"; end
if isempty(LHS);   	LHS    	= mean(xdat(:)) - abs(0.25*range(xdat(:))); end
if isempty(RHS);    RHS    	= mean(xdat(:)) + abs(0.25*range(xdat(:))); end
if isempty(BGR);	BGR     = 0; end
if isempty(ARGS);	ARGS    = {1}; end
if isempty(plot_result); plot_result = 0; end
if isempty(WIN);	WIN    = abs(0.01*range(xdat(:))); end

%% Ensuring the data is in the form of 1D column vectors
if size(xdat, 2) > 1; xdat = xdat'; end
if size(ydat, 2) > 1; ydat = ydat'; end

%% - 1 - Determination of the PES background
% - EXTRACTING POLYNOMIAL BACKGROUND
if strcmpi(bTYPE,"Poly")
    if length(ARGS) == 1;   [roi_xdat, roi_ydat, roi_bgrnd]     = BgrndPoly(xdat, ydat, ARGS{1}, LHS, RHS, WIN, WIN, plot_result);
    else;                   [roi_xdat, roi_ydat, roi_bgrnd]     = BgrndPoly(xdat, ydat, [], LHS, RHS, WIN, WIN, plot_result);
    end
% - EXTRACTING SHIRLEY BACKGROUND 
elseif strcmpi(bTYPE,"Shir")
    [roi_xdat, roi_ydat, roi_bgrnd]     = BgrndShirley(xdat, ydat, LHS, RHS, [], [], plot_result);
% - EXTRACTING LINEAR-SHIRLEY BACKGROUND
elseif strcmpi(bTYPE,"LinShir")
    if length(ARGS) == 1;       [roi_xdat, roi_ydat, roi_bgrnd]     = BgrndShirleyOffset(xdat, ydat, ARGS{1}, [], LHS, RHS, WIN, WIN, plot_result);
    elseif length(ARGS) == 2;   [roi_xdat, roi_ydat, roi_bgrnd]     = BgrndShirleyOffset(xdat, ydat, ARGS{1}, ARGS{2}, LHS, RHS, WIN, WIN, plot_result);
    else;                       [roi_xdat, roi_ydat, roi_bgrnd]     = BgrndShirleyOffset(xdat, ydat, [], [], LHS, RHS, WIN, WIN, plot_result);
    end
% - EXTRACTING STEP-FDD-PRODUCT BACKGROUND
elseif strcmpi(bTYPE,"StepFDDGpL")
    if length(ARGS) == 1;       [roi_xdat, roi_ydat, roi_bgrnd]     = BgrndStepFDDGpL(xdat, ydat, ARGS{1}, [], [], LHS, RHS, WIN, WIN, plot_result);
    elseif length(ARGS) == 2;   [roi_xdat, roi_ydat, roi_bgrnd]     = BgrndStepFDDGpL(xdat, ydat, ARGS{1}, ARGS{2}, [], LHS, RHS, WIN, WIN, plot_result);
    elseif length(ARGS) == 3;   [roi_xdat, roi_ydat, roi_bgrnd]     = BgrndStepFDDGpL(xdat, ydat, ARGS{1}, ARGS{2}, ARGS{3}, LHS, RHS, WIN, WIN, plot_result);
    else;                       [roi_xdat, roi_ydat, roi_bgrnd]     = BgrndStepFDDGpL(xdat, ydat, [], [], [], LHS, RHS, WIN, WIN, plot_result);
    end
% - EXTRACTING STEP-FDD-SUM BACKGROUND
elseif strcmpi(bTYPE,"StepFDDGsL")
    if length(ARGS) == 1;       [roi_xdat, roi_ydat, roi_bgrnd]     = BgrndStepFDDGsL(xdat, ydat, ARGS{1}, [], [], LHS, RHS, WIN, WIN, plot_result);
    elseif length(ARGS) == 2;   [roi_xdat, roi_ydat, roi_bgrnd]     = BgrndStepFDDGsL(xdat, ydat, ARGS{1}, ARGS{2}, [], LHS, RHS, WIN, WIN, plot_result);
    elseif length(ARGS) == 3;   [roi_xdat, roi_ydat, roi_bgrnd]     = BgrndStepFDDGsL(xdat, ydat, ARGS{1}, ARGS{2}, ARGS{3}, LHS, RHS, WIN, WIN, plot_result);
    else;                       [roi_xdat, roi_ydat, roi_bgrnd]     = BgrndStepFDDGsL(xdat, ydat, [], [], [], LHS, RHS, WIN, WIN, plot_result);
    end
% - EXTRACTING NO BACKGROUND
elseif strcmpi(bTYPE,"none")
    [~, lhsIndx]	= min(abs(xdat - LHS));
    [~, rhsIndx]	= min(abs(xdat - RHS));
    roi_xdat        = xdat(lhsIndx:rhsIndx);
    roi_ydat        = ydat(lhsIndx:rhsIndx);
    roi_bgrnd       = 0.*roi_ydat;
end
% - Adding a constant background 
roi_bgrnd = roi_bgrnd - BGR;

end