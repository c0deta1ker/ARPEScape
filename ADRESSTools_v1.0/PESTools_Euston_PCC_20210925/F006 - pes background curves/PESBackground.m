function [roi_xdat, roi_int, roi_bgrnd] = PESBackground(xdat, int, bTYPE, LHS, RHS, ORD, LAM, DEL, BGR)
% [roi_xdat, roi_int, roi_bgrnd] = PESBackground(xdat, int, bTYPE, LHS, RHS, ORD, LAM, DEL, BGR)
%   Function that evaluates the background associated with photoelectron 
%   spectroscopy (PES) data.
%
%   REQ. FUNCTIONS: none
%
%   IN:
%   -   xdat:    	Nx1 column vector of the input domain (binding energy for XPS)
%   -   int:        Nx1 column vector of the intensity range (intensity for XPS)
%   -   bTYPE:    	type of background to use for fitting. Default: "Poly" ("none", "Shir", "LinShir")
%   -   LHS:      	1x2 vector or scalar of the [START, END] or START point of the LHS background
%   -   RHS:    	1x2 vector or scalar of the [START, END] or END point of the RHS background
%   -   ORD:     	positive integer of the Polynomial Order for "Poly" background.
%   -   LAM:        scalar for "LinShir" mixing ratio: lambda = 0 is Pure Shirley; lambda = 1 is Pure Linear.
%   -   DEL:     	scalar for the "LinShir" curve offset in binding energy.
%   -   BGR:     	scalar for a constant background to be included.
%
%   OUT:
%   -   int:     	Nx1 column vector of the output background intensity

%% Default parameters
% Default based on inputs
if nargin < 9; BGR	= 0;  end
if nargin < 8; DEL	= 0;  end
if nargin < 7; LAM	= 0;  end
if nargin < 6; ORD	= 1; end
if nargin < 5; RHS	= mean(xdat(:)) + 2; end
if nargin < 4; LHS	= mean(xdat(:)) - 2; end
if nargin < 3; bTYPE = "Poly"; end
% Default based on empty inputs
if isempty(BGR);	BGR     = 0; end
if isempty(DEL);	DEL     = 0; end
if isempty(LAM);	LAM     = 0; end
if isempty(ORD);	ORD  	= 1; end
if isempty(RHS);    RHS    	= mean(xdat(:)) + 2; end
if isempty(LHS);   	LHS    	= mean(xdat(:)) - 2; end
if isempty(bTYPE);  bTYPE 	= "Poly"; end

%% - 1 - Determination of the photoemission spectrum background
% Ensuring xdat is a column vector
if size(xdat, 2) >1; xdat = xdat'; end
if size(int, 2) >1; int = int'; end
ORD = round(ceil(ORD));    % making sure the order is an integer
% - Extracting the background
if bTYPE == "Poly"
    [roi_xdat, roi_int, roi_bgrnd]      = BgrndPoly(xdat, int, LHS, RHS, ORD);
elseif bTYPE == "Shir"
    [roi_xdat, roi_int, roi_bgrnd]      = BgrndShirley(xdat, int, LHS, RHS);
elseif bTYPE == "LinShir"
    [roi_xdat, roi_int, roi_bgrnd]      = BgrndShirleyOffset(xdat, int, LHS, RHS, LAM, DEL);
elseif bTYPE == "U2Tougaard"
    'UNDER DEVELOPMENT'
elseif bTYPE == "U4Tougaard"
    'UNDER DEVELOPMENT'
elseif bTYPE == "none"
    [~, lhsIndx]	= min(abs(xdat - LHS));
    [~, rhsIndx]	= min(abs(xdat - RHS));
    roi_xdat    = xdat(lhsIndx:rhsIndx);
    roi_int     = int(lhsIndx:rhsIndx);
    roi_bgrnd	= 0.*roi_int;
end
% - Adding a constant background 
roi_bgrnd = roi_bgrnd - BGR;

end