function fitStr = xps_solver(xpsStr, cTYPE, solve_type)
% fitStr = xps_solver(xpsStr, cTYPE, iPESCurves, dPESCurves, bTYPE, iPESBgrnd, dPESBgrnd, solve_type)
%   Function that runs a a Global Optimisation Algorithm based on Simulated
%   Annealing to determine the best curve fit to the input XPS data and
%   model (initial conditions and bounds). The temperature of the Simulated
%   Annealing is set to walk a large span of the parameter space to ensure
%   the global minimum is obtained. The minimisation function is based on
%   minimising the standard deviation of the sum of the squared residuals.
%
%   REQ. FUNCTIONS: none
%
%   IN:
%   -   xpsStr:         MATLAB data-structure that contains the XPS data.
%   -   cTYPE:          1xN vector of the type of curve to use for the nth state.
%   -   iPESCurves:   	Nx8 array: the n'th peak parameters [BE,INT,FWHM,MR,LSE,LSI,LSW,ASY]
%   -   dPESCurves: 	Nx8 array: the n'th peak parameters uncertainties [BE,INT,FWHM,MR,LSE,LSI,LSW,ASY]
%   -   bTYPE:          string of the type of background to use for fitting. Default: "Poly" ("none", "Shir", "LinShir")
%   -   iPESBgrnd:      1x5 vector of the background parameters: [LHS,RHS,ORD,LAM,DEL,BGR]
%   -   dPESBgrnd:     	1x5 vector of the background parameter uncertainties: [LHS,RHS,ORD,LAM,DEL,BGR]
%   -   solve_type:     String of either "lsqnonlin", "simulannealbnd"
%
%   OUT:
%   -   fitStr:         MATLAB data-structure that contains all the fit parameters / variables / information


