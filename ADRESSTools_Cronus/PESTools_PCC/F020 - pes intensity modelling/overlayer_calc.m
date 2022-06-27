function d_o = overlayer_calc(theta, C_m, imfp_m, I_m, C_o, imfp_o, I_o)
% d_o = overlayer_calc(C_m, imfp_m, I_m, C_o, imfp_o, I_o)
%   Function that determines the thickness of a uniform oxide-film on an
%   oxidised metal surface. This is calculated from the relative
%   photoelectron intensity of the oxidic or the metallic component, with
%   respect to that of the corresponding oxide or clean metal, respectively.
%   This model assumes that the underlying metallic layer is an infinitely
%   thick slab, with a uniform oxide layer of thickness d_o on top.
%
%   IN:
%   -   theta:      scalar of the photoelectron take-off angle.
%   -   C_m:        scalar of the volume density of metal atoms in the metal [mole / cc].
%   -   imfp_m:  	scalar of the imfp of electrons in the metal [nm].
%   -   I_m:        scalar of the total peak area of the metal photoelectron peak.
%   -   C_o:        scalar of the volume density of metal atoms in the oxide [mole / cc].
%   -   imfp_o:   	scalar of the imfp of electrons in the oxide [nm].
%   -   I_o:        scalar of the total peak area of the oxide photoelectron peak.
%
%   OUT:
%   -   d_o:      	thickness of the oxide in nm.

%% 1 : Determination of the oxide thickness
d_o = imfp_o .* sin(deg2rad(theta)) .* log(1 + (C_m.*imfp_m.*I_o) ./ (C_o.*imfp_o.*I_m));

end