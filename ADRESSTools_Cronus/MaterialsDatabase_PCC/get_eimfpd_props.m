function [ke_dat, imfp] = get_eimfpd_props(element)
% [ke_dat, imfp] = get_eimfpd_props(element)
%   This is a function that extracts all of the experimental IMFPs from the
%   Electron Inelastic Mean Free Path Database (eIMFPD) for the element 
%   defined here by the user. This database contains all the experimental
%   IMFP values for elements 1 - 92, for kinetic energies in the approximate
%   range of 10 - 10,000 eV.
%
%   IN:
%   -   element:	char/string of the element; e.g. "H", "He", "Si", "In"...
%
%   OUT:
%   -   ke_dat:     N×M array of the N'th electron kinetic energies taken from M'th file.
%   -   imfp:   	N×M array of the N'th imfp taken from M'th file.

%% Default parameters
if nargin < 1; element = []; end
if isempty(element); element = []; end

%% - 1 - Loading the eIMFP properties from our MATLAB database
eIMFPD_PCC = load('eIMFPD_PCC.mat'); eIMFPD_PCC = eIMFPD_PCC.eIMFPD_PCC;

%% - 2 - Find the index of the element defined by the user
% -- Find the element
if ~isempty(element)
    % -- Find the row-number / index for the element of choice (case-insensitive)
    element  	= char(element);
    ele_indx 	= find(strcmpi([eIMFPD_PCC.ATOM_SYMB], element));
    if isempty(ele_indx)
        msg = 'Element could not be identified. Only use atomic-symbols for elements 1 - 92; H, He, Li, Be..., Pa, U'; error(msg)
    end
% -- If no element is defined, return an error
else
    msg = 'Element could not be identified. Only use atomic-symbols for elements 1 - 92; H, He, Li, Be..., Pa, U'; error(msg)
end

%% - 3 - Extracting the IMFP data
data  	= eIMFPD_PCC.DataTable{ele_indx};
ke_dat 	= [];
imfp    = [];
for i = 1:eIMFPD_PCC.Length(ele_indx)
    ke_dat(:,i)	= table2array(data(:,2*i-1));
    imfp(:,i)	= table2array(data(:,2*i));
end
end