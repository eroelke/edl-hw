function rho = atm_exponential(h,rho0,h0,H)
%%ATM_EXPONENTIAL Exponential atmospheric density profile
%
% Density 'rho' depending on altitude 'h' exclusively through the
% exponential atmospheric density profile defined by 'rho0', the density at
% altitude 'h0', and the scale height 'H'.
% Units of 'h', 'h0', 'H', and of 'rho', 'rho0' have to be consistent.
%
% Author:
%   Davide Amato
%   University of Colorado Boulder
%   davide.amato@colorado.edu
%

rho = rho0 .* exp(-(h - h0)./H );

end