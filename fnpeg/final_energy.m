function [root, isTerminal, direction] = final_energy(~, xSph, eF)
%%%FINAL_ENERGY Event function specifying the final condition for the
%%%integration of the longitudinal dynamics.
%
% Terminates the integration when the current negative energy is equal to
% the final value. All units are dimensionless.
%
% AUTHOR:
%   Davide Amato, CU Boulder, davide.amato@colorado.edu
% 

%% Initializations
isTerminal = 1;
direction  = 0;

% Unpack
r = xSph(1); v = xSph(4);

%% Condition 1: final energy
e = 1/r - 0.5 * v^2;
root = e - eF;

end