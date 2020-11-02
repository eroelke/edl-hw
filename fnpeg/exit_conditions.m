function [roots, isTerminal, direction] = exit_conditions(t, x, lRef, tRef, hf, ef, mDry, planet, integ)
%FINAL_FNPEG Event function specifying final conditions for the
%integration of the FNPEG system of equations.
%%TODO
%
% The following events terminate the integration:
% roots(1) = 0: reached final energy.
% roots(2) = 0: collision with the planetary surface.
% roots(3) = 0: reached final time.
%
% Author: Davide Amato, CU Boulder, davide.amato@colorado.edu
% 
% Last change: 11 July 2019.
%
%% Initialization
roots = ones(4,1);
isTerminal = ones(4,1);
direction  = zeros(4,1);
direction(4) = -1;

% Unpack, auxiliary variables
rMag = norm(x(1:3));
vMag = norm(x(4:6));
mass = x(7);

%% Condition 1: final altitude
roots(1) = (rMag - planet.r / lRef) - hf;

%% Condition 2: final energy
e = -rv2energy(1, rMag, vMag);
roots(2) = ef - e;

%% Condition 3: final time
roots(3) = integ.tF / tRef - t;

%% Condition 4: negative mass
roots(4) = mass - mDry;

end