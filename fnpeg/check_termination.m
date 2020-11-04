%FINAL_FNPEG Event function specifying final conditions for the
%integration of the FNPEG system of equations.
% 
function [roots, isTerminal, direction] = check_termination(t, x, lRef, tRef, hf, ef, planet, integ)
roots = ones(3,1);
isTerminal = ones(3,1);
direction  = zeros(3,1);

rMag = norm(x(1:3));
vMag = norm(x(4:6));

% Condition 1: final altitude
roots(1) = (rMag - planet.r / lRef) - hf;

% Condition 2: final energy
e = -rv2xi(1, rMag, vMag);
roots(2) = ef - e;

% Condition 3: final time
roots(3) = integ.tF / tRef - t;
end