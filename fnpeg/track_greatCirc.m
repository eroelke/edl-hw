%% 
% Given initial lat/lon,
% the great circle angle "s", and the azimuth "psi" (measured clockwise
% from the North), return the latitude and longitude of the final point
% "phi", "lambda".
% Inputs:
%   long0: initial longitude, deg
%   lat0: initial latitude, deg
%   s: great circle angle, deg (aka range)
%   psi: azimuth, deg
% Outputs:
%   lambda: longitude of final point
%   phi:    latitude of final point
% 
function [lambda, phi] = track_greatCirc(long0, lat0, s, az)
% put vars on the stack
sphi0 = sind(lat0);
cphi0 = cosd(lat0);
ss = sind(s);
cs = cosd(s);
spsi = sind(az);
cpsi = cosd(az);

% compute return vals
lambda = long0 + atan2d(ss*spsi, cphi0*cs - sphi0*ss*cpsi); %longitude
phi = 90 - acosd(sphi0*cs + cphi0*ss*cpsi); %latitude

end